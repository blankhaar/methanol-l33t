"""Build the separate A/E :sup:`13`CH3OH ``vt=0,1,2`` LAMDA files.

The spectroscopic model is isotope-specific: Xu et al. (2014) Table 2 is
loaded into the vendored BELGI Hamiltonian and the quoted constant permanent
dipoles are evaluated with BELGI ``INTEN``.  Collision dynamics are explicitly
labelled proxies obtained by exact quantum-identity mapping of the normalized
main-isotopologue calculations followed by the documented reduced-mass
prefactor.  Nothing in this module performs a fit.
"""

from __future__ import annotations

from collections import Counter
import csv
from dataclasses import dataclass, replace
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Iterable

from ...build import (
    MasterCollision,
    MasterLevel,
    SelectedCollision,
    _balanced_levels,
    _calculated_collisions,
    _collision_coverage,
    _write_sidecars,
)
from ...collision_integrity import (
    assert_collision_selection,
    audit_basecol_collision_master,
)
from ...hamiltonian import build_frozen_levels, sha256_file
from ...lamda import read_lamda, validate_lamda, write_lamda
from ...models import (
    CollisionPartner,
    CollisionTransition,
    LamdaData,
    Level,
    RadiativeTransition,
)
from ...paths import portable_path
from ...quantum import LevelKey
from ...radiative_integrity import audit_radiative_integrity
from ...strengths import build_frozen_strengths
from .cdms import (
    CDMS_DATASET_ID,
    catalog_file_to_overrides,
    rescale_einstein_a_nu3,
)
from .collisions import (
    ISOTOPOLOGUE_PROXY_STATUS,
    scale_rate_tuple,
    scaling_factor,
    scaling_metadata,
)
from .profile import (
    DIPOLES_DEBYE,
    HAMILTONIAN_DATASET_ID,
    LAMDA_FILENAMES,
    MOLECULAR_MASS_U,
    RADIATIVE_DATASET_ID,
    SPECIES,
    SPECIES_NAMES,
    load_assignment_sets,
)


SUBPROJECT_RELATIVE_PATH = Path("isotopologues/13ch3oh")
DEFAULT_CONFIG_RELATIVE_PATH = SUBPROJECT_RELATIVE_PATH / "config/build.json"
DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH = (
    SUBPROJECT_RELATIVE_PATH / "config/observed-frequencies.json"
)
MAIN_COLLISION_MASTER_RELATIVE_PATH = Path(
    "data/normalized/master/collisions.csv"
)
LEVEL_OUTPUT_RELATIVE_PATH = Path(
    "data/normalized/belgi-xu2014-frozen"
)
STRENGTH_OUTPUT_RELATIVE_PATH = Path(
    "data/normalized/belgi-xu2014-constant-dipole"
)
STRENGTH_DIPOLE_MODEL_ID = (
    "Xu2014_13CH3OH_constant_mu_a_mu_b_cos3gamma_zero"
)
MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE = 1.0e-2


class CH3OH13BuildError(RuntimeError):
    """Raised when a product cannot be built without guessing."""


@dataclass(frozen=True, slots=True)
class FrequencyOverride:
    frequency_mhz: float
    uncertainty_mhz: float | None
    source: str
    candidate_count: int
    source_row: int | None


def subproject_root(project_root: str | Path) -> Path:
    return Path(project_root).resolve() / SUBPROJECT_RELATIVE_PATH


def _load_config(root: Path) -> dict[str, object]:
    path = root / DEFAULT_CONFIG_RELATIVE_PATH
    if not path.is_file():
        raise CH3OH13BuildError(f"missing isotope build configuration: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("torsional_states") != [0, 1, 2]:
        raise CH3OH13BuildError("13CH3OH build requires torsional states [0, 1, 2]")
    try:
        configured_mass = float(payload.get("molecular_mass_u", math.nan))
    except (TypeError, ValueError) as exc:
        raise CH3OH13BuildError("configured molecular mass must be numeric") from exc
    if not math.isclose(
        configured_mass, MOLECULAR_MASS_U, rel_tol=0.0, abs_tol=1.0e-12
    ):
        raise CH3OH13BuildError(
            f"configured mass {configured_mass} does not match {MOLECULAR_MASS_U} u"
        )
    configured_dipoles = payload.get("dipoles_debye")
    try:
        dipoles_match = (
            isinstance(configured_dipoles, dict)
            and set(configured_dipoles) == set(DIPOLES_DEBYE)
            and all(
                math.isclose(
                    float(configured_dipoles[name]),
                    value,
                    rel_tol=0.0,
                    abs_tol=1.0e-15,
                )
                for name, value in DIPOLES_DEBYE.items()
            )
        )
    except (TypeError, ValueError):
        dipoles_match = False
    if not dipoles_match:
        raise CH3OH13BuildError("configured dipoles do not match the Xu 2014 prescription")

    j_max = payload.get("j_max")
    abs_k_max = payload.get("abs_k_max")
    if (
        isinstance(j_max, bool)
        or not isinstance(j_max, int)
        or isinstance(abs_k_max, bool)
        or not isinstance(abs_k_max, int)
        or not (1 <= j_max <= 22)
        or not (0 <= abs_k_max <= min(j_max, 14))
    ):
        raise CH3OH13BuildError(
            "13CH3OH strength scope requires 1 <= j_max <= 22 and "
            "0 <= abs_k_max <= min(j_max, 14)"
        )
    try:
        ground_max_cm1 = float(payload["ground_manifold_max_energy_cm1"])
        effective_dipole = float(
            payload["minimum_effective_transition_dipole_debye"]
        )
        generation_smu2 = float(payload["strength_generation_smu2_min_debye2"])
        temperatures = tuple(
            float(value) for value in payload["collision_temperature_grid_k"]
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise CH3OH13BuildError(
            "13CH3OH numeric build policy is missing or malformed"
        ) from exc
    if not math.isfinite(ground_max_cm1) or ground_max_cm1 <= 0.0:
        raise CH3OH13BuildError(
            "ground-manifold maximum energy must be finite and positive"
        )
    if not math.isclose(
        effective_dipole,
        MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise CH3OH13BuildError("13CH3OH products require the declared 0.01-D cutoff")
    if not math.isfinite(generation_smu2) or not math.isclose(
        generation_smu2,
        effective_dipole**2,
        rel_tol=0.0,
        abs_tol=1.0e-18,
    ):
        raise CH3OH13BuildError(
            "strength generation floor must equal the squared effective-dipole cutoff"
        )
    if (
        not temperatures
        or any(not math.isfinite(value) or value <= 0.0 for value in temperatures)
        or any(right <= left for left, right in zip(temperatures, temperatures[1:]))
    ):
        raise CH3OH13BuildError(
            "collision temperature grid must be positive and increase strictly"
        )
    if payload.get("calculated_rate_extrapolation") != "forbidden":
        raise CH3OH13BuildError("calculated collision-rate extrapolation must be forbidden")
    return payload


def _optional_float(value: str | None) -> float | None:
    if value is None or not value.strip():
        return None
    return float(value)


def _split_floats(value: str) -> tuple[float, ...]:
    return tuple(float(token) for token in value.split(";") if token)


def _load_main_collision_master(root: Path) -> list[MasterCollision]:
    """Load only the collision master, avoiding the unrelated main line list."""
    path = root / MAIN_COLLISION_MASTER_RELATIVE_PATH
    if not path.is_file():
        raise CH3OH13BuildError(
            f"missing normalized main-isotopologue collision master: {path}; "
            "run `python -m methanol_lamda normalize` first"
        )
    result: list[MasterCollision] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            result.append(
                MasterCollision(
                    species=row["species"],
                    collider_id=int(row["collider_id"]),
                    upper_key=LevelKey.decode(row["upper_key"]),
                    lower_key=LevelKey.decode(row["lower_key"]),
                    temperatures=_split_floats(row["temperatures_k"]),
                    rates=_split_floats(row["rates_cm3_s"]),
                    source=row["source"],
                    status=row["status"],
                    raw_sha256=row.get("raw_sha256", ""),
                    source_upper_id=row.get("source_upper_id", ""),
                    source_lower_id=row.get("source_lower_id", ""),
                    source_upper_energy_cm1=_optional_float(
                        row.get("source_upper_energy_cm1")
                    ),
                    source_lower_energy_cm1=_optional_float(
                        row.get("source_lower_energy_cm1")
                    ),
                    mapping_confidence=row.get("mapping_confidence", ""),
                )
            )
    if not result:
        raise CH3OH13BuildError("main-isotopologue collision master is empty")
    return result


def _load_frozen_levels(path: Path, constants_sha256: str) -> list[MasterLevel]:
    result: list[MasterLevel] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            result.append(
                MasterLevel(
                    key=LevelKey.decode(row["level_key"]),
                    energy_cm1=float(row["energy_cm1"]),
                    weight=float(row["weight"]),
                    source=row["source"],
                    source_index=row["source_index"],
                    mapping_confidence=row["mapping_confidence"],
                    raw_sha256=constants_sha256,
                )
            )
    if not result:
        raise CH3OH13BuildError(f"frozen level table is empty: {path}")
    if any(level.source != HAMILTONIAN_DATASET_ID for level in result):
        raise CH3OH13BuildError("frozen level table has an unexpected dataset identity")
    return result


def _configured_observed_overrides(
    root: Path,
) -> tuple[dict[frozenset[LevelKey], FrequencyOverride], dict[str, object]]:
    """Load optional CDMS laboratory rows and explicit reviewed overrides."""
    path = root / DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH
    payload = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    result: dict[frozenset[LevelKey], FrequencyOverride] = {}
    catalog_relative = Path(
        str(payload.get("catalog", "data/raw/official/cdms-033502/c033502.cat"))
    )
    catalog_path = subproject_root(root) / catalog_relative
    catalog_summary: dict[str, object] = {
        "path": portable_path(root, catalog_path),
        "present": catalog_path.is_file(),
        "dataset_id": CDMS_DATASET_ID,
        "selected_experimental_rows": 0,
    }
    if catalog_path.is_file():
        selected = catalog_file_to_overrides(catalog_path)
        for override in selected.values():
            result[override.identity.unordered_endpoints] = FrequencyOverride(
                frequency_mhz=override.frequency_mhz,
                uncertainty_mhz=override.uncertainty_mhz,
                source=override.source_id,
                candidate_count=override.candidate_count,
                source_row=override.source_row,
            )
        catalog_summary.update(
            {
                "sha256": sha256_file(catalog_path),
                "selected_experimental_rows": len(selected),
                "selection": "negative tag -033502 and exact endpoint identity",
            }
        )

    # Explicit reviewed records have final precedence.  They use the same
    # exact endpoint schema as the generated CDMS crosswalk.
    manual_count = 0
    for row in payload.get("overrides", []):
        if not isinstance(row, dict):
            raise CH3OH13BuildError("observed-frequency override must be an object")
        upper = LevelKey.decode(str(row["upper_key"]))
        lower = LevelKey.decode(str(row["lower_key"]))
        if upper.species != lower.species:
            raise CH3OH13BuildError("observed-frequency override changes A/E symmetry")
        frequency = float(row["frequency_mhz"])
        if not math.isfinite(frequency) or frequency <= 0.0:
            raise CH3OH13BuildError("observed-frequency override must be positive")
        result[frozenset((upper, lower))] = FrequencyOverride(
            frequency_mhz=frequency,
            uncertainty_mhz=(
                None
                if row.get("uncertainty_mhz") is None
                else float(row["uncertainty_mhz"])
            ),
            source=str(row["source"]),
            candidate_count=int(row.get("candidate_count", 1)),
            source_row=(
                None if row.get("source_row") is None else int(row["source_row"])
            ),
        )
        manual_count += 1
    catalog_summary["explicit_reviewed_override_count"] = manual_count
    catalog_summary["available_override_count"] = len(result)
    return result, catalog_summary


def _select_radiative(
    path: Path,
    species: str,
    selected_levels: list[MasterLevel],
    key_to_index: dict[LevelKey, int],
    overrides: dict[frozenset[LevelKey], FrequencyOverride],
    *,
    minimum_smu2_debye2: float,
    constants_sha256: str,
) -> tuple[list[RadiativeTransition], list[dict[str, object]], dict[str, object]]:
    if not math.isfinite(minimum_smu2_debye2) or minimum_smu2_debye2 <= 0.0:
        raise CH3OH13BuildError(
            "minimum radiative S mu^2 must be finite and positive"
        )
    energies = {level.key: level.energy_cm1 for level in selected_levels}
    candidates: list[tuple[LevelKey, LevelKey, dict[str, str]]] = []
    candidate_pairs: set[frozenset[LevelKey]] = set()
    rejected = 0
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["species"] != species:
                continue
            upper = LevelKey.decode(row["upper_key"])
            lower = LevelKey.decode(row["lower_key"])
            if upper not in key_to_index or lower not in key_to_index:
                continue
            smu2 = float(row["smu2_debye2"])
            if not math.isfinite(smu2) or smu2 <= 0.0:
                raise CH3OH13BuildError(
                    "BELGI radiative row contains non-positive or non-finite S mu^2"
                )
            if smu2 < minimum_smu2_debye2:
                rejected += 1
                continue
            pair = frozenset((upper, lower))
            if pair in candidate_pairs:
                raise CH3OH13BuildError(
                    f"duplicate BELGI radiative endpoint pair {upper.encode()} / {lower.encode()}"
                )
            candidate_pairs.add(pair)
            candidates.append((upper, lower, row))

    transitions: list[RadiativeTransition] = []
    provenance: list[dict[str, object]] = []
    used_override_pairs: set[frozenset[LevelKey]] = set()
    selected_keys = frozenset(key_to_index)
    species_override_pairs = {
        pair
        for pair in overrides
        if len(pair) == 2
        and pair <= selected_keys
        and all(key.species == species for key in pair)
    }
    for upper, lower, row in candidates:
        if energies[upper] <= energies[lower]:
            raise CH3OH13BuildError(
                "BELGI radiative upper endpoint is not above the selected lower endpoint: "
                f"{upper.encode()} -> {lower.encode()}"
            )
        calculated_frequency_mhz = float(row["frequency_mhz"])
        calculated_a_s = float(row["a_s"])
        smu2 = float(row["smu2_debye2"])
        if (
            not math.isfinite(calculated_frequency_mhz)
            or not math.isfinite(calculated_a_s)
            or calculated_frequency_mhz <= 0.0
            or calculated_a_s <= 0.0
        ):
            raise CH3OH13BuildError("BELGI radiative row contains a non-positive value")
        pair = frozenset((upper, lower))
        override = overrides.get(pair)
        if override is None:
            adopted_frequency_mhz = calculated_frequency_mhz
            adopted_a_s = calculated_a_s
            frequency_source = HAMILTONIAN_DATASET_ID
            status = "calculated_xu2014_table2_frequency"
            uncertainty_mhz = None
            candidate_count = 0
            source_row = None
        else:
            adopted_frequency_mhz = override.frequency_mhz
            if (
                not math.isfinite(adopted_frequency_mhz)
                or adopted_frequency_mhz <= 0.0
            ):
                raise CH3OH13BuildError(
                    "observed-frequency override must be finite and positive"
                )
            adopted_a_s = rescale_einstein_a_nu3(
                calculated_a_s,
                calculated_frequency_mhz,
                adopted_frequency_mhz,
            )
            frequency_source = override.source
            status = "calculated_strength_experimental_cdms_frequency"
            uncertainty_mhz = override.uncertainty_mhz
            candidate_count = override.candidate_count
            source_row = override.source_row
            used_override_pairs.add(pair)
        transitions.append(
            RadiativeTransition(
                0,
                key_to_index[upper],
                key_to_index[lower],
                adopted_a_s,
                adopted_frequency_mhz / 1000.0,
            )
        )
        provenance.append(
            {
                "upper_key": upper.encode(),
                "lower_key": lower.encode(),
                "a_s": adopted_a_s,
                "frequency_ghz": adopted_frequency_mhz / 1000.0,
                "smu2_debye2": smu2,
                "effective_transition_dipole_debye": math.sqrt(smu2),
                "a_source": RADIATIVE_DATASET_ID,
                "frequency_source": frequency_source,
                "status": status,
                "raw_sha256": constants_sha256,
                "source_upper_id": upper.encode(),
                "source_lower_id": lower.encode(),
                "mapping_confidence": "exact_quantum_identity",
                "duplicate_candidates": 1,
                "calculated_frequency_mhz": calculated_frequency_mhz,
                "adopted_frequency_mhz": adopted_frequency_mhz,
                "frequency_uncertainty_mhz": uncertainty_mhz,
                "frequency_candidate_count": candidate_count,
                "frequency_source_row": source_row,
            }
        )

    order = sorted(
        range(len(transitions)),
        key=lambda index: (transitions[index].upper, transitions[index].lower),
    )
    transitions = [
        replace(transitions[index], index=number)
        for number, index in enumerate(order, start=1)
    ]
    provenance = [
        {**provenance[index], "transition_index": number}
        for number, index in enumerate(order, start=1)
    ]
    summary = {
        "retained_transition_count": len(transitions),
        "rejected_below_cutoff": rejected,
        "minimum_smu2_debye2": minimum_smu2_debye2,
        "minimum_effective_transition_dipole_debye": math.sqrt(
            minimum_smu2_debye2
        ),
        "available_observed_override_count": len(species_override_pairs),
        "matched_observed_override_count": len(used_override_pairs),
        "unmatched_observed_override_count": len(
            species_override_pairs - used_override_pairs
        ),
    }
    return transitions, provenance, summary


def _scale_calculated_collisions(
    rows: Iterable[SelectedCollision],
) -> list[SelectedCollision]:
    scaled: list[SelectedCollision] = []
    for row in rows:
        factor = scaling_factor(row.collider_id)
        source_status = row.status
        scaled.append(
            replace(
                row,
                rates=scale_rate_tuple(row.rates, row.collider_id),
                status=ISOTOPOLOGUE_PROXY_STATUS,
                multiplier=row.multiplier * factor,
                transform=(
                    f"{row.transform}; source_status={source_status}; "
                    "k13=sqrt(mu12/mu13)*k12 at same temperature"
                ),
                # Mapping confidence describes endpoint identity only.  The
                # dynamics approximation belongs in status and transform.
                mapping_confidence=row.mapping_confidence,
            )
        )
    return scaled


def _partners_13ch3oh(
    rows: list[SelectedCollision], target_temperatures: tuple[float, ...]
) -> list[CollisionPartner]:
    descriptions = {
        2: "13CH3OH-pH2 [reduced-mass-scaled 12CH3OH proxy]",
        3: "13CH3OH-oH2 [reduced-mass-scaled 12CH3OH proxy]",
        6: "13CH3OH-He [reduced-mass-scaled 12CH3OH proxy]",
    }
    partners: list[CollisionPartner] = []
    for collider_id in (2, 3, 6):
        chosen = sorted(
            (row for row in rows if row.collider_id == collider_id),
            key=lambda row: (row.upper, row.lower),
        )
        if not chosen:
            continue
        partners.append(
            CollisionPartner(
                collider_id,
                descriptions[collider_id],
                target_temperatures,
                [
                    CollisionTransition(index, row.upper, row.lower, row.rates)
                    for index, row in enumerate(chosen, start=1)
                ],
            )
        )
    return partners


def _write_isotope_sidecars(
    stage_root: Path,
    filename: str,
    radiative: list[dict[str, object]],
    collisions: list[SelectedCollision],
    target_temperatures: tuple[float, ...],
) -> None:
    """Write standard collision provenance plus the extended frequency audit."""
    standard_radiative_fields = (
        "transition_index",
        "upper_key",
        "lower_key",
        "a_s",
        "frequency_ghz",
        "smu2_debye2",
        "effective_transition_dipole_debye",
        "a_source",
        "frequency_source",
        "status",
        "raw_sha256",
        "source_upper_id",
        "source_lower_id",
        "mapping_confidence",
        "duplicate_candidates",
    )
    # The shared writer provides the collision contract.  Give its radiative
    # writer only the fields it declares, then replace that smaller file with
    # the isotope-specific superset below.
    core_rows = [
        {field: row[field] for field in standard_radiative_fields}
        for row in radiative
    ]
    _write_sidecars(
        stage_root,
        filename,
        core_rows,
        collisions,
        target_temperatures,
    )
    path = stage_root / "outputs/provenance" / f"{filename}.radiative.csv"
    extra_fields = (
        "calculated_frequency_mhz",
        "adopted_frequency_mhz",
        "frequency_uncertainty_mhz",
        "frequency_candidate_count",
        "frequency_source_row",
    )
    fields = ("file", *standard_radiative_fields, *extra_fields)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in radiative:
            writer.writerow({"file": filename, **row})


def _write_level_sidecar(
    stage_root: Path,
    filename: str,
    levels: list[MasterLevel],
) -> Path:
    path = stage_root / "outputs/provenance" / f"{filename}.levels.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "file",
        "level_index",
        "level_key",
        "energy_cm1",
        "weight",
        "source",
        "source_index",
        "mapping_confidence",
        "raw_sha256",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for index, level in enumerate(levels, start=1):
            writer.writerow(
                {
                    "file": filename,
                    "level_index": index,
                    "level_key": level.key.encode(),
                    "energy_cm1": f"{level.energy_cm1:.17g}",
                    "weight": f"{level.weight:.17g}",
                    "source": level.source,
                    "source_index": level.source_index,
                    "mapping_confidence": level.mapping_confidence,
                    "raw_sha256": level.raw_sha256,
                }
            )
    return path


def _radiative_components(level_count: int, rows: list[RadiativeTransition]) -> int:
    parent = list(range(level_count + 1))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    for row in rows:
        left, right = find(row.upper), find(row.lower)
        if left != right:
            parent[right] = left
    return len({find(index) for index in range(1, level_count + 1)})


def build_products(
    project_root: str | Path,
    *,
    compiler: str | Path | None = None,
) -> dict[str, object]:
    """Generate the two isotope-specific proxy LAMDA products and sidecars."""
    root = Path(project_root).resolve()
    isotope_root = subproject_root(root)
    collision_source_audit = audit_basecol_collision_master(root)
    config = _load_config(root)
    assignments = load_assignment_sets(root)
    constants_path = assignments.constants_path
    constants_sha256 = sha256_file(constants_path)
    j_max = int(config["j_max"])
    abs_k_max = int(config["abs_k_max"])
    vt_max = 2
    ground_max_cm1 = float(config["ground_manifold_max_energy_cm1"])
    minimum_smu2 = float(config["strength_generation_smu2_min_debye2"])
    declared_minimum_smu2 = (
        float(config["minimum_effective_transition_dipole_debye"]) ** 2
    )
    if not math.isclose(
        minimum_smu2,
        declared_minimum_smu2,
        rel_tol=0.0,
        abs_tol=1.0e-18,
    ):
        raise CH3OH13BuildError(
            "strength generation floor must equal the squared effective-dipole cutoff"
        )
    target_temperatures = tuple(
        float(value) for value in config["collision_temperature_grid_k"]
    )
    if not target_temperatures or any(
        right <= left
        for left, right in zip(target_temperatures, target_temperatures[1:])
    ):
        raise CH3OH13BuildError("collision temperature grid must increase strictly")

    level_dir = isotope_root / LEVEL_OUTPUT_RELATIVE_PATH
    strength_dir = isotope_root / STRENGTH_OUTPUT_RELATIVE_PATH
    level_manifest = build_frozen_levels(
        root,
        constants_path,
        root / "vendor/belgi",
        level_dir,
        j_max=j_max,
        abs_k_max=abs_k_max,
        vt_max=vt_max,
        compiler=compiler,
        require_xu_final=False,
        dataset_id=HAMILTONIAN_DATASET_ID,
    )
    strength_manifest = build_frozen_strengths(
        root,
        constants_path,
        root / "vendor/belgi",
        strength_dir,
        j_max=j_max,
        abs_k_max=abs_k_max,
        vt_max=vt_max,
        smu2_min_debye2=minimum_smu2,
        compiler=compiler,
        require_xu_final=False,
        dataset_id=RADIATIVE_DATASET_ID,
        reference_levels_path=level_dir / "levels.csv",
        dipole_model_name="Xu et al. (2014) 13CH3OH constant permanent dipoles",
        dipole_model_id=STRENGTH_DIPOLE_MODEL_ID,
        dipole_provenance_path=(
            isotope_root / "data/raw/user-supplied/xu-et-al-2014.pdf"
        ),
        completeness_claim=(
            "selection-complete above the declared threshold for the Xu 2014 "
            "13CH3OH constant-permanent-dipole operator"
        ),
    )

    all_levels = _load_frozen_levels(level_dir / "levels.csv", constants_sha256)
    collision_master = _load_main_collision_master(root)
    observed_overrides, observed_summary = _configured_observed_overrides(root)
    prepared: list[
        tuple[
            str,
            LamdaData,
            list[MasterLevel],
            list[dict[str, object]],
            list[SelectedCollision],
            dict[str, object],
            list[dict[str, object]],
            list[dict[str, object]],
        ]
    ] = []
    file_records: list[dict[str, object]] = []
    for species in SPECIES:
        selected_levels, template, tie_breaks = _balanced_levels(
            all_levels,
            species,
            j_max=j_max,
            k_max=abs_k_max,
            ground_max_cm1=ground_max_cm1,
        )
        key_to_index = {
            level.key: index for index, level in enumerate(selected_levels, start=1)
        }
        lamda_levels = [
            Level(index, level.energy_cm1, level.weight, level.key.encode())
            for index, level in enumerate(selected_levels, start=1)
        ]
        radiative, radiative_provenance, radiative_summary = _select_radiative(
            strength_dir / "lines.csv",
            species,
            selected_levels,
            key_to_index,
            observed_overrides,
            minimum_smu2_debye2=minimum_smu2,
            constants_sha256=constants_sha256,
        )
        calculated = _calculated_collisions(
            collision_master,
            species,
            key_to_index,
            selected_levels,
            target_temperatures,
        )
        scaled = _scale_calculated_collisions(calculated)
        if not radiative or not scaled:
            raise CH3OH13BuildError(
                f"{species}-13CH3OH produced an empty radiative or collision network"
            )
        collider_counts = Counter(row.collider_id for row in scaled)
        if set(collider_counts) != {2, 3, 6}:
            raise CH3OH13BuildError(
                f"{species}-13CH3OH is missing a p-H2, o-H2, or He partner"
            )
        filename = LAMDA_FILENAMES[species]
        data = LamdaData(
            molecule=(
                f"{SPECIES_NAMES[species]} vt=0,1,2 "
                "[12CH3OH reduced-mass collision proxy]"
            ),
            mass_amu=MOLECULAR_MASS_U,
            levels=lamda_levels,
            radiative=radiative,
            collision_partners=_partners_13ch3oh(scaled, target_temperatures),
        )
        assert_collision_selection(
            data, scaled, target_temperatures, context=filename
        )
        radiative_audit = audit_radiative_integrity(
            data,
            context=filename,
            source_rows=radiative_provenance,
            minimum_smu2_debye2=minimum_smu2,
        )
        coverage = _collision_coverage(
            filename, species, selected_levels, scaled
        )
        prepared.append(
            (
                filename,
                data,
                selected_levels,
                radiative_provenance,
                scaled,
                radiative_summary,
                tie_breaks,
                coverage,
            )
        )
        file_records.append(
            {
                "file": filename,
                "species": species,
                "levels": len(lamda_levels),
                "rotational_template_size": len(template),
                "radiative_transitions": len(radiative),
                "radiative_integrity": radiative_audit,
                "radiative_components": _radiative_components(
                    len(lamda_levels), radiative
                ),
                "radiative_selection": radiative_summary,
                "collision_rows": len(scaled),
                "collision_rows_by_collider": {
                    str(collider_id): collider_counts[collider_id]
                    for collider_id in (2, 3, 6)
                },
                "collision_coverage": coverage,
                "level_energy_serialization_tie_breaks": tie_breaks,
            }
        )

    output_dir = isotope_root / "outputs"
    manifest: dict[str, object] = {
        "schema_version": 1,
        "product_status": "validated isotope-specific spectroscopy with explicitly labelled collision proxy",
        "isotopologue": "13CH3OH",
        "separate_symmetry_species": list(SPECIES),
        "molecular_mass_u": MOLECULAR_MASS_U,
        "configuration": config,
        "hamiltonian": {
            "source": "Xu et al. (2014), Table 2",
            "doi": "10.1016/j.jms.2014.06.005",
            "reported_parameter_count": 103,
            "explicit_zero_slot_count": 102,
            "constants_sha256": constants_sha256,
            "fit_performed": False,
            "rounded_published_constants": True,
            "level_manifest_sha256": sha256_file(level_dir / "manifest.json"),
        },
        "radiative_model": {
            "dataset_id": RADIATIVE_DATASET_ID,
            "dipoles_debye": dict(DIPOLES_DEBYE),
            "cos_3gamma_terms_zero": True,
            "minimum_effective_transition_dipole_debye": math.sqrt(minimum_smu2),
            "minimum_smu2_debye2": minimum_smu2,
            "strength_manifest_sha256": sha256_file(
                strength_dir / "manifest.json"
            ),
            "observed_frequencies": observed_summary,
        },
        "collision_proxy": {
            "source_isotopologue": "12CH3OH",
            "mapping": "exact canonical quantum identity",
            "dynamics_status": ISOTOPOLOGUE_PROXY_STATUS,
            "same_temperature": True,
            "temperature_extrapolation": False,
            "scaling_by_collider": {
                str(collider_id): {
                    "name": scaling_metadata(collider_id).collider_name,
                    "factor": scaling_metadata(collider_id).factor,
                    "source_reduced_mass_u": scaling_metadata(
                        collider_id
                    ).source_reduced_mass_u,
                    "target_reduced_mass_u": scaling_metadata(
                        collider_id
                    ).target_reduced_mass_u,
                    "formula": scaling_metadata(collider_id).formula,
                    "assumption": scaling_metadata(collider_id).assumption,
                }
                for collider_id in (2, 3, 6)
            },
            "source_master_sha256": sha256_file(
                root / MAIN_COLLISION_MASTER_RELATIVE_PATH
            ),
        },
        "files": file_records,
        "collision_source_level_audit": collision_source_audit,
        "limitations": [
            "Table 2 contains rounded published constants rather than the full internal fit vector.",
            "The current CDMS tag 033502 is an older vt=0,1 entry; vt=2 frequencies remain BELGI predictions.",
            "If the optional CDMS catalog is absent, all frequencies remain BELGI predictions.",
            "Collision state identity is exact, but isotope dynamics and resonance structure are only reduced-mass proxies.",
            "No missing source collision rows are invented and no rate is extrapolated beyond its source grid.",
        ],
        "build_inputs_sha256": {
            portable_path(root, root / DEFAULT_CONFIG_RELATIVE_PATH): sha256_file(
                root / DEFAULT_CONFIG_RELATIVE_PATH
            ),
            portable_path(
                root, root / DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH
            ): sha256_file(root / DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH),
            portable_path(root, constants_path): constants_sha256,
            portable_path(
                root, isotope_root / "data/raw/user-supplied/xu-et-al-2014.pdf"
            ): sha256_file(
                isotope_root / "data/raw/user-supplied/xu-et-al-2014.pdf"
            ),
            portable_path(
                root, root / MAIN_COLLISION_MASTER_RELATIVE_PATH
            ): sha256_file(root / MAIN_COLLISION_MASTER_RELATIVE_PATH),
        },
        "upstream_run_manifests": {
            "levels": level_manifest,
            "strengths": strength_manifest,
        },
    }

    with tempfile.TemporaryDirectory(
        prefix=".ch3oh13-stage-", dir=isotope_root
    ) as temporary:
        stage_root = Path(temporary)
        stage_output = stage_root / "outputs"
        stage_output.mkdir(parents=True)
        records = {record["file"]: record for record in file_records}
        for (
            filename,
            data,
            levels,
            radiative_provenance,
            collisions,
            _,
            _,
            _,
        ) in prepared:
            data_path = stage_output / filename
            write_lamda(data_path, data)
            _write_isotope_sidecars(
                stage_root,
                filename,
                radiative_provenance,
                collisions,
                target_temperatures,
            )
            level_sidecar = _write_level_sidecar(
                stage_root, filename, levels
            )
            rad_sidecar = (
                stage_output / "provenance" / f"{filename}.radiative.csv"
            )
            collision_sidecar = (
                stage_output / "provenance" / f"{filename}.collisions.csv"
            )
            record = records[filename]
            record["sha256"] = sha256_file(data_path)
            record["levels_sidecar_sha256"] = sha256_file(level_sidecar)
            record["radiative_sidecar_sha256"] = sha256_file(rad_sidecar)
            record["collision_sidecar_sha256"] = sha256_file(collision_sidecar)

        validation = _validate_staged_products(stage_output, manifest)
        manifest["validation"] = validation
        manifest_path = stage_output / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        (output_dir / "provenance").mkdir(parents=True, exist_ok=True)
        for path in sorted((stage_output / "provenance").iterdir()):
            os.replace(path, output_dir / "provenance" / path.name)
        for path in sorted(stage_output.glob("*.dat")):
            os.replace(path, output_dir / path.name)
        os.replace(manifest_path, output_dir / "manifest.json")

    report_dir = isotope_root / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "status": "pass",
        "manifest_sha256": sha256_file(output_dir / "manifest.json"),
        "files": manifest["validation"]["files"],  # type: ignore[index]
        "scientific_status": {
            "spectroscopy": "Xu 2014 Table-2 BELGI candidate",
            "frequencies": (
                "CDMS experimental overrides where exact and available; BELGI otherwise"
            ),
            "collisions": ISOTOPOLOGUE_PROXY_STATUS,
        },
    }
    (report_dir / "validation.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def _validate_staged_products(
    output_dir: Path, manifest: dict[str, object]
) -> dict[str, object]:
    results: list[dict[str, object]] = []
    for record in manifest["files"]:  # type: ignore[union-attr]
        filename = str(record["file"])
        data = read_lamda(output_dir / filename)
        errors = validate_lamda(data)
        if not math.isclose(
            data.mass_amu, MOLECULAR_MASS_U, rel_tol=0.0, abs_tol=5.0e-7
        ):
            errors.append("molecular mass does not match exact 13CH3OH mass")
        partner_ids = [partner.collider_id for partner in data.collision_partners]
        if partner_ids != [2, 3, 6]:
            errors.append(f"collision partner IDs are {partner_ids}, expected [2, 3, 6]")
        species = str(record["species"])
        level_keys = [LevelKey.decode(level.quantum_label) for level in data.levels]
        if any(key.species != species for key in level_keys):
            errors.append("LAMDA file mixes A/E symmetry species")
        by_vt = Counter(key.vt for key in level_keys)
        if len(set(by_vt.values())) != 1 or set(by_vt) != {0, 1, 2}:
            errors.append(f"unbalanced torsional manifolds: {dict(by_vt)}")
        rad_path = output_dir / "provenance" / f"{filename}.radiative.csv"
        with rad_path.open(encoding="utf-8", newline="") as handle:
            effective_dipoles = [
                float(row["effective_transition_dipole_debye"])
                for row in csv.DictReader(handle)
            ]
        if not effective_dipoles or min(effective_dipoles) < 0.01 - 1.0e-14:
            errors.append("radiative sidecar violates the 0.01 D cutoff")
        if errors:
            raise CH3OH13BuildError(
                f"staged {filename} validation failed: {'; '.join(errors)}"
            )
        results.append(
            {
                "file": filename,
                "status": "pass",
                "levels": len(data.levels),
                "levels_by_vt": {str(vt): by_vt[vt] for vt in (0, 1, 2)},
                "radiative_transitions": len(data.radiative),
                "minimum_effective_transition_dipole_debye": min(
                    effective_dipoles
                ),
                "collision_partners": {
                    str(partner.collider_id): len(partner.transitions)
                    for partner in data.collision_partners
                },
                "independent_parse": "methanol_lamda.lamda.read_lamda",
                "format_validation": "pass",
            }
        )
    return {"status": "pass", "files": results}


def validate_products(project_root: str | Path) -> dict[str, object]:
    """Re-validate installed isotope outputs against their manifest hashes."""
    root = Path(project_root).resolve()
    output_dir = subproject_root(root) / "outputs"
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.is_file():
        raise CH3OH13BuildError("13CH3OH output manifest is missing; run the build")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validation = _validate_staged_products(output_dir, manifest)
    for record in manifest["files"]:
        filename = record["file"]
        checks = {
            "sha256": output_dir / filename,
            "levels_sidecar_sha256": output_dir
            / "provenance"
            / f"{filename}.levels.csv",
            "radiative_sidecar_sha256": output_dir
            / "provenance"
            / f"{filename}.radiative.csv",
            "collision_sidecar_sha256": output_dir
            / "provenance"
            / f"{filename}.collisions.csv",
        }
        for field, path in checks.items():
            if sha256_file(path) != record[field]:
                raise CH3OH13BuildError(
                    f"installed {path.name} does not match manifest field {field}"
                )
    return validation


__all__ = [
    "CH3OH13BuildError",
    "FrequencyOverride",
    "build_products",
    "subproject_root",
    "validate_products",
]
