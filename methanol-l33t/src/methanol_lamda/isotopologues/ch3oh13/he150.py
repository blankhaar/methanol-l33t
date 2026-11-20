"""Build the exact RF2011 150-state He-basis products for :sup:`13`CH3OH.

The state set is the complete canonical LevelKey set carried by BASECOL
collision IDs 476 (A) and 477 (E).  It is deliberately *not* an energy-rank
truncation of the isotope spectrum.  Xu et al. (2014) supplies the isotope
Hamiltonian and radiative strengths; every He coefficient is reconstructed
from the normalized main-isotopologue source row, reoriented by detailed
balance when the isotope energy ordering requires it, and finally multiplied
by the documented reduced-mass factor.
"""

from __future__ import annotations

from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

from ...build import (
    RF_TORSION_HE_SOURCE,
    MasterCollision,
    MasterLevel,
    SelectedCollision,
    _calculated_collisions,
)
from ...hamiltonian import build_frozen_levels, sha256_file
from ...independent import read_summary
from ...integrity import RawBindingError, verify_dataset_record
from ...lamda import read_lamda, validate_lamda, write_lamda
from ...models import (
    CollisionPartner,
    CollisionTransition,
    LamdaData,
    Level,
)
from ...paths import portable_path
from ...physics import einstein_a_from_smu2
from ...quantum import LevelKey
from ...radiative_integrity import audit_radiative_integrity
from ...strengths import build_frozen_strengths
from .build import (
    DEFAULT_CONFIG_RELATIVE_PATH,
    DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH,
    LEVEL_OUTPUT_RELATIVE_PATH,
    MAIN_COLLISION_MASTER_RELATIVE_PATH,
    MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE,
    STRENGTH_DIPOLE_MODEL_ID,
    STRENGTH_OUTPUT_RELATIVE_PATH,
    CH3OH13BuildError,
    _configured_observed_overrides,
    _load_config,
    _load_frozen_levels,
    _load_main_collision_master,
    _radiative_components,
    _scale_calculated_collisions,
    _select_radiative,
    _write_isotope_sidecars,
    _write_level_sidecar,
    subproject_root,
)
from .collisions import (
    ISOTOPOLOGUE_PROXY_STATUS,
    scaling_factor,
    scaling_metadata,
)
from .profile import (
    HAMILTONIAN_DATASET_ID,
    MOLECULAR_MASS_U,
    RADIATIVE_DATASET_ID,
    SPECIES,
    SPECIES_NAMES,
    load_assignment_sets,
)


PROFILE_CONFIG_RELATIVE_PATH = Path(
    "isotopologues/13ch3oh/config/rf2011-he-150.json"
)
PROFILE_ID = "rf2011_he_150_reduced_mass_proxy"
MANIFEST_FILENAME = "rf2011-he-150-reduced-mass-proxy-manifest.json"
REPORT_RELATIVE_PATH = Path(
    "isotopologues/13ch3oh/reports/rf2011-he-150-validation.json"
)
EXPECTED_LEVEL_COUNT = 150
EXPECTED_LEVELS_BY_VT = {0: 100, 1: 25, 2: 25}
EXPECTED_J_MAX_BY_VT = {0: 9, 1: 4, 2: 4}
EXPECTED_COLLISION_ROWS = EXPECTED_LEVEL_COUNT * (EXPECTED_LEVEL_COUNT - 1) // 2
EXPECTED_RADIATIVE_ROWS = {"A": 863, "E": 981}
EXPECTED_REORIENTATION_ROWS = {"A": 22, "E": 5}
EXPECTED_REORIENTATION_METHODS = {
    "A": {"source_degenerate": 14, "source_nondegenerate": 8},
    "E": {"source_degenerate": 0, "source_nondegenerate": 5},
}
EXPECTED_ALL_ZERO_ROWS = {"A": 210, "E": 216}
EXPECTED_PARTIALLY_ZERO_ROWS = {"A": 4, "E": 0}
EXPECTED_COLLISION_ROWS_BY_VT_PAIR = {
    "0-0": 4950,
    "0-1": 2500,
    "0-2": 2500,
    "1-1": 300,
    "1-2": 625,
    "2-2": 300,
}
HE_TEMPERATURES = tuple(float(value) for value in range(10, 401, 10))
FILENAMES = {
    "A": "a-13ch3oh-vt012-rf2011-he-reduced-mass-proxy.dat",
    "E": "e-13ch3oh-vt012-rf2011-he-reduced-mass-proxy.dat",
}


@dataclass(frozen=True, slots=True)
class SpectroscopyArtifacts:
    constants_path: Path
    constants_sha256: str
    level_dir: Path
    strength_dir: Path
    level_manifest: dict[str, object]
    strength_manifest: dict[str, object]


def _load_profile_config(root: Path) -> dict[str, object]:
    path = root / PROFILE_CONFIG_RELATIVE_PATH
    if not path.is_file():
        raise CH3OH13BuildError(f"missing 150-state profile configuration: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    expected_scalars = {
        "profile_id": PROFILE_ID,
        "source_dataset_id": RF_TORSION_HE_SOURCE,
        "collider_id": 6,
        "level_count": EXPECTED_LEVEL_COUNT,
        "collision_row_count": EXPECTED_COLLISION_ROWS,
    }
    for field, expected in expected_scalars.items():
        if payload.get(field) != expected:
            raise CH3OH13BuildError(
                f"150-state profile field {field!r} must equal {expected!r}"
            )
    if payload.get("source_collision_ids") != {"A": 476, "E": 477}:
        raise CH3OH13BuildError("150-state profile must bind RF2011 IDs 476 and 477")
    if payload.get("levels_by_vt") != {"0": 100, "1": 25, "2": 25}:
        raise CH3OH13BuildError("150-state profile must declare the 100/25/25 basis")
    if payload.get("j_max_by_vt") != {"0": 9, "1": 4, "2": 4}:
        raise CH3OH13BuildError("150-state profile must declare Jmax=9/4/4")
    try:
        temperatures = tuple(float(value) for value in payload["temperature_grid_k"])
        cutoff = float(payload["minimum_effective_transition_dipole_debye"])
    except (KeyError, TypeError, ValueError) as exc:
        raise CH3OH13BuildError("malformed 150-state numeric profile") from exc
    if temperatures != HE_TEMPERATURES:
        raise CH3OH13BuildError("150-state profile must retain the native 10--400 K grid")
    if not math.isclose(
        cutoff,
        MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise CH3OH13BuildError("150-state profile must use the 0.01-D cutoff")
    if payload.get("filenames") != FILENAMES:
        raise CH3OH13BuildError("150-state profile filenames do not match the contract")
    return payload


def _prepare_spectroscopy(
    root: Path, compiler: str | Path | None
) -> SpectroscopyArtifacts:
    """Regenerate isotope spectroscopy without requiring a prior full build."""
    config = _load_config(root)
    assignments = load_assignment_sets(root)
    constants_path = assignments.constants_path
    constants_sha256 = sha256_file(constants_path)
    isotope_root = subproject_root(root)
    level_dir = isotope_root / LEVEL_OUTPUT_RELATIVE_PATH
    strength_dir = isotope_root / STRENGTH_OUTPUT_RELATIVE_PATH
    j_max = int(config["j_max"])
    abs_k_max = int(config["abs_k_max"])
    minimum_smu2 = float(config["strength_generation_smu2_min_debye2"])
    level_manifest = build_frozen_levels(
        root,
        constants_path,
        root / "vendor/belgi",
        level_dir,
        j_max=j_max,
        abs_k_max=abs_k_max,
        vt_max=2,
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
        vt_max=2,
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
    return SpectroscopyArtifacts(
        constants_path=constants_path,
        constants_sha256=constants_sha256,
        level_dir=level_dir,
        strength_dir=strength_dir,
        level_manifest=level_manifest,
        strength_manifest=strength_manifest,
    )


def _read_collision_completeness(root: Path) -> dict[str, object]:
    path = root / "data/normalized/master/completeness.json"
    if not path.is_file():
        raise CH3OH13BuildError(
            "missing normalized collision completeness; run the root normalize command"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def _verified_source_rows(
    root: Path,
    species: str,
    rows: list[MasterCollision],
) -> tuple[list[MasterCollision], set[LevelKey], dict[str, object]]:
    chosen = [
        row
        for row in rows
        if row.species == species and row.source == RF_TORSION_HE_SOURCE
    ]
    record = (
        _read_collision_completeness(root)
        .get("collisions", {})
        .get(species, {})
        .get(RF_TORSION_HE_SOURCE, {})
    )
    if not isinstance(record, dict):
        raise CH3OH13BuildError(f"missing RF2011 {species} completeness record")
    try:
        verify_dataset_record(
            root,
            record,
            dataset_id=RF_TORSION_HE_SOURCE,
            normalized_transition_count=len(chosen),
        )
    except RawBindingError as exc:
        raise CH3OH13BuildError(
            f"RF2011 {species} collision source is not raw-bound: {exc}"
        ) from exc
    expected_collision_id = 476 if species == "A" else 477
    coverage = record.get("coverage")
    datasets = coverage.get("datasets") if isinstance(coverage, dict) else None
    if (
        record.get("basecol_collision_ids") != [expected_collision_id]
        or record.get("selection_complete") is not True
        or not isinstance(coverage, dict)
        or coverage.get("complete_pair_triangle_in_each_raw_dataset") is not True
        or not isinstance(datasets, list)
        or len(datasets) != 1
        or not isinstance(datasets[0], dict)
        or datasets[0].get("basecol_collision_id") != expected_collision_id
        or datasets[0].get("mapped_level_count") != EXPECTED_LEVEL_COUNT
        or datasets[0].get("normalized_transition_count")
        != EXPECTED_COLLISION_ROWS
        or tuple(float(value) for value in datasets[0].get("usable_temperature_grid_k", []))
        != HE_TEMPERATURES
    ):
        raise CH3OH13BuildError(
            f"RF2011 {species} completeness record does not bind the exact "
            f"BASECOL ID{expected_collision_id} 150-state native-grid triangle"
        )
    bound_hashes = {
        str(artifact["sha256"])
        for artifact in record.get("raw_artifacts", [])
        if isinstance(artifact, dict) and "sha256" in artifact
    }
    if len(chosen) != EXPECTED_COLLISION_ROWS:
        raise CH3OH13BuildError(
            f"RF2011 {species} source has {len(chosen)} rows; expected 11175"
        )
    for row in chosen:
        if (
            row.collider_id != 6
            or row.status != "calculated"
            or row.temperatures != HE_TEMPERATURES
            or len(row.rates) != len(HE_TEMPERATURES)
            or any(rate < 0.0 or not math.isfinite(rate) for rate in row.rates)
            or row.raw_sha256 not in bound_hashes
            or not row.source_upper_id
            or not row.source_lower_id
            or row.mapping_confidence != "exact_quantum_identity"
            or row.source_upper_energy_cm1 is None
            or row.source_lower_energy_cm1 is None
            or row.source_upper_energy_cm1 < row.source_lower_energy_cm1
        ):
            raise CH3OH13BuildError(
                f"RF2011 {species} source row violates the exact-source contract"
            )
    keys = {row.upper_key for row in chosen} | {row.lower_key for row in chosen}
    pairs = {frozenset((row.upper_key, row.lower_key)) for row in chosen}
    if (
        len(keys) != EXPECTED_LEVEL_COUNT
        or len(pairs) != EXPECTED_COLLISION_ROWS
        or any(len(pair) != 2 for pair in pairs)
    ):
        raise CH3OH13BuildError(
            f"RF2011 {species} source is not the complete 150-state pair triangle"
        )
    vt_counts = Counter(key.vt for key in keys)
    j_max_by_vt = {
        vt: max(key.j for key in keys if key.vt == vt) for vt in EXPECTED_LEVELS_BY_VT
    }
    if dict(vt_counts) != EXPECTED_LEVELS_BY_VT or j_max_by_vt != EXPECTED_J_MAX_BY_VT:
        raise CH3OH13BuildError(
            f"RF2011 {species} keys do not match the audited 100/25/25 J=9/4/4 basis"
        )
    return chosen, keys, record


def _select_isotope_levels(
    levels: list[MasterLevel], species: str, keys: set[LevelKey]
) -> list[MasterLevel]:
    selected = [level for level in levels if level.key.species == species and level.key in keys]
    if len(selected) != EXPECTED_LEVEL_COUNT or {level.key for level in selected} != keys:
        raise CH3OH13BuildError(
            f"isotope Hamiltonian does not contain all 150 RF2011 {species} identities"
        )
    selected.sort(key=lambda row: (row.energy_cm1, row.key.encode()))
    if any(
        right.energy_cm1 <= left.energy_cm1
        for left, right in zip(selected, selected[1:])
    ):
        raise CH3OH13BuildError(
            f"RF2011-basis {species}-13CH3OH energies are not strictly ordered"
        )
    if not math.isclose(selected[0].energy_cm1, 0.0, rel_tol=0.0, abs_tol=1.0e-9):
        raise CH3OH13BuildError(f"{species}-13CH3OH basis is not zeroed at ground")
    for level in selected:
        if not math.isclose(
            level.weight, 2 * level.key.j + 1, rel_tol=0.0, abs_tol=1.0e-9
        ):
            raise CH3OH13BuildError(
                f"unexpected weight for isotope level {level.key.encode()}"
            )
    return selected


def _keyset_sha256(keys: set[LevelKey]) -> str:
    text = "".join(
        f"{key.encode()}\n" for key in sorted(keys, key=lambda key: key.encode())
    )
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _vt_pair_counts(rows: list[SelectedCollision]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        left, right = sorted((row.upper_key.vt, row.lower_key.vt))
        counts[f"{left}-{right}"] += 1
    return dict(sorted(counts.items()))


def _he_partner(species: str, rows: list[SelectedCollision]) -> CollisionPartner:
    ordered = sorted(rows, key=lambda row: (row.upper, row.lower))
    source_id = 476 if species == "A" else 477
    return CollisionPartner(
        6,
        (
            f"13CH3OH-He [RF2011 BASECOL-ID{source_id}; "
            "reduced-mass-scaled 12CH3OH proxy]"
        ),
        HE_TEMPERATURES,
        [
            CollisionTransition(index, row.upper, row.lower, row.rates)
            for index, row in enumerate(ordered, start=1)
        ],
    )


def _validate_staged_products(
    output_dir: Path, manifest: dict[str, object]
) -> dict[str, object]:
    results: list[dict[str, object]] = []
    for record in manifest["files"]:  # type: ignore[union-attr]
        filename = str(record["file"])
        species = str(record["species"])
        data = read_lamda(output_dir / filename)
        errors = validate_lamda(data)
        if not math.isclose(
            data.mass_amu, MOLECULAR_MASS_U, rel_tol=0.0, abs_tol=5.0e-7
        ):
            errors.append("molecular mass does not match 13CH3OH")
        keys = [LevelKey.decode(level.quantum_label) for level in data.levels]
        level_by_index = {level.index: level for level in data.levels}
        vt_counts = Counter(key.vt for key in keys)
        if len(keys) != EXPECTED_LEVEL_COUNT or len(set(keys)) != EXPECTED_LEVEL_COUNT:
            errors.append("level set is not 150 unique identities")
        elif _keyset_sha256(set(keys)) != record.get("basis_level_key_sha256"):
            errors.append("level identities do not match the manifested RF2011 basis")
        if dict(vt_counts) != EXPECTED_LEVELS_BY_VT:
            errors.append(f"torsional counts are {dict(vt_counts)}")
        if any(key.species != species for key in keys):
            errors.append("LAMDA file mixes symmetry species")
        if data.levels and not math.isclose(
            data.levels[0].energy_cm1, 0.0, rel_tol=0.0, abs_tol=1.0e-9
        ):
            errors.append("ground-state energy is not zero")
        if any(
            not math.isclose(level.weight, 2 * key.j + 1, rel_tol=0.0, abs_tol=1.0e-9)
            for level, key in zip(data.levels, keys)
        ):
            errors.append("one or more level weights differ from 2J+1")

        partner: CollisionPartner | None = None
        zero_rows = -1
        partially_zero_rows = -1
        if [partner.collider_id for partner in data.collision_partners] != [6]:
            errors.append("150-state product must contain only collider ID 6")
        else:
            partner = data.collision_partners[0]
            if partner.temperatures_k != HE_TEMPERATURES:
                errors.append("He temperature grid is not the native 10--400 K grid")
            if len(partner.transitions) != EXPECTED_COLLISION_ROWS:
                errors.append("He collision block is not the complete 11175-row triangle")
            endpoint_pairs = {
                frozenset((row.upper, row.lower)) for row in partner.transitions
            }
            if len(endpoint_pairs) != EXPECTED_COLLISION_ROWS:
                errors.append("He collision block contains duplicate/missing endpoint pairs")
            zero_rows = sum(
                all(rate == 0.0 for rate in row.rates_cm3_s)
                for row in partner.transitions
            )
            partially_zero_rows = sum(
                any(rate == 0.0 for rate in row.rates_cm3_s)
                and not all(rate == 0.0 for rate in row.rates_cm3_s)
                for row in partner.transitions
            )
            if zero_rows != EXPECTED_ALL_ZERO_ROWS[species]:
                errors.append(f"all-zero collision count is {zero_rows}")
            if partially_zero_rows != EXPECTED_PARTIALLY_ZERO_ROWS[species]:
                errors.append(
                    f"partially-zero collision count is {partially_zero_rows}"
                )
        radiative_components = _radiative_components(len(data.levels), data.radiative)
        if len(data.radiative) != EXPECTED_RADIATIVE_ROWS[species]:
            errors.append(
                f"radiative transition count is {len(data.radiative)}, expected "
                f"{EXPECTED_RADIATIVE_ROWS[species]}"
            )
        if radiative_components != 1:
            errors.append("radiative graph is disconnected")

        rad_path = output_dir / "provenance" / f"{filename}.radiative.csv"
        with rad_path.open(encoding="utf-8", newline="") as handle:
            rad_rows = list(csv.DictReader(handle))
        radiative_by_index = {row.index: row for row in data.radiative}
        effective: list[float] = []
        seen_radiative: set[int] = set()
        if len(rad_rows) != len(data.radiative):
            errors.append("radiative provenance row count differs from the LAMDA file")
        for row in rad_rows:
            try:
                index = int(row["transition_index"])
                a_s = float(row["a_s"])
                frequency_ghz = float(row["frequency_ghz"])
                smu2 = float(row["smu2_debye2"])
                effective_dipole = float(row["effective_transition_dipole_debye"])
            except (KeyError, TypeError, ValueError):
                errors.append("radiative provenance has malformed numeric data")
                continue
            effective.append(effective_dipole)
            transition = radiative_by_index.get(index)
            if index in seen_radiative or transition is None:
                errors.append("radiative provenance has duplicate/unknown indices")
                continue
            seen_radiative.add(index)
            upper = level_by_index[transition.upper]
            lower = level_by_index[transition.lower]
            try:
                reconstructed_a = einstein_a_from_smu2(
                    frequency_ghz * 1000.0, smu2, upper.weight
                )
            except ValueError:
                errors.append("radiative provenance has non-positive physical values")
                continue
            if (
                row.get("file") != filename
                or row.get("upper_key") != upper.quantum_label
                or row.get("lower_key") != lower.quantum_label
                or row.get("a_source") != RADIATIVE_DATASET_ID
                or row.get("mapping_confidence") != "exact_quantum_identity"
                or not math.isclose(a_s, transition.a_s, rel_tol=5.0e-10)
                or not math.isclose(
                    frequency_ghz, transition.frequency_ghz, rel_tol=5.0e-10
                )
                or not math.isclose(
                    effective_dipole, math.sqrt(smu2), rel_tol=5.0e-13
                )
                or not math.isclose(a_s, reconstructed_a, rel_tol=5.0e-12)
            ):
                errors.append("radiative provenance is inconsistent with S mu^2/file")
                break
        if (
            not effective
            or min(effective)
            < MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE - 1.0e-14
        ):
            errors.append("radiative provenance violates the 0.01-D contract")

        collision_path = output_dir / "provenance" / f"{filename}.collisions.csv"
        with collision_path.open(encoding="utf-8", newline="") as handle:
            collision_rows = list(csv.DictReader(handle))
        factor = scaling_factor(6)
        actual_collisions = (
            {row.index: row for row in partner.transitions} if partner is not None else {}
        )
        try:
            source_record = manifest["collision_proxy"]["source_records"][species]  # type: ignore[index]
            bound_hashes = {
                str(artifact["sha256"])
                for artifact in source_record["raw_artifacts"]
            }
        except (KeyError, TypeError):
            bound_hashes = set()
            errors.append("manifest lacks RF2011 raw-artifact bindings")
        if len(collision_rows) != EXPECTED_COLLISION_ROWS:
            errors.append("collision provenance row count is not 11175")
        seen_collisions: set[int] = set()
        reorientation_rows = 0
        reorientation_methods = {
            "source_degenerate": 0,
            "source_nondegenerate": 0,
        }
        for row in collision_rows:
            try:
                index = int(row["transition_index"])
                upper_index = int(row["upper_index"])
                lower_index = int(row["lower_index"])
                multiplier = float(row["multiplier"])
                temperatures = tuple(
                    float(value) for value in row["temperatures_k"].split(";")
                )
                rates = tuple(float(value) for value in row["rates_cm3_s"].split(";"))
                source_temperatures = tuple(
                    float(value) for value in row["source_temperatures_k"].split(";")
                )
                source_rates = tuple(
                    float(value) for value in row["source_rates_cm3_s"].split(";")
                )
                source_upper_energy = float(row["source_upper_energy_cm1"])
                source_lower_energy = float(row["source_lower_energy_cm1"])
                source_upper_key = LevelKey.decode(row["source_upper_key"])
                source_lower_key = LevelKey.decode(row["source_lower_key"])
            except (KeyError, TypeError, ValueError):
                errors.append("collision provenance has malformed data")
                continue
            transition = actual_collisions.get(index)
            if index in seen_collisions or transition is None:
                errors.append("collision provenance has duplicate/unknown indices")
                continue
            seen_collisions.add(index)
            transform = row.get("transform", "")
            reorientation_rows += "reorientation" in transform
            reorientation_methods["source_degenerate"] += (
                "source_degenerate_detailed_balance_reorientation" in transform
            )
            reorientation_methods["source_nondegenerate"] += (
                "source_downward_to_upward_detailed_balance_reorientation"
                in transform
            )
            if (
                row.get("file") != filename
                or row.get("collider_id") != "6"
                or upper_index != transition.upper
                or lower_index != transition.lower
                or row.get("upper_key")
                != level_by_index[transition.upper].quantum_label
                or row.get("lower_key")
                != level_by_index[transition.lower].quantum_label
                or temperatures != HE_TEMPERATURES
                or source_temperatures != HE_TEMPERATURES
                or len(rates) != len(HE_TEMPERATURES)
                or len(source_rates) != len(HE_TEMPERATURES)
                or any(not math.isfinite(value) or value < 0.0 for value in rates)
                or any(
                    not math.isfinite(value) or value < 0.0 for value in source_rates
                )
                or any(
                    not math.isclose(expected, actual, rel_tol=5.0e-9, abs_tol=1.0e-99)
                    for expected, actual in zip(rates, transition.rates_cm3_s)
                )
                or row.get("source") != RF_TORSION_HE_SOURCE
                or row.get("status") != ISOTOPOLOGUE_PROXY_STATUS
                or row.get("mapping_confidence") != "exact_quantum_identity"
                or not math.isclose(multiplier, factor, rel_tol=0.0, abs_tol=1.0e-15)
                or "sqrt(mu12/mu13)" not in transform
                or row.get("raw_sha256") not in bound_hashes
                or not row.get("source_upper_id")
                or not row.get("source_lower_id")
                or source_upper_energy < source_lower_energy
                or frozenset((source_upper_key, source_lower_key))
                != frozenset(
                    (
                        LevelKey.decode(row["upper_key"]),
                        LevelKey.decode(row["lower_key"]),
                    )
                )
            ):
                errors.append("collision provenance violates the proxy contract")
                break
        if reorientation_rows != EXPECTED_REORIENTATION_ROWS[species]:
            errors.append(
                f"source-basis reorientation count is {reorientation_rows}"
            )
        if reorientation_methods != EXPECTED_REORIENTATION_METHODS[species]:
            errors.append(
                f"source-basis reorientation methods are {reorientation_methods}"
            )
        if (
            record.get("radiative_transitions") != EXPECTED_RADIATIVE_ROWS[species]
            or record.get("radiative_components") != 1
            or record.get("collision_rows") != EXPECTED_COLLISION_ROWS
            or record.get("collision_rows_by_vt_pair")
            != EXPECTED_COLLISION_ROWS_BY_VT_PAIR
            or record.get("source_basis_reorientation_rows")
            != EXPECTED_REORIENTATION_ROWS[species]
            or record.get("all_zero_collision_rows")
            != EXPECTED_ALL_ZERO_ROWS[species]
            or record.get("partially_zero_collision_rows")
            != EXPECTED_PARTIALLY_ZERO_ROWS[species]
        ):
            errors.append("manifested file statistics violate the fixed profile")
        if errors:
            raise CH3OH13BuildError(
                f"staged {filename} validation failed: {'; '.join(errors)}"
            )
        summary = read_summary(output_dir / filename)
        if summary.collider_counts != ((6, EXPECTED_COLLISION_ROWS, len(HE_TEMPERATURES)),):
            raise CH3OH13BuildError(f"{filename} failed independent-reader validation")
        results.append(
            {
                "file": filename,
                "status": "pass",
                "levels": len(data.levels),
                "levels_by_vt": {str(vt): vt_counts[vt] for vt in (0, 1, 2)},
                "radiative_transitions": len(data.radiative),
                "radiative_components": 1,
                "minimum_effective_transition_dipole_debye": min(effective),
                "collision_partner": 6,
                "collision_rows": EXPECTED_COLLISION_ROWS,
                "temperature_count": len(HE_TEMPERATURES),
                "source_basis_reorientation_rows": reorientation_rows,
                "all_zero_collision_rows": zero_rows,
                "partially_zero_collision_rows": partially_zero_rows,
                "radiative_einstein_a_reconstruction": "pass",
                "collision_sidecar_consistency": "pass",
                "independent_parse": "methanol_lamda.independent.read_summary",
                "format_validation": "pass",
            }
        )
    return {"status": "pass", "files": results}


def build_he150_products(
    project_root: str | Path,
    *,
    compiler: str | Path | None = None,
) -> dict[str, object]:
    """Generate the A/E 150-state isotope He-proxy pair."""
    root = Path(project_root).resolve()
    isotope_root = subproject_root(root)
    profile = _load_profile_config(root)
    spectroscopy = _prepare_spectroscopy(root, compiler)
    all_levels = _load_frozen_levels(
        spectroscopy.level_dir / "levels.csv", spectroscopy.constants_sha256
    )
    collision_master = _load_main_collision_master(root)
    overrides, observed_summary = _configured_observed_overrides(root)
    minimum_smu2 = MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE**2
    prepared: list[
        tuple[
            str,
            LamdaData,
            list[MasterLevel],
            list[dict[str, object]],
            list[SelectedCollision],
            dict[str, object],
        ]
    ] = []
    file_records: list[dict[str, object]] = []
    source_records: dict[str, object] = {}

    for species in SPECIES:
        source_rows, source_keys, source_record = _verified_source_rows(
            root, species, collision_master
        )
        source_records[species] = {
            "dataset_id": source_record["dataset_id"],
            "normalized_transition_count": source_record[
                "normalized_transition_count"
            ],
            "raw_artifacts": source_record["raw_artifacts"],
        }
        selected_levels = _select_isotope_levels(all_levels, species, source_keys)
        key_to_index = {
            level.key: index for index, level in enumerate(selected_levels, start=1)
        }
        lamda_levels = [
            Level(index, level.energy_cm1, level.weight, level.key.encode())
            for index, level in enumerate(selected_levels, start=1)
        ]
        radiative, radiative_provenance, radiative_summary = _select_radiative(
            spectroscopy.strength_dir / "lines.csv",
            species,
            selected_levels,
            key_to_index,
            overrides,
            minimum_smu2_debye2=minimum_smu2,
            constants_sha256=spectroscopy.constants_sha256,
        )
        calculated = _calculated_collisions(
            source_rows,
            species,
            key_to_index,
            selected_levels,
            HE_TEMPERATURES,
        )
        scaled = _scale_calculated_collisions(calculated)
        if len(scaled) != EXPECTED_COLLISION_ROWS:
            raise CH3OH13BuildError(
                f"{species}-13CH3OH retained {len(scaled)} RF2011 rows; expected 11175"
            )
        if any(row.collider_id != 6 for row in scaled):
            raise CH3OH13BuildError("150-state profile contains a non-He collision")
        filename = FILENAMES[species]
        data = LamdaData(
            molecule=(
                f"{SPECIES_NAMES[species]} vt=0,1,2 "
                "[RF2011-He-150; reduced-mass collision proxy]"
            ),
            mass_amu=MOLECULAR_MASS_U,
            levels=lamda_levels,
            radiative=radiative,
            collision_partners=[_he_partner(species, scaled)],
        )
        radiative_audit = audit_radiative_integrity(
            data,
            context=filename,
            source_rows=radiative_provenance,
            minimum_smu2_debye2=minimum_smu2,
        )
        errors = validate_lamda(data)
        if errors:
            raise CH3OH13BuildError(
                f"{filename} failed preflight: {'; '.join(errors[:5])}"
            )
        status_counts = Counter(row.status for row in scaled)
        transform_counts = Counter(row.transform for row in scaled)
        record: dict[str, object] = {
            "file": filename,
            "species": species,
            "product_profile": PROFILE_ID,
            "levels": len(lamda_levels),
            "levels_by_vt": {
                str(vt): sum(level.key.vt == vt for level in selected_levels)
                for vt in (0, 1, 2)
            },
            "basis_level_key_sha256": _keyset_sha256(source_keys),
            "radiative_transitions": len(radiative),
            "radiative_integrity": radiative_audit,
            "radiative_components": _radiative_components(len(lamda_levels), radiative),
            "radiative_selection": radiative_summary,
            "collision_rows": len(scaled),
            "collision_rows_by_vt_pair": _vt_pair_counts(scaled),
            "collision_rows_by_status": dict(sorted(status_counts.items())),
            "collision_transforms": dict(sorted(transform_counts.items())),
            "source_basis_reorientation_rows": sum(
                "reorientation" in row.transform for row in scaled
            ),
            "all_zero_collision_rows": sum(
                all(rate == 0.0 for rate in row.rates) for row in scaled
            ),
            "partially_zero_collision_rows": sum(
                any(rate == 0.0 for rate in row.rates)
                and not all(rate == 0.0 for rate in row.rates)
                for row in scaled
            ),
            "collision_partner_ids": [6],
            "temperature_count": len(HE_TEMPERATURES),
            "temperature_min_k": HE_TEMPERATURES[0],
            "temperature_max_k": HE_TEMPERATURES[-1],
        }
        if (
            record["radiative_transitions"] != EXPECTED_RADIATIVE_ROWS[species]
            or record["radiative_components"] != 1
            or record["collision_rows_by_vt_pair"]
            != EXPECTED_COLLISION_ROWS_BY_VT_PAIR
            or record["source_basis_reorientation_rows"]
            != EXPECTED_REORIENTATION_ROWS[species]
            or record["all_zero_collision_rows"]
            != EXPECTED_ALL_ZERO_ROWS[species]
            or record["partially_zero_collision_rows"]
            != EXPECTED_PARTIALLY_ZERO_ROWS[species]
            or {
                "source_degenerate": sum(
                    count
                    for transform, count in transform_counts.items()
                    if "source_degenerate_detailed_balance_reorientation"
                    in transform
                ),
                "source_nondegenerate": sum(
                    count
                    for transform, count in transform_counts.items()
                    if "source_downward_to_upward_detailed_balance_reorientation"
                    in transform
                ),
            }
            != EXPECTED_REORIENTATION_METHODS[species]
        ):
            raise CH3OH13BuildError(
                f"{filename} statistics differ from the audited RF2011-150 profile"
            )
        file_records.append(record)
        prepared.append(
            (
                filename,
                data,
                selected_levels,
                radiative_provenance,
                scaled,
                record,
            )
        )

    factor = scaling_metadata(6)
    completeness_path = root / "data/normalized/master/completeness.json"
    profile_path = root / PROFILE_CONFIG_RELATIVE_PATH
    manifest: dict[str, object] = {
        "schema_version": 1,
        "product_profile": PROFILE_ID,
        "isotopologue": "13CH3OH",
        "scope": (
            "exact unbalanced RF2011 ID476/477 scattering basis; 100/25/25 "
            "levels in vt=0/1/2; diagnostic and not a replacement for the "
            "balanced 1089-level Class II products"
        ),
        "profile_configuration": profile,
        "molecular_mass_u": MOLECULAR_MASS_U,
        "level_energy_policy": (
            "Xu et al. 2014 Table-2 BELGI isotope term values; RF2011 source "
            "energies retained only for detailed-balance provenance"
        ),
        "radiative_model": {
            "dataset_id": RADIATIVE_DATASET_ID,
            "minimum_effective_transition_dipole_debye": (
                MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE
            ),
            "minimum_smu2_debye2": minimum_smu2,
            "observed_frequencies": observed_summary,
            "strength_manifest_sha256": sha256_file(
                spectroscopy.strength_dir / "manifest.json"
            ),
        },
        "collision_proxy": {
            "source_isotopologue": "12CH3OH",
            "source_dataset_id": RF_TORSION_HE_SOURCE,
            "source_collision_ids": {"A": 476, "E": 477},
            "basis_mapping": "exact canonical quantum identity",
            "dynamics_status": ISOTOPOLOGUE_PROXY_STATUS,
            "collider_id": 6,
            "native_temperature_grid_k": list(HE_TEMPERATURES),
            "same_temperature": True,
            "temperature_extrapolation": False,
            "source_reorientation_before_scaling": True,
            "scaling": {
                "factor": factor.factor,
                "formula": factor.formula,
                "assumption": factor.assumption,
                "source_reduced_mass_u": factor.source_reduced_mass_u,
                "target_reduced_mass_u": factor.target_reduced_mass_u,
            },
            "source_records": source_records,
            "source_master_sha256": sha256_file(
                root / MAIN_COLLISION_MASTER_RELATIVE_PATH
            ),
            "source_completeness_sha256": sha256_file(completeness_path),
        },
        "hamiltonian": {
            "dataset_id": HAMILTONIAN_DATASET_ID,
            "constants_sha256": spectroscopy.constants_sha256,
            "level_manifest_sha256": sha256_file(
                spectroscopy.level_dir / "manifest.json"
            ),
        },
        "files": file_records,
        "limitations": [
            "The 100/25/25 source basis is deliberately unbalanced and excludes the 23.1-GHz J=10 system.",
            "Table 2 contains rounded published constants rather than the authors' internal fit vector.",
            "When the optional CDMS catalog is absent, all transition frequencies remain BELGI predictions.",
            "Collision identities are exact, but isotope dynamics and resonances are reduced-mass proxies.",
        ],
        "build_inputs_sha256": {
            portable_path(root, profile_path): sha256_file(profile_path),
            portable_path(root, root / DEFAULT_CONFIG_RELATIVE_PATH): sha256_file(
                root / DEFAULT_CONFIG_RELATIVE_PATH
            ),
            portable_path(
                root, root / DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH
            ): sha256_file(root / DEFAULT_OBSERVED_CONFIG_RELATIVE_PATH),
            portable_path(root, spectroscopy.constants_path): spectroscopy.constants_sha256,
            portable_path(
                root,
                isotope_root / "data/raw/user-supplied/xu-et-al-2014.pdf",
            ): sha256_file(
                isotope_root / "data/raw/user-supplied/xu-et-al-2014.pdf"
            ),
            portable_path(
                root, root / MAIN_COLLISION_MASTER_RELATIVE_PATH
            ): sha256_file(root / MAIN_COLLISION_MASTER_RELATIVE_PATH),
            portable_path(root, completeness_path): sha256_file(completeness_path),
        },
        "upstream_run_manifests": {
            "levels": spectroscopy.level_manifest,
            "strengths": spectroscopy.strength_manifest,
        },
    }

    output_dir = isotope_root / "outputs"
    with tempfile.TemporaryDirectory(
        prefix=".ch3oh13-stage-he150-", dir=isotope_root
    ) as temporary:
        stage_output = Path(temporary) / "outputs"
        stage_output.mkdir(parents=True)
        for filename, data, levels, radiative_rows, collisions, record in prepared:
            data_path = stage_output / filename
            write_lamda(data_path, data)
            stage_root = stage_output.parent
            _write_isotope_sidecars(
                stage_root,
                filename,
                radiative_rows,
                collisions,
                HE_TEMPERATURES,
            )
            level_path = _write_level_sidecar(stage_root, filename, levels)
            radiative_path = stage_output / "provenance" / f"{filename}.radiative.csv"
            collision_path = stage_output / "provenance" / f"{filename}.collisions.csv"
            record.update(
                {
                    "sha256": sha256_file(data_path),
                    "levels_sidecar_sha256": sha256_file(level_path),
                    "radiative_sidecar_sha256": sha256_file(radiative_path),
                    "collision_sidecar_sha256": sha256_file(collision_path),
                }
            )
        validation = _validate_staged_products(stage_output, manifest)
        manifest["validation"] = validation
        manifest_path = stage_output / MANIFEST_FILENAME
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (output_dir / "provenance").mkdir(parents=True, exist_ok=True)
        for path in sorted((stage_output / "provenance").iterdir()):
            os.replace(path, output_dir / "provenance" / path.name)
        for path in sorted(stage_output.glob("*.dat")):
            os.replace(path, output_dir / path.name)
        os.replace(manifest_path, output_dir / MANIFEST_FILENAME)

    report_path = root / REPORT_RELATIVE_PATH
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "status": "pass",
        "product_profile": PROFILE_ID,
        "manifest_sha256": sha256_file(output_dir / MANIFEST_FILENAME),
        "files": manifest["validation"]["files"],  # type: ignore[index]
        "scientific_status": {
            "spectroscopy": "Xu 2014 Table-2 BELGI candidate",
            "frequencies": (
                "CDMS experimental overrides where exact and available; BELGI otherwise"
            ),
            "collisions": ISOTOPOLOGUE_PROXY_STATUS,
            "basis": "exact RF2011 ID476/477 100/25/25 scattering basis",
        },
    }
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _reconstruct_and_compare_collisions(
    root: Path, manifest: dict[str, object]
) -> None:
    """Independently reconstruct every installed rate from normalized sources."""
    output_dir = subproject_root(root) / "outputs"
    collision_master = _load_main_collision_master(root)
    for record in manifest["files"]:  # type: ignore[union-attr]
        species = str(record["species"])
        filename = str(record["file"])
        data = read_lamda(output_dir / filename)
        source_rows, source_keys, _ = _verified_source_rows(
            root, species, collision_master
        )
        file_keys = {LevelKey.decode(level.quantum_label) for level in data.levels}
        if file_keys != source_keys:
            raise CH3OH13BuildError(
                f"installed {filename} does not use the exact RF2011 key set"
            )
        levels = [
            MasterLevel(
                key=LevelKey.decode(level.quantum_label),
                energy_cm1=level.energy_cm1,
                weight=level.weight,
                source=HAMILTONIAN_DATASET_ID,
                source_index=str(level.index),
                mapping_confidence="exact_quantum_identity",
            )
            for level in data.levels
        ]
        key_to_index = {
            LevelKey.decode(level.quantum_label): level.index for level in data.levels
        }
        expected = _scale_calculated_collisions(
            _calculated_collisions(
                source_rows,
                species,
                key_to_index,
                levels,
                HE_TEMPERATURES,
            )
        )
        expected.sort(key=lambda row: (row.upper, row.lower))
        actual = data.collision_partners[0].transitions
        if len(expected) != len(actual):
            raise CH3OH13BuildError(f"installed {filename} collision count changed")
        for expected_row, actual_row in zip(expected, actual):
            if (
                expected_row.upper != actual_row.upper
                or expected_row.lower != actual_row.lower
                or any(
                    not math.isclose(
                        expected_rate,
                        actual_rate,
                        rel_tol=5.0e-9,
                        abs_tol=1.0e-99,
                    )
                    for expected_rate, actual_rate in zip(
                        expected_row.rates, actual_row.rates_cm3_s
                    )
                )
            ):
                raise CH3OH13BuildError(
                    f"installed {filename} collision reconstruction mismatch"
                )


def validate_he150_products(project_root: str | Path) -> dict[str, object]:
    """Validate hashes, structure, exact basis, and every installed rate."""
    root = Path(project_root).resolve()
    _load_profile_config(root)
    output_dir = subproject_root(root) / "outputs"
    manifest_path = output_dir / MANIFEST_FILENAME
    if not manifest_path.is_file():
        raise CH3OH13BuildError(
            "150-state isotope manifest is missing; run the build-150 command"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("product_profile") != PROFILE_ID:
        raise CH3OH13BuildError("installed 150-state manifest has the wrong profile")
    validation = _validate_staged_products(output_dir, manifest)
    for record in manifest["files"]:
        filename = str(record["file"])
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
    _reconstruct_and_compare_collisions(root, manifest)
    validation["collision_reconstruction"] = (
        "all 22350 rows reproduced from raw-bound source grids, detailed-balance "
        "orientation, and reduced-mass scaling"
    )
    return validation


__all__ = [
    "FILENAMES",
    "HE_TEMPERATURES",
    "MANIFEST_FILENAME",
    "PROFILE_ID",
    "build_he150_products",
    "validate_he150_products",
]
