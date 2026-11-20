"""Build ``vt=0..4 + vCO=1`` methanol LAMDA inputs.

The CO-stretch transition operator is evaluated by BELGI with a 0.27-D
body-fixed ``a``-axis moment and the ground-state rotational eigenfunctions.
No nominal-K selection rule is imposed.  Main-isotopologue upper-state terms
come from the Xu/HITRAN-format line compilation; the carbon-13 analogue uses
the published Mukhopadhyay et al. subband origins and K=0 B/D differences.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
import csv
import io
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import tempfile

from .build import MasterLevel, SelectedCollision, _calculated_collisions, _partners, load_master
from .collision_integrity import (
    assert_collision_selection,
    assert_parent_collision_inheritance,
    audit_basecol_collision_master,
)
from .extended import (
    ExtendedBuildError,
    _assert_manifest_binds_file,
    _copy_vt0_intraband,
    _inherit_radiative_subset,
    _load_extended_levels,
    _prefer_native_calculated_collisions,
    _select_levels,
    _select_radiative,
)
from .hamiltonian import (
    _compile_driver,
    _find_compiler,
    parse_final_constants,
    render_complete_cst_vector,
    sha256_file,
)
from .independent import read_summary
from .lamda import read_lamda, validate_lamda, write_lamda
from .models import LamdaData, Level, RadiativeTransition
from .physics import einstein_a_from_smu2
from .quantum import LevelKey
from .radiative_integrity import audit_radiative_integrity
from .strengths import _parse_strength_summary


CM1_TO_MHZ = 29979.2458
MAIN_HITRAN_RELATIVE_PATH = Path(
    "data/raw/third-party-mirror/hitran2004-ch3oh-950-1300/"
    "ch3oh-950-1300.hitran2004"
)
MAIN_HITRAN_METADATA_RELATIVE_PATH = MAIN_HITRAN_RELATIVE_PATH.with_name(
    "source.json"
)
ISOTOPE_THESIS_RELATIVE_PATH = Path(
    "data/raw/official/library-and-archives-canada/nq23873.pdf"
)
MAIN_LEVELS = Path("data/normalized/belgi-xu2008-vt04-frozen/levels.csv")
MAIN_LINES = Path("data/normalized/belgi-xu2008-vt04-constant-dipole/lines.csv")
ISOTOPE_LEVELS = Path(
    "isotopologues/13ch3oh/data/normalized/belgi-xu2014-vt04-frozen/levels.csv"
)
ISOTOPE_LINES = Path(
    "isotopologues/13ch3oh/data/normalized/"
    "belgi-xu2014-vt04-constant-dipole/lines.csv"
)
ISOTOPE_CONSTANTS = Path(
    "isotopologues/13ch3oh/constants/table2-2014-constants.txt"
)
ISOTOPE_CO_MODEL = Path(
    "isotopologues/13ch3oh/data/raw/transcribed/"
    "mukhopadhyay1995-co-stretch.json"
)
GENERATION_FLOOR_SMU2 = 1.0e-6
PRODUCT_FLOOR_SMU2 = 1.0e-4
VCO_TRANSITION_MOMENT_DEBYE = 0.27
HITRAN_E_GROUND_ZERO_CM1 = 5.4897
HITRAN_LOWER_ENERGY_TOLERANCE_CM1 = 2.5e-4


class COStretchError(ExtendedBuildError):
    """Raised when a CO-stretch source or generated product is inconsistent."""


@dataclass(frozen=True, slots=True, order=True)
class COStretchKey:
    species: str
    j: int
    k: int
    parity: str = ""
    vt: int = 5

    def __post_init__(self) -> None:
        ground = LevelKey(self.species, 0, self.j, self.k, self.parity)
        object.__setattr__(self, "species", ground.species)

    @property
    def rotational_identity(self) -> tuple[int, int, str]:
        return self.j, self.k, self.parity

    def encode(self) -> str:
        if self.species == "A":
            return f"A:vCO1:vt0:J{self.j}:K{self.k}:{self.parity}"
        return f"E:vCO1:vt0:J{self.j}:K{self.k:+d}"

    @classmethod
    def from_ground(cls, key: LevelKey) -> "COStretchKey":
        if key.vt != 0:
            raise COStretchError("CO-stretch identity must be copied from vt=0")
        return cls(key.species, key.j, key.k, key.parity)


@dataclass(frozen=True, slots=True)
class RotationalStrength:
    left: LevelKey
    right: LevelKey
    smu2_debye2: float


def _extract_csv(stdout: str, begin: str, end: str) -> str:
    marker = begin + "\n"
    if marker not in stdout or end not in stdout:
        raise COStretchError("BELGI CO-stretch output is missing CSV markers")
    return stdout.split(marker, 1)[1].split("\n" + end, 1)[0]


def _run_transition_operator(
    root: Path,
    constants_path: Path,
    *,
    compiler: str | Path | None,
) -> tuple[list[RotationalStrength], dict[str, int | float]]:
    assignments = parse_final_constants(constants_path, require_xu_final=False)
    vector = render_complete_cst_vector(assignments)
    compiler_path = _find_compiler(compiler)
    with tempfile.TemporaryDirectory(prefix="belgi-vco-strengths-") as temporary:
        executable, _, _, _ = _compile_driver(
            root,
            root / "vendor/belgi",
            Path(temporary),
            compiler_path,
            driver_filename="belgi_frozen_strengths.f",
            executable_name="belgi_frozen_strengths_vco",
        )
        run = subprocess.run(
            [str(executable), "22", "14", "0", f"{GENERATION_FLOOR_SMU2:.17g}", "VCO_A"],
            input=vector,
            text=True,
            capture_output=True,
            cwd=temporary,
            check=False,
        )
    if run.returncode:
        raise COStretchError(
            f"BELGI CO-stretch operator failed ({run.returncode}): {run.stderr}"
        )
    rows = csv.DictReader(
        io.StringIO(
            _extract_csv(
                run.stdout,
                "BEGIN_FROZEN_STRENGTHS_CSV",
                "END_FROZEN_STRENGTHS_CSV",
            )
        )
    )
    result: list[RotationalStrength] = []
    seen: set[frozenset[LevelKey]] = set()
    for row in rows:
        species = row["species"].strip().upper()
        left = LevelKey(
            species,
            0,
            int(row["upper_j"]),
            int(row["upper_k"]),
            row["upper_parity"].strip(),
        )
        right = LevelKey(
            species,
            0,
            int(row["lower_j"]),
            int(row["lower_k"]),
            row["lower_parity"].strip(),
        )
        identity = frozenset((left, right))
        if identity in seen:
            raise COStretchError("duplicate rotational pair from BELGI CO operator")
        seen.add(identity)
        result.append(RotationalStrength(left, right, float(row["smu2_debye2"])))
    return result, _parse_strength_summary(run.stdout)


def _parse_hitran_local(label: str) -> LevelKey:
    fields = label.replace('"', " ").split()
    if len(fields) != 3:
        raise COStretchError(f"unrecognized HITRAN methanol label {label!r}")
    j, k, symmetry = int(fields[0]), int(fields[1]), fields[2].upper()
    if symmetry in {"A+", "A-"}:
        return LevelKey("A", 0, j, k, symmetry[-1])
    if symmetry == "E1":
        return LevelKey("E", 0, j, k)
    if symmetry == "E2":
        return LevelKey("E", 0, j, -k)
    raise COStretchError(f"unsupported HITRAN methanol symmetry {symmetry!r}")


def _main_ground_energies(path: Path) -> dict[LevelKey, float]:
    """Load the model's canonical ground-state terms used for source mapping."""

    result = {
        row.key: row.energy_cm1
        for row in _load_extended_levels(path)
        if isinstance(row.key, LevelKey) and row.key.vt == 0
    }
    if not result:
        raise COStretchError("main ground-state reference is empty")
    return result


def _resolve_hitran_ground_identity(
    raw: LevelKey,
    lower_energy_cm1: float,
    ground_energies: dict[LevelKey, float],
) -> tuple[LevelKey | None, float | None, bool]:
    """Resolve a HITRAN lower label against the frozen Hamiltonian energy.

    The archived high-|K| V8 rows use ``E1`` for both signed E ladders.  The
    lower-state term printed on the same row uniquely identifies the intended
    sign.  A parity is never changed by this resolver: its explicit source
    label is retained even when a near-degenerate partner has a similar term.
    """

    candidates = [raw]
    if raw.species == "E" and raw.k:
        candidates.append(LevelKey("E", 0, raw.j, -raw.k))
    candidates = [key for key in candidates if key in ground_energies]
    if not candidates:
        return None, None, False
    source_relative = lower_energy_cm1 - (
        HITRAN_E_GROUND_ZERO_CM1 if raw.species == "E" else 0.0
    )
    matches = sorted(
        (
            abs(ground_energies[key] - source_relative),
            key,
        )
        for key in candidates
        if abs(ground_energies[key] - source_relative)
        <= HITRAN_LOWER_ENERGY_TOLERANCE_CM1
    )
    if len(matches) != 1:
        detail = ", ".join(
            f"{key.encode()}:{abs(ground_energies[key] - source_relative):.9g}"
            for key in candidates
        )
        raise COStretchError(
            f"HITRAN lower identity {raw.encode()} is not uniquely energy-bound "
            f"within {HITRAN_LOWER_ENERGY_TOLERANCE_CM1:g} cm-1 ({detail})"
        )
    residual, resolved = matches[0]
    return resolved, residual, resolved != raw


def _resolve_hitran_upper_identity(
    raw_upper: LevelKey,
    raw_lower: LevelKey,
    resolved_lower: LevelKey,
) -> LevelKey:
    """Carry an energy-resolved E sign across a same-|K| source transition."""

    if (
        raw_upper.species == "E"
        and abs(raw_upper.k) == abs(raw_lower.k)
        and raw_upper.k
    ):
        sign = -1 if resolved_lower.k < 0 else 1
        return LevelKey("E", 0, raw_upper.j, sign * abs(raw_upper.k))
    return raw_upper


def _main_co_terms(
    path: Path,
    ground_levels_path: Path,
) -> tuple[dict[LevelKey, float], dict[str, object]]:
    ground_energies = _main_ground_energies(ground_levels_path)
    candidates: dict[LevelKey, list[tuple[float, float, int]]] = defaultdict(list)
    source_rows_total = 0
    source_rows_mapped = 0
    source_rows_outside_model = 0
    e_sign_corrections = 0
    lower_residuals: list[float] = []
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        if line[67:82].strip() != "V8" or line[82:97].strip() != "GROUND":
            continue
        raw_upper = _parse_hitran_local(line[97:112])
        raw_lower = _parse_hitran_local(line[112:127])
        source_rows_total += 1
        if raw_upper.species != raw_lower.species:
            raise COStretchError("HITRAN V8 row changes A/E nuclear-spin species")
        wavenumber = float(line[3:15])
        intensity = float(line[15:25])
        einstein_a = float(line[25:35])
        lower_energy = float(line[45:55])
        lower, residual, corrected = _resolve_hitran_ground_identity(
            raw_lower, lower_energy, ground_energies
        )
        if lower is None:
            source_rows_outside_model += 1
            continue
        upper = _resolve_hitran_upper_identity(raw_upper, raw_lower, lower)
        if upper not in ground_energies:
            source_rows_outside_model += 1
            continue
        if abs(upper.j - lower.j) > 1 or (upper.j == 0 and lower.j == 0):
            raise COStretchError(
                f"V8 source row {number} violates the rank-1 rotational rule"
            )
        if residual is None:
            raise RuntimeError('Invariant failed in co_stretch.py: residual is not None')
        lower_residuals.append(residual)
        e_sign_corrections += int(corrected)
        score = einstein_a if einstein_a > 0.0 else intensity
        candidates[upper].append(
            (ground_energies[lower] + wavenumber, score, number)
        )
        source_rows_mapped += 1
    if not candidates:
        raise COStretchError("no V8 -> GROUND lines found in HITRAN-format source")
    terms: dict[LevelKey, float] = {}
    conflict_rows = 0
    maximum_spread = 0.0
    for key, values in candidates.items():
        ordered = sorted(values)
        clusters: list[list[tuple[float, float, int]]] = []
        for value in ordered:
            if not clusters or value[0] - clusters[-1][-1][0] > 0.005:
                clusters.append([value])
            else:
                clusters[-1].append(value)
        chosen = max(
            clusters,
            key=lambda cluster: (sum(value[1] for value in cluster), len(cluster)),
        )
        representative = max(chosen, key=lambda value: value[1])
        terms[key] = representative[0]
        conflict_rows += len(values) - len(chosen)
        maximum_spread = max(maximum_spread, max(v[0] for v in values) - min(v[0] for v in values))
    return terms, {
        "source_rows_total": source_rows_total,
        "source_rows_mapped_to_frozen_ground_model": source_rows_mapped,
        "source_rows_outside_frozen_ground_model": source_rows_outside_model,
        "unique_upper_levels": len(terms),
        "A_upper_levels": sum(key.species == "A" for key in terms),
        "E_upper_levels": sum(key.species == "E" for key in terms),
        "rows_outside_selected_term_cluster": conflict_rows,
        "maximum_duplicate_term_spread_cm1": maximum_spread,
        "cluster_tolerance_cm1": 0.005,
        "term_representative": "largest-A row in highest-summed-A cluster",
        "source_energy_reference": (
            "source lower identity validated by printed lower term; upper term "
            "formed as frozen-ground lower energy plus source wavenumber"
        ),
        "hitran_E_ground_zero_cm1": HITRAN_E_GROUND_ZERO_CM1,
        "lower_energy_tolerance_cm1": HITRAN_LOWER_ENERGY_TOLERANCE_CM1,
        "maximum_lower_energy_residual_cm1": max(lower_residuals),
        "E_signed_K_rows_corrected_from_lower_energy": e_sign_corrections,
    }


def _load_c13_model(root: Path) -> tuple[dict[str, dict[int, float]], dict[str, float], dict[str, float]]:
    payload = json.loads((root / ISOTOPE_CO_MODEL).read_text(encoding="utf-8"))
    origins = {
        species: {int(k): float(value) for k, value in payload["origins_cm1"][species].items()}
        for species in ("A", "E")
    }
    table = payload["table_IV_k0_power_series"]
    delta_b = {
        species: float(table[species]["B_co_cm1"]) - float(table[species]["B_ground_cm1"])
        for species in ("A", "E")
    }
    delta_d = {
        species: (
            float(table[species]["D_co_1e-6_cm1"])
            - float(table[species]["D_ground_1e-6_cm1"])
        ) * 1.0e-6
        for species in ("A", "E")
    }
    return origins, delta_b, delta_d


def _c13_co_term(
    ground: MasterLevel,
    model: tuple[dict[str, dict[int, float]], dict[str, float], dict[str, float]],
) -> tuple[float, str]:
    key = ground.key
    all_origins, delta_b, delta_d = model
    origins = all_origins[key.species]
    lookup = key.k
    if lookup in origins:
        origin = origins[lookup]
        status = "Mukhopadhyay1995_TableI_origin_plus_TableIV_K0_BD"
    else:
        origin = statistics.median(origins.values())
        status = "extrapolated_median_origin_outside_published_K_range_plus_TableIV_K0_BD"
    x = key.j * (key.j + 1)
    shift = origin + delta_b[key.species] * x - delta_d[key.species] * x * x
    return ground.energy_cm1 + shift, status


def _co_levels(
    selected: list[MasterLevel],
    *,
    isotopologue: str,
    main_terms: dict[LevelKey, float] | None,
    c13_model: tuple[dict[str, dict[int, float]], dict[str, float], dict[str, float]] | None,
) -> tuple[list[MasterLevel], Counter[str]]:
    ground = {row.key: row for row in selected if row.key.vt == 0}
    status_counts: Counter[str] = Counter()
    offsets: dict[str, list[float]] = defaultdict(list)
    if main_terms is not None:
        for key, term in main_terms.items():
            if key in ground:
                offsets[key.species].append(term - ground[key].energy_cm1)
    result: list[MasterLevel] = []
    for key, row in ground.items():
        if isotopologue == "12CH3OH":
            if main_terms is None:
                raise RuntimeError('Invariant failed in co_stretch.py: main_terms is not None')
            if key in main_terms:
                energy = main_terms[key]
                status = "Xu2004_HITRAN_format_direct_upper_term"
            else:
                energy = row.energy_cm1 + statistics.median(offsets[key.species])
                status = "extrapolated_species_median_vibrational_offset"
        else:
            if c13_model is None:
                raise COStretchError("carbon-13 CO-stretch source model is missing")
            energy, status = _c13_co_term(row, c13_model)
        status_counts[status] += 1
        result.append(
            MasterLevel(
                key=COStretchKey.from_ground(key),  # type: ignore[arg-type]
                energy_cm1=energy,
                weight=row.weight,
                source=status,
                source_index=key.encode(),
                mapping_confidence=(
                    "exact_source_quantum_identity"
                    if "direct" in status or "TableI" in status
                    else "source_model_extrapolation"
                ),
            )
        )
    return result, status_counts


def _copy_co_collisions(
    ground_rows: list[SelectedCollision],
    levels: list[MasterLevel],
    key_to_index: dict[object, int],
    temperatures: tuple[float, ...],
) -> list[SelectedCollision]:
    by_key = {row.key: row for row in levels}
    result: list[SelectedCollision] = []
    for ground in ground_rows:
        if ground.collider_id not in {2, 3, 6}:
            continue
        upper_key = COStretchKey.from_ground(ground.upper_key)
        lower_key = COStretchKey.from_ground(ground.lower_key)
        rates = ground.rates
        transform = ground.transform + "; vt0_exact_rotational_identity_vCO1_intraband_copy"
        if by_key[upper_key].energy_cm1 < by_key[lower_key].energy_cm1:
            source_upper = by_key[ground.upper_key]
            source_lower = by_key[ground.lower_key]
            gap = source_upper.energy_cm1 - source_lower.energy_cm1
            upper_key, lower_key = lower_key, upper_key
            rates = tuple(
                rate
                * (source_upper.weight / source_lower.weight)
                * (math.exp(-gap * 1.4387768775039338 / temperature) if gap else 1.0)
                for rate, temperature in zip(ground.rates, temperatures)
            )
            transform += "; detailed_balance_for_vCO1_energy_reversal"
        result.append(
            replace(
                ground,
                upper=key_to_index[upper_key],
                lower=key_to_index[lower_key],
                upper_key=upper_key,  # type: ignore[arg-type]
                lower_key=lower_key,  # type: ignore[arg-type]
                rates=rates,
                status="vt0-intraband-identity-copy-to-vCO1",
                transform=transform,
            )
        )
    return result


def _co_radiative(
    base_lines_path: Path,
    species: str,
    strengths: list[RotationalStrength],
    levels: list[MasterLevel],
    key_to_index: dict[object, int],
) -> tuple[list[RadiativeTransition], Counter[str]]:
    energy = {row.key: row.energy_cm1 for row in levels}
    weight = {row.key: row.weight for row in levels}
    rows: list[tuple[object, object, float, str]] = []
    provenance_counts: Counter[str] = Counter()
    for strength in strengths:
        if strength.left.species != species or strength.smu2_debye2 < PRODUCT_FLOOR_SMU2:
            continue
        orientations = [(strength.left, strength.right)]
        if strength.left != strength.right:
            orientations.append((strength.right, strength.left))
        for co_identity, ground_identity in orientations:
            co = COStretchKey.from_ground(co_identity)
            if co not in key_to_index or ground_identity not in key_to_index:
                continue
            if energy[co] <= energy[ground_identity]:
                raise COStretchError("CO-stretch emission endpoint is not above vt=0")
            rows.append((co, ground_identity, strength.smu2_debye2, "vCO1_to_vt0_0.27D_a_axis_BELGI"))
            provenance_counts["cross_vibrational"] += 1

    with base_lines_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["species"].upper() != species:
                continue
            upper = LevelKey.decode(row["upper_key"])
            lower = LevelKey.decode(row["lower_key"])
            if upper.vt != 0 or lower.vt != 0:
                continue
            smu2 = float(row["smu2_debye2"])
            if smu2 < PRODUCT_FLOOR_SMU2:
                continue
            co_upper = COStretchKey.from_ground(upper)
            co_lower = COStretchKey.from_ground(lower)
            if co_upper not in key_to_index or co_lower not in key_to_index:
                continue
            if energy[co_upper] < energy[co_lower]:
                co_upper, co_lower = co_lower, co_upper
            rows.append((co_upper, co_lower, smu2, "vCO1_intraband_vt0_permanent_operator_copy"))
            provenance_counts["co_intraband"] += 1

    result: list[RadiativeTransition] = []
    seen: set[frozenset[object]] = set()
    for upper, lower, smu2, _ in rows:
        identity = frozenset((upper, lower))
        if identity in seen:
            raise COStretchError("duplicate CO-stretch radiative endpoint pair")
        seen.add(identity)
        frequency_mhz = (energy[upper] - energy[lower]) * CM1_TO_MHZ
        result.append(
            RadiativeTransition(
                0,
                key_to_index[upper],
                key_to_index[lower],
                einstein_a_from_smu2(frequency_mhz, smu2, weight[upper]),
                frequency_mhz / 1000.0,
            )
        )
    return result, provenance_counts


def _regularize(levels: list[MasterLevel]) -> tuple[list[MasterLevel], int]:
    ordered = sorted(levels, key=lambda row: (row.energy_cm1, row.key.encode()))
    result: list[MasterLevel] = []
    changes = 0
    for row in ordered:
        if result and row.energy_cm1 <= result[-1].energy_cm1:
            if row.energy_cm1 != result[-1].energy_cm1:
                raise COStretchError("source levels cross under deterministic energy ordering")
            row = replace(row, energy_cm1=math.nextafter(result[-1].energy_cm1, math.inf))
            changes += 1
        result.append(row)
    return result, changes


def _build_one_isotopologue(
    root: Path,
    *,
    isotopologue: str,
    levels_path: Path,
    lines_path: Path,
    collision_rows: list[object],
    strengths: list[RotationalStrength],
    main_terms: dict[LevelKey, float] | None,
    c13_model: tuple[dict[str, dict[int, float]], dict[str, float], dict[str, float]] | None,
) -> dict[str, object]:
    policy = json.loads((root / "config/extended-build.json").read_text(encoding="utf-8"))
    temperatures = tuple(float(value) for value in policy["collision_temperature_grid_k"])
    all_levels = _load_extended_levels(levels_path)
    if isotopologue == "12CH3OH":
        mass = 32.04186
    else:
        from .isotopologues.ch3oh13.profile import MOLECULAR_MASS_U

        mass = MOLECULAR_MASS_U
    output_root = (
        root / "outputs/co-stretch"
        if isotopologue == "12CH3OH"
        else root / "isotopologues/13ch3oh/outputs/co-stretch"
    )
    if isotopologue == "12CH3OH":
        upstream_root = root / "outputs/extended"
        stem = "ch3oh"
    else:
        upstream_root = root / "isotopologues/13ch3oh/outputs/extended"
        stem = "13ch3oh"
    upstream_manifest_path = upstream_root / "manifest.json"
    if not upstream_manifest_path.is_file():
        raise COStretchError(
            f"missing upstream extended manifest {upstream_manifest_path.relative_to(root)}"
        )
    records: list[dict[str, object]] = []
    prepared: list[tuple[str, LamdaData]] = []
    for species in ("A", "E"):
        torsional, template, torsional_ties = _select_levels(
            all_levels,
            species,
            j_max=int(policy["j_max"]),
            k_max=int(policy["abs_k_max"]),
            ground_max_cm1=float(policy["ground_manifold_max_energy_cm1"]),
        )
        co, co_status = _co_levels(
            torsional,
            isotopologue=isotopologue,
            main_terms=main_terms,
            c13_model=c13_model,
        )
        combined, combined_ties = _regularize(torsional + co)
        key_to_index = {row.key: index for index, row in enumerate(combined, 1)}
        torsional_index = {
            key: index
            for key, index in key_to_index.items()
            if isinstance(key, LevelKey)
        }
        native_rows = _calculated_collisions(
            collision_rows,
            species,
            torsional_index,
            torsional,
            temperatures,
            allow_sqrt_endpoint_extrapolation=True,
        )
        ground_rows = [
            row
            for row in native_rows
            if row.upper_key.vt == row.lower_key.vt == 0
        ]
        copied_rows = _copy_vt0_intraband(
            ground_rows, torsional, key_to_index, temperatures
        )
        collisions, collision_precedence = _prefer_native_calculated_collisions(
            native_rows,
            copied_rows,
        )
        collisions.extend(
            _copy_co_collisions(ground_rows, combined, key_to_index, temperatures)
        )
        if isotopologue == "13CH3OH":
            from .isotopologues.ch3oh13.build import _partners_13ch3oh, _scale_calculated_collisions

            collisions = _scale_calculated_collisions(collisions)
            partners = _partners_13ch3oh(collisions, temperatures)
        else:
            partners = _partners(collisions, temperatures)
        radiative = _select_radiative(
            lines_path,
            species,
            key_to_index,
            minimum_smu2_debye2=PRODUCT_FLOOR_SMU2,
        )
        upstream_filename = f"{species.lower()}-{stem}-vt01234-intraband-copy"
        upstream_filename += (
            "-reduced-mass-proxy.dat"
            if isotopologue == "13CH3OH"
            else ".dat"
        )
        upstream_path = upstream_root / upstream_filename
        if not upstream_path.is_file():
            raise COStretchError(
                f"missing upstream extended product {upstream_path.relative_to(root)}"
            )
        _assert_manifest_binds_file(upstream_manifest_path, upstream_path)
        upstream = read_lamda(upstream_path)
        radiative, inheritance = _inherit_radiative_subset(
            radiative,
            combined,
            upstream,
        )
        co_radiative, radiative_counts = _co_radiative(
            lines_path, species, strengths, combined, key_to_index
        )
        radiative.extend(co_radiative)
        radiative.sort(key=lambda row: (row.upper, row.lower))
        radiative = [replace(row, index=index) for index, row in enumerate(radiative, 1)]
        data = LamdaData(
            molecule=(
                f"{species}-{isotopologue} vt=0..4 + vCO=1 "
                "[exploratory intraband collision copies]"
            ),
            mass_amu=mass,
            levels=[
                Level(index, row.energy_cm1, row.weight, row.key.encode())
                for index, row in enumerate(combined, 1)
            ],
            radiative=radiative,
            collision_partners=partners,
        )
        assert_collision_selection(
            data, collisions, temperatures, context=upstream_path.name
        )
        collision_inheritance = assert_parent_collision_inheritance(
            upstream, data, context=upstream_path.name
        )
        stem = "ch3oh" if isotopologue == "12CH3OH" else "13ch3oh"
        filename = f"{species.lower()}-{stem}-vt01234-vco1-intraband-copy.dat"
        radiative_audit = audit_radiative_integrity(
            data,
            context=filename,
            parent=upstream,
            minimum_smu2_debye2=PRODUCT_FLOOR_SMU2,
        )
        errors = validate_lamda(data)
        if errors:
            raise COStretchError(
                f"{species}-{isotopologue} CO-stretch validation failed: {errors[:10]}"
            )
        prepared.append((filename, data))
        records.append(
            {
                "file": filename,
                "species": species,
                "levels": len(combined),
                "levels_by_manifold": {
                    **{f"vt{vt}": sum(isinstance(row.key, LevelKey) and row.key.vt == vt for row in combined) for vt in range(5)},
                    "vCO1_vt0": len(co),
                },
                "rotational_template_size": len(template),
                "co_energy_status_counts": dict(sorted(co_status.items())),
                "radiative_transitions": len(radiative),
                "radiative_integrity": radiative_audit,
                "radiative_inheritance": {
                    **inheritance,
                    "upstream_path": str(upstream_path.relative_to(root)),
                    "upstream_sha256": sha256_file(upstream_path),
                    "upstream_manifest_path": str(
                        upstream_manifest_path.relative_to(root)
                    ),
                    "upstream_manifest_sha256": sha256_file(
                        upstream_manifest_path
                    ),
                },
                "collision_inheritance": {
                    **collision_inheritance,
                    "upstream_path": str(upstream_path.relative_to(root)),
                    "upstream_sha256": sha256_file(upstream_path),
                    "upstream_manifest_path": str(
                        upstream_manifest_path.relative_to(root)
                    ),
                    "upstream_manifest_sha256": sha256_file(
                        upstream_manifest_path
                    ),
                },
                "co_radiative_counts": dict(radiative_counts),
                "collision_rows": len(collisions),
                "collision_precedence": collision_precedence,
                "collision_rows_by_collider": dict(Counter(row.collider_id for row in collisions)),
                "energy_tie_breaks": torsional_ties,
                "combined_energy_tie_break_count": combined_ties,
                "structural_validation": "passed",
            }
        )
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "exploratory-source-bound-input-candidate",
        "isotopologue": isotopologue,
        "transition_operator": {
            "vibrational_transition_moment_debye": VCO_TRANSITION_MOMENT_DEBYE,
            "body_fixed_direction": "approximate CO bond / BELGI a axis",
            "rotational_wavefunctions": "ground-vibrational BELGI eigenfunctions reused",
            "nominal_delta_k_filter": False,
            "rank_rule": "BELGI INTEN; |Delta J| <= 1; J=0 to J=0 excluded",
        },
        "selection": {
            "operator_generation_minimum_smu2_debye2": GENERATION_FLOOR_SMU2,
            "product_minimum_smu2_debye2": PRODUCT_FLOOR_SMU2,
        },
        "collisions": {
            "torsional_intermanifold_rows": (
                "available He vt0-vt1, vt0-vt2, and vt1-vt2 rows retained"
            ),
            "vibrational_intermanifold_rows": 0,
            "torsional_intraband_precedence": (
                "available native vt1/vt2 rows before vt0 exact-identity copies"
            ),
            "vCO1_intraband": "vt0 exact rotational-identity copy",
        },
        "spectroscopy": (
            {
                "frequency_source": "Xu et al. 2004 V8 upper terms in HITRAN-2004-format mirror file",
                "source_path": str(MAIN_HITRAN_RELATIVE_PATH),
                "source_sha256": sha256_file(root / MAIN_HITRAN_RELATIVE_PATH),
                "source_metadata_path": str(MAIN_HITRAN_METADATA_RELATIVE_PATH),
                "source_metadata_sha256": sha256_file(
                    root / MAIN_HITRAN_METADATA_RELATIVE_PATH
                ),
                "uncovered_identity_policy": "labelled species-median vibrational-offset extrapolation",
            }
            if isotopologue == "12CH3OH"
            else {
                "frequency_source": "Mukhopadhyay et al. 1995 Tables I and IV transcription",
                "source_path": str(ISOTOPE_CO_MODEL),
                "source_sha256": sha256_file(root / ISOTOPE_CO_MODEL),
                "archival_context_path": str(ISOTOPE_THESIS_RELATIVE_PATH),
                "archival_context_sha256": sha256_file(root / ISOTOPE_THESIS_RELATIVE_PATH),
                "uncovered_identity_policy": "labelled median-origin extrapolation outside published K range",
                "precision_upgrade_gap": "AIP PAPS JCPSA-102-6444-32 is access-controlled",
            }
        ),
        "build_inputs_sha256": {
            str(levels_path.relative_to(root)): sha256_file(levels_path),
            str(lines_path.relative_to(root)): sha256_file(lines_path),
            "config/extended-build.json": sha256_file(root / "config/extended-build.json"),
            "fortran/belgi_frozen_strengths.f": sha256_file(root / "fortran/belgi_frozen_strengths.f"),
        },
        "files": records,
    }
    with tempfile.TemporaryDirectory(prefix=".co-stretch-stage-", dir=output_root.parent) as temporary:
        stage = Path(temporary)
        by_name = {row["file"]: row for row in records}
        for filename, data in prepared:
            path = stage / filename
            write_lamda(path, data)
            reread = read_lamda(path)
            independent = read_summary(path)
            if (
                len(reread.levels) != len(data.levels)
                or len(reread.radiative) != len(data.radiative)
                or independent.level_count != len(data.levels)
                or independent.radiative_count != len(data.radiative)
                or {row[0] for row in independent.collider_counts} != {2, 3, 6}
            ):
                raise COStretchError(f"{filename} failed write/read cross-check")
            by_name[filename]["sha256"] = sha256_file(path)
            by_name[filename]["validation"] = {
                "primary_parser": "pass",
                "independent_parser": "pass",
                "format_validation": "pass",
            }
        output_root.mkdir(parents=True, exist_ok=True)
        for path in stage.iterdir():
            os.replace(path, output_root / path.name)
        manifest_path = output_root / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return manifest


def build_co_stretch_products(
    project_root: str | Path,
    *,
    compiler: str | Path | None = None,
) -> dict[str, object]:
    """Build and validate main- and carbon-13 CO-stretch LAMDA pairs."""
    root = Path(project_root).resolve()
    collision_source_audit = audit_basecol_collision_master(root)
    required = [
        root / MAIN_HITRAN_RELATIVE_PATH,
        root / MAIN_HITRAN_METADATA_RELATIVE_PATH,
        root / ISOTOPE_THESIS_RELATIVE_PATH,
        root / MAIN_LEVELS,
        root / MAIN_LINES,
        root / ISOTOPE_LEVELS,
        root / ISOTOPE_LINES,
        root / ISOTOPE_CONSTANTS,
        root / ISOTOPE_CO_MODEL,
    ]
    missing = [str(path.relative_to(root)) for path in required if not path.is_file()]
    if missing:
        raise COStretchError("missing CO-stretch build inputs: " + ", ".join(missing))
    strengths_main, strength_summary_main = _run_transition_operator(
        root,
        root / "data/raw/local-cache/xu2008-final-constants.txt",
        compiler=compiler,
    )
    strengths_c13, strength_summary_c13 = _run_transition_operator(
        root,
        root / ISOTOPE_CONSTANTS,
        compiler=compiler,
    )
    terms, term_audit = _main_co_terms(
        root / MAIN_HITRAN_RELATIVE_PATH,
        root / MAIN_LEVELS,
    )
    _, _, main_collisions = load_master(root)
    from .isotopologues.ch3oh13.build import _load_main_collision_master

    main = _build_one_isotopologue(
        root,
        isotopologue="12CH3OH",
        levels_path=root / MAIN_LEVELS,
        lines_path=root / MAIN_LINES,
        collision_rows=main_collisions,
        strengths=strengths_main,
        main_terms=terms,
        c13_model=None,
    )
    carbon13 = _build_one_isotopologue(
        root,
        isotopologue="13CH3OH",
        levels_path=root / ISOTOPE_LEVELS,
        lines_path=root / ISOTOPE_LINES,
        collision_rows=_load_main_collision_master(root),
        strengths=strengths_c13,
        main_terms=None,
        c13_model=_load_c13_model(root),
    )
    result = {
        "schema_version": 1,
        "main": main,
        "carbon13": carbon13,
        "operator_enumeration": {
            "main": strength_summary_main,
            "carbon13": strength_summary_c13,
        },
        "main_frequency_source_audit": term_audit,
        "collision_source_level_audit": collision_source_audit,
        "source_sha256": {
            str(MAIN_HITRAN_RELATIVE_PATH): sha256_file(root / MAIN_HITRAN_RELATIVE_PATH),
            str(MAIN_HITRAN_METADATA_RELATIVE_PATH): sha256_file(
                root / MAIN_HITRAN_METADATA_RELATIVE_PATH
            ),
            str(ISOTOPE_THESIS_RELATIVE_PATH): sha256_file(root / ISOTOPE_THESIS_RELATIVE_PATH),
            str(ISOTOPE_CO_MODEL): sha256_file(root / ISOTOPE_CO_MODEL),
        },
    }
    (root / "outputs/co-stretch/build-summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


__all__ = ["COStretchError", "COStretchKey", "build_co_stretch_products"]
