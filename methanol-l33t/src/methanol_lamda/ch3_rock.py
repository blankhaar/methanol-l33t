"""Build uncoupled ``vt=0..4 + vCO=1 + vrock=1`` methanol inputs.

This is an excitation-model input, not a coupled vibrational spectroscopy
fit.  Direct V7 data are retained where available; missing level origins are
explicit source-model extrapolations.  The carbon-13 radiative strengths copy
the main-isotopologue line strengths by exact rotational identity and rescale
Einstein A with frequency cubed.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import tempfile

from .co_stretch import (
    CM1_TO_MHZ,
    HITRAN_E_GROUND_ZERO_CM1,
    MAIN_HITRAN_RELATIVE_PATH,
    MAIN_HITRAN_METADATA_RELATIVE_PATH,
    MAIN_LEVELS,
    PRODUCT_FLOOR_SMU2,
    _main_ground_energies,
    _parse_hitran_local,
    _resolve_hitran_ground_identity,
    _resolve_hitran_upper_identity,
)
from .collision_integrity import assert_parent_collision_inheritance
from .hamiltonian import sha256_file
from .extended import _assert_manifest_binds_file
from .independent import read_summary
from .lamda import read_lamda, validate_lamda, write_lamda
from .models import (
    CollisionPartner,
    CollisionTransition,
    LamdaData,
    Level,
    RadiativeTransition,
)
from .physics import CM1_TO_K, SMU2_TO_A, einstein_a_from_smu2
from .quantum import LevelKey
from .radiative_integrity import audit_radiative_integrity


C13_ROCK_MODEL = Path(
    "isotopologues/13ch3oh/data/raw/transcribed/"
    "predoi-cross1997-ch3-rock-origins.json"
)
C13_THESIS = Path("data/raw/official/library-and-archives-canada/nq23873.pdf")
ROCK_TERM_CLUSTER_TOLERANCE_CM1 = 0.005
V7_LABEL_CORRECTION_TOLERANCE_CM1 = 0.001


class CH3RockError(RuntimeError):
    """Raised when the uncoupled CH3-rock product is inconsistent."""


@dataclass(frozen=True, slots=True, order=True)
class RockKey:
    species: str
    j: int
    k: int
    parity: str = ""

    def __post_init__(self) -> None:
        ground = LevelKey(self.species, 0, self.j, self.k, self.parity)
        object.__setattr__(self, "species", ground.species)

    @classmethod
    def from_ground(cls, key: LevelKey) -> "RockKey":
        if key.vt != 0:
            raise CH3RockError("CH3-rock identity must be copied from vt=0")
        return cls(key.species, key.j, key.k, key.parity)

    def encode(self) -> str:
        if self.species == "A":
            return f"A:vrock1:vt0:J{self.j}:K{self.k}:{self.parity}"
        return f"E:vrock1:vt0:J{self.j}:K{self.k:+d}"


@dataclass(frozen=True, slots=True)
class SourceLine:
    upper: LevelKey
    lower: LevelKey
    source_upper_j: int
    source_upper_weight: float
    source_lower_weight: float
    wavenumber_cm1: float
    einstein_a_s: float
    source_line: int
    inferred_upper_energy_cm1: float
    identity_basis: str
    raw_record_sha256: str
    lower_mapping_residual_cm1: float
    label_correction_anchor_residual_cm1: float | None

    @property
    def smu2_debye2(self) -> float:
        # The 26 corrected V7 assignments carry the source statistical weight
        # for the printed (incorrect) upper J.  Recover the physical line
        # strength from that source convention first, then let the product
        # recompute A with the corrected upper-state weight.  Otherwise the
        # label repair would silently change S mu^2 by (2J+1).
        frequency_mhz = self.wavenumber_cm1 * CM1_TO_MHZ
        return self.einstein_a_s * (2.0 * self.source_upper_j + 1.0) / (
            SMU2_TO_A * frequency_mhz**3
        )


def _rank1_allowed(upper: LevelKey, lower: LevelKey) -> bool:
    """Return whether rotational endpoints can couple through a rank-1 dipole."""

    return abs(upper.j - lower.j) <= 1 and not (
        upper.j == 0 and lower.j == 0
    )


def _main_rock_source(
    path: Path,
    ground_levels_path: Path,
) -> tuple[dict[LevelKey, float], list[SourceLine], dict[str, object]]:
    ground_energies = _main_ground_energies(ground_levels_path)
    candidates: dict[LevelKey, list[tuple[float, float, int]]] = defaultdict(list)
    candidate_lines: dict[int, SourceLine] = {}
    source_rows_total = 0
    source_rows_label_corrected = 0
    source_rows_outside_model = 0
    e_sign_corrections = 0
    lower_residuals: list[float] = []
    correction_anchor_residuals: list[float] = []
    source_lines = [
        (number, line)
        for number, line in enumerate(
            path.read_text(encoding="ascii").splitlines(), 1
        )
        if line
        and not line.startswith("#")
        and line[67:82].strip() == "V7"
        and line[82:97].strip() == "GROUND"
    ]
    valid_term_anchors: dict[LevelKey, list[float]] = defaultdict(list)
    for _, line in source_lines:
        upper = _parse_hitran_local(line[97:112])
        lower = _parse_hitran_local(line[112:127])
        if _rank1_allowed(upper, lower):
            lower_relative = float(line[45:55]) - (
                HITRAN_E_GROUND_ZERO_CM1 if lower.species == "E" else 0.0
            )
            valid_term_anchors[upper].append(
                lower_relative + float(line[3:15])
            )

    for number, line in source_lines:
        raw_upper = _parse_hitran_local(line[97:112])
        raw_lower = _parse_hitran_local(line[112:127])
        source_upper_j = raw_upper.j
        source_rows_total += 1
        if raw_upper.species != raw_lower.species:
            raise CH3RockError("V7 source row changes A/E species")
        wavenumber = float(line[3:15])
        einstein_a = float(line[25:35])
        lower_energy = float(line[45:55])
        source_upper_weight = float(line[146:153])
        source_lower_weight = float(line[153:160])
        if not math.isclose(
            source_upper_weight,
            8.0 * (2.0 * source_upper_j + 1.0),
            rel_tol=0.0,
            abs_tol=1.0e-9,
        ) or not math.isclose(
            source_lower_weight,
            8.0 * (2.0 * raw_lower.j + 1.0),
            rel_tol=0.0,
            abs_tol=1.0e-9,
        ):
            raise CH3RockError(
                f"V7 source row {number} has an unexpected statistical weight"
            )
        if not math.isfinite(einstein_a) or einstein_a <= 0.0:
            raise CH3RockError(
                f"V7 source row {number} has no positive Einstein A for S mu^2"
            )
        label_corrected = False
        anchor_residual: float | None = None
        if not _rank1_allowed(raw_upper, raw_lower):
            proposed = (
                LevelKey("E", 0, raw_upper.j - 1, raw_upper.k)
                if raw_upper.species == "E"
                and raw_lower.species == "E"
                and raw_upper.j - raw_lower.j == 2
                and raw_upper.k == -5
                and raw_lower.k == -6
                else None
            )
            source_term = (
                lower_energy - HITRAN_E_GROUND_ZERO_CM1 + wavenumber
            )
            anchor_residual = (
                math.inf
                if proposed is None or proposed not in valid_term_anchors
                else min(
                    abs(anchor - source_term)
                    for anchor in valid_term_anchors[proposed]
                )
            )
            if anchor_residual > V7_LABEL_CORRECTION_TOLERANCE_CM1:
                raise CH3RockError(
                    f"V7 source row {number} has an unresolved non-rank-1 identity"
                )
            if proposed is None:
                raise RuntimeError('Invariant failed in ch3_rock.py: proposed is not None')
            raw_upper = proposed
            label_corrected = True
            source_rows_label_corrected += 1
            correction_anchor_residuals.append(anchor_residual)
        lower, residual, corrected = _resolve_hitran_ground_identity(
            raw_lower, lower_energy, ground_energies
        )
        if lower is None:
            source_rows_outside_model += 1
            continue
        identity_basis = "frozen_ground_lower_energy_crosswalk"
        if label_corrected:
            identity_basis += "; source_upper_J_minus_one_energy_validated"
        upper = _resolve_hitran_upper_identity(raw_upper, raw_lower, lower)
        if upper not in ground_energies:
            source_rows_outside_model += 1
            continue
        if not _rank1_allowed(upper, lower):
            raise CH3RockError(
                f"V7 source row {number} violates rank-1 after identity resolution"
            )
        if residual is not None:
            lower_residuals.append(residual)
        e_sign_corrections += int(corrected)
        frequency_mhz = wavenumber * CM1_TO_MHZ
        score = (
            einstein_a * (2.0 * source_upper_j + 1.0)
            / (SMU2_TO_A * frequency_mhz**3)
        )
        inferred_upper = lower_energy - (
            HITRAN_E_GROUND_ZERO_CM1 if lower.species == "E" else 0.0
        ) + wavenumber
        candidates[upper].append((inferred_upper, score, number))
        candidate_lines[number] = SourceLine(
            upper,
            lower,
            source_upper_j,
            source_upper_weight,
            source_lower_weight,
            wavenumber,
            einstein_a,
            number,
            inferred_upper,
            identity_basis,
            hashlib.sha256(line.encode("ascii")).hexdigest(),
            float(residual),
            anchor_residual,
        )
    if not candidates or not candidate_lines:
        raise CH3RockError("no V7 -> GROUND rows found")
    terms: dict[LevelKey, float] = {}
    selected_source_rows: set[int] = set()
    rejected = 0
    for key, values in candidates.items():
        clusters: list[list[tuple[float, float, int]]] = []
        for value in sorted(values):
            if (
                not clusters
                or value[0] - clusters[-1][0][0]
                > ROCK_TERM_CLUSTER_TOLERANCE_CM1
            ):
                clusters.append([value])
            else:
                clusters[-1].append(value)
        chosen = max(
            clusters,
            key=lambda group: (
                sum(x[1] for x in group),
                len(group),
            ),
        )
        terms[key] = max(chosen, key=lambda x: x[1])[0]
        within_representative = [
            value
            for value in chosen
            if abs(value[0] - terms[key]) <= ROCK_TERM_CLUSTER_TOLERANCE_CM1
        ]
        selected_source_rows.update(value[2] for value in within_representative)
        rejected += len(values) - len(within_representative)
    direct = [
        candidate_lines[number]
        for number in sorted(selected_source_rows)
        if number in candidate_lines
    ]
    corrected_direct = [
        row
        for row in direct
        if "source_upper_J_minus_one" in row.identity_basis
    ]
    return terms, direct, {
        "source_rows_mapped_before_term_cluster_selection": sum(
            len(rows) for rows in candidates.values()
        ),
        "source_rows_retained_in_selected_term_clusters": len(direct),
        "source_rows_total": source_rows_total,
        "source_rows_outside_frozen_ground_model": source_rows_outside_model,
        "source_rows_upper_J_minus_one_energy_validated": (
            source_rows_label_corrected
        ),
        "source_label_correction_pattern": (
            "E2 K=5<-6, printed DeltaJ=2; adopt upper J-1 only when "
            "the inferred upper term matches independent rank-1 anchors"
        ),
        "source_label_correction_tolerance_cm1": (
            V7_LABEL_CORRECTION_TOLERANCE_CM1
        ),
        "source_label_correction_maximum_anchor_residual_cm1": max(
            correction_anchor_residuals
        ),
        "corrected_rows_retained_in_selected_term_clusters": [
            row.source_line for row in corrected_direct
        ],
        "corrected_rows_above_product_strength_floor": [
            row.source_line
            for row in corrected_direct
            if row.smu2_debye2 >= PRODUCT_FLOOR_SMU2
        ],
        "unique_upper_levels": len(terms),
        "A_upper_levels": sum(key.species == "A" for key in terms),
        "E_upper_levels": sum(key.species == "E" for key in terms),
        "rows_outside_selected_term_cluster": rejected,
        "cluster_tolerance_cm1": ROCK_TERM_CLUSTER_TOLERANCE_CM1,
        "cluster_definition": "maximum diameter, not single-link adjacency",
        "term_representative": (
            "largest-Smu2 row in the highest-summed-Smu2 source-term cluster"
        ),
        "source_energy_reference": (
            "source upper term = printed lower term minus the A/E species "
            "zero plus source wavenumber; in-model lower identities are "
            "validated against the frozen Hamiltonian"
        ),
        "hitran_E_ground_zero_cm1": HITRAN_E_GROUND_ZERO_CM1,
        "maximum_lower_energy_residual_cm1": max(lower_residuals),
        "E_signed_K_rows_corrected_from_lower_energy": e_sign_corrections,
    }


def _c13_origins(path: Path) -> tuple[dict[str, dict[int, float]], dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    origins = {
        species: {int(k): float(value) for k, value in payload["origins_cm1"][species].items()}
        for species in ("A", "E")
    }
    return origins, payload["source"]


def _ground_key(label: str) -> LevelKey | None:
    try:
        key = LevelKey.decode(label)
    except ValueError:
        return None
    return key if key.vt == 0 else None


def _rock_levels(
    base: LamdaData,
    *,
    isotopologue: str,
    main_terms: dict[LevelKey, float],
    c13_origins: dict[str, dict[int, float]],
) -> tuple[list[Level], Counter[str]]:
    offsets: dict[str, list[float]] = defaultdict(list)
    ground: dict[LevelKey, Level] = {}
    for level in base.levels:
        key = _ground_key(level.quantum_label)
        if key is not None:
            ground[key] = level
            if key in main_terms:
                offsets[key.species].append(main_terms[key] - level.energy_cm1)
    counts: Counter[str] = Counter()
    result: list[Level] = []
    for key, level in ground.items():
        if isotopologue == "12CH3OH":
            if key in main_terms:
                energy = main_terms[key]
                status = "Xu2004_V7_direct_upper_term"
            else:
                energy = level.energy_cm1 + statistics.median(offsets[key.species])
                status = "extrapolated_species_median_vibrational_offset"
        else:
            lookup = abs(key.k) if key.species == "A" else key.k
            if lookup in c13_origins[key.species]:
                origin = c13_origins[key.species][lookup]
                status = "PredoiCross1997_Table4.1_observed_subband_origin"
            else:
                origin = statistics.median(c13_origins[key.species].values())
                status = "extrapolated_species_median_origin_outside_published_K"
            energy = level.energy_cm1 + origin
        counts[status] += 1
        result.append(Level(0, energy, level.weight, RockKey.from_ground(key).encode()))
    return result, counts


def _merge_levels(base: LamdaData, rock: list[Level]) -> tuple[list[Level], dict[int, int], dict[str, int], int]:
    tagged = [(level.energy_cm1, level.quantum_label, level, level.index) for level in base.levels]
    tagged.extend((level.energy_cm1, level.quantum_label, level, None) for level in rock)
    tagged.sort(key=lambda row: (row[0], row[1]))
    levels: list[Level] = []
    old_to_new: dict[int, int] = {}
    label_to_new: dict[str, int] = {}
    ties = 0
    for _, _, level, old_index in tagged:
        energy = level.energy_cm1
        if levels and energy <= levels[-1].energy_cm1:
            if energy != levels[-1].energy_cm1:
                raise CH3RockError("levels cross under deterministic ordering")
            energy = math.nextafter(levels[-1].energy_cm1, math.inf)
            ties += 1
        index = len(levels) + 1
        levels.append(Level(index, energy, level.weight, level.quantum_label))
        label_to_new[level.quantum_label] = index
        if old_index is not None:
            old_to_new[old_index] = index
    return levels, old_to_new, label_to_new, ties


def _radiative(
    base: LamdaData,
    levels: list[Level],
    old_to_new: dict[int, int],
    label_to_new: dict[str, int],
    direct: list[SourceLine],
    *,
    species: str,
    isotopologue: str,
) -> tuple[list[RadiativeTransition], Counter[str]]:
    energy = {row.index: row.energy_cm1 for row in levels}
    weight = {row.index: row.weight for row in levels}
    old_level = {row.index: row for row in base.levels}
    rows: dict[tuple[int, int], tuple[float, float, str]] = {}
    for row in base.radiative:
        upper = old_to_new[row.upper]
        lower = old_to_new[row.lower]
        rows[(upper, lower)] = (
            row.a_s,
            row.frequency_ghz,
            "base_A_and_frequency_exact_copy",
        )

    for row in base.radiative:
        upper_key = _ground_key(old_level[row.upper].quantum_label)
        lower_key = _ground_key(old_level[row.lower].quantum_label)
        if upper_key is None or lower_key is None:
            continue
        rock_upper = label_to_new[RockKey.from_ground(upper_key).encode()]
        rock_lower = label_to_new[RockKey.from_ground(lower_key).encode()]
        if energy[rock_upper] < energy[rock_lower]:
            rock_upper, rock_lower = rock_lower, rock_upper
        frequency_mhz = (energy[rock_upper] - energy[rock_lower]) * CM1_TO_MHZ
        source_smu2 = row.a_s * old_level[row.upper].weight / (
            SMU2_TO_A * (row.frequency_ghz * 1.0e3) ** 3
        )
        rows[(rock_upper, rock_lower)] = (
            einstein_a_from_smu2(frequency_mhz, source_smu2, weight[rock_upper]),
            frequency_mhz / 1.0e3,
            "rock_intraband_vt0_line_strength_copy",
        )

    direct_endpoints: set[tuple[int, int]] = set()
    for source in direct:
        if (
            source.upper.species != species
            or not _rank1_allowed(source.upper, source.lower)
            or source.einstein_a_s <= 0.0
            or source.smu2_debye2 < PRODUCT_FLOOR_SMU2
        ):
            continue
        rock_label = RockKey.from_ground(source.upper).encode()
        if rock_label not in label_to_new or source.lower.encode() not in label_to_new:
            raise CH3RockError(
                f"above-floor V7 source row {source.source_line} has a missing "
                "selected product endpoint"
            )
        upper = label_to_new[rock_label]
        lower = label_to_new[source.lower.encode()]
        if energy[upper] <= energy[lower]:
            raise CH3RockError(
                f"V7 source row {source.source_line} has a nonpositive product gap"
            )
        if (
            isotopologue == "12CH3OH"
            and abs(energy[upper] - source.inferred_upper_energy_cm1)
            > ROCK_TERM_CLUSTER_TOLERANCE_CM1 + 1.0e-12
        ):
            raise CH3RockError(
                f"V7 source row {source.source_line} is outside its adopted "
                "upper-term cluster"
            )
        frequency_mhz = (energy[upper] - energy[lower]) * CM1_TO_MHZ
        if isotopologue == "12CH3OH":
            a_value = einstein_a_from_smu2(
                frequency_mhz,
                source.smu2_debye2,
                weight[upper],
            )
            provenance = "Xu2004_V7_direct_line_strength_frequency_rescaled"
        else:
            a_value = einstein_a_from_smu2(
                frequency_mhz, source.smu2_debye2, weight[upper]
            )
            provenance = "main_V7_line_strength_identity_copy_to_13CH3OH"
        endpoint = (upper, lower)
        if endpoint in direct_endpoints or endpoint in rows:
            raise CH3RockError(
                f"duplicate V7 product endpoint from source row {source.source_line}"
            )
        direct_endpoints.add(endpoint)
        rows[endpoint] = (a_value, frequency_mhz / 1.0e3, provenance)

    counts = Counter(status for _, _, status in rows.values())
    result = []
    for index, ((upper, lower), (a_value, frequency_ghz, _)) in enumerate(
        sorted(rows.items()), 1
    ):
        result.append(RadiativeTransition(index, upper, lower, a_value, frequency_ghz))
    return result, counts


def _collisions(
    base: LamdaData,
    levels: list[Level],
    old_to_new: dict[int, int],
    label_to_new: dict[str, int],
) -> list[CollisionPartner]:
    energy = {row.index: row.energy_cm1 for row in levels}
    old_level = {row.index: row for row in base.levels}
    result: list[CollisionPartner] = []
    for partner in base.collision_partners:
        rows: list[CollisionTransition] = [
            CollisionTransition(
                0,
                old_to_new[row.upper],
                old_to_new[row.lower],
                row.rates_cm3_s,
            )
            for row in partner.transitions
        ]
        for row in partner.transitions:
            upper_key = _ground_key(old_level[row.upper].quantum_label)
            lower_key = _ground_key(old_level[row.lower].quantum_label)
            if upper_key is None or lower_key is None:
                continue
            upper = label_to_new[RockKey.from_ground(upper_key).encode()]
            lower = label_to_new[RockKey.from_ground(lower_key).encode()]
            rates = row.rates_cm3_s
            if energy[upper] < energy[lower]:
                upper, lower = lower, upper
                gap = old_level[row.upper].energy_cm1 - old_level[row.lower].energy_cm1
                rates = tuple(
                    rate
                    * (old_level[row.upper].weight / old_level[row.lower].weight)
                    * math.exp(-gap * CM1_TO_K / temperature)
                    for rate, temperature in zip(rates, partner.temperatures_k)
                )
            rows.append(CollisionTransition(0, upper, lower, rates))
        rows.sort(key=lambda row: (row.upper, row.lower))
        rows = [replace(row, index=index) for index, row in enumerate(rows, 1)]
        result.append(
            CollisionPartner(
                partner.collider_id,
                partner.description,
                partner.temperatures_k,
                rows,
            )
        )
    return result


def _build_one(
    root: Path,
    base_path: Path,
    base_manifest_path: Path,
    output_root: Path,
    *,
    species: str,
    isotopologue: str,
    main_terms: dict[LevelKey, float],
    direct: list[SourceLine],
    c13_origins: dict[str, dict[int, float]],
) -> dict[str, object]:
    if not base_manifest_path.is_file():
        raise CH3RockError(
            f"missing prerequisite manifest {base_manifest_path.relative_to(root)}"
        )
    _assert_manifest_binds_file(base_manifest_path, base_path)
    base = read_lamda(base_path)
    rock, status = _rock_levels(
        base,
        isotopologue=isotopologue,
        main_terms=main_terms,
        c13_origins=c13_origins,
    )
    levels, old_to_new, label_to_new, ties = _merge_levels(base, rock)
    radiative, radiative_counts = _radiative(
        base,
        levels,
        old_to_new,
        label_to_new,
        direct,
        species=species,
        isotopologue=isotopologue,
    )
    if radiative_counts["base_A_and_frequency_exact_copy"] != len(base.radiative):
        raise CH3RockError(
            f"{species}-{isotopologue} did not inherit every base radiative row"
        )
    weight = {row.index: row.weight for row in levels}
    invalid_rank1 = [
        row.index
        for row in radiative
        if abs(weight[row.upper] - weight[row.lower]) > 2.0
        or (weight[row.upper] == 1.0 and weight[row.lower] == 1.0)
    ]
    if invalid_rank1:
        raise CH3RockError(
            f"{species}-{isotopologue} contains non-rank-1 radiative rows: "
            f"{invalid_rank1[:10]}"
        )
    partners = _collisions(base, levels, old_to_new, label_to_new)
    data = LamdaData(
        f"{species}-{isotopologue} vt=0..4 + vCO=1 + vrock=1 [uncoupled exploratory]",
        base.mass_amu,
        levels,
        radiative,
        partners,
    )
    collision_inheritance = assert_parent_collision_inheritance(
        base, data, context=base_path.name
    )
    stem = "ch3oh" if isotopologue == "12CH3OH" else "13ch3oh"
    filename = f"{species.lower()}-{stem}-vt01234-vco1-vrock1-uncoupled.dat"
    radiative_audit = audit_radiative_integrity(
        data,
        context=filename,
        parent=base,
        minimum_smu2_debye2=PRODUCT_FLOOR_SMU2,
    )
    errors = validate_lamda(data)
    if errors:
        raise CH3RockError(f"{species}-{isotopologue} validation failed: {errors[:10]}")
    output_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ch3-rock-stage-", dir=output_root.parent) as temporary:
        staged = Path(temporary) / filename
        write_lamda(staged, data)
        reread = read_lamda(staged)
        independent = read_summary(staged)
        if (
            len(reread.levels) != len(levels)
            or len(reread.radiative) != len(radiative)
            or independent.level_count != len(levels)
            or independent.radiative_count != len(radiative)
        ):
            raise CH3RockError(f"{filename} failed dual-reader validation")
        digest = sha256_file(staged)
        os.replace(staged, output_root / filename)
    return {
        "file": filename,
        "sha256": digest,
        "species": species,
        "levels": len(levels),
        "rock_levels": len(rock),
        "rock_energy_status_counts": dict(sorted(status.items())),
        "radiative_transitions": len(radiative),
        "radiative_integrity": radiative_audit,
        "radiative_provenance_counts": dict(sorted(radiative_counts.items())),
        "radiative_inheritance": {
            "upstream_path": str(base_path.relative_to(root)),
            "upstream_sha256": sha256_file(base_path),
            "upstream_manifest_path": str(base_manifest_path.relative_to(root)),
            "upstream_manifest_sha256": sha256_file(base_manifest_path),
            "inherited_radiative_rows": len(base.radiative),
            "endpoint_policy": "exact ordered level identity",
            "inherited_fields": ["a_s", "frequency_ghz"],
        },
        "collision_inheritance": {
            **collision_inheritance,
            "upstream_path": str(base_path.relative_to(root)),
            "upstream_sha256": sha256_file(base_path),
            "upstream_manifest_path": str(base_manifest_path.relative_to(root)),
            "upstream_manifest_sha256": sha256_file(base_manifest_path),
        },
        "collision_rows_by_partner": {
            str(row.collider_id): len(row.transitions) for row in partners
        },
        "energy_tie_breaks": ties,
        "validation": "primary-independent-structural-pass",
    }


def build_ch3_rock_products(project_root: str | Path) -> dict[str, object]:
    """Build all four independent CH3-rock exploratory products."""
    root = Path(project_root).resolve()
    source_path = root / MAIN_HITRAN_RELATIVE_PATH
    c13_path = root / C13_ROCK_MODEL
    required = [source_path, c13_path, root / C13_THESIS, root / MAIN_LEVELS]
    required.append(root / MAIN_HITRAN_METADATA_RELATIVE_PATH)
    missing = [str(path.relative_to(root)) for path in required if not path.is_file()]
    if missing:
        raise CH3RockError("missing CH3-rock inputs: " + ", ".join(missing))
    main_terms, direct, source_audit = _main_rock_source(
        source_path,
        root / MAIN_LEVELS,
    )
    c13_origins, c13_source = _c13_origins(c13_path)
    records: list[dict[str, object]] = []
    for isotopologue, base_root, output_root in (
        ("12CH3OH", root / "outputs/co-stretch", root / "outputs/vibrational-full"),
        (
            "13CH3OH",
            root / "isotopologues/13ch3oh/outputs/co-stretch",
            root / "isotopologues/13ch3oh/outputs/vibrational-full",
        ),
    ):
        base_manifest_path = base_root / "manifest.json"
        stem = "ch3oh" if isotopologue == "12CH3OH" else "13ch3oh"
        for species in ("A", "E"):
            base = base_root / f"{species.lower()}-{stem}-vt01234-vco1-intraband-copy.dat"
            if not base.is_file():
                raise CH3RockError(
                    f"missing prerequisite {base.relative_to(root)}; run build-co-stretch first"
                )
            records.append(
                _build_one(
                    root,
                    base,
                    base_manifest_path,
                    output_root,
                    species=species,
                    isotopologue=isotopologue,
                    main_terms=main_terms,
                    direct=direct,
                    c13_origins=c13_origins,
                )
            )
    manifest = {
        "schema_version": 1,
        "status": "exploratory-source-bound-uncoupled-input",
        "model_scope": {
            "stretch_rock_mixing": "neglected by owner-reviewed excitation-model approximation",
            "rock_rotational_functions": "ground-state rotational identities",
            "rock_intraband_collisions": "vt0 exact-identity copy",
            "torsional_interband_collisions": (
                "available He vt0-vt1, vt0-vt2, and vt1-vt2 rows inherited "
                "from the CO-tier prerequisite"
            ),
            "stretch_or_rock_interband_collisions": 0,
            "carbon13_cross_band_strengths": "main V7 line-strength copy by exact rotational endpoints",
        },
        "main_source_audit": source_audit,
        "carbon13_source": c13_source,
        "source_sha256": {
            str(MAIN_HITRAN_RELATIVE_PATH): sha256_file(source_path),
            str(MAIN_HITRAN_METADATA_RELATIVE_PATH): sha256_file(
                root / MAIN_HITRAN_METADATA_RELATIVE_PATH
            ),
            str(C13_ROCK_MODEL): sha256_file(c13_path),
            str(C13_THESIS): sha256_file(root / C13_THESIS),
        },
        "files": records,
    }
    for output_root in (
        root / "outputs/vibrational-full",
        root / "isotopologues/13ch3oh/outputs/vibrational-full",
    ):
        own = [
            row for row in records
            if ("13ch3oh" in row["file"]) == ("isotopologues" in str(output_root))
        ]
        payload = {**manifest, "files": own}
        (output_root / "manifest.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return manifest


__all__ = [
    "C13_ROCK_MODEL",
    "CH3RockError",
    "RockKey",
    "build_ch3_rock_products",
]
