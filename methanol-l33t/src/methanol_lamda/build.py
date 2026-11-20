"""Construct balanced vt=0,1,2 LAMDA products from normalized master tables.

This module is intentionally strict.  It can build an explicitly named
``radiative-incomplete`` fallback, but it cannot turn a missing source into a
zero transition or a generic selection-rule line strength.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
import csv
import itertools
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Iterable

from .collision_integrity import (
    assert_collision_selection,
    audit_basecol_collision_master,
)
from .lamda import write_lamda
from .integrity import (
    RawBindingError,
    sha256_path,
    verify_dataset_record,
    verify_pinned_sha256,
)
from .models import CollisionPartner, CollisionTransition, LamdaData, Level, RadiativeTransition
from .physics import (
    CM1_TO_K,
    interpolate_or_sqrt_extrapolate_rate,
    interpolate_rate,
    upward_from_downward,
)
from .quantum import LevelKey
from .radiative_integrity import audit_radiative_integrity


COMMON_TEMPERATURES = (10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0)
DEFAULT_MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE = 1.0e-2
PRODUCTS = ("calculated", "rf2011-proxy", "cragg2005")
MEKHTIEV_SOURCE = "mekhtiev1999-linestrength-tables"
CDMS_SOURCE = "cdms-032504-v3"
XU_FROZEN_LEVEL_SOURCE = "local-belgi-xu2008-frozen-hamiltonian"
XU_CONSTANT_DIPOLE_SOURCE = "local-belgi-xu2008-constant-dipole"
DAGDIGIAN_SOURCE = "dagdigian2024-h2"
RF_H2_SOURCE = "rabli-flower2010-h2-vt012"
RF_HE_SOURCE = "rabli-flower2010-he-vt012"
RF_TORSION_HE_SOURCE = "rabli-flower2011-he-torsion-changing"


class MissingSourcesError(RuntimeError):
    """Raised when a production build lacks essential normalized inputs."""

    def __init__(self, missing: Iterable[str]) -> None:
        self.missing = tuple(missing)
        super().__init__("missing essential sources: " + ", ".join(self.missing))


class MappingError(ValueError):
    """Raised when mapped source endpoints contradict their physical roles."""


@dataclass(frozen=True, slots=True)
class MasterLevel:
    key: LevelKey
    energy_cm1: float
    weight: float
    source: str
    source_index: str
    mapping_confidence: str
    raw_sha256: str = ""


@dataclass(frozen=True, slots=True)
class MasterRadiative:
    species: str
    upper_key: LevelKey
    lower_key: LevelKey
    a_s: float
    frequency_ghz: float
    source: str
    status: str
    frequency_source: str
    raw_sha256: str = ""
    source_upper_id: str = ""
    source_lower_id: str = ""
    mapping_confidence: str = ""
    smu2_debye2: float = math.nan


@dataclass(frozen=True, slots=True)
class MasterCollision:
    species: str
    collider_id: int
    upper_key: LevelKey
    lower_key: LevelKey
    temperatures: tuple[float, ...]
    rates: tuple[float, ...]
    source: str
    status: str
    raw_sha256: str = ""
    source_upper_id: str = ""
    source_lower_id: str = ""
    source_upper_energy_cm1: float | None = None
    source_lower_energy_cm1: float | None = None
    mapping_confidence: str = ""


@dataclass(frozen=True, slots=True)
class SelectedCollision:
    collider_id: int
    upper: int
    lower: int
    upper_key: LevelKey
    lower_key: LevelKey
    rates: tuple[float, ...]
    source: str
    status: str
    multiplier: float = 1.0
    transform: str = "loglog_interpolation_no_extrapolation"
    raw_sha256: str = ""
    source_upper_id: str = ""
    source_lower_id: str = ""
    source_upper_key: str = ""
    source_lower_key: str = ""
    source_upper_energy_cm1: float | None = None
    source_lower_energy_cm1: float | None = None
    mapping_confidence: str = ""
    source_temperatures: tuple[float, ...] = ()
    source_rates: tuple[float, ...] = ()


def _split_floats(value: str) -> tuple[float, ...]:
    return tuple(float(part) for part in value.replace(",", ";").split(";") if part.strip())


def _optional_float(value: str | None) -> float | None:
    return float(value) if value not in (None, "") else None


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_master(root: str | Path) -> tuple[list[MasterLevel], list[MasterRadiative], list[MasterCollision]]:
    root = Path(root)
    directory = root / "data" / "normalized" / "master"
    paths = {name: directory / f"{name}.csv" for name in ("levels", "radiative", "collisions")}
    missing = [str(path.relative_to(root)) for path in paths.values() if not path.exists()]
    if missing:
        raise MissingSourcesError(missing)
    levels = [
        MasterLevel(
            key=LevelKey.decode(row["level_key"]), energy_cm1=float(row["energy_cm1"]),
            weight=float(row["weight"]), source=row["source"], source_index=row.get("source_index", ""),
            mapping_confidence=row.get("mapping_confidence", "unrecorded"),
            raw_sha256=row.get("raw_sha256", ""),
        )
        for row in _read_rows(paths["levels"])
    ]
    radiative = [
        MasterRadiative(
            species=row["species"].upper(), upper_key=LevelKey.decode(row["upper_key"]),
            lower_key=LevelKey.decode(row["lower_key"]), a_s=float(row["a_s"]),
            frequency_ghz=float(row["frequency_ghz"]), source=row["source"],
            status=row.get("status", "calculated"), frequency_source=row.get("frequency_source", row["source"]),
            raw_sha256=row.get("raw_sha256", ""),
            source_upper_id=row.get("source_upper_id", ""),
            source_lower_id=row.get("source_lower_id", ""),
            mapping_confidence=row.get("mapping_confidence", ""),
            smu2_debye2=float(row["smu2_debye2"]),
        )
        for row in _read_rows(paths["radiative"])
    ]
    collisions = [
        MasterCollision(
            species=row["species"].upper(), collider_id=int(row["collider_id"]),
            upper_key=LevelKey.decode(row["upper_key"]), lower_key=LevelKey.decode(row["lower_key"]),
            temperatures=_split_floats(row["temperatures_k"]), rates=_split_floats(row["rates_cm3_s"]),
            source=row["source"], status=row.get("status", "calculated"),
            raw_sha256=row.get("raw_sha256", ""),
            source_upper_id=row.get("source_upper_id", ""),
            source_lower_id=row.get("source_lower_id", ""),
            source_upper_energy_cm1=_optional_float(row.get("source_upper_energy_cm1")),
            source_lower_energy_cm1=_optional_float(row.get("source_lower_energy_cm1")),
            mapping_confidence=row.get("mapping_confidence", ""),
        )
        for row in _read_rows(paths["collisions"])
    ]
    return levels, radiative, collisions


def _balanced_levels(
    levels: list[MasterLevel], species: str, *, j_max: int, k_max: int, ground_max_cm1: float
) -> tuple[
    list[MasterLevel], set[tuple[int, int, str]], list[dict[str, object]]
]:
    by_vt: dict[int, dict[tuple[int, int, str], MasterLevel]] = {0: {}, 1: {}, 2: {}}
    for level in levels:
        if level.key.species == species and level.key.vt in by_vt:
            identity = level.key.rotational_identity
            if identity in by_vt[level.key.vt]:
                raise ValueError(f"duplicate canonical level {level.key.encode()}")
            by_vt[level.key.vt][identity] = level
    # The scientific template is selected from vt=0.  Intersecting all three
    # manifolds would silently hide a missing pump counterpart by shrinking the
    # ground-state model.
    template = {
        identity for identity in by_vt[0]
        if identity[0] <= j_max and abs(identity[1]) <= k_max
        and by_vt[0][identity].energy_cm1 <= ground_max_cm1
    }
    if not template:
        raise MissingSourcesError([f"balanced {species} vt=0,1,2 level template"])
    for vt in (1, 2):
        missing = template - set(by_vt[vt])
        if missing:
            examples = ", ".join(map(str, sorted(missing)[:5]))
            raise MissingSourcesError(
                [f"{len(missing)} {species} vt={vt} counterparts to selected vt=0 labels ({examples})"]
            )
    selected = [by_vt[vt][identity] for vt in (0, 1, 2) for identity in template]
    for level in selected:
        expected_weight = 2 * level.key.j + 1
        if not math.isclose(level.weight, expected_weight, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(
                f"{level.key.encode()} weight {level.weight} violates the declared 2J+1 convention; "
                "normalize any common nuclear-spin factor before building"
            )
    selected.sort(
        key=lambda row: (
            row.energy_cm1,
            int(row.source_index) if row.source_index.isdigit() else math.inf,
            row.key.encode(),
        )
    )
    # The final Xu Hamiltonian has one exact high-K A+/A- degeneracy in the
    # configured basis (vt=1,J=9,K=9).  Keep the exact value in the normalized
    # master, but make the second serialized level one IEEE-754 step higher so
    # LAMDA/RADEX readers that require strictly ordered indices can represent
    # both parity states.  This is a numerical encoding tie-break (~1e-13
    # cm-1), not a claimed spectroscopic splitting.
    regularized: list[MasterLevel] = []
    adjustments: list[dict[str, object]] = []
    for row in selected:
        if regularized and row.energy_cm1 <= regularized[-1].energy_cm1:
            if row.energy_cm1 != regularized[-1].energy_cm1:
                raise ValueError(
                    f"{species} canonical levels cannot be ordered without "
                    "crossing a distinct energy"
                )
            original = row.energy_cm1
            adjusted = math.nextafter(regularized[-1].energy_cm1, math.inf)
            row = replace(row, energy_cm1=adjusted)
            adjustments.append(
                {
                    "level_key": row.key.encode(),
                    "source_index": row.source_index,
                    "original_energy_cm1": original,
                    "serialized_energy_cm1": adjusted,
                    "offset_cm1": adjusted - original,
                    "method": "IEEE-754 nextafter toward +infinity",
                    "physical_claim": "none; exact Hamiltonian degeneracy retained in master",
                }
            )
        regularized.append(row)
    if any(
        b.energy_cm1 <= a.energy_cm1 for a, b in zip(regularized, regularized[1:])
    ):
        raise ValueError(f"{species} canonical energies are not strictly increasing")
    return regularized, template, adjustments


def _source_priority(source: str) -> int:
    return {
        MEKHTIEV_SOURCE: 0,
        XU_CONSTANT_DIPOLE_SOURCE: 0,
        DAGDIGIAN_SOURCE: 0,
        RF_H2_SOURCE: 1,
        RF_HE_SOURCE: 1,
        RF_TORSION_HE_SOURCE: 1,
        "cdms-032504-v3": 2,
        "jpl-032003-v3": 3,
    }.get(source, 10)


def _level_completeness(root: Path, species: str, rows: list[MasterLevel]) -> bool:
    """Require one complete, raw-bound canonical level dataset."""
    path = root / "data" / "normalized" / "master" / "completeness.json"
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        record = payload.get("levels", {}).get(species, {})
        dataset_id = record.get("dataset_id")
        if dataset_id not in {CDMS_SOURCE, XU_FROZEN_LEVEL_SOURCE}:
            return False
        source_rows = [
            row for row in rows if row.key.species == species and row.source == dataset_id
        ]
        species_rows = [row for row in rows if row.key.species == species]
        verify_dataset_record(
            root, record, dataset_id=dataset_id,
            normalized_transition_count=len(source_rows),
        )
        bound_hashes = {
            artifact["sha256"] for artifact in record.get("raw_artifacts", [])
            if isinstance(artifact, dict) and "sha256" in artifact
        }
        coverage = record.get("coverage", {})
        if dataset_id == CDMS_SOURCE:
            scope_complete = (
                int(coverage.get("j_max", -1)) >= 40
                and int(coverage.get("vt_max", -1)) >= 2
            )
        else:
            scope_complete = (
                int(coverage.get("j_max", -1)) >= 22
                and int(coverage.get("abs_k_max", -1)) >= 14
                and int(coverage.get("vt_max", -1)) >= 2
            )
        return (
            bool(source_rows)
            and len(source_rows) == len(species_rows)
            and scope_complete
            and not any(
            row.raw_sha256 not in bound_hashes
            or not row.source_index
            or row.mapping_confidence not in {"exact", "exact_quantum_identity"}
            for row in source_rows
            )
        )
    except (OSError, ValueError, TypeError, RawBindingError):
        return False


def _select_radiative(
    rows: list[MasterRadiative], species: str, key_to_index: dict[LevelKey, int], levels: list[MasterLevel],
    frequency_overrides: dict[frozenset[LevelKey], tuple[float, str]],
    *, minimum_smu2_debye2: float,
) -> tuple[list[RadiativeTransition], list[dict[str, object]], dict[str, object]]:
    if not math.isfinite(minimum_smu2_debye2) or minimum_smu2_debye2 < 0:
        raise ValueError("minimum radiative S mu^2 must be finite and nonnegative")
    energy = {row.key: row.energy_cm1 for row in levels}
    grouped: dict[frozenset[LevelKey], list[MasterRadiative]] = {}
    for row in rows:
        if row.species != species:
            continue
        if row.source not in {
            MEKHTIEV_SOURCE, XU_CONSTANT_DIPOLE_SOURCE, CDMS_SOURCE,
            "jpl-032003-v3",
        }:
            raise ValueError(f"unsupported radiative dataset id {row.source}")
        if row.upper_key not in key_to_index or row.lower_key not in key_to_index:
            continue
        if row.a_s <= 0 or row.frequency_ghz <= 0:
            raise ValueError("non-positive master radiative value")
        if row.smu2_debye2 <= 0 or not math.isfinite(row.smu2_debye2):
            raise ValueError("non-positive or non-finite master radiative S mu^2")
        if energy[row.upper_key] <= energy[row.lower_key]:
            raise MappingError(
                "radiative mapping error: source upper endpoint "
                f"{row.upper_key.encode()} is not above source lower endpoint "
                f"{row.lower_key.encode()}; Einstein A_ul cannot be reused after endpoint reversal"
            )
        grouped.setdefault(frozenset((row.upper_key, row.lower_key)), []).append(row)
    transitions: list[RadiativeTransition] = []
    provenance: list[dict[str, object]] = []
    rejected_smu2: list[float] = []
    retained_smu2: list[float] = []
    for candidates in grouped.values():
        chosen = min(candidates, key=lambda row: (_source_priority(row.source), row.source))
        if chosen.smu2_debye2 < minimum_smu2_debye2:
            rejected_smu2.append(chosen.smu2_debye2)
            continue
        retained_smu2.append(chosen.smu2_debye2)
        upper, lower = chosen.upper_key, chosen.lower_key
        override = frequency_overrides.get(frozenset((upper, lower)))
        frequency_ghz = override[0] if override else chosen.frequency_ghz
        frequency_source = override[1] if override else chosen.frequency_source
        a_s = chosen.a_s
        if override:
            # The strength and upper-state degeneracy are unchanged; A_ul is
            # therefore rescaled exactly for the adopted observed frequency.
            a_s *= (frequency_ghz / chosen.frequency_ghz) ** 3
        transitions.append(
            RadiativeTransition(0, key_to_index[upper], key_to_index[lower], a_s, frequency_ghz)
        )
        provenance.append(
            {
                "upper_key": upper.encode(), "lower_key": lower.encode(), "a_source": chosen.source,
                "frequency_source": frequency_source, "status": chosen.status,
                "raw_sha256": chosen.raw_sha256,
                "source_upper_id": chosen.source_upper_id,
                "source_lower_id": chosen.source_lower_id,
                "mapping_confidence": chosen.mapping_confidence,
                "a_s": a_s, "frequency_ghz": frequency_ghz,
                "smu2_debye2": chosen.smu2_debye2,
                "effective_transition_dipole_debye": math.sqrt(chosen.smu2_debye2),
                "duplicate_candidates": len(candidates),
            }
        )
    order = sorted(range(len(transitions)), key=lambda idx: (transitions[idx].upper, transitions[idx].lower))
    transitions = [replace(transitions[idx], index=number) for number, idx in enumerate(order, start=1)]
    provenance = [dict(provenance[idx], transition_index=number) for number, idx in enumerate(order, start=1)]
    selection = {
        "candidate_pairs_in_selected_level_set": len(grouped),
        "retained_transitions": len(transitions),
        "rejected_below_cutoff": len(rejected_smu2),
        "minimum_smu2_debye2": minimum_smu2_debye2,
        "minimum_effective_transition_dipole_debye": math.sqrt(minimum_smu2_debye2),
        "minimum_retained_smu2_debye2": min(retained_smu2) if retained_smu2 else None,
        "maximum_rejected_smu2_debye2": max(rejected_smu2) if rejected_smu2 else None,
    }
    return transitions, provenance, selection


def _radiative_strength_cutoff(policy: dict[str, object]) -> tuple[float, float]:
    """Return the declared effective-dipole and equivalent S mu^2 cutoffs."""
    effective_dipole = float(
        policy.get(
            "minimum_effective_transition_dipole_debye",
            DEFAULT_MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE,
        )
    )
    if not math.isfinite(effective_dipole) or effective_dipole <= 0:
        raise ValueError(
            "minimum_effective_transition_dipole_debye must be finite and positive"
        )
    return effective_dipole, effective_dipole**2


def _frequency_overrides(root: Path, species: str) -> dict[frozenset[LevelKey], tuple[float, str]]:
    path = root / "config" / "observed-frequencies.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: dict[frozenset[LevelKey], tuple[float, str]] = {}
    for row in payload.get("overrides", []):
        if row["species"].upper() != species:
            continue
        keys = frozenset((LevelKey.decode(row["upper_key"]), LevelKey.decode(row["lower_key"])))
        result[keys] = (float(row["frequency_mhz"]) / 1000.0, row["source"])
    return result


def _radiative_completeness(
    root: Path, species: str, rows: list[MasterRadiative]
) -> str | None:
    """Return the complete raw-bound radiative dataset selected for a species."""
    path = root / "data" / "normalized" / "master" / "completeness.json"
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        record = payload.get("radiative", {}).get(species, {})
        dataset_id = record.get("dataset_id")
        if dataset_id not in {MEKHTIEV_SOURCE, XU_CONSTANT_DIPOLE_SOURCE}:
            return None
        count = sum(
            row.species == species and row.source == dataset_id for row in rows
        )
        verify_dataset_record(
            root, record, dataset_id=dataset_id,
            normalized_transition_count=count,
        )
        bound_hashes = {
            artifact["sha256"] for artifact in record.get("raw_artifacts", [])
            if isinstance(artifact, dict) and "sha256" in artifact
        }
        source_rows = [
            row for row in rows if row.species == species and row.source == dataset_id
        ]
        if any(
            row.raw_sha256 not in bound_hashes
            or not row.source_upper_id or not row.source_lower_id
            or row.mapping_confidence not in {"exact", "exact_quantum_identity"}
            for row in source_rows
        ):
            return None
        coverage = record.get("coverage", {})
        scope_complete = bool(
            int(coverage.get("j_max", -1)) >= 22
            and int(coverage.get("abs_k_max", -1)) >= 14
            and int(coverage.get("vt_max", -1)) >= 2
        )
        if not scope_complete:
            return None
        if dataset_id == MEKHTIEV_SOURCE:
            if not (
                float(coverage.get("smu2_min_debye2", math.inf)) <= 1.0e-6
                and float(coverage.get("frequency_min_khz", math.inf)) <= 1.0
            ):
                return None
        else:
            expected_operator = {
                "mu_a": 0.899,
                "mu_a_cos_3gamma": 0.0,
                "mu_b": -1.44,
                "mu_b_cos_3gamma": 0.0,
                "mu_c_cos_3gamma": 0.0,
            }
            operator = record.get("operator", {})
            if not isinstance(operator, dict) or any(
                not math.isclose(
                    float(operator.get(name, math.nan)), value,
                    rel_tol=0.0, abs_tol=1.0e-15,
                )
                for name, value in expected_operator.items()
            ):
                return None
            if float(coverage.get("smu2_min_debye2", math.inf)) > 1.0e-14:
                return None
        return dataset_id
    except (OSError, ValueError, TypeError, RawBindingError):
        return None


def _collision_completeness(
    root: Path, species: str, rows: list[MasterCollision]
) -> bool:
    """Bind every required calculated collision family to pinned raw files."""
    path = root / "data" / "normalized" / "master" / "completeness.json"
    if not path.exists():
        return False
    required = {
        DAGDIGIAN_SOURCE: {(2, (0, 0)), (3, (0, 0))},
        RF_H2_SOURCE: {(2, (0, 0)), (2, (1, 1)), (2, (2, 2))},
        RF_HE_SOURCE: {(6, (0, 0)), (6, (1, 1)), (6, (2, 2))},
        RF_TORSION_HE_SOURCE: {(6, (0, 1)), (6, (1, 2))},
    }
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = payload.get("collisions", {}).get(species, {})
        for source, conditions in required.items():
            source_rows = [
                row for row in rows
                if row.species == species and row.source == source
                and row.status.startswith("calculated")
            ]
            verify_dataset_record(
                root, records.get(source, {}), dataset_id=source,
                normalized_transition_count=len(source_rows),
            )
            bound_hashes = {
                artifact["sha256"] for artifact in records[source].get("raw_artifacts", [])
                if isinstance(artifact, dict) and "sha256" in artifact
            }
            if any(
                row.raw_sha256 not in bound_hashes
                or not row.source_upper_id or not row.source_lower_id
                or row.mapping_confidence not in {"exact", "exact_quantum_identity"}
                or row.source_upper_energy_cm1 is None
                or row.source_lower_energy_cm1 is None
                or row.source_upper_energy_cm1 < row.source_lower_energy_cm1
                for row in source_rows
            ):
                return False
            actual = {
                (row.collider_id, tuple(sorted((row.upper_key.vt, row.lower_key.vt))))
                for row in source_rows
            }
            if not conditions <= actual:
                return False
        return True
    except (OSError, ValueError, TypeError, RawBindingError):
        return False


def _collision_choice_rank(row: MasterCollision) -> tuple[int, int, str]:
    source = row.source
    vt0 = row.upper_key.vt == row.lower_key.vt == 0
    if row.collider_id in {2, 3} and vt0 and source == DAGDIGIAN_SOURCE:
        return 0, 0, source
    if source in {RF_H2_SOURCE, RF_HE_SOURCE, RF_TORSION_HE_SOURCE}:
        return 1, 0, source
    return 5, 0, source


def _selected_collision_coverage(rows: list[SelectedCollision]) -> bool:
    """Ensure truncation/mapping did not erase a required partner or family."""
    present = {
        (row.source, row.collider_id, tuple(sorted((row.upper_key.vt, row.lower_key.vt))))
        for row in rows
    }
    required = {
        (DAGDIGIAN_SOURCE, 2, (0, 0)),
        (DAGDIGIAN_SOURCE, 3, (0, 0)),
        (RF_H2_SOURCE, 2, (1, 1)),
        (RF_H2_SOURCE, 2, (2, 2)),
        (RF_HE_SOURCE, 6, (0, 0)),
        (RF_HE_SOURCE, 6, (1, 1)),
        (RF_HE_SOURCE, 6, (2, 2)),
        (RF_TORSION_HE_SOURCE, 6, (0, 1)),
        (RF_TORSION_HE_SOURCE, 6, (1, 2)),
    }
    return required <= present


def _calculated_collisions(
    rows: list[MasterCollision], species: str, key_to_index: dict[LevelKey, int], levels: list[MasterLevel],
    target_temperatures: tuple[float, ...],
    *,
    allow_sqrt_endpoint_extrapolation: bool = False,
) -> list[SelectedCollision]:
    by_key = {level.key: level for level in levels}
    groups: dict[tuple[int, frozenset[LevelKey]], list[MasterCollision]] = {}
    for row in rows:
        if row.species != species or not row.status.startswith("calculated"):
            continue
        if row.source not in {
            DAGDIGIAN_SOURCE, RF_H2_SOURCE, RF_HE_SOURCE, RF_TORSION_HE_SOURCE
        }:
            raise ValueError(f"unsupported calculated collision dataset id {row.source}")
        if row.collider_id not in {2, 3, 6}:
            raise ValueError(f"unsupported collider id {row.collider_id}; expected 2, 3, or 6")
        if row.upper_key not in key_to_index or row.lower_key not in key_to_index:
            continue
        if len(row.temperatures) != len(row.rates) or any(rate < 0 or not math.isfinite(rate) for rate in row.rates):
            raise ValueError("malformed normalized collision grid")
        groups.setdefault((row.collider_id, frozenset((row.upper_key, row.lower_key))), []).append(row)
    selected: list[SelectedCollision] = []
    for (collider_id, _), candidates in groups.items():
        chosen = min(candidates, key=_collision_choice_rank)
        first, second = chosen.upper_key, chosen.lower_key
        e_first, e_second = by_key[first].energy_cm1, by_key[second].energy_cm1
        reoriented = e_first < e_second
        upper, lower = (second, first) if reoriented else (first, second)
        source_rates = chosen.rates
        transform = (
            "loglog_interpolation_with_nearest_endpoint_sqrt_temperature_extrapolation"
            if allow_sqrt_endpoint_extrapolation
            else "loglog_interpolation_no_extrapolation"
        )
        if reoriented:
            if (
                chosen.source_upper_energy_cm1 is None
                or chosen.source_lower_energy_cm1 is None
            ):
                raise MappingError(
                    "collision endpoint reversal requires both source term values"
                )
            source_gap = (
                chosen.source_upper_energy_cm1 - chosen.source_lower_energy_cm1
            )
            if source_gap < 0:
                raise MappingError(
                    "normalized source collision is not downward in its source energy basis"
                )
            # The downloaded coefficient is source-downward (first -> second),
            # but after adopting the Xu energy order the desired downward row
            # is the reverse process (second -> first).  Recover that source-
            # upward coefficient by detailed balance using the source gap.
            if source_gap == 0:
                source_rates = tuple(
                    rate * by_key[first].weight / by_key[second].weight
                    for rate in chosen.rates
                )
                transform = (
                    "source_degenerate_detailed_balance_reorientation_rate_ratio_"
                    "then_loglog_interpolation"
                )
            else:
                source_rates = tuple(
                    upward_from_downward(
                        rate,
                        by_key[first].weight,
                        by_key[second].weight,
                        source_gap,
                        temperature,
                    )
                    for rate, temperature in zip(chosen.rates, chosen.temperatures)
                )
                transform = (
                    "source_downward_to_upward_detailed_balance_reorientation_"
                    "then_loglog_interpolation"
                )
        rate_transform = (
            interpolate_or_sqrt_extrapolate_rate
            if allow_sqrt_endpoint_extrapolation else interpolate_rate
        )
        try:
            rates = tuple(
                rate_transform(chosen.temperatures, source_rates, t)
                for t in target_temperatures
            )
        except ValueError as exc:
            if "extrapolation" in str(exc):
                continue
            raise
        if allow_sqrt_endpoint_extrapolation and "sqrt_temperature" not in transform:
            transform += "; nearest_endpoint_sqrt_temperature_extrapolation"
        selected.append(
            SelectedCollision(
                collider_id, key_to_index[upper], key_to_index[lower], upper, lower, rates,
                chosen.source, chosen.status, transform=transform,
                raw_sha256=chosen.raw_sha256,
                source_upper_id=chosen.source_upper_id,
                source_lower_id=chosen.source_lower_id,
                source_upper_key=chosen.upper_key.encode(),
                source_lower_key=chosen.lower_key.encode(),
                source_upper_energy_cm1=chosen.source_upper_energy_cm1,
                source_lower_energy_cm1=chosen.source_lower_energy_cm1,
                mapping_confidence=chosen.mapping_confidence,
                source_temperatures=chosen.temperatures,
                source_rates=chosen.rates,
            )
        )
    return selected


def _rf2011_proxy(calculated: list[SelectedCollision]) -> list[SelectedCollision]:
    result = list(calculated)
    existing = {(row.collider_id, row.upper_key, row.lower_key) for row in result}
    for he in calculated:
        if he.collider_id != 6 or he.upper_key.vt == he.lower_key.vt:
            continue
        for collider_id, multiplier, label in ((2, 1.0, "p-H2=1xHe"), (3, 3.0, "o-H2=3xHe")):
            identity = (collider_id, he.upper_key, he.lower_key)
            if identity in existing:
                continue
            result.append(
                replace(
                    he, collider_id=collider_id, rates=tuple(multiplier * value for value in he.rates),
                    source="Rabli-Flower-2011-He-prescription", status="proxy",
                    multiplier=multiplier, transform=f"{he.transform}; {label}",
                )
            )
            existing.add(identity)
    return result


def _cragg_rate(upper: MasterLevel, lower: MasterLevel, temperature: float) -> tuple[float, str] | None:
    if upper.key.vt != lower.key.vt:
        return 5.0e-14, "constant_torsion-changing"
    if upper.key.vt == 0:
        # Appendix B treats missing high-J ground-state rows by inheritance
        # from a corresponding J=9 row; do not substitute the excited-state
        # closed-form propensities when no such analogue exists.
        return None
    delta_j = abs(upper.key.j - lower.key.j)
    delta_k = abs(upper.key.k - lower.key.k)
    is_a_doublet = (
        upper.key.species == "A" and delta_j == 0 and delta_k == 0
        and upper.key.parity != lower.key.parity
    )
    if is_a_doublet:
        base = 3.0e-11
        method = "propensity_A_asymmetry_doublet"
    elif delta_k == 0 and delta_j >= 1:
        base = 3.0e-11 / delta_j
        method = "propensity_deltaK0_inverse_deltaJ"
    elif delta_k > 0:
        base = 1.0e-11 / max(1, delta_j)
        method = "propensity_deltaK_nonzero_inverse_max1_deltaJ"
    else:
        # Cragg Appendix B does not give an unambiguous closed rule here.
        return None
    return base * math.sqrt(temperature / 20.0), method


def _cragg2005(
    calculated: list[SelectedCollision], levels: list[MasterLevel], key_to_index: dict[LevelKey, int],
    target_temperatures: tuple[float, ...],
) -> list[SelectedCollision]:
    """Add the Appendix-B mixture coefficients to p-H2 and He blocks.

    Cragg used one 80% p-H2/20% He density mixture and assigned the same
    coefficient to both partners.  It published no separate o-H2 propensity
    law, so collider 3 is never filled here.
    """
    result = list(calculated)
    existing = {(row.collider_id, row.upper_key, row.lower_key) for row in result}
    calculated_ph2 = {
        frozenset((row.upper_key, row.lower_key)): row
        for row in calculated if row.collider_id == 2
    }
    level_by_key = {level.key: level for level in levels}
    source = "Cragg-Sobolev-Godfrey-2005-Appendix-B"
    mixture_note = "assumed 80% p-H2/20% He density mixture; identical coefficient in collider 2 and 6"
    for first, second in itertools.combinations(levels, 2):
        upper, lower = (first, second) if first.energy_cm1 > second.energy_cm1 else (second, first)
        if upper.key.vt == lower.key.vt == 0:
            identity = (2, upper.key, lower.key)
            if identity in existing:
                continue
            max_j = max(upper.key.j, lower.key.j)
            analogue: SelectedCollision | None = None
            if max_j > 9:
                shift = max_j - 9
                try:
                    first_key = LevelKey(
                        upper.key.species, 0, upper.key.j - shift, upper.key.k, upper.key.parity
                    )
                    second_key = LevelKey(
                        lower.key.species, 0, lower.key.j - shift, lower.key.k, lower.key.parity
                    )
                except ValueError:
                    pass
                else:
                    if first_key in level_by_key and second_key in level_by_key:
                        analogue = calculated_ph2.get(frozenset((first_key, second_key)))
            is_a_doublet = (
                upper.key.species == "A" and upper.key.j == lower.key.j
                and upper.key.k == lower.key.k and upper.key.parity != lower.key.parity
            )
            if analogue is None and is_a_doublet:
                try:
                    doublet = frozenset((LevelKey("A", 0, 9, 4, "+"), LevelKey("A", 0, 9, 4, "-")))
                except ValueError:  # pragma: no cover - fixed valid quantum numbers
                    doublet = frozenset()
                analogue = calculated_ph2.get(doublet)
            if analogue is not None:
                result.append(
                    SelectedCollision(
                        2, key_to_index[upper.key], key_to_index[lower.key], upper.key, lower.key,
                        analogue.rates, source, "propensity",
                        transform="ground-vt J=9 analogue from calculated p-H2 row",
                    )
                )
                existing.add(identity)
            continue
        values = [_cragg_rate(upper, lower, temperature) for temperature in target_temperatures]
        if any(value is None for value in values):
            continue
        rates = tuple(value[0] for value in values if value is not None)
        method = values[0][1]  # type: ignore[index]
        status = "constant" if upper.key.vt != lower.key.vt else "propensity"
        temperature_transform = (
            "temperature-independent" if status == "constant" else "sqrt(T/20K)"
        )
        for collider_id in (2, 6):
            identity = (collider_id, upper.key, lower.key)
            if identity in existing:
                continue
            result.append(
                SelectedCollision(
                    collider_id, key_to_index[upper.key], key_to_index[lower.key], upper.key, lower.key,
                    rates, source, status,
                    transform=f"{method}; {temperature_transform}; {mixture_note}",
                )
            )
            existing.add(identity)
    return result


def _partners(rows: list[SelectedCollision], target_temperatures: tuple[float, ...]) -> list[CollisionPartner]:
    descriptions = {2: "CH3OH-pH2", 3: "CH3OH-oH2", 6: "CH3OH-He"}
    partners: list[CollisionPartner] = []
    for collider_id in (2, 3, 6):
        chosen = sorted((row for row in rows if row.collider_id == collider_id), key=lambda row: (row.upper, row.lower))
        if not chosen:
            continue
        transitions = [
            CollisionTransition(index, row.upper, row.lower, row.rates)
            for index, row in enumerate(chosen, start=1)
        ]
        partners.append(CollisionPartner(collider_id, descriptions[collider_id], target_temperatures, transitions))
    return partners


def _collision_coverage(
    filename: str, species: str, levels: list[MasterLevel], rows: list[SelectedCollision]
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    possible: dict[tuple[int, int], int] = {}
    for lower_vt in (0, 1, 2):
        for upper_vt in range(lower_vt, 3):
            left = sum(level.key.vt == lower_vt for level in levels)
            right = sum(level.key.vt == upper_vt for level in levels)
            possible[(lower_vt, upper_vt)] = left * (left - 1) // 2 if lower_vt == upper_vt else left * right
    present: Counter[tuple[int, int, int]] = Counter()
    for row in rows:
        pair = tuple(sorted((row.upper_key.vt, row.lower_key.vt)))
        present[(row.collider_id, pair[0], pair[1])] += 1
    j_values = [level.key.j for level in levels]
    for collider_id in (2, 3, 6):
        for (lower_vt, upper_vt), expected in possible.items():
            count = present[(collider_id, lower_vt, upper_vt)]
            result.append(
                {
                    "file": filename, "species": species, "vt_pair": f"{lower_vt}-{upper_vt}",
                    "collider_id": collider_id, "expected_pair_rows": expected,
                    "present_rows": count, "missing_rows": expected - count,
                    "j_min": min(j_values), "j_max": max(j_values),
                }
            )
    return result


def _partition_capture(
    selected: list[MasterLevel], all_levels: list[MasterLevel], species: str
) -> dict[str, float]:
    available = [level for level in all_levels if level.key.species == species]
    zero = min(level.energy_cm1 for level in available)
    result: dict[str, float] = {}
    for temperature in (100.0, 200.0, 300.0):
        total = sum(
            level.weight * math.exp(-(level.energy_cm1 - zero) * CM1_TO_K / temperature)
            for level in available
        )
        kept = sum(
            level.weight * math.exp(-(level.energy_cm1 - zero) * CM1_TO_K / temperature)
            for level in selected
        )
        result[f"{temperature:g}K"] = kept / total if total else 0.0
    return result


def _write_sidecars(
    root: Path, filename: str, radiative: list[dict[str, object]], collisions: list[SelectedCollision],
    target_temperatures: tuple[float, ...],
) -> None:
    sidecar = root / "outputs" / "provenance"
    sidecar.mkdir(parents=True, exist_ok=True)
    rad_fields = [
        "file", "transition_index", "upper_key", "lower_key", "a_s", "frequency_ghz",
        "smu2_debye2", "effective_transition_dipole_debye",
        "a_source", "frequency_source", "status", "raw_sha256", "source_upper_id",
        "source_lower_id", "mapping_confidence", "duplicate_candidates",
    ]
    with (sidecar / f"{filename}.radiative.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rad_fields, lineterminator="\n")
        writer.writeheader()
        for row in radiative:
            writer.writerow({"file": filename, **row})
    collision_fields = [
        "file", "transition_index", "collider_id", "upper_index", "lower_index",
        "upper_key", "lower_key", "temperatures_k", "rates_cm3_s",
        "source", "status", "multiplier", "transform",
        "raw_sha256", "source_upper_id", "source_lower_id", "source_upper_energy_cm1",
        "source_lower_energy_cm1", "source_upper_key", "source_lower_key",
        "mapping_confidence", "source_temperatures_k",
        "source_rates_cm3_s",
    ]
    with (sidecar / f"{filename}.collisions.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=collision_fields, lineterminator="\n")
        writer.writeheader()
        for collider_id in (2, 3, 6):
            ordered = sorted(
                (row for row in collisions if row.collider_id == collider_id),
                key=lambda row: (row.upper, row.lower),
            )
            for transition_index, row in enumerate(ordered, start=1):
                writer.writerow(
                    {
                        "file": filename, "transition_index": transition_index,
                        "collider_id": row.collider_id, "upper_index": row.upper,
                        "lower_index": row.lower, "upper_key": row.upper_key.encode(),
                        "lower_key": row.lower_key.encode(),
                        "temperatures_k": ";".join(
                            f"{value:.9g}" for value in target_temperatures
                        ),
                        "rates_cm3_s": ";".join(
                            f"{value:.12g}" for value in row.rates
                        ),
                        "source": row.source, "status": row.status,
                        "multiplier": row.multiplier, "transform": row.transform,
                        "raw_sha256": row.raw_sha256,
                        "source_upper_id": row.source_upper_id,
                        "source_lower_id": row.source_lower_id,
                        "source_upper_key": row.source_upper_key,
                        "source_lower_key": row.source_lower_key,
                        "source_upper_energy_cm1": row.source_upper_energy_cm1,
                        "source_lower_energy_cm1": row.source_lower_energy_cm1,
                        "mapping_confidence": row.mapping_confidence,
                        "source_temperatures_k": ";".join(f"{value:.9g}" for value in row.source_temperatures),
                        "source_rates_cm3_s": ";".join(f"{value:.12g}" for value in row.source_rates),
                    }
                )


def build_products(
    root: str | Path, *, accept_incomplete_radiative: bool = False,
    j_max: int | None = None, k_max: int | None = None, ground_max_cm1: float | None = None,
    target_temperatures: tuple[float, ...] | None = None,
) -> dict[str, object]:
    root = Path(root)
    policy_path = root / "config" / "build.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8")) if policy_path.exists() else {}
    j_max = int(policy.get("j_max", 18)) if j_max is None else j_max
    k_max = int(policy.get("abs_k_max", 9)) if k_max is None else k_max
    ground_max_cm1 = (
        float(policy.get("ground_manifold_max_energy_cm1", 417.0))
        if ground_max_cm1 is None else ground_max_cm1
    )
    target_temperatures = (
        tuple(float(value) for value in policy.get("collision_temperature_grid_k", COMMON_TEMPERATURES))
        if target_temperatures is None else target_temperatures
    )
    if not target_temperatures or any(
        right <= left for left, right in zip(target_temperatures, target_temperatures[1:])
    ):
        raise ValueError("configured collision temperature grid must increase strictly")
    minimum_effective_dipole, minimum_smu2 = _radiative_strength_cutoff(policy)
    levels, radiative_rows, collision_rows = load_master(root)
    collision_source_audit = audit_basecol_collision_master(root)
    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "configuration": {
            "vt": [0, 1, 2], "j_max": j_max, "k_max": k_max,
            "ground_max_cm1": ground_max_cm1, "temperatures_k": target_temperatures,
            "nuclear_spin_factor": "common factor stripped; weights 2J+1",
        },
        "radiative_strength_cutoff": {
            "interpretation": "mu_eff = sqrt(S mu^2)",
            "retention_rule": "mu_eff >= minimum; equivalently S mu^2 >= minimum squared",
            "minimum_effective_transition_dipole_debye": minimum_effective_dipole,
            "minimum_smu2_debye2": minimum_smu2,
        },
        "files": [],
        "level_energy_serialization_tie_breaks": {},
        "collision_source_level_audit": collision_source_audit,
    }
    prepared: list[
        tuple[str, LamdaData, list[dict[str, object]], list[SelectedCollision]]
    ] = []
    for species in ("A", "E"):
        if not _level_completeness(root, species, levels):
            raise MissingSourcesError([f"raw-bound canonical {species} level table"])
        selected_levels, template, energy_tie_breaks = _balanced_levels(
            levels, species, j_max=j_max, k_max=k_max, ground_max_cm1=ground_max_cm1
        )
        manifest["level_energy_serialization_tie_breaks"][species] = energy_tie_breaks  # type: ignore[index]
        key_to_index = {level.key: index for index, level in enumerate(selected_levels, start=1)}
        lamda_levels = [
            Level(index, level.energy_cm1, level.weight, level.key.encode())
            for index, level in enumerate(selected_levels, start=1)
        ]
        complete_radiative_source = _radiative_completeness(
            root, species, radiative_rows
        )
        radiative_complete = complete_radiative_source is not None
        strength_rows = (
            [row for row in radiative_rows if row.source == complete_radiative_source]
            if radiative_complete else radiative_rows
        )
        radiative, rad_provenance, rad_selection = _select_radiative(
            strength_rows, species, key_to_index, selected_levels,
            _frequency_overrides(root, species),
            minimum_smu2_debye2=minimum_smu2,
        )
        try:
            for row in rad_provenance:
                verify_pinned_sha256(root, str(row["raw_sha256"]))
                if (
                    not row["source_upper_id"] or not row["source_lower_id"]
                    or row["mapping_confidence"] not in {"exact", "exact_quantum_identity"}
                ):
                    raise RawBindingError("radiative mapping provenance is incomplete")
        except RawBindingError as exc:
            raise MissingSourcesError([f"raw-bound mapped {species} radiative rows: {exc}"]) from exc
        if not radiative:
            raise MissingSourcesError([f"mapped {species} radiative transitions"])
        if not radiative_complete and not accept_incomplete_radiative:
            raise MissingSourcesError([f"selection-complete {species} radiative table"])
        if not _collision_completeness(root, species, collision_rows):
            raise MissingSourcesError([f"raw-bound calculated {species} collision source families"])
        calculated = _calculated_collisions(
            collision_rows, species, key_to_index, selected_levels, target_temperatures
        )
        if not _selected_collision_coverage(calculated):
            raise MissingSourcesError([f"mapped calculated {species} collision source coverage"])
        variants = {
            "calculated": calculated,
            "rf2011-proxy": _rf2011_proxy(calculated),
            "cragg2005": _cragg2005(
                calculated, selected_levels, key_to_index, target_temperatures
            ),
        }
        for product, selected_collisions in variants.items():
            suffix = "" if radiative_complete else "-radiative-incomplete"
            filename = f"{species.lower()}-ch3oh-vt012-{product}{suffix}.dat"
            data = LamdaData(
                molecule=f"{species}-CH3OH vt=0,1,2 [{product}{suffix}]", mass_amu=32.04186,
                levels=lamda_levels, radiative=radiative,
                collision_partners=_partners(selected_collisions, target_temperatures),
            )
            assert_collision_selection(
                data,
                selected_collisions,
                target_temperatures,
                context=filename,
            )
            radiative_audit = audit_radiative_integrity(
                data,
                context=filename,
                source_rows=rad_provenance,
                minimum_smu2_debye2=minimum_smu2,
            )
            prepared.append((filename, data, rad_provenance, selected_collisions))
            status_counts: dict[str, int] = {}
            for row in selected_collisions:
                status_counts[row.status] = status_counts.get(row.status, 0) + 1
            manifest["files"].append(  # type: ignore[union-attr]
                {
                    "file": filename, "species": species, "product": product,
                    "radiative_status": (
                        f"complete_{complete_radiative_source}"
                        if radiative_complete else "incomplete_fallback"
                    ),
                    "levels": len(lamda_levels), "rotational_template_size": len(template),
                    "radiative_transitions": len(radiative),
                    "radiative_integrity": radiative_audit,
                    "radiative_strength_selection": rad_selection,
                    "collision_rows_by_status": status_counts,
                    "collision_coverage": _collision_coverage(
                        filename, species, selected_levels, selected_collisions
                    ),
                    "partition_capture_fraction_of_available_master": _partition_capture(
                        selected_levels, levels, species
                    ),
                    "partition_capture_warning": "A diagnostic against the available master only; not proof of convergence beyond 200 K gas or 300 K dust",
                }
            )
    # Nothing is installed until both A and E have passed all source, mapping,
    # balancing, and collision-transformation preflight checks.  The manifest
    # is promoted last and is the coherence marker for the staged six-file set.
    with tempfile.TemporaryDirectory(prefix=".methanol-lamda-stage-", dir=root) as directory:
        stage_root = Path(directory)
        stage_output = stage_root / "outputs"
        stage_output.mkdir(parents=True)
        for filename, data, rad_provenance, selected_collisions in prepared:
            write_lamda(stage_output / filename, data)
            _write_sidecars(
                stage_root, filename, rad_provenance, selected_collisions, target_temperatures
            )
        record_by_name = {row["file"]: row for row in manifest["files"]}  # type: ignore[index]
        for filename in record_by_name:
            data_path = stage_output / filename
            rad_path = stage_output / "provenance" / f"{filename}.radiative.csv"
            col_path = stage_output / "provenance" / f"{filename}.collisions.csv"
            record_by_name[filename]["sha256"] = sha256_path(data_path)
            record_by_name[filename]["radiative_sidecar_sha256"] = sha256_path(rad_path)
            record_by_name[filename]["collision_sidecar_sha256"] = sha256_path(col_path)
        build_inputs = {}
        for relative in (
            "config/build.json", "config/observed-frequencies.json",
            "data/raw/source-ledger.json",
            "data/normalized/master/levels.csv",
            "data/normalized/master/radiative.csv",
            "data/normalized/master/collisions.csv",
            "data/normalized/master/completeness.json",
        ):
            path = root / relative
            if path.exists():
                build_inputs[relative] = sha256_path(path)
        manifest["build_inputs_sha256"] = build_inputs
        manifest["validation_status"] = "candidate_pending_validation"
        staged_manifest = stage_output / "manifest.json"
        staged_manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (output_dir / "provenance").mkdir(parents=True, exist_ok=True)
        for path in sorted((stage_output / "provenance").iterdir()):
            os.replace(path, output_dir / "provenance" / path.name)
        for path in sorted(stage_output.glob("*.dat")):
            os.replace(path, output_dir / path.name)
        os.replace(staged_manifest, output_dir / "manifest.json")
    (root / "reports" / "build-blocker.json").unlink(missing_ok=True)
    return manifest
