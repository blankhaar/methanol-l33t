"""Enumerate Xu-style methanol strengths from the frozen BELGI eigenvectors.

The Xu et al. (2008) line list used a constant body-fixed dipole model with
``mu_a = +0.899 D`` and ``mu_b = -1.44 D``.  The three cos(3 gamma)
coefficients are deliberately zero here.  This module keeps that radiative
operator separate from the fitted Hamiltonian constants, drives the local
BELGI enumerator, validates its output, and writes a provenance-rich line
table with Einstein A coefficients.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
import csv
import io
import json
import math
import platform
from pathlib import Path
import subprocess
import sys
import tempfile

from .hamiltonian import (
    CST_PARAMETER_ORDER,
    FrozenLevel,
    HamiltonianError,
    _belgi_sources,
    _compile_driver,
    _expected_level_count,
    _find_compiler,
    _parse_level_stdout,
    _portable_compile_command,
    belgi_source_inventory,
    parse_final_constants,
    render_complete_cst_vector,
    sha256_file,
)
from .paths import portable_path
from .physics import einstein_a_from_smu2
from .quantum import LevelKey
from .xu import CM1_TO_MHZ, read_linecal_rows


DATASET_ID = "local-belgi-xu2008-constant-dipole"
DEFAULT_SMU2_MIN_DEBYE2 = 1.0e-14
DEFAULT_DIPOLE_MODEL_NAME = "Xu2008 constant-dipole line-strength model"
DEFAULT_LINE_DIPOLE_MODEL_ID = "Xu2008_constant_mu_a_mu_b_cos3gamma_zero"
DRIVER_FILENAME = "belgi_frozen_strengths.f"
DRIVER_EXECUTABLE = "belgi_frozen_strengths"

# Decimal spellings are passed directly to Fortran.  Slots 103, 105, and 106
# are the cos(3 gamma) terms and stay zero by construction.
XU2008_DIPOLE_LITERALS = {
    "AMUA0": "0.899D0",
    "AMUA3": "0.0D0",
    "AMUB0": "-1.44D0",
    "AMUB3": "0.0D0",
    "AMUC3": "0.0D0",
}
XU2008_DIPOLES_DEBYE = {
    "mu_a": 0.899,
    "mu_a_cos_3gamma": 0.0,
    "mu_b": -1.44,
    "mu_b_cos_3gamma": 0.0,
    "mu_c_cos_3gamma": 0.0,
}

_STRENGTH_BASE_FIELDS = (
    "species",
    "upper_vt", "upper_j", "upper_k", "upper_parity",
    "lower_vt", "lower_j", "lower_k", "lower_parity",
    "upper_energy_raw_cm1", "lower_energy_raw_cm1", "frequency_mhz",
    "smu2_debye2",
)
_AMPLITUDE_DIAGNOSTIC_FIELDS = frozenset(
    ("amplitude_a_debye", "amplitude_b_debye", "component_cancellation_factor")
)
_NATIVE_DIAGNOSTIC_FIELDS = frozenset(
    (
        "dipole_components", "contribution_a_percent", "contribution_b_percent",
        "cancellation_factor",
    )
)


class StrengthError(HamiltonianError):
    """Raised when constant-dipole enumeration or validation fails."""


@dataclass(frozen=True, slots=True)
class FrozenStrength:
    species: str
    upper_key: LevelKey
    lower_key: LevelKey
    upper_energy_raw_cm1: float
    lower_energy_raw_cm1: float
    frequency_mhz: float
    smu2_debye2: float
    amplitude_a_debye: float | None
    amplitude_b_debye: float | None
    dipole_components: str
    contribution_a_percent: int | None
    contribution_b_percent: int | None
    component_cancellation_factor: float

    @property
    def g_upper(self) -> int:
        return 2 * self.upper_key.j + 1

    @property
    def einstein_a_s(self) -> float:
        return einstein_a_from_smu2(
            self.frequency_mhz, self.smu2_debye2, self.g_upper
        )


def apply_xu2008_constant_dipoles(assignments: dict[str, str]) -> dict[str, str]:
    """Return a copy with BELGI dipole slots 102--106 fixed to Xu values."""
    unknown = set(assignments) - set(CST_PARAMETER_ORDER)
    if unknown:
        raise StrengthError(f"unknown BELGI constants: {sorted(unknown)}")
    result = dict(assignments)
    result.update(XU2008_DIPOLE_LITERALS)
    return result


def _parse_level_key(row: dict[str, str], prefix: str, species: str) -> LevelKey:
    return LevelKey(
        species,
        int(row[f"{prefix}_vt"]),
        int(row[f"{prefix}_j"]),
        int(row[f"{prefix}_k"]),
        row[f"{prefix}_parity"].strip(),
    )


def _parse_strength_table(stdout: str) -> list[FrozenStrength]:
    begin = "BEGIN_FROZEN_STRENGTHS_CSV\n"
    end = "END_FROZEN_STRENGTHS_CSV"
    if begin not in stdout or end not in stdout:
        raise StrengthError("BELGI output does not contain the frozen-strength CSV markers")
    table = stdout.split(begin, 1)[1].split(end, 1)[0]
    reader = csv.DictReader(io.StringIO(table))
    fieldnames = set(reader.fieldnames or ())
    missing = set(_STRENGTH_BASE_FIELDS) - fieldnames
    if missing:
        raise StrengthError(f"BELGI strength table is missing fields: {sorted(missing)}")
    has_amplitudes = _AMPLITUDE_DIAGNOSTIC_FIELDS <= fieldnames
    has_native_diagnostics = _NATIVE_DIAGNOSTIC_FIELDS <= fieldnames
    if not has_amplitudes and not has_native_diagnostics:
        raise StrengthError(
            "BELGI strength table needs either amplitude or native INTEN diagnostics"
        )

    strengths: list[FrozenStrength] = []
    seen: set[tuple[LevelKey, LevelKey]] = set()
    for source_row, row in enumerate(reader, 2):
        try:
            species = row["species"].strip().upper()
            upper = _parse_level_key(row, "upper", species)
            lower = _parse_level_key(row, "lower", species)
            base_values = tuple(
                float(row[name])
                for name in (
                    "upper_energy_raw_cm1", "lower_energy_raw_cm1", "frequency_mhz",
                    "smu2_debye2",
                )
            )
            if has_amplitudes:
                amp_a = float(row["amplitude_a_debye"])
                amp_b = float(row["amplitude_b_debye"])
                cancellation = float(row["component_cancellation_factor"])
            else:
                amp_a = None
                amp_b = None
                cancellation = float(row["cancellation_factor"])
            if has_native_diagnostics:
                dipole_components = row["dipole_components"].strip()
                contribution_a = int(row["contribution_a_percent"])
                contribution_b = int(row["contribution_b_percent"])
            else:
                dipole_components = ""
                contribution_a = None
                contribution_b = None
        except (KeyError, TypeError, ValueError) as exc:
            raise StrengthError(
                f"invalid BELGI strength row {source_row}: {exc}"
            ) from exc
        if upper.species != lower.species or upper.species != species:
            raise StrengthError(f"mixed-species strength row {source_row}")
        diagnostic_numbers = (cancellation,) if amp_a is None else (amp_a, amp_b, cancellation)
        if not all(math.isfinite(value) for value in (*base_values, *diagnostic_numbers)):
            raise StrengthError(f"non-finite value in BELGI strength row {source_row}")
        upper_energy, lower_energy, frequency, smu2 = base_values
        if upper_energy <= lower_energy or frequency <= 0 or smu2 <= 0:
            raise StrengthError(
                f"non-positive energy ordering, frequency, or strength at row {source_row}"
            )
        if abs(upper.j - lower.j) > 1 or (upper.j == 0 and lower.j == 0):
            raise StrengthError(f"row {source_row} violates rank-1 rotational selection")
        expected_frequency = (upper_energy - lower_energy) * CM1_TO_MHZ
        frequency_tolerance = max(1.0e-6, abs(expected_frequency) * 1.0e-10)
        if abs(frequency - expected_frequency) > frequency_tolerance:
            raise StrengthError(
                f"row {source_row} frequency is inconsistent with its raw energies"
            )
        # Native INTEN initializes CF=2 and leaves that sentinel untouched when
        # its internal unweighted S threshold is not crossed.  Preserve that
        # diagnostic for very weak rows rather than misrepresenting it as a
        # physical cancellation factor.
        cancellation_is_sentinel = abs(cancellation - 2.0) <= 1.0e-12
        if (
            not cancellation_is_sentinel
            and (cancellation < -1.0e-12 or cancellation > 1.0 + 1.0e-10)
        ):
            raise StrengthError(f"invalid cancellation factor at row {source_row}")
        if set(dipole_components) - {"A", "B"}:
            raise StrengthError(f"unexpected dipole component at row {source_row}")
        if contribution_a is not None and not (-100 <= contribution_a <= 100):
            raise StrengthError(f"invalid a-dipole contribution at row {source_row}")
        if contribution_b is not None and not (-100 <= contribution_b <= 100):
            raise StrengthError(f"invalid b-dipole contribution at row {source_row}")
        pair = (upper, lower)
        if pair in seen:
            raise StrengthError(
                f"duplicate BELGI strength identity {upper.encode()} -> {lower.encode()}"
            )
        seen.add(pair)
        strengths.append(
            FrozenStrength(
                species, upper, lower, upper_energy, lower_energy, frequency,
                smu2, amp_a, amp_b, dipole_components, contribution_a,
                contribution_b, cancellation,
            )
        )
    return strengths


def _summary_integer(value: str, name: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise StrengthError(f"invalid {name} in BELGI strength summary: {value!r}") from exc
    if parsed < 0:
        raise StrengthError(f"negative {name} in BELGI strength summary")
    return parsed


def _parse_strength_summary(stdout: str) -> dict[str, int | float]:
    lines = [
        line.strip() for line in stdout.splitlines()
        if line.strip().startswith("FROZEN_STRENGTH_SUMMARY,")
    ]
    if len(lines) != 1:
        raise StrengthError(
            f"expected one FROZEN_STRENGTH_SUMMARY line; found {len(lines)}"
        )
    tokens = next(csv.reader([lines[0]]))[1:]
    aliases = {
        "candidates": "candidate_count",
        "candidate_count": "candidate_count",
        "retained": "retained_count",
        "retained_count": "retained_count",
        "legacy_checks": "legacy_check_count",
        "legacy_check_count": "legacy_check_count",
        "maxdiff": "max_abs_legacy_difference_debye2",
        "max_abs_difference": "max_abs_legacy_difference_debye2",
        "max_abs_legacy_difference": "max_abs_legacy_difference_debye2",
        "max_abs_legacy_difference_debye2": "max_abs_legacy_difference_debye2",
    }
    raw: dict[str, str]
    if len(tokens) == 4 and not any("=" in token for token in tokens):
        raw = dict(
            zip(
                (
                    "candidate_count", "retained_count", "legacy_check_count",
                    "max_abs_legacy_difference_debye2",
                ),
                tokens,
                strict=True,
            )
        )
    elif tokens and all("=" in token for token in tokens):
        raw = {}
        for token in tokens:
            name, value = token.split("=", 1)
            normalized = aliases.get(name.strip().lower())
            if normalized is None or normalized in raw:
                raise StrengthError(f"invalid BELGI strength summary token {token!r}")
            raw[normalized] = value.strip()
    elif len(tokens) == 8:
        raw = {}
        for name, value in zip(tokens[::2], tokens[1::2], strict=True):
            normalized = aliases.get(name.strip().lower())
            if normalized is None or normalized in raw:
                raise StrengthError(f"invalid BELGI strength summary field {name!r}")
            raw[normalized] = value.strip()
    else:
        raise StrengthError("unrecognized BELGI strength summary format")

    required = {
        "candidate_count", "retained_count", "legacy_check_count",
        "max_abs_legacy_difference_debye2",
    }
    if set(raw) != required:
        raise StrengthError(
            f"BELGI strength summary fields differ from required fields: {sorted(raw)}"
        )
    try:
        maximum_difference = float(raw["max_abs_legacy_difference_debye2"])
    except ValueError as exc:
        raise StrengthError("invalid maximum legacy difference in strength summary") from exc
    if not math.isfinite(maximum_difference) or maximum_difference < 0:
        raise StrengthError("maximum legacy difference must be finite and non-negative")
    summary: dict[str, int | float] = {
        "candidate_count": _summary_integer(raw["candidate_count"], "candidate count"),
        "retained_count": _summary_integer(raw["retained_count"], "retained count"),
        "legacy_check_count": _summary_integer(
            raw["legacy_check_count"], "legacy check count"
        ),
        "max_abs_legacy_difference_debye2": maximum_difference,
    }
    floor_lines = [
        line.strip() for line in stdout.splitlines()
        if line.strip().startswith("FROZEN_STRENGTH_FLOOR_COUNTS,")
    ]
    if len(floor_lines) != 1:
        raise StrengthError(
            f"expected one FROZEN_STRENGTH_FLOOR_COUNTS line; found {len(floor_lines)}"
        )
    floor_tokens = next(csv.reader([floor_lines[0]]))[1:]
    if len(floor_tokens) != 6:
        raise StrengthError("invalid BELGI strength-floor summary")
    for name, value in zip(
        (
            "positive_strength_count", "strength_count_ge_1e-16",
            "strength_count_ge_1e-14", "strength_count_ge_1e-12",
            "strength_count_ge_1e-10", "strength_count_ge_1e-8",
        ),
        floor_tokens,
        strict=True,
    ):
        summary[name] = _summary_integer(value, name.replace("_", " "))
    floor_counts = [
        int(summary[name])
        for name in (
            "positive_strength_count", "strength_count_ge_1e-16",
            "strength_count_ge_1e-14", "strength_count_ge_1e-12",
            "strength_count_ge_1e-10", "strength_count_ge_1e-8",
        )
    ]
    if any(right > left for left, right in zip(floor_counts, floor_counts[1:])):
        raise StrengthError("strength-floor counts are not monotonic")
    if floor_counts[0] > int(summary["candidate_count"]):
        raise StrengthError("positive-strength count exceeds candidate count")
    return summary


def _parse_strength_stdout(
    stdout: str,
) -> tuple[list[FrozenLevel], list[FrozenStrength], dict[str, int | float]]:
    return (
        _parse_level_stdout(stdout),
        _parse_strength_table(stdout),
        _parse_strength_summary(stdout),
    )


def _expected_rank1_candidate_count(j_max: int, abs_k_max: int, vt_max: int) -> int:
    """Return unordered same-species pairs sent to the rank-1 evaluator.

    The exact rank-1 prefilter applies ``abs(delta J) <= 1`` and excludes
    J=0/J=0 pairs before invoking the matrix-element routine.
    """
    counts = [
        (vt_max + 1) * (1 + 2 * min(j, abs_k_max))
        for j in range(j_max + 1)
    ]
    per_species = sum(count * (count - 1) // 2 for count in counts)
    per_species += sum(left * right for left, right in zip(counts, counts[1:]))
    per_species -= counts[0] * (counts[0] - 1) // 2
    return 2 * per_species


def _validate_enumeration(
    levels: list[FrozenLevel],
    strengths: list[FrozenStrength],
    summary: dict[str, int | float],
    *,
    j_max: int,
    abs_k_max: int,
    vt_max: int,
    smu2_min_debye2: float,
) -> None:
    expected_levels = _expected_level_count(j_max, abs_k_max, vt_max)
    if len(levels) != expected_levels:
        raise StrengthError(f"expected {expected_levels} levels, received {len(levels)}")
    expected_candidates = _expected_rank1_candidate_count(j_max, abs_k_max, vt_max)
    if summary["candidate_count"] != expected_candidates:
        raise StrengthError(
            f"expected {expected_candidates} rank-1 candidates, received "
            f"{summary['candidate_count']}"
        )
    if summary["retained_count"] != len(strengths):
        raise StrengthError("BELGI retained-count summary disagrees with the strength table")

    energies = {level.key: level.energy_raw_cm1 for level in levels}
    for strength in strengths:
        if strength.upper_key not in energies or strength.lower_key not in energies:
            raise StrengthError("strength table references a level absent from the level table")
        for reported, expected, key in (
            (strength.upper_energy_raw_cm1, energies[strength.upper_key], strength.upper_key),
            (strength.lower_energy_raw_cm1, energies[strength.lower_key], strength.lower_key),
        ):
            if abs(reported - expected) > max(1.0e-11, abs(expected) * 1.0e-12):
                raise StrengthError(f"strength energy disagrees for {key.encode()}")
        if strength.smu2_debye2 < smu2_min_debye2 * (1.0 - 1.0e-12):
            raise StrengthError("strength table contains a row below its declared threshold")


def _compare_reference_levels(
    levels: list[FrozenLevel], reference_path: Path, *, project_root: Path
) -> dict[str, object]:
    """Compare strength-run energies with the separately frozen level artifact."""
    if not reference_path.is_file():
        return {
            "status": "reference_not_present",
            "path": portable_path(project_root, reference_path),
        }
    reference: dict[LevelKey, float] = {}
    with reference_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = LevelKey.decode(row["level_key"])
            reference[key] = float(row["energy_raw_cm1"])
    mismatches: list[dict[str, object]] = []
    for level in levels:
        expected = reference.get(level.key)
        if expected is None or abs(level.energy_raw_cm1 - expected) > 1.0e-11:
            mismatches.append(
                {
                    "level_key": level.key.encode(),
                    "strength_run_raw_cm1": level.energy_raw_cm1,
                    "reference_raw_cm1": expected,
                }
            )
    if mismatches:
        raise StrengthError(
            f"strength-run energies disagree with frozen levels: {mismatches[:5]}"
        )
    return {
        "status": "passed",
        "path": portable_path(project_root, reference_path),
        "sha256": sha256_file(reference_path),
        "matched_level_count": len(levels),
    }


def _radiative_graph_stats(
    levels: list[FrozenLevel], strengths: list[FrozenStrength]
) -> dict[str, dict[str, object]]:
    """Return connected-component and torsional-coupling counts by species."""
    result: dict[str, dict[str, object]] = {}
    for species in ("A", "E"):
        nodes = {level.key for level in levels if level.key.species == species}
        adjacency = {key: set() for key in nodes}
        vt_pairs: Counter[tuple[int, int]] = Counter()
        edge_count = 0
        for row in strengths:
            if row.species != species:
                continue
            adjacency[row.upper_key].add(row.lower_key)
            adjacency[row.lower_key].add(row.upper_key)
            vt_pairs[tuple(sorted((row.upper_key.vt, row.lower_key.vt)))] += 1
            edge_count += 1
        isolated = sum(not neighbors for neighbors in adjacency.values())
        unseen = set(nodes)
        component_sizes: list[int] = []
        while unseen:
            seed = unseen.pop()
            stack = [seed]
            size = 0
            while stack:
                key = stack.pop()
                size += 1
                new = adjacency[key] & unseen
                unseen.difference_update(new)
                stack.extend(new)
            component_sizes.append(size)
        component_sizes.sort(reverse=True)
        result[species] = {
            "node_count": len(nodes),
            "edge_count": edge_count,
            "isolated_node_count": isolated,
            "connected_component_count": len(component_sizes),
            "component_sizes": component_sizes,
            "vt_pair_edge_counts": {
                f"{left}-{right}": count
                for (left, right), count in sorted(vt_pairs.items())
            },
        }
    return result


def _write_strengths_csv(
    path: Path,
    strengths: list[FrozenStrength],
    *,
    source_id: str,
    dipole_model_id: str = DEFAULT_LINE_DIPOLE_MODEL_ID,
) -> None:
    fields = [
        "species", "transition_index", "upper_key", "lower_key",
        "upper_vt", "upper_j", "upper_k", "upper_parity",
        "lower_vt", "lower_j", "lower_k", "lower_parity",
        "upper_energy_raw_cm1", "lower_energy_raw_cm1", "frequency_mhz",
        "smu2_debye2", "a_s", "g_upper", "amplitude_a_debye",
        "amplitude_b_debye", "dipole_components", "contribution_a_percent",
        "contribution_b_percent", "component_cancellation_factor",
        "cancellation_factor_available", "source", "dipole_model",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for index, strength in enumerate(strengths, 1):
            writer.writerow(
                {
                    "species": strength.species,
                    "transition_index": index,
                    "upper_key": strength.upper_key.encode(),
                    "lower_key": strength.lower_key.encode(),
                    "upper_vt": strength.upper_key.vt,
                    "upper_j": strength.upper_key.j,
                    "upper_k": strength.upper_key.k,
                    "upper_parity": strength.upper_key.parity,
                    "lower_vt": strength.lower_key.vt,
                    "lower_j": strength.lower_key.j,
                    "lower_k": strength.lower_key.k,
                    "lower_parity": strength.lower_key.parity,
                    "upper_energy_raw_cm1": f"{strength.upper_energy_raw_cm1:.17g}",
                    "lower_energy_raw_cm1": f"{strength.lower_energy_raw_cm1:.17g}",
                    "frequency_mhz": f"{strength.frequency_mhz:.17g}",
                    "smu2_debye2": f"{strength.smu2_debye2:.17g}",
                    "a_s": f"{strength.einstein_a_s:.17g}",
                    "g_upper": strength.g_upper,
                    "amplitude_a_debye": (
                        "" if strength.amplitude_a_debye is None
                        else f"{strength.amplitude_a_debye:.17g}"
                    ),
                    "amplitude_b_debye": (
                        "" if strength.amplitude_b_debye is None
                        else f"{strength.amplitude_b_debye:.17g}"
                    ),
                    "dipole_components": strength.dipole_components,
                    "contribution_a_percent": (
                        "" if strength.contribution_a_percent is None
                        else strength.contribution_a_percent
                    ),
                    "contribution_b_percent": (
                        "" if strength.contribution_b_percent is None
                        else strength.contribution_b_percent
                    ),
                    "component_cancellation_factor": (
                        f"{strength.component_cancellation_factor:.17g}"
                    ),
                    "cancellation_factor_available": (
                        abs(strength.component_cancellation_factor - 2.0) > 1.0e-12
                    ),
                    "source": source_id,
                    "dipole_model": dipole_model_id,
                }
            )


def _write_xu_observation_crosswalk(
    path: Path,
    xu_output: Path,
    levels: list[FrozenLevel],
    strengths: list[FrozenStrength],
) -> dict[str, object]:
    """Retain every Xu measurement and map only exact endpoint identities.

    Multiple measurements of one transition remain separate rows.  The 142
    locally spectrum-equivalent A rows are deliberately left unresolved: a
    rowwise parity flip is not a globally valid level mapping.
    """
    energies = {level.key: level.energy_raw_cm1 for level in levels}
    line_indices = {
        (strength.upper_key, strength.lower_key): index
        for index, strength in enumerate(strengths, 1)
    }
    line_by_index = {index: strength for index, strength in enumerate(strengths, 1)}
    fields = [
        "source_line", "species", "xu_upper_key", "xu_lower_key",
        "observation", "unit", "observed_mhz", "xu_calculation",
        "xu_calculated_mhz", "obs_minus_calc", "fit_status", "weight_code",
        "comment", "mapping_status", "line_transition_index",
        "generated_upper_key", "generated_lower_key", "generated_frequency_mhz",
        "direct_calculation_error", "direct_identity_tolerance",
    ]
    status_counts: Counter[str] = Counter()
    fit_counts: Counter[str] = Counter()
    mapped_line_indices: set[int] = set()
    source_rows = read_linecal_rows(xu_output)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in source_rows:
            fit_counts[row.fit_status] += 1
            line_index: int | None = None
            generated: FrozenStrength | None = None
            direct_error: float | None = None
            tolerance: float | None = None
            if row.upper_key not in energies or row.lower_key not in energies:
                status = "outside_frozen_level_scope"
            else:
                difference_cm1 = energies[row.upper_key] - energies[row.lower_key]
                predicted = difference_cm1 * CM1_TO_MHZ if row.unit == "MHz" else difference_cm1
                direct_error = predicted - row.calculation
                tolerance = 0.0005001 if row.unit == "MHz" else 0.000005001
                if abs(direct_error) > tolerance:
                    status = "quantum_identity_unresolved"
                else:
                    line_index = line_indices.get((row.upper_key, row.lower_key))
                    if line_index is None:
                        line_index = line_indices.get((row.lower_key, row.upper_key))
                    if line_index is None:
                        status = "direct_identity_below_strength_floor_or_zero"
                    else:
                        status = "mapped_to_retained_line"
                        generated = line_by_index[line_index]
                        mapped_line_indices.add(line_index)
            status_counts[status] += 1
            writer.writerow(
                {
                    "source_line": row.source_line,
                    "species": row.species,
                    "xu_upper_key": row.upper_key.encode(),
                    "xu_lower_key": row.lower_key.encode(),
                    "observation": f"{row.observation:.17g}",
                    "unit": row.unit,
                    "observed_mhz": f"{row.observed_mhz:.17g}",
                    "xu_calculation": f"{row.calculation:.17g}",
                    "xu_calculated_mhz": f"{row.calculated_mhz:.17g}",
                    "obs_minus_calc": (
                        "" if row.obs_minus_calc is None else f"{row.obs_minus_calc:.17g}"
                    ),
                    "fit_status": row.fit_status,
                    "weight_code": "" if row.weight_code is None else row.weight_code,
                    "comment": row.comment,
                    "mapping_status": status,
                    "line_transition_index": "" if line_index is None else line_index,
                    "generated_upper_key": "" if generated is None else generated.upper_key.encode(),
                    "generated_lower_key": "" if generated is None else generated.lower_key.encode(),
                    "generated_frequency_mhz": (
                        "" if generated is None else f"{generated.frequency_mhz:.17g}"
                    ),
                    "direct_calculation_error": (
                        "" if direct_error is None else f"{direct_error:.17g}"
                    ),
                    "direct_identity_tolerance": (
                        "" if tolerance is None else f"{tolerance:.17g}"
                    ),
                }
            )
    mapped_rows = status_counts["mapped_to_retained_line"]
    return {
        "source_path": str(xu_output),
        "source_sha256": sha256_file(xu_output),
        "source_row_count": len(source_rows),
        "fit_status_counts": dict(sorted(fit_counts.items())),
        "mapping_status_counts": dict(sorted(status_counts.items())),
        "mapped_measurement_rows": mapped_rows,
        "unique_retained_lines_with_measurements": len(mapped_line_indices),
        "additional_measurements_of_already_mapped_lines": (
            mapped_rows - len(mapped_line_indices)
        ),
        "frequency_replacement_applied": False,
        "frequency_policy": (
            "all measurements are preserved; choosing among repeated observations "
            "and replacing predicted frequencies is a separate reviewed step"
        ),
    }


def build_frozen_strengths(
    project_root: str | Path,
    constants_path: str | Path,
    belgi_source_dir: str | Path,
    output_dir: str | Path,
    *,
    j_max: int = 22,
    abs_k_max: int = 14,
    vt_max: int = 2,
    smu2_min_debye2: float = DEFAULT_SMU2_MIN_DEBYE2,
    compiler: str | Path | None = None,
    require_xu_final: bool = True,
    xu_output: str | Path | None = None,
    dataset_id: str = DATASET_ID,
    reference_levels_path: str | Path | None = None,
    dipole_model_name: str = DEFAULT_DIPOLE_MODEL_NAME,
    dipole_model_id: str = DEFAULT_LINE_DIPOLE_MODEL_ID,
    dipole_provenance_path: str | Path | None = None,
    completeness_claim: str | None = None,
    parity_repair: bool = False,
    cached_dipoles: bool = False,
) -> dict[str, object]:
    """Compile, enumerate, validate, and write Xu constant-dipole strengths."""
    if parity_repair and dataset_id == DATASET_ID:
        dataset_id += "-parity-full9-candidate-v1"
    root = Path(project_root).resolve()
    constants = Path(constants_path).expanduser().resolve()
    source_dir = Path(belgi_source_dir).expanduser().resolve()
    target = Path(output_dir).expanduser().resolve()
    if not (0 <= vt_max <= 4 and 0 <= abs_k_max <= min(j_max, 14) and 1 <= j_max <= 22):
        raise StrengthError(
            "require 0 <= vt_max <= 4, 0 <= abs_k_max <= min(j_max, 14), "
            "and 1 <= j_max <= 22"
        )
    if not math.isfinite(smu2_min_debye2) or smu2_min_debye2 <= 0:
        raise StrengthError("S mu^2 threshold must be finite and positive")

    hamiltonian_assignments = parse_final_constants(
        constants, require_xu_final=require_xu_final
    )
    assignments = apply_xu2008_constant_dipoles(hamiltonian_assignments)
    cst_vector = render_complete_cst_vector(assignments)
    compiler_path = _find_compiler(compiler)
    implementations, includes = _belgi_sources(source_dir)

    with tempfile.TemporaryDirectory(prefix="belgi-frozen-strengths-") as temporary:
        build_dir = Path(temporary)
        executable, command, compiler_version, source_patches = _compile_driver(
            root,
            source_dir,
            build_dir,
            compiler_path,
            parity_repair=parity_repair,
            cached_dipoles=cached_dipoles,
            driver_filename=DRIVER_FILENAME,
            executable_name=DRIVER_EXECUTABLE,
        )
        compile_command_template = _portable_compile_command(
            command, root, source_dir, build_dir
        )
        completed = subprocess.run(
            [
                str(executable), str(j_max), str(abs_k_max), str(vt_max),
                f"{smu2_min_debye2:.17g}",
            ],
            input=cst_vector,
            text=True,
            capture_output=True,
            cwd=build_dir,
            check=False,
        )
        if completed.returncode != 0:
            raise StrengthError(
                f"BELGI frozen-strength run failed ({completed.returncode})\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        levels, strengths, summary = _parse_strength_stdout(completed.stdout)

    _validate_enumeration(
        levels,
        strengths,
        summary,
        j_max=j_max,
        abs_k_max=abs_k_max,
        vt_max=vt_max,
        smu2_min_debye2=smu2_min_debye2,
    )
    reference_path = (
        root / "data" / "normalized" / "belgi-xu2008-frozen" / "levels.csv"
        if reference_levels_path is None
        else Path(reference_levels_path).expanduser()
    )
    if not reference_path.is_absolute():
        reference_path = root / reference_path
    reference_path = reference_path.resolve()
    reference_level_validation = _compare_reference_levels(
        levels, reference_path, project_root=root
    )
    radiative_graph = _radiative_graph_stats(levels, strengths)

    target.mkdir(parents=True, exist_ok=True)
    lines_path = target / "lines.csv"
    _write_strengths_csv(
        lines_path,
        strengths,
        source_id=dataset_id,
        dipole_model_id=dipole_model_id,
    )
    observation_path: Path | None = None
    observation_validation: dict[str, object] | None = None
    if xu_output is not None:
        xu_path = Path(xu_output).expanduser().resolve()
        observation_path = target / "observations.csv"
        observation_validation = _write_xu_observation_crosswalk(
            observation_path, xu_path, levels, strengths
        )
        observation_validation["source_path"] = portable_path(root, xu_path)
    source_inventory = belgi_source_inventory(source_dir)
    driver_path = root / "fortran" / DRIVER_FILENAME
    constants_sha256 = sha256_file(constants)
    dipole_note_path = (
        root / "data" / "raw" / "user-supplied" / "xu2008-dipole-model.txt"
        if dipole_provenance_path is None
        else Path(dipole_provenance_path).expanduser()
    )
    if not dipole_note_path.is_absolute():
        dipole_note_path = root / dipole_note_path
    dipole_note_path = dipole_note_path.resolve()
    dipole_provenance = {
        "kind": "project-owner-supplied Xu et al. quotation and explicit zero-term instruction",
        "path": portable_path(root, dipole_note_path),
        "sha256": sha256_file(dipole_note_path) if dipole_note_path.is_file() else None,
    }
    manifest: dict[str, object] = {
        "schema_version": 2,
        "dataset_id": dataset_id,
        "method": (
            "complete same-species rank-1 enumeration from frozen BELGI eigenvectors; "
            "no observations read and no fit performed"
        ),
        "hamiltonian": {
            "constants_path": portable_path(root, constants),
            "constants_sha256": constants_sha256,
            "named_assignments": len(hamiltonian_assignments),
            "total_cst_slots": len(CST_PARAMETER_ORDER),
        },
        "dipole_model": {
            "name": dipole_model_name,
            "units": "Debye",
            "coefficients": XU2008_DIPOLES_DEBYE,
            "belgi_slots_102_through_106": XU2008_DIPOLE_LITERALS,
            "overrides_hamiltonian_file": True,
            "cos_3gamma_terms_forced_zero": True,
            "provenance": dipole_provenance,
        },
        "selection": {
            "same_species_only": True,
            "rank": 1,
            "delta_j_absolute_max": 1,
            "j0_to_j0_excluded": True,
            "nominal_delta_k_prefilter": False,
            "smu2_min_debye2": smu2_min_debye2,
            "selection_complete_above_threshold": True,
            "floor_sensitivity_counts": {
                "positive": summary["positive_strength_count"],
                "ge_1e-16_debye2": summary["strength_count_ge_1e-16"],
                "ge_1e-14_debye2": summary["strength_count_ge_1e-14"],
                "ge_1e-12_debye2": summary["strength_count_ge_1e-12"],
                "ge_1e-10_debye2": summary["strength_count_ge_1e-10"],
                "ge_1e-8_debye2": summary["strength_count_ge_1e-8"],
            },
        },
        "scope": {
            "j_max": j_max,
            "abs_k_max": abs_k_max,
            "vt_max": vt_max,
            "level_count": len(levels),
            "candidate_count": summary["candidate_count"],
            "retained_count": len(strengths),
            "a_retained_count": sum(row.species == "A" for row in strengths),
            "e_retained_count": sum(row.species == "E" for row in strengths),
        },
        "validation": {
            "calculation_backend": ("exact full-operator BLAS contraction" if cached_dipoles else "native BELGI INTEN"),
            "parity_repair_enabled": parity_repair,
            "cached_complete_operator": cached_dipoles,
            "native_inten_evaluation_count": summary["legacy_check_count"],
            "surrogate_contraction_used": False,
            "energies_reproduced_in_same_driver_run": True,
            "separate_frozen_level_artifact": reference_level_validation,
            "radiative_graph": radiative_graph,
            "einstein_a_formula": (
                "A_ul = 1.16395e-20 * frequency_MHz^3 * "
                "Smu2_Debye2 / (2*J_upper+1) s^-1"
            ),
        },
        "program": {
            "driver": portable_path(root, driver_path),
            "driver_sha256": sha256_file(driver_path),
            "belgi_source_dir": portable_path(root, source_dir),
            "belgi_source_inventory": source_inventory,
            "compiled_belgi_source_sha256": {
                path.name: sha256_file(path)
                for path in sorted((*implementations, *includes), key=lambda path: path.name)
            },
            "compiler": compiler_path.name,
            "compiler_version": compiler_version,
            "compile_command_template": compile_command_template,
            "runtime_environment": {
                "platform": f"{platform.system()} {platform.release()}",
                "machine": platform.machine(),
                "lapack_linkage": (
                    "Apple Accelerate" if sys.platform == "darwin" else "system LAPACK/BLAS"
                ),
            },
            "source_patches": source_patches,
        },
        "outputs": {
            "lines_csv": portable_path(root, lines_path),
            "lines_sha256": sha256_file(lines_path),
        },
    }
    if observation_validation is not None and observation_path is not None:
        manifest["observed_frequencies"] = observation_validation
        manifest["outputs"]["observations_csv"] = portable_path(root, observation_path)  # type: ignore[index]
        manifest["outputs"]["observations_sha256"] = sha256_file(observation_path)  # type: ignore[index]
    completeness_path = target / "completeness.json"
    per_species_candidates = _expected_rank1_candidate_count(
        j_max, abs_k_max, vt_max
    ) // 2
    completeness_limitations = [
        "all cos(3 gamma) dipole coefficients are zero by project instruction",
        "predicted Hamiltonian frequencies remain in lines.csv",
    ]
    if xu_output is not None:
        completeness_limitations.extend(
            [
                "observations.csv preserves repeated Xu measurements without choosing a replacement",
                "142 A-species Xu rows lack a globally consistent parity-identity mapping",
            ]
        )
    completeness = {
        "schema_version": 2,
        "dataset_id": dataset_id,
        "claim": (
            "selection-complete for the declared Xu2008 constant-dipole operator"
            if completeness_claim is None else completeness_claim
        ),
        "not_claimed": "complete torsion-dependent Mekhtiev dipole model",
        "scope": {"j_max": j_max, "abs_k_max": abs_k_max, "vt_max": vt_max},
        "operator": XU2008_DIPOLES_DEBYE,
        "smu2_min_debye2": smu2_min_debye2,
        "species": {
            "A": {
                "candidate_count": per_species_candidates,
                "retained_count": sum(row.species == "A" for row in strengths),
                "selection_complete_above_threshold": True,
                "radiative_graph": radiative_graph["A"],
            },
            "E": {
                "candidate_count": per_species_candidates,
                "retained_count": sum(row.species == "E" for row in strengths),
                "selection_complete_above_threshold": True,
                "radiative_graph": radiative_graph["E"],
            },
        },
        "source_bindings": {
            "hamiltonian_constants_sha256": constants_sha256,
            "dipole_instruction_sha256": dipole_provenance["sha256"],
            "belgi_driver_sha256": sha256_file(driver_path),
            "belgi_source_tree_sha256": source_inventory["tree_sha256"],
            "lines_sha256": sha256_file(lines_path),
        },
        "limitations": completeness_limitations,
    }
    completeness_path.write_text(
        json.dumps(completeness, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest["outputs"]["completeness_json"] = portable_path(root, completeness_path)  # type: ignore[index]
    manifest["outputs"]["completeness_sha256"] = sha256_file(completeness_path)  # type: ignore[index]
    manifest_path = target / "manifest.json"
    manifest["outputs"]["manifest_json"] = portable_path(root, manifest_path)  # type: ignore[index]
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest
