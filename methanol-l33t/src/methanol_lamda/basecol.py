"""Strict reader and quantum-number mapper for BASECOL 2023 text tables.

The Rabli & Flower CH3OH--He downloads comprise one ``.lev`` level table and
one ``.kij`` rate table per symmetry species.  Source energies are retained for
provenance and ordering checks only; callers map the exact quantum identities
onto the project's canonical spectroscopic energies.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable

from .quantum import LevelKey


class BasecolFormatError(ValueError):
    """Raised when a BASECOL file is malformed or internally inconsistent."""


@dataclass(frozen=True, slots=True)
class BasecolLevel:
    source_index: int
    source_energy_cm1: float
    key: LevelKey
    rovibrational_symmetry: str


@dataclass(frozen=True, slots=True)
class BasecolRate:
    initial_level: int
    final_level: int
    rates_cm3_s: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class BasecolRateTable:
    temperatures_k: tuple[float, ...]
    transitions: tuple[BasecolRate, ...]
    maximum_initial_level: int
    maximum_final_level: int


def _finite_float(value: str, *, context: str) -> float:
    try:
        result = float(value.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise BasecolFormatError(f"{context}: invalid floating-point value {value!r}") from exc
    if not math.isfinite(result):
        raise BasecolFormatError(f"{context}: non-finite value {value!r}")
    return result


def _next_nonempty(lines: list[str], start: int) -> tuple[int, str]:
    for index in range(start, len(lines)):
        if lines[index].strip():
            return index, lines[index].strip()
    raise BasecolFormatError("unexpected end of BASECOL header")


def _value_after(lines: list[str], label_index: int, *, context: str) -> int:
    _, value = _next_nonempty(lines, label_index + 1)
    try:
        return int(value)
    except ValueError as exc:
        raise BasecolFormatError(f"{context}: expected an integer after header label") from exc


def _canonical_key(
    species: str, *, vt: int, j: int, k: int, symmetry: str
) -> LevelKey:
    species = species.upper()
    if species == "A":
        if symmetry not in {"A+", "A-"}:
            raise BasecolFormatError(
                f"A-CH3OH level has unexpected rovibrational symmetry {symmetry!r}"
            )
        if k == 0 and symmetry != "A+":
            raise BasecolFormatError("A-CH3OH K=0 source level must have A+ symmetry")
        return LevelKey("A", vt, j, k, symmetry[-1])
    if species == "E":
        if symmetry not in {"E1", "E2"}:
            raise BasecolFormatError(
                f"E-CH3OH level has unexpected rovibrational symmetry {symmetry!r}"
            )
        if k == 0 and symmetry != "E1":
            raise BasecolFormatError("E-CH3OH K=0 source level must have E1 symmetry")
        # Rabli/Flower use the E1/E2 convention with a non-negative printed K.
        # Comparison with the source term values fixes E1 -> +K, E2 -> -K.
        signed_k = -k if symmetry == "E2" else k
        return LevelKey("E", vt, j, signed_k)
    raise BasecolFormatError(f"species must be A or E, found {species!r}")


def read_basecol_levels(path: str | Path, species: str) -> tuple[BasecolLevel, ...]:
    """Read one BASECOL ``.lev`` file and map all rows to canonical keys."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    count_labels = [
        index for index, line in enumerate(lines)
        if line.strip().casefold() == "number of lines :"
    ]
    if len(count_labels) != 1:
        raise BasecolFormatError("level file must contain one 'Number of Lines :' label")
    expected_count = _value_after(lines, count_labels[0], context="level count")
    header_indices = [
        index for index, line in enumerate(lines)
        if line.split() == ["Level", "Energy", "J", "K", "rovibSym", "vt"]
    ]
    if len(header_indices) != 1:
        raise BasecolFormatError("cannot locate the six-column BASECOL level header")

    rows: list[BasecolLevel] = []
    for line_number, line in enumerate(lines[header_indices[0] + 1 :], start=header_indices[0] + 2):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 6:
            raise BasecolFormatError(
                f"{path.name}:{line_number}: expected 6 level fields, found {len(fields)}"
            )
        try:
            source_index = int(fields[0])
            j = int(fields[2])
            k = int(fields[3])
            vt = int(fields[5])
        except ValueError as exc:
            raise BasecolFormatError(
                f"{path.name}:{line_number}: malformed integer level field"
            ) from exc
        energy = _finite_float(fields[1], context=f"{path.name}:{line_number} energy")
        if source_index <= 0 or energy < 0 or j < 0 or k < 0:
            raise BasecolFormatError(f"{path.name}:{line_number}: invalid level value")
        try:
            key = _canonical_key(
                species, vt=vt, j=j, k=k, symmetry=fields[4]
            )
        except ValueError as exc:
            raise BasecolFormatError(f"{path.name}:{line_number}: {exc}") from exc
        rows.append(BasecolLevel(source_index, energy, key, fields[4]))

    if len(rows) != expected_count:
        raise BasecolFormatError(
            f"level row count {len(rows)} does not match declared {expected_count}"
        )
    indices = [row.source_index for row in rows]
    if indices != list(range(1, expected_count + 1)):
        raise BasecolFormatError("level indices are not consecutive from 1")
    if len({row.key for row in rows}) != len(rows):
        raise BasecolFormatError("two BASECOL rows map to the same canonical level identity")
    if any(right.source_energy_cm1 < left.source_energy_cm1 for left, right in zip(rows, rows[1:])):
        raise BasecolFormatError("source level energies decrease with BASECOL index")
    return tuple(rows)


def _indices_matching(lines: Iterable[str], label: str) -> list[int]:
    normalized = label.casefold()
    return [
        index for index, line in enumerate(lines)
        if line.strip().casefold() == normalized
    ]


def read_basecol_rates(path: str | Path) -> BasecolRateTable:
    """Read and fully validate one BASECOL ``.kij`` state-to-state table."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    if not any(line.strip().startswith("BASECOL V1") for line in lines[:5]):
        raise BasecolFormatError("rate file does not declare BASECOL V1")

    temperature_labels = _indices_matching(lines, "number of temperatures :")
    transition_labels = _indices_matching(lines, "number of transitions :")
    initial_labels = _indices_matching(lines, "maximum  initial level value :")
    final_labels = _indices_matching(lines, "maximum  final level value :")
    if not (
        len(temperature_labels) == len(transition_labels) == 1
        and len(initial_labels) == len(final_labels) == 2
    ):
        raise BasecolFormatError("rate file has an unexpected BASECOL header layout")
    temperature_count = _value_after(
        lines, temperature_labels[0], context="temperature count"
    )
    transition_count = _value_after(
        lines, transition_labels[0], context="transition count"
    )
    target_max_initial = _value_after(
        lines, initial_labels[0], context="target maximum initial level"
    )
    target_max_final = _value_after(
        lines, final_labels[0], context="target maximum final level"
    )
    collider_max_initial = _value_after(
        lines, initial_labels[1], context="collider maximum initial level"
    )
    collider_max_final = _value_after(
        lines, final_labels[1], context="collider maximum final level"
    )
    if collider_max_initial != 1 or collider_max_final != 1:
        raise BasecolFormatError("expected the structureless He collider at level 1 -> 1")

    order_indices = [
        index for index, line in enumerate(lines)
        if line.strip().casefold().startswith("order :")
    ]
    if len(order_indices) != 1:
        raise BasecolFormatError("rate file must contain one 'order :' declaration")
    temperature_line_index, temperature_line = _next_nonempty(lines, order_indices[0] + 1)
    temperatures = tuple(
        _finite_float(value, context="collision temperature")
        for value in temperature_line.split()
    )
    if len(temperatures) != temperature_count:
        raise BasecolFormatError(
            f"temperature grid has {len(temperatures)} values, expected {temperature_count}"
        )
    if any(value <= 0 for value in temperatures) or any(
        right <= left for left, right in zip(temperatures, temperatures[1:])
    ):
        raise BasecolFormatError("collision temperatures must be positive and strictly increasing")

    rows: list[BasecolRate] = []
    seen_pairs: set[tuple[int, int]] = set()
    for line_number, line in enumerate(lines[temperature_line_index + 1 :], start=temperature_line_index + 2):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 4 + temperature_count:
            raise BasecolFormatError(
                f"{path.name}:{line_number}: expected {4 + temperature_count} fields, found {len(fields)}"
            )
        try:
            initial, final, collider_initial, collider_final = map(int, fields[:4])
        except ValueError as exc:
            raise BasecolFormatError(
                f"{path.name}:{line_number}: malformed state index"
            ) from exc
        if collider_initial != 1 or collider_final != 1:
            raise BasecolFormatError(
                f"{path.name}:{line_number}: collider transition is not 1 -> 1"
            )
        if not (1 <= final < initial <= target_max_initial):
            raise BasecolFormatError(
                f"{path.name}:{line_number}: target transition {initial} -> {final} is not downward by source index"
            )
        rates = tuple(
            _finite_float(value, context=f"{path.name}:{line_number} rate")
            for value in fields[4:]
        )
        if any(value < 0 for value in rates):
            raise BasecolFormatError(f"{path.name}:{line_number}: negative rate coefficient")
        pair = (initial, final)
        if pair in seen_pairs:
            raise BasecolFormatError(f"duplicate BASECOL transition {initial} -> {final}")
        seen_pairs.add(pair)
        rows.append(BasecolRate(initial, final, rates))

    if len(rows) != transition_count:
        raise BasecolFormatError(
            f"rate row count {len(rows)} does not match declared {transition_count}"
        )
    if target_max_final != target_max_initial - 1:
        raise BasecolFormatError("target maximum final level is inconsistent with a complete pair table")
    expected_pairs = {
        (initial, final)
        for initial in range(2, target_max_initial + 1)
        for final in range(1, initial)
    }
    if transition_count != len(expected_pairs) or seen_pairs != expected_pairs:
        missing = len(expected_pairs - seen_pairs)
        extra = len(seen_pairs - expected_pairs)
        raise BasecolFormatError(
            f"state-to-state table is not complete: {missing} pairs missing, {extra} extra"
        )
    return BasecolRateTable(
        temperatures, tuple(rows), target_max_initial, target_max_final
    )


def validate_level_rate_pair(
    levels: tuple[BasecolLevel, ...], rates: BasecolRateTable
) -> None:
    """Assert that a ``.lev``/``.kij`` pair spans the same target basis."""
    if len(levels) != rates.maximum_initial_level:
        raise BasecolFormatError(
            f"level table has {len(levels)} rows but rate table addresses "
            f"{rates.maximum_initial_level} levels"
        )
