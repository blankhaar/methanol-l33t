"""A deliberately separate, minimal LAMDA reader used for cross-validation.

It does not import or call :mod:`methanol_lamda.lamda`; agreement therefore
catches a useful class of writer/parser state-machine mistakes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class IndependentSummary:
    molecule: str
    mass_amu: float
    level_count: int
    radiative_count: int
    collider_counts: tuple[tuple[int, int, int], ...]


def read_summary(path: str | Path) -> IndependentSummary:
    raw = Path(path).read_text(encoding="utf-8")
    rows = [
        line.strip()
        for line in raw.splitlines()
        if line.strip() and not line.lstrip().startswith(("!", "#"))
    ]
    cursor = 0

    def take() -> str:
        nonlocal cursor
        if cursor >= len(rows):
            raise ValueError("independent reader reached unexpected EOF")
        value = rows[cursor]
        cursor += 1
        return value

    molecule = take()
    mass = float(take().split()[0].replace("D", "E"))
    level_count = int(take().split()[0])
    level_ids = []
    for _ in range(level_count):
        fields = take().split()
        if len(fields) < 4:
            raise ValueError("malformed level in independent reader")
        level_ids.append(int(fields[0]))
    radiative_count = int(take().split()[0])
    for expected in range(1, radiative_count + 1):
        fields = take().split()
        if len(fields) < 5 or int(fields[0]) != expected:
            raise ValueError("malformed/non-consecutive radiative row")
        if int(fields[1]) not in level_ids or int(fields[2]) not in level_ids:
            raise ValueError("unresolved radiative level")
    partner_count = int(take().split()[0])
    colliders: list[tuple[int, int, int]] = []
    for _ in range(partner_count):
        collider_id = int(take().split()[0])
        transition_count = int(take().split()[0])
        temperature_count = int(take().split()[0])
        temperatures: list[str] = []
        while len(temperatures) < temperature_count:
            temperatures.extend(take().split())
        if len(temperatures) != temperature_count:
            raise ValueError("collision temperature count mismatch")
        for expected in range(1, transition_count + 1):
            fields = take().split()
            if len(fields) != 3 + temperature_count or int(fields[0]) != expected:
                raise ValueError("malformed/non-consecutive collision row")
            if int(fields[1]) not in level_ids or int(fields[2]) not in level_ids:
                raise ValueError("unresolved collision level")
        colliders.append((collider_id, transition_count, temperature_count))
    if cursor != len(rows):
        raise ValueError(f"independent reader found {len(rows) - cursor} trailing data rows")
    return IndependentSummary(molecule, mass, level_count, radiative_count, tuple(colliders))
