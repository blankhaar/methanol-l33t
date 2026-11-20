"""LAMDA/RADEX text parser, writer, and structural validation.

The parser deliberately does not depend on astroquery so it can also serve as
one of the two independent readers used by the validation workflow.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

from .models import CollisionPartner, CollisionTransition, LamdaData, Level, RadiativeTransition


class LamdaFormatError(ValueError):
    """Raised when a LAMDA file is structurally invalid."""


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


class _Cursor:
    def __init__(self, lines: Iterable[str]) -> None:
        self.lines = list(lines)
        self.position = 0

    def next_data(self) -> str:
        while self.position < len(self.lines):
            raw = self.lines[self.position]
            self.position += 1
            line = raw.strip()
            if line and not line.startswith("!") and not line.startswith("#"):
                return line
        raise LamdaFormatError("unexpected end of file")

    def next_values(self, count: int) -> list[str]:
        values: list[str] = []
        while len(values) < count:
            values.extend(self.next_data().split())
        if len(values) != count:
            raise LamdaFormatError(f"expected {count} values, found {len(values)}")
        return values


def parse_lamda_text(text: str) -> LamdaData:
    cursor = _Cursor(text.splitlines())
    molecule = cursor.next_data()
    mass = _float(cursor.next_data().split()[0])
    nlevels = int(cursor.next_data().split()[0])
    levels: list[Level] = []
    for _ in range(nlevels):
        fields = cursor.next_data().split(maxsplit=3)
        if len(fields) < 4:
            raise LamdaFormatError(f"malformed level row: {fields!r}")
        levels.append(Level(int(fields[0]), _float(fields[1]), _float(fields[2]), fields[3].strip()))

    nrad = int(cursor.next_data().split()[0])
    radiative: list[RadiativeTransition] = []
    for _ in range(nrad):
        fields = cursor.next_data().split()
        if len(fields) < 5:
            raise LamdaFormatError(f"malformed radiative row: {fields!r}")
        eu = _float(fields[5]) if len(fields) >= 6 else None
        radiative.append(
            RadiativeTransition(
                index=int(fields[0]),
                upper=int(fields[1]),
                lower=int(fields[2]),
                a_s=_float(fields[3]),
                frequency_ghz=_float(fields[4]),
                upper_energy_k=eu,
            )
        )

    npartners = int(cursor.next_data().split()[0])
    partners: list[CollisionPartner] = []
    for _ in range(npartners):
        identity = cursor.next_data().split(maxsplit=1)
        collider_id = int(identity[0])
        description = identity[1] if len(identity) == 2 else ""
        ntrans = int(cursor.next_data().split()[0])
        ntemps = int(cursor.next_data().split()[0])
        temperatures = tuple(_float(value) for value in cursor.next_values(ntemps))
        transitions: list[CollisionTransition] = []
        for _ in range(ntrans):
            fields = cursor.next_data().split()
            if len(fields) != 3 + ntemps:
                raise LamdaFormatError(
                    f"collision row has {len(fields) - 3} rates, expected {ntemps}: {fields[:4]!r}"
                )
            transitions.append(
                CollisionTransition(
                    index=int(fields[0]),
                    upper=int(fields[1]),
                    lower=int(fields[2]),
                    rates_cm3_s=tuple(_float(value) for value in fields[3:]),
                )
            )
        partners.append(CollisionPartner(collider_id, description, temperatures, transitions))
    return LamdaData(molecule, mass, levels, radiative, partners)


def read_lamda(path: str | Path) -> LamdaData:
    return parse_lamda_text(Path(path).read_text(encoding="utf-8", errors="strict"))


def _fmt_rate(value: float) -> str:
    return f"{value:.8e}"


def lamda_text(data: LamdaData) -> str:
    out = [
        "!MOLECULE",
        data.molecule,
        "!MOLECULAR WEIGHT",
        f"{data.mass_amu:.8g}",
        "!NUMBER OF ENERGY LEVELS",
        str(len(data.levels)),
        "!LEVEL + ENERGY(CM-1) + WEIGHT + QUANTUM NOS.",
    ]
    for level in data.levels:
        # BELGI resolves some high-K A+/A- doublets below 1e-9 cm-1.  Retain
        # enough digits that serialization cannot collapse or reverse them.
        out.append(
            f"{level.index:6d} {level.energy_cm1:.17g} "
            f"{level.weight:12.6f} {level.quantum_label}"
        )
    out.extend(
        [
            "!NUMBER OF RADIATIVE TRANSITIONS",
            str(len(data.radiative)),
            "!TRANS + UP + LOW + EINSTEIN_A(S-1) + FREQ(GHZ) + E_UPPER(K)",
        ]
    )
    level_by_index = {level.index: level for level in data.levels}
    cm1_to_k = 1.438776877
    for transition in data.radiative:
        eu = transition.upper_energy_k
        if eu is None:
            eu = level_by_index[transition.upper].energy_cm1 * cm1_to_k
        out.append(
            f"{transition.index:6d} {transition.upper:6d} {transition.lower:6d} "
            f"{transition.a_s:.10e} {transition.frequency_ghz:.15g} {eu:.12g}"
        )
    out.extend(["!NUMBER OF COLLISION PARTNERS", str(len(data.collision_partners))])
    for partner in data.collision_partners:
        out.extend(
            [
                "!COLLISIONS BETWEEN",
                f"{partner.collider_id} {partner.description}".rstrip(),
                "!NUMBER OF COLLISIONAL TRANSITIONS",
                str(len(partner.transitions)),
                "!NUMBER OF COLLISION TEMPERATURES",
                str(len(partner.temperatures_k)),
                "!COLLISION TEMPERATURES",
                " ".join(f"{temperature:.8g}" for temperature in partner.temperatures_k),
                "!TRANS + UP + LOW + RATE COEFFS(CM3 S-1)",
            ]
        )
        for transition in partner.transitions:
            rates = " ".join(_fmt_rate(rate) for rate in transition.rates_cm3_s)
            out.append(f"{transition.index:7d} {transition.upper:6d} {transition.lower:6d} {rates}")
    return "\n".join(out) + "\n"


def write_lamda(path: str | Path, data: LamdaData) -> None:
    Path(path).write_text(lamda_text(data), encoding="utf-8")


def validate_lamda(data: LamdaData, *, require_downward: bool = True) -> list[str]:
    """Return validation errors; an empty list means the structure passed."""
    errors: list[str] = []
    if not data.molecule.strip():
        errors.append("molecule name is empty")
    if not math.isfinite(data.mass_amu) or data.mass_amu <= 0:
        errors.append("molecular mass is not finite and positive")
    indices = [level.index for level in data.levels]
    if indices != list(range(1, len(indices) + 1)):
        errors.append("level indices are not consecutive from 1")
    energies = [level.energy_cm1 for level in data.levels]
    if any(not math.isfinite(value) for value in energies):
        errors.append("one or more level energies are not finite")
    if any(right <= left for left, right in zip(energies, energies[1:])):
        errors.append("level energies do not increase strictly with file index")
    if any(not math.isfinite(level.weight) or level.weight <= 0 for level in data.levels):
        errors.append("one or more statistical weights are not finite and positive")
    level_by_index = {level.index: level for level in data.levels}

    rad_indices = [transition.index for transition in data.radiative]
    if rad_indices != list(range(1, len(rad_indices) + 1)):
        errors.append("radiative transition indices are not consecutive from 1")
    seen_radiative: set[tuple[int, int]] = set()
    for transition in data.radiative:
        if transition.upper not in level_by_index or transition.lower not in level_by_index:
            errors.append(f"radiative transition {transition.index} has unresolved level index")
            continue
        if not math.isfinite(transition.a_s) or transition.a_s <= 0:
            errors.append(f"radiative transition {transition.index} has non-positive Einstein A")
        if not math.isfinite(transition.frequency_ghz) or transition.frequency_ghz <= 0:
            errors.append(f"radiative transition {transition.index} has non-positive frequency")
        if require_downward and level_by_index[transition.upper].energy_cm1 <= level_by_index[transition.lower].energy_cm1:
            errors.append(f"radiative transition {transition.index} is not downward")
        identity = (transition.upper, transition.lower)
        if identity in seen_radiative:
            errors.append(f"radiative transition {transition.index} duplicates endpoints {identity}")
        seen_radiative.add(identity)

    seen_colliders: set[int] = set()
    for partner in data.collision_partners:
        if partner.collider_id in seen_colliders:
            errors.append(f"collider id {partner.collider_id} occurs more than once")
        seen_colliders.add(partner.collider_id)
        if any(not math.isfinite(t) or t <= 0 for t in partner.temperatures_k):
            errors.append(f"collider {partner.collider_id} has invalid temperatures")
        if any(right <= left for left, right in zip(partner.temperatures_k, partner.temperatures_k[1:])):
            errors.append(f"collider {partner.collider_id} temperature grid is not strictly increasing")
        collision_indices = [transition.index for transition in partner.transitions]
        if collision_indices != list(range(1, len(collision_indices) + 1)):
            errors.append(f"collider {partner.collider_id} transition indices are not consecutive from 1")
        seen_collisions: set[tuple[int, int]] = set()
        for transition in partner.transitions:
            if transition.upper not in level_by_index or transition.lower not in level_by_index:
                errors.append(
                    f"collider {partner.collider_id} transition {transition.index} has unresolved level index"
                )
                continue
            if len(transition.rates_cm3_s) != len(partner.temperatures_k):
                errors.append(f"collider {partner.collider_id} transition {transition.index} grid mismatch")
            if any(not math.isfinite(rate) or rate < 0 for rate in transition.rates_cm3_s):
                errors.append(f"collider {partner.collider_id} transition {transition.index} has invalid rate")
            if require_downward and level_by_index[transition.upper].energy_cm1 <= level_by_index[transition.lower].energy_cm1:
                errors.append(f"collider {partner.collider_id} transition {transition.index} is not downward")
            identity = (transition.upper, transition.lower)
            if identity in seen_collisions:
                errors.append(
                    f"collider {partner.collider_id} transition {transition.index} duplicates endpoints {identity}"
                )
            seen_collisions.add(identity)
    return errors
