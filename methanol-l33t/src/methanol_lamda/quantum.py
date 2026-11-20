"""Canonical methanol level identities.

E levels retain signed K.  A levels retain non-negative K plus an explicit
parity label.  Energies and source row numbers never participate in identity.
"""

from __future__ import annotations

from dataclasses import dataclass
import re


@dataclass(frozen=True, slots=True, order=True)
class LevelKey:
    species: str
    vt: int
    j: int
    k: int
    parity: str = ""

    def __post_init__(self) -> None:
        species = self.species.upper()
        if species not in {"A", "E"}:
            raise ValueError("species must be A or E")
        object.__setattr__(self, "species", species)
        if self.vt not in {0, 1, 2, 3, 4} or self.j < 0 or abs(self.k) > self.j:
            raise ValueError("invalid methanol quantum numbers")
        if species == "A":
            if self.k < 0:
                raise ValueError("canonical A-state K must be non-negative")
            if self.parity not in {"+", "-"}:
                raise ValueError("A state requires explicit + or - parity")
        elif self.parity:
            raise ValueError("E state uses signed K and no parity field")

    @property
    def rotational_identity(self) -> tuple[int, int, str]:
        return self.j, self.k, self.parity

    def encode(self) -> str:
        if self.species == "A":
            return f"A:vt{self.vt}:J{self.j}:K{self.k}:{self.parity}"
        return f"E:vt{self.vt}:J{self.j}:K{self.k:+d}"

    @classmethod
    def decode(cls, value: str) -> "LevelKey":
        match = re.fullmatch(r"A:vt([0-4]):J(\d+):K(\d+):([+-])", value)
        if match:
            return cls("A", int(match[1]), int(match[2]), int(match[3]), match[4])
        match = re.fullmatch(r"E:vt([0-4]):J(\d+):K([+-]\d+)", value)
        if match:
            return cls("E", int(match[1]), int(match[2]), int(match[3]))
        raise ValueError(f"invalid canonical level key: {value!r}")


def from_legacy_lamda_label(species: str, label: str, vt: int = 0) -> LevelKey:
    """Map the historical EMAA ``J_K`` label to an explicit canonical key.

    The legacy A file uses signed K to distinguish the two A parity ladders.
    For K=0 the historical convention is A+.
    """
    match = re.fullmatch(r"\s*(\d+)_([+-]?\d+)\s*", label)
    if not match:
        raise ValueError(f"unsupported legacy quantum label {label!r}")
    j, signed_k = int(match[1]), int(match[2])
    if species.upper() == "A":
        parity = "-" if signed_k < 0 else "+"
        return LevelKey("A", vt, j, abs(signed_k), parity)
    return LevelKey("E", vt, j, signed_k)
