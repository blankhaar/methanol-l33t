"""Small immutable-ish data models used by the parsers and writer."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True, frozen=True)
class Level:
    index: int
    energy_cm1: float
    weight: float
    quantum_label: str


@dataclass(slots=True, frozen=True)
class RadiativeTransition:
    index: int
    upper: int
    lower: int
    a_s: float
    frequency_ghz: float
    upper_energy_k: float | None = None


@dataclass(slots=True, frozen=True)
class CollisionTransition:
    index: int
    upper: int
    lower: int
    rates_cm3_s: tuple[float, ...]


@dataclass(slots=True)
class CollisionPartner:
    collider_id: int
    description: str
    temperatures_k: tuple[float, ...]
    transitions: list[CollisionTransition] = field(default_factory=list)


@dataclass(slots=True)
class LamdaData:
    molecule: str
    mass_amu: float
    levels: list[Level] = field(default_factory=list)
    radiative: list[RadiativeTransition] = field(default_factory=list)
    collision_partners: list[CollisionPartner] = field(default_factory=list)

