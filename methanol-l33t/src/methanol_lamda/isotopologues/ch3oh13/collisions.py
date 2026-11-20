"""Reduced-mass collision-rate proxies for :sup:`13`CH3OH.

The scaling implemented here changes only the reduced-mass prefactor in the
thermal rate integral.  It assumes that the state-to-state cross-section
integral is unchanged by carbon-isotope substitution and compares rates at the
same kinetic temperature::

    k_13,c(T) = sqrt(mu_12,c / mu_13,c) k_12,c(T)

This is an explicit proxy, not a new scattering calculation.  Quantum-state
mapping, detailed-balance reorientation, interpolation, and source precedence
belong to the isotope-specific product builder; this module only supplies the
mass transform and immutable provenance describing it.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, fields, is_dataclass, replace
import math
from types import MappingProxyType
from typing import TypeVar, cast


# Neutral-isotope mass sums in unified atomic mass units.  Keeping the numeric
# inputs here, instead of deriving them from rounded LAMDA display masses,
# makes every generated multiplier reproducible.
MASS_12CH3OH_U = 32.026214748490
MASS_13CH3OH_U = 33.029569583560
MASS_H2_U = 2.015650064460
MASS_HE4_U = 4.002603254130

ISOTOPOLOGUE_PROXY_STATUS = "isotopologue_proxy_reduced_mass_scaled"
SCALING_FORMULA = "k13(T) = sqrt(mu12 / mu13) * k12(T)"
SCALING_ASSUMPTION = (
    "same state-to-state cross-section integral at the same kinetic temperature"
)

ColliderReference = int | str
CollisionRowT = TypeVar("CollisionRowT")


@dataclass(frozen=True, slots=True)
class CollisionScalingMetadata:
    """Complete immutable provenance for one collider's isotope transform."""

    collider_id: int
    collider_name: str
    collider_mass_u: float
    source_isotopologue: str
    target_isotopologue: str
    source_target_mass_u: float
    target_target_mass_u: float
    source_reduced_mass_u: float
    target_reduced_mass_u: float
    factor: float
    formula: str
    assumption: str
    temperature_mapping: str
    status: str


def reduced_mass(first_mass_u: float, second_mass_u: float) -> float:
    """Return the two-body reduced mass in unified atomic mass units."""
    try:
        first = float(first_mass_u)
        second = float(second_mass_u)
    except (TypeError, ValueError) as exc:
        raise TypeError("reduced-mass inputs must be real numbers") from exc
    if not math.isfinite(first) or not math.isfinite(second):
        raise ValueError("reduced-mass inputs must be finite")
    if first <= 0.0 or second <= 0.0:
        raise ValueError("reduced-mass inputs must be positive")
    return first * second / (first + second)


_COLLIDERS: Mapping[int, tuple[str, float]] = MappingProxyType(
    {
        2: ("p-H2", MASS_H2_U),
        3: ("o-H2", MASS_H2_U),
        6: ("He", MASS_HE4_U),
    }
)

_COLLIDER_ALIASES: Mapping[str, int] = MappingProxyType(
    {
        "ph2": 2,
        "parah2": 2,
        "h2para": 2,
        "oh2": 3,
        "orthoh2": 3,
        "h2ortho": 3,
        "he": 6,
        "helium": 6,
        "4he": 6,
        "he4": 6,
    }
)


def _collider_id(collider: ColliderReference) -> int:
    if isinstance(collider, bool):
        raise TypeError("collider must be an integer ID or name")
    if isinstance(collider, int):
        if collider not in _COLLIDERS:
            raise ValueError(f"unsupported collider ID {collider}; expected 2, 3, or 6")
        return collider
    if isinstance(collider, str):
        normalized = "".join(character for character in collider.casefold() if character.isalnum())
        try:
            return _COLLIDER_ALIASES[normalized]
        except KeyError as exc:
            raise ValueError(
                f"unsupported collider name {collider!r}; expected p-H2, o-H2, or He"
            ) from exc
    raise TypeError("collider must be an integer ID or name")


def _metadata(collider_id: int) -> CollisionScalingMetadata:
    collider_name, collider_mass = _COLLIDERS[collider_id]
    source_mu = reduced_mass(MASS_12CH3OH_U, collider_mass)
    target_mu = reduced_mass(MASS_13CH3OH_U, collider_mass)
    return CollisionScalingMetadata(
        collider_id=collider_id,
        collider_name=collider_name,
        collider_mass_u=collider_mass,
        source_isotopologue="12CH3OH",
        target_isotopologue="13CH3OH",
        source_target_mass_u=MASS_12CH3OH_U,
        target_target_mass_u=MASS_13CH3OH_U,
        source_reduced_mass_u=source_mu,
        target_reduced_mass_u=target_mu,
        factor=math.sqrt(source_mu / target_mu),
        formula=SCALING_FORMULA,
        assumption=SCALING_ASSUMPTION,
        temperature_mapping="same_temperature_no_rescaling",
        status=ISOTOPOLOGUE_PROXY_STATUS,
    )


COLLISION_SCALING_METADATA: Mapping[int, CollisionScalingMetadata] = MappingProxyType(
    {collider_id: _metadata(collider_id) for collider_id in _COLLIDERS}
)


def scaling_metadata(collider: ColliderReference) -> CollisionScalingMetadata:
    """Return immutable mass-scaling provenance for a collider ID or name."""
    return COLLISION_SCALING_METADATA[_collider_id(collider)]


def scaling_factor(collider: ColliderReference) -> float:
    """Return ``sqrt(mu_12 / mu_13)`` for p-H2, o-H2, or He."""
    return scaling_metadata(collider).factor


def scale_rate_tuple(
    rates_cm3_s: Iterable[float], collider: ColliderReference
) -> tuple[float, ...]:
    """Scale a nonempty rate sequence, preserving every exact zero as zero."""
    factor = scaling_factor(collider)
    scaled: list[float] = []
    for index, value in enumerate(rates_cm3_s):
        try:
            rate = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"collision rate {index} must be a real number") from exc
        if not math.isfinite(rate):
            raise ValueError(f"collision rate {index} must be finite")
        if rate < 0.0:
            raise ValueError(f"collision rate {index} must be nonnegative")
        scaled.append(rate if rate == 0.0 else factor * rate)
    if not scaled:
        raise ValueError("collision-rate tuple must not be empty")
    return tuple(scaled)


def scale_collision_row(
    row: CollisionRowT, collider: ColliderReference | None = None
) -> CollisionRowT:
    """Return a rate-scaled copy of a frozen dataclass collision row.

    The row may call its rate field either ``rates`` (the normalized/product
    convention) or ``rates_cm3_s`` (the LAMDA model convention).  If it has a
    ``collider_id`` field, that ID is used by default and an explicitly passed
    collider must agree with it.  All non-rate fields are preserved verbatim.
    """
    if not is_dataclass(row) or isinstance(row, type):
        raise TypeError("collision row must be a frozen dataclass instance")
    parameters = getattr(type(row), "__dataclass_params__", None)
    if parameters is None or not parameters.frozen:
        raise TypeError("collision row must be a frozen dataclass instance")

    field_names = {field.name for field in fields(row)}
    rate_field = (
        "rates" if "rates" in field_names
        else "rates_cm3_s" if "rates_cm3_s" in field_names
        else None
    )
    if rate_field is None:
        raise TypeError("collision row must define a rates or rates_cm3_s field")

    row_collider = getattr(row, "collider_id", None)
    if collider is None:
        if row_collider is None:
            raise ValueError("collider is required when the row has no collider_id field")
        selected_collider = _collider_id(row_collider)
    else:
        selected_collider = _collider_id(collider)
        if row_collider is not None and _collider_id(row_collider) != selected_collider:
            raise ValueError(
                f"explicit collider {selected_collider} conflicts with row collider {row_collider}"
            )

    scaled_rates = scale_rate_tuple(getattr(row, rate_field), selected_collider)
    return cast(CollisionRowT, replace(row, **{rate_field: scaled_rates}))
