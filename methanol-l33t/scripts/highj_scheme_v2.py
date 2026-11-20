#!/usr/bin/env python3
"""Version 2 exact-branch high-J collision analogues, never a fitted rate set.

The caller supplies a masked donor network. This module cannot access withheld
rates. It keeps signed E K and A parity, never changes species/manifold/collider,
and never introduces an intermanifold edge. A donor coefficient is transferred
in the mapped direction; when source energies reverse the mapped direction its
upward coefficient is obtained by source detailed balance first. The consumer
must derive the target reverse rate using TARGET energies and weights.

For nearest_shift, a nonzero target delta J requires a nonzero donor delta J
with the same sign. If the exact jump is unavailable, an edge-anchored branch
donor uses the closest available smaller nonzero jump and its exact jump ratio.
The cragg_edge comparator retains the frozen v1 convention if explicitly used.
These are candidate assumptions requiring held-out and excitation tests.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

HC_OVER_K_CM_K = 1.4387768775039338


@dataclass(frozen=True, order=True)
class Key:
    species: str
    manifold: str
    j: int
    k: int
    parity: str = ""

    def with_j(self, j: int) -> "Key | None":
        if j < abs(self.k) or j < 0:
            return None
        return Key(self.species, self.manifold, j, self.k, self.parity)

    def label(self) -> str:
        return f"{self.species}:{self.manifold}:J{self.j}:K{self.k:+d}:{self.parity}"


@dataclass(frozen=True)
class State:
    key: Key
    energy_cm1: float
    weight: float


@dataclass(frozen=True)
class Rate:
    upper: Key
    lower: Key
    values: tuple[float, ...]


@dataclass(frozen=True)
class DonorChoice:
    mapped_upper: Key
    mapped_lower: Key
    scale: float
    rule: str
    reverse_source_direction: bool
    donor_upper: Key
    donor_lower: Key


@dataclass
class SourceNetwork:
    states: Mapping[Key, State]
    rates: Mapping[frozenset[Key], Rate]
    temperatures: tuple[float, ...]
    collider: str


def _valid_choice(
    upper: Key, lower: Key, a: Key | None, b: Key | None,
    network: SourceNetwork, edge: int, rule: str, scale_delta_j: bool,
    nonzero_same_signed_delta_j: bool = False,
) -> DonorChoice | None:
    if a is None or b is None or a == b or max(a.j, b.j) > edge:
        return None
    if a not in network.states or b not in network.states:
        return None
    donor = network.rates.get(frozenset((a, b)))
    if donor is None:
        return None
    target_dj, donor_dj = abs(upper.j - lower.j), abs(a.j - b.j)
    if nonzero_same_signed_delta_j and target_dj > 0:
        if donor_dj == 0 or (upper.j-lower.j)*(a.j-b.j) <= 0:
            return None
    # delta-J=0 is a separate branch rule, not division by zero.
    if target_dj == 0:
        if donor_dj != 0:
            return None
        scale = 1.0
    else:
        scale = max(1, donor_dj) / target_dj if scale_delta_j else 1.0
    return DonorChoice(a, b, scale, rule, donor.upper != a,
                       donor.upper, donor.lower)


def select_donor(
    upper: Key, lower: Key, network: SourceNetwork, edge: int,
    method: str = "nearest_shift",
) -> DonorChoice | None:
    """Select an analogue without using target coefficients or target energies.

    Inputs upper/lower are ordered by TARGET energy, not by J. The network must
    already be masked to j<=edge; a hard assertion prevents donor leakage.
    This function is for absent pairs only. Existing direct coefficients should
    always be preserved by the calling builder.
    """
    if method not in {"nearest_shift", "cragg_edge"}:
        raise ValueError(f"unknown method: {method}")
    if (upper.species, upper.manifold) != (lower.species, lower.manifold):
        raise ValueError("cross-species/intermanifold extrapolation prohibited")
    if upper == lower or max(upper.j, lower.j) <= edge:
        return None
    if abs(upper.k) > edge or abs(lower.k) > edge:
        return None  # No same-K ladder exists within the masked source.
    shift = max(upper.j, lower.j) - edge
    shifted = (upper.with_j(upper.j - shift), lower.with_j(lower.j - shift))
    clamped = (upper.with_j(min(upper.j, edge)), lower.with_j(min(lower.j, edge)))
    exactly_one_high = (upper.j > edge) != (lower.j > edge)
    proposals = []
    if method == "cragg_edge" and exactly_one_high:
        proposals.append((*clamped, "fixed_low_endpoint_propensity", True))
    proposals.append((*shifted, "same_delta_j_shift", False))
    proposals.append((*clamped, "clamped_branch_propensity", True))
    for a, b, rule, scaled in proposals:
        choice = _valid_choice(upper, lower, a, b, network, edge, rule, scaled,
                               nonzero_same_signed_delta_j=method == "nearest_shift")
        if choice is not None:
            if not (max(choice.donor_upper.j, choice.donor_lower.j) <= edge):
                raise RuntimeError('Invariant failed in highj_scheme_v2.py: max(choice.donor_upper.j, choice.donor_lower.j) <= edge')
            if not (choice.mapped_upper.k == upper.k):
                raise RuntimeError('Invariant failed in highj_scheme_v2.py: choice.mapped_upper.k == upper.k')
            if not (choice.mapped_lower.k == lower.k):
                raise RuntimeError('Invariant failed in highj_scheme_v2.py: choice.mapped_lower.k == lower.k')
            if not (choice.mapped_upper.parity == upper.parity):
                raise RuntimeError('Invariant failed in highj_scheme_v2.py: choice.mapped_upper.parity == upper.parity')
            if not (choice.mapped_lower.parity == lower.parity):
                raise RuntimeError('Invariant failed in highj_scheme_v2.py: choice.mapped_lower.parity == lower.parity')
            return choice
    if method == "nearest_shift" and upper.j != lower.j:
        signed_jump = upper.j-lower.j
        # The largest donor J is strictly the specified edge. No additional
        # lower-anchor assumption or coefficient-based search is introduced.
        for jump in range(min(abs(signed_jump), edge), 0, -1):
            ja, jb = ((edge, edge-jump) if signed_jump > 0
                      else (edge-jump, edge))
            a, b = upper.with_j(ja), lower.with_j(jb)
            choice = _valid_choice(upper, lower, a, b, network, edge,
                "compressed_nonzero_delta_j", True,
                nonzero_same_signed_delta_j=True)
            if choice is not None:
                if not (0 < abs(a.j-b.j) <= abs(signed_jump)):
                    raise RuntimeError('Invariant failed in highj_scheme_v2.py: 0 < abs(a.j-b.j) <= abs(signed_jump)')
                if not ((a.j-b.j)*signed_jump > 0):
                    raise RuntimeError('Invariant failed in highj_scheme_v2.py: (a.j-b.j)*signed_jump > 0')
                if not (max(a.j,b.j) == edge):
                    raise RuntimeError('Invariant failed in highj_scheme_v2.py: max(a.j,b.j) == edge')
                return choice
    return None


def mapped_rates(choice: DonorChoice, network: SourceNetwork) -> tuple[float, ...]:
    """Transfer the complete donor curve with explicit source-direction handling."""
    rate = network.rates[frozenset((choice.donor_upper, choice.donor_lower))]
    if len(rate.values) != len(network.temperatures):
        raise ValueError("temperature/rate length mismatch")
    factors = [choice.scale] * len(network.temperatures)
    if choice.reverse_source_direction:
        high = network.states[rate.upper]
        low = network.states[rate.lower]
        gap = high.energy_cm1 - low.energy_cm1
        if gap < -1e-10:
            raise ValueError("source downward rate is not energy ordered")
        factors = [choice.scale * high.weight / low.weight *
                   math.exp(-HC_OVER_K_CM_K * max(0.0, gap) / t)
                   for t in network.temperatures]
    values = tuple(value * factor for value, factor in zip(rate.values, factors))
    if any(not math.isfinite(x) or x < 0 for x in values):
        raise ValueError("invalid extrapolated coefficient")
    return values


def target_upward_rates(
    downward: Sequence[float], temperatures: Sequence[float],
    upper: State, lower: State,
) -> tuple[float, ...]:
    gap = upper.energy_cm1 - lower.energy_cm1
    if gap < -1e-10 or upper.weight <= 0 or lower.weight <= 0:
        raise ValueError("invalid target energy ordering or weights")
    return tuple(k * upper.weight / lower.weight *
                 math.exp(-HC_OVER_K_CM_K * max(0.0, gap) / t)
                 for k, t in zip(downward, temperatures))
