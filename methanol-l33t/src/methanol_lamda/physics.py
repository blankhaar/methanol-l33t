"""Unit conversions and rate transformations used by the build.

Constants are written out here so every numerical convention is reviewable.
"""

from __future__ import annotations

import bisect
import math

CM1_TO_K = 1.4387768775039338
MHZ_TO_CM1 = 1.0 / 29979.2458
SMU2_TO_A = 1.16395e-20


def einstein_a_from_smu2(frequency_mhz: float, smu2_debye2: float, g_upper: float) -> float:
    """Return A_ul/s^-1 for frequency/MHz and S*mu^2/Debye^2.

    ``g_upper`` must use the same statistical-weight convention as the supplied
    line strength.  In particular, normalize any CDMS special-catalog nuclear-
    spin factor before calling this function; those products are not assumed to
    share the ordinary catalog convention.
    """
    if frequency_mhz <= 0 or smu2_debye2 <= 0 or g_upper <= 0:
        raise ValueError("frequency, S mu^2, and upper degeneracy must be positive")
    return SMU2_TO_A * frequency_mhz**3 * smu2_debye2 / g_upper


def einstein_a_from_catalog_intensity(
    frequency_mhz: float,
    intensity_nm2_mhz: float,
    lower_energy_cm1: float,
    g_upper: float,
    partition_function: float,
    temperature_k: float = 300.0,
) -> float:
    """Convert ordinary JPL/CDMS integrated intensity to Einstein A.

    Implements the exact JPL catalog equation.  It is intentionally not used
    for starred/special-strength CDMS catalogs whose third column has different
    semantics.
    """
    values = (frequency_mhz, intensity_nm2_mhz, g_upper, partition_function, temperature_k)
    if any(value <= 0 or not math.isfinite(value) for value in values):
        raise ValueError("catalog conversion inputs must be finite and positive")
    upper_energy_cm1 = lower_energy_cm1 + frequency_mhz * MHZ_TO_CM1
    denominator = math.exp(-lower_energy_cm1 * CM1_TO_K / temperature_k) - math.exp(
        -upper_energy_cm1 * CM1_TO_K / temperature_k
    )
    if denominator <= 0:
        raise ValueError("non-positive Boltzmann population difference")
    return (
        2.7964e-16
        * intensity_nm2_mhz
        * frequency_mhz**2
        * (partition_function / g_upper)
        / denominator
    )


def upward_from_downward(
    downward_rate: float,
    g_upper: float,
    g_lower: float,
    delta_energy_cm1: float,
    temperature_k: float,
) -> float:
    """Detailed balance: k_lu from canonical downward k_ul."""
    _check_rate_balance_inputs(downward_rate, g_upper, g_lower, delta_energy_cm1, temperature_k)
    return downward_rate * (g_upper / g_lower) * math.exp(
        -delta_energy_cm1 * CM1_TO_K / temperature_k
    )


def downward_from_upward(
    upward_rate: float,
    g_upper: float,
    g_lower: float,
    delta_energy_cm1: float,
    temperature_k: float,
) -> float:
    """Detailed balance: recover canonical downward k_ul from k_lu."""
    _check_rate_balance_inputs(upward_rate, g_upper, g_lower, delta_energy_cm1, temperature_k)
    return upward_rate * (g_lower / g_upper) * math.exp(
        delta_energy_cm1 * CM1_TO_K / temperature_k
    )


def _check_rate_balance_inputs(
    rate: float, g_upper: float, g_lower: float, delta_energy_cm1: float, temperature_k: float
) -> None:
    if rate < 0 or not math.isfinite(rate):
        raise ValueError("rate must be finite and non-negative")
    if g_upper <= 0 or g_lower <= 0 or delta_energy_cm1 <= 0 or temperature_k <= 0:
        raise ValueError("weights, energy gap, and temperature must be positive")


def interpolate_rate(
    source_temperatures: tuple[float, ...],
    source_rates: tuple[float, ...],
    target_temperature: float,
) -> float:
    """Interpolate a rate without extrapolation.

    Positive rates use log(rate)-log(T), as is conventional for collision-rate
    tables.  An interval touching zero falls back to linear interpolation so
    an actual calculated zero is preserved.
    """
    if len(source_temperatures) != len(source_rates) or not source_temperatures:
        raise ValueError("temperature/rate grid mismatch")
    if any(b <= a for a, b in zip(source_temperatures, source_temperatures[1:])):
        raise ValueError("source temperatures must increase strictly")
    if target_temperature < source_temperatures[0] or target_temperature > source_temperatures[-1]:
        raise ValueError("refusing collision-rate extrapolation")
    pos = bisect.bisect_left(source_temperatures, target_temperature)
    if pos < len(source_temperatures) and source_temperatures[pos] == target_temperature:
        return source_rates[pos]
    left = pos - 1
    right = pos
    t0, t1 = source_temperatures[left], source_temperatures[right]
    r0, r1 = source_rates[left], source_rates[right]
    if r0 < 0 or r1 < 0:
        raise ValueError("rates must be non-negative")
    if r0 == 0 or r1 == 0:
        fraction = (target_temperature - t0) / (t1 - t0)
        return r0 + fraction * (r1 - r0)
    fraction = math.log(target_temperature / t0) / math.log(t1 / t0)
    return math.exp(math.log(r0) + fraction * math.log(r1 / r0))


def interpolate_or_sqrt_extrapolate_rate(
    source_temperatures: tuple[float, ...],
    source_rates: tuple[float, ...],
    target_temperature: float,
) -> float:
    """Interpolate in range and extrapolate from an endpoint as ``sqrt(T)``.

    This is the explicit high-temperature policy for exploratory maser grids.
    It preserves the calculated value at the nearest tabulated endpoint and
    applies only the weak kinetic scaling ``k(T) = k(T_edge) sqrt(T/T_edge)``
    outside the source grid.  It does not manufacture missing state-to-state
    transitions.
    """
    if len(source_temperatures) != len(source_rates) or not source_temperatures:
        raise ValueError("temperature/rate grid mismatch")
    if any(b <= a for a, b in zip(source_temperatures, source_temperatures[1:])):
        raise ValueError("source temperatures must increase strictly")
    if target_temperature <= 0 or not math.isfinite(target_temperature):
        raise ValueError("target temperature must be finite and positive")
    if any(rate < 0 or not math.isfinite(rate) for rate in source_rates):
        raise ValueError("rates must be finite and non-negative")
    if source_temperatures[0] <= target_temperature <= source_temperatures[-1]:
        return interpolate_rate(source_temperatures, source_rates, target_temperature)
    edge = 0 if target_temperature < source_temperatures[0] else -1
    return source_rates[edge] * math.sqrt(
        target_temperature / source_temperatures[edge]
    )
