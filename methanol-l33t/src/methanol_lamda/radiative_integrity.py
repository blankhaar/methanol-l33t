"""Fail-closed radiative-transition checks for generated LAMDA products.

The audit uses canonical quantum labels rather than product-local level
indices.  It enforces the electric-dipole rank-1 rotational rule, binds base
products to their row-level source provenance, and requires exact radiative
inheritance between successive molecular tiers.  A laboratory frequency may
differ from the adopted Hamiltonian energy gap only when an explicit source
row accounts for it; descendants may inherit that exact source-bound value.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Iterable, Mapping

from .models import LamdaData, RadiativeTransition
from .physics import SMU2_TO_A


CM1_TO_GHZ = 29.9792458
LAMDA_CM1_TO_K = 1.438776877
DEFAULT_MINIMUM_SMU2_DEBYE2 = 1.0e-4
CALCULATED_FREQUENCY_TOLERANCE_GHZ = 1.0e-7
MAXIMUM_SOURCE_FREQUENCY_OFFSET_GHZ = 0.1
SMU2_RELATIVE_TOLERANCE = 1.0e-9

_A_LABEL = re.compile(
    r"^A:(?P<manifold>vt[0-4]|vCO1:vt0|vrock1:vt0):"
    r"J(?P<j>0|[1-9][0-9]*):K(?P<k>0|[1-9][0-9]*):(?P<parity>[+-])$"
)
_E_LABEL = re.compile(
    r"^E:(?P<manifold>vt[0-4]|vCO1:vt0|vrock1:vt0):"
    r"J(?P<j>0|[1-9][0-9]*):K(?P<sign>[+-])(?P<k>0|[1-9][0-9]*)$"
)


class RadiativeIntegrityError(ValueError):
    """A radiative row is forbidden, unbound, or internally inconsistent."""


@dataclass(frozen=True, slots=True)
class _QuantumIdentity:
    raw: str
    species: str
    manifold: str
    j: int
    k: int
    parity: str


def _parse_label(raw: str) -> _QuantumIdentity:
    match = _A_LABEL.fullmatch(raw)
    if match:
        j = int(match.group("j"))
        k = int(match.group("k"))
        parity = match.group("parity")
        if k > j or (k == 0 and parity != "+"):
            raise RadiativeIntegrityError(f"noncanonical A level label {raw!r}")
        return _QuantumIdentity(
            raw, "A", match.group("manifold"), j, k, parity
        )
    match = _E_LABEL.fullmatch(raw)
    if match:
        j = int(match.group("j"))
        magnitude = int(match.group("k"))
        if magnitude > j or (magnitude == 0 and match.group("sign") != "+"):
            raise RadiativeIntegrityError(f"noncanonical E level label {raw!r}")
        k = magnitude if match.group("sign") == "+" else -magnitude
        return _QuantumIdentity(raw, "E", match.group("manifold"), j, k, "")
    raise RadiativeIntegrityError(f"noncanonical methanol level label {raw!r}")


def _source_ledger(
    source_rows: Iterable[Mapping[str, object]],
    *,
    context: str,
) -> dict[tuple[str, str], Mapping[str, object]]:
    ledger: dict[tuple[str, str], Mapping[str, object]] = {}
    for row in source_rows:
        try:
            identity = (str(row["upper_key"]), str(row["lower_key"]))
        except KeyError as exc:
            raise RadiativeIntegrityError(
                f"{context}: radiative source row lacks endpoint identities"
            ) from exc
        reverse = (identity[1], identity[0])
        if identity in ledger or reverse in ledger:
            raise RadiativeIntegrityError(
                f"{context}: duplicate or reversed radiative source row {identity}"
            )
        ledger[identity] = row
    return ledger


def _product_ledger(
    data: LamdaData,
    labels: Mapping[int, _QuantumIdentity],
    *,
    context: str,
) -> dict[tuple[str, str], RadiativeTransition]:
    ledger: dict[tuple[str, str], RadiativeTransition] = {}
    for row in data.radiative:
        try:
            identity = (labels[row.upper].raw, labels[row.lower].raw)
        except KeyError as exc:
            raise RadiativeIntegrityError(
                f"{context}: radiative row {row.index} has an unresolved level index"
            ) from exc
        reverse = (identity[1], identity[0])
        if identity in ledger or reverse in ledger:
            raise RadiativeIntegrityError(
                f"{context}: duplicate or reversed radiative endpoint {identity}"
            )
        ledger[identity] = row
    return ledger


def _inferred_smu2(upper_weight: float, row: RadiativeTransition) -> float:
    return (
        row.a_s
        * upper_weight
        / (SMU2_TO_A * (row.frequency_ghz * 1.0e3) ** 3)
    )


def _radiative_components(data: LamdaData) -> int:
    parent = list(range(len(data.levels) + 1))

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for row in data.radiative:
        upper = root(row.upper)
        lower = root(row.lower)
        if upper != lower:
            parent[lower] = upper
    return len({root(level.index) for level in data.levels})


def audit_radiative_integrity(
    data: LamdaData,
    *,
    context: str,
    source_rows: Iterable[Mapping[str, object]] | None = None,
    parent: LamdaData | None = None,
    minimum_smu2_debye2: float = DEFAULT_MINIMUM_SMU2_DEBYE2,
) -> dict[str, object]:
    """Audit one generated radiative network and raise on any hard failure.

    Exactly one of ``source_rows`` and ``parent`` may be supplied.  Source rows
    must reproduce the complete product graph, its A/frequency values, and its
    S mu^2 values.  A parent must be an exact labelled subset; its complete
    radiative graph and A/frequency tokens must survive unchanged.
    """

    if source_rows is not None and parent is not None:
        raise RadiativeIntegrityError(
            f"{context}: source rows and parent inheritance are mutually exclusive"
        )
    if (
        not math.isfinite(minimum_smu2_debye2)
        or minimum_smu2_debye2 <= 0.0
    ):
        raise RadiativeIntegrityError(
            f"{context}: minimum S mu^2 must be finite and positive"
        )
    if not data.levels or not data.radiative:
        raise RadiativeIntegrityError(
            f"{context}: radiative integrity requires nonempty levels and transitions"
        )
    if not math.isfinite(data.mass_amu) or data.mass_amu <= 0.0:
        raise RadiativeIntegrityError(
            f"{context}: molecular mass is not finite and positive"
        )
    level_indices = [level.index for level in data.levels]
    if level_indices != list(range(1, len(data.levels) + 1)):
        raise RadiativeIntegrityError(
            f"{context}: level indices are not consecutive from one"
        )
    energies = [level.energy_cm1 for level in data.levels]
    if any(not math.isfinite(value) for value in energies) or any(
        right <= left for left, right in zip(energies, energies[1:])
    ):
        raise RadiativeIntegrityError(
            f"{context}: level energies are non-finite or not strictly increasing"
        )
    transition_indices = [row.index for row in data.radiative]
    if transition_indices != list(range(1, len(data.radiative) + 1)):
        raise RadiativeIntegrityError(
            f"{context}: radiative indices are not consecutive from one"
        )

    molecule_match = re.match(r"^([AE])-", data.molecule.strip())
    if molecule_match is None:
        raise RadiativeIntegrityError(
            f"{context}: molecule name does not declare A or E species"
        )
    declared_species = molecule_match.group(1)

    labels: dict[int, _QuantumIdentity] = {}
    seen_labels: set[str] = set()
    for level in data.levels:
        if level.index in labels:
            raise RadiativeIntegrityError(
                f"{context}: duplicate level index {level.index}"
            )
        label = _parse_label(level.quantum_label)
        if label.raw in seen_labels:
            raise RadiativeIntegrityError(
                f"{context}: duplicate level identity {label.raw}"
            )
        if label.species != declared_species:
            raise RadiativeIntegrityError(
                f"{context}: {declared_species} product contains {label.species} level {label.raw}"
            )
        seen_labels.add(label.raw)
        if not math.isclose(
            level.weight, 2.0 * label.j + 1.0, rel_tol=0.0, abs_tol=1.0e-12
        ):
            raise RadiativeIntegrityError(
                f"{context}: level {label.raw} has weight {level.weight}, "
                f"expected {2 * label.j + 1}"
            )
        labels[level.index] = label

    product = _product_ledger(data, labels, context=context)
    source = (
        _source_ledger(source_rows, context=context)
        if source_rows is not None
        else {}
    )
    if source_rows is not None and set(source) != set(product):
        raise RadiativeIntegrityError(
            f"{context}: source/product radiative graph mismatch "
            f"(missing={len(set(source) - set(product))}, "
            f"extra={len(set(product) - set(source))})"
        )

    inherited: dict[tuple[str, str], RadiativeTransition] = {}
    parent_level_labels: set[str] = set()
    if parent is not None:
        parent_labels = {
            level.index: _parse_label(level.quantum_label)
            for level in parent.levels
        }
        parent_level_labels = {label.raw for label in parent_labels.values()}
        missing_levels = parent_level_labels - seen_labels
        if missing_levels:
            raise RadiativeIntegrityError(
                f"{context}: child omits {len(missing_levels)} parent levels"
            )
        inherited = _product_ledger(parent, parent_labels, context=f"{context}:parent")
        restricted_child = {
            identity: row
            for identity, row in product.items()
            if identity[0] in parent_level_labels
            and identity[1] in parent_level_labels
        }
        if set(restricted_child) != set(inherited):
            raise RadiativeIntegrityError(
                f"{context}: parent/child radiative graph mismatch "
                f"(missing={len(set(inherited) - set(restricted_child))}, "
                f"extra={len(set(restricted_child) - set(inherited))})"
            )
        changed = [
            identity
            for identity, parent_row in inherited.items()
            if restricted_child[identity].a_s != parent_row.a_s
            or restricted_child[identity].frequency_ghz
            != parent_row.frequency_ghz
        ]
        if changed:
            raise RadiativeIntegrityError(
                f"{context}: {len(changed)} inherited A/frequency rows changed; "
                f"first={changed[0]}"
            )

    source_frequency_rows = 0
    maximum_frequency_offset = 0.0
    minimum_smu2 = math.inf
    maximum_smu2 = 0.0
    delta_j_counts = {0: 0, 1: 0}
    level_by_index = {level.index: level for level in data.levels}
    for identity, row in product.items():
        upper_label = labels[row.upper]
        lower_label = labels[row.lower]
        upper_level = level_by_index[row.upper]
        lower_level = level_by_index[row.lower]
        if upper_label.species != lower_label.species:
            raise RadiativeIntegrityError(
                f"{context}: radiative row {row.index} mixes A and E species"
            )
        delta_j = abs(upper_label.j - lower_label.j)
        if delta_j > 1 or (upper_label.j == 0 and lower_label.j == 0):
            raise RadiativeIntegrityError(
                f"{context}: forbidden rank-1 radiative row {row.index}: "
                f"{identity[0]} -> {identity[1]}"
            )
        delta_j_counts[delta_j] += 1
        if (
            not math.isfinite(row.a_s)
            or not math.isfinite(row.frequency_ghz)
            or row.a_s <= 0.0
            or row.frequency_ghz <= 0.0
        ):
            raise RadiativeIntegrityError(
                f"{context}: radiative row {row.index} has non-positive/non-finite A or frequency"
            )
        if upper_level.energy_cm1 <= lower_level.energy_cm1:
            raise RadiativeIntegrityError(
                f"{context}: radiative row {row.index} is not downward in the adopted levels"
            )
        if row.upper_energy_k is not None and not math.isclose(
            row.upper_energy_k,
            upper_level.energy_cm1 * LAMDA_CM1_TO_K,
            rel_tol=0.0,
            abs_tol=1.0e-6,
        ):
            raise RadiativeIntegrityError(
                f"{context}: radiative row {row.index} has inconsistent upper-state energy"
            )

        smu2 = _inferred_smu2(upper_level.weight, row)
        if not math.isfinite(smu2) or smu2 < minimum_smu2_debye2 * (1.0 - 1.0e-9):
            raise RadiativeIntegrityError(
                f"{context}: radiative row {row.index} has implausible/ineligible "
                f"S mu^2={smu2:.12g} D^2"
            )
        minimum_smu2 = min(minimum_smu2, smu2)
        maximum_smu2 = max(maximum_smu2, smu2)

        expected_frequency = (
            upper_level.energy_cm1 - lower_level.energy_cm1
        ) * CM1_TO_GHZ
        offset = abs(row.frequency_ghz - expected_frequency)
        maximum_frequency_offset = max(maximum_frequency_offset, offset)
        explicitly_bound = identity in source or identity in inherited
        if offset > CALCULATED_FREQUENCY_TOLERANCE_GHZ:
            if not explicitly_bound:
                raise RadiativeIntegrityError(
                    f"{context}: radiative row {row.index} frequency differs from "
                    f"the level gap by {offset:.12g} GHz without source/inheritance provenance"
                )
            if offset > MAXIMUM_SOURCE_FREQUENCY_OFFSET_GHZ:
                raise RadiativeIntegrityError(
                    f"{context}: radiative row {row.index} has suspicious source-bound "
                    f"frequency offset {offset:.12g} GHz"
                )
            source_frequency_rows += 1

        source_row = source.get(identity)
        if source_row is not None:
            try:
                expected_a = float(source_row["a_s"])
                expected_nu = float(source_row["frequency_ghz"])
                expected_smu2 = float(source_row["smu2_debye2"])
            except (KeyError, TypeError, ValueError) as exc:
                raise RadiativeIntegrityError(
                    f"{context}: radiative source row {identity} has malformed physics fields"
                ) from exc
            if not math.isclose(
                row.a_s, expected_a, rel_tol=5.0e-11, abs_tol=0.0
            ) or not math.isclose(
                row.frequency_ghz,
                expected_nu,
                rel_tol=5.0e-13,
                abs_tol=1.0e-12,
            ):
                raise RadiativeIntegrityError(
                    f"{context}: radiative row {row.index} differs from its source A/frequency"
                )
            if not math.isclose(
                smu2,
                expected_smu2,
                rel_tol=SMU2_RELATIVE_TOLERANCE,
                abs_tol=0.0,
            ):
                raise RadiativeIntegrityError(
                    f"{context}: radiative row {row.index} violates A-nu^3-S mu^2 consistency"
                )

    components = _radiative_components(data)
    if components != 1:
        raise RadiativeIntegrityError(
            f"{context}: radiative graph has {components} disconnected components"
        )

    return {
        "status": "pass",
        "levels": len(data.levels),
        "radiative_rows": len(data.radiative),
        "rank_rule": "abs(Delta J) <= 1; J=0 to J=0 excluded",
        "forbidden_rows": 0,
        "duplicate_or_reversed_endpoints": 0,
        "radiative_components": components,
        "delta_j_rows": {str(key): value for key, value in delta_j_counts.items()},
        "source_bound_rows": len(source),
        "source_frequency_rows": source_frequency_rows,
        "inherited_rows": len(inherited),
        "minimum_inferred_smu2_debye2": minimum_smu2,
        "maximum_inferred_smu2_debye2": maximum_smu2,
        "maximum_frequency_energy_offset_ghz": maximum_frequency_offset,
    }
