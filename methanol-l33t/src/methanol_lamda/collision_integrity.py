"""Fail-closed collision provenance and level-assignment checks.

These checks deliberately operate on quantum-label identities, never on the
numerical level indices written to a LAMDA file.  Level indices are local to a
particular product and can change whenever a manifold is added or energies are
reordered.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Protocol

from .basecol import (
    read_basecol_levels,
    read_basecol_rates,
    validate_level_rate_pair,
)
from .integrity import sha256_path
from .models import LamdaData
from .normalize import BASECOL_INPUTS


class CollisionIntegrityError(ValueError):
    """A collision row lost its source or product level identity."""


class SelectedCollisionLike(Protocol):
    collider_id: int
    upper: int
    lower: int
    upper_key: object
    lower_key: object
    rates: tuple[float, ...]


_MASTER_FIELDS = (
    "species",
    "collider_id",
    "upper_key",
    "lower_key",
    "temperatures_k",
    "rates_cm3_s",
    "source",
    "status",
    "raw_sha256",
    "source_upper_id",
    "source_lower_id",
    "source_upper_energy_cm1",
    "source_lower_energy_cm1",
    "mapping_confidence",
)


def _encoded(key: object) -> str:
    encode = getattr(key, "encode", None)
    if not callable(encode):
        raise CollisionIntegrityError(f"collision key {key!r} is not encodable")
    return str(encode())


def audit_basecol_collision_master(root: str | Path) -> dict[str, int]:
    """Reconstruct every normalized BASECOL row from the raw level/rate pair.

    The comparison covers source level IDs, their exact quantum identities,
    source energies, raw-file hashes, usable temperature columns, and every
    serialized rate coefficient.  It therefore detects a wrong level-table
    pairing as well as a later CSV mutation or endpoint reassignment.
    """

    root = Path(root)
    master = root / "data/normalized/master"
    with (master / "levels.csv").open(encoding="utf-8", newline="") as stream:
        canonical_keys = {row["level_key"] for row in csv.DictReader(stream)}

    families = {str(record["family"]) for record in BASECOL_INPUTS}
    actual: dict[tuple[int, int, int], dict[str, str]] = {}
    with (master / "collisions.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("source") not in families:
                continue
            try:
                collision_id = int(row["source_upper_id"].split(":")[1][3:])
                upper_id = int(row["source_upper_id"].rsplit(":", 1)[1])
                lower_id = int(row["source_lower_id"].rsplit(":", 1)[1])
            except (KeyError, ValueError, IndexError) as exc:
                raise CollisionIntegrityError(
                    "normalized BASECOL row has malformed source level IDs"
                ) from exc
            identity = (collision_id, upper_id, lower_id)
            if identity in actual:
                raise CollisionIntegrityError(
                    f"duplicate normalized BASECOL source identity {identity}"
                )
            actual[identity] = row

    expected: dict[tuple[int, int, int], dict[str, str]] = {}
    raw_level_rows = 0
    mapped_level_rows = 0
    raw_collision_rows = 0
    excluded_collision_rows = 0
    for record in BASECOL_INPUTS:
        species = str(record["species"])
        family = str(record["family"])
        collision_id = int(record["collision_id"])
        collider_id = int(record["collider_id"])
        raw_root = root / "data/raw/official" / str(record["directory"])
        level_path = raw_root / str(record["levels"])
        rate_path = raw_root / str(record["rates"])
        levels = read_basecol_levels(level_path, species)
        rates = read_basecol_rates(rate_path)
        validate_level_rate_pair(levels, rates)
        level_by_index = {level.source_index: level for level in levels}
        in_scope = {
            level.source_index
            for level in levels
            if level.key.encode() in canonical_keys
        }
        maximum_usable = record.get("maximum_usable_temperature_k")
        usable_positions = tuple(
            index
            for index, temperature in enumerate(rates.temperatures_k)
            if maximum_usable is None or temperature <= float(maximum_usable)
        )
        status = (
            "calculated"
            if maximum_usable is None
            else "calculated_source_grid_quarantined_180_200K"
        )
        rate_sha256 = sha256_path(rate_path)
        raw_level_rows += len(levels)
        mapped_level_rows += len(in_scope)
        raw_collision_rows += len(rates.transitions)
        for transition in rates.transitions:
            if (
                transition.initial_level not in in_scope
                or transition.final_level not in in_scope
            ):
                excluded_collision_rows += 1
                continue
            upper = level_by_index[transition.initial_level]
            lower = level_by_index[transition.final_level]
            identity = (collision_id, upper.source_index, lower.source_index)
            if identity in expected:
                raise CollisionIntegrityError(
                    f"duplicate raw BASECOL source identity {identity}"
                )
            expected[identity] = {
                "species": species,
                "collider_id": str(collider_id),
                "upper_key": upper.key.encode(),
                "lower_key": lower.key.encode(),
                "temperatures_k": ";".join(
                    f"{rates.temperatures_k[index]:.9g}"
                    for index in usable_positions
                ),
                "rates_cm3_s": ";".join(
                    f"{transition.rates_cm3_s[index]:.12g}"
                    for index in usable_positions
                ),
                "source": family,
                "status": status,
                "raw_sha256": rate_sha256,
                "source_upper_id": (
                    f"BASECOL:IDC{collision_id}:level:{upper.source_index}"
                ),
                "source_lower_id": (
                    f"BASECOL:IDC{collision_id}:level:{lower.source_index}"
                ),
                "source_upper_energy_cm1": f"{upper.source_energy_cm1:.12g}",
                "source_lower_energy_cm1": f"{lower.source_energy_cm1:.12g}",
                "mapping_confidence": "exact_quantum_identity",
            }

    missing = set(expected) - set(actual)
    extra = set(actual) - set(expected)
    mismatches: list[tuple[tuple[int, int, int], str]] = []
    for identity in set(expected) & set(actual):
        for field in _MASTER_FIELDS:
            if expected[identity][field] != actual[identity].get(field, ""):
                mismatches.append((identity, field))
                break
    if missing or extra or mismatches:
        raise CollisionIntegrityError(
            "normalized BASECOL collision master does not reconstruct from raw "
            f"sources (missing={len(missing)}, extra={len(extra)}, "
            f"mismatched={len(mismatches)}; samples="
            f"{sorted(missing)[:1] + sorted(extra)[:1] + mismatches[:1]})"
        )
    return {
        "source_datasets": len(BASECOL_INPUTS),
        "raw_level_rows": raw_level_rows,
        "mapped_level_rows": mapped_level_rows,
        "raw_collision_rows": raw_collision_rows,
        "out_of_canonical_scope_collision_rows": excluded_collision_rows,
        "normalized_collision_rows": len(actual),
        "missing_rows": 0,
        "extra_rows": 0,
        "mismatched_rows": 0,
    }


def assert_collision_selection(
    data: LamdaData,
    selected_rows: Iterable[SelectedCollisionLike],
    temperatures_k: tuple[float, ...],
    *,
    context: str,
) -> dict[str, int]:
    """Bind each generated collision row to its selected quantum identities."""

    labels = {level.index: level.quantum_label for level in data.levels}
    expected: dict[tuple[int, str, str], tuple[float, ...]] = {}
    for row in selected_rows:
        try:
            actual_upper_label = labels[row.upper]
            actual_lower_label = labels[row.lower]
        except KeyError as exc:
            raise CollisionIntegrityError(
                f"{context}: selected collision has unresolved level index"
            ) from exc
        upper_label = _encoded(row.upper_key)
        lower_label = _encoded(row.lower_key)
        if (actual_upper_label, actual_lower_label) != (
            upper_label,
            lower_label,
        ):
            raise CollisionIntegrityError(
                f"{context}: collider {row.collider_id} indices "
                f"({row.upper}, {row.lower}) name "
                f"{actual_upper_label}->{actual_lower_label}, expected "
                f"{upper_label}->{lower_label}"
            )
        identity = (row.collider_id, upper_label, lower_label)
        reverse = (row.collider_id, lower_label, upper_label)
        if identity in expected or reverse in expected:
            raise CollisionIntegrityError(
                f"{context}: duplicate or reversed selected collision {identity}"
            )
        expected[identity] = tuple(row.rates)

    actual: dict[tuple[int, str, str], tuple[float, ...]] = {}
    for partner in data.collision_partners:
        if tuple(partner.temperatures_k) != tuple(temperatures_k):
            raise CollisionIntegrityError(
                f"{context}: collider {partner.collider_id} temperature grid drifted"
            )
        for row in partner.transitions:
            try:
                identity = (
                    partner.collider_id,
                    labels[row.upper],
                    labels[row.lower],
                )
            except KeyError as exc:
                raise CollisionIntegrityError(
                    f"{context}: serialized collision has unresolved level index"
                ) from exc
            reverse = (identity[0], identity[2], identity[1])
            if identity in actual or reverse in actual:
                raise CollisionIntegrityError(
                    f"{context}: duplicate or reversed product collision {identity}"
                )
            actual[identity] = tuple(row.rates_cm3_s)

    missing = set(expected) - set(actual)
    extra = set(actual) - set(expected)
    changed = {
        identity
        for identity in set(expected) & set(actual)
        if expected[identity] != actual[identity]
    }
    if missing or extra or changed:
        raise CollisionIntegrityError(
            f"{context}: collision selection/product mismatch "
            f"(missing={len(missing)}, extra={len(extra)}, changed={len(changed)})"
        )
    return {"collision_rows": len(actual), "colliders": len(data.collision_partners)}


def _collision_ledger(
    data: LamdaData,
) -> dict[
    tuple[int, str, str],
    tuple[str, str, tuple[float, ...], tuple[float, ...]],
]:
    labels = {level.index: level.quantum_label for level in data.levels}
    ledger: dict[
        tuple[int, str, str],
        tuple[str, str, tuple[float, ...], tuple[float, ...]],
    ] = {}
    for partner in data.collision_partners:
        temperatures = tuple(partner.temperatures_k)
        for row in partner.transitions:
            try:
                identity = (
                    partner.collider_id,
                    labels[row.upper],
                    labels[row.lower],
                )
            except KeyError as exc:
                raise CollisionIntegrityError(
                    "collision product contains an unresolved level index"
                ) from exc
            unordered = (
                identity[0],
                min(identity[1], identity[2]),
                max(identity[1], identity[2]),
            )
            if unordered in ledger:
                raise CollisionIntegrityError(
                    f"duplicate or reversed collision identity {identity}"
                )
            ledger[unordered] = (
                identity[1],
                identity[2],
                temperatures,
                tuple(row.rates_cm3_s),
            )
    return ledger


def assert_parent_collision_inheritance(
    parent: LamdaData,
    child: LamdaData,
    *,
    context: str,
) -> dict[str, int | str | list[str]]:
    """Require exact preservation of the parent's labelled collision graph."""

    parent_labels = {level.quantum_label for level in parent.levels}
    child_labels = {level.quantum_label for level in child.levels}
    missing_levels = parent_labels - child_labels
    if missing_levels:
        raise CollisionIntegrityError(
            f"{context}: child omits {len(missing_levels)} parent collision levels"
        )
    parent_ledger = _collision_ledger(parent)
    child_ledger = _collision_ledger(child)
    restricted_child = {
        identity: value
        for identity, value in child_ledger.items()
        if identity[1] in parent_labels and identity[2] in parent_labels
    }
    missing = set(parent_ledger) - set(restricted_child)
    extra = set(restricted_child) - set(parent_ledger)
    changed: set[tuple[int, str, str]] = set()
    reoriented = 0
    parent_level_by_label = {
        level.quantum_label: level for level in parent.levels
    }
    for identity in set(parent_ledger) & set(restricted_child):
        parent_upper, parent_lower, parent_temperatures, parent_rates = (
            parent_ledger[identity]
        )
        child_upper, child_lower, child_temperatures, child_rates = (
            restricted_child[identity]
        )
        common_count = len(parent_temperatures)
        if child_temperatures[:common_count] != parent_temperatures:
            changed.add(identity)
            continue
        if (child_upper, child_lower) == (parent_upper, parent_lower):
            expected_rates = parent_rates
        elif (child_upper, child_lower) == (parent_lower, parent_upper):
            reoriented += 1
            source_upper = parent_level_by_label[parent_upper]
            source_lower = parent_level_by_label[parent_lower]
            gap = source_upper.energy_cm1 - source_lower.energy_cm1
            ratio = source_upper.weight / source_lower.weight
            expected_rates = tuple(
                rate * ratio
                if gap == 0.0
                else rate
                * ratio
                * math.exp(-gap * 1.4387768775039338 / temperature)
                for rate, temperature in zip(parent_rates, parent_temperatures)
            )
        else:  # pragma: no cover - unordered identity makes this unreachable
            changed.add(identity)
            continue
        if tuple(f"{value:.8e}" for value in expected_rates) != tuple(
            f"{value:.8e}" for value in child_rates[:common_count]
        ):
            changed.add(identity)
    if missing or changed:
        raise CollisionIntegrityError(
            f"{context}: parent collision inheritance failed "
            f"(missing={len(missing)}, extra={len(extra)}, changed={len(changed)})"
        )
    return {
        "upstream_level_count": len(parent_labels),
        "inherited_collision_rows": len(parent_ledger),
        "endpoint_policy": "exact ordered quantum-label identity",
        "inherited_fields": [
            "common_temperatures_k",
            "rates_cm3_s_on_common_grid",
        ],
        "reoriented_by_detailed_balance": reoriented,
        "generated_target_only_rows": len(child_ledger) - len(parent_ledger),
        "generated_rows_between_parent_levels": len(extra),
    }
