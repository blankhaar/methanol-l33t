"""Independent, fail-closed audit of the final methanol LAMDA products.

This module intentionally uses only the Python standard library and does not
import the builders, their data models, or either existing LAMDA parser.  It
therefore provides an independent check of the serialized products.  LAMDA
files are consumed line by line.  Collision rows are retained only as compact
quantum-label ledgers for the current parent/child tier edge, so inheritance
can be checked without trusting product-local level indices.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterator, Mapping, Sequence


_A_LABEL = re.compile(
    r"^A:(?P<manifold>vt[0-4]|vCO1:vt0|vrock1:vt0):"
    r"J(?P<j>0|[1-9][0-9]*):K(?P<k>0|[1-9][0-9]*):(?P<parity>[+-])$"
)
_E_LABEL = re.compile(
    r"^E:(?P<manifold>vt[0-4]|vCO1:vt0|vrock1:vt0):"
    r"J(?P<j>0|[1-9][0-9]*):K(?P<sign>[+-])(?P<k>0|[1-9][0-9]*)$"
)


class AuditParseError(ValueError):
    """A product could not be parsed without guessing."""


@dataclass(frozen=True, slots=True)
class QuantumLabel:
    raw: str
    species: str
    manifold: str
    j: int
    k: int
    parity: str | None


@dataclass(frozen=True, slots=True)
class LevelRow:
    index: int
    energy_cm1: float
    weight: float
    label: QuantumLabel
    energy_token: str
    weight_token: str


@dataclass(frozen=True, slots=True)
class RadiativeRow:
    index: int
    upper_index: int
    lower_index: int
    upper_label: str
    lower_label: str
    a_s: float
    frequency_ghz: float
    a_token: str
    frequency_token: str

    @property
    def endpoints(self) -> tuple[str, str]:
        return self.upper_label, self.lower_label


@dataclass(frozen=True, slots=True)
class CollisionRow:
    upper_label: str
    lower_label: str
    rates_cm3_s: tuple[float, ...]
    rate_tokens: tuple[str, ...]


@dataclass(slots=True)
class ParsedProduct:
    molecule: str
    mass_amu: float
    levels_by_index: dict[int, LevelRow]
    levels_by_label: dict[str, LevelRow]
    radiative_by_index: dict[int, RadiativeRow]
    radiative_by_endpoints: dict[tuple[str, str], RadiativeRow]
    collision_rows_by_collider: dict[int, int]
    collision_temperatures_by_collider: dict[int, list[float]]
    collision_by_identity: dict[tuple[int, str, str], CollisionRow]


class _DataCursor:
    """Streaming cursor over non-comment, nonblank LAMDA records."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self._stream = path.open("r", encoding="utf-8", errors="strict")
        self.line_number = 0

    def close(self) -> None:
        self._stream.close()

    def __enter__(self) -> "_DataCursor":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def next(self) -> tuple[int, str]:
        for raw in self._stream:
            self.line_number += 1
            value = raw.strip()
            if value and not value.startswith(("!", "#")):
                return self.line_number, value
        raise AuditParseError(
            f"{self.path}: unexpected EOF after physical line {self.line_number}"
        )

    def optional(self) -> tuple[int, str] | None:
        try:
            return self.next()
        except AuditParseError:
            return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _number(token: str, *, context: str) -> float:
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise AuditParseError(f"{context}: invalid floating-point token {token!r}") from exc
    if not math.isfinite(value):
        raise AuditParseError(f"{context}: non-finite value {token!r}")
    return value


def _integer(token: str, *, context: str) -> int:
    try:
        return int(token)
    except ValueError as exc:
        raise AuditParseError(f"{context}: invalid integer token {token!r}") from exc


def _parse_label(raw: str) -> QuantumLabel:
    match = _A_LABEL.fullmatch(raw)
    if match:
        j = int(match.group("j"))
        k = int(match.group("k"))
        parity = match.group("parity")
        if k == 0 and parity != "+":
            raise AuditParseError(f"noncanonical A K=0 parity in {raw!r}")
        return QuantumLabel(raw, "A", match.group("manifold"), j, k, parity)
    match = _E_LABEL.fullmatch(raw)
    if match:
        j = int(match.group("j"))
        magnitude = int(match.group("k"))
        sign = match.group("sign")
        if magnitude == 0 and sign != "+":
            raise AuditParseError(f"noncanonical E K=-0 label {raw!r}")
        k = magnitude if sign == "+" else -magnitude
        return QuantumLabel(raw, "E", match.group("manifold"), j, k, None)
    raise AuditParseError(f"noncanonical methanol quantum label {raw!r}")


def _failure(
    failures: list[dict[str, object]],
    gate: str,
    product: str,
    message: str,
    *,
    line: int | None = None,
) -> None:
    record: dict[str, object] = {"gate": gate, "product": product, "message": message}
    if line is not None:
        record["line"] = line
    failures.append(record)


def _parse_product(
    path: Path,
    *,
    product_id: str,
    species: str,
    isotopologue: str,
    expected_mass: float,
    tier: Mapping[str, Any],
    failures: list[dict[str, object]],
) -> ParsedProduct:
    allowed_manifolds = set(tier["allowed_manifolds"])
    expected_colliders = {int(key): int(value) for key, value in tier["collision_rows_by_collider"].items()}
    expected_temperatures = [float(value) for value in tier["collision_temperatures_k"]]

    with _DataCursor(path) as cursor:
        _, molecule = cursor.next()
        mass_line, mass_record = cursor.next()
        mass_token = mass_record.split()[0]
        mass_amu = _number(mass_token, context=f"{path}:{mass_line} mass")
        if abs(mass_amu - expected_mass) > 5e-7:
            _failure(
                failures,
                "isotopologue_mass",
                product_id,
                f"mass {mass_amu:.12g} differs from configured {expected_mass:.12g}",
                line=mass_line,
            )
        if not molecule.startswith(f"{species}-"):
            _failure(
                failures,
                "molecule_identity",
                product_id,
                f"molecule name {molecule!r} does not start with {species!r}",
            )
        if isotopologue == "13CH3OH" and "13CH3OH" not in molecule:
            _failure(failures, "molecule_identity", product_id, "13CH3OH marker is absent")
        if isotopologue == "12CH3OH" and "13CH3OH" in molecule:
            _failure(failures, "molecule_identity", product_id, "12C product names 13CH3OH")

        count_line, count_record = cursor.next()
        level_count = _integer(count_record.split()[0], context=f"{path}:{count_line} level count")
        if level_count < 1 or level_count > 100_000:
            raise AuditParseError(f"{path}:{count_line}: implausible level count {level_count}")
        if level_count != int(tier["levels"]):
            _failure(
                failures,
                "expected_counts",
                product_id,
                f"level count {level_count} != expected {tier['levels']}",
                line=count_line,
            )

        levels_by_index: dict[int, LevelRow] = {}
        levels_by_label: dict[str, LevelRow] = {}
        previous_energy = -math.inf
        manifolds: dict[str, int] = {}
        for expected_index in range(1, level_count + 1):
            line_number, record = cursor.next()
            fields = record.split(maxsplit=3)
            if len(fields) != 4:
                raise AuditParseError(f"{path}:{line_number}: malformed level row")
            index = _integer(fields[0], context=f"{path}:{line_number} level index")
            energy = _number(fields[1], context=f"{path}:{line_number} level energy")
            weight = _number(fields[2], context=f"{path}:{line_number} level weight")
            label = _parse_label(fields[3])
            if index != expected_index:
                _failure(
                    failures,
                    "level_indices",
                    product_id,
                    f"level index {index} occurs where {expected_index} is required",
                    line=line_number,
                )
            if index in levels_by_index:
                _failure(failures, "level_indices", product_id, f"duplicate level index {index}", line=line_number)
            if label.raw in levels_by_label:
                _failure(failures, "level_labels", product_id, f"duplicate level label {label.raw}", line=line_number)
            if energy <= previous_energy:
                _failure(
                    failures,
                    "level_order",
                    product_id,
                    f"energy {energy:.17g} does not strictly exceed {previous_energy:.17g}",
                    line=line_number,
                )
            previous_energy = energy
            if label.species != species:
                _failure(failures, "quantum_labels", product_id, f"wrong species in {label.raw}", line=line_number)
            if label.manifold not in allowed_manifolds:
                _failure(failures, "quantum_labels", product_id, f"disallowed manifold in {label.raw}", line=line_number)
            if not (0 <= label.j <= 22):
                _failure(failures, "quantum_labels", product_id, f"J outside 0..22 in {label.raw}", line=line_number)
            if abs(label.k) > min(14, label.j):
                _failure(failures, "quantum_labels", product_id, f"|K| invalid in {label.raw}", line=line_number)
            expected_weight = float(2 * label.j + 1)
            if weight != expected_weight:
                _failure(
                    failures,
                    "statistical_weights",
                    product_id,
                    f"weight {weight:.12g} != 2J+1={expected_weight:.12g} for {label.raw}",
                    line=line_number,
                )
            row = LevelRow(index, energy, weight, label, fields[1], fields[2])
            levels_by_index[index] = row
            levels_by_label[label.raw] = row
            manifolds[label.manifold] = manifolds.get(label.manifold, 0) + 1

        expected_per_manifold = tier.get("levels_per_manifold")
        if expected_per_manifold is not None:
            expected_manifold_counts = {
                manifold: int(expected_per_manifold) for manifold in allowed_manifolds
            }
            if manifolds != expected_manifold_counts:
                _failure(
                    failures,
                    "quantum_labels",
                    product_id,
                    f"manifold level census {dict(sorted(manifolds.items()))!r} != expected {dict(sorted(expected_manifold_counts.items()))!r}",
                )

        rad_count_line, rad_count_record = cursor.next()
        radiative_count = _integer(
            rad_count_record.split()[0], context=f"{path}:{rad_count_line} radiative count"
        )
        if radiative_count < 0 or radiative_count > 2_000_000:
            raise AuditParseError(f"{path}:{rad_count_line}: implausible radiative count {radiative_count}")
        radiative_by_index: dict[int, RadiativeRow] = {}
        radiative_by_endpoints: dict[tuple[str, str], RadiativeRow] = {}
        seen_unordered_radiative: set[tuple[int, int]] = set()
        for expected_index in range(1, radiative_count + 1):
            line_number, record = cursor.next()
            fields = record.split()
            if len(fields) < 5:
                raise AuditParseError(f"{path}:{line_number}: malformed radiative row")
            index = _integer(fields[0], context=f"{path}:{line_number} radiative index")
            upper = _integer(fields[1], context=f"{path}:{line_number} radiative upper")
            lower = _integer(fields[2], context=f"{path}:{line_number} radiative lower")
            a_s = _number(fields[3], context=f"{path}:{line_number} Einstein A")
            frequency = _number(fields[4], context=f"{path}:{line_number} frequency")
            if index != expected_index:
                _failure(
                    failures,
                    "radiative_indices",
                    product_id,
                    f"radiative index {index} occurs where {expected_index} is required",
                    line=line_number,
                )
            if upper not in levels_by_index or lower not in levels_by_index:
                _failure(
                    failures,
                    "radiative_endpoints",
                    product_id,
                    f"unresolved endpoints ({upper}, {lower})",
                    line=line_number,
                )
                continue
            upper_level = levels_by_index[upper]
            lower_level = levels_by_index[lower]
            endpoints = (upper_level.label.raw, lower_level.label.raw)
            if endpoints in radiative_by_endpoints:
                _failure(failures, "radiative_endpoints", product_id, f"duplicate endpoints {endpoints}", line=line_number)
            unordered = (min(upper, lower), max(upper, lower))
            if unordered in seen_unordered_radiative:
                _failure(failures, "radiative_endpoints", product_id, f"reverse/duplicate endpoint pair {unordered}", line=line_number)
            seen_unordered_radiative.add(unordered)
            if a_s <= 0.0 or frequency <= 0.0:
                _failure(failures, "radiative_physics", product_id, "A and frequency must be positive", line=line_number)
            if upper_level.energy_cm1 <= lower_level.energy_cm1:
                _failure(failures, "radiative_downward", product_id, f"transition {index} is not downward", line=line_number)
            delta_j = abs(upper_level.label.j - lower_level.label.j)
            if delta_j > 1 or (upper_level.label.j == 0 and lower_level.label.j == 0):
                _failure(
                    failures,
                    "rank1_selection",
                    product_id,
                    f"transition {index} violates rank-1 J rule: {upper_level.label.j}->{lower_level.label.j}",
                    line=line_number,
                )
            row = RadiativeRow(
                index,
                upper,
                lower,
                upper_level.label.raw,
                lower_level.label.raw,
                a_s,
                frequency,
                fields[3],
                fields[4],
            )
            radiative_by_index[index] = row
            radiative_by_endpoints[endpoints] = row

        partner_count_line, partner_count_record = cursor.next()
        partner_count = _integer(
            partner_count_record.split()[0], context=f"{path}:{partner_count_line} partner count"
        )
        if partner_count != len(expected_colliders):
            _failure(
                failures,
                "collision_topology",
                product_id,
                f"partner count {partner_count} != expected {len(expected_colliders)}",
                line=partner_count_line,
            )
        collision_rows_by_collider: dict[int, int] = {}
        collision_temperatures_by_collider: dict[int, list[float]] = {}
        collision_by_identity: dict[tuple[int, str, str], CollisionRow] = {}
        for _ in range(partner_count):
            identity_line, identity_record = cursor.next()
            collider = _integer(identity_record.split()[0], context=f"{path}:{identity_line} collider")
            transition_count_line, transition_count_record = cursor.next()
            transition_count = _integer(
                transition_count_record.split()[0],
                context=f"{path}:{transition_count_line} collision count",
            )
            temperature_count_line, temperature_count_record = cursor.next()
            temperature_count = _integer(
                temperature_count_record.split()[0],
                context=f"{path}:{temperature_count_line} collision temperature count",
            )
            temperature_tokens: list[str] = []
            while len(temperature_tokens) < temperature_count:
                _, temperature_record = cursor.next()
                temperature_tokens.extend(temperature_record.split())
            if len(temperature_tokens) != temperature_count:
                raise AuditParseError(
                    f"{path}: collider {collider} supplies {len(temperature_tokens)} temperature tokens, expected {temperature_count}"
                )
            temperatures = [
                _number(token, context=f"{path} collider {collider} temperature")
                for token in temperature_tokens
            ]
            if collider in collision_rows_by_collider:
                _failure(failures, "collision_topology", product_id, f"duplicate collider {collider}", line=identity_line)
            collision_rows_by_collider[collider] = transition_count
            collision_temperatures_by_collider[collider] = temperatures
            expected_rows = expected_colliders.get(collider)
            if expected_rows is None:
                _failure(failures, "collision_topology", product_id, f"unexpected collider {collider}", line=identity_line)
            elif transition_count != expected_rows:
                _failure(
                    failures,
                    "collision_counts",
                    product_id,
                    f"collider {collider} has {transition_count} rows, expected {expected_rows}",
                    line=transition_count_line,
                )
            if temperatures != expected_temperatures:
                _failure(
                    failures,
                    "collision_temperatures",
                    product_id,
                    f"collider {collider} grid {temperatures!r} != expected {expected_temperatures!r}",
                    line=temperature_count_line,
                )
            seen_collision_pairs: set[tuple[int, int]] = set()
            for expected_index in range(1, transition_count + 1):
                line_number, record = cursor.next()
                fields = record.split()
                if len(fields) != 3 + temperature_count:
                    raise AuditParseError(
                        f"{path}:{line_number}: collision row has {len(fields) - 3} rates, expected {temperature_count}"
                    )
                index = _integer(fields[0], context=f"{path}:{line_number} collision index")
                upper = _integer(fields[1], context=f"{path}:{line_number} collision upper")
                lower = _integer(fields[2], context=f"{path}:{line_number} collision lower")
                rates = [
                    _number(token, context=f"{path}:{line_number} collision rate")
                    for token in fields[3:]
                ]
                if index != expected_index:
                    _failure(
                        failures,
                        "collision_indices",
                        product_id,
                        f"collider {collider} index {index} occurs where {expected_index} is required",
                        line=line_number,
                    )
                if upper not in levels_by_index or lower not in levels_by_index:
                    _failure(
                        failures,
                        "collision_endpoints",
                        product_id,
                        f"collider {collider} has unresolved endpoints ({upper}, {lower})",
                        line=line_number,
                    )
                    continue
                pair = (min(upper, lower), max(upper, lower))
                if pair in seen_collision_pairs:
                    _failure(
                        failures,
                        "collision_endpoints",
                        product_id,
                        f"collider {collider} repeats or reverses pair {pair}",
                        line=line_number,
                    )
                seen_collision_pairs.add(pair)
                upper_level = levels_by_index[upper]
                lower_level = levels_by_index[lower]
                unordered_identity = (
                    collider,
                    min(upper_level.label.raw, lower_level.label.raw),
                    max(upper_level.label.raw, lower_level.label.raw),
                )
                if unordered_identity in collision_by_identity:
                    _failure(
                        failures,
                        "collision_endpoints",
                        product_id,
                        "collision quantum-label identity is duplicated or reversed",
                        line=line_number,
                    )
                else:
                    collision_by_identity[unordered_identity] = CollisionRow(
                        upper_level.label.raw,
                        lower_level.label.raw,
                        tuple(rates),
                        tuple(fields[3:]),
                    )
                if upper_level.energy_cm1 <= lower_level.energy_cm1:
                    _failure(
                        failures,
                        "collision_downward",
                        product_id,
                        f"collider {collider} transition {index} is not downward",
                        line=line_number,
                    )
                if any(rate < 0.0 for rate in rates):
                    _failure(
                        failures,
                        "collision_rates",
                        product_id,
                        f"collider {collider} transition {index} has a negative rate",
                        line=line_number,
                    )
                upper_manifold = upper_level.label.manifold
                lower_manifold = lower_level.label.manifold
                if collider in {2, 3} and upper_manifold != lower_manifold:
                    _failure(
                        failures,
                        "collision_topology",
                        product_id,
                        f"H2 collider {collider} crosses manifolds {upper_manifold}->{lower_manifold}",
                        line=line_number,
                    )
                if collider == 6 and upper_manifold != lower_manifold:
                    if not (upper_manifold.startswith("vt") and lower_manifold.startswith("vt")):
                        _failure(
                            failures,
                            "collision_topology",
                            product_id,
                            f"He row invents a vibrational cross-manifold pair {upper_manifold}->{lower_manifold}",
                            line=line_number,
                        )

        missing_colliders = sorted(set(expected_colliders) - set(collision_rows_by_collider))
        if missing_colliders:
            _failure(
                failures,
                "collision_topology",
                product_id,
                f"missing colliders {missing_colliders}",
            )
        trailing = cursor.optional()
        if trailing is not None:
            _failure(
                failures,
                "trailing_records",
                product_id,
                f"unexpected trailing data record {trailing[1]!r}",
                line=trailing[0],
            )

    return ParsedProduct(
        molecule,
        mass_amu,
        levels_by_index,
        levels_by_label,
        radiative_by_index,
        radiative_by_endpoints,
        collision_rows_by_collider,
        collision_temperatures_by_collider,
        collision_by_identity,
    )


def _manifest_entry(
    root: Path,
    manifest_path: Path,
    product_path: Path,
    *,
    product_id: str,
    species: str,
    file_sha256: str,
    failures: list[dict[str, object]],
) -> tuple[dict[str, Any], str]:
    try:
        manifest_sha256 = _sha256(manifest_path)
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        _failure(failures, "manifest_binding", product_id, f"cannot read manifest {manifest_path}: {exc}")
        return {}, ""
    matches = [
        entry
        for entry in payload.get("files", [])
        if isinstance(entry, dict) and entry.get("file") == product_path.name
    ]
    if len(matches) != 1:
        _failure(
            failures,
            "manifest_binding",
            product_id,
            f"manifest has {len(matches)} entries for {product_path.name!r}; expected one",
        )
        return {}, manifest_sha256
    entry = matches[0]
    if entry.get("sha256") != file_sha256:
        _failure(
            failures,
            "manifest_binding",
            product_id,
            f"file SHA-256 {file_sha256} != manifest {entry.get('sha256')!r}",
        )
    if entry.get("species") != species:
        _failure(
            failures,
            "manifest_binding",
            product_id,
            f"manifest species {entry.get('species')!r} != {species!r}",
        )
    return entry, manifest_sha256


def _manifest_declared_input_hash_gate(
    root: Path,
    manifest_path: Path,
    *,
    product_id: str,
    failures: list[dict[str, object]],
) -> int:
    """Verify every path-addressable source/build input declared by a manifest."""

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    declared: set[tuple[str, str, str]] = set()

    def collect(node: object, context: str) -> None:
        if isinstance(node, dict):
            for field in ("source_sha256", "build_inputs_sha256"):
                values = node.get(field)
                if isinstance(values, dict):
                    for relative, expected in values.items():
                        if isinstance(relative, str) and isinstance(expected, str):
                            declared.add((relative, expected, f"{context}/{field}"))
            for key, value in node.items():
                if key.endswith("_path") and isinstance(value, str):
                    expected = node.get(f"{key[:-5]}_sha256")
                    if isinstance(expected, str):
                        declared.add((value, expected, f"{context}/{key}"))
                collect(value, f"{context}/{key}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                collect(value, f"{context}/{index}")

    collect(payload, "manifest")
    checked = 0
    for relative, expected, context in sorted(declared):
        candidate = (root / relative).resolve()
        try:
            candidate.relative_to(root)
        except ValueError:
            _failure(
                failures,
                "manifest_input_binding",
                product_id,
                f"{context} escapes repository root: {relative!r}",
            )
            continue
        if not candidate.is_file():
            _failure(
                failures,
                "manifest_input_binding",
                product_id,
                f"{context} is missing: {relative!r}",
            )
            continue
        actual = _sha256(candidate)
        checked += 1
        if actual != expected:
            _failure(
                failures,
                "manifest_input_binding",
                product_id,
                f"{context} SHA-256 {actual} != {expected} for {relative!r}",
            )
    return checked


def _read_overlay_sidecar(
    path: Path,
    *,
    expected_sha256: object,
    expected_count: int,
    product_id: str,
    parsed: ParsedProduct,
    failures: list[dict[str, object]],
) -> tuple[set[tuple[str, str]], dict[str, object]]:
    sha256 = _sha256(path)
    if sha256 != expected_sha256:
        _failure(
            failures,
            "sidecar_binding",
            product_id,
            f"radiative sidecar SHA-256 {sha256} != manifest {expected_sha256!r}",
        )
    overlay_endpoints: set[tuple[str, str]] = set()
    seen_indices: set[int] = set()
    row_count = 0
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"transition_index", "upper_key", "lower_key", "a_s", "frequency_ghz", "status"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            _failure(failures, "sidecar_binding", product_id, "radiative sidecar lacks required columns")
            return set(), {"path": str(path), "sha256": sha256, "rows": 0, "experimental_rows": 0}
        for row_count, row in enumerate(reader, start=1):
            try:
                index = int(row["transition_index"])
                a_s = float(row["a_s"])
                frequency = float(row["frequency_ghz"])
            except (TypeError, ValueError) as exc:
                _failure(failures, "sidecar_binding", product_id, f"sidecar row {row_count} is malformed: {exc}")
                continue
            if index in seen_indices:
                _failure(failures, "sidecar_binding", product_id, f"sidecar repeats transition index {index}")
                continue
            seen_indices.add(index)
            transition = parsed.radiative_by_index.get(index)
            if transition is None:
                _failure(failures, "sidecar_binding", product_id, f"sidecar index {index} is absent from LAMDA")
                continue
            endpoints = (row["upper_key"], row["lower_key"])
            if endpoints != transition.endpoints:
                _failure(
                    failures,
                    "sidecar_binding",
                    product_id,
                    f"sidecar index {index} endpoints {endpoints} != LAMDA {transition.endpoints}",
                )
            if not math.isclose(a_s, transition.a_s, rel_tol=1e-9, abs_tol=0.0):
                _failure(failures, "sidecar_binding", product_id, f"sidecar index {index} A differs from LAMDA")
            if not math.isclose(frequency, transition.frequency_ghz, rel_tol=5e-13, abs_tol=1e-12):
                _failure(failures, "sidecar_binding", product_id, f"sidecar index {index} frequency differs from LAMDA")
            status = row["status"].lower()
            if "observed_frequency" in status or "experimental" in status:
                overlay_endpoints.add(transition.endpoints)
    if row_count != len(parsed.radiative_by_index):
        _failure(
            failures,
            "sidecar_binding",
            product_id,
            f"sidecar has {row_count} rows, LAMDA has {len(parsed.radiative_by_index)}",
        )
    if seen_indices != set(parsed.radiative_by_index):
        _failure(failures, "sidecar_binding", product_id, "sidecar transition-index census is incomplete")
    if len(overlay_endpoints) != expected_count:
        _failure(
            failures,
            "experimental_overlays",
            product_id,
            f"experimental overlay count {len(overlay_endpoints)} != expected {expected_count}",
        )
    return overlay_endpoints, {
        "path": str(path),
        "sha256": sha256,
        "rows": row_count,
        "experimental_rows": len(overlay_endpoints),
    }


def _frequency_gap_gate(
    parsed: ParsedProduct,
    *,
    exceptions: set[tuple[str, str]],
    tolerance_ghz: float,
    cm1_to_ghz: float,
    product_id: str,
    failures: list[dict[str, object]],
) -> dict[str, object]:
    checked = 0
    maximum_checked = 0.0
    maximum_exception = 0.0
    for endpoints, transition in parsed.radiative_by_endpoints.items():
        upper = parsed.levels_by_index[transition.upper_index]
        lower = parsed.levels_by_index[transition.lower_index]
        residual = abs(
            transition.frequency_ghz
            - (upper.energy_cm1 - lower.energy_cm1) * cm1_to_ghz
        )
        if endpoints in exceptions:
            maximum_exception = max(maximum_exception, residual)
            continue
        checked += 1
        maximum_checked = max(maximum_checked, residual)
        if residual > tolerance_ghz:
            _failure(
                failures,
                "frequency_energy_consistency",
                product_id,
                f"transition {transition.index} residual {residual:.12g} GHz exceeds {tolerance_ghz:.12g}",
            )
    return {
        "checked_rows": checked,
        "exact_exceptions": len(exceptions),
        "maximum_checked_residual_ghz": maximum_checked,
        "maximum_exception_residual_ghz": maximum_exception,
        "tolerance_ghz": tolerance_ghz,
    }


def _edge_gate(
    root: Path,
    parent_spec: Mapping[str, Any],
    child_spec: Mapping[str, Any],
    parent: ParsedProduct,
    child: ParsedProduct,
    *,
    child_manifest_entry: Mapping[str, Any],
    parent_file_sha256: str,
    parent_manifest_sha256: str,
    energy_tolerance_cm1: float,
    overlay_endpoints: set[tuple[str, str]],
    failures: list[dict[str, object]],
) -> dict[str, object]:
    child_id = str(child_spec["id"])
    missing_levels = 0
    energy_mismatches = 0
    weight_mismatches = 0
    max_energy_delta = 0.0
    for label, parent_level in parent.levels_by_label.items():
        child_level = child.levels_by_label.get(label)
        if child_level is None:
            missing_levels += 1
            continue
        delta = abs(parent_level.energy_cm1 - child_level.energy_cm1)
        max_energy_delta = max(max_energy_delta, delta)
        if delta > energy_tolerance_cm1:
            energy_mismatches += 1
        if parent_level.weight_token != child_level.weight_token:
            weight_mismatches += 1
    if missing_levels or energy_mismatches or weight_mismatches:
        _failure(
            failures,
            "tier_level_nesting",
            child_id,
            f"parent nesting failures: missing={missing_levels}, energy={energy_mismatches}, weight={weight_mismatches}",
        )

    missing_radiative = 0
    token_mismatches = 0
    for endpoints, parent_row in parent.radiative_by_endpoints.items():
        child_row = child.radiative_by_endpoints.get(endpoints)
        if child_row is None:
            missing_radiative += 1
            continue
        if (parent_row.a_token, parent_row.frequency_token) != (
            child_row.a_token,
            child_row.frequency_token,
        ):
            token_mismatches += 1
    if missing_radiative or token_mismatches:
        _failure(
            failures,
            "tier_radiative_nesting",
            child_id,
            f"parent radiative failures: missing={missing_radiative}, exact-token mismatches={token_mismatches}",
        )
    missing_overlays = len(overlay_endpoints - set(child.radiative_by_endpoints))
    if missing_overlays:
        _failure(
            failures,
            "experimental_overlays",
            child_id,
            f"{missing_overlays} base experimental overlay endpoints are absent",
        )

    parent_collision_keys = set(parent.collision_by_identity)
    child_collision_keys = {
        identity
        for identity in child.collision_by_identity
        if identity[1] in parent.levels_by_label
        and identity[2] in parent.levels_by_label
    }
    missing_collisions = parent_collision_keys - child_collision_keys
    added_collisions_between_parent_levels = (
        child_collision_keys - parent_collision_keys
    )
    collision_token_mismatches = 0
    reoriented_collisions = 0
    for identity in parent_collision_keys & child_collision_keys:
        collider = identity[0]
        parent_row = parent.collision_by_identity[identity]
        child_row = child.collision_by_identity[identity]
        parent_temperatures = parent.collision_temperatures_by_collider[collider]
        child_temperatures = child.collision_temperatures_by_collider[collider]
        common_count = len(parent_temperatures)
        if child_temperatures[:common_count] != parent_temperatures:
            collision_token_mismatches += 1
            continue
        if (
            child_row.upper_label,
            child_row.lower_label,
        ) == (parent_row.upper_label, parent_row.lower_label):
            expected_tokens = parent_row.rate_tokens
        elif (
            child_row.upper_label,
            child_row.lower_label,
        ) == (parent_row.lower_label, parent_row.upper_label):
            reoriented_collisions += 1
            source_upper = parent.levels_by_label[parent_row.upper_label]
            source_lower = parent.levels_by_label[parent_row.lower_label]
            source_gap = source_upper.energy_cm1 - source_lower.energy_cm1
            ratio = source_upper.weight / source_lower.weight
            expected_tokens = tuple(
                f"{(
                    rate * ratio
                    if source_gap == 0.0
                    else rate
                    * ratio
                    * math.exp(
                        -source_gap * 1.4387768775039338 / temperature
                    )
                ):.8e}"
                for rate, temperature in zip(
                    parent_row.rates_cm3_s, parent_temperatures
                )
            )
        else:  # pragma: no cover - unordered endpoint key makes this unreachable
            collision_token_mismatches += 1
            continue
        if child_row.rate_tokens[:common_count] != expected_tokens:
            collision_token_mismatches += 1
    if missing_collisions or collision_token_mismatches:
        _failure(
            failures,
            "tier_collision_nesting",
            child_id,
            "parent collision failures: "
            f"missing={len(missing_collisions)}, common-grid token mismatches="
            f"{collision_token_mismatches}",
        )

    inheritance = child_manifest_entry.get("radiative_inheritance")
    if not isinstance(inheritance, dict):
        _failure(failures, "upstream_manifest_binding", child_id, "radiative_inheritance record is absent")
    else:
        checks = {
            "upstream_path": str(parent_spec["path"]),
            "upstream_sha256": parent_file_sha256,
            "upstream_manifest_path": str(parent_spec["manifest"]),
            "upstream_manifest_sha256": parent_manifest_sha256,
            "inherited_radiative_rows": len(parent.radiative_by_endpoints),
        }
        for field, expected in checks.items():
            if inheritance.get(field) != expected:
                _failure(
                    failures,
                    "upstream_manifest_binding",
                    child_id,
                    f"{field}={inheritance.get(field)!r}, expected {expected!r}",
                )
        inherited_fields = set(inheritance.get("inherited_fields", []))
        if inherited_fields != {"a_s", "frequency_ghz"}:
            _failure(
                failures,
                "upstream_manifest_binding",
                child_id,
                f"inherited_fields={sorted(inherited_fields)!r}, expected ['a_s', 'frequency_ghz']",
            )

    collision_inheritance = child_manifest_entry.get("collision_inheritance")
    if isinstance(collision_inheritance, dict):
        collision_checks = {
            "upstream_path": str(parent_spec["path"]),
            "upstream_sha256": parent_file_sha256,
            "upstream_manifest_path": str(parent_spec["manifest"]),
            "upstream_manifest_sha256": parent_manifest_sha256,
            "inherited_collision_rows": len(parent.collision_by_identity),
        }
        for field, expected in collision_checks.items():
            if collision_inheritance.get(field) != expected:
                _failure(
                    failures,
                    "upstream_collision_manifest_binding",
                    child_id,
                    f"{field}={collision_inheritance.get(field)!r}, expected {expected!r}",
                )
        inherited_fields = set(
            collision_inheritance.get("inherited_fields", [])
        )
        if inherited_fields != {
            "common_temperatures_k",
            "rates_cm3_s_on_common_grid",
        }:
            _failure(
                failures,
                "upstream_collision_manifest_binding",
                child_id,
                f"collision inherited_fields={sorted(inherited_fields)!r}",
            )

    return {
        "parent": str(parent_spec["id"]),
        "child": child_id,
        "parent_levels": len(parent.levels_by_label),
        "missing_parent_levels": missing_levels,
        "level_energy_mismatches": energy_mismatches,
        "level_weight_token_mismatches": weight_mismatches,
        "maximum_level_energy_delta_cm1": max_energy_delta,
        "parent_radiative_rows": len(parent.radiative_by_endpoints),
        "missing_parent_radiative_rows": missing_radiative,
        "a_frequency_token_mismatches": token_mismatches,
        "experimental_overlay_rows_present": len(overlay_endpoints) - missing_overlays,
        "parent_collision_rows": len(parent.collision_by_identity),
        "missing_parent_collision_rows": len(missing_collisions),
        "collision_common_grid_token_mismatches": collision_token_mismatches,
        "collision_rows_added_between_parent_levels": len(
            added_collisions_between_parent_levels
        ),
        "collision_rows_reoriented_by_detailed_balance": reoriented_collisions,
        "collision_manifest_inheritance_record_present": isinstance(
            collision_inheritance, dict
        ),
    }


def _product_metrics(
    product_id: str,
    path: Path,
    manifest_path: Path,
    parsed: ParsedProduct,
    *,
    file_sha256: str,
    manifest_sha256: str,
    frequency_metrics: Mapping[str, object],
    sidecar_metrics: Mapping[str, object] | None,
) -> dict[str, object]:
    manifolds: dict[str, int] = {}
    for level in parsed.levels_by_label.values():
        manifolds[level.label.manifold] = manifolds.get(level.label.manifold, 0) + 1
    record: dict[str, object] = {
        "id": product_id,
        "path": str(path),
        "sha256": file_sha256,
        "manifest": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "molecule": parsed.molecule,
        "mass_amu": parsed.mass_amu,
        "levels": len(parsed.levels_by_label),
        "levels_by_manifold": dict(sorted(manifolds.items())),
        "radiative_rows": len(parsed.radiative_by_endpoints),
        "collision_rows_by_collider": {
            str(key): value for key, value in sorted(parsed.collision_rows_by_collider.items())
        },
        "collision_temperatures_k_by_collider": {
            str(key): value
            for key, value in sorted(parsed.collision_temperatures_by_collider.items())
        },
        "frequency_energy_gate": dict(frequency_metrics),
    }
    if sidecar_metrics is not None:
        record["radiative_sidecar"] = dict(sidecar_metrics)
    return record


def audit_repository(
    root: str | Path,
    config_path: str | Path = "config/postbuild-audit.json",
) -> dict[str, object]:
    """Audit every explicitly configured product and return a JSON-ready report."""
    root_path = Path(root).resolve()
    config_file = Path(config_path)
    if not config_file.is_absolute():
        config_file = root_path / config_file
    config_file = config_file.resolve()
    config = json.loads(config_file.read_text(encoding="utf-8"))
    failures: list[dict[str, object]] = []
    product_metrics: list[dict[str, object]] = []
    edge_metrics: list[dict[str, object]] = []
    audited_manifest_inputs: set[Path] = set()
    manifest_declared_input_hashes_checked = 0
    manifest_radiative_integrity_records_checked = 0
    configured_products = sum(len(chain.get("products", [])) for chain in config.get("chains", []))

    for chain in config.get("chains", []):
        species = str(chain["species"])
        isotopologue = str(chain["isotopologue"])
        expected_mass = float(config["isotopologue_masses_amu"][isotopologue])
        overlay_expected = int(chain.get("overlay_expected", 0))
        overlay_endpoints: set[tuple[str, str]] = set()
        previous: ParsedProduct | None = None
        previous_spec: Mapping[str, Any] | None = None
        previous_file_sha256 = ""
        previous_manifest_sha256 = ""
        chain_masses: list[float] = []

        for position, product_spec in enumerate(chain.get("products", [])):
            product_id = str(product_spec["id"])
            path = root_path / str(product_spec["path"])
            manifest_path = root_path / str(product_spec["manifest"])
            try:
                file_sha256 = _sha256(path)
                entry, manifest_sha256 = _manifest_entry(
                    root_path,
                    manifest_path,
                    path,
                    product_id=product_id,
                    species=species,
                    file_sha256=file_sha256,
                    failures=failures,
                )
                if manifest_path not in audited_manifest_inputs:
                    manifest_declared_input_hashes_checked += (
                        _manifest_declared_input_hash_gate(
                            root_path,
                            manifest_path,
                            product_id=product_id,
                            failures=failures,
                        )
                    )
                    audited_manifest_inputs.add(manifest_path)
                tier = config["tiers"][product_spec["tier"]]
                parsed = _parse_product(
                    path,
                    product_id=product_id,
                    species=species,
                    isotopologue=isotopologue,
                    expected_mass=expected_mass,
                    tier=tier,
                    failures=failures,
                )
            except (OSError, UnicodeError, AuditParseError, KeyError, TypeError, ValueError) as exc:
                _failure(failures, "product_parse", product_id, str(exc))
                previous = None
                previous_spec = None
                continue

            chain_masses.append(parsed.mass_amu)
            expected_levels = int(tier["levels"])
            expected_radiative = int(product_spec["radiative"])
            manifest_counts = {
                "levels": entry.get("levels"),
                "radiative_transitions": entry.get("radiative_transitions"),
            }
            actual_counts = {
                "levels": len(parsed.levels_by_label),
                "radiative_transitions": len(parsed.radiative_by_endpoints),
            }
            configured_counts = {
                "levels": expected_levels,
                "radiative_transitions": expected_radiative,
            }
            for field in ("levels", "radiative_transitions"):
                if actual_counts[field] != configured_counts[field]:
                    _failure(
                        failures,
                        "expected_counts",
                        product_id,
                        f"actual {field}={actual_counts[field]}, expected {configured_counts[field]}",
                    )
                if manifest_counts[field] != actual_counts[field]:
                    _failure(
                        failures,
                        "manifest_binding",
                        product_id,
                        f"manifest {field}={manifest_counts[field]!r}, actual {actual_counts[field]}",
                    )
            radiative_integrity = entry.get("radiative_integrity")
            if not isinstance(radiative_integrity, dict):
                _failure(
                    failures,
                    "manifest_radiative_integrity",
                    product_id,
                    "radiative_integrity record is absent",
                )
            else:
                manifest_radiative_integrity_records_checked += 1
                expected_integrity = {
                    "status": "pass",
                    "levels": len(parsed.levels_by_label),
                    "radiative_rows": len(parsed.radiative_by_endpoints),
                    "forbidden_rows": 0,
                    "duplicate_or_reversed_endpoints": 0,
                    "radiative_components": 1,
                }
                mismatched_integrity = {
                    field: (radiative_integrity.get(field), expected)
                    for field, expected in expected_integrity.items()
                    if radiative_integrity.get(field) != expected
                }
                if mismatched_integrity:
                    _failure(
                        failures,
                        "manifest_radiative_integrity",
                        product_id,
                        f"radiative_integrity mismatch {mismatched_integrity}",
                    )
            actual_collision_counts = {
                str(key): value for key, value in sorted(parsed.collision_rows_by_collider.items())
            }
            actual_collision_total = sum(actual_collision_counts.values())
            if "collision_rows" in entry and entry["collision_rows"] != actual_collision_total:
                _failure(
                    failures,
                    "manifest_binding",
                    product_id,
                    f"manifest collision_rows={entry['collision_rows']!r}, actual {actual_collision_total}",
                )
            for field in ("collision_rows_by_collider", "collision_rows_by_partner"):
                if field in entry:
                    manifest_collision_counts = {
                        str(key): int(value) for key, value in entry[field].items()
                    }
                    if manifest_collision_counts != actual_collision_counts:
                        _failure(
                            failures,
                            "manifest_binding",
                            product_id,
                            f"manifest {field}={manifest_collision_counts!r}, actual {actual_collision_counts!r}",
                        )

            sidecar_metrics: dict[str, object] | None = None
            if position == 0:
                sidecar_value = product_spec.get("radiative_sidecar")
                if not sidecar_value:
                    _failure(failures, "sidecar_binding", product_id, "base product has no radiative sidecar")
                else:
                    sidecar_path = root_path / str(sidecar_value)
                    try:
                        overlay_endpoints, sidecar_metrics = _read_overlay_sidecar(
                            sidecar_path,
                            expected_sha256=entry.get("radiative_sidecar_sha256"),
                            expected_count=overlay_expected,
                            product_id=product_id,
                            parsed=parsed,
                            failures=failures,
                        )
                        sidecar_metrics["path"] = str(
                            sidecar_path.relative_to(root_path)
                        )
                    except (OSError, UnicodeError, csv.Error) as exc:
                        _failure(failures, "sidecar_binding", product_id, str(exc))
                frequency_exceptions = overlay_endpoints
            elif previous is not None and previous_spec is not None:
                edge_metrics.append(
                    _edge_gate(
                        root_path,
                        previous_spec,
                        product_spec,
                        previous,
                        parsed,
                        child_manifest_entry=entry,
                        parent_file_sha256=previous_file_sha256,
                        parent_manifest_sha256=previous_manifest_sha256,
                        energy_tolerance_cm1=float(config["edge_energy_tolerance_cm1"]),
                        overlay_endpoints=overlay_endpoints,
                        failures=failures,
                    )
                )
                frequency_exceptions = set(previous.radiative_by_endpoints)
            else:
                _failure(
                    failures,
                    "tier_nesting",
                    product_id,
                    "parent product was unavailable; inherited exceptions cannot be established",
                )
                frequency_exceptions = set()

            frequency_metrics = _frequency_gap_gate(
                parsed,
                exceptions=frequency_exceptions,
                tolerance_ghz=float(config["frequency_gap_tolerance_ghz"]),
                cm1_to_ghz=float(config["cm1_to_ghz"]),
                product_id=product_id,
                failures=failures,
            )
            product_metrics.append(
                _product_metrics(
                    product_id,
                    Path(str(product_spec["path"])),
                    Path(str(product_spec["manifest"])),
                    parsed,
                    file_sha256=file_sha256,
                    manifest_sha256=manifest_sha256,
                    frequency_metrics=frequency_metrics,
                    sidecar_metrics=sidecar_metrics,
                )
            )
            previous = parsed
            previous_spec = product_spec
            previous_file_sha256 = file_sha256
            previous_manifest_sha256 = manifest_sha256

        if chain_masses and any(value != chain_masses[0] for value in chain_masses[1:]):
            _failure(
                failures,
                "isotopologue_mass",
                str(chain["id"]),
                f"tier masses are inconsistent: {chain_masses!r}",
            )

    gate_failure_counts: dict[str, int] = {}
    for failure in failures:
        gate = str(failure["gate"])
        gate_failure_counts[gate] = gate_failure_counts.get(gate, 0) + 1
    expected_edges = sum(max(0, len(chain.get("products", [])) - 1) for chain in config.get("chains", []))
    if len(edge_metrics) != expected_edges:
        _failure(
            failures,
            "tier_nesting",
            "repository",
            f"audited {len(edge_metrics)} tier edges, expected {expected_edges}",
        )
        gate_failure_counts["tier_nesting"] = gate_failure_counts.get("tier_nesting", 0) + 1

    return {
        "schema_version": 1,
        "status": "pass" if not failures else "fail",
        "root": ".",
        "config": str(config_file.relative_to(root_path)),
        "config_sha256": _sha256(config_file),
        "configured_products": configured_products,
        "audited_products": len(product_metrics),
        "configured_tier_edges": expected_edges,
        "audited_tier_edges": len(edge_metrics),
        "manifest_declared_input_hashes_checked": (
            manifest_declared_input_hashes_checked
        ),
        "manifest_radiative_integrity_records_checked": (
            manifest_radiative_integrity_records_checked
        ),
        "products": product_metrics,
        "tier_edges": edge_metrics,
        "gate_failure_counts": dict(sorted(gate_failure_counts.items())),
        "failures": failures,
    }


def _default_root() -> Path:
    return Path(__file__).resolve().parents[2]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=_default_root())
    parser.add_argument("--config", default="config/postbuild-audit.json")
    parser.add_argument("--output", type=Path, help="write the complete JSON report here")
    arguments = parser.parse_args(argv)
    try:
        report = audit_repository(arguments.root, arguments.config)
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        report = {
            "schema_version": 1,
            "status": "fail",
            "failures": [{"gate": "configuration", "product": "repository", "message": str(exc)}],
        }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        sys.stdout.write(rendered)
    else:
        output = arguments.output
        if not output.is_absolute():
            output = Path(arguments.root) / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    return 0 if report.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
