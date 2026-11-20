"""CDMS 033502 frequency imports for :sup:`13`CH3OH.

The current CDMS entry is a mixed catalog: calculated rows have a positive
species tag, while laboratory-frequency replacements have a negative tag.
Only rows tagged ``-033502`` are eligible here.  The module deliberately does
not import CDMS intensities or level energies into the isotope builder.

CDMS reserves six two-character quantum-number slots for each state even when
``QNFMT`` declares fewer.  In the current QNFMT=1404 catalog the four fields
are ``J, signed Ka, Kc, parity+vt``; older fixtures may spell parity and vt
separately.  The parity prefix is absent for E symmetry, so removing positional
blanks would change the identity of following fields.
"""

from __future__ import annotations

from collections.abc import Hashable, Iterable, Mapping
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import TypeVar

from ...quantum import LevelKey


CDMS_TAG = 33502
CDMS_TAG_TEXT = "033502"
CDMS_DATASET_ID = "cdms-13ch3oh-033502-v2"
SELECTION_POLICY = (
    "negative CDMS tag -033502 only; exact documented endpoint identity; "
    "minimum uncertainty; source-row tie break"
)

CrosswalkId = TypeVar("CrosswalkId", bound=Hashable)


class CH3OH13CDMSError(ValueError):
    """Raised when a CDMS row cannot be interpreted without guessing."""


@dataclass(frozen=True, slots=True)
class CDMS033502Row:
    """One fixed-width CDMS catalog row with positional QN blanks retained."""

    frequency_mhz: float
    uncertainty_mhz: float
    log10_intensity: float
    degrees_of_freedom: int
    lower_energy_cm1: float
    upper_degeneracy: int
    tag: int
    qn_format: int
    upper_qns: tuple[str, ...]
    lower_qns: tuple[str, ...]
    source_row: int
    raw: str

    @property
    def quantum_numbers_per_state(self) -> int:
        return abs(self.qn_format) % 10

    @property
    def is_experimental_033502(self) -> bool:
        return self.tag == -CDMS_TAG


@dataclass(frozen=True, slots=True, order=True)
class CDMSTransitionIdentity:
    """Ordered catalog endpoints for one radiative transition."""

    upper_key: LevelKey
    lower_key: LevelKey

    @property
    def unordered_endpoints(self) -> frozenset[LevelKey]:
        return frozenset((self.upper_key, self.lower_key))


@dataclass(frozen=True, slots=True)
class CDMSFrequencyOverride:
    """The selected laboratory frequency for one exact endpoint pair."""

    identity: CDMSTransitionIdentity
    frequency_mhz: float
    uncertainty_mhz: float
    source_row: int
    candidate_count: int
    tag: int = -CDMS_TAG

    @property
    def source_id(self) -> str:
        return f"{CDMS_DATASET_ID}:experimental-row-{self.source_row}"

    def as_json_record(self, *, transition_id: Hashable | None = None) -> dict[str, object]:
        record: dict[str, object] = {
            "species": self.identity.upper_key.species,
            "upper_key": self.identity.upper_key.encode(),
            "lower_key": self.identity.lower_key.encode(),
            "frequency_mhz": self.frequency_mhz,
            "uncertainty_mhz": self.uncertainty_mhz,
            "source_row": self.source_row,
            "source": self.source_id,
            "source_id": self.source_id,
            "candidate_count": self.candidate_count,
            "tag": self.tag,
        }
        if transition_id is not None:
            record["transition_id"] = str(transition_id)
        return record


def _qn_count(qn_format: int) -> int:
    count = abs(qn_format) % 10
    if not 1 <= count <= 6:
        raise CH3OH13CDMSError(
            f"QNFMT {qn_format} declares {count} quantum numbers per state; expected 1-6"
        )
    return count


def _fixed_qns(raw: str, qn_format: int) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split the two six-slot states without collapsing internal blanks."""
    count = _qn_count(qn_format)
    field = raw[55:79].ljust(24)
    upper_field = field[:12]
    lower_field = field[12:24]

    def state_values(value: str) -> tuple[str, ...]:
        return tuple(value[offset : offset + 2].strip() for offset in range(0, 2 * count, 2))

    return state_values(upper_field), state_values(lower_field)


def parse_cdms_033502_line(line: str, *, source_row: int = 1) -> CDMS033502Row:
    """Parse one ordinary JPL/CDMS fixed-width catalog line.

    Trailing unused quantum-number blanks may have been stripped by a text
    transport; they are restored before the two six-slot state fields are
    split.  Internal blanks, including the E-symmetry parity slot, are never
    removed.
    """
    raw = line.rstrip("\r\n")
    if len(raw) < 55:
        raise CH3OH13CDMSError("catalog line shorter than the 55-column header")
    if source_row <= 0:
        raise CH3OH13CDMSError("source row must be positive")
    try:
        qn_format = int(raw[51:55])
        upper_qns, lower_qns = _fixed_qns(raw, qn_format)
        row = CDMS033502Row(
            frequency_mhz=float(raw[0:13]),
            uncertainty_mhz=float(raw[13:21]),
            log10_intensity=float(raw[21:29]),
            degrees_of_freedom=int(raw[29:31]),
            lower_energy_cm1=float(raw[31:41]),
            upper_degeneracy=int(raw[41:44]),
            tag=int(raw[44:51]),
            qn_format=qn_format,
            upper_qns=upper_qns,
            lower_qns=lower_qns,
            source_row=source_row,
            raw=raw,
        )
    except (TypeError, ValueError) as exc:
        raise CH3OH13CDMSError(f"malformed fixed-width catalog row: {exc}") from exc
    numeric_values = (
        row.frequency_mhz,
        row.uncertainty_mhz,
        row.log10_intensity,
        row.lower_energy_cm1,
    )
    if not all(math.isfinite(value) for value in numeric_values):
        raise CH3OH13CDMSError("catalog row contains a non-finite numeric field")
    if row.frequency_mhz <= 0.0:
        raise CH3OH13CDMSError("catalog frequency must be positive")
    if row.uncertainty_mhz < 0.0:
        raise CH3OH13CDMSError("catalog uncertainty must be nonnegative")
    return row


def parse_cdms_033502(lines: Iterable[str]) -> list[CDMS033502Row]:
    """Parse non-comment rows while retaining their physical source numbers."""
    parsed: list[CDMS033502Row] = []
    for source_row, line in enumerate(lines, start=1):
        if not line.strip() or line.lstrip().startswith(("#", "!")):
            continue
        try:
            parsed.append(parse_cdms_033502_line(line, source_row=source_row))
        except CH3OH13CDMSError as exc:
            raise CH3OH13CDMSError(f"catalog source row {source_row}: {exc}") from exc
    return parsed


def read_cdms_033502(path: str | Path) -> list[CDMS033502Row]:
    return parse_cdms_033502(Path(path).read_text(encoding="ascii").splitlines())


def _integer_quantum_number(token: str, *, field: str) -> int:
    """Decode the ordinary integer and CDMS two-character overflow forms."""
    try:
        return int(token)
    except ValueError:
        pass
    if len(token) == 2 and token[1].isdigit():
        digit = int(token[1])
        if "a" <= token[0] <= "z":
            return -(10 * (ord(token[0]) - ord("a") + 1) + digit)
        if "A" <= token[0] <= "Z":
            return 100 + 10 * (ord(token[0]) - ord("A")) + digit
    raise CH3OH13CDMSError(f"unsupported {field} quantum number {token!r}")


def state_qns_to_level_key(qns: tuple[str, ...]) -> LevelKey:
    """Map documented 033502 labels to one canonical level identity."""
    if len(qns) == 4:
        j_token, ka_token, kc_token, parity_vt = qns
        if not parity_vt:
            raise CH3OH13CDMSError("CDMS 033502 combined parity/vt field is blank")
        parity = parity_vt[:-1]
        vt_token = parity_vt[-1]
    elif len(qns) == 5:
        j_token, ka_token, kc_token, parity, vt_token = qns
    else:
        raise CH3OH13CDMSError(
            "CDMS 033502 state needs four combined or five expanded quantum "
            f"numbers; found {len(qns)}"
        )
    j = _integer_quantum_number(j_token, field="J")
    signed_ka = _integer_quantum_number(ka_token, field="Ka")
    kc = _integer_quantum_number(kc_token, field="Kc")
    vt = _integer_quantum_number(vt_token, field="vt")
    if vt not in {0, 1}:
        raise CH3OH13CDMSError(
            f"CDMS tag {CDMS_TAG_TEXT} documents vt=0,1 only; found vt={vt}"
        )
    if j < 0 or not 0 <= kc <= j or abs(signed_ka) > j:
        raise CH3OH13CDMSError(
            f"invalid 033502 rotational labels J={j}, Ka={signed_ka}, Kc={kc}"
        )
    try:
        if parity in {"+", "-"}:
            return LevelKey("A", vt, j, abs(signed_ka), parity)
        if parity == "":
            return LevelKey("E", vt, j, signed_ka)
    except ValueError as exc:
        raise CH3OH13CDMSError(f"invalid 033502 level labels {qns!r}: {exc}") from exc
    raise CH3OH13CDMSError(
        f"033502 parity must be '+', '-', or blank; found {parity!r}"
    )


def row_transition_identity(row: CDMS033502Row) -> CDMSTransitionIdentity:
    """Return the exact canonical endpoints printed on a catalog row."""
    if row.quantum_numbers_per_state not in {4, 5}:
        raise CH3OH13CDMSError(
            f"CDMS 033502 row {row.source_row} has incompatible QNFMT {row.qn_format}"
        )
    upper = state_qns_to_level_key(row.upper_qns)
    lower = state_qns_to_level_key(row.lower_qns)
    if upper.species != lower.species:
        raise CH3OH13CDMSError(
            f"CDMS 033502 row {row.source_row} changes A/E symmetry"
        )
    if upper == lower:
        raise CH3OH13CDMSError(
            f"CDMS 033502 row {row.source_row} has identical endpoints"
        )
    return CDMSTransitionIdentity(upper, lower)


def select_experimental_overrides(
    rows: Iterable[CDMS033502Row],
) -> dict[CDMSTransitionIdentity, CDMSFrequencyOverride]:
    """Select one negative-tag laboratory frequency per exact transition.

    Positive 033502 predictions and rows belonging to every other tag are
    ignored.  Repeated laboratory rows are ranked first by the reported
    uncertainty and then by physical source-row number.  Frequency and raw
    spelling are tertiary keys only to make pathological same-row inputs
    independent of iterable order.
    """
    grouped: dict[CDMSTransitionIdentity, list[CDMS033502Row]] = {}
    for row in rows:
        if not row.is_experimental_033502:
            continue
        identity = row_transition_identity(row)
        grouped.setdefault(identity, []).append(row)

    selected: dict[CDMSTransitionIdentity, CDMSFrequencyOverride] = {}
    for identity, candidates in grouped.items():
        chosen = min(
            candidates,
            key=lambda row: (
                row.uncertainty_mhz,
                row.source_row,
                row.frequency_mhz,
                row.raw,
            ),
        )
        selected[identity] = CDMSFrequencyOverride(
            identity=identity,
            frequency_mhz=chosen.frequency_mhz,
            uncertainty_mhz=chosen.uncertainty_mhz,
            source_row=chosen.source_row,
            candidate_count=len(candidates),
        )
    return selected


def catalog_file_to_overrides(
    path: str | Path,
) -> dict[CDMSTransitionIdentity, CDMSFrequencyOverride]:
    """Read a local c033502-style file and select its safe overrides."""
    return select_experimental_overrides(read_cdms_033502(path))


def crosswalk_experimental_overrides(
    overrides: Mapping[CDMSTransitionIdentity, CDMSFrequencyOverride],
    transition_endpoints: Mapping[CrosswalkId, tuple[LevelKey, LevelKey]],
) -> dict[CrosswalkId, CDMSFrequencyOverride]:
    """Crosswalk product transition IDs by endpoint identity.

    Frequencies are direction-independent, so the product endpoint order may
    be reversed.  Ambiguous duplicate unordered CDMS identities are rejected.
    """
    by_endpoints: dict[frozenset[LevelKey], CDMSFrequencyOverride] = {}
    for identity, override in overrides.items():
        endpoint_set = identity.unordered_endpoints
        if len(endpoint_set) != 2:
            raise CH3OH13CDMSError("CDMS override must have two distinct endpoints")
        if endpoint_set in by_endpoints and by_endpoints[endpoint_set] != override:
            raise CH3OH13CDMSError("ambiguous CDMS overrides for one endpoint pair")
        by_endpoints[endpoint_set] = override

    matched: dict[CrosswalkId, CDMSFrequencyOverride] = {}
    for transition_id, endpoints in transition_endpoints.items():
        if len(endpoints) != 2:
            raise CH3OH13CDMSError(
                f"crosswalk transition {transition_id!r} must have two endpoints"
            )
        endpoint_set = frozenset(endpoints)
        override = by_endpoints.get(endpoint_set)
        if override is not None:
            matched[transition_id] = override
    return matched


def overrides_to_json(
    overrides: Mapping[CDMSTransitionIdentity, CDMSFrequencyOverride],
    *,
    transition_endpoints: Mapping[CrosswalkId, tuple[LevelKey, LevelKey]] | None = None,
) -> str:
    """Serialize selected overrides, optionally restricted by a crosswalk."""
    records: list[dict[str, object]]
    if transition_endpoints is None:
        records = [
            override.as_json_record()
            for _, override in sorted(
                overrides.items(),
                key=lambda item: (
                    item[0].upper_key.encode(),
                    item[0].lower_key.encode(),
                ),
            )
        ]
    else:
        matched = crosswalk_experimental_overrides(overrides, transition_endpoints)
        records = [
            matched[transition_id].as_json_record(transition_id=transition_id)
            for transition_id in sorted(matched, key=lambda value: str(value))
        ]
    payload = {
        "dataset_id": CDMS_DATASET_ID,
        "catalog_tag": CDMS_TAG_TEXT,
        "selection_policy": SELECTION_POLICY,
        "unit": "MHz",
        "override_count": len(records),
        "overrides": records,
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def catalog_file_to_override_json(
    catalog_path: str | Path,
    *,
    output_path: str | Path | None = None,
    transition_endpoints: Mapping[CrosswalkId, tuple[LevelKey, LevelKey]] | None = None,
) -> str:
    """Convert a local catalog into deterministic override JSON.

    The JSON text is always returned.  Supplying ``output_path`` additionally
    writes that exact text, allowing callers to retain a reviewed crosswalk
    artifact without coupling this module to the full product builder.
    """
    text = overrides_to_json(
        catalog_file_to_overrides(catalog_path),
        transition_endpoints=transition_endpoints,
    )
    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text, encoding="utf-8")
    return text


def rescale_einstein_a_nu3(
    einstein_a_s: float,
    calculated_frequency: float,
    adopted_frequency: float,
) -> float:
    """Rescale ``A_ul`` when only the adopted frequency changes.

    The two frequencies may use any common unit.  With the line strength and
    upper-state degeneracy fixed, ``A_ul`` scales exactly as frequency cubed.
    """
    try:
        a_value = float(einstein_a_s)
        calculated = float(calculated_frequency)
        adopted = float(adopted_frequency)
    except (TypeError, ValueError) as exc:
        raise TypeError("Einstein A and frequencies must be real numbers") from exc
    if not all(math.isfinite(value) for value in (a_value, calculated, adopted)):
        raise ValueError("Einstein A and frequencies must be finite")
    if a_value < 0.0:
        raise ValueError("Einstein A must be nonnegative")
    if calculated <= 0.0 or adopted <= 0.0:
        raise ValueError("frequencies must be positive")
    return a_value * (adopted / calculated) ** 3


__all__ = [
    "CDMS_DATASET_ID",
    "CDMS_TAG",
    "CDMS_TAG_TEXT",
    "SELECTION_POLICY",
    "CH3OH13CDMSError",
    "CDMS033502Row",
    "CDMSFrequencyOverride",
    "CDMSTransitionIdentity",
    "catalog_file_to_override_json",
    "catalog_file_to_overrides",
    "crosswalk_experimental_overrides",
    "overrides_to_json",
    "parse_cdms_033502",
    "parse_cdms_033502_line",
    "read_cdms_033502",
    "rescale_einstein_a_nu3",
    "row_transition_identity",
    "select_experimental_overrides",
    "state_qns_to_level_key",
]
