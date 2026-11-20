"""Reviewed selection of observed frequencies from the Xu LINECAL crosswalk."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Mapping

from .xu import CM1_TO_MHZ


SELECTION_POLICY = (
    "exact endpoint identity; fitted Xu rows only; minimum explicit LINECAL UNC; "
    "conflicting equal-uncertainty measurements require review"
)


class ObservedFrequencySelectionError(ValueError):
    """Raised rather than making an arbitrary observed-frequency choice."""


@dataclass(frozen=True, slots=True)
class XuObservedFrequency:
    line_transition_index: int
    species: str
    upper_key: str
    lower_key: str
    frequency_mhz: float
    uncertainty_mhz: float
    source_line: int
    candidate_count: int

    @property
    def source_id(self) -> str:
        return f"xu2008-globalfit-observation:line-{self.source_line}"


def _uncertainty_mhz(unit: str, code: str) -> float:
    try:
        value = int(code)
    except (TypeError, ValueError) as exc:
        raise ObservedFrequencySelectionError(
            f"mapped fitted Xu observation has invalid UNC field {code!r}"
        ) from exc
    if value <= 0:
        raise ObservedFrequencySelectionError(
            f"mapped fitted Xu observation has non-positive UNC field {code!r}"
        )
    if unit == "MHz":
        return value * 0.001
    if unit == "cm-1":
        return value * 0.00001 * CM1_TO_MHZ
    raise ObservedFrequencySelectionError(f"unsupported Xu observation unit {unit!r}")


def select_xu_observed_frequencies(
    rows: Iterable[Mapping[str, str]],
) -> tuple[dict[int, XuObservedFrequency], dict[str, object]]:
    """Choose one safe observed frequency for each exactly mapped Xu line.

    Excluded measurements and unresolved identities remain provenance only.
    Repeated fitted measurements are ranked by the explicit ``UNC`` printed by
    GLOBALFIT.  A conflicting tie is an error, because source order is not a
    scientific frequency-selection rule.
    """
    groups: dict[int, list[XuObservedFrequency]] = defaultdict(list)
    row_counts: Counter[str] = Counter()
    for row in rows:
        row_counts["source_rows"] += 1
        mapping_status = row.get("mapping_status", "")
        fit_status = row.get("fit_status", "")
        if mapping_status != "mapped_to_retained_line":
            row_counts[f"ignored_mapping_{mapping_status or 'missing'}"] += 1
            continue
        if fit_status != "fitted":
            row_counts[f"ignored_fit_{fit_status or 'missing'}"] += 1
            continue
        try:
            line_index = int(row.get("line_transition_index", ""))
            source_line = int(row.get("source_line", ""))
            frequency_mhz = float(row.get("observed_mhz", ""))
        except (TypeError, ValueError) as exc:
            raise ObservedFrequencySelectionError(
                "mapped fitted Xu observation has malformed line or frequency identity"
            ) from exc
        if line_index <= 0 or source_line <= 0 or not math.isfinite(frequency_mhz) or frequency_mhz <= 0:
            raise ObservedFrequencySelectionError(
                "mapped fitted Xu observation has non-positive or non-finite identity/frequency"
            )
        species = row.get("species", "").upper()
        upper_key = row.get("generated_upper_key", "")
        lower_key = row.get("generated_lower_key", "")
        if species not in {"A", "E"} or not upper_key or not lower_key:
            raise ObservedFrequencySelectionError(
                f"mapped fitted Xu observation {source_line} lacks canonical endpoints"
            )
        groups[line_index].append(
            XuObservedFrequency(
                line_transition_index=line_index,
                species=species,
                upper_key=upper_key,
                lower_key=lower_key,
                frequency_mhz=frequency_mhz,
                uncertainty_mhz=_uncertainty_mhz(
                    row.get("unit", ""), row.get("weight_code", "")
                ),
                source_line=source_line,
                candidate_count=0,
            )
        )
        row_counts["fitted_exact_measurements"] += 1

    selected: dict[int, XuObservedFrequency] = {}
    single_count = 0
    repeated_count = 0
    for line_index, candidates in groups.items():
        identities = {
            (candidate.species, candidate.upper_key, candidate.lower_key)
            for candidate in candidates
        }
        if len(identities) != 1:
            raise ObservedFrequencySelectionError(
                f"Xu line index {line_index} maps to inconsistent canonical endpoints"
            )
        minimum_uncertainty = min(candidate.uncertainty_mhz for candidate in candidates)
        best = [
            candidate for candidate in candidates
            if math.isclose(
                candidate.uncertainty_mhz,
                minimum_uncertainty,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
        ]
        frequencies = {candidate.frequency_mhz for candidate in best}
        if len(frequencies) != 1:
            details = ", ".join(
                f"source line {candidate.source_line}: {candidate.frequency_mhz:.12g} MHz"
                for candidate in sorted(best, key=lambda item: item.source_line)
            )
            raise ObservedFrequencySelectionError(
                f"Xu line index {line_index} has conflicting equally precise observations ({details})"
            )
        chosen = min(best, key=lambda candidate: candidate.source_line)
        selected[line_index] = XuObservedFrequency(
            line_transition_index=chosen.line_transition_index,
            species=chosen.species,
            upper_key=chosen.upper_key,
            lower_key=chosen.lower_key,
            frequency_mhz=chosen.frequency_mhz,
            uncertainty_mhz=chosen.uncertainty_mhz,
            source_line=chosen.source_line,
            candidate_count=len(candidates),
        )
        if len(candidates) == 1:
            single_count += 1
        else:
            repeated_count += 1

    summary: dict[str, object] = {
        "policy": SELECTION_POLICY,
        "source_rows": row_counts["source_rows"],
        "fitted_exact_measurements": row_counts["fitted_exact_measurements"],
        "selected_line_count": len(selected),
        "single_measurement_lines": single_count,
        "repeated_measurement_lines": repeated_count,
        "lower_precision_measurements_not_selected": (
            row_counts["fitted_exact_measurements"] - len(selected)
        ),
        "ignored_rows_by_reason": {
            key.removeprefix("ignored_"): value
            for key, value in sorted(row_counts.items())
            if key.startswith("ignored_")
        },
    }
    return selected, summary


def read_xu_observed_frequencies(
    path: str | Path,
) -> tuple[dict[int, XuObservedFrequency], dict[str, object]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return select_xu_observed_frequencies(csv.DictReader(handle))
