"""Safe extraction of A/E microwave and infrared rows from Xu/BELGI output."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
from pathlib import Path
import re

from .quantum import LevelKey


CM1_TO_MHZ = 29979.2458


@dataclass(frozen=True, slots=True)
class XuLinecalRow:
    source_line: int
    species: str
    upper_vt: int
    upper_j: int
    upper_k: int
    upper_kc: int
    upper_parity: str
    lower_vt: int
    lower_j: int
    lower_k: int
    lower_kc: int
    lower_parity: str
    observation: float
    calculation: float
    obs_minus_calc: float | None
    unit: str
    fit_status: str
    weight_code: int | None
    comment: str

    @property
    def upper_key(self) -> LevelKey:
        return LevelKey(
            self.species, self.upper_vt, self.upper_j, self.upper_k, self.upper_parity
        )

    @property
    def lower_key(self) -> LevelKey:
        return LevelKey(
            self.species, self.lower_vt, self.lower_j, self.lower_k, self.lower_parity
        )

    @property
    def observed_mhz(self) -> float:
        return self.observation if self.unit == "MHz" else self.observation * CM1_TO_MHZ

    @property
    def calculated_mhz(self) -> float:
        return self.calculation if self.unit == "MHz" else self.calculation * CM1_TO_MHZ

    @property
    def uncertainty_mhz(self) -> float | None:
        """Return Xu's parenthetical ``UNC`` field in MHz.

        GLOBALFIT prints microwave observations to 0.001 MHz and infrared
        observations to 0.00001 cm-1.  The parenthetical integer is an
        uncertainty in those last printed digits.  Formally excluded rows do
        not carry the field and therefore return ``None``.
        """
        if self.weight_code is None:
            return None
        if self.unit == "MHz":
            return self.weight_code * 0.001
        return self.weight_code * 0.00001 * CM1_TO_MHZ


_NUMBER = r"[+-]?\d+(?:\.\d+)?"
_RESIDUAL = rf"(?:{_NUMBER}|\*{{3,}})"
_A_LINE = re.compile(
    rf"^\s*([012])\s+(\d+)\s+(\d+)\s+(\d+)\s+([+-])\s+"
    rf"([012])\s+(\d+)\s+(\d+)\s+(\d+)\s+([+-])\s+"
    rf"({_NUMBER})\s*(?:\(\s*(\d+)\s*\))?\s+"
    rf"({_NUMBER})\s+({_RESIDUAL})(.*)$"
)
_E_LINE = re.compile(
    rf"^\s*([012])\s+(\d+)\s+(-?\d+)\s+(\d+)\s+"
    rf"([012])\s+(\d+)\s+(-?\d+)\s+(\d+)\s+"
    rf"({_NUMBER})\s*(?:\(\s*(\d+)\s*\))?\s+"
    rf"({_NUMBER})\s+({_RESIDUAL})(.*)$"
)


def _residual(token: str) -> float | None:
    """Return None where the fixed-width residual field overflowed as stars."""
    return None if token.startswith("*") else float(token)


def output_has_intensities(text: str) -> bool:
    match = re.search(r"^\s*Intensity\?\s+([01])\s*$", text, flags=re.MULTILINE | re.IGNORECASE)
    return bool(match and match.group(1) == "1")


def _unit_from_printed_observation(token: str, source_line: int) -> str:
    decimals = len(token.rsplit(".", 1)[1]) if "." in token else 0
    if decimals == 3:
        return "MHz"
    if decimals == 5:
        return "cm-1"
    raise ValueError(
        f"Xu LINECAL row {source_line} has unsupported observation precision {token!r}"
    )


def parse_linecal_rows(text: str) -> list[XuLinecalRow]:
    """Parse all fitted and formally excluded rows from a GLOBALFIT listing.

    A rows carry explicit parity; E rows carry signed K and no parity.  The
    program's fixed output precision distinguishes MHz (three decimals) from
    infrared wavenumber (five decimals).  Parenthetical values are retained as
    raw fit-weight codes rather than mislabelled as universal uncertainties.
    """
    rows: list[XuLinecalRow] = []
    has_marker = "BEGIN LINECAL" in text
    inside = not has_marker
    for source_line, line in enumerate(text.splitlines(), 1):
        if "BEGIN LINECAL" in line:
            inside = True
            continue
        if inside and line.lstrip().startswith("STD.DEV.(UNITLESS)"):
            break
        if not inside:
            continue

        match = _A_LINE.match(line)
        if match:
            values = match.groups()
            observation_token = values[10]
            weight_code = int(values[11]) if values[11] is not None else None
            rows.append(
                XuLinecalRow(
                    source_line=source_line,
                    species="A",
                    upper_vt=int(values[0]), upper_j=int(values[1]),
                    upper_k=int(values[2]), upper_kc=int(values[3]), upper_parity=values[4],
                    lower_vt=int(values[5]), lower_j=int(values[6]),
                    lower_k=int(values[7]), lower_kc=int(values[8]), lower_parity=values[9],
                    observation=float(observation_token), calculation=float(values[12]),
                    obs_minus_calc=_residual(values[13]),
                    unit=_unit_from_printed_observation(observation_token, source_line),
                    fit_status="fitted" if weight_code is not None else "excluded",
                    weight_code=weight_code, comment=values[14].strip(),
                )
            )
            continue

        match = _E_LINE.match(line)
        if match:
            values = match.groups()
            observation_token = values[8]
            weight_code = int(values[9]) if values[9] is not None else None
            rows.append(
                XuLinecalRow(
                    source_line=source_line,
                    species="E",
                    upper_vt=int(values[0]), upper_j=int(values[1]),
                    upper_k=int(values[2]), upper_kc=int(values[3]), upper_parity="",
                    lower_vt=int(values[4]), lower_j=int(values[5]),
                    lower_k=int(values[6]), lower_kc=int(values[7]), lower_parity="",
                    observation=float(observation_token), calculation=float(values[10]),
                    obs_minus_calc=_residual(values[11]),
                    unit=_unit_from_printed_observation(observation_token, source_line),
                    fit_status="fitted" if weight_code is not None else "excluded",
                    weight_code=weight_code, comment=values[12].strip(),
                )
            )
    return rows


def read_linecal_rows(path: str | Path) -> list[XuLinecalRow]:
    return parse_linecal_rows(Path(path).read_text(encoding="utf-8", errors="replace"))


def parse_observed_lines(text: str) -> list[XuLinecalRow]:
    """Return the fitted A and E rows, excluding formal ``*`` rejections."""
    return [row for row in parse_linecal_rows(text) if row.fit_status == "fitted"]


def read_observed_lines(path: str | Path) -> list[XuLinecalRow]:
    return parse_observed_lines(Path(path).read_text(encoding="utf-8", errors="replace"))


def write_observed_csv(path: str | Path, rows: list[XuLinecalRow]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fields = list(XuLinecalRow.__dataclass_fields__)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
