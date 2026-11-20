"""Parsers for fixed-width JPL/CDMS catalog records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True, slots=True)
class CatalogLine:
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
    raw: str


def _catalog_qns(field: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    chunks = tuple(field[pos : pos + 2].strip() for pos in range(0, len(field), 2))
    chunks = tuple(value for value in chunks if value)
    if len(chunks) % 2:
        raise ValueError("odd number of quantum-number fields")
    half = len(chunks) // 2
    return chunks[:half], chunks[half:]


def parse_catalog_line(line: str) -> CatalogLine:
    """Parse one standard 80+ column JPL/CDMS catalog line.

    Special CDMS ``Smu2`` catalogs share much of this layout but do not have an
    ordinary log-intensity third field; callers must tag their semantics before
    doing any intensity conversion.
    """
    if len(line.rstrip("\n")) < 55:
        raise ValueError("catalog line shorter than fixed-width header")
    raw = line.rstrip("\n")
    upper_qns, lower_qns = _catalog_qns(raw[55:])
    return CatalogLine(
        frequency_mhz=float(raw[0:13]),
        uncertainty_mhz=float(raw[13:21]),
        log10_intensity=float(raw[21:29]),
        degrees_of_freedom=int(raw[29:31]),
        lower_energy_cm1=float(raw[31:41]),
        upper_degeneracy=int(raw[41:44]),
        tag=int(raw[44:51]),
        qn_format=int(raw[51:55]),
        upper_qns=upper_qns,
        lower_qns=lower_qns,
        raw=raw,
    )


def parse_catalog(lines: Iterable[str]) -> list[CatalogLine]:
    parsed: list[CatalogLine] = []
    for number, line in enumerate(lines, start=1):
        if not line.strip() or line.lstrip().startswith(("#", "!")):
            continue
        try:
            parsed.append(parse_catalog_line(line))
        except Exception as exc:
            raise ValueError(f"catalog line {number}: {exc}") from exc
    return parsed


def read_catalog(path: str | Path) -> list[CatalogLine]:
    return parse_catalog(Path(path).read_text(encoding="ascii").splitlines())
