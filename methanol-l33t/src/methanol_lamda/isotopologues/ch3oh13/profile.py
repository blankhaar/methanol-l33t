"""Immutable setup for the :sup:`13`CH3OH BELGI/LAMDA pipeline.

The published Xu et al. (2014) Table-2 constants describe the Hamiltonian
only.  BELGI's five dipole slots therefore stay out of the Hamiltonian
assignment set and are overlaid only for a line-strength calculation.  This
separation prevents a radiative-model choice from silently becoming part of
the fitted spectroscopic parameter vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from ...hamiltonian import CST_PARAMETER_ORDER, parse_final_constants
from .collisions import MASS_13CH3OH_U


ISOTOPOLOGUE = "13CH3OH"
MOLECULAR_MASS_U = MASS_13CH3OH_U
SPECIES = ("A", "E")

HAMILTONIAN_DATASET_ID = "local-belgi-xu2014-13ch3oh-frozen-hamiltonian"
RADIATIVE_DATASET_ID = "local-belgi-xu2014-13ch3oh-constant-dipole"
DEFAULT_CONSTANTS_RELATIVE_PATH = Path(
    "isotopologues/13ch3oh/constants/table2-2014-constants.txt"
)

SPECIES_NAMES: Mapping[str, str] = MappingProxyType(
    {
        "A": "A-13CH3OH",
        "E": "E-13CH3OH",
    }
)
LAMDA_FILENAMES: Mapping[str, str] = MappingProxyType(
    {
        "A": "a-13ch3oh-vt012-reduced-mass-proxy.dat",
        "E": "e-13ch3oh-vt012-reduced-mass-proxy.dat",
    }
)

# Permanent IAM dipoles quoted by Xu et al. (2014), in Debye.  The three
# torsion-dependent cos(3 gamma) terms are deliberately zero.
DIPOLES_DEBYE: Mapping[str, float] = MappingProxyType(
    {
        "mu_a": 0.899,
        "mu_a_cos_3gamma": 0.0,
        "mu_b": -1.44,
        "mu_b_cos_3gamma": 0.0,
        "mu_c_cos_3gamma": 0.0,
    }
)
DIPOLE_ASSIGNMENTS: Mapping[str, str] = MappingProxyType(
    {
        "AMUA0": "0.899D0",
        "AMUA3": "0.0D0",
        "AMUB0": "-1.44D0",
        "AMUB3": "0.0D0",
        "AMUC3": "0.0D0",
    }
)

# Independent transcription of the names having a nonblank 13CH3OH value in
# Xu et al. (2014), Table 2.  Validating the exact set, rather than only the
# count, catches accidental ingestion of the adjacent parent-isotopologue
# column.  The order here is the paper order and has no vector semantics.
TABLE2_PARAMETER_NAMES = frozenset(
    """
    FPARA V3 RHORHO OA B C DAB AK4 V6 AK3 GV AK2 DELTA C1 FV AK5
    ODAB C2 DAC DBC ALV AK1 ODELTA C4 DJ DJK DK ODELN ODELK DABJ
    DABK AK4B V9 AK3B AMV BK1 DELTAB C3 AK7J ANV BK2 ODAB6 C11
    DAC6 AK3J AK3K ODELTB C12 AK6J AK6K GVJ AK2J DELTAJ C1J AK2K
    DELTAK C1K OFV AK5J ODABJ C2J AK5K ODABK C2K DACJ DBCJ DACK
    DBCK OLV AK1J DAGJ C4J AK1K C4K HJ HJK HKJ HK OHJK OHK AK4C
    AK3C AK4BJ AK4BK V9J ODAB9 AK3BJ AK3BK G4J2K2 G4K4 C6J2K2
    CABJ C6BCJ C6K4 C6BCK GAJ2K2 GAK4 AK6JK GK6 FJ2K4 DACJK
    AGJ2K4 AGK6
    """.split()
)

_CST_PARAMETER_NAMES = frozenset(CST_PARAMETER_ORDER)
TABLE2_ZERO_PARAMETER_NAMES = _CST_PARAMETER_NAMES - TABLE2_PARAMETER_NAMES
if len(TABLE2_PARAMETER_NAMES) != 103:  # pragma: no cover - source invariant
    raise RuntimeError("13CH3OH Table-2 profile must name exactly 103 parameters")
if not TABLE2_PARAMETER_NAMES <= _CST_PARAMETER_NAMES:  # pragma: no cover
    raise RuntimeError("13CH3OH Table-2 profile contains a non-BELGI parameter")
if len(TABLE2_ZERO_PARAMETER_NAMES) != 102:  # pragma: no cover
    raise RuntimeError("13CH3OH Hamiltonian must leave exactly 102 CST slots at zero")


class CH3OH13ProfileError(ValueError):
    """Raised when the isotope setup is inconsistent with the audited profile."""


@dataclass(frozen=True, slots=True)
class CH3OH13AssignmentSets:
    """Separate immutable Hamiltonian and radiative-operator assignments."""

    constants_path: Path
    hamiltonian_assignments: Mapping[str, str]
    strength_assignments: Mapping[str, str]
    hamiltonian_zero_parameter_names: frozenset[str]


def constants_path(project_root: str | Path) -> Path:
    """Return the repository-local Xu et al. (2014) constants artifact."""
    return Path(project_root).resolve() / DEFAULT_CONSTANTS_RELATIVE_PATH


def load_assignment_sets(
    project_root: str | Path,
    *,
    source: str | Path | None = None,
) -> CH3OH13AssignmentSets:
    """Load and audit the Table-2 Hamiltonian, then overlay strength dipoles.

    ``hamiltonian_assignments`` contains exactly the 103 published Table-2
    parameters.  The complete 205-slot Hamiltonian vector obtains its other
    102 values from BELGI's explicit-zero serializer.  ``strength_assignments``
    adds the five dipole slots, including explicit zeros for all three
    cos(3 gamma) components.
    """
    if source is None:
        path = constants_path(project_root)
    else:
        candidate = Path(source).expanduser()
        path = (
            candidate.resolve()
            if candidate.is_absolute()
            else (Path(project_root).resolve() / candidate).resolve()
        )
    assignments = parse_final_constants(path, require_xu_final=False)
    actual_names = frozenset(assignments)
    if len(assignments) != 103 or actual_names != TABLE2_PARAMETER_NAMES:
        missing = sorted(TABLE2_PARAMETER_NAMES - actual_names)
        extra = sorted(actual_names - TABLE2_PARAMETER_NAMES)
        raise CH3OH13ProfileError(
            "13CH3OH Table-2 constants must contain the audited exact 103-name "
            f"set; found {len(assignments)} assignments, missing={missing}, extra={extra}"
        )

    zero_names = _CST_PARAMETER_NAMES - actual_names
    if len(zero_names) != 102:
        raise CH3OH13ProfileError(
            f"13CH3OH Hamiltonian complement must contain 102 slots; found {len(zero_names)}"
        )
    if actual_names & frozenset(DIPOLE_ASSIGNMENTS):
        raise CH3OH13ProfileError(
            "Table-2 Hamiltonian constants must not contain BELGI dipole slots"
        )

    hamiltonian = MappingProxyType(dict(assignments))
    strength_values = dict(assignments)
    strength_values.update(DIPOLE_ASSIGNMENTS)
    strength = MappingProxyType(strength_values)
    return CH3OH13AssignmentSets(
        constants_path=path,
        hamiltonian_assignments=hamiltonian,
        strength_assignments=strength,
        hamiltonian_zero_parameter_names=frozenset(zero_names),
    )


__all__ = [
    "CH3OH13AssignmentSets",
    "CH3OH13ProfileError",
    "DEFAULT_CONSTANTS_RELATIVE_PATH",
    "DIPOLE_ASSIGNMENTS",
    "DIPOLES_DEBYE",
    "HAMILTONIAN_DATASET_ID",
    "ISOTOPOLOGUE",
    "LAMDA_FILENAMES",
    "MOLECULAR_MASS_U",
    "RADIATIVE_DATASET_ID",
    "SPECIES",
    "SPECIES_NAMES",
    "TABLE2_PARAMETER_NAMES",
    "TABLE2_ZERO_PARAMETER_NAMES",
    "constants_path",
    "load_assignment_sets",
]
