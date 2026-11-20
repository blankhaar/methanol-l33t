"""Prediction-only access to the frozen Xu/BELGI methanol Hamiltonian.

The historical ``globalfit`` program normally reads observations and performs
a least-squares fit.  For molecular-data production we instead load the final
119 named coefficients, set every other BELGI coefficient explicitly to zero,
and diagonalize the Hamiltonian without reading or fitting a single line.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import io
import json
import math
import platform
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from .quantum import LevelKey
from .paths import portable_path
from .xu import CM1_TO_MHZ, read_linecal_rows


DATASET_ID = "local-belgi-xu2008-frozen-hamiltonian"
XU_FINAL_CONSTANTS_SHA256 = "11353561e41ef266477e76898282d46b58f7fe92f228cd1f41061b9c95716f3d"
XU_FINAL_NONZERO_COUNT = 119
KTRONC = 10
STATE_TRACKING_LABEL = 1
BELGI_VENDORED_SOURCE_TREE_SHA256 = (
    "6ed8a782fd0dd2cbbe7ef01ca285ae700f8a50e91a0ad2d944a8ca04321a0e8e"
)

# Canonical CST order from the WRITE list in BELGI data11.f.  The older local
# const_order.txt/const_clean.txt files are incomplete and contain a displaced
# AK4BJ value, so positional data from those files must never be used.
CST_PARAMETER_ORDER = tuple(
    """
    OA B C DJ DJK DK ODELN ODELK HJ HJK
    HKJ HK OHJ OHJK OHK FV OFV GV ALV AK1
    AK2 AK3 AK4 AK5 AK6 AK7 DAB ODAB DABJ DABK
    OLV C1 C2 C3 DUMMY1 AK5J C2J GVJ AK2J C1J
    ANV BK2 AMV DAC DBC C4 ODELTA DELTA AK1J AK5K
    AK2K AK1K AK3J AK3K BK1 AK7J C11 ODAB6 AK6K AK7K
    AK66 AK6J AK3B AK4B C2K C1K C4J C4K DELTAB ODELTB
    C12 FPARA RHO V3 V6 V9 RHORHO LJ LK LJKKK
    LJJKK LJJJK OLKKK OLJKK OLJJK OLJJJ HABK HABJ HABJK GJ6
    GJ4K2 GJ2K4 GK6 FK6 AK6KK AK6JK AK6JJ AK7KK AK7JK AK7JJ
    AK3BK AMUA0 AMUA3 AMUB0 AMUB3 AMUC3 AGJ6 AGJ4K2 AGJ2K4 AGK6
    FJ6 FJ4K2 FJ2K4 AK7BC V9J V9K ODABK C6J4 C6J2K2 C6K4
    G4J4 G4J2K2 G4K4 GAJ4 GAJ2K2 GAK4 AK3BJ AK4BJ AK4BK ODABJ
    ODAC ODBC DELTAJ DELTAK DAGJ DAGK CABK CABJ CBCK CBCJ
    DG4J DG4K DACK DBCJ DACJ V12 V12J AK4C AK3C DAC6
    DBC6 GB AK8 AK9 AK8J AK8K AK9J AK9K DELTJK DELTKK
    ODAB9 C9BC C6BCJ C6BCK AG3BCJ AG3BCK G4BCJ G4BCK AK4BC AK3BC
    AK7AB BC4 CORREC AK76 OBC4 AK3BK6 AK4BK4 AK3CK4 AK4D AK3D
    AK4CK AK3CK AK3BK4 AK4BK6 AK4CK4 DACJJ DACJK DACKK DBCK DBCJK
    ODABKK C2JK C1JK C4JK C6BCJK G3BCJK G4BCJK C2KK C1KK C4KK
    C6BCKK G3BCKK G4BCKK DG3J DG3K
    """.split()
)
if len(CST_PARAMETER_ORDER) != 205:  # pragma: no cover - source-code invariant
    raise RuntimeError("BELGI CST parameter order must contain exactly 205 entries")

DIPOLE_PARAMETERS = frozenset({"AMUA0", "AMUA3", "AMUB0", "AMUB3", "AMUC3"})

# Source files needed by the prediction-only driver.  Include files are found
# through -I and are checksum-bound separately in the run manifest.
BELGI_IMPLEMENTATION_FILES = (
    "setup.f", "rottor.f", "htorse.f", "aset.f", "hrtset.f",
    "h.f", "encal.f", "vset.f", "vrtdi.f", "resum.f", "matout.f",
    "iposa.f", "ipose.f", "parity.f", "irep.f", "encal2.f", "vtors.f",
    "rule.f", "printener.f", "printderivl.f", "expectation.f", "factt.f",
    "ad3j.f", "inten.f", "pythag.f", "blockdata.f",
)

_SETUP_BOUNDS_BUG = "      max = 2*nn2+1"
_SETUP_BOUNDS_FIX = "      max = kdmx"

_ASSIGNMENT = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9]*)\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?)\s*,?\s*$"
)


class HamiltonianError(RuntimeError):
    """Raised when the frozen-Hamiltonian build or validation fails."""


@dataclass(frozen=True, slots=True)
class FrozenLevel:
    key: LevelKey
    energy_raw_cm1: float
    energy_cm1: float


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_final_constants(path: str | Path, *, require_xu_final: bool = True) -> dict[str, str]:
    """Parse named Fortran assignments without relying on positional files.

    Returned values retain their decimal spelling.  The Fortran namelist read
    therefore performs the same conversion to double precision as the
    historical calculation.
    """
    source = Path(path)
    assignments: dict[str, str] = {}
    for number, line in enumerate(source.read_text(encoding="ascii").splitlines(), 1):
        if not line.strip():
            continue
        match = _ASSIGNMENT.fullmatch(line)
        if not match:
            raise HamiltonianError(f"invalid constants assignment at {source}:{number}: {line!r}")
        name, literal = match.groups()
        name = name.upper()
        if name not in CST_PARAMETER_ORDER:
            raise HamiltonianError(f"unknown BELGI constant {name!r} at {source}:{number}")
        if name in assignments:
            raise HamiltonianError(f"duplicate BELGI constant {name!r} at {source}:{number}")
        value = float(literal.replace("D", "E").replace("d", "e"))
        if not math.isfinite(value):
            raise HamiltonianError(f"non-finite BELGI constant {name!r}")
        assignments[name] = literal

    if require_xu_final:
        digest = sha256_file(source)
        if digest != XU_FINAL_CONSTANTS_SHA256:
            raise HamiltonianError(
                f"constants SHA-256 {digest} is not the audited Xu final vector "
                f"{XU_FINAL_CONSTANTS_SHA256}"
            )
        if len(assignments) != XU_FINAL_NONZERO_COUNT:
            raise HamiltonianError(
                f"Xu final vector needs {XU_FINAL_NONZERO_COUNT} named assignments; "
                f"found {len(assignments)}"
            )
    return assignments


def render_complete_cste_namelist(assignments: dict[str, str]) -> str:
    """Serialize all 205 CST slots, including explicit zeros.

    Explicit serialization also neutralizes a historical ``data11.f`` default-
    initialization typo that leaves ``AMUB3`` unset when it is omitted.
    """
    unknown = set(assignments) - set(CST_PARAMETER_ORDER)
    if unknown:
        raise HamiltonianError(f"unknown BELGI constants: {sorted(unknown)}")
    lines = ["CM-1", "&CSTE"]
    for name in CST_PARAMETER_ORDER:
        lines.append(f" {name} = {assignments.get(name, '0.0D0')},")
    lines.append("/")
    return "\n".join(lines) + "\n"


def render_complete_cst_vector(assignments: dict[str, str]) -> str:
    """Serialize the canonical 205-slot vector for the prediction driver."""
    unknown = set(assignments) - set(CST_PARAMETER_ORDER)
    if unknown:
        raise HamiltonianError(f"unknown BELGI constants: {sorted(unknown)}")
    values = [assignments.get(name, "0.0D0") for name in CST_PARAMETER_ORDER]
    return "\n".join(("BELGI_FROZEN_CST_205", *values, ""))


def _find_compiler(explicit: str | Path | None) -> Path:
    if explicit is not None:
        candidate = Path(explicit).expanduser().resolve()
        if not candidate.is_file():
            raise HamiltonianError(f"Fortran compiler not found: {candidate}")
        return candidate
    for name in ("gfortran", "gfortran-15", "gfortran-14", "gfortran-13"):
        found = shutil.which(name)
        if found:
            return Path(found).resolve()
    raise HamiltonianError("no gfortran compiler found on PATH; use --compiler PATH")


def _belgi_sources(source_dir: Path) -> tuple[list[Path], list[Path]]:
    implementations = [source_dir / name for name in BELGI_IMPLEMENTATION_FILES]
    includes = sorted(source_dir.glob("common_*.f"))
    missing = [path for path in implementations if not path.is_file()]
    if missing:
        raise HamiltonianError(f"missing BELGI source files: {[path.name for path in missing]}")
    if not includes:
        raise HamiltonianError(f"no BELGI common include files found under {source_dir}")
    return implementations, includes


def belgi_source_inventory(source_dir: Path) -> dict[str, object]:
    """Checksum every vendored Fortran source and its deterministic inventory."""
    hashes = {
        path.name: sha256_file(path)
        for path in sorted(source_dir.glob("*.f"), key=lambda path: path.name)
    }
    payload = "".join(f"{digest}  {name}\n" for name, digest in hashes.items())
    tree_sha256 = hashlib.sha256(payload.encode("ascii")).hexdigest()
    return {
        "file_count": len(hashes),
        "file_sha256": hashes,
        "tree_sha256": tree_sha256,
        "audited_bundled_tree_sha256": BELGI_VENDORED_SOURCE_TREE_SHA256,
        "matches_audited_bundled_tree": tree_sha256 == BELGI_VENDORED_SOURCE_TREE_SHA256,
    }


def _portable_compile_command(
    command: list[str], root: Path, source_dir: Path, build_dir: Path
) -> list[str]:
    """Replace checkout/build locations and compiler executable with replay tokens."""
    replacements = sorted(
        (
            (str(root), "{PROJECT_ROOT}"),
            (str(source_dir), "{BELGI_SOURCE_DIR}"),
            (str(build_dir), "{BUILD_DIR}"),
        ),
        key=lambda item: len(item[0]),
        reverse=True,
    )
    result: list[str] = []
    for index, token in enumerate(command):
        if index == 0:
            result.append("{FORTRAN_COMPILER}")
            continue
        for source, replacement in replacements:
            token = token.replace(source, replacement)
        result.append(token)
    return result


def _compile_driver(
    project_root: Path,
    source_dir: Path,
    build_dir: Path,
    compiler: Path,
    *,
    driver_filename: str = "belgi_frozen_levels.f",
    executable_name: str = "belgi_frozen_levels",
    parity_repair: bool = False,
    cached_dipoles: bool = False,
) -> tuple[Path, list[str], str, dict[str, object]]:
    implementations, _ = _belgi_sources(source_dir)
    compiled_implementations: list[Path] = []
    source_patches: dict[str, object] = {}
    for source in implementations:
        if source.name != "setup.f":
            compiled_implementations.append(source)
            continue
        setup_text = source.read_text(encoding="ascii")
        if setup_text.count(_SETUP_BOUNDS_BUG) != 1:
            raise HamiltonianError(
                "setup.f does not contain the audited VECTPR bounds defect exactly once"
            )
        compiled_setup = build_dir / "setup.f"
        compiled_setup.write_text(
            setup_text.replace(_SETUP_BOUNDS_BUG, _SETUP_BOUNDS_FIX),
            encoding="ascii",
        )
        compiled_implementations.append(compiled_setup)
        source_patches[source.name] = {
            "reason": (
                "prevent SETUP from clearing 95 entries in VECTPR's "
                "93-element KDMX dimension and corrupting ENCAL2 state tracking"
            ),
            "original_sha256": sha256_file(source),
            "compiled_sha256": sha256_file(compiled_setup),
            "replacement": f"{_SETUP_BOUNDS_BUG.strip()} -> {_SETUP_BOUNDS_FIX.strip()}",
        }
    repair_modules: list[Path] = []
    if parity_repair:
        from .parity_repair import materialize_belgi_repair_sources
        compiled_implementations, repair_modules, repair_patches = materialize_belgi_repair_sources(
            project_root, compiled_implementations, build_dir
        )
        source_patches.update(repair_patches)
    driver = project_root / "fortran" / driver_filename
    if not driver.is_file():
        raise HamiltonianError(f"prediction driver is missing: {driver}")
    if cached_dipoles:
        if not parity_repair:
            raise HamiltonianError("cached full-vector dipoles require explicit parity_repair=True")
        from .parity_repair import materialize_cached_strength_driver
        driver, cache_sources, cache_patches = materialize_cached_strength_driver(project_root, driver, build_dir)
        repair_modules.extend(cache_sources)
        source_patches.update(cache_patches)
    executable = build_dir / executable_name
    if sys.platform == "darwin":
        lapack_link_flags = ["-Wl,-framework", "-Wl,Accelerate"]
    else:
        lapack_link_flags = ["-llapack", "-lblas"]
    command = [
        str(compiler), "-std=legacy", "-ffixed-line-length-none", "-ffree-line-length-none",
        "-fallow-argument-mismatch", "-O2", "-I", str(source_dir),
        *(str(path) for path in repair_modules), str(driver), *(str(path) for path in compiled_implementations),
        *lapack_link_flags, "-o", str(executable),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False, cwd=build_dir)
    if completed.returncode != 0:
        raise HamiltonianError(
            "BELGI driver compilation failed\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    version = subprocess.run(
        [str(compiler), "--version"], text=True, capture_output=True, check=True
    ).stdout.splitlines()[0]
    return executable, command, version, source_patches


def _parse_level_stdout(stdout: str) -> list[FrozenLevel]:
    begin = "BEGIN_FROZEN_LEVELS_CSV\n"
    end = "END_FROZEN_LEVELS_CSV"
    if begin not in stdout or end not in stdout:
        raise HamiltonianError("BELGI output does not contain the frozen-level CSV markers")
    table = stdout.split(begin, 1)[1].split(end, 1)[0]
    levels: list[FrozenLevel] = []
    seen: set[LevelKey] = set()
    for row in csv.DictReader(io.StringIO(table)):
        key = LevelKey(
            row["species"], int(row["vt"]), int(row["j"]), int(row["k"]),
            row["parity"].strip(),
        )
        if key in seen:
            raise HamiltonianError(f"duplicate BELGI level identity {key.encode()}")
        seen.add(key)
        raw = float(row["energy_raw_cm1"])
        relative = float(row["energy_species_cm1"])
        if not math.isfinite(raw) or not math.isfinite(relative):
            raise HamiltonianError(f"non-finite energy for {key.encode()}")
        levels.append(FrozenLevel(key, raw, relative))
    return levels


def _expected_level_count(j_max: int, abs_k_max: int, vt_max: int) -> int:
    per_torsional_state = sum(1 + 2 * min(j, abs_k_max) for j in range(j_max + 1))
    return 2 * (vt_max + 1) * per_torsional_state


def validate_against_xu_calculations(
    levels: list[FrozenLevel], xu_output: str | Path, *, project_root: str | Path | None = None
) -> dict[str, object]:
    """Check energy differences against every Xu CALCULATION row in scope.

    Xu prints microwave calculations to 0.001 MHz and infrared calculations to
    0.00001 cm-1.  Half of the last printed digit is therefore the strongest
    meaningful reproduction tolerance.
    """
    energies = {level.key: level.energy_raw_cm1 for level in levels}
    matched: list[dict[str, object]] = []
    direct_identity_failures: list[dict[str, object]] = []
    rowwise_spectrum_equivalences: list[dict[str, object]] = []
    spectrum_failures: list[dict[str, object]] = []
    mapping_targets: dict[LevelKey, set[LevelKey]] = {}
    tolerances = {"MHz": 0.0005001, "cm-1": 0.000005001}
    for row in read_linecal_rows(xu_output):
        if row.upper_key not in energies or row.lower_key not in energies:
            continue
        predicted_cm1 = energies[row.upper_key] - energies[row.lower_key]
        predicted = predicted_cm1 * CM1_TO_MHZ if row.unit == "MHz" else predicted_cm1
        direct_error = predicted - row.calculation
        spectrum_error = direct_error
        match_mode = "direct_quantum_labels"
        continuity_upper_key = row.upper_key
        continuity_lower_key = row.lower_key
        direct_ok = abs(direct_error) <= tolerances[row.unit]
        if not direct_ok:
            direct_identity_failures.append(
                {
                    "source_line": row.source_line,
                    "species": row.species,
                    "unit": row.unit,
                    "upper_key": row.upper_key.encode(),
                    "lower_key": row.lower_key.encode(),
                    "predicted": predicted,
                    "xu_calculation": row.calculation,
                    "direct_error": direct_error,
                    "tolerance": tolerances[row.unit],
                }
            )

        # ENCAL2 follows eigenvectors continuously through near-degenerate A
        # parity pairs.  Xu's legacy LINECAL labels use the older local IPOSA
        # assignment at some of those crossings.  The Hamiltonian spectrum is
        # identical as an unordered local spectrum.  A rowwise parity exchange
        # is useful evidence for that spectrum, but it is not accepted as an
        # exact quantum-identity regression: the substitutions must form one
        # global bijection before they can be used as a label mapping.
        equivalent: tuple[float, LevelKey, LevelKey, float] | None = None
        if not direct_ok and row.species == "A":
            upper_variants = [row.upper_key]
            lower_variants = [row.lower_key]
            if row.upper_key.k > 0:
                upper_variants.append(
                    LevelKey(
                        "A", row.upper_key.vt, row.upper_key.j, row.upper_key.k,
                        "-" if row.upper_key.parity == "+" else "+",
                    )
                )
            if row.lower_key.k > 0:
                lower_variants.append(
                    LevelKey(
                        "A", row.lower_key.vt, row.lower_key.j, row.lower_key.k,
                        "-" if row.lower_key.parity == "+" else "+",
                    )
                )
            for upper_key in upper_variants:
                for lower_key in lower_variants:
                    if upper_key not in energies or lower_key not in energies:
                        continue
                    alternative_cm1 = energies[upper_key] - energies[lower_key]
                    alternative = (
                        alternative_cm1 * CM1_TO_MHZ
                        if row.unit == "MHz" else alternative_cm1
                    )
                    alternative_error = alternative - row.calculation
                    candidate = (
                        abs(alternative_error), upper_key, lower_key, alternative_error
                    )
                    if equivalent is None or candidate[0] < equivalent[0]:
                        equivalent = candidate
            if equivalent is not None and equivalent[0] <= tolerances[row.unit]:
                spectrum_error = equivalent[3]
                continuity_upper_key = equivalent[1]
                continuity_lower_key = equivalent[2]
                match_mode = "rowwise_a_parity_spectrum_equivalent"
                rowwise_spectrum_equivalences.append(
                    {
                        "source_line": row.source_line,
                        "unit": row.unit,
                        "xu_upper_key": row.upper_key.encode(),
                        "xu_lower_key": row.lower_key.encode(),
                        "continuity_upper_key": continuity_upper_key.encode(),
                        "continuity_lower_key": continuity_lower_key.encode(),
                        "direct_error": direct_error,
                        "spectrum_equivalent_error": spectrum_error,
                    }
                )

        matched.append(
            {
                "species": row.species,
                "unit": row.unit,
                "direct_error": direct_error,
                "spectrum_error": spectrum_error,
                "match_mode": match_mode,
            }
        )
        for source_key, target_key in (
            (row.upper_key, continuity_upper_key),
            (row.lower_key, continuity_lower_key),
        ):
            mapping_targets.setdefault(source_key, set()).add(target_key)
        if abs(spectrum_error) > tolerances[row.unit]:
            spectrum_failures.append(
                {
                    "source_line": row.source_line,
                    "species": row.species,
                    "unit": row.unit,
                    "upper_key": row.upper_key.encode(),
                    "lower_key": row.lower_key.encode(),
                    "predicted": predicted,
                    "xu_calculation": row.calculation,
                    "direct_error": direct_error,
                    "best_rowwise_spectrum_error": spectrum_error,
                    "tolerance": tolerances[row.unit],
                }
            )
    if not matched:
        raise HamiltonianError("no Xu CALCULATION rows matched the generated levels")
    if spectrum_failures:
        raise HamiltonianError(
            "frozen Hamiltonian failed "
            f"{len(spectrum_failures)}/{len(matched)} Xu spectrum checks; "
            f"first failures: {spectrum_failures[:5]}"
        )
    groups: dict[str, dict[str, float | int]] = {}
    for species in ("A", "E"):
        for unit in ("MHz", "cm-1"):
            selected = [row for row in matched if (row["species"], row["unit"]) == (species, unit)]
            if selected:
                direct_errors = [float(row["direct_error"]) for row in selected]
                spectrum_errors = [float(row["spectrum_error"]) for row in selected]
                groups[f"{species}_{unit}"] = {
                    "rows": len(selected),
                    "direct_quantum_identity_rows": sum(
                        abs(error) <= tolerances[unit] for error in direct_errors
                    ),
                    "direct_quantum_identity_failure_rows": sum(
                        abs(error) > tolerances[unit] for error in direct_errors
                    ),
                    "rowwise_a_parity_spectrum_equivalent_rows": sum(
                        row["match_mode"] == "rowwise_a_parity_spectrum_equivalent"
                        for row in selected
                    ),
                    "max_abs_direct_identity_error": max(abs(error) for error in direct_errors),
                    "max_abs_spectrum_error": max(abs(error) for error in spectrum_errors),
                    "spectrum_error_rms": math.sqrt(
                        sum(error * error for error in spectrum_errors) / len(spectrum_errors)
                    ),
                    "tolerance": tolerances[unit],
                }
    ambiguous_keys = sorted(
        key.encode() for key, targets in mapping_targets.items() if len(targets) > 1
    )
    return {
        "energy_spectrum_passed": True,
        "direct_quantum_identity_passed": not direct_identity_failures,
        "xu_output": (
            portable_path(project_root, xu_output)
            if project_root is not None
            else str(Path(xu_output).expanduser().resolve())
        ),
        "xu_sha256": sha256_file(xu_output),
        "matched_rows": len(matched),
        "direct_quantum_identity_rows": len(matched) - len(direct_identity_failures),
        "direct_quantum_identity_failure_rows": len(direct_identity_failures),
        "direct_quantum_identity_failures": direct_identity_failures,
        "rowwise_a_parity_spectrum_equivalent_rows": len(rowwise_spectrum_equivalences),
        "rowwise_a_parity_spectrum_equivalences": rowwise_spectrum_equivalences,
        "global_bijective_parity_relabeling_established": False,
        "ambiguous_source_level_count": len(ambiguous_keys),
        "ambiguous_source_levels": ambiguous_keys,
        "interpretation": (
            "The energy spectrum reproduces every Xu calculation at print precision. "
            "The rowwise A parity substitutions are not a globally consistent level "
            "mapping and therefore do not validate those 142 quantum identities."
        ),
        "groups": groups,
    }


def build_frozen_levels(
    project_root: str | Path,
    constants_path: str | Path,
    belgi_source_dir: str | Path,
    output_dir: str | Path,
    *,
    j_max: int = 22,
    abs_k_max: int = 14,
    vt_max: int = 2,
    compiler: str | Path | None = None,
    require_xu_final: bool = True,
    xu_output: str | Path | None = None,
    dataset_id: str = DATASET_ID,
    parity_repair: bool = False,
) -> dict[str, object]:
    """Compile the prediction driver, diagonalize, and write levels + manifest."""
    if parity_repair and dataset_id == DATASET_ID:
        dataset_id += "-parity-full9-candidate-v1"
    root = Path(project_root).resolve()
    constants = Path(constants_path).expanduser().resolve()
    source_dir = Path(belgi_source_dir).expanduser().resolve()
    target = Path(output_dir).expanduser().resolve()
    if not (0 <= vt_max <= 4 and 0 <= abs_k_max <= j_max <= 46 and j_max >= 1):
        raise HamiltonianError("require 0 <= vt_max <= 4, 0 <= abs_k_max <= j_max, and 1 <= j_max <= 46")
    assignments = parse_final_constants(constants, require_xu_final=require_xu_final)
    cst_vector = render_complete_cst_vector(assignments)
    compiler_path = _find_compiler(compiler)
    implementations, includes = _belgi_sources(source_dir)

    with tempfile.TemporaryDirectory(prefix="belgi-frozen-hamiltonian-") as temporary:
        build_dir = Path(temporary)
        executable, command, compiler_version, source_patches = _compile_driver(
            root, source_dir, build_dir, compiler_path, parity_repair=parity_repair
        )
        compile_command_template = _portable_compile_command(
            command, root, source_dir, build_dir
        )
        completed = subprocess.run(
            [str(executable), str(j_max), str(abs_k_max), str(vt_max)],
            input=cst_vector,
            text=True,
            capture_output=True,
            cwd=build_dir,
            check=False,
        )
        if completed.returncode != 0:
            raise HamiltonianError(
                f"BELGI frozen-Hamiltonian run failed ({completed.returncode})\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        levels = _parse_level_stdout(completed.stdout)

    expected = _expected_level_count(j_max, abs_k_max, vt_max)
    if len(levels) != expected:
        raise HamiltonianError(f"expected {expected} levels, received {len(levels)}")
    xu_validation = (
        validate_against_xu_calculations(levels, xu_output, project_root=root)
        if xu_output is not None else None
    )

    target.mkdir(parents=True, exist_ok=True)
    levels_path = target / "levels.csv"
    fields = [
        "species", "level_key", "vt", "j", "k", "parity", "energy_cm1",
        "energy_raw_cm1", "weight", "source", "source_index",
        "mapping_confidence", "constants_sha256",
    ]
    with levels_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for index, level in enumerate(levels, 1):
            writer.writerow(
                {
                    "species": level.key.species,
                    "level_key": level.key.encode(),
                    "vt": level.key.vt,
                    "j": level.key.j,
                    "k": level.key.k,
                    "parity": level.key.parity,
                    "energy_cm1": f"{level.energy_cm1:.17g}",
                    "energy_raw_cm1": f"{level.energy_raw_cm1:.17g}",
                    "weight": 2 * level.key.j + 1,
                    "source": dataset_id,
                    "source_index": index,
                    "mapping_confidence": (
                        "global_parity_overlap_candidate" if parity_repair and level.key.species == "A"
                        else "legacy_ENCAL2_full9" if parity_repair
                        else "exact_internal_quantum_identity"
                    ),
                    "constants_sha256": sha256_file(constants),
                }
            )

    source_inventory = belgi_source_inventory(source_dir)
    driver_path = root / "fortran" / "belgi_frozen_levels.f"
    manifest: dict[str, object] = {
        "schema_version": 2,
        "dataset_id": dataset_id,
        "method": "prediction-only diagonalization; no observations read and no fit performed",
        "constants_loader": "named Python validation followed by direct canonical 205-slot vector input",
        "constants": {
            "path": portable_path(root, constants),
            "sha256": sha256_file(constants),
            "units": (
                "Hamiltonian coefficients are cm-1 except dimensionless RHO/RHORHO; "
                "AMU slots are Debye and are zero in this energy-only vector"
            ),
            "named_assignments": len(assignments),
            "explicit_zero_slots": len(CST_PARAMETER_ORDER) - len(assignments),
            "total_cst_slots": len(CST_PARAMETER_ORDER),
            "dipole_parameters": {name: assignments.get(name, "0.0D0") for name in sorted(DIPOLE_PARAMETERS)},
        },
        "basis": {
            "ktronc": KTRONC,
            "state_tracking_label": STATE_TRACKING_LABEL,
            "state_tracking": (
                "A: full9 gauge-invariant global parity-constrained overlap assignment; E: ENCAL2 identity with full9 export"
                if parity_repair else "ENCAL2 eigenvector continuity"
            ),
            "parity_repair_enabled": parity_repair,
            "torsional_states": list(range(vt_max + 1)),
        },
        "scope": {
            "j_max": j_max,
            "abs_k_max": abs_k_max,
            "vt_max": vt_max,
            "level_count": len(levels),
            "a_level_count": sum(level.key.species == "A" for level in levels),
            "e_level_count": sum(level.key.species == "E" for level in levels),
            "energy_zero": {
                "A": {
                    "level_key": "A:vt0:J0:K0:+",
                    "raw_cm1": next(
                        level.energy_raw_cm1 - level.energy_cm1
                        for level in levels if level.key.species == "A"
                    ),
                },
                "E": {
                    "level_key": "E:vt0:J1:K-1",
                    "raw_cm1": next(
                        level.energy_raw_cm1 - level.energy_cm1
                        for level in levels if level.key.species == "E"
                    ),
                },
            },
        },
        "program": {
            "driver": portable_path(root, driver_path),
            "driver_sha256": sha256_file(driver_path),
            "belgi_source_dir": portable_path(root, source_dir),
            "belgi_source_inventory": source_inventory,
            "compiled_belgi_source_sha256": {
                path.name: sha256_file(path)
                for path in sorted((*implementations, *includes), key=lambda path: path.name)
            },
            "compiler": compiler_path.name,
            "compiler_version": compiler_version,
            "compile_command_template": compile_command_template,
            "runtime_environment": {
                "platform": f"{platform.system()} {platform.release()}",
                "machine": platform.machine(),
                "lapack_linkage": (
                    "Apple Accelerate" if sys.platform == "darwin" else "system LAPACK/BLAS"
                ),
            },
            "source_patches": source_patches,
        },
        "outputs": {
            "levels_csv": portable_path(root, levels_path),
            "levels_sha256": sha256_file(levels_path),
        },
    }
    if xu_validation is not None:
        manifest["validation"] = {"xu_calculations": xu_validation}
    manifest_path = target / "manifest.json"
    manifest["outputs"]["manifest_json"] = portable_path(root, manifest_path)  # type: ignore[index]
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
