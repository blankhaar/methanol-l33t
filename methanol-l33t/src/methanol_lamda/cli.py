"""Command-line workflow for source capture, normalization, build, and QA."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .build import MissingSourcesError, build_products
from .ch3_rock import build_ch3_rock_products
from .co_stretch import build_co_stretch_products
from .fetch import FetchError, fetch_sources
from .hamiltonian import HamiltonianError, build_frozen_levels
from .extended import build_extended_products
from .he_build import build_he_products
from .normalize import normalize_workspace
from .paths import resolve_from_root
from .report import generate_reports, write_blocker
from .strengths import DEFAULT_SMU2_MIN_DEBYE2, build_frozen_strengths
from .validate import validate_workspace
from .vibrational_revision import build_vibrational_revision


DEFAULT_CONSTANTS = Path("data/raw/local-cache/xu2008-final-constants.txt")
DEFAULT_BELGI_SOURCE = Path("vendor/belgi")
DEFAULT_XU_OUTPUT = Path("data/raw/local-cache/xu2008-globalfit-output.txt")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="methanol-lamda")
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="project root (default: cwd)")
    sub = parser.add_subparsers(dest="command", required=True)
    fetch = sub.add_parser("fetch", help="capture configured primary sources immutably")
    fetch.add_argument("--offline", action="store_true")
    sub.add_parser("normalize", help="normalize captured sources/reference files")
    build = sub.add_parser("build", help="build six A/E data products")
    build.add_argument("--accept-incomplete-radiative", action="store_true")
    sub.add_parser(
        "build-he",
        help="build the exact A/E 150-state RF2011 BASECOL He-only products",
    )
    validate = sub.add_parser("validate", help="run dual-parser and science validation")
    validate.add_argument("paths", nargs="*")
    sub.add_parser("report", help="generate comparison and coverage reports")
    hamiltonian = sub.add_parser(
        "hamiltonian", help="diagonalize the frozen Xu/BELGI Hamiltonian without fitting"
    )
    hamiltonian.add_argument("--constants", type=Path, default=DEFAULT_CONSTANTS)
    hamiltonian.add_argument("--belgi-source", type=Path, default=DEFAULT_BELGI_SOURCE)
    hamiltonian.add_argument(
        "--output-dir", type=Path,
        default=Path("data/normalized/belgi-xu2008-frozen"),
    )
    hamiltonian.add_argument("--j-max", type=int, default=22)
    hamiltonian.add_argument("--abs-k-max", type=int, default=14)
    hamiltonian.add_argument("--vt-max", type=int, default=2)
    hamiltonian.add_argument("--compiler", type=Path)
    hamiltonian.add_argument(
        "--xu-output", type=Path, default=DEFAULT_XU_OUTPUT,
        help="Xu GLOBALFIT listing for exact CALCULATION-column regression",
    )
    strengths = sub.add_parser(
        "strengths",
        help="enumerate the Xu2008 constant-dipole line list from frozen BELGI states",
    )
    strengths.add_argument("--constants", type=Path, default=DEFAULT_CONSTANTS)
    strengths.add_argument("--belgi-source", type=Path, default=DEFAULT_BELGI_SOURCE)
    strengths.add_argument(
        "--output-dir", type=Path,
        default=Path("data/normalized/belgi-xu2008-constant-dipole"),
    )
    sub.add_parser(
        "build-extended",
        help="build exploratory A/E vt=0..4 files with vt0 intraband collision copies",
    )
    co_stretch = sub.add_parser(
        "build-co-stretch",
        help="build main- and carbon-13 vt=0..4 + vCO=1 exploratory files",
    )
    co_stretch.add_argument("--compiler", type=Path)
    sub.add_parser(
        "build-ch3-rock",
        help="add uncoupled source-bound vrock=1 manifolds to the CO-stretch files",
    )
    revision = sub.add_parser(
        "build-vibrational-revision",
        help="build source-corrected 2026-10-16 candidates without replacing originals",
    )
    revision.add_argument("--output-dir", type=Path)
    revision.add_argument("--compiler", type=Path)
    strengths.add_argument("--j-max", type=int, default=22)
    strengths.add_argument("--abs-k-max", type=int, default=14)
    strengths.add_argument("--vt-max", type=int, default=2)
    strengths.add_argument(
        "--smu2-min", type=float, default=DEFAULT_SMU2_MIN_DEBYE2,
        help="minimum retained S mu^2 in Debye^2",
    )
    strengths.add_argument("--compiler", type=Path)
    strengths.add_argument(
        "--xu-output", type=Path, default=DEFAULT_XU_OUTPUT,
        help="Xu GLOBALFIT listing to retain and identity-map observed frequencies",
    )
    strengths.add_argument(
        "--reference-levels", type=Path,
        help="separately generated levels.csv used to verify the strength-driver energies",
    )
    all_command = sub.add_parser("all", help="run fetch, normalize, build, validate, report")
    all_command.add_argument("--offline", action="store_true")
    all_command.add_argument(
        "--allow-incomplete", action="store_true",
        help="finish audit/report steps and return success when essential inputs are absent",
    )
    return parser


def _print(value: object) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    root: Path = args.root.resolve()
    try:
        if args.command == "fetch":
            _print(fetch_sources(root, offline=args.offline))
        elif args.command == "normalize":
            _print(normalize_workspace(root))
        elif args.command == "build":
            manifest = build_products(root, accept_incomplete_radiative=args.accept_incomplete_radiative)
            paths = [root / "outputs" / row["file"] for row in manifest["files"]]
            validation = validate_workspace(root, paths)
            status = generate_reports(root)
            _print({"manifest": manifest, "validation": validation, "status": status})
            return 0 if validation["production_usable"] else 1
        elif args.command == "build-he":
            manifest = build_he_products(root)
            paths = [root / "outputs" / row["file"] for row in manifest["files"]]
            validation = validate_workspace(root, paths)
            status = generate_reports(root)
            _print({"manifest": manifest, "validation": validation, "status": status})
            return 0 if validation["passed"] else 1
        elif args.command == "validate":
            report = validate_workspace(root, args.paths or None)
            _print(report)
            return 0 if report["passed"] else 1
        elif args.command == "report":
            _print(generate_reports(root))
        elif args.command == "hamiltonian":
            output_dir = resolve_from_root(root, args.output_dir)
            _print(
                build_frozen_levels(
                    root,
                    resolve_from_root(root, args.constants),
                    resolve_from_root(root, args.belgi_source),
                    output_dir,
                    j_max=args.j_max,
                    abs_k_max=args.abs_k_max,
                    vt_max=args.vt_max,
                    compiler=(
                        resolve_from_root(root, args.compiler)
                        if args.compiler is not None else None
                    ),
                    require_xu_final=True,
                    xu_output=resolve_from_root(root, args.xu_output),
                )
            )
        elif args.command == "strengths":
            output_dir = resolve_from_root(root, args.output_dir)
            _print(
                build_frozen_strengths(
                    root,
                    resolve_from_root(root, args.constants),
                    resolve_from_root(root, args.belgi_source),
                    output_dir,
                    j_max=args.j_max,
                    abs_k_max=args.abs_k_max,
                    vt_max=args.vt_max,
                    smu2_min_debye2=args.smu2_min,
                    compiler=(
                        resolve_from_root(root, args.compiler)
                        if args.compiler is not None else None
                    ),
                    require_xu_final=True,
                    xu_output=resolve_from_root(root, args.xu_output),
                    reference_levels_path=(
                        resolve_from_root(root, args.reference_levels)
                        if args.reference_levels is not None else None
                    ),
                )
            )
        elif args.command == "build-extended":
            _print(build_extended_products(root))
        elif args.command == "build-co-stretch":
            _print(
                build_co_stretch_products(
                    root,
                    compiler=(
                        resolve_from_root(root, args.compiler)
                        if args.compiler is not None else None
                    ),
                )
            )
        elif args.command == "build-ch3-rock":
            _print(build_ch3_rock_products(root))
        elif args.command == "build-vibrational-revision":
            _print(build_vibrational_revision(root, output_dir=args.output_dir,
                                               compiler=args.compiler))
        elif args.command == "all":
            fetch_sources(root, offline=args.offline)
            normalize_workspace(root)
            build_succeeded = True
            try:
                build_products(root)
            except MissingSourcesError as exc:
                build_succeeded = False
                write_blocker(root, exc.missing)
                if not args.allow_incomplete:
                    raise
            validation = validate_workspace(root)
            status = generate_reports(root)
            _print({"validation": validation, "status": status})
            if build_succeeded and not validation["production_usable"]:
                return 1
        return 0
    except MissingSourcesError as exc:
        write_blocker(root, exc.missing)
        print(str(exc), file=sys.stderr)
        return 2
    except FetchError as exc:
        print(str(exc), file=sys.stderr)
        return 3
    except HamiltonianError as exc:
        print(str(exc), file=sys.stderr)
        return 4


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
