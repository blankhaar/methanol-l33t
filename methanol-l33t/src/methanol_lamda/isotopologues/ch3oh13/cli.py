"""Command-line interface for the separate 13CH3OH build."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .build import build_products, validate_products
from .extended import build_extended_products
from .he150 import build_he150_products, validate_he150_products


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="methanol-lamda-13ch3oh",
        description="Build or validate the A/E 13CH3OH vt=0,1,2 LAMDA files.",
    )
    subparsers = result.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build", help="run BELGI and generate both files")
    build.add_argument("--root", type=Path, default=Path.cwd())
    build.add_argument("--compiler", type=Path)
    validate = subparsers.add_parser("validate", help="validate installed products")
    validate.add_argument("--root", type=Path, default=Path.cwd())
    build_150 = subparsers.add_parser(
        "build-150",
        help="generate the exact RF2011 150-state He-proxy files",
    )
    build_150.add_argument("--root", type=Path, default=Path.cwd())
    build_150.add_argument("--compiler", type=Path)
    build_extended = subparsers.add_parser(
        "build-extended",
        help="generate the exploratory vt=0..4 A/E isotope candidates",
    )
    build_extended.add_argument("--root", type=Path, default=Path.cwd())
    build_extended.add_argument("--compiler", type=Path)
    validate_150 = subparsers.add_parser(
        "validate-150",
        help="validate the installed RF2011 150-state He-proxy files",
    )
    validate_150.add_argument("--root", type=Path, default=Path.cwd())
    return result


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    if arguments.command in {"build", "build-150", "build-extended"}:
        if arguments.command == "build":
            payload = build_products(arguments.root, compiler=arguments.compiler)
        elif arguments.command == "build-150":
            payload = build_he150_products(
                arguments.root, compiler=arguments.compiler
            )
        else:
            payload = build_extended_products(
                arguments.root, compiler=arguments.compiler
            )
        validation_status = (
            payload["validation"]["status"]
            if "validation" in payload
            else "pass"
        )
        summary = {
            "status": validation_status,
            "files": [
                {
                    "file": row["file"],
                    "levels": row["levels"],
                    "radiative_transitions": row["radiative_transitions"],
                    "collision_rows": row["collision_rows"],
                }
                for row in payload["files"]
            ],
        }
    elif arguments.command == "validate":
        summary = validate_products(arguments.root)
    else:
        summary = validate_he150_products(arguments.root)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


__all__ = ["main", "parser"]
