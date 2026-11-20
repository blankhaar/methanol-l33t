"""Workspace-level audit and comparison reports."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from .integrity import sha256_path
from .lamda import read_lamda
from .paths import portable_path, resolve_from_root
from .physics import upward_from_downward
from .validate import provenance_coverage


DAGDIGIAN_H2_DATASET = "dagdigian2024-h2"
RF2010_H2_DATASET = "rabli-flower2010-h2-vt012"
COLLISION_OVERLAP_TEMPERATURES_K = (10.0, 20.0, 30.0, 50.0, 100.0, 150.0, 200.0)


def _semicolon_floats(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(";") if item)


def _percentile(values: list[float], fraction: float) -> float | None:
    """Return a linearly interpolated sample percentile without NumPy."""
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def _canonical_collision_rates(
    row: dict[str, str],
    levels: dict[str, tuple[float, float]],
) -> tuple[tuple[str, str], dict[float, float], bool, bool]:
    """Orient one source-downward row to the canonical-energy downward process.

    The normalized BASECOL rows retain the direction and term values of their
    source calculation.  A pair whose Xu/BELGI energy ordering is reversed
    therefore needs the inverse rate coefficient, obtained by detailed
    balance in the *source* energy basis.  This is the same convention used by
    the product builder.
    """
    first = row["upper_key"]
    second = row["lower_key"]
    try:
        first_energy, first_weight = levels[first]
        second_energy, second_weight = levels[second]
    except KeyError as exc:
        raise ValueError(f"collision overlap row references missing level {exc.args[0]}") from exc

    canonical_degenerate = first_energy == second_energy
    if first_energy > second_energy:
        canonical_upper, canonical_lower = first, second
    elif second_energy > first_energy:
        canonical_upper, canonical_lower = second, first
    else:
        # There is no physical up/down distinction at exact degeneracy.  A
        # lexical tie-break nevertheless puts both source families onto the
        # same directed coefficient for comparison.
        canonical_upper, canonical_lower = sorted((first, second), reverse=True)

    temperatures = _semicolon_floats(row["temperatures_k"])
    rates = _semicolon_floats(row["rates_cm3_s"])
    if not temperatures or len(temperatures) != len(rates):
        raise ValueError("malformed normalized collision temperature/rate grid")
    if any(
        not math.isfinite(value) or value <= 0
        for value in temperatures
    ) or any(not math.isfinite(rate) or rate < 0 for rate in rates):
        raise ValueError("non-finite, negative, or non-positive collision grid value")
    if any(right <= left for left, right in zip(temperatures, temperatures[1:])):
        raise ValueError("collision temperatures must increase strictly")

    reoriented = canonical_upper != first
    if reoriented:
        source_upper_energy = float(row["source_upper_energy_cm1"])
        source_lower_energy = float(row["source_lower_energy_cm1"])
        source_gap = source_upper_energy - source_lower_energy
        if source_gap < 0:
            raise ValueError("normalized collision row is upward in its source energy basis")
        if source_gap == 0:
            rates = tuple(rate * first_weight / second_weight for rate in rates)
        else:
            rates = tuple(
                upward_from_downward(
                    rate, first_weight, second_weight, source_gap, temperature
                )
                for temperature, rate in zip(temperatures, rates)
            )
    return (
        (canonical_upper, canonical_lower),
        {
            temperature: rate
            for temperature, rate in zip(temperatures, rates)
            if temperature in COLLISION_OVERLAP_TEMPERATURES_K
        },
        reoriented,
        canonical_degenerate,
    )


def write_collision_overlap_report(root: str | Path) -> dict[str, object] | None:
    """Compare overlapping Dagdigian and Rabli/Flower para-H2 vt=0 rates.

    Only exact shared quantum-state pairs and temperatures are compared.  Zero
    rates are classified separately rather than replaced by a floor, so ratio
    statistics contain only pairs for which both calculations are positive.
    """
    root = Path(root)
    level_path = root / "data" / "normalized" / "master" / "levels.csv"
    collision_path = root / "data" / "normalized" / "master" / "collisions.csv"
    if not level_path.is_file() or not collision_path.is_file():
        return None

    levels: dict[str, tuple[float, float]] = {}
    with level_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            levels[row["level_key"]] = (float(row["energy_cm1"]), float(row["weight"]))

    sources = (DAGDIGIAN_H2_DATASET, RF2010_H2_DATASET)
    rates_by_source: dict[str, dict[str, dict[tuple[str, str], dict[float, float]]]] = {
        source: {"A": {}, "E": {}} for source in sources
    }
    diagnostics: dict[str, dict[str, dict[str, int]]] = {
        source: {
            species: {
                "selected_vt0_para_h2_rows": 0,
                "detailed_balance_reoriented_rows": 0,
                "canonical_energy_degenerate_rows": 0,
            }
            for species in ("A", "E")
        }
        for source in sources
    }
    with collision_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            source = row["source"]
            if source not in rates_by_source or row["collider_id"] != "2":
                continue
            if ":vt0:" not in row["upper_key"] or ":vt0:" not in row["lower_key"]:
                continue
            species = row["species"].upper()
            if species not in {"A", "E"}:
                raise ValueError(f"unexpected methanol species {species!r}")
            pair, rates, reoriented, degenerate = _canonical_collision_rates(row, levels)
            if pair in rates_by_source[source][species]:
                raise ValueError(
                    f"duplicate {source} p-H2 vt=0 collision pair {species} {pair}"
                )
            rates_by_source[source][species][pair] = rates
            diagnostic = diagnostics[source][species]
            diagnostic["selected_vt0_para_h2_rows"] += 1
            diagnostic["detailed_balance_reoriented_rows"] += int(reoriented)
            diagnostic["canonical_energy_degenerate_rows"] += int(degenerate)

    result_rows: list[dict[str, object]] = []
    pair_counts: dict[str, int] = {}
    for species in ("A", "E"):
        dagdigian = rates_by_source[DAGDIGIAN_H2_DATASET][species]
        rabli_flower = rates_by_source[RF2010_H2_DATASET][species]
        overlap = sorted(set(dagdigian) & set(rabli_flower))
        pair_counts[species] = len(overlap)
        for temperature in COLLISION_OVERLAP_TEMPERATURES_K:
            both_positive: list[float] = []
            both_zero = 0
            dagdigian_zero = 0
            rabli_flower_zero = 0
            missing_temperature = 0
            for pair in overlap:
                dagdigian_rates = dagdigian[pair]
                rabli_flower_rates = rabli_flower[pair]
                if temperature not in dagdigian_rates or temperature not in rabli_flower_rates:
                    missing_temperature += 1
                    continue
                dagdigian_rate = dagdigian_rates[temperature]
                rabli_flower_rate = rabli_flower_rates[temperature]
                if dagdigian_rate > 0 and rabli_flower_rate > 0:
                    both_positive.append(dagdigian_rate / rabli_flower_rate)
                elif dagdigian_rate == 0 and rabli_flower_rate == 0:
                    both_zero += 1
                elif dagdigian_rate == 0:
                    dagdigian_zero += 1
                else:
                    rabli_flower_zero += 1
            compared = len(overlap) - missing_temperature
            outside_factor_2 = sum(ratio > 2.0 or ratio < 0.5 for ratio in both_positive)
            outside_factor_10 = sum(ratio > 10.0 or ratio < 0.1 for ratio in both_positive)
            ratio_count = len(both_positive)
            result_rows.append(
                {
                    "species": species,
                    "temperature_k": temperature,
                    "overlapping_quantum_pair_count": len(overlap),
                    "compared_pair_count": compared,
                    "missing_temperature_pair_count": missing_temperature,
                    "both_positive_count": ratio_count,
                    "both_zero_count": both_zero,
                    "dagdigian_zero_rabli_positive_count": dagdigian_zero,
                    "dagdigian_positive_rabli_zero_count": rabli_flower_zero,
                    "zero_mismatch_count": dagdigian_zero + rabli_flower_zero,
                    "ratio_dagdigian_over_rabli_median": _percentile(both_positive, 0.5),
                    "ratio_dagdigian_over_rabli_p10": _percentile(both_positive, 0.1),
                    "ratio_dagdigian_over_rabli_p90": _percentile(both_positive, 0.9),
                    "ratio_outside_factor_2_count": outside_factor_2,
                    "ratio_outside_factor_2_fraction": (
                        outside_factor_2 / ratio_count if ratio_count else None
                    ),
                    "ratio_outside_factor_10_count": outside_factor_10,
                    "ratio_outside_factor_10_fraction": (
                        outside_factor_10 / ratio_count if ratio_count else None
                    ),
                }
            )

    reports = root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    csv_path = reports / "collision-overlap.csv"
    fields = list(result_rows[0]) if result_rows else [
        "species", "temperature_k", "overlapping_quantum_pair_count",
        "compared_pair_count", "missing_temperature_pair_count", "both_positive_count",
        "both_zero_count", "dagdigian_zero_rabli_positive_count",
        "dagdigian_positive_rabli_zero_count", "zero_mismatch_count",
        "ratio_dagdigian_over_rabli_median", "ratio_dagdigian_over_rabli_p10",
        "ratio_dagdigian_over_rabli_p90", "ratio_outside_factor_2_count",
        "ratio_outside_factor_2_fraction", "ratio_outside_factor_10_count",
        "ratio_outside_factor_10_fraction",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(result_rows)

    payload: dict[str, object] = {
        "schema_version": 1,
        "collider": {"id": 2, "name": "para-H2"},
        "torsional_scope": "vt=0 to vt=0",
        "sources": {
            "numerator": DAGDIGIAN_H2_DATASET,
            "denominator": RF2010_H2_DATASET,
        },
        "comparison": (
            "Exact shared canonical quantum pairs and exact shared temperatures; "
            "each source-downward coefficient is first oriented to the Xu/BELGI "
            "canonical-energy downward process using detailed balance in its source "
            "energy basis when necessary"
        ),
        "zero_policy": (
            "No floor or imputation: both-zero and one-source-zero pairs are counted "
            "separately and excluded from rate-ratio summaries"
        ),
        "temperatures_k": list(COLLISION_OVERLAP_TEMPERATURES_K),
        "overlapping_quantum_pair_count_by_species": pair_counts,
        "source_diagnostics": diagnostics,
        "rows": result_rows,
    }
    json_path = reports / "collision-overlap.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _xu_radiative_variant(root: Path) -> dict[str, object] | None:
    path = (
        root / "data" / "normalized" / "belgi-xu2008-constant-dipole"
        / "manifest.json"
    )
    completeness = path.with_name("completeness.json")
    if not path.is_file() or not completeness.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        model = json.loads(completeness.read_text(encoding="utf-8"))
        if payload.get("dataset_id") != "local-belgi-xu2008-constant-dipole":
            return None
        coefficients = payload["dipole_model"]["coefficients"]
        scope = payload["scope"]
        graphs = payload["validation"]["radiative_graph"]
        return {
            "dataset_id": payload["dataset_id"],
            "status": "complete_with_release_caveats",
            "selection_complete_for_declared_operator": (
                model.get("claim")
                == "selection-complete for the declared Xu2008 constant-dipole operator"
            ),
            "candidate_count": int(scope["candidate_count"]),
            "retained_line_count": int(scope["retained_count"]),
            "mu_a_debye": float(coefficients["mu_a"]),
            "mu_b_debye": float(coefficients["mu_b"]),
            "cos_3gamma_terms_debye": {
                "mu_a3": float(coefficients["mu_a_cos_3gamma"]),
                "mu_b3": float(coefficients["mu_b_cos_3gamma"]),
                "mu_c3": float(coefficients["mu_c_cos_3gamma"]),
            },
            "graph_connected_within_each_species": all(
                int(graphs[species]["connected_component_count"]) == 1
                for species in ("A", "E")
            ),
            "isolated_level_count": sum(
                int(graphs[species]["isolated_node_count"])
                for species in ("A", "E")
            ),
        }
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def write_blocker(root: str | Path, missing: list[str] | tuple[str, ...]) -> Path:
    root = Path(root)
    radiative_variant = _xu_radiative_variant(root)
    payload = {
        "status": "blocked_scientific_integrity",
        "missing_artifacts": list(missing),
        "consequence": (
            "The requested balanced Class II product set was not emitted in "
            "this run; absent or unverified inputs were not fabricated."
        ),
        "minimum_unblock": [
            "Supply or regenerate every artifact listed in missing_artifacts",
            "Bind each recovered input to data/raw/SHA256SUMS and the appropriate completeness record",
            "Rerun normalize, build, validate, and report without weakening the declared source or release gates",
        ],
        "available_intermediate_product": "A/E RF2011 150-state He-only LAMDA pair",
        "primary_blocker": missing[0] if missing else "unspecified required artifact",
        "mekhtiev_torsion_dependent_model": "optional alternate radiative model",
    }
    if radiative_variant is not None:
        payload["available_radiative_variant"] = radiative_variant
    path = root / "reports" / "build-blocker.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def generate_reports(root: str | Path) -> dict[str, object]:
    root = Path(root)
    manifest_path = root / "outputs" / "manifest.json"
    manifest: dict[str, object] = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    declared_names = [row["file"] for row in manifest.get("files", [])]  # type: ignore[union-attr]
    output_files = (
        [root / "outputs" / name for name in declared_names if (root / "outputs" / name).exists()]
        if declared_names else sorted((root / "outputs").glob("*.dat"))
    )
    comparison: list[dict[str, object]] = []
    for path in output_files:
        data = read_lamda(path)
        sidecar = root / "outputs" / "provenance" / f"{path.name}.collisions.csv"
        status: dict[str, set[tuple[str, str, str]]] = {}
        if sidecar.exists():
            with sidecar.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    status.setdefault(row["status"], set()).add(
                        (row["collider_id"], row["upper_index"], row["lower_index"])
                    )
        comparison.append(
            {
                "file": path.name, "levels": len(data.levels), "radiative": len(data.radiative),
                **{f"collision_rows_{name}": len(keys) for name, keys in sorted(status.items())},
            }
        )
    reports = root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    collision_overlap = write_collision_overlap_report(root)
    if comparison:
        fields = sorted({key for row in comparison for key in row}, key=lambda key: (key != "file", key))
        with (reports / "product-comparison.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(comparison)
    coverage = provenance_coverage(root)
    if manifest:
        for file_record in manifest.get("files", []):
            for row in file_record.get("collision_coverage", []):
                if not row["missing_rows"]:
                    continue
                coverage.append(
                    {
                        "file": row["file"], "species": row["species"], "vt_pair": row["vt_pair"],
                        "collider_id": row["collider_id"], "source": "none", "status": "missing",
                        "j_min": row["j_min"], "j_max": row["j_max"],
                        "temperature_min_k": "", "temperature_max_k": "",
                        "transition_rows": row["missing_rows"], "rate_values": 0,
                    }
                )
    if coverage:
        fields = [
            "file", "species", "vt_pair", "collider_id", "j_min", "j_max",
            "temperature_min_k", "temperature_max_k", "source", "status",
            "transition_rows", "rate_values",
        ]
        with (reports / "rate-coverage.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(coverage)
    validation_path = reports / "validation.json"
    validation_current = False
    validation_production_usable = False
    if validation_path.exists():
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        expected_paths = {portable_path(root, path) for path in output_files}
        validated_paths = set(validation.get("files", {}))
        recorded_hashes = validation.get("validated_artifacts_sha256", {})
        hashes_current = bool(recorded_hashes)
        for raw_path, expected_hash in recorded_hashes.items():
            artifact = resolve_from_root(root, raw_path)
            if not artifact.exists() or sha256_path(artifact) != expected_hash:
                hashes_current = False
                break
        validation_current = (
            bool(output_files) and validated_paths == expected_paths and hashes_current
        )
        validation_production_usable = bool(validation.get("production_usable"))
    build_blocked = (reports / "build-blocker.json").exists()
    radiative_variant = _xu_radiative_variant(root)
    release_caveats = [
        "resolve the 142 A-species direct parity identities with a global level mapping",
        "broaden A-species strength validation and investigate the factor-two legacy discrepancy at 107.013831 GHz",
    ]
    status = {
        "production_outputs": len(output_files), "comparison_rows": len(comparison),
        "coverage_rows": len(coverage),
        "collision_overlap_available": collision_overlap is not None,
        "collision_overlap_rows": (
            len(collision_overlap["rows"]) if collision_overlap is not None else 0
        ),
        "collision_overlap_pair_count_by_species": (
            collision_overlap["overlapping_quantum_pair_count_by_species"]
            if collision_overlap is not None else {}
        ),
        "build_blocked": build_blocked,
        "build_blocker_primary": (
            "missing H2 and larger within-manifold He rates for the intended combined vt=0,1,2 products"
            if build_blocked else None
        ),
        "validation_current": validation_current,
        "validation_production_usable": validation_production_usable,
        "status_scope": "final A/E vt=0,1,2 LAMDA production set",
        "radiative_variant": radiative_variant,
        "release_caveats": release_caveats if radiative_variant is not None else [],
        "mekhtiev_torsion_dependent_model": {
            "blocks_selected_radiative_variant": False,
            "role": "optional alternate dipole model",
        },
        "scientifically_complete": (
            bool(output_files)
            and not build_blocked
            and all("radiative-incomplete" not in path.name for path in output_files)
            and validation_current
            and validation_production_usable
        ),
    }
    (reports / "status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return status
