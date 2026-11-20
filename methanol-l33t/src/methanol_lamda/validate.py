"""Structural and scientific validation for generated and reference files."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Iterable

from .build import DEFAULT_MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE
from .integrity import (
    RawBindingError,
    sha256_path,
    verify_pinned_sha256,
    verify_raw_binding,
)
from .independent import read_summary
from .lamda import read_lamda, validate_lamda
from .physics import (
    downward_from_upward,
    einstein_a_from_catalog_intensity,
    einstein_a_from_smu2,
    interpolate_rate,
    upward_from_downward,
)
from .paths import portable_path, resolve_from_root
from .quantum import LevelKey


MASER_BENCHMARKS = (
    {
        "name": "6.7 GHz", "species": "A", "frequency_ghz": 6.6685192,
        "keys": frozenset(("A:vt0:J5:K1:+", "A:vt0:J6:K0:+")), "uncertainty_mhz": 0.0008,
    },
    {
        "name": "12.2 GHz", "species": "E", "frequency_ghz": 12.178596415,
        "keys": frozenset(("E:vt0:J2:K+0", "E:vt0:J3:K-1")), "uncertainty_mhz": 0.000135,
    },
    {
        "name": "23.1 GHz", "species": "A", "frequency_ghz": 23.1210242,
        "keys": frozenset(("A:vt0:J9:K2:+", "A:vt0:J10:K1:+")), "uncertainty_mhz": 0.0005,
    },
    {
        "name": "107.0 GHz", "species": "A", "frequency_ghz": 107.013803,
        "keys": frozenset(("A:vt0:J3:K1:+", "A:vt0:J4:K0:+")), "uncertainty_mhz": 0.005,
    },
)


def _einstein_benchmark_check(
    root: Path, policy: dict[str, object]
) -> tuple[list[str], dict[str, object]]:
    """Recompute raw-bound trusted A fixtures through both conversion paths."""
    path = root / "data" / "normalized" / "master" / "einstein-a-benchmarks.csv"
    errors: list[str] = []
    counts: Counter[tuple[str, str]] = Counter()
    maximum_relative_error = 0.0
    if not path.exists():
        return ["Einstein-A benchmark table is missing"], {"rows": 0}
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    tolerance = float(policy.get("relative_tolerance", 1.0e-5))
    minimum = int(policy.get("minimum_per_species_per_mode", 2))
    for number, row in enumerate(rows, start=2):
        species = row.get("species", "").upper()
        mode = row.get("mode", "")
        if species not in {"A", "E"} or mode not in {"smu2", "ordinary_catalog_intensity"}:
            errors.append(f"Einstein-A benchmark row {number} has invalid species/mode")
            continue
        try:
            verify_raw_binding(
                root, {"path": row["raw_path"], "sha256": row["raw_sha256"]}
            )
            frequency = float(row["frequency_mhz"])
            g_upper = float(row["g_upper"])
            trusted = float(row["trusted_a_s"])
            if mode == "smu2":
                computed = einstein_a_from_smu2(
                    frequency, float(row["smu2_debye2"]), g_upper
                )
            else:
                computed = einstein_a_from_catalog_intensity(
                    frequency,
                    float(row["catalog_intensity_nm2_mhz"]),
                    float(row["lower_energy_cm1"]),
                    g_upper,
                    float(row["partition_function"]),
                    float(row.get("temperature_k") or 300.0),
                )
        except (KeyError, ValueError, RawBindingError) as exc:
            errors.append(f"Einstein-A benchmark row {number} cannot be verified: {exc}")
            continue
        if trusted <= 0:
            errors.append(f"Einstein-A benchmark row {number} has non-positive trusted A")
            continue
        relative_error = abs(computed / trusted - 1.0)
        maximum_relative_error = max(maximum_relative_error, relative_error)
        if relative_error > tolerance:
            errors.append(
                f"Einstein-A benchmark row {number} relative error {relative_error:.6g} "
                f"exceeds {tolerance:.6g}"
            )
        if not row.get("source_id") or not row.get("source_row_id"):
            errors.append(f"Einstein-A benchmark row {number} lacks source identity")
        counts[(species, mode)] += 1
    for species in ("A", "E"):
        for mode in ("smu2", "ordinary_catalog_intensity"):
            if counts[(species, mode)] < minimum:
                errors.append(
                    f"Einstein-A benchmarks need {minimum} {species}/{mode} rows; "
                    f"found {counts[(species, mode)]}"
                )
    return errors, {
        "rows": len(rows),
        "counts": {f"{species}:{mode}": count for (species, mode), count in sorted(counts.items())},
        "maximum_relative_error": maximum_relative_error,
        "relative_tolerance": tolerance,
    }


def _scientific_checks(path: Path) -> dict[str, object]:
    data = read_lamda(path)
    result: dict[str, object] = {"errors": [], "warnings": [], "metrics": {}}
    errors: list[str] = result["errors"]  # type: ignore[assignment]
    warnings: list[str] = result["warnings"]  # type: ignore[assignment]
    errors.extend(validate_lamda(data))
    summary = read_summary(path)
    if (summary.level_count, summary.radiative_count) != (len(data.levels), len(data.radiative)):
        errors.append("independent reader disagrees with primary reader")
    independent_colliders = {row[0]: row[1:] for row in summary.collider_counts}
    primary_colliders = {
        partner.collider_id: (len(partner.transitions), len(partner.temperatures_k))
        for partner in data.collision_partners
    }
    if independent_colliders != primary_colliders:
        errors.append("independent reader disagrees on collision blocks")

    decoded: dict[int, LevelKey] = {}
    for level in data.levels:
        try:
            decoded[level.index] = LevelKey.decode(level.quantum_label)
        except ValueError:
            pass
    if len(decoded) != len(data.levels):
        warnings.append("legacy/non-canonical labels: vt balance and maser QN checks unavailable")
        status = "legacy_vt0_reference_not_classII" if "v=0" in data.molecule.casefold() else "reference_or_unclassified"
        if status == "legacy_vt0_reference_not_classII":
            warnings.append("file is vt=0 only and is structurally valid but scientifically unsuitable for Class II pumping")
        result["metrics"] = {
            "levels": len(data.levels), "radiative": len(data.radiative),
            "colliders": primary_colliders, "science_status": status,
        }
        return result

    he_only_basis = "[rf2011-he-only]" in data.molecule.casefold()
    species = {key.species for key in decoded.values()}
    if len(species) != 1:
        errors.append("A and E levels occur in the same file")
    species_name = next(iter(species)) if species else "?"
    templates: dict[int, set[tuple[int, int, str]]] = defaultdict(set)
    for key in decoded.values():
        templates[key.vt].add(key.rotational_identity)
    if set(templates) != {0, 1, 2}:
        errors.append(f"torsional manifolds present are {sorted(templates)}, expected [0, 1, 2]")
    elif he_only_basis:
        expected_counts = {0: 100, 1: 25, 2: 25}
        actual_counts = {vt: len(templates[vt]) for vt in (0, 1, 2)}
        if actual_counts != expected_counts:
            errors.append(
                f"RF2011 He basis has torsional counts {actual_counts}, expected {expected_counts}"
            )
        expected_jmax = {0: 9, 1: 4, 2: 4}
        actual_jmax = {
            vt: max(identity[0] for identity in templates[vt])
            for vt in (0, 1, 2)
        }
        if actual_jmax != expected_jmax:
            errors.append(
                f"RF2011 He basis has J maxima {actual_jmax}, expected {expected_jmax}"
            )
    elif not (templates[0] == templates[1] == templates[2]):
        for vt in (1, 2):
            missing = sorted(templates[0] - templates[vt])
            extra = sorted(templates[vt] - templates[0])
            if missing or extra:
                errors.append(f"vt={vt} template differs: {len(missing)} missing, {len(extra)} extra")

    edge_counts: Counter[tuple[int, int]] = Counter()
    degree: Counter[int] = Counter()
    adjacency: dict[int, set[int]] = defaultdict(set)
    transition_by_keys: dict[frozenset[str], float] = {}
    for transition in data.radiative:
        up, low = decoded[transition.upper], decoded[transition.lower]
        edge_counts[tuple(sorted((up.vt, low.vt)))] += 1
        degree[transition.upper] += 1
        degree[transition.lower] += 1
        adjacency[transition.upper].add(transition.lower)
        adjacency[transition.lower].add(transition.upper)
        transition_by_keys[frozenset((up.encode(), low.encode()))] = transition.frequency_ghz
    # Operationalize "substantial" pump coupling as at least one adjacent-
    # torsion edge per two rotational labels.  This is deliberately a floor,
    # not a completeness claim; the selection-complete source assertion is a
    # separate build gate.
    minimum_cross_edges = max(1, math.ceil(len(templates.get(0, ())) / 2))
    for pair in ((0, 1), (1, 2)):
        if edge_counts[pair] < minimum_cross_edges:
            errors.append(
                f"radiative graph has only {edge_counts[pair]} vt={pair[0]}<->{pair[1]} links; "
                f"minimum substantial-coupling floor is {minimum_cross_edges}"
            )
    isolated = [index for index in decoded if degree[index] == 0]
    if isolated:
        errors.append(f"{len(isolated)} levels have no radiative edge")

    def reachable(start: set[int]) -> set[int]:
        visited = set(start)
        frontier = list(start)
        while frontier:
            current = frontier.pop()
            for neighbor in adjacency[current] - visited:
                visited.add(neighbor)
                frontier.append(neighbor)
        return visited

    component_by_level: dict[int, frozenset[int]] = {}
    for index in decoded:
        if index not in component_by_level:
            component = frozenset(reachable({index}))
            for member in component:
                component_by_level[member] = component
    unpumped_ground = [
        index for index, key in decoded.items()
        if key.vt == 0
        and not {0, 1, 2} <= {decoded[member].vt for member in component_by_level[index]}
    ]
    if unpumped_ground:
        errors.append(
            f"{len(unpumped_ground)} vt=0 levels lack a radiative path through both vt=1 and vt=2"
        )

    collider_ids = {partner.collider_id for partner in data.collision_partners}
    required_colliders = {6} if he_only_basis else {2, 3, 6}
    missing_colliders = required_colliders - collider_ids
    extra_colliders = collider_ids - required_colliders if he_only_basis else set()
    if missing_colliders:
        errors.append(f"missing required collision partners {sorted(missing_colliders)}")
    if extra_colliders:
        errors.append(f"RF2011 He-only file has unexpected collision partners {sorted(extra_colliders)}")
    empty_colliders = [
        partner.collider_id for partner in data.collision_partners if not partner.transitions
    ]
    if empty_colliders:
        errors.append(f"collision partners with no state-to-state rows {sorted(empty_colliders)}")
    if he_only_basis:
        he_partner = next(
            (partner for partner in data.collision_partners if partner.collider_id == 6),
            None,
        )
        if he_partner is not None:
            if len(he_partner.transitions) != 11_175:
                errors.append(
                    f"RF2011 He block has {len(he_partner.transitions)} rows, expected 11175"
                )
            if he_partner.temperatures_k != tuple(float(value) for value in range(10, 401, 10)):
                errors.append("RF2011 He temperature grid is not 10..400 K in 10 K steps")

    masers: list[dict[str, object]] = []
    available_key_strings = {key.encode() for key in decoded.values()}
    for benchmark in MASER_BENCHMARKS:
        if benchmark["species"] != species_name:
            continue
        if he_only_basis and not set(benchmark["keys"]) <= available_key_strings:
            entry = dict(benchmark)
            entry["keys"] = sorted(benchmark["keys"])
            entry["found"] = False
            entry["status"] = "outside_RF2011_150_level_basis"
            masers.append(entry)
            continue
        actual = transition_by_keys.get(benchmark["keys"])
        entry = dict(benchmark)
        entry["keys"] = sorted(benchmark["keys"])
        entry["found"] = actual is not None
        if actual is None:
            errors.append(f"missing principal {benchmark['name']} transition by quantum identity")
        else:
            difference_khz = abs(actual - float(benchmark["frequency_ghz"])) * 1.0e6
            entry.update({"actual_ghz": actual, "difference_khz": difference_khz})
            tolerance_khz = max(2.0, 3.0 * float(benchmark["uncertainty_mhz"]) * 1000.0)
            if difference_khz > tolerance_khz:
                errors.append(f"{benchmark['name']} frequency differs by {difference_khz:.3f} kHz")
            encoded_to_index = {key.encode(): index for index, key in decoded.items()}
            start = {
                encoded_to_index[key] for key in benchmark["keys"]
                if key in encoded_to_index
            }
            reached_vt = {decoded[index].vt for index in reachable(start)} if start else set()
            entry["radiatively_reachable_vt"] = sorted(reached_vt)
            if not {0, 1, 2} <= reached_vt:
                errors.append(
                    f"principal {benchmark['name']} system has no radiative path through vt=1 and vt=2"
                )
        masers.append(entry)

    # A sampled numerical round trip tests the detailed-balance implementation
    # used when energy ordering changes during mapping.
    balance_max_relative = 0.0
    level_by_id = {level.index: level for level in data.levels}
    for partner in data.collision_partners:
        for row in partner.transitions[:: max(1, len(partner.transitions) // 25)]:
            upper, lower = level_by_id[row.upper], level_by_id[row.lower]
            gap = upper.energy_cm1 - lower.energy_cm1
            for temperature, down in zip(partner.temperatures_k, row.rates_cm3_s):
                if down == 0:
                    continue
                up = upward_from_downward(down, upper.weight, lower.weight, gap, temperature)
                recovered = downward_from_upward(up, upper.weight, lower.weight, gap, temperature)
                balance_max_relative = max(balance_max_relative, abs(recovered / down - 1.0))
    if balance_max_relative > 1e-12:
        errors.append(f"detailed-balance round trip relative error {balance_max_relative:.3g}")

    result["metrics"] = {
        "levels": len(data.levels), "radiative": len(data.radiative), "colliders": primary_colliders,
        "template_sizes": {str(vt): len(value) for vt, value in sorted(templates.items())},
        "radiative_vt_pairs": {f"{a}-{b}": count for (a, b), count in sorted(edge_counts.items())},
        "minimum_adjacent_torsion_edges": minimum_cross_edges,
        "isolated_levels": len(isolated), "maser_benchmarks": masers,
        "radiative_components": len(set(component_by_level.values())),
        "ground_levels_without_vt1_vt2_path": len(unpumped_ground),
        "detailed_balance_max_relative_error": balance_max_relative,
        "product_profile": "rf2011_he_only" if he_only_basis else "full_class_ii",
        "science_status": (
            "radiative_incomplete_candidate"
            if not errors and "radiative-incomplete" in path.name
            else "rf2011_he_only_candidate" if not errors and he_only_basis
            else "vt012_candidate" if not errors else "failed"
        ),
    }
    return result


def validate_paths(paths: Iterable[str | Path]) -> dict[str, object]:
    files: dict[str, object] = {}
    for raw_path in paths:
        path = Path(raw_path)
        try:
            files[str(path)] = _scientific_checks(path)
        except Exception as exc:
            files[str(path)] = {"errors": [f"validator exception: {type(exc).__name__}: {exc}"], "warnings": [], "metrics": {}}
    structural_passed = all(not result["errors"] for result in files.values())  # type: ignore[index,union-attr]
    production_usable = bool(files) and structural_passed and all(
        result.get("metrics", {}).get("science_status") == "vt012_candidate"  # type: ignore[union-attr]
        for result in files.values()
    )
    return {
        "files": files,
        "passed": structural_passed,
        "structural_passed": structural_passed,
        "production_usable": production_usable,
    }


def _check_sidecars(root: Path, path: Path, result: dict[str, object]) -> None:
    try:
        path.resolve().relative_to((root / "outputs").resolve())
    except ValueError:
        return
    if path.parent.resolve() != (root / "outputs").resolve():
        return
    errors: list[str] = result["errors"]  # type: ignore[assignment]
    data = read_lamda(path)
    policy_path = root / "config" / "build.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8")) if policy_path.exists() else {}
    minimum_effective_dipole = float(
        policy.get(
            "minimum_effective_transition_dipole_debye",
            DEFAULT_MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE,
        )
    )
    minimum_smu2 = minimum_effective_dipole**2
    directory = root / "outputs" / "provenance"
    rad_path = directory / f"{path.name}.radiative.csv"
    col_path = directory / f"{path.name}.collisions.csv"
    if not rad_path.exists():
        errors.append("radiative provenance sidecar is missing")
    else:
        with rad_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != len(data.radiative):
            errors.append("radiative provenance row count does not match LAMDA transitions")
        if any(row.get("file") != path.name for row in rows):
            errors.append("radiative provenance contains a different filename key")
        transitions = {row.index: row for row in data.radiative}
        level_keys = {level.index: level.quantum_label for level in data.levels}
        seen_indices: set[int] = set()
        for row in rows:
            try:
                index = int(row["transition_index"])
                a_s = float(row["a_s"])
                frequency = float(row["frequency_ghz"])
                smu2 = float(row["smu2_debye2"])
                effective_dipole = float(row["effective_transition_dipole_debye"])
            except (KeyError, ValueError):
                errors.append("radiative provenance contains malformed numeric data")
                continue
            if index in seen_indices:
                errors.append(f"radiative provenance duplicates transition {index}")
            seen_indices.add(index)
            transition = transitions.get(index)
            if transition is None:
                errors.append(f"radiative provenance references absent transition {index}")
                continue
            if row.get("upper_key") != level_keys[transition.upper] or row.get("lower_key") != level_keys[transition.lower]:
                errors.append(f"radiative provenance endpoint mismatch at transition {index}")
            if not math.isclose(a_s, transition.a_s, rel_tol=5e-11, abs_tol=0.0):
                errors.append(f"radiative provenance Einstein A mismatch at transition {index}")
            if not math.isclose(frequency, transition.frequency_ghz, rel_tol=5e-10, abs_tol=0.0):
                errors.append(f"radiative provenance frequency mismatch at transition {index}")
            if not (
                math.isfinite(smu2) and smu2 > 0
                and math.isfinite(effective_dipole) and effective_dipole > 0
            ):
                errors.append(
                    f"radiative provenance has non-positive/non-finite strength at transition {index}"
                )
            else:
                if smu2 < minimum_smu2:
                    errors.append(
                        f"radiative transition {index} is below the configured effective-dipole cutoff"
                    )
                if not math.isclose(
                    effective_dipole, math.sqrt(smu2), rel_tol=5e-13, abs_tol=0.0
                ):
                    errors.append(
                        f"radiative provenance effective dipole mismatch at transition {index}"
                    )
                try:
                    upper_key = LevelKey.decode(row["upper_key"])
                    expected_a_s = einstein_a_from_smu2(
                        frequency * 1000.0, smu2, 2 * upper_key.j + 1
                    )
                except (KeyError, ValueError):
                    errors.append(
                        f"radiative provenance cannot reconstruct Einstein A at transition {index}"
                    )
                else:
                    if not math.isclose(a_s, expected_a_s, rel_tol=5e-12, abs_tol=0.0):
                        errors.append(
                            f"radiative provenance S mu^2/Einstein A mismatch at transition {index}"
                        )
            if (
                len(row.get("raw_sha256", "")) != 64
                or not row.get("a_source")
                or not row.get("source_upper_id")
                or not row.get("source_lower_id")
                or row.get("mapping_confidence") not in {"exact", "exact_quantum_identity"}
            ):
                errors.append(f"radiative provenance is not raw-bound/exact at transition {index}")
    if not col_path.exists():
        errors.append("collision provenance sidecar is missing")
        return
    with col_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        compact_collision_sidecar = "rates_cm3_s" in (reader.fieldnames or ())
        rows = list(reader)
    if compact_collision_sidecar:
        _check_compact_collision_sidecar(root, path, data, rows, errors)
        return
    expected_values = sum(
        len(partner.transitions) * len(partner.temperatures_k) for partner in data.collision_partners
    )
    if len(rows) != expected_values:
        errors.append("collision provenance value count does not match LAMDA collision grids")
    if any(row.get("file") != path.name for row in rows):
        errors.append("collision provenance contains a different filename key")
    level_keys: dict[int, str] = {}
    for level in data.levels:
        try:
            level_keys[level.index] = LevelKey.decode(level.quantum_label).encode()
        except ValueError:
            break
    expected_values: dict[tuple[int, int, int, float], float] = {}
    for partner in data.collision_partners:
        for transition in partner.transitions:
            for temperature, rate in zip(partner.temperatures_k, transition.rates_cm3_s):
                expected_values[(partner.collider_id, transition.upper, transition.lower, temperature)] = rate
    sidecar_values: dict[tuple[int, int, int, float], float] = {}
    for row in rows:
        try:
            key = (
                int(row["collider_id"]), int(row["upper_index"]),
                int(row["lower_index"]), float(row["temperature_k"]),
            )
            rate = float(row["rate_cm3_s"])
        except (KeyError, ValueError):
            errors.append("collision provenance contains a malformed numeric identity or value")
            continue
        if key in sidecar_values:
            errors.append(f"collision provenance duplicates value identity {key}")
        sidecar_values[key] = rate
        if level_keys:
            if row.get("upper_key") != level_keys.get(key[1]) or row.get("lower_key") != level_keys.get(key[2]):
                errors.append(f"collision provenance quantum identity disagrees at {key}")
    if set(sidecar_values) != set(expected_values):
        missing = len(set(expected_values) - set(sidecar_values))
        extra = len(set(sidecar_values) - set(expected_values))
        errors.append(f"collision provenance/LAMDA identities differ: {missing} missing, {extra} extra")
    for key in set(sidecar_values) & set(expected_values):
        if not math.isclose(sidecar_values[key], expected_values[key], rel_tol=5e-8, abs_tol=0.0):
            errors.append(f"collision provenance rate disagrees with LAMDA value at {key}")
    if any(not row.get("source") or not row.get("transform") for row in rows):
        errors.append("collision provenance omits source or transformation method")
    level_by_index = {level.index: level for level in data.levels}
    for row in rows:
        status = row.get("status", "")
        if status not in {"calculated", "proxy"}:
            continue
        try:
            upper_index = int(row["upper_index"])
            lower_index = int(row["lower_index"])
            target_temperature = float(row["temperature_k"])
            target_rate = float(row["rate_cm3_s"])
            source_temperatures = tuple(
                float(value) for value in row["source_temperatures_k"].split(";") if value
            )
            source_rates = tuple(
                float(value) for value in row["source_rates_cm3_s"].split(";") if value
            )
            source_upper_energy = float(row["source_upper_energy_cm1"])
            source_lower_energy = float(row["source_lower_energy_cm1"])
            multiplier = float(row["multiplier"])
        except (KeyError, ValueError):
            errors.append("calculated/proxy collision provenance is numerically incomplete")
            continue
        if (
            len(row.get("raw_sha256", "")) != 64
            or not row.get("source_upper_id") or not row.get("source_lower_id")
            or row.get("mapping_confidence") not in {"exact", "exact_quantum_identity"}
            or not source_temperatures or len(source_temperatures) != len(source_rates)
            or source_upper_energy < source_lower_energy
        ):
            errors.append("calculated/proxy collision provenance is not raw-bound and exact")
            continue
        transform = row["transform"]
        canonical_upper = row["upper_key"]
        canonical_lower = row["lower_key"]
        reoriented = "detailed_balance_reorientation" in transform
        if reoriented:
            endpoints_match = (
                row.get("source_upper_key") == canonical_lower
                and row.get("source_lower_key") == canonical_upper
            )
        else:
            endpoints_match = (
                row.get("source_upper_key") == canonical_upper
                and row.get("source_lower_key") == canonical_lower
            )
        if not endpoints_match:
            errors.append("calculated/proxy source and canonical collision endpoints disagree")
            continue
        working_rates = source_rates
        if reoriented:
            upper = level_by_index[upper_index]
            lower = level_by_index[lower_index]
            source_gap = source_upper_energy - source_lower_energy
            if source_gap == 0:
                # Source upper is canonical lower after endpoint reversal.
                working_rates = tuple(
                    rate * lower.weight / upper.weight for rate in source_rates
                )
            else:
                working_rates = tuple(
                    upward_from_downward(
                        rate, lower.weight, upper.weight, source_gap, temperature
                    )
                    for rate, temperature in zip(source_rates, source_temperatures)
                )
        try:
            reconstructed = interpolate_rate(
                source_temperatures, working_rates, target_temperature
            ) * multiplier
        except ValueError:
            errors.append("calculated/proxy collision provenance cannot reconstruct output rate")
            continue
        if not math.isclose(reconstructed, target_rate, rel_tol=5e-8, abs_tol=0.0):
            errors.append("calculated/proxy collision provenance does not reconstruct output rate")
    statuses = {row.get("status", "") for row in rows}
    allowed = {"calculated", "calculated_zeroed_negative_artifact", "proxy", "propensity", "constant"}
    if not statuses <= allowed:
        errors.append(f"collision provenance has unsupported statuses {sorted(statuses - allowed)}")
    if "-calculated" in path.name and any(not status.startswith("calculated") for status in statuses):
        errors.append("calculated-only filename contains approximate collision provenance")
    if "-rf2011-proxy" in path.name and not statuses <= {
        "calculated", "calculated_zeroed_negative_artifact", "proxy"
    }:
        errors.append("RF2011 filename contains non-RF approximation provenance")
    if "-cragg2005" in path.name and not statuses <= {
        "calculated", "calculated_zeroed_negative_artifact", "propensity", "constant"
    }:
        errors.append("Cragg filename contains an unrelated proxy status")


def _semicolon_floats(value: str) -> tuple[float, ...]:
    return tuple(float(part) for part in value.split(";") if part)


def _check_compact_collision_sidecar(
    root: Path,
    path: Path,
    data: object,
    rows: list[dict[str, str]],
    errors: list[str],
) -> None:
    """Validate the one-row-per-transition provenance used by dense matrices."""
    partners = {partner.collider_id: partner for partner in data.collision_partners}  # type: ignore[attr-defined]
    expected_count = sum(len(partner.transitions) for partner in partners.values())
    if len(rows) != expected_count:
        errors.append(
            f"compact collision provenance has {len(rows)} rows, expected {expected_count}"
        )
    level_by_index = {level.index: level for level in data.levels}  # type: ignore[attr-defined]
    expected = {
        (collider_id, transition.index): (partner, transition)
        for collider_id, partner in partners.items()
        for transition in partner.transitions
    }
    seen: set[tuple[int, int]] = set()
    raw_hashes: set[str] = set()
    for row in rows:
        try:
            collider_id = int(row["collider_id"])
            transition_index = int(row["transition_index"])
            upper_index = int(row["upper_index"])
            lower_index = int(row["lower_index"])
            temperatures = _semicolon_floats(row["temperatures_k"])
            rates = _semicolon_floats(row["rates_cm3_s"])
        except (KeyError, ValueError):
            errors.append("compact collision provenance contains malformed numeric data")
            continue
        identity = (collider_id, transition_index)
        if identity in seen:
            errors.append(f"compact collision provenance duplicates {identity}")
            continue
        seen.add(identity)
        target = expected.get(identity)
        if target is None:
            errors.append(f"compact collision provenance references absent row {identity}")
            continue
        partner, transition = target
        if (
            row.get("file") != path.name
            or upper_index != transition.upper
            or lower_index != transition.lower
            or row.get("upper_key") != level_by_index[upper_index].quantum_label
            or row.get("lower_key") != level_by_index[lower_index].quantum_label
        ):
            errors.append(f"compact collision provenance endpoint mismatch at {identity}")
        if temperatures != partner.temperatures_k or len(rates) != len(transition.rates_cm3_s):
            errors.append(f"compact collision provenance grid mismatch at {identity}")
            continue
        if any(
            not math.isclose(left, right, rel_tol=5e-8, abs_tol=0.0)
            for left, right in zip(rates, transition.rates_cm3_s)
        ):
            errors.append(f"compact collision provenance rate mismatch at {identity}")
        status = row.get("status", "")
        raw_derived = status.startswith("calculated") or status == "proxy"
        if (
            not row.get("source")
            or not status
            or not row.get("transform")
            or (not raw_derived and status not in {"propensity", "constant"})
        ):
            errors.append(f"compact collision provenance is incomplete at {identity}")
            continue
        if not raw_derived:
            # Analytic Cragg variants have no state-to-state raw rate row to
            # reconstruct.  Their complete formula/method label is retained in
            # ``transform``; equality to the LAMDA grid was checked above.
            continue

        try:
            source_temperatures = _semicolon_floats(row["source_temperatures_k"])
            source_rates = _semicolon_floats(row["source_rates_cm3_s"])
            source_upper_energy = float(row["source_upper_energy_cm1"])
            source_lower_energy = float(row["source_lower_energy_cm1"])
            multiplier = float(row.get("multiplier", "1"))
        except (KeyError, ValueError):
            errors.append(f"compact collision provenance is numerically incomplete at {identity}")
            continue
        raw_sha = row.get("raw_sha256", "")
        if (
            len(raw_sha) != 64
            or not row.get("source_upper_id")
            or not row.get("source_lower_id")
            or row.get("mapping_confidence") not in {"exact", "exact_quantum_identity"}
            or not source_temperatures
            or len(source_temperatures) != len(source_rates)
            or source_upper_energy < source_lower_energy
        ):
            errors.append(f"compact collision provenance is incomplete at {identity}")
            continue
        raw_hashes.add(raw_sha)
        transform = row["transform"]
        reoriented = "detailed_balance_reorientation" in transform
        if reoriented:
            endpoints_match = (
                row.get("source_upper_key") == row.get("lower_key")
                and row.get("source_lower_key") == row.get("upper_key")
            )
        else:
            endpoints_match = (
                row.get("source_upper_key") == row.get("upper_key")
                and row.get("source_lower_key") == row.get("lower_key")
            )
        if not endpoints_match:
            errors.append(f"compact collision source endpoint mismatch at {identity}")
            continue
        working_rates = source_rates
        if reoriented:
            source_gap = source_upper_energy - source_lower_energy
            canonical_upper = level_by_index[upper_index]
            canonical_lower = level_by_index[lower_index]
            if source_gap == 0:
                working_rates = tuple(
                    rate * canonical_lower.weight / canonical_upper.weight
                    for rate in source_rates
                )
            else:
                working_rates = tuple(
                    upward_from_downward(
                        rate,
                        canonical_lower.weight,
                        canonical_upper.weight,
                        source_gap,
                        temperature,
                    )
                    for rate, temperature in zip(source_rates, source_temperatures)
                )
        try:
            reconstructed = tuple(
                interpolate_rate(source_temperatures, working_rates, temperature)
                for temperature in temperatures
            )
        except ValueError:
            errors.append(f"compact collision provenance cannot reconstruct {identity}")
            continue
        reconstructed = tuple(multiplier * rate for rate in reconstructed)
        if any(
            not math.isclose(left, right, rel_tol=5e-8, abs_tol=0.0)
            for left, right in zip(reconstructed, rates)
        ):
            errors.append(f"compact collision transform does not reconstruct {identity}")
    if seen != set(expected):
        errors.append("compact collision provenance does not cover every LAMDA row")
    for digest in raw_hashes:
        try:
            verify_pinned_sha256(root, digest)
        except RawBindingError as exc:
            errors.append(f"compact collision provenance raw binding failed: {exc}")


def validate_workspace(root: str | Path, paths: Iterable[str | Path] | None = None) -> dict[str, object]:
    root = Path(root).resolve()
    manifest_path = root / "outputs" / "manifest.json"
    manifest: dict[str, object] = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if paths is None:
        declared = [
            root / "outputs" / row["file"]
            for row in manifest.get("files", [])  # type: ignore[union-attr]
            if (root / "outputs" / row["file"]).exists()
        ]
        generated = declared or sorted((root / "outputs").glob("*.dat"))
        paths = generated or sorted((root / "data" / "raw" / "local-cache").glob("lamda-*.dat"))
    else:
        paths = [resolve_from_root(root, path) for path in paths]
    report = validate_paths(paths)
    for raw_path, result in report["files"].items():  # type: ignore[union-attr]
        _check_sidecars(root, Path(raw_path), result)  # type: ignore[arg-type]
    report["files"] = {
        portable_path(root, raw_path): result
        for raw_path, result in report["files"].items()  # type: ignore[union-attr]
    }
    manifest_errors: list[str] = []
    if manifest:
        he_only_manifest = manifest.get("product_profile") == "rf2011_he_only"
        policy_path = root / "config" / "build.json"
        policy = json.loads(policy_path.read_text(encoding="utf-8")) if policy_path.exists() else {}
        expected_effective_dipole = float(
            policy.get(
                "minimum_effective_transition_dipole_debye",
                DEFAULT_MINIMUM_EFFECTIVE_TRANSITION_DIPOLE_DEBYE,
            )
        )
        cutoff = manifest.get("radiative_strength_cutoff", {})
        try:
            declared_effective_dipole = float(
                cutoff["minimum_effective_transition_dipole_debye"]  # type: ignore[index]
            )
            declared_smu2 = float(cutoff["minimum_smu2_debye2"])  # type: ignore[index]
        except (KeyError, TypeError, ValueError):
            manifest_errors.append("radiative strength cutoff is missing or malformed")
        else:
            if not math.isclose(
                declared_effective_dipole, expected_effective_dipole,
                rel_tol=0.0, abs_tol=1.0e-15,
            ) or not math.isclose(
                declared_smu2, expected_effective_dipole**2,
                rel_tol=0.0, abs_tol=1.0e-18,
            ):
                manifest_errors.append("manifest radiative strength cutoff differs from build policy")
        partition_thresholds = policy.get("minimum_partition_capture_fraction", {})
        benchmark_policy = policy.get("einstein_a_benchmark_policy", {})
        if benchmark_policy.get("required") is True and not he_only_manifest:
            benchmark_errors, benchmark_metrics = _einstein_benchmark_check(
                root, benchmark_policy
            )
            manifest_errors.extend(benchmark_errors)
            report["einstein_a_benchmarks"] = benchmark_metrics
        declared_names = {row["file"] for row in manifest.get("files", [])}  # type: ignore[union-attr]
        validated_names = {Path(path).name for path in report["files"]}
        if declared_names != validated_names:
            manifest_errors.append("validated file set differs from manifest")
        for row in manifest.get("files", []):  # type: ignore[union-attr]
            filename = row["file"]
            selection = row.get("radiative_strength_selection", {})
            try:
                candidate_count = int(selection["candidate_pairs_in_selected_level_set"])
                retained_count = int(selection["retained_transitions"])
                rejected_count = int(selection["rejected_below_cutoff"])
                selection_smu2 = float(selection["minimum_smu2_debye2"])
                declared_radiative_count = int(row["radiative_transitions"])
            except (KeyError, TypeError, ValueError):
                manifest_errors.append(
                    f"radiative strength selection is missing or malformed for {filename}"
                )
            else:
                matching_result = next(
                    (
                        result for validated_path, result in report["files"].items()  # type: ignore[union-attr]
                        if Path(validated_path).name == filename
                    ),
                    None,
                )
                actual_radiative_count = (
                    matching_result.get("metrics", {}).get("radiative")
                    if isinstance(matching_result, dict) else None
                )
                if candidate_count != retained_count + rejected_count:
                    manifest_errors.append(
                        f"radiative selection arithmetic is inconsistent for {filename}"
                    )
                if retained_count != declared_radiative_count or retained_count != actual_radiative_count:
                    manifest_errors.append(
                        f"radiative retained count differs from emitted rows for {filename}"
                    )
                if not math.isclose(
                    selection_smu2, expected_effective_dipole**2,
                    rel_tol=0.0, abs_tol=1.0e-18,
                ):
                    manifest_errors.append(
                        f"radiative selection cutoff differs from build policy for {filename}"
                    )
            capture = row.get("partition_capture_fraction_of_available_master", {})
            for temperature, threshold in (
                {} if he_only_manifest else partition_thresholds
            ).items():
                try:
                    actual_capture = float(capture[temperature])
                    required_capture = float(threshold)
                except (KeyError, TypeError, ValueError):
                    manifest_errors.append(
                        f"partition-capture metric missing for {filename} at {temperature}"
                    )
                    continue
                if actual_capture < required_capture:
                    manifest_errors.append(
                        f"partition capture {actual_capture:.6g} below {required_capture:.6g} "
                        f"for {filename} at {temperature}"
                    )
            artifacts = (
                (root / "outputs" / filename, row.get("sha256")),
                (
                    root / "outputs" / "provenance" / f"{filename}.radiative.csv",
                    row.get("radiative_sidecar_sha256"),
                ),
                (
                    root / "outputs" / "provenance" / f"{filename}.collisions.csv",
                    row.get("collision_sidecar_sha256"),
                ),
            )
            if row.get("level_sidecar_sha256"):
                artifacts += (
                    (
                        root / "outputs" / "provenance" / f"{filename}.levels.csv",
                        row.get("level_sidecar_sha256"),
                    ),
                )
            for artifact, expected_hash in artifacts:
                if not artifact.exists() or len(str(expected_hash or "")) != 64:
                    manifest_errors.append(f"manifest artifact missing or unpinned: {artifact.name}")
                elif sha256_path(artifact) != expected_hash:
                    manifest_errors.append(f"manifest checksum mismatch: {artifact.name}")
        for relative, expected_hash in manifest.get("build_inputs_sha256", {}).items():  # type: ignore[union-attr]
            artifact = root / relative
            if not artifact.exists() or sha256_path(artifact) != expected_hash:
                manifest_errors.append(f"build input changed or disappeared: {relative}")
    report["manifest_errors"] = manifest_errors
    structural_passed = all(
        not result["errors"] for result in report["files"].values()  # type: ignore[index,union-attr]
    )
    release_gate_passed = structural_passed and not manifest_errors
    report["passed"] = release_gate_passed
    report["structural_passed"] = structural_passed
    report["release_gate_passed"] = release_gate_passed
    expected_names = {
        f"{species}-ch3oh-vt012-{product}.dat"
        for species in ("a", "e")
        for product in ("calculated", "rf2011-proxy", "cragg2005")
    }
    actual_names = {Path(path).name for path in report["files"]}
    complete_product_set = actual_names == expected_names
    report["complete_product_set"] = complete_product_set
    report["production_usable"] = (
        complete_product_set and release_gate_passed and all(  # type: ignore[arg-type]
        result.get("metrics", {}).get("science_status") == "vt012_candidate"
        for result in report["files"].values()  # type: ignore[union-attr]
        )
    )
    validated_artifacts: dict[str, str] = {}
    for raw_path in report["files"]:
        path = resolve_from_root(root, raw_path)
        if path.exists():
            validated_artifacts[portable_path(root, path)] = sha256_path(path)
        for suffix in ("levels.csv", "radiative.csv", "collisions.csv"):
            sidecar = root / "outputs" / "provenance" / f"{path.name}.{suffix}"
            if sidecar.exists():
                validated_artifacts[portable_path(root, sidecar)] = sha256_path(sidecar)
    if manifest_path.exists():
        validated_artifacts[portable_path(root, manifest_path)] = sha256_path(manifest_path)
    report["validated_artifacts_sha256"] = validated_artifacts
    target = root / "reports" / "validation.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def provenance_coverage(root: str | Path) -> list[dict[str, object]]:
    """Aggregate collision sidecars into the requested coverage matrix."""
    root = Path(root)
    manifest_path = root / "outputs" / "manifest.json"
    allowed_files: set[str] | None = None
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        allowed_files = {row["file"] for row in manifest.get("files", [])}
    grouped: dict[tuple[str, str, int, int, int, str, str], dict[str, object]] = {}
    for path in sorted((root / "outputs" / "provenance").glob("*.collisions.csv")):
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if allowed_files is not None and row["file"] not in allowed_files:
                    continue
                upper, lower = LevelKey.decode(row["upper_key"]), LevelKey.decode(row["lower_key"])
                if row.get("temperatures_k"):
                    row_temperatures = _semicolon_floats(row["temperatures_k"])
                    rate_value_count = len(_semicolon_floats(row["rates_cm3_s"]))
                else:
                    row_temperatures = (float(row["temperature_k"]),)
                    rate_value_count = 1
                key = (
                    row["file"], upper.species, min(upper.vt, lower.vt), max(upper.vt, lower.vt),
                    int(row["collider_id"]), row["source"], row["status"],
                )
                entry = grouped.setdefault(
                    key,
                    {
                        "file": row["file"], "species": upper.species,
                        "vt_pair": f"{min(upper.vt, lower.vt)}-{max(upper.vt, lower.vt)}",
                        "collider_id": int(row["collider_id"]), "source": row["source"],
                        "status": row["status"], "j_min": min(upper.j, lower.j),
                        "j_max": max(upper.j, lower.j), "temperature_min_k": min(row_temperatures),
                        "temperature_max_k": max(row_temperatures), "rate_values": 0,
                        "_transitions": set(),
                    },
                )
                entry["j_min"] = min(int(entry["j_min"]), upper.j, lower.j)
                entry["j_max"] = max(int(entry["j_max"]), upper.j, lower.j)
                entry["temperature_min_k"] = min(float(entry["temperature_min_k"]), min(row_temperatures))
                entry["temperature_max_k"] = max(float(entry["temperature_max_k"]), max(row_temperatures))
                entry["rate_values"] = int(entry["rate_values"]) + rate_value_count
                entry["_transitions"].add((row["upper_index"], row["lower_index"]))  # type: ignore[union-attr]
    result = []
    for entry in grouped.values():
        transitions = entry.pop("_transitions")
        entry["transition_rows"] = len(transitions)  # type: ignore[arg-type]
        result.append(entry)
    return result
