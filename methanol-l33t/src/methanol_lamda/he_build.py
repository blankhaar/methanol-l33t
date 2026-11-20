"""Build the exact 150-state Rabli--Flower 2011 CH3OH--He products."""

from __future__ import annotations

from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import tempfile

from .build import (
    RF_TORSION_HE_SOURCE,
    XU_FROZEN_LEVEL_SOURCE,
    MasterCollision,
    MasterLevel,
    SelectedCollision,
    MissingSourcesError,
    _calculated_collisions,
    _frequency_overrides,
    _level_completeness,
    _radiative_completeness,
    _radiative_strength_cutoff,
    _select_radiative,
    load_master,
)
from .independent import read_summary
from .integrity import (
    RawBindingError,
    sha256_path,
    verify_dataset_record,
    verify_pinned_sha256,
)
from .lamda import read_lamda, validate_lamda, write_lamda
from .models import CollisionPartner, CollisionTransition, LamdaData, Level
from .quantum import LevelKey
from .radiative_integrity import audit_radiative_integrity


HE_TEMPERATURES = tuple(float(value) for value in range(10, 401, 10))


def _read_completeness(root: Path) -> dict[str, object]:
    path = root / "data" / "normalized" / "master" / "completeness.json"
    if not path.is_file():
        raise MissingSourcesError(["data/normalized/master/completeness.json"])
    return json.loads(path.read_text(encoding="utf-8"))


def _verified_he_source_rows(
    root: Path, species: str, rows: list[MasterCollision]
) -> tuple[list[MasterCollision], set[LevelKey]]:
    chosen = [
        row for row in rows
        if row.species == species and row.source == RF_TORSION_HE_SOURCE
    ]
    record = (
        _read_completeness(root)
        .get("collisions", {})
        .get(species, {})
        .get(RF_TORSION_HE_SOURCE, {})
    )
    try:
        verify_dataset_record(
            root, record, dataset_id=RF_TORSION_HE_SOURCE,
            normalized_transition_count=len(chosen),
        )
    except (AttributeError, RawBindingError) as exc:
        raise MissingSourcesError(
            [f"raw-bound BASECOL RF2011 {species}-CH3OH--He matrix: {exc}"]
        ) from exc
    bound_hashes = {
        artifact["sha256"] for artifact in record.get("raw_artifacts", [])
        if isinstance(artifact, dict) and "sha256" in artifact
    }
    if len(chosen) != 11_175:
        raise MissingSourcesError(
            [f"complete 11175-row BASECOL RF2011 {species} He matrix"]
        )
    for row in chosen:
        if (
            row.collider_id != 6
            or row.status != "calculated"
            or row.temperatures != HE_TEMPERATURES
            or len(row.rates) != len(HE_TEMPERATURES)
            or any(rate < 0 or not math.isfinite(rate) for rate in row.rates)
            or row.raw_sha256 not in bound_hashes
            or not row.source_upper_id or not row.source_lower_id
            or row.mapping_confidence != "exact_quantum_identity"
            or row.source_upper_energy_cm1 is None
            or row.source_lower_energy_cm1 is None
            or row.source_upper_energy_cm1 < row.source_lower_energy_cm1
        ):
            raise MissingSourcesError(
                [f"complete, exact, nonnegative {species} He collision provenance"]
            )
    keys = {row.upper_key for row in chosen} | {row.lower_key for row in chosen}
    pairs = {frozenset((row.upper_key, row.lower_key)) for row in chosen}
    if len(keys) != 150 or len(pairs) != len(chosen):
        raise MissingSourcesError(
            [f"complete 150-level pair triangle for {species}-CH3OH--He"]
        )
    expected = len(keys) * (len(keys) - 1) // 2
    if len(pairs) != expected:
        raise MissingSourcesError(
            [f"all {expected} unordered {species} He state pairs"]
        )
    return chosen, keys


def _select_he_levels(
    all_levels: list[MasterLevel], species: str, keys: set[LevelKey]
) -> list[MasterLevel]:
    selected = [
        row for row in all_levels if row.key.species == species and row.key in keys
    ]
    if len(selected) != 150 or {row.key for row in selected} != keys:
        raise MissingSourcesError([f"150 canonical Xu/BELGI {species} collision levels"])
    selected.sort(key=lambda row: (row.energy_cm1, row.key.encode()))
    if any(right.energy_cm1 <= left.energy_cm1 for left, right in zip(selected, selected[1:])):
        raise ValueError(f"{species} RF2011-basis canonical term values are not strictly unique")
    if not math.isclose(selected[0].energy_cm1, 0.0, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(f"{species} standalone term values are not zeroed at their ground state")
    for row in selected:
        if not math.isclose(row.weight, 2 * row.key.j + 1, rel_tol=0.0, abs_tol=1.0e-9):
            raise ValueError(f"unexpected statistical weight for {row.key.encode()}")
    return selected


def _write_level_sidecar(path: Path, filename: str, levels: list[MasterLevel]) -> None:
    fields = [
        "file", "level_index", "level_key", "energy_cm1", "weight",
        "energy_source", "energy_status", "source_index", "raw_sha256",
        "mapping_confidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(levels, start=1):
            writer.writerow(
                {
                    "file": filename,
                    "level_index": index,
                    "level_key": row.key.encode(),
                    "energy_cm1": f"{row.energy_cm1:.17g}",
                    "weight": f"{row.weight:.9g}",
                    "energy_source": row.source,
                    "energy_status": (
                        "Xu-2008 fitted-Hamiltonian term value constrained by observed lines; "
                        "BASECOL theoretical energy rejected as the output term value"
                    ),
                    "source_index": row.source_index,
                    "raw_sha256": row.raw_sha256,
                    "mapping_confidence": row.mapping_confidence,
                }
            )


def _write_radiative_sidecar(
    path: Path, filename: str, rows: list[dict[str, object]]
) -> None:
    fields = [
        "file", "transition_index", "upper_key", "lower_key", "a_s",
        "frequency_ghz", "smu2_debye2", "effective_transition_dipole_debye",
        "a_source", "frequency_source", "status",
        "raw_sha256", "source_upper_id", "source_lower_id",
        "mapping_confidence", "duplicate_candidates",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({"file": filename, **row})


def _write_compact_collision_sidecar(
    path: Path, filename: str, rows: list[SelectedCollision]
) -> None:
    fields = [
        "file", "transition_index", "collider_id", "upper_index", "lower_index",
        "upper_key", "lower_key", "temperatures_k", "rates_cm3_s", "source",
        "status", "transform", "raw_sha256", "source_upper_id",
        "source_lower_id", "source_upper_key", "source_lower_key",
        "source_upper_energy_cm1", "source_lower_energy_cm1",
        "mapping_confidence", "source_temperatures_k", "source_rates_cm3_s",
    ]
    ordered = sorted(rows, key=lambda row: (row.upper, row.lower))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(ordered, start=1):
            writer.writerow(
                {
                    "file": filename,
                    "transition_index": index,
                    "collider_id": row.collider_id,
                    "upper_index": row.upper,
                    "lower_index": row.lower,
                    "upper_key": row.upper_key.encode(),
                    "lower_key": row.lower_key.encode(),
                    "temperatures_k": ";".join(f"{value:.9g}" for value in HE_TEMPERATURES),
                    "rates_cm3_s": ";".join(f"{value:.12g}" for value in row.rates),
                    "source": row.source,
                    "status": row.status,
                    "transform": row.transform,
                    "raw_sha256": row.raw_sha256,
                    "source_upper_id": row.source_upper_id,
                    "source_lower_id": row.source_lower_id,
                    "source_upper_key": row.source_upper_key,
                    "source_lower_key": row.source_lower_key,
                    "source_upper_energy_cm1": row.source_upper_energy_cm1,
                    "source_lower_energy_cm1": row.source_lower_energy_cm1,
                    "mapping_confidence": row.mapping_confidence,
                    "source_temperatures_k": ";".join(
                        f"{value:.9g}" for value in row.source_temperatures
                    ),
                    "source_rates_cm3_s": ";".join(
                        f"{value:.12g}" for value in row.source_rates
                    ),
                }
            )


def _vt_pair_counts(rows: list[SelectedCollision]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        pair = sorted((row.upper_key.vt, row.lower_key.vt))
        counts[f"{pair[0]}-{pair[1]}"] += 1
    return dict(sorted(counts.items()))


def build_he_products(root: str | Path) -> dict[str, object]:
    """Emit A/E LAMDA files for the exact BASECOL RF2011 150-state basis."""
    root = Path(root).resolve()
    policy_path = root / "config" / "build.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8")) if policy_path.exists() else {}
    minimum_effective_dipole, minimum_smu2 = _radiative_strength_cutoff(policy)
    levels, radiative_rows, collision_rows = load_master(root)
    prepared: list[
        tuple[
            str, LamdaData, list[MasterLevel], list[dict[str, object]],
            list[SelectedCollision], dict[str, object],
        ]
    ] = []
    manifest: dict[str, object] = {
        "schema_version": 1,
        "product_profile": "rf2011_he_only",
        "scope": (
            "Exact BASECOL 150-state CH3OH--He collision basis; not the still-planned "
            "combined p-H2/o-H2/He Class II production set"
        ),
        "temperature_grid_k": list(HE_TEMPERATURES),
        "level_energy_policy": (
            "Xu et al. 2008 fitted-Hamiltonian term values are used; BASECOL theoretical "
            "term values are retained only in collision provenance"
        ),
        "radiative_frequency_policy": (
            "fitted exact-identity Xu observations first, curated overrides last, "
            "Xu fitted-Hamiltonian prediction otherwise"
        ),
        "radiative_strength_cutoff": {
            "interpretation": "mu_eff = sqrt(S mu^2)",
            "retention_rule": "mu_eff >= minimum; equivalently S mu^2 >= minimum squared",
            "minimum_effective_transition_dipole_debye": minimum_effective_dipole,
            "minimum_smu2_debye2": minimum_smu2,
        },
        "files": [],
    }

    for species in ("A", "E"):
        if not _level_completeness(root, species, levels):
            raise MissingSourcesError([f"raw-bound canonical {species} level table"])
        source_collisions, collision_keys = _verified_he_source_rows(
            root, species, collision_rows
        )
        selected_levels = _select_he_levels(levels, species, collision_keys)
        key_to_index = {
            level.key: index for index, level in enumerate(selected_levels, start=1)
        }
        lamda_levels = [
            Level(index, level.energy_cm1, level.weight, level.key.encode())
            for index, level in enumerate(selected_levels, start=1)
        ]

        complete_radiative_source = _radiative_completeness(
            root, species, radiative_rows
        )
        if complete_radiative_source is None:
            raise MissingSourcesError([f"selection-complete {species} radiative table"])
        strength_rows = [
            row for row in radiative_rows if row.source == complete_radiative_source
        ]
        radiative, radiative_provenance, radiative_selection = _select_radiative(
            strength_rows, species, key_to_index, selected_levels,
            _frequency_overrides(root, species),
            minimum_smu2_debye2=minimum_smu2,
        )
        for row in radiative_provenance:
            try:
                verify_pinned_sha256(root, str(row["raw_sha256"]))
            except RawBindingError as exc:
                raise MissingSourcesError(
                    [f"raw-bound {species} radiative provenance: {exc}"]
                ) from exc
        selected_collisions = _calculated_collisions(
            source_collisions, species, key_to_index, selected_levels, HE_TEMPERATURES
        )
        if len(selected_collisions) != 11_175:
            raise MissingSourcesError(
                [f"all 11175 mapped {species}-CH3OH--He collision rows"]
            )
        partner = CollisionPartner(
            6,
            f"CH3OH-He Rabli-Flower-2011 BASECOL-ID{476 if species == 'A' else 477}",
            HE_TEMPERATURES,
            [
                CollisionTransition(index, row.upper, row.lower, row.rates)
                for index, row in enumerate(
                    sorted(selected_collisions, key=lambda item: (item.upper, item.lower)),
                    start=1,
                )
            ],
        )
        filename = f"{species.lower()}-ch3oh-vt012-rf2011-he.dat"
        data = LamdaData(
            molecule=f"{species}-CH3OH vt=0,1,2 [RF2011-He-only]",
            mass_amu=32.04186,
            levels=lamda_levels,
            radiative=radiative,
            collision_partners=[partner],
        )
        radiative_audit = audit_radiative_integrity(
            data,
            context=filename,
            source_rows=radiative_provenance,
            minimum_smu2_debye2=minimum_smu2,
        )
        structural_errors = validate_lamda(data)
        if structural_errors:
            raise ValueError(
                f"{filename} failed preflight: " + "; ".join(structural_errors[:5])
            )
        status_counts = Counter(row.status for row in selected_collisions)
        transform_counts = Counter(row.transform for row in selected_collisions)
        file_record: dict[str, object] = {
            "file": filename,
            "species": species,
            "product": "rf2011_he_only",
            "levels": len(lamda_levels),
            "levels_by_vt": {
                str(vt): sum(level.key.vt == vt for level in selected_levels)
                for vt in (0, 1, 2)
            },
            "radiative_transitions": len(radiative),
            "radiative_integrity": radiative_audit,
            "radiative_strength_selection": radiative_selection,
            "radiative_observed_frequency_rows": sum(
                str(row["frequency_source"]).startswith(
                    "xu2008-globalfit-observation:"
                )
                for row in radiative_provenance
            ),
            "collision_rows": len(selected_collisions),
            "collision_rows_by_vt_pair": _vt_pair_counts(selected_collisions),
            "collision_rows_by_status": dict(sorted(status_counts.items())),
            "collision_transforms": dict(sorted(transform_counts.items())),
            "all_zero_collision_rows": sum(
                all(rate == 0.0 for rate in row.rates)
                for row in selected_collisions
            ),
            "collision_partner_ids": [6],
            "temperature_count": len(HE_TEMPERATURES),
            "temperature_min_k": HE_TEMPERATURES[0],
            "temperature_max_k": HE_TEMPERATURES[-1],
        }
        manifest["files"].append(file_record)  # type: ignore[union-attr]
        prepared.append(
            (
                filename, data, selected_levels, radiative_provenance,
                selected_collisions, file_record,
            )
        )

    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".methanol-he-stage-", dir=root) as directory:
        stage = Path(directory) / "outputs"
        provenance = stage / "provenance"
        provenance.mkdir(parents=True)
        for filename, data, selected_levels, rad_rows, col_rows, record in prepared:
            data_path = stage / filename
            level_path = provenance / f"{filename}.levels.csv"
            rad_path = provenance / f"{filename}.radiative.csv"
            col_path = provenance / f"{filename}.collisions.csv"
            write_lamda(data_path, data)
            _write_level_sidecar(level_path, filename, selected_levels)
            _write_radiative_sidecar(rad_path, filename, rad_rows)
            _write_compact_collision_sidecar(col_path, filename, col_rows)
            parsed = read_lamda(data_path)
            round_trip_errors = validate_lamda(parsed)
            if round_trip_errors:
                raise ValueError(
                    f"{filename} failed writer/parser round trip: "
                    + "; ".join(round_trip_errors[:5])
                )
            independent = read_summary(data_path)
            if (
                independent.level_count != len(data.levels)
                or independent.radiative_count != len(data.radiative)
                or independent.collider_counts != ((6, 11_175, 40),)
            ):
                raise ValueError(f"{filename} failed independent-reader cross-check")
            record.update(
                {
                    "sha256": sha256_path(data_path),
                    "level_sidecar_sha256": sha256_path(level_path),
                    "radiative_sidecar_sha256": sha256_path(rad_path),
                    "collision_sidecar_sha256": sha256_path(col_path),
                }
            )

        inputs: dict[str, str] = {}
        for relative in (
            "config/build.json",
            "config/observed-frequencies.json",
            "data/raw/SHA256SUMS",
            "data/normalized/master/levels.csv",
            "data/normalized/master/radiative.csv",
            "data/normalized/master/collisions.csv",
            "data/normalized/master/completeness.json",
        ):
            path = root / relative
            if path.is_file():
                inputs[relative] = sha256_path(path)
        manifest["build_inputs_sha256"] = inputs
        manifest["validation_status"] = "dual_parser_preflight_passed"
        manifest_path = stage / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (output_dir / "provenance").mkdir(parents=True, exist_ok=True)
        for path in sorted(provenance.iterdir()):
            os.replace(path, output_dir / "provenance" / path.name)
        for path in sorted(stage.glob("*.dat")):
            os.replace(path, output_dir / path.name)
        # Keep a profile-specific manifest even when a later full build makes
        # its own manifest.json the active six-product coherence marker.
        dedicated_manifest = output_dir / "rf2011-he-manifest.json"
        dedicated_manifest.write_text(
            manifest_path.read_text(encoding="utf-8"), encoding="utf-8"
        )
        os.replace(manifest_path, output_dir / "manifest.json")
    return manifest
