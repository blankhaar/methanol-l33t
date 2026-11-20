"""Normalize captured reference inputs without upgrading their scientific status."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .basecol import (
    read_basecol_levels,
    read_basecol_rates,
    validate_level_rate_pair,
)
from .integrity import sha256_path, verify_raw_binding
from .lamda import read_lamda
from .observed import read_xu_observed_frequencies
from .quantum import from_legacy_lamda_label
from .xu import output_has_intensities, read_observed_lines, write_observed_csv


XU_FROZEN_LEVEL_DATASET = "local-belgi-xu2008-frozen-hamiltonian"
XU_CONSTANT_DIPOLE_DATASET = "local-belgi-xu2008-constant-dipole"
DAGDIGIAN_H2_DATASET = "dagdigian2024-h2"
RF2010_H2_DATASET = "rabli-flower2010-h2-vt012"
RF2010_HE_DATASET = "rabli-flower2010-he-vt012"
RF2011_HE_DATASET = "rabli-flower2011-he-torsion-changing"

# One entry per downloaded rate matrix.  BASECOL reuses a level table between
# the He and p-H2 calculations for each RF2010 torsional manifold, and between
# the p-/o-H2 Dagdigian calculations.  Keeping that relationship explicit
# avoids inventing duplicate level artifacts while preserving every collision
# ID in the normalized provenance.
BASECOL_INPUTS = (
    # Rabli & Flower (2010), torsion-preserving He.
    {"species": "A", "family": RF2010_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 463, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-A_IDT371_vt1_IDC463_vc1.lev", "rates": "CH3OH-A_He_statetostate_ID463_v1.kij"},
    {"species": "A", "family": RF2010_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 464, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-A_IDT372_vt1_IDC464_vc1.lev", "rates": "CH3OH-A_He_statetostate_ID464_v1.kij"},
    {"species": "A", "family": RF2010_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 465, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-A_IDT373_vt1_IDC465_vc1.lev", "rates": "CH3OH-A_He_statetostate_ID465_v1.kij"},
    {"species": "E", "family": RF2010_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 466, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-E_IDT374_vt1_IDC466_vc1.lev", "rates": "CH3OH-E_He_statetostate_ID466_v1.kij"},
    {"species": "E", "family": RF2010_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 467, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-E_IDT377_vt1_IDC467_vc1.lev", "rates": "CH3OH-E_He_statetostate_ID467_v1.kij"},
    {"species": "E", "family": RF2010_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 468, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-E_IDT378_vt1_IDC468_vc1.lev", "rates": "CH3OH-E_He_statetostate_ID468_v1.kij"},
    # Rabli & Flower (2010), torsion-preserving para-H2.
    {"species": "A", "family": RF2010_H2_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 469, "collider_id": 2, "directory": "basecol2023", "levels": "CH3OH-A_IDT371_vt1_IDC463_vc1.lev", "rates": "CH3OH-A_H2-para_statetostate_ID469_v1.kij"},
    {"species": "A", "family": RF2010_H2_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 470, "collider_id": 2, "directory": "basecol2023", "levels": "CH3OH-A_IDT372_vt1_IDC464_vc1.lev", "rates": "CH3OH-A_H2-para_statetostate_ID470_v1.kij"},
    {"species": "A", "family": RF2010_H2_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 472, "collider_id": 2, "directory": "basecol2023", "levels": "CH3OH-A_IDT373_vt1_IDC465_vc1.lev", "rates": "CH3OH-A_H2-para_statetostate_ID472_v1.kij"},
    {"species": "E", "family": RF2010_H2_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 473, "collider_id": 2, "directory": "basecol2023", "levels": "CH3OH-E_IDT374_vt1_IDC466_vc1.lev", "rates": "CH3OH-E_H2-para_statetostate_ID473_v1.kij"},
    {"species": "E", "family": RF2010_H2_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 474, "collider_id": 2, "directory": "basecol2023", "levels": "CH3OH-E_IDT377_vt1_IDC467_vc1.lev", "rates": "CH3OH-E_H2-para_statetostate_ID474_v1.kij", "maximum_usable_temperature_k": 170.0},
    {"species": "E", "family": RF2010_H2_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 475, "collider_id": 2, "directory": "basecol2023", "levels": "CH3OH-E_IDT378_vt1_IDC468_vc1.lev", "rates": "CH3OH-E_H2-para_statetostate_ID475_v1.kij"},
    # Rabli & Flower (2011), rotational and torsion-changing He.
    {"species": "A", "family": RF2011_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 476, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-A_IDT379_vt1_IDC476_vc1.lev", "rates": "CH3OH-A_He_statetostate_ID476_v1.kij"},
    {"species": "E", "family": RF2011_HE_DATASET, "release": "BASECOL V1 2023-04", "collision_id": 477, "collider_id": 6, "directory": "basecol2023", "levels": "CH3OH-E_IDT380_vt1_IDC477_vc1.lev", "rates": "CH3OH-E_He_statetostate_ID477_v1.kij"},
    # Dagdigian low-level vt=0 p-/o-H2 matrices.
    {"species": "A", "family": DAGDIGIAN_H2_DATASET, "release": "BASECOL V1 2025-07", "collision_id": 580, "collider_id": 3, "directory": "basecol2025", "levels": "CH3OH-A_IDT464_vt1_IDC580_vc1.lev", "rates": "CH3OH-A_H2-ortho_statetostate_ID580_v1.kij"},
    {"species": "A", "family": DAGDIGIAN_H2_DATASET, "release": "BASECOL V1 2025-07", "collision_id": 581, "collider_id": 2, "directory": "basecol2025", "levels": "CH3OH-A_IDT464_vt1_IDC580_vc1.lev", "rates": "CH3OH-A_H2-para_statetostate_ID581_v1.kij"},
    {"species": "E", "family": DAGDIGIAN_H2_DATASET, "release": "BASECOL V1 2025-07", "collision_id": 576, "collider_id": 3, "directory": "basecol2025", "levels": "CH3OH-E_IDT462_vt1_IDC576_vc1.lev", "rates": "CH3OH-E_H2-ortho_statetostate_ID576_v1.kij"},
    {"species": "E", "family": DAGDIGIAN_H2_DATASET, "release": "BASECOL V1 2025-07", "collision_id": 578, "collider_id": 2, "directory": "basecol2025", "levels": "CH3OH-E_IDT462_vt1_IDC576_vc1.lev", "rates": "CH3OH-E_H2-para_statetostate_ID578_v1.kij"},
)


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalize_legacy_lamda(source: Path, target_dir: Path, species: str) -> dict[str, int]:
    """Normalize a historical vt=0 file as a *reference*, not a vt=0..2 input."""
    data = read_lamda(source)
    key_by_index = {
        level.index: from_legacy_lamda_label(species, level.quantum_label, vt=0)
        for level in data.levels
    }
    stem = f"{species.lower()}-ch3oh-vt0-rf2010"
    level_rows = [
        {
            "species": species.upper(), "level_key": key_by_index[level.index].encode(),
            "source_index": level.index, "vt": 0, "j": key_by_index[level.index].j,
            "k": key_by_index[level.index].k, "parity": key_by_index[level.index].parity,
            "energy_cm1": f"{level.energy_cm1:.9g}", "weight": f"{level.weight:.9g}",
            "quantum_label_raw": level.quantum_label, "source": "LAMDA-RF2010-local-cache",
            "mapping_confidence": "exact_label",
        }
        for level in data.levels
    ]
    level_fields = list(level_rows[0])
    _write_csv(target_dir / f"{stem}-levels.csv", level_fields, level_rows)

    radiative_rows = [
        {
            "species": species.upper(), "transition_index": row.index,
            "upper_key": key_by_index[row.upper].encode(), "lower_key": key_by_index[row.lower].encode(),
            "a_s": f"{row.a_s:.12g}", "frequency_ghz": f"{row.frequency_ghz:.12g}",
            "source": "LAMDA-RF2010-local-cache", "status": "calculated_legacy_reference",
        }
        for row in data.radiative
    ]
    _write_csv(target_dir / f"{stem}-radiative.csv", list(radiative_rows[0]), radiative_rows)

    collision_rows: list[dict[str, object]] = []
    for partner in data.collision_partners:
        for row in partner.transitions:
            collision_rows.append(
                {
                    "species": species.upper(), "collider_id": partner.collider_id,
                    "collider_description": partner.description, "transition_index": row.index,
                    "upper_key": key_by_index[row.upper].encode(),
                    "lower_key": key_by_index[row.lower].encode(),
                    "temperatures_k": ";".join(f"{x:.9g}" for x in partner.temperatures_k),
                    "rates_cm3_s": ";".join(f"{x:.12g}" for x in row.rates_cm3_s),
                    "source": "Rabli-Flower-2010-pH2-vt0", "status": "calculated_legacy_reference",
                }
            )
    if collision_rows:
        _write_csv(target_dir / f"{stem}-collisions.csv", list(collision_rows[0]), collision_rows)
    return {
        "levels": len(level_rows), "radiative": len(radiative_rows),
        "collision_rows": len(collision_rows), "colliders": len(data.collision_partners),
    }


def _raw_binding(root: Path, relative: str, digest: str) -> dict[str, str]:
    artifact = {"path": relative, "sha256": digest}
    verify_raw_binding(root, artifact)
    return artifact


def normalize_xu_constant_dipole_master(root: Path) -> dict[str, object] | None:
    """Install the frozen Xu levels and radiative list into the master contract.

    Collision data remain a separate source family.  This function deliberately
    does not create an empty ``collisions.csv``: the build then reports that
    precise missing artifact instead of treating absent rates as zeros.
    """
    frozen_dir = root / "data" / "normalized" / "belgi-xu2008-frozen"
    strength_dir = root / "data" / "normalized" / "belgi-xu2008-constant-dipole"
    level_input = frozen_dir / "levels.csv"
    level_manifest_path = frozen_dir / "manifest.json"
    line_input = strength_dir / "lines.csv"
    observation_input = strength_dir / "observations.csv"
    strength_manifest_path = strength_dir / "manifest.json"
    model_completeness_path = strength_dir / "completeness.json"
    required = (
        level_input, level_manifest_path, line_input, strength_manifest_path,
        model_completeness_path,
    )
    if not all(path.is_file() for path in required):
        return None

    level_manifest = json.loads(level_manifest_path.read_text(encoding="utf-8"))
    strength_manifest = json.loads(strength_manifest_path.read_text(encoding="utf-8"))
    model_completeness = json.loads(model_completeness_path.read_text(encoding="utf-8"))
    if level_manifest.get("dataset_id") != XU_FROZEN_LEVEL_DATASET:
        raise ValueError("unexpected frozen Xu level dataset identity")
    if strength_manifest.get("dataset_id") != XU_CONSTANT_DIPOLE_DATASET:
        raise ValueError("unexpected Xu constant-dipole dataset identity")
    if model_completeness.get("dataset_id") != XU_CONSTANT_DIPOLE_DATASET:
        raise ValueError("unexpected Xu radiative completeness identity")

    constants_sha = str(strength_manifest["hamiltonian"]["constants_sha256"])
    dipole_sha = str(
        strength_manifest["dipole_model"]["provenance"]["sha256"]
    )
    xu_sha = str(
        level_manifest.get("validation", {})
        .get("xu_calculations", {})
        .get("xu_sha256", "")
    )
    constants_binding = _raw_binding(
        root, "data/raw/local-cache/xu2008-final-constants.txt", constants_sha
    )
    dipole_binding = _raw_binding(
        root, "data/raw/user-supplied/xu2008-dipole-model.txt", dipole_sha
    )
    level_bindings = [constants_binding]
    if xu_sha:
        level_bindings.append(
            _raw_binding(
                root, "data/raw/local-cache/xu2008-globalfit-output.txt", xu_sha
            )
        )

    observed_by_line = {}
    observed_summary: dict[str, object] | None = None
    radiative_bindings = [constants_binding, dipole_binding]
    if observation_input.is_file():
        expected_observation_sha = str(
            strength_manifest.get("outputs", {}).get("observations_sha256", "")
        )
        actual_observation_sha = sha256_path(observation_input)
        if expected_observation_sha and actual_observation_sha != expected_observation_sha:
            raise ValueError("Xu observation crosswalk does not match its strength manifest")
        observed_by_line, observed_summary = read_xu_observed_frequencies(
            observation_input
        )
        observation_manifest = strength_manifest.get("observed_frequencies", {})
        observation_source_sha = str(observation_manifest.get("source_sha256", ""))
        if not observation_source_sha:
            raise ValueError("Xu observation crosswalk lacks a raw source binding")
        if xu_sha and observation_source_sha != xu_sha:
            raise ValueError("Xu level and observation manifests bind different listings")
        observation_source_path = str(
            observation_manifest.get(
                "source_path", "data/raw/local-cache/xu2008-globalfit-output.txt"
            )
        )
        radiative_bindings.append(
            _raw_binding(root, observation_source_path, observation_source_sha)
        )

    master = root / "data" / "normalized" / "master"
    master.mkdir(parents=True, exist_ok=True)
    level_output = master / "levels.csv"
    level_fields = [
        "species", "level_key", "energy_cm1", "weight", "source",
        "source_index", "mapping_confidence", "raw_sha256",
    ]
    level_counts = {"A": 0, "E": 0}
    with level_input.open(encoding="utf-8", newline="") as source, level_output.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(target, fieldnames=level_fields, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            species = row["species"].upper()
            level_counts[species] += 1
            writer.writerow(
                {
                    "species": species,
                    "level_key": row["level_key"],
                    "energy_cm1": row["energy_cm1"],
                    "weight": row["weight"],
                    "source": XU_FROZEN_LEVEL_DATASET,
                    "source_index": row["source_index"],
                    "mapping_confidence": "exact_quantum_identity",
                    "raw_sha256": constants_sha,
                }
            )

    radiative_output = master / "radiative.csv"
    radiative_fields = [
        "species", "upper_key", "lower_key", "a_s", "frequency_ghz",
        "smu2_debye2", "source", "status", "frequency_source", "raw_sha256",
        "source_upper_id", "source_lower_id", "mapping_confidence",
    ]
    radiative_counts = {"A": 0, "E": 0}
    observed_frequency_counts = {"A": 0, "E": 0}
    with line_input.open(encoding="utf-8", newline="") as source, radiative_output.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(target, fieldnames=radiative_fields, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            species = row["species"].upper()
            radiative_counts[species] += 1
            predicted_frequency_mhz = float(row["frequency_mhz"])
            a_s = float(row["a_s"])
            observed = observed_by_line.get(int(row["transition_index"]))
            if observed is None:
                frequency_mhz = predicted_frequency_mhz
                frequency_source = XU_FROZEN_LEVEL_DATASET
                status = "calculated_constant_dipole"
            else:
                if (
                    observed.species != species
                    or observed.upper_key != row["upper_key"]
                    or observed.lower_key != row["lower_key"]
                ):
                    raise ValueError(
                        f"Xu observed frequency endpoint mismatch at line {row['transition_index']}"
                    )
                frequency_mhz = observed.frequency_mhz
                # S mu^2 and g_u are unchanged, so A_ul follows nu^3 exactly.
                a_s *= (frequency_mhz / predicted_frequency_mhz) ** 3
                frequency_source = observed.source_id
                status = "calculated_constant_dipole_observed_frequency"
                observed_frequency_counts[species] += 1
            writer.writerow(
                {
                    "species": species,
                    "upper_key": row["upper_key"],
                    "lower_key": row["lower_key"],
                    "a_s": f"{a_s:.17g}",
                    "frequency_ghz": f"{frequency_mhz / 1000.0:.17g}",
                    "smu2_debye2": row["smu2_debye2"],
                    "source": XU_CONSTANT_DIPOLE_DATASET,
                    "status": status,
                    "frequency_source": frequency_source,
                    "raw_sha256": dipole_sha,
                    "source_upper_id": row["upper_key"],
                    "source_lower_id": row["lower_key"],
                    "mapping_confidence": "exact_quantum_identity",
                }
            )

    master_completeness_path = master / "completeness.json"
    master_completeness: dict[str, object] = {}
    if master_completeness_path.is_file():
        master_completeness = json.loads(
            master_completeness_path.read_text(encoding="utf-8")
        )
    level_records: dict[str, object] = {}
    radiative_records: dict[str, object] = {}
    for species in ("A", "E"):
        level_records[species] = {
            "dataset_id": XU_FROZEN_LEVEL_DATASET,
            "selection_complete": True,
            "normalized_transition_count": level_counts[species],
            "raw_artifacts": level_bindings,
            "coverage": {"j_max": 22, "abs_k_max": 14, "vt_max": 2},
        }
        radiative_records[species] = {
            "dataset_id": XU_CONSTANT_DIPOLE_DATASET,
            "selection_complete": True,
            "normalized_transition_count": radiative_counts[species],
            "raw_artifacts": radiative_bindings,
            "coverage": {
                "j_max": 22,
                "abs_k_max": 14,
                "vt_max": 2,
                "rank": 1,
                "smu2_min_debye2": float(model_completeness["smu2_min_debye2"]),
            },
            "operator": model_completeness["operator"],
            "claim": model_completeness["claim"],
            "observed_frequency_selection": {
                "applied_line_count": observed_frequency_counts[species],
                "raw_listing_sha256": (
                    radiative_bindings[-1]["sha256"] if observed_summary else None
                ),
                "policy": observed_summary.get("policy") if observed_summary else None,
            },
        }
    master_completeness["levels"] = level_records
    master_completeness["radiative"] = radiative_records
    master_completeness["normalized_outputs"] = {
        "levels": {
            "path": str(level_output.relative_to(root)),
            "sha256": sha256_path(level_output),
        },
        "radiative": {
            "path": str(radiative_output.relative_to(root)),
            "sha256": sha256_path(radiative_output),
        },
    }
    master_completeness_path.write_text(
        json.dumps(master_completeness, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "level_rows": sum(level_counts.values()),
        "radiative_rows": sum(radiative_counts.values()),
        "level_rows_by_species": level_counts,
        "radiative_rows_by_species": radiative_counts,
        "observed_frequency_rows_by_species": observed_frequency_counts,
        "observed_frequency_selection": observed_summary,
        "collision_rows_installed": False,
        "collision_master_status": "missing_source_artifact",
    }


def normalize_basecol_collisions_master(root: Path) -> dict[str, object] | None:
    """Install every captured BASECOL matrix into the collision master.

    BASECOL source energies establish the direction of each downloaded rate and
    are retained verbatim for detailed-balance reorientation.  Output level
    energies always come from the canonical Xu/BELGI Hamiltonian.  The RF2010
    source basis contains two J=15, |K|=15 states per torsional manifold beyond
    the frozen canonical |K|<=14 scope; transitions touching those states are
    counted and explicitly excluded instead of expanding the radiative model.
    """
    official = root / "data" / "raw" / "official"
    all_paths = [
        official / str(record["directory"]) / str(record[name])
        for record in BASECOL_INPUTS for name in ("levels", "rates")
    ]
    if not any(path.is_file() for path in all_paths):
        return None
    missing = sorted({
        str(path.relative_to(root)) for path in all_paths if not path.is_file()
    })
    if missing:
        raise ValueError("incomplete BASECOL capture: " + ", ".join(missing))

    master = root / "data" / "normalized" / "master"
    level_master_path = master / "levels.csv"
    if not level_master_path.is_file():
        raise ValueError("canonical master levels must be normalized before BASECOL rates")
    with level_master_path.open(encoding="utf-8", newline="") as handle:
        canonical_keys = {row["level_key"] for row in csv.DictReader(handle)}

    collision_fields = [
        "species", "collider_id", "upper_key", "lower_key",
        "temperatures_k", "rates_cm3_s", "source", "status", "raw_sha256",
        "source_upper_id", "source_lower_id", "source_upper_energy_cm1",
        "source_lower_energy_cm1", "mapping_confidence",
    ]
    collision_output = master / "collisions.csv"
    basecol_families = {
        DAGDIGIAN_H2_DATASET, RF2010_H2_DATASET,
        RF2010_HE_DATASET, RF2011_HE_DATASET,
    }
    preserved_rows: list[dict[str, object]] = []
    if collision_output.is_file():
        with collision_output.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("source") not in basecol_families:
                    preserved_rows.append(
                        {field: row.get(field, "") for field in collision_fields}
                    )

    imported_rows: list[dict[str, object]] = []
    groups: dict[tuple[str, str], dict[str, object]] = {}
    for input_record in BASECOL_INPUTS:
        species = str(input_record["species"])
        family = str(input_record["family"])
        collision_id = int(input_record["collision_id"])
        collider_id = int(input_record["collider_id"])
        raw_dir = official / str(input_record["directory"])
        level_path = raw_dir / str(input_record["levels"])
        rate_path = raw_dir / str(input_record["rates"])
        source_levels = read_basecol_levels(level_path, species)
        source_rates = read_basecol_rates(rate_path)
        validate_level_rate_pair(source_levels, source_rates)
        source_by_index = {row.source_index: row for row in source_levels}
        in_scope_indices = {
            row.source_index for row in source_levels
            if row.key.encode() in canonical_keys
        }
        excluded_levels = [
            row for row in source_levels if row.source_index not in in_scope_indices
        ]

        maximum_usable = input_record.get("maximum_usable_temperature_k")
        if maximum_usable is None:
            usable_positions = tuple(range(len(source_rates.temperatures_k)))
            status = "calculated"
        else:
            usable_positions = tuple(
                index for index, temperature in enumerate(source_rates.temperatures_k)
                if temperature <= float(maximum_usable)
            )
            status = "calculated_source_grid_quarantined_180_200K"
        usable_temperatures = tuple(
            source_rates.temperatures_k[index] for index in usable_positions
        )
        if not usable_temperatures:
            raise ValueError(f"BASECOL ID{collision_id} has no usable temperature columns")

        level_relative = str(level_path.relative_to(root))
        rate_relative = str(rate_path.relative_to(root))
        level_sha = sha256_path(level_path)
        rate_sha = sha256_path(rate_path)
        level_binding = _raw_binding(root, level_relative, level_sha)
        rate_binding = _raw_binding(root, rate_relative, rate_sha)

        imported_count = 0
        excluded_transition_count = 0
        all_zero = 0
        partial_zero = 0
        vt_pairs: set[tuple[int, int]] = set()
        for transition in source_rates.transitions:
            if (
                transition.initial_level not in in_scope_indices
                or transition.final_level not in in_scope_indices
            ):
                excluded_transition_count += 1
                continue
            upper = source_by_index[transition.initial_level]
            lower = source_by_index[transition.final_level]
            if upper.source_energy_cm1 < lower.source_energy_cm1:
                raise ValueError(
                    f"BASECOL ID{collision_id} source row {upper.source_index} -> "
                    f"{lower.source_index} is upward in the source energy basis"
                )
            usable_rates = tuple(
                transition.rates_cm3_s[index] for index in usable_positions
            )
            zero_count = sum(rate == 0.0 for rate in usable_rates)
            all_zero += zero_count == len(usable_rates)
            partial_zero += 0 < zero_count < len(usable_rates)
            vt_pairs.add(tuple(sorted((upper.key.vt, lower.key.vt))))
            imported_rows.append(
                {
                    "species": species,
                    "collider_id": collider_id,
                    "upper_key": upper.key.encode(),
                    "lower_key": lower.key.encode(),
                    "temperatures_k": ";".join(
                        f"{value:.9g}" for value in usable_temperatures
                    ),
                    "rates_cm3_s": ";".join(
                        f"{value:.12g}" for value in usable_rates
                    ),
                    "source": family,
                    "status": status,
                    "raw_sha256": rate_sha,
                    "source_upper_id": (
                        f"BASECOL:IDC{collision_id}:level:{upper.source_index}"
                    ),
                    "source_lower_id": (
                        f"BASECOL:IDC{collision_id}:level:{lower.source_index}"
                    ),
                    "source_upper_energy_cm1": f"{upper.source_energy_cm1:.12g}",
                    "source_lower_energy_cm1": f"{lower.source_energy_cm1:.12g}",
                    "mapping_confidence": "exact_quantum_identity",
                }
            )
            imported_count += 1

        group = groups.setdefault(
            (species, family),
            {
                "normalized_transition_count": 0,
                "raw_transition_count": 0,
                "out_of_canonical_scope_transition_count": 0,
                "out_of_canonical_scope_levels": set(),
                "all_zero_transition_rows": 0,
                "partially_zero_transition_rows": 0,
                "raw_artifacts": {},
                "datasets": [],
                "vt_pairs": set(),
                "releases": set(),
            },
        )
        group["normalized_transition_count"] = int(group["normalized_transition_count"]) + imported_count
        group["raw_transition_count"] = int(group["raw_transition_count"]) + len(source_rates.transitions)
        group["out_of_canonical_scope_transition_count"] = int(group["out_of_canonical_scope_transition_count"]) + excluded_transition_count
        group["all_zero_transition_rows"] = int(group["all_zero_transition_rows"]) + all_zero
        group["partially_zero_transition_rows"] = int(group["partially_zero_transition_rows"]) + partial_zero
        group["out_of_canonical_scope_levels"].update(
            row.key.encode() for row in excluded_levels
        )
        group["vt_pairs"].update(vt_pairs)
        group["releases"].add(str(input_record["release"]))
        group["raw_artifacts"][level_relative] = level_binding
        group["raw_artifacts"][rate_relative] = rate_binding
        dataset_summary: dict[str, object] = {
            "basecol_collision_id": collision_id,
            "collider_id": collider_id,
            "raw_level_count": len(source_levels),
            "mapped_level_count": len(in_scope_indices),
            "raw_transition_count": len(source_rates.transitions),
            "normalized_transition_count": imported_count,
            "out_of_canonical_scope_transition_count": excluded_transition_count,
            "temperature_grid_k": list(source_rates.temperatures_k),
            "usable_temperature_grid_k": list(usable_temperatures),
            "all_zero_transition_rows_on_usable_grid": all_zero,
            "partially_zero_transition_rows_on_usable_grid": partial_zero,
        }
        if maximum_usable is not None:
            dataset_summary["quarantined_temperature_columns_k"] = [
                value for value in source_rates.temperatures_k
                if value > float(maximum_usable)
            ]
            dataset_summary["quarantine_reason"] = (
                "ID474 has a discontinuous three-column splice: 26,151/32,640 "
                "rows change zero status or by >10x at 170->180 K; preserve raw "
                "columns but do not normalize, interpolate, or repair them"
            )
        group["datasets"].append(dataset_summary)

    imported_rows.sort(
        key=lambda row: (
            str(row["species"]), str(row["source"]),
            int(str(row["source_upper_id"]).split(":")[1].removeprefix("IDC")),
            int(str(row["source_upper_id"]).rsplit(":", 1)[1]),
            int(str(row["source_lower_id"]).rsplit(":", 1)[1]),
        )
    )
    _write_csv(collision_output, collision_fields, preserved_rows + imported_rows)

    completeness_path = master / "completeness.json"
    completeness = (
        json.loads(completeness_path.read_text(encoding="utf-8"))
        if completeness_path.is_file() else {}
    )
    collision_records = completeness.setdefault("collisions", {})
    result_groups: dict[str, object] = {"A": {}, "E": {}}
    for (species, family), group in sorted(groups.items()):
        species_records = collision_records.setdefault(species, {})
        record = {
            "dataset_id": family,
            "database_releases": sorted(group["releases"]),
            "basecol_collision_ids": sorted(
                dataset["basecol_collision_id"] for dataset in group["datasets"]
            ),
            "selection_complete": True,
            "normalized_transition_count": group["normalized_transition_count"],
            "raw_artifacts": [
                group["raw_artifacts"][path]
                for path in sorted(group["raw_artifacts"])
            ],
            "coverage": {
                "raw_transition_count": group["raw_transition_count"],
                "out_of_canonical_scope_transition_count": group["out_of_canonical_scope_transition_count"],
                "out_of_canonical_scope_levels": sorted(group["out_of_canonical_scope_levels"]),
                "vt_pairs": [list(pair) for pair in sorted(group["vt_pairs"])],
                "complete_pair_triangle_in_each_raw_dataset": True,
                "all_zero_transition_rows_on_usable_grids": group["all_zero_transition_rows"],
                "partially_zero_transition_rows_on_usable_grids": group["partially_zero_transition_rows"],
                "negative_or_nonfinite_rate_values": 0,
                "datasets": sorted(
                    group["datasets"], key=lambda row: row["basecol_collision_id"]
                ),
            },
            "mapping": (
                "exact (vt,J,K,A parity/E signed-K) quantum identity; source "
                "energies retained for rate-direction provenance and detailed "
                "balance only; canonical Xu/BELGI energies used in products"
            ),
        }
        species_records[family] = record
        result_groups[species][family] = {
            "collision_rows": group["normalized_transition_count"],
            "raw_collision_rows": group["raw_transition_count"],
            "out_of_scope_rows": group["out_of_canonical_scope_transition_count"],
            "basecol_collision_ids": record["basecol_collision_ids"],
        }
    outputs = completeness.setdefault("normalized_outputs", {})
    outputs["collisions"] = {
        "path": str(collision_output.relative_to(root)),
        "sha256": sha256_path(collision_output),
    }
    completeness_path.write_text(
        json.dumps(completeness, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "dataset_families": sorted(basecol_families),
        "collision_rows": len(imported_rows),
        "collision_rows_by_species_and_family": result_groups,
        "canonical_energy_source": XU_FROZEN_LEVEL_DATASET,
        "id474_temperature_quarantine_k": [180.0, 190.0, 200.0],
    }


def normalize_basecol_rf2011_he_master(root: Path) -> dict[str, object] | None:
    """Backward-compatible entry point; all captured BASECOL families are normalized."""
    return normalize_basecol_collisions_master(root)


def normalize_workspace(root: str | Path) -> dict[str, object]:
    root = Path(root)
    raw = root / "data" / "raw" / "local-cache"
    target = root / "data" / "normalized" / "reference-vt0"
    result: dict[str, object] = {"reference_vt0": {}}
    for species in ("A", "E"):
        source = raw / f"lamda-{species.lower()}-ch3oh-vt0-rf2010.dat"
        if source.exists():
            result["reference_vt0"][species] = normalize_legacy_lamda(source, target, species)  # type: ignore[index]
    xu_path = raw / "xu2008-globalfit-output.txt"
    if xu_path.exists():
        text = xu_path.read_text(encoding="utf-8", errors="replace")
        lines = read_observed_lines(xu_path)
        write_observed_csv(root / "data" / "normalized" / "xu2008-observed-lines.csv", lines)
        result["xu2008"] = {
            "fitted_observed_lines_exported": len(lines),
            "intensities_enabled_in_captured_output": output_has_intensities(text),
            "captured_listing_contains_strengths": False,
            "reason": "captured BELGI globalfit output declares Intensity? 0",
        }
    master_result = normalize_xu_constant_dipole_master(root)
    if master_result is not None:
        result["xu2008_constant_dipole_master"] = master_result
    basecol_result = normalize_basecol_collisions_master(root)
    if basecol_result is not None:
        result["basecol_collision_master"] = basecol_result
        if master_result is not None:
            master_result["collision_rows_installed"] = True
            master_result["collision_master_status"] = "all_captured_basecol_families_complete"
    report = root / "reports" / "normalization.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result
