"""Build the exploratory ``vt=0..4`` :sup:`13`CH3OH analogue.

The isotope-specific Xu et al. (2014) Hamiltonian is evaluated directly with
BELGI.  The collision construction is deliberately parallel to the main
isotopologue candidate: available p-H2, o-H2, and He rates are
temperature-extended under the owner-reviewed policy, with native vt1/vt2
and He torsion-changing rows preferred before exact-identity vt0 copies,
and only then reduced-mass scaled.  Absent broad o-H2 coverage is not
synthesized from p-H2.
"""

from __future__ import annotations

from collections import Counter
import json
import math
import os
from pathlib import Path
import tempfile

from ...build import _calculated_collisions
from ...collision_integrity import (
    assert_collision_selection,
    assert_parent_collision_inheritance,
    audit_basecol_collision_master,
)
from ...extended import (
    ExtendedBuildError,
    _assert_manifest_binds_file,
    _copy_vt0_intraband,
    _inherit_radiative_subset,
    _load_extended_levels,
    _prefer_native_calculated_collisions,
    _select_levels,
    _select_radiative,
)
from ...hamiltonian import build_frozen_levels, sha256_file
from ...independent import read_summary
from ...lamda import read_lamda, validate_lamda, write_lamda
from ...models import LamdaData, Level
from ...paths import portable_path
from ...radiative_integrity import audit_radiative_integrity
from ...strengths import build_frozen_strengths
from .build import (
    MAIN_COLLISION_MASTER_RELATIVE_PATH,
    STRENGTH_DIPOLE_MODEL_ID,
    _load_main_collision_master,
    _partners_13ch3oh,
    _scale_calculated_collisions,
    subproject_root,
)
from .collisions import ISOTOPOLOGUE_PROXY_STATUS, scaling_metadata
from .profile import (
    HAMILTONIAN_DATASET_ID,
    MOLECULAR_MASS_U,
    RADIATIVE_DATASET_ID,
    SPECIES,
    load_assignment_sets,
)


EXTENDED_CONFIG_RELATIVE_PATH = Path(
    "isotopologues/13ch3oh/config/extended-build.json"
)
LEVEL_OUTPUT_RELATIVE_PATH = Path(
    "data/normalized/belgi-xu2014-vt04-frozen"
)
STRENGTH_OUTPUT_RELATIVE_PATH = Path(
    "data/normalized/belgi-xu2014-vt04-constant-dipole"
)


def _load_policies(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    parent_path = root / "config/extended-build.json"
    isotope_path = root / EXTENDED_CONFIG_RELATIVE_PATH
    try:
        parent = json.loads(parent_path.read_text(encoding="utf-8"))
        isotope = json.loads(isotope_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ExtendedBuildError("missing or malformed extended-build policy") from exc
    if parent.get("torsional_states") != [0, 1, 2, 3, 4]:
        raise ExtendedBuildError("parent extended policy must select vt=0..4")
    if isotope.get("torsional_states") != [0, 1, 2, 3, 4]:
        raise ExtendedBuildError("13CH3OH extended policy must select vt=0..4")
    spectroscopy = isotope.get("spectroscopy")
    if not isinstance(spectroscopy, dict):
        raise ExtendedBuildError("13CH3OH extended spectroscopy policy is missing")
    generation_floor = float(spectroscopy["strength_generation_smu2_min_debye2"])
    product_floor = float(
        spectroscopy["final_minimum_effective_transition_dipole_debye"]
    ) ** 2
    parent_generation_floor = float(
        parent["strength_generation_smu2_min_debye2"]
    )
    parent_product_floor = (
        float(parent["minimum_effective_transition_dipole_debye"]) ** 2
    )
    if not math.isclose(
        generation_floor,
        parent_generation_floor,
        rel_tol=0.0,
        abs_tol=1.0e-18,
    ):
        raise ExtendedBuildError(
            "13CH3OH generation floor must match the parent extended policy"
        )
    if not math.isclose(
        product_floor, parent_product_floor, rel_tol=0.0, abs_tol=1.0e-18
    ):
        raise ExtendedBuildError(
            "13CH3OH product floor must match the parent extended policy"
        )
    if generation_floor > product_floor:
        raise ExtendedBuildError("generation floor cannot exceed the LAMDA product floor")
    return parent, isotope


def _generate_spectroscopy(
    root: Path,
    isotope_root: Path,
    parent_policy: dict[str, object],
    isotope_policy: dict[str, object],
    *,
    compiler: str | Path | None,
) -> tuple[Path, Path, dict[str, object], dict[str, object]]:
    assignments = load_assignment_sets(root)
    spectroscopy = isotope_policy["spectroscopy"]
    assert isinstance(spectroscopy, dict)
    level_dir = isotope_root / LEVEL_OUTPUT_RELATIVE_PATH
    strength_dir = isotope_root / STRENGTH_OUTPUT_RELATIVE_PATH
    level_manifest = build_frozen_levels(
        root,
        assignments.constants_path,
        root / "vendor/belgi",
        level_dir,
        j_max=int(parent_policy["j_max"]),
        abs_k_max=int(parent_policy["abs_k_max"]),
        vt_max=4,
        compiler=compiler,
        require_xu_final=False,
        dataset_id=HAMILTONIAN_DATASET_ID,
    )
    strength_manifest = build_frozen_strengths(
        root,
        assignments.constants_path,
        root / "vendor/belgi",
        strength_dir,
        j_max=int(parent_policy["j_max"]),
        abs_k_max=int(parent_policy["abs_k_max"]),
        vt_max=4,
        smu2_min_debye2=float(
            spectroscopy["strength_generation_smu2_min_debye2"]
        ),
        compiler=compiler,
        require_xu_final=False,
        dataset_id=RADIATIVE_DATASET_ID,
        reference_levels_path=level_dir / "levels.csv",
        dipole_model_name="Xu et al. (2014) 13CH3OH constant permanent dipoles",
        dipole_model_id=STRENGTH_DIPOLE_MODEL_ID,
        dipole_provenance_path=(
            isotope_root / "data/raw/user-supplied/xu-et-al-2014.pdf"
        ),
        completeness_claim=(
            "selection-complete above S mu^2=1e-6 D^2 for the Xu 2014 "
            "13CH3OH constant-permanent-dipole operator"
        ),
    )
    return level_dir, strength_dir, level_manifest, strength_manifest


def build_extended_products(
    project_root: str | Path,
    *,
    compiler: str | Path | None = None,
) -> dict[str, object]:
    """Generate and validate the A/E ``vt=0..4`` isotope candidates."""
    root = Path(project_root).resolve()
    isotope_root = subproject_root(root)
    collision_source_audit = audit_basecol_collision_master(root)
    parent_policy, isotope_policy = _load_policies(root)
    level_dir, strength_dir, level_manifest, strength_manifest = (
        _generate_spectroscopy(
            root,
            isotope_root,
            parent_policy,
            isotope_policy,
            compiler=compiler,
        )
    )
    all_levels = _load_extended_levels(level_dir / "levels.csv")
    collision_master = _load_main_collision_master(root)
    temperatures = tuple(
        float(value) for value in parent_policy["collision_temperature_grid_k"]
    )
    minimum_smu2 = (
        float(parent_policy["minimum_effective_transition_dipole_debye"]) ** 2
    )
    output_names = isotope_policy["output_files"]
    if not isinstance(output_names, dict):
        raise ExtendedBuildError("13CH3OH extended output names are missing")

    prepared: list[tuple[str, LamdaData]] = []
    records: list[dict[str, object]] = []
    upstream_manifest_path = isotope_root / "outputs/manifest.json"
    if not upstream_manifest_path.is_file():
        raise ExtendedBuildError(
            "missing canonical 13CH3OH vt=0..2 manifest; run build first"
        )
    for species in SPECIES:
        selected_levels, template, tie_breaks = _select_levels(
            all_levels,
            species,
            j_max=int(parent_policy["j_max"]),
            k_max=int(parent_policy["abs_k_max"]),
            ground_max_cm1=float(
                parent_policy["ground_manifold_max_energy_cm1"]
            ),
        )
        key_to_index = {
            level.key: index
            for index, level in enumerate(selected_levels, start=1)
        }
        native_rows = _calculated_collisions(
            collision_master,
            species,
            key_to_index,
            selected_levels,
            temperatures,
            allow_sqrt_endpoint_extrapolation=True,
        )
        ground_rows = [
            row
            for row in native_rows
            if row.upper_key.vt == row.lower_key.vt == 0
        ]
        copied_rows = _copy_vt0_intraband(
            ground_rows, selected_levels, key_to_index, temperatures
        )
        merged_rows, collision_precedence = _prefer_native_calculated_collisions(
            native_rows,
            copied_rows,
        )
        collisions = _scale_calculated_collisions(merged_rows)
        radiative = _select_radiative(
            strength_dir / "lines.csv",
            species,
            key_to_index,
            minimum_smu2_debye2=minimum_smu2,
        )
        upstream_path = (
            isotope_root
            / "outputs"
            / f"{species.lower()}-13ch3oh-vt012-reduced-mass-proxy.dat"
        )
        if not upstream_path.is_file():
            raise ExtendedBuildError(
                f"missing canonical upstream product {portable_path(root, upstream_path)}"
            )
        _assert_manifest_binds_file(upstream_manifest_path, upstream_path)
        upstream = read_lamda(upstream_path)
        radiative, inheritance = _inherit_radiative_subset(
            radiative,
            selected_levels,
            upstream,
        )
        upstream_sidecar = (
            isotope_root
            / "outputs/provenance"
            / f"{upstream_path.name}.radiative.csv"
        )
        if not upstream_sidecar.is_file():
            raise ExtendedBuildError(
                f"missing upstream radiative sidecar {portable_path(root, upstream_sidecar)}"
            )
        data = LamdaData(
            molecule=(
                f"{species}-13CH3OH vt=0..4 "
                "[exploratory intraband-copy reduced-mass proxy]"
            ),
            mass_amu=MOLECULAR_MASS_U,
            levels=[
                Level(index, row.energy_cm1, row.weight, row.key.encode())
                for index, row in enumerate(selected_levels, start=1)
            ],
            radiative=radiative,
            collision_partners=_partners_13ch3oh(collisions, temperatures),
        )
        assert_collision_selection(
            data, collisions, temperatures, context=upstream_path.name
        )
        collision_inheritance = assert_parent_collision_inheritance(
            upstream, data, context=upstream_path.name
        )
        radiative_audit = audit_radiative_integrity(
            data,
            context=str(output_names[species]),
            parent=upstream,
            minimum_smu2_debye2=minimum_smu2,
        )
        errors = validate_lamda(data)
        if errors:
            raise ExtendedBuildError(
                f"{species}-13CH3OH extended validation failed: {errors[:10]}"
            )
        filename = str(output_names[species])
        prepared.append((filename, data))
        collision_counts = Counter(
            (row.collider_id, row.upper_key.vt) for row in collisions
        )
        records.append(
            {
                "file": filename,
                "species": species,
                "levels": len(data.levels),
                "levels_by_vt": {
                    str(vt): sum(
                        row.key.vt == vt for row in selected_levels
                    )
                    for vt in range(5)
                },
                "rotational_template_size": len(template),
                "radiative_transitions": len(radiative),
                "radiative_integrity": radiative_audit,
                "radiative_inheritance": {
                    **inheritance,
                    "upstream_path": portable_path(root, upstream_path),
                    "upstream_sha256": sha256_file(upstream_path),
                    "upstream_manifest_path": portable_path(
                        root, upstream_manifest_path
                    ),
                    "upstream_manifest_sha256": sha256_file(
                        upstream_manifest_path
                    ),
                    "upstream_radiative_sidecar_path": portable_path(
                        root, upstream_sidecar
                    ),
                    "upstream_radiative_sidecar_sha256": sha256_file(
                        upstream_sidecar
                    ),
                },
                "collision_inheritance": {
                    **collision_inheritance,
                    "upstream_path": str(upstream_path.relative_to(root)),
                    "upstream_sha256": sha256_file(upstream_path),
                    "upstream_manifest_path": str(
                        upstream_manifest_path.relative_to(root)
                    ),
                    "upstream_manifest_sha256": sha256_file(
                        upstream_manifest_path
                    ),
                },
                "collision_rows": len(collisions),
                "collision_precedence": collision_precedence,
                "collision_rows_by_collider_and_vt": {
                    f"collider{collider}_vt{vt}": count
                    for (collider, vt), count in sorted(collision_counts.items())
                },
                "energy_tie_breaks": tie_breaks,
                "scientific_limit": (
                    "vt=3,4 are Xu-2014 Hamiltonian extrapolations; CO stretch absent"
                ),
            }
        )

    output_dir = isotope_root / "outputs/extended"
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": isotope_policy["status"],
        "isotopologue": "13CH3OH",
        "configuration": {
            "parent": parent_policy,
            "isotope": isotope_policy,
        },
        "spectroscopy": {
            "hamiltonian_source": "Xu et al. (2014), Table 2",
            "doi": "10.1016/j.jms.2014.06.005",
            "fit_performed": False,
            "vt3_vt4_status": "extrapolative prediction outside fitted data domain",
            "level_manifest": level_manifest,
            "strength_manifest": strength_manifest,
        },
        "collisions": {
            "status": ISOTOPOLOGUE_PROXY_STATUS,
            "source_isotopologue": "12CH3OH",
            "temperature_extrapolation": parent_policy[
                "collision_temperature_transform"
            ],
            "native_precedence": (
                "available vt1/vt2 and He torsion-changing rows precede any "
                "vt0 exact-identity intraband copy"
            ),
            "intraband_copy": (
                "vt0 exact rotational endpoint identity fills uncovered "
                "intramanifold pairs, including vt3/vt4"
            ),
            "interband_rows": "available He vt0-vt1, vt0-vt2, and vt1-vt2 only",
            "scaling_by_collider": {
                str(collider_id): {
                    "name": scaling_metadata(collider_id).collider_name,
                    "factor": scaling_metadata(collider_id).factor,
                    "formula": scaling_metadata(collider_id).formula,
                }
                for collider_id in (2, 3, 6)
            },
        },
        "build_inputs_sha256": {
            portable_path(root, root / "config/extended-build.json"): sha256_file(
                root / "config/extended-build.json"
            ),
            portable_path(root, root / EXTENDED_CONFIG_RELATIVE_PATH): sha256_file(
                root / EXTENDED_CONFIG_RELATIVE_PATH
            ),
            portable_path(
                root, root / MAIN_COLLISION_MASTER_RELATIVE_PATH
            ): sha256_file(root / MAIN_COLLISION_MASTER_RELATIVE_PATH),
            portable_path(root, level_dir / "levels.csv"): sha256_file(
                level_dir / "levels.csv"
            ),
            portable_path(root, strength_dir / "lines.csv"): sha256_file(
                strength_dir / "lines.csv"
            ),
        },
        "files": records,
        "collision_source_level_audit": collision_source_audit,
        "limitations": [
            "The Xu-2014 constants are rounded published values.",
            "vt=3 and vt=4 are predictions beyond the fitted torsional domain.",
            "Collision dynamics are reduced-mass-scaled main-isotopologue proxies.",
            "Only the available He torsion-changing collision rates are included; missing interband identities remain absent.",
            "The CO-stretch and CH3-rocking manifolds are absent.",
        ],
    }

    with tempfile.TemporaryDirectory(
        prefix=".ch3oh13-extended-stage-", dir=isotope_root
    ) as temporary:
        stage = Path(temporary)
        by_name = {str(record["file"]): record for record in records}
        for filename, data in prepared:
            path = stage / filename
            write_lamda(path, data)
            primary = read_lamda(path)
            independent = read_summary(path)
            if (
                independent.level_count != len(primary.levels)
                or independent.radiative_count != len(primary.radiative)
                or {row[0] for row in independent.collider_counts} != {2, 3, 6}
            ):
                raise ExtendedBuildError(
                    f"{filename} failed independent-reader cross-check"
                )
            record = by_name[filename]
            record["sha256"] = sha256_file(path)
            record["validation"] = {
                "primary_parser": "pass",
                "independent_parser": "pass",
                "format_validation": "pass",
            }
        manifest_path = stage / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        for path in stage.iterdir():
            os.replace(path, output_dir / path.name)
    return manifest


__all__ = ["build_extended_products"]
