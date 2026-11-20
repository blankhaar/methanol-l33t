"""Build the owner-reviewed exploratory ``vt=0..4`` LAMDA pair.

Available native ``vt=1,2`` p-H2/He rows take precedence, native He
torsion-changing rows are retained, and exact ``vt=0`` identity copies fill
only missing intraband identities (including the added ``vt=3,4`` manifolds).
The much smaller o-H2 coverage remains visible rather than being filled from
p-H2.  CO-stretch and CH3-rocking successors are built by separate modules.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

from .build import MasterLevel, SelectedCollision, _calculated_collisions, _partners, load_master
from .collision_integrity import (
    assert_collision_selection,
    assert_parent_collision_inheritance,
    audit_basecol_collision_master,
)
from .hamiltonian import HamiltonianError
from .lamda import read_lamda, validate_lamda, write_lamda
from .models import LamdaData, Level, RadiativeTransition
from .quantum import LevelKey
from .radiative_integrity import audit_radiative_integrity


class ExtendedBuildError(HamiltonianError):
    """Raised when an exploratory extended input cannot be built safely."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _assert_manifest_binds_file(manifest_path: Path, product_path: Path) -> dict[str, object]:
    """Return the unique manifest row that checksum-binds ``product_path``."""

    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ExtendedBuildError(
            f"missing or malformed prerequisite manifest {manifest_path}"
        ) from exc
    files = payload.get("files")
    if not isinstance(files, list):
        raise ExtendedBuildError(
            f"prerequisite manifest {manifest_path} has no file ledger"
        )
    matches = [
        row
        for row in files
        if isinstance(row, dict) and row.get("file") == product_path.name
    ]
    if len(matches) != 1:
        raise ExtendedBuildError(
            f"prerequisite manifest does not bind exactly one {product_path.name}"
        )
    expected = str(matches[0].get("sha256", ""))
    actual = _sha256(product_path)
    if expected != actual:
        raise ExtendedBuildError(
            f"prerequisite manifest hash mismatch for {product_path.name}: "
            f"expected {expected or '<missing>'}, found {actual}"
        )
    return matches[0]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_extended_levels(path: Path) -> list[MasterLevel]:
    levels: list[MasterLevel] = []
    for row in _read_csv(path):
        key = LevelKey.decode(row["level_key"])
        levels.append(
            MasterLevel(
                key=key,
                energy_cm1=float(row["energy_cm1"]),
                weight=float(row["weight"]),
                source=row.get("source", "local-belgi-vt04"),
                source_index=row.get("source_index", ""),
                mapping_confidence="exact_generated_identity",
            )
        )
    return levels


def _select_levels(
    all_levels: list[MasterLevel],
    species: str,
    *,
    j_max: int,
    k_max: int,
    ground_max_cm1: float,
) -> tuple[list[MasterLevel], set[tuple[int, int, str]], list[dict[str, object]]]:
    by_vt: dict[int, dict[tuple[int, int, str], MasterLevel]] = {
        vt: {} for vt in range(5)
    }
    for level in all_levels:
        if level.key.species != species or level.key.vt not in by_vt:
            continue
        identity = level.key.rotational_identity
        if identity in by_vt[level.key.vt]:
            raise ExtendedBuildError(f"duplicate level {level.key.encode()}")
        by_vt[level.key.vt][identity] = level
    template = {
        identity
        for identity, level in by_vt[0].items()
        if level.key.j <= j_max
        and abs(level.key.k) <= k_max
        and level.energy_cm1 <= ground_max_cm1
    }
    if not template:
        raise ExtendedBuildError(f"empty {species} vt=0 rotational template")
    for vt in range(1, 5):
        missing = template - set(by_vt[vt])
        if missing:
            raise ExtendedBuildError(
                f"{species} vt={vt} lacks {len(missing)} selected rotational counterparts"
            )
    selected = [by_vt[vt][identity] for vt in range(5) for identity in template]
    selected.sort(key=lambda row: (row.energy_cm1, row.key.encode()))
    regularized: list[MasterLevel] = []
    adjustments: list[dict[str, object]] = []
    for row in selected:
        if regularized and row.energy_cm1 <= regularized[-1].energy_cm1:
            if row.energy_cm1 != regularized[-1].energy_cm1:
                raise ExtendedBuildError(
                    f"{species} levels cross under deterministic energy ordering"
                )
            original = row.energy_cm1
            row = replace(
                row, energy_cm1=math.nextafter(regularized[-1].energy_cm1, math.inf)
            )
            adjustments.append(
                {
                    "level_key": row.key.encode(),
                    "original_energy_cm1": original,
                    "serialized_energy_cm1": row.energy_cm1,
                    "method": "IEEE-754 nextafter; no physical splitting claimed",
                }
            )
        regularized.append(row)
    return regularized, template, adjustments


def _select_radiative(
    path: Path,
    species: str,
    key_to_index: dict[LevelKey, int],
    *,
    minimum_smu2_debye2: float,
) -> list[RadiativeTransition]:
    selected: list[RadiativeTransition] = []
    seen: set[frozenset[LevelKey]] = set()
    for row in _read_csv(path):
        if row["species"].upper() != species:
            continue
        upper = LevelKey.decode(row["upper_key"])
        lower = LevelKey.decode(row["lower_key"])
        if upper not in key_to_index or lower not in key_to_index:
            continue
        if float(row["smu2_debye2"]) < minimum_smu2_debye2:
            continue
        identity = frozenset((upper, lower))
        if identity in seen:
            raise ExtendedBuildError(
                f"duplicate extended radiative endpoints {upper.encode()} / {lower.encode()}"
            )
        seen.add(identity)
        selected.append(
            RadiativeTransition(
                index=0,
                upper=key_to_index[upper],
                lower=key_to_index[lower],
                a_s=float(row["a_s"]),
                frequency_ghz=float(row["frequency_mhz"]) / 1000.0,
            )
        )
    selected.sort(key=lambda row: (row.upper, row.lower))
    return [replace(row, index=index) for index, row in enumerate(selected, 1)]


def _inherit_radiative_subset(
    generated: list[RadiativeTransition],
    target_levels: list[MasterLevel],
    source: LamdaData,
) -> tuple[list[RadiativeTransition], dict[str, object]]:
    """Inherit an upstream radiative subnetwork by exact quantum identity.

    The generated target must contain exactly the same endpoint graph when it
    is restricted to the upstream level set.  Only Einstein A and frequency
    are inherited; target indices and the target level ordering are retained.
    This preserves reviewed experimental-frequency overlays through every
    larger molecular tier without using frequency proximity.
    """

    target_labels = {
        index: level.key.encode()
        for index, level in enumerate(target_levels, 1)
    }
    if len(set(target_labels.values())) != len(target_labels):
        raise ExtendedBuildError("target radiative inheritance labels are not unique")
    source_labels = {level.index: level.quantum_label for level in source.levels}
    if len(set(source_labels.values())) != len(source_labels):
        raise ExtendedBuildError("source radiative inheritance labels are not unique")
    target_label_set = set(target_labels.values())
    source_label_set = set(source_labels.values())
    missing_levels = source_label_set - target_label_set
    if missing_levels:
        raise ExtendedBuildError(
            f"target tier omits {len(missing_levels)} upstream level identities"
        )

    source_by_pair: dict[tuple[str, str], RadiativeTransition] = {}
    for row in source.radiative:
        try:
            pair = (source_labels[row.upper], source_labels[row.lower])
        except KeyError as exc:
            raise ExtendedBuildError(
                f"upstream radiative row {row.index} has an unresolved level index"
            ) from exc
        if pair in source_by_pair:
            raise ExtendedBuildError(
                f"duplicate upstream radiative identity {pair[0]} -> {pair[1]}"
            )
        source_by_pair[pair] = row

    generated_by_pair: dict[tuple[str, str], RadiativeTransition] = {}
    for row in generated:
        pair = (target_labels[row.upper], target_labels[row.lower])
        if pair in generated_by_pair:
            raise ExtendedBuildError(
                f"duplicate target radiative identity {pair[0]} -> {pair[1]}"
            )
        if pair[0] in source_label_set and pair[1] in source_label_set:
            generated_by_pair[pair] = row

    source_pairs = set(source_by_pair)
    target_pairs = set(generated_by_pair)
    if source_pairs != target_pairs:
        missing = sorted(source_pairs - target_pairs)
        extra = sorted(target_pairs - source_pairs)
        raise ExtendedBuildError(
            "upstream radiative endpoint graph does not match the generated "
            f"target subset (missing={missing[:3]}, extra={extra[:3]})"
        )

    inherited: list[RadiativeTransition] = []
    for row in generated:
        pair = (target_labels[row.upper], target_labels[row.lower])
        upstream = source_by_pair.get(pair)
        if upstream is None:
            inherited.append(row)
        else:
            inherited.append(
                replace(
                    row,
                    a_s=upstream.a_s,
                    frequency_ghz=upstream.frequency_ghz,
                    upper_energy_k=None,
                )
            )
    return inherited, {
        "upstream_level_count": len(source.levels),
        "inherited_radiative_rows": len(source_by_pair),
        "endpoint_policy": "exact ordered quantum-label identity; no frequency matching",
        "inherited_fields": ["a_s", "frequency_ghz"],
        "generated_target_only_rows": len(generated) - len(source_by_pair),
    }


def _copy_vt0_intraband(
    ground_rows: list[SelectedCollision],
    selected_levels: list[MasterLevel],
    key_to_index: dict[LevelKey, int],
    target_temperatures: tuple[float, ...],
) -> list[SelectedCollision]:
    level_by_key = {level.key: level for level in selected_levels}
    result: list[SelectedCollision] = []
    for ground in ground_rows:
        if ground.collider_id not in {2, 3, 6}:
            continue
        if ground.upper_key.vt != 0 or ground.lower_key.vt != 0:
            raise ExtendedBuildError("intraband-copy source must have two vt=0 endpoints")
        source_upper = level_by_key[ground.upper_key]
        source_lower = level_by_key[ground.lower_key]
        source_gap = source_upper.energy_cm1 - source_lower.energy_cm1
        if source_gap < 0:
            raise ExtendedBuildError("vt=0 collision source is not downward")
        for vt in range(5):
            upper_key = LevelKey(
                ground.upper_key.species,
                vt,
                ground.upper_key.j,
                ground.upper_key.k,
                ground.upper_key.parity,
            )
            lower_key = LevelKey(
                ground.lower_key.species,
                vt,
                ground.lower_key.j,
                ground.lower_key.k,
                ground.lower_key.parity,
            )
            if upper_key not in key_to_index or lower_key not in key_to_index:
                raise ExtendedBuildError("selected manifold lacks a collision-copy endpoint")
            rates = ground.rates
            transform = ground.transform
            if level_by_key[upper_key].energy_cm1 < level_by_key[lower_key].energy_cm1:
                upper_key, lower_key = lower_key, upper_key
                if source_gap == 0:
                    rates = tuple(
                        rate * source_upper.weight / source_lower.weight
                        for rate in ground.rates
                    )
                else:
                    rates = tuple(
                        rate
                        * (source_upper.weight / source_lower.weight)
                        * math.exp(-source_gap * 1.4387768775039338 / temperature)
                        for rate, temperature in zip(ground.rates, target_temperatures)
                    )
                transform += "; detailed_balance_for_excited_manifold_energy_reversal"
            if vt:
                transform += "; vt0_exact_rotational_identity_intraband_copy"
            result.append(
                replace(
                    ground,
                    upper=key_to_index[upper_key],
                    lower=key_to_index[lower_key],
                    upper_key=upper_key,
                    lower_key=lower_key,
                    rates=rates,
                    status=(
                        ground.status if vt == 0 else "vt0-intraband-identity-copy"
                    ),
                    transform=transform,
                )
            )
    identities = {(row.collider_id, row.upper_key, row.lower_key) for row in result}
    if len(identities) != len(result):
        raise ExtendedBuildError("duplicate collision identity after intraband copying")
    return result


def _prefer_native_calculated_collisions(
    native: list[SelectedCollision],
    intraband_copies: list[SelectedCollision],
) -> tuple[list[SelectedCollision], dict[str, int]]:
    """Merge collision rows with native calculated coefficients first.

    Exact-identity vt0 intraband copies fill only identities absent from the
    calculated network.  This keeps the available vt1/vt2 p-H2 and He rows,
    the available He torsion-changing rows, and the small o-H2 source subset;
    copies remain an explicit approximation for genuinely uncovered added
    manifolds.
    """

    def identity(row: SelectedCollision) -> tuple[int, frozenset[LevelKey]]:
        return row.collider_id, frozenset((row.upper_key, row.lower_key))

    merged = {identity(row): row for row in intraband_copies}
    if len(merged) != len(intraband_copies):
        raise ExtendedBuildError("duplicate intraband collision-copy identity")
    native_by_identity = {identity(row): row for row in native}
    if len(native_by_identity) != len(native):
        raise ExtendedBuildError("duplicate native calculated collision identity")
    overwritten = len(set(merged) & set(native_by_identity))
    merged.update(native_by_identity)
    result = sorted(
        merged.values(),
        key=lambda row: (row.collider_id, row.upper, row.lower),
    )
    return result, {
        "native_calculated_rows": len(native),
        "candidate_vt0_identity_copies": len(intraband_copies),
        "copies_suppressed_by_native_precedence": overwritten,
        "retained_identity_copy_rows": len(result) - len(native),
        "total_rows": len(result),
    }


def build_extended_products(root: str | Path) -> dict[str, object]:
    """Build and structurally validate exploratory A/E vt=0..4 files."""
    root = Path(root).resolve()
    policy_path = root / "config" / "extended-build.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    levels_path = root / "data" / "normalized" / "belgi-xu2008-vt04-frozen" / "levels.csv"
    lines_path = (
        root / "data" / "normalized" / "belgi-xu2008-vt04-constant-dipole" / "lines.csv"
    )
    strength_manifest_path = lines_path.parent / "manifest.json"
    missing = [
        str(path.relative_to(root))
        for path in (levels_path, lines_path, strength_manifest_path)
        if not path.is_file()
    ]
    if missing:
        raise ExtendedBuildError("missing generated extended spectroscopy: " + ", ".join(missing))
    strength_manifest = json.loads(strength_manifest_path.read_text(encoding="utf-8"))
    generated_floor = float(strength_manifest["selection"]["smu2_min_debye2"])
    policy_floor = float(policy["strength_generation_smu2_min_debye2"])
    if not math.isclose(generated_floor, policy_floor, rel_tol=0.0, abs_tol=1.0e-18):
        raise ExtendedBuildError(
            "extended strength artifact does not implement the configured generation floor"
        )
    all_levels = _load_extended_levels(levels_path)
    collision_source_audit = audit_basecol_collision_master(root)
    _, _, collision_rows = load_master(root)
    temperatures = tuple(float(value) for value in policy["collision_temperature_grid_k"])
    minimum_smu2 = float(policy["minimum_effective_transition_dipole_debye"]) ** 2
    output_dir = root / "outputs" / "extended"
    upstream_manifest_path = root / "outputs" / "manifest.json"
    if not upstream_manifest_path.is_file():
        raise ExtendedBuildError(
            "missing canonical vt=0..2 manifest; run the balanced build first"
        )
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": policy["status"],
        "configuration": policy,
        "input_sha256": {
            "policy": _sha256(policy_path),
            "levels": _sha256(levels_path),
            "radiative_lines": _sha256(lines_path),
            "radiative_manifest": _sha256(strength_manifest_path),
            "collision_master": _sha256(root / "data" / "normalized" / "master" / "collisions.csv"),
        },
        "files": [],
        "collision_source_level_audit": collision_source_audit,
    }
    with tempfile.TemporaryDirectory(prefix=".methanol-extended-stage-", dir=root) as temporary:
        stage = Path(temporary)
        for species in ("A", "E"):
            selected_levels, template, tie_breaks = _select_levels(
                all_levels,
                species,
                j_max=int(policy["j_max"]),
                k_max=int(policy["abs_k_max"]),
                ground_max_cm1=float(policy["ground_manifold_max_energy_cm1"]),
            )
            key_to_index = {
                level.key: index for index, level in enumerate(selected_levels, 1)
            }
            native_rows = _calculated_collisions(
                collision_rows,
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
            collisions, collision_precedence = _prefer_native_calculated_collisions(
                native_rows,
                copied_rows,
            )
            radiative = _select_radiative(
                lines_path,
                species,
                key_to_index,
                minimum_smu2_debye2=minimum_smu2,
            )
            upstream_path = (
                root
                / "outputs"
                / f"{species.lower()}-ch3oh-vt012-calculated.dat"
            )
            if not upstream_path.is_file():
                raise ExtendedBuildError(
                    f"missing canonical upstream product {upstream_path.relative_to(root)}"
                )
            _assert_manifest_binds_file(upstream_manifest_path, upstream_path)
            upstream = read_lamda(upstream_path)
            radiative, inheritance = _inherit_radiative_subset(
                radiative,
                selected_levels,
                upstream,
            )
            upstream_sidecar = (
                root
                / "outputs/provenance"
                / f"{upstream_path.name}.radiative.csv"
            )
            if not upstream_sidecar.is_file():
                raise ExtendedBuildError(
                    f"missing upstream radiative sidecar {upstream_sidecar.relative_to(root)}"
                )
            data = LamdaData(
                molecule=f"{species}-CH3OH vt=0..4 [exploratory intraband-copy]",
                mass_amu=32.04186,
                levels=[
                    Level(index, row.energy_cm1, row.weight, row.key.encode())
                    for index, row in enumerate(selected_levels, 1)
                ],
                radiative=radiative,
                collision_partners=_partners(collisions, temperatures),
            )
            assert_collision_selection(
                data, collisions, temperatures, context=upstream_path.name
            )
            collision_inheritance = assert_parent_collision_inheritance(
                upstream, data, context=upstream_path.name
            )
            filename = f"{species.lower()}-ch3oh-vt01234-intraband-copy.dat"
            radiative_audit = audit_radiative_integrity(
                data,
                context=filename,
                parent=upstream,
                minimum_smu2_debye2=minimum_smu2,
            )
            errors = validate_lamda(data)
            if errors:
                raise ExtendedBuildError(
                    f"{species} extended LAMDA validation failed: {errors[:10]}"
                )
            staged_path = stage / filename
            write_lamda(staged_path, data)
            collision_counts = Counter(
                (row.collider_id, row.upper_key.vt) for row in collisions
            )
            manifest["files"].append(  # type: ignore[union-attr]
                {
                    "file": filename,
                    "species": species,
                    "sha256": _sha256(staged_path),
                    "levels": len(data.levels),
                    "rotational_template_size": len(template),
                    "radiative_transitions": len(data.radiative),
                    "radiative_integrity": radiative_audit,
                    "radiative_inheritance": {
                        **inheritance,
                        "upstream_path": str(upstream_path.relative_to(root)),
                        "upstream_sha256": _sha256(upstream_path),
                        "upstream_manifest_path": str(
                            upstream_manifest_path.relative_to(root)
                        ),
                        "upstream_manifest_sha256": _sha256(
                            upstream_manifest_path
                        ),
                        "upstream_radiative_sidecar_path": str(
                            upstream_sidecar.relative_to(root)
                        ),
                        "upstream_radiative_sidecar_sha256": _sha256(
                            upstream_sidecar
                        ),
                    },
                    "collision_inheritance": {
                        **collision_inheritance,
                        "upstream_path": str(upstream_path.relative_to(root)),
                        "upstream_sha256": _sha256(upstream_path),
                        "upstream_manifest_path": str(
                            upstream_manifest_path.relative_to(root)
                        ),
                        "upstream_manifest_sha256": _sha256(
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
                    "structural_validation": "passed",
                    "scientific_limit": "CO stretch absent; not a complete Td=400-500 K input",
                }
            )
        output_dir.mkdir(parents=True, exist_ok=True)
        for staged_path in stage.iterdir():
            staged_path.replace(output_dir / staged_path.name)
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
