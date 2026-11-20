"""Source-adjudicated 2026-10-16 vibrational candidates; originals are immutable.

This revision is intentionally separate from the historical builders. It binds
its parent bytes, calibrates emitted CO absorption, reconstructs corrected rock
terms, and adds a theoretical parallel band. See the accompanying source profile
and scientific report for the limits of the uncoupled effective operator.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import statistics
import tempfile

from .ch3_rock import (
    C13_ROCK_MODEL, RockKey, SourceLine, _c13_origins, _collisions,
    _ground_key, _main_rock_source, _merge_levels, _radiative,
)
from .co_stretch import (
    CM1_TO_MHZ, HITRAN_E_GROUND_ZERO_CM1, ISOTOPE_CONSTANTS,
    MAIN_HITRAN_RELATIVE_PATH, MAIN_LEVELS, PRODUCT_FLOOR_SMU2,
    VCO_TRANSITION_MOMENT_DEBYE, COStretchKey, _main_ground_energies,
    _parse_hitran_local, _resolve_hitran_ground_identity,
    _resolve_hitran_upper_identity, _run_transition_operator,
)
from .collision_integrity import assert_parent_collision_inheritance
from .hamiltonian import sha256_file
from .independent import read_summary
from .lamda import read_lamda, validate_lamda, write_lamda
from .models import LamdaData, Level, RadiativeTransition
from .physics import CM1_TO_K, SMU2_TO_A, einstein_a_from_smu2
from .quantum import LevelKey
from .radiative_integrity import audit_radiative_integrity

PROFILE = Path('config/vibrational-revision-2026-10-16')
DEFAULT_OUTPUT = Path('outputs/vibrational-revised-2026-10-16')
C_CM_S = 2.99792458e10


def _json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def absorption(a: float, nu: float, elo: float, gu: float,
               temperature: float, q_spin: float) -> float:
    """HITRAN-style integrated absorption in cm/molecule."""
    return (a * gu * math.exp(-CM1_TO_K * elo / temperature)
            * -math.expm1(-CM1_TO_K * nu / temperature)
            / (8 * math.pi * C_CM_S * nu**2 * q_spin))


def a_from_absorption(s: float, nu: float, elo: float, gu: float,
                      temperature: float, q_spin: float) -> float:
    return s / absorption(1.0, nu, elo, gu, temperature, q_spin)


def _parent_path(root: Path, isotope: str, species: str) -> Path:
    home = root if isotope == '12' else root / 'isotopologues/13ch3oh'
    stem = 'ch3oh' if isotope == '12' else '13ch3oh'
    return home / f'outputs/co-stretch/{species.lower()}-{stem}-vt01234-vco1-intraband-copy.dat'


def _calibrate_co(root: Path, profile: dict) -> dict:
    settings = profile['co_calibration']
    raw = (root / MAIN_HITRAN_RELATIVE_PATH).read_text().splitlines()
    ground = _main_ground_energies(root / MAIN_LEVELS)
    records = []
    for species in ('A', 'E'):
        data = read_lamda(_parent_path(root, '12', species))
        levels = {v.index: v for v in data.levels}
        rows = {(levels[r.upper].quantum_label, levels[r.lower].quantum_label): r
                for r in data.radiative}
        by_label = {v.quantum_label: v for v in data.levels}
        for source in settings['rows']:
            if source['wavenumber_cm1'] == settings['excluded_feature_cm1']:
                continue
            line = raw[source['raw_row'] - 1]
            u0, l0 = _parse_hitran_local(line[97:112]), _parse_hitran_local(line[112:127])
            if u0.species != species:
                continue
            if not (abs(float(line[3:15]) - source['wavenumber_cm1']) < 3e-4):
                raise RuntimeError("Invariant failed in vibrational_revision.py: abs(float(line[3:15]) - source['wavenumber_cm1']) < 3e-4")
            l, _, _ = _resolve_hitran_ground_identity(l0, float(line[45:55]), ground)
            if l is None:
                raise ValueError('calibration lower state is unmapped')
            u = _resolve_hitran_upper_identity(u0, l0, l)
            ul, ll = COStretchKey.from_ground(u).encode(), l.encode()
            row = rows[(ul, ll)]
            elo = by_label[ll].energy_cm1 + (HITRAN_E_GROUND_ZERO_CM1 if species == 'E' else 0)
            # Use the published feature frequency in the laboratory comparison.
            predicted = absorption(row.a_s, source['wavenumber_cm1'], elo, 2*u.j+1,
                                   settings['temperature_K'], settings['Q_over_common_spin'])
            predicted *= 2.68675e19 * 273.15 / settings['temperature_K']
            records.append(dict(source, upper_key=ul, lower_key=ll,
                                original_predicted=predicted,
                                original_over_xu=predicted/source['xu_cm2_inverse_atm_inverse']))
        del data
    if len(records) != 13:
        raise ValueError('expected 13 isolated laboratory calibration features')
    scale = 1 / statistics.median(r['original_over_xu'] for r in records)
    for row in records:
        row['revised_predicted'] = row['original_predicted'] * scale
        row['revised_over_xu'] = row['original_over_xu'] * scale
        row['revised_over_earlier'] = row['revised_predicted'] / row['earlier_cm2_inverse_atm_inverse']
    return dict(rate_scale=scale, effective_moment_debye=.27*math.sqrt(scale),
                rows=sorted(records, key=lambda r: r['wavenumber_cm1']),
                status='in-sample robust normalization, not independent accuracy certification')


def _corrected_rock_source(root: Path, profile: dict):
    raw = (root / MAIN_HITRAN_RELATIVE_PATH).read_text().splitlines()
    ground = _main_ground_energies(root / MAIN_LEVELS)
    old_terms, old_rows, old_audit = _main_rock_source(root / MAIN_HITRAN_RELATIVE_PATH,
                                                     root / MAIN_LEVELS)
    # Preserve the previously screened A and E-5 series. E+4 belongs to E+3.
    terms = {k: v for k, v in old_terms.items() if not (k.species == 'E' and k.k == 4)}
    selected = {r.source_line: r for r in old_rows if not (r.upper.species == 'E' and r.upper.k == 4)}
    groups = defaultdict(list)
    quarantined = set(profile['source_assignment']['quarantine_raw_rows'])
    for number, line in enumerate(raw, 1):
        if not line.startswith('391') or line[67:82].strip() != 'V7':
            continue
        u0, l0 = _parse_hitran_local(line[97:112]), _parse_hitran_local(line[112:127])
        if not (u0.species == 'E' and u0.k == 4 and l0.k == 3):
            continue
        if number in quarantined:
            continue
        u, l = LevelKey('E', 0, u0.j, 3), LevelKey('E', 0, l0.j, 4)
        if u not in ground or l not in ground:
            continue
        nu = float(line[3:15])
        term = ground[l] + nu
        selected[number] = SourceLine(
            u, l, u0.j, float(line[146:153]), float(line[153:160]), nu,
            float(line[25:35]), number, term,
            'Xu2004_Table3_K3_from_K4; frozen_corrected_lower_energy; multibranch_term',
            hashlib.sha256(line.encode('ascii')).hexdigest(),
            ground[l] + HITRAN_E_GROUND_ZERO_CM1 - float(line[45:55]), None)
        groups[u].append(term)
    spreads = {k.encode(): max(v)-min(v) for k, v in groups.items()}
    if not groups or max(spreads.values()) > .003:
        raise ValueError('corrected E+3 source fails combination-difference gate')
    terms.update({k: statistics.median(v) for k, v in groups.items()})
    direct, ledger = [], []
    for number, source in sorted(selected.items()):
        line = raw[number-1]
        s = float(line[15:25])
        # Even unchanged rows are consistently reconverted from source S.
        elo = ground[source.lower] + (HITRAN_E_GROUND_ZERO_CM1 if source.lower.species == 'E' else 0)
        a = a_from_absorption(s, source.wavenumber_cm1, elo, 2*source.upper.j+1,
                             296, profile['empirical_intensity_Q296_over_common_spin'])
        fixed = replace(source, source_upper_j=source.upper.j, einstein_a_s=a)
        status = 'accepted_empirical_strength' if s > 1.0001e-26 else 'position_only_placeholder_intensity'
        ledger.append(dict(raw_row=number, raw_sha256=source.raw_record_sha256,
                           upper_key=source.upper.encode(), lower_key=source.lower.encode(),
                           basis=source.identity_basis, status=status,
                           raw_s_cm_per_molecule=s, raw_a_s=float(line[25:35]),
                           raw_lower_energy_cm1=float(line[45:55]),
                           corrected_lower_energy_absolute_cm1=elo,
                           corrected_a_s=a, corrected_smu2_debye2=fixed.smu2_debye2,
                           wavenumber_cm1=source.wavenumber_cm1,
                           inferred_upper_energy_cm1=source.inferred_upper_energy_cm1,
                           adopted_upper_energy_cm1=terms[source.upper]))
        if status == 'accepted_empirical_strength':
            direct.append(fixed)
    selected_numbers = set(selected)
    for number, line in enumerate(raw, 1):
        if line.startswith('391') and line[67:82].strip() == 'V7' and number not in selected_numbers:
            ledger.append(dict(raw_row=number,
                               raw_sha256=hashlib.sha256(line.encode('ascii')).hexdigest(),
                               status=('quarantined_additional_J_ambiguity' if number in quarantined
                                       else 'outside_basis_or_rejected_historical_term_cluster')))
    return terms, direct, sorted(ledger, key=lambda x: x['raw_row']), dict(
        corrected_E3_term_count=len(groups), corrected_E3_branch_spreads_cm1=spreads,
        max_corrected_E3_branch_spread_cm1=max(spreads.values()),
        quarantined_rows=sorted(quarantined), historical_selection_audit=old_audit)


def _origin(species: str, k: int, isotope: str, main: dict, c13: dict):
    """K-resolved source or explicit isotope/family proxy; never a global median."""
    table = main[species] if isotope == '12' else c13[species]
    if k in table:
        return table[k], ('Xu2004_Table3_origin' if isotope == '12' else 'PredoiCross1997_Table4.1_origin')
    if isotope == '12' and k in c13[species]:
        shared = main[species].keys() & c13[species].keys()
        shift = statistics.median(main[species][q]-c13[species][q] for q in shared)
        return c13[species][k]+shift, f'isotope_pattern_proxy_median_12minus13_shift_{shift:.9g}'
    # Signed K modulo 3 retains the torsional family; tie break is deterministic.
    candidates = [q for q in table if (q-k) % 3 == 0]
    if not candidates:
        raise ValueError(f'no source origin in K family for {species} {k}')
    nearest = min(candidates, key=lambda q: (abs(q-k), abs(q), q))
    return table[nearest], f'constant_nearest_K_mod3_origin_proxy_from_K{nearest:+d}'


def _levels(base: LamdaData, isotope: str, terms: dict, main: dict, c13: dict):
    ground = {_ground_key(v.quantum_label): v for v in base.levels
              if _ground_key(v.quantum_label) is not None}
    offsets = defaultdict(list)
    for key, v in ground.items():
        if key in terms:
            offsets[(key.species, key.k, key.parity)].append((key.j*(key.j+1), terms[key]-v.energy_cm1))
    result, ledger = [], []
    for key, v in ground.items():
        origin, origin_basis = _origin(key.species, key.k, isotope, main, c13)
        if isotope == '12' and key in terms:
            energy, basis = terms[key], 'source_selected_direct_term'
        elif isotope == '12' and offsets[(key.species, key.k, key.parity)]:
            points = sorted([(0, origin)] + offsets[(key.species, key.k, key.parity)])
            x = key.j * (key.j+1)
            right = next((i for i, p in enumerate(points) if p[0] >= x), None)
            if right is None:
                offset = points[-1][1]
                basis = 'constant_offset_above_last_measured_J'
            elif right == 0:
                offset = points[0][1]
                basis = 'published_J0_origin'
            else:
                (x0, y0), (x1, y1) = points[right-1:right+1]
                offset = y0+(y1-y0)*(x-x0)/(x1-x0)
                basis = 'offset_interpolation_in_JJ1_between_origin_and_terms'
            energy = v.energy_cm1 + offset
        else:
            energy, basis = v.energy_cm1+origin, 'ground_rotational_spacing_plus_K_origin'
        label = RockKey.from_ground(key).encode()
        result.append(Level(0, energy, v.weight, label))
        ledger.append(dict(quantum_label=label, energy_cm1=energy,
                           ground_energy_cm1=v.energy_cm1, offset_cm1=energy-v.energy_cm1,
                           basis=basis, origin_cm1=origin, origin_basis=origin_basis))
    return result, ledger


def _scaled_co(base: LamdaData, scale: float):
    levels = {v.index: v for v in base.levels}
    rows, counts = [], Counter()
    for r in base.radiative:
        u, l = levels[r.upper], levels[r.lower]
        if ':vCO1:' in u.quantum_label and _ground_key(l.quantum_label) is not None:
            strength = r.a_s * u.weight / (SMU2_TO_A * (r.frequency_ghz*1000)**3) * scale
            if strength < PRODUCT_FLOOR_SMU2:
                counts['CO_dropped_below_revised_strength_floor'] += 1
                continue
            rows.append(replace(r, index=len(rows)+1, a_s=r.a_s*scale))
            counts['CO_crossband_rescaled'] += 1
        else:
            rows.append(replace(r, index=len(rows)+1))
            counts['parent_other_radiative_unchanged'] += 1
    return LamdaData(base.molecule, base.mass_amu, base.levels, rows, base.collision_partners), counts


def _parallel(data: LamdaData, strengths: list, mu2: float):
    by_label = {v.quantum_label: v for v in data.levels}
    existing = {(r.upper, r.lower) for r in data.radiative}
    new, ledger = [], []
    duplicates = 0
    for strength in strengths:
        s = strength.smu2_debye2 * mu2 / VCO_TRANSITION_MOMENT_DEBYE**2
        if s < PRODUCT_FLOOR_SMU2:
            continue
        pairs = [(strength.left, strength.right)]
        if strength.left != strength.right:
            pairs.append((strength.right, strength.left))
        for ru, gl in pairs:
            ul, ll = RockKey.from_ground(ru).encode(), gl.encode()
            if ul not in by_label or ll not in by_label:
                continue
            u, l = by_label[ul], by_label[ll]
            if (u.index, l.index) in existing:
                duplicates += 1  # empirical row supersedes the approximate operator
                continue
            nu = (u.energy_cm1-l.energy_cm1)*CM1_TO_MHZ
            a = einstein_a_from_smu2(nu, s, u.weight)
            new.append(RadiativeTransition(0, u.index, l.index, a, nu/1000))
            existing.add((u.index, l.index))
            ledger.append(dict(upper_key=ul, lower_key=ll, smu2_debye2=s,
                               basis='Sunaga2026_total_vibrational_strength_to_pure_a_BELGI_proxy'))
    data.radiative = [replace(r, index=i) for i, r in enumerate(
        sorted(data.radiative+new, key=lambda r: (r.upper, r.lower)), 1)]
    return ledger, duplicates


def build_vibrational_revision(project_root: str | Path, *, output_dir: str | Path | None = None,
                               compiler: str | Path | None = None) -> dict:
    root = Path(project_root).resolve()
    out = Path(output_dir) if output_dir is not None else root / DEFAULT_OUTPUT
    if not out.is_absolute():
        out = root / out
    out = out.resolve()
    if out.exists():
        raise FileExistsError(f'refusing to replace any existing output directory: {out}')
    profile = json.loads((root/PROFILE/'sources.json').read_text())
    inputs = json.loads((root/PROFILE/'inputs.json').read_text())
    for path, expected in inputs.items():
        if sha256_file(root/path) != expected:
            raise ValueError(f'bound input changed: {path}')
    calibration = _calibrate_co(root, profile)
    terms, direct, source_ledger, source_audit = _corrected_rock_source(root, profile)
    c13, _ = _c13_origins(root/C13_ROCK_MODEL)
    main = {s: {int(k): v for k, v in rows.items()}
            for s, rows in profile['main_parallel_rock_origins_cm1'].items()}
    manifest = dict(revision=profile['revision'], status=profile['status'], inputs=inputs,
                    profile_sha256=sha256_file(root/PROFILE/'sources.json'),
                    co_calibration=calibration, rock_source_audit=source_audit, files=[])
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.vibrational-revision-', dir=out.parent) as tmp:
        stage = Path(tmp)/'products'
        stage.mkdir()
        _json(stage/'source-corrections.json', source_ledger)
        for isotope in ('12', '13'):
            constants = (root/'data/raw/local-cache/xu2008-final-constants.txt'
                         if isotope == '12' else root/ISOTOPE_CONSTANTS)
            strengths, operator_summary = _run_transition_operator(root, constants, compiler=compiler)
            for species in ('A', 'E'):
                parent_path = _parent_path(root, isotope, species)
                parent = read_lamda(parent_path)
                base, co_counts = _scaled_co(parent, calibration['rate_scale'])
                rock, energy_ledger = _levels(base, isotope, terms, main, c13)
                levels, oldmap, labelmap, ties = _merge_levels(base, rock)
                rad, rad_counts = _radiative(base, levels, oldmap, labelmap, direct,
                                             species=species, isotopologue=f'{isotope}CH3OH')
                data = LamdaData(
                    f'{species}-{isotope}CH3OH vt=0..4+vCO1+vrock1 [revised 2026-10-16 uncoupled proxies]',
                    base.mass_amu, levels, rad, _collisions(base, levels, oldmap, labelmap))
                mu2 = (profile['rock_strength']['squared_transition_moment_au'][species]
                       * profile['rock_strength']['atomic_dipole_to_debye']**2)
                parallel, duplicates = _parallel(data, strengths, mu2)
                filename = f'{species.lower()}-{isotope}ch3oh-vt01234-vco1-vrock1-revised.dat'
                # The scaled parent is exact for this tier. The independent audit
                # separately proves that only CO crossband rates changed upstream.
                audit = audit_radiative_integrity(data, context=filename, parent=base)
                collision_audit = assert_parent_collision_inheritance(parent, data, context=filename)
                errors = validate_lamda(data)
                if errors:
                    raise ValueError(errors[:10])
                rock_indices = {v.index for v in levels if ':vrock1:' in v.quantum_label}
                decay = {r.upper for r in data.radiative if r.upper in rock_indices and r.lower not in rock_indices}
                if decay != rock_indices:
                    raise ValueError(f'{filename}: rock levels lack direct vibrational decay')
                write_lamda(stage/filename, data)
                check = read_summary(stage/filename)
                if check.level_count != 2541 or check.radiative_count != len(data.radiative):
                    raise ValueError('independent reader count mismatch')
                _json(stage/f'{species.lower()}-{isotope}-energy-ledger.json', energy_ledger)
                _json(stage/f'{species.lower()}-{isotope}-parallel-ledger.json', parallel)
                manifest['files'].append(dict(
                    file=filename, sha256=sha256_file(stage/filename),
                    parent=str(parent_path.relative_to(root)), levels=len(levels),
                    radiative_transitions=len(data.radiative),
                    collision_rows_by_partner={str(p.collider_id):len(p.transitions) for p in data.collision_partners},
                    rock_energy_counts=dict(Counter(r['basis'] for r in energy_ledger)),
                    rock_origin_counts=dict(Counter(r['origin_basis'] for r in energy_ledger)),
                    co_counts=dict(co_counts), inherited_and_empirical_radiative_counts=dict(rad_counts),
                    parallel_rows_added=len(parallel), empirical_duplicates_superseding_operator=duplicates,
                    rock_moment_debye=math.sqrt(mu2), rock_direct_decay_levels=len(decay),
                    radiative_audit=audit, collision_inheritance=collision_audit,
                    energy_ties=ties, operator_summary=operator_summary))
                del parent, base, data
        for path, expected in inputs.items():
            if sha256_file(root/path) != expected:
                raise ValueError(f'original input changed during build: {path}')
        manifest['originals_unchanged'] = True
        manifest['ledgers'] = {p.name: sha256_file(p) for p in sorted(stage.glob('*.json'))}
        manifest['implementation_sha256'] = sha256_file(Path(__file__))
        _json(stage/'manifest.json', manifest)
        stage.rename(out)
    return manifest
