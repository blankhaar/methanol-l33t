#!/usr/bin/env python3
"""Reproduce and extract the frozen L33T v3 molecular data, with explicit inputs."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE/'src'))


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')


def lock():
    return json.loads((PACKAGE/'provenance/inputs.lock.json').read_text())


def safe_relative(value):
    p = Path(value)
    if p.is_absolute() or '..' in p.parts:
        raise ValueError(f'unsafe manifest path: {value}')
    return p


def check_inputs(root):
    failures = []
    for record in lock()['files']:
        p = root/safe_relative(record['path'])
        if not p.is_file():
            failures.append(f"missing: {record['path']}")
        elif p.stat().st_size != record['bytes'] or sha(p) != record['sha256']:
            failures.append(f"checksum mismatch: {record['path']}")
    if failures:
        raise ValueError('\n'.join(failures))
    return len(lock()['files'])


def import_inputs(args):
    """Import only exact named files; no recursive copy, download, or source edits."""
    if args.dest.exists():
        raise FileExistsError('input destination must be new')
    roots = {'resource': args.resource_root.resolve(), 'repair': args.repair_root.resolve()}
    selected = []
    for record in lock()['files']:
        src = roots[record['source_root']]/safe_relative(record['source_path'])
        if not src.is_file() or sha(src) != record['sha256']:
            raise ValueError(f"source absent or changed: {record['source_root']}:{record['source_path']}")
        selected.append((src, safe_relative(record['path'])))
    args.dest.mkdir(parents=True)
    for src, relative in selected:
        dst = args.dest/relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
    count = check_inputs(args.dest)
    print(f'Imported and verified {count} files. Third-party redistribution remains separate.')


def normalize(inputs):
    check_inputs(inputs)
    # Compatibility ledger for the original reader; derived from the same lock,
    # so this is not independent evidence beyond check_inputs.
    sums = ''.join(f"{r['sha256']}  {r['path']}\n" for r in lock()['files']
                   if r['path'].startswith('data/raw/'))
    (inputs/'data/raw/SHA256SUMS').write_text(sums)
    master = inputs/'data/normalized/master/collisions.csv'
    expected = lock()['normalized_collision_sha256']
    if master.exists():
        if sha(master) != expected:
            raise ValueError('existing normalized collision master differs from the frozen reference')
        print('Normalized collisions already match the frozen reference.', flush=True)
        return
    from methanol_lamda.normalize import normalize_basecol_collisions_master
    report = normalize_basecol_collisions_master(inputs)
    if not report or sha(master) != expected:
        raise ValueError('normalized BASECOL data do not reproduce the frozen collision master')
    from review_support import portable_paths
    dump(inputs/'NORMALIZATION.json', portable_paths(report, {inputs.resolve(): '@inputs'}))
    print(f"Normalized {report['collision_rows']} collision rows; exact checksum passes.", flush=True)


def run(script, *args):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1',
               PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')
    # Keep invariant checks enabled even when the caller uses optimized Python.
    env.pop('PYTHONOPTIMIZE', None)
    command = [sys.executable, '-B', str(PACKAGE/'scripts'/script), *map(str, args)]
    print('Running '+script, flush=True)
    subprocess.run(command, env=env, cwd=PACKAGE, check=True)


def spectroscopy(inputs, output, isotope):
    check_inputs(inputs)
    run('spectroscopy.py', '--inputs', inputs.resolve(), '--output', output.resolve(),
        '--isotope', 'A'+isotope, '--modes', 'PERMANENT', 'VCO_A', 'VCO_C')


def build(inputs, spec12, spec13, output):
    normalize(inputs)
    if output.exists():
        raise FileExistsError('product destination must be new')
    output.mkdir(parents=True)
    for iso, spec in [('12', spec12), ('13', spec13)]:
        run('build_l33t.py', '--inputs', inputs.resolve(), '--isotope', iso,
            '--spectroscopy', spec.resolve(), '--main-spectroscopy', spec12.resolve(),
            '--vco-c-spectroscopy', spec.resolve(), '--output', (output/iso).resolve())
    records = verify(output)
    for iso in ('12', '13'):
        path = output/iso/'BUILD-MANIFEST.json'
        manifest = json.loads(path.read_text())
        manifest['status'] = 'verified-exact-v3-DAT-and-expanded-provenance'
        manifest['verification'] = dict(structural='passed', DAT='exact frozen v3 bytes',
                                       sidecars='exact decompressed content')
        dump(path, manifest)
    return records


def find_product(root, name):
    matches = [p for p in (root/name, root/('13' if '13ch3oh' in name else '12')/name) if p.is_file()]
    if len(matches) != 1:
        raise ValueError(f'expected one unambiguous product {name}; found {len(matches)}')
    return matches[0]


def verify(root, archive_bytes=False):
    from validate_products import inspect_lamda
    from review_support import expanded_hash
    records = []
    for reference in json.loads((PACKAGE/'provenance/products.json').read_text())['products']:
        path = find_product(root, reference['file'])
        digest = sha(path)
        if digest != reference['sha256']:
            raise ValueError(f'product differs from frozen v3 bytes: {path.name}; retained for diagnosis')
        for name, info in reference['sidecars'].items():
            sidecar = path.parent/name
            if not sidecar.is_file() or expanded_hash(sidecar) != {
                k: info[k] for k in ('uncompressed_sha256', 'uncompressed_bytes')
            }:
                raise ValueError(f'missing or changed provenance sidecar: {name}')
            if archive_bytes and sha(sidecar) != info['sha256']:
                raise ValueError(f'compressed archive bytes differ: {name}')
        result = inspect_lamda(path)
        if result['levels'] != reference['levels'] or result['radiative'] != reference['radiative']:
            raise ValueError(f'reference count mismatch: {path.name}')
        records.append(dict(file=path.name, sha256=digest, bytes=path.stat().st_size, **result))
        print(f'PASS {path.name}: exact v3 DAT, expanded sidecars and structural audit'
              + ('; compressed bytes exact' if archive_bytes else ''), flush=True)
    return records


def extract(args):
    """Copy a chosen isotope/species subset, with provenance, into a new directory."""
    verify(args.products)
    if args.output.exists():
        raise FileExistsError('extraction destination must be new')
    selected = []
    for ref in json.loads((PACKAGE/'provenance/products.json').read_text())['products']:
        name = ref['file']
        if args.isotope != 'both' and f'-{args.isotope}ch3oh-' not in name:
            continue
        if args.species != 'both' and not name.startswith(args.species.lower()+'-'):
            continue
        src = find_product(args.products, name)
        tag = name.split('-l33t-')[0]
        siblings = [src.parent/f'{tag}.{suffix}.csv.gz' for suffix in
                    ('levels', 'radiative-ledger', 'collision-provenance')]
        if not all(p.is_file() for p in siblings):
            raise ValueError('extraction requires all three aligned provenance sidecars')
        selected.append((src, siblings))
    args.output.mkdir(parents=True)
    inventory = []
    for src, siblings in selected:
        dst = args.output/(src.name+'.gz' if args.gzip else src.name)
        if args.gzip:
            with src.open('rb') as source, dst.open('xb') as raw:
                with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as target:
                    shutil.copyfileobj(source, target)
        else:
            shutil.copyfile(src, dst)
        inventory.append(dict(file=dst.name, sha256=sha(dst), uncompressed_sha256=sha(src)))
        for p in siblings:
            shutil.copyfile(p, args.output/p.name)
            inventory.append(dict(file=p.name, sha256=sha(args.output/p.name)))
    dump(args.output/'EXTRACTION.json', dict(version='l33t-v3', release_version='1.0', files=inventory,
        status='adopted for numerical experiments, with model approximations documented in the v1.0 companion'))
    print(f'Extracted {len(selected)} LAMDA files with provenance to {args.output}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('inputs', help='import exact inputs from existing authorized local trees')
    p.add_argument('--resource-root', type=Path, required=True)
    p.add_argument('--repair-root', type=Path, required=True)
    p.add_argument('--dest', type=Path, required=True)
    p = sub.add_parser('check-inputs')
    p.add_argument('--inputs', type=Path, required=True)
    p = sub.add_parser('normalize')
    p.add_argument('--inputs', type=Path, required=True)
    p = sub.add_parser('spectroscopy')
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--isotope', choices=('12', '13'), required=True)
    p = sub.add_parser('build')
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--spectroscopy12', type=Path, required=True)
    p.add_argument('--spectroscopy13', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('rebuild', help='run normalization, both Hamiltonians, radiation, collisions and validation')
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('verify')
    p.add_argument('--products', type=Path, required=True)
    p.add_argument('--archive-bytes', action='store_true', help='also require original gzip bytes')
    p = sub.add_parser('coverage', help='measure per-collider endpoint coverage')
    p.add_argument('--products', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p = sub.add_parser('compare', help='label-keyed numerical DAT comparison; does not change release locks')
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--products', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p = sub.add_parser('extract')
    p.add_argument('--products', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--isotope', choices=('12','13','both'), default='both')
    p.add_argument('--species', choices=('A','E','both'), default='both')
    p.add_argument('--gzip', action='store_true')
    args = parser.parse_args()
    if args.command == 'inputs':
        import_inputs(args)
    elif args.command == 'check-inputs':
        print(f'PASS: {check_inputs(args.inputs)} exact inputs')
    elif args.command == 'normalize':
        normalize(args.inputs)
    elif args.command == 'spectroscopy':
        spectroscopy(args.inputs, args.output, args.isotope)
    elif args.command == 'build':
        build(args.inputs, args.spectroscopy12, args.spectroscopy13, args.output)
    elif args.command == 'verify':
        verify(args.products, archive_bytes=args.archive_bytes)
    elif args.command == 'coverage':
        from review_support import coverage
        verify(args.products)
        records = [coverage(find_product(args.products, r['file'])) for r in
                   json.loads((PACKAGE/'provenance/products.json').read_text())['products']]
        dump(args.report, dict(products=records))
    elif args.command == 'compare':
        from review_support import compare_dat
        from validate_products import inspect_lamda
        records = []
        for r in json.loads((PACKAGE/'provenance/products.json').read_text())['products']:
            baseline = find_product(args.reference, r['file'])
            candidate = find_product(args.products, r['file'])
            inspect_lamda(baseline); inspect_lamda(candidate)
            records.append(compare_dat(baseline, candidate))
            print(r['file'], records[-1]['passed'], flush=True)
        dump(args.report, dict(products=records, passed=all(r['passed'] for r in records),
                               scope='Numerical comparison only; not release identity or physical validation.'))
        if not all(r['passed'] for r in records):
            raise SystemExit(1)
    elif args.command == 'extract':
        extract(args)
    elif args.command == 'rebuild':
        if args.output.exists():
            raise FileExistsError('rebuild output must be new')
        args.output.mkdir(parents=True)
        start = time.monotonic()
        normalize(args.inputs)
        for iso in ('12', '13'):
            spectroscopy(args.inputs, args.output/f'spectroscopy-{iso}', iso)
        records = build(args.inputs, args.output/'spectroscopy-12', args.output/'spectroscopy-13', args.output/'products')
        dump(args.output/'REPRODUCTION.json', dict(status='exact-v3-DAT-and-expanded-provenance-pass',
            release_version='1.0', product_version='l33t-v3',
            wall_seconds=time.monotonic()-start, python=sys.version,
            verification='exact DAT and decompressed sidecars; structural audit', products=records))


if __name__ == '__main__':
    main()
