"""Portable content checks, measured coverage and label-keyed DAT comparison."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import tempfile


def portable_paths(value, roots):
    """Replace known local roots in shareable run records, retaining role and hashes."""
    if isinstance(value, dict):
        return {k: portable_paths(v, roots) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [portable_paths(v, roots) for v in value]
    if isinstance(value, str):
        for root, role in sorted(roots.items(), key=lambda x: len(str(x[0])), reverse=True):
            value = value.replace(str(root), role)
    return value


def expanded_hash(path):
    digest = hashlib.sha256()
    size = 0
    with gzip.open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    return dict(uncompressed_sha256=digest.hexdigest(), uncompressed_bytes=size)


def coverage(path):
    """Count actual endpoints, including zero-only rows; do not infer completeness."""
    tag = path.name.split('-l33t-')[0]
    with gzip.open(path.parent / (tag + '.levels.csv.gz'), 'rt') as f:
        levels = {r['label']: r for r in csv.DictReader(f)}
    result = {}
    with gzip.open(path.parent / (tag + '.collision-provenance.csv.gz'), 'rt') as f:
        for row in csv.DictReader(f):
            c = row['collider']
            item = result.setdefault(c, dict(rows=0, labels=set(), extrapolated_rows=0))
            item['rows'] += 1
            item['labels'].update((row['upper'], row['lower']))
            item['extrapolated_rows'] += row['status'] == 'l33t-JK-extrapolated'
    for item in result.values():
        touched = item.pop('labels')
        if not touched <= levels.keys():
            raise ValueError('collision provenance references unknown levels')
        item.update(levels_touched=len(touched), levels_without_rows=len(levels)-len(touched))
        item['by_manifold'] = {}
        for label, level in levels.items():
            man = level['manifold']
            counts = item['by_manifold'].setdefault(man, dict(total=0, touched=0))
            counts['total'] += 1
            counts['touched'] += label in touched
    return dict(file=path.name, levels=len(levels), colliders=result,
                definition='Endpoints present in rows; a present zero rate is not positive coupling.')


def dat_records(path):
    """Yield identities and numerical fields, independent of transition row order."""
    with open(path) as f:
        def line():
            for raw in f:
                if raw.strip() and not raw.lstrip().startswith(('!', '#')):
                    return raw.split()
            raise ValueError('truncated LAMDA file')
        line()  # descriptive title is not a numerical field
        yield 'mass', 'mass', [float(line()[0])]
        labels = {}
        for _ in range(int(line()[0])):
            r = line()
            labels[int(r[0])] = r[3]
            yield 'level', r[3], [float(r[1]), float(r[2])]
        for _ in range(int(line()[0])):
            r = line()
            key = json.dumps([labels[int(r[1])], labels[int(r[2])]])
            yield 'radiative', key, list(map(float, r[3:]))
        for _ in range(int(line()[0])):
            collider = line()[0]
            n = int(line()[0]); nt = int(line()[0])
            grid = list(map(float, line()))
            if len(grid) != nt:
                raise ValueError('temperature width mismatch')
            yield 'grid', collider, grid
            for _ in range(n):
                r = line()
                key = json.dumps([collider, labels[int(r[1])], labels[int(r[2])]])
                yield 'collision', key, list(map(float, r[3:]))
        if any(r.strip() and not r.lstrip().startswith(('!', '#')) for r in f):
            raise ValueError('unexpected trailing records')


def tolerances(kind, index):
    # Zeros in strengths/rates remain exact. These are numerical, not physical errors.
    if kind in ('mass', 'grid') or (kind == 'level' and index == 1):
        return 0.0, 0.0
    if kind == 'level':
        return 1e-10, 1e-8  # cm^-1
    if kind == 'radiative' and index in (1, 2):
        return 1e-9, 1e-9  # GHz, K
    return 1e-7, 0.0  # Einstein A and collisional rates


def compare_dat(reference, candidate):
    """Disk-backed matching keeps memory bounded for millions of collision rows."""
    result = dict(reference=reference.name, candidate=candidate.name, matched={},
                  missing=0, extra=0, numerical_failures=0, examples=[], maxima={})
    with tempfile.TemporaryDirectory(prefix='l33t-compare-') as temp:
        db = sqlite3.connect(str(Path(temp) / 'rows.sqlite'))
        try:
            db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
            db.execute('CREATE TABLE rows(kind TEXT, identity TEXT, vals TEXT, seen INTEGER DEFAULT 0, PRIMARY KEY(kind,identity)) WITHOUT ROWID')
            db.executemany('INSERT INTO rows(kind,identity,vals) VALUES (?,?,?)',
                           ((k, key, json.dumps(v, allow_nan=False)) for k, key, v in dat_records(reference)))
            for kind, key, actual in dat_records(candidate):
                row = db.execute('SELECT vals,seen FROM rows WHERE kind=? AND identity=?', (kind,key)).fetchone()
                if row is None:
                    result['extra'] += 1
                    if len(result['examples']) < 10: result['examples'].append(['extra',kind,key])
                    continue
                if row[1]:
                    raise ValueError('duplicate candidate identity: '+key)
                db.execute('UPDATE rows SET seen=1 WHERE kind=? AND identity=?', (kind,key))
                expected = json.loads(row[0])
                result['matched'][kind] = result['matched'].get(kind,0)+1
                failed = len(actual) != len(expected)
                for i, (a,b) in enumerate(zip(actual,expected)):
                    rt, at = tolerances(kind,i)
                    finite = math.isfinite(a) and math.isfinite(b)
                    error = abs(a-b) if finite else math.inf
                    relative = error/max(abs(a),abs(b)) if max(abs(a),abs(b)) else 0.0
                    stat = result['maxima'].setdefault(f'{kind}:{i}', dict(absolute=0.0,relative=0.0,rtol=rt,atol=at))
                    stat['absolute'] = max(stat['absolute'],error)
                    stat['relative'] = max(stat['relative'],relative)
                    failed |= not finite or error > at + rt*abs(b)
                if failed:
                    result['numerical_failures'] += 1
                    if len(result['examples']) < 10: result['examples'].append(['numerical',kind,key])
            result['missing'] = db.execute('SELECT count(*) FROM rows WHERE seen=0').fetchone()[0]
        finally:
            db.close()
    result['passed'] = not any(result[k] for k in ('missing','extra','numerical_failures'))
    result['scope'] = 'Numerical DAT fields and exact identities; sidecar provenance is checked separately by verify.'
    return result
