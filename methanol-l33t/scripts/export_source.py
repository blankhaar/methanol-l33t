#!/usr/bin/env python3
"""Create a source ZIP, optionally with all locked inputs for local review."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
DIRS=('src','fortran','scripts','tests','docs','provenance')
TOP=('README.md','LICENSE','LICENSE-NOTICE.md','CITATION.cff','requirements.txt','.gitignore')


def source_files():
    files=[ROOT/x for x in TOP]
    for directory in DIRS:
        for p in (ROOT/directory).rglob('*'):
            if not p.is_file() or '__pycache__' in p.parts or 'build' in p.parts:
                continue
            if p.suffix in ('.pyc','.pyo','.aux','.log','.out','.bbl','.blg','.fls','.fdb_latexmk'):
                continue
            files.append(p)
    for p in files:
        if p.is_symlink() or not p.is_file():
            raise ValueError(f'unsupported/missing source member: {p}')
        if p.stat().st_size > 10*1024**2:
            raise ValueError(f'unexpectedly large source member: {p}')
        if p.suffix in ('.py','.f','.f90','.json','.md','.tex','.txt','.yml','.cff'):
            txt=p.read_text()
            if ('/'+'Users'+'/') in txt or ('/'+'home'+'/') in txt:
                raise ValueError(f'personal absolute path in source member: {p}')
    return sorted(set(files))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--inputs',type=Path,help='include precisely the 97 locked inputs in a local review bundle')
    args=ap.parse_args()
    files=source_files()
    input_files=[]
    if args.inputs:
        from l33t import check_inputs, lock, safe_relative
        check_inputs(args.inputs)
        input_files=[(args.inputs/safe_relative(r['path']), 'inputs/'+r['path']) for r in lock()['files']]
        if any(p.is_symlink() for p,_ in input_files):
            raise ValueError('input bundle requires regular files')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    manifest=[]
    with zipfile.ZipFile(args.output,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p, relative in [(p,str(p.relative_to(ROOT))) for p in files]+input_files:
            data=p.read_bytes();name='methanol-l33t/'+relative
            info=zipfile.ZipInfo(name,date_time=(2026,11,20,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            z.writestr(info,data)
            manifest.append(dict(path=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
        info=zipfile.ZipInfo('methanol-l33t/SOURCE-MANIFEST.json',date_time=(2026,11,20,0,0,0))
        info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
        z.writestr(info,json.dumps(dict(version='1.0',product='l33t-v3',
            locked_inputs_included=bool(args.inputs), scope=('private reproduction bundle' if args.inputs else 'public source distribution'),files=manifest),indent=2)+'\n')
    with zipfile.ZipFile(args.output) as z:
        if z.testzip() is not None:
            raise ValueError('ZIP integrity failure')
    print(json.dumps(dict(file=str(args.output),members=len(manifest)+1,bytes=args.output.stat().st_size,
                         sha256=hashlib.sha256(args.output.read_bytes()).hexdigest())))


if __name__=='__main__':
    main()
