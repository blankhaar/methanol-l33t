#!/usr/bin/env python3
"""Bounded native unpruned INTEN versus cached PERMANENT/VCO_A comparison."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys

PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE/'src'))
from methanol_lamda.hamiltonian import _compile_driver, _find_compiler, parse_final_constants, render_complete_cst_vector
from l33t import check_inputs, dump, sha
from review_support import portable_paths


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--jmax',type=int,default=12)
    args=ap.parse_args()
    if not 1<=args.jmax<=12:ap.error('bounded check requires 1 <= Jmax <= 12')
    inputs=args.inputs.resolve();out=args.output.resolve()
    check_inputs(inputs);out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    executables={};compilers={}
    for cached in (False,True):
        build=out/('cached' if cached else 'native');build.mkdir()
        exe,cmd,version,patches=_compile_driver(PACKAGE,inputs/'vendor/belgi',build,_find_compiler(None),
            driver_filename='belgi_frozen_strengths.f',executable_name='strengths',parity_repair=True,cached_dipoles=cached)
        executables[cached]=exe
        compilers[str(cached)]=dict(command=cmd,version=version,patches=patches)
    records=[]
    for iso,constant in [('12','data/raw/local-cache/xu2008-final-constants.txt'),
                         ('13','isotopologues/13ch3oh/constants/table2-2014-constants.txt')]:
        raw=render_complete_cst_vector(parse_final_constants(inputs/constant,require_xu_final=iso=='12'))
        for mode,vt in [('PERMANENT',4),('VCO_A',0)]:
            results={}
            for cached,exe in executables.items():
                run_dir=out/f'{iso}-{mode}-{cached}';run_dir.mkdir();(run_dir/'assignment').mkdir()
                run=subprocess.run([str(exe),str(args.jmax),str(args.jmax),str(vt),'1e-10',mode],
                    input=raw,text=True,capture_output=True,cwd=run_dir,check=True,timeout=1800,
                    env=dict(env,BELGI_REPAIR_AUDIT_DIR=str(run_dir/'assignment')))
                block=run.stdout.split('BEGIN_FROZEN_STRENGTHS_CSV\n',1)[1].split('END_FROZEN_STRENGTHS_CSV',1)[0]
                rows={}
                for row in csv.DictReader(io.StringIO(block)):
                    key=tuple(row[x].strip() for x in ('species','upper_vt','upper_j','upper_k','upper_parity','lower_vt','lower_j','lower_k','lower_parity'))
                    if key in rows:raise ValueError('duplicate native transition')
                    rows[key]=float(row['smu2_debye2'])
                results[cached]=rows
                dump(run_dir/'SUMMARY.json',dict(rows=len(rows),stdout_sha256=hashlib.sha256(run.stdout.encode()).hexdigest(),
                    summary=next(x for x in run.stdout.splitlines() if x.startswith('FROZEN_STRENGTH_SUMMARY,'))))
            a,b=results[False],results[True]
            common=a.keys()&b.keys()
            errors=[abs(a[k]-b[k])/max(abs(a[k]),abs(b[k])) for k in common]
            # Near-cancelled lines use a small absolute allowance as well.
            failures=sum(abs(a[k]-b[k])>1e-18+1e-10*abs(a[k]) for k in common)
            rec=dict(isotope=iso,mode=mode,rows=len(a),identities_equal=a.keys()==b.keys(),
                maximum_relative=max(errors,default=0.),failures=failures,rtol=1e-10,atol_D2=1e-18)
            records.append(rec);print(json.dumps(rec),flush=True)
    result=dict(Jmax=args.jmax,passed=all(x['identities_equal'] and not x['failures'] for x in records),
        checks=records,compilation=compilers,scope='Shared repaired eigensystem; tests contraction equivalence, not independent Hamiltonian physics or the c operator.')
    dump(out/'NATIVE-CACHE.json',portable_paths(result,{PACKAGE:'@package',inputs:'@inputs',out:'@output'}))
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
