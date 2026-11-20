#!/usr/bin/env python3
"""Stage 1 of the l33t build: parity-repaired BELGI spectroscopy to J<=46.

Runs the codex parity/full-vector repair (block solve, global bijective A
assignment, full nine-block export, cached exact dipole contraction; pinned in
D45/D46) with a patched copy of the strengths driver:
  * level capacity 2,500 -> 20,000; output box J<=46 (BELGI's array limit), |K|<=46;
  * a VCO_B mode: the cross-vibrational operator along b instead of a (v1 A'' proxy);
  * a VCO_C mode (v2): the operator along c. In the a-quantized basis b and c differ only
    by the relative sign of their q = +-1 components, b ~ (T-1 - T+1)/sqrt2 and
    c ~ i(T-1 + T+1)/sqrt2; the common phase i drops out of |amplitude|^2. The
    cached contraction therefore omits its delta = +-1 factor. Cross-vibrational pairs
    are restricted to those with at least one vt0 endpoint (the ground function that
    carries the band), so VCO_C with vt<=4 also yields the Delta vt != 0 branches.
The codex worktree is read-only; the patches are applied to a staged copy.

Modes (Xu 2008 12CH3OH constants for --isotope A12; Xu et al. 2014 Table 2 13CH3OH
constants for --isotope A13):
  PERMANENT  vt<=4, J<=JMAX, |K|<=KMAX, permanent Xu dipoles, floor 1e-6 D^2
  VCO_A      vt=0, pure-a 0.27-D operator between ground functions (CO, rock, A' bands)
  VCO_B      vt=0, pure-b 0.27-D operator between ground functions (v1 A'' proxy)
  VCO_C      vt0 <-> vt<=4, pure-c 0.27-D operator between ground functions (v2 A'' bands)

Usage: run_spectroscopy_l33t.py --inputs DIR --output DIR [--isotope A12|A13] [--jmax 46 --kmax 22]
       [--modes PERMANENT VCO_A VCO_B VCO_C]
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

DRIVER_PATCHES = [
    ("      parameter (maxlev=2500)", "      parameter (maxlev=20000)"),
    ("      if (output_jmax .lt. 1 .or. output_jmax .gt. 22) stop 2",
     "      if (output_jmax .lt. 1 .or. output_jmax .gt. 46) stop 2"),
    ("      if (output_kmax .lt. 0 .or. output_kmax .gt. 14 .or.",
     "      if (output_kmax .lt. 0 .or. output_kmax .gt. 46 .or."),
    ("      cross_vibrational_mode = dipole_mode .eq. 'VCO_A'",
     "      cross_vibrational_mode = dipole_mode .eq. 'VCO_A' .or.\n     &    dipole_mode .eq. 'VCO_B' .or.\n"
     "     &    dipole_mode .eq. 'VCO_C'"),
    ("         cst(1,102) = 0.27d0\n         cst(1,103) = 0.d0\n         cst(1,104) = 0.d0\n",
     "         cst(1,102) = 0.27d0\n         cst(1,103) = 0.d0\n         cst(1,104) = 0.d0\n"
     "         if (dipole_mode .eq. 'VCO_B' .or. dipole_mode .eq. 'VCO_C') then\n"
     "C           l33t: the same 0.27-D reference moment along b (perpendicular bands);\n"
     "C           VCO_C turns the cached b contraction into the c contraction.\n"
     "            cst(1,102) = 0.d0\n            cst(1,104) = 0.27d0\n         endif\n"
     "         if (dipole_mode .eq. 'VCO_C') call repair_set_c_axis()\n"),
    ("            if (ia .eq. ib .and. .not. cross_vibrational_mode)\n     &         goto 20\n",
     "            if (ia .eq. ib .and. .not. cross_vibrational_mode)\n     &         goto 20\n"
     "C           l33t: a band is carried by a vt0 ground function (no-op for vt=0 runs).\n"
     "            if (cross_vibrational_mode .and. svt(ia) .ne. 0 .and.\n"
     "     &          svt(ib) .ne. 0) goto 20\n"),
]
CACHE_PATCHES = [
    ("  real(8), save :: cached_mua=huge(1d0),cached_mub=huge(1d0)\n",
     "  real(8), save :: cached_mua=huge(1d0),cached_mub=huge(1d0)\n"
     "  logical, save :: c_axis=.false.  ! l33t VCO_C: perpendicular c operator\n"),
    ("        if(delta/=0) angular=angular*delta/sqrt(2d0)\n",
     "        if(delta/=0) then\n"
     "          ! b = (T(-1)-T(+1))/sqrt2 carries the sign delta; c = i(T(-1)+T(+1))/sqrt2 does not\n"
     "          ! (its common phase i drops out of the squared amplitude).\n"
     "          if(c_axis) then\n"
     "            angular=angular/sqrt(2d0)\n"
     "          else\n"
     "            angular=angular*delta/sqrt(2d0)\n"
     "          end if\n"
     "        end if\n"),
    ("end module\n",
     "end module\n\nsubroutine repair_set_c_axis()\n  use methanol_dipole_cache\n  implicit none\n  c_axis=.true.\nend subroutine\n"),
]
MODE_VT = {"PERMANENT": 4, "VCO_A": 0, "VCO_B": 0, "VCO_C": 4}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=(__doc__ or 'Generate native L33T torsion-rotation spectroscopy.').split("\n")[0])
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--jmax", type=int, default=46)
    ap.add_argument("--kmax", type=int, default=22)
    ap.add_argument("--floor", type=float, default=1e-6)
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--isotope", choices=("A12", "A13"), default="A12")
    ap.add_argument("--modes", nargs="+", choices=tuple(MODE_VT), default=["PERMANENT", "VCO_A", "VCO_C"])
    a = ap.parse_args()
    inputs = a.inputs.resolve()
    a.output = a.output.resolve()  # the compiler runs with cwd=build; paths must be absolute
    repo = Path(__file__).resolve().parents[1]
    if a.output.exists():
        raise FileExistsError("choose a new output directory")
    stage = a.output / "staged-project"
    (stage / "fortran").mkdir(parents=True)
    for f in (repo / "fortran").iterdir():
        shutil.copy2(f, stage / "fortran" / f.name)
    os.symlink(repo / "src", stage / "src")
    os.symlink(inputs / "vendor", stage / "vendor")
    driver = stage / "fortran" / "belgi_frozen_strengths.f"
    original_sha = sha(driver)
    text = driver.read_text()
    for old, new in DRIVER_PATCHES:
        if text.count(old) != 1:
            raise RuntimeError(f"driver patch anchor not unique: {old[:60]}")
        text = text.replace(old, new)
    driver.write_text(text)
    cache = stage / "fortran" / "cached_dipole_repair.f90"
    cache_original_sha = sha(cache)
    text = cache.read_text()
    for old, new in CACHE_PATCHES:
        if text.count(old) != 1:
            raise RuntimeError(f"cache patch anchor not unique: {old[:60]}")
        text = text.replace(old, new)
    cache.write_text(text)
    sys.path.insert(0, str(repo / "src"))
    from methanol_lamda.hamiltonian import (_compile_driver, _find_compiler, parse_final_constants,
                                            render_complete_cst_vector, _parse_level_stdout)
    from methanol_lamda.strengths import _parse_strength_stdout, _write_strengths_csv, _validate_enumeration
    build = a.output / "build"
    build.mkdir()
    exe, cmd, version, patches = _compile_driver(stage, stage / "vendor/belgi", build, _find_compiler(None),
                                                 driver_filename="belgi_frozen_strengths.f",
                                                 executable_name="l33t_strengths", parity_repair=True, cached_dipoles=True)
    constants = inputs / ("data/raw/local-cache/xu2008-final-constants.txt" if a.isotope == "A12" else "isotopologues/13ch3oh/constants/table2-2014-constants.txt")
    raw = render_complete_cst_vector(parse_final_constants(constants, require_xu_final=a.isotope == "A12"))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    runs = []
    for mode in a.modes:
        vt = MODE_VT[mode]
        d = a.output / a.isotope / mode
        (d / "assignment").mkdir(parents=True)
        args = [str(exe), str(a.jmax), str(a.kmax), str(vt), repr(a.floor), mode]
        t0 = time.monotonic()
        with (d / "native.stdout").open("w") as so, (d / "native.stderr").open("w") as se:
            run = subprocess.run(args, input=raw, text=True, stdout=so, stderr=se, cwd=d,
                                 env=dict(env, BELGI_REPAIR_AUDIT_DIR=str(d / "assignment")), timeout=a.timeout)
        seconds = time.monotonic() - t0
        if run.returncode:
            raise RuntimeError(f"{mode} failed ({run.returncode}): {(d / 'native.stderr').read_text()[-3000:]}")
        out = (d / "native.stdout").read_text()
        levels = _parse_level_stdout(out)
        with (d / "levels.csv").open("w", newline="") as h:
            w = csv.writer(h)
            w.writerow(["species", "level_key", "vt", "j", "k", "parity", "energy_cm1", "energy_raw_cm1", "weight"])
            for lv in levels:
                w.writerow([lv.key.species, lv.key.encode(), lv.key.vt, lv.key.j, lv.key.k, lv.key.parity,
                            format(lv.energy_cm1, ".17g"), format(lv.energy_raw_cm1, ".17g"), 2 * lv.key.j + 1])
        rawcsv = out.split("BEGIN_FROZEN_STRENGTHS_CSV\n", 1)[1].split("END_FROZEN_STRENGTHS_CSV", 1)[0]
        (d / "native-lines.csv").write_text(rawcsv)
        summary = next(x for x in out.splitlines() if x.startswith("FROZEN_STRENGTH_SUMMARY,"))
        rec = dict(mode=mode, vt_max=vt, seconds=seconds, levels=len(levels), command=args, summary_line=summary)
        if mode == "PERMANENT":
            lev, rows, summ = _parse_strength_stdout(out)
            _validate_enumeration(lev, rows, summ, j_max=a.jmax, abs_k_max=a.kmax, vt_max=4, smu2_min_debye2=a.floor)
            _write_strengths_csv(d / "lines.csv", rows, source_id=f"l33t-parity-full9-{a.isotope}",
                                 dipole_model_id="Xu2008_constant_mu_a_mu_b_cos3gamma_zero")
            rec["strength_rows"] = len(rows)
            del rows, lev
        (d / "native.stdout").unlink()  # large; levels/lines/native-lines retain its content
        rec["files"] = {p.name: dict(bytes=p.stat().st_size, sha256=sha(p)) for p in d.iterdir() if p.is_file()}
        runs.append(rec)
        print(json.dumps({k: rec[k] for k in ("mode", "seconds", "levels", "summary_line")}), flush=True)
    manifest = dict(status="l33t-stage1-spectroscopy", isotope=a.isotope, inputs=str(inputs), driver_original_sha256=original_sha,
                    driver_patched_sha256=sha(driver), driver_patches=[o for o, _ in DRIVER_PATCHES],
                    cache_original_sha256=cache_original_sha, cache_patched_sha256=sha(cache),
                    cache_patches=[o for o, _ in CACHE_PATCHES], modes=a.modes, compiler=version,
                    compile_command=cmd, source_patches=patches, constants=str(constants), constants_sha256=sha(constants),
                    jmax=a.jmax, kmax=a.kmax, floor_debye2=a.floor, runs=runs, script_sha256=sha(Path(__file__)))
    from review_support import portable_paths
    manifest["native_cache_comparison"] = "not performed by generation; run check_native_cache.py separately"
    manifest = portable_paths(manifest, {inputs: "@inputs", a.output: "@output", repo: "@package"})
    (a.output / "SPECTROSCOPY-MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
