"""Explicit, checksum-described prediction overlays for the parity repair.

The historical BELGI source tree remains frozen. Opt-in consumers compile these
small reviewed transformations plus the new physical-parity module. No patch
modifies a pre-existing molecular input or removes radiative rows by their labels.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import re

REPAIR_ID = "native-PER23-blocksolve-global-all9-overlap-v1"


def require(condition, message="BELGI patch anchor changed"):
    if not condition:
        raise RuntimeError(message)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def materialize_belgi_repair_sources(
    project_root: Path, implementations: list[Path], build_dir: Path
) -> tuple[list[Path], list[Path], dict[str, object]]:
    """Write compilation-only overlays; return sources, extra module and pins."""
    descriptions = {
        "rottor.f": "A: explicit physical PER23 eigenspaces, H block solves and full-state bijective reference assignment; E: full unchanged symmetric eigenproblem",
        "encal.f": "A: consume globally assigned parity-consistent column; both species retain full nine torsional blocks",
        "inten.f": "remove state-dependent coefficient-product pruning; use double-precision sqrt(2) in complete signed dipole contractions",
    }
    result: list[Path] = []
    patches: dict[str, object] = {}
    for source in implementations:
        if source.name not in descriptions:
            result.append(source)
            continue
        text = source.read_text(encoding="ascii")
        if source.name == "rottor.f":
            old = "      call dsyev('V', 'U', nrotor, EGVC, ndmx, egvl, \n     *     work, lwork, ierr)"
            require(text.count(old) == 1, "audited full rotation-torsion solve changed")
            text = text.replace(old, "      call repair_diagonalize(N,ISIG,NDIMTO,H,EGVL,EGVC,A)\n      ierr=0")
        elif source.name == "encal.f":
            old = "      integer encal2"
            require(text.count(old) == 1)
            text = text.replace(old, "      integer encal2, repair_state_index")
            old = "      if ( (label .eq. 1 .or. label .eq. 2) .and. (N .gt. 5) .and."
            require(text.count(old) == 1)
            text = text.replace(old,
                "      if (isig .eq. 0) then\n"
                "         nbr(lo)=repair_state_index(N,IV,KK,PR)\n"
                "      else if ( (label .eq. 1 .or. label .eq. 2) .and. (N .gt. 5) .and.")
            old = "         IF(IVTORP.GE.9) EGVC3(II,IVTORP)=0.D0"
            require(text.count(old) == 1)
            text = text.replace(old, "C        Retain the complete nine-block eigenvector without projection.")
        elif source.name == "inten.f":
            text, count = re.subn(r"^\s*IF\(ABS\(CPROD(?:0|1|M1)\)\.LT\.EPS\) GO TO \d+[^\n]*$", "C         Complete contraction: no state-dependent coefficient pruning.", text, flags=re.MULTILINE)
            require(count == 9, f"expected nine coefficient-product pruning sites, got {count}")
            require(text.count("SQRT(2.)") == 6)
            text = text.replace("SQRT(2.)", "SQRT(2.D0)")
        dest = build_dir / ("repaired-" + source.name)
        dest.write_text(text, encoding="ascii")
        result.append(dest)
        patches[source.name] = {
            "repair_id": REPAIR_ID,
            "reason": descriptions[source.name],
            "original_sha256": _sha(source),
            "compiled_sha256": _sha(dest),
        }
    module = project_root / "fortran/parity_state_repair.f90"
    patches["parity_state_repair.f90"] = {
        "repair_id": REPAIR_ID, "source_sha256": _sha(module),
        "reason": "native physical parity block diagonalization and gauge-invariant full-state Hungarian assignment",
    }
    return result, [module], patches


def materialize_cached_strength_driver(
    project_root: Path, driver: Path, build_dir: Path
) -> tuple[Path, list[Path], dict[str, object]]:
    """Opt-in exact matrix contraction; independent direct INTEN remains built."""
    require(driver.name == 'belgi_frozen_strengths.f')
    text = driver.read_text(encoding='ascii')
    old = "      double precision se(maxlev),smu2_min,rint,cf,frequency"
    require(text.count(old)==1)
    text=text.replace(old, old+'\n      double precision repair_a,repair_b,repair_component_cf\n'
        '      common/REPAIR_INTEN_COMPONENTS/repair_a,repair_b,\n'
        '     &     repair_component_cf')
    old='            call inten(rint,cst,1,ivvpair,inkpair,muabc,ictt,cf)'
    require(text.count(old)==1)
    text=text.replace(old,old.replace('call inten(', 'call inten_cached('))
    old="     & 'contribution_b_percent,cancellation_factor'"
    require(text.count(old)==1)
    text=text.replace(old,"     & 'contribution_b_percent,cancellation_factor,'//\n"
        "     & 'amplitude_a_debye,amplitude_b_debye,'//\n"
        "     & 'component_cancellation_factor'")
    old='     &           se(ilow),frequency,ssq,components,ictt(1),ictt(2),cf'
    require(text.count(old)==1)
    text=text.replace(old,old+',\n     &           repair_a,repair_b,repair_component_cf')
    old="     &       I0,',',I0,',',ES25.17E3)"
    require(text.count(old)==1)
    text=text.replace(old,"     &       I0,',',I0,',',ES25.17E3,3(',',ES25.17E3))")
    old='     &             candidates,0.d0'
    require(text.count(old)==1)
    text=text.replace(old,'     &             0,0.d0')
    text=text.replace('C     The fourth field counts direct native-INTEN evaluations.  No surrogate\nC     contraction is used, hence its difference from native INTEN is zero.',
        'C     Exact cached full-operator contractions; no per-row native INTEN\nC     calls claimed. Separate bounded native comparison validates the cache.')
    dest=build_dir/'repaired-cached-strengths.f';dest.write_text(text,encoding='ascii')
    extras=[project_root/'fortran/cached_dipole_repair.f90', project_root/'fortran/cached_inten.f']
    patches={'cached_strength_driver':{'repair_id':REPAIR_ID,'original_sha256':_sha(driver),'compiled_sha256':_sha(dest),
        'reason':'Exact full-operator BLAS contractions with signed components; no forbidden-edge filtering; direct INTEN validation is separate'},
        **{x.name:{'source_sha256':_sha(x)} for x in extras}}
    return dest,extras,patches
