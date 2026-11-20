"""Audited v3 source helpers; only functions needed for this product are retained."""
from contextlib import contextmanager
import csv
import math
import statistics
from pathlib import Path
from methanol_lamda.quantum import LevelKey
from methanol_lamda.co_stretch import (COStretchKey, CM1_TO_MHZ, MAIN_HITRAN_RELATIVE_PATH,
    HITRAN_E_GROUND_ZERO_CM1, _parse_hitran_local, _resolve_hitran_ground_identity,
    _resolve_hitran_upper_identity)
from methanol_lamda.ch3_rock import C13_ROCK_MODEL
from methanol_lamda.physics import einstein_a_from_smu2
from methanol_lamda import vibrational_revision as vr

def physical_parity(k):
    if k.species!='A':return None
    # CO and parallel rock fundamental identities inherit the vt0 symmetry.
    vt=k.vt if isinstance(k,LevelKey) else 0
    return (1 if k.parity=='+' else -1)*(-1)**(k.j+vt)


def read_strengths(path):
    result=[];seen=set();intern={}
    def canonical(k):
        if k not in intern:intern[k]=k
        return intern[k]
    with Path(path).open() as f:
        for row_number,row in enumerate(csv.DictReader(f),1):
            if row.get('upper_key'):
                a,b=LevelKey.decode(row['upper_key']),LevelKey.decode(row['lower_key'])
            else:
                sp=row['species'].strip().upper()
                a=LevelKey(sp,int(row.get('upper_vt',0)),int(row['upper_j']),int(row['upper_k']),row.get('upper_parity','').strip())
                b=LevelKey(sp,int(row.get('lower_vt',0)),int(row['lower_j']),int(row['lower_k']),row.get('lower_parity','').strip())
            a,b=canonical(a),canonical(b)
            pair=frozenset((a,b))
            if pair in seen:raise ValueError('duplicate native strength endpoint pair')
            s=float(row['smu2_debye2'])
            if s<=0 or not math.isfinite(s):raise ValueError('invalid strength in native admitted CSV')
            if a.species=='A' and physical_parity(a)==physical_parity(b):raise ValueError('native A parity violation before product cutoff')
            if abs(a.j-b.j)>1 or a.j==b.j==0:raise ValueError('native rank-one violation before product cutoff')
            # The immutable native CSV retains amplitudes and every diagnostic;
            # keep only its row identity in memory for the product crosswalk.
            seen.add(pair);result.append((a,b,s,{'transition_index':row.get('transition_index',row.get('row_index',row_number))}))
    return result


@contextmanager
def new_main_ground(path):
    old=vr.MAIN_LEVELS;vr.MAIN_LEVELS=Path(path).resolve()
    try:yield
    finally:vr.MAIN_LEVELS=old


def calibrate_repaired_co(root,profile,new_ground,vco,co_by_species):
    """Reapply the SAME thirteen-feature normalization to the repaired operator."""
    raw=(root/MAIN_HITRAN_RELATIVE_PATH).read_text().splitlines()
    strengths={frozenset((a,b)):s for a,b,s,_ in vco};rows=[]
    settings=profile['co_calibration']
    for source in settings['rows']:
        if source['wavenumber_cm1']==settings['excluded_feature_cm1']:continue
        line=raw[source['raw_row']-1]
        if abs(float(line[3:15])-source['wavenumber_cm1'])>=3e-4:raise ValueError('CO calibration source row/frequency mismatch')
        u0,l0=_parse_hitran_local(line[97:112]),_parse_hitran_local(line[112:127])
        l,_,_=_resolve_hitran_ground_identity(l0,float(line[45:55]),new_ground)
        if l is None:raise ValueError('CO calibration lower identity unmapped with repaired ground levels')
        u=_resolve_hitran_upper_identity(u0,l0,l);nu=source['wavenumber_cm1']
        s=strengths[frozenset((u,l))]
        if s<=0:raise ValueError('calibration feature absent from repaired pure-a operator')
        cu=COStretchKey.from_ground(u)
        if cu not in co_by_species[u.species]:raise ValueError('CO calibration outside selected basis')
        elo=new_ground[l]+(HITRAN_E_GROUND_ZERO_CM1 if l.species=='E' else 0)
        # Match the prior policy: A uses emitted model frequency, the absorption
        # comparison uses the published feature frequency.
        modelnu=(co_by_species[u.species][cu].energy_cm1-new_ground[l])*CM1_TO_MHZ
        a=einstein_a_from_smu2(modelnu,s,2*u.j+1)
        predicted=vr.absorption(a,nu,elo,2*u.j+1,settings['temperature_K'],settings['Q_over_common_spin'])*2.68675e19*273.15/settings['temperature_K']
        rows.append(dict(source,upper_key=cu.encode(),lower_key=l.encode(),unscaled_smu2_debye2=s,original_predicted=predicted,original_over_xu=predicted/source['xu_cm2_inverse_atm_inverse']))
    if len(rows)!=13:raise ValueError('CO normalization must retain exactly 13 source features')
    scale=1/statistics.median(r['original_over_xu'] for r in rows)
    for r in rows:r['revised_over_xu']=r['original_over_xu']*scale
    return dict(rate_scale=scale,effective_moment_debye=.27*math.sqrt(scale),rows=rows,status='Same source-bound in-sample normalization; no new absolute-accuracy claim')

