"""Synthetic builder fixtures and portable-verification failure cases."""
import csv
import gzip
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'scripts'),str(ROOT/'src')]
from build_l33t import copy_orientation, orient_serialized_tie, missing_pair_donor, phys_parity, read_c_strengths
from highj_scheme_v2 import Key, State, Rate, SourceNetwork
from hybrid_e_lowj import build_factors
from reapply_source_frequency_policy import frequency_maps
from methanol_lamda.physics import einstein_a_from_smu2, interpolate_or_sqrt_extrapolate_rate
from methanol_lamda.isotopologues.ch3oh13.collisions import scaling_factor
from methanol_lamda.quantum import LevelKey
from review_support import expanded_hash, compare_dat, portable_paths
from test_release import small_file


class BuilderRules(unittest.TestCase):
    def test_absolute_einstein_prefactor(self):
        self.assertAlmostEqual(einstein_a_from_smu2(100000.,2.,5.)/4.6558e-6,1.,places=14)

    def test_identity_copy_uses_source_gap_when_order_flips(self):
        phys=dict(u=40.,l=10.,tu=1.,tl=2.); weights=dict(u=5.,l=3.)
        u,l,r,flipped=copy_orientation('u','l','tu','tl',(2e-10,),phys,weights,(100.,))
        self.assertEqual((u,l,flipped),('tl','tu',True))
        self.assertAlmostEqual(r[0]/(2e-10*5/3*math.exp(-1.4387768775039338*30/100)),1.)

    def test_identity_copy_without_reversal_is_unchanged(self):
        self.assertEqual(copy_orientation('u','l','tu','tl',(2.,),dict(u=9,l=0,tu=2,tl=1),{},(100.,)),('tu','tl',(2.,),False))

    def test_only_exact_physical_ties_may_reorient(self):
        u,l,r,flip=orient_serialized_tie('u','l',(3.,),dict(u=1.,l=1.),dict(u=1.,l=2.),dict(u=5.,l=3.))
        self.assertEqual((u,l,r,flip),('l','u',(5.,),True))
        with self.assertRaisesRegex(ValueError,'physically downward'):
            orient_serialized_tie('u','l',(3.,),dict(u=0.,l=1.),dict(u=0.,l=1.),dict(u=5.,l=3.))

    def network(self,u,l):
        return SourceNetwork({u:State(u,10.,2*u.j+1),l:State(l,0.,2*l.j+1)},
                             {frozenset((u,l)):Rate(u,l,(1e-10,))},(100.,),'p-H2')

    def test_signed_k_translation(self):
        for sign in (-1,1):
            u,l=Key('E','vt0',20,sign*15),Key('E','vt0',19,sign*14)
            a,b=Key('E','vt0',12,sign*11),Key('E','vt0',11,sign*10)
            choice,prefix=missing_pair_donor(u,l,self.network(a,b))
            self.assertEqual((choice.mapped_upper,choice.mapped_lower),(a,b))
            self.assertEqual(prefix,'k_translation_s4+')

    def test_a_translation_canonicalizes_zero_k(self):
        u,l=Key('A','vt0',13,3,'+'),Key('A','vt0',3,1,'-')
        a,b=Key('A','vt0',12,2,'+'),Key('A','vt0',2,0,'+')
        choice,_=missing_pair_donor(u,l,self.network(a,b))
        self.assertEqual(choice.mapped_lower,b)

    def test_nu10_overall_parity_flip(self):
        self.assertEqual(phys_parity('A:vt0:J3:K1:+'),-1)
        self.assertEqual(phys_parity('A:vt1:J3:K1:+'),1)
        self.assertEqual(phys_parity('A:vnu10:vt0:J3:K1:+'),1)
        self.assertEqual(phys_parity('A:vnu6:vt0:J3:K1:+'),-1)

    def test_c_operator_rejects_b_parity(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'c.csv'
            header='species,upper_vt,upper_j,upper_k,upper_parity,lower_vt,lower_j,lower_k,lower_parity,smu2_debye2\n'
            p.write_text(header+'A,0,1,1,-,0,0,0,+,1.0\n')
            self.assertEqual(len(read_c_strengths(p,LevelKey)),1)
            p.write_text(header+'A,0,1,1,+,0,0,0,+,1.0\n')
            with self.assertRaisesRegex(ValueError,'opposite-parity'):
                read_c_strengths(p,LevelKey)

    def test_isotope_scaling_direction_and_partner(self):
        self.assertAlmostEqual(scaling_factor(2),scaling_factor(3))
        self.assertLess(scaling_factor(6),scaling_factor(2))
        m12,m13=32.026214748490,33.029569583560
        for collider,mc in [(2,2.015650064460),(6,4.002603254130)]:
            expected=math.sqrt((m12/(m12+mc))/(m13/(m13+mc)))
            self.assertAlmostEqual(scaling_factor(collider),expected,places=14)

    def master(self, numerator=2., rf_e=3., missing_quartet=False):
        master={}
        for tag,k in [('Dp',numerator),('Do',numerator),('RFp',1.)]:
            u,l='A:vt0:J1:K0:+','A:vt0:J0:K0:+'
            master[tag,'A',frozenset((u,l))]=dict(u=u,l=l,eu=10.,el=0.,ts=(100.,),ks=(k,))
        j=2 if missing_quartet else 1
        u,l=f'E:vt0:J{j}:K+0',f'E:vt0:J{j-1}:K+0'
        for tag,k in [('Dp',10.),('Do',10.),('RFp',rf_e)]:
            master[tag,'E',frozenset((u,l))]=dict(u=u,l=l,eu=10.,el=0.,ts=(100.,),ks=(k,))
        return master

    def test_hybrid_quartet_magnitude_and_fallback(self):
        for missing in (False,True):
            fac,_=build_factors(self.master(missing_quartet=missing),(100.,),interpolate_or_sqrt_extrapolate_rate)
            for r in fac.values():
                self.assertAlmostEqual(r['f'][0],0.6)
                self.assertEqual(r['ratio_source'],['class-median' if missing else 'quartet'])

    def test_hybrid_honours_zero_ratio_and_rf_zero(self):
        for numerator,rf_e in [(0.,3.),(2.,0.)]:
            fac,_=build_factors(self.master(numerator,rf_e),(100.,),interpolate_or_sqrt_extrapolate_rate)
            for r in fac.values():self.assertEqual(r['f'],[0.])

    def test_observed_frequency_override_precedence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);home=root/'isotopologues/13ch3oh'
            (root/'config').mkdir();(home/'config').mkdir(parents=True)
            row=dict(upper_key='E:vt0:J1:K+0',lower_key='E:vt0:J0:K+0',frequency_mhz=10,source='fixture')
            xu=root/'xu.json';xu.write_text(json.dumps([row]))
            (root/'config/observed-frequencies.json').write_text(json.dumps({'overrides':[dict(row,frequency_mhz=20)]}))
            (home/'config/observed-frequencies.json').write_text(json.dumps({'catalog':'fixture','overrides':[]}))
            with patch('reapply_source_frequency_policy.catalog_file_to_overrides',return_value={}):
                maps=frequency_maps(root,xu)
                self.assertEqual(next(iter(maps['12'].values()))['frequency_mhz'],20)
                xu.write_text(json.dumps([row,row]))
                with self.assertRaisesRegex(ValueError,'duplicate'):frequency_maps(root,xu)


class PortableChecks(unittest.TestCase):
    def test_gzip_representation_is_not_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a.gz',Path(tmp)/'b.gz'
            data=b'row,data\n'*1000
            a.write_bytes(gzip.compress(data,compresslevel=1,mtime=0))
            b.write_bytes(gzip.compress(data,compresslevel=9,mtime=1))
            self.assertNotEqual(a.read_bytes(),b.read_bytes())
            self.assertEqual(expanded_hash(a),expanded_hash(b))
            b.write_bytes(gzip.compress(data+b'changed'))
            self.assertNotEqual(expanded_hash(a),expanded_hash(b))

    def compare(self,other):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a.dat',Path(tmp)/'b.dat'
            a.write_text(small_file());b.write_text(other)
            return compare_dat(a,b)

    def test_numerical_comparison_accepts_tolerance(self):
        self.assertTrue(self.compare(small_file().replace('1e-10','1.00000001e-10',1))['passed'])

    def test_numerical_comparison_rejects_rate_error(self):
        self.assertFalse(self.compare(small_file().replace('1e-10','1.001e-10',1))['passed'])

    def test_numerical_comparison_rejects_changed_identity(self):
        r=self.compare(small_file().replace('J1:K+0','J1:K+1'))
        self.assertFalse(r['passed']);self.assertGreater(r['missing'],0);self.assertGreater(r['extra'],0)

    def test_numerical_comparison_rejects_grid_change(self):
        self.assertFalse(self.compare(small_file().replace('10 20 30','11 20 30'))['passed'])

    def test_transition_row_order_is_not_identity(self):
        rows=small_file().splitlines()
        rows[2]='3';rows.insert(5,'3 2 5 E:vt0:J2:K+0')
        rows[6]='2';rows.insert(8,'2 3 2 2e-6 29.9792458 2.877553755')
        original='\n'.join(rows)+'\n'
        rows[7],rows[8]='1 '+rows[8].split(' ',1)[1],'2 '+rows[7].split(' ',1)[1]
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a',Path(tmp)/'b';a.write_text(original);b.write_text('\n'.join(rows)+'\n')
            self.assertTrue(compare_dat(a,b)['passed'])

    def test_zero_rate_has_no_absolute_tolerance(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a',Path(tmp)/'b'
            a.write_text(small_file().replace('1e-10','0',1))
            b.write_text(small_file().replace('1e-10','1e-300',1))
            self.assertFalse(compare_dat(a,b)['passed'])

    def test_patch_guard_survives_optimized_python(self):
        command='import sys; sys.path.insert(0,sys.argv[1]); from methanol_lamda.parity_repair import require; require(False)'
        r=subprocess.run([sys.executable,'-B','-O','-c',command,str(ROOT/'src')],capture_output=True,text=True)
        self.assertNotEqual(r.returncode,0);self.assertIn('RuntimeError',r.stderr)

    def test_recursive_provenance_paths(self):
        self.assertEqual(portable_paths({'command':['/sample/source/file','/sample/out/run']},
                                       {Path('/sample/source'):'@inputs',Path('/sample/out'):'@output'}),
                         {'command':['@inputs/file','@output/run']})
