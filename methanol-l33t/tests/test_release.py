"""Scientific invariants and malformed-file tests for the portable workflow."""
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'src')]
from methanol_lamda.physics import einstein_a_from_smu2, upward_from_downward, interpolate_or_sqrt_extrapolate_rate
from highj_scheme_v2 import Key, State, Rate, DonorChoice, SourceNetwork, select_donor, mapped_rates, target_upward_rates
from validate_products import inspect_lamda
from l33t import safe_relative


def small_file():
    grid='10 20 30 50 75 100 150 200 300 400'
    lines=['test E molecule','32.04186','2','1 0 1 E:vt0:J0:K+0',
           '2 1 3 E:vt0:J1:K+0','1','1 2 1 1e-6 29.9792458 1.4387768775','3']
    for c in (2,3,6):
        lines.extend([f'{c} test','1','10',grid,'1 2 1 '+' '.join(['1e-10']*10)])
    return '\n'.join(lines)+'\n'


class ScientificInvariants(unittest.TestCase):
    def test_frequency_override_preserves_strength_via_nu_cubed(self):
        self.assertAlmostEqual(einstein_a_from_smu2(2000,0.3,5)/einstein_a_from_smu2(1000,0.3,5),8)

    def test_lte_reverse_uses_target_energy(self):
        u,l=Key('E','vt0',15,1),Key('E','vt0',14,1)
        du,dl=Key('E','vt0',12,1),Key('E','vt0',11,1)
        # Target upper maps to the source LOWER level: reverse source direction.
        net=SourceNetwork({du:State(du,30,25),dl:State(dl,5,23)},
            {frozenset((du,dl)):Rate(du,dl,(2e-10,))},(100.,),'p-H2')
        choice=DonorChoice(dl,du,0.5,'test',True,du,dl)
        down=mapped_rates(choice,net)[0]
        source_expected=1e-10*25/23*math.exp(-1.4387768775039338*25/100)
        self.assertAlmostEqual(down/source_expected,1)
        upper,lower=State(u,200,31),State(l,80,29)
        up=target_upward_rates((down,),(100.,),upper,lower)[0]
        self.assertAlmostEqual(up/down,31/29*math.exp(-1.4387768775039338*120/100))
        self.assertAlmostEqual(up,upward_from_downward(down,31,29,120,100))

    def test_nonzero_jump_never_uses_zero_jump_donor(self):
        u,l=Key('E','vt0',19,1),Key('E','vt0',18,-1)
        a,b=Key('E','vt0',12,1),Key('E','vt0',12,-1)
        net=SourceNetwork({a:State(a,5,25),b:State(b,0,25)},
            {frozenset((a,b)):Rate(a,b,(1e-10,))},(100.,),'p-H2')
        self.assertIsNone(select_donor(u,l,net,12))

    def test_exact_jump_and_signed_k_preserved(self):
        u,l=Key('E','vt0',20,-1),Key('E','vt0',18,2)
        a,b=Key('E','vt0',12,-1),Key('E','vt0',10,2)
        net=SourceNetwork({a:State(a,5,25),b:State(b,0,21)},
            {frozenset((a,b)):Rate(a,b,(1e-10,))},(100.,),'p-H2')
        choice=select_donor(u,l,net,12)
        self.assertEqual((choice.mapped_upper,choice.mapped_lower,choice.scale),(a,b,1))
        with self.assertRaises(ValueError):
            select_donor(u,Key('E','vt1',18,2),net,12)

    def test_zero_rate_and_temperature_endpoint(self):
        f=interpolate_or_sqrt_extrapolate_rate
        self.assertEqual(f((10.,20.),(0.,2.),15),1)
        self.assertEqual(f((10.,20.),(1.,2.),80),4)


class SerializationAndCustody(unittest.TestCase):
    def inspect(self,text):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'test.dat';p.write_text(text)
            return inspect_lamda(p)

    def test_complete_structure(self):
        result=self.inspect(small_file())
        self.assertEqual((result['levels'],result['radiative'],result['collision_rows']),(2,1,3))

    def test_reject_truncation(self):
        with self.assertRaisesRegex(ValueError,'truncated'):
            self.inspect(small_file().rsplit('\n',2)[0]+'\n')

    def test_reject_reversed_endpoints(self):
        with self.assertRaisesRegex(ValueError,'endpoints'):
            self.inspect(small_file().replace('1 2 1 1e-6','1 1 2 1e-6'))

    def test_reject_nan_negative_and_mixed_species(self):
        for text in [small_file().replace('1e-10','nan',1),
                     small_file().replace('1e-10','-1e-10',1),
                     small_file().replace('E:vt0:J1:K+0','A:vt0:J1:K0:+')]:
            with self.subTest(text=text[-30:]),self.assertRaises(ValueError):
                self.inspect(text)

    def test_manifest_escape_rejected(self):
        for p in ('../outside','/absolute'):
            with self.assertRaises(ValueError):
                safe_relative(p)


if __name__=='__main__':
    unittest.main()
