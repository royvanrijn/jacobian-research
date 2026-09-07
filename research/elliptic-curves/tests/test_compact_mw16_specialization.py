"""Exact rational coordinate and failure-boundary regressions for MW16 inputs."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import compact_mw16_specialization as spec


class CompactMW16SpecializationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.families=json.loads(spec.ATLAS.read_text())['families']
        cls.native={r['presentation_id']:r for r in json.loads((ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json').read_text())['presentations']}

    def test_rational_parameter_transports_each_equation(self):
        t=F(2,3)
        for family in self.families:
            with self.subTest(family=family['fibration_id']):
                model,points=spec.specialize(family,t)
                self.assertEqual(len(points),16)
                old=spec.native_parameter(family,t)
                self.assertIsNotNone(old)
                a,b,c,d=map(F,family['base_matrix_a_b_c_d'])
                scale=F(family['scale_u'])/(c*t+d)**2/t.denominator**2
                native=self.native[family['presentation_id']]['pencil']
                self.assertEqual(model[3]*scale**4,spec.polynomial(native['A_coefficients_low_to_high'],old))
                self.assertEqual(model[4]*scale**6,spec.polynomial(native['B_coefficients_low_to_high'],old))

    def test_corrupt_point_is_rejected(self):
        family=copy.deepcopy(self.families[0])
        values=family['sections'][0]['Y']['numerator_coefficients_low_to_high']
        values[0]=str(F(values[0])+1)
        with self.assertRaisesRegex(ArithmeticError,'section equation'):
            spec.specialize(family,1)

    def test_missing_section_is_rejected(self):
        family=copy.deepcopy(self.families[0]);family['sections'].pop()
        with self.assertRaisesRegex(ArithmeticError,'section roster'):
            spec.specialize(family,1)

    def test_singular_fibre_is_rejected(self):
        family=copy.deepcopy(self.families[0])
        family['A_coefficients_low_to_high']=['0'];family['B_coefficients_low_to_high']=['0']
        with self.assertRaisesRegex(ArithmeticError,'singular fibre'):
            spec.specialize(family,1)

    def test_section_pole_is_rejected(self):
        family=copy.deepcopy(self.families[0])
        family['sections'][0]['X']['denominator_coefficients_low_to_high']=['-1','1']
        with self.assertRaisesRegex(ArithmeticError,'pole'):
            spec.specialize(family,1)


if __name__=='__main__': unittest.main()
