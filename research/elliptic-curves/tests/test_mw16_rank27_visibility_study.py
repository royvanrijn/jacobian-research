"""Adversarial checks for the new pointwise visibility certificate."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'elliptic-curves/cas'))
import study_mw16_rank27_visibility as study


class VisibilityWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(study.INPUT.read_text())
        cls.record = cls.data['arms']['adaptive'][85]['search']
        cls.point = tuple(map(study.cert.F, cls.data['point_proof']['discovery_points'][-1]))

    def test_winning_signed_point_has_raw_square_witness(self):
        result = study.locate(self.record, self.point)
        self.assertEqual(result['coordinate'], ['-73571', '79466'])
        self.assertTrue(result['returned'])
        self.assertTrue(result['square_hit_recorded'])
        negative = study.locate(self.record, (self.point[0], -self.point[1]))
        self.assertGreater(negative['minimum_affine_height'], 100000)

    def test_quartic_coefficient_corruption_fails(self):
        record = copy.deepcopy(self.record)
        record['coefficients'][0] = str(int(record['coefficients'][0]) + 1)
        with self.assertRaises(ArithmeticError):
            study.locate(record, self.point)

    def test_off_curve_oracle_fails(self):
        with self.assertRaises(ValueError):
            study.locate(self.record, (self.point[0], self.point[1] + 1))

    def test_square_hit_and_point_output_are_separate_witnesses(self):
        record = copy.deepcopy(self.record)
        record['primitive_square_hits'] = []
        result = study.locate(record, self.point)
        self.assertTrue(result['returned'])
        self.assertFalse(result['square_hit_recorded'])

    def test_wrong_specialization_fails_before_geometry(self):
        data = copy.deepcopy(self.data)
        data['point_proof']['parameter'] = '1867/270'
        with self.assertRaises(ArithmeticError):
            study.expected(data)


if __name__ == '__main__':
    unittest.main()
