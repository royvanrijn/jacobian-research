"""Proof-tamper regressions for the retained blind recovery, no new searches."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
import verify_constructed_class_blind as exact
import replay_constructed_class_quartic as quartic


class BlindRecoveryTamperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = exact.EXPERIMENT
        cls.covers = json.loads((cls.base/'covers.json').read_text())
        cls.witness = json.loads((cls.base/'cover-worker-witness.json').read_text())
        cls.f = list(map(F, cls.covers['cubic_ascending']))
        cls.model = list(map(F, cls.covers['original_curve']))

    def verify(self, proposed):
        return exact.verify_case(self.covers['cases'][0], proposed, self.f, self.model)

    def test_retained_points_replay(self):
        for source, proposed in zip(self.covers['cases'], self.witness['cases']):
            result = exact.verify_case(source, proposed, self.f, self.model)
            self.assertTrue(result['exact_square_identity'])

    def test_bad_projective_transport_rejected(self):
        p = copy.deepcopy(self.witness['cases'][0])
        p['coordinate_matrix'][0][0] = str(F(p['coordinate_matrix'][0][0])+1)
        with self.assertRaises(AssertionError):
            self.verify(p)

    def test_wrong_primitive_point_rejected(self):
        p = copy.deepcopy(self.witness['cases'][0])
        p['primitive_cover_point'][0] = str(int(p['primitive_cover_point'][0])+1)
        with self.assertRaises(AssertionError):
            self.verify(p)

    def test_swapped_classes_rejected(self):
        with self.assertRaises(AssertionError):
            self.verify(self.witness['cases'][1])

    def test_wrong_quartic_root_rejected(self):
        p = self.witness['cases'][0]
        q = json.loads((self.base/'cover-worker-quartic-c6.json').read_text())
        quartic.replay_case(p, q)
        q['reduced_quartic_point'][1] = str(F(q['reduced_quartic_point'][1])+1)
        with self.assertRaises(AssertionError):
            quartic.replay_case(p, q)


if __name__ == '__main__':
    unittest.main()
