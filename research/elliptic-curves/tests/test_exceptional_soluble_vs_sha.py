"""Arithmetic corruption and incomplete-pairing regressions; Sage Python."""
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import runpy
import unittest

from sage.all import EllipticCurve, GF, Matrix, QQ

ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(str(ROOT/'elliptic-curves/cas/compare_exceptional_soluble_vs_sha.sage'))


class ComparisonTests(unittest.TestCase):
    def test_restricted_pairing_certifies_obstruction_only(self):
        b = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]
        r = M['pairing_profile'](b, [[0, 0, 1]])
        self.assertEqual(r['soluble_intersection_dimension_interval'], [1, 1])
        self.assertEqual(r['sha_image_dimension_lower_bound'], 2)
        self.assertEqual(r['provably_insoluble_class_count'], 6)
        self.assertIsNone(r['whole_curve_rank_upper_bound'])
        self.assertIsNone(r['full_radical_dimension'])
        self.assertEqual(M['pairing_profile']([[0]], [])['soluble_intersection_dimension_interval'], [0, 1])

    def test_unknown_and_false_witness_rejected(self):
        for b in ([[0, None], [None, 0]], [[1]], [[0, 1], [0, 0]], [[0, 2], [2, 0]]):
            with self.assertRaises(ArithmeticError):
                M['pairing_profile'](b, [])
        with self.assertRaises(ArithmeticError):
            M['pairing_profile']([[0, 1], [1, 0]], [[1, 0]])

    def test_profile_invariant_under_change_of_basis(self):
        b = Matrix(GF(2), [[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        c = Matrix(GF(2), [[1, 0, 1], [1, 1, 0], [0, 0, 1]])
        q = Matrix(GF(2), [[0, 0, 1]])
        rows = lambda x: [list(map(int, r)) for r in x.rows()]
        self.assertEqual(M['pairing_profile'](rows(b), rows(q)),
                         M['pairing_profile'](rows(c*b*c.transpose()), rows(q*c.inverse())))

    def test_marked_cover_reduction_and_corruption(self):
        e = EllipticCurve(QQ, [0, 0, 0, -1, 1])
        r = M['reduce_pointed'](e, e(0, 1), 'P')
        c = M['verify_pointed'](r, e)
        self.assertGreater(c['total_presentation_slots_bits'], c['coefficient_slots_bits'])
        for field in ('quartic', 'witness_X_Z_Y', 'short_point'):
            bad = deepcopy(r)
            bad[field][0] = str(QQ(bad[field][0])+1)
            with self.assertRaises(ArithmeticError):
                M['verify_pointed'](bad, e)
        bad = deepcopy(r)
        bad['raw_parameter_from_reduced'] = [['0', '0'], ['0', '0']]
        with self.assertRaises(ArithmeticError):
            M['verify_pointed'](bad, e)
        bad = deepcopy(r)
        bad['cover_map'] = 'phi'
        with self.assertRaises(ArithmeticError):
            M['verify_pointed'](bad, e)

    def test_four_torsion_sha_model_is_invisible_on_two_torsion(self):
        # Abstract symplectic module, not an assertion of an elliptic curve.
        pairing = lambda x, y: Fraction(x[0]*y[1]-x[1]*y[0], 4) % 1
        two_torsion = [(0, 0), (2, 0), (0, 2), (2, 2)]
        self.assertTrue(all(pairing(x, y) == 0 for x in two_torsion for y in two_torsion))
        self.assertEqual(pairing((1, 0), (0, 1)), Fraction(1, 4))


if __name__ == '__main__':
    unittest.main()
