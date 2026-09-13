"""Explicit small Sage group-law controls; excluded from navigation checks."""

from copy import deepcopy
from pathlib import Path
import sys
import runpy
import unittest
from unittest.mock import patch

from sage.all import EllipticCurve, QQ

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'elkies-k3/scripts'))
import replay_r17_norm12_native_icarm_quotient_audit as replay


class RetainedRelationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.helper = runpy.run_path(str(replay.PRODUCER))

    def setUp(self):
        self.curve = EllipticCurve(QQ, [-7, 10])
        self.points = [self.curve(1, 2), self.curve(2, 2)]
        first, second = self.points
        self.targets = [2 * first - 3 * second, -first + 4 * second, self.curve(0)]
        self.rows = [[2, -1, 0], [-3, 4, 0]]

    def test_exact_relations_need_no_numerical_heights(self):
        with patch.object(type(self.curve), 'height_pairing_matrix', side_effect=AssertionError('height discovery forbidden')):
            replay.verify_integer_relations(self.curve, self.points, self.targets, self.rows)

    def test_wrong_integer_relation_is_rejected(self):
        rows = deepcopy(self.rows)
        rows[0][0] += 1
        with self.assertRaisesRegex(ArithmeticError, 'relation 1 failed'):
            replay.verify_integer_relations(self.curve, self.points, self.targets, rows)

    def test_point_order_is_part_of_the_witness(self):
        with self.assertRaises(ArithmeticError):
            replay.verify_integer_relations(self.curve, self.points[::-1], self.targets, self.rows)

    def test_missing_target_is_not_silently_zipped_away(self):
        with self.assertRaisesRegex(ArithmeticError, 'exact dimensions'):
            replay.verify_integer_relations(self.curve, self.points, self.targets[:-1], self.rows)

    def test_retained_primes_replay_without_prime_selection(self):
        record = self.helper['finite_reduction_certificate'](self.curve, self.points)
        with patch('mod2_reduction_independence.find_mod2_reduction_certificate', side_effect=AssertionError('prime selection forbidden')):
            self.assertEqual(replay.verify_retained_independence(self.curve, self.points, record, self.helper), record)

    def test_changed_finite_signature_is_rejected(self):
        record = self.helper['finite_reduction_certificate'](self.curve, self.points)
        record['signatures'][0]['rows'][0][0] ^= 1
        with self.assertRaisesRegex(ArithmeticError, 'signature differs'):
            replay.verify_retained_independence(self.curve, self.points, record, self.helper)

    def test_duplicate_finite_prime_is_rejected(self):
        record = self.helper['finite_reduction_certificate'](self.curve, self.points)
        record['certificate_primes'].append(record['certificate_primes'][0])
        record['signatures'].append(deepcopy(record['signatures'][0]))
        with self.assertRaisesRegex(ArithmeticError, 'prime roster'):
            replay.verify_retained_independence(self.curve, self.points, record, self.helper)

    def test_unverified_independence_rank_is_rejected(self):
        record = self.helper['finite_reduction_certificate'](self.curve, self.points)
        record['combined_exact_rank_over_F2'] += 1
        with self.assertRaisesRegex(ArithmeticError, 'do not prove'):
            replay.verify_retained_independence(self.curve, self.points, record, self.helper)


if __name__ == '__main__':
    unittest.main()
