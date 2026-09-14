"""Reject incomplete counts and unsupported rank conclusions in this proof."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_common_quartic_control_rank_gate.py'
SPEC = importlib.util.spec_from_file_location('control_rank_gate', PATH)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class ControlParentRankTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        packet = V.read(V.DEFAULT/'input.json')
        cls.rows = [(control, V.read(V.DEFAULT/(control['name']+'-checkpoint.json')), p)
                    for control, p in zip(packet['controls'], packet['primes'])]

    def test_both_complete_surface_counts(self):
        bounds = [V.check_control(*row)['any_rational_fibration_mw_rank_upper']
                  for row in self.rows]
        self.assertEqual(bounds, [10, 8])

    def test_missing_infinity_fibre_is_rejected(self):
        control, record, p = deepcopy(self.rows[0])
        record['fibre_counts_finite_then_infinity'].pop()
        with self.assertRaisesRegex(ValueError, 'including infinity'):
            V.check_control(control, record, p)

    def test_changed_fibre_count_is_rejected(self):
        control, record, p = deepcopy(self.rows[1])
        record['fibre_counts_finite_then_infinity'][0] += 1
        with self.assertRaisesRegex(ValueError, 'complete fibre counts'):
            V.check_control(control, record, p)

    def test_bad_reduction_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'good infinity'):
            V.squarefree_degree24([1]+[0]*24, 5)
        with self.assertRaisesRegex(ValueError, 'non-squarefree'):
            V.squarefree_degree24([0]*24+[1], 5)

    def test_false_upper_bound_is_rejected(self):
        control, record, p = deepcopy(self.rows[0])
        record['any_rational_fibration_mw_rank_upper'] = 9
        with self.assertRaisesRegex(ValueError, 'trace or rank bound'):
            V.check_control(control, record, p)

    def test_exact_or_geometric_promotion_is_rejected(self):
        for field in ['exact_inherited_rank', 'geometric_picard_rank']:
            with self.subTest(field=field):
                control, record, p = deepcopy(self.rows[0])
                record[field] = record['rational_picard_rank_upper']
                with self.assertRaisesRegex(ValueError, 'promoted to exact or geometric'):
                    V.check_control(control, record, p)


if __name__ == '__main__':
    unittest.main()
