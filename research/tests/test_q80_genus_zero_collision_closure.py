"""Boundary controls for the finite arithmetic used in genus-zero closure."""
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_q80_genus_zero_collision_closure.py'
SPEC=importlib.util.spec_from_file_location('q80_collision_closure',PATH)
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)
COUNT_PATH=PATH.with_name('verify_q80_genus_zero_closure.py')
COUNT_SPEC=importlib.util.spec_from_file_location('q80_exact_count',COUNT_PATH)
COUNT=importlib.util.module_from_spec(COUNT_SPEC);COUNT_SPEC.loader.exec_module(COUNT)

class CollisionClosureTests(unittest.TestCase):
    def test_simple_contact_has_unique_hensel_direction(self):
        self.assertNotEqual(V.hensel_determinant([1],[0,1,0,0,0,0,1],0),0)
    def test_repeated_contact_is_not_claimed_rigid(self):
        self.assertEqual(V.hensel_determinant([1],[0,0,1,0,0,0,1],0),0)
    def test_infinity_contact_is_in_the_degree_bounds(self):
        self.assertNotEqual(V.hensel_determinant([1,0,0,0,0,0,0,0,1],[0,1,0,0,0,1],0),0)
    def test_actual_identity_required(self):
        V.point_identity([], [1], [1], [1])
        with self.assertRaises(ValueError):V.point_identity([], [2], [1], [1])
    def test_high_degree_point_does_not_enter_census(self):
        with self.assertRaises(ValueError):V.point_identity([0,0,0,0,0,1],[1],[1],[1])
    def test_distinct_repeated_roots_prevent_closure(self):
        with self.assertRaises(ValueError):V.repeated_labels([{'b':2,'root':19},{'b':2,'root':58}])
    def test_equal_repeated_labels_remain_single_root(self):
        self.assertEqual(V.repeated_labels([{'b':2,'root':58},{'b':2,'root':58}]),{'2':[58]})
    def test_exact_intervals_include_both_signs_and_fractional_centres(self):
        self.assertEqual(COUNT.exact_short_count([[2,-1],[-1,2]],2)['signed_norm_histogram'],{'0':1,'2':6})
    def test_nonpositive_lattice_is_rejected(self):
        with self.assertRaises(ValueError):COUNT.exact_short_count([[0]],4)

if __name__=='__main__':unittest.main()
