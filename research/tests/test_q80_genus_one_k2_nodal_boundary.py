"""Failure controls for the finite nodal-boundary arithmetic."""
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py'
SPEC=importlib.util.spec_from_file_location('q80_k2_boundary_test',PATH)
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)

class NodalBoundaryTests(unittest.TestCase):
    def test_scalar_is_retained_in_square_test(self):
        self.assertTrue(V.square_up_to_scalar([3,6,3]))
    def test_changed_coefficient_is_not_a_square(self):
        self.assertFalse(V.square_up_to_scalar([3,6,4]))
    def test_odd_order_at_zero_is_rejected(self):
        self.assertFalse(V.square_up_to_scalar([0,5,10,5]))
    def test_even_zero_and_degree_drop_are_kept(self):
        self.assertTrue(V.square_up_to_scalar([0,0,5]))
    def test_nodal_quadratic_has_no_rational_root(self):
        self.assertTrue(V.irreducible([88,62,1]))
    def test_split_quadratic_is_not_irreducible(self):
        self.assertFalse(V.irreducible(V.mul([-35%131,1],[-75%131,1])))
    def test_augmented_rank_detects_impossible_lift(self):
        self.assertEqual(V.rank([[1,0],[0,1],[1,1]]),2)
        self.assertEqual(V.rank([[1,0,0],[0,1,0],[1,1,1]]),3)

if __name__=='__main__':unittest.main()
