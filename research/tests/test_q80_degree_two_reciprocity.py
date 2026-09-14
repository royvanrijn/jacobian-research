"""Failure controls for the quadratic norm and degeneration boundary."""
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_q80_degree_two_reciprocity.py'
SPEC=importlib.util.spec_from_file_location('degree_two_reciprocity',PATH)
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)

class DegreeTwoReciprocityTests(unittest.TestCase):
    def test_quadratic_norm_character(self):
        # i^2=-1; Norm(1+i)=2 is nonsquare at131.
        self.assertEqual(V.mul(131,131),130)
        self.assertFalse(V.square(132))
        self.assertTrue(V.square(2)) # A base-field element becomes square.
        self.assertEqual(V.mul(132,V.inv(132)),1)
    def test_missing_roots_fail(self):
        V.check_roster(130,0,[0,1,130])
        for roots in [[],[0],[0,1],[0,1,129]]:
            with self.assertRaises(ValueError):V.check_roster(130,0,roots)
    def test_nodal_fibre_is_retained(self):
        V.check_roster(128,2,[1,1,129])
        self.assertEqual(V.character_value(128,129,129),9)
        with self.assertRaises(ValueError):V.character_value(128,1,1)
    def test_two_torsion_uses_regularized_value(self):
        self.assertEqual(V.characters(130,0,[0],[(0,0)],1),[1])
        self.assertEqual(V.characters(130,0,[0],[(0,0)],2),[0])
    def test_pole_and_cancellation(self):
        self.assertEqual(V.character_value(0,1,None),1)
        self.assertEqual(V.rational_at([1,129,1],[1,129,1],1,4),1)
        self.assertIsNone(V.rational_at([1],[1,129,1],1,4))
    def test_false_section_rejected(self):
        with self.assertRaises(ValueError):V.characters(130,0,[0],[(2,1)],1)

if __name__=='__main__':unittest.main()
