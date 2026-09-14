"""Arithmetic boundaries required for a branch-field exclusion witness."""
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_q80_rational_branch_torsion.py'
SPEC = importlib.util.spec_from_file_location('branch_torsion',PATH)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class BranchFieldWitnessTests(unittest.TestCase):
    def test_irreducible_field_with_split_good_reduction(self):
        # Q(i), prime above101 with i=10; x^3+x+1 is irreducible mod101.
        V.check_witness([1,0,1],[1],[1],101,10)

    def test_actual_torsion_cannot_be_excluded(self):
        with self.assertRaisesRegex(ValueError,'residue root'):
            V.check_witness([1,0,1],[1],[0],101,10)

    def test_wrong_branch_residue_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'branch root'):
            V.check_witness([1,0,1],[1],[1],101,11)

    def test_ramified_or_nonmonic_reduction_is_deferred(self):
        for q,t in [([101,0,1],0),([1,1,101],100)]:
            with self.assertRaisesRegex(ValueError,'leading unit'):
                V.check_witness(q,[1],[1],101,t)

    def test_singular_fibre_is_not_a_witness(self):
        with self.assertRaisesRegex(ValueError,'singular parent'):
            V.check_witness([1,0,1],[0],[0],101,10)


if __name__ == '__main__':
    unittest.main()
