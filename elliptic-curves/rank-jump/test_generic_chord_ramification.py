"""Regressions for support cancellation and composite parity witnesses."""
import unittest
from math import gcd
from generic_chord_ramification import exponents,nullspace
from generic_chord_ramification_v2 import insert


class ChordRamificationTests(unittest.TestCase):
    def test_three_distinct_roots_cancel_despite_odd_norm_parity(self):
        # Local control only: f=X^3-2 at p=31, not a global soluble block.
        roots=[4,7,20];p=31
        norms=[2-x**3 for x in roots]
        self.assertTrue(all(n%p==0 and n%(p*p)!=0 for n in norms))
        self.assertTrue(all((x**3-2)%p==0 for x in roots))
        self.assertEqual(1^1^1,1)
        self.assertEqual((7^1)^(7^2)^(7^4),0)
        self.assertEqual(nullspace([6,5,3]),[7])

    def test_coprime_refinement_preserves_integer_norms(self):
        basis=[]
        for n in [12,18,72,50,49,343]:insert(basis,n)
        self.assertTrue(all(gcd(a,b)==1 for i,a in enumerate(basis) for b in basis[i+1:]))
        for n in [12,18,72,50,49,343]:
            product=1
            for j,e in exponents(n,basis).items():product*=basis[j]**e
            self.assertEqual(product,n)

    def test_private_coordinate_blocks_all_cancellation(self):
        # Add a different private coordinate to each of the three local rows.
        self.assertEqual(nullspace([6|8,5|16,3|32]),[])


if __name__=='__main__':unittest.main()
