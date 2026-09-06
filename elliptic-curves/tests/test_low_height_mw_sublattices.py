"""Small arithmetic failure regressions, not research-success assertions."""
import json
from math import comb
import unittest
from search_low_height_mw_sublattices import shell, parity, run_gp


class CombinationTests(unittest.TestCase):
    def test_shell_is_complete_unoriented_small_box(self):
        for n in (3,5):
            vectors=list(shell(n))
            self.assertEqual(len(vectors),sum(comb(n,k)*2**(k-1) for k in range(1,4)))
            self.assertEqual(len(set(map(tuple,vectors))),len(vectors))
            self.assertTrue(all(next(a for a in v if a)==1 for v in vectors))
            self.assertEqual(parity(vectors)['largest_coset'],4)

    def test_common_cover_does_not_mean_primitive_subgroup(self):
        # Three independent points in one parity class generate index four.
        lines=run_gp('B=[1,0,0;1,2,0;1,0,2];print(abs(vecprod(matsnf(B))));')
        self.assertEqual(lines,['4'])

    def test_bounded_qfminim_second_field_is_not_the_minimum(self):
        lines=run_gp('G=[2.,0;0,3.];Q=qfminim(G,3,100,2);print(Q[2]);print(vecmin(vector(matsize(Q[3])[2],j,Q[3][,j]~*G*Q[3][,j])));')
        self.assertEqual(float(lines[0]),3.)
        self.assertEqual(float(lines[1]),2.)

    def test_full_rank_integer_kernel_is_empty(self):
        self.assertEqual(run_gp('print(matkerint(matid(3)));'),['[;]'])


if __name__=='__main__':
    unittest.main()
