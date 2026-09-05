import unittest
from fractions import Fraction as F
import retrospective as r
import local_collision as lc
import affine_selmer as af
from affine_selmer_analysis import extension_rank


class AffineSelmerTests(unittest.TestCase):
    def test_joint_inconsistency_with_each_equation_soluble(self):
        # x=0 and x=1: local existence need not give one common correction.
        equations=[1,3]
        self.assertTrue(all(af.affine_solve([row],1)["consistent"] for row in equations))
        result=af.affine_solve(equations,1)
        self.assertFalse(result["consistent"])
        self.assertEqual(lc.lift(result["inconsistent_row_combination"],equations),2)

    def test_affine_solution(self):
        result=af.affine_solve([0b1101,0b0110],3)
        self.assertTrue(result["consistent"])
        x=result["particular_anchor_mask"]
        self.assertEqual((x&0b101).bit_count()%2,1)
        self.assertEqual((x&0b110).bit_count()%2,0)

    def test_local_square(self):
        self.assertTrue(af.local_square(F(1,9),2))
        self.assertFalse(af.local_square(F(5),2))
        self.assertTrue(af.local_square(F(25,9),7))
        self.assertFalse(af.local_square(F(7),7))

    def test_alternating_extension_rank_from_radical_only(self):
        # All alternating 4x4 forms and all possible new columns.
        for mask in range(64):
            matrix=[[0]*4 for _ in range(4)]
            k=0
            for i in range(4):
                for j in range(i+1,4):
                    matrix[i][j]=matrix[j][i]=(mask>>k)&1;k+=1
            radical=lc.orthogonal(map(r.pack,matrix),4)
            for column in range(16):
                pairings=[(column&v).bit_count()%2 for v in radical]
                extended=[r.pack(row)|(((column>>i)&1)<<4) for i,row in enumerate(matrix)]+[column]
                self.assertEqual(extension_rank(matrix,pairings),r.rank(extended))


if __name__=="__main__":
    unittest.main()
