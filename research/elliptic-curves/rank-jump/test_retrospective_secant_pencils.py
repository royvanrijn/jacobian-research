import unittest
from fractions import Fraction as F
import retrospective as r
import retrospective_secant_pencils as pencils
import verify_retrospective_secant_pencils as verify


class SecantPencilTests(unittest.TestCase):
    def test_square_decision_is_exact(self):
        for value in [F(0),F(4,9),F(1),F(25,49)]:self.assertTrue(pencils.rational_square(value))
        for value in [F(-1),F(2),F(4,3),F(17,19)]:self.assertFalse(pencils.rational_square(value))

    def test_original_anchor_and_generic_pair_only(self):
        inp=r.read(pencils.INPUT);out=r.read(pencils.OUTPUT)
        for row in out['rows']:
            self.assertEqual(row['selected_original_generic_indices'],[0,1])
            source=next(x for x in inp['rows'] if x['case_index']==row['case_index'])
            A,B=map(F,source['short_model'][3:]);s=F(row['secant_slope']);b=F(row['secant_intercept'])
            for point in source['generic_pair']:
                x,y=map(F,point)
                self.assertEqual(y*y,x**3+A*x+B)
                self.assertEqual(y,s*x+b)
                # At t=0 the same abscissae are roots; at t=1 the original point returns.
                self.assertEqual(x**3+A*x+B-(s*x+b)**2,0)

    def test_apparent_jump_accounting_does_not_change_original_quotient(self):
        for row in r.read(pencils.OUTPUT)['rows']:
            n=row['original_known_independent_rank'];m=row['original_generic_rank'];g=row['base_arithmetic_generic_rank']
            self.assertEqual(n-g,(m-g)+row['original_observed_quotient_rank'])
            self.assertEqual(row['original_generic_quotient_contribution_of_selected_pair'],0)

    def test_polynomial_arithmetic_cancellation(self):
        self.assertEqual(verify.mul([F(-1),F(1)],[F(1),F(1)]),[-1,0,1])
        self.assertEqual(verify.add([F(1),F(2)],[-1,-2]),[0])


if __name__=='__main__':unittest.main()
