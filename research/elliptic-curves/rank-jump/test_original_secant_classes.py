import unittest
from fractions import Fraction as F
import retrospective as r
import original_secant_classes as s


class OriginalSecantTests(unittest.TestCase):
    def test_known_rational_intercept_is_detected(self):
        # The prior oblique positive control: all three points lie on y=5*x.
        model,pts=r.short([0,F(47,2),0,F(-3,2),1],[[-1,-5],[2,10],[F(1,2),F(5,2)]])
        A,B=map(F,model[3:])
        for i,j in [(0,1),(0,2),(1,2)]:
            row=s.intercept(A,B,pts[i],pts[j],1)
            self.assertEqual(row['C'],'1')
            self.assertTrue(row['rationally_soluble'])

    def test_signature_preserves_squareclasses_and_multiplies(self):
        primes=[3,5,7,11]
        a,b=F(-45,14),F(77,125)
        self.assertEqual(s.squareclass_signature(a,primes),s.squareclass_signature(a*F(35,33)**2,primes))
        self.assertEqual(s.squareclass_signature(a*b,primes),s.squareclass_signature(a,primes)^s.squareclass_signature(b,primes))

    def test_local_alias_is_not_global_equality(self):
        # All primes in this deliberately tiny dictionary miss the obstruction.
        self.assertEqual(s.squareclass_signature(F(1),[3]),s.squareclass_signature(F(13),[3]))
        self.assertFalse(s.rational_square(F(13)))

    def test_degenerate_secant_is_not_a_negative_solubility_claim(self):
        self.assertEqual(s.intercept(F(0),F(1),[0,1],[0,-1],1)['status'],'DEGENERATE')
        row=s.intercept(F(-1),F(1),[0,1],[1,1],1)
        self.assertEqual(row,{'status':'DEGENERATE','reason':'horizontal secant'})


if __name__=='__main__':unittest.main()
