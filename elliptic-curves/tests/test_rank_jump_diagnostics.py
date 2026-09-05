"""Independent small arithmetic checks for the diagnostic-only additions."""
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from compare_bounded_prime_selectors import short_mod,trace,correlation,contributions
from audit_mw18_retained_visibility import transport


class DiagnosticTests(unittest.TestCase):
    def test_general_model_count_agrees_with_completed_square(self):
        # Direct enumeration on the general equation independently tests both
        # the x/y translations and the character-sum count used by selectors.
        model=(1,-1,1,-7,10)
        for p in (5,7,11,13,17,19):
            a,b=short_mod(model,p)
            if (4*a**3+27*b*b)%p==0:continue
            chi=[0]+[-1]*(p-1)
            for y in range(1,p):chi[y*y%p]=1
            count=1+sum((y*y+model[0]*x*y+model[2]*y-x**3-model[1]*x*x-model[3]*x-model[4])%p==0
                        for x in range(p) for y in range(p))
            self.assertEqual(trace(a,b,p,chi),p+1-count)

    def test_bad_primes_do_not_become_good_prime_scores(self):
        self.assertIsNone(trace(0,0,7,[0,1,1,-1,1,-1,-1]))
        self.assertEqual(contributions(None,7),(0,0))

    def test_ties_and_constant_outcomes(self):
        self.assertEqual(correlation([2,2,4],[1,1,5]),1)
        self.assertEqual(correlation([2,2,4],[5,5,1]),-1)
        self.assertIsNone(correlation([1,2,3],[0,0,0]))

    def test_public_to_native_scale_and_reject_twist(self):
        self.assertEqual(transport([0,0,0,-1,1],[(0,1)],[0,0,0,'-1/16','1/64']),[{'x':'0','y':'1/8'}])
        with self.assertRaises(ArithmeticError):
            transport([0,0,0,-1,1],[(0,1)],[0,0,0,'-1/16','1/32'])


if __name__=='__main__':unittest.main()
