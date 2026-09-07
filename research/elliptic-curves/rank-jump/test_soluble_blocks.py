import unittest
from fractions import Fraction as F
import retrospective as r
from cubic_bridge import Cubic
import quadric_rulings as qr
import linear_twist_blocks as lb


class SolubleBlockTests(unittest.TestCase):
    def test_quadratic_arithmetic(self):
        qr.Quad.B=F(2)
        h=qr.Quad(0,1)
        self.assertEqual(h*h,2)
        self.assertEqual((1+h)/(1-h),-3-2*h)

    def test_small_ruling_example(self):
        # (1,1) on y^2=x^3-2x+2; no high-rank data.
        qr.Quad.B=F(2)
        K=Cubic(-2,2)
        for sign in (-1,1):
            result=qr.section(K,F(1),F(1),sign)
            self.assertTrue(result["point_identity"])
            self.assertTrue(result["non_invariant_section"])
        self.assertEqual(len(qr.quadric_isomorphism(K,F(1),F(1))),4)

    def test_linear_twist_identity(self):
        for A,B,a in [(F(-2),F(2),F(1)),(F(0),F(3),F(2)),
                      (F(5,7),F(-11,3),F(4,9))]:
            lb.twist_identity(A,B,a)

    def test_common_shift_loses_at_most_one_direction(self):
        # Include shifts outside the k-dimensional span.
        for k in range(1,7):
            for eta in range(1<<(k+1)):
                rows=[(1<<i)^eta for i in range(k)]
                expected=k-int(eta<(1<<k) and eta.bit_count()%2==1)
                self.assertEqual(r.rank(rows),expected)


if __name__=="__main__":
    unittest.main()
