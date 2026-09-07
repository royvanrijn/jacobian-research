from itertools import product
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from exact_parity_ellipsoid import enumerate_coset


class ExactParityEllipsoidTests(unittest.TestCase):
    def test_skew_gram_against_complete_integer_box(self):
        G=[[4,2],[2,4]]
        for r in product(range(2),repeat=2):
            expected={}
            # Minimum eigenvalue2 bounds every coordinate by sqrt(10/2)<3.
            for z in product(range(-3,4),repeat=2):
                if any((z[i]-r[i])%2 for i in range(2)):continue
                q=sum(z[i]*G[i][j]*z[j] for i in range(2) for j in range(2))
                if q<=10:expected[str(q)]=expected.get(str(q),0)+1
            self.assertEqual(enumerate_coset(G,r,10)['norm_counts'],expected)

    def test_empty_coset_and_nonpositive_input(self):
        self.assertEqual(enumerate_coset([[4]],[1],2)['norm_counts'],{})
        with self.assertRaises(ValueError):enumerate_coset([[0]],[0],1)


if __name__=='__main__':unittest.main()
