"""Independent checks of field-degree hypotheses and halving quartics."""
import unittest
import retrospective as r
import local_collision as lc
import halving_fields as hf


class HalvingFieldTests(unittest.TestCase):
    def test_subspace_generator_is_complete_in_dimension_four(self):
        actual={tuple(x) for x in hf.subspaces(4)}
        generated={tuple(lc.canonical(x)) for x in __import__('itertools').product(range(16),repeat=4)}
        self.assertEqual(actual,generated)
        self.assertEqual(len(actual),67)

    def test_degree_depends_on_mod2_dimension_not_free_rank(self):
        # A formally independent list P,2Q contributes only one mod-2
        # direction; an integer rank count cannot replace the hypothesis.
        classes=[1,0];d=r.rank(classes)
        self.assertEqual(4**d,4)
        self.assertNotEqual(4**d,4**len(classes))

    def test_independent_sage_discriminant_and_frobenius_checks(self):
        try:from sage.all import QQ,GF,PolynomialRing
        except ImportError:self.skipTest('optional independent Sage polynomial checks')
        R=PolynomialRing(QQ,'A,B,p');A,B,p=R.gens();S=PolynomialRing(R,'z');z=S.gen()
        h=z**4-4*p*z**3-2*A*z*z-(4*A*p+8*B)*z+A*A-4*B*p
        self.assertEqual(h.discriminant(),4096*(p**3+A*p+B)**2*(-4*A**3-27*B*B))
        out=r.read(hf.OUTPUT)
        for row in out['paired_quartics']:
            for pattern,w in row['frobenius_witnesses'].items():
                F=PolynomialRing(GF(w['prime']),'x');poly=F(w['reduction'])
                degrees=sorted(int(q.degree()) for q,e in poly.factor() for _ in range(e))
                self.assertEqual(degrees,[4] if pattern=='4' else [1,3])
        # Recover y_Q and verify both the curve and duplication identities
        # in a retained quartic quotient, independently of modular tests.
        row=out['paired_quartics'][0];T=PolynomialRing(QQ,'z');z=T.gen()
        H=T(list(map(QQ,row['halving_quartic_ascending'])));Q=T.quotient(H);zz=Q.gen()
        A,B=map(QQ,row['short_model'][3:]);p,q=map(QQ,row['point'])
        yy=(zz**3-3*p*zz**2-A*zz-A*p-2*B)/(2*q)
        self.assertEqual(yy**2,zz**3+A*zz+B)
        self.assertEqual((3*zz**2+A)**2-4*yy**2*(p+2*zz),0)


if __name__=='__main__':unittest.main()
