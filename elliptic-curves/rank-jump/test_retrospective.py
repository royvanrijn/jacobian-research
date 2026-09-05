"""Independent small-model checks for mathematical failure modes, not searches."""
import itertools
import sys
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/'cas'))
import retrospective as r
import blocks
from mod2_reduction_independence import finite_add,finite_curve_points


class ArithmeticTests(unittest.TestCase):
    def test_characters_are_homomorphisms_including_two_torsion(self):
        for p in (5,7,11,13,17):
            for A,B in ((-1,0),(-2,1),(1,1)):
                model=['0','0','0',str(A),str(B)]
                roots=r.roots_at(str(A),str(B),p)
                if not roots:continue
                pts=finite_curve_points(A%p,B%p,p)
                def character(P):
                    return 0 if P is None else r.point_signature(model,P,[(p,roots)])
                for P in pts:
                    for Q in pts:
                        self.assertEqual(character(finite_add(P,Q,A%p,p)),character(P)^character(Q))

    def test_nonintegral_point_reduces_to_identity(self):
        # Twice (3,5) on y^2=x^3-2 is (129/100,-383/1000).
        # At the good prime 5 it reduces to O although the model has rational 2-torsion mod 5.
        model=['0','0','0','0','-2'];P=['129/100','-383/1000']
        roots=r.roots_at(model[3],model[4],5)
        self.assertEqual(roots,(3,))
        r.short(model,[P])
        self.assertEqual(r.point_signature(model,P,[(5,roots)]),0)

    def test_smith_forms_against_determinantal_divisors(self):
        # For 2x2 full-rank matrices, d1=gcd(entries), d1*d2=abs(det).
        from math import gcd
        for a,b,c,d in itertools.product(range(-2,3),repeat=4):
            det=a*d-b*c
            if det:
                first=gcd(gcd(a,b),gcd(c,d))
                self.assertEqual(blocks.smith_diagonal([[a,b],[c,d]]),[first,abs(det)//first])
        self.assertEqual(blocks.smith_diagonal([[2,0],[0,12],[0,0]]),[2,12])

    def test_alternating_blocks_exhaustive_through_dimension_four(self):
        for n in range(1,5):
            pairs=[(i,j) for i in range(n) for j in range(i)]
            for mask in range(1<<len(pairs)):
                matrix=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(pairs):matrix[i][j]=matrix[j][i]=(mask>>k)&1
                out=blocks.ct_normal_form(matrix)
                null_count=sum(all((r.pack(row)&v).bit_count()%2==0 for row in matrix) for v in range(1<<n))
                self.assertEqual(null_count,1<<out['restricted_radical_dimension'])
                self.assertEqual(out['pairing_rank'],r.rank(map(r.pack,matrix)))

    def test_finite_alias_does_not_prove_global_dependence(self):
        # Projection F2^3 -> F2^2 aliases two genuinely independent ambient vectors.
        ambient=[0b001,0b101]
        self.assertEqual(r.rank(ambient),2)
        self.assertEqual(r.rank([x&3 for x in ambient]),1)


class CoverSquareTests(unittest.TestCase):
    def test_exact_polynomial_squareclass_gate(self):
        import cover_experiment as c
        for coeff in itertools.product(range(-1,2),repeat=4):
            p=c.trim(list(map(r.F,coeff)));square=c.mul(p,p)
            root=c.polynomial_square_root(square)
            self.assertIsNotNone(root)
            self.assertEqual(c.mul(root,root),square)
        self.assertIsNone(c.polynomial_square_root(list(map(r.F,[1,1,1]))))
        # Equal squareclasses can have different degrees and nonconstant square ratios.
        f=list(map(r.F,[1,0,1]));g=c.mul(f,list(map(r.F,[1,2,1])))
        self.assertIsNotNone(c.polynomial_square_root(c.mul(f,g)))


if __name__=='__main__':unittest.main()
