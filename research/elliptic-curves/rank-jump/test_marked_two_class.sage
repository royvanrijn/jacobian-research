#!/usr/bin/env sage-python
"""Small exact theorem regressions; no high-rank discovery or shell search."""
import json
from pathlib import Path
import sys
import unittest
from sage.all import QQ, ZZ, PolynomialRing, matrix, GF, prime_range

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from test_shared_value_block import SharedValueTests
from test_oblique_split_cubic import ObliqueBlockTests

R=PolynomialRing(QQ,'t');t=R.gen()


def common_genus(d,e,control):
    d,e=R(d),R(e)
    for f in [d,e]:
        if f.degree() not in [1,2] or f.gcd(f.derivative()).degree()!=0:
            raise ValueError('not a squarefree genus-zero double cover')
        if not f(control) or not f(control).is_square():
            raise ValueError('missing nonzero rational lift')
    if d.monic()==e.monic():
        if not (d.leading_coefficient()/e.leading_coefficient()).is_square():
            raise ValueError('lost constant twist')
        return 0,'equal'
    b=d.degree()+e.degree()-d.gcd(e).degree()+int(d.degree()%2 or e.degree()%2)
    assert b in [3,4]
    return int(b-3),'independent'


def finite_rank(c,n):
    m=(n*n-c*c-2)//2
    P=[(0,c),(-1,n),(2,n)]
    f=t**3+m*t*t-(m+3)*t+c*c
    rows=[];torsion_prime=None
    for p in prime_range(3,200):
        if f.discriminant()%p==0:continue
        roots=[a for a in range(int(p)) if f(a)%p==0]
        if not roots and torsion_prime is None:torsion_prime=int(p)
        for a in roots:
            bits=[]
            for x,y in P:
                v=(x-a)%p
                if not v:v=f.derivative()(a)%p
                bits.append(int(pow(int(v),(int(p)-1)//2,int(p))==p-1))
            rows.append(bits)
    assert torsion_prime is not None
    return matrix(GF(2),rows).rank()


class PropagationTests(unittest.TestCase):
    def test_equal_extension_keeps_square_constant(self):
        self.assertEqual(common_genus(t,4*t,1),(0,'equal'))

    def test_nonsquare_constant_twist_cannot_be_ignored(self):
        with self.assertRaises(ValueError):common_genus(t,2*t,1)

    def test_shared_branch_infinity(self):
        self.assertEqual(common_genus(t,t+3,1),(0,'independent'))

    def test_shared_finite_branch(self):
        self.assertEqual(common_genus(t*(t+3),t*(t+8),1),(0,'independent'))

    def test_disjoint_branches_give_genus_one_only(self):
        self.assertEqual(common_genus(t*(t+3),(t+1)*(t+7),1),(1,'independent'))
        # No infinitude assertion is returned merely from these rational lifts.

    def test_repeated_branch_rejected(self):
        with self.assertRaises(ValueError):common_genus((t-1)**2,t+3,0)

    def test_four_prescribed_anchors(self):
        self.assertEqual([finite_rank(c,n) for c,n in [(1,5),(3,7),(3,5),(1,1)]],[3,3,2,1])

    def test_joint_local_corrections(self):
        # G has joint local image span(1,1). Both separate local projections
        # are full, but (1,0) cannot be corrected by one global generic word.
        G=matrix(GF(2),[[1],[1]])
        H=matrix(GF(2),[[1,1,0],[1,0,0]])
        self.assertEqual(2-(H.rank()-G.rank()),1)
        self.assertEqual((3-H.rank())-(1-G.rank()),1)
        self.assertEqual(G[0:1,:].rank(),H[0:1,:].rank())
        self.assertEqual(G[1:2,:].rank(),H[1:2,:].rank())

    def test_K3_shell_formula(self):
        for h in [0,1]:
            norm=4*h+10
            D_square=4*(-2)+4*(h+4)-norm
            D_O=2*(-2)+(h+4)
            self.assertEqual((D_square,D_O),(-2,h))


if __name__=='__main__':unittest.main()
