#!/usr/bin/env python3
"""Universal invariants and independent checks of the norm-square rank gate."""
from math import isqrt
import retrospective as r
import norm_square_rank as gate


def verify():
    from sage.all import QQ,PolynomialRing,EllipticCurve,RealIntervalField
    R=PolynomialRing(QQ,["a","b"]);F=R.fraction_field();A,B=map(F,R.gens())
    U=PolynomialRing(F,"u");u=U.gen()
    D=1+A*u*u+B*u**3;delta=-4*A**3-27*B*B
    a2=2*A*u*D
    a4=(A+3*B*u+A*A*u*u)*D**2
    a6=(B+A*B*u*u-B*B*u**3)*D**3
    b2,b4,b6=4*a2,2*a4,4*a6;b8=4*a2*a6-a4*a4
    discriminant=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    c4=b2*b2-24*b4
    assert discriminant==16*delta*D**8
    assert c4==16*D**2*(A*A*u*u-9*B*u-3*A)
    assert D.discriminant()==delta
    assert D.gcd(c4//D**2)==1
    assert [a2.degree(),a4.degree(),a6.degree()]==[4,8,12]
    assert discriminant.degree()==24
    assert discriminant[24]==16*delta*B**8
    assert (c4//D**2).resultant(D)==4096*delta**2
    # Degenerate leading c4 coefficient A=0: it still has no common root with D.
    V=PolynomialRing(QQ,"v");v=V.gen()
    assert (1+v**3).gcd(-9*v)==1

    data=r.read(gate.OUTPUT);control=data["control"]
    a,b=QQ(control["A"]),QQ(control["B"]);q=QQ(control["u"])
    assert q==-a/b and q!=0 and 1+a*q*q+b*q**3==1
    assert isqrt(int(b))**2!=b
    f=V([b,a,0,1])
    assert f.is_irreducible()
    disc=f.discriminant()
    assert disc>0 and isqrt(int(disc))**2!=disc
    roots=f.roots(RealIntervalField(128),multiplicities=False)
    assert len(roots)==3
    assert sum(x<0 for x in roots)==1 and sum(x>0 for x in roots)==2
    signs=[1-q*x for x in roots]
    assert sum(x<0 for x in signs)==2 and sum(x>0 for x in signs)==1
    base=EllipticCurve([0,a,0,0,b*b]);Q=base([0,b])
    assert list((2*Q)[:2])==[-a,-b]
    assert list(map(str,base.a_invariants()))==control["auxiliary_curve_ainvariants"]
    print("PASS universal D-twist invariants, good infinity, and coefficient-only S4 control")


if __name__=="__main__":
    verify()
