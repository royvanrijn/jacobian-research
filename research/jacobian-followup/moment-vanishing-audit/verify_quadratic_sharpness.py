#!/usr/bin/env python3
"""Exact small-rank Gaussian witnesses for the credited all-rank theorem.

No source producer code is imported. The general existence proof is written.
"""
from pathlib import Path
from math import factorial
from fractions import Fraction
import hashlib,json
import sympy as sp


def compositions(total,parts):
    if parts==1:
        yield (total,);return
    for i in range(total+1):
        for tail in compositions(total-i,parts-1):yield (i,)+tail


def real_moment(e):
    return 0 if e%2 else factorial(e)//(2**(e//2)*factorial(e//2))


def main():
    root=sp.sqrt(3);K=sp.QQ.algebraic_field(root);d=K.from_sympy((3+root)/6)
    q=lambda n,m=1:K.from_sympy(sp.Rational(n,m))
    h0=d-q(1,2);h1=q(1,3)-d;h2=d/q(2)-q(1,6)
    determinant=h0*h2-h1*h1
    assert determinant==q(1,18)-d/q(4) and determinant!=K.zero
    alpha=h1/h0;beta2=determinant/(h0*h0);gamma2=h0
    # Ordinary symmetric Ttilde has trace -1, determinant delta and
    # l^T adj(I-t*Ttilde) l = h0-t/6.
    assert -alpha*(1+alpha)-beta2==d
    assert gamma2*(1+alpha)==-q(1,6)
    # HT=T^T H for the cyclic Hankel realization.
    assert -d*h0-h1==h2
    rows=[]
    for m in range(1,7):
        value=K.zero
        for n1,nb,n2,ng,n0 in compositions(m,5):
            e1=2*n1+nb+ng;e2=nb+2*n2
            if e1%2 or e2%2:continue
            assert nb%2==ng%2==0
            multinomial=factorial(m)
            for n in (n1,nb,n2,ng,n0):multinomial//=factorial(n)
            term=q(multinomial*real_moment(e1)*real_moment(e2))
            term*= (alpha/q(2))**n1 * beta2**(nb//2)
            term*= (-(1+alpha)/q(2))**n2 * gamma2**(ng//2)*q(1,2)**n0
            value+=term
        if m<=4:assert value==K.zero
        if m==5:assert value==q(2)*(5*d-1) and value!=K.zero
        rows.append(str(sp.expand(K.to_sympy(value))))
    rankone=[]
    for m in range(1,5):
        value=0
        for n2,ng,n0 in compositions(m,3):
            e=2*n2+ng
            if e%2:continue
            coefficient=factorial(m)//(factorial(n2)*factorial(ng)*factorial(n0))
            value+=coefficient*(-2)**(ng//2)*(-1)**n0*real_moment(e)
        rankone.append(value)
    assert rankone==[0,0,-4,-48]
    rec={'scope':'Exact ranks one and two. All-rank sharpness uses the written projective-gap and symmetric-realization argument.',
         'rank_one_moments':rankone,'rank_two_delta':'(3+sqrt(3))/6',
         'rank_two_moments_through_six':rows,'first_four_zero':True,
         'fifth_equals':'2*(5*delta-1), nonzero in Q(delta)',
         'ordinary_symmetric_realization':'trace, determinant, linear transfer numerator and nonsingular Hankel form verified',
         'sympy_version':sp.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('quadratic_sharpness_verification.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
