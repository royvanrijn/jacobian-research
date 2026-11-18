#!/usr/bin/env python3
"""Exact spin, real-spectrum and symmetric-normal-form controls.

Uses sympy==1.14.0 for polynomial identities and exact Sturm sequences.
The uniform quantifiers are supplied by the accompanying written proofs.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

def clifford_multiply(left,right):
    out={}
    for a,x in left.items():
        for b,y in right.items():
            inversions=sum((b & ((1<<i)-1)).bit_count() for i in range(4) if (a>>i)&1)
            mask=a^b
            out[mask]=out.get(mask,Fraction(0))+(-1)**inversions*x*y
    return {mask:value for mask,value in out.items() if value}
def clifford_scale(value,scalar):
    return {mask:scalar*coefficient for mask,coefficient in value.items()}
def sturm_real_count(poly,z):
    seq=sp.sturm(poly,z)
    def variations(infinity_sign):
        signs=[]
        for value in seq:
            q=sp.Poly(value,z)
            sign=sp.sign(q.LC())*infinity_sign**q.degree()
            if sign:signs.append(sign)
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    return variations(-1)-variations(1)
def parity(permutation):
    return sum(permutation[i]>permutation[j]
               for i in range(len(permutation)) for j in range(i+1,len(permutation)))%2

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    if sp.__version__!="1.14.0":raise RuntimeError("Requires pinned sympy==1.14.0")
    basis=[{1<<i:Fraction(1)} for i in range(4)]
    def difference(i,j):return {1<<i:Fraction(1),1<<j:Fraction(-1)}
    d=clifford_scale(clifford_multiply(difference(0,1),difference(2,3)),Fraction(1,2))
    h=clifford_scale(clifford_multiply(difference(0,2),difference(1,3)),Fraction(1,2))
    minus_one={0:Fraction(-1)}
    assert clifford_multiply(d,d)==minus_one
    assert clifford_multiply(h,h)==minus_one
    assert clifford_multiply(h,d)==clifford_scale(clifford_multiply(d,h),-1)
    assert clifford_multiply(clifford_multiply(h,d),clifford_scale(h,-1))==clifford_scale(d,-1)
    for lift,permutation in [(d,[1,0,3,2]),(h,[2,3,0,1])]:
        inverse=clifford_scale(lift,-1)
        for i in range(4):
            assert clifford_multiply(clifford_multiply(lift,basis[i]),inverse)==basis[permutation[i]]
    # Admissible long cycles have odd-square length. The block argument in
    # the note proves the all-m assertion; these are independent finite controls.
    for m in (3,5,7,9,11,13):
        r=m*m
        for a in range(r):
            if sp.gcd(a,r)==1:
                assert parity([(a*j)%r for j in range(r)])==0
    # A nonadmissible eleven-cycle has a power conjugator of odd parity.
    assert parity([(2*j)%11 for j in range(11)])==1

    z,u,s,x,k=sp.symbols("z u s x k")
    dcoef=(47+7*x**2)/168
    ecoef=(379+357*x**2+196*x**3)/21168
    a=-(7*x**2-3)**2/87808
    hreal=z**8-z**6+dcoef*z**4-ecoef*z**2
    breal=-(6912-14112*x**2+5488*x**3+8575*x**4
             -7203*x**5+31213*x**6+16807*x**7)/533655864
    constant=((7*x**2-3)**2*(4-7*x+7*x**2)**2
              *(48+56*x+35*x**2)**3/379718108242115328)
    qreal=s**4+breal*s**2+constant
    critical=7*u**4-5*u**3+3*dcoef*u**2-ecoef*u-a
    assert sp.expand(sp.discriminant(hreal-s*z+a,z)+7**7*qreal**2)==0
    expected=x**2*(x-2)**6*(x+1)**2*(48+56*x+35*x**2)/1259712
    assert sp.factor(sp.discriminant(critical,u)-expected)==0
    assert sp.factor(breal**2-4*constant
                     +x**6*(x-2)**5*(x+1)**3/10806531246)==0
    seed_h=sp.expand(hreal.subs(x,1));seed_a=a.subs(x,1);seed_q=qreal.subs(x,1)
    test_u=[0,sp.Rational(1,100),sp.Rational(1,10),sp.Rational(2,5),sp.Rational(4,5)]
    assert [sp.sign(critical.subs({x:1,u:v}))for v in test_u]==[1,-1,1,-1,1]
    samples8=[sp.Integer(0),sp.Rational(1,150),sp.Integer(1)]
    assert [sturm_real_count(seed_h-v*z+seed_a,z)for v in samples8]==[2,6,2]
    assert [sp.sign(seed_q.subs(s,v))for v in samples8]==[1,-1,1]
    assert sturm_real_count(seed_q,s)==4

    # The infinite sequence is x=1+11/(1+121k), k>=1. Every denominator
    # is a unit modulo11 and x reduces to1, preserving the original seed's
    # arithmetic Frobenius certificate before the complex coordinate twist.
    seq=1+sp.Rational(11,1)/(1+121*k)
    assert sp.factor(seq-1)==11/(121*k+1)
    for coeff in (dcoef,ecoef,a):
        num,den=sp.fraction(sp.cancel(coeff.subs(x,seq)-coeff.subs(x,1)))
        assert all(int(c)%11==0 for c in sp.Poly(num,k,domain=sp.QQ).all_coeffs())
        assert int(sp.Poly(den,k).eval(0))%11!=0
        assert all(int(c)%11==0 for c in sp.Poly(den,k).all_coeffs()[:-1])
    assert sp.diff(ecoef,x)==(714*x+588*x**2)/21168

    h10=9*z**10-864*z**8+28161*z**6-350150*z**4+1127500*z**2
    w10=59049*s**4-126938174857200000*s**2+68130645548641000000000000000
    assert [sturm_real_count(h10-v*z+250000,z)for v in (1,1040000,1100000)]==[6,6,2]
    # Positive branch roots lie in the three intervals represented by these samples.
    assert [sp.Poly(w10,s).count_roots(0,v) for v in (1,1040000,1100000)]==[0,1,2]

    # Uniform symmetric normal form, tested symbolically in all eleven
    # free G coefficients of the first new admissible degree n=26.
    coeff=sp.symbols("a0:11")
    G=u**11+sum(coeff[j]*u**j for j in range(11))
    N=(u-1)**2*G
    S=(3*u+1)*G+2*u*(u-1)*sp.diff(G,u)
    assert sp.expand(2*u*sp.diff(N,u)-N-(u-1)*S)==0
    assert sp.Poly(S,u).degree()==12 and sp.Poly(S,u).LC()==25
    assert S.subs(u,0)==coeff[0]
    assert sp.expand(S.subs(u,1)-4*G.subs(u,1))==0

    # Negative controls distinguish a false commuting-lift assumption and
    # a degenerate boundary member of the real family.
    assert clifford_multiply(h,d)!=clifford_multiply(d,h)
    assert sp.factor(expected.subs(x,2))==0
    assert sturm_real_count(seed_h+seed_a,z)!=8
    report={
        "schema_version":1,"kind":"internal_exact_controls","python":platform.python_version(),
        "sympy":sp.__version__,"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic":"Q Clifford algebra; Q[x,s]; exact rational Sturm sequences",
        "spin":{"double_transposition_lift_square":-1,"commuting_elements_have_anticommuting_lifts":True,
                "conjugation_changes_lift_sign":True},
        "odd_square_cycle_controls":[3,5,7,9,11,13],
        "real_degree8_seed":{"H":str(seed_h),"a":str(seed_a),"Q":str(seed_q),
                            "critical_signs":[1,-1,1,-1,1],
                            "regular_real_root_count_samples":[2,6,2]},
        "infinite_real_family":"c=-1,t=-(1+11/(1+121k)), integer k>=1",
        "degree10_real_root_count_samples":[6,6,2],
        "degree26_shape":{"degree_G":11,"degree_S":12,"LC_S":25,
                         "generic_coefficient_derivative_identity":True},
        "uniform_proof_limits":"All-degree spin, real-root-count and normal-form quantifiers use the written proof. The checker does not produce a degree26 paired-critical model or prove uniform braid transitivity.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "result":"PASS"
    }
    out=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(out)
    print(out,end="")
if __name__=="__main__":main()
