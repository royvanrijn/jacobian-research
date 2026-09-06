#!/usr/bin/env python3
"""Replay the small scalar-cup control without class groups or norm solvers."""
from fractions import Fraction
from math import lcm
import retrospective as r
import scalar_cup as scalar


def valuation(n,p):
    assert n
    v=0
    while n%p==0:n//=p;v+=1
    return v


def hensel(coefficients,root,p,precision):
    def evaluate(cs,t):return sum(c*t**i for i,c in enumerate(cs))
    derivative=[i*c for i,c in enumerate(coefficients)][1:]
    modulus=p
    assert evaluate(coefficients,root)%p==0
    assert evaluate(derivative,root)%p
    for _ in range(1,precision):
        root=(root-evaluate(coefficients,root)*pow(evaluate(derivative,root),-1,modulus*p))%(modulus*p)
        modulus*=p
    assert evaluate(coefficients,root)%modulus==0
    return root


def local_value(coefficients,root,p,precision):
    cs=list(map(Fraction,coefficients));d=lcm(*(x.denominator for x in cs))
    modulus=p**precision
    value=sum(int(x*d)*pow(root,i,modulus) for i,x in enumerate(cs))%modulus
    assert value, "precision cannot determine valuation"
    v=valuation(value,p);vd=valuation(d,p)
    unit=(value//p**v)*pow(d//p**vd,-1,p)%p
    return v-vd,unit


def norm_value(a,b,root,imaginary,p,precision):
    cs=list(map(Fraction,a))+list(map(Fraction,b))
    d=lcm(*(x.denominator for x in cs));modulus=p**precision
    value=sum(int(cs[i]*d)*pow(root,i,modulus) for i in range(3))
    value+=imaginary*sum(int(cs[i+3]*d)*pow(root,i,modulus) for i in range(3))
    value%=modulus
    assert value,"precision cannot determine norm-witness valuation"
    return valuation(value,p)-valuation(d,p)


def verify(data=None):
    from sage.all import QQ,PolynomialRing,NumberField,RealIntervalField,pari
    if data is None:data=r.read(scalar.CONTROL)
    assert data["status"]=="PASS"
    assert data["bindings"]==scalar.bindings()
    fcoeff=data["polynomial_ascending"];R=PolynomialRing(QQ,"t");f=R(fcoeff)
    assert f.is_irreducible() and f.discriminant()==163**2
    K=R.quotient(f,"theta")
    field=NumberField(f,"a");nf=pari.nfinit(pari(f))
    for witness in data["norm_witnesses"]:
        beta=field(list(map(QQ,witness["beta"])))
        assert all(e(beta)>0 for e in field.embeddings(RealIntervalField(128)))
        for prime in data["S_finite"]:
            for P in pari.idealprimedec(nf,prime):
                assert pari.nfislocalpower(nf,P,pari.Mod(pari(R(witness["beta"])),pari(f)),2)==1
    p=5;precision=64
    roots=[hensel(fcoeff,t,p,precision) for t in range(p) if sum(c*t**i for i,c in enumerate(fcoeff))%p==0]
    imag=[hensel([1,0,1],t,p,precision) for t in (2,3)]
    assert len(roots)==3
    betas=[w["beta"] for w in data["norm_witnesses"]]
    local=[[local_value(beta,t,p,precision) for t in roots] for beta in betas]
    artin=[];observed=[]
    for beta in betas:
        row=[]
        for other in local:
            assert all(v%2==0 for v,u in other)
            row.append(sum((v//2)*int(pow(local_value(beta,t,p,precision)[1],2,p)==p-1)
                           for t,(v,u) in zip(roots,other))%2)
        artin.append(row)
    for witness in data["norm_witnesses"]:
        beta,a,b=[K(list(map(QQ,witness[key]))) for key in ("beta","norm_a","norm_b")]
        assert a*a+b*b==beta
        # All coordinates are integral away from 5; the norm is a 5-unit.
        # Thus no prime outside 5 can occur in the norm-witness ideal.
        assert abs(R(witness["beta"]).resultant(f))==625
        for key in ("beta","norm_a","norm_b"):
            for x in map(Fraction,witness[key]):
                assert x.denominator==5**valuation(x.denominator,5)
        parities=[]
        for t in roots:
            exponents=[norm_value(witness["norm_a"],witness["norm_b"],t,i,p,precision) for i in imag]
            assert exponents[0]%2==exponents[1]%2
            assert sum(exponents)==local_value(witness["beta"],t,p,precision)[0]
            parities.append(exponents[0]%2)
        observed.append([sum(e*int(pow(unit,2,p)==p-1) for e,(v,unit) in zip(parities,row))%2
                         for row in local])
    assert artin==data["Artin_matrix"]
    assert observed==data["independent_norm_cup_matrix"]==scalar.symmetrize(artin)
    print("PASS independent rational norm identities and 5-adic parity/Artin evaluations")
    return observed


if __name__=="__main__":
    verify()
