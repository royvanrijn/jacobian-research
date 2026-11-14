#!/usr/bin/env python3
"""Independent exact regressions; all-order proofs are in the adjacent notes.

No external producer code is imported. Uses rational polynomial arithmetic.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import factorial, comb
import hashlib
import json
from pathlib import Path


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for mon, c in poly.items(): out[mon] += c
    return {m:c for m,c in out.items() if c}


def scale(poly, c): return {m:c*v for m,v in poly.items() if c*v}


def mul(a, b):
    out = defaultdict(F)
    for ma, ca in a.items():
        for mb, cb in b.items(): out[tuple(x+y for x,y in zip(ma,mb))] += ca*cb
    return {m:c for m,c in out.items() if c}


def var(i, n): return {tuple(int(j == i) for j in range(n)):F(1)}


def sphere(poly, n):
    # First two complex coordinates of the unit sphere in C^n.
    return sum((c*F(factorial(a)*factorial(b)*factorial(n-1),
                    factorial(a+b+n-1))
                for (a,aa,b,bb),c in poly.items() if a == aa and b == bb), F(0))


def gaussian(poly):
    out=F(0)
    for (z,w,t),c in poly.items():
        if z == w and t % 2 == 0:
            out += c*factorial(z)*F(factorial(t),2**(t//2)*factorial(t//2))
    return out


def main():
    z0,b0,z1,b1=[var(i,4) for i in range(4)]
    a=add(mul(z0,b0),mul(z1,b1))
    u=scale(mul(z0,b1),2); v=scale(mul(z1,b0),2)
    tau=add(mul(z0,b0),scale(mul(z1,b1),-1))
    assert add(mul(a,a),scale(mul(u,v),-1),scale(mul(tau,tau),-1)) == {}
    P=mul(add(a,u), add(mul(mul(a,a),v),scale(mul(add(scale(a,2),u),mul(tau,tau)),-1)))
    assert {sum(m) for m in P} == {8}
    power={(0,0,0,0):F(1)}; hopf=[]
    for m in range(1,7):
        power=mul(power,P); marked=power
        for s in range(m+2):
            if s: marked=mul(marked,u)
            for n in (2,3,4,6):
                actual=sphere(marked,n)
                cm=F(4**m*factorial(m)**2,factorial(2*m+1))
                angular=cm*comb(m-1,s-1) if 1 <= s <= m else F(0)
                k=4*m+s
                radial=F(factorial(k+1)*factorial(n-1),factorial(n+k-1))
                assert actual == angular*radial,(m,s,n,actual,angular*radial)
                hopf.append({'m':m,'s':s,'complex_sphere_dimension':n,'moment':str(actual)})
    z,w,t=[var(i,3) for i in range(3)]
    one={(0,0,0):F(1)}
    P3=mul(z,add(scale(z,2),scale(t,3),scale(one,5)))
    powers=[one]
    for m in range(1,8): powers.append(mul(powers[-1],P3))
    checks=0
    for degree in range(5):
        for i in range(degree+1):
            for j in range(degree-i+1):
                k=degree-i-j; Q={(i,j,k):F(1)}
                for m in range(degree+1,degree+3):
                    assert gaussian(mul(Q,powers[m])) == 0
                    checks += 1
    sharp=[]
    for d in range(1,7):
        # P=Z, Q=W^d attains the last possible mixed index in rank three.
        value=gaussian({(d,d,0):F(1)})
        assert value == factorial(d)
        sharp.append({'degree_Q':d,'last_nonzero_m':d,'moment':str(value)})
    chain=[]
    for h in range(1,5):
        # h circular pairs; one additional unused real coordinate is allowed.
        vs=[var(i,2*h) for i in range(2*h)]
        chainP=add(vs[2*(h-1)],*(mul(vs[2*i],vs[2*(i+1)+1]) for i in range(h-1)))
        chainpower={tuple([0]*(2*h)):F(1)}
        for m in range(1,3*h+2):
            chainpower=mul(chainpower,chainP)
            for d in range(4):
                Q={tuple(d if i == 1 else 0 for i in range(2*h)):F(1)}
                value=F(0)
                for mon,c in mul(Q,chainpower).items():
                    if all(mon[2*i] == mon[2*i+1] for i in range(h)):
                        term=c
                        for i in range(h): term *= factorial(mon[2*i])
                        value += term
                assert value == (factorial(m) if m == h*d else 0),(h,d,m,value)
                if m == h*d: chain.append({'circular_pairs':h,'degree_Q':d,'m':m,'moment':str(value)})
    # Rank-one quadratic prefix trap: l^2=-2, A=1, centered constant -1.
    # kappa_j=2^(j-1)(j-1)!+j! 2^(j-3) l^2 for j>=2.
    cumulants=[F(0),F(0)]+[F(2**(j-1)*factorial(j-1))-F(2)*F(2)**(j-3)*factorial(j) for j in range(2,6)]
    moments=[F(1)]
    for j in range(1,6):
        moments.append(sum((F(comb(j-1,k-1))*cumulants[k]*moments[j-k] for k in range(1,j+1)),F(0)))
    assert moments[1:4] == [0,0,-4]
    component_checks=0
    plus={(1,):F(1),(2,):F(1)}; minus={(-1,):F(1),(-3,):F(2)}
    pp={(0,):F(1)}; pm={(0,):F(1)}
    for m in range(1,12):
        pp=mul(pp,plus); pm=mul(pm,minus)
        assert pp.get((0,),0)+pm.get((0,),0) == 0
        mixed=F(1,2)*(mul(pp,{(-7,):F(1)}).get((0,),0)+mul(pm,{(8,):F(1)}).get((0,),0))
        if m>8: assert mixed == 0
        component_checks += 1
    receipt={'coverage':'Bounded independent rational checks; no full Lie-theoretic or all-rank sharpness replay.',
             'opposite_component_separator_checks':component_checks,
             'hopf_sphere_marker_checks':len(hopf),'hopf_rows':hopf,
             'ternary_quadratic_mixed_checks':checks,'cutoff_equality_witnesses':sharp,
             'quadratic_chain_equality_witnesses':chain,
             'rank_one_prefix_moments':[str(x) for x in moments[1:]],
             'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).with_name('verification.json')
    dest.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k != 'hopf_rows'},indent=2))


if __name__ == '__main__': main()
