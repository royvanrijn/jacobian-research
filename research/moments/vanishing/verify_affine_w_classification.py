#!/usr/bin/env python3
"""Exact bounded regressions for the written affine-W classification.

Uses only rational arithmetic. No CAS output or finite prefix is promoted
into an all-order theorem. The adjacent proof supplies every universal step.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
from verify_moment_transfers import add, mul, scale, gaussian


def smul(a,b,n):
    out=[F(0)]*(n+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:n+1-i]): out[i+j]+=x*y
    return out


def sadd(*args): return [sum(x,F(0)) for x in zip(*args)]


def sc(a,c): return [x*c for x in a]


def unit(n): return [F(1)]+[F(0)]*n


def sinv(a,n):
    assert a[0]
    out=[1/a[0]]+[F(0)]*n
    for i in range(1,n+1): out[i]=-sum((a[j]*out[i-j] for j in range(1,i+1)),F(0))/a[0]
    return out


def sexp(a,n):
    assert a[0] == 0
    out=unit(n)
    for i in range(1,n+1): out[i]=sum((j*a[j]*out[i-j] for j in range(1,i+1)),F(0))/i
    return out


def sroot_inverse(a,n):
    assert a[0] == 1
    b=a[:]; b[0]=0; out=unit(n); power=unit(n); coefficient=F(1)
    for j in range(1,n+1):
        power=smul(power,b,n); coefficient*=F(-1,2)-(j-1); coefficient/=j
        out=sadd(out,sc(power,coefficient))
    return out


def compose(poly,g,n):
    out=[F(0)]*(n+1); power=unit(n)
    for c in poly:
        out=sadd(out,sc(power,F(c))); power=smul(power,g,n)
    return out


def tilt(h,v,u,r,n):
    t=[F(0),F(1)]+[F(0)]*(n-1); g=[F(0)]*(n+1)
    for _ in range(n+1): g=smul(t,compose(h,g,n),n)
    derivative=[j*F(h[j]) for j in range(1,len(h))]
    J=sadd(unit(n),sc(smul(t,compose(derivative,g,n),n),-1))
    K=sadd(unit(n),sc(smul(t,compose(v,g,n),n),-2))
    S=sadd(smul(t,compose(r,g,n),n),sc(smul(smul(smul(t,t,n),smul(compose(u,g,n),compose(u,g,n),n),n),sinv(K,n),n),F(1,2)))
    return smul(smul(sexp(S,n),sinv(J,n),n),sroot_inverse(K,n),n)


def polynomial(h,v,u,r):
    out={}
    for coeffs,w,t in ((h,1,0),(v,0,2),(u,0,1),(r,0,0)):
        out=add(out,{(i,w,t):F(c) for i,c in enumerate(coeffs) if c})
    return out


def ordinary_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=F(x)*F(y)
    return out


def main():
    # Compare (4) to direct independent Wick expansion, including both h charts
    # and arbitrary polynomials that do not satisfy the vanishing premise.
    tests=[([2,3,1],[1,-2],[3,1,2],[-1,2,0,1]),
           ([0,2,1],[1,2],[3,1,2],[-1,0,2,1]),
           ([0,0,2],[0,3],[0,1,-2],[0,4,1,3]),
           ([1,0,0,2],[1,-1],[2,1],[0,3,0,1]),
           ([2],[0],[1,3],[0,F(-1,4),F(-3,2),F(-9,4)])]
    rows=[]
    for h,v,u,r in tests:
        P=polynomial(h,v,u,r); power={(0,0,0):F(1)}; generated=tilt(h,v,u,r,8); values=[]
        for m in range(1,9):
            power=mul(power,P); value=gaussian(power)
            assert value == factorial(m)*generated[m],(h,m,value,generated[m])
            values.append(str(value))
        rows.append({'h':list(map(str,h)),'moments_through_8':values})
    # Infinite linear-branch k freedom is tested with nonconstant repairs.
    failures=[]
    for a,b,k in ((F(1),F(2),[F(0)]),(F(2),F(-1),[F(1),F(2)])):
        h=[a,b]; h2=ordinary_mul(h,h); h3=ordinary_mul(h2,h)
        v=[-b,-F(3)*b*b/(2*a),-b*b*b/(2*a*a)]
        u=ordinary_mul(h2,k); r=[F(0)]+[-a*a*x/2 for x in ordinary_mul(h,ordinary_mul(k,k))]
        P=polynomial(h,v,u,r); power={(0,0,0):F(1)}; mixed=[]
        for m in range(1,7):
            power=mul(power,P); assert gaussian(power) == 0
            value=gaussian(mul({(1,0,0):F(1)},power))
            assert value == factorial(m)*a*b**(m-1)
            mixed.append(str(value))
        failures.append({'a':str(a),'b':str(b),'repair_degree':len(k)-1,'Z_marked_moments':mixed})
    # Sharp equality in the safe branch, including one quintic beyond cubics.
    sharp=[]
    for a,u,j in ((F(1),[0,1],1),(F(2),[1,3],1),(F(1),[0,0,1],2)):
        r=[F(0)]+[-x/(2*a) for x in ordinary_mul(u,u)]
        P=polynomial([a],[0],u,r); powers=[{(0,0,0):F(1)}]
        maxd=3 if j==1 else 2
        for m in range(1,(2*j+1)*maxd+2):powers.append(mul(powers[-1],P))
        for d in range(1,maxd+1):
            last=(2*j+1)*d; Q={(0,d,0):F(1)}
            value=gaussian(mul(Q,powers[last])); expected=F(factorial(last))*(-u[-1]**2*a**(2*j-1)/2)**d
            assert value == expected and value != 0
            assert gaussian(mul(Q,powers[last+1])) == 0
            sharp.append({'a':str(a),'deg_u':j,'deg_Q':d,'last_nonzero_m':last,'moment':str(value)})
    # A rational specialization checks C^2 B=A^3 independently.
    # The universal coefficient identity is established in the written proof.
    a,b=F(2),F(-3); h=[a,b]; h3=ordinary_mul(ordinary_mul(h,h),h)
    v=[-b,-3*b*b/(2*a),-b**3/(2*a*a)]
    B=h+[F(0)]*(len(v)+1-len(h))
    for i,c in enumerate(v):B[i+1]-=2*c
    assert [a*a*x for x in B] == h3
    receipt={'scope':'Exact bounded checks for all-degree written classification, not an exhaustive cubic search.',
             'formula_vs_direct_Wick_rows':rows,'repaired_failure_rows':failures,
             'sharp_safe_rows':sharp,'linear_branch_polynomial_identity_at_a2_bminus3':'PASS',
             'finite_gate':'Written numerator degree <=5D-4 gives 5D-3 sufficient pure moments for total degree D>=2. Cubics:12; quartics:17. The h(0)=0 numerator has degree <=3. No optimality claim.',
             'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'arithmetic_helper_sha256':hashlib.sha256(Path(__file__).with_name('verify_moment_transfers.py').read_bytes()).hexdigest()}
    Path(__file__).with_name('affine_w_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
