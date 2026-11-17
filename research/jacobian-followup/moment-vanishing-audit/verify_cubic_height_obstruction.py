#!/usr/bin/env python3
"""Exact algebra and finite regressions for the formal cubic-height proof.

The universal rational identities are checked with indeterminates. Bounded
Wick and degree tests are separate; none asserts an exhaustive cubic census.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import sympy as sp
from verify_moment_transfers import add, mul, gaussian
from verify_affine_w_classification import smul, sadd, sc, unit, sinv, sexp, compose, polynomial


def system(A,B,C,Ap,Bp,Cp):
    p2=-2*B/(3*A); q2=-C/(3*A)
    p3=(4*B**2-3*A*C)/(9*A**2); q3=(2*B*C-3*A)/(9*A**2)
    p4=(-8*B**3+12*A*B*C-18*A**2)/(27*A**3)
    q4=(6*A*B-4*B**2*C+3*A*C**2)/(27*A**3)
    U=Ap*q3+Bp*q2; V=Ap*p3+Bp*p2+Cp
    U2=Ap*q4+Bp*q3+Cp*q2; V2=Ap*p4+Bp*p3+Cp*p2
    return U,V,U2,V2


def universal_identities():
    A,B,C,Ap,Bp,Cp,ff,a=sp.symbols('A B C Ap Bp Cp F a')
    U,V,U2,V2=system(A,B,C,Ap,Bp,Cp)
    Fp=U2+V2*ff-ff*(U+V*ff)
    X=(B**2-3*A*C)/a**4; ap=a*Ap/(3*A)
    Xp=(2*B*Bp-3*Ap*C-3*A*Cp)/a**4-4*ap*X/a
    f=-a*ff-B/a**2
    fp=-ap*ff-a*Fp-Bp/a**2+2*B*ap/a**3
    residual=sp.cancel((fp-Xp*(X-f**2)).subs(A,a**3/3))
    pullback=sp.cancel((V+a*Xp).subs(A,a**3/3))
    assert residual == pullback == 0
    # Verify the C-source recursion that produces p2,...,q4 directly.
    p=-2*B/(3*A);q=-C/(3*A)
    out=[]
    for n in (3,4):
        p,q=sp.cancel(sp.diff(p,C)+p*(-2*B/(3*A))+q),sp.cancel(sp.diff(q,C)+p*(-C/(3*A)))
        if n==3:
            assert sp.cancel(p-(4*B**2-3*A*C)/(9*A**2))==0
            assert sp.cancel(q-(2*B*C-3*A)/(9*A**2))==0
        else:
            assert sp.cancel(p-(-8*B**3+12*A*B*C-18*A**2)/(27*A**3))==0
            assert sp.cancel(q-(6*A*B-4*B**2*C+3*A*C**2)/(27*A**3))==0
        out.append(n)
    return {'Riccati_residual':str(residual),'V_plus_a_Xprime':str(pullback),'source_recursions_checked':out}


def cubic_integral(A,V,C,n):
    S=add({(i,3):x for i,x in enumerate(A) if x},
          {(i,2):x for i,x in enumerate(V) if x},
          {(i,1):x for i,x in enumerate(C) if x})
    out=unit(n); power={(0,0):F(1)}
    for k in range(1,n+1):
        power={m:c for m,c in mul(power,S).items() if m[0]<=n}
        for (i,t),c in power.items():
            if t%2==0:out[i]+=c*F(factorial(t),2**(t//2)*factorial(t//2)*factorial(k))
    return out


def tilted(h,s,v,u,r,n):
    t=[F(0),F(1)]+[F(0)]*(n-1);g=[F(0)]*(n+1)
    for _ in range(n+1):g=smul(t,compose(h,g,n),n)
    J=sadd(unit(n),sc(smul(t,compose([j*F(h[j]) for j in range(1,len(h))],g,n),n),-1))
    A=smul(t,compose(s,g,n),n);V=smul(t,compose(v,g,n),n);C=smul(t,compose(u,g,n),n)
    R=smul(t,compose(r,g,n),n)
    return smul(smul(sexp(R,n),sinv(J,n),n),cubic_integral(A,V,C,n),n)


def wick_rows():
    inputs=[([2,3,1],[1],[1,-2],[3,1,2],[-1,2,0,1]),
            ([0,2,1],[1],[1,2],[3,1,2],[-1,0,2,1]),
            ([0,0,2],[0,1],[0,3],[0,1,-2],[0,4,1,3]),
            ([1,2],[0,1],[2,-1],[1,3],[0,2,0,1])]
    rows=[]
    for h,s,v,u,r in inputs:
        P=add(polynomial(h,v,u,r),{(i,0,3):F(c) for i,c in enumerate(s) if c})
        predicted=tilted(h,s,v,u,r,8);power={(0,0,0):F(1)};values=[]
        for m in range(1,9):
            power=mul(power,P);actual=gaussian(power)
            assert actual==factorial(m)*predicted[m],(h,s,m,actual,predicted[m])
            values.append(str(actual))
        rows.append({'h':list(map(str,h)),'s':list(map(str,s)),'pure_moments_through_8':values})
    return rows



def wt_face_rows():
    # G is given in the (Z,T) monomial basis; h has genuine T coupling.
    inputs=[(F(2),F(1),F(3),{(3,0):F(1),(2,1):F(-2),(0,3):F(1),(1,0):F(2),(0,0):F(-1)}),
            (F(0),F(0),F(2),{}),
            (F(1),F(2),F(0),{(1,1):F(1),(3,0):F(-1,2)})]
    rows=[];n=8;t=[F(0),F(1)]+[F(0)]*(n-1)
    for aa,bb,cc,G in inputs:
        q=sadd(unit(n),sc(t,-bb));alpha=sc(smul(t,sinv(q,n),n),aa);beta=sc(smul(t,sinv(q,n),n),cc)
        g=add({(i,0):x for i,x in enumerate(alpha) if x},{(i,1):x for i,x in enumerate(beta) if x})
        gp={(0,0):F(1)};polys=[gp]
        for i in range(1,4):polys.append({m:x for m,x in mul(polys[-1],g).items() if m[0]<=n})
        exponent={}
        for (i,j),x in G.items():
            exponent=add(exponent,{(m+1,k+j):x*v for (m,k),v in polys[i].items() if m+1<=n})
        coeffs=[[F(0)]*(n+1) for _ in range(4)]
        for (i,k),x in exponent.items():coeffs[k][i]+=x
        R,C,V,A=coeffs
        predicted=smul(smul(sexp(R,n),sinv(q,n),n),cubic_integral(A,V,C,n),n)
        P=add({(0,1,0):aa,(1,1,0):bb,(0,1,1):cc},{(i,0,j):x for (i,j),x in G.items()})
        P={m:x for m,x in P.items() if x};power={(0,0,0):F(1)};values=[]
        for m in range(1,n+1):
            power=mul(power,P);v=gaussian(power);assert v==factorial(m)*predicted[m]
            values.append(str(v))
        rows.append({'a':str(aa),'b':str(bb),'c':str(cc),'pure_moments_through_8':values})
    return rows


def degree_rows():
    z=sp.Symbol('z');out=[]
    for D in (3,4):
        h=1+2*z+3*z**(D-1);s=2+z**(D-3) if D>3 else sp.Integer(2)
        v=4+z**(D-2);u=1+2*z+z**(D-1);r=1-z+z**D
        pp=z*s;qq=-h/2+z*v;cc=z*u;j=h-z*sp.diff(h,z)
        A=pp/h;B=qq/h;C=cc/h
        U,V,U2,V2=map(sp.cancel,system(A,B,C,sp.diff(A,z),sp.diff(B,z),sp.diff(C,z)))
        clear=h**6*pp**7
        q2=sp.cancel(clear*V)
        q1=sp.cancel(clear*(-U*V-sp.diff(V,z)-V*V2))
        q0=sp.cancel(clear*(-V*sp.diff(U,z)-V**2*U2+U*sp.diff(V,z)+U*V*V2))
        polys=[]
        for expr,bound in zip((q2,q1,q0),(13*D-20,13*D-19,13*D-18)):
            n,d=sp.fraction(expr);assert not d.has(z)
            deg=sp.degree(n,z);assert deg<=bound;polys.append(int(deg))
        H=j/h;R=z*r/h;L=sp.diff(H,z)/H-sp.diff(R,z)
        residual=sp.cancel((h**2*j)**2*(q2*(sp.diff(L,z)+L**2)+q1*L+q0))
        n,d=sp.fraction(residual);assert not d.has(z) and n!=0
        degree=int(sp.degree(n,z));assert degree<=19*D-24
        out.append({'D':D,'Q2_Q1_Q0_degrees':polys,'residual_numerator_degree':degree,'bound':19*D-24,'sufficient_prefix':19*D-22})
    return out


def main():
    rec={'universal_identities':universal_identities(),'Wick_vs_GL1_cubic_height':wick_rows(),
         'WT_face_Wick_vs_GL1':wt_face_rows(),'degree_specializations':degree_rows(),'sympy_version':sp.__version__,
         'scope':'Universal rational algebra plus stated finite regressions. The all-coefficient and no-algebraic-Riccati proofs are written, not inferred from samples.',
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'helper_hashes':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('verify_moment_transfers.py','verify_affine_w_classification.py')}}
    Path(__file__).with_name('cubic_height_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
