#!/usr/bin/env python3
"""Exact identities for QCC1. All-degree invariant and stability arguments are written."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
x,q,p,z,h=sp.symbols("x q p z h")
coordinates=(x,q,p,z)
a=1-3*x*q
c=(1+3*x*q)/2
R=x*(2-3*x*q)
D0=c*p-3*q**2*z
B=3*x**2*p+2*a*z

def pb(f,g):
    return sp.expand(sp.diff(f,p)*sp.diff(g,x)-sp.diff(f,x)*sp.diff(g,p)+
                     sp.diff(f,z)*sp.diff(g,q)-sp.diff(f,q)*sp.diff(g,z))

def normal(f,g):
    out=0
    if f==0 or g==0:return sp.Integer(0)
    for i in range(sp.Poly(f,p).degree()+1):
        for j in range(sp.Poly(f,z).degree()+1):
            out+=h**(i+j)*sp.diff(f,p,i,z,j)*sp.diff(g,x,i,q,j)/(factorial(i)*factorial(j))
    return sp.expand(out)

def zero(f):
    assert sp.cancel(f)==0,sp.cancel(f)

def V(f):
    return sp.expand(c*sp.diff(f,x)-3*q**2*sp.diff(f,q))

def main():
    zero(pb(D0,R)-1)
    for generator in (x,q,B):
        zero(pb(R,generator))
    zero(2*a*c+9*x**2*q**2-1)
    zero(2*a*D0+3*q**2*B-p)
    zero(-3*x**2*D0+c*B-z)
    zero(pb(D0,a)+sp.Rational(3,2)*q*a)
    zero(pb(D0,q)+3*q**2)
    zero(pb(D0,B)-sp.Rational(9,2)*q*B)

    U,Vmid,W=a**3*B,a*q*B,q**3*B**2
    for f in (U,Vmid,W):
        zero(pb(R,f));zero(pb(D0,f))
    zero(U*W-Vmid**3)
    zero(pb(U,Vmid)-2*U)
    zero(pb(Vmid,W)-2*W)
    zero(pb(U,W)-6*Vmid**2)

    # Arbitrary monomial controls for the written semigroup factorization.
    factorizations=0
    aa,bb,qq=sp.symbols("aa bb qq")
    for k in range(8):
        for j in range((3*k)//2+1):
            i=3*k-2*j
            lhs=aa**i*qq**j*bb**k
            if j<=k:
                rhs=(aa**3*bb)**(k-j)*(aa*qq*bb)**j
            else:
                rhs=(aa*qq*bb)**i*(qq**3*bb**2)**(j-k)
            zero(lhs-rhs);factorizations+=1

    # A noncanonical polynomial fiber shift, distinct from a canonical R shift.
    H=B
    image_p=p-2*a*H
    image_z=z+3*x**2*H
    substitution={x:x,q:q,p:image_p,z:image_z}
    zero(D0.subs(substitution,simultaneous=True)-(D0-H))
    zero(B.subs(substitution,simultaneous=True)-B)
    zero(pb(image_p,x)-1+6*a*x**2)
    assert pb(image_p,x)!=1
    Hgood=R**2+R**4
    gp,gz=p-2*a*Hgood,z+3*x**2*Hgood
    zero(pb(gp,x)-1);zero(pb(gz,q)-1);zero(pb(gp,gz))

    b=sp.symbols("b")
    beta=b-9*q**2
    y=q-x*beta/3
    u=x*y
    Hlong=(y**4*(18*u**2+78*u+125)/20+
           3*beta*y**2*(u**3+5*u**2+10*u-5)/10-
           beta**2*(9*u+2)/6-x**2*beta**3/6)
    zero(Hlong.subs(x,0)-(-b**2/3+sp.Rational(9,2)*q**2*b-sp.Rational(29,4)*q**4))
    assert sp.diff(Hlong.subs(x,0),b)!=0

    zero(normal(D0,R)-normal(R,D0)-h)
    tested=[]
    for m in (0,1,2,5):
        Dm=D0+sp.Rational(3*m,2)*h*q
        zero(normal(Dm,R)-normal(R,Dm)-h)
        basis=[a**(m-2*j)*q**j for j in range(m//2+1)]
        for f in basis:
            zero(V(f)+sp.Rational(3*m,2)*q*f)
        tested.append({"m":m,"dimension":len(basis)})
    divergence=sp.diff(c,x)+sp.diff(-3*q**2,q)
    zero(divergence+sp.Rational(9,2)*q)
    zero(D0+h*divergence/2-(D0-sp.Rational(9,4)*h*q))

    receipt={
        "claim":"QCC1","status":"PASS",
        "classical_core_and_coordinate_inverse":True,
        "toric_centralizer_generator_identities":True,
        "toric_Poisson_character":True,
        "semigroup_factorization_controls":factorizations,
        "noncanonical_and_canonical_shear_controls":True,
        "Long_correction_B_dependence":True,
        "normal_order_exact_Weyl_core":True,
        "polynomial_eigenbasis_controls":tested,
        "symmetric_ordering_shift":str(h*divergence/2),
        "coverage_boundary":"Exact symbolic identities and listed finite controls only. Full invariant-ring and eigenbasis exhaustion, all-degree shear classification, the finite-eigenspace centralizer criterion and stable noncompletion use written proofs. Quantum obstruction applies to the specified normal-ordered resonant family, not every ordering or Long's corrected D.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("CUBIC_CORE_NONCOMPLETION.md","check_cubic_core.py")}
    }
    (HERE/"cubic_core_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
