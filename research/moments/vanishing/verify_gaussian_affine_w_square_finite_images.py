#!/usr/bin/env python3
"""Bounded exact controls for affine W-square finite actual images."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from math import factorial
from pathlib import Path
import sympy as s

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

z,w,y,t,x=s.symbols('z w y t x')
la,mu,c,a,h,u,k=s.symbols('lambda mu c a h u k', nonzero=True)
ell=la*z+mu
d=1-2*a*t
q=c*c+2*d*ell/t
A=t*ell+t*t*c*c/(2*d)
B=t*h+t*t*c*u/d
phi=t*k+t*t*u*u/(2*d)-(z-B)**2/(4*A)
N=t*k*q-d*(z-t*h)**2/(2*t*t)+c*u*(z-t*h)+t*ell*u*u
assert s.cancel(phi-N/q)==0
Ncoeff=(-z*z/(2*t*t)+(a*z*z+z*h)/t
        +2*ell*k-2*a*z*h-h*h/2+c*u*z
        +t*((c*c-4*a*ell)*k+a*h*h-c*u*h+ell*u*u))
assert s.expand(N-Ncoeff)==0

# Independently expand U after the exact moving-root coordinate.  The
# samples only control the general derivative/divisibility proof in text.
divisibility=[]
for cs in (0,3):
    hs=2+z+3*z**3
    us=1-2*z+z**2
    ks=3+2*z+z**4
    ells=2*z+1
    Ns=s.expand(t*ks*(cs**2+2*ells/t)-(z-t*hs)**2/(2*t*t)
                +cs*us*(z-t*hs)+t*ells*us**2)
    beta=-s.Rational(1,2)-t*cs**2/4
    scale=t/4
    expr=s.Poly(s.expand(Ns.subs(z,beta+scale*x*x)),x)
    K=s.degree(Ns,z)-1
    for j in range(2,K+1):
        coeff=expr.coeff_monomial(x**(2*j+2))
        if coeff==0:
            continue
        deriv=s.expand(scale**(j+1)*s.diff(Ns,z,j+1).subs(z,beta)/factorial(j+1))
        assert s.expand(coeff-deriv)==0
        quotient=s.cancel(coeff/t**j)
        assert quotient.is_polynomial(t)
        assert s.degree(coeff,t)>=j
        divisibility.append({'c':cs,'j':j,'time_degree':int(s.degree(coeff,t))})

ORDER=4
def trunc(f):
    f=s.expand(f)
    return s.expand(sum(co*t**int(exp[0]) for exp,co in s.Poly(f,t).terms()
                        if int(exp[0])<=ORDER))
def mul(f,g):
    return trunc(f*g)
def ep(f):
    out=s.Integer(1)
    power=s.Integer(1)
    for j in range(1,ORDER+1):
        power=mul(power,f)
        out+=power/s.factorial(j)
    return trunc(out)
def wick(f):
    out=s.Integer(0)
    for (i,j,b),co in s.Poly(s.expand(f),z,w,y).terms():
        if i==j and b%2==0:
            out+=co*s.factorial(i)*s.factorial2(b-1)
    return out

def vacuum_control(lv,mv,cv,h0,h1,u0,k0,k1,av):
    P=(lv*z+mv)*w*w+(h0+h1*z+cv*y)*w+av*y*y+u0*y+k0+k1*z
    moment=s.Integer(1)
    actual=s.Integer(1)
    for n in range(1,ORDER+1):
        moment=s.expand(moment*P)
        actual+=wick(moment)*t**n/s.factorial(n)
    r=s.Integer(0)
    for _ in range(ORDER+1):
        r=trunc(t*(lv*r*r+h1*r+k1))
    J=trunc(1-t*(2*lv*r+h1))
    Jinv=s.Integer(0)
    jdiff=1-J
    pwr=s.Integer(1)
    for _ in range(ORDER+1):
        Jinv+=pwr
        pwr=mul(pwr,jdiff)
    dinv=sum((2*av*t)**n for n in range(ORDER+1))
    dsqrt=sum(s.binomial(2*n,n)*av**n*t**n/2**n for n in range(ORDER+1))
    R=trunc(t*(mv*r*r+h0*r+k0)+mul(t*t*(cv*r+u0)**2,dinv)/2)
    predicted=mul(mul(Jinv,dsqrt),ep(R))
    assert s.expand(actual-predicted)==0
    return {'lambda':lv,'mu':mv,'c':cv,'a':av,
            'h0':h0,'h1':h1,'u0':u0,'k0':k0,'k1':k1,
            'actual_vacuum_jet':str(s.expand(actual))}

vacua=[vacuum_control(1,2,3,4,5,6,7,8,9),
        vacuum_control(1,2,3,4,0,5,6,0,0),
        vacuum_control(2,0,0,3,1,2,4,1,0)]

# Exact low-degree discriminant relation follows from the selected root.
r,h1,k1=s.symbols('r h1 k1')
root_relation=r-t*(la*r*r+h1*r+k1)
J=1-t*(2*la*r+h1)
assert s.rem(s.expand(J*J-((1-h1*t)**2-4*la*k1*t*t)),root_relation,r)==0

result={'status':'PASS','scope':'bounded algebraic controls only; unrestricted proof is written',
        'generic_schur_identities':2,'a0_high_coefficient_controls':divisibility,
        'original_wick_vacuum_controls':vacua,'moment_orders':[0,ORDER],
        'low_degree_discriminant_identity':True,
        'formal_verification':False,'independent_replay':False,'external_review':False}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,
                        help="Write the exact JSON receipt to this path")
    args=parser.parse_args()
    here=Path(__file__).resolve()
    dependencies=[
        "GAUSSIAN_AFFINE_W_SQUARE_FINITE_IMAGES.md",
        "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
        "GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md",
        "GAUSSIAN_MOVING_COMMON_ROOT_SCALAR_EXCLUSION.md",
        "GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md",
        "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md",
        "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
    ]
    paths=[here]+[here.with_name(name) for name in dependencies]
    receipt={
        "status":"passed", "checked_on":"2026-11-17",
        "environment":{"python":platform.python_version(),"sympy":s.__version__},
        "sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "checks":result,
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "unrestricted_result_from_finite_prefix":False,
        "scope":"Bounded exact algebraic and original Wick controls; the unrestricted theorem is the pinned written proof with retained comparison imports.",
    }
    encoded=json.dumps(receipt,indent=2)+"\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding="utf-8")
    print(encoded,end="")


if __name__=="__main__":
    main()
