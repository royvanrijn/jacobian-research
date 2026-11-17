#!/usr/bin/env python3
"""Bounded exact original-moment controls for the algebraic angular theorem."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as s

if not __debug__:
    raise RuntimeError('Run without -O: exact assertions must remain enabled')

z,w,y,t=s.symbols('z w y t')
ORDER=5

def wick(expr):
    return s.expand(sum(co*s.factorial(i)*s.factorial2(k-1)
                        for (i,j,k),co in s.Poly(s.expand(expr),z,w,y).terms()
                        if i==j and k%2==0))

def truncate(expr):
    return s.expand(sum(co*t**int(ex[0])
                        for ex,co in s.Poly(s.expand(expr),t).terms()
                        if int(ex[0])<=ORDER))

def periods(phase):
    powers=[s.Integer(1)]
    for _ in range(ORDER):
        powers.append(s.expand(powers[-1]*phase))
    def evaluate(marker):
        return truncate(sum(wick(marker*power)*t**m/s.factorial(m)
                            for m,power in enumerate(powers)))
    return evaluate

def project_nonpositive(expr):
    return s.expand(sum(co*z**i*w**j*y**k*t**e
                        for (i,j,k,e),co in s.Poly(s.expand(expr),z,w,y,t).terms()
                        if i<=j))

def min_weight(expr):
    return min(i-j for i,j,k,e in s.Poly(s.expand(expr),z,w,y,t).monoms())

def transport(expr,R,sign=1):
    expr=project_nonpositive(expr)
    if expr==0:
        return s.Integer(0)
    gamma=min(i-j for i,j,k in s.Poly(R,z,w,y).monoms())
    count=max(0,-min_weight(expr)//gamma)
    return project_nonpositive(sum(sign**n*t**n*expr*R**n/s.factorial(n)
                                   for n in range(count+1)))

core=(z*w)**2+z*w*y+y**2
remainder=z+z**2*w*y**2+z**2*y
full=core+remainder
phi0,phi=periods(core),periods(full)
markers=[z**i*w**j*y**k
         for i in range(4) for j in range(4-i) for k in range(4-i-j)]
for Q in markers:
    weight=min_weight(Q)
    N=max(0,-weight)
    assert truncate(phi(Q)-sum(t**n*phi0(Q*remainder**n)/s.factorial(n)
                               for n in range(N+1)))==0
    assert truncate(phi0(Q)-sum((-t)**n*phi(Q*remainder**n)/s.factorial(n)
                                for n in range(N+1)))==0
    Q0=sum(co*z**i*w**j*y**k
           for (i,j,k),co in s.Poly(Q,z,w,y).terms() if i==j)
    assert phi0(Q)==phi(Q0)

inputs=[1-2*w,w**3+y*t*w**2+z,z*w**2+t**2*y*w**3]
for Q in inputs:
    assert transport(transport(Q,remainder,1),remainder,-1)==project_nonpositive(Q)
    assert transport(transport(Q,remainder,-1),remainder,1)==project_nonpositive(Q)
    TQ=transport(Q,remainder)
    lhs=project_nonpositive(s.diff(TQ,t)+core*TQ)
    rhs=transport(s.diff(Q,t)+full*Q,remainder)
    assert s.expand(lhs-rhs)==0

# Deliberately invalid ideal quotient: discarding Z would incorrectly kill ZW.
assert wick(z)==0 and wick(z*w)==1

phiS=periods((z*w)**2)
phiSZ=periods((z*w)**2+z)
f=phiS(1)
assert truncate(f-1-2*t*phiS(z*w))==0
assert truncate(phiSZ(w)-t*phiS(z*w))==0
assert phiSZ(1-2*w)==1
for j in range(6):
    for n in range(ORDER+1):
        assert phiS((z*w)**j).coeff(t,n)==s.Rational(factorial(2*n+j),factorial(n))

# Safe high-degree remainder: marker time degree depends on angular weight,
# not phase degree.  Here the minimum positive weight is two.
safeR=z**2+z**3*w*y**3+z**4*w**2*y
safe=periods(safeR)
safe_checks=0
for Q in markers:
    series=safe(Q)
    D=s.Poly(Q,z,w,y).total_degree()
    assert series==0 or s.degree(series,t)<=D//2
    safe_checks+=1
assert safe(w**4)==40692*t**2

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,
                    default=Path(__file__).with_name('gaussian_angular_core_full_image_verification.json'))
args=parser.parse_args()
here=Path(__file__).resolve().parent
pins=['GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md',Path(__file__).name,
      'GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md',
      'MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md']
result=dict(status='passed',scope='bounded exact controls; all-degree algebraic proof is written separately',
            environment=dict(python=platform.python_version(),sympy=s.__version__),
            original_moment_order=[0,ORDER],original_marker_count=len(markers),
            finite_transport_identities=2*len(markers),
            exact_amplitude_inverse_and_connection_cases=len(inputs),
            strict_literal_image_example=True,positive_subspace_not_ideal=True,
            safe_minimum_weight=2,safe_cutoff_marker_count=safe_checks,
            sha256={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in pins},
            formal_verification=False,independent_replay=False,external_review=False)
args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))
