#!/usr/bin/env python3
"""Exact controls for one fixed A4 Keller map of a rational coefficient family.

The written proof covers all ranks and parameters. This replay tests a genuine
rational degree-five family with varying linear/cubic coefficients, including
the full identity, polynomial denominator clearing and entire quotient fibre.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path)
    args=parser.parse_args()
    assert sp.__version__=='1.14.0'
    w,T,S,V,Y,Z,Q=sp.symbols('w T S V Y Z Q')
    P=T**5+w/(w-1)*T**4+(w+1)*T**3+(w**2+1)*T**2+(w-2)*T+(w**3+1)/(w+2)
    c=[sp.expand(P).coeff(T,j) for j in range(6)]
    r=sp.cancel(c[1]/c[3])
    L=(w-1)*(w+1)**5/sp.Integer(2048)
    h={j:sp.cancel(c[j]/c[1]*r**j*L**(j+1)) for j in (4,5)}
    for value in h.values():
        assert sp.Poly(value,w,domain=sp.QQ).degree()>0
        assert sp.fraction(value)[1].is_number
    G=S+V*S**3+sum(h[j]*V**j*S**j for j in h)
    E=G+Y*S**2-Z/2
    gamma={V:c[3]/(c[1]*L**2),Y:c[2]/(c[1]*L),Z:-2*L*c[0]/c[1]}
    assert sp.cancel(E.subs(gamma)-L/c[1]*P.subs(T,S/L))==0
    for j in (4,5):
        assert sp.cancel(h[j]*gamma[V]**j-c[j]/c[1]*L**(1-j))==0
    assert sp.cancel(E.subs(gamma)-L/c[1]*P.subs(T,S/L)+S**3)!=0
    beta=sp.cancel((sp.diff(G,S)-1-V*S**2)/S)
    B=Q+beta
    C=2*G-B*S**2
    D=1-S*Q+V*S**2
    jac=sp.det(sp.Matrix([[sp.diff(f,z) for z in (S,Q)] for f in (B,C)]))
    assert sp.expand(jac+2*D)==0
    assert sp.expand(D.subs(Q,-2*Y-beta)-sp.diff(E,S))==0
    # Preserved w makes all coefficient derivatives disappear from the minor.
    jac4=sp.det(sp.Matrix([[sp.diff(f,z) for z in (w,V,S,Q)] for f in (w,V,-B/2,C)]))
    assert sp.expand(jac4-D)==0
    x,y,q,z=sp.symbols('x y q z')
    t=1+x*y
    chart=(t*q,x/t,y+x*q)
    chart_jac=sp.det(sp.Matrix([[sp.diff(f,v) for v in (x,y,q)] for f in chart]))
    assert sp.cancel(chart_jac-1/t)==0
    assert sp.cancel(D.subs(dict(zip((V,S,Q),chart)))-1/t)==0
    assert sp.cancel(chart_jac*t**2/t-1)==0

    # Complete fibre at w=3: L=1, v=4 and all Taylor coefficients are nonzero.
    w0=sp.Integer(3)
    assert L.subs(w,w0)==1
    v0=sp.cancel(gamma[V].subs(w,w0))
    Y0=sp.cancel(gamma[Y].subs(w,w0))
    Z0=sp.cancel(gamma[Z].subs(w,w0))
    assert (v0,Y0,Z0)==(4,10,-sp.Rational(56,5))
    h0={j:value.subs(w,w0) for j,value in h.items()}
    modulus=sp.Poly(E.subs({w:w0,V:v0,Y:Y0,Z:Z0}),S,domain=sp.QQ)
    assert modulus.degree()==5
    assert sp.expand(modulus.as_expr()-P.subs({w:w0,T:S}))==0
    assert sp.gcd(modulus,modulus.diff()).degree()==0
    def red(f):return sp.rem(sp.Poly(f,S,domain=sp.QQ),modulus).as_expr()
    def mul(*factors):
        out=sp.Integer(1)
        for factor in factors:out=red(out*factor)
        return out
    def power(f,n):
        out=sp.Integer(1)
        for _ in range(n):out=mul(out,f)
        return out
    d=red(modulus.diff().as_expr())
    dinv=sp.invert(d,modulus.as_expr(),S)
    assert mul(d,dinv)==1
    root_Q=red(-2*Y0-beta.subs({w:w0,V:v0}))
    root_x=mul(S,dinv)
    root_y=red(root_Q-v0*S)
    root_z=mul(d,d,red(v0*d-mul(root_y,root_y,1+3*dinv)))
    root_t=red(1+mul(root_x,root_y))
    root_q=red(mul(root_t,root_t,root_z)+mul(root_y,root_y,1+3*root_t))
    assert red(root_t-dinv)==0
    assert red(root_q-v0*d)==0
    assert mul(root_t,root_q)==v0
    Bvalue=red(root_y+3*mul(root_x,root_q))
    Cvalue=red(mul(root_x,5-3*root_t)-mul(power(root_x,3),root_z))
    for j,value in h0.items():
        Bvalue=red(Bvalue+j*value*mul(power(root_t,2),power(root_x,j-2),power(root_q,j)))
        Cvalue=red(Cvalue-(j-2)*value*power(mul(root_x,root_q),j))
    assert red(-Bvalue/2-Y0)==0
    assert red(Cvalue-Z0)==0
    assert mul(root_x,d)==S
    assert red(Cvalue-Z0-1)!=0

    root=Path(__file__).resolve().parents[3]
    inputs=[
        'research/jacobian/universal-keller-structure/PROOF.md',
        'unpublished/common-arithmetic-fibers/sections/02-quadratic-gauge.tex',
        'unpublished/common-arithmetic-fibers/sections/03-realization.tex',
        'research/arithmetic/ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md',
    ]
    report=dict(schema_version=1,kind='internal_exact_replay',result='PASS',
        command='.venv/bin/python research/jacobian/monodromy-theory/verify_fixed_keller_finite_etale_family.py',
        python=platform.python_version(),sympy=sp.__version__,
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs={path:hashlib.sha256((root/path).read_bytes()).hexdigest() for path in inputs},
        family_polynomial=str(P),denominator_clearer=str(L),
        polynomial_profile_degrees={str(j):sp.Poly(value,w).degree() for j,value in h.items()},
        checks=dict(full_parameter_identity=True,polynomial_profiles=True,
            four_dimensional_compact_jacobian='D',full_polynomial_jacobian=1,
            complete_quotient_reconstruction=True,quotient_degree=5,
            target_at_w_three=['3',str(v0),str(Y0),str(Z0)],invalid_controls=2),
        scope='Exact rational degree-five family and complete quotient controls. The written proof supplies uniform degree, arbitrary finite-etale families and the infinite rank301 Hasse application. Dimension four; arbitrary fixed A3 product-family realization remains open.',
        independent_replay=False,formal_verification=False,external_review=False)
    output=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.write_receipt:args.write_receipt.write_text(output)
    print(output,end='')


if __name__=='__main__':main()
