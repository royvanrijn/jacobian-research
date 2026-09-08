"""Sage-free specialization and point transport for the compact six-family atlas."""
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT=Path(__file__).resolve().parents[2]
ATLAS=ROOT/'artifacts/generated-results/elliptic-curves/compact_six_r17_atlas_v1.json'


def polynomial(coefficients,t):
    value=F(0)
    for c in reversed(coefficients):value=value*t+F(c)
    return value


def rational(record,t):
    return polynomial(record['numerator_coefficients_low_to_high'],t)/polynomial(record['denominator_coefficients_low_to_high'],t)


def specialize(family,parameter):
    t=F(parameter);d=t.denominator
    A=polynomial(family['A_coefficients_low_to_high'],t)*d**8
    B=polynomial(family['B_coefficients_low_to_high'],t)*d**12
    model=(F(0),F(0),F(0),A,B)
    if 4*A**3+27*B**2==0:raise ArithmeticError('singular fibre')
    points=tuple((rational(r['X'],t)*d**4,rational(r['Y'],t)*d**6) for r in family['sections'])
    if len(points)!=17 or any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('specialized section failed')
    return model,points


def family_check(family,parameter,model,points):
    original,generic=specialize(family,parameter)
    if not cert.isomorphic(original,model):raise ArithmeticError('specialized model is not isomorphic')
    a,b=cert.weierstrass_invariants(original),cert.weierstrass_invariants(model)
    u=cert.square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
    if u is None:raise ArithmeticError('missing rational scale')
    for sign in (1,-1):
        v=sign*u
        if tuple((x/v**2,y/v**3) for x,y in generic)==tuple(points[:17]):return str(v)
    raise ArithmeticError('generic-section transport failed')
