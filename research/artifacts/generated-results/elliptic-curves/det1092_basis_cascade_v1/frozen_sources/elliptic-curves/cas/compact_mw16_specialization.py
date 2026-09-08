"""Exact finite-parameter specialization of the portable compact MW16 atlas.

Coordinates use t=n/d and x=d^4 X(t), y=d^6 Y(t). Membership is checked;
specialized independence must still be certified on each candidate.
"""
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json'


def polynomial(coefficients, parameter):
    result=F(0)
    for value in reversed(coefficients):
        result=result*parameter+F(value)
    return result


def rational(record, parameter):
    numerator=polynomial(record['numerator_coefficients_low_to_high'],parameter)
    denominator=polynomial(record['denominator_coefficients_low_to_high'],parameter)
    if not denominator:
        raise ArithmeticError('retained section has a pole at this parameter')
    return numerator/denominator


def specialize(family, parameter):
    t=F(parameter);d=t.denominator
    A=polynomial(family['A_coefficients_low_to_high'],t)*d**8
    B=polynomial(family['B_coefficients_low_to_high'],t)*d**12
    if not 4*A**3+27*B**2:
        raise ArithmeticError('singular fibre')
    model=(F(0),F(0),F(0),A,B)
    sections=family['sections']
    if len(sections)!=16 or [p['basis_index'] for p in sections]!=list(range(16)):
        raise ArithmeticError('incomplete generic MW16 section roster')
    points=tuple((rational(s['X'],t)*d**4,rational(s['Y'],t)*d**6) for s in sections)
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):
        raise ArithmeticError('specialized section equation failed')
    return model,points


def native_parameter(family, parameter):
    t=F(parameter);a,b,c,d=map(F,family['base_matrix_a_b_c_d'])
    if a*d==b*c:
        raise ArithmeticError('singular base map')
    return None if c*t+d==0 else (a*t+b)/(c*t+d)
