#!/usr/bin/env python3
"""Complete SageMath verifier: two trivial-torsion generic rank >= 18 families.

Usage: sage -python verify.py (or Miniforge sage environment Python)
This script is intentionally stand-alone: it imports neither Rank Hunter nor
Rank Ledger's original verify.py, and does not use the genus-one rank-19 pair.
"""
from pathlib import Path
import hashlib
import json

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix

ROOT = Path(__file__).resolve().parent
HASHES = {
    'parent-rank17.json': '18d5ef159dd44621a520a5095ffc4ea27e5ea2cc0a97d7e0bb7e7dc100ea7db3',
    'cover-4339.json': 'dcdd3403e0ed5b07c1a5cc00ffd068e7de0b272bf96d198cc73cea7eed134dc8',
    'cover-36773.json': '7b4246abbaf6bb846a0445518272f924ffc71f85dfb9aad7271feb21709baff5',
}


def load(name):
    file = ROOT / 'data' / name
    if not file.is_file():
        raise RuntimeError('Missing '+str(file)+'. Restore the checksum-pinned proof data.')
    data = file.read_bytes()
    assert hashlib.sha256(data).hexdigest() == HASHES[name], name+' hash mismatch'
    return json.loads(data)


def build(data):
    ring = PolynomialRing(QQ, 'U')
    field = ring.fraction_field()
    a4 = ring(data['jacobian_model']['a4_low_to_high'])
    a6 = ring(data['jacobian_model']['a6_low_to_high'])
    curve = EllipticCurve(field, [0,0,0,a4,a6])
    sections = []
    for row in data['neighbor_sections']:
        x = field(ring(row['x_numerator'].replace('^','**'))) / ring(row['x_denominator'].replace('^','**'))
        y = field(ring(row['y_numerator'].replace('^','**'))) / ring(row['y_denominator'].replace('^','**'))
        sections.append(curve(x,y))  # exact defining identity, or raises
    return ring, field, curve, sections


def shioda_height(p, chi):
    if not p:
        return QQ(0)
    x = p[0]
    return QQ(max(x.numerator().degree(), x.denominator().degree()+2*chi))


def check_model(label, data, curve, sections, ring, chi, determinant):
    delta = ring(curve.discriminant())
    assert ring(curve.a4()).degree() <= 4*chi, label+' a4 degree'
    assert ring(curve.a6()).degree() <= 6*chi, label+' a6 degree'
    assert delta.degree() == 12*chi, label+' delta degree'
    assert delta.gcd(delta.derivative()).degree() == 0, label+' non-I1 singular fibers'
    # With all singular fibers I1, no reducible-component terms contribute.
    # Our minimal Weierstrass models admit this pole-degree Shioda height formula.
    hs = [shioda_height(p, chi) for p in sections]
    gram = matrix(QQ, len(sections))
    for i, p in enumerate(sections):
        gram[i,i] = hs[i]
        for j in range(i):
            gram[i,j] = gram[j,i] = (shioda_height(p+sections[j],chi)-hs[i]-hs[j])/2
    expected = matrix(QQ, data['exact_generic_section_basis']['height_matrix'])
    assert gram == expected, label+' Gram matrix changed'
    assert gram.is_positive_definite(), label+' Gram not positive definite'
    assert gram.det() == determinant, label+' unexpected determinant'
    assert len(sections) == len(expected.rows())
    print(f'{label}: {len(sections)} exact rational sections, squarefree discriminant, '
          f'positive Gram, determinant {determinant}: PASS')
    return gram


def check_cover_identity(parent_data, cover_data, ring):
    p4 = ring(parent_data['jacobian_model']['a4_low_to_high'])
    p6 = ring(parent_data['jacobian_model']['a6_low_to_high'])
    c4 = ring(cover_data['jacobian_model']['a4_low_to_high'])
    c6 = ring(cover_data['jacobian_model']['a6_low_to_high'])
    base = cover_data['quadratic_basechange']
    n = ring(base['parameter_numerator'])
    d = ring(base['parameter_denominator'])
    t = n/d
    assert d**8*p4(t) == c4, 'a4 pullback mismatch'
    assert d**12*p6(t) == c6, 'a6 pullback mismatch'
    q = ring(base['branch_coefficients'])
    w = ring(base['branch_root_numerator'])
    assert q[2]*n*n + q[1]*n*d + q[0]*d*d == w*w, 'quadratic cover identity mismatch'
    assert t(0) == QQ(15122)/33609, 'incorrect common base fiber'
    assert t.derivative() != 0, 'parameter map constant'
    return t


def torsion_fiber(parent_data):
    a4 = PolynomialRing(QQ,'T')(parent_data['jacobian_model']['a4_low_to_high'])
    a6 = PolynomialRing(QQ,'T')(parent_data['jacobian_model']['a6_low_to_high'])
    t0 = QQ(15122)/33609
    A, B = a4(t0), a6(t0)
    curve = EllipticCurve(QQ, [0,0,0,A,B])
    assert curve.discriminant() != 0
    result = []
    for p, expected in [(31,42),(37,46),(59,71)]:
        assert (4*GF(p)(A)**3 + 27*GF(p)(B)**2) != 0
        reduced = EllipticCurve(GF(p), [GF(p)(0),GF(p)(0),GF(p)(0),GF(p)(A),GF(p)(B)])
        n = reduced.cardinality()
        assert n == expected, (p,n,expected)
        result.append(n)
    from math import gcd
    assert gcd(*result) == 1
    print('Torsion at t=15122/33609: good reduction orders 42,46,71, gcd 1: PASS')
    return curve


def main():
    pdata = load('parent-rank17.json')
    pR,pK,pE,pP = build(pdata)
    assert len(pP)==17
    pH = check_model('K3 parent',pdata,pE,pP,pR,2,1152)
    for index in (4339,36773):
        data = load(f'cover-{index}.json')
        R,K,E,P = build(data)
        assert len(P)==18
        check_cover_identity(pdata,data,R)
        G = check_model(f'rational cover {index}',data,E,P,R,4,452984832)
        assert G.submatrix(0,0,17,17) == 2*pH, 'inherited-height block mismatch'
        assert G.submatrix(0,0,17,17).det() == 150994944
        assert G.det()/G.submatrix(0,0,17,17).det() == 3
        print(f'  cover {index}: inherited block = 2*parent; extra Schur complement 3: PASS')
    torsion_fiber(pdata)
    print('CERTIFIED: two nonconstant families over Q(U) with trivial torsion and generic rank >=18.')
    print('NOTE: This is a lower bound, not a proof that the exact generic rank is 18.')


if __name__ == '__main__':
    main()
