#!/usr/bin/env sage-python
"""Construct the exact branch-first incidence curve for one frozen pencil.

All finite slopes through its first basepoint are retained symbolically.
Rational branch points require Y^2=R8(m); a rational residual point on302
also requires H^2=D4(m). Their normalized fibre product has genus nine.
This is a construction of the incidence curve, not a rational-point search.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, prod

ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT/'elliptic-curves/cas/icarm_curve302.py'
SOURCE = ROOT/'artifacts/generated-results/elkies-k3-curve302-cubic-pencil-overlap-v1.json'
OUT = ROOT/'artifacts/generated-results/elkies-k3-curve302-rational-branch-carrier-v1.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return list(map(str, poly.list()))


def compute():
    source = json.loads(SOURCE.read_text())
    pencil = source['pencils'][0]
    path = ROOT/pencil['checkpoint']
    assert digest(path) == pencil['checkpoint_sha256']
    assert digest(PUBLIC) == source['inputs'][str(PUBLIC.relative_to(ROOT))]
    checkpoint = json.loads(path.read_text())
    assert checkpoint['anchors_one_based'] == pencil['anchors_one_based'] == list(range(1, 9))
    data = runpy.run_path(str(PUBLIC))
    target = EllipticCurve(QQ, list(map(lambda c: QQ(str(c)), data['GENERAL_WEIERSTRASS_COEFFICIENTS'])))
    c4, c6 = target.c_invariants()
    a, b = -27*c4, -54*c6
    curve = EllipticCurve(QQ, [a, b])
    points = [curve(36*QQ(str(x))+15, 108*(2*QQ(str(y))+QQ(str(x))+1)) for x, y in data['POINTS']]
    base = points[:8]+[-sum(points[:8], curve(0))]
    assert len(set(base)) == 9 and all(q[2] == 1 for q in base)
    point = base[0]
    terms = checkpoint['F1']

    def f1(x, y):
        return sum(QQ(t['coefficient'])*x**t['exponents'][0]*y**t['exponents'][1] for t in terms)

    assert all(f1(q[0], q[1]) == 0 for q in base)
    mr = PolynomialRing(QQ, 'm')
    m = mr.gen()
    sr = PolynomialRing(mr, 's')
    s = sr.gen()
    x, y = point[0]+s, point[1]+m*s
    g0, r0 = (y*y-x**3-a*x-b).quo_rem(s)
    g1, r1 = sr(f1(x, y)).quo_rem(s)
    assert not r0 and not r1 and g0.degree() == g1.degree() == 2
    dc = g0[1]**2-4*g0[2]*g0[0]
    db = 2*g0[1]*g1[1]-4*(g0[2]*g1[0]+g1[2]*g0[0])
    da = g1[1]**2-4*g1[2]*g1[0]
    result = g0.resultant(g1)
    assert db*db-4*da*dc == 16*result
    assert result.degree() == 8 and dc.degree() == 4
    assert dc == m**4-6*point[0]*m*m+8*point[1]*m-3*point[0]**2-4*a
    assert result.gcd(result.derivative()).degree() == 0
    assert dc.gcd(dc.derivative()).degree() == 0
    assert result.gcd(dc).degree() == 0
    assert da != 0
    exceptional_checks = {}
    for label, polynomial in [('da', da), ('D4', dc)]:
        reduction = polynomial.change_ring(GF(17))
        assert reduction.degree() == polynomial.degree() == 4
        assert all(reduction(z) != 0 for z in GF(17))
        exceptional_checks[label] = {
            'prime': 17, 'degree': 4,
            'coefficients_low_to_high': list(map(int, reduction.list())),
            'projective_root_count': 0}
    slopes = [(q[1]-point[1])/(q[0]-point[0]) for q in base[1:]]
    assert len(set(slopes)) == 8
    assert result.monic() == prod(m-v for v in slopes)
    known = []
    for q, slope in zip(base[1:], slopes):
        h = 2*g0[2](slope)*(q[0]-point[0])+g0[1](slope)
        assert h*h == dc(slope) and result(slope) == 0
        assert point[0]+(-g0[1](slope)+h)/(2*g0[2](slope)) == q[0]
        known.append({'m': str(slope), 'Y': '0', 'H': str(h),
                      'point_on_short302': list(map(str, q)),
                      'excluded_reason': 'Line contains another pencil basepoint; resultant zero.'})

    # Check the missing projective slope directly, rather than inferring its
    # arithmetic from a leading coefficient in a singular affine model.
    vr = PolynomialRing(QQ, 'v')
    v = vr.gen()
    vertical0, rem0 = vr((point[1]+v)**2-point[0]**3-a*point[0]-b).quo_rem(v)
    vertical1, rem1 = vr(f1(point[0], point[1]+v)).quo_rem(v)
    assert not rem0 and not rem1
    vc = vertical0[1]**2-4*vertical0[2]*vertical0[0]
    vb = 2*vertical0[1]*vertical1[1]-4*(vertical0[2]*vertical1[0]+vertical1[2]*vertical0[0])
    va = vertical1[1]**2-4*vertical1[2]*vertical1[0]
    vd = vb*vb-4*va*vc
    assert va != 0 and vd != 0 and not vd.is_square()
    return {'schema': 'curve302.rational-branch-carrier.v1',
            'status': 'EXPLICIT_GENUS_NINE_INCIDENCE_CURVE_NO_POINT_SEARCH',
            'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), PUBLIC, SOURCE, path]},
            'source_pencil_zero_based': 0, 'basepoint_index_zero_based': 0,
            'short302_coefficients_a_b': [str(a), str(b)], 'basepoint': list(map(str, point)),
            'F1': terms,
            'residual_g0_coefficients_in_s': [coefficients(c) for c in g0.list()],
            'residual_g1_coefficients_in_s': [coefficients(c) for c in g1.list()],
            'branch_discriminant_in_u': {'constant': coefficients(dc), 'linear': coefficients(db), 'quadratic': coefficients(da)},
            'R8_low_to_high': coefficients(result), 'D4_low_to_high': coefficients(dc),
            'R8_leading_coefficient': str(result.leading_coefficient()),
            'R8_rational_roots': list(map(str, slopes)),
            'squarefree_degrees': [8, 4], 'gcd_degree': 0,
            'genus_Y_squared_R8': 3, 'genus_H_squared_D4': 1,
            'geometric_branch_union_size': 12, 'normalized_compositum_genus': 9,
            'known_degenerate_points': known,
            'vertical_branch_discriminant': str(vd), 'vertical_rational_branch_points': False,
            'no_rational_roots_of_da_or_D4': exceptional_checks,
            'identity': 'disc_u(disc_s(g0+u*g1))=16*Res_s(g0,g1)=16*R8(m)',
            'point_recovery': 's=(-g0[1](m)+H)/(2*g0[2](m)); Q=(xP+s,yP+m*s)',
            'branch_recovery': 'u=(-db(m)+/-4*Y)/(2*da(m))',
            'open_conditions': ['R8(m) != 0', 'D4(m) != 0', 'da(m) != 0'],
            'equivalence_scope': 'On this open set, rational points (m,Y,H) correspond to finite pointed lines with two rational branch points and rational residual points on302, with choices of branch order and residual root. Both da=0 and D4=0 have no rational slopes by the projective mod17 certificates. The vertical line has nonrational branch points. Thus only the eight known degenerate slopes R8=0 are removed from the rational-line incidence problem.',
            'boundary': 'All slopes for this one basepoint and pencil are accounted for symbolically. Rational points off the eight degenerate slopes remain UNKNOWN. This constructs an incidence curve, not a new elliptic surface, full MW basis, or parent.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    signal.alarm(120)
    rendered = json.dumps(compute(), indent=2, sort_keys=True)+'\n'
    if args.check:
        assert OUT.read_text() == rendered, 'Certificate changed'
    else:
        assert not OUT.exists(), 'Preserve existing certificate; use --check'
        OUT.write_text(rendered)
    print('PASS_BRANCH_CARRIER_GENUS_3_1_9', flush=True)
