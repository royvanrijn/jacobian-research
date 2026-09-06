#!/usr/bin/env sage-python
"""Exact index-three enlargement and saturation of the constructed 302 span.

No parameter search. Replay the pinned baseline first. The only enumeration
is the 512 parity vectors of a fixed rank-nine lattice; short vectors supply
witnesses, and completeness of short-vector enumeration is not required.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

from sage.all import (EllipticCurve, PolynomialRing, QQ, ZZ, identity_matrix,
                      matrix, pari, vector)
from sage.env import SAGE_VERSION

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/generated-results/elkies-k3-curve302-nine-direction-k3-v1.json'
OUTPUT = ROOT / 'artifacts/generated-results/elkies-k3-curve302-section-saturation-v1.json'
BASE_HASH = '51c9ef1f87dcb31effb83ab71f5b50a18ee56c351b8d4764f973143e43b64a77'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def rows(m):
    return [list(map(int, r)) for r in m.rows()]


def rational(f):
    return {'numerator': list(map(str, f.numerator().list())),
            'denominator': list(map(str, f.denominator().list()))}


def build():
    if digest(BASE) != BASE_HASH:
        raise ArithmeticError('baseline changed; replay and review dependency')
    old = json.loads(BASE.read_text())
    G = matrix(ZZ, old['generic_section_lattice']['height_gram'])
    T = identity_matrix(QQ, 9)
    # Q is the third intersection of the line P1 P2 with the source cubic.
    # Q = (C1+...+C8)/3-C1 in the pointed cubic group.
    T[7] = vector(QQ, [-QQ(2)/3]+[QQ(1)/3]*7+[0])
    H = matrix(ZZ, T*G*T.T)
    source_H = matrix(ZZ, H[:8, :8]/2)
    assert T.det() == QQ(1)/3 and H.det() == 512
    assert source_H.det() == 1 and source_H.is_positive_definite()
    assert all(source_H[i, i] % 2 == 0 for i in range(8))
    image = matrix(ZZ, T*matrix(ZZ, old['specialization_overlap']['matrix_9_by_31_rows']))
    assert image.rank() == 9
    assert image.elementary_divisors() == [1]*8+[3]
    assert list(image[7]) == [-2, -1]+[0]*29

    # Every proper even integral overlattice has an order-two first step:
    # det H is a power of two. Enumerate every possible such coset exactly.
    candidates = []
    for bits in product(range(2), repeat=9):
        v = vector(ZZ, bits)
        if any(v) and all(c % 2 == 0 for c in H*v) and v*H*v % 8 == 0:
            candidates.append(tuple(bits))
    short = matrix(ZZ, pari(H).qfminim(8)[2]).columns()
    witnesses = {}
    for w in short:
        key = tuple(int(c % 2) for c in w)
        if w*H*w == 8 and key in candidates:
            ww = tuple(map(int, w))
            witnesses[key] = min(witnesses.get(key, ww), ww)
    assert len(candidates) == 71 and set(witnesses) == set(candidates)
    # A half of each witness has height 2, impossible on a 24I1 K3, where
    # every nonzero section has height 4+2(P.O). This excludes every coset.
    for key, w in witnesses.items():
        w = vector(ZZ, w)
        assert tuple(int(c % 2) for c in w) == key and w*H*w == 8
    print('SAT302|index=3|det=512|all_71_overlattice_cosets_excluded', flush=True)

    U = PolynomialRing(QQ, 'u'); u = U.gen(); K = U.fraction_field()
    R = PolynomialRing(K, names=('X', 'Y', 'Z')); x, y, z = R.gens()
    c = old['pointed_cubic_family']
    def polynomial(records):
        return sum(K(e['coefficient'])*x**e['exponents'][0]*y**e['exponents'][1]*z**e['exponents'][2]
                   for e in records)
    F0, F1 = polynomial(c['F0']), polynomial(c['F1'])
    F = F0+u*F1
    P = vector(K, c['zero']); P /= P[2]
    constants = [vector(K, p) for p in c['constant_sections']]
    V = PolynomialRing(K, 'v'); v = V.gen()
    line = P+v*(constants[0]-P)
    restriction = V(F(*line))
    assert restriction(0) == restriction(1) == 0 and restriction.degree() == 3
    lam = -restriction[2]/restriction[3]-1
    Q = P+lam*(constants[0]-P)
    assert F(*Q) == 0
    basis = constants[:7]+[Q]
    print('SAT302|residual_line_section=exact', flush=True)

    # Explicit non-flex cubic conversion. P2 is the tangent residual at P;
    # P3 need only be a point on the tangent at P2, not on the cubic.
    dx, dy = (F.derivative(w)(*P) for w in (x, y))
    tangent_direction = vector(K, [1, -dx/dy, 0])
    restriction = V(F(*(P+v*tangent_direction)))
    assert restriction[0] == restriction[1] == 0 and restriction[3]
    P2 = P-restriction[2]/restriction[3]*tangent_direction
    assert F(*P2) == 0 and P2 != P
    dx, dy = (F.derivative(w)(*P2) for w in (x, y))
    P3 = vector(K, [1, -dx/dy, 0])
    M = matrix(K, [P, P2, P3]).T
    assert M.det()
    F2 = R(M.act_on_polynomial(F))
    pulled = R(F2(x*x, y*z, x*z))
    assert all(e[0] >= 2 and e[2] >= 1 for e in pulled.dict())
    transformed = R({(e[0]-2, e[1], e[2]-1): v for e, v in pulled.dict().items()})
    assert transformed*x*x*z == pulled
    a = K(transformed.monomial_coefficient(x**3))
    b = K(transformed.monomial_coefficient(y*y*z))
    ab = a*b
    W = R(transformed(-x, y/b, ab*z)/a)
    E = EllipticCurve(W(x, y, 1))
    assert W == E.defining_polynomial()(x, y, z)
    MI = M.inverse()
    source_a = U(old['source_jacobian']['A_coefficients_low_to_high'])
    source_b = U(old['source_jacobian']['B_coefficients_low_to_high'])
    J = EllipticCurve(K, [source_a, source_b])
    compact_iso = E.isomorphism_to(J)
    # Forward coordinates in the M frame are (-xz, bxy, z^2/(ab)).
    # The verified transformations are invertible on a dense open set.
    def on_weierstrass(p):
        px, py, pz = MI*p
        xx, yy = -ab*px/pz, a*b*b*px*py/pz**2
        # Verify the final compact point equation directly; the polynomial
        # transformation above already verifies the intermediate model.
        us, rs, ss, ts = compact_iso.tuple()
        return J([(xx-rs)/us**2, (yy-ss*(xx-rs)-ts)/us**3])
    wbasis = []
    for i, p in enumerate(basis):
        assert F(*p) == 0
        wbasis.append(on_weierstrass(p))
        print('SAT302|generic_Weierstrass_section=%s' % (i+1), flush=True)

    # Choose the fibre isomorphism by the known translation, resolving its
    # possible sign with all eight points, rather than by j equality alone.
    n = list(map(QQ, c['N_coefficients_low_to_high']))
    d = list(map(QQ, c['D_coefficients_low_to_high']))
    n += [QQ(0)]*(3-len(n)); d += [QQ(0)]*(3-len(d))
    branch = U((n[1]-u*d[1])**2-4*(n[2]-u*d[2])*(n[0]-u*d[0]))
    assert branch.degree() == 2 and branch.gcd(branch.derivative()).degree() == 0
    source_delta = -16*(4*source_a**3+27*source_b**2)
    assert branch.gcd(source_delta).degree() == 0
    E0 = EllipticCurve(QQ, [f(0) for f in J.a_invariants()])
    f0 = F0.change_ring(QQ)
    fx, fy, fz = f0.parent().gens()
    short = EllipticCurve(QQ, [-f0.monomial_coefficient(fx*fz*fz), -f0.monomial_coefficient(fz**3)])
    p0 = short([v(0) for v in P])
    at_zero = [E0([v(0) for v in p]) for p in wbasis]
    expected = [short([v(0) for v in p])-p0 for p in basis]
    matches = [iso for iso in E0.isomorphisms(short)
               if all(iso(q) == e for q, e in zip(at_zero, expected))]
    assert len(matches) == 1
    target_iso = matches[0]
    assert expected[-1] == -2*p0-short([v(0) for v in constants[0]])
    public9 = vector(K, [sum(map(QQ, cs)) for cs in c['moving_section_coordinates_low_to_high']])
    # public9 is only on the u=0 fibre: evaluate the map before point testing.
    px, py, pz = [f(0) for f in MI*public9]
    xx, yy = -ab(0)*px/pz, (a*b*b)(0)*px*py/pz**2
    us, rs, ss, ts = [f(0) for f in compact_iso.tuple()]
    ninth_at_zero = E0([(xx-rs)/us**2, (yy-ss*(xx-rs)-ts)/us**3])
    assert target_iso(ninth_at_zero) == short([f(0) for f in public9])-p0

    # The K3 uses this same explicit gauge with u=N(t)/D(t). Store the
    # straight-line evaluation map, avoiding needless high-degree expansion.
    # Its ninth cubic section was verified identically in the pinned baseline.
    result = {
        'schema': 'elkies-k3.curve302-section-saturation.v1',
        'status': 'PASS_EXACT_SOURCE_MW8_BASIS_AND_SATURATED_K3_RANK9_SPAN',
        'input_sha256': {str(BASE.relative_to(ROOT)): digest(BASE)},
        'source_sha256': digest(Path(__file__)),
        'software': {'sage': SAGE_VERSION},
        'source_basis': {
            'rank': 8, 'full_geometric_and_arithmetic_MW_basis': True,
            'torsion_order': 1, 'height_gram': rows(source_H), 'determinant': 1,
            'cubic_sections': ['P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'Q12(u)'],
            'Q12_coordinates': [rational(p) for p in Q],
            'Weierstrass_a_invariants': [rational(f) for f in J.a_invariants()],
            'Weierstrass_basis_xy': [[rational(f) for f in p.xy()] for p in wbasis],
            'specialization_parameter': '0',
            'fibre_to_short_u_r_s_t': list(map(str, target_iso.tuple())),
            'short_to_public': 'x=(X-15)/36; y=(Y/108-x-1)/2',
        },
        'cubic_to_Weierstrass': {
            'inverse_frame_matrix': [[rational(f) for f in row] for row in MI.rows()],
            'a': rational(a), 'b': rational(b),
            'intermediate_to_compact_u_r_s_t': [rational(f) for f in compact_iso.tuple()],
            'rule': 'Set (v,w,z)=M_inverse(u)*(X,Y,Z); x=-a*b*v/z; y=a*b^2*v*w/z^2. Then apply stored (k,r,s,q): x_c=(x-r)/k^2, y_c=(y-s*(x-r)-q)/k^3. Extend across removable points; P1 maps to O.',
            'polynomial_transformation_identity': True,
        },
        'K3': {
            'base_change': 'u=N(t)/D(t), with N,D and ninth cubic section L(t) from pinned input',
            'Weierstrass_equation': 'the polynomial Weierstrass family A_K3,B_K3 in the pinned baseline',
            'sections': 'Pull back the eight source basis points and cubic_to_Weierstrass(L(t)) along u=N(t)/D(t), then multiply x by D(t)^2 and y by D(t)^3',
            'specialization_parameter': '1',
            'specialization_to_source_fibre': 'Divide x,y by D(1)^2,D(1)^3, then apply the stored source fibre isomorphism and short-to-public map',
            'ninth_specialization_transport_checked': True,
            'height_gram': rows(H), 'determinant': 512,
            'old_in_new_index': 3,
            'new_in_old_rational_matrix': [list(map(str, row)) for row in T.rows()],
            'saturated_in_full_geometric_MW': True,
            'rank': 9, 'full_generic_rank': 'UNKNOWN',
            'full_generic_arithmetic_rank_bounds': [9, 17],
            'specialization_matrix': rows(image),
            'specialization_smith_factors': list(map(int, image.elementary_divisors())),
            'displayed_public_quotient': 'Z^22 + Z/3',
            'torsion_order': 1,
        },
        'saturation_certificate': {
            'candidate_cosets': len(candidates),
            'height_two_witnesses': [{'parity': list(key), 'twice_vector': list(witnesses[key])}
                                     for key in sorted(witnesses)],
            'proof': 'det=2^9, so any proper even integral overlattice has an order-two first step. Every admissible coset contains a vector of height 2. All fibres are I1, so a nonzero K3 section has height 4+2(P.O)>=4. No such step exists.',
        },
        'remaining_rank_gate': {
            'branch_polynomial_coefficients_low_to_high': list(map(str, branch.list())),
            'branch_degree': 2, 'branch_squarefree': True,
            'branch_disjoint_from_source_singular_fibres': True,
            'twist_equation': 'y^2=x^3+branch(u)^2*A_source(u)*x+branch(u)^3*B_source(u)',
            'rank_identity': 'rank K3(Q(t)) = 8 + rank source_twist(Q(u))',
            'known_anti_invariant': '2*S-h, where h=Q+C1 and S is the moving ninth section',
            'known_anti_invariant_K3_height': 8,
            'full_twist_rank': 'UNKNOWN; at least 1',
        },
        'proof_boundary': 'This completes a generic MW8 basis for an imposed rational pencil and saturates the known K3 rank-nine span. It neither proves full K3 rank nine nor recovers hidden provenance or a high-rank parent. No new rational curve rank or parameter search.',
        'replay': 'sage -python elkies-k3/scripts/construct_curve302_nine_direction_k3.sage --check && sage -python elkies-k3/scripts/certify_curve302_section_saturation.sage --check',
    }
    print('SAT302|source_MW_basis=8|K3_saturated_span=9|full_K3_rank=UNKNOWN', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    rendered = json.dumps(build(), indent=2, sort_keys=True)+'\n'
    if args.check:
        assert args.output.read_text() == rendered, 'certificate differs from exact replay'
    else:
        args.output.write_text(rendered)


if __name__ == '__main__':
    main()
