#!/usr/bin/env sage-python
"""Independent point identity, generic-source transport and divisor replay."""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
CERT = ART / 'det1092_rr_inherited_point_v1.json'
OUT = ART / 'det1092_rr_inherited_point_replay_v1.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify():
    d = json.loads(CERT.read_text())
    for path, digest in d['inputs'].items():
        assert sha(ROOT/path) == digest
    pp = ART/'curve302_recovered_mw17_parent_v1.json'
    np = ART/'det1092_first_centre_rr_net_v1.json'
    up = ART/'det1092_universal_rr_descent_preflight_v1.json'
    rp = ART/'det1092_rr_net_reducible_locus_v2.json'
    parent, net, uni, red = [json.loads(p.read_text()) for p in [pp, np, up, rp]]
    uni = uni['universal_genus2']
    R = PolynomialRing(QQ, 'T')
    T = R.gen()
    K = R.fraction_field()
    def rat(row):
        return K(R(row['numerator']))/R(row['denominator'])
    E = EllipticCurve(K, [rat(row) for row in parent['a_invariants']])
    basis = [E([rat(row) for row in point]) for point in parent['basis_weierstrass_coordinates']]
    w = -vector(QQ, net['trace_word'])
    C = sum((n*P for n, P in zip(w, basis)), E(0))
    gram = matrix(QQ, parent['generic_height_gram'])
    assert w*gram*w == 10
    h = R(uni['h'])
    assert h.degree() == 3 and h**2 == C[0].denominator()
    A, B = [[R(row) for row in net[key]] for key in ['A', 'B']]
    p = d['point_formula']
    alpha, beta0, beta1 = [QQ(p[key]) for key in ['alpha', 'beta0', 'beta1']]
    assert alpha and A[2] == alpha*h and B[2] == (beta0+beta1*T)*h
    ue, ve = -beta0/alpha, -beta1/alpha
    assert list(map(QQ, p['unique_base_parameter'])) == [ue, ve]
    co = list(map(QQ, red['canonical_O_plus_C_member']['coefficients_in_A_tA_B']))
    assert co[2] and [co[0]/co[2], co[1]/co[2]] == [ue, ve]
    vertical = [B[i]+(ue+ve*T)*A[i] for i in range(3)]
    g = K(vertical[1])/h**2
    assert g in QQ and g != 0 and vertical[2] == 0
    assert vertical[0]+vertical[1]*C[0] == 0
    # Rebuild the source sextic identity, without calling the constructor.
    UV = PolynomialRing(QQ, ['u', 'v'])
    u, v = UV.gens()
    S = PolynomialRing(UV, 'T')
    t = S.gen()
    f1, f2 = [S(B[i])+(u+v*t)*S(A[i]) for i in [1, 2]]
    hh = S(h)
    cx, cy = C[0]+E.b2()/12, C[1]+(E.a1()*C[0]+E.a3())/2
    nx, ny, aa = S(R(cx*h**2)), S(R(cy*h**3)), S(R(-E.c4()/48))
    m = -f1+S(R(E.a1()))*f2/2
    q = sum(QQ(row['coefficient'])*u**row['u']*v**row['v']*t**row['T']
            for row in uni['sparse_q'])
    c = QQ(uni['scale'])
    assert (hh**4*m**4-6*nx*hh**2*m**2*f2**2-8*ny*hh*m*f2**3
            -(3*nx**2+4*aa*hh**4)*f2**4) == hh**10*c*q
    # Homogeneous evaluation handles rational points over T=infinity too.
    tt, zz = -beta0-alpha*u, beta1+alpha*v
    def read_sparse(rows):
        return sum(QQ(r['coefficient'])*u**r['u']*v**r['v'] for r in rows)
    def at_point(poly, degree):
        return sum(UV(poly[i])*tt**i*zz**(degree-i) for i in range(degree+1))
    H = at_point(S(h), 3)
    F = at_point(f1, 6)
    W = H*QQ(g)**2
    assert F == H**2*QQ(g)
    assert W**2 == c*at_point(q, 6)
    for key, value in [('t0', tt), ('z0', zz), ('H', H), ('G', UV(g)), ('W', W)]:
        assert read_sparse(p[key]) == value
    assert tt(ue, ve) == zz(ue, ve) == 0
    assert tt.derivative(u) != 0 and zz.derivative(v) != 0
    fixed = R([a(0, 0) for a in q.list()])
    x, y = map(QQ, p['fixed_0_0_point'])
    assert fixed.degree() == 6 and fixed.gcd(fixed.derivative()) == 1
    assert y != 0 and y*y == c*fixed(x)
    assert x == tt(0, 0)/zz(0, 0) and y == W(0, 0)/zz(0, 0)**3
    assert net['generic_class']['D_dot_O'] == 2*(-2)+5 == 1
    # Reuse the certified pairs; no new enumeration or point search.
    words = [vector(QQ, word) for pair in red['section_pairs'] for word in pair['words']]
    assert len(words) == len(set(tuple(word) for word in words)) == 44
    assert all(1+x*gram*x-w*gram*x == 1 for x in words)
    return {
        'classification': 'verified application and new deduction',
        'status': 'PASS_INDEPENDENT_UNIFORM_INHERITED_RR_POINT',
        'generic_source_sextic_identity': True,
        'uniform_homogeneous_point_identity': True,
        'unique_formula_basepoint_is_old_reducible_member': True,
        'fixed_0_0_rational_point_verified': True,
        'generic_sections_with_degree_one_intersection': 44,
        'distinctness_boundary': '44 distinct points generically; special intersections can coincide.',
        'rational_and_everywhere_local_solubility': 'PROVED on every smooth integral rational RR member',
        'curve_fake_2_Selmer_set_nonempty': 'PROVED by the inherited rational point, not a Selmer computation',
        'Jacobian_2_Selmer_and_Cassels_Tate': 'NOT_COMPUTED',
        'new_MW_direction': 'NOT_CLAIMED',
        'limits': {'wall_seconds': 25, 'point_searches': 0, 'Selmer_runs': 0,
                   'class_group_runs': 0, 'control_sweeps': 0, 'pilot_changes': 0},
        'inputs': {str(path.relative_to(ROOT)): sha(path) for path in [CERT, pp, np, up, rp]},
        'checker_sha256': sha(Path(__file__))
    }


if __name__ == '__main__':
    signal.alarm(25)
    result = verify()
    payload = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if OUT.exists():
        assert OUT.read_text() == payload
    else:
        OUT.write_text(payload)
    print(result['status'], flush=True)
