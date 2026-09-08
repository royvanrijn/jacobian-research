#!/usr/bin/env sage-python
"""Independent parent-point, residue, and bielliptic-map replay.

No constructor import; local root counts checked using Frobenius gcds.
"""
import hashlib
import json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, vector, matrix, power_mod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
PREFIX = 'det1092_rr_fibre_cover_gate'
CERT = ART / (PREFIX + '_v1.json')
PROTO = ART / (PREFIX + '_protocol_v1.json')
OUT = ART / (PREFIX + '_replay_v1.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def check():
    d, p = read(CERT), read(PROTO)
    assert d['protocol_sha256'] == sha(PROTO)
    for path, h in p['inputs'].items():
        assert sha(ROOT/path) == h
    parent_path = ART / 'curve302_recovered_mw17_parent_v1.json'
    net_path = ART / 'det1092_first_centre_rr_net_v1.json'
    parent, net = read(parent_path), read(net_path)
    R = PolynomialRing(QQ, 't')
    K = R.fraction_field()
    def dec(v):
        return K(R(v['numerator']))/R(v['denominator'])
    E = EllipticCurve(K, [dec(v) for v in parent['a_invariants']])
    base = [E([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    word = -vector(ZZ, net['trace_word'])
    assert word*matrix(QQ, parent['generic_height_gram'])*word == 10
    C = sum((n*P for n, P in zip(word, base)), E(0))
    cx = C[0]+E.b2()/12
    cy = C[1]+(E.a1()*C[0]+E.a3())/2
    aa, bb = -E.c4()/48, -E.c6()/864
    # Independently recover the unchanged roster and check fibre j-invariants.
    roster_path = ROOT / 'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v1/roster.json'
    transport_path = ART / 'det1092_reduced_parameter_chart_v1/generic-proof.json'
    a, b, c, dd = map(QQ, read(transport_path)['parameter_matrix'])
    expected = {}
    for row in read(roster_path):
        s = QQ(row['parameter'])
        t = s if row['presentation'] == 'normalized' else (a*s+b)/(c*s+dd)
        key = 'case-'+hashlib.sha256(str(t).encode()).hexdigest()[:12]
        expected[key] = (t, row['id'])
        if row['presentation'] == 'reduced':
            assert E.j_invariant()(t) == EllipticCurve(QQ, list(map(QQ, row['model']))).j_invariant()
    assert len(expected) == len(d['rows']) == 11
    assert p['cases'] == [{'case': r['case'], 'parameter': r['parameter']} for r in d['rows']]
    summary = []
    for row in d['rows']:
        t, label = expected[row['case']]
        assert t == QQ(row['parameter'])
        x, y, a0, b0 = cx(t), cy(t), aa(t), bb(t)
        assert list(map(QQ, row['centre'])) == [x, y]
        assert list(map(QQ, row['short_curve'])) == [a0, b0]
        assert y*y == x**3+a0*x+b0 and y and 4*a0**3+27*b0*b0
        coeff = [-3*x*x-4*a0, -8*y, -6*x, QQ(0), QQ(1)]
        assert list(map(QQ, row['halving_coefficients'])) == coeff
        assert [r['prime'] for r in row['trials']] == p['primes']
        obstruction_primes = []
        for trial in row['trials']:
            prime = ZZ(trial['prime'])
            assert prime.is_prime() and prime <= 197
            excluded = any(v.denominator() % prime == 0 for v in coeff)
            if excluded:
                assert trial['status'] == 'EXCLUDED_NONINTEGRAL_CHART'
                continue
            S = PolynomialRing(GF(prime), 'X')
            X = S.gen()
            f = S(coeff)
            roots_polynomial = f.gcd(power_mod(X, prime, f)-X)
            assert roots_polynomial.degree() == len(trial['roots'])
            assert all(f(z) == 0 for z in trial['roots'])
            assert len(set(trial['roots'])) == len(trial['roots'])
            assert list(map(int, f.list())) == trial['coefficients']
            assert bool(f.discriminant()) == trial['squarefree']
            no_root = roots_polynomial.degree() == 0
            assert (trial['status'] == 'NO_LOCAL_HALF') == no_root
            if no_root:
                obstruction_primes.append(int(prime))
        assert obstruction_primes == row['obstruction_primes']
        assert row['status'] == ('Q_BRANCH_COVER_OBSTRUCTED' if obstruction_primes else 'UNKNOWN')
        summary.append({'case': row['case'], 'roster_id': label,
                        'status': row['status'], 'obstruction_primes': obstruction_primes})
    # Prove quartic irreducibility mod17 without calling a factoring routine.
    S = PolynomialRing(GF(17), 'X')
    X = S.gen()
    f = X**4+2*X**2+8*X+14
    assert f == S([-3*cx(0)**2-4*aa(0), -8*cy(0), -6*cx(0), 0, 1])
    assert f.gcd(power_mod(X, 17**2, f)-X) == 1
    assert power_mod(X, 17**4, f)-X == 0
    # Independent coefficient proof over a general short elliptic curve.
    B = PolynomialRing(QQ, ['a', 'cx', 'cy', 'n'])
    a, xC, yC, n = B.gens()
    F = B.fraction_field()
    M = PolynomialRing(F, 'm')
    m = M.gen()
    H = m**4-6*xC*m*m-8*yC*m-3*xC*xC-4*a
    qx = (m*m-xC)/2
    qy = m*(qx-xC)-yC
    b = yC*yC-xC**3-a*xC
    assert (qy*qy-qx**3-a*qx-b) % H == 0
    assert (2*m*qy-3*qx*qx-a) % H == 0
    assert m*m-2*qx == xC and m*(qx-xC)-qy == yC
    L = n-m
    numerator_e = qx*L+2*qy
    branch_identity = L*(numerator_e**3+a*numerator_e*L*L+b*L**3)-qy*qy*(n**4-6*xC*n*n-8*yC*n-3*xC*xC-4*a)
    assert branch_identity % H == 0
    # Independent general even-sextic smoothness and second quotient identity.
    T = PolynomialRing(QQ, ['a', 'b', 'e'])
    a, b, e = T.gens()
    Z = PolynomialRing(T, 'z')
    z = Z.gen()
    fe = e**3+a*e+b
    sextic = (z*z+e)**3+a*(z*z+e)+b
    assert sextic.discriminant() == -64*fe*(-4*a**3-27*b*b)**2
    assert fe**2*sextic == fe**3+(3*e*e+a)*fe**2*z*z+3*e*fe**2*z**4+fe**2*z**6
    # Formula for the chord map after translating by Q, by cleared addition.
    V = PolynomialRing(QQ, ['m', 'qx', 'qy', 'x', 'y'])
    m, qx, qy, x, y = V.gens()
    a = 2*m*qy-3*qx*qx
    b = qy*qy-qx**3-a*qx
    xC, yC = m*m-2*qx, m*(3*qx-m*m)-qy
    dx = x-qx
    px_num = (y-qy)**2-(x+qx)*dx**2
    py_num = -qy*dx**3+(y-qy)*(qx*dx**2-px_num)
    chord = py_num+yC*dx**3-(m*dx+2*qy)*(px_num-xC*dx**2)
    assert chord == (-y+qy-m*dx)*(y*y-x**3-a*x-b)
    return {'classification': 'verified application and new deduction',
            'status': 'PASS_INDEPENDENT_RR_FIBRE_COVER_GATE',
            'inputs': {str(path.relative_to(ROOT)): sha(path) for path in
                       [CERT, PROTO, parent_path, net_path, roster_path, transport_path]},
            'checker_sha256': sha(Path(__file__)),
            'summary_after_arithmetic': summary,
            'verified_identities': ['parent centre', 'roster parameter transports and fibre j',
                                    'all modular root counts via Frobenius gcd',
                                    'quartic irreducibility mod17', 'doubling modulo H',
                                    'branch-value identity modulo H', 'sextic discriminant',
                                    'two elliptic quotient maps', 'translated chord map'],
            'zero_minimum_half_field_degree': 4,
            'point_searches': 0, 'Selmer_runs': 0, 'class_group_runs': 0,
            'boundary': 'Excludes only the one-RR-pair branch construction over Q; it is not a genus2 Selmer panel or a rank-incidence theorem.'}


if __name__ == '__main__':
    data = check()
    text = json.dumps(data, indent=2, sort_keys=True)+'\n'
    if OUT.exists():
        assert OUT.read_text() == text
    else:
        OUT.write_text(text)
    print(data['status'], flush=True)
