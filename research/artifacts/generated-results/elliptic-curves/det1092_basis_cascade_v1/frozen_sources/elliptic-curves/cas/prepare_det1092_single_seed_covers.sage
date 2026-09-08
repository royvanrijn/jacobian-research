#!/usr/bin/env sage-python
"""Freeze exactly two existing covers and nine existing parameter addresses.

The second cover is calibrated by the historical first-unlock RR pencil.
No V3 artifacts, later points, control outcomes, or new search are read.
"""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve, matrix, vector, version

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_single_seed_covers_v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def retain(p, data):
    text = json.dumps(data, indent=2, sort_keys=True)+'\n'
    if p.exists():
        assert p.read_text() == text
    else:
        p.write_text(text)


def build():
    paths = [ART/'curve302_recovered_mw17_parent_v1.json',
             ART/'det1092_first_centre_rr_net_v1.json',
             ART/'det1092_first_witness_pencil_genus_gate_v1.json',
             ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent, net, pencil, roster = [json.loads(p.read_text()) for p in paths]
    R = PolynomialRing(QQ, 't')
    K = R.fraction_field()

    def rat(r):
        return K(R(r['numerator']))/R(r['denominator'])

    def rec(f):
        f = K(f)
        return {'numerator': list(map(str, f.numerator().list())),
                'denominator': list(map(str, f.denominator().list()))}

    E = EllipticCurve(K, [rat(r) for r in parent['a_invariants']])
    basis = [E([rat(r) for r in p]) for p in parent['basis_weierstrass_coordinates']]
    w = -vector(QQ, net['trace_word'])
    Z = sum((n*P for n, P in zip(w, basis)), E(0))
    G = matrix(QQ, parent['generic_height_gram'])
    assert w*G*w == 10 and G.det() == 1092
    h = R(Z[0].denominator().sqrt())
    assert h.degree() == 3
    A, B = [[R(r) for r in net[k]] for k in ['A', 'B']]
    cx = Z[0]+E.b2()/12
    cy = Z[1]+(E.a1()*Z[0]+E.a3())/2
    a = -E.c4()/48
    cases = [{'index': i, 'parameter': r['parameter']}
             for i, r in enumerate(roster['cases'])]
    assert len(cases) == 9 and cases[-1]['parameter'] == '0'
    # No fitting or new member selection: these two members already exist.
    members = [QQ(0), QQ(pencil['u0'])]
    for i, u in enumerate(members):
        f0, f1, f2 = [B[j]+u*A[j] for j in range(3)]
        m = -f1/f2+E.a1()/2
        H = m**4-6*cx*m*m-8*cy*m-3*cx*cx-4*a
        F = R(f2**4*H/h**6)
        assert F.degree() == 6 and F.gcd(F.derivative()) == 1
        delta = R(E.discriminant())
        assert delta.degree() == 24 and delta.gcd(delta.derivative()) == 1
        assert F.gcd(delta) == 1
        x0 = (m*m-cx)/2-E.b2()/12
        x1 = h**3/(2*f2**2)
        y0 = m*(x0+E.b2()/12-cx)-cy-(E.a1()*x0+E.a3())/2
        y1 = (m-E.a1()/2)*x1
        t0 = -R(f2/h)[0]/R(f2/h)[1]
        assert h(t0) and F(t0) and F(t0).is_square()
        data = {
            'schema': 'det1092.equation-only-single-seed-cover.v1',
            'a_invariants': parent['a_invariants'],
            'generic_sections': parent['basis_weierstrass_coordinates'],
            'generic_height_gram': parent['generic_height_gram'],
            'trace_word': list(map(int, w)), 'trace_point': [rec(v) for v in Z[:2]],
            'line': [list(map(str, f.list())) for f in [f0, f1, f2]],
            'h': list(map(str, h.list())),
            'cover_polynomial': list(map(str, F.list())),
            'point_map': {'x': [rec(x0), rec(x1)], 'y': [rec(y0), rec(y1)]},
            'inherited_parameter': str(t0),
            'exclusion_factors': [list(map(str, delta.list())),
                                  list(map(str, F.list())), list(map(str, f2.list())),
                                  list(map(str, h.list()))],
            'cases': cases,
        }
        retain(DIR/f'cover-{i:02d}-input.json', data)
    protocol = {
        'classification': 'verified application / frozen calibrated replay protocol',
        'members': [
            {'index': 0, 'u': '0', 'v': '0', 'selection': 'existing generic-input member B'},
            {'index': 1, 'u': str(members[1]), 'v': '0',
             'selection': 'existing historical first-unlock member; retrospectively calibrated'}],
        'cases': [{'index': i, 'label': r['label'],
                   'parameter': r['parameter'], 'reduced_parameter': r['reduced_parameter']}
                  for i, r in enumerate(roster['cases'])],
        'limits': {'wall_seconds_per_command': 25, 'existing_members': 2,
                   'existing_parameters': 9, 'new_parameters': 0,
                   'modular_obstruction_primes': [3,5,7,11,13,17,19,23,29,31,37,41,43,47],
                   'point_searches': 0, 'class_groups': 0, 'Selmer_runs': 0,
                   'V3_artifacts_read': 0, 'pilot_changes': 0},
        'boundary': 'The execution payloads contain equations and generic sections, not exceptional points or chart hits. Cover01 still encodes historical calibration; replay on302 is not a held-out discovery. Controls are not used to choose either cover.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths+[Path(__file__)]},
        'payload_sha256': {f'cover-{i:02d}-input.json': sha(DIR/f'cover-{i:02d}-input.json') for i in range(2)},
        'software': version(),
    }
    retain(DIR/'protocol.json', protocol)
    print('PASS_TWO_EXISTING_SEED_COVERS_FROZEN', flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    DIR.mkdir(parents=True, exist_ok=True)
    build()
