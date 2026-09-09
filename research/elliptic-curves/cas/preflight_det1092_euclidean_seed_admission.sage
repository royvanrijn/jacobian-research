#!/usr/bin/env sage-python
"""Cold generic-word -> conic -> seed admission on three fixed old addresses.

Each action is checkpointed and limited to 25 seconds. No saved conic/seed,
later point, search output or catalogue is a construction input. Addresses
are retrospective regression choices, NOT a prospective parameter selector.
"""
import argparse
import csv
import hashlib
import json
import signal
import time
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector
from sage.env import SAGE_VERSION
from split_seed_descent import build_frame, Classifier, point_record

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_euclidean_seed_admission_v1'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
TABLE = ART / 'curve302_parent_degree2_multisection_orbits_v1.tsv'
CHART = ART / 'det1092_reduced_parameter_chart_v1/generic-proof.json'
ROSTER = [dict(label='302', original_parameter='0'),
          dict(label='old-positive', reduced_parameter='1926/2699'),
          dict(label='old-dependent', reduced_parameter='-528/3635')]
R = PolynomialRing(QQ, 't')
K = R.fraction_field()


def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def coefficients(f): return list(map(str, f.list()))
def decode(v): return K(R(v['numerator'])) / R(v['denominator'])
def save(name, data):
    p = OUT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        assert read(p) == data, 'immutable checkpoint changed: ' + str(p)
    else:
        with p.open('x') as stream:
            json.dump(data, stream, indent=2, sort_keys=True)
            stream.write('\n')


def freeze():
    sources = [PARENT, TABLE, CHART, Path(__file__),
               Path(__file__).with_name('split_seed_descent.py'),
               Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage')]
    save('protocol.json', dict(
        classification='verified application preflight, not a new theorem or rank result',
        mask=8044, roster=ROSTER, sage_version=SAGE_VERSION,
        limits=dict(seconds_per_action=25, orbits=1, old_addresses=3,
                    prime_cap=1009, halving_steps=8, point_searches=0,
                    new_parameters=0, full_atlas_runs=0),
        rule='Rebuild the negative trace from generic sections and the one fixed table word. '
             'Construct its Euclidean conic without saved trace, conic or point inputs. '
             'At each old address, first seal a generic-only finite footprint; then test q. '
             'For a nonzero rational square use its nonnegative root and classify that branch. '
             'No substitution after a miss or a cap.',
        boundary='Both orbit and evaluation addresses are fixed historical regressions. '
                 'The constructor/admission procedure is equation-only, but this roster '
                 'does not demonstrate prospective parameter selection or a302 seed.',
        inputs={str(p.relative_to(ROOT)): sha(p) for p in sources}))


def protocol():
    p = read(OUT / 'protocol.json')
    for name, digest in p['inputs'].items():
        assert sha(ROOT / name) == digest, 'input changed: ' + name
    return p


def generic():
    parent = read(PARENT)
    original = EllipticCurve(K, [decode(a) for a in parent['a_invariants']])
    A, B = R(-original.c4()/48), R(-original.c6()/864)
    E = EllipticCurve(K, [A, B])
    basis = []
    for row in parent['basis_weierstrass_coordinates']:
        x, y = map(decode, row)
        basis.append(E([x+original.b2()/12, y+(original.a1()*x+original.a3())/2]))
    return parent, original, E, A, B, basis


def construct():
    p = protocol()
    parent, original, E, A, B, basis = generic()
    with TABLE.open() as stream:
        matches = [row for row in csv.DictReader(stream, delimiter='\t')
                   if int(row['orbit_mask']) == p['mask']]
    assert len(matches) == 1 and matches[0]['category'] == 'rational'
    w = vector(ZZ, matches[0]['parent_MW17_w'].split())
    assert w*matrix(ZZ, parent['generic_height_gram'])*w == 10
    trace = -sum((n*P for n, P in zip(w, basis)), E(0))
    h = R(trace[0].denominator()).sqrt().monic()
    N, V = R(trace[0]*h*h), R(trace[1]*h**3)
    assert h.degree() == 3 and N.gcd(h) == 1
    assert V*V == N**3 + A*N*h**4 + B*h**6
    def divide(a, b):
        quotient, remainder = a.quo_rem(b)
        assert not remainder
        return quotient
    m = (-V*N.inverse_mod(h*h)) % (h*h)
    g = divide(m*N+V, h*h)
    b = divide(m*m-N, h*h)
    k = divide(m*b-2*g, h*h)
    q = divide(4*m*k-3*b*b-4*A, h*h)
    assert 1 <= q.degree() <= 2 and q.gcd(q.derivative()) == 1
    x0, x1, y0, y1 = b/2, h/2, -h*k/2, -m/2
    assert y0*y0+y1*y1*q == x0**3+3*x0*x1*x1*q+A*x0+B
    assert 2*y0*y1 == 3*x0*x0*x1+x1**3*q+A*x1
    save('construction.json', dict(status='PASS_COLD_GENERIC_EUCLIDEAN_CONSTRUCTION',
        mask=p['mask'], word=list(map(int, w)),
        h=coefficients(h), nx=coefficients(N), ny=coefficients(V),
        m=coefficients(m), g=coefficients(g), b=coefficients(b), k=coefficients(k),
        q=coefficients(q), maps=[coefficients(f) for f in [x0, x1, y0, y1]],
        convention='W^2=q(t); X=x0+x1*W,Y=y0+y1*W in short original parent',
        protocol_sha256=sha(OUT/'protocol.json')))
    print('PASS_COLD_GENERIC_EUCLIDEAN_CONSTRUCTION', flush=True)


def run(label):
    p = protocol()
    entry = next(row for row in p['roster'] if row['label'] == label)
    if 'original_parameter' in entry:
        tau = QQ(entry['original_parameter'])
    else:
        s = QQ(entry['reduced_parameter'])
        aa, bb, cc, dd = map(QQ, read(CHART)['parameter_matrix'])
        assert aa*dd-bb*cc and cc*s+dd
        tau = (aa*s+bb)/(cc*s+dd)
    _, original, _, A, B, basis = generic()
    E = EllipticCurve(QQ, [A(tau), B(tau)])
    assert E.discriminant()
    specialized = [E([v(tau) for v in P.xy()]) for P in basis]
    frame = build_frame(E, specialized, p['limits']['prime_cap'])
    save(label+'/generic-frame.json', frame)
    # No candidate coordinates or square evaluation before this checkpoint.
    c = read(OUT/'construction.json')
    q = R(c['q'])(tau)
    row = dict(label=label, original_parameter=str(tau), q_value=str(q),
        protocol_sha256=sha(OUT/'protocol.json'),
        construction_sha256=sha(OUT/'construction.json'),
        frame_sha256=sha(OUT/label/'generic-frame.json'))
    if q < 0:
        row.update(status='EXACT_NONSPLIT', square_certificate=dict(negative=True))
    else:
        n, d = ZZ(q.numerator()), ZZ(q.denominator())
        a, b = n.isqrt(), d.isqrt()
        row['square_certificate'] = dict(numerator=str(n), denominator=str(d),
            numerator_floor_sqrt=str(a), denominator_floor_sqrt=str(b))
        if a*a != n or b*b != d:
            row['status'] = 'EXACT_NONSPLIT'
        else:
            W = QQ(a)/b
            x0, x1, y0, y1 = [R(v)(tau) for v in c['maps']]
            P, companion = E([x0+x1*W, y0+y1*W]), E([x0-x1*W, y0-y1*W])
            trace = sum((n*Q for n, Q in zip(c['word'], specialized)), E(0))
            assert P+companion == trace
            row.update(root=str(W), point=point_record(P), companion=point_record(companion))
            if not W:
                assert 2*P == trace
                row['status'] = 'RAMIFIED_INHERITED'
            else:
                result = Classifier(frame).classify(P, p['limits']['halving_steps'])
                row.update(status=result['status'], admission=result)
    save(label+'/result.json', row)
    print(label, row['status'], row.get('admission', {}).get('reason'), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze', 'construct', 'run'])
    parser.add_argument('--case', choices=[r['label'] for r in ROSTER])
    args = parser.parse_args()
    signal.alarm(25)
    began = time.monotonic()
    if args.action == 'freeze': freeze()
    elif args.action == 'construct': construct()
    else:
        assert args.case is not None
        run(args.case)
    print('seconds', round(time.monotonic()-began, 3), flush=True)
