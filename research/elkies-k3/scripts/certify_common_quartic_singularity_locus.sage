#!/usr/bin/env sage -python
"""Exact low-genus controls and a constant-abscissa discriminant gate.

This does not search the fixed parents' full quartic-abscissa systems.
The prime was selected after the small discovery probes retained in execution.json.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import GF, PolynomialRing, QQ, EllipticCurve, gcd
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1'
SOURCE = 'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
NAMES = ['published-r17', 'alternate-q80', 'curve302-parent', 'x1092-class1']
PRIME = 1009


def write_new(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write('\n')


def coeffs(f):
    return [str(c) for c in f.list()]


def freeze(out):
    source = ROOT / SOURCE
    packet = json.loads(source.read_text())
    assert [p['name'] for p in packet['parents']] == NAMES
    helpers = ['elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py',
               'elkies-k3/scripts/verify_r17_mestre_shared_twist.py']
    write_new(out / 'input.json', {
        'schema': 'common-quartic-singularities-input-v1',
        'source': SOURCE, 'source_sha256': sha256(source.read_bytes()).hexdigest(),
        'parents': NAMES, 'prime': PRIME,
        'producer_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic_helpers': {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in helpers},
        'controls': [
            {'name': 'four-nodes-genus-one', 'D': [1, 1, 0, 0, 1],
             'square_factor': [1], 'd': [1, 1, 0, 0, 1], 's': [1, 0, 0, 0, 1]},
            {'name': 'five-nodes-genus-zero', 'D': [0, 0, 1, 0, 1],
             'square_factor': [0, 1], 'd': [1, 0, 1], 's': [1, 0, 0, 0, 1]},
        ],
        'limits': {'cpu_seconds': 40, 'address_space_gib': 4},
        'scope': 'Symbolic identities, two explicitly constructed control parents, and four complete constant-x parameter lines. No fixed-parent quartic search, point search, descent, or Mestre auxiliary-function variation.',
        'selection': 'Controls are derived algebraically. Prime1009 was chosen after the documented modular probes; it is a proof witness, not a predeclared search sample.',
        'sage_version': version,
    })
    print('Frozen four constant-abscissa gates and two new-parent controls', flush=True)


def constant_gate(parent, p):
    R = PolynomialRing(GF(p), 't')
    A, B = [R([GF(p)(QQ(v)) for v in parent[k]]) for k in ('A', 'B')]
    assert (A.degree(), B.degree()) == (8, 12)
    Ap, Bp = A.derivative(), B.derivative()
    H = B * Ap**3 - A * Bp * Ap**2 - Bp**3
    assert H.degree() == 33 and gcd(H, Ap) == 1
    K = PolynomialRing(GF(p), 'c'); c = K.gen()
    S = PolynomialRing(K, 't')
    N = K(S(H).resultant(c * S(Ap) + S(Bp))).monic()
    assert N.degree() == 33 and gcd(N, N.derivative()) == 1
    assert gcd(H, H.derivative()) == 1
    return {'name': parent['name'], 'prime': p,
            'critical_t_polynomial_mod_p': [int(x) for x in H.list()],
            'critical_value_norm_monic_mod_p': [int(x) for x in N.list()],
            'constant_x_genus_lower_bound': 4,
            'critical_values_over_Qbar': 33,
            'scope': 'Every algebraic constant x=c, with no bound on c; not nonconstant quartic x(t).'}


def control(row, p):
    R = PolynomialRing(QQ, 't'); t = R.gen()
    D, e, d, s = [R(row[k]) for k in ('D', 'square_factor', 'd', 's')]
    assert D == e**2 * d
    A, B = D * (2*s + 1) - 1, D * s**2
    delta = 4*A**3 + 27*B**2
    assert (A.degree(), B.degree(), delta.degree()) == (8, 12, 24)
    assert all(gcd(f, f.derivative()) == 1 for f in (delta, d, s, s+1))
    assert all(gcd(D, f) == 1 for f in (delta, s, s+1))
    pts = [(R(0), e*s), (R(1), e*(s+1)), (D-1, e*(1-D-s))]
    for x, r in pts:
        assert x**3 + A*x + B == d*r**2
        assert x.degree() <= 4 and 2*r.degree() + d.degree() <= 12
    R0 = PolynomialRing(GF(p), 't'); disc0 = R0(delta)
    assert disc0.degree() == 24 and gcd(disc0, disc0.derivative()) == 1
    assert gcd(disc0, R0(D)) == 1
    g = (int(d.degree())-2)//2
    ans = {'name': row['name'], 'A': coeffs(A), 'B': coeffs(B),
           'D': coeffs(D), 'square_factor': coeffs(e), 'd': coeffs(d), 's': coeffs(s),
           'points_and_sum': [{'x': coeffs(x), 'y_over_w': coeffs(r)} for x, r in pts],
           'branch_degree': int(d.degree()), 'normalization_genus': g,
           'image_arithmetic_genus': 5, 'nodes_per_image': 5-g,
           'height_matrix': [[8, -4], [-4, 8]], 'height_determinant': 48,
           'new_directions_modulo_full_inherited_group': 2,
           'smoothness_prime': p, 'inherited_parent_rank': 'UNKNOWN',
           'is_a_retained_mw17_parent': False}
    if g == 1:
        E = EllipticCurve(QQ, [-4, 1]); P = E(0, 1)
        ans['rational_base_certificate'] = {
            'type': 'elliptic', 'a_invariants': [0, 0, 0, -4, 1],
            'forward': ['X=2*(w+t^2)', 'Y=2*t*X+1'],
            'inverse': ['t=(Y-1)/(2*X)', 'w=X/2-t^2'],
            'point': ['0', '1'],
            'twice_point': [str(c) for c in list(2*P)[:2]],
            'three_times_point': [str(c) for c in list(3*P)[:2]],
            'nontorsion_proof': '3P has a nonintegral x-coordinate on the integral short model; Lutz-Nagell.'}
    else:
        ans['rational_base_certificate'] = {
            'type': 'rational', 't': '2*z/(1-z^2)', 'w': '(1+z^2)/(1-z^2)'}
    return ans


def run(out):
    packet = json.loads((out / 'input.json').read_text())
    assert sha256((ROOT / packet['source']).read_bytes()).hexdigest() == packet['source_sha256']
    resource.setrlimit(resource.RLIMIT_CPU, (40, 45))
    resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
    start = time.process_time()
    parents = json.loads((ROOT / packet['source']).read_text())['parents']
    gates, controls = [], []
    for parent in parents:
        gate = constant_gate(parent, packet['prime'])
        write_new(out / (parent['name'] + '-checkpoint.json'), gate)
        gates.append(gate)
        print(parent['name'] + ': all constant x have normalization genus >=4', flush=True)
    for row in packet['controls']:
        result = control(row, packet['prime'])
        write_new(out / (row['name'] + '-checkpoint.json'), result)
        controls.append(result)
        print(row['name'] + ': two independent gains on an infinite rational base', flush=True)
    write_new(out / 'result.json', {
        'schema': 'common-quartic-singularities-result-v1',
        'input_sha256': sha256((out / 'input.json').read_bytes()).hexdigest(),
        'parents': gates, 'controls': controls,
        'fixed_parent_full_quartic_chart': 'UNKNOWN', 'mw17_endpoint_complete': False,
        'component_cpu_seconds': round(time.process_time()-start, 6)})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('freeze', 'run'))
    p.add_argument('--output', type=Path, default=DEFAULT)
    args = p.parse_args()
    (freeze if args.mode == 'freeze' else run)(args.output)


if __name__ == '__main__':
    main()
