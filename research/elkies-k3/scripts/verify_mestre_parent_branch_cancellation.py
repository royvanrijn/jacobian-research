#!/usr/bin/env python3
"""Portable replay of four generic Mestre branch gates; no Sage import.

Reconstructs short coefficients with b-invariant formulas, verifies section
and addition identities over Q[t], and rechecks finite witnesses with Rabin
tests, Euclidean gcds and Horner evaluation. Written valuation and surface
height proofs in the canonical note remain part of the theorem.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1'
SPEC = importlib.util.spec_from_file_location(
    'mestre_replay_arithmetic', ROOT/'elkies-k3/scripts/verify_r17_mestre_shared_twist.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
V, require = M.V, M.require
SOURCES = {
    'published-r17': 'artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/input.json',
    'alternate-q80': 'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json',
    'curve302-parent': 'artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json',
    'x1092-class1': 'artifacts/generated-results/elliptic-curves/x1092_class1_realization_compact_parent_v1.json',
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def source_coefficients(name, data):
    if name == 'published-r17':
        return V.poly(data['A']), V.poly(data['B'])
    if name == 'alternate-q80':
        w = data['weierstrass_model']
        return V.poly(w['A_coefficients_low_to_high']), V.poly(w['B_coefficients_low_to_high'])
    coefficients = []
    for a in data['a_invariants']:
        den = V.poly(a['denominator'])
        require(len(den) == 1, 'nonpolynomial source is outside this packet')
        coefficients.append(V.scale(V.poly(a['numerator']), 1/den[0]))
    a1, a2, a3, a4, a6 = coefficients
    b2 = V.add(V.power(a1, 2), V.scale(a2, 4))
    b4 = V.add(V.mul(a1, a3), V.scale(a4, 2))
    b6 = V.add(V.power(a3, 2), V.scale(a6, 4))
    A = V.add(V.scale(V.power(b2, 2), Q(-1, 48)), V.scale(b4, Q(1, 2)))
    B = V.add(V.sub(V.scale(V.power(b2, 3), Q(1, 864)),
                    V.scale(V.mul(b2, b4), Q(1, 24))), V.scale(b6, Q(1, 4)))
    return A, B


def verify_field(f, record, label):
    p = record['irreducible_prime']
    require(5 <= p < 2000, 'witness outside frozen prime range')
    require(M.irreducible(M.reduction(f, p), p), 'irreducible root field witness')
    place = record['local_place']
    ell, r = place['prime'], place['root']
    require(5 <= ell < 2000 and 0 <= r < ell, 'local place outside frozen range')
    require(ell % 4 == 3 and (label != 'A' or ell % 3 == 2), 'cyclotomic congruence')
    g = M.reduction(f, ell)
    dg = M.trim([i*c % ell for i, c in enumerate(g)][1:])
    require(M.evaluate(g, r, ell) == 0 and M.evaluate(dg, r, ell) != 0,
            'simple degree-one place')


def verify_parent(parent, row):
    name = parent['name']
    require(row['name'] == name and row['status'] == 'PASS', 'parent coverage or incomplete row')
    require(parent['source'] == SOURCES[name], 'source allowlist')
    path = ROOT/parent['source']
    require(parent['source_sha256'] == digest(path), 'retained source hash')
    A, B = source_coefficients(name, json.loads(path.read_text()))
    require((A, B) == (V.poly(parent['A']), V.poly(parent['B'])), 'short-model projection')
    core = V.add(V.scale(V.power(B, 2), 9261), V.scale(V.power(A, 3), 400))
    D = V.scale(V.mul(V.mul(A, B), core), -5)
    delta = V.add(V.scale(V.power(A, 3), 4), V.scale(V.power(B, 2), 27))
    degrees = [len(f)-1 for f in (A, B, core, D, delta)]
    require(degrees == row['degrees'] == [8, 12, 24, 44, 24], 'degree and infinity gate')
    p = row['smooth_disjoint_prime']
    require(p is not None and 5 <= p < 2000, 'missing smoothness witness')
    d, e = M.reduction(D, p), M.reduction(delta, p)
    for f in (d, e):
        df = M.trim([i*c % p for i, c in enumerate(f)][1:])
        require(len(M.ff_gcd(f, df, p)) == 1, 'squarefreeness')
    require(len(M.ff_gcd(d, e, p)) == 1, 'branch-discriminant disjointness')
    require(set(row['fields']) == {'A', 'B'}, 'coefficient divisor coverage')
    for label, f in [('A', A), ('B', B)]:
        verify_field(f, row['fields'][label], label)

    # The two ordinate factors multiply sqrt(D). Test after clearing every
    # rational denominator, independently of the producer's finite gate.
    points = [((V.scale(B, -21), V.scale(A, 5)), V.scale(V.power(A, 2), 25)),
              ((V.scale(B, -21), V.scale(A, 20)), V.scale(V.power(A, 2), 200))]
    for (xn, xd), yd in points:
        rhs = V.add(V.add(V.power(xn, 3), V.mul(A, V.mul(xn, V.power(xd, 2)))),
                    V.mul(B, V.power(xd, 3)))
        require(V.mul(D, V.power(xd, 3)) == V.mul(V.power(yd, 2), rhs),
                'Mestre section equation')
    # P2 has x2=x1/4 and y2=y1/8. Compute x(P1+P2) by the chord law.
    (xn, xd), yd = points[0]
    slope_num = V.scale(V.mul(D, V.power(xd, 2)), 49)
    slope_den = V.scale(V.mul(V.power(yd, 2), V.power(xn, 2)), 36)
    sum_num = V.sub(V.mul(slope_num, xd), V.scale(V.mul(xn, slope_den), Q(5, 4)))
    sum_den = V.mul(slope_den, xd)
    plus_num = V.add(V.scale(V.power(B, 2), Q(-7, 15)),
                     V.scale(V.power(A, 3), Q(-20, 81)))
    plus_den = V.mul(A, B)
    require(M.same_fraction((sum_num, sum_den), (plus_num, plus_den)), 'sum identity')
    require(len(M.ff_gcd(M.reduction(plus_num, p), M.reduction(plus_den, p), p)) == 1,
            'sum pole cancellation')
    require(len(plus_num)-len(plus_den) <= 4 and len(B)-len(A) <= 4, 'infinity poles')
    single_poles, sum_poles = len(A)-1, len(plus_den)-1
    single_height = 8 + 2*single_poles
    cross = Q(8+2*sum_poles-2*single_height, 2)
    require((single_height, cross) == (24, 0), 'height inputs')
    return {'name': name, 'status': 'PASS', 'fixed_u_cover_genus': 21,
            'height_matrix': [[single_height, int(cross)], [int(cross), single_height]],
            'rank_lower_bound_with_inherited_mw17': 19,
            'forced_branch_degree_for_every_rational_u': 20,
            'genus_lower_bound_for_every_rational_u': 9}


def verify(directory):
    packet = json.loads((directory/'input.json').read_text())
    result = json.loads((directory/'result.json').read_text())
    require(packet['schema'] == 'mestre-parent-branch-input-v1', 'input schema')
    require(result['schema'] == 'mestre-parent-branch-result-v1', 'result schema')
    require(result['input_sha256'] == digest(directory/'input.json'), 'input binding')
    require(packet['u'] == '2' and (packet['prime_min'], packet['prime_max_exclusive']) == (5, 2000),
            'frozen gate')
    names = list(SOURCES)
    require([p['name'] for p in packet['parents']] == names, 'input parent coverage')
    require([p['name'] for p in result['parents']] == names, 'result parent coverage')
    require(result['status'] == 'PASS' and result['positive_low_genus_cover'] is False, 'result scope')
    rows = [verify_parent(p, r) for p, r in zip(packet['parents'], result['parents'])]
    return {'status': 'PASS', 'parents': rows, 'positive_endpoint_complete': False,
            'parent_mw17_theorems': 'inherited; no rank or basis reconstruction',
            'proof_boundary': 'Exact arithmetic replay plus the written valuation and height proofs; no formal verification or exclusion of other identities.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    answer = verify(args.input_dir)
    answer['elapsed_seconds'] = round(time.monotonic()-start, 6)
    if args.output:
        with args.output.open('x') as stream:
            json.dump(answer, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(json.dumps(answer, sort_keys=True))
