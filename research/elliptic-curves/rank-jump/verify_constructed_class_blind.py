#!/usr/bin/env python3
"""Independent exact-fraction replay of blind degree-four cover recovery.

No Sage, Magma, elliptic-curve implementation, minimiser or point search is
used. The producer supplies proposed matrices and coordinates only.
"""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import permutations
import json
from math import gcd
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT = ROOT/'artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1'
PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def det(m):
    n = len(m)
    result = F(0)
    for p in permutations(range(n)):
        term = F((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
        for i in range(n):
            term *= m[i][p[i]]
        result += term
    return result


def mul(f, a, b):
    c = [F(0)]*5
    for i in range(3):
        for j in range(3):
            c[i+j] += a[i]*b[j]
    for k in [4, 3]:
        for j in range(3):
            c[k-3+j] -= c[k]*f[j]
    return c[:3]


def norm(f, a):
    columns = [mul(f, a, [F(i == j) for i in range(3)]) for j in range(3)]
    return det([[columns[j][i] for j in range(3)] for i in range(3)])


def quadrics(f, beta):
    result = [[F(0)]*10 for _ in range(3)]
    for k, (i, j) in enumerate(PAIRS):
        if j == 3:
            continue
        product = mul(f, beta, mul(f, [F(i == q) for q in range(3)], [F(j == q) for q in range(3)]))
        for l in range(3):
            result[l][k] = product[l]*(1 if i == j else 2)
    result[1][-1] += 1
    return result


def evaluate(q, p):
    return sum(c*p[i]*p[j] for c, (i, j) in zip(q, PAIRS))


def substitute(q, m):
    # Original column coordinates = m * reduced column coordinates.
    out = {pair: F(0) for pair in PAIRS}
    for c, (i, j) in zip(q, PAIRS):
        for k in range(4):
            for l in range(4):
                out[tuple(sorted((k, l)))] += c*m[i][k]*m[j][l]
    return [out[p] for p in PAIRS]


def verify_case(source, proposal, f, model):
    beta = list(map(F, source['beta_ascending']))
    q0, q1, q2 = quadrics(f, beta)
    ns = F(source['positive_norm_square_root'])
    assert ns > 0 and norm(f, beta) == ns*ns == F(source['norm'])
    m = [list(map(F, row)) for row in proposal['coordinate_matrix']]
    a = [list(map(F, row)) for row in proposal['equation_matrix']]
    assert len(m) == 4 and all(len(row) == 4 for row in m) and det(m)
    assert len(a) == 2 and all(len(row) == 2 for row in a) and det(a)
    reduced = [list(map(F, row)) for row in proposal['reduced_quadrics']]
    assert len(reduced) == 2 and all(len(row) == 10 for row in reduced)
    moved = [substitute(q, m) for q in [q2, q1]]
    assert [[sum(a[i][j]*moved[j][k] for j in range(2)) for k in range(10)] for i in range(2)] == reduced
    rp = list(map(F, proposal['reduced_point']))
    assert len(rp) == 4 and any(rp)
    assert all(evaluate(q, rp) == 0 for q in reduced)
    transported = [sum(m[i][j]*rp[j] for j in range(4)) for i in range(4)]
    cp = list(map(int, proposal['primitive_cover_point']))
    assert len(cp) == 4 and gcd(*cp) == 1 and cp[3]
    ratio = next(transported[i]/cp[i] for i in range(4) if cp[i])
    assert ratio and transported == [ratio*x for x in cp]
    p = list(map(F, cp))
    assert evaluate(q1, p) == 0 and evaluate(q2, p) == 0
    xi = [v/p[3] for v in p[:3]]
    X = evaluate(q0, p)/p[3]**2
    Y = ns*norm(f, p[:3])/p[3]**3
    assert mul(f, beta, mul(f, xi, xi)) == [X, F(-1), F(0)]
    assert Y*Y == sum(c*X**i for i, c in enumerate(f))
    x, y = X/4, (Y-X)/8
    a1, a2, a3, a4, a6 = model
    assert y*y+a1*x*y+a3*y == x**3+a2*x*x+a4*x+a6
    short = [36*x+3, 216*y+108*x]
    return {'column': source['column'], 'status': 'PASS_INDEPENDENT_FRACTION_REPLAY',
        'point_on_original_model': list(map(str, [x, y])),
        'point_on_cubic_model': list(map(str, [X, Y])),
        'point_on_short_model': list(map(str, short)),
        'primitive_cover_point': list(map(str, cp)),
        'exact_matrix_transport': True, 'exact_square_identity': True,
        'coefficientwise_quadric_identity': True}


def verify(witness_path, output):
    started = time.monotonic()
    source_path = EXPERIMENT/'covers.json'
    source = json.loads(source_path.read_text())
    proposed = json.loads(witness_path.read_text())
    assert proposed['input_sha256'] == sha(source_path)
    f = list(map(F, source['cubic_ascending']))
    assert f[3] == 1
    model = list(map(F, source['original_curve']))
    assert model[:3] == [1, 0, 0] and model[3] == f[1]/16 and model[4] == f[0]/64
    by_column = {c['column']: c for c in source['cases']}
    assert len(proposed['cases']) == 2 and {c['column'] for c in proposed['cases']} == set(by_column)
    rows = [verify_case(by_column[c['column']], c, f, model) for c in proposed['cases']]
    result = {'schema': 'constructed-class-blind-independent-replay.v1', 'status': 'PASS',
        'cases': rows, 'wall_seconds': time.monotonic()-started,
        'bindings': {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__).resolve(), source_path, witness_path]},
        'boundary': 'Independent rational replay proves both fixed covers have rational points with exact reduction and elliptic transport. Provenance of discovery is separately bound to the frozen worker transcript. No full rank or parameter criterion is asserted.'}
    with output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print('PASS exact reduction, cover points, field identities and elliptic transport', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--witness', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.witness.resolve(), args.output.resolve())
