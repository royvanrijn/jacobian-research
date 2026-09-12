#!/usr/bin/env python3
"""Independent rational replay of conic, binary-quartic and point transports."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
import verify_constructed_class_blind as exact


def add(a, b):
    return [(a[i] if i < len(a) else F(0))+(b[i] if i < len(b) else F(0)) for i in range(max(len(a), len(b)))]


def scale(a, c):
    return [c*x for x in a]


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def power(a, n):
    out = [F(1)]
    for _ in range(n):
        out = mul(out, a)
    return out


def substitute_quadric(q, vector):
    out = [F(0)]*5
    for c, (i, j) in zip(q, exact.PAIRS):
        out = add(out, scale(mul(vector[i], vector[j]), c))
    return out


def replay_case(witness, record):
    qs = [list(map(F, q)) for q in witness['reduced_quadrics']]
    assert qs[0][:4] == [0]*4 and qs[1][0] == 1
    pm = [list(map(F, row)) for row in record['conic_parametrization_matrix']]
    assert exact.det(pm)
    param = [[F(0)]]+[[pm[2][j], pm[1][j], pm[0][j]] for j in range(3)]
    assert substitute_quadric(qs[0], param) == [0]*5
    linear = [F(0)]*3
    for j in range(1, 4):
        linear = add(linear, scale(param[j], qs[1][j]))
    raw = add(mul(linear, linear), scale(substitute_quadric(qs[1], param), -4))
    assert raw == list(reversed(list(map(F, record['raw_quartic_descending']))))
    assert linear == list(reversed(list(map(F, record['linear_term_on_conic_param_descending']))))
    T = [list(map(F, row)) for row in record['degree2_binary_matrix']]
    k = F(record['degree2_scalar'])
    assert k and exact.det(T)
    a, b = [T[1][0], T[0][0]], [T[1][1], T[0][1]]
    moved = [F(0)]*5
    for i, c in enumerate(raw):
        moved = add(moved, scale(mul(power(a, i), power(b, 4-i)), c*k*k))
    reduced = list(reversed(list(map(F, record['reduced_quartic_descending']))))
    assert moved == reduced
    r, y, t = map(F, record['reduced_quartic_point'])
    assert t and y*y == sum(c*r**i*t**(4-i) for i, c in enumerate(reduced))
    assert max(abs((r/t).numerator), (r/t).denominator) <= record['search_bound_abscissa_height']
    a, b = r*T[0][0]+t*T[1][0], r*T[0][1]+t*T[1][1]
    conic = [a*a*pm[0][j]+a*b*pm[1][j]+b*b*pm[2][j] for j in range(3)]
    q0 = (y/k-sum(qs[1][j]*conic[j-1] for j in range(1, 4)))/2
    point = [q0]+conic
    assert point == list(map(F, witness['reduced_point'])) == list(map(F, record['reduced_degree4_point']))
    assert all(exact.evaluate(q, point) == 0 for q in qs)
    return {'column': witness['column'], 'status': 'PASS_FULL_QUARTIC_TRANSPORT',
            'quartic_parameter': str(r/t), 'bounded_height': record['search_bound_abscissa_height']}


def replay(witness_path, output):
    start = time.monotonic()
    witness = json.loads(witness_path.read_text())
    paths = [witness_path, Path(__file__).resolve(), Path(exact.__file__).resolve()]
    rows = []
    for case in witness['cases']:
        p = witness_path.parent/case['quartic_search_certificate']
        paths.append(p)
        rows.append(replay_case(case, json.loads(p.read_text())))
    result = {'status': 'PASS', 'cases': rows, 'wall_seconds': time.monotonic()-start,
        'bindings': {str(p.relative_to(exact.ROOT)): exact.sha(p) for p in paths},
        'boundary': 'Independent Fraction polynomial arithmetic checks the full conic parametrisation, degree-two transformation, square hit and map to the reduced degree-four cover. Original-cover transport is checked separately.'}
    with output.open('x') as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write('\n')
    print('PASS full conic-to-quartic-to-cover transport')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--witness', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    replay(a.witness.resolve(), a.output.resolve())
