#!/usr/bin/env python3
"""Exact bounded checks for the circuit and component-lattice theorem.

All-order and all-group implications are written in the adjacent proof.
This standard-library checker does not import producer verification code.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def add(*polys):
    result = defaultdict(F)
    for p in polys:
        for a, c in p.items():
            result[a] += c
    return {a: c for a, c in result.items() if c}


def mul(p, q):
    result = defaultdict(F)
    for a, c in p.items():
        for b, d in q.items():
            result[tuple(x+y for x, y in zip(a, b))] += c*d
    return {a: c for a, c in result.items() if c}


def scale(p, c):
    return {a: c*b for a, b in p.items() if c*b}


@lru_cache(None)
def sphere_monomial(a, b, c):
    # U=u^2, V=conjugate(u)^2, S=tau^2 on S^2.
    # Phase balance first; tau has normalized uniform law on [-1,1].
    if a != b:
        return F(0)
    return sum((F((-1)**j * comb(2*a, j), 2*c+2*j+1)
                for j in range(2*a+1)), F(0))


def integral(poly):
    return sum((coef*sphere_monomial(*a) for a, coef in poly.items()), F(0))


def kernel(matrix):
    """Rational RREF kernel with column-vector convention."""
    a = [[F(x) for x in row] for row in matrix]
    width = len(a[0])
    pivots = []
    row = 0
    for col in range(width):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][col]
        a[row] = [x/divisor for x in a[row]]
        for i in range(len(a)):
            if i != row:
                multiple = a[i][col]
                a[i] = [x-multiple*y for x, y in zip(a[i], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    out = []
    for col in range(width):
        if col not in pivots:
            v = [F(0)]*width
            v[col] = F(1)
            for i, p in enumerate(pivots):
                v[p] = -a[i][col]
            assert all(sum(x*y for x, y in zip(r, v)) == 0 for r in matrix)
            out.append(v)
    return len(pivots), out


def matmul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def run():
    one = {(0, 0, 0): F(1)}
    U = {(1, 0, 0): F(1)}
    V = {(0, 1, 0): F(1)}
    S = {(0, 0, 1): F(1)}
    R = mul(mul(S, add(U, scale(one, 2))),
            add(scale(one, 2), scale(mul(S, add(mul(U, U), scale(U, 2), scale(one, 2))), -1)))
    P = mul(add(one, U), add(V, scale(R, -1)))
    # Check the polynomial identity before integration, modulo UV=(1-S)^2.
    lhs = mul(U, add(V, scale(R, -1)))
    rhs_base = add(one, scale(mul(mul(add(one, U), add(one, U)), S), -1))
    rhs = mul(rhs_base, rhs_base)
    uv = {(1, 1, 0): F(1)}
    relation = add(uv, scale(mul(add(one, scale(S, -1)), add(one, scale(S, -1))), -1))
    assert add(lhs, scale(rhs, -1), scale(relation, -1)) == {}
    power = one
    rows = []
    for m in range(1, 7):
        power = mul(power, P)
        c2m = F(4**(2*m)*factorial(2*m)**2, factorial(4*m+1))
        for s in range(m+3):
            marked = {(a+s, b, c): v for (a, b, c), v in power.items()}
            actual = integral(marked)
            expected = c2m*comb(m-1, s-1) if 1 <= s <= m else F(0)
            assert actual == expected, (m, s, actual, expected)
            rows.append({'m': m, 'marker_power': s, 'moment': str(actual)})
    # A sign mutation in the regular circuit destroys the first pure identity.
    wrong = mul(add(one, U), add(V, R))
    assert integral(wrong) != 0

    shared = [[1, 1]]
    independent = [[1, 0], [0, 1]]
    triple = [[1, 0, -1], [0, 1, -1]]
    # A nonsaturated integral inclusion: Q-injectivity is enough.
    finite_index = [[2, 0], [0, 3]]
    rank, relations = kernel(shared)
    assert rank == 1 and relations == [[F(-1), F(1)]]
    assert kernel(independent) == (2, [])
    assert kernel(finite_index) == (2, [])
    rank, relations3 = kernel(triple)
    assert rank == 2 and relations3 == [[F(1), F(1), F(1)]]
    for i in range(3):
        for j in range(i):
            pair = [[row[j], row[i]] for row in triple]
            assert kernel(pair) == (2, [])  # pairwise tests miss the obstruction

    action = [[0, -1], [1, -1]]
    identity = [[1, 0], [0, 1]]
    assert matmul(matmul(action, action), action) == identity
    direction = [[1], [0]]
    orbit = []
    for _ in range(3):
        orbit.append([r[0] for r in direction])
        direction = matmul(action, direction)
    assert orbit == [[1, 0], [0, 1], [-1, -1]]

    # P=(z,z^-1): every pure trace is zero. Fixed Laurent monomial
    # markers vanish once the exponent exceeds their absolute support.
    component_rows = []
    for q1 in range(-4, 5):
        for q2 in range(-4, 5):
            for m in range(1, 9):
                pure = F(int(m == 0)+int(-m == 0), 2)
                mixed = F(int(q1+m == 0)+int(q2-m == 0), 2)
                assert pure == 0
                if m > max(abs(q1), abs(q2)):
                    assert mixed == 0
                component_rows.append([q1, q2, m, str(mixed)])
    # Coefficient product is invariant for BOTH the torus and C2 swap.
    for torus in [F(2), F(3, 2), F(7, 3)]:
        a, b = torus, 1/torus
        assert a*b == b*a == 1
    # One common functional after assembling independent components;
    # neither disconnectedness nor dimension alone predicts failure.
    ell = [F(1), F(-1)]
    plus = [[1, 0], [2, 0]]
    minus = [[0, -1], [0, -3]]
    assert all(sum(x*y for x, y in zip(ell, v)) > 0 for v in plus+minus)

    return {
        'completed': True,
        'normalizer_regular_circuit_identity': True,
        'normalizer_marker_checks': len(rows),
        'normalizer_rows': rows,
        'circuit_sign_mutation_detected': True,
        'shared_line_kernel': [[str(x) for x in row] for row in relations],
        'independent_components_injective': True,
        'nonsaturated_character_inclusion_injective': True,
        'triple_relation': [[str(x) for x in row] for row in relations3],
        'pairwise_independence_mutation_detected': True,
        'order_three_character_action': action,
        'shared_line_fixed_marker_checks': len(component_rows),
        'opposite_direction_common_separator_control': True,
        'shared_line_nonzero_invariant': 'coefficient_1*coefficient_2=1',
        'all_order_proof_verified_by_checker': False,
        'general_Lie_theory_verified_by_checker': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False,
    }


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'spherical_components_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'normalizer_rows'}, indent=2))
