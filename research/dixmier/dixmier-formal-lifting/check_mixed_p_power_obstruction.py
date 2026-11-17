#!/usr/bin/env python3
"""Exact finite controls for QPN1; universal coverage is the written proof."""
from math import comb
from pathlib import Path
import hashlib
import json
import tempfile

HERE = Path(__file__).resolve().parent


def mul(left, right, p, h=1):
    out = {}
    for (i, e), x in left.items():
        for (j, f), y in right.items():
            derivative = 1
            for k in range(i + 1):
                key = (i + j - k, e + f - 2 * k)
                value = x * y * comb(i, k) * derivative * pow(h, k, p)
                out[key] = (out.get(key, 0) + value) % p
                derivative = derivative * (f - 2 * k) % p
                if derivative == 0:
                    break
    return {key: value for key, value in out.items() if value}


def power(op, exponent, p, h=1):
    out = {(0, 0): 1}
    for _ in range(exponent):
        out = mul(out, op, p, h)
    return out


def admissible(op, n):
    if n < 2 or not op.get((n, 0)) or not op.get((0, 1)):
        return False
    for (j, e), value in op.items():
        if not value or (j, e) in [(n, 0), (0, 1)]:
            continue
        if j < 0 or e < 0:
            return False
        if e % 2 == 0:
            if j - e // 2 >= n:
                return False
        elif j >= e // 2:
            return False
    return True


cases = [
    ('n2_even_mixed', 2, 7, 1,
     {(2, 0): 1, (0, 1): 1, (2, 2): 2}),
    ('n2_odd_mixed_and_even', 2, 7, 1,
     {(2, 0): 1, (0, 1): 1, (1, 5): 3, (0, 3): 2, (1, 2): 4}),
    ('n2_generic_order_above_anchor', 2, 7, 1,
     {(2, 0): 1, (0, 1): 1, (4, 6): 2}),
    ('n2_nonmonic_h_and_lower_terms', 2, 11, 3,
     {(2, 0): 2, (0, 1): 4, (1, 0): 6, (0, 0): 5, (2, 2): 3}),
    ('n3_mixed_even_and_odd', 3, 11, 2,
     {(3, 0): 2, (0, 1): 3, (1, 5): 4, (3, 2): 2}),
]
controls = []
for label, n, p, h, op in cases:
    assert admissible(op, n), label
    assert p > 2 * n + 1
    target = ((n - 1) * p - (2 * n + 1), -2 * p)
    expected = (
        (n * n - 1) * pow(24 * n, -1, p)
        * pow(op[(n, 0)], p - 2, p) * op[(0, 1)]**2
        * pow(h, p + 1, p)
    ) % p
    result = power(op, p, p, h)
    actual = result.get(target, 0)
    assert actual == expected != 0, (label, actual, expected)
    assert target[0] > 0 and target[0] % p
    assert mul(result, {(0, 1): 1}, p, h) != mul({(0, 1): 1}, result, p, h)
    controls.append({
        'case': label, 'n': n, 'p': p, 'h': h,
        'S_terms': [[j, e, c] for (j, e), c in sorted(op.items())],
        'generic_a_order': max(j for j, _ in op),
        'target_a_order': target[0], 'target_v_exponent': target[1],
        'actual_coefficient': actual, 'predicted_coefficient': expected,
        'normal_order_terms': len(result), 'commutator_with_v_nonzero': True,
    })

p = 7
boundary = {(2, 0): 1, (0, 1): 1, (1, 3): 3}
assert not admissible(boundary, 2)
base_value = power({(2, 0): 1, (0, 1): 1}, p, p).get((2, -14), 0)
boundary_value = power(boundary, p, p).get((2, -14), 0)
assert base_value == 4 and boundary_value == 0

p = 11
P = {(2, 0): 1, (0, 2): pow(2, -1, p)}
Q = {(1, 0): -1}
left, right = mul(P, Q, p), mul(Q, P, p)
commutator = {key: (left.get(key, 0)-right.get(key, 0)) % p
              for key in set(left) | set(right)}
commutator = {key: val for key, val in commutator.items() if val}
assert commutator == {(0, 0): 1}
assert not admissible(P, 2)
assert all(j % p == 0 and e % p == 0 for j, e in power(P, p, p))

names = [
    'MIXED_NEWTON_P_POWER_OBSTRUCTION.md',
    'check_mixed_p_power_obstruction.py',
    'UNIVERSAL_P_POWER_OBSTRUCTION.md',
    'p_power_controls.json',
]
receipt = {
    'claim': 'QPN1', 'status': 'PASS',
    'full_sparse_Ore_power_controls': controls,
    'equal_deficit_boundary': {
        'p': 7, 'n': 2, 'added_term': '3*v^3*a',
        'support_gate_accepts': False,
        'base_target_coefficient': base_value,
        'boundary_target_coefficient': boundary_value,
        'conclusion': 'Witness preservation fails; no mate is claimed.',
    },
    'positive_pair': {
        'P': 'a^2+v^2/2', 'Q': '-a', 'p': 11,
        'commutator': '1', 'support_gate_accepts': False,
        'p_power_central': True,
    },
    'hashes': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
               for name in names},
    'coverage_boundary': (
        'Full finite sparse powers for declared mixed inputs; exact support '
        'gate, target coefficient and noncentrality; an equality-deficit '
        'cancellation and a positive canonical pair. Universal support '
        'selection, all finite orders, good reduction and arbitrary '
        'coefficient-cover descent use QPN1 and QPC1 written proofs. '
        'No universal Long normalization or scalarization is established.'
    ),
    'independent_replay': False, 'formal_verification': False,
    'external_review': False,
}
out = Path(tempfile.mkdtemp(prefix='qpn1-check-')) / 'mixed_p_power_controls.json'
out.write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'claim': 'QPN1', 'status': 'PASS',
                  'full_power_cases': len(controls), 'receipt': str(out)}, indent=2))
