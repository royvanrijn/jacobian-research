#!/usr/bin/env python3
"""Exact finite parameter-order controls for the QPN2 written theorem."""
from math import comb
from pathlib import Path
import hashlib
import json
import tempfile

HERE = Path(__file__).resolve().parent


def mul(left, right, p, ramification=1, eta_unit=1):
    """Normal key (a order, v exponent, t exponent), eta=eta_unit*t^e."""
    out = {}
    for (i, x, s), c in left.items():
        for (j, y, z), d in right.items():
            derivative = 1
            for k in range(i + 1):
                key = (i+j-k, x+y-2*k, s+z+ramification*k)
                value = c*d*comb(i, k)*derivative*pow(eta_unit, k, p)
                out[key] = (out.get(key, 0) + value) % p
                derivative = derivative*(y-2*k) % p
                if derivative == 0:
                    break
    return {key: value for key, value in out.items() if value}


def power(op, exponent, p, ramification=1, eta_unit=1):
    out = {(0, 0, 0): 1}
    for _ in range(exponent):
        out = mul(out, op, p, ramification, eta_unit)
    return out


def coefficient_series(op, a_order, v_exponent):
    return {t: value for (i, x, t), value in op.items()
            if i == a_order and x == v_exponent}


def support_gate(op, n, strict_even=False):
    if n < 2 or any(i < 0 or x < 0 or t < 0
                    for (i, x, t), value in op.items() if value):
        return False
    g = coefficient_series(op, n, 0)
    mu = coefficient_series(op, 0, 1)
    if not g or not mu.get(0) or (not strict_even and not g.get(0)):
        return False
    for (i, x, _), value in op.items():
        if not value or (i, x) in [(n, 0), (0, 1)]:
            continue
        if x % 2 == 0:
            if i-x//2 > n or (strict_even and i-x//2 == n):
                return False
        elif i > x//2:
            return False
    return True


cases = [
    ('closed_even_and_odd_boundary', 2, 7, 1, 1, False,
     {(2, 0, 0): 1, (0, 1, 0): 1, (3, 2, 0): 2, (1, 3, 0): 3}),
    ('Long_classical_and_small_boundary_corrections', 2, 7, 1, 1, False,
     {(2, 0, 0): -3, (0, 1, 0): 2, (1, 0, 0): 4,
      (3, 2, 1): 2, (1, 3, 1): 3, (2, 5, 1): 1}),
    ('ramified_closed_cone', 2, 7, 2, 3, False,
     {(2, 0, 0): 2, (0, 1, 0): 3, (4, 4, 0): 2, (2, 5, 0): 4}),
    ('strict_even_small_anchor', 2, 7, 1, 1, True,
     {(2, 0, 2): 3, (0, 1, 0): 2, (2, 2, 0): 3, (1, 3, 0): 4}),
    ('closed_n3_boundary', 3, 11, 1, 2, False,
     {(3, 0, 0): 2, (0, 1, 0): 3, (4, 2, 0): 1, (1, 3, 0): 2}),
]
controls = []
for label, n, p, ram, unit, strict, S in cases:
    assert p > 2*n+1
    assert support_gate(S, n, strict), label
    g = coefficient_series(S, n, 0)
    g_order = min(g)
    g_lead = g[g_order] % p
    mu0 = coefficient_series(S, 0, 1)[0] % p
    J = (n-1)*p-(2*n+1)
    target = coefficient_series(power(S, p, p, ram, unit), J, -2*p)
    expected_order = ram*(p+1) + (g_order*(p-2) if strict else 0)
    expected_lead = ((n*n-1)*pow(24*n, -1, p)*pow(g_lead, p-2, p)
                     *mu0**2*pow(unit, p+1, p)) % p
    assert target and min(target) == expected_order, (label, target)
    assert target[expected_order] == expected_lead != 0, label
    controls.append({
        'case': label, 'n': n, 'p': p, 'ramification': ram,
        'eta_unit': unit, 'strict_even_version': strict,
        'input_terms': [[i, x, t, c] for (i, x, t), c in sorted(S.items())],
        'generic_a_order': max(i for i, _, _ in S),
        'target_a_order': J, 'target_v_exponent': -2*p,
        'target_t_series': [[t, value] for t, value in sorted(target.items())],
        'expected_t_order': expected_order, 'actual_t_order': min(target),
        'expected_leading_coefficient': expected_lead,
        'actual_leading_coefficient': target[expected_order],
    })

# With a parameter pole, the old equality-boundary cancellation returns.
p = 7
pole = {(2, 0, 0): 1, (0, 1, 0): 1, (1, 3, -1): 2*pow(3, -1, p) % p}
assert not support_gate(pole, 2)
pole_target = coefficient_series(power(pole, p, p), 2, -14)
assert not pole_target, pole_target

# Positive first-order pair, including its t-dependent normal-order term.
p = 7
Slinear = {(1, 0, 0): 1, (0, 1, 0): 1}
Tlinear = {(2, 0, 0): -pow(2, -1, p), (1, 1, 0): -1,
           (0, -1, 1): -pow(2, -1, p)}
A, B = mul(Slinear, Tlinear, p), mul(Tlinear, Slinear, p)
comm = {key: (A.get(key, 0)-B.get(key, 0)) % p for key in set(A) | set(B)}
assert {key: value for key, value in comm.items() if value} == {(0, 0, 1): 1}
assert not support_gate(Slinear, 1)

# The elementary equality-boundary momentum shift is exact, not a limit.
c = 2
E = {(1, 0, 0): 1, (0, 3, 0): c*pow(2, -1, p) % p}
shifted = mul(E, E, p)
assert shifted.get((1, 3, 0)) == c
assert shifted.get((0, 1, 1)) == 3*c*pow(2, -1, p) % p
assert shifted.get((0, 6, 0)) == c*c*pow(4, -1, p) % p

names = [
    'CLOSED_NEWTON_CONE_OBSTRUCTION.md', 'check_closed_newton_obstruction.py',
    'closed_newton_source_audit.json', 'UNIVERSAL_P_POWER_OBSTRUCTION.md',
    'p_power_controls.json', 'MIXED_NEWTON_P_POWER_OBSTRUCTION.md',
    'mixed_p_power_controls.json',
]
receipt = {
    'claim': 'QPN2', 'status': 'PASS', 'parameter_power_controls': controls,
    'parameter_pole_boundary': {'p': 7, 'coefficient': '2/(3*t)',
        'support_gate_accepts': False, 'target_coefficient_zero': True,
        'boundary': 'This is cancellation of one witness; no mate is asserted.'},
    'positive_n1_pair': {'p': 7, 'commutator': 't', 'support_gate_accepts': False},
    'exact_momentum_shift_identity': True,
    'hashes': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
               for name in names},
    'coverage_boundary': (
        'Finite sparse parameter-polynomial powers for declared closed-cone '
        'inputs, strict-even small anchor, a ramified parameter, parameter-pole '
        'cancellation, a positive n1 pair and an exact momentum shift. '
        'Universal deficit selection, t-order separation, Artin approximation, '
        'spreading and lifting in good characteristic use written proofs and '
        'the inspected source audit; no finite test proves those steps.'
    ),
    'independent_replay': False, 'formal_verification': False,
    'external_review': False,
}
output = Path(tempfile.mkdtemp(prefix='qpn2-check-')) / 'closed_newton_controls.json'
output.write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'claim': 'QPN2', 'status': 'PASS',
                  'power_cases': len(controls), 'receipt': str(output)}, indent=2))
