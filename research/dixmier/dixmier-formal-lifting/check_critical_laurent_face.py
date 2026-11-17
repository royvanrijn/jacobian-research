#!/usr/bin/env python3
"""Exact finite controls for QCL1; universal coverage uses its written proof."""
from pathlib import Path
from math import comb
import hashlib
import json
import tempfile

import sympy as sp

HERE = Path(__file__).resolve().parent
V, s = sp.symbols('V s', nonzero=True)


def clean(op):
    return {j: sp.cancel(c) for j, c in op.items() if sp.cancel(c) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for j, c in op.items():
            out[j] = out.get(j, 0) + c
    return clean(out)


def scale(op, scalar):
    return clean({j: scalar*c for j, c in op.items()})


def delta(f):
    return sp.diff(f, V)/V


def mul(left, right, eta=1):
    out = {}
    for i, x in left.items():
        for j, y in right.items():
            derivative = y
            for k in range(i+1):
                key = i+j-k
                out[key] = out.get(key, 0) + comb(i, k)*x*derivative
                derivative = eta*delta(derivative)
    return clean(out)


def power(op, n, eta=1):
    out = {0: sp.Integer(1)}
    for _ in range(n):
        out = mul(out, op, eta)
    return out


def comm(left, right, eta=1):
    return add(mul(left, right, eta), scale(mul(right, left, eta), -1))


def face_gate(terms, e=1):
    """Input (a order, v exponent, parameter order, nonzero scalar)."""
    if any(j < 0 or r < 0 for j, _, r, c in terms if c):
        return False
    if not any(j == 2 and m == 0 and r == 0 and c for j,m,r,c in terms):
        return False
    if not any(j == 0 and m == 1 and r == 0 and c for j,m,r,c in terms):
        return False
    values = [(5*r+e*(j+2*m-2), j) for j,m,r,c in terms if c]
    return all(F >= 0 for F, _ in values) and max(j for F,j in values if F == 0) == 2


def blowup(terms, e=1):
    out = {}
    for j, m, r, c in terms:
        F = 5*r+e*(j+2*m-2)
        out[j] = out.get(j, 0) + c*s**F*V**m
    return clean(out)


# A quadratic Laurent face, with both linear and quadratic negative tails.
terms = [
    (2, 0, 0, 1), (0, 1, 0, 1),
    (2, -5, 2, 2), (2, -10, 4, 1),
    (1, -2, 1, 1), (0, -4, 2, sp.Rational(1, 6)),
    (7, 0, 2, 3),
]
assert face_gate(terms)
rescaled = blowup(terms)
P = {j: sp.cancel(c.subs(s, 0)) for j,c in rescaled.items()
     if sp.cancel(c.subs(s, 0)) != 0}
assert max(rescaled) == 7 and max(P) == 2
assert rescaled[7] == 3*s**15
assert P == {2: (V**10+2*V**5+1)/V**10,
             1: 1/V**2, 0: (6*V**5+1)/(6*V**4)}

# First inverse-square equality survives, while an order-four equality does not.
assert face_gate([(2, 0, 0, 1), (0, 1, 0, 1), (0, -4, 2, 7)])
assert not face_gate([(2, 0, 0, 1), (0, 1, 0, 1), (4, -6, 2, 1)])
assert not face_gate([(2, 0, 0, 1), (0, 1, 0, 1), (0, -5, 2, 1)])
assert face_gate([(2, 0, 0, 2), (0, 1, 0, 3), (0, -4, 4, 5)], e=2)

# Polynomial mixed shapes far outside QPN2 are all above the critical face.
arbitrary_polynomial = [
    (2, 0, 0, 1), (0, 1, 0, 1),
    (10, 1, 1, 2), (20, 2, 2, 3), (100, 0, 1, 4),
]
assert face_gate(arbitrary_polynomial)
polynomial_scaled = blowup(arbitrary_polynomial)
assert max(polynomial_scaled) == 100
assert clean({j: c.subs(s, 0) for j,c in polynomial_scaled.items()}) == {2: 1, 0: V}

# General shape of the exact quadratic normalization, in a rational root example.
rho = 1+V**-5
H, J = P[1], P[0]
ell = sp.cancel(H/rho-delta(rho))
E1 = {1: rho, 0: ell/2}
eps = lambda f: sp.cancel(rho*delta(f))
U = sp.cancel(J-ell**2/4-eps(ell)/2)
assert add(power(E1, 2), {0: U}) == P
epsU = eps(U)
assert sp.limit(U/V, V, sp.oo) == 1
assert sp.limit(epsU, V, sp.oo) == 0
assert sp.limit(V*epsU, V, sp.oo) == 1

# Declared finite Lenard controls; the all-order dominance proof is written.
derivatives = [U]
for _ in range(4):
    derivatives.append(eps(derivatives[-1]))
R1 = U/4
R2 = (derivatives[2]+3*U**2)/16
R3 = (derivatives[4]+10*U*derivatives[2]+5*derivatives[1]**2+10*U**3)/64
lenard_leads = []
for j, Rj in enumerate([R1, R2, R3], 1):
    actual = sp.limit(Rj/V**j, V, sp.oo)
    beta = sp.Rational(comb(2*j, j), 2*4**j)
    assert actual == beta != 0
    lenard_leads.append({'j': j, 'leading_coefficient': str(actual)})

# Pole removal in a positive Airy pair. Its potential grows like V^2,
# so it correctly lies outside QCL1's linear-growth infinity hypothesis.
Airy = {2: sp.Integer(1), 0: V**2/2}
Tregular = {1: -sp.Integer(1)}
assert comm(Airy, Tregular) == {0: sp.Integer(1)}
Tpolar = add(Tregular, scale(power(Airy, 2), s**-3), scale(Airy, s**-1))
assert comm(Airy, Tpolar) == {0: sp.Integer(1)}
leading3 = clean({j: sp.limit(s**3*c, s, 0) for j,c in Tpolar.items()})
assert leading3 == power(Airy, 2)
after3 = add(Tpolar, scale(power(Airy, 2), -s**-3))
leading1 = clean({j: sp.limit(s*c, s, 0) for j,c in after3.items()})
assert leading1 == Airy
after1 = add(after3, scale(Airy, -s**-1))
assert after1 == Tregular

# Positive order-one control: the quadratic-face requirement is substantive.
Plinear = {1: sp.Integer(1), 0: V}
Tlinear = {2: -sp.Rational(1, 2), 1: -V, 0: -1/(2*V)}
assert comm(Plinear, Tlinear) == {0: sp.Integer(1)}

names = [
    'CRITICAL_LAURENT_FACE_OBSTRUCTION.md', 'check_critical_laurent_face.py',
    'QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md',
    'EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md',
    'UNIVERSAL_P_POWER_OBSTRUCTION.md',
]
receipt = {
    'claim': 'QCL1', 'status': 'PASS',
    'rescaling': {
        't': 's^5', 'v': 's^(2e)*V', 'a': 's^e*A',
        'S_scale': 's^(-2e)', 'T_scale': 's^(-3e)',
        'input_terms': [[j,m,r,str(c)] for j,m,r,c in terms],
        'generic_a_order': 7, 'critical_a_order': 2,
        'higher_order_correction_after_rescaling': '3*s^15*A^7',
        'ramified_e2_control': True,
    },
    'support_boundaries': {
        'inverse_square_t2_vminus4_accepted': True,
        'critical_t2_vminus6_a4_rejected': True,
        'negative_weight_t2_vminus5_rejected': True,
        'unrestricted_polynomial_mixed_shapes_accepted': True,
        'polynomial_generic_order100_critical_order2': True,
    },
    'quadratic_normalization': {
        'rho': str(rho), 'P_equals_E1_squared_plus_U': True,
        'U_over_V_limit': 1, 'epsilon_U_limit': 0,
        'V_times_epsilon_U_limit': 1,
    },
    'finite_Lenard_leading_controls': lenard_leads,
    'finite_pole_removal': {
        'positive_control': 'P=A^2+V^2/2, T=-A',
        'initial_poles': [3, 1], 'subtractions': ['s^-3*P^2', 's^-1*P'],
        'final_mate': '-A', 'commutator': '1',
        'boundary': 'The positive potential has quadratic V growth, outside QCL1.',
    },
    'positive_first_order_pair': {'P': 'A+V', 'commutator': '1'},
    'hashes': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
               for name in names},
    'coverage_boundary': (
        'Exact finite rescaling, support boundaries, quadratic normalization, '
        'infinity limits, Lenard terms1,2,3, and positive canonical pairs/pole '
        'subtractions. All-order centralizer exhaustion and quadratic mate '
        'exclusion depend on QSO2; finite termination of pole descent and '
        'critical-root normalization are written proofs, not finite scans.'
    ),
    'independent_replay': False, 'formal_verification': False,
    'external_review': False,
}
output = Path(tempfile.mkdtemp(prefix='qcl1-check-')) / 'critical_laurent_face_controls.json'
output.write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'claim': 'QCL1', 'status': 'PASS', 'receipt': str(output)}, indent=2))
