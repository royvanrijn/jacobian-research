#!/usr/bin/env python3
"""Exact certificate for the minimal-boundary cubic branch collapse.

This checks the algebra after a map has entered the weighted or cancellation
suspension branch.  It intentionally does not claim to verify the open
intrinsic-marking extraction or positive-chart straightening theorem.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402
from master_cancellation import hensel_jet, parameter_polynomial  # noqa: E402


# A marked cubic primitive with a double zero at 0 and a zero at 1 is unique
# up to its nonzero scalar.
a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3")
H_general = a3 * w**3 + a2 * w**2 + a1 * w + a0
constraints = (
    H_general.subs(w, 0),
    sp.diff(H_general, w).subs(w, 0),
    H_general.subs(w, 1),
)
solution = sp.solve(constraints, (a0, a1, a2), dict=True)
assert solution == [{a0: 0, a1: 0, a2: -a3}]
assert sp.factor(H_general.subs(solution[0]) - a3 * w**2 * (w - 1)) == 0


# Normalize the weighted scalar and vertical coefficient by diagonal source
# and target automorphisms, then compare the normalized lift with F_0.
c, b = sp.symbols("c b", nonzero=True)
weighted_general = WeightedSeedModel(
    sp.diff(c * w**2 * (1 - w), w), c=c, b=b
).mapping()
weighted_normal = WeightedSeedModel(
    sp.diff(w**2 * (1 - w), w), c=1, b=1
).mapping()

weighted_scaled_source = tuple(
    sp.cancel(component.subs(z, z / b)) for component in weighted_general
)
for actual, expected in zip(
    weighted_scaled_source,
    (weighted_normal[0], c * weighted_normal[1], weighted_normal[2]),
):
    assert sp.factor(actual - expected) == 0

u = 1 + x * y
foundational = (
    sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y)),
    sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)),
    sp.expand(2 * x - 3 * x**2 * y - x**3 * z),
)
assert sp.factor(
    sp.Matrix(foundational).jacobian((x, y, z)).det() + 2
) == 0
weighted_foundational_source = tuple(
    sp.expand(component.subs(z, -z / 2)) for component in weighted_normal
)
for actual, expected in zip(
    foundational,
    (
        weighted_foundational_source[0],
        2 * weighted_foundational_source[1],
        2 * weighted_foundational_source[2],
    ),
):
    assert sp.expand(actual - expected) == 0


# The cancellation degree equation r(m+1)+1=3 has only (m,r)=(1,1).
cancellation_types = [
    (m, r)
    for r in range(1, 3)
    if 2 % r == 0
    for m in (2 // r - 1,)
    if m >= 1
]
assert cancellation_types == [(1, 1)]

A_symbol, q = sp.symbols("A q")
assert sp.expand(parameter_polynomial(1, 1, q) - (q - 3)) == 0
jet = hensel_jet(1, 1, A_symbol, q)
assert sp.rem(
    sp.Poly(jet - (3 + 9 * A_symbol), q),
    sp.Poly(q - 3, q),
).is_zero


# Build the (1,1) cancellation map and verify its explicit LR equivalence
# with F_0 for an arbitrary nonzero integral scale C.
C = sp.symbols("C", nonzero=True)
A = 1 + x * y
h = 3 + 9 * A
B = sp.expand(A**2 * z + y**2 * h)
P = sp.expand(A * B)
Q = sp.expand(y + x * B)
t = sp.symbols("t")
s = x / A
R = sp.cancel(C * sp.integrate(1 - t * (Q - P * t), (t, 0, s)))
assert sp.denom(R) == 1
assert sp.factor(sp.Matrix((P, Q, R)).jacobian((x, y, z)).det() + C) == 0

cancellation_scaled_source = tuple(
    sp.expand(component.subs(z, 3 * z)) for component in (P, Q, R)
)
for actual, expected in zip(
    cancellation_scaled_source,
    (3 * foundational[0], foundational[1], C * foundational[2] / 2),
):
    assert sp.expand(actual - expected) == 0

assert sp.Poly(
    sp.integrate((1 - t * (sp.Symbol("Q") - sp.Symbol("P") * t)), t),
    t,
).degree() == 3
assert sp.Poly(c * w**2 * (1 - w), w).degree() == 3

print("PASS: every marked weighted cubic seed is c*W^2*(1-W)")
print("PASS: arbitrary weighted cubic parameters normalize to the foundational map")
print("PASS: degree three forces cancellation type (m,r)=(1,1) and jet h=3+9*A")
print("PASS: the foundational map is the explicit LR normal form of both cubic branches")
