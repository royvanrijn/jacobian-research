#!/usr/bin/env python3
"""Exact audit of the oriented cubic chart's affine descent."""

from __future__ import annotations

import sympy as sp


X = sp.symbols("X")
u1, v1, u2, v2, u3, v3 = sp.symbols(
    "u1 v1 u2 v2 u3 v3"
)
variables = (u1, v1, u2, v2, u3, v3)
L1 = u1 * X + v1
L2 = u2 * X + v2
L3 = u3 * X + v3
Q = sp.expand(L1 * L2)
P = sp.expand(Q * L3)
p0, p1, p2, p3 = sp.Poly(P, X).all_coeffs()
r12 = u1 * v2 - v1 * u2
r13 = u1 * v3 - v1 * u3
r23 = u2 * v3 - v2 * u3


# The linear--quadratic quotient lands on the normalized slice.
resultant_L3_Q = sp.factor(sp.resultant(L3, Q, X))
assert sp.factor(resultant_L3_Q - r13 * r23) == 0
assert sp.Poly(P, X).nth(2) == p1
print("PASS: r13=r23=p1=1 descends to Res(L3,Q)=p1=1")


# The discriminant double cover is exactly ordering the roots of Q.
disc_P = sp.factor(sp.discriminant(P, X))
disc_Q = sp.factor(sp.discriminant(Q, X))
assert sp.factor(disc_Q - r12**2) == 0
assert sp.factor(disc_P - resultant_L3_Q**2 * disc_Q) == 0
print("PASS: the source and target orientation double covers are Cartesian")


# Swapping L1,L2 fixes Q and P, negates D=r12, and swaps the two selected
# collision resultants.
swap = {
    u1: u2,
    v1: v2,
    u2: u1,
    v2: v1,
}
simultaneous_swap = {
    expression: sp.expand(expression.xreplace(swap))
    for expression in (Q, P, r12, r13, r23)
}
assert simultaneous_swap[Q] == Q
assert simultaneous_swap[P] == P
assert simultaneous_swap[r12] == -r12
assert simultaneous_swap[r13] == r23
assert simultaneous_swap[r23] == r13
print("PASS: the involution negates orientation and exchanges the dicritical primes")


# Independent linear--quadratic ambient determinant.
a, b, c, d, e = sp.symbols("a b c d e")
linear = a * X + b
quadratic = c * X**2 + d * X + e
linear_quadratic = sp.Poly(sp.expand(linear * quadratic), X).all_coeffs()
resultant = sp.factor(sp.resultant(linear, quadratic, X))
ambient_outputs = (*linear_quadratic, resultant)
ambient_jacobian = sp.factor(
    sp.Matrix(ambient_outputs).jacobian((a, b, c, d, e)).det()
)
assert sp.factor(ambient_jacobian + resultant**2) == 0
print("PASS: the quotient slice has residue Jacobian -1")


# Reuse the explicit A^3 coordinates of the normalized factorization slice.
y, z = sp.symbols("y z")
b_forward = 1 + a * y
c_forward = 1 - sp.Rational(3, 2) * a * y + a**2 * z
d_forward = (
    sp.Rational(1, 2) * y
    - a * z
    + sp.Rational(3, 2) * a * y**2
    - a**2 * y * z
)
e_forward = (
    -2 * z
    + 4 * y**2
    - 4 * a * y * z
    + 3 * a * y**3
    - 2 * a**2 * y**2 * z
)
forward = {
    b: b_forward,
    c: c_forward,
    d: d_forward,
    e: e_forward,
}
middle = a * d + b * c
assert sp.expand((middle - 1).subs(forward)) == 0
assert sp.expand((resultant - 1).subs(forward)) == 0

G = tuple(
    sp.expand(component.subs(forward))
    for component in (
        linear_quadratic[0],
        linear_quadratic[2],
        linear_quadratic[3],
    )
)
assert sp.factor(sp.Matrix(G).jacobian((a, y, z)).det()) == -1


# Exact linear equivalence with the foundational map.
z1, z2, z3 = sp.symbols("z1 z2 z3")
unit = 1 + z1 * z2
foundational = (
    unit**3 * z3 + z2**2 * unit * (4 + 3 * z1 * z2),
    z2 + 3 * z1 * unit**2 * z3 + 3 * z1 * z2**2 * (4 + 3 * z1 * z2),
    2 * z1 - 3 * z1**2 * z2 - z1**3 * z3,
)
G_sub = tuple(
    sp.expand(component.subs({a: z1, y: z2, z: -z3 / 2}))
    for component in G
)
assert sp.expand(foundational[0] - G_sub[2]) == 0
assert sp.expand(foundational[1] - 2 * G_sub[1]) == 0
assert sp.expand(foundational[2] - 2 * G_sub[0]) == 0
print("PASS: affine descent is exactly the foundational cubic Keller map")
print("PASS oriented cubic affine descent")
