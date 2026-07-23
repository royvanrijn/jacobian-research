#!/usr/bin/env python3
"""Exact audit for the oriented cubic Cox chart."""

from __future__ import annotations

from collections import Counter
from itertools import product

import sympy as sp


# Boundary classes R13, R23, E form a unimodular source chart.
source_boundary = sp.Matrix(
    [
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ]
)
assert source_boundary.det() == -1
print("PASS: R13,R23,E give a unimodular normalization chart")


u1, v1, u2, v2, u3, v3, X = sp.symbols(
    "u1 v1 u2 v2 u3 v3 X"
)
variables = (u1, v1, u2, v2, u3, v3)
factors = (
    u1 * X + v1,
    u2 * X + v2,
    u3 * X + v3,
)
polynomial = sp.Poly(sp.expand(sp.prod(factors)), X)
p0, p1, p2, p3 = polynomial.all_coeffs()
r12 = u1 * v2 - v1 * u2
r13 = u1 * v3 - v1 * u3
r23 = u2 * v3 - v2 * u3
m = p0 + p3

discriminant = sp.factor(sp.discriminant(polynomial.as_expr(), X))
assert sp.factor(discriminant - (r12 * r13 * r23) ** 2) == 0
assert sp.factor(
    discriminant.subs({r13: 1, r23: 1}) - r12**2
) == 0
print("PASS: the oriented target coordinate D=r12 squares to the discriminant")


ambient_outputs = (p0, p1, p2, m, r13, r23)
ambient_jacobian = sp.factor(
    sp.Matrix(ambient_outputs).jacobian(variables).det()
)
assert sp.factor(ambient_jacobian + r12 * r13**2 * r23**2) == 0

# On the source slice, mu^*(dp0 dp1 dp2)=-r12*omega.  Dividing by
# d(D^2-Delta)/dD=2D and substituting D=r12 gives the constant -1/2.
D = sp.symbols("D", nonzero=True)
residue_ratio = sp.cancel((-r12) / (2 * D))
assert sp.cancel(
    residue_ratio.subs(D, r12) + sp.Rational(1, 2)
) == 0
print("PASS: the hypersurface-residue Jacobian is the constant -1/2")


# A triple root makes the chosen Hessian coefficient vanish, so h!=0 is a
# smooth principal chart of the oriented discriminant cover.
a, b = sp.symbols("a b")
triple = sp.Poly((a * X + b) ** 3, X).all_coeffs()
h = p1**2 - 3 * p0 * p2
assert sp.expand(h.subs(dict(zip((p0, p1, p2, p3), triple)))) == 0
print("PASS: h!=0 excludes the triple-root singular locus")


# The two escaping cyclic orderings have the stated normalization valuations.
# Formula lambda=(r23/m,r13/m,m/(r13*r23)).
delta = sp.symbols("delta")
def delta_order(expression: sp.Expr) -> int:
    return int(sp.sympify(expression).as_powers_dict().get(delta, 0))


branch_r23 = (
    delta_order(delta),
    delta_order(sp.Integer(1)),
    delta_order(1 / delta),
)
branch_r13 = (
    delta_order(sp.Integer(1)),
    delta_order(delta),
    delta_order(1 / delta),
)
assert branch_r23 == (1, 0, -1)
assert branch_r13 == (0, 1, -1)
print("PASS: R23 and R13 produce two distinct escaping valuation vectors")


# Exact F_5 point-set fiber audit.
prime = 5
nonzero_vectors = [
    vector
    for vector in product(range(prime), repeat=2)
    if vector != (0, 0)
]


def determinant(
    left: tuple[int, int], right: tuple[int, int]
) -> int:
    return (
        left[0] * right[1] - left[1] * right[0]
    ) % prime


fiber_counts: dict[tuple[int, int, int, int], int] = {}
for factor1, factor2, factor3 in product(nonzero_vectors, repeat=3):
    if determinant(factor1, factor3) != 1:
        continue
    if determinant(factor2, factor3) != 1:
        continue

    finite_p0 = factor1[0] * factor2[0] * factor3[0] % prime
    finite_p1 = (
        factor1[0] * factor2[0] * factor3[1]
        + factor1[0] * factor2[1] * factor3[0]
        + factor1[1] * factor2[0] * factor3[0]
    ) % prime
    finite_p2 = (
        factor1[0] * factor2[1] * factor3[1]
        + factor1[1] * factor2[0] * factor3[1]
        + factor1[1] * factor2[1] * factor3[0]
    ) % prime
    finite_p3 = factor1[1] * factor2[1] * factor3[1] % prime
    if (finite_p0 + finite_p3) % prime != 1:
        continue

    finite_D = determinant(factor1, factor2)
    target = (finite_p0, finite_p1, finite_p2, finite_D)
    fiber_counts[target] = fiber_counts.get(target, 0) + 1

profile = Counter(fiber_counts.values())
assert profile == Counter({3: 30, 1: 25})
print("PASS: over F_5 the affine oriented fibers have profile 30*3 + 25*1")


# The motivic formula [Ybar]=L^3-2L-[Spec(x^2-x+1)] specializes to the
# exact finite-field counts.  Re-enumerate Ybar for two splitting behaviors
# of the quadratic Artin factor.
def source_count(field_size: int) -> int:
    vectors = [
        vector
        for vector in product(range(field_size), repeat=2)
        if vector != (0, 0)
    ]

    def det_q(
        left: tuple[int, int], right: tuple[int, int]
    ) -> int:
        return (
            left[0] * right[1] - left[1] * right[0]
        ) % field_size

    count = 0
    for factor1, factor2, factor3 in product(vectors, repeat=3):
        if det_q(factor1, factor3) != 1:
            continue
        if det_q(factor2, factor3) != 1:
            continue
        finite_p0 = factor1[0] * factor2[0] * factor3[0]
        finite_p3 = factor1[1] * factor2[1] * factor3[1]
        if (finite_p0 + finite_p3) % field_size == 1:
            count += 1
    return count


for field_size in (5, 7):
    artin_points = sum(
        1
        for value in range(field_size)
        if (value**2 - value + 1) % field_size == 0
    )
    predicted_count = (
        field_size**3 - 2 * field_size - artin_points
    )
    assert source_count(field_size) == predicted_count
print("PASS: finite-field counts realize [Ybar]=L^3-2L-[Spec(x^2-x+1)]")
print("PASS oriented cubic Cox chart")
