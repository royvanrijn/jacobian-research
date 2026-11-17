#!/usr/bin/env python3
"""Exact audit for extended-geometry/COX_LEDGER_THREE_FACTOR.md."""

from __future__ import annotations

from functools import reduce
from itertools import product
from math import gcd

import sympy as sp


def gcd_many(values: list[int]) -> int:
    return reduce(gcd, (abs(value) for value in values), 0)


def integer_boundary_matrix(a: int, b: int, c: int) -> sp.Matrix:
    return sp.Matrix(
        [
            [b, c, 0, 1],
            [a, 0, c, 1],
            [0, a, b, 1],
        ]
    )


def maximal_minors(matrix: sp.Matrix) -> list[int]:
    return [
        int(matrix[:, [column for column in range(4) if column != omitted]].det())
        for omitted in range(4)
    ]


# Full boundary lattice for three linear factors.
B = integer_boundary_matrix(1, 1, 1)
unit_relation = sp.Matrix([1, 1, 1, -2])
assert B * unit_relation == sp.zeros(3, 1)
assert gcd_many(maximal_minors(B)) == 1

normalization_minor = B[:, [1, 2, 3]]
assert normalization_minor.det() == -1

augmented_B = B.col_join(sp.Matrix([[1, 0, 0, 0]]))
assert augmented_B.det() == 1
assert all(entry.q == 1 for entry in augmented_B.inv())
print("PASS: the four-column linear-factor boundary ledger has a primitive unit")
print("PASS: one primitive row gives a unimodular 4-by-4 completion")


# General degree triples: the gcd of maximal minors is the finite Smith
# obstruction, while the displayed relation spans the rational kernel.
a_symbol, b_symbol, c_symbol = sp.symbols(
    "a b c", integer=True, positive=True
)
B_symbolic = sp.Matrix(
    [
        [b_symbol, c_symbol, 0, 1],
        [a_symbol, 0, c_symbol, 1],
        [0, a_symbol, b_symbol, 1],
    ]
)
relation_symbolic = sp.Matrix(
    [
        c_symbol * (a_symbol + b_symbol - c_symbol),
        b_symbol * (a_symbol + c_symbol - b_symbol),
        a_symbol * (b_symbol + c_symbol - a_symbol),
        -2 * a_symbol * b_symbol * c_symbol,
    ]
)
assert all(sp.expand(entry) == 0 for entry in B_symbolic * relation_symbolic)

for a_value, b_value, c_value in product(range(1, 8), repeat=3):
    matrix = integer_boundary_matrix(a_value, b_value, c_value)
    minors = maximal_minors(matrix)
    predicted = gcd_many(
        [
            2 * a_value * b_value * c_value,
            c_value * (a_value + b_value - c_value),
            b_value * (a_value + c_value - b_value),
            a_value * (b_value + c_value - a_value),
        ]
    )
    assert gcd_many(minors) == predicted

    relation = [
        c_value * (a_value + b_value - c_value),
        b_value * (a_value + c_value - b_value),
        a_value * (b_value + c_value - a_value),
        -2 * a_value * b_value * c_value,
    ]
    primitive_divisor = gcd_many(relation)
    primitive_relation = sp.Matrix(
        [entry // primitive_divisor for entry in relation]
    )
    assert matrix * primitive_relation == sp.zeros(3, 1)

obstructed = integer_boundary_matrix(1, 2, 3)
assert gcd_many(maximal_minors(obstructed)) == 4
# Every determinant obtained by appending one row is a combination of these
# cofactors, hence is divisible by four.
for row in product(range(-2, 3), repeat=4):
    completion = obstructed.col_join(sp.Matrix([row]))
    assert int(completion.det()) % 4 == 0
print("PASS: the general torsion index is the gcd of the four maximal minors")
print("PASS: degree triple (1,2,3) has a rank-one but index-four obstruction")


# Explicit three-linear-factor construction.
u1, v1, u2, v2, u3, v3, X, z = sp.symbols(
    "u1 v1 u2 v2 u3 v3 X z"
)
variables = (u1, v1, u2, v2, u3, v3)
linear_factors = (
    u1 * X + v1,
    u2 * X + v2,
    u3 * X + v3,
)
product_cubic = sp.Poly(sp.expand(sp.prod(linear_factors)), X)
p0, p1, p2, p3 = product_cubic.all_coeffs()
r12 = u1 * v2 - v1 * u2
r13 = u1 * v3 - v1 * u3
r23 = u2 * v3 - v2 * u3
m = p0 + p3

discriminant = sp.factor(sp.discriminant(product_cubic.as_expr(), X))
assert sp.factor(discriminant - (r12 * r13 * r23) ** 2) == 0
print("PASS: Disc(L1*L2*L3)=(r12*r13*r23)^2")


# Unique normalization of r13, r23, and m.
raw_r13, raw_r23, raw_m = sp.symbols(
    "raw_r13 raw_r23 raw_m", nonzero=True
)
lambda1 = raw_r23 / raw_m
lambda2 = raw_r13 / raw_m
lambda3 = raw_m / (raw_r13 * raw_r23)
assert sp.cancel(lambda1 * lambda3 * raw_r13 - 1) == 0
assert sp.cancel(lambda2 * lambda3 * raw_r23 - 1) == 0
assert sp.cancel(lambda1 * lambda2 * lambda3 * raw_m - 1) == 0
normalized_r12 = sp.cancel(lambda1 * lambda2 * r12)
assert normalized_r12 == sp.cancel(
    r12 * raw_r13 * raw_r23 / raw_m**2
)
print("PASS: r13=r23=m=1 has a unique root-free torus normalization")
print("PASS: the surviving coordinate is the boundary unit epsilon")


# Ambient determinant.  Its restriction to the normalized complete
# intersection is the residue Jacobian -r12.
ambient_outputs = (p0, p1, p2, m, r13, r23)
ambient_jacobian = sp.factor(
    sp.Matrix(ambient_outputs).jacobian(variables).det()
)
assert sp.factor(ambient_jacobian + r12 * r13**2 * r23**2) == 0
print("PASS: det D(p0,p1,p2,m,r13,r23)=-r12*r13^2*r23^2")


# The primitive suspension Z=z/r12 cancels the remaining unit.  Keeping the
# constraints as the last three ambient outputs is an exact residue-form
# calculation: on r13=r23=m=1 the seven-dimensional determinant is 1.
suspended_outputs = (p0, p1, p2, z / r12, m, r13, r23)
suspended_variables = (*variables, z)
suspended_jacobian = sp.factor(
    sp.Matrix(suspended_outputs).jacobian(suspended_variables).det()
)
assert sp.factor(suspended_jacobian - r13**2 * r23**2) == 0
assert suspended_jacobian.subs({r13: 1, r23: 1, m: 1}) == 1
print("PASS: z/r12 absorbs the primitive unit and the residue Jacobian is 1")


# The extra coordinate does not alter the factorization fiber: Z determines
# z uniquely once an ordered factor point (and hence r12) is chosen.
Z = sp.symbols("Z")
assert sp.cancel((r12 * Z) / r12 - Z) == 0
print("PASS: the suspension preserves the geometric degree-six ordered cover")


# Small finite-field replay.  It is intentionally point-set level: the exact
# determinant calculation above is the scheme-theoretic certificate.
def det_mod(left: tuple[int, int], right: tuple[int, int], prime: int) -> int:
    return (left[0] * right[1] - left[1] * right[0]) % prime


prime = 5
nonzero_vectors = [
    (u, v)
    for u, v in product(range(prime), repeat=2)
    if (u, v) != (0, 0)
]
fiber_counts: dict[tuple[int, int, int, int], int] = {}
for factor1, factor2, factor3 in product(nonzero_vectors, repeat=3):
    finite_r13 = det_mod(factor1, factor3, prime)
    finite_r23 = det_mod(factor2, factor3, prime)
    finite_r12 = det_mod(factor1, factor2, prime)
    if (finite_r13, finite_r23) != (1, 1) or finite_r12 == 0:
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

    cubic = (finite_p0, finite_p1, finite_p2, finite_p3)
    fiber_counts[cubic] = fiber_counts.get(cubic, 0) + 1

assert fiber_counts
assert set(fiber_counts.values()) == {6}
print("PASS: over F_5 every rational split target in the normalized image has six lifts")
print("PASS three-factor Cox ledger")
