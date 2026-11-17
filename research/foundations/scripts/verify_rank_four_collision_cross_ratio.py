#!/usr/bin/env python3
"""Exact rank-four collision-frame and cross-ratio obstruction algebra.

This checker proves the scoped projective-root statement used by
RANK_FOUR_COLLISION_CROSS_RATIO.md.  It does not assert that every Keller
incidence equivalence must arise from PGL_2.
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp


def assert_zero(expression: sp.Expr) -> None:
    """Assert an exact rational identity."""

    assert sp.cancel(expression) == 0


def mobius_matrix(
    source_roots: tuple[sp.Expr, sp.Expr, sp.Expr],
    target_roots: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> sp.Matrix:
    """Signed maximal minors for the unique three-point interpolation."""

    interpolation = sp.Matrix(
        [
            [
                source_roots[index],
                1,
                -target_roots[index] * source_roots[index],
                -target_roots[index],
            ]
            for index in range(3)
        ]
    )
    coefficients = []
    for column in range(4):
        remaining = [index for index in range(4) if index != column]
        coefficients.append(
            (-1) ** column * interpolation[:, remaining].det()
        )
    return sp.Matrix(
        [
            [coefficients[0], coefficients[1]],
            [coefficients[2], coefficients[3]],
        ]
    )


# ---------------------------------------------------------------------------
# 1. Ordered pairs are not frames in rank four, while ordered triples are.
# ---------------------------------------------------------------------------

labels = tuple(range(4))
ordered_pairs = {
    (first, second)
    for first in labels
    for second in labels
    if first != second
}
ordered_triples = {
    (first, second, third)
    for first in labels
    for second in labels
    for third in labels
    if len({first, second, third}) == 3
}
assert len(ordered_pairs) == 12
assert len(ordered_triples) == 24
assert all(len(set(labels) - set(pair)) == 2 for pair in ordered_pairs)

frames_from_triples = {
    triple + tuple(set(labels) - set(triple))
    for triple in ordered_triples
}
assert frames_from_triples == set(permutations(labels))

# Pointwise stabilizers in S_4 have the expected sizes: an ordered pair
# leaves a transposition, while an ordered triple leaves only the identity.
symmetric_group = tuple(permutations(labels))
pair_stabilizer = {
    permutation
    for permutation in symmetric_group
    if (permutation[0], permutation[1]) == (0, 1)
}
triple_stabilizer = {
    permutation
    for permutation in symmetric_group
    if (permutation[0], permutation[1], permutation[2]) == (0, 1, 2)
}
assert len(pair_stabilizer) == 2
assert len(triple_stabilizer) == 1


# ---------------------------------------------------------------------------
# 2. The fourth-point residual and cross-ratio defect.
# ---------------------------------------------------------------------------

r = sp.symbols("r1:5")
q0, q1, q2, q3 = sp.symbols("q0 q1 q2 q3")
u = tuple(
    q0 + q1 * root + q2 * root**2 + q3 * root**3
    for root in r
)

first_three_matrix = mobius_matrix(r[:3], u[:3])
matrix_a, matrix_b, matrix_c, matrix_d = tuple(first_three_matrix)

vandermonde_three_r = (
    (r[0] - r[1]) * (r[0] - r[2]) * (r[1] - r[2])
)
vandermonde_three_u = (
    (u[0] - u[1]) * (u[0] - u[2]) * (u[1] - u[2])
)
assert_zero(
    first_three_matrix.det()
    - vandermonde_three_r * vandermonde_three_u
)

e1 = sum(r)
e2 = sum(
    r[first] * r[second]
    for first in range(4)
    for second in range(first + 1, 4)
)
e3 = sum(
    r[first] * r[second] * r[third]
    for first in range(4)
    for second in range(first + 1, 4)
    for third in range(second + 1, 4)
)
e4 = sp.prod(r)

projective_defect = (
    q2**2
    - q1 * q3
    + q2 * q3 * e1
    + q3**2 * e2
)
vandermonde_four = sp.prod(
    r[first] - r[second]
    for first in range(4)
    for second in range(first + 1, 4)
)

fourth_residual = (
    matrix_a * r[3]
    + matrix_b
    - u[3] * (matrix_c * r[3] + matrix_d)
)
assert_zero(
    fourth_residual + vandermonde_four * projective_defect
)

# The same defect is the numerator of equality of the two labeled cross
# ratios.  Denominators are products of root differences and divided
# differences, all units on the framed primitive overlap.
cross_ratio_equality_numerator = (
    (u[0] - u[2])
    * (u[1] - u[3])
    * (r[0] - r[3])
    * (r[1] - r[2])
    - (r[0] - r[2])
    * (r[1] - r[3])
    * (u[0] - u[3])
    * (u[1] - u[2])
)
assert_zero(
    cross_ratio_equality_numerator
    + vandermonde_four * projective_defect
)


# ---------------------------------------------------------------------------
# 3. Primitive-element boundary and the universal quartic target form.
# ---------------------------------------------------------------------------


def divided_difference(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    """(q(first)-q(second))/(first-second) without division."""

    return (
        q1
        + q2 * (first + second)
        + q3 * (first**2 + first * second + second**2)
    )


primitive_determinant = sp.prod(
    divided_difference(r[first], r[second])
    for first in range(4)
    for second in range(first + 1, 4)
)
transformed_vandermonde_four = sp.prod(
    u[first] - u[second]
    for first in range(4)
    for second in range(first + 1, 4)
)
assert_zero(
    transformed_vandermonde_four
    - vandermonde_four * primitive_determinant
)

# Normalize a quartic relation to have linear coefficient one:
# E=S+b*S^2+pi*S^3+u4*pi^4*S^4-c/2.
quartic_parameter = -e3**3 / e1**4
target_pi = e1 / e3
target_b = -e2 / e3
target_c = 2 * e4 / e3
leading_coefficient = -1 / e3
assert_zero(
    quartic_parameter * target_pi**4 - leading_coefficient
)

S = sp.symbols("S")
normalized_quartic = (
    S
    + target_b * S**2
    + target_pi * S**3
    + quartic_parameter * target_pi**4 * S**4
    - target_c / 2
)
for root in r:
    assert_zero(normalized_quartic.subs(S, root))

# Clearing the invertible universal-quartic coefficient gives a polynomial
# equation on the actual (u4,pi,b)-target chart.
cleared_projective_defect = (
    quartic_parameter
    * target_pi**4
    * (q2**2 - q1 * q3)
    - q2 * q3 * target_pi
    + target_b * q3**2
)
assert_zero(
    cleared_projective_defect
    - quartic_parameter * target_pi**4 * projective_defect
)

# A genuinely quadratic change has defect q2^2.  Over a reduced
# characteristic-zero base its projective locus is therefore q2=0.
assert_zero(projective_defect.subs(q3, 0) - q2**2)


# ---------------------------------------------------------------------------
# 4. Exact witness cards on one normalized quartic root configuration.
# ---------------------------------------------------------------------------

root_witness = {r[index]: index + 1 for index in range(4)}

# q(r)=r+r^2 is primitive on {1,2,3,4} but is not projective.
quadratic_witness = {
    q0: 0,
    q1: 1,
    q2: 1,
    q3: 0,
    **root_witness,
}
quadratic_values = tuple(
    sp.expand(value.subs(quadratic_witness)) for value in u
)
assert len(set(quadratic_values)) == 4
assert projective_defect.subs(quadratic_witness) == 1
assert fourth_residual.subs(quadratic_witness) != 0

# On the same roots, q(r)=35*r+r^3 lies on the non-affine projective conic
# and remains primitive.
cubic_projective_witness = {
    q0: 0,
    q1: 35,
    q2: 0,
    q3: 1,
    **root_witness,
}
cubic_values = tuple(
    sp.expand(value.subs(cubic_projective_witness)) for value in u
)
assert len(set(cubic_values)) == 4
assert projective_defect.subs(cubic_projective_witness) == 0
assert primitive_determinant.subs(cubic_projective_witness) != 0
assert fourth_residual.subs(cubic_projective_witness) == 0

print("PASS: rank-four ordered triples, not pairs, are full S4 frames")
print("PASS: the fourth-root interpolation residual is -Vandermonde*Psi")
print("PASS: Psi is exactly the labeled cross-ratio defect")
print("PASS: the primitive and projective boundaries are separated")
print("PASS: the universal quartic target equation for Psi is polynomial")
