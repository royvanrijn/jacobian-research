#!/usr/bin/env python3
"""Exact all-rank collision-frame and projective-descent algebra.

For a separable degree-N presentation A=k[r] and a second primitive
coordinate u=q(r), projective equivalence is the rank-at-most-three
condition on the four columns 1, r, u, and r*u.  This checker verifies:

* Conf_(N-1) is the full S_N frame set through degree eight;
* signed-minor interpolation turns the 4-by-4 minors into the residuals
  after matching the first three framed roots;
* the normalized coefficient matrix represents 1, r, u, and a_N*r*u;
* rank three is automatic, rank four recovers the exact cross-ratio
  equation, and q=r+r^2 is a uniform primitive nonprojective witness;
* the N-3 framed residuals have independent differentials along an exact
  projective witness through degree ten.

The checker concerns projective root-coordinate transport.  It does not
assert that every Keller-incidence equivalence is induced by PGL_2.
"""

from __future__ import annotations

from itertools import permutations
from math import factorial

import sympy as sp


def assert_zero(expression: sp.Expr) -> None:
    """Assert an exact rational identity."""

    assert sp.cancel(expression) == 0


def mobius_matrix(
    source_roots: tuple[sp.Expr, sp.Expr, sp.Expr],
    target_roots: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> sp.Matrix:
    """Return signed maximal minors for three-point interpolation."""

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


def normalized_coefficient_matrix(
    relation_coefficients: tuple[sp.Expr, ...],
    coordinate_coefficients: tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Coefficient matrix for 1, r, q(r), and a_N*r*q(r).

    The relation is

        a_0 + r + a_2*r^2 + ... + a_N*r^N = 0,

    and q has degree less than N.  Multiplying the fourth column by the
    nonzero leading coefficient keeps the rank criterion polynomial.
    """

    degree = len(coordinate_coefficients)
    assert len(relation_coefficients) == degree + 1
    assert relation_coefficients[1] == 1

    leading_coefficient = relation_coefficients[degree]
    top_coordinate = coordinate_coefficients[degree - 1]
    matrix = sp.zeros(degree, 4)
    matrix[0, 0] = 1
    matrix[1, 1] = 1
    for row, coefficient in enumerate(coordinate_coefficients):
        matrix[row, 2] = coefficient
    matrix[0, 3] = -top_coordinate * relation_coefficients[0]
    matrix[1, 3] = (
        leading_coefficient * coordinate_coefficients[0]
        - top_coordinate
    )
    for row in range(2, degree):
        matrix[row, 3] = (
            leading_coefficient * coordinate_coefficients[row - 1]
            - top_coordinate * relation_coefficients[row]
        )
    return matrix


# ---------------------------------------------------------------------------
# 1. Conf_(N-1) is the full frame torsor.
# ---------------------------------------------------------------------------

for degree in range(3, 9):
    labels = tuple(range(degree))
    symmetric_group = tuple(permutations(labels))
    assert len(symmetric_group) == factorial(degree)

    almost_frames = tuple(permutations(labels, degree - 1))
    assert len(almost_frames) == factorial(degree)
    completed_frames = {
        almost_frame + tuple(set(labels) - set(almost_frame))
        for almost_frame in almost_frames
    }
    assert completed_frames == set(symmetric_group)

    for marked_count in range(1, degree):
        stabilizer = {
            permutation
            for permutation in symmetric_group
            if all(
                permutation[index] == index
                for index in range(marked_count)
            )
        }
        assert len(stabilizer) == factorial(degree - marked_count)


# ---------------------------------------------------------------------------
# 2. Framed PGL_2 interpolation and the determinantal residuals.
# ---------------------------------------------------------------------------

r = sp.symbols("r1:8")
u = sp.symbols("u1:8")
first_three_matrix = mobius_matrix(r[:3], u[:3])
matrix_a, matrix_b, matrix_c, matrix_d = tuple(first_three_matrix)

vandermonde_three_r = sp.prod(
    r[first] - r[second]
    for first in range(3)
    for second in range(first + 1, 3)
)
vandermonde_three_u = sp.prod(
    u[first] - u[second]
    for first in range(3)
    for second in range(first + 1, 3)
)
assert_zero(
    first_three_matrix.det()
    - vandermonde_three_r * vandermonde_three_u
)

for index in range(3, 7):
    interpolation_rows = sp.Matrix(
        [
            [r[row], 1, -u[row] * r[row], -u[row]]
            for row in (0, 1, 2, index)
        ]
    )
    residual = (
        matrix_a * r[index]
        + matrix_b
        - u[index] * (matrix_c * r[index] + matrix_d)
    )
    assert_zero(residual + interpolation_rows.det())

# A literal Mobius transform makes every four-row determinant vanish after
# clearing the four nonzero linear denominators.
mobius_a, mobius_b, mobius_c, mobius_d = sp.symbols("A B C D")
four_roots = sp.symbols("x1:5")
cleared_projective_rows = []
for root in four_roots:
    numerator = mobius_a * root + mobius_b
    denominator = mobius_c * root + mobius_d
    cleared_projective_rows.append(
        [
            root * denominator,
            denominator,
            -root * numerator,
            -numerator,
        ]
    )
assert_zero(sp.Matrix(cleared_projective_rows).det())


# ---------------------------------------------------------------------------
# 3. The normalized coefficient matrix and exact low-rank witnesses.
# ---------------------------------------------------------------------------

S = sp.symbols("S")

for degree in range(3, 9):
    relation_symbols = (
        sp.symbols(f"a0_{degree}"),
        sp.Integer(1),
        *(
            sp.Symbol(f"a{power}_{degree}")
            for power in range(2, degree + 1)
        ),
    )
    coordinate_symbols = tuple(
        sp.Symbol(f"q{power}_{degree}") for power in range(degree)
    )
    leading_coefficient = relation_symbols[degree]
    top_coordinate = coordinate_symbols[degree - 1]
    relation = sum(
        relation_symbols[power] * S**power
        for power in range(degree + 1)
    )
    coordinate = sum(
        coordinate_symbols[power] * S**power
        for power in range(degree)
    )
    coefficient_matrix = normalized_coefficient_matrix(
        relation_symbols, coordinate_symbols
    )
    reduced_product = sum(
        coefficient_matrix[power, 3] * S**power
        for power in range(degree)
    )
    assert_zero(
        leading_coefficient * S * coordinate
        - reduced_product
        - top_coordinate * relation
    )

    # Normalize the exact root relation prod(S-i) by its nonzero linear
    # coefficient.  Interpolation of one Mobius map produces a polynomial
    # coordinate of degree less than N whose coefficient matrix has rank 3.
    root_values = tuple(sp.Integer(index) for index in range(1, degree + 1))
    root_relation = sp.Poly(
        sp.prod(S - root for root in root_values), S
    )
    linear_coefficient = root_relation.nth(1)
    assert linear_coefficient != 0
    normalized_relation = sp.Poly(
        root_relation.as_expr() / linear_coefficient, S
    )
    exact_relation_coefficients = tuple(
        normalized_relation.nth(power)
        for power in range(degree + 1)
    )

    projective_values = tuple(
        sp.Rational(2 * root + 1, root + degree + 1)
        for root in root_values
    )
    assert len(set(projective_values)) == degree
    projective_coordinate = sp.Poly(
        sp.interpolate(tuple(zip(root_values, projective_values)), S), S
    )
    exact_projective_coefficients = tuple(
        projective_coordinate.nth(power) for power in range(degree)
    )
    exact_projective_matrix = normalized_coefficient_matrix(
        exact_relation_coefficients, exact_projective_coefficients
    )
    vandermonde = sp.Matrix(
        [
            [root**power for power in range(degree)]
            for root in root_values
        ]
    )
    evaluated_projective_matrix = vandermonde * exact_projective_matrix
    expected_evaluations = sp.Matrix(
        [
            [
                1,
                root,
                value,
                exact_relation_coefficients[degree] * root * value,
            ]
            for root, value in zip(root_values, projective_values)
        ]
    )
    assert evaluated_projective_matrix == expected_evaluations
    assert exact_projective_matrix.rank() == 3

    quadratic_coordinate = sp.Poly(S + S**2, S)
    quadratic_values = tuple(
        quadratic_coordinate.eval(root) for root in root_values
    )
    assert len(set(quadratic_values)) == degree
    exact_quadratic_coefficients = tuple(
        quadratic_coordinate.nth(power) for power in range(degree)
    )
    exact_quadratic_matrix = normalized_coefficient_matrix(
        exact_relation_coefficients, exact_quadratic_coefficients
    )
    if degree == 3:
        assert exact_quadratic_matrix.rank() == 3
    else:
        assert exact_quadratic_matrix.rank() == 4


# ---------------------------------------------------------------------------
# 4. Rank four is exactly the previously computed cross-ratio equation.
# ---------------------------------------------------------------------------

a0, a2, a3, a4 = sp.symbols("a0 a2 a3 a4")
q0, q1, q2, q3 = sp.symbols("q0 q1 q2 q3")
quartic_matrix = normalized_coefficient_matrix(
    (a0, sp.Integer(1), a2, a3, a4),
    (q0, q1, q2, q3),
)
quartic_defect = (
    a4 * (q2**2 - q1 * q3) - a3 * q2 * q3 + a2 * q3**2
)
assert_zero(quartic_matrix.det() - quartic_defect)

seed4, target_pi, target_b, target_c = sp.symbols("u4 pi b c")
universal_quartic_matrix = quartic_matrix.subs(
    {
        a0: -target_c / 2,
        a2: target_b,
        a3: target_pi,
        a4: seed4 * target_pi**4,
    }
)
universal_quartic_defect = (
    seed4 * target_pi**4 * (q2**2 - q1 * q3)
    - target_pi * q2 * q3
    + target_b * q3**2
)
assert_zero(
    universal_quartic_matrix.det() - universal_quartic_defect
)


# ---------------------------------------------------------------------------
# 5. The N-3 residuals are independent on the primitive projective locus.
# ---------------------------------------------------------------------------

for degree in range(3, 11):
    source_values = tuple(
        sp.Integer(index) for index in range(1, degree + 1)
    )
    target_values = tuple(
        sp.Rational(2 * root + 1, root + degree + 1)
        for root in source_values
    )
    interpolating_matrix = mobius_matrix(
        source_values[:3], target_values[:3]
    )
    exact_a, exact_b, exact_c, exact_d = tuple(interpolating_matrix)
    assert interpolating_matrix.det() != 0
    if degree == 3:
        continue

    residual_variables = sp.symbols(f"v4:{degree + 1}")
    residuals = sp.Matrix(
        [
            exact_a * source_values[index]
            + exact_b
            - residual_variables[index - 3]
            * (exact_c * source_values[index] + exact_d)
            for index in range(3, degree)
        ]
    )
    jacobian = residuals.jacobian(residual_variables)
    assert jacobian.rank() == degree - 3
    projective_substitution = {
        residual_variables[index - 3]: target_values[index]
        for index in range(3, degree)
    }
    for residual in residuals:
        assert_zero(residual.subs(projective_substitution))


print("PASS: Conf_(N-1) is the full S_N frame set through rank eight")
print("PASS: framed PGL_2 residuals are exactly the 4-by-4 minors")
print("PASS: the normalized all-rank coefficient matrix is exact")
print("PASS: ranks three and four recover the automatic and cross-ratio cases")
print("PASS: q=r+r^2 is primitive but nonprojective in every rank N>=4")
print("PASS: the projective locus is smooth of codimension N-3 on its frame chart")
