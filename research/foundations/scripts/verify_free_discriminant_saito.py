#!/usr/bin/env python3
"""Exact Saito-matrix checks for the first marked-root constructions.

The three-dimensional divisor used here is the reduced *ledger divisor*: the
reduced branch discriminant together with the additional divisorial
reconstruction boundary when that boundary is present.  This distinction is
essential for the weighted and cancellation quartics.
"""

from __future__ import annotations

import sympy as sp


T = sp.symbols("T")


def discriminant_without_unit(polynomial: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Return the primitive squarefree part of a discriminant."""

    discriminant = sp.factor(sp.discriminant(polynomial, variable))
    _, factors = sp.factor_list(discriminant)
    return sp.factor(sp.prod(factor for factor, _ in factors))


def assert_saito(
    divisor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    matrix: sp.Matrix,
    determinant_unit: sp.Rational | sp.Integer,
) -> None:
    """Verify Saito's determinant and logarithmic-derivation criteria."""

    assert matrix.rows == matrix.cols == len(variables)
    assert sp.expand(matrix.det() - determinant_unit * divisor) == 0
    for column in range(matrix.cols):
        derivative = sum(
            matrix[row, column] * sp.diff(divisor, variables[row])
            for row in range(matrix.rows)
        )
        quotient = sp.cancel(derivative / divisor)
        assert sp.denom(quotient) == 1


def root_lift(
    inverse_equation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    coefficients: tuple[sp.Expr, ...],
) -> sp.Expr:
    """Lift a logarithmic target derivation to the marked-root incidence.

    The returned polynomial tau is the unique representative of degree below
    deg_T(E) for which delta(E) + tau E_T vanishes modulo E.
    """

    fraction_field = sp.QQ.frac_field(*variables)
    equation = sp.Poly(inverse_equation, T, domain=fraction_field)
    derivative = sp.Poly(sp.diff(inverse_equation, T), T, domain=fraction_field)
    derivative_inverse = sp.invert(derivative, equation).as_expr()
    target_derivative = sum(
        coefficient * sp.diff(inverse_equation, variable)
        for coefficient, variable in zip(coefficients, variables, strict=True)
    )
    lift = sp.rem(
        sp.Poly(
            sp.cancel(-target_derivative * derivative_inverse),
            T,
            domain=fraction_field,
        ),
        equation,
    ).as_expr()
    assert sp.denom(sp.cancel(lift)) == 1
    remainder = sp.rem(
        sp.Poly(
            sp.expand(target_derivative + lift * sp.diff(inverse_equation, T)),
            T,
            domain=fraction_field,
        ),
        equation,
    )
    assert remainder.is_zero
    return sp.factor(lift)


# Foundational cubic, in normalized quadratic-gauge target coordinates.
P, B, C = sp.symbols("P B C")
cubic_equation = P * T**3 - B * T**2 / 2 + T - C / 2
cubic_branch = (
    B**3 * C - B**2 - 18 * B * C * P + 27 * C**2 * P**2 + 16 * P
)
assert sp.expand(sp.discriminant(cubic_equation, T) + cubic_branch / 4) == 0
cubic_saito = sp.Matrix(
    [
        [2 * P, -2 * B, 0],
        [B, 3 * B * C - 16, B**2 - 12 * P],
        [-C, 9 * C**2, 3 * B * C - 4],
    ]
)
assert_saito(cubic_branch, (P, B, C), cubic_saito, 8)
cubic_lifts = [
    root_lift(cubic_equation, (P, B, C), tuple(cubic_saito[:, column]))
    for column in range(3)
]
assert cubic_lifts == [-T, T * (3 * C + 4 * T), B * T - 2]


# First noncubic weighted construction: H(W)=W^3(1-W).
A, B, C = sp.symbols("A B C")
weighted_equation = T**3 * (1 - T) - B * C * T + A * C**2
weighted_branch = (
    256 * A**3 * C**3
    - 192 * A**2 * B * C**2
    + 27 * A**2 * C
    - 6 * A * B**2 * C
    + 27 * B**4 * C
    - 4 * B**3
)
weighted_ledger = C * weighted_branch
assert sp.expand(
    sp.discriminant(weighted_equation, T) + C**3 * weighted_branch
) == 0
weighted_saito = sp.Matrix(
    [
        [2 * A, -2 * B**2, -14 * B**2],
        [
            B,
            8 * B**2 * C - 8 * A * C - B,
            128 * A * B * C**2 - 64 * A * C + B,
        ],
        [
            -C,
            16 * B * C**2 - 3 * C,
            256 * A * C**3 - 48 * B * C**2 + 3 * C,
        ],
    ]
)
assert_saito(weighted_ledger, (A, B, C), weighted_saito, -16)
weighted_lifts = [
    root_lift(weighted_equation, (A, B, C), tuple(weighted_saito[:, column]))
    for column in range(3)
]
assert weighted_lifts[0] == 0


# First noncubic cancellation construction: (m,r)=(2,1), with C=1.
P, Q, R = sp.symbols("P Q R")
cancellation_equation = (
    T
    - Q**2 * T**2 / 2
    + 2 * P * Q * T**3 / 3
    - P**2 * T**4 / 4
    - R
)
cancellation_branch = (
    1728 * P**4 * R**3
    - 3456 * P**3 * Q * R**2
    + 288 * P**2 * Q**4 * R**2
    + 1656 * P**2 * Q**2 * R
    - 729 * P**2
    - 288 * P * Q**5 * R
    + 136 * P * Q**3
    + 12 * Q**8 * R
    - 6 * Q**6
)
cancellation_ledger = P * cancellation_branch
assert sp.expand(
    sp.discriminant(cancellation_equation, T)
    - P**2 * cancellation_branch / 432
) == 0
cancellation_a = (
    -384 * P**2 * Q * R**2 - 8 * P * Q**2 * R - 3 * Q**3 + 54 * P
)
cancellation_b = (
    144 * Q**4 * R**2
    + 192 * P**2 * R**3
    - 608 * P * Q * R**2
    + 24 * Q**2 * R
    - 45
)
cancellation_c = -96 * P**2 * Q**2 * R - 2 * P * Q**3 + 81 * P**2
cancellation_d = (
    36 * Q**5 * R
    + 48 * P**2 * Q * R**2
    - 152 * P * Q**2 * R
    - 18 * Q**3
    + 54 * P
)
cancellation_saito = sp.Matrix(
    [
        [3 * P, 0, 0],
        [Q, cancellation_a, cancellation_c],
        [-2 * R, cancellation_b, cancellation_d],
    ]
)
assert_saito(cancellation_ledger, (P, Q, R), cancellation_saito, -27)
cancellation_lifts = [
    root_lift(
        cancellation_equation,
        (P, Q, R),
        tuple(cancellation_saito[:, column]),
    )
    for column in range(3)
]
assert cancellation_lifts[0] == -2 * T


# First noncubic quadratic gauge.  The full three-dimensional branch and
# ledger divisors fail the perfect-Jacobian-ideal freeness test in the
# companion Singular certificate.
P, B, C = sp.symbols("P B C")
quadratic_full_equation = (
    P**4 * T**4 - 2 * P * T**3 + (-P - B) * T**2 + 2 * T - C
)
quadratic_full_branch = (
    B**4 * C * P**2
    + 4 * B**3 * C * P**3
    + B**3 * C
    - B**3 * P**2
    + 8 * B**2 * C**2 * P**6
    + 6 * B**2 * C * P**4
    + 20 * B**2 * C * P**3
    + 3 * B**2 * C * P
    - 3 * B**2 * P**3
    - B**2
    + 16 * B * C**2 * P**7
    + 36 * B * C**2 * P**4
    - 36 * B * C * P**6
    + 4 * B * C * P**5
    + 40 * B * C * P**4
    + 3 * B * C * P**2
    + 18 * B * C * P
    - 3 * B * P**4
    - 18 * B * P**3
    - 2 * B * P
    + 16 * C**3 * P**10
    + 8 * C**2 * P**8
    - 48 * C**2 * P**7
    + 36 * C**2 * P**5
    + 27 * C**2 * P**2
    - 36 * C * P**7
    + C * P**6
    + 20 * C * P**5
    - 6 * C * P**4
    + C * P**3
    + 18 * C * P**2
    + 27 * P**6
    - P**5
    - 18 * P**4
    - P**2
    - 16 * P
)
assert sp.expand(
    sp.discriminant(quadratic_full_equation, T)
    + 16 * P**2 * quadratic_full_branch
) == 0
assert sp.factor_list(quadratic_full_branch)[1] == [(quadratic_full_branch, 1)]

# Every fixed-P plane section is free; here is an exact rank-minimal basis on
# the construction's canonical P=1 section.
x, y = sp.symbols("x y")
quadratic_equation = T**4 - 2 * T**3 + (-1 - x) * T**2 + 2 * T - y
quadratic_section = (
    x**4 * y
    + 5 * x**3 * y
    - x**3
    + 8 * x**2 * y**2
    + 29 * x**2 * y
    - 4 * x**2
    + 52 * x * y**2
    + 29 * x * y
    - 23 * x
    + 16 * y**3
    + 23 * y**2
    - 2 * y
    - 9
)
assert sp.expand(
    sp.discriminant(quadratic_equation, T) + 16 * quadratic_section
) == 0
v0 = sp.Matrix(
    [
        x**3 + 4 * x**2 + 4 * x * y + 19 * x + 10 * y + 10,
        2 * x**2 * y + 7 * x * y + 8 * y**2 - 3 * x + y - 7,
    ]
)
v1 = sp.Matrix(
    [
        4 * x**3 * y
        + 31 * x**2 * y
        + 16 * x * y**2
        + 15 * x**2
        + 201 * x * y
        + 20 * y**2
        + 25 * x
        + 100 * y
        + 80,
        8 * x**2 * y**2
        + 68 * x * y**2
        + 32 * y**3
        + 128 * x * y
        + 249 * y**2
        + 62 * y
        - 155,
    ]
)
v2 = sp.Matrix(
    [
        2 * x**4 + 8 * x**2 * y + x**2 - 52 * x * y - 147 * x - 105 * y - 105,
        4 * x**3 * y
        + 3 * x**2 * y
        + 16 * x * y**2
        - 6 * x**2
        - 84 * x * y
        - 124 * y**2
        + 5 * x
        - 28 * y
        + 96,
    ]
)
quadratic_section_saito = sp.Matrix.hstack(v1 - 8 * v2, 11 * v0 + 2 * v2)
assert_saito(quadratic_section, (x, y), quadratic_section_saito, -1500)
for column in range(2):
    root_lift(
        quadratic_equation,
        (x, y),
        tuple(quadratic_section_saito[:, column]),
    )


# First external control: the type-A3 reflection discriminant, written as the
# discriminant of the depressed quartic T^4+pT^2+qT+r.  This passes the
# finite-flat, free-divisor, explicit-normalization, and regular-lift gates;
# the missing datum is a reciprocal affine reconstruction chart.
p, q, r = sp.symbols("p q r")
reflection_equation = T**4 + p * T**2 + q * T + r
reflection_discriminant = (
    16 * p**4 * r
    - 4 * p**3 * q**2
    - 128 * p**2 * r**2
    + 144 * p * q**2 * r
    - 27 * q**4
    + 256 * r**3
)
assert sp.expand(
    sp.discriminant(reflection_equation, T) - reflection_discriminant
) == 0
reflection_saito = sp.Matrix(
    [
        [2 * p, -16 * r, -6 * q],
        [3 * q, 2 * p * q, 2 * p**2 - 8 * r],
        [4 * r, 3 * q**2 - 8 * p * r, p * q],
    ]
)
assert_saito(reflection_discriminant, (p, q, r), reflection_saito, 2)
reflection_lifts = [
    root_lift(
        reflection_equation,
        (p, q, r),
        tuple(reflection_saito[:, column]),
    )
    for column in range(3)
]
assert reflection_lifts == [T, -4 * T**3 - 4 * T * p - 3 * q, -2 * T**2 - p]


# The transverse last-target-coordinate direction always reconstructs the
# marked root with reciprocal derivative.  These identities are the exact
# bridge from the Saito calculation to the determinant ledgers.
assert sp.diff(cubic_equation, C) == -sp.Rational(1, 2)
assert sp.diff(cancellation_equation, R) == -1
assert sp.diff(quadratic_equation, y) == -1
assert sp.diff(weighted_equation, A) == C**2

print("PASS: four inverse discriminants and their reduced ledger factors")
print("PASS: full-target Saito bases for cubic, weighted quartic, cancellation quartic")
print("PASS: fixed-P Saito basis for the quadratic-gauge quartic")
print("PASS: logarithmic columns lift regularly to the marked-root incidences")
print("PASS: type-A3 reflection-discriminant control and regular root lifts")
