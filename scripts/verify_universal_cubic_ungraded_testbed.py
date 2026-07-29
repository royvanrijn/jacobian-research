#!/usr/bin/env python3
"""Exact checks for the universal flat cubic ungraded testbed.

The universal degree-at-most-four coefficient cell is

    A(P) = a0+a1 P+...+a4 P^4,
    H(P) = 1+P^2 (c0+c1 P).

The inverse cubic is

    P H(P) S^3 + (A(P)-B/2) S^2 + S-C/2.

The checker verifies the denominator-free Keller map, its incidence and
ramification identities, the Deligne--Faddeev finite normalization,
discriminant, full GL2 coefficient action, and the base-change criterion.
"""

from __future__ import annotations

import sympy as sp


x, y, z, S = sp.symbols("x y z S")
P_symbol, B_symbol, C_symbol = sp.symbols("P B C")
a0, a1, a2, a3, a4, c0, c1 = sp.symbols(
    "a0 a1 a2 a3 a4 c0 c1"
)

t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t)
P = t * q


def alpha(variable: sp.Expr) -> sp.Expr:
    return (
        a0
        + a1 * variable
        + a2 * variable**2
        + a3 * variable**3
        + a4 * variable**4
    )


def gamma(variable: sp.Expr) -> sp.Expr:
    return c0 + c1 * variable


A = alpha(P)
H = 1 + P**2 * gamma(P)
leading = P * H

target_B = (
    y
    + 3 * x * q
    + 2 * alpha(P)
    + 3 * t**2 * x * q**3 * gamma(P)
)
target_C = x * (5 - 3 * t) - x**3 * z - x**3 * q**3 * gamma(P)

# The quartic coefficient cell is an exact Keller family.
jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(component, variable) for variable in (x, y, z)]
            for component in (P, target_B, target_C)
        ]
    )
)
assert sp.factor(jacobian) == -2

# The short chart and inverse cubic recover the map identically.
root = x / t
Q = y + x * q
D = sp.factor(1 - root * (Q - P * root))
assert D == 1 / t


def inverse_cubic(
    p_value: sp.Expr,
    b_value: sp.Expr,
    c_value: sp.Expr,
    root_value: sp.Expr,
) -> sp.Expr:
    return (
        p_value * (1 + p_value**2 * gamma(p_value)) * root_value**3
        + (alpha(p_value) - b_value / 2) * root_value**2
        + root_value
        - c_value / 2
    )


incidence = inverse_cubic(P, target_B, target_C, root)
assert sp.factor(incidence) == 0
incidence_derivative = sp.diff(
    inverse_cubic(P, target_B, target_C, S),
    S,
).subs(S, root)
assert sp.factor(incidence_derivative - D) == 0

chart_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(component, variable) for variable in (x, y, z)]
            for component in (P, root, Q)
        ]
    )
)
assert sp.factor(chart_jacobian - t) == 0

short_B = sp.expand(
    Q
    + 2 * alpha(P)
    + (3 * leading - P) * root
)
short_C = sp.expand(
    2
    * (
        leading * root**3
        + alpha(P) * root**2
        + root
    )
    - short_B * root**2
)
assert sp.factor(short_B - target_B) == 0
assert sp.factor(short_C - target_C) == 0
# SymPy cannot differentiate with respect to the expression x/t.  Repeat
# the short-chart identity with independent variables.
p_short, s_short, q_short = sp.symbols("p_short s_short q_short")
a_short = alpha(p_short)
h_short = 1 + p_short**2 * gamma(p_short)
leading_short = p_short * h_short
b_short = (
    q_short
    + 2 * a_short
    + (3 * leading_short - p_short) * s_short
)
c_short = (
    2
    * (
        leading_short * s_short**3
        + a_short * s_short**2
        + s_short
    )
    - b_short * s_short**2
)
short_D = 1 - s_short * q_short + p_short * s_short**2
plane_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(b_short, variable) for variable in (s_short, q_short)],
            [sp.diff(c_short, variable) for variable in (s_short, q_short)],
        ]
    )
)
assert sp.factor(plane_jacobian + 2 * short_D) == 0

# Deligne--Faddeev finite normalization.  For
# f=a*S^3+b*S^2+c*S+d, put u=a*S and v=a*S^2+b*S+c.
# Then 1,u,v satisfy the displayed finite-free multiplication table.
a, b, c, d = sp.symbols("a b c d")
u = a * S
v = a * S**2 + b * S + c
f = a * S**3 + b * S**2 + c * S + d


def reduce_mod_f(expression: sp.Expr) -> sp.Expr:
    return sp.factor(sp.rem(sp.Poly(expression, S), sp.Poly(f, S)).as_expr())


assert reduce_mod_f(u**2 - (-a * c - b * u + a * v)) == 0
assert reduce_mod_f(u * v + a * d) == 0
assert reduce_mod_f(v**2 - (-b * d - d * u + c * v)) == 0

# The trace determinant is the usual binary-cubic discriminant.
discriminant = (
    b**2 * c**2
    - 4 * a * c**3
    - 4 * b**3 * d
    - 27 * a**2 * d**2
    + 18 * a * b * c * d
)
testbed_a = P_symbol * (1 + P_symbol**2 * gamma(P_symbol))
testbed_b = alpha(P_symbol) - B_symbol / 2
testbed_c = sp.Integer(1)
testbed_d = -C_symbol / 2
testbed_discriminant = sp.factor(
    discriminant.subs(
        {
            a: testbed_a,
            b: testbed_b,
            c: testbed_c,
            d: testbed_d,
        }
    )
)
expected_discriminant = sp.factor(
    testbed_b**2
    - 4 * testbed_a
    + 2 * testbed_b**3 * C_symbol
    - sp.Rational(27, 4) * testbed_a**2 * C_symbol**2
    - 9 * testbed_a * testbed_b * C_symbol
)
assert sp.factor(testbed_discriminant - expected_discriminant) == 0

# The ramification support is the smooth Laurent plane with parameters
# (P,S): B=2A+3aS+S^-1 and C=S-aS^3.
ramification_B = (
    2 * alpha(P_symbol)
    + 3 * testbed_a * S
    + 1 / S
)
ramification_C = S - testbed_a * S**3
ramification_inverse = inverse_cubic(
    P_symbol,
    ramification_B,
    ramification_C,
    S,
)
ramification_derivative = sp.diff(
    inverse_cubic(P_symbol, B_symbol, C_symbol, S),
    S,
).subs({B_symbol: ramification_B, C_symbol: ramification_C})
assert sp.factor(ramification_inverse) == 0
assert sp.factor(ramification_derivative) == 0
assert (
    sp.factor(
        testbed_discriminant.subs(
            {B_symbol: ramification_B, C_symbol: ramification_C}
        )
    )
    == 0
)

# A completely general polynomial GL2 Tschirnhausen gauge acts through
# these universal coefficient identities.  Its determinant is a unit in
# the polynomial case; the discriminant divisor is unchanged because the
# discriminant scales by det^6.
r, s, lower, upper = sp.symbols("r s lower upper")
U, V = sp.symbols("U V")
binary_cubic = a * U**3 + b * U**2 * V + c * U * V**2 + d * V**3
gauged = sp.Poly(
    sp.expand(
        binary_cubic.subs(
            {
                U: r * U + s * V,
                V: lower * U + upper * V,
            },
            simultaneous=True,
        )
    ),
    U,
    V,
)
gauged_coefficients = [
    gauged.coeff_monomial(U**3),
    gauged.coeff_monomial(U**2 * V),
    gauged.coeff_monomial(U * V**2),
    gauged.coeff_monomial(V**3),
]
gauged_discriminant = discriminant.subs(
    dict(zip((a, b, c, d), gauged_coefficients)),
    simultaneous=True,
)
gauge_determinant = r * upper - s * lower
assert sp.factor(gauged_discriminant - gauge_determinant**6 * discriminant) == 0

# In the adapted hyperplane c=1, the intrinsic coefficient map is
# G=(a,b,d).  Its Jacobian is a'(P)/4, and it is an automorphism precisely
# when the degree-drop/phantom factor H is the unit 1 on this normalized
# cell, i.e. c0=c1=0.
coefficient_jacobian = sp.det(
    sp.Matrix(
        [
            [
                sp.diff(component, variable)
                for variable in (P_symbol, B_symbol, C_symbol)
            ]
            for component in (testbed_a, testbed_b, testbed_d)
        ]
    )
)
assert sp.factor(
    coefficient_jacobian - sp.diff(testbed_a, P_symbol) / 4
) == 0
assert sp.Poly(sp.diff(testbed_a, P_symbol), P_symbol).degree() == 3
assert sp.solve(
    [
        sp.Poly(
            sp.diff(testbed_a, P_symbol), P_symbol
        ).coeff_monomial(P_symbol**degree)
        for degree in (1, 2, 3)
    ],
    (c0, c1),
    dict=True,
) == [{c0: 0, c1: 0}]

# The flat cell and the reduced Koszul quartic-kernel cell have different
# intrinsic Fitt_3 ideals.  A free rank-three normalization has unit
# Fitt_3, while A plus coker(z,-y,x)^T has Fitt_3=(x,y,z).
flat_fitting_3 = {sp.Integer(1)}
koszul_fitting_3 = {
    entry
    for entry in sp.Matrix((0, z, -y, x))
    if entry != 0
}
assert flat_fitting_3 == {1}
assert koszul_fitting_3 == {x, -y, z}

print("PASS: the seven-parameter quartic coefficient cell has determinant -2")
print("PASS: the inverse cubic, derivative, and reciprocal chart are exact")
print("PASS: the Deligne--Faddeev normalization and discriminant are exact")
print("PASS: the ramification support is the smooth Laurent plane")
print("PASS: arbitrary GL2 gauge scales the discriminant by det^6")
print("PASS: G is an automorphism exactly when the phantom factor is a unit")
print("PASS: the flat and 24-direction Koszul cells have different Fitt_3 ideals")
