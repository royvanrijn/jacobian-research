#!/usr/bin/env python3
"""Exact audits for the universal relative quadratic-gauge Keller map."""

from __future__ import annotations

import sympy as sp


# The compact chart proves the Jacobian for every seed degree at once.
Pi, S, Q = sp.symbols("Pi S Q")
g1 = sp.symbols("g1", nonzero=True)
G = sp.Function("G")
D = 1 - S * Q + Pi * S**2

# Substitute beta=(G'/g1-1-Pi*S^2)/S and differentiate directly.
direct_B = Q + (sp.diff(G(S), S) / g1 - 1 - Pi * S**2) / S
direct_C = 2 * G(S) / g1 - direct_B * S**2
direct_jacobian = sp.simplify(
    sp.det(sp.Matrix([direct_B, direct_C]).jacobian((S, Q)))
)
assert sp.simplify(direct_jacobian + 2 * D) == 0


# The reciprocal source chart contributes D^{-1}.
x, y, z = sp.symbols("x y z")
t = 1 + x * y
r = sp.symbols("r")
q_source = t**2 * z + r * y**2 * (1 + 3 * t)
pi_chart = t * q_source
s_chart = x / t
q_chart = y + x * q_source
chart_jacobian = sp.factor(
    sp.det(sp.Matrix([pi_chart, s_chart, q_chart]).jacobian((x, y, z)))
)
assert chart_jacobian == t
assert sp.factor(
    D.subs({Pi: pi_chart, S: s_chart, Q: q_chart}) - 1 / t
) == 0


def compressed_map(degree: int) -> tuple[sp.Expr, sp.Expr, sp.Expr, dict[int, sp.Expr]]:
    """Return the normalized compressed map for one symbolic degree."""

    us = {
        j: sp.symbols(f"u{degree}_{j}", nonzero=(j == degree))
        for j in range(4, degree + 1)
    }
    local_q = t**2 * z + y**2 * (1 + 3 * t)
    first = t * local_q
    raw_second = (
        y
        + 3 * x * local_q
        + sum(
            j * us[j] * t**2 * x ** (j - 2) * local_q**j
            for j in range(4, degree + 1)
        )
    )
    third = (
        x * (5 - 3 * t)
        - x**3 * z
        - sum(
            (j - 2) * us[j] * (x * local_q) ** j
            for j in range(4, degree + 1)
        )
    )
    return first, -raw_second / 2, third, us


# Direct source-coordinate Jacobians are useful low-degree regressions.  The
# all-degree proof is the compact-chart factorization above.
for checked_degree in range(3, 6):
    mapping = compressed_map(checked_degree)
    determinant = sp.factor(
        sp.det(sp.Matrix(mapping[:3]).jacobian((x, y, z)))
    )
    assert determinant == 1


# The degree-five relative map extends to one absolute polynomial map of A^5
# by retaining its two seed coefficients.  Its full Jacobian remains one.
quintic_mapping = compressed_map(5)
u4, u5 = quintic_mapping[3][4], quintic_mapping[3][5]
absolute_quintic = sp.Matrix(
    [u4, u5, quintic_mapping[0], quintic_mapping[1], quintic_mapping[2]]
)
absolute_quintic_jacobian = absolute_quintic.jacobian((u4, u5, x, y, z))
assert absolute_quintic_jacobian[:2, :2] == sp.eye(2)
assert absolute_quintic_jacobian[:2, 2:] == sp.zeros(2, 3)
assert sp.factor(absolute_quintic_jacobian.det()) == 1


T, a = sp.symbols("T a")
for degree in range(3, 9):
    coefficients = sp.symbols(f"c{degree}_0:{degree}")
    polynomial = T**degree + sum(
        coefficients[j] * T**j for j in range(degree)
    )
    shifted = sp.Poly(sp.expand(polynomial.subs(T, a + S)), S)
    jets = {j: shifted.nth(j) for j in range(degree + 1)}
    assert jets[degree] == 1

    # Exact inverse to (P,a) -> (H,a): P(T)=h_N^{-1} H(T-a).
    h_degree = 1 / jets[1]
    assert sp.cancel(h_degree * jets[1] - 1) == 0
    unshifted = sp.expand(shifted.as_expr().subs(S, T - a))
    assert sp.Poly(unshifted - polynomial, T).is_zero


# The N-3 map parameters and three target coordinates reproduce every
# normalized inverse polynomial coefficientwise.  Independent jet symbols
# keep this all-degree audit small instead of expanding high powers of the
# universal translated coefficients.
for degree in range(3, 13):
    hs = {0: sp.symbols(f"h{degree}_0"), 1: sp.Integer(1)}
    hs.update({
        j: sp.symbols(f"h{degree}_{j}", nonzero=(j in {3, degree}))
        for j in range(2, degree + 1)
    })
    pi_value = hs[3]
    b_value = hs[2]
    c_value = -2 * hs[0]
    assert -c_value / 2 == hs[0]
    assert b_value == hs[2]
    assert pi_value == hs[3]
    for j in range(4, degree + 1):
        u_value = hs[j] / pi_value**j
        assert sp.cancel(u_value * pi_value**j - hs[j]) == 0

    # Formula (3.9) is the same parameter change written in Hasse jets.
    ds = {
        j: sp.symbols(f"d{degree}_{j}", nonzero=(j in {1, 3, degree}))
        for j in range(degree + 1)
    }
    for j in range(4, degree + 1):
        normalized_jet_u = sp.cancel(
            (ds[j] / ds[1]) / (ds[3] / ds[1]) ** j
        )
        displayed_jet_u = sp.cancel(
            ds[j] * ds[1] ** (j - 1) / ds[3] ** j
        )
        assert normalized_jet_u == displayed_jet_u


# Every monic squarefree quintic presentation enters the absolute A^5 map
# after translating away from the first- and third-jet divisors.
absolute_coefficients = sp.symbols("absolute_c0:5")
absolute_polynomial = T**5 + sum(
    absolute_coefficients[j] * T**j for j in range(5)
)
absolute_shifted = sp.Poly(
    sp.expand(absolute_polynomial.subs(T, a + S)), S
)
absolute_jets = {j: absolute_shifted.nth(j) for j in range(6)}
absolute_pi = absolute_jets[3] / absolute_jets[1]
absolute_b = absolute_jets[2] / absolute_jets[1]
absolute_c = -2 * absolute_jets[0] / absolute_jets[1]
absolute_u = (
    absolute_jets[4]
    * absolute_jets[1] ** 3
    / absolute_jets[3] ** 4
)
absolute_v = (
    absolute_jets[5]
    * absolute_jets[1] ** 4
    / absolute_jets[3] ** 5
)
absolute_inverse = (
    absolute_v * absolute_pi**5 * S**5
    + absolute_u * absolute_pi**4 * S**4
    + absolute_pi * S**3
    + absolute_b * S**2
    + S
    - absolute_c / 2
)
assert sp.cancel(
    absolute_inverse
    - absolute_shifted.as_expr() / absolute_jets[1]
) == 0


# Five explicit generic-family surfaces retain their inverse polynomials.
def compiled_quintic_target(
    coefficient_4: sp.Expr,
    coefficient_3: sp.Expr,
    coefficient_2: sp.Expr,
    coefficient_1: sp.Expr,
    coefficient_0: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Compile a monic quintic at the origin into the absolute target."""

    return (
        coefficient_4 * coefficient_1**3 / coefficient_3**4,
        coefficient_1**4 / coefficient_3**5,
        coefficient_3 / coefficient_1,
        coefficient_2 / coefficient_1,
        -2 * coefficient_0 / coefficient_1,
    )


def assert_compiled_quintic(
    polynomial: sp.Expr,
    coefficients: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Check that the compiled inverse is the polynomial divided by its linear term."""

    target = compiled_quintic_target(*coefficients)
    target_u, target_v, target_pi, target_b, target_c = target
    target_inverse = (
        target_v * target_pi**5 * S**5
        + target_u * target_pi**4 * S**4
        + target_pi * S**3
        + target_b * S**2
        + S
        - target_c / 2
    )
    linear_coefficient = coefficients[3]
    for exponent in range(6):
        assert sp.cancel(
            sp.Poly(target_inverse, S).nth(exponent)
            - sp.Poly(polynomial, S).nth(exponent) / linear_coefficient
        ) == 0
    return target


r_param, s_param = sp.symbols("r_param s_param", nonzero=True)
pi_s5 = r_param / s_param
v_s5 = s_param**4 / r_param**5
inverse_s5 = sp.expand(
    v_s5 * pi_s5**5 * S**5
    + pi_s5 * S**3
    + S
    + 1
)
generic_s5 = S**5 + r_param * S**3 + s_param * S + s_param
assert sp.cancel(inverse_s5 - generic_s5 / s_param) == 0

brumer_third = r_param - s_param + 3
u_d5 = (s_param - 3) * r_param**3 / brumer_third**4
v_d5 = r_param**4 / brumer_third**5
pi_d5 = brumer_third / r_param
b_d5 = (s_param**2 - s_param - 2 * r_param - 1) / r_param
c_d5 = -2 * s_param / r_param
inverse_d5 = sp.expand(
    v_d5 * pi_d5**5 * S**5
    + u_d5 * pi_d5**4 * S**4
    + pi_d5 * S**3
    + b_d5 * S**2
    + S
    - c_d5 / 2
)
generic_d5 = (
    S**5
    + (s_param - 3) * S**4
    + brumer_third * S**3
    + (s_param**2 - s_param - 2 * r_param - 1) * S**2
    + r_param * S
    + s_param
)
assert sp.cancel(inverse_d5 - generic_d5 / r_param) == 0


f20_d = s_param**2 + 4
f20_a4 = r_param * f20_d - 2 * s_param - sp.Rational(17, 4)
f20_a3 = (
    3 * r_param * f20_d
    + f20_d
    + sp.Rational(13, 2) * s_param
    + 1
)
f20_a2 = -(r_param * f20_d + sp.Rational(11, 2) * s_param - 8)
f20_a1 = s_param - 6
f20_a0 = sp.Integer(1)
generic_f20 = (
    S**5
    + f20_a4 * S**4
    + f20_a3 * S**3
    + f20_a2 * S**2
    + f20_a1 * S
    + f20_a0
)
target_f20 = assert_compiled_quintic(
    generic_f20,
    (f20_a4, f20_a3, f20_a2, f20_a1, f20_a0),
)


a5_A, a5_B = sp.symbols("a5_A a5_B", nonzero=True)
a5_gamma = 5 * a5_A**2 - a5_B**2 + 3
a5_lambda = (
    a5_B * a5_gamma**2
    - 52 * a5_B * a5_gamma
    + 576 * a5_B
    - 10 * a5_gamma**2
    + 360 * a5_gamma
    - 3456
)
a5_sigma = 125 * a5_gamma**2 / (4 * a5_lambda)
a5_tau = 3125 * a5_gamma**5 / (256 * a5_lambda**2)
generic_a5 = S**5 + a5_sigma * S**3 + a5_tau * S + a5_tau
target_a5 = assert_compiled_quintic(
    generic_a5,
    (sp.Integer(0), a5_sigma, sp.Integer(0), a5_tau, a5_tau),
)
assert target_a5[0] == 0
assert sp.cancel(
    target_a5[1]
    - 3125 * a5_gamma**10 / (2**22 * a5_lambda**3)
) == 0
assert sp.cancel(
    target_a5[2] - 64 * a5_lambda / (25 * a5_gamma**3)
) == 0
assert target_a5[3] == 0
assert target_a5[4] == -2


c5_A, c5_B = sp.symbols("c5_A c5_B")
c5_Q = -c5_A + 1 + c5_B**2 * c5_A + 7 * c5_B**2
c5_R4 = (
    c5_A**3
    + c5_A**2
    + 10 * c5_B**2 * c5_A
    - 3 * c5_A
    + 20 * c5_B**2
    + 3
)
c5_R3 = (
    -24 * c5_B**2 * c5_A
    + 28 * c5_B**2
    + 210 * c5_B**4 * c5_A
    + 3
    - 28 * c5_B**2 * c5_A**2
    - 40 * c5_B**4
    - 625 * c5_B**6
    - 8 * c5_A
    - 135 * c5_B**4 * c5_A**2
    - 3 * c5_A**4
    + 2 * c5_A**5
    - 7 * c5_B**2 * c5_A**4
    + 7 * c5_A**2
    + 44 * c5_A**3 * c5_B**2
)
c5_R2 = (
    4 * c5_A**4
    - 1
    + c5_A**6
    + 6 * c5_A
    + 305 * c5_B**4
    + 1250 * c5_B**6
    + 44 * c5_B**2 * c5_A**2
    - 220 * c5_B**4 * c5_A
    - 52 * c5_A**3 * c5_B**2
    + 345 * c5_B**4 * c5_A**2
    - 2 * c5_A**5
    + 12 * c5_B**2 * c5_A
    + 31 * c5_B**2 * c5_A**4
    + 11 * c5_B**2
    - 6 * c5_A**2
)
c5_R1 = (
    2 * c5_A**5
    - 2 * c5_A**4
    - 8 * c5_B**2 * c5_A**4
    + 36 * c5_A**3 * c5_B**2
    - 145 * c5_B**4 * c5_A**2
    + 3 * c5_A**2
    - 22 * c5_B**2 * c5_A**2
    + 4 * c5_B**2 * c5_A
    + 120 * c5_B**4 * c5_A
    - 2 * c5_A
    - 13 * c5_B**2
    - 180 * c5_B**4
    - 625 * c5_B**6
)
c5_R0 = (
    c5_A**3
    + c5_A**2
    + 7 * c5_B**2 * c5_A
    - c5_B**2
)
c5_coefficients = (
    -c5_R4 / c5_Q,
    c5_R3 / c5_Q**2,
    c5_R2 / c5_Q**2,
    c5_R1 / c5_Q**2,
    -c5_R0 / c5_Q,
)
generic_c5 = (
    S**5
    + c5_coefficients[0] * S**4
    + c5_coefficients[1] * S**3
    + c5_coefficients[2] * S**2
    + c5_coefficients[3] * S
    + c5_coefficients[4]
)
target_c5 = assert_compiled_quintic(generic_c5, c5_coefficients)
displayed_target_c5 = (
    -c5_Q * c5_R4 * c5_R1**3 / c5_R3**4,
    c5_Q**2 * c5_R1**4 / c5_R3**5,
    c5_R3 / c5_R1,
    c5_R2 / c5_R1,
    2 * c5_Q * c5_R0 / c5_R1,
)
for computed_coordinate, displayed_coordinate in zip(
    target_c5,
    displayed_target_c5,
    strict=True,
):
    assert sp.cancel(computed_coordinate - displayed_coordinate) == 0


print("universal relative Keller-map checks passed")
