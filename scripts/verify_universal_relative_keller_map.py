#!/usr/bin/env python3
"""Exact audits for the universal relative and promoted absolute Keller maps."""

from __future__ import annotations

from itertools import permutations
from math import factorial

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
# all-degree vertical proof is the compact-chart factorization above.
for checked_degree in range(3, 6):
    mapping = compressed_map(checked_degree)
    determinant = sp.factor(
        sp.det(sp.Matrix(mapping[:3]).jacobian((x, y, z)))
    )
    assert determinant == 1


# In every rank, retaining all N-3 seed coefficients gives one polynomial
# self-map of A^N.  Its full Jacobian has an identity upper block, a zero
# upper-right block, and the same vertical Jacobian.  The compact identity
# above then proves determinant one for arbitrary N.
for checked_degree in range(3, 9):
    promoted = compressed_map(checked_degree)
    promoted_parameters = tuple(
        promoted[3][j] for j in range(4, checked_degree + 1)
    )
    promoted_source = (*promoted_parameters, x, y, z)
    promoted_mapping = (*promoted_parameters, *promoted[:3])
    assert len(promoted_source) == checked_degree
    assert len(promoted_mapping) == checked_degree

    promoted_jacobian = sp.Matrix(promoted_mapping).jacobian(promoted_source)
    parameter_count = checked_degree - 3
    assert (
        promoted_jacobian[:parameter_count, :parameter_count]
        == sp.eye(parameter_count)
    )
    assert (
        promoted_jacobian[:parameter_count, parameter_count:]
        == sp.zeros(parameter_count, 3)
    )
    assert (
        promoted_jacobian[parameter_count:, parameter_count:]
        == sp.Matrix(promoted[:3]).jacobian((x, y, z))
    )


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

    # Promote all N-3 normalized seed coefficients and compile this arbitrary
    # monic presentation into one target of the fixed A^N map.
    pi_value = jets[3] / jets[1]
    b_value = jets[2] / jets[1]
    c_value = -2 * jets[0] / jets[1]
    compiled_inverse = (
        S
        + b_value * S**2
        + pi_value * S**3
        + sum(
            (
                jets[j]
                * jets[1] ** (j - 1)
                / jets[3] ** j
            )
            * pi_value**j
            * S**j
            for j in range(4, degree + 1)
        )
        - c_value / 2
    )
    for exponent in range(degree + 1):
        assert sp.cancel(
            sp.Poly(compiled_inverse, S).nth(exponent)
            - jets[exponent] / jets[1]
        ) == 0


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
    promoted_inverse = S + b_value * S**2 + pi_value * S**3 - c_value / 2
    for j in range(4, degree + 1):
        u_value = hs[j] / pi_value**j
        assert sp.cancel(u_value * pi_value**j - hs[j]) == 0
        promoted_inverse += u_value * pi_value**j * S**j

    promoted_inverse = sp.Poly(promoted_inverse, S)
    assert promoted_inverse.degree() == degree
    for exponent in range(degree + 1):
        assert sp.cancel(
            promoted_inverse.nth(exponent) - hs[exponent]
        ) == 0

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


def compile_promoted_target(
    polynomial: sp.Expr,
    degree: int,
    translation: sp.Expr,
) -> tuple[tuple[sp.Expr, ...], sp.Poly]:
    """Compile a supplied presentation into one target of U_degree."""

    shifted = sp.Poly(
        sp.expand(polynomial.subs(T, translation + S)),
        S,
    )
    assert shifted.degree() == degree
    jets = {j: shifted.nth(j) for j in range(degree + 1)}
    assert jets[1] != 0
    assert jets[3] != 0
    normalized = sp.Poly(
        sp.cancel(shifted.as_expr() / jets[1]),
        S,
    )
    pi_value = sp.cancel(jets[3] / jets[1])
    target = tuple(
        sp.cancel(
            jets[j] * jets[1] ** (j - 1) / jets[3] ** j
        )
        for j in range(4, degree + 1)
    ) + (
        pi_value,
        sp.cancel(jets[2] / jets[1]),
        sp.cancel(-2 * jets[0] / jets[1]),
    )
    return target, normalized


def promoted_inverse_from_target(
    degree: int,
    target: tuple[sp.Expr, ...],
) -> sp.Poly:
    """Rebuild the selected inverse polynomial from an absolute target."""

    parameter_count = degree - 3
    parameters = target[:parameter_count]
    pi_value, b_value, c_value = target[parameter_count:]
    inverse = (
        S
        + b_value * S**2
        + pi_value * S**3
        + sum(
            parameters[j - 4] * pi_value**j * S**j
            for j in range(4, degree + 1)
        )
        - c_value / 2
    )
    return sp.Poly(sp.cancel(inverse), S)


collision_root_1, collision_root_2 = sp.symbols(
    "collision_root_1 collision_root_2"
)


def leading_exponents(basis: sp.GroebnerBasis) -> set[tuple[int, ...]]:
    """Return the leading exponent vectors of an exact Groebner basis."""

    return {
        tuple(polynomial.LM(order=basis.order).exponents)
        for polynomial in basis.polys
    }


def check_relative_collision_witness(
    polynomial: sp.Expr,
    degree: int,
) -> None:
    """Audit the diagonal/off-diagonal CRT for one separable root cover."""

    polynomial_over_q = sp.Poly(polynomial, T, domain=sp.QQ)
    assert polynomial_over_q.degree() == degree
    polynomial_1 = sp.expand(polynomial.subs(T, collision_root_1))
    polynomial_2 = sp.expand(polynomial.subs(T, collision_root_2))
    divided_difference = sp.Poly(
        sp.cancel(
            (polynomial_2 - polynomial_1)
            / (collision_root_2 - collision_root_1)
        ),
        collision_root_2,
        collision_root_1,
        domain=sp.QQ,
    ).as_expr()
    assert sp.expand(
        polynomial_2
        - polynomial_1
        - (collision_root_2 - collision_root_1) * divided_difference
    ) == 0

    derivative_inverse = sp.invert(
        sp.Poly(sp.diff(polynomial, T), T, domain=sp.QQ),
        polynomial_over_q,
    ).as_expr()
    inverse_remainder = sp.rem(
        sp.Poly(
            sp.expand(
                derivative_inverse * sp.diff(polynomial, T) - 1
            ),
            T,
            domain=sp.QQ,
        ),
        polynomial_over_q,
    )
    assert inverse_remainder.is_zero

    diagonal_idempotent = sp.expand(
        derivative_inverse.subs(T, collision_root_1)
        * divided_difference
    )
    collision_basis = sp.groebner(
        [polynomial_1, polynomial_2],
        collision_root_2,
        collision_root_1,
        order="lex",
        domain=sp.QQ,
    )
    diagonal_basis = sp.groebner(
        [polynomial_1, collision_root_2 - collision_root_1],
        collision_root_2,
        collision_root_1,
        order="lex",
        domain=sp.QQ,
    )
    off_diagonal_basis = sp.groebner(
        [polynomial_1, divided_difference],
        collision_root_2,
        collision_root_1,
        order="lex",
        domain=sp.QQ,
    )
    comaximal_basis = sp.groebner(
        [
            polynomial_1,
            collision_root_2 - collision_root_1,
            divided_difference,
        ],
        collision_root_2,
        collision_root_1,
        order="lex",
        domain=sp.QQ,
    )

    assert collision_basis.reduce(
        diagonal_idempotent**2 - diagonal_idempotent
    )[1] == 0
    assert diagonal_basis.reduce(diagonal_idempotent - 1)[1] == 0
    assert off_diagonal_basis.reduce(diagonal_idempotent)[1] == 0
    assert off_diagonal_basis.reduce(polynomial_2)[1] == 0
    assert len(comaximal_basis.polys) == 1
    assert comaximal_basis.polys[0].as_expr() == 1

    # With variable order root_2 > root_1, these rectangular standard
    # monomial sets have the displayed ranks N^2, N, and N(N-1).
    assert leading_exponents(collision_basis) == {
        (degree, 0),
        (0, degree),
    }
    assert leading_exponents(diagonal_basis) == {
        (1, 0),
        (0, degree),
    }
    assert leading_exponents(off_diagonal_basis) == {
        (degree - 1, 0),
        (0, degree),
    }


def check_ordered_configuration_action(degree: int) -> None:
    """Audit every ordered-distinct-root orbit of S_degree."""

    symmetric_group = tuple(permutations(range(degree)))
    assert len(symmetric_group) == factorial(degree)
    for tuple_length in range(1, degree + 1):
        base_tuple = tuple(range(tuple_length))
        orbit = {
            tuple(permutation[index] for index in base_tuple)
            for permutation in symmetric_group
        }
        stabilizer_size = sum(
            tuple(permutation[index] for index in base_tuple) == base_tuple
            for permutation in symmetric_group
        )
        expected_rank = factorial(degree) // factorial(
            degree - tuple_length
        )
        assert len(orbit) == expected_rank
        assert stabilizer_size == factorial(degree - tuple_length)
        assert len(orbit) * stabilizer_size == factorial(degree)


# Adversarial witness cards cover connected fields, the split algebra, and
# disconnected products.  Each target is pinned independently of the general
# coefficientwise identity above.
expected_connected_targets = {
    3: (
        sp.Rational(1, 3),
        sp.Integer(1),
        sp.Rational(2, 3),
    ),
    4: (
        sp.Rational(1, 4),
        sp.Integer(1),
        sp.Rational(3, 2),
        sp.Rational(1, 2),
    ),
    5: (
        sp.Rational(1, 16),
        sp.Rational(1, 160),
        sp.Integer(2),
        sp.Integer(2),
        sp.Rational(2, 5),
    ),
    6: (
        sp.Rational(81, 4000),
        sp.Rational(243, 100000),
        sp.Rational(243, 2000000),
        sp.Rational(10, 3),
        sp.Rational(5, 2),
        sp.Rational(1, 3),
    ),
}
for witness_degree, expected_target in expected_connected_targets.items():
    witness_polynomial = T**witness_degree - 2
    target, normalized = compile_promoted_target(
        witness_polynomial,
        witness_degree,
        sp.Integer(1),
    )
    assert target == expected_target
    assert promoted_inverse_from_target(witness_degree, target) == normalized
    assert sp.discriminant(witness_polynomial, T) != 0

    # Eisenstein at two certifies that Q[T]/(T^N-2) is a field.
    witness_coefficients = sp.Poly(witness_polynomial, T).all_coeffs()
    assert witness_coefficients[0] == 1
    assert all(coefficient % 2 == 0 for coefficient in witness_coefficients[1:])
    assert witness_coefficients[-1] % 4 != 0


# The all-rank arithmetic stress family P_N=T^N-T-1 has Galois group S_N
# over Q by Osada's theorem.  The checker does not re-prove that theorem; it
# certifies the closed-form target and inverse identity for the tested ranks.
for witness_degree in range(3, 13):
    witness_polynomial = T**witness_degree - T - 1
    target, normalized = compile_promoted_target(
        witness_polynomial,
        witness_degree,
        sp.Integer(1),
    )
    pi_value = sp.Rational(
        witness_degree * (witness_degree - 2),
        6,
    )
    expected_target = tuple(
        sp.cancel(
            sp.binomial(witness_degree, j)
            / (witness_degree - 1)
            / pi_value**j
        )
        for j in range(4, witness_degree + 1)
    ) + (
        pi_value,
        sp.Rational(witness_degree, 2),
        sp.Rational(2, witness_degree - 1),
    )
    assert target == expected_target
    assert promoted_inverse_from_target(witness_degree, target) == normalized
    assert sp.discriminant(witness_polynomial, T) != 0
    if witness_degree <= 8:
        check_relative_collision_witness(
            witness_polynomial,
            witness_degree,
        )
        check_ordered_configuration_action(witness_degree)


split_quartic = sp.prod(T - root for root in range(1, 5))
split_quartic_target, split_quartic_inverse = compile_promoted_target(
    split_quartic,
    4,
    sp.Integer(0),
)
assert split_quartic_target == (
    sp.Rational(-25, 2),
    sp.Rational(1, 5),
    sp.Rational(-7, 10),
    sp.Rational(24, 25),
)
assert promoted_inverse_from_target(4, split_quartic_target) == split_quartic_inverse
assert sp.discriminant(split_quartic, T) == 144


product_witnesses = {
    4: (
        (T**2 - 2) * (T**2 - 3),
        (
            sp.Rational(-27, 32),
            sp.Rational(-2, 3),
            sp.Rational(-1, 6),
            sp.Rational(2, 3),
        ),
    ),
    5: (
        (T**2 - 2) * (T**3 - 3),
        (
            sp.Rational(-1715, 4096),
            sp.Rational(2401, 32768),
            sp.Rational(-8, 7),
            sp.Rational(-1, 7),
            sp.Rational(4, 7),
        ),
    ),
    6: (
        (T**2 - 2) * (T**4 - 3),
        (
            sp.Rational(-26, 81),
            sp.Rational(8, 81),
            sp.Rational(-8, 729),
            sp.Rational(-3, 2),
            sp.Integer(0),
            sp.Rational(1, 2),
        ),
    ),
}
for witness_degree, (witness_polynomial, expected_target) in product_witnesses.items():
    target, normalized = compile_promoted_target(
        witness_polynomial,
        witness_degree,
        sp.Integer(1),
    )
    assert target == expected_target
    assert promoted_inverse_from_target(witness_degree, target) == normalized
    assert sp.discriminant(witness_polynomial, T) != 0


# Boundary attacks: these show why every open condition in the theorem is
# needed without threatening polynomiality of the ambient fixed map.
bad_linear_translation = sp.Poly((T**6 - 2).subs(T, S), S)
assert bad_linear_translation.nth(1) == 0

bad_cubic_translation = sp.Poly((T**4 + T).subs(T, S), S)
assert bad_cubic_translation.nth(1) == 1
assert bad_cubic_translation.nth(3) == 0
assert sp.discriminant(T**4 + T, T) != 0

quartic_connected_target = expected_connected_targets[4]
zero_pi_target = (
    quartic_connected_target[0],
    sp.Integer(0),
    quartic_connected_target[2],
    quartic_connected_target[3],
)
assert promoted_inverse_from_target(4, zero_pi_target).degree() < 4
zero_top_parameter_target = (
    sp.Integer(0),
    *quartic_connected_target[1:],
)
assert promoted_inverse_from_target(4, zero_top_parameter_target).degree() < 4

repeated_root_polynomial = (T - 1) ** 2 * (T**2 + 1)
repeated_target, repeated_inverse = compile_promoted_target(
    repeated_root_polynomial,
    4,
    sp.Integer(-3),
)
assert promoted_inverse_from_target(4, repeated_target) == repeated_inverse
assert sp.discriminant(repeated_root_polynomial, T) == 0


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


print(
    "universal relative and absolute Keller-map algebraic checks passed "
    "with exact ordered-collision regressions "
    "(S_N monodromy and atomicity are theorem-level imports)"
)
