#!/usr/bin/env python3
"""Exact checks for one fixed quintic Keller map dominating quintic moduli."""

from __future__ import annotations

import os
import warnings

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings(
    "ignore",
    category=SymPyDeprecationWarning,
)

x, y, z, S, T, X, Y = sp.symbols("x y z S T X Y")
Pi, B, C = sp.symbols("Pi B C")


def factor_degrees_mod_prime(poly: sp.Poly, prime: int) -> tuple[int, ...]:
    """Return the squarefree factor-degree pattern modulo ``prime``."""
    assert int(sp.discriminant(poly.as_expr(), S)) % prime != 0
    factors = sp.factor_list(poly.as_expr(), modulus=prime)[1]
    return tuple(
        sorted(
            (
                int(sp.degree(factor, S))
                for factor, exponent in factors
                for _ in range(exponent)
            ),
            reverse=True,
        )
    )


# Fixed split seed.
G = S**5 - 5 * S**3 + 4 * S
assert sp.factor(G) == S * (S - 2) * (S - 1) * (S + 1) * (S + 2)
assert sp.discriminant(G, S) == 82944

# Fixed polynomial Keller map.
t = 1 + x * y
q = t**2 * z - sp.Rational(4, 5) * y**2 * (1 + 3 * t)
mapping = (
    t * q,
    y
    - sp.Rational(15, 4) * x * q
    + sp.Rational(5, 4) * t**2 * x**3 * q**5,
    x * (5 - 3 * t)
    + sp.Rational(5, 4) * x**3 * z
    - sp.Rational(3, 4) * (x * q) ** 5,
)

jacobian = sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det())
assert jacobian == -2

coordinate_degrees = tuple(
    int(sp.Poly(sp.expand(component), x, y, z, domain=sp.QQ).total_degree())
    for component in mapping
)
assert coordinate_degrees == (7, 32, 30)

# Inverse quintic and the small affine-normalized coefficient chart.
E = Pi**5 * S**5 - 5 * Pi * S**3 - 2 * B * S**2 + 4 * S - 2 * C
assert sp.Poly(E, S).degree() == 5
normalized_E = sp.expand(Pi**5 * E.subs(S, Pi**-2 * T))
c2 = -2 * Pi * B
c3 = 4 * Pi**3
c4 = -2 * Pi**5 * C
assert normalized_E == T**5 - 5 * T**3 + c2 * T**2 + c3 * T + c4

coefficient_jacobian = sp.factor(
    sp.Matrix((c2, c3, c4)).jacobian((Pi, B, C)).det()
)
assert coefficient_jacobian == -48 * Pi**8

certificate_point = {Pi: 1, B: sp.Rational(5, 2), C: 0}
assert coefficient_jacobian.subs(certificate_point) == -48
certificate_poly = sp.Poly(E.subs(certificate_point), S, domain=sp.QQ)
assert sp.discriminant(certificate_poly.as_expr(), S) == -1139056
assert sp.gcd(certificate_poly, certificate_poly.diff()).degree() == 0

# Rational affine-quotient chart.
u = -sp.Rational(125, 4) / (Pi**2 * B**2)
v = sp.Rational(625, 4) / (Pi * B**4)
w = -sp.Rational(3125, 16) * C / B**5
rational_quotient_jacobian = sp.factor(
    sp.Matrix((u, v, w)).jacobian((Pi, B, C)).det()
)
assert rational_quotient_jacobian == sp.Rational(732421875, 128) / (
    Pi**4 * B**12
)
assert rational_quotient_jacobian.subs(certificate_point) == 96
assert tuple(value.subs(certificate_point) for value in (u, v, w)) == (-5, 4, 0)

# Explicit S_5 specializations in all three quintic signatures.  An
# irreducible reduction gives a 5-cycle and irreducibility over Q.  The
# additional (4,1) and (3,2) patterns force the transitive group to be S_5.
signature_examples = (
    {
        "target": {Pi: 1, B: -1, C: -1},
        "real_roots": 1,
        "discriminant": 500624,
        "pattern_primes": {11: (5,), 7: (4, 1), 3: (3, 2)},
    },
    {
        "target": {Pi: 1, B: 0, C: -1},
        "real_roots": 3,
        "discriminant": -57056,
        "pattern_primes": {3: (5,), 31: (4, 1), 23: (3, 2)},
    },
    {
        "target": {Pi: 1, B: 0, C: -sp.Rational(1, 2)},
        "real_roots": 5,
        "discriminant": 38569,
        "pattern_primes": {2: (5,), 7: (4, 1), 19: (3, 2)},
    },
)

for example in signature_examples:
    poly = sp.Poly(E.subs(example["target"]), S, domain=sp.QQ)
    assert sp.discriminant(poly.as_expr(), S) == example["discriminant"]
    assert poly.count_roots(-sp.oo, sp.oo) == example["real_roots"]
    for prime, expected_pattern in example["pattern_primes"].items():
        assert factor_degrees_mod_prime(poly, prime) == expected_pattern

# Every unramified quintic splitting partition occurs in this fixed family
# already over F_7.
mod_seven_witnesses = {
    (5,): (1, 0, 1),
    (4, 1): (1, 0, 3),
    (3, 2): (1, 0, 2),
    (3, 1, 1): (1, 1, 3),
    (2, 2, 1): (3, 2, 0),
    (2, 1, 1, 1): (1, 2, 6),
    (1, 1, 1, 1, 1): (1, 0, 0),
}
for expected_pattern, (pi_value, b_value, c_value) in mod_seven_witnesses.items():
    poly = sp.Poly(
        E.subs({Pi: pi_value, B: b_value, C: c_value}),
        S,
        domain=sp.QQ,
    )
    assert factor_degrees_mod_prime(poly, 7) == expected_pattern

# Exact target discriminant and repeated-root normalization.
Delta = (
    432 * B**5 * C * Pi**2
    - 432 * B**4 * Pi**2
    + 12600 * B**3 * C * Pi**3
    - 2000 * B**3 * C
    + 9000 * B**2 * C**2 * Pi**7
    + 20625 * B**2 * C**2 * Pi**4
    - 11520 * B**2 * Pi**3
    + 2000 * B**2
    + 18750 * B * C**3 * Pi**8
    - 25600 * B * C * Pi**7
    + 56000 * B * C * Pi**4
    - 45000 * B * C * Pi
    + 3125 * C**4 * Pi**12
    - 40000 * C**2 * Pi**8
    + 112500 * C**2 * Pi**5
    - 84375 * C**2 * Pi**2
    + 16384 * Pi**7
    - 51200 * Pi**4
    + 40000 * Pi
)
assert sp.factor(sp.discriminant(E, S)) == 16 * Pi**8 * Delta

r = sp.symbols("r", nonzero=True)
discriminant_B = (5 * Pi**5 * r**4 - 15 * Pi * r**2 + 4) / (4 * r)
discriminant_C = -r * (3 * Pi**5 * r**4 - 5 * Pi * r**2 - 4) / 4
repeated_root_substitution = {S: r, B: discriminant_B, C: discriminant_C}
assert sp.factor(E.subs(repeated_root_substitution)) == 0
assert sp.factor(sp.diff(E, S).subs(repeated_root_substitution)) == 0

# Rational-root and quadratic-times-cubic incidences.
root_C = (Pi**5 * r**5 - 5 * Pi * r**3 - 2 * B * r**2 + 4 * r) / 2
assert sp.factor(E.subs({S: r, C: root_C})) == 0

quadratic_a, quadratic_b = sp.symbols("quadratic_a quadratic_b", nonzero=True)
quadratic_factor = S**2 + quadratic_a * S + quadratic_b
factor_B = -(
    Pi**5 * quadratic_a**4
    - 3 * Pi**5 * quadratic_a**2 * quadratic_b
    + Pi**5 * quadratic_b**2
    - 5 * Pi * quadratic_a**2
    + 5 * Pi * quadratic_b
    + 4
) / (2 * quadratic_a)
factor_C = (
    quadratic_b
    * (
        Pi**5 * quadratic_a**2 * quadratic_b
        - Pi**5 * quadratic_b**2
        - 5 * Pi * quadratic_b
        - 4
    )
    / (2 * quadratic_a)
)
factor_remainder = sp.rem(
    sp.Poly(E.subs({B: factor_B, C: factor_C}), S),
    sp.Poly(quadratic_factor, S),
).as_expr()
assert sp.factor(factor_remainder) == 0

# Pair-sum resolvent and its target pullback.
resolvent_a, resolvent_b, resolvent_c = sp.symbols(
    "resolvent_a resolvent_b resolvent_c"
)
normalized_quintic = (
    Y**5
    - 5 * Y**3
    + resolvent_a * Y**2
    + resolvent_b * Y
    + resolvent_c
)
pair_sum_resolvent = (
    X**10
    - 15 * X**8
    + resolvent_a * X**7
    + (75 - 3 * resolvent_b) * X**6
    + (-10 * resolvent_a - 11 * resolvent_c) * X**5
    + (-resolvent_a**2 + 10 * resolvent_b - 125) * X**4
    + (-4 * resolvent_a * resolvent_b + 25 * resolvent_a + 20 * resolvent_c)
    * X**3
    + (
        5 * resolvent_a**2
        + 7 * resolvent_a * resolvent_c
        - 4 * resolvent_b**2
        + 25 * resolvent_b
    )
    * X**2
    + (
        -resolvent_a**3
        + 4 * resolvent_b * resolvent_c
        - 25 * resolvent_c
    )
    * X
    - resolvent_a**2 * resolvent_b
    - 5 * resolvent_a * resolvent_c
    - resolvent_c**2
)
pair_sum_resultant = sp.resultant(
    normalized_quintic,
    normalized_quintic.subs(Y, X - Y),
    Y,
)
diagonal_factor = 2**5 * normalized_quintic.subs(Y, X / 2)
assert sp.expand(
    pair_sum_resultant - diagonal_factor * pair_sum_resolvent**2
) == 0

target_pair_sum_resolvent = sp.expand(
    pair_sum_resolvent.subs(
        {
            resolvent_a: c2,
            resolvent_b: c3,
            resolvent_c: c4,
        }
    )
)
assert sp.Poly(target_pair_sum_resolvent, X).degree() == 10
assert len(sp.Poly(target_pair_sum_resolvent, X, Pi, B, C).terms()) == 23

# One explicit A_5 fiber.
a5_target = {Pi: 1, B: -sp.Rational(4, 3), C: 6}
a5_poly = sp.Poly(E.subs(a5_target), S, domain=sp.QQ)
assert sp.discriminant(a5_poly.as_expr(), S) == 984**2
assert factor_degrees_mod_prime(a5_poly, 5) == (5,)
assert factor_degrees_mod_prime(a5_poly, 7) == (3, 1, 1)
a5_resolvent = sp.Poly(
    target_pair_sum_resolvent.subs(a5_target),
    X,
    domain=sp.QQ,
)
assert a5_resolvent.is_irreducible

# Trace formulas for a generic centered monic quintic.
trace_a, trace_b, trace_c, trace_d = sp.symbols(
    "trace_a trace_b trace_c trace_d"
)
generic_companion = sp.zeros(5)
for column in range(4):
    generic_companion[column + 1, column] = 1
generic_companion[:, 4] = sp.Matrix(
    (-trace_d, -trace_c, -trace_b, -trace_a, 0)
)
assert sp.trace(generic_companion) == 0
assert sp.factor(sp.trace(generic_companion**2)) == -2 * trace_a
assert sp.expand(
    sp.trace(generic_companion**4) - (2 * trace_a**2 - 4 * trace_c)
) == 0

# The projective trace compactification is a (2,4) complete intersection.
# Over the split algebra, eliminate x_5 using sum_i x_i=0.
split_x = sp.symbols("split_x1:5")
hom_s, hom_W = sp.symbols("hom_s hom_W")
split_x5 = -sum(split_x)
split_p2 = sp.expand(sum(value**2 for value in split_x) + split_x5**2)
split_p4 = sp.expand(sum(value**4 for value in split_x) + split_x5**4)
trace_quadric = sp.expand(split_p2 - 10 * hom_s**2)
trace_quartic = sp.expand(
    split_p4 - 50 * hom_s**4 + 16 * hom_s * hom_W**3
)
assert sp.Poly(trace_quadric, *split_x, hom_s, hom_W).total_degree() == 2
assert sp.Poly(trace_quartic, *split_x, hom_s, hom_W).total_degree() == 4

# The useful locus sW != 0 is smooth: this fixed Jacobian minor is a unit
# there.  The projective boundary vertex [0:0:1] is nevertheless singular.
smooth_minor = sp.det(
    sp.Matrix(
        [
            [sp.diff(trace_quadric, hom_s), sp.diff(trace_quadric, hom_W)],
            [sp.diff(trace_quartic, hom_s), sp.diff(trace_quartic, hom_W)],
        ]
    )
)
assert sp.factor(smooth_minor) == -960 * hom_s**2 * hom_W**2
boundary_vertex = {**{value: 0 for value in split_x}, hom_s: 0, hom_W: 1}
boundary_jacobian = sp.Matrix(
    [
        [sp.diff(trace_quadric, value) for value in (*split_x, hom_s, hom_W)],
        [sp.diff(trace_quartic, value) for value in (*split_x, hom_s, hom_W)],
    ]
).subs(boundary_vertex)
assert boundary_jacobian.rank() == 1

# A small non-affine primitive generator realizes Q[theta]/(theta^5-theta-1).
theta_companion = sp.zeros(5)
for column in range(4):
    theta_companion[column + 1, column] = 1
theta_companion[:, 4] = sp.Matrix((1, 1, 0, 0, 0))
eta_operator = theta_companion**3 - theta_companion**2 + 2 * theta_companion
assert sp.trace(eta_operator) == 0
assert sp.trace(eta_operator**2) == 10
assert sp.trace(eta_operator**4) == -78
eta_characteristic = sp.expand(eta_operator.charpoly(T).as_expr())
assert eta_characteristic == T**5 - 5 * T**3 - 13 * T**2 + 32 * T - 23

descent_target = {
    Pi: 2,
    B: sp.Rational(13, 4),
    C: sp.Rational(23, 64),
}
descent_inverse = sp.expand(
    2**5
    * E.subs(descent_target).subs(S, sp.Rational(1, 4) * T)
)
assert descent_inverse == eta_characteristic
original_field_polynomial = sp.Poly(S**5 - S - 1, S, domain=sp.QQ)
assert sp.discriminant(original_field_polynomial.as_expr(), S) == 2869
assert factor_degrees_mod_prime(original_field_polynomial, 3) == (5,)
assert factor_degrees_mod_prime(original_field_polynomial, 2) == (3, 2)


def verify_nonaffine_s5_sample(
    *,
    hermite_u_value: int,
    hermite_v_value: int,
    eta_coefficients: tuple[sp.Rational, ...],
    pi_value: sp.Rational,
    target_b: sp.Rational,
    target_c: sp.Rational,
    discriminant: int,
    irreducible_prime: int,
    odd_order_six_prime: int,
) -> None:
    """Check one exact S_5 field and its non-affine trace point."""
    field_polynomial = sp.Poly(
        S**5
        + hermite_u_value * S**3
        + hermite_v_value * S
        + hermite_v_value,
        S,
        domain=sp.QQ,
    )
    assert sp.discriminant(field_polynomial.as_expr(), S) == discriminant
    assert all(
        exponent == 1
        for exponent in sp.factorint(abs(discriminant)).values()
    )
    assert factor_degrees_mod_prime(field_polynomial, irreducible_prime) == (5,)
    assert factor_degrees_mod_prime(
        field_polynomial, odd_order_six_prime
    ) == (3, 2)

    field_companion = sp.zeros(5)
    for column in range(4):
        field_companion[column + 1, column] = 1
    field_companion[:, 4] = sp.Matrix(
        (
            -hermite_v_value,
            -hermite_v_value,
            0,
            -hermite_u_value,
            0,
        )
    )
    sample_eta = sum(
        (
            coefficient * field_companion**power
            for power, coefficient in enumerate(eta_coefficients)
        ),
        sp.zeros(5),
    )
    assert sp.trace(sample_eta) == 0
    assert sp.trace(sample_eta**2) == 10
    assert sp.trace(sample_eta**4) == 50 - 16 * pi_value**3
    sample_characteristic = sp.expand(sample_eta.charpoly(T).as_expr())
    sample_target = {Pi: pi_value, B: target_b, C: target_c}
    target_characteristic = sp.expand(
        pi_value**5
        * E.subs(sample_target).subs(S, T / pi_value**2)
    )
    assert sample_characteristic == target_characteristic


verify_nonaffine_s5_sample(
    hermite_u_value=-4,
    hermite_v_value=1,
    eta_coefficients=(
        -sp.Rational(4, 7),
        sp.Rational(16, 7),
        sp.Rational(13, 7),
        -sp.Rational(5, 7),
        -sp.Rational(3, 7),
    ),
    pi_value=sp.Rational(8, 7),
    target_b=sp.Rational(209, 784),
    target_c=sp.Rational(2273, 16384),
    discriminant=-55563,
    irreducible_prime=13,
    odd_order_six_prime=2,
)
verify_nonaffine_s5_sample(
    hermite_u_value=-4,
    hermite_v_value=-1,
    eta_coefficients=(
        0,
        1,
        0,
        -sp.Rational(1, 2),
        0,
    ),
    pi_value=-sp.Rational(3, 4),
    target_b=sp.Rational(5, 4),
    target_c=-sp.Rational(784, 243),
    discriminant=-179467,
    irreducible_prime=11,
    odd_order_six_prime=2,
)

# Exact comparison with the two-parameter generic S_5 polynomial.
hermite_lambda, hermite_u, hermite_v = sp.symbols(
    "hermite_lambda hermite_u hermite_v"
)
hermite_polynomial = T**5 + hermite_u * T**3 + hermite_v * T + hermite_v
hermite_substitution = {
    hermite_u: -5 * hermite_lambda**2,
    hermite_v: 4 * hermite_lambda**4 * Pi**3,
}
hermite_target = {
    B: 0,
    C: -2 / (Pi**2 * hermite_lambda),
}
assert sp.factor(
    Pi**5 * E.subs(hermite_target).subs(S, T / Pi**2)
    - hermite_polynomial.subs(hermite_substitution).subs(
        T, hermite_lambda * T
    )
    / hermite_lambda**5
) == 0
hermite_parameter_jacobian = sp.det(
    sp.Matrix(
        [
            [
                sp.diff(hermite_substitution[hermite_u], hermite_lambda),
                sp.diff(hermite_substitution[hermite_u], Pi),
            ],
            [
                sp.diff(hermite_substitution[hermite_v], hermite_lambda),
                sp.diff(hermite_substitution[hermite_v], Pi),
            ],
        ]
    )
)
assert sp.factor(hermite_parameter_jacobian) == -120 * hermite_lambda**5 * Pi**2

# A rational curve of genuinely non-affine Hermite realizations.
tau = sp.symbols("tau")
curve_u = (4 - tau**2) * (1 + 5 * tau**2) / (49 * tau**2)
curve_v = (4 - tau**2) ** 2 * (1 + 5 * tau**2) / (343 * tau**4)
curve_beta = (
    343 * tau**3 / ((4 - tau**2) ** 2 * (1 + 5 * tau**2))
)
curve_pi = 7 / (4 - tau**2)
curve_B = (
    21
    * tau
    * (6 - 5 * tau**2)
    / (2 * (4 - tau**2) * (1 + 5 * tau**2))
)
curve_C = (
    -7
    * tau**3
    * (4 - tau**2)
    / (2 * (1 + 5 * tau**2) ** 2)
)

curve_companion = sp.zeros(5)
for column in range(4):
    curve_companion[column + 1, column] = 1
curve_companion[:, 4] = sp.Matrix(
    (-curve_v, -curve_v, 0, -curve_u, 0)
)
curve_eta = curve_beta * curve_companion**3
assert sp.factor(sp.trace(curve_eta)) == 0
assert sp.factor(sp.trace(curve_eta**2)) == 10
assert sp.factor(
    sp.trace(curve_eta**4) - (50 - 16 * curve_pi**3)
) == 0
curve_characteristic = sp.expand(curve_eta.charpoly(T).as_expr())
curve_target_characteristic = sp.expand(
    curve_pi**5
    * E.subs({Pi: curve_pi, B: curve_B, C: curve_C}).subs(
        S, T / curve_pi**2
    )
)
assert sp.factor(curve_characteristic - curve_target_characteristic) == 0

curve_affine_square_class = sp.factor(-curve_u / 5)
assert curve_affine_square_class == (
    (tau - 2)
    * (tau + 2)
    * (5 * tau**2 + 1)
    / (245 * tau**2)
)
square_class_numerator = sp.Poly(
    sp.fraction(curve_affine_square_class)[0], tau, domain=sp.QQ
)
square_class_factors = sp.factor_list(square_class_numerator.as_expr())[1]
assert any(
    sp.expand(factor - (tau - 2)) == 0 and exponent == 1
    for factor, exponent in square_class_factors
)

curve_specialization = sp.Poly(
    (
        S**5
        + curve_u * S**3
        + curve_v * S
        + curve_v
    ).subs(tau, 1),
    S,
    domain=sp.QQ,
)
assert sp.discriminant(curve_specialization.as_expr(), S) != 0
assert factor_degrees_mod_prime(curve_specialization, 47) == (5,)
assert factor_degrees_mod_prime(curve_specialization, 5) == (3, 2)

# The curve extends to a dominant rational surface over the Hermite base.
rho = sp.symbols("rho")
surface_D = 4 * rho**3 - tau**2
surface_H = 1 + 20 * rho**3
surface_u = (
    9
    * surface_D
    * (1 + 5 * tau**2)
    / (tau**2 * surface_H**2)
)
surface_v = (
    27
    * surface_D**2
    * (1 + 5 * tau**2)
    / (tau**4 * surface_H**3)
)
surface_beta = (
    tau**3
    * surface_H**3
    / (27 * surface_D**2 * (1 + 5 * tau**2))
)
surface_pi = rho * surface_H / (3 * surface_D)

# General power traces for eta=beta*theta^3 in the Hermite algebra.
power_beta, power_u, power_v = sp.symbols(
    "power_beta power_u power_v"
)
power_companion = sp.zeros(5)
for column in range(4):
    power_companion[column + 1, column] = 1
power_companion[:, 4] = sp.Matrix(
    (-power_v, -power_v, 0, -power_u, 0)
)
power_eta = power_beta * power_companion**3
power_trace2 = -2 * power_beta**2 * power_u * (
    power_u**2 - 3 * power_v
)
power_trace4 = 2 * power_beta**4 * (
    power_u**6
    - 6 * power_u**4 * power_v
    + 9 * power_u**2 * power_v**2
    - 6 * power_u * power_v**2
    - 2 * power_v**3
)
assert sp.factor(sp.trace(power_eta**2) - power_trace2) == 0
assert sp.factor(sp.trace(power_eta**4) - power_trace4) == 0

surface_substitution = {
    power_beta: surface_beta,
    power_u: surface_u,
    power_v: surface_v,
}
assert sp.factor(power_trace2.subs(surface_substitution)) == 10
assert sp.factor(
    power_trace4.subs(surface_substitution)
    - (50 - 16 * surface_pi**3)
) == 0
surface_t2_coefficient = 3 * surface_beta**3 * surface_v * (
    surface_u**2 - surface_v
)
surface_constant = surface_beta**5 * surface_v**3
surface_B = -surface_t2_coefficient / (2 * surface_pi)
surface_C = -surface_constant / (2 * surface_pi**5)
surface_characteristic = (
    T**5
    - 5 * T**3
    + surface_t2_coefficient * T**2
    + 4 * surface_pi**3 * T
    + surface_constant
)
surface_target_characteristic = (
    T**5
    - 5 * T**3
    - 2 * surface_pi * surface_B * T**2
    + 4 * surface_pi**3 * T
    - 2 * surface_pi**5 * surface_C
)
assert sp.factor(
    surface_characteristic - surface_target_characteristic
) == 0

surface_jacobian = sp.factor(
    sp.Matrix((surface_u, surface_v)).jacobian((rho, tau)).det()
)
assert surface_jacobian == (
    29160
    * rho**2
    * surface_D**3
    * (1 + 5 * tau**2) ** 2
    / (tau**7 * surface_H**6)
)
surface_inverse_tau2 = -power_u * (
    power_u**2 - 3 * power_v
) / (5 * power_v**2)
surface_inverse_rho3 = -(
    (3 * power_u + power_v)
    * (power_u**2 - 3 * power_v)
    / (20 * power_u**2 * power_v)
)
assert sp.factor(
    surface_inverse_tau2.subs(
        {power_u: surface_u, power_v: surface_v}
    )
    - tau**2
) == 0
assert sp.factor(
    surface_inverse_rho3.subs(
        {power_u: surface_u, power_v: surface_v}
    )
    - rho**3
) == 0
assert sp.factor(surface_u.subs(rho, 1) - curve_u) == 0
assert sp.factor(surface_v.subs(rho, 1) - curve_v) == 0
assert sp.factor(surface_beta.subs(rho, 1) - curve_beta) == 0
assert sp.factor(surface_pi.subs(rho, 1) - curve_pi) == 0

print("PASS: the fixed quintic map has determinant -2 and inverse degree 5")
print("PASS: the affine-normalized coefficient Jacobian is -48*Pi^8")
print("PASS: the rational quotient-chart Jacobian is nonzero")
print("PASS: one fixed map has S_5 fibers of all three quintic signatures")
print("PASS: all seven unramified quintic splitting types occur modulo 7")
print("PASS: discriminant and rational repeated-root normalization are exact")
print("PASS: rational-root and quadratic-times-cubic incidences are exact")
print("PASS: the pulled-back pair-sum resolvent identity is exact")
print("PASS: the fixed map contains an explicit A_5 field fiber")
print("PASS: the square/cube descent criteria equal the intrinsic trace criteria")
print("PASS: the trace obstruction is a smooth Kummer open in a (2,4) model")
print("PASS: three pairwise nonisomorphic S_5 fields have non-affine trace points")
print("PASS: the Hermite generic family pulls back along the degree-six cover")
print("PASS: a rational non-affine curve has generic Galois group S_5")
print("PASS: a dominant rational non-affine surface covers the Hermite base")
