#!/usr/bin/env python3
"""Exact certificates for an arithmetic zoo in one quintic Keller map."""

from __future__ import annotations

import os
import warnings
from math import isqrt

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


x, y, z, S, T, X = sp.symbols("x y z S T X")
Pi, B, C = sp.symbols("Pi B C")


def factor_degrees_mod_prime(poly: sp.Poly, prime: int) -> tuple[int, ...]:
    """Return the squarefree factor-degree pattern modulo ``prime``."""

    discriminant = sp.discriminant(poly.as_expr(), S)
    assert int(poly.LC()) % prime != 0
    assert int(discriminant) % prime != 0
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


def inverse_polynomial(pi: sp.Rational, b: sp.Rational, c: sp.Rational) -> sp.Poly:
    return sp.Poly(
        pi**5 * S**5 - 5 * pi * S**3 - 2 * b * S**2 + 4 * S - 2 * c,
        S,
        domain=sp.QQ,
    )


def normalized_polynomial(
    pi: sp.Rational, b: sp.Rational, c: sp.Rational
) -> sp.Poly:
    """Return pi^5 E(pi^-2 T), the centered monic presentation."""

    return sp.Poly(
        T**5
        - 5 * T**3
        - 2 * pi * b * T**2
        + 4 * pi**3 * T
        - 2 * pi**5 * c,
        T,
        domain=sp.QQ,
    )


def pair_sum_resolvent(a: sp.Rational, b: sp.Rational, c: sp.Rational) -> sp.Poly:
    """Pair-sum resolvent of T^5-5T^3+aT^2+bT+c."""

    return sp.Poly(
        X**10
        - 15 * X**8
        + a * X**7
        + (75 - 3 * b) * X**6
        + (-10 * a - 11 * c) * X**5
        + (-a**2 + 10 * b - 125) * X**4
        + (-4 * a * b + 25 * a + 20 * c) * X**3
        + (5 * a**2 + 7 * a * c - 4 * b**2 + 25 * b) * X**2
        + (-a**3 + 4 * b * c - 25 * c) * X
        - a**2 * b
        - 5 * a * c
        - c**2,
        X,
        domain=sp.QQ,
    )


def companion_matrix(poly: sp.Poly) -> sp.Matrix:
    """Multiplication by the residue class of the polynomial variable."""

    assert poly.is_monic
    degree = poly.degree()
    matrix = sp.zeros(degree)
    for column in range(degree - 1):
        matrix[column + 1, column] = 1
    matrix[:, degree - 1] = sp.Matrix(
        [-coefficient for coefficient in reversed(poly.all_coeffs()[1:])]
    )
    return matrix


def valuation(value: int, prime: int) -> int:
    """Return the exact ``prime``-adic valuation of a nonzero integer."""

    value = abs(int(value))
    assert value != 0
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


# The one fixed split-seed quadratic-gauge Keller map.
t = 1 + x * y
q_source = t**2 * z - sp.Rational(4, 5) * y**2 * (1 + 3 * t)
mapping = (
    t * q_source,
    y
    - sp.Rational(15, 4) * x * q_source
    + sp.Rational(5, 4) * t**2 * x**3 * q_source**5,
    x * (5 - 3 * t)
    + sp.Rational(5, 4) * x**3 * z
    - sp.Rational(3, 4) * (x * q_source) ** 5,
)
assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -2

E = Pi**5 * S**5 - 5 * Pi * S**3 - 2 * B * S**2 + 4 * S - 2 * C


# Row 1: the completely split fiber.
split_target = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
split_poly = inverse_polynomial(*split_target)
split_factorization = S * (S - 1) * (S + 1) * (S - 2) * (S + 2)
assert sp.expand(split_poly.as_expr() - split_factorization) == 0
assert sp.discriminant(split_poly.as_expr(), S) == 82944


# Row 2: an irreducible S_5 field.  Modulo 2 it is irreducible, and
# modulo 79 it has the cycle type of a transposition.
s5_target = (sp.Integer(1), sp.Integer(0), -sp.Rational(1, 2))
s5_poly = inverse_polynomial(*s5_target)
s5_discriminant = int(sp.discriminant(s5_poly.as_expr(), S))
assert s5_poly.as_expr() == S**5 - 5 * S**3 + 4 * S + 1
assert s5_discriminant == 38569
assert isqrt(s5_discriminant) ** 2 != s5_discriminant
assert factor_degrees_mod_prime(s5_poly, 2) == (5,)
assert factor_degrees_mod_prime(s5_poly, 79) == (2, 1, 1, 1)


# Row 3: an irreducible A_5 field.  Clearing the denominator does not
# change the roots or the modular factor degrees at 5 and 7.
a5_target = (sp.Integer(1), -sp.Rational(4, 3), sp.Integer(6))
a5_poly = inverse_polynomial(*a5_target)
a5_integral = sp.Poly(3 * a5_poly.as_expr(), S, domain=sp.QQ)
assert sp.discriminant(a5_poly.as_expr(), S) == 984**2
assert factor_degrees_mod_prime(a5_integral, 5) == (5,)
assert factor_degrees_mod_prime(a5_integral, 7) == (3, 1, 1)


# Row 4: an explicit product K_2 x K_3.
product_target = (
    sp.Integer(1),
    -sp.Rational(3, 2),
    -sp.Rational(9, 2),
)
product_poly = inverse_polynomial(*product_target)
quadratic = sp.Poly(S**2 + S + 1, S, domain=sp.QQ)
cubic = sp.Poly(S**3 - S**2 - 5 * S + 9, S, domain=sp.QQ)
assert sp.expand(product_poly.as_expr() - quadratic.as_expr() * cubic.as_expr()) == 0
assert sp.discriminant(quadratic.as_expr(), S) == -3
assert factor_degrees_mod_prime(cubic, 5) == (3,)
assert sp.resultant(quadratic.as_expr(), cubic.as_expr(), S) == 181
assert sp.discriminant(product_poly.as_expr(), S) == 80198928


# Row 5: a Hasse failure.  The irreducible quadratic and cubic have the
# same quadratic resolvent field Q(sqrt(-3)).
hasse_target = (
    sp.Integer(4),
    -sp.Rational(335, 27),
    sp.Rational(4807, 20736),
)
hasse_poly = inverse_polynomial(*hasse_target)
hasse_quadratic = sp.Poly(192 * S**2 - 72 * S + 19, S, domain=sp.QQ)
hasse_cubic = sp.Poly(
    55296 * S**3 + 20736 * S**2 + 1224 * S - 253,
    S,
    domain=sp.QQ,
)
assert sp.expand(
    hasse_poly.as_expr()
    - hasse_quadratic.as_expr() * hasse_cubic.as_expr() / 10368
) == 0
assert sp.discriminant(hasse_quadratic.as_expr(), S) == -3 * 56**2
assert factor_degrees_mod_prime(hasse_cubic, 7) == (3,)
assert sp.discriminant(hasse_cubic.as_expr(), S) == -3 * 28366848**2
assert sp.resultant(
    hasse_quadratic.as_expr(),
    hasse_cubic.as_expr(),
    S,
) == 93138374098944

# The displayed integral models can be bad only at 2, 3, 7, and 19.
possible_ramified_primes = set(
    sp.factorint(
        abs(
            int(hasse_quadratic.LC())
            * int(hasse_cubic.LC())
            * int(sp.discriminant(hasse_quadratic.as_expr(), S))
            * int(sp.discriminant(hasse_cubic.as_expr(), S))
        )
    )
)
assert possible_ramified_primes == {2, 3, 7, 19}

# At 2, put U=8S.  The transformed cubic has the simple root U=1 mod 2.
U = sp.symbols("U")
cubic_at_2 = sp.Poly(
    sp.expand(hasse_cubic.as_expr().subs(S, U / 8)),
    U,
    domain=sp.QQ,
)
assert cubic_at_2.as_expr() == 108 * U**3 + 324 * U**2 + 153 * U - 253
assert int(cubic_at_2.eval(1)) % 2 == 0
assert int(cubic_at_2.diff().eval(1)) % 2 != 0

# At 3, put U=3S and apply strong Hensel at U=-7.
cubic_at_3 = sp.Poly(
    sp.expand(hasse_cubic.as_expr().subs(S, U / 3)),
    U,
    domain=sp.QQ,
)
assert cubic_at_3.as_expr() == 2048 * U**3 + 2304 * U**2 + 408 * U - 253
assert valuation(cubic_at_3.eval(-7), 3) == 7
assert valuation(cubic_at_3.diff().eval(-7), 3) == 2

# Strong Hensel at S=5 handles 7, and ordinary Hensel at S=0 handles 19.
assert valuation(hasse_quadratic.eval(5), 7) == 3
assert valuation(hasse_quadratic.diff().eval(5), 7) == 1
assert int(hasse_quadratic.eval(0)) % 19 == 0
assert int(hasse_quadratic.diff().eval(0)) % 19 != 0

hasse_discriminant = sp.factor(sp.discriminant(hasse_poly.as_expr(), S))
assert hasse_discriminant == sp.Rational(28002486064, 729) ** 2

# The improved Hasse target has primitive projective coordinates
# [20736:82944:-257280:4807], of height 257280.
assert max(abs(value) for value in (20736, 82944, -257280, 4807)) == 257280


# The remaining transitive quintic groups occur in the same inverse pencil.
#
# C_5.  Let theta be a root of the real 11th-cyclotomic subfield polynomial.
# The substitution theta -> theta^2-2 has exact order five.  The displayed
# eta is primitive and has the normalized target polynomial as its
# characteristic polynomial, so its degree-five field has five automorphisms.
c5_target = (sp.Integer(1), -sp.Rational(15, 11), sp.Rational(331, 242))
c5_poly = inverse_polynomial(*c5_target)
c5_integral = sp.Poly(121 * c5_poly.as_expr(), S, domain=sp.QQ)
assert factor_degrees_mod_prime(c5_integral, 2) == (5,)
assert sp.discriminant(c5_poly.as_expr(), S) == sp.Rational(
    109 * 2663, 11**4
) ** 2

theta = sp.symbols("theta")
cyclotomic_subfield = sp.Poly(
    theta**5 + theta**4 - 4 * theta**3 - 3 * theta**2 + 3 * theta + 1,
    theta,
    domain=sp.QQ,
)
theta_image = theta**2 - 2
assert sp.rem(
    cyclotomic_subfield.as_expr().subs(theta, theta_image),
    cyclotomic_subfield.as_expr(),
    theta,
) == 0
iterate = theta
for exponent in range(1, 6):
    iterate = sp.rem(
        iterate.subs(theta, theta_image),
        cyclotomic_subfield.as_expr(),
        theta,
    )
    if exponent < 5:
        assert sp.expand(iterate - theta) != 0
assert sp.expand(iterate - theta) == 0

theta_operator = companion_matrix(cyclotomic_subfield)
c5_eta = (
    -3 * theta_operator**4
    + theta_operator**3
    - 4 * theta_operator
    + 15 * sp.eye(5)
) / 11
assert sp.Poly(c5_eta.charpoly(T).as_expr(), T, domain=sp.QQ) == normalized_polynomial(
    *c5_target
)
c5_resolvent = pair_sum_resolvent(
    sp.Rational(30, 11), sp.Integer(4), -sp.Rational(331, 121)
)
c5_resolvent_factors = (
    sp.Poly(
        121 * X**5 - 1210 * X**3 + 715 * X**2 + 1364 * X + 23,
        X,
        domain=sp.QQ,
    ),
    sp.Poly(
        121 * X**5 - 605 * X**3 - 385 * X**2 + 209 * X + 43,
        X,
        domain=sp.QQ,
    ),
)
assert all(factor.is_irreducible for factor in c5_resolvent_factors)
assert sp.expand(
    c5_resolvent.as_expr()
    - c5_resolvent_factors[0].as_expr()
    * c5_resolvent_factors[1].as_expr()
    / 11**4
) == 0

# D_5.  Irreducibility gives transitivity.  The pair-sum resolvent has two
# irreducible quintic factors, leaving C_5 or D_5; the (2,2,1) Frobenius
# pattern modulo 7 excludes C_5.
d5_target = (
    sp.Rational(5, 2),
    -sp.Rational(27, 8),
    -sp.Rational(738, 3125),
)
d5_poly = inverse_polynomial(*d5_target)
d5_normalized = normalized_polynomial(*d5_target)
d5_integral = sp.Poly(8 * d5_normalized.as_expr(), T, domain=sp.QQ)
assert factor_degrees_mod_prime(
    sp.Poly(d5_integral.as_expr().subs(T, S), S, domain=sp.QQ), 11
) == (5,)
assert factor_degrees_mod_prime(
    sp.Poly(d5_integral.as_expr().subs(T, S), S, domain=sp.QQ), 7
) == (2, 2, 1)
assert sp.discriminant(d5_normalized.as_expr(), T) == sp.Rational(
    2048625, 256
) ** 2
d5_resolvent = pair_sum_resolvent(
    sp.Rational(135, 8), sp.Rational(125, 2), sp.Rational(369, 8)
)
d5_resolvent_factors = (
    sp.Poly(
        16 * X**5 + 60 * X**3 + 500 * X**2 - 585 * X + 2196,
        X,
        domain=sp.QQ,
    ),
    sp.Poly(
        32 * X**5 - 600 * X**3 - 460 * X**2 - 180 * X - 5553,
        X,
        domain=sp.QQ,
    ),
)
assert all(factor.is_irreducible for factor in d5_resolvent_factors)
assert sp.expand(
    d5_resolvent.as_expr()
    - d5_resolvent_factors[0].as_expr()
    * d5_resolvent_factors[1].as_expr()
    / 512
) == 0

# F_20.  Work first in the sparse De Moivre presentation
# H=T^5-10T^3+20T+20.  Cayley's solvability sextic
# P(Z)^2-2^10 Disc(H) Z has a rational root, so the transitive group lies
# in F_20.  Its nonsquare discriminant excludes C_5 and D_5.  The short
# element eta below transports this field to the fixed inverse pencil.
f20_source = sp.Poly(T**5 - 10 * T**3 + 20 * T + 20, T, domain=sp.QQ)
assert factor_degrees_mod_prime(
    sp.Poly(f20_source.as_expr().subs(T, S), S, domain=sp.QQ), 29
) == (5,)
assert sp.discriminant(f20_source.as_expr(), T) == 231200000
cayley_P = X**3 - 700 * X**2 + 110000 * X - 15880000
cayley_sextic = sp.Poly(
    cayley_P**2 - 2**10 * 231200000 * X,
    X,
    domain=sp.QQ,
)
assert sp.factor(cayley_sextic.as_expr()) == (X - 500) * (
    X**5
    - 900 * X**4
    + 260000 * X**3
    - 55760000 * X**2
    + 6452000000 * X
    - 504348800000
)

f20_target = (
    sp.Rational(31, 5),
    sp.Rational(5229, 310),
    sp.Rational(9618099, 114516604),
)
f20_operator = companion_matrix(f20_source)
f20_eta = (
    sp.Rational(22, 5) * sp.eye(5)
    + sp.Rational(9, 5) * f20_operator
    - sp.Rational(11, 10) * f20_operator**2
    - sp.Rational(3, 5) * f20_operator**3
)
assert sp.trace(f20_eta) == 0
assert sp.trace(f20_eta**2) == 10
assert sp.trace(f20_eta**4) == 50 - 16 * sp.Rational(31, 5) ** 3
assert sp.Poly(
    f20_eta.charpoly(T).as_expr(), T, domain=sp.QQ
) == normalized_polynomial(*f20_target)
f20_integral = sp.Poly(
    6250 * normalized_polynomial(*f20_target).as_expr(), T, domain=sp.QQ
)
assert factor_degrees_mod_prime(
    sp.Poly(f20_integral.as_expr().subs(T, S), S, domain=sp.QQ), 29
) == (5,)
f20_discriminant = sp.factor(
    sp.discriminant(normalized_polynomial(*f20_target).as_expr(), T)
)
assert f20_discriminant == sp.Rational(
    17 * 6665679406891, 4 * 5**8
) ** 2 * 5

# A tempting route to infinitely many Hasse failures is the classical
# pure-cubic family Q(sqrt(-3)) x Q(cuberoot(m)).  It cannot meet this
# normalized trace chart even after independent affine changes on the two
# factors.  If eta_2=u*sqrt(-3)+v and eta_3=w*cuberoot(m)+s, trace zero
# gives s=-2v/3 and trace(eta^2)=10 would force
#     5 v^2 - 9 u^2 = 15.
# The projective conic 5V^2-9U^2=15W^2 has no Q_5-point: modulo 5 first
# 5|U, and after division by 5 one gets V^2=3W^2 mod 5.
u_trace, v_trace = sp.symbols("u_trace v_trace")
s_trace = -sp.Rational(2, 3) * v_trace
pure_trace2 = 2 * (-3 * u_trace**2 + v_trace**2) + 3 * s_trace**2
assert sp.expand(
    pure_trace2
    - 10
    - sp.Rational(2, 3) * (5 * v_trace**2 - 9 * u_trace**2 - 15)
) == 0


# Every displayed row is on the complete reconstruction chart and is
# squarefree of degree five.
targets_and_polynomials = (
    (split_target, split_poly),
    (s5_target, s5_poly),
    (a5_target, a5_poly),
    (product_target, product_poly),
    (hasse_target, hasse_poly),
    (c5_target, c5_poly),
    (d5_target, d5_poly),
    (f20_target, inverse_polynomial(*f20_target)),
)
for target, polynomial in targets_and_polynomials:
    assert target[0] != 0
    assert polynomial.degree() == 5
    assert sp.discriminant(polynomial.as_expr(), S) != 0
    assert sp.expand(E.subs(dict(zip((Pi, B, C), target))) - polynomial.as_expr()) == 0


print("PASS: one fixed degree-five map has determinant -2")
print("PASS: its split, S_5, A_5, and K_2 x K_3 fibers are exact")
print("PASS: the same map realizes C_5, D_5, F_20, A_5, and S_5")
print("PASS: cyclic automorphism, pair-resolvent, and Cayley certificates are exact")
print("PASS: the Hasse row has common quadratic resolvent Q(sqrt(-3))")
print("PASS: exact local witnesses cover 2, 3, 7, and 19")
print("PASS: the improved Hasse target has projective height 257280")
print("PASS: the pure-cubic infinitude ansatz has a Q_5 trace obstruction")
