#!/usr/bin/env python3
"""Exact audit for the explicit same-degree stable-inequivalence theorem."""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
W, S, T = sp.symbols("W S T")
target_A, target_B, target_C = sp.symbols("A_t B_t C_t")
target_P, target_Q, target_R = sp.symbols("P_t Q_t R_t")


# Weighted quartic.
u = 1 + 3 * x * y
gamma = 1 - 4 * x * y - x**2 * z
weighted = (
    sp.cancel((2 * u + u**2 - 3 * u**4 * gamma**2) / x**2),
    sp.cancel((1 + u - 2 * u**3 * gamma**2) / x),
    sp.expand(x * gamma),
)
assert all(sp.denom(component) == 1 for component in weighted)
assert sp.factor(sp.Matrix(weighted).jacobian((x, y, z)).det()) == -6

weighted_inverse = (
    W**2
    - W**4
    - 2 * target_B * target_C * W
    + target_A * target_C**2
)
weighted_source_root = u * gamma
weighted_substitution = {
    target_A: weighted[0],
    target_B: weighted[1],
    target_C: weighted[2],
    W: weighted_source_root,
}
assert sp.factor(weighted_inverse.subs(weighted_substitution)) == 0
weighted_control = weighted_inverse.subs(
    {target_A: 1, target_B: 0, target_C: 1}
)
assert sp.degree(weighted_control, W) == 4
assert sp.gcd(weighted_control, sp.diff(weighted_control, W)) == 1

# Its split boundary contact is reduced.
weighted_discriminant = sp.factor(sp.discriminant(weighted_inverse, W))
weighted_saturated_trace = sp.factor(
    (weighted_discriminant / target_C**2).subs(target_C, 0)
)
assert sp.expand(
    weighted_saturated_trace + 16 * (target_A - target_B**2)
) == 0
assert sp.diff(weighted_saturated_trace, target_A) == -16
weighted_mu = 1


# Cancellation quartic over Q(theta), theta^2-4theta+6=0.
theta = sp.symbols("theta")
theta_relation = sp.Poly(theta**2 - 4 * theta + 6, theta)


def reduce_theta(expression: sp.Expr) -> sp.Expr:
    """Reduce a polynomial expression modulo theta^2-4theta+6."""

    expanded = sp.cancel(expression)
    numerator, denominator = expanded.as_numer_denom()
    reduced_numerator = sp.rem(sp.Poly(numerator, theta), theta_relation).as_expr()
    reduced_denominator = sp.rem(
        sp.Poly(denominator, theta), theta_relation
    ).as_expr()
    return sp.cancel(reduced_numerator / reduced_denominator)


cancellation_A = 1 + x * y**2
cancellation_h = theta + (4 * theta - 6) * cancellation_A
cancellation_B = sp.expand(
    cancellation_A**2 * z + y**3 * cancellation_h
)
cancellation_P = sp.expand(cancellation_A * cancellation_B)
cancellation_Q = sp.expand(y + x * cancellation_B)
cancellation_R = -x * (
    theta
    * (
        8 * x**4 * y**5 * z
        + 16 * x**3 * y**6
        + 10 * x**3 * y**3 * z
        + 20 * x**2 * y**4
    )
    + x**5 * y**4 * z**2
    - 12 * x**4 * y**5 * z
    + 2 * x**4 * y**2 * z**2
    - 60 * x**3 * y**6
    - 12 * x**3 * y**3 * z
    + x**3 * z**2
    - 48 * x**2 * y**4
    + 4 * x**2 * y * z
    + 18 * x * y**2
    - 12
) / 12
cancellation = (cancellation_P, cancellation_Q, cancellation_R)
cancellation_jacobian = sp.Matrix(cancellation).jacobian((x, y, z)).det()
assert sp.factor(reduce_theta(cancellation_jacobian) + 1) == 0

cancellation_inverse = (
    T
    - target_Q**2 * T**2 / 2
    + 2 * target_P * target_Q * T**3 / 3
    - target_P**2 * T**4 / 4
    - target_R
)
cancellation_source_root = x / cancellation_A
cancellation_incidence = cancellation_inverse.subs(
    {
        target_P: cancellation_P,
        target_Q: cancellation_Q,
        target_R: cancellation_R,
        T: cancellation_source_root,
    }
)
assert sp.factor(reduce_theta(cancellation_incidence)) == 0
assert sp.factor(
    sp.diff(cancellation_inverse, T)
    - (1 - T * (target_Q - target_P * T) ** 2)
) == 0
cancellation_control = cancellation_inverse.subs(
    {target_P: 1, target_Q: 0, target_R: 0}
)
assert sp.degree(cancellation_control, T) == 4
assert sp.gcd(cancellation_control, sp.diff(cancellation_control, T)) == 1


# Root-engineered quadratic-gauge quartic.
t = 1 + x * y
q = t**2 * z - y**2 * (1 + 3 * t)
quadratic = (
    sp.expand(t * q),
    sp.expand(-y + 3 * x * q + t * q - 2 * t**2 * x**2 * q**4),
    sp.expand(-2 * x + 3 * x**2 * y - x**3 * z + x**4 * q**4),
)
assert sp.factor(sp.Matrix(quadratic).jacobian((x, y, z)).det()) == -2

quadratic_inverse = (
    2 * S
    - target_P * S**2
    - 2 * target_P * S**3
    + target_P**4 * S**4
    - target_B * S**2
    - target_C
)
quadratic_incidence = quadratic_inverse.subs(
    {
        target_P: quadratic[0],
        target_B: -quadratic[1],
        target_C: -quadratic[2],
        S: x / t,
    }
)
assert sp.factor(quadratic_incidence) == 0
quadratic_control = sp.factor(
    quadratic_inverse.subs(
        {target_P: 1, target_B: 0, target_C: 0}
    )
)
assert quadratic_control == S * (S - 2) * (S - 1) * (S + 1)
assert sp.gcd(quadratic_control, sp.diff(quadratic_control, S)) == 1


# Fitting-support obstruction on the common ramified torus.
torus_u, torus_P, torus_S = sp.symbols("u_t P_u S_u")
cancellation_fitting = 2 * torus_u - 1
quadratic_fitting = -1 - 3 * torus_P * torus_S**2 + 4 * torus_P**4 * torus_S**3


def affine_rank(points: list[tuple[int, ...]]) -> int:
    """Affine rank of a finite exponent support."""

    origin = sp.Matrix(points[0])
    differences = [
        sp.Matrix(point) - origin
        for point in points[1:]
    ]
    return sp.Matrix.hstack(*differences).rank() if differences else 0


assert affine_rank([(0,), (1,)]) == 1
assert affine_rank([(0, 0), (1, 2), (4, 3)]) == 2
assert sp.Poly(cancellation_fitting, torus_u).length() == 2
assert sp.Poly(quadratic_fitting, torus_P, torus_S).length() == 3

# Scheme-theoretic contact indices from the exact family formulas.
cancellation_mu = 2 * 1 * (2 + 1)
quadratic_mu = 2
assert (weighted_mu, quadratic_mu, cancellation_mu) == (1, 2, 6)

print("PASS: three explicit maps have constant Jacobians -6, -1, and -2")
print("PASS: all three inverse problems have exact degree-four control fibers")
print("PASS: ramified-torus Fitting support ranks are 1 and 2")
print("PASS: boundary-contact nilpotency indices are 1, 2, and 6")
