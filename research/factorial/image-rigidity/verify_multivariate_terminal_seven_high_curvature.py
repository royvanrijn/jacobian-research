#!/usr/bin/env python3
"""Certify terminal-seven angular separation for curvature at least 73/10."""

from fractions import Fraction as Q
from math import comb, factorial

from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, radial_replay
from verify_multivariate_terminal_six_private_angles import require


CURVATURE = Q(73, 10)
COSINES = (Q(2167, 2500), Q(2969, 5000), Q(14, 125), Q(-3303, 10000), Q(-379, 500))


def exponential_cap(x):
    return 1 / sum(x**k / factorial(k) for k in range(33))


def verify_anchor():
    knots = (Q(1), *COSINES, Q(-1))
    require(all(left > right for left, right in zip(knots, knots[1:])), "ray order")
    sextic = (8*x**3 - 6*x for x in knots)
    septic = (8*x**3 - 4*x**2 - 4*x + 1 for x in knots)
    require(all(value * sign > 0 for value, sign in
                zip(sextic, (1, 1, -1, -1, 1, 1, -1), strict=True)), "sextic anchor signs")
    require(all(value * sign > 0 for value, sign in
                zip(septic, (1, -1, -1, 1, 1, -1, -1), strict=True)), "septic anchor signs")


def verify():
    verify_anchor()
    total = 0
    for order in (6, 7):
        coefficients = [(Q(1), Q(1))]
        coefficients.extend(
            (Q(0), comb(order, k) * exponential_cap(CURVATURE * k * (order-k) / (2*order)))
            for k in range(1, order)
        )
        coefficients.append((Q(1), Q(1)))
        require(sum(hi * INNER**k for k, (_, hi) in enumerate(coefficients) if k) < 1,
                ("inner radius", order))
        require(sum(hi * OUTER**k for k, (_, hi) in enumerate(coefficients[:-1])) < OUTER**order,
                ("outer radius", order))
        for cosine in COSINES:
            count = radial_replay(order, coefficients, (cosine, cosine))
            require(count is not None, ("unresolved ray", order, cosine))
            total += count
    print(f"Certified both terminal-seven arrays on {total} radial cells.")


if __name__ == "__main__":
    verify()
