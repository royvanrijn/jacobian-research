#!/usr/bin/env python3
"""Exact fixed-ray separation for terminal nine at large total curvature."""

from fractions import Fraction as Q
from math import comb

from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, radial_replay
from verify_multivariate_terminal_six_private_angles import exponential_lower, rational_cap, require


CURVATURE = Q(1493, 200)
COSINES = (
    Q(23097, 25000), Q(3741, 5000), Q(4147, 10000), Q(199, 2000),
    Q(-377, 1250), Q(-317, 500), Q(-8729, 10000),
)


def coefficient_bounds(order, curvature):
    """Outward coefficient caps using two equal exponential factors."""
    interior = []
    for k in range(1, order):
        half_argument = curvature * k * (order - k) / 4
        inverse_lower = exponential_lower(half_argument)**2
        power_bound = Q(comb(order, k)**order) / inverse_lower
        interior.append(rational_cap(power_bound, order))
    return (Q(1), *interior, Q(1))


def anchor_value(order, x):
    """Sign polynomial of the binomial anchor on an upper ray."""
    degree = order // 2
    if degree == 0:
        return Q(1)
    previous, current = Q(1), x if order % 2 == 0 else 2 * x - 1
    for _ in range(1, degree):
        previous, current = current, 2 * x * current - previous
    return current


def verify_fixed_threshold(terminal, curvature, cosines):
    """Check alternating upper root sectors for consecutive orders."""
    lower = terminal - 1
    require(len(cosines) == terminal - 2, "number of rays")
    knots = (Q(1), *cosines, Q(-1))
    require(all(left > right for left, right in zip(knots, knots[1:])), "ray order")

    total = 0
    for order in (lower, terminal):
        for index, x in enumerate(knots):
            previous_roots = (index + (order == terminal)) // 2
            require((-1) ** previous_roots * anchor_value(order, x) > 0,
                    ("anchor signs", order, index))

        caps = coefficient_bounds(order, curvature)
        require(sum(cap * INNER**k for k, cap in enumerate(caps) if k) < 1,
                ("inner radius", order))
        require(sum(cap * OUTER**k for k, cap in enumerate(caps[:-1])) < OUTER**order,
                ("outer radius", order))
        boxes = [(Q(1), Q(1)) if k in (0, order) else (Q(0), cap)
                 for k, cap in enumerate(caps)]
        for cosine in cosines:
            cells = radial_replay(order, boxes, (cosine, cosine))
            require(cells is not None, ("unresolved ray", order, cosine))
            total += cells
    return total


def main():
    count = verify_fixed_threshold(9, CURVATURE, COSINES)
    print(f"Verified {count} terminal-nine radial cells.")


if __name__ == "__main__":
    main()
