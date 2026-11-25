#!/usr/bin/env python3
"""Exact terminal-seven paired replay with cap 1/8 on .10 <= u <= .48."""

from fractions import Fraction as Q

from verify_multivariate_terminal_seven_paired_light import (
    RAY_KNOTS,
    check_anchor,
    coefficient_model,
    envelope,
    radial,
)
from verify_multivariate_terminal_seven_private_linear import verify_skew_convexity
from verify_multivariate_terminal_six_intermediate_cayley import interval, log, rational, scalar
from verify_multivariate_terminal_six_private_angles import require


CAP = Q(1, 8)
FIRST_RAYS = {360: 362, 380: 320, 400: 283, 420: 248, 440: 211, 460: 161, 480: 109}
RAYS = {k: (FIRST_RAYS.get(k, row[0]), *row[1:])
        for k, row in RAY_KNOTS.items() if 100 <= k <= 480}


def tight_envelope(order):
    rho, groups = envelope(CAP, order)
    if order == 6:
        return rho, groups
    x = rational(6 * CAP / (7 - CAP))
    g = tuple(log(scalar(1) + scalar(h) * x) for h in range(7))
    inner = (scalar(6) * g[3] - g[1] - g[2] - g[4] - g[5] - g[6]) / (scalar(7) * g[6])
    return rho, ((groups[0][0], inner.symmetric_error()),) + groups[1:]


def cosines(left, right, index, u):
    a, b = RAYS[left][index], RAYS[right][index]
    values = tuple(Q(a, 1000) + Q(b - a, right - left) * (v - Q(left, 1000)) for v in u)
    return interval(min(values), max(values))


def verify_ray(order, index, left, right, env):
    pending = [((Q(left, 1000), Q(right, 1000)), 0)]
    count = 0
    while pending:
        u, depth = pending.pop()
        used = radial(order, cosines(left, right, index, u), coefficient_model(order, u, env))
        if used is not None:
            count += used
            continue
        midpoint = sum(u) / 2
        require(depth < 15, ("unresolved light-cap ray", order, index, left, right, u))
        pending.extend((((u[0], midpoint), depth + 1), ((midpoint, u[1]), depth + 1)))
    return count


def main():
    verify_skew_convexity()
    check_anchor()
    require(RAYS[100] == RAY_KNOTS[100], "matching Gaussian anchor")
    for row in RAYS.values():
        require(all(1000 > a > b > -1000 for a, b in zip(row, row[1:])), "ordered rays")
    bounds = {j: tight_envelope(j) for j in (6, 7)}
    knots = tuple(RAYS)
    for left, right in zip(knots, knots[1:]):
        counts = [verify_ray(j, i, left, right, bounds[j]) for j in (6, 7) for i in range(5)]
        print(left, right, counts, flush=True)


if __name__ == "__main__":
    main()
