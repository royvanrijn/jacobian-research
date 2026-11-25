#!/usr/bin/env python3
"""Prove linked terminal-seven separation at fully shared curvature at least 17/4."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod

import sympy as sym

from verify_multivariate_terminal_seven_high_curvature import (
    CURVATURE as EXISTING_CURVATURE,
)
from verify_multivariate_terminal_six_joint_factor_links import check_nonnegative_on_interval, x
from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER,
    OUTER,
    coefficient_box,
    log_bounds,
    radial_replay,
    rounded_interval,
)
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds,
    exponential_lower,
    multiply_interval,
    require,
)


CURVATURE = Q(17, 4)
FLOOR_CAP = Q(1, 6)
RAYS = {
    0: (88552, 57226, 11736, -37095, -75204),
    24: (86639, 55636, 5417, -40814, -78573),
    48: (84332, 46762, 592, -42237, -80429),
    72: (81763, 44625, -3102, -42593, -81914),
    96: (78847, 42437, -6619, -43789, -82954),
    120: (75429, 39627, -10504, -44058, -84344),
    144: (72249, 36208, -13400, -46006, -84665),
    168: (69296, 33944, -13411, -47066, -84571),
    171: (69118, 33440, -13498, -47348, -84527),
}


def possible_rows():
    """Rows that can raise the order-six cap above its floor."""
    for gap in range(1, int(EXISTING_CURVATURE * (1 + 1 / (6 * FLOOR_CAP))) + 1):
        for lower in range(1, int(gap / (6 * FLOOR_CAP) - Q(1, 6)) + 1):
            curvature = Q(gap**2, lower + gap)
            x6 = Q(gap, 6 * lower + 1)
            if x6 <= FLOOR_CAP or curvature >= EXISTING_CURVATURE:
                continue
            u_squared = prod(Q(6 * lower + h, 6 * lower + 5 * gap + h)
                             for h in range(1, gap + 1))
            other_loss = exponential_lower(Q(5, 6) * max(CURVATURE - curvature, 0))
            yield x6, Q(gap, 7 * lower + 1), u_squared / other_loss


def reciprocal_caps(rows, outer_lower):
    lower = FLOOR_CAP
    upper = 6 * FLOOR_CAP / (7 - FLOOR_CAP)
    for x6, x7, u_squared_cap in rows:
        if outer_lower**2 <= u_squared_cap:
            lower = max(lower, x6)
            upper = max(upper, x7)
    return lower, upper


def skew_weights(order, k):
    return tuple((Q(k <= m < order-k) - Q(order - 2*k, order)) / 2
                 for m in range(order))


def verify_skew_convexity():
    for order in (6, 7):
        for k in range(1, (order + 1) // 2):
            first = sum(w * m / (1 + m*x) for m, w in enumerate(skew_weights(order, k)))
            derivative_scale = 1 + (order - 1)*x
            second = derivative_scale * sym.diff(derivative_scale * first, x)
            check_nonnegative_on_interval(second)


def weighted_log_bounds(weights, cap):
    intervals = [multiply_interval((weight, weight), log_bounds(1 + m*cap))
                 for m, weight in enumerate(weights) if weight and m]
    if not intervals:
        return Q(0), Q(0)
    return tuple(map(sum, zip(*intervals, strict=True)))


@lru_cache(maxsize=None)
def exponent_bounds(cap6, cap7):
    denominator_lower, denominator_upper = log_bounds(1 + 5*cap6)
    rho_lower = log_bounds(1 + Q(36, 7)*cap6)[0] / denominator_upper
    rho_upper = 1 + log_bounds(1 / (1-cap6/7))[1] / denominator_lower
    require(rho_lower <= Q(36, 35) <= rho_upper, "Gaussian ratio")
    result = {}
    for order, cap, ratio in ((6, cap6, (1, 1)), (7, cap7, (rho_lower, rho_upper))):
        local_outer = log_bounds(1 + (order-1)*cap)[1] / 2
        skew_outer = log_bounds(1 + (order-1)*cap)[0] / 2
        bounds = []
        for k in range(1, order):
            distance = min(k, order-k)
            mean_lower = sum(log_bounds((1+(order-1-m)*cap)/(1+m*cap))[0]
                             for m in range(distance)) / (2 * local_outer)
            skew_upper = weighted_log_bounds(skew_weights(order, distance), cap)[1] / skew_outer
            mean_upper = Q(k * (order-k), order-1)
            bounds.append(rounded_interval(ratio[0] * (mean_lower-skew_upper),
                                           ratio[1] * (mean_upper+skew_upper)))
            require(0 < bounds[-1][0] <= mean_upper * (1 if order == 6 else Q(36, 35))
                    <= bounds[-1][1], "Gaussian anchor envelope")
        result[order] = tuple(bounds)
    return result


def boxes(u, bounds, total_caps):
    sextic = [(lo, min(hi, cap)) for (lo, hi), cap in
              zip(coefficient_box(6, u, bounds[6]), total_caps[6], strict=True)]
    septic = [(lo, min(hi, cap)) for (lo, hi), cap in
              zip(coefficient_box(7, u, bounds[7]), total_caps[7], strict=True)]
    previous = [hi / comb(6, k) for k, (_, hi) in enumerate(sextic)]
    for k in range(1, 7):
        lo, hi = septic[k]
        left = previous[k] / exponential_lower(CURVATURE * k*k * Q(1, 91))
        right = previous[k-1] / exponential_lower(CURVATURE * (7-k)**2 * Q(1, 91))
        septic[k] = lo, min(hi, comb(7, k) * min(left, right))
    return {6: sextic, 7: septic}


def ray_at(index, left, right, u):
    first, last = Q(RAYS[left][index], 100000), Q(RAYS[right][index], 100000)
    return first + (last-first) * (1000*u-left) / (right-left)


def verify_knot(order, index, left, right, rows, total_caps):
    stack = [((Q(left, 1000), Q(right, 1000)), 0)]
    count = 0
    while stack:
        u, depth = stack.pop()
        caps = reciprocal_caps(rows, u[0])
        coefficients = boxes(u, exponent_bounds(*caps), total_caps)
        if any(lo > hi for array in coefficients.values() for lo, hi in array):
            continue
        cosines = tuple(sorted(ray_at(index, left, right, end) for end in u))
        cells = radial_replay(order, coefficients[order], cosines)
        if cells is not None:
            count += cells
            continue
        require(depth < 16, ("unresolved ray", order, index, u, caps))
        midpoint = sum(u) / 2
        stack.extend((((u[0], midpoint), depth+1), ((midpoint, u[1]), depth+1)))
    return count


def verify_anchor():
    for values in RAYS.values():
        cosines = (100000, *values, -100000)
        require(all(lo > hi for lo, hi in zip(cosines, cosines[1:])), "ray order")
    knots = (Q(1), *(Q(v, 100000) for v in RAYS[0]), Q(-1))
    for values, signs in (
        ((8*c**3-6*c for c in knots), (1, 1, -1, -1, 1, 1, -1)),
        ((8*c**3-4*c**2-4*c+1 for c in knots), (1, -1, -1, 1, 1, -1, -1)),
    ):
        require(all(value*sign > 0 for value, sign in zip(values, signs, strict=True)),
                "binomial anchor signs")


def verify():
    verify_skew_convexity()
    verify_anchor()
    knots = sorted(RAYS)
    require(Q(1, 84) > Q(1, 91), "Gaussian endpoint attenuation")
    require(exponential_lower(Q(5, 12) * CURVATURE) * Q(knots[-1], 1000) > 1,
            "outer mean cutoff")
    rows = tuple(possible_rows())
    require(all(x6 <= 2 and x7 <= 2 for x6, x7, _ in rows), "convexity range")
    total_caps = {j: coefficient_bounds(j, (), CURVATURE) for j in (6, 7)}
    for j in (6, 7):
        caps = total_caps[j]
        require(sum(a * INNER**k for k, a in enumerate(caps) if k) < 1,
                ("small circle", j))
        require(sum(a * OUTER**k for k, a in enumerate(caps[:-1])) < OUTER**j,
                ("large circle", j))
    total = 0
    for left, right in zip(knots, knots[1:]):
        cells = sum(verify_knot(j, i, left, right, rows, total_caps)
                    for j in (6, 7) for i in range(5))
        total += cells
        print(f"[{left}/1000,{right}/1000]: {cells} radial cells.", flush=True)
    print(f"Verified {total} radial cells.", flush=True)


if __name__ == "__main__":
    verify()
