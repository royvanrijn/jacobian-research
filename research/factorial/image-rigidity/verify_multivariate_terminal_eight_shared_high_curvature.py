#!/usr/bin/env python3
"""Prove linked terminal-eight separation at fully shared curvature at least 19/4."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod

import sympy as sym

from verify_multivariate_terminal_eight_high_curvature import (
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


CURVATURE = Q(19, 4)
FLOOR_CAP = Q(1, 7)
MARKED_START = 92
MARKED_CURVATURE = Q(4, 3)
MARKED_OUTER_SQUARE = Q(6, 35)
RAYS = {
    0: (91360, 65613, 30940, -13225, -49942, -80972),
    24: (89867, 64295, 25308, -15948, -54468, -83956),
    48: (87985, 62124, 20956, -17557, -56837, -85224),
    72: (85758, 60162, 17086, -19013, -58650, -85669),
    92: (83574, 57717, 13768, -21135, -60275, -86136),
    96: (83137, 57228, 13104, -21559, -60600, -86230),
    120: (80197, 55062, 10876, -22823, -61286, -85478),
    132: (79343, 53954, 10654, -24011, -62032, -85512),
}


def possible_rows(floor):
    """Return each row above the given reciprocal floor in the curvature band."""
    gap_bound = int(EXISTING_CURVATURE * (1 + 1 / (7 * floor)))
    for gap in range(1, gap_bound + 1):
        lower_bound = int(gap / (7 * floor) - Q(1, 7))
        for lower in range(1, lower_bound + 1):
            curvature = Q(gap**2, lower + gap)
            if Q(gap, 7 * lower + 1) <= floor or curvature >= EXISTING_CURVATURE:
                continue
            outer_square = prod(Q(7 * lower + h, 7 * lower + 6 * gap + h)
                                for h in range(1, gap + 1))
            yield lower, gap, curvature, outer_square


def profile_outer_square(core_square, core_curvature):
    # The individual loss is 3/7 for the log mean, or 6/7 for its square.
    other_loss = exponential_lower(Q(6, 7) * max(CURVATURE - core_curvature, 0))
    return core_square / other_loss


def marked_exclusion():
    threshold = Q(MARKED_START, 1000)**2
    for _, _, curvature, outer_square in possible_rows(Q(1, 8)):
        pair_cap = profile_outer_square(MARKED_OUTER_SQUARE * outer_square,
                                        MARKED_CURVATURE + curvature)
        require(pair_cap < threshold, ("second retained row after marking", curvature))
    require(threshold * exponential_lower(Q(45, 14)) > MARKED_OUTER_SQUARE,
            "marked remainder curvature below 15/4")


def reciprocal_caps(rows, outer_lower):
    lower = FLOOR_CAP
    upper = 7 * FLOOR_CAP / (8 - FLOOR_CAP)
    for x7, x8, u_squared_cap in rows:
        if outer_lower**2 <= u_squared_cap:
            lower = max(lower, x7)
            upper = max(upper, x8)
    return lower, upper


def skew_weights(order, k):
    return tuple((Q(k <= m < order-k) - Q(order - 2*k, order)) / 2
                 for m in range(order))


def verify_skew_convexity():
    for order in (7, 8):
        for k in range(1, (order + 1) // 2):
            first = sum(w * m / (1 + m*x) for m, w in enumerate(skew_weights(order, k)))
            derivative_scale = 1 + (order - 1)*x
            second = derivative_scale * sym.diff(derivative_scale * first, x)
            check_nonnegative_on_interval(second.subs(x, x/2))


def weighted_log_bounds(weights, cap):
    intervals = [multiply_interval((weight, weight), log_bounds(1 + m*cap))
                 for m, weight in enumerate(weights) if weight and m]
    if not intervals:
        return Q(0), Q(0)
    return tuple(map(sum, zip(*intervals, strict=True)))


@lru_cache(maxsize=None)
def exponent_bounds(cap7, cap8):
    denominator_lower, denominator_upper = log_bounds(1 + 6*cap7)
    rho_lower = log_bounds(1 + Q(49, 8)*cap7)[0] / denominator_upper
    rho_upper = 1 + log_bounds(1 / (1-cap7/8))[1] / denominator_lower
    require(rho_lower <= Q(49, 48) <= rho_upper, "Gaussian ratio")
    result = {}
    for order, cap, ratio in ((7, cap7, (1, 1)), (8, cap8, (rho_lower, rho_upper))):
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
            require(0 < bounds[-1][0] <= mean_upper * (1 if order == 7 else Q(49, 48))
                    <= bounds[-1][1], "Gaussian anchor envelope")
        result[order] = tuple(bounds)
    return result


def boxes(u, bounds, total_caps):
    put = {
        order: [(lo, min(hi, cap)) for (lo, hi), cap in
                zip(coefficient_box(order, u, bounds[order]), total_caps[order], strict=True)]
        for order in (7, 8)
    }
    lower = [hi / comb(7, k) for k, (_, hi) in enumerate(put[7])]
    for k in range(1, 8):
        lo, hi = put[8][k]
        left = lower[k] / exponential_lower(CURVATURE * k*k * Q(1, 120))
        right = lower[k-1] / exponential_lower(CURVATURE * (8-k)**2 * Q(1, 120))
        put[8][k] = lo, min(hi, comb(8, k) * min(left, right))
    return put


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
        ((8*c**4-8*c**2+1 for c in knots), (1, -1, -1, 1, 1, -1, -1, 1)),
        ((8*c**3-4*c**2-4*c+1 for c in knots), (1, 1, -1, -1, 1, 1, -1, -1)),
    ):
        require(all(value*sign > 0 for value, sign in zip(values, signs, strict=True)),
                "binomial anchor signs")


def verify():
    verify_skew_convexity()
    verify_anchor()
    knots = sorted(RAYS)
    require(Q(1, 112) > Q(1, 120), "Gaussian endpoint attenuation")
    require(exponential_lower(Q(3, 7) * CURVATURE) * Q(knots[-1], 1000) > 1,
            "outer mean cutoff")
    marked_exclusion()
    rows = tuple((Q(d, 7*a+1), Q(d, 8*a+1), profile_outer_square(v2, s))
                 for a, d, s, v2 in possible_rows(FLOOR_CAP))
    marked = (Q(1, 4), Q(2, 9),
              profile_outer_square(MARKED_OUTER_SQUARE, MARKED_CURVATURE))
    require(marked in rows, "marked row identity")
    unmarked = tuple(row for row in rows if row != marked)
    require(all(x7 <= 1 and x8 <= 1 for x7, x8, _ in rows), "convexity range")
    total_caps = {j: coefficient_bounds(j, (), CURVATURE) for j in (7, 8)}
    for j in (7, 8):
        caps = total_caps[j]
        require(sum(a * INNER**k for k, a in enumerate(caps) if k) < 1,
                ("small circle", j))
        require(sum(a * OUTER**k for k, a in enumerate(caps[:-1])) < OUTER**j,
                ("large circle", j))
    total = 0
    for left, right in zip(knots, knots[1:]):
        candidates = unmarked if left >= MARKED_START else rows
        cells = sum(verify_knot(j, i, left, right, candidates, total_caps)
                    for j in (7, 8) for i in range(6))
        total += cells
        print(f"[{left}/1000,{right}/1000]: {cells} radial cells.", flush=True)
    print(f"Verified {total} radial cells.", flush=True)


if __name__ == "__main__":
    verify()
