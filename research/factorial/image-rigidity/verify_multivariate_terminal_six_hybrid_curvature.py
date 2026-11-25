#!/usr/bin/env python3
"""Exact high-curvature exclusion combining the shared outer mean and total loss."""

from fractions import Fraction as Q
from math import comb, prod

from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER, OUTER, check_defect, check_skew, coefficient_box, log_bounds,
    radial_replay, rounded_interval,
)
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds, exponential_lower, require,
)


STRENGTH = 4
ATTENUATION_BASE = Q(65, 66)
RAY_KNOTS = {
    0: (8412, -1433), 10: (8311, -1753), 20: (8197, -1934),
    30: (8099, -2116), 40: (7981, -2289), 45: (7914, -2351),
    50: (7868, -2423), 60: (7748, -2590), 70: (7624, -2752),
    80: (7482, -2911), 90: (7348, -3024), 100: (7210, -3177),
    101: (7196, -3230), 110: (7101, -3264), 120: (6944, -3412),
    125: (6856, -3465), 130: (6768, -3558), 135: (6679, -3610),
    140: (6588, -3702), 145: (6497, -3700), 150: (6405, -3735),
    155: (6314, -3710), 160: (6259, -3740), 165: (6189, -3740),
    170: (6134, -3740), 175: (6072, -3730), 180: (6024, -3720),
    185: (5982, -3677), 190: (5956, -3660), 195: (5930, -3684),
    200: (5890, -3710), 205: (5850, -3740),
}
PARAMETER_CAPS = (
    (0, 20, Q(1), Q(6, 7), Q(55, 48), (18986, 8768, 6933, 20993, 6808, 7806, 3518, 102037)),
    (20, 45, Q(5, 6), Q(5, 7), Q(25, 22), (17465, 8158, 5902, 19338, 5401, 6644, 3036, 102151)),
    (45, 101, Q(2, 3), Q(4, 7), Q(9, 8), (15634, 7396, 4761, 17344, 3959, 5361, 2489, 102297)),
    (101, 205, Q(1, 2), Q(3, 7), Q(10, 9), (13362, 6413, 3503, 14865, 2537, 3950, 1868, 102493)),
)


def check_exponents(x5, x6, upper_ratio, data):
    alpha, chi, defect, beta, link, defect_y, defect_z, lower_ratio = (
        Q(x, 100000) for x in data
    )
    a, b, c, d = (1 + k * x5 for k in range(1, 5))
    require(log_bounds(1 + Q(25, 24) * 4 * x5)[0]
            > lower_ratio * log_bounds(d)[1], "outer mean transfer")
    check_skew((a * b * c)**2 / d**3, d, 10, alpha)
    check_skew(b**4 / (a * c * d), d, 10, chi)
    check_defect(Q(1, 2), d, defect)
    a, b, c, d, e = (1 + k * x6 for k in range(1, 6))
    check_skew(a * b * c * d / e**2, e, 6, beta)
    check_skew((a * d)**3 / ((b * c)**2 * e), e, 2, link)
    check_defect(Q(3, 5), e, defect_y)
    check_defect(Q(1, 5), e, defect_z)
    eta = (4 * beta + link) / 5
    outer_five = rounded_interval(1 - alpha, 1 + alpha)
    inner_five = rounded_interval(Q(3, 2) - defect - chi, Q(3, 2) + chi)
    outer_six = rounded_interval(lower_ratio * (1 - beta), upper_ratio * (1 + beta))
    inner_six = rounded_interval(lower_ratio * (Q(8, 5) - defect_y - eta),
                                 upper_ratio * (Q(8, 5) + eta))
    central = rounded_interval(lower_ratio * (Q(9, 5) - defect_y - defect_z),
                               upper_ratio * Q(9, 5))
    return {
        5: (outer_five, inner_five, inner_five, outer_five),
        6: (outer_six, inner_six, central, inner_six, outer_six),
    }


def rays_at(knot):
    first, third = RAY_KNOTS[knot]
    return Q(first, 10000), Q(2, 5), Q(third, 10000), Q(-7, 10)


def interpolate_ray(left, right, index, value):
    t = (1000 * value - left) / (right - left)
    return (1 - t) * rays_at(left)[index] + t * rays_at(right)[index]


def check_knot(order, index, left, right, exponent_bounds, total_caps):
    stack = [((Q(left, 1000), Q(right, 1000)), 0)]
    accepted = 0
    while stack:
        u, depth = stack.pop()
        coefficients = [(lo, min(hi, cap)) for (lo, hi), cap in zip(
            coefficient_box(order, u, exponent_bounds[order]), total_caps[order], strict=True
        )]
        if order == 6:
            # Apply the fixed endpoint-distance attenuation from both ends.
            previous = [
                min(hi, cap) / comb(5, k)
                for k, ((_, hi), cap) in enumerate(zip(
                    coefficient_box(5, u, exponent_bounds[5]), total_caps[5], strict=True
                ))
            ]
            for k in range(1, 6):
                lower, upper = coefficients[k]
                left_cap = previous[k] * ATTENUATION_BASE ** (STRENGTH * k * k)
                right_cap = previous[k - 1] * ATTENUATION_BASE ** (STRENGTH * (6 - k) ** 2)
                coefficients[k] = lower, min(upper, comb(6, k) * min(left_cap, right_cap))
        if any(lo > hi for lo, hi in coefficients):
            continue
        cosine = tuple(sorted(interpolate_ray(left, right, index, end) for end in u))
        radial = radial_replay(order, coefficients, cosine)
        if radial is not None:
            accepted += radial
            continue
        require(depth < 16, ("unresolved hybrid ray", order, index, u))
        midpoint = (u[0] + u[1]) / 2
        stack.append(((u[0], midpoint), depth + 1))
        stack.append(((midpoint, u[1]), depth + 1))
    return accepted


def verify():
    for a, d, cutoff in ((1, 6, 20), (1, 5, 45), (1, 4, 101), (2, 6, 101)):
        require(prod(Q(5 * a + h, 5 * a + 4 * d + h) for h in range(1, d + 1))
                < Q(cutoff, 1000)**2, ("integer row exclusion", a, d, cutoff))
    bands = [(a, b, check_exponents(x5, x6, ratio, data))
             for a, b, x5, x6, ratio, data in PARAMETER_CAPS]
    last = max(RAY_KNOTS)
    require(exponential_lower(Q(2, 5) * STRENGTH) * Q(last, 1000) > 1,
            "rigorous outer mean cutoff")
    for knot in RAY_KNOTS:
        a, b, c, d = rays_at(knot)
        require(1 > a > b > c > d > -1, ("ray order", knot))
    a, b, c, d = rays_at(0)
    require((4 * a - 1)**2 > 5 and a*a < Q(3, 4) and 0 < b < Q(4, 5)
            and c < 0 and (1 - 4 * c)**2 < 5
            and (1 - 4 * d)**2 > 5 and d*d < Q(3, 4), "binomial anchor")
    caps = {order: coefficient_bounds(order, (), STRENGTH) for order in (5, 6)}
    for order in (5, 6):
        cap = caps[order]
        require(sum(c * INNER**k for k, c in enumerate(cap) if k) < 1,
                ("small circle", order))
        require(sum(c * OUTER**k for k, c in enumerate(cap[:-1])) < OUTER**order,
                ("large circle", order))
    knots = sorted(RAY_KNOTS)
    for order in (5, 6):
        counts = []
        for ray in range(4):
            counts.append(sum(
                check_knot(order, ray, a, b, exponent_bounds, caps)
                for band_left, band_right, exponent_bounds in bands
                for a, b in zip(knots, knots[1:])
                if band_left <= a < b <= band_right
            ))
        print(f"Order {order}: {counts} accepted radial intervals.", flush=True)


if __name__ == "__main__":
    verify()
