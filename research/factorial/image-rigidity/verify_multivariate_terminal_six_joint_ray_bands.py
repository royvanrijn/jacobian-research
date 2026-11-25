#!/usr/bin/env python3
"""Exact fully shared terminal-six bands tied by one common outer mean."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb

from verify_multivariate_terminal_six_diffuse_bands import (
    bernstein_coefficients, bisect_bernstein,
)
from verify_multivariate_terminal_six_private_angles import (
    multiply_interval, require, sine_ratio_intervals,
)


# t, followed by bounds for quintic outer and inner skew, quintic defect,
# sextic outer and linked skew, sextic defects, and the lower order log ratio.
BANDS = (
    (0, 28, Q(1, 8), (5050, 2510, 509, 5666, 147, 578, 287, 103401)),
    (28, 36, Q(2, 25), (3465, 1728, 240, 3874, 47, 271, 135, 103620)),
    (36, 44, Q(11, 200), (2484, 1240, 124, 2772, 18, 139, 69, 103764)),
    (44, 50, Q(7, 200), (1638, 819, 54, 1825, 5, 60, 30, 103895)),
    (50, 55, Q(23, 1000), (1100, 550, 25, 1225, 2, 27, 14, 103981)),
    (55, 60, Q(3, 200), (729, 365, 11, 811, 1, 12, 6, 104042)),
    (60, 65, Q(9, 1000), (443, 222, 4, 492, 1, 5, 3, 104090)),
    (65, 70, Q(1, 200), (248, 124, 2, 276, 1, 2, 1, 104123)),
    (70, 72, Q(7, 2000), (174, 87, 1, 194, 1, 1, 1, 104136)),
)

# Values of the continuous affine rays at their rational knots.
RAYS = {
    0: (840, 400, -150, -730), 5: (795, 385, -222, -739),
    10: (731, 356, -291, -743), 15: (658, 317, -356, -702),
    20: (573, 251, -415, -724), 25: (484, 176, -472, -741),
    28: (429, 130, -504, -744), 30: (400, 100, -516, -754),
    35: (308, -55, -567, -762), 36: (290, -70, -577, -763),
    40: (221, -108, -607, -776), 44: (146, -148, -647, -777),
    45: (131, -163, -646, -794), 50: (36, -221, -691, -796),
    55: (-63, -293, -728, -816), 60: (-165, -352, -767, -836),
    65: (-267, -437, -801, -857), 70: (-370, -539, -831, -882),
    72: (-411, -566, -843, -893),
}
INNER, OUTER = Q(1, 8), Q(8)


def log_bounds(value):
    require(value >= 1, ("log bound domain", value))
    argument = (value - 1) / (value + 1)
    lower = 2 * sum(argument ** (2 * k + 1) / (2 * k + 1) for k in range(16))
    upper = lower + 2 * argument**33 / (33 * (1 - argument**2))
    return lower, upper


def check_skew(value, base, factor, limit):
    require(2 * log_bounds(value)[1] <= limit * factor * log_bounds(base)[0],
            ("local skew", value, base, factor, limit))


def check_defect(parameter, base, limit):
    ratio = ((1 + parameter) * base + 1 - parameter) / ((1 - parameter) * base + 1 + parameter)
    require(log_bounds(ratio)[0] >= (parameter - limit) * log_bounds(base)[1],
            ("local defect", parameter, base, limit))


def rounded_interval(lower, upper):
    size = 10000
    return Q(lower.numerator * size // lower.denominator, size), Q(
        -((-upper.numerator * size) // upper.denominator), size
    )


def exponents(cap, values):
    alpha, chi, defect, beta, linked, defect_y, defect_z, log_lower = (
        Q(value, 100000) for value in values
    )
    upper_x6 = 5 * cap / (6 - cap)
    log_upper = Q(25, 4) / (6 - cap)
    a, b, c, d = (1 + k * cap for k in range(1, 5))
    require(log_bounds(1 + Q(25, 24) * 4 * cap)[0]
            >= log_lower * log_bounds(d)[1], ("order ratio", cap))
    check_skew((a * b * c)**2 / d**3, d, 10, alpha)
    check_skew(b**4 / (a * c * d), d, 10, chi)
    check_defect(Q(1, 2), d, defect)
    a, b, c, d, e = (1 + k * upper_x6 for k in range(1, 6))
    check_skew(a * b * c * d / e**2, e, 6, beta)
    check_skew((a * d)**3 / ((b * c)**2 * e), e, 2, linked)
    check_defect(Q(3, 5), e, defect_y)
    check_defect(Q(1, 5), e, defect_z)
    eta = (4 * beta + linked) / 5
    outer_five = rounded_interval(1 - alpha, 1 + alpha)
    inner_five = rounded_interval(Q(3, 2) - defect - chi, Q(3, 2) + chi)
    outer_six = rounded_interval(log_lower * (1 - beta), log_upper * (1 + beta))
    inner_six = rounded_interval(log_lower * (Q(8, 5) - defect_y - eta),
                                 log_upper * (Q(8, 5) + eta))
    central = rounded_interval(log_lower * (Q(9, 5) - defect_y - defect_z),
                               log_upper * Q(9, 5))
    return {
        5: (outer_five, inner_five, inner_five, outer_five),
        6: (outer_six, inner_six, central, inner_six, outer_six),
    }


@lru_cache(maxsize=None)
def power_bounds(value, exponent):
    if value in (0, 1):
        return value, value
    scale = 10**8
    p, q = exponent.numerator, exponent.denominator
    rhs = scale**q * value.numerator**p
    factor = value.denominator**p
    # This proposal only speeds up the exact integer comparisons below.
    approx = int(float(value)**float(exponent) * scale)
    lo, hi = max(0, approx - 1), min(scale, approx + 2)
    if lo**q * factor <= rhs <= hi**q * factor:
        return Q(lo, scale), Q(hi, scale)
    lo, hi = 0, scale
    while lo + 1 < hi:
        middle = (lo + hi) // 2
        if middle**q * factor <= rhs:
            lo = middle
        else:
            hi = middle
    return Q(lo, scale), Q(hi, scale)


def coefficient_box(order, u, bounds):
    coefficients = [(Q(1), Q(1))]
    for k, (low_exponent, high_exponent) in enumerate(bounds, 1):
        coefficients.append((comb(order, k) * power_bounds(u[0], high_exponent)[0],
                             comb(order, k) * power_bounds(u[1], low_exponent)[1]))
    return coefficients + [(Q(1), Q(1))]


def cosine_at(left, right, index, value):
    x, y = RAYS[left][index], RAYS[right][index]
    return Q(x, 1000) + Q(y - x, 1000) * (100 * value - left) / (right - left)


def radial_replay(order, coefficients, cosines):
    ratios = sine_ratio_intervals(order, cosines)
    projections = []
    for k in range(order + 1):
        directions = [ratios[j - k] if j >= k else (-ratios[k - j][1], -ratios[k - j][0])
                      for j in range(order + 1)]
        products = [multiply_interval(a, b)
                    for a, b in zip(coefficients, directions, strict=True)]
        for p in ([x[0] for x in products], [-x[1] for x in products]):
            projections.append(bernstein_coefficients(p, INNER, OUTER))
    stack = [(projections, 0)]
    count = 0
    while stack:
        candidates, depth = stack.pop()
        count += 1
        if any(min(item) > 0 for item in candidates):
            continue
        possible = [item for item in candidates if max(item) > 0]
        if not possible or depth == 14:
            return None
        halves = [bisect_bernstein(item) for item in possible]
        stack.extend((child, depth + 1) for child in zip(*halves, strict=True))
    return count


def verify_knot_interval(order, index, left, right, bounds):
    u = Q(left, 100), Q(right, 100)
    coefficients = coefficient_box(order, u, bounds)
    require(sum(b * INNER**k for k, (_, b) in enumerate(coefficients) if k) < 1,
            ("small circle", order, left, right))
    require(sum(b * OUTER**k for k, (_, b) in enumerate(coefficients[:-1])) < OUTER**order,
            ("large circle", order, left, right))
    count = 0
    stack = [(u, 0)]
    while stack:
        u, depth = stack.pop()
        cosines = tuple(sorted(cosine_at(left, right, index, end) for end in u))
        coefficients = coefficient_box(order, u, bounds)
        radial = radial_replay(order, coefficients, cosines)
        if radial is not None:
            count += radial
            continue
        require(depth < 15, ("unresolved outer-mean interval", order, index, left, right, u))
        midpoint = (u[0] + u[1]) / 2
        stack.append(((u[0], midpoint), depth + 1))
        stack.append(((midpoint, u[1]), depth + 1))
    return count


def check_anchor():
    for knots in RAYS.values():
        require(1000 > knots[0] > knots[1] > knots[2] > knots[3] > -1000, "ray order")
    x1, x2, x3, x4 = (Q(value, 1000) for value in RAYS[0])
    # cos(pi/5)=(1+sqrt(5))/4, cos(3pi/5)=(1-sqrt(5))/4.
    require((4 * x1 - 1)**2 > 5 and x1**2 < Q(3, 4) and 0 < x2 < Q(4, 5)
            and x3 < 0 and (1 - 4 * x3)**2 < 5
            and (1 - 4 * x4)**2 > 5 and x4**2 < Q(3, 4), "binomial anchor")


if __name__ == "__main__":
    check_anchor()
    for left, right, cap, values in BANDS:
        bounds = exponents(cap, values)
        knots = [k for k in RAYS if left <= k <= right]
        boxes = sum(verify_knot_interval(order, index, a, b, bounds[order])
                    for a, b in zip(knots, knots[1:])
                    for order in (5, 6) for index in range(4))
        print(f"[{left}/100,{right}/100], Gauss cap {cap}: {boxes} radial intervals.",
              flush=True)
