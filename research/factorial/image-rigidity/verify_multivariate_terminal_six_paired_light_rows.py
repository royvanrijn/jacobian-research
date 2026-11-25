#!/usr/bin/env python3
"""Certify terminal-six paired projections for shared rows with cap 1/8."""

from fractions import Fraction as Q
from math import comb

from verify_multivariate_terminal_six_diffuse_bands import (
    bernstein_coefficients,
    bisect_bernstein,
)
from verify_multivariate_terminal_six_joint_ray_bands import (
    check_defect,
    check_skew,
    log_bounds,
    power_bounds,
)
from verify_multivariate_terminal_six_private_angles import (
    multiply_interval,
    require,
    sine_ratio_intervals,
)


RAYS = {
    28: (439, 13, -481, -757), 30: (404, -15, -499, -766),
    32: (369, -44, -517, -775), 34: (333, -72, -534, -783),
    36: (296, -100, -551, -791), 38: (259, -129, -568, -799),
    40: (222, -157, -585, -807), 42: (185, -185, -602, -815),
    44: (147, -213, -618, -823), 46: (122, -242, -633, -830),
    48: (92, -270, -649, -838), 50: (56, -298, -665, -845),
    52: (23, -326, -680, -852), 54: (-20, -354, -695, -858),
    56: (-69, -382, -710, -865), 58: (-115, -410, -726, -872),
    60: (-166, -438, -740, -878),
}
INNER, OUTER = Q(1, 8), Q(8)
RATIO_LO, RATIO_HI, RATIO_ANCHOR = Q(1034, 1000), Q(1052, 1000), Q(25, 24)
QUINTIC_OUTER_SKEW = Q(5050, 100000)
QUINTIC_INNER_SKEW = Q(2510, 100000)
QUINTIC_LINK = Q(8, 10000)
QUINTIC_DEFECT = Q(509, 100000)
SEXTIC_OUTER_SKEW = Q(5666, 100000)
SEXTIC_INNER_SKEW = Q(4563, 100000)
SEXTIC_LINK = Q(147, 100000)
SEXTIC_DEFECT = Q(578, 100000)
SEXTIC_CENTRAL_LINK = Q(9, 100000)


def validate_envelopes():
    t5, t6 = Q(1, 8), Q(5, 47)
    a, b, c, d = (1 + k * t5 for k in range(1, 5))
    low_log5, high_log5 = log_bounds(d)
    require(log_bounds(1 + Q(25, 6) * t5)[0] >= RATIO_LO * high_log5,
            "lower outer-log ratio")
    require(log_bounds(1 / (1 - t5 / 6))[1] <= (RATIO_HI - 1) * low_log5,
            "upper outer-log ratio")
    require(RATIO_LO < RATIO_ANCHOR < RATIO_HI, "anchor outer-log ratio")
    check_skew((a * b * c)**2 / d**3, d, 10, QUINTIC_OUTER_SKEW)
    check_skew(b**4 / (a * c * d), d, 10, QUINTIC_INNER_SKEW)
    check_defect(Q(1, 2), d, QUINTIC_DEFECT)
    require(Q(6, 5) * t5**4 <= QUINTIC_LINK * low_log5, "quintic linked skew")

    a, b, c, d, e = (1 + k * t6 for k in range(1, 6))
    check_skew(a * b * c * d / e**2, e, 6, SEXTIC_OUTER_SKEW)
    check_skew((a * d)**3 / ((b * c)**2 * e), e, 2, SEXTIC_LINK)
    check_defect(Q(3, 5), e, SEXTIC_DEFECT)
    require((4 * SEXTIC_OUTER_SKEW + SEXTIC_LINK) / 5 <= SEXTIC_INNER_SKEW,
            "sextic inner skew")
    require(Q(12, 5) * t6**5 <= SEXTIC_CENTRAL_LINK * log_bounds(e)[0],
            "sextic central link")

    U = log_bounds(Q(25, 7))[1]
    require(max(QUINTIC_OUTER_SKEW, SEXTIC_OUTER_SKEW * RATIO_HI) * U < Q(2, 25),
            "outer skew bound for the root annulus")
    require(Q(15, 23) * ((1 + INNER)**6 - 1) + INNER**5 < 1, "root annulus")


def gaussian_exponents(order):
    return tuple(Q(k * (order - k), order - 1) for k in range(order + 1))


def coefficient_model(order, u):
    reference_ratio = Q(1) if order == 5 else RATIO_ANCHOR
    midpoint = (u[0] + u[1]) / 2
    exponents = gaussian_exponents(order)
    nominal = tuple(
        (comb(order, k) * lo, comb(order, k) * hi)
        for k, exponent in enumerate(exponents)
        for lo, hi in [power_bounds(midpoint, reference_ratio * exponent)
                       if exponent else (Q(1), Q(1))]
    )

    upper_log = log_bounds(1 / u[0])[1]
    mean_lo = -reference_ratio * log_bounds(midpoint / u[0])[1]
    mean_hi = reference_ratio * log_bounds(u[1] / midpoint)[1]
    if order == 5:
        variables = (
            ((0, 2, 1, -1, -2, 0), -QUINTIC_INNER_SKEW, QUINTIC_INNER_SKEW),
            ((0, 1, 0, 0, -1, 0), -QUINTIC_LINK, QUINTIC_LINK),
            ((0, 0, 1, 1, 0, 0), 0, QUINTIC_DEFECT),
        )
    else:
        mean_lo -= (RATIO_HI - reference_ratio) * upper_log
        mean_hi += (reference_ratio - RATIO_LO) * upper_log
        upper_log *= RATIO_HI
        variables = (
            ((0, Q(5, 4), 1, 0, -1, Q(-5, 4), 0), -SEXTIC_INNER_SKEW, SEXTIC_INNER_SKEW),
            ((0, 1, 0, 0, 0, -1, 0), -SEXTIC_LINK / 4, SEXTIC_LINK / 4),
            ((0, 0, 1, Q(3, 2), 1, 0, 0), 0, SEXTIC_DEFECT),
            ((0, 0, 0, -1, 0, 0, 0), 0, SEXTIC_CENTRAL_LINK),
        )

    pieces = tuple((vector, (lo * upper_log, hi * upper_log)) for vector, lo, hi in variables)
    pieces += ((exponents, (mean_lo, mean_hi)),)
    radii = tuple(
        sum(abs(vector[k]) * max(-lo, hi) for vector, (lo, hi) in pieces)
        for k in range(order + 1)
    )
    require(max(radii) < 1, ("exponential remainder", order, u))
    remainder = tuple(
        multiply_interval(coefficient, (Q(1), 1 + radius**2 / (2 * (1 - radius))))
        for coefficient, radius in zip(nominal, radii, strict=True)
    )
    basis = [((Q(1), Q(1)), remainder)]
    for vector, interval in pieces:
        basis.append((interval, tuple(
            multiply_interval(coefficient, (value, value))
            for coefficient, value in zip(nominal, vector, strict=True)
        )))
    return basis


def projected_bernstein(order, cosines, model):
    ratios = sine_ratio_intervals(order, cosines)
    projections = []
    for rotation in range(order + 1):
        directions = tuple(
            ratios[k - rotation] if k >= rotation
            else (-ratios[rotation - k][1], -ratios[rotation - k][0])
            for k in range(order + 1)
        )
        groups = []
        for _, basis in model:
            pairs = tuple(
                multiply_interval(coefficient, direction)
                for coefficient, direction in zip(basis, directions, strict=True)
            )
            groups.append(tuple(
                bernstein_coefficients(values, INNER, OUTER)
                for values in zip(*pairs, strict=True)
            ))
        projections.append(tuple(groups))
    return tuple(projections)


def projection_is_signed(projection, model):
    positive = negative = True
    group_entries = (tuple(zip(*group, strict=True)) for group in projection)
    for entries in zip(*group_entries, strict=True):
        bounds = (
            multiply_interval(interval, coefficient)
            for (interval, _), coefficient in zip(model, entries, strict=True)
        )
        lo, hi = (sum(side) for side in zip(*bounds, strict=True))
        positive = positive and lo > 0
        negative = negative and hi < 0
        if not (positive or negative):
            return False
    return True


def split_projection(projection):
    groups = []
    for lower, upper in projection:
        lower_left, lower_right = bisect_bernstein(lower)
        upper_left, upper_right = bisect_bernstein(upper)
        groups.append(((lower_left, upper_left), (lower_right, upper_right)))
    return tuple(zip(*groups, strict=True))


def radial_replay(order, cosines, model):
    stack = [(projected_bernstein(order, cosines, model), 0)]
    count = 0
    while stack:
        projections, depth = stack.pop()
        if any(projection_is_signed(p, model) for p in projections):
            count += 1
            continue
        if depth >= 13:
            return None
        children = zip(*(split_projection(p) for p in projections), strict=True)
        stack.extend((child, depth + 1) for child in children)
    return count


def ray_cosine(left, right, index, u):
    start, finish = RAYS[left][index], RAYS[right][index]
    return Q(start, 1000) + Q(finish - start, 1000) * (100 * u - left) / (right - left)


def verify_knot_interval(order, index, left, right):
    pending = [((Q(left, 100), Q(right, 100)), 0)]
    count = 0
    while pending:
        u, depth = pending.pop()
        cosines = tuple(sorted(ray_cosine(left, right, index, endpoint) for endpoint in u))
        used = radial_replay(order, cosines, coefficient_model(order, u))
        if used is not None:
            count += used
            continue
        require(depth < 10, ("unresolved ray interval", order, index, left, right, u))
        middle = (u[0] + u[1]) / 2
        pending.extend((((u[0], middle), depth + 1), ((middle, u[1]), depth + 1)))
    return count


def affine_interval(constant, terms):
    bounds = [multiply_interval((factor, factor), interval) for factor, interval in terms]
    return tuple(constant + sum(side) for side in zip(*bounds, strict=True))


def verify_anchor():
    for knots in RAYS.values():
        require(1000 > knots[0] > knots[1] > knots[2] > knots[3] > -1000, "ray order")
    require(gaussian_exponents(5) == (0, 1, Q(3, 2), Q(3, 2), 1, 0),
            "quintic reference exponents")
    require(gaussian_exponents(6) == (0, 1, Q(8, 5), Q(9, 5), Q(8, 5), 1, 0),
            "sextic reference exponents")
    u = Q(7, 25)
    v = power_bounds(u, Q(3, 2))
    w = power_bounds(u, RATIO_ANCHOR)
    y = power_bounds(u, Q(8, 5) * RATIO_ANCHOR)
    z = power_bounds(u, Q(9, 5) * RATIO_ANCHOR)

    def quintic_quotient(x):
        return affine_interval(x**2 + (5 * u - 1) * x - 5 * u - 1, ((10, v),))

    def sextic_quotient(x):
        return affine_interval(x**3 - 3 * x, ((6 * x**2 - 12, w), (15 * x, y), (20, z)))

    abscissae = tuple(Q(c, 500) for c in RAYS[28])
    quintic = tuple(quintic_quotient(x) for x in abscissae)
    sextic = tuple(sextic_quotient(x) for x in abscissae)
    require(quintic[0][0] > 0 > quintic[1][1]
            and quintic[2][1] < 0 < quintic[3][0], "quintic anchor sectors")
    require(sextic_quotient(2)[0] > 0 > sextic[0][1]
            and sextic[1][1] < 0 < sextic[2][0]
            and sextic[3][0] > 0 > sextic_quotient(-2)[1], "sextic anchor sectors")


def main():
    validate_envelopes()
    verify_anchor()
    knots = tuple(sorted(RAYS))
    total = 0
    for left, right in zip(knots, knots[1:]):
        count = sum(verify_knot_interval(order, index, left, right)
                    for order in (5, 6) for index in range(4))
        total += count
        print(f"[{left}/100,{right}/100]: {count} radial cells.", flush=True)
    print(f"All {total} radial cells certified.", flush=True)


if __name__ == "__main__":
    main()
