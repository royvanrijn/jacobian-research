#!/usr/bin/env python3
"""A rational moving-ray certificate for terminal six with one private side."""

from fractions import Fraction as Q
from functools import lru_cache

from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds, multiply_interval, require, sine_ratio_intervals,
)
from verify_multivariate_terminal_six_small_shared import CASES, check_anchor


P2 = (
    Q(1, 2), Q(3, 7), Q(10, 9),
    Q(7, 50), Q(7, 100), Q(1, 25),
    Q(3, 20), Q(3, 100), Q(13, 100), Q(1, 25), Q(1, 50),
)
P1_LOW = (
    Q(5, 6), Q(5, 7), Q(25, 22),
    Q(9, 50), Q(17, 200), Q(3, 50),
    Q(1, 5), Q(11, 200), Q(9, 50), Q(7, 100), Q(1, 25),
)
P1_HIGH = (
    Q(1, 3), Q(1, 3), Q(35, 32),
    Q(11, 100), Q(11, 200), Q(1, 40),
    Q(13, 100), Q(17, 1000), Q(11, 100), Q(3, 100), Q(3, 200),
)


def ratio_bound(argument, base, root_degree, limit):
    require(
        argument ** (2 * limit.denominator) <= base ** (root_degree * limit.numerator),
        ("skew constant", argument, base, root_degree, limit),
    )


def defect_bound(parameter, base, limit):
    required = parameter - limit
    ratio = ((1 + parameter) * base + 1 - parameter) / ((1 - parameter) * base + 1 + parameter)
    require(ratio**required.denominator >= base**required.numerator,
            ("defect constant", parameter, base, limit))


def coefficient_exponents(profile):
    x5, x6, order_ratio, alpha, chi, quintic_defect, beta, skew_link, eta, first_defect, second_defect = profile
    a, b, c, d = (1 + k * x5 for k in range(1, 5))
    ratio_bound((a * b * c)**2 / d**3, d, 10, alpha)
    ratio_bound(b**4 / (a * c * d), d, 10, chi)
    defect_bound(Q(1, 2), d, quintic_defect)
    a, b, c, d, e = (1 + k * x6 for k in range(1, 6))
    ratio_bound(a * b * c * d / e**2, e, 6, beta)
    ratio_bound((a * d)**3 / ((b * c)**2 * e), e, 2, skew_link)
    require((4 * beta + skew_link) / 5 <= eta, ("inner sextic skew", profile))
    defect_bound(Q(3, 5), e, first_defect)
    defect_bound(Q(1, 5), e, second_defect)
    outer_five = (1 - alpha, 1 + alpha)
    inner_five = (Q(3, 2) - quintic_defect - chi, Q(3, 2) + chi)
    outer_six = (1 - beta, order_ratio * (1 + beta))
    inner_six = (Q(8, 5) - first_defect - eta, order_ratio * (Q(8, 5) + eta))
    central_six = (Q(9, 5) - first_defect - second_defect, order_ratio * Q(9, 5))
    return {
        5: (outer_five, inner_five, inner_five, outer_five),
        6: (outer_six, inner_six, central_six, inner_six, outer_six),
    }


def ray_cosines(private_rows, u):
    if private_rows == (1,):
        if u <= Q(1, 5):
            first = Q(84, 100) - Q(525, 1000) * u
        elif u <= Q(3, 8):
            first = Q(735, 1000) - Q(26, 35) * (u - Q(1, 5))
        else:
            first = Q(605, 1000) - Q(736, 1000) * (u - Q(3, 8))
        third = (Q(-15, 100) - Q(64, 100) * u if u <= Q(3, 8)
                 else Q(-39, 100) - Q(352, 1000) * (u - Q(3, 8)))
        return first, Q(55, 100) - Q(82, 100) * u, third, Q(-7, 10) - Q(16, 100) * u**2
    if private_rows == (1, 1):
        return (Q(84, 100) - Q(26, 100) * u, Q(6, 10) - Q(1, 2) * u,
                Q(-15, 100) - Q(26, 100) * u, Q(-73, 100) - Q(7, 100) * u)
    return (Q(84, 100) - Q(19, 100) * u, Q(6, 10) - Q(49, 100) * u,
            Q(-15, 100) - Q(21, 100) * u, Q(-73, 100) - Q(63, 1000) * u)


@lru_cache(maxsize=None)
def power_interval(argument, exponent):
    if argument in (0, 1):
        return argument, argument
    scale = 10000
    rhs = scale**exponent.denominator * argument.numerator**exponent.numerator
    factor = argument.denominator**exponent.numerator
    lower, upper = 0, scale
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if middle**exponent.denominator * factor <= rhs:
            lower = middle
        else:
            upper = middle
    return Q(lower, scale), Q(upper, scale)


def verify_moving_ray(order, private_rows, index, exponents, u_range):
    reference = coefficient_bounds(order, private_rows, Q(0))
    low_reference = (Q(1),) + tuple(c - Q(1, 10000) for c in reference[1:-1]) + (Q(1),)
    require(sum(c / 5**k for k, c in enumerate(reference) if k) < 1,
            ("inner-circle domination", private_rows, order))
    require(sum(c * 5**k for k, c in enumerate(reference[:-1])) < 5**order,
            ("outer-circle domination", private_rows, order))
    radius_range = (Q(1, 5), Q(5))
    stack = [(u_range, radius_range, 0)]
    boxes = 0
    while stack:
        u, r, depth = stack.pop()
        boxes += 1
        cosines = (ray_cosines(private_rows, u[1])[index],
                   ray_cosines(private_rows, u[0])[index])
        ratios = sine_ratio_intervals(order, cosines)
        powers = [(r[0] ** k, r[1] ** k) for k in range(order + 1)]
        multipliers = [(Q(1), Q(1))]
        for k, (lo, hi) in enumerate(exponents, 1):
            multipliers.append((
                low_reference[k] * power_interval(u[0], hi)[0],
                reference[k] * power_interval(u[1], lo)[1],
            ))
        multipliers.append((Q(1), Q(1)))
        coefficients = [multiply_interval(value, power)
                        for value, power in zip(multipliers, powers, strict=True)]
        for rotation in range(order + 1):
            lower = upper = Q(0)
            for k, coefficient in enumerate(coefficients):
                if k >= rotation:
                    sign = ratios[k - rotation]
                else:
                    lo, hi = ratios[rotation - k]
                    sign = (-hi, -lo)
                value = multiply_interval(sign, coefficient)
                lower += value[0]
                upper += value[1]
            if lower > 0 or upper < 0:
                break
        else:
            require(depth < 28, ("unresolved moving ray", private_rows, order, index, u, r))
            if (u[1] - u[0]) * (radius_range[1] - radius_range[0]) > (
                r[1] - r[0]
            ) * (u_range[1] - u_range[0]):
                midpoint = (u[0] + u[1]) / 2
                stack.append(((u[0], midpoint), r, depth + 1))
                stack.append(((midpoint, u[1]), r, depth + 1))
            else:
                midpoint = (r[0] + r[1]) / 2
                stack.append((u, (r[0], midpoint), depth + 1))
                stack.append((u, (midpoint, r[1]), depth + 1))
    return boxes


def verify_private_rows(private_rows, slabs, exponents):
    cosines = ray_cosines(private_rows, Q(1))
    fifth_centers, sixth_centers = CASES[private_rows][2:]
    check_anchor(5, private_rows, cosines, fifth_centers)
    check_anchor(6, private_rows, cosines, sixth_centers)
    boxes = 0
    for lo, hi, name in slabs:
        for u in (lo, hi):
            first, second, third, fourth = ray_cosines(private_rows, u)
            require(1 > first > second > third > fourth > -1, ("angle ordering", private_rows, u))
        for order in (5, 6):
            for index in range(4):
                boxes += verify_moving_ray(order, private_rows, index, exponents[name][order], (lo, hi))
    return boxes


if __name__ == "__main__":
    exponents = {name: coefficient_exponents(data)
                 for name, data in (("P2", P2), ("P1 low", P1_LOW), ("P1 high", P1_HIGH))}
    for private_rows, slabs in (
        ((1,), (
            (Q(0), Q(1, 5), "P1 low"),
            (Q(1, 5), Q(3, 8), "P1 low"),
            (Q(3, 8), Q(1), "P1 high"),
        )),
        ((1, 1), ((Q(0), Q(1), "P2"),)),
        ((2,), ((Q(0), Q(1), "P2"),)),
    ):
        print(f"Private rows {private_rows}: {verify_private_rows(private_rows, slabs, exponents)} rectangles.",
              flush=True)
