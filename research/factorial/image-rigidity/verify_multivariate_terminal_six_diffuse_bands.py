#!/usr/bin/env python3
"""Exact angular certificates for fully shared profiles at terminal six."""

from fractions import Fraction
from math import comb, factorial

from verify_multivariate_terminal_six_private_angles import (
    exponential_lower, rational_cap, require, sine_ratios,
)


ANCHOR_COSINES = {
    Fraction(5, 4): (-181, -429, -760, -861),
    Fraction(3, 2): (-60, -350, -708, -850),
    Fraction(2): (135, -210, -640, -807),
    Fraction(3): (402, -20, -502, -768),
    Fraction(4): (577, 120, -403, -721),
    Fraction(5): (674, 210, -330, -680),
    Fraction(6): (740, 278, -280, -648),
}

BANDS = (
    (Fraction(5, 4), Fraction(3, 2), Fraction(1, 500)),
    (Fraction(3, 2), Fraction(2), Fraction(1, 100)),
    (Fraction(2), Fraction(3), Fraction(1, 50)),
    (Fraction(3), Fraction(4), Fraction(1, 25)),
    (Fraction(4), Fraction(5), Fraction(1, 25)),
    (Fraction(5), Fraction(6), Fraction(1, 25)),
)

INNER_RADIUS = Fraction(1, 8)
OUTER_RADIUS = Fraction(8)


def bernstein_coefficients(polynomial, left, right):
    degree = len(polynomial) - 1
    length = right - left
    transformed = [
        length**k * sum(polynomial[j] * comb(j, k) * left**(j - k)
                        for j in range(k, degree + 1))
        for k in range(degree + 1)
    ]
    return tuple(
        sum(transformed[k] * Fraction(comb(i, k), comb(degree, k))
            for k in range(i + 1))
        for i in range(degree + 1)
    )


def bisect_bernstein(coefficients):
    row = coefficients
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(x + y) / 2 for x, y in zip(row, row[1:])]
        left.append(row[0])
        right.append(row[-1])
    return tuple(left), tuple(reversed(right))


def verify_ray(order, lower, upper, cosine):
    ratios = sine_ratios(order, cosine)
    projections = []
    for rotation in range(order + 1):
        signs = [
            ratios[k - rotation] if k >= rotation else -ratios[rotation - k]
            for k in range(order + 1)
        ]
        for overall in (1, -1):
            coefficients = [
                min(overall * direction * lo, overall * direction * hi)
                for direction, lo, hi in zip(signs, lower, upper, strict=True)
            ]
            projections.append(bernstein_coefficients(coefficients, INNER_RADIUS, OUTER_RADIUS))

    stack = [(projections, 0)]
    boxes = 0
    while stack:
        candidates, depth = stack.pop()
        boxes += 1
        if any(min(item) > 0 for item in candidates):
            continue
        possible = [item for item in candidates if max(item) > 0]
        require(possible and depth < 18, ("unresolved ray", order, cosine, depth))
        halves = [bisect_bernstein(item) for item in possible]
        stack.extend((children, depth + 1) for children in zip(*halves, strict=True))
    return boxes


def exp_negative_bounds(value):
    require(0 <= value < 66, ("exponential tail", value))
    lower = exponential_lower(value)
    upper = lower + value**65 / factorial(65) / (1 - value / 66)
    return 1 / upper, 1 / lower


def coefficient_box(order, left, right, relative_gap):
    lower, upper = [Fraction(1)], [Fraction(1)]
    for k in range(1, order):
        binomial_power = Fraction(comb(order, k)**order)
        weight = Fraction(k * (order - k), 2)
        lo = exp_negative_bounds(weight * right / (1 - relative_gap))[0]
        hi = exp_negative_bounds(weight * left)[1]
        lower.append(max(Fraction(), rational_cap(binomial_power * lo, order) - Fraction(1, 10000)))
        upper.append(rational_cap(binomial_power * hi, order))
    return lower + [Fraction(1)], upper + [Fraction(1)]


def check_anchor(order, strength, cosines):
    coefficients = []
    for k in range(1, order // 2 + 1):
        lo, hi = exp_negative_bounds(strength * Fraction(k * (order - k), 2 * order))
        coefficients.append((comb(order, k) * lo, comb(order, k) * hi))
    signs = (1, -1, -1, 1) if order == 5 else (-1, -1, 1, 1)
    for cosine, sign in zip(cosines, signs, strict=True):
        x = 2 * cosine
        if order == 5:
            constant, weights = x * x - x - 1, (x - 1, 1)
        else:
            constant, weights = x**3 - 3 * x, (x * x - 2, x, 1)
        bound = sign * constant + sum(
            min(sign * w * lo, sign * w * hi)
            for w, (lo, hi) in zip(weights, coefficients, strict=True)
        )
        require(bound > 0, ("anchor", order, strength, cosine))
    if order == 6:
        require(-2 + 2 * coefficients[0][1] - 2 * coefficients[1][0] + coefficients[2][1] < 0,
                ("negative real endpoint", strength))


def verify_band(left, right, relative_gap):
    count = (right - left) * 200
    require(count.denominator == 1, "subinterval cover")
    cosines_left, cosines_right = ANCHOR_COSINES[left], ANCHOR_COSINES[right]
    boxes = 0
    for k in range(count.numerator):
        lo = left + Fraction(k, 200)
        hi = lo + Fraction(1, 200)
        midpoint = (lo + hi) / 2
        t = (midpoint - left) / (right - left)
        cosines = tuple(
            Fraction(x, 1000) + Fraction(y - x, 1000) * t
            for x, y in zip(cosines_left, cosines_right, strict=True)
        )
        require(1 > cosines[0] > cosines[1] > cosines[2] > cosines[3] > -1,
                ("separating ray order", lo, cosines))
        for order in (5, 6):
            check_anchor(order, midpoint, cosines)
            lower, upper = coefficient_box(order, lo, hi, relative_gap)
            require(sum(c * INNER_RADIUS**i for i, c in enumerate(upper) if i) < 1,
                    ("inner-circle domination", order, lo))
            require(sum(c * OUTER_RADIUS**i for i, c in enumerate(upper[:-1])) < OUTER_RADIUS**order,
                    ("outer-circle domination", order, lo))
            boxes += sum(verify_ray(order, lower, upper, cosine) for cosine in cosines)
    return boxes


if __name__ == "__main__":
    for left, right, relative_gap in BANDS:
        count = verify_band(left, right, relative_gap)
        print(f"[{left}, {right}], relative gap at most {relative_gap}: {count} radial intervals")
