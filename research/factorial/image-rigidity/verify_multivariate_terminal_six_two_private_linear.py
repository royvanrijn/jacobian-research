#!/usr/bin/env python3
"""Exact terminal-six angle separation with one private linear row per side."""

from fractions import Fraction
from math import factorial

from verify_multivariate_terminal_six_private_angles import exponential_lower, require
from verify_multivariate_terminal_six_small_shared import verify_ray


SHARED_LIMIT = Fraction(19, 5)
LOSS_FACTOR = Fraction(13, 6)
SIZE_LIMIT = 50
BAND_COUNT = 152


def verify_row_loss():
    # If d/l >= 1/3 and l >= 51, then d^2/(l+d) >= l/12.
    require(Fraction(SIZE_LIMIT + 1, 12) > SHARED_LIMIT, "unretained row bound")
    # d^2/(50+d) increases with d and excludes all d >= 16.
    require(Fraction(16**2, SIZE_LIMIT + 16) > SHARED_LIMIT, "row gap cutoff")
    counts = {5: 0, 6: 0}
    for base in range(1, SIZE_LIMIT + 1):
        for gap in range(1, 16):
            strength = Fraction(gap**2, base + gap)
            if strength > SHARED_LIMIT:
                continue
            for order in counts:
                const = factorial(order * base)
                lead = factorial(order * (base + gap))
                for k in range(1, order):
                    reciprocal_power = Fraction(
                        const ** (order - k) * lead**k,
                        factorial(order * base + k * gap) ** order,
                    )
                    exponent = LOSS_FACTOR * strength * k * (order - k) / 2
                    require(reciprocal_power <= exponential_lower(exponent),
                            ("row loss", base, gap, order, k))
                    counts[order] += 1
    return counts


def exp_negative_bounds(value):
    require(0 <= value < 66, ("exponential tail", value))
    lower = exponential_lower(value)
    upper = lower + value**65 / factorial(65) / (1 - value / 66)
    return 1 / upper, 1 / lower


def coefficient_box(order, strength_lo, strength_hi):
    lower, upper = [], []
    for k in range(order + 1):
        weight = Fraction(k * (order - k), 2 * order)
        lower.append(exp_negative_bounds(LOSS_FACTOR * strength_hi * weight)[0])
        upper.append(exp_negative_bounds(strength_lo * weight)[1])
    return lower, upper


def cosines_at(strength):
    first = Fraction(14, 25) + Fraction(3, 8) * strength / (strength + Fraction(3, 2))
    return (first, Fraction(1, 5), Fraction(-27, 100), Fraction(-13, 20))


def check_anchor(order, strength, cosines):
    coefficients = [
        exp_negative_bounds(strength * Fraction(k * (order - k), 2 * order))
        for k in range(1, order // 2 + 1)
    ]
    signs = (1, -1, -1, 1) if order == 5 else (-1, -1, 1, 1)
    for sign, cosine in zip(signs, cosines, strict=True):
        x = 2 * cosine
        if order == 5:
            constant, factors = x**2 - x - 1, (x - 1, 1)
        else:
            constant, factors = x**3 - 3 * x, (x**2 - 2, x, 1)
        lower = sign * constant + sum(
            min(sign * a * lo, sign * a * hi)
            for a, (lo, hi) in zip(factors, coefficients, strict=True)
        )
        require(lower > 0, ("anchor sign", order, strength, cosine))


def verify_all_bands():
    radial_boxes = {5: 0, 6: 0}
    require(SHARED_LIMIT / BAND_COUNT == Fraction(1, 40), "curvature cover")
    for band in range(BAND_COUNT):
        lo, hi = Fraction(band, 40), Fraction(band + 1, 40)
        midpoint = (lo + hi) / 2
        cosines = cosines_at(midpoint)
        for order in radial_boxes:
            check_anchor(order, midpoint, cosines)
            lower, upper = coefficient_box(order, lo, hi)
            radial_boxes[order] += sum(verify_ray(order, lower, upper, x) for x in cosines)
    return radial_boxes


if __name__ == "__main__":
    print("Retained-row coefficient comparisons:", verify_row_loss())
    print("Certified radial intervals:", verify_all_bands())
