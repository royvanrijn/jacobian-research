#!/usr/bin/env python3
"""Enclose the square-root divisor mark and a longer-modulus correction."""

from fractions import Fraction

from verify_eisenstein_factor_moments import (
    factor_constant,
    linear_interval,
    show_interval,
)

ETA = Fraction(1, 5)
HALF = Fraction(1, 2)
THREE_FIFTHS = Fraction(3, 5)
SCALE = 10**12


def log_bounds(ratio):
    x = (ratio - 1) / (ratio + 1)
    terms = 18
    partial = 2 * sum(x**(2*k + 1) / (2*k + 1) for k in range(terms))
    tail = 2 * x**(2*terms + 1) / ((2*terms + 1) * (1 - x*x))
    return partial, partial + tail


def scaled_round_down(value):
    return value.numerator * SCALE // value.denominator


def scaled_round_up(value):
    return -scaled_round_down(-value)


def integrate(left, right, log_argument, grid=512):
    """Enclose the integral for a monotone logarithm argument."""
    if left == right:
        return Fraction(0), Fraction(0)
    assert ETA <= left <= right <= THREE_FIFTHS
    step = (right - left) / grid
    lower = upper = 0
    for index in range(grid):
        a = left + index * step
        b = a + step
        arguments = sorted((log_argument(a), log_argument(b)))
        logs = log_bounds(arguments[0])[0], log_bounds(arguments[1])[1]
        endpoints = a * (1 - a), b * (1 - b)
        divisor_low = min(endpoints)
        divisor_high = HALF**2 if a <= HALF <= b else max(endpoints)
        lower += scaled_round_down(logs[0] / divisor_high)
        upper += scaled_round_up(logs[1] / divisor_low)
    return step * lower / SCALE, step * upper / SCALE


def marginal_integral(left, right):
    return integrate(left, right, lambda s: (1 - s - ETA) / ETA)


def certificate(theta, constants, large, bonus=False):
    c1, c2, c3, c4 = constants
    left = 1 - theta
    assert Fraction(3, 10) <= left <= HALF
    selected = marginal_integral(left, HALF)
    if theta >= THREE_FIFTHS:
        bounds = [c1, c2, c3, c4, large, selected]
        weights = [1, -HALF, -Fraction(5, 4), -Fraction(1, 8), HALF, HALF]
        if bonus:
            bounds.append(integrate(left, Fraction(2, 5), lambda s: left / (1-s-left)))
            weights.append(HALF)
        return linear_interval(bounds, weights)
    partial_large = marginal_integral(HALF, theta)
    return linear_interval(
        [c1, c2, c3, c4, large, partial_large, selected],
        [1, -Fraction(1, 2), -Fraction(5, 4), -Fraction(1, 8), 1, -HALF, HALF],
    )


def divide_interval(numerator, denominator):
    quotients = [n / d for n in numerator for d in denominator]
    return min(quotients), max(quotients)


def main():
    constants = [factor_constant(ETA, k) for k in range(1, 5)]
    total = linear_interval(constants, [1, HALF, Fraction(1, 4), Fraction(1, 8)])
    large = marginal_integral(HALF, THREE_FIFTHS)
    print("triple constant", show_interval(constants[2]),
          "triple constant with a large prime", show_interval(large))
    endpoints = ((HALF, -1), (Fraction(629, 1000), -1), (Fraction(631, 1000), 1),
                 (Fraction(16, 25), 1), (Fraction(33, 50), 1))
    for theta, sign in endpoints:
        numerator = certificate(theta, constants, large)
        assert numerator[0] > 0 if sign > 0 else numerator[1] < 0
        margin = divide_interval(numerator, total)
        if theta == HALF:
            assert Fraction(-65, 1000) < margin[0] <= margin[1] < Fraction(-63, 1000)
        if theta == Fraction(16, 25):
            assert Fraction(9, 1000) < margin[0] <= margin[1] < Fraction(1, 100)
        print("pair-product exponent", theta, "margin/total", show_interval(margin))

    improved = ((Fraction(617, 1000), -1), (Fraction(31, 50), 1), (Fraction(5, 8), 1),
                (Fraction(16, 25), 1))
    for theta, sign in improved:
        numerator = certificate(theta, constants, large, bonus=True)
        assert numerator[0] > 0 if sign > 0 else numerator[1] < 0
        margin = divide_interval(numerator, total)
        if theta == Fraction(31, 50):
            assert Fraction(25, 10000) < margin[0] <= margin[1] < Fraction(34, 10000)
        print("improved exponent", theta, "margin/total", show_interval(margin))


if __name__ == "__main__":
    main()
