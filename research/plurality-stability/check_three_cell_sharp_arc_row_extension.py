"""Directed extension of the sharp scalar arc bound above 11/20.

The radial moments use piecewise linear interpolation with a rigorous
second derivative remainder. The finite Gram norm uses weighted row
sums. Floating point arithmetic chooses weights, and Arb proves the
bound for those exact dyadic weights on the entire correlation slab.
"""

import json
from fractions import Fraction
from math import factorial

from flint import arb, ctx

from check_three_cell_global_scalar_through_half import radial_coefficients, upper_float
from check_three_cell_sharp_positive_arc_envelope import sharp_moments


MODE_ORDER = 112
ROW_ORDER = 24
RADIAL_LIMIT = 12
RADIAL_PANELS = 2048


def enclose_rational(value):
    return arb(value.numerator) / value.denominator


def radial_moments(d, kappa, beta, row):
    exponents = (beta, arb(1), 1 - kappa)
    orders = (2 * MODE_ORDER, MODE_ORDER, 0)
    output = [[arb(0) for _ in range(order + 1)] for order in orders]
    width = arb(RADIAL_LIMIT) / RADIAL_PANELS

    def endpoint(x):
        polynomial, first, second = row[-1], arb(0), arb(0)
        for coefficient in reversed(row[:-1]):
            second = second * x + 2 * first
            first = first * x + polynomial
            polynomial = polynomial * x + coefficient
        assert polynomial > 0
        factors = [8 * d * (-exponent * x**2).exp() / polynomial for exponent in exponents]
        return polynomial, first, second, factors

    left = arb(0)
    left_polynomial, _, _, left_factors = endpoint(left)
    for panel in range(RADIAL_PANELS):
        right = (panel + 1) * width
        right_polynomial, first, second, right_factors = endpoint(right)
        log_slope_bound = first / left_polynomial
        curvature_bound = second / left_polynomial
        errors = []
        for exponent, left_factor in zip(exponents, left_factors):
            assert exponent > 0
            derivative_bound = left_factor * (
                (2 * exponent * right + log_slope_bound)**2
                + 2 * exponent + curvature_bound + log_slope_bound**2
            )
            error = (derivative_bound * width**2 / 8).upper()
            errors.append((-error).union(error))

        left_power, right_power = left, right
        mass = right - left
        for n in range(2 * MODE_ORDER + 1):
            left_power *= left
            right_power *= right
            first_moment = (right_power - left_power) / (n + 2)
            left_weight = (right * mass - first_moment) / width
            right_weight = (first_moment - left * mass) / width
            for moments, at_left, at_right, error in zip(output, left_factors, right_factors, errors):
                if n < len(moments):
                    moments[n] += at_left * left_weight + at_right * right_weight + error * mass
            mass = first_moment
        left, left_polynomial, left_factors = right, right_polynomial, right_factors

    radius = arb(RADIAL_LIMIT)
    for moments, exponent in zip(output, exponents):
        for n in range(len(moments)):
            decay = exponent - arb(max(0, n - 1)) / (2 * radius**2)
            assert decay > 0
            tail = 4 * d * radius**(n - 1) * (-exponent * radius**2).exp() / (row[0] * decay)
            moments[n] += arb(0).union(tail)
    tt, uu, (vv,) = output
    assert vv > 0
    return tt, uu, vv


def gram_norm(modes, tt, uu, vv):
    size = len(modes)
    roots = [mode.sqrt() for mode in modes]
    entries = [[arb(0) for _ in modes] for _ in modes]
    for i in range(size):
        for j in range(i + 1):
            covariance = tt[i + j] - uu[i] * uu[j] / vv
            entries[i][j] = entries[j][i] = abs(covariance) * roots[i] * roots[j]

    approximate = [[upper_float(entry) for entry in line] for line in entries]
    trial = [1.0] * size
    for _ in range(24):
        image = [sum(a * w for a, w in zip(line, trial)) for line in approximate]
        scale = max(image)
        trial = [a / scale for a in image]

    weights = [enclose_rational(Fraction(max(value, 2**-55))) for value in trial]
    bound = arb(0)
    for i, line in enumerate(entries):
        row_sum = sum((entry * weight for entry, weight in zip(line, weights)), arb(0))
        bound = bound.union((row_sum / weights[i]).upper()).upper()
    return bound


def slab_bound(rho, angular):
    d = 1 - rho**2
    kappa = rho**2 / d
    beta = 1 + kappa
    rate = 2 * rho / d
    row, modes = radial_coefficients(d, beta, rate, *angular, ROW_ORDER)
    finite = gram_norm(modes, *radial_moments(d, kappa, beta, row))

    first = MODE_ORDER + 1
    small_ratio = 4 * rate / (first + 1)
    eta = 1 + 19 * kappa / 10
    large_ratio = rate / eta
    assert 0 < small_ratio < 1 and 0 < large_ratio < 1
    small_tail = 16 * (4 * rate)**first / (row[0] * factorial(first) * (1 - small_ratio))
    large_tail = 3200 * large_ratio**first / (31 * d * eta * (1 - large_ratio))
    return finite + small_tail + large_tail


def certificate():
    ctx.prec = 100
    assert (arb.pi() / 10).cos()**2 > arb(9) / 10
    assert (2 * arb.pi() / 9).cos()**2 < arb(3) / 5
    angular = sharp_moments(MODE_ORDER)
    target = arb(1999) / 1000
    upper_bounds = []
    widths = set()

    def certify(left, right, depth):
        enclosed = enclose_rational(left).union(enclose_rational(right))
        bound = slab_bound(enclosed, angular)
        if bound < target:
            upper_bounds.append(upper_float(bound))
            widths.add(right - left)
            return
        assert depth < 8, (str(left), str(right), upper_float(bound))
        midpoint = (left + right) / 2
        certify(left, midpoint, depth + 1)
        certify(midpoint, right, depth + 1)

    for start in range(550, 580):
        certify(Fraction(start, 1000), Fraction(start + 1, 1000), 0)

    return {
        "status": "pass",
        "precision_bits": ctx.prec,
        "radial_panels_per_correlation_ball": RADIAL_PANELS,
        "radial_limit": RADIAL_LIMIT,
        "mode_order": MODE_ORDER,
        "row_order": ROW_ORDER,
        "interval": {
            "rho_bounds": ["11/20", "29/50"],
            "correlation_slabs": len(upper_bounds),
            "slab_widths": [str(width) for width in sorted(widths)],
            "total_upper": max(upper_bounds),
        },
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
