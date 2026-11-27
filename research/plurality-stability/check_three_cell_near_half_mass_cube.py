"""Directed scalar bound for the planar three-cell mass cube near rho=1/2.

Run with python-flint. The mathematical reduction and tail estimate are in
THREE_CELL_NEAR_HALF_MASS_CUBE.md.
"""

import argparse
from fractions import Fraction
import json
from math import factorial, inf, nextafter

from flint import arb, ctx


def upper_float(value):
    return nextafter(float(value.upper()), inf)


def ball(value):
    return arb(value.numerator) / value.denominator


def angular_mixture():
    displacement = arb(11) / 25
    slope = arb(76) / 25
    tangent_angle = 5 * arb.pi() / 12
    tangent_slope = arb.pi() * tangent_angle.sin()
    tangent_at_sixth = tangent_angle.cos() + arb(1) / 2 - tangent_slope / 12
    assert slope > 3 and slope > tangent_slope
    assert tangent_at_sixth < slope / 6

    first = -arb(1) / 2
    last = arb(419) / 500
    lower = arb(17) / 20
    assert first + slope * displacement < last
    angle = arb.pi() * displacement
    assert 2 * (1 - angle.cos()) / angle**2 > lower
    upper = [(first, arb(5) / 12), ((first + last) / 2, arb(1) / 2), (last, arb(1) / 12)]
    assert all(abs(node) < lower for node, _ in upper)
    return upper, lower


def trace_bound(rho, panels, upper, lower):
    d = 1 - rho * rho
    kappa = rho * rho / d
    rate = 2 * rho / d
    beta = 1 + kappa
    center = arb(14) / 25
    radius = 6
    order = 12

    gamma = [arb(1), arb.pi().sqrt() / 2]
    for n in range(2, order + 1):
        gamma.append(n * gamma[n - 2] / 2)
    coefficients = []
    for n in range(1, order + 1):
        power = beta ** (1 + n // 2)
        if n % 2:
            power *= beta.sqrt()
        moment = lower**n - sum((weight * node**n for node, weight in upper), arb(0))
        coefficients.append(
            rate**n * moment * gamma[n] / (power * factorial(n))
        )
    assert all(coefficient > 0 for coefficient in coefficients)

    exponents = [(rate * node / 2 - kappa, rate * node / 2, weight)
                 for node, weight in [*upper, (lower, arb(1))]]
    center_factors = [(s * center**2).exp() for s, _, _ in exponents]
    center_diagonal = sum((weight * value**2
                           for (_, _, weight), value in zip(exponents, center_factors)), arb(0)) / (2 * d)
    assert center_diagonal < 2

    width = arb(radius) / panels
    core = arb(0)
    for index in range(panels):
        x = (index * width).union((index + 1) * width)
        polynomial = sum((coefficient * x**i for i, coefficient in enumerate(coefficients)), arb(0))
        assert polynomial > 0
        diagonal = arb(0)
        for (s, v, weight), center_factor in zip(exponents, center_factors):
            factor = (s * x**2).exp()
            diagonal += weight * ((factor - center_factor)**2
                                  + 2 * factor * center_factor * (1 - (-v * (x - center)**2).exp()))
        value = 4 * (-(1 - kappa) * x**2).exp() * diagonal / polynomial
        core += width * value.upper()

    tail_rate = 1 + kappa - rate * lower
    assert tail_rate > 0 and rate * lower - 2 * kappa > 0
    tail = 12 * (-tail_rate * radius**2).exp() / (tail_rate * coefficients[-1] * radius**order)
    return core, tail


def certificate(panels=4096, point_only=False):
    ctx.prec = 100
    upper, lower = angular_mixture()

    point_core, point_tail = trace_bound(arb(1) / 2, panels, upper, lower)
    assert point_core + point_tail < arb(44) / 25
    result = {
        "status": "pass",
        "mass_radius": "11/50",
        "panels_per_correlation_ball": panels,
        "precision_bits": ctx.prec,
        "point": {
            "rho": "1/2",
            "core_upper": upper_float(point_core),
            "tail_upper": upper_float(point_tail),
            "total_upper": upper_float(point_core + point_tail),
        },
    }
    if point_only:
        return result

    maximum_upper = 0.0
    for index in range(440, 530):
        rho = ball(Fraction(index, 1000)).union(ball(Fraction(index + 1, 1000)))
        core, tail = trace_bound(rho, panels, upper, lower)
        assert core + tail < arb(48) / 25, index
        maximum_upper = max(maximum_upper, upper_float(core + tail))
    result["interval"] = {
        "rho_bounds": ["11/25", "53/100"],
        "slabs": 90,
        "slab_width": "1/1000",
        "total_upper": maximum_upper,
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps(certificate(point_only=args.point_only), indent=2))
