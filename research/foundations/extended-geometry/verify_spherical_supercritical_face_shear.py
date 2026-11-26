#!/usr/bin/env python3
"""Replay supercritical weighted-face moment examples in the literal chart."""

from fractions import Fraction
from math import comb

from verify_spherical_critical_face_p_power_moments import multiply, power


def moment(phase, order):
    return sum(
        Fraction(coefficient, height + 1)
        for (angle, height), coefficient in power(phase, order).items()
        if angle == 0 and height % 2 == 0
    )


def add(left, right):
    result = left.copy()
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, 0) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def regular_and_in_strip(phase, minimum, maximum, d=2):
    for angle, height in phase:
        assert height >= 0
        assert minimum <= height + Fraction(angle, d) <= maximum
    for angle in {angle for angle, _ in phase if angle < 0}:
        coefficients = {
            height: coefficient for (term_angle, height), coefficient in phase.items()
            if term_angle == angle
        }
        # Check the requisite order at both roots of 1-t^2.
        for _ in range(-angle):
            assert all(sum(value * root**height for height, value in coefficients.items())
                       == 0 for root in (-1, 1))
            coefficients = {height - 1: height * value
                            for height, value in coefficients.items() if height > 0}


def cluster_phase(clusters, d=2, separation=1):
    phase = {(-2, 0): 1}
    for slope, negative_count, positive_count in clusters:
        w = {(d, 0): 1, (0, 1): -1, (separation, 0): -slope}
        phase = multiply(phase, power(add(w, {(0, 0): -1}), negative_count))
        phase = multiply(phase, power(add(w, {(0, 0): 1}), positive_count))
    regular_and_in_strip(phase, Fraction(-2, d), 4 - Fraction(2, d), d)
    return phase


def verify():
    w = {(2, 0): 1, (0, 1): -1}
    quadratic = add(power(w, 2), {(0, 0): -1})
    f = multiply(multiply(power(w, 30), power(quadratic, 4)), {(-2, 0): 1, (-4, 1): -2})
    perturbation = power({exponent: -value for exponent, value in w.items()}, 33)
    g = add(f, perturbation)
    for phase in (f, g):
        regular_and_in_strip(phase, 29, 37)
    assert moment(f, 2) == 0
    assert moment(g, 2) == Fraction(1, 67)
    assert (75 * moment(g, 2)).numerator % 5 == 0

    w[(1, 0)] = -1
    quadratic = add(power(w, 2), {(0, 0): -1})
    f = add(multiply(multiply(power(w, 16), power(quadratic, 2)), {(-2, 0): 1}),
            power(w, 17))
    regular_and_in_strip(f, 7, 19)
    assert moment(f, 2) == Fraction(1, 35)
    assert (39 * moment(f, 2)).numerator % 13 == 0
    print("Order-two moments verified: 0, 1/67, and 1/35")

    pure = cluster_phase([(0, 1, 2), (1, 1, 0)])
    critical = cluster_phase([(0, 1, 1), (1, 1, 1)])
    separated = cluster_phase([(0, 1, 0), (1, 1, 0), (2, 0, 1), (3, 0, 1)])
    assert all(moment(pure, order) == 0 for order in range(1, 4))
    assert moment(critical, 1) == Fraction(-2, 3)
    assert moment(critical, 2) == Fraction(8, 15)
    assert moment(separated, 1) == Fraction(14, 3)
    print("Three cluster regimes checked in the literal Laurent chart")

    phase = cluster_phase([(0, 1, 2), (1, 1, 0)], d=3, separation=2)
    assert moment(phase, 1) == Fraction(2, 3)
    assert moment(phase, 2) == Fraction(64, 105)
    p, n, large_cluster, negative_count = 53, 27, 21, 13
    order = (p**2 - 1) // n
    assert p**2 - 1 == n * order
    denominator = (large_cluster * order + 1) * comb(large_cluster * order, negative_count * order)
    assert denominator % p == 0 and denominator % p**2 != 0
    print("Quadratic separation checked, including exact depth one at p=53")


if __name__ == "__main__":
    verify()
