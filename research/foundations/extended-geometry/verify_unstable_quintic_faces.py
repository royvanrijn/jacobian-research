#!/usr/bin/env python3
"""Replay the lattice census and cyclic-cover signatures using exact integers."""

from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import gcd


EXPECTED_FACES = {
    (5, 5, 3): (-5, -2, 1),
    (2, 2, 1): (-3, -2, -1, 0, 1),
    (3, 5, 2): (-5, -3, -1, 1),
    (3, 6, 2): (-4, -2, 0, 2),
    (2, 3, 1): (-2, -1, 0, 1),
    (3, 3, 1): (-1, 0, 1),
    (4, 10, 3): (-5, -2, 1),
    (3, 7, 2): (-3, -1, 1),
    (4, 11, 3): (-4, -1, 2),
    (5, 15, 4): (-5, -1, 3),
    (4, 12, 3): (-3, 0, 3),
    (3, 8, 2): (-2, 0, 2),
    (2, 4, 1): (-1, 0, 1, 2),
    (3, 9, 2): (-1, 1, 3),
}


def enumerate_faces():
    nodes = [
        (angle, height)
        for angle in range(-5, 5)
        for height in range(max(0, -2 * angle), (5 if angle < 0 else 4) - angle + 1)
    ]
    negative = [(angle, height) for angle, height in nodes if angle < 0]
    positive = [(angle, height) for angle, height in nodes if angle > 0]
    top = [(angle, height) for angle, height in negative if angle + height == 5]
    faces = defaultdict(set)
    for (left, high), (right, low) in product(negative, positive):
        rho = Fraction(high - low, right - left)
        intercept = high + rho * left
        if rho <= 1 or intercept >= 5:
            continue
        if all(height + rho * angle > intercept for angle, height in top):
            continue
        faces[rho, intercept].update(
            angle for angle, height in nodes if height + rho * angle == intercept
        )
    return {key: tuple(sorted(angles)) for key, angles in faces.items() if len(angles) >= 3}


def signature(order, exponents):
    assert sum(exponents) == 2 * order
    assert gcd(order, *exponents) == 1
    return {
        index: sum(index * exponent % order for exponent in exponents) // order - 1
        for index in range(1, order)
    }


def genus(order, exponents):
    return 1 - order + sum(order - gcd(order, exponent) for exponent in exponents) // 2


def main():
    faces = enumerate_faces()
    assert faces == {
        (Fraction(numerator, denominator), Fraction(intercept, denominator)): angles
        for (numerator, intercept, denominator), angles in EXPECTED_FACES.items()
    }

    for (rho, intercept), angles in faces.items():
        if intercept.denominator == 1:
            continue
        step = rho.denominator
        order = intercept.numerator
        left, right = angles[0], angles[-1]
        assert step == intercept.denominator
        assert all((angle - left) % step == 0 for angle in angles)
        exponents = (order + left,) + (step,) * ((right - left) // step) + (order - right,)
        if order in (7, 9, 10, 11, 15):
            dimensions = signature(order, exponents)
            primitive = {i: d for i, d in dimensions.items() if gcd(i, order) == 1}
            assert all(d + primitive[order - i] == 2 for i, d in primitive.items())
            assert sum(d == 1 for d in primitive.values()) <= 4

    for order, exponents, expected in (
        (7, (4, 2, 2, 6), (1, 1, 2, 0, 1, 1)),
        (9, (8, 2, 2, 6), (1, 1, 2, 0, 1, 1)),
        (10, (5, 3, 3, 9), (1, 2, 0, 1)),
        (11, (7, 3, 3, 9), (1, 1, 2, 0, 0, 2, 2, 0, 1, 1)),
        (15, (10, 4, 4, 12), (1, 1, 0, 2, 0, 2, 1, 1)),
    ):
        dims = signature(order, exponents)
        assert tuple(dims[i] for i in dims if gcd(i, order) == 1) == expected
        assert sum(dims.values()) == genus(order, exponents)

        for divisor in range(2, order):
            if order % divisor == 0:
                assert sum(d for i, d in dims.items() if order // gcd(i, order) == divisor) <= 2

        outer_left, root, _, outer_right = exponents
        repeated = (outer_left, 2 * root, outer_right)
        dims = signature(order, repeated)
        assert all(d + dims[order - i] <= 1 for i, d in dims.items() if gcd(i, order) == 1)

    for order, exponents, expected_genus in (
        (3, (2, 1, 1, 2), 2),
        (3, (2, 2, 2), 1),
        (3, (1, 1, 1, 1, 2), 3),
        (3, (1, 2, 1, 2), 2),
        (4, (3, 1, 1, 3), 3),
        (4, (3, 2, 3), 1),
        (4, (2, 1, 1, 1, 3), 4),
        (4, (2, 2, 1, 3), 2),
        (4, (2, 3, 3), 1),
        (5, (3, 3, 4), 2),
        (5, (4, 2, 2, 2), 4),
        (5, (4, 4, 2), 2),
    ):
        assert genus(order, exponents) == expected_genus == sum(signature(order, exponents).values())

    print("Verified 14 unstable quintic Newton faces and the cyclotomic signatures.")


if __name__ == "__main__":
    main()
