#!/usr/bin/env python3
"""Exact joint modular shape exclusion at orders five and six."""

from collections import Counter
from fractions import Fraction
from math import comb, factorial, gcd, lcm, prod


def primitive_parent(order, rows):
    coefficients = [
        comb(order, k) * prod(
            factorial((order - k) * a + k * b) for a, b in rows
        )
        for k in range(order + 1)
    ]
    content = gcd(*coefficients)
    return [value // content for value in coefficients]


def valuation(value, prime):
    depth = 0
    while value % prime == 0:
        depth += 1
        value //= prime
    return depth


def newton_slopes(coefficients, prime):
    hull = []
    for index, value in enumerate(coefficients):
        depth = valuation(value, prime)
        while len(hull) > 1:
            left, left_depth = hull[-2]
            middle, middle_depth = hull[-1]
            if (middle_depth - left_depth) * (index - middle) < (
                (depth - middle_depth) * (middle - left)
            ):
                break
            hull.pop()
        hull.append((index, depth))
    slopes = Counter()
    for (left, low), (right, high) in zip(hull, hull[1:]):
        slopes[Fraction(high - low, right - left)] += right - left
    return slopes


def trim(polynomial):
    while polynomial and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def residue_gcd(first, second, prime):
    first = trim([value % prime for value in first])
    second = trim([value % prime for value in second])
    while second:
        inverse = pow(second[-1], -1, prime)
        while len(first) >= len(second):
            shift = len(first) - len(second)
            multiplier = first[-1] * inverse % prime
            for index, value in enumerate(second):
                first[shift + index] = (first[shift + index] - multiplier * value) % prime
            trim(first)
        first, second = second, first
    inverse = pow(first[-1], -1, prime)
    return [value * inverse % prime for value in first]


def main():
    rows = [(0, 1)] * 3 + [(1, 0)] * 3
    parents = [primitive_parent(order, rows) for order in (5, 6)]
    assert parents == [[100, 4, 1, 1, 4, 100], [3600, 100, 16, 9, 16, 100, 3600]]
    endpoint_product = prod(parent[0] * parent[-1] for parent in parents)
    assert endpoint_product == 2**12 * 3**4 * 5**8
    shared_slopes = {}
    for prime in (2, 3, 5):
        first, second = (newton_slopes(parent, prime) for parent in parents)
        shared_slopes[prime] = first & second
        assert sum(shared_slopes[prime].values()) >= 2
    assert shared_slopes == {2: Counter({-1: 2, 1: 2}), 3: Counter({0: 4}), 5: Counter({0: 2})}
    divisor = prod(prime**int(max(slopes) - min(slopes)) for prime, slopes in shared_slopes.items())
    assert divisor == 4
    joint_modulus, joint_residue = 1, 0
    separate = {}
    for prime, expected in ((7, [1, 5, 1]), (31, [1, 3, 1])):
        assert all(parent[0] * parent[-1] % prime for parent in parents)
        quadratic = residue_gcd(*parents, prime)
        assert quadratic == expected
        constant, linear, _ = quadratic
        local_shape = linear**2 * pow(constant, -1, prime) % prime
        separate[prime] = [
            lift for lift in range(4 * divisor) if (lift - divisor * local_shape) % prime == 0
        ]
        joint_residue += joint_modulus * ((local_shape - joint_residue)
                                          * pow(joint_modulus, -1, prime) % prime)
        joint_modulus = lcm(joint_modulus, prime)
    assert separate == {7: [2, 9], 31: [5]}
    assert (joint_modulus, joint_residue) == (217, 102)
    integer_residue = divisor * joint_residue % joint_modulus
    assert integer_residue == 191
    assert not any((lift - integer_residue) % joint_modulus == 0 for lift in range(4 * divisor))
    print("Verified shared widths, two separate compatible traces, and joint trace exclusion.")


if __name__ == "__main__":
    main()
