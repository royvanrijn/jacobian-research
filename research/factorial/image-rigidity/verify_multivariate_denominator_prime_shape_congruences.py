#!/usr/bin/env python3
"""Exact shape congruences from two unequal p-adic root valuations."""

from fractions import Fraction

from verify_multivariate_joint_native_shape_reconstruction import (
    newton_slopes,
    primitive_parent,
    valuation,
)


def unit_roots(coefficients, prime, root_valuation, precision=1):
    modulus = prime**precision
    depths = [valuation(value, prime) + index * root_valuation
              for index, value in enumerate(coefficients)]
    bottom = min(depths)
    scaled = [
        value // prime**valuation(value, prime) * prime**(depth - bottom) % modulus
        for value, depth in zip(coefficients, depths, strict=True)
    ]
    return {
        root for root in range(1, modulus) if root % prime
        and sum(value * pow(root, index, modulus) for index, value in enumerate(scaled))
        % modulus == 0
    }


def split_congruence(prime, budget_depth, spread, budget_unit, units, precision):
    modulus = prime**precision
    u, v = units
    first = u * pow(v, -1, modulus) % modulus
    second = v * pow(u, -1, modulus) % modulus
    scaled_shape = (first + 2 * prime**spread + prime**(2 * spread) * second) % modulus
    factor = prime**(budget_depth - spread)
    return (factor * budget_unit * scaled_shape % (factor * modulus), factor * modulus)


def ordinary_profile():
    rows = [(0, 1)] * 3 + [(1, 0)] * 3
    parents = [primitive_parent(order, rows) for order in (5, 6)]
    shared = newton_slopes(parents[0], 2) & newton_slopes(parents[1], 2)
    assert shared == {-1: 2, 1: 2}
    units = [
        set.intersection(*(unit_roots(parent, 2, slope) for parent in parents))
        for slope in (-1, 1)
    ]
    assert units == [{1}, {1}]
    assert split_congruence(2, 2, 2, 1, [1, 1], 1) == (1, 2)
    candidates = [lift for lift in range(16) if lift % 7 == 2]
    assert candidates == [2, 9]
    assert [lift for lift in candidates if lift % 4 == 0 or lift % 2 == 1] == [9]


def quadratic_controls():
    counts = 0
    for prime, s, t, linear_unit, constant_unit, extra in [
        (2, -2, 0, 1, 3, 0),
        (3, 0, 1, 2, 5, 3),
        (5, -1, 2, 3, 2, 2),
    ]:
        spread = t - s
        depth = spread + extra
        linear = Fraction(prime)**s * linear_unit
        constant = Fraction(prime)**(s + t) * constant_unit
        assert linear**2 < 4 * constant
        exact_shape = linear**2 / constant
        exact_budget = prime**depth * constant_unit
        integer = exact_budget * exact_shape
        assert integer.denominator == 1
        for precision in (1, 2, 3):
            modulus = prime**precision
            units = [
                root for root in range(1, modulus) if root % prime
                and (root**2 + linear_unit * root + prime**spread * constant_unit)
                % modulus == 0
            ]
            other_units = [
                root for root in range(1, modulus) if root % prime
                and (prime**spread * root**2 + linear_unit * root + constant_unit)
                % modulus == 0
            ]
            assert len(units) == len(other_units) == 1
            residue, joint_modulus = split_congruence(
                prime, depth, spread, constant_unit, (units[0], other_units[0]), precision
            )
            assert integer % joint_modulus == residue
            assert valuation(integer.numerator, prime) == extra
            counts += 1
    return counts


if __name__ == "__main__":
    ordinary_profile()
    print(f"Verified ordinary branch residues and {quadratic_controls()} valued shape congruences.")
