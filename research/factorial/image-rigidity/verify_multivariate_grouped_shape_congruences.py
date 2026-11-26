#!/usr/bin/env python3
"""Exact coefficientwise projective congruences for a two-slope quartic."""

from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, lcm

from verify_multivariate_joint_native_shape_reconstruction import (
    newton_slopes,
    primitive_parent,
    valuation,
)


def rational_valuation(value, prime):
    return valuation(value.numerator, prime) - valuation(value.denominator, prime)


def power_sums(monic, maximum):
    degree = len(monic) - 1
    out = [Fraction(degree)]
    for i in range(1, maximum + 1):
        value = sum(monic[degree - j] * out[i - j]
                    for j in range(1, min(i - 1, degree) + 1))
        if i <= degree:
            value += i * monic[degree - i]
        out.append(-value)
    return out


def elementary_of_powers(monic, power):
    degree = len(monic) - 1
    sums = power_sums(monic, power * degree)
    elementary = [Fraction(1)]
    for k in range(1, degree + 1):
        elementary.append(sum(
            (-1)**(i - 1) * elementary[k - i] * sums[i * power]
            for i in range(1, k + 1)
        ) / k)
    return elementary


def remainder(polynomial, monic, modulus):
    polynomial = [value % modulus for value in polynomial]
    for index in range(len(polynomial) - 1, 1, -1):
        factor = polynomial[index]
        for j, value in enumerate(monic):
            polynomial[index - 2 + j] = (
                polynomial[index - 2 + j] - factor * value
            ) % modulus
    return polynomial[:2]


def lift_positive_group(normalized, depth):
    factor = [1, 0, 1]
    modulus = 3
    assert remainder(normalized, factor, modulus) == [0, 0]
    for _ in range(1, depth):
        next_modulus = 3 * modulus
        candidates = [[factor[0] + modulus * constant, factor[1] + modulus * linear, 1]
                      for constant, linear in product(range(3), repeat=2)]
        factors = [candidate for candidate in candidates
                   if remainder(normalized, candidate, next_modulus) == [0, 0]]
        assert len(factors) == 1
        factor = factors[0]
        modulus = next_modulus
    return factor


def fourth_power_elementary(factor, modulus):
    constant, linear, lead = factor
    assert lead == 1
    return [1, (linear**4 - 4 * linear**2 * constant + 2 * constant**2) % modulus,
            constant**4 % modulus]


def ordinary_quartic():
    rows = [(0, 1)] * 3 + [(1, 0)] * 3
    parent = primitive_parent(4, rows)
    assert parent == [144, 9, 4, 9, 144]
    assert newton_slopes(parent, 3) == {-1: 2, 1: 2}
    monic = [Fraction(value, parent[-1]) for value in parent]
    assert monic[0] == 1
    shape = elementary_of_powers(monic, 4)
    depths = [0, 4, 8, 4, 0]
    assert [-rational_valuation(value, 3) for value in shape] == depths
    normalized = [value * 3**index // 9 for index, value in enumerate(parent)]
    assert normalized == [16, 3, 4, 27, 1296]
    calculations = 0
    for precision in range(1, 9):
        modulus = 3**precision
        upper = lift_positive_group(normalized, precision)
        inverse = pow(upper[0], -1, modulus)
        lower = [inverse, upper[1] * inverse % modulus, 1]
        unit_elementaries = [fourth_power_elementary(factor, modulus)
                             for factor in (lower, upper)]
        unit_norm = lower[0] * upper[0] % modulus
        assert unit_norm == 1
        for k, depth in enumerate(depths):
            calculated = sum(
                3**(depth - 4 * i + 4 * j)
                * unit_elementaries[0][i] * unit_elementaries[1][j]
                for i, j in product(range(3), repeat=2) if i + j == k
            ) % modulus
            actual = 3**depth * shape[k]
            actual_residue = actual.numerator * pow(actual.denominator, -1, modulus) % modulus
            assert calculated == actual_residue
            calculations += 1
    return calculations


def fractional_control():
    monic = list(map(Fraction, [3, 0, 3, 0, 1, 0, 1]))
    shape = [value / 3**k for k, value in enumerate(elementary_of_powers(monic, 6))]
    assert shape == [
        1, Fraction(-2, 3), Fraction(55, 9), -4, Fraction(29, 3), -6, 1
    ]
    valuation_list = [Fraction(-1)] * 2 + [Fraction(1, 2)] * 4
    minimum = sum(valuation_list[:3])
    assert minimum == Fraction(-3, 2)
    assert (-minimum).__floor__() == 1
    for k, value in enumerate(shape):
        minimum = sum(valuation_list[:k])
        assert rational_valuation(value, 3) >= minimum.__ceil__()




def frobenius_gain_control():
    def middle_shape(factors, modulus):
        elementary = [fourth_power_elementary(factor, modulus) for factor in factors]
        inverse_norm = pow(factors[0][0] * factors[1][0], -2, modulus)
        return inverse_norm * sum(
            2**(8 - 4 * i + 4 * j) * elementary[0][i] * elementary[1][j]
            for i, j in product(range(3), repeat=2) if i + j == 2
        ) % modulus

    base = [[1, 0, 1], [3, 1, 1]]
    for precision in range(1, 5):
        deeper = 2**(precision + 1)
        expected = middle_shape(base, deeper)
        for digits in product(range(2), repeat=4):
            moved = [[base[group][j] + 2**precision * digits[2 * group + j]
                      for j in range(2)] + [1] for group in range(2)]
            assert middle_shape(moved, deeper) == expected

    parent = [4, 4, 4, 1, 1]
    assert newton_slopes(parent, 2) == {Fraction(-2, 3): 3, 0: 1}
    shape = elementary_of_powers(list(map(Fraction, parent)), 4)
    middle = shape[2] / 4**2
    first_possible = Fraction(-4, 3).__ceil__()
    assert first_possible == -1 and first_possible % 2
    assert rational_valuation(middle, 2) >= first_possible + 1


def ambiguity_control():
    denominators = (6, 10)
    common_scale = lcm(*denominators)
    scales = [common_scale // item for item in denominators]
    depths = {2: (2, 2), 3: (2, 1), 5: (1, 2)}
    original = [(0, 0), (1, 1)]
    scaled = [tuple(scale * value for scale, value in zip(scales, pair, strict=True))
              for pair in original]
    levels = {}
    for prime, coordinate_depths in depths.items():
        depth = min(value + valuation(scale, prime)
                    for value, scale in zip(coordinate_depths, scales, strict=True))
        assert depth == 2
        for level in range(1, depth + 1):
            modulus = prime**level
            levels[prime, level] = {tuple(value % modulus for value in pair)
                                    for pair in scaled}
            assert len(levels[prime, level]) == 2
    possible = [
        vector for vector in product(range(7), repeat=2)
        if all(tuple(value % prime**level for value in vector) in allowed
               for (prime, level), allowed in levels.items())
    ]
    assert possible == scaled
    differences = [gcd(*(abs(a - b) for a, b in zip(left, right, strict=True)))
                   for left, right in combinations(possible, 2)]
    assert differences == [1]
    collision_lower_bound_for_three = 1
    for (prime, _), allowed in levels.items():
        classes = len(allowed)
        quotient, extra = divmod(3, classes)
        collision_lower_bound_for_three *= prime**(classes * comb(quotient, 2) + extra * quotient)
    assert collision_lower_bound_for_three == 900 > 6**comb(3, 2)


if __name__ == "__main__":
    checks = ordinary_quartic()
    fractional_control()
    ambiguity_control()
    frobenius_gain_control()
    print(f"Verified {checks} lifted quartic coefficient congruences, fractional rounding, Frobenius precision, and the correlated larger-sieve bound.")
