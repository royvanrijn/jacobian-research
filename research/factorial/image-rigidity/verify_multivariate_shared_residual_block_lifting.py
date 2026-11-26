#!/usr/bin/env python3
"""Replay exact shared-block lifts in an ordinary paired factorial profile."""

from collections import Counter
from itertools import product

from verify_multivariate_height_capped_native_unit_interpolation import coefficient_data
from verify_multivariate_joint_native_shape_reconstruction import (
    newton_slopes,
    primitive_parent,
    residue_gcd,
    trim,
)


def divide_mod(polynomial, monic, modulus):
    polynomial = trim([value % modulus for value in polynomial])
    quotient = [0] * max(0, len(polynomial) - len(monic) + 1)
    for index in range(len(quotient) - 1, -1, -1):
        leading = polynomial[index + len(monic) - 1]
        quotient[index] = leading
        for offset, coefficient in enumerate(monic):
            polynomial[index + offset] = (
                polynomial[index + offset] - leading * coefficient
            ) % modulus
    return trim(quotient), trim(polynomial)


def mixed_polynomial(first, second, left, right):
    output = [0] * max(len(first) + len(left) - 1, len(second) + len(right) - 1)
    for parent, multiplier in ((first, left), (second, right)):
        for i, coefficient in enumerate(parent):
            for j, factor in enumerate(multiplier):
                output[i + j] += coefficient * factor
    return output


def lift_block(mixture, residual, maximum_depth):
    found = residual
    modulus = 3
    answers = {}
    assert not divide_mod(mixture, found, modulus)[1]
    for depth in range(2, maximum_depth + 1):
        candidates = [
            [value + modulus * digit for value, digit in zip(found[:-1], digits, strict=True)]
            + [1]
            for digits in product(range(3), repeat=len(found) - 1)
        ]
        modulus *= 3
        lifts = [candidate for candidate in candidates
                 if not divide_mod(mixture, candidate, modulus)[1]]
        assert len(lifts) == 1
        found = lifts[0]
        answers[depth] = found
    return answers


def main():
    rows = Counter({(2, 5): 4, (5, 2): 4})
    primitive = [primitive_parent(order, list(rows.elements())) for order in (8, 9)]
    first_slopes, second_slopes = (newton_slopes(parent, 3) for parent in primitive)
    assert set(first_slopes) & set(second_slopes) == {0}
    parents = []
    for order in (8, 9):
        data = coefficient_data(rows, order, 3, 27)
        minimum = min(value for value, _ in data)
        parents.append([3**(value - minimum) * unit % 27 for value, unit in data])
    assert parents == [
        [16, 0, 22, 23, 0, 23, 22, 0, 16],
        [9, 1, 1, 0, 0, 0, 0, 1, 1, 9],
    ]
    residual = [trim([value % 3 for value in parent]) for parent in parents]
    block = residue_gcd(*residual, 3)
    assert block == [1, 1, 1, 1]
    assert residue_gcd(block, [1, 0, 0], 3) == [1]
    quotient = [divide_mod(parent, block, 3) for parent in residual]
    assert all(not remainder for _, remainder in quotient)
    assert residue_gcd(block, quotient[0][0], 3) == [1, 1]
    assert residue_gcd(block, quotient[1][0], 3) == [1, 0, 1]
    mixture = mixed_polynomial(*parents, [1, 0, -1, 1, -1], [1, 0, 1, 1, 1])
    assert trim([coefficient % 3 for coefficient in mixture]) == block
    linear = lift_block(mixture, [1, 1], 3)
    quadratic = lift_block(mixture, [1, 0, 1], 3)
    full = lift_block(mixture, block, 3)
    assert linear == {2: [7, 1], 3: [16, 1]}
    assert quadratic == {2: [1, 3, 1], 3: [19, 12, 1]}
    assert full == {2: [7, 4, 1, 1], 3: [7, 22, 1, 1]}
    assert [divide_mod(parent, linear[2], 9)[1] for parent in parents] == [[3], [3]]
    assert all(not divide_mod(parent, quadratic[2], 9)[1] for parent in parents)
    assert [divide_mod(parent, quadratic[3], 27)[1] for parent in parents] == [
        [9, 9], [18, 18],
    ]
    assert [divide_mod(parent, full[2], 9)[1] for parent in parents] == [
        [6, 0, 6], [6, 0, 6],
    ]
    print("Verified crossed residual cofactors, six unique lifts, and no shared unit block.")



def fractional_block():
    rows = Counter({(0, 3): 8, (3, 0): 5, (3, 4): 1, (5, 1): 1})
    parents = []
    for order in (5, 6):
        data = coefficient_data(rows, order, 2, 2)
        depths = [2 * valuation + 9 * k for k, (valuation, _) in enumerate(data)]
        minimum = min(depths)
        parents.append([1 << (depth - minimum) if depth - minimum < 2 else 0
                        for depth in depths])
    assert parents == [[1, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1, 0]]

    def multiply_digits(first, second):
        constant = (first & 1) * (second & 1)
        linear = ((first & 1) * (second >> 1) + (first >> 1) * (second & 1)) % 2
        return constant | (linear << 1)

    def remainder_digits(parent, factor):
        coefficients = parent[:]
        for index in range(len(coefficients) - 1, 1, -1):
            leading = coefficients[index]
            for offset, value in enumerate(factor):
                coefficients[index - 2 + offset] ^= multiply_digits(leading, value)
        return coefficients[:2]

    residues = []
    for constant, linear in product((0, 1), repeat=2):
        factor = [1 + 2 * constant, 2 * linear, 1]
        if all(remainder_digits(parent, factor) == [0, 0] for parent in parents):
            residues.append((constant, linear))
    assert residues == [(0, 0), (0, 1), (1, 0), (1, 1)]
    allowed = [digits for digits in residues if all(
        digit == 0 or 1 % 2 == (index - 2) * 9 % 2
        for index, digit in enumerate(digits)
    )]
    assert allowed == [(0, 0), (0, 1)]


if __name__ == "__main__":
    main()
    fractional_block()
    print("Verified four arbitrary first digits and two rational digits for a partial block.")
