#!/usr/bin/env python3
"""Exact normal directions for a multiple Sylvester contact."""

from fractions import Fraction


def multiply(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def subtract(left, right):
    length = max(len(left), len(right))
    left = left + [0] * (length - len(left))
    right = right + [0] * (length - len(right))
    return [a - b for a, b in zip(left, right, strict=True)]


def remainder(polynomial, monic):
    degree = len(monic) - 1
    output = list(map(Fraction, polynomial)) + [0] * max(0, degree - len(polynomial))
    for i in range(len(output) - 1, degree - 1, -1):
        leading = output[i]
        for j, value in enumerate(monic):
            output[i - degree + j] -= leading * value
    return output[:degree]


def rank(matrix):
    matrix = [list(map(Fraction, row)) for row in matrix]
    nonzero = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(nonzero, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[nonzero], matrix[pivot] = matrix[pivot], matrix[nonzero]
        basis_row = matrix[nonzero]
        for i in range(nonzero + 1, len(matrix)):
            factor = matrix[i][column] / basis_row[column]
            matrix[i] = [a - factor * b for a, b in zip(matrix[i], basis_row, strict=True)]
        nonzero += 1
    return nonzero


def weight_factor(coefficient_array, weight):
    degree = len(coefficient_array) - 1
    return [value * weight(degree, k) for k, value in enumerate(coefficient_array)]


def strength(j, k):
    return Fraction(k * (k - j), 2 * j)


def skew(j, k):
    return (Fraction(k, j) - Fraction(1, 2)) * Fraction(k * (j - k), 6 * j)


def inverse_height(j, k):
    return -strength(j, k) / (2 * j)


def quartic(j, k):
    return Fraction(j, 12) * ((Fraction(k, j) - Fraction(1, 2))**4 - Fraction(1, 16))


def check():
    first_pair = [1, 1, 1]
    second_pair = [2, 1, 1]
    square = multiply(first_pair, first_pair)
    repeated_gcd = multiply(square, second_pair)
    cofactor = [1, 1]
    first = repeated_gcd
    second = multiply(repeated_gcd, cofactor)
    parent_derivatives = [
        (weight_factor(first, weight), weight_factor(second, weight))
        for weight in (strength, skew, inverse_height, quartic)
    ]
    parent_derivatives.insert(
        1, ([0], weight_factor(second, lambda j, k: Fraction(-k, 84)))
    )
    vectors = [subtract(multiply(cofactor, lower), upper)
               for lower, upper in parent_derivatives]

    certificates = [
        (square, [Fraction(493, 130), Fraction(726, 65), Fraction(-131, 10), 1],
         Fraction(-121, 1820), Fraction(-498713, 4815720)),
        (multiply(first_pair, second_pair),
         [Fraction(698, 1061), Fraction(363, 1061), Fraction(-1303, 1061), 1],
         Fraction(-121, 12732), Fraction(-2048203, 78607368)),
    ]
    for candidate, normal, constant, slope in certificates:
        columns = [remainder(vector, candidate) for vector in vectors]
        assert rank([row[:3] for row in zip(*columns, strict=True)]) == 3
        dot_products = [sum(a * b for a, b in zip(normal, column, strict=True))
                        for column in columns]
        assert dot_products == [0, 0, 0, constant, slope]
        assert constant < 0 and slope < 0

    for multiplier, expected_rank in (
        ([1], 6), (first_pair, 4), (square, 2),
        (multiply(first_pair, second_pair), 2), (repeated_gcd, 0),
    ):
        columns = [
            remainder([0] * k + multiplier, repeated_gcd)
            for k in range(len(repeated_gcd) - 1)
        ]
        assert rank(list(zip(*columns, strict=True))) == expected_rank
    for (lower, upper), vector in zip(parent_derivatives, vectors, strict=True):
        for k in range(6):
            phi = [0] * k + [1]
            kernel_projection = subtract(
                multiply(lower, multiply(cofactor, phi)), multiply(upper, phi)
            )
            assert remainder(kernel_projection, repeated_gcd) == remainder(
               multiply(vector, phi), repeated_gcd
            )

    bezout_first = subtract([1], first_pair)
    assert remainder(multiply(bezout_first, second_pair), square) == [1, 0, 0, 0]
    common_delta = [Fraction(1, 13), Fraction(2, 13)]
    factor = subtract(square, [-value for value in common_delta])
    cofactor_delta = [Fraction(1, 13), Fraction(1, 13)]
    lower_moved = multiply(factor, subtract(second_pair, [-value for value in cofactor_delta]))
    observed = remainder(multiply(bezout_first, subtract(lower_moved, first)), square)
    exact = remainder(multiply(common_delta, multiply(
        bezout_first, subtract(second_pair, [-value for value in cofactor_delta])
    )), square)
    assert observed == exact
    print("Verified the Sylvester kernel projection, both repeated-divisor gates, and the whole-block response.")



if __name__ == "__main__":
    check()
