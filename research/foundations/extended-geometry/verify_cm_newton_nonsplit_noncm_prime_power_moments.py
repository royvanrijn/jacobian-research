"""Replay coherent returns, a primitive rank, and rational sphere moments."""

from collections import defaultdict
from fractions import Fraction
from math import comb, gcd


def multiply(left, right):
    result = defaultdict(Fraction)
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            result[i + k, j + ell] += a * b
    return {degree: value for degree, value in result.items() if value}


def coherent_dimension(refined_order, character, exponents):
    residue_sum = sum(character * exponent % refined_order for exponent in exponents)
    quotient, remainder = divmod(residue_sum, refined_order)
    assert remainder == 0
    return quotient - 1


def frobenius_entry(prime, branch_floors, parameter, row=1, column=1):
    g, a, b = branch_floors
    height = prime * column - row - g
    coefficient = sum(
        comb(a, i) * (-1) ** (a - i)
        * comb(b, height - i) * (-parameter) ** (b - height + i)
        for i in range(max(0, height - b), min(a, height) + 1)
    )
    return coefficient % prime


def return_coefficients(prime, iterate, order, parameter):
    coefficients = []
    character = (order + 1) // 2
    first_character = character
    for _ in range(iterate):
        branch_floors = tuple(
            prime * (character * exponent % order) // order
            for exponent in (order - 2, 2, 2)
        )
        coefficients.append(frobenius_entry(prime, branch_floors, parameter))
        character = prime * character % order
    assert character == first_character
    return coefficients


def check_example(prime, iterate, parameter):
    order = (prime**iterate - 1) // (prime - 1)
    assert gcd(order, 6 * prime) == 1 and parameter % 9 in (3, 6)
    exponents = order - 2, 2, 2, order - 2
    odd_orbit = [
        coherent_dimension(2 * order, 2 * prime**j + order, exponents)
        for j in range(iterate)
    ]
    assert odd_orbit[0] == 0 and odd_orbit[-1] == 2
    assert return_coefficients(prime, iterate, order, parameter) == [1] * iterate

    # The regular face lift, expanded in literal angles and heights,
    # with 3y+7/2 added at smaller weight.
    phase = {
        (-2, order - 2): parameter,
        (-2, order): -2 * parameter,
        (-2, order + 2): parameter,
        (0, order): -1 - parameter,
        (2, order - 2): 1,
        (-1, 0): Fraction(3, 2),
        (-1, 2): Fraction(-3, 2),
        (0, 0): Fraction(7, 2),
    }
    powered = {(0, 0): Fraction(1)}
    for _ in range(prime - 1):
        powered = multiply(powered, phase)
    value = sum(
        coefficient / (height + 1)
        for (angle, height), coefficient in powered.items()
        if angle == 0 and height % 2 == 0
    )
    difference = prime**iterate * value - 1
    assert difference.denominator % prime != 0 and difference.numerator % prime == 0
    print(f"p={prime}, N={order}: odd dimensions {odd_orbit}, return and moment residue 1")


def check_rank_at_31():
    prime, iterate, order, parameter = 5, 3, 31, 3
    refined_order = 2 * order
    exponents = (order - 2, 2, 2, order - 2)
    primitive_labels = {
        label for label in range(1, refined_order) if gcd(label, order) == 1
    }
    rank_by_parity = [0, 0]
    while primitive_labels:
        first = min(primitive_labels)
        orbit = tuple(first * prime**j % refined_order for j in range(iterate))
        assert len(set(orbit)) == iterate and set(orbit) <= primitive_labels
        assert prime * orbit[-1] % refined_order == first
        primitive_labels.difference_update(orbit)

        dimensions = [coherent_dimension(refined_order, label, exponents) for label in orbit]
        if 0 in dimensions:
            continue
        dimension = dimensions[0]
        assert dimension == (2 if first % 2 else 1)
        assert dimensions == [dimension] * iterate
        for label in orbit:
            branch_floors = tuple(
                prime * (label * exponent % refined_order) // refined_order
                for exponent in exponents[:-1]
            )
            indices = range(1, dimension + 1)
            step = [
                [frobenius_entry(prime, branch_floors, parameter, row, column) for column in indices]
                for row in indices
            ]
            determinant = step[0][0]
            if dimension == 2:
                determinant = step[0][0] * step[1][1] - step[0][1] * step[1][0]
            assert determinant % prime != 0
        rank_by_parity[first % 2] += iterate * dimension

    digit_count = sum(
        gcd(i, order) == 1 and all(4 * i // prime**j % prime >= 2 for j in range(iterate))
        for i in range(1, order)
    )
    assert digit_count == 6 and rank_by_parity == [30, 2 * digit_count]
    print(f"p={prime}, N={order}: 5-rank {sum(rank_by_parity)} from coherent orbits, "
          f"digit count {digit_count}")


def check_iterated_return_at_31():
    prime, order, parameter = 5, 31, 3
    moment_order = (prime**6 - 1) // order
    # The top coefficient is [u^m](u-1)^m(u-c)^m. Expanding both
    # factors directly gives the central binomial convolution.
    coefficient = sum(
        comb(moment_order, j)**2 * pow(parameter, j, prime)
        for j in range(moment_order + 1)
    ) % prime
    assert moment_order == 504 and coefficient == 1
    print(f"p={prime}, N={order}: order {moment_order} top coefficient {coefficient}")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    for prime, iterate, parameter in ((5, 3, 3), (5, 5, 3), (11, 3, 6), (17, 3, 3)):
        check_example(prime, iterate, parameter)
    check_rank_at_31()
    check_iterated_return_at_31()


if __name__ == "__main__":
    main()
