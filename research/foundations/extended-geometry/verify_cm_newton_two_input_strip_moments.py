#!/usr/bin/env python3
"""Replay the rank-two and rank-three direct strip congruences."""

from collections import defaultdict


def multiply(left, right, modulus):
    product = defaultdict(int)
    for (left_angle, left_height), left_coefficient in left.items():
        for (right_angle, right_height), right_coefficient in right.items():
            exponent = (left_angle + right_angle, left_height + right_height)
            product[exponent] += left_coefficient * right_coefficient
    return {
        exponent: coefficient % modulus
        for exponent, coefficient in product.items()
        if coefficient % modulus
    }


def power(polynomial, exponent, modulus):
    result = {(0, 0): 1}
    for _ in range(exponent):
        result = multiply(result, polynomial, modulus)
    return result


def normalized_moment_residue(literal_power, prime):
    residue = 0
    for (angle, height), coefficient in literal_power.items():
        if angle != 0 or height % 2:
            continue
        denominator = height + 1
        assert denominator < prime**2
        if denominator % prime == 0:
            residue += coefficient * pow(denominator // prime, -1, prime)
    return residue % prime


def check_parameter(parameter, expected_inputs, expected_residues):
    n, a, s, prime = 10, 5, 3, 17
    order = (s * prime - 1) // n
    first_factor = {(0, 0): -1, (1, 0): 1}
    second_factor = {(0, 0): -parameter, (1, 0): 1}
    angular_polynomial = multiply(
        power(first_factor, 3, prime), power(second_factor, 3, prime), prime
    )
    angular_power = power(angular_polynomial, order, prime)
    assert tuple(angular_power.get((index, 0), 0) for index in
                 (order, order + prime)) == expected_inputs

    phase = {
        (index - 1, n - 2 * (index - 1)): coefficient
        for (index, _), coefficient in angular_polynomial.items()
    }
    phase[-1, n] = -angular_polynomial[0, 0]
    # The lift is regular and these strict lower terms have odd height.
    phase[0, n - 1] = 7
    phase[1, n - 3] = 5
    for angle, height in phase:
        assert height >= 0 and height + 2 * angle <= n
        assert (a + 1) * height + 2 * a * angle >= (a + 1) * n - 2 * a

    residues = tuple(
        normalized_moment_residue(power(phase, exponent, prime), prime)
        for exponent in (order, order + prime)
    )
    assert residues == expected_residues
    print(f"c={parameter}: inputs {expected_inputs}, residues {residues}")


def check_third_input():
    n, a, s, prime = 18, 8, 5, 83
    order = (s * prime - 1) // n
    angular_polynomial = power({(3, 0): 1, (0, 0): -1}, 3, prime)
    angular_power = power(angular_polynomial, order, prime)
    inputs = tuple(angular_power.get((order + k * prime, 0), 0) for k in range(3))
    assert inputs == (0, 0, 74)

    phase = {
        (index - 1, n - 2 * (index - 1)): coefficient
        for (index, _), coefficient in angular_polynomial.items()
    }
    phase[-1, n] = -angular_polynomial[0, 0]
    phase[0, n - 1] = 7
    for angle, height in phase:
        assert height >= 0 and height + 2 * angle <= n
        assert (a + 1) * height + 2 * a * angle >= (a + 1) * n - 2 * a

    base_power = power(phase, order, prime)
    residues = []
    for k in range(3):
        frobenius_power = {
            (prime * angle, prime * height): coefficient
            for (angle, height), coefficient in power(phase, k, prime).items()
        }
        literal_power = multiply(base_power, frobenius_power, prime)
        residues.append(normalized_moment_residue(literal_power, prime))
    assert residues == [0, 0, 46]
    print(f"rank three: inputs {inputs}, residues {tuple(residues)}")


if __name__ == "__main__":
    check_parameter(2, (14, 7), (16, 8))
    check_parameter(-1, (0, 5), (0, 16))
    check_third_input()
