#!/usr/bin/env python3
"""Replay generic-fiber pole orders and low moments of the shear word."""

from fractions import Fraction

from verify_spherical_critical_face_p_power_moments import multiply, power


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, 0) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def scale(polynomial, scalar):
    return {exponent: scalar*coefficient for exponent, coefficient in polynomial.items()}


def shear_x(point, direction):
    x, y, t = point
    return (x, add(y, scale(multiply(x, t), -direction), scale(power(x, 3), Fraction(-1, 2))),
            add(t, scale(power(x, 2), direction)))


def shear_y(point, direction):
    x, y, t = point
    return (add(x, scale(multiply(y, t), -direction), scale(power(y, 3), Fraction(-1, 2))), y,
            add(t, scale(power(y, 2), direction)))


def word(point, direction):
    return shear_x(shear_y(shear_x(point, direction), direction), direction)


def angular_coefficient(polynomial, angle):
    return {(0, degree): coefficient for (term_angle, degree), coefficient in polynomial.items()
            if term_angle == angle}


def check_ends():
    # Exponents here are those of z and of the symbolic fiber value lambda.
    x = {(1, 0): 1}
    y = {(-1, 0): Fraction(1, 2), (-1, 2): Fraction(-1, 2)}
    t = {(0, 1): 1}
    inverse = word((x, y, t), -1)
    largest = tuple(max(angle for angle, _ in polynomial) for polynomial in inverse)
    smallest = tuple(min(angle for angle, _ in polynomial) for polynomial in inverse)
    assert largest == (9, 27, 18)
    assert smallest == (-3, -9, -6)

    constant = angular_coefficient(y, -1)
    infinity_coefficients = (Fraction(1, 16), Fraction(-1, 8192), Fraction(-1, 256))
    zero_coefficients = (scale(power(constant, 3), Fraction(-1, 2)),
                         scale(power(constant, 9), Fraction(1, 16)),
                         scale(power(constant, 6), Fraction(-1, 4)))
    for polynomial, high, low, infinity, zero in zip(inverse, largest, smallest,
                                                     infinity_coefficients, zero_coefficients):
        assert angular_coefficient(polynomial, high) == {(0, 0): infinity}
        assert angular_coefficient(polynomial, low) == zero
    print("Pole orders at infinity: (9, 27, 18); at zero: (3, 9, 6)")


def check_forward():
    # Now the symbols are the original x and t in the literal chart.
    x = {(1, 0): 1}
    y = {(-1, 0): Fraction(1, 2), (-1, 2): Fraction(-1, 2)}
    t = {(0, 1): 1}
    _, _, phase = word((x, y, t), 1)
    restricted = {angle: coefficient for (angle, height), coefficient in phase.items()
                  if height == 0}
    assert max(restricted) == 18 and restricted[18] == Fraction(1, 256)
    expected = (0, Fraction(1, 3))
    for order, target in enumerate(expected, 1):
        moment = sum(coefficient/Fraction(height+1)
                     for (angle, height), coefficient in power(phase, order).items()
                     if angle == 0 and height % 2 == 0)
        assert moment == target
    print("Degree 18 and literal moments 0, 1/3 verified")


if __name__ == "__main__":
    check_ends()
    check_forward()
