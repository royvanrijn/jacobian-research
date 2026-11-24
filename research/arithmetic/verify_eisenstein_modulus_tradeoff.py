#!/usr/bin/env python3
"""Certified limiting signs for the Eisenstein roughness/modulus trade-off."""

from fractions import Fraction

from verify_eisenstein_divisor_modulus import HALF, divide_interval, integrate
from verify_eisenstein_factor_moments import factor_constant, linear_interval, show_interval


def profile_constants(eta):
    constants = [factor_constant(eta, j) for j in range(1, 5)]
    total = linear_interval(constants, [1, HALF, Fraction(1, 4), Fraction(1, 8)])
    return constants, total


def triple_marginal(eta, lower, upper):
    if upper <= lower:
        return Fraction(0), Fraction(0)
    return integrate(lower, upper, lambda s: (1 - eta - s) / eta)


def square_root_numerator(eta, constants):
    large = triple_marginal(eta, HALF, 1 - 2 * eta)
    return linear_interval(
        [*constants, large],
        [1, -HALF, -Fraction(5, 4), -Fraction(1, 8), 1],
    )


def signed_numerator(eta, theta, constants):
    baseline = square_root_numerator(eta, constants)
    selected = triple_marginal(eta, 1 - theta, HALF)
    charged = triple_marginal(eta, HALF, theta)
    return linear_interval([baseline, selected, charged], [1, 1, -1])


def verify_interval(label, parameter, margin, sign):
    if sign < 0:
        assert margin[1] < 0
    else:
        assert margin[0] > 0
    print(label, parameter, "margin/total", show_interval(margin))


def main():
    for eta, sign in ((Fraction(107, 500), -1), (Fraction(43, 200), 1),
                      (Fraction(11, 50), 1)):
        constants, total = profile_constants(eta)
        margin = divide_interval(square_root_numerator(eta, constants), total)
        verify_interval("square-root eta", eta, margin, sign)

    eta = Fraction(1, 5)
    constants, total = profile_constants(eta)
    for theta, sign in ((Fraction(591, 1000), -1), (Fraction(74, 125), 1),
                        (Fraction(119, 200), 1), (Fraction(3, 5), 1)):
        margin = divide_interval(signed_numerator(eta, theta, constants), total)
        verify_interval("signed modulus", theta, margin, sign)

    # A balanced triple saving its single high pair at theta=3/5.
    balanced = 7 * 13 * 37
    assert 37**2 < balanced and 91**2 > balanced and 91**5 < balanced**3
    assert 7**5 > balanced + 31
    # A triple paying for a low pair at the same modulus.
    unbalanced = 7 * 13 * 97
    assert 97**2 > unbalanced and 91**2 < unbalanced and 91**5 > unbalanced**2
    assert 7**5 > unbalanced + 31


if __name__ == "__main__":
    main()
