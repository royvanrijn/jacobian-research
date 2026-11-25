#!/usr/bin/env python3
"""Exact certificate for an ordinary three-column Gaussian correlation path."""

from functools import cache
from math import comb, prod, factorial

from flint import fmpq, fmpq_poly, fmpz_poly


ROWS = ((5, 2, 3), (2, 5, 3), (5, 4, 4), (4, 5, 4), (5, 5, 4))


def require(condition, label):
    if not condition:
        raise AssertionError(label)


@cache
def raw_moment(i, j, k):
    """Polynomial in the squared correlation for E[U^i V^j W^k]."""
    powers = [i * a + j * b + k * c for a, b, c in ROWS]
    p, q = powers[:2]
    factorial_factor = prod(factorial(power) for power in powers)
    wick = fmpz_poly([comb(p, ell) * comb(q, ell) for ell in range(min(p, q) + 1)])
    return factorial_factor * wick


@cache
def symmetric_moment(i, j, k):
    """Polynomial in the squared correlation for E[(U+V)^i (U-V)^j W^k]."""
    value = fmpz_poly()
    for a in range(i + 1):
        for b in range(j + 1):
            multiplier = (-1) ** (j - b) * comb(i, a) * comb(j, b)
            value += multiplier * raw_moment(a + b, i + j - a - b, k)
    return value


def moments_with_field(i, j):
    """Coefficients in t for E[(U+V+tW)^i (U-V)^j]."""
    return [comb(i, k) * symmetric_moment(i - k, j, k) for k in range(i + 1)]


def multiply(left, right):
    coefficients = [fmpz_poly() for _ in range(len(left) + len(right) - 1)]
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            coefficients[i + j] += first * second
    return coefficients


def weighted_sum(*terms):
    coefficients = [fmpz_poly() for _ in range(max(len(group) for _, group in terms))]
    for factor, group in terms:
        for i, coefficient in enumerate(group):
            coefficients[i] += factor * coefficient
    return coefficients


def at_correlation(polynomial, correlation):
    return fmpq_poly([coefficient(correlation) for coefficient in polynomial])


def at_field(polynomial, field):
    return sum(
        (fmpq_poly(coefficient) * field**i for i, coefficient in enumerate(polynomial)),
        fmpq_poly(),
    )


def sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def common_sign(values):
    signs = {sign(value) for value in values}
    return signs.pop() if len(signs) == 1 else 0


def bernstein_coefficients(polynomial, left, right):
    """Bernstein coefficients after mapping [left, right] to [0, 1]."""
    power = polynomial(fmpq_poly([left, right - left]))
    degree = max(0, power.degree())
    return [
        sum(power[j] * fmpq(comb(k, j), comb(degree, j)) for j in range(k + 1))
        for k in range(degree + 1)
    ]


def rectangle_sign(polynomial, q_left, q_right, t_left, t_right):
    """Prove one strict sign by converting in t and then in correlation."""
    degree = len(polynomial) - 1
    power = [fmpq_poly() for _ in range(degree + 1)]
    for i, coefficient in enumerate(polynomial):
        for j in range(i + 1):
            multiplier = comb(i, j) * t_left ** (i - j) * (t_right - t_left) ** j
            power[j] += fmpq_poly(coefficient) * multiplier
    signs = []
    for k in range(degree + 1):
        beta = sum(
            (power[j] * fmpq(comb(k, j), comb(degree, j)) for j in range(k + 1)),
            fmpq_poly(),
        )
        signs.append(common_sign(bernstein_coefficients(beta, q_left, q_right)))
    return common_sign(signs)


def isolate(polynomial, left, right, steps=22):
    left_sign = sign(polynomial(left))
    require(left_sign * sign(polynomial(right)) == -1, "Root bisection bracket")
    for _ in range(steps):
        midpoint = (left + right) / 2
        midpoint_sign = sign(polynomial(midpoint))
        if midpoint_sign == 0:
            return midpoint, midpoint
        if midpoint_sign == left_sign:
            left = midpoint
        else:
            right = midpoint
    return left, right


def main():
    difference_variance = symmetric_moment(0, 2, 0)
    difference_fourth = symmetric_moment(0, 4, 0)
    second = moments_with_field(2, 0)
    cubic = weighted_sum(
        (difference_variance, moments_with_field(3, 0)),
        (-3, multiply(second, moments_with_field(1, 2))),
    )
    quartic = weighted_sum(
        (difference_variance**2, moments_with_field(4, 0)),
        (-6 * difference_variance, multiply(second, moments_with_field(2, 2))),
        (difference_fourth, multiply(second, second)),
    )

    q_left, q_right = fmpq(21, 25), fmpq(17, 20)
    t_left, t_right = fmpq(-2540), fmpq(-2520)
    require(
        common_sign(bernstein_coefficients(at_field(cubic, t_left), q_left, q_right)) == -1,
        "Cubic negative on left field edge",
    )
    require(
        common_sign(bernstein_coefficients(at_field(cubic, t_right), q_left, q_right)) == 1,
        "Cubic positive on right field edge",
    )
    derivative = [(i + 1) * cubic[i + 1] for i in range(len(cubic) - 1)]
    require(
        rectangle_sign(derivative, q_left, q_right, t_left, t_right) == 1,
        "Cubic increasing on the entire rectangle",
    )

    for correlation, expected in [(q_left, -1), (q_right, 1)]:
        cubic_at_q = at_correlation(cubic, correlation)
        quartic_at_q = at_correlation(quartic, correlation)
        left, right = isolate(cubic_at_q, t_left, t_right)
        require(
            common_sign(bernstein_coefficients(quartic_at_q, left, right)) == expected,
            "Quartic sign throughout the cubic root bracket",
        )
        print(f"Correlation {correlation}, cubic root in [{left}, {right}], quartic sign {expected:+d}")

    for endpoint in (fmpq(0), fmpq(1)):
        require(difference_variance(endpoint) > 0, "Positive endpoint difference variance")
        cubic_at_q = at_correlation(cubic, endpoint)
        quartic_at_q = at_correlation(quartic, endpoint)
        require(cubic_at_q.degree() == 3, "No zero in the r=0 projective chart")
        require(cubic_at_q.gcd(quartic_at_q).degree() == 0, "Endpoint coprimality")
        print(f"Correlation {endpoint}: endpoint forms are projectively regular")

    print("PASS: exact interior physical collision with regular ordinary endpoints")


if __name__ == "__main__":
    main()
