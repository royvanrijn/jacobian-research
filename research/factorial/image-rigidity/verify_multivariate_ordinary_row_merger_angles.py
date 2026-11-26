"""Exact arithmetic checks for the ordinary row-merger angle comparison."""

from fractions import Fraction
from math import comb, factorial, prod


def multiply(left, right, length):
    result = [Fraction(0)] * length
    for i, a in enumerate(left):
        for j, b in enumerate(right[: length - i]):
            result[i + j] += a * b
    return result


def inverse_power(a, b, c, power, length):
    # Expand (a + b*s + c*s**2)**(-power) by the binomial series.
    shift = [Fraction(0), b / a, c / a]
    expansion = [Fraction(0)] * length
    term = [Fraction(1)]
    for k in range(length):
        weight = (-1) ** k * comb(power + k - 1, k) / a**power
        for i, coefficient in enumerate(term):
            expansion[i] += weight * coefficient
        term = multiply(term, shift, length)
    return expansion


def lower_exponents(order, t):
    radius = 3 * t + Fraction(7, 4)
    center = (
        Fraction(8 * order, 3) * radius**2
        + Fraction(2, 3) * radius
        - Fraction(1, 12)
    )
    left = center + radius - Fraction(1, 2)
    right = 3 * center - radius - Fraction(3, 2)
    assert left.denominator == right.denominator == 1
    assert left >= 0 and right >= 0
    return int(left), int(right)


def check_lattice_series(order):
    a = Fraction(8 * order**2, 3)
    u = Fraction(2 * order, 3)
    g = Fraction(-7 * order, 12)
    split = [(u + order, g + 1)] + [
        (u - Fraction(order, 3), g + Fraction(h, 3))
        for h in range(1, 4)
    ]
    merged = [(u, g + Fraction(h, 4)) for h in range(1, 5)]
    leading = {
        1: [0, 0, 0, 0, (-g / 2 - Fraction(67, 144)) / a**3],
        2: [0, 0, 1 / (2 * a**3)],
        3: [0, 0, 3 / (2 * a**4)],
        4: [0, 0, 3 / a**5],
    }
    for power, expected in leading.items():
        length = len(expected)
        delta = [Fraction(0)] * length
        for sign, pairs in ((1, split), (-1, merged)):
            for b, c in pairs:
                for i, term in enumerate(inverse_power(a, b, c, power, length)):
                    delta[i] += sign * term
        assert delta == expected, (order, power, delta)
    b2 = leading[1][-1] - Fraction(order - 1, 2) * leading[2][-1]
    b4 = leading[3][-1] / 6
    theta = 2 * order - 3 + 4 * (2 / a) * b2 / b4
    assert b2 == Fraction(6 * order - 31, 144) / a**3
    assert b4 == 1 / (4 * a**4)
    assert theta == Fraction(10 * order, 3) - Fraction(89, 9)
    if order == 4:
        assert (theta - 3) ** 2 < 6
        assert b2 < 0
    else:
        assert order - 1 < theta < Fraction(361, 100) * (order - 20)
        assert b2 > 0
    for t in (0, 1, 100):
        lower_exponents(order, t)


def quartic_coefficients(rows):
    endpoint = prod(factorial(4 * a) for a, _ in rows)
    return tuple(
        Fraction(comb(4, k) * prod(factorial((4 - k) * a + k * b) for a, b in rows),
                 endpoint)
        for k in (1, 2)
    )


def check_quartic():
    left, right = lower_exponents(4, 0)
    assert (left, right) == (35, 98)
    split = [(left, left + 1), (left + 1, left), (right, right + 3), (right + 3, right)]
    merged = [(left + right, left + right + 4), (left + right + 4, left + right)]
    a_split, b_split = quartic_coefficients(split)
    a_merged, b_merged = quartic_coefficients(merged)
    # For z^(-2) F(z), the real variable is x=z+1/z and the
    # quadratic is x**2+a*x+b-2. The vertex of the merged quadratic
    # lies strictly between -2 and 2, where its endpoint signs are
    # positive. The split quadratic is negative at that vertex.
    vertex = -a_merged / 2
    assert -2 < vertex < 2
    assert 4 - 2 * a_merged + b_merged - 2 > 0
    assert 4 + 2 * a_merged + b_merged - 2 > 0
    assert vertex**2 + a_merged * vertex + b_merged - 2 < 0
    assert vertex**2 + a_split * vertex + b_split - 2 < 0
    # The split-minus-merged quadratic is positive for x >= -2.
    # Thus each of its roots is strictly inside the corresponding
    # half of the merged root interval.
    assert a_split > a_merged
    assert b_split - b_merged > 2 * (a_split - a_merged)
    ratio_a = a_split / a_merged
    ratio_b = b_split / b_merged
    assert 1 < ratio_a < ratio_b < ratio_a**2


if __name__ == "__main__":
    check_quartic()
    orders = (4, 304, 310, 1000)
    for order in orders:
        check_lattice_series(order)
    print("Exact quartic fixture and four lattice series at each of four orders passed.")

