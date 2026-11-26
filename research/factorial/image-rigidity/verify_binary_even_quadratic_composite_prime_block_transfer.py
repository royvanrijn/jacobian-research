"""Replay block congruences using exact factorials and ramified valuations."""

from fractions import Fraction
from math import comb, factorial


def multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def homogeneous_powers(q, degree):
    powers = [[1]]
    for _ in range(degree):
        powers.append(multiply(powers[-1], q))
    return powers


def moment(coefficients, x_shift=0, y_shift=0):
    degree = len(coefficients) - 1
    return sum(
        value * factorial(i + x_shift) * factorial(degree - i + y_shift)
        for i, value in enumerate(coefficients)
    )


def remove_power(value, p):
    valuation = 0
    while value % p == 0:
        value //= p
        valuation += 1
    return valuation


def ramified_weight(coefficient, k_power, p):
    if coefficient == 0:
        return float("inf")
    return p * (
        remove_power(coefficient.numerator, p) - remove_power(coefficient.denominator, p)
    ) + 2 * k_power


def ramified_terms(order, homogeneous, p, strength, m):
    terms = [Fraction(0)] * p
    for i in range(order + 1):
        blocks, k_power = divmod(order - i, p)
        # k^p = p^2 * strength in this test field.
        terms[k_power] += Fraction(
            comb(order, i) * homogeneous[i] * (p**2 * strength) ** blocks, p ** (2 * m)
        )
    return terms


def check_homogeneous_blocks(q, p, m):
    h = (p - 1) // 2
    powers = homogeneous_powers(q, m * p + 3)
    values = [moment(coefficients) for coefficients in powers]
    block_powers = homogeneous_powers([a**p for a in q], m)
    block_moments = [moment(coefficients) for coefficients in block_powers]
    lift = powers[p][:]
    for i, value in enumerate(q):
        lift[i * p] -= value**p
    assert all(value % p == 0 for value in lift)
    lift = [value // p for value in lift]
    corrections = []
    for r in range(h):
        product = multiply(lift, powers[r])
        weights = [
            (u, factorial(u) * factorial(p + 2 * r - u))
            for u in range(2 * r + 1, p)
        ]
        corrections.append((
            sum(value * product[p + u] for u, value in weights),
            sum(value * product[u] for u, value in weights),
        ))

    for i, observed in enumerate(values):
        high, low = divmod(i, p)
        expected = block_moments[high] * values[low]
        if high and low < h:
            left, right = corrections[low]
            expected -= high * (
                left * moment(block_powers[high - 1], x_shift=1)
                + right * moment(block_powers[high - 1], y_shift=1)
            )
        assert (observed - p ** (2 * high) * expected) % p ** (2 * high + 1) == 0
    return values, block_moments


def check_transfers(p, m, values, block_moments):
    h = (p - 1) // 2
    for strength in (-2, 1, 5):
        lower = sum(
            comb(m - 1, i) * block_moments[i] * strength ** (m - 1 - i)
            for i in range(m)
        )
        preceding = ramified_terms(m * p - 1, values, p, strength, m - 1)
        for r in range(h, p):
            preceding[p - 1 - r] -= lower * comb(p - 1, r) * values[r]
        assert all(ramified_weight(value, i, p) >= p for i, value in enumerate(preceding))
        for d in range(4):
            order = m * p + d
            interior = ramified_terms(order, values, p, strength, m)
            for r in range(h + 1, p):
                k_power = p + d - r
                interior[k_power] -= Fraction(
                    -m * factorial(d) * values[r] * lower, p * factorial(r) * factorial(k_power)
                )
            interior[0] -= sum(
                comb(m, i) * strength ** (m - i) * Fraction(values[i * p + d], p ** (2 * i))
                for i in range(m + 1)
            )
            assert all(ramified_weight(value, i, p) > 0 for i, value in enumerate(interior))
            for axis in (1, 2, p):
                axis_values = [axis**i * factorial(2 * i) for i in range(order + 1)]
                boundary = ramified_terms(order, axis_values, p, strength, m)
                boundary[0] -= factorial(2 * d) * axis**d * sum(
                    comb(m, i) * factorial(2 * i) * axis ** (i * p) * strength ** (m - i)
                    for i in range(m + 1)
                )
                assert all(ramified_weight(value, i, p) > 0 for i, value in enumerate(boundary))


if __name__ == "__main__":
    total = 0
    for p, m in ((11, 2), (13, 2), (13, 3), (17, 3), (17, 5), (19, 6)):
        for q in ((0, 1, 0), (0, 2, 1), (3, 2, 1), (1, -2, 1), (0, 3, -1)):
            values, block_moments = check_homogeneous_blocks(q, p, m)
            check_transfers(p, m, values, block_moments)
            total += len(values)
    print(f"Checked {total:,} homogeneous blocks and 1,530 ramified moment and boundary transfers.")
