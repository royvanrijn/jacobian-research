"""Exact checks of the factorial normalization and the Cartier matrix pairing."""

from collections import defaultdict
from math import factorial, prod


def compositions(total, width):
    if width == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, width - 1):
            yield (first,) + tail


def indices(total, width, *, cumulative=False):
    sums = range(width, total + 1) if cumulative else (total,)
    return [
        tuple(value + 1 for value in index)
        for size in sums
        for index in compositions(size - width, width)
    ]


def poly_power(polynomial, exponent, modulus=None):
    width = len(next(iter(polynomial)))
    out = {(0,) * width: 1}
    for _ in range(exponent):
        new = defaultdict(int)
        for x, a in out.items():
            for y, b in polynomial.items():
                position = tuple(i + j for i, j in zip(x, y))
                new[position] += a * b
        if modulus is None:
            out = {key: value for key, value in new.items() if value}
        else:
            out = {key: value % modulus for key, value in new.items() if value % modulus}
    return out


def coefficient_matrix(polynomial, power, rows, columns):
    return [
        [
            polynomial.get(tuple(power * i - j for i, j in zip(row, column)), 0)
            for column in columns
        ]
        for row in rows
    ]


def multiply(left, right, prime):
    return [
        [sum(a * b for a, b in zip(row, column)) % prime for column in zip(*right)]
        for row in left
    ]


def unit(value, prime):
    valuation = 0
    while value and value % prime == 0:
        value //= prime
        valuation += 1
    return valuation, value % prime


def factorial_weight(index):
    return prod(factorial(value - 1) for value in index)


def factorial_moment_check(polynomial, prime, depth, offset, sigma):
    width = len(sigma)
    degree = min(sum(index) for index in polynomial)
    assert sum(sigma) == degree and min(sigma) > 0
    assert prime > max(offset * degree, max(map(sum, polynomial)))
    power = prime ** depth
    order = power - offset
    coefficients = poly_power(polynomial, order)
    normalization = prod(factorial(order * value) for value in sigma)
    baseline_valuation, baseline_unit = unit(normalization, prime)
    relative_scale = prod(
        factorial(offset * value - 1)
        * pow(factorial(value - 1), -1, prime)
        for value in sigma
    ) % prime
    target = indices(degree, width)
    inputs = indices(offset * degree, width, cumulative=True)
    surviving = {
        tuple(power * i - j for i, j in zip(t, z)): (
            relative_scale * (-1) ** (offset * degree - sum(z))
            * factorial_weight(t) * pow(factorial_weight(z), -1, prime)
        ) % prime
        for t in target
        for z in inputs
    }
    factorials = {i: factorial(i) for i in range(order * max(map(sum, polynomial)) + 1)}
    literal = 0
    term_count = 0
    for index, coefficient in coefficients.items():
        exponent_factorial = prod(factorials[i] for i in index)
        valuation, coefficient_unit = unit(exponent_factorial, prime)
        assert valuation >= baseline_valuation
        ratio_residue = 0
        if valuation == baseline_valuation:
            ratio_residue = coefficient_unit * pow(baseline_unit, -1, prime) % prime
        assert ratio_residue == surviving.get(index, 0)
        literal += coefficient * ratio_residue
        term_count += 1

    observed = coefficient_matrix(coefficients, power, target, inputs)
    lowest = {x: c for x, c in polynomial.items() if sum(x) == degree}
    A = coefficient_matrix(poly_power(lowest, prime - 1, prime), prime, target, target)
    folded = coefficient_matrix(
        poly_power(polynomial, prime - offset, prime), prime, target, inputs
    )
    for _ in range(depth - 1):
        folded = multiply(A, folded, prime)
    assert [[c % prime for c in row] for row in observed] == folded
    weighted = sum(
        relative_scale * factorial_weight(t) * folded[i][j]
        * (-1) ** (offset * degree - sum(z))
        * pow(factorial_weight(z), -1, prime)
        for i, t in enumerate(target)
        for j, z in enumerate(inputs)
    )
    assert literal % prime == weighted % prime

    moment = sum(coefficient * prod(factorials[i] for i in index)
                 for index, coefficient in coefficients.items())
    modulus_power = prime ** baseline_valuation
    quotient, remainder = divmod(moment, modulus_power)
    assert not remainder
    assert quotient * pow(baseline_unit, -1, prime) % prime == weighted % prime
    return term_count, weighted % prime


def point_count(polynomial, prime):
    count = 0
    width = len(next(iter(polynomial)))
    for first in range(width):
        for tail in product_residues(prime, width - first - 1):
            point = (0,) * first + (1,) + tail
            value = sum(
                c * prod(pow(x, e, prime) for x, e in zip(point, exponents))
                for exponents, c in polynomial.items()
            )
            count += value % prime == 0
    return count


def product_residues(prime, width):
    if width == 0:
        yield ()
        return
    for x in range(prime):
        for tail in product_residues(prime, width - 1):
            yield (x,) + tail


def main():
    cubic = {(3, 0, 0): 1, (0, 3, 0): 1, (0, 0, 3): 1}
    cubic_tail = cubic | {(3, 1, 0): 2}
    quartic = {(4, 0, 0, 0): 1, (0, 4, 0, 0): 1, (0, 0, 4, 0): 1, (0, 0, 0, 4): 1}
    quartic_tail = quartic | {(0, 2, 2, 1): -1}
    mixed = {(2, 1): 1, (1, 2): 1, (1, 3): -2}
    mixed_four = {(4, 0): 1, (2, 2): -2, (0, 4): 3, (4, 1): 1, (2, 5): 2}
    profiles = (
        (cubic, 5, 1, 1, (1, 1, 1)),
        (cubic, 7, 2, 1, (1, 1, 1)),
        (cubic_tail, 5, 2, 1, (1, 1, 1)),
        (cubic_tail, 7, 1, 2, (1, 1, 1)),
        (mixed, 5, 2, 1, (1, 2)),
        (mixed, 7, 2, 1, (2, 1)),
        (mixed, 7, 2, 2, (2, 1)),
        (mixed_four, 11, 1, 2, (1, 3)),
        (mixed_four, 11, 1, 2, (2, 2)),
        (mixed_four, 11, 2, 1, (3, 1)),
        (quartic_tail, 7, 1, 1, (1, 1, 1, 1)),
    )
    count = 0
    residues = []
    for polynomial, prime, depth, offset, sigma in profiles:
        tested, residue = factorial_moment_check(polynomial, prime, depth, offset, sigma)
        count += tested
        residues.append(residue)
    assert residues[4] == residues[5] == 0
    for prime in (5, 7, 11, 13, 17, 19):
        scalar = poly_power(cubic, prime - 1, prime).get((prime - 1,) * 3, 0)
        assert scalar == (1 - point_count(cubic, prime)) % prime
    for prime in (5, 7, 11, 13):
        scalar = poly_power(quartic, prime - 1, prime).get((prime - 1,) * 4, 0)
        assert scalar == (point_count(quartic, prime) - 1) % prime
    print(f"{count} exact factorial weights and {len(profiles)} normalized matrix pairings")
    print("Ten cubic and quartic point counts and the nonhomogeneous cancellation witness checked")


if __name__ == "__main__":
    main()
