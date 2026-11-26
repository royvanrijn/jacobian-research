#!/usr/bin/env python3
"""Reproduce the large-prime first-position quadratic control.

Requires python-flint (replayed with version 0.9.0). Both block sizes
compute the complete finite sum by exact arithmetic modulo p. The
product tree and multipoint remainder tree give soft-square-root
complexity in p when the block size is of order sqrt(p).
"""

from math import comb, factorial, isqrt

from flint import nmod_poly


def consecutive_sum(prime, argument, start, stop, term=1):
    total = 0
    for index in range(start, stop):
        total = (total + term) % prime
        multiplier = argument * (1 - 4 * index * index) % prime
        term = term * multiplier % prime
    return term, total


def block_polynomials(prime, argument, size):
    """Return the final term and the sum of incoming terms of a block."""
    one = nmod_poly([1], prime)
    blocks = [
        (nmod_poly([
            argument * (1 - 4 * index * index) % prime,
            -8 * argument * index % prime,
            -4 * argument % prime,
        ], prime), one)
        for index in range(size)
    ]
    while len(blocks) > 1:
        merged = []
        for i in range(0, len(blocks) - 1, 2):
            left_product, left_sum = blocks[i]
            right_product, right_sum = blocks[i + 1]
            merged.append((
                left_product * right_product,
                left_sum + left_product * right_sum,
            ))
        if len(blocks) % 2:
            merged.append(blocks[-1])
        blocks = merged
    return blocks[0]


def point_tree(prime, step, count):
    levels = [[nmod_poly([-step * k % prime, 1], prime) for k in range(count)]]
    while len(levels[-1]) > 1:
        previous = levels[-1]
        merged = [previous[i] * previous[i + 1] for i in range(0, len(previous) - 1, 2)]
        if len(previous) % 2:
            merged.append(previous[-1])
        levels.append(merged)
    return levels


def evaluate(polynomial, levels):
    remainders = [polynomial % levels[-1][0]]
    for level in reversed(levels[:-1]):
        remainders = [remainders[i // 2] % divisor for i, divisor in enumerate(level)]
    return [int(value[0]) if value else 0 for value in remainders]


def factorial_moment(prime, argument, block_size):
    half = (prime - 1) // 2
    # If a_i = binom(half + 1, i) (2i)! argument^i, then
    # a_{i+1} = argument (1 - 4i^2) a_i modulo prime.
    # The multiplier at i = half is zero, so a_{half + 1} vanishes.
    count = (half + 1) // block_size
    if count == 0:
        return consecutive_sum(prime, argument, 0, half + 1)[1]

    block_product, block_sum = block_polynomials(prime, argument, block_size)
    levels = point_tree(prime, block_size, count)
    products = evaluate(block_product, levels)
    sums = evaluate(block_sum, levels)

    term, total = 1, 0
    for index, (product, subtotal) in enumerate(zip(products, sums, strict=True)):
        if index in {0, count // 2, count - 1}:
            start = index * block_size
            assert (product, subtotal) == consecutive_sum(
                prime, argument, start, start + block_size
            )
        total = (total + term * subtotal) % prime
        term = term * product % prime
    terminal, remainder = consecutive_sum(prime, argument, count * block_size, half + 1, term)
    assert terminal == 0
    return (total + remainder) % prime


def main():
    for prime in (7, 11, 13, 101, 1009, 10007):
        half = (prime - 1) // 2
        argument = 2
        direct = sum(
            comb(half + 1, i) * factorial(2 * i) * pow(argument, i, prime)
            for i in range(half + 2)
        ) % prime
        assert factorial_moment(prime, argument, isqrt(half) + 1) == direct

    prime = 5_944_184_339
    assert all(prime % divisor for divisor in range(3, isqrt(prime) + 1, 2))
    w = -3_369_377_549 * pow(81_133_029_162, -1, prime) % prime
    g = 8 * pow(15 * (1 + 4 * w + 24 * w * w), -1, prime) % prime
    argument = (w + 2 * g * g) % prime
    assert (w, g, argument) == (43_405_095, 4_250_342_174, 1_951_318_268)
    for block_size in (65_536, 45_017):
        result = factorial_moment(prime, argument, block_size)
        assert result == 1_059_228_468
        print(f"Block {block_size} verified g_{{(p+1)/2}}({argument}) = {result} mod {prime}")


if __name__ == "__main__":
    main()
