#!/usr/bin/env python3
"""Exact finite checks for the private-prime certificate and counting identity."""

from fractions import Fraction
from itertools import product
from math import comb, factorial, isqrt, prod
import json


WINDOWS = (
    range(200, 207),
    range(202, 209),
    range(205, 210),
    range(198, 204),
)
SIGNED_GAPS = (1, -2, 3, -1)
CUTOFF = 9
MIN_PRIME = CUTOFF**2 * max(map(abs, SIGNED_GAPS)) + CUTOFF + 1


def is_prime(value):
    return value > 1 and all(value % divisor for divisor in range(2, isqrt(value) + 1))


def ramp(base, signed_gap, orders):
    return tuple(
        order * base + offset
        for order in orders
        for offset in range(1, order * abs(signed_gap) + 1)
    )


def valuation(value, prime):
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def factorial_valuation(value, prime):
    exponent = 0
    while value:
        value //= prime
        exponent += value
    return exponent


def ordinates(rows, order, prime):
    result = []
    for index in range(order + 1):
        ordinate = valuation(comb(order, index), prime)
        for base, signed_gap in rows:
            draw = index if signed_gap > 0 else order - index
            start = order * base
            stop = start + draw * abs(signed_gap)
            ordinate += factorial_valuation(stop, prime) - factorial_valuation(start, prime)
        result.append(ordinate)
    return result


def slopes(array):
    hull = []
    for index, value in enumerate(array):
        while len(hull) >= 2:
            left, middle = hull[-2:]
            old = Fraction(middle[1] - left[1], middle[0] - left[0])
            new = Fraction(value - middle[1], index - middle[0])
            if old < new:
                break
            hull.pop()
        hull.append((index, value))
    return {
        Fraction(right[1] - left[1], right[0] - left[0])
        for left, right in zip(hull, hull[1:])
    }


def run():
    counts = dict(tuples=0, certifying_primes=0, fixed_prime_removals=0,
                  exact_expectations=0, factorial_integral_replays=0)
    for terminal in range(2, CUTOFF):
        orders = (terminal - 1, terminal)
        candidates = tuple(
            {
                base: terminal * base + 1
                for base in window
                if terminal * base + 1 >= MIN_PRIME and is_prime(terminal * base + 1)
            }
            for window in WINDOWS
        )
        first_offer = next(
            base for offers in candidates for base in offers
        )
        fixed_rows = ((13, 2), (19, -3), (first_offer, 1), (first_offer, 1))
        fixed_factors = tuple(
            factor for row in fixed_rows for factor in ramp(*row, orders)
        )
        private_frequencies = {False: 0, True: 0}
        for bases in product(*WINDOWS):
            rows = tuple(zip(bases, SIGNED_GAPS))
            factors = tuple(ramp(*row, orders) for row in rows)
            counts["tuples"] += 1
            for index, (base, signed_gap) in enumerate(rows):
                if base not in candidates[index]:
                    continue
                prime = candidates[index][base]
                if any(value % prime == 0 for i, strip in enumerate(factors)
                       if i != index for value in strip):
                    continue
                private_frequencies[False] += 1
                if any(value % prime == 0 for value in fixed_factors):
                    counts["fixed_prime_removals"] += 1
                    continue
                private_frequencies[True] += 1
                counts["certifying_primes"] += 1
                left = ordinates(rows + fixed_rows, orders[0], prime)
                right = ordinates(rows + fixed_rows, terminal, prime)
                endpoint = 0 if signed_gap > 0 else terminal
                assert left == [0] * terminal
                assert right == [int(index != endpoint) for index in range(terminal + 1)]
                assert slopes(left) == {0}
                assert slopes(right) == {Fraction(1 if signed_gap > 0 else -1, terminal)}
                if counts["factorial_integral_replays"] < 12:
                    for index, expected in enumerate(right):
                        value = comb(terminal, index)
                        for a, d in rows + fixed_rows:
                            offset = index if d > 0 else terminal - index
                            value *= factorial(terminal * a + offset * abs(d)) // factorial(terminal * a)
                        assert valuation(value, prime) == expected
                    counts["factorial_integral_replays"] += 1
        tuple_count = prod(map(len, WINDOWS))
        for use_fixed in (False, True):
            mean = Fraction()
            for index, offers in enumerate(candidates):
                for prime in offers.values():
                    if use_fixed and any(value % prime == 0 for value in fixed_factors):
                        continue
                    probability = Fraction(1, len(WINDOWS[index]))
                    for other, (window, gap) in enumerate(zip(WINDOWS, SIGNED_GAPS)):
                        if other == index:
                            continue
                        avoid = sum(
                            all(value % prime != 0 for value in ramp(base, gap, orders))
                            for base in window
                        )
                        probability *= Fraction(avoid, len(window))
                    mean += probability
            assert mean == Fraction(private_frequencies[use_fixed], tuple_count)
            counts["exact_expectations"] += 1
    assert counts["certifying_primes"] > 0
    assert counts["fixed_prime_removals"] > 0
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
