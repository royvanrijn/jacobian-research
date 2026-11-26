#!/usr/bin/env python3
"""Exact arithmetic checks for the populated upper height certificates."""

from collections import Counter
from fractions import Fraction
from math import comb, isqrt, prod
import json


J = 12
T = 90


def prime(value):
    return value > 1 and all(value % d for d in range(2, isqrt(value) + 1))


def valuation(value, p):
    exponent = 0
    while value % p == 0:
        value //= p
        exponent += 1
    return exponent


def factorial_valuation(value, p):
    total = 0
    while value:
        value //= p
        total += value
    return total


def favorable(rows, height):
    starts = Counter(sign for a, _, sign in rows if a == height)
    straddling = Counter(sign for a, d, sign in rows if a <= height < a+d)
    return next((sign for sign in (-1, 1) if starts[sign] > straddling[-sign]), 0)


def ordinates(rows, order, p, direction):
    return [
        valuation(comb(order, k), p) + sum(
            factorial_valuation(order*a+d*(k if sign == direction else order-k), p)
            - factorial_valuation(order*a, p)
            for a, d, sign in rows
        )
        for k in range(order+1)
    ]


def replay(rows, heights, counts):
    maximum = max(a+d for a, d, _ in rows)
    assert all(2*J*height >= (2*J-1)*maximum for height in heights)
    for terminal in range(2, J):
        height = next(
            a for a in heights if favorable(rows, a) and prime(terminal*a+1)
        )
        p = terminal*height+1
        direction = favorable(rows, height)
        assert p > (terminal-1)*maximum and 2*p > terminal*maximum
        upper = ordinates(rows, terminal, p, direction)
        lower = ordinates(rows, terminal-1, p, direction)
        assert not any(lower)
        assert min(upper[1:]) > upper[0]
        counts["private_height_certificates"] += 1
        counts["valuation_entries"] += len(upper)+len(lower)


def multiplier(a, b, order, index):
    rising = prod(range(order*a+1, order*a+index+1))
    falling = prod(range(order*b+order-index+1, order*b+order+1))
    return Fraction(rising, falling)


def replay_pair_coefficient(a, b, order, index):
    increasing = prod(range(order*a+1, order*a+index+1))
    dec_short = prod(range(order*b+1, order*b+order-index+1))
    dec_full = prod(range(order*b+1, order*b+order+1))
    assert multiplier(a, b, order, index) == Fraction(increasing*dec_short, dec_full)


def fixtures():
    for maximum in (7000, 10000, 18000):
        original = [
            (0, maximum, 1),
            (17, maximum-17, -1),
            (maximum//3, maximum//4, -1),
            (maximum-2, 2, 1),
        ]
        yield maximum, original


def run():
    counts = Counter()
    for M, original in fixtures():
        selected_heights = [
            next(a for a in range(M+T, M+2*T) if prime(order*a+1))
            for order in range(2, J)
        ]
        selected = [(a, 1, 1) for a in selected_heights] + [(M, 1, -1)]*len(selected_heights)
        replay(original+selected, selected_heights, counts)
        assert sum(d*sign for _, d, sign in selected) == 0

        first_even = 2*((M+2)//2)
        predetermined_heights = [first_even+2*r for r in range(T)]
        predetermined = (
            [(a, 1, 1) for a in predetermined_heights] +
            [(a-1, 1, -1) for a in predetermined_heights]
        )
        replay(original+predetermined, predetermined_heights, counts)
        assert sum(d*sign for _, d, sign in predetermined) == 0

        # A moving native profile has alternating favorable directions
        # despite wide oppositely oriented ramps straddling its high sites.
        mixed_heights = list(range(M-2*T, M-T))
        mixed = original + [
            (a, 1, 1 if a % 2 else -1) for a in mixed_heights for _ in range(3)
        ]
        assert {favorable(mixed, a) for a in mixed_heights} == {-1, 1}
        replay(mixed, mixed_heights, counts)
        counts["moving_mixed_profiles"] += 1

        for order in range(2, J):
            for k in range(order+1):
                for a in selected_heights:
                    replay_pair_coefficient(a, M, order, k)
                    assert 1 <= multiplier(a, M, order, k) <= Fraction(M+2*T, M)**k
                    counts["bounded_extension_multiplier_checks"] += 1
                for a in (predetermined_heights[0], predetermined_heights[-1]):
                    replay_pair_coefficient(a, a-1, order, k)
                    assert 1 <= multiplier(a, a-1, order, k) <= (
                        1 + Fraction(k, order*a-k+1)
                    )**k
                    counts["predetermined_multiplier_checks"] += 1
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
