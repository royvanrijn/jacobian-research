#!/usr/bin/env python3
"""Exact checks for the coefficient envelope and its native prime tests."""

from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb, factorial, gcd, lcm, prod
import json


def envelope(rows, order):
    return tuple(
        max(order * a + d * (k if sign > 0 else order - k) for a, d, sign in rows)
        for k in range(order + 1)
    )


def window(rows, order):
    height = envelope(rows, order)
    floor = min(height)
    ceiling = min(max(height[k], height[k + 2]) for k in range(order - 1))
    return height, floor, ceiling


def row_corner(rows, order, sign):
    base = max(a for a, _, _ in rows)
    opposite = [2 * d - order * (base + 1 - a) for a, d, s in rows if s != sign]
    return (
        (base, 1, sign) in rows
        and all(a + d <= base + 1 for a, d, s in rows if s == sign)
        and max(opposite, default=-2) in (-1, 0)
    )


def pair_corner(rows, lower, sign):
    base = max(a for a, _, _ in rows)
    return (
        lower % 2 == 0
        and (base, 1, sign) in rows
        and (base, lower // 2, -sign) in rows
        and all(a + d <= base + 1 for a, d, s in rows if s == sign)
        and all(2 * d <= lower * (base + 1 - a) for a, d, s in rows if s != sign)
    )


def literal_counts(rows, order, mark):
    return tuple(
        sum(
            order * a < mark <= order * a + d * (k if s > 0 else order - k)
            for a, d, s in rows
        )
        for k in range(order + 1)
    )


def coefficients(rows, order):
    return tuple(
        comb(order, k) * prod(
            prod(range(order * a + 1, order * a + d * (k if s > 0 else order-k) + 1))
            for a, d, s in rows
        )
        for k in range(order + 1)
    )


def prime_divisors(value):
    result = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            result.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        result.append(value)
    return result


def valuation(value, prime):
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def slope_widths(ordinates):
    hull = []
    for k, y in enumerate(ordinates):
        while len(hull) >= 2:
            left, right = hull[-2:]
            slope = Fraction(right[1] - left[1], right[0] - left[0])
            if slope < Fraction(y - right[1], k - right[0]):
                break
            hull.pop()
        hull.append((k, y))
    return {
        Fraction(right[1] - left[1], right[0] - left[0]): right[0] - left[0]
        for left, right in zip(hull, hull[1:])
    }


def replay_bound(marked, unmarked, prime):
    a = slope_widths([valuation(c, prime) for c in marked])
    b = slope_widths([valuation(c, prime) for c in unmarked])
    assert sum(min(a[slope], b[slope]) for slope in a.keys() & b.keys()) <= 1


def prime_tests(rows, marked_order, other_order, mark, counts):
    literal = literal_counts(rows, marked_order, mark)
    missing = tuple(k for k, count in enumerate(literal) if count == 0)
    outside = tuple(k for k, count in enumerate(literal) if count)
    marked = coefficients(rows, marked_order)
    other = coefficients(rows, other_order)
    greatest_distance = 0
    numerator_gcd = 0
    lcms = [
        lcm(
            *(
                order * a + h
                for order in (marked_order, other_order)
                for h in range(1, order * d + 1)
                if order != marked_order or order * a + h != mark
            )
        )
        for a, d, _ in rows
    ]
    for k in outside:
        q = min(missing, key=lambda entry: abs(k - entry))
        distance = abs(k - q)
        greatest_distance = max(greatest_distance, distance)
        for length in range(1, other_order + 1):
            if k > q:
                suffix = Fraction(other[other_order - length], other[-1])
            else:
                suffix = Fraction(other[length], other[0])
            ratio = Fraction(marked[k], marked[q])**length * suffix**distance
            numerator_gcd = gcd(numerator_gcd, ratio.numerator)

    divisor_conditions_fail = False
    for prime in prime_divisors(mark):
        depth = valuation(mark, prime)
        caps = [valuation(integer, prime) for integer in lcms]
        total = sum(caps)
        b_marked = valuation(factorial(marked_order), prime)
        b_other = valuation(factorial(other_order), prime)
        covered = True
        for k in outside:
            excess = sum(
                max(0, depth - cap)
                for (a, d, s), cap in zip(rows, caps)
                if marked_order * a < mark <= marked_order * a + d * (k if s > 0 else marked_order - k)
            )
            distance = min(abs(k - q) for q in missing)
            covered &= excess > (1 + distance) * total + b_marked + distance * b_other
        if covered:
            replay_bound(marked, other, prime)
            counts["layer_certificates"] += 1
        if min(literal[k] for k in outside) * depth > (
            (greatest_distance + 2) * total + b_marked + greatest_distance * b_other
        ):
            assert covered
            counts["divisibility_certificates"] += 1
            divisor_conditions_fail = True

    if numerator_gcd > 1:
        for prime in prime_divisors(numerator_gcd):
            replay_bound(marked, other, prime)
        counts["gcd_certificates"] += 1
        if not divisor_conditions_fail:
            counts["gcd_improves_divisibility"] += 1


def run():
    counts = Counter()
    types = tuple(product(range(4), range(1, 4), (-1, 1)))
    for size in (1, 2, 3):
        for index, rows in enumerate(combinations_with_replacement(types, size)):
            greatest_gap = max(d for _, d, _ in rows)
            smallest_gap = min(d for _, d, _ in rows)
            top = max(a + d for a, d, _ in rows)
            legal = {}
            for order in range(3, 10):
                height, floor, ceiling = window(rows, order)
                assert smallest_gap <= ceiling - floor <= 2 * greatest_gap
                allowed = []
                for mark in range(floor + 1, order * top + 3):
                    literal = literal_counts(rows, order, mark)
                    missing = tuple(k for k, count in enumerate(literal) if count == 0)
                    assert missing == tuple(k for k, value in enumerate(height) if value < mark)
                    assert (missing[-1] - missing[0] <= 1) == (mark <= ceiling)
                    if mark <= ceiling and mark % order:
                        allowed.append(mark)
                legal[order] = allowed
                assert (not allowed) == any(row_corner(rows, order, sign) for sign in (-1, 1))
                counts["order_windows"] += 1
            for lower in range(3, 9):
                upper = lower + 1
                failure = not legal[lower] and not legal[upper]
                assert failure == any(pair_corner(rows, lower, sign) for sign in (-1, 1))
                counts["consecutive_classifications"] += 1
                counts["simultaneous_corners"] += failure
                if index % 31 == 0:
                    for marked_order, other_order in ((lower, upper), (upper, lower)):
                        for mark in legal[marked_order]:
                            prime_tests(rows, marked_order, other_order, mark, counts)

    for even in range(4, 17, 2):
        gap = even // 2
        for copies in (1, 3):
            a = 20
            b = a + gap + 1
            rows = ((a, 1, 1), (a, gap, -1), (b, 1, -1), (b, gap, 1)) * copies
            assert sum(g * s for _, g, s in rows) == 0
            for order in (even, even + 1):
                for mark in range(order * a + 1, order * (b + gap) + 1):
                    if mark % order:
                        array = literal_counts(rows, order, mark)
                        minimum = [k for k, count in enumerate(array) if count == min(array)]
                        assert minimum[-1] - minimum[0] >= 2
                        counts["sharp_odd_literal_failures"] += 1

    for gap in (6, 12, 25, 40, 80):
        baseline = 8 * gap**2
        rows = ((baseline, gap, 1), (baseline, gap, -1), (baseline - 2 * gap, 1, 1))
        for order in range(3, 10):
            _, low, high = window(rows, order)
            width = high - low
            marks = tuple(mark for mark in range(low + 1, high + 1) if mark % order == 1)
            if len(marks) < 2:
                continue
            assert width >= gap and min(d for _, d, _ in rows) == 1
            primes = {p for mark in marks for p in prime_divisors(mark)}
            for prime in primes:
                vals = [valuation(mark, prime) for mark in marks]
                if prime > width:
                    assert sum(value > 0 for value in vals) == 1
                else:
                    assert sum(vals) - max(vals) <= valuation(factorial(len(marks)-1), prime)
                counts["mark_progression_prime_checks"] += 1

    assert counts["simultaneous_corners"] > 0
    assert counts["layer_certificates"] > 0
    assert counts["gcd_improves_divisibility"] > 0
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
