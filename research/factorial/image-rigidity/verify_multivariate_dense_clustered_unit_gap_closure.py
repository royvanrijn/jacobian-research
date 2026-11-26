#!/usr/bin/env python3
"""Exact checks for the gap-one short-cluster prime certificate."""

from collections import Counter
from fractions import Fraction
from math import gcd
import json


CUTOFF = 12
ORIENTATIONS = ((2, 0), (1, 3), (2, 2), (0, 3), (5, 5), (3, 1))


def factorization(value):
    factors = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            value //= divisor
        divisor += 1
    if value > 1:
        factors[value] = factors.get(value, 0) + 1
    return factors


def factorial_valuation(value, prime):
    result = 0
    while value:
        value //= prime
        result += value
    return result


def slopes(ordinates):
    hull = []
    for index, value in enumerate(ordinates):
        while len(hull) >= 2:
            a, b = hull[-2:]
            previous = Fraction(b[1] - a[1], b[0] - a[0])
            joining = Fraction(value - b[1], index - b[0])
            if previous < joining:
                break
            hull.pop()
        hull.append((index, value))
    return {
        Fraction(b[1] - a[1], b[0] - a[0])
        for a, b in zip(hull, hull[1:])
    }


def ordinates(cluster, order, prime):
    output = []
    for index in range(order + 1):
        result = 0
        for base, (plus, minus) in cluster.items():
            start = order * base
            baseline = factorial_valuation(start, prime)
            forward = factorial_valuation(start + index, prime) - baseline
            backward = factorial_valuation(start + order - index, prime) - baseline
            result += plus * forward + minus * backward
        output.append(result)
    return output


def mark(base, counts, even):
    plus, minus = counts
    offset = 1 if plus != minus else even // 2 + 1
    raw = even * base + offset
    step = even // gcd(even, offset)
    intercept = offset // gcd(even, offset)
    assert gcd(step, intercept) == 1
    return raw, offset, raw // gcd(even, offset)


def replay_certificate(cluster, even, odd, base, prime, depth):
    plus, minus = cluster[base]
    raw, offset, _ = mark(base, cluster[base], even)
    hits = [
        (order, height, h)
        for order in (even, odd)
        for height in cluster
        for h in range(1, order + 1)
        if (order * height + h) % prime == 0
    ]
    assert hits == [(even, base, offset)]
    assert raw % prime**depth == 0 and raw % prime**(depth + 1) != 0
    even_array = ordinates(cluster, even, prime)
    odd_array = ordinates(cluster, odd, prime)
    if offset == 1:
        expected = [minus * depth] + [(plus + minus) * depth] * (even - 1) + [plus * depth]
        expected_slopes = {Fraction((plus - minus) * depth, even)}
    else:
        half = even // 2
        expected = [minus * depth] * half + [0] + [plus * depth] * half
        expected_slopes = {Fraction(-minus * depth, half), Fraction(plus * depth, half)}
    assert even_array == expected
    assert odd_array == [0] * (odd + 1)
    assert slopes(even_array) == expected_slopes
    assert slopes(odd_array) == {0}
    assert 0 not in expected_slopes


def first_prime_at_least(lower, residue, modulus):
    return next(
        candidate
        for candidate in range(lower, lower + 1000)
        if candidate % modulus == residue and factorization(candidate) == {candidate: 1}
    )


def run():
    counts = Counter()
    for even in range(4, CUTOFF, 2):
        for odd in (even - 1, even + 1):
            for width in (1, 3, 7):
                threshold = 2 * CUTOFF**2 * (width + 1)
                for base in range(1800, 1870, 3):
                    cluster = {
                        base + index: ORIENTATIONS[(base + index) % len(ORIENTATIONS)]
                        for index in range(width + 1)
                        if (base + index) % 5
                    }
                    if not cluster:
                        continue
                    counts["order_cluster_pairs"] += 1
                    branch_hits = Counter()
                    prime_to_heights = {}
                    for height, kinds in cluster.items():
                        raw, offset, normalized = mark(height, kinds, even)
                        if gcd(even, offset) > 1:
                            counts["nontrivial_gcd_normalizations"] += 1
                        for prime, exponent in factorization(normalized).items():
                            for depth in range(1, exponent + 1):
                                branch_hits[offset, prime**depth] += 1
                            if prime <= threshold:
                                continue
                            counts["large_prime_certificates"] += 1
                            assert prime not in prime_to_heights
                            prime_to_heights[prime] = height
                            replay_certificate(cluster, even, odd, height, prime, exponent)
                    for (_, prime_power), number in branch_hits.items():
                        assert (number - 1) * prime_power <= width
                        counts["residue_counts"] += 1
                    if prime_to_heights:
                        # A retained ramp starting at the marked height can hit its factor.
                        prime, height = next(iter(prime_to_heights.items()))
                        fixed_ramp_length = even * even
                        forbidden = {
                            test
                            for test in prime_to_heights
                            if any((even * height + h) % test == 0
                                   for h in range(1, fixed_ramp_length + 1))
                        }
                        blocked = {prime_to_heights[test] for test in forbidden}
                        assert prime in forbidden and height in blocked
                        assert len(blocked) <= len(forbidden)
                        counts["retained_overlap_fixtures"] += 1

    width = 3
    threshold = 2 * CUTOFF**2 * (width + 1)
    for even, depth, kinds in ((4, 2, (7, 2)), (6, 3, (4, 4))):
        offset = 1 if kinds[0] != kinds[1] else even // 2 + 1
        common = gcd(even, offset)
        prime = first_prime_at_least(threshold + 1, 2, 3)
        raw = common * prime**depth
        assert (raw - offset) % even == 0
        height = (raw - offset) // even
        cluster = {height - 1: (5, 2), height: kinds, height + 2: (3, 3)}
        for odd in (even - 1, even + 1):
            replay_certificate(cluster, even, odd, height, prime, depth)
            counts["deep_prime_certificates"] += 1

    assert counts["large_prime_certificates"] > 0
    assert counts["nontrivial_gcd_normalizations"] > 0
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
