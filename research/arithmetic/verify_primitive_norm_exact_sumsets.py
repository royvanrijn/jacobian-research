#!/usr/bin/env python3
"""Exact bounded primitive binary-norm sumset enumeration; no extrapolation.

Python standard library only. Every marked value is generated from an integral
pair with gcd 1. The short-summand bitset sweep supplies positive witnesses;
any remaining locally admissible target is checked against *all* summands.
"""
import argparse
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path


def norm_values(discriminant, bound):
    flags = bytearray(bound + 1)
    if discriminant == -8:
        for y in range(isqrt(bound // 2) + 1):
            for x in range(isqrt(bound - 2*y*y) + 1):
                if gcd(x, y) == 1:
                    flags[x*x + 2*y*y] = 1
    else:
        d = -discriminant
        c = (d + 1) // 4
        for y in range(isqrt(4*bound // d) + 1):
            # Reflection (x,y)->(-x-y,y) permits 2*x+y >= 0.
            upper = (isqrt(4*bound - d*y*y) - y) // 2
            for x in range(-(y//2), upper + 1):
                if gcd(x, y) == 1:
                    flags[x*x + x*y + c*y*y] = 1
    flags[0] = 0
    return flags


def admissible(discriminant, n):
    if discriminant == -8:
        return n % 16 != 0
    d = -discriminant
    return n % (d*d) != 0 and (d == 7 or n % 2 == 0)


def scan(discriminant, bound, sweep_bound):
    flags = norm_values(discriminant, bound)
    values = [n for n in range(1, bound + 1) if flags[n]]
    packed = bytearray((bound + 8) // 8)
    for n in values:
        packed[n//8] |= 1 << (n % 8)
    bits = int.from_bytes(packed, 'little')
    covered = 0
    for a in values:
        if a > sweep_bound:
            break
        covered |= bits << a
    covered &= (1 << (bound + 1)) - 1
    packed_covered = covered.to_bytes((bound + 8)//8, 'little')
    exceptions = []
    remaining = 0
    forbidden_hits = []
    for n in range(2, bound + 1):
        marked = bool(packed_covered[n//8] & (1 << (n % 8)))
        local = admissible(discriminant, n)
        if marked and not local:
            forbidden_hits.append(n)
        if local and not marked:
            remaining += 1
            for a in values:
                if 2*a > n:
                    exceptions.append(n)
                    break
                if flags[n-a]:
                    break
            else:
                exceptions.append(n)
    assert not forbidden_hits, forbidden_hits[:10]
    return dict(discriminant=discriminant, bound=bound,
                primitive_value_count=len(values), sweep_bound=sweep_bound,
                exhaustive_residual_targets=remaining,
                locally_admissible_global_exceptions=exceptions,
                scope='Exact only for 2 <= N <= bound; no unrestricted claim.')


def in_split_semigroup(n):
    """Independent trial-factor membership in C, with primes 1,3 mod 8."""
    if n < 1:
        return False
    p = 2
    while p*p <= n:
        if n % p == 0:
            if p % 8 not in (1, 3):
                return False
            while n % p == 0:
                n //= p
        p += 1
    return n == 1 or n % 8 in (1, 3)


def structural_checks():
    # Gcd-coordinate enumeration is independent of the trial-factor test.
    flags = norm_values(-8, 10000)
    for n in range(1, 10001):
        support = in_split_semigroup(n) or (
            n % 2 == 0 and in_split_semigroup(n//2))
        assert bool(flags[n]) == support, n
    for t in range(1, 101):
        in_ternary = any(flags[t*t-x*x] for x in range(t))
        in_c_sum = any(in_split_semigroup(a) and
                       in_split_semigroup(2*t-a)
                       for a in range(1, 2*t))
        assert in_ternary == in_c_sum, t
        if t % 2:
            in_primitive_sum = any(flags[a] and flags[2*t-a]
                                   for a in range(1, 2*t))
            assert in_primitive_sum == in_c_sum, t
    return dict(split_support_bound=10000, square_factorization_targets=100,
                scope='Finite cross-checks of the proved splitting-support '
                'description and ternary square-factorization equivalence.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bound', type=int, default=1000000)
    ap.add_argument('--sweep-bound', type=int, default=5000)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    assert args.bound >= 2 and args.sweep_bound >= 1
    rows = [scan(d, args.bound, args.sweep_bound)
            for d in (-8, -11, -7, -19, -43, -67, -163)]
    # Independent hand checks exercise nonrepresentation, boundaries and gcd.
    assert norm_values(-8, 16) == bytearray(
        int(n in {1, 2, 3, 6, 9, 11}) for n in range(17))
    assert norm_values(-11, 16) == bytearray(
        int(n in {1, 3, 5, 9, 11, 15}) for n in range(17))
    result = dict(algorithm='Integral-pair enumeration, exact integer bitset '
                  'witness sweep, exhaustive residual scan.',
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scans=rows, structural_checks=structural_checks())
    out = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(out)
    print(out, end='')


if __name__ == '__main__':
    main()
