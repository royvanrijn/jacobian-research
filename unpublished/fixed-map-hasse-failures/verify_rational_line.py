#!/usr/bin/env python3
"""Exact finite audits for the rational target line (Python standard library).

Run from the repository root:
    python3 unpublished/fixed-map-hasse-failures/verify_rational_line.py

The checks compare independent finite descriptions of projective content,
height, local cube obstructions and bounded-height enumeration.  They also
check rational polygon areas and Euler-factor algebra.  No finite calculation
here proves the local classification for all primes or the asymptotic: those
are the written arguments in sections/03-local-global.tex and
sections/04-rational-counting.tex.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
from fractions import Fraction as Q
from functools import lru_cache, reduce
from hashlib import sha256
import json
from math import gcd, lcm
from pathlib import Path
import platform


@lru_cache(maxsize=None)
def factor(n: int) -> tuple[tuple[int, int], ...]:
    assert n >= 1
    factors = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            factors.append((p, e))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def strip(n: int, p: int) -> tuple[int, int]:
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e, n


def is_cube(n: int) -> bool:
    return all(e % 3 == 0 for _, e in factor(abs(n)))


def support_allowed(n: int) -> bool:
    return all(p % 3 == 1 or e % 3 == 0 for p, e in factor(n))


def local_criterion(m: int, n: int) -> bool:
    """The claimed criterion, including cubes, for nonzero reduced m/n."""
    assert m and n > 0 and gcd(m, n) == 1
    if not support_allowed(abs(m)) or not support_allowed(n):
        return False
    em, um = strip(m, 3)
    en, un = strip(n, 3)
    return (em - en) % 3 == 0 and (um * pow(un, -1, 9)) % 9 in (1, 8)


@lru_cache(maxsize=None)
def unit_cube_residues(modulus: int) -> frozenset[int]:
    return frozenset(pow(x, 3, modulus) for x in range(modulus)
                     if gcd(x, modulus) == 1)


def modular_local_test(m: int, n: int, p: int) -> bool:
    """Independent finite unit-residue test after extracting valuations.

    At 3 use modulus 9; elsewhere use p**2.  Testing these finite residues
    is evidence for, rather than a replacement for, the Hensel arguments.
    """
    quadratic_modulus = 9 if p == 3 else p
    if any((x * x + x + 1) % quadratic_modulus == 0
           for x in range(quadratic_modulus)):
        return True
    em, um = strip(m, p)
    en, un = strip(n, p)
    modulus = 9 if p == 3 else p * p
    return ((em - en) % 3 == 0
            and um * pow(un, -1, modulus) % modulus
            in unit_cube_residues(modulus))


def raw_coordinates(m: int, n: int) -> tuple[int, ...]:
    return 9 * n, -9 * n, 32 * m, 24 * m + 3 * n


def direct_content(m: int, n: int) -> int:
    return reduce(gcd, raw_coordinates(m, n))


def content_formula(m: int, n: int) -> int:
    return gcd(n, 8) * gcd(m, 3)


def height_from_rational_coordinates(m: int, n: int) -> int:
    """Clear denominators of the actual target, without using raw formula."""
    a = Q(m, n)
    coordinates = (Q(1), Q(-1), 32 * a / 9, (8 * a + 1) / 3)
    denominator = lcm(*(c.denominator for c in coordinates))
    integers = tuple(int(c * denominator) for c in coordinates)
    return max(map(abs, integers)) // reduce(gcd, integers)


def height_formula(m: int, n: int) -> int:
    return max(map(abs, raw_coordinates(m, n))) // content_formula(m, n)


def polygon_area(vertices: tuple[tuple[Q, Q], ...]) -> Q:
    return abs(sum(x * v - y * u
                   for (x, y), (u, v) in zip(vertices, vertices[1:] + vertices[:1]))) / 2


def finite_pair_checks(bound: int = 128) -> dict:
    primes = tuple(p for p in range(2, bound + 1) if factor(p) == ((p, 1),))
    coprime_pairs = 0
    local_pairs = 0
    for n in range(1, bound + 1):
        for m in range(-bound, bound + 1):
            if not m or gcd(m, n) != 1:
                continue
            coprime_pairs += 1
            assert direct_content(m, n) == content_formula(m, n)
            assert height_from_rational_coordinates(m, n) == height_formula(m, n)
            claimed = local_criterion(m, n)
            modular = all(modular_local_test(m, n, p) for p in primes)
            assert claimed == modular, (m, n, claimed, modular)
            if claimed:
                local_pairs += 1
                assert content_formula(m, n) == (
                    (8 if n % 2 == 0 else 1) * (3 if m % 3 == 0 else 1))
    return {"numerator_absolute_bound": bound, "denominator_bound": bound,
            "coprime_pairs": coprime_pairs, "locally_admissible_including_cubes": local_pairs,
            "tested_primes": list(primes)}


def enumerate_optimized(maximum_height: int) -> list[int]:
    """Complete bounded count using proved d<=24 and support prefiltering."""
    # H<=B implies 32|m|<=24B and 9n<=24B, including both signs.
    max_m = 3 * maximum_height // 4
    max_n = 8 * maximum_height // 3
    numerators = [m for m in range(1, max_m + 1) if support_allowed(m)]
    denominators = [n for n in range(1, max_n + 1) if support_allowed(n)]
    heights = []
    for n in denominators:
        for m in numerators:
            # Even the largest possible content cannot meet the height bound.
            if 32 * m > 24 * maximum_height:
                break
            if gcd(m, n) != 1 or not local_criterion(m, n):
                continue
            if is_cube(m) and is_cube(n):
                continue
            for signed_m in (m, -m):
                h = height_from_rational_coordinates(signed_m, n)
                assert h == height_formula(signed_m, n)
                if h <= maximum_height:
                    heights.append(h)
    return sorted(heights)


def enumerate_naive(maximum_height: int) -> list[int]:
    """Independent tiny-domain oracle, with no support prefilter."""
    heights = []
    for n in range(1, 8 * maximum_height // 3 + 1):
        for m in range(-3 * maximum_height // 4, 3 * maximum_height // 4 + 1):
            if not m or gcd(m, n) != 1 or (is_cube(m) and is_cube(n)):
                continue
            h = height_from_rational_coordinates(m, n)
            if h > maximum_height:
                continue
            primes = {3} | {p for p, _ in factor(abs(m))} | {p for p, _ in factor(n)}
            if all(modular_local_test(m, n, p) for p in primes):
                heights.append(h)
    return sorted(heights)


def algebra_checks() -> dict:
    positive = polygon_area(((Q(0), Q(0)), (Q(1, 32), Q(0)),
                             (Q(1, 32), Q(1, 12)),
                             (Q(1, 36), Q(1, 9)), (Q(0), Q(1, 9))))
    negative = polygon_area(((Q(0), Q(0)), (Q(1, 32), Q(0)),
                             (Q(1, 32), Q(1, 9)), (Q(0), Q(1, 9))))
    assert positive == Q(71, 20736)
    assert negative == Q(72, 20736)
    # The sums are geometric series with exponents 3j; the coefficients
    # are the squared content on the numerator and denominator branches.
    weight_2 = 1 + (1 + 8**2) * Q(1, 2**3 - 1)
    weight_3 = 1 + (3**2 + 1) * Q(1, 3**3 - 1)
    assert weight_2 == Q(72, 7) and weight_3 == Q(18, 13)
    # Relative to a baseline product already including the inert prime 2,
    # remove its ordinary coprime factor (1+1/8)/(1-1/8) = 9/7.
    constant_rational_factor = ((positive + negative) * Q(1, 3)
                                * weight_2 * weight_3 * Q(7, 9) * Q(1, 4))
    assert constant_rational_factor == Q(11, 1728)
    for p in (2, 5, 11, 17, 23):
        x = Q(1, p)
        # Inert baseline: univariate unit/cube factors squared times
        # coprimality correction, then zeta(2) normalization.
        assert (1 - x*x) * (1 + x**3) / (1 - x**3) == (
            (1 - x*x) / (1 - x**3)**2 * (1 - x**6))
    return {"positive_area": str(positive), "negative_area": str(negative),
            "signed_area": str(positive + negative),
            "content_weight_at_2": str(weight_2), "content_weight_at_3": str(weight_3),
            "constant": "11/(1728*sqrt(3)*zeta(2)) * product_{p=2 mod3}(1+p^-3)/(1-p^-3)"}


def example_checks() -> list[dict]:
    examples = [(152, 1, True), (19, 37, True), (19, 8, True),
                (513, 1, True), (1, 513, True), (513, 8, True),
                (-19, 37, True), (1, 1, False),
                (8, 1, False), (27, 1, False), (2, 1, False),
                (7, 1, False), (19, 2, False), (19, 3, False), (57, 1, False)]
    result = []
    for m, n, expected_hasse in examples:
        local = local_criterion(m, n)
        hasse = local and not (is_cube(m) and is_cube(n))
        assert hasse == expected_hasse, (m, n, local, hasse)
        result.append({"a": str(Q(m, n)), "everywhere_local": local,
                       "hasse_failure": hasse,
                       "content": direct_content(m, n),
                       "height": height_from_rational_coordinates(m, n)})
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-height", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert args.maximum_height >= 80
    result = {"scope": "Internal exact finite checks, not a proof of the unrestricted classification or asymptotic",
              "python": platform.python_version(),
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "pair_checks": finite_pair_checks(), "algebra": algebra_checks(),
              "examples": example_checks()}
    tiny_naive = enumerate_naive(80)
    tiny_optimized = enumerate_optimized(80)
    assert tiny_naive == tiny_optimized
    heights = enumerate_optimized(args.maximum_height)
    bounds = sorted({80, 100, 300, 1000, 3000, 10000, args.maximum_height})
    bounds = [b for b in bounds if b <= args.maximum_height]
    result["height_counts"] = {str(b): bisect_right(heights, b) for b in bounds}
    result["enumeration"] = {"maximum_height": args.maximum_height,
                             "numerator_absolute_bound": 3 * args.maximum_height // 4,
                             "denominator_bound": 8 * args.maximum_height // 3,
                             "naive_and_optimized_agree_through_height": 80}
    result["outcome"] = "pass"
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
