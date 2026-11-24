#!/usr/bin/env python3
"""Exact constants and a finite CRT control for prime-discriminant core sumsets."""

import hashlib
import json
import platform
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, prod
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run this certificate checker with Python assertions active.")

from verify_primitive_norm_minus11 import prime_flags


def divisor_bound(h, dyadic_exponent):
    limit = 2**h
    assert (2*h+2)**h < 2*(2*h+1)**h
    flags = prime_flags(limit-1)
    factors = [
        max(Fraction((a+1)**h, q**a) for a in range(2*h+1))
        for q, prime in enumerate(flags) if prime
    ]
    numerator = prod(f.numerator for f in factors)
    denominator = prod(f.denominator for f in factors)
    assert numerator < denominator << (h*dyadic_exponent)
    assert numerator >= denominator << (h*(dyadic_exponent-1))
    # 2^k exceeds (256*2^dyadic_exponent)^(4h/(h-4)), or equals it.
    power = Fraction(4*h*(dyadic_exponent+8), h-4)
    k = -(-power.numerator//power.denominator)
    return {
        "h": h, "checked_primes": len(factors),
        "divisor_bound_dyadic_exponent": dyadic_exponent,
        "cutoff_dyadic_exponent": k,
        "discriminant_exponent": str(Fraction(3*h,h-4)),
    }


def example_gap_set(odd_primes, bits, bound):
    allowed = [q for q, selected in zip(odd_primes, bits, strict=True) if selected]
    core = [1]
    for n in range(3, bound, 2):
        remaining = n
        for q in allowed:
            while remaining % q == 0:
                remaining //= q
        if remaining == 1:
            core.append(n)
    sums = {a+b for a in core for b in core if a+b <= bound}
    return tuple(n for n in range(4, bound+1, 2) if n not in sums)


def crt_control():
    bound = 12
    odd_primes = [q for q, prime in enumerate(prime_flags(bound-1)) if q > 2 and prime]
    odd_modulus = prod(odd_primes)
    modulus = 8*odd_modulus
    classes = Counter(
        tuple(pow(a, (q-1)//2, q) == 1 for q in odd_primes)
        for a in range(3, modulus, 8) if gcd(a, odd_modulus) == 1
    )
    all_runs = tuple(range(4, bound+1, 2))
    patterns = Counter()
    bits_list = tuple(product((False, True), repeat=len(odd_primes)))
    for bits in bits_list:
        gaps = example_gap_set(odd_primes, bits, bound)
        assert classes[bits] == prod((q-1)//2 for q in odd_primes)
        assert (gaps == all_runs) == (not any(bits))
        patterns[gaps] += 1
    assert sum(patterns.values()) == 2**len(odd_primes)
    assert Fraction(sum(weight for gaps, weight in patterns.items() if 4 in gaps),
                    sum(patterns.values())) == Fraction(1,2)
    assert Fraction(sum(weight for gaps, weight in patterns.items() if 6 in gaps),
                    sum(patterns.values())) == Fraction(1,4)
    assert all(4 in gaps for gaps in patterns if 6 in gaps)
    return {
        "target_bound": bound, "modulus": modulus, "odd_primes": odd_primes,
        "classes_per_splitting_pattern": classes[bits_list[0]],
        "initial_run_probability": str(Fraction(1, len(bits_list))),
        "gap_distribution_numerators": [
            {"gaps": gaps, "weight": weight} for gaps, weight in sorted(patterns.items())
        ],
        "gap_distribution_denominator": len(bits_list),
    }


def main():
    assert Fraction(56,81) > Fraction(11,16)
    assert Fraction(37,6*256) < Fraction(1,24)
    source = Path(__file__).resolve().parent
    names = ("verify_prime_discriminant_core_uniform.py", "verify_primitive_norm_minus11.py",
             "../foundations/extended-geometry/sources/soprounov_pnt_arithmetic_progressions.pdf")
    print(json.dumps({
        "python": platform.python_version(),
        "inputs_sha256": {
            name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in names
        },
        "uniform_divisor_bounds": [divisor_bound(h,b) for h,b in ((8,16),(10,36),(16,685))],
        "crt_control": crt_control(),
    }, indent=2))


if __name__ == "__main__":
    main()
