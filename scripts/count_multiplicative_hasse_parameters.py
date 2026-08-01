#!/usr/bin/env python3
"""Enumerate multiplicative Hasse parameters on the fixed target line.

The clean family consists of noncubes a > 1 all of whose prime divisors are
1 modulo 9.  The broader admissible family allows prime divisors 1 modulo 3
but retains a = 1 modulo 9.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
from pathlib import Path


def prime_sieve(bound: int) -> bytearray:
    sieve = bytearray(b"\x01") * (bound + 1)
    if bound >= 0:
        sieve[0] = 0
    if bound >= 1:
        sieve[1] = 0
    for prime in range(2, math.isqrt(bound) + 1):
        if sieve[prime]:
            start = prime * prime
            count = (bound - start) // prime + 1
            sieve[start : bound + 1 : prime] = b"\x00" * count
    return sieve


def support_sieve(bound: int, primes: list[int], modulus: int) -> bytearray:
    supported = bytearray(b"\x01") * (bound + 1)
    supported[0] = 0
    for prime in primes:
        if prime % modulus != 1:
            count = bound // prime
            supported[prime : bound + 1 : prime] = b"\x00" * count
    return supported


def cube_values(bound: int, supported: bytearray) -> set[int]:
    root_bound = round(bound ** (1 / 3))
    while (root_bound + 1) ** 3 <= bound:
        root_bound += 1
    while root_bound**3 > bound:
        root_bound -= 1
    return {
        root**3
        for root in range(2, root_bound + 1)
        if supported[root]
    }


def digest(values: list[int]) -> str:
    payload = ",".join(map(str, values)).encode("ascii")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def semiprime_count(bound: int, allowed_primes: list[int]) -> int:
    total = 0
    for index, prime in enumerate(allowed_primes):
        if prime * prime > bound:
            break
        total += bisect.bisect_right(allowed_primes, bound // prime) - index
    return total


def scaled(count: int, bound: int, exponent: float) -> float:
    return count * math.log(bound) ** exponent / bound


def enumerate_parameters(bound: int) -> dict[str, object]:
    prime_flags = prime_sieve(bound)
    primes = [number for number in range(2, bound + 1) if prime_flags[number]]
    clean_support = support_sieve(bound, primes, 9)
    broad_support = support_sieve(bound, primes, 3)
    clean_cubes = cube_values(bound, clean_support)
    broad_cubes = cube_values(bound, broad_support)

    clean = [
        number
        for number in range(2, bound + 1)
        if clean_support[number] and number not in clean_cubes
    ]
    broad = [
        number
        for number in range(2, bound + 1)
        if (
            broad_support[number]
            and number % 9 == 1
            and number not in broad_cubes
        )
    ]

    assert set(clean).issubset(broad)
    for parameter in broad:
        coordinates = (9, -9, 32 * parameter, 24 * parameter + 3)
        assert math.gcd(*coordinates) == 1
        assert max(map(abs, coordinates)) == 32 * parameter

    allowed_primes = [prime for prime in primes if prime % 9 == 1]
    truncated_product = math.prod(
        1 - prime**-2 for prime in primes if prime % 3 == 2
    )
    l_one_chi_minus_three = math.pi / (3 * math.sqrt(3))
    full_constant_truncated = math.sqrt(
        (2 / 3) * l_one_chi_minus_three * truncated_product
    ) / (3 * math.sqrt(math.pi))
    checkpoints = sorted(
        {
            min(bound, checkpoint)
            for checkpoint in (1_000, 10_000, 100_000, 1_000_000, bound)
        }
    )
    rows = []
    for checkpoint in checkpoints:
        clean_count = bisect.bisect_right(clean, checkpoint)
        broad_count = bisect.bisect_right(broad, checkpoint)
        prime_count = bisect.bisect_right(allowed_primes, checkpoint)
        two_prime_count = semiprime_count(
            checkpoint,
            allowed_primes[:prime_count],
        )
        rows.append(
            {
                "bound": checkpoint,
                "clean_count": clean_count,
                "broad_count": broad_count,
                "prime_count": prime_count,
                "two_prime_product_count": two_prime_count,
                "clean_scaled_log_power_5_over_6": scaled(
                    clean_count, checkpoint, 5 / 6
                ),
                "broad_scaled_log_power_1_over_2": scaled(
                    broad_count, checkpoint, 1 / 2
                ),
                "two_prime_to_prime_ratio": (
                    two_prime_count / prime_count if prime_count else 0.0
                ),
            }
        )

    return {
        "generator": "scripts/count_multiplicative_hasse_parameters.py",
        # Keep the pinned mathematical artifact independent of the interpreter
        # patch release used to reproduce it.
        "software": "Python standard library",
        "bound": bound,
        "clean_definition": (
            "a>1, every prime divisor is 1 mod 9, and a is not a cube"
        ),
        "broad_definition": (
            "a>1, a=1 mod 9, every prime divisor is 1 mod 3, "
            "and a is not a cube"
        ),
        "target": "(-1,32*a/9,(8*a+1)/3)",
        "primitive_coordinates": "[9:-9:32*a:24*a+3]",
        "height": "32*a",
        "clean_count": len(clean),
        "broad_count": len(broad),
        "broad_asymptotic_constant_euler_product_cutoff": bound,
        "broad_asymptotic_constant_truncated": full_constant_truncated,
        "clean_digest": digest(clean),
        "broad_digest": digest(broad),
        "clean_first_25": clean[:25],
        "clean_last_25": clean[-25:],
        "broad_first_25": broad[:25],
        "broad_last_25": broad[-25:],
        "checkpoints": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.bound < 19:
        raise SystemExit("--bound must be at least 19")

    result = enumerate_parameters(args.bound)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
