#!/usr/bin/env python3
"""Verify the exponent divisor certificate for an Erdős #364 power-of-two center."""

import argparse
import csv
import hashlib
from math import gcd, isqrt, prod
from pathlib import Path
from sys import set_int_max_str_digits


SWEEP_SHA256 = "def14e2946677f0b942ee093ec1039489a73ce0a71abb50f8c8cd3b9b69344f5"
# Each row gives a forced index, a simple prime divisor at twice that index,
# a Lucas primality witness, and the complete factorization of prime - 1.
ADDITIONAL = (
    (43, 2932031007403, 3, {2: 1, 3: 1, 7: 2, 43: 1, 127: 1, 337: 1, 5419: 1}),
    (
        127, 170141183460469231731687303715884105727, 43,
        {2: 1, 3: 3, 7: 2, 19: 1, 43: 1, 73: 1, 127: 1, 337: 1, 5419: 1,
         92737: 1, 649657: 1, 77158673929: 1},
    ),
    (
        127, 56713727820156410577229101238628035243, 3,
        {2: 1, 3: 2, 7: 2, 19: 1, 43: 1, 73: 1, 127: 1, 337: 1, 5419: 1,
         92737: 1, 649657: 1, 77158673929: 1},
    ),
    (
        129, 11053036065049294753459639, 13,
        {2: 1, 3: 2, 7: 1, 11: 1, 43: 1, 53: 1, 97: 1, 127: 1, 337: 1,
         5419: 1, 155542661: 1},
    ),
    (129, 1591582393, 5, {2: 3, 3: 3, 43: 1, 349: 1, 491: 1}),
    (129, 15686603697451, 3, {2: 1, 3: 1, 5: 2, 7: 1, 43: 1, 347433083: 1}),
)


def trial_prime(number):
    if number <= 2 or number % 2 == 0:
        return number == 2
    return all(number % divisor for divisor in range(3, isqrt(number) + 1, 2))


def verify_prime(prime, witness, factors):
    if prod(q**e for q, e in factors.items()) != prime - 1:
        raise ValueError(f"incomplete primality certificate for {prime}")
    if not all(trial_prime(q) for q in factors):
        raise ValueError(f"composite factor in certificate for {prime}")
    if pow(witness, prime - 1, prime) != 1:
        raise ValueError(f"invalid order certificate for {prime}")
    for q in factors:
        if gcd(pow(witness, (prime - 1) // q, prime) - 1, prime) != 1:
            raise ValueError(f"invalid {q}-part of order for {prime}")


def verify(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SWEEP_SHA256:
        raise ValueError("unexpected SHA-256 for the base-two prime sweep")
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    forced = {int(row["prime"]) for row in rows}
    if len(rows) != 789 or len(forced) != 789 or max(forced) != 999375343:
        raise ValueError("unexpected prime sweep totals")
    cutoff_divisor = prod(forced)
    for index, prime, witness, factors in ADDITIONAL:
        verify_prime(prime, witness, factors)
        if cutoff_divisor % index or gcd(index, prime) != 1 or prime in forced:
            raise ValueError(f"invalid forcing index for {prime}")
        residue = pow(2, 2 * index, prime**2)
        if residue % prime != 1 or residue == 1:
            raise ValueError(f"expected a simple divisor at index {index}")
        forced.add(prime)

    required = prod(forced)
    if len(forced) != 795 or not 10**6283 < required < 10**6284:
        raise ValueError("incorrect exponent certificate bound")
    print(f"{len(forced)} forced odd primes")
    print(f"Exponent divisor has {len(str(required))} decimal digits")
    print(f"Exponent divisor exceeds 10^6283")
    print(f"Base-two Wieferich prime count required: at least 2^{len(forced)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sweep", type=Path)
    args = parser.parse_args()
    set_int_max_str_digits(0)
    verify(args.sweep)


if __name__ == "__main__":
    main()
