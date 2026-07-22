#!/usr/bin/env python3
"""Exact Frobenius--Jordan certificates for the odd square family through 1057.

This is intentionally separate from the degree-at-most-thirty regression:
factoring the degree-721 through degree-1057 reductions is exact but slow.
"""
from __future__ import annotations

import math

import sympy as sp

from experiment_fixed_r_modular_factorization import (
    X,
    block_size_compatible,
    factor_degrees,
    isolated_prime_cycle,
)
from master_cancellation import (
    parameter_discriminant,
    parameter_discriminant_is_square,
    parameter_polynomial,
)


# Each sieve maps an unramified auxiliary prime to the complete tuple of
# irreducible factor degrees.  Their subset-sum intersection is {0, degree}.
# The Jordan pair is (auxiliary prime, isolated prime-cycle length).  In a
# composite degree, one displayed type also excludes every proper block size
# by the exact wreath-product compatibility criterion.
CERTIFICATES = {
    161: {
        "sieve": {
            11: (81, 41, 26, 12, 1),
            13: (131, 22, 8),
        },
        "jordan": (13, 131),
        "blocks": (13, (7, 23)),
    },
    241: {
        "sieve": {
            5: (149, 50, 36, 5, 1),
            13: (101, 66, 60, 12, 2),
        },
        "jordan": (5, 149),
    },
    337: {
        "sieve": {
            5: (161, 160, 11, 3, 2),
            17: (203, 109, 25),
        },
        "jordan": (17, 109),
    },
    449: {
        "sieve": {
            19: (320, 126, 3),
            47: (133, 130, 127, 38, 9, 8, 4),
        },
        "jordan": (47, 127),
    },
    577: {
        "sieve": {
            5: (258, 159, 90, 41, 23, 4, 2),
            7: (529, 29, 11, 6, 2),
            19: (374, 104, 53, 34, 5, 4, 3),
        },
        "jordan": (19, 53),
    },
    721: {
        "sieve": {
            5: (481, 128, 102, 9, 1),
            7: (402, 226, 62, 18, 6, 5, 2),
            11: (355, 305, 47, 5, 4, 3, 2),
        },
        "jordan": (11, 47),
        "blocks": (5, (7, 103)),
    },
    881: {
        "sieve": {
            5: (428, 259, 189, 4, 1),
            17: (532, 139, 83, 71, 56),
        },
        "jordan": (17, 139),
    },
    1057: {
        "sieve": {
            7: (974, 41, 32, 8, 2),
            17: (350, 345, 195, 106, 46, 14, 1),
        },
        "jordan": (7, 41),
        "blocks": (7, (7, 151)),
    },
}


def subset_sums(degrees: tuple[int, ...]) -> set[int]:
    sums = {0}
    for degree in degrees:
        sums |= {value + degree for value in tuple(sums)}
    return sums


for degree, certificate in CERTIFICATES.items():
    # These are exactly the consecutive odd parameters a=9,11,...,23 in
    # degree m=2*a^2-1, following the existing 17,49,97 certificates.
    square_argument = (degree + 1) // 2
    a = math.isqrt(square_argument)
    assert a * a == square_argument and a % 2 == 1

    polynomial = sp.Poly(parameter_polynomial(degree, 1, X), X, domain=sp.ZZ)
    discriminant = int(parameter_discriminant(degree, 1))
    possible = set(range(degree + 1))

    for prime, expected_type in certificate["sieve"].items():
        assert discriminant % prime != 0
        actual_type = factor_degrees(polynomial, prime)
        assert actual_type == expected_type
        possible &= subset_sums(actual_type)
    assert possible == {0, degree}

    jordan_prime, cycle_prime = certificate["jordan"]
    jordan_type = certificate["sieve"][jordan_prime]
    assert sp.isprime(cycle_prime) and cycle_prime <= degree - 3
    assert isolated_prime_cycle(jordan_type, cycle_prime)
    other_lengths = list(jordan_type)
    other_lengths.remove(cycle_prime)
    assert math.lcm(*other_lengths) % cycle_prime != 0

    if "blocks" in certificate:
        block_prime, block_sizes = certificate["blocks"]
        assert set(block_sizes) == {
            size for size in range(2, degree // 2 + 1) if degree % size == 0
        }
        block_type = certificate["sieve"][block_prime]
        assert all(
            not block_size_compatible(block_type, block_size)
            for block_size in block_sizes
        )
        primitive_description = (
            f"type mod {block_prime} excludes block sizes {block_sizes}"
        )
    else:
        assert sp.isprime(degree)
        primitive_description = "prime-degree transitivity"

    assert parameter_discriminant_is_square(degree, 1)
    discriminant_root = math.isqrt(discriminant)
    assert discriminant_root * discriminant_root == discriminant
    print(
        f"PASS A{degree}: primes {tuple(certificate['sieve'])}; "
        f"{jordan_type} mod {jordan_prime} isolates a {cycle_prime}-cycle; "
        f"{primitive_description}",
        flush=True,
    )

assert tuple(CERTIFICATES) == (161, 241, 337, 449, 577, 721, 881, 1057)
print(
    "PASS: consecutive odd square-family witnesses Gal=A_m through m=1057",
    flush=True,
)
