#!/usr/bin/env python3
"""Exact Frobenius--Jordan certificates, including two odd square-family cases."""
from __future__ import annotations

import math
import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

from master_cancellation import (
    parameter_discriminant,
    parameter_discriminant_is_square,
    parameter_polynomial,
)


q = sp.symbols("q")

# pair: (prime giving (n-1,1), prime giving the extractable cycle, cycle).
# The first entry is None when prime-degree transitivity supplies primitivity.
CERTIFICATES = {
    (1, 7): (7, 19, 2),
    (7, 1): (11, 7, 3),
    (1, 8): (43, 7, 2),
    (2, 4): (29, 2, 3),
    (4, 2): (17, 2, 3),
    (8, 1): (41, 13, 3),
    (1, 9): (43, 31, 5),
    (3, 3): (7, 29, 2),
    (9, 1): (23, 3, 2),
    (1, 10): (37, 31, 5),
    (2, 5): (29, 23, 2),
    (5, 2): (43, 5, 3),
    (10, 1): (31, 7, 3),
    (1, 11): (103, 41, 3),
    (11, 1): (23, 17, 5),
    (1, 12): (43, 11, 5),
    (2, 6): (23, 29, 3),
    (3, 4): (107, 29, 2),
    (4, 3): (3, 11, 5),
    (6, 2): (79, 23, 7),
    (12, 1): (59, 5, 2),
    (1, 13): (13, 37, 3),
    (13, 1): (97, 11, 2),
    (1, 14): (83, 13, 2),
    (2, 7): (47, 13, 5),
    (7, 2): (41, 13, 7),
    (14, 1): (59, 7, 5),
    (1, 15): (37, 11, 3),
    (3, 5): (29, 11, 3),
    (5, 3): (37, 13, 3),
    (15, 1): (7, 3, 3),
    (1, 16): (101, 7, 3),
    (2, 8): (41, 2, 5),
    (4, 4): (103, 29, 3),
    (8, 2): (139, 2, 3),
    (16, 1): (37, 5, 3),
    (1, 17): (331, 13, 2),
    (17, 1): (None, 7, 7),
    (1, 18): (73, 43, 2),
    (2, 9): (59, 41, 3),
    (3, 6): (113, 13, 2),
    (6, 3): (31, 13, 11),
    (9, 2): (71, 11, 5),
    (18, 1): (7, 17, 5),
    (1, 19): (43, 53, 2),
    (19, 1): (43, 17, 13),
    (1, 20): (157, 17, 5),
    (2, 10): (67, 17, 7),
    (4, 5): (41, 17, 2),
    (5, 4): (109, 13, 13),
    (10, 2): (137, 2, 5),
    (20, 1): (19, 5, 3),
    (1, 21): (137, 7, 17),
    (3, 7): (67, 17, 5),
    (7, 3): (127, 17, 7),
    (21, 1): (13, 7, 3),
    (1, 22): (67, 17, 3),
    (2, 11): (127, 19, 13),
    (11, 2): (229, 11, 2),
    (22, 1): (79, 7, 3),
    (1, 23): (None, 17, 7),
    (23, 1): (17, 23, 7),
    (1, 24): (101, 53, 13),
    (2, 12): (127, 11, 2),
    (3, 8): (17, 7, 3),
    (4, 6): (367, 17, 5),
    (6, 4): (37, 2, 5),
    (8, 3): (113, 17, 3),
    (12, 2): (131, 2, 5),
    (24, 1): (17, 7, 2),
    (1, 25): (89, 11, 3),
    (5, 5): (73, 5, 5),
    (25, 1): (11, 5, 7),
    (1, 26): (83, 5, 7),
    (2, 13): (73, 11, 23),
    (13, 2): (223, 11, 3),
    (26, 1): (13, 17, 5),
    (1, 27): (79, 13, 2),
    (3, 9): (73, 3, 13),
    (9, 3): (541, 3, 13),
    (27, 1): (47, 3, 7),
    (1, 28): (307, 23, 7),
    (2, 14): (373, 23, 5),
    (4, 7): (191, 13, 2),
    (7, 4): (101, 7, 7),
    (14, 2): (11, 7, 5),
    (28, 1): (53, 7, 7),
    (1, 29): (29, 23, 23),
    (29, 1): (None, 7, 17),
    (1, 30): (149, 5, 7),
    (2, 15): (107, 5, 23),
    (3, 10): (443, 3, 2),
    (5, 6): (59, 13, 19),
    (6, 5): (83, 19, 2),
    (10, 3): (37, 3, 11),
    (15, 2): (37, 5, 17),
    (30, 1): (13, 3, 3),
}


def modular_factor_degrees(polynomial: sp.Poly, prime: int) -> tuple[int, ...]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SymPyDeprecationWarning)
        factors = sp.factor_list(polynomial.as_expr(), modulus=prime)[1]
    degrees: list[int] = []
    for factor, multiplicity in factors:
        degrees.extend([int(sp.degree(factor, q))] * multiplicity)
    return tuple(sorted(degrees, reverse=True))


def subset_sums(degrees: tuple[int, ...]) -> set[int]:
    sums = {0}
    for degree in degrees:
        sums |= {value + degree for value in tuple(sums)}
    return sums


def assert_irreducible_by_degree_sieve(polynomial: sp.Poly) -> None:
    possible = set(range(polynomial.degree() + 1))
    for prime in sp.primerange(2, 100):
        possible &= subset_sums(modular_factor_degrees(polynomial, int(prime)))
        if possible == {0, polynomial.degree()}:
            return
    raise AssertionError(f"degree sieve left possible factor degrees {possible}")


for pair, (primitive_prime, jordan_prime, cycle_prime) in CERTIFICATES.items():
    m, r = pair
    degree = m * r
    polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
    discriminant = int(parameter_discriminant(m, r))

    assert 7 <= degree <= 30
    assert_irreducible_by_degree_sieve(polynomial)

    # Dedekind's theorem applies because both reductions are unramified.
    assert discriminant % jordan_prime != 0
    jordan_type = modular_factor_degrees(polynomial, jordan_prime)
    if primitive_prime is None:
        # A transitive permutation group of prime degree is primitive.
        assert sp.isprime(degree)
        primitive_description = "prime-degree primitivity"
    else:
        assert discriminant % primitive_prime != 0
        primitive_type = modular_factor_degrees(polynomial, primitive_prime)
        assert primitive_type == (degree - 1, 1)
        primitive_description = f"({degree - 1},1) mod {primitive_prime}"

    # Taking the lcm of the other cycle lengths kills those cycles but leaves
    # the unique cycle_prime-cycle nontrivial, hence isolates that prime cycle.
    assert sp.isprime(cycle_prime) and cycle_prime <= degree - 3
    assert jordan_type.count(cycle_prime) == 1
    other_lengths = list(jordan_type)
    other_lengths.remove(cycle_prime)
    exponent = math.lcm(*other_lengths)
    assert exponent % cycle_prime != 0

    square_discriminant = parameter_discriminant_is_square(m, r)
    expected_alternating = pair in {
        (4, 3),
        (2, 8),
        (16, 1),
        (17, 1),
        (1, 24),
        (12, 2),
    }
    assert square_discriminant == expected_alternating
    group = f"A{degree}" if square_discriminant else f"S{degree}"
    print(
        f"PASS (m,r)={pair}: {primitive_description}; "
        f"{jordan_type} mod {jordan_prime} isolates a {cycle_prime}-cycle; "
        f"Gal={group}"
    )

assert set(CERTIFICATES) == {
    (m, r)
    for m in range(1, 31)
    for r in range(1, 31)
    if 7 <= m * r <= 30
}
print("PASS: complete Frobenius--Jordan classification in degrees 7--30")

# The next two members of the odd square-discriminant family
# (m,r)=(2*a^2-1,1), after (17,1).  The listed unramified factorizations give
# exact factor-degree sieves and isolate the displayed prime cycles.
ODD_SQUARE_CERTIFICATES = {
    (49, 1): {
        "sieve": {
            7: (22, 13, 9, 3, 2),
            11: (22, 10, 8, 6, 3),
            13: (24, 15, 10),
            19: (34, 12, 3),
        },
        "jordan": (7, 13),
    },
    (97, 1): {
        "sieve": {
            5: (89, 6, 2),
            19: (91, 5, 1),
            23: (83, 12, 2),
        },
        "jordan": (5, 89),
    },
}

for pair, certificate in ODD_SQUARE_CERTIFICATES.items():
    m, r = pair
    degree = m * r
    polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
    discriminant = int(parameter_discriminant(m, r))
    possible = set(range(degree + 1))

    for prime, expected_type in certificate["sieve"].items():
        assert discriminant % prime != 0
        actual_type = modular_factor_degrees(polynomial, prime)
        assert actual_type == expected_type
        possible &= subset_sums(actual_type)
    assert possible == {0, degree}

    jordan_prime, cycle_prime = certificate["jordan"]
    jordan_type = certificate["sieve"][jordan_prime]
    assert sp.isprime(cycle_prime) and cycle_prime <= degree - 3
    assert jordan_type.count(cycle_prime) == 1
    other_lengths = list(jordan_type)
    other_lengths.remove(cycle_prime)
    assert math.lcm(*other_lengths) % cycle_prime != 0

    if degree == 49:
        # A transitive imprimitive group of degree 7^2 embeds in S_7 wr S_7,
        # whose order has no prime divisor 13.  The isolated 13-cycle therefore
        # rules out the only possible nontrivial block size.
        assert cycle_prime == 13 and degree == 7 * 7 and cycle_prime > 7
        primitive_description = "13-cycle excludes the 7-by-7 block system"
    else:
        assert sp.isprime(degree)
        primitive_description = "prime-degree transitivity"

    assert parameter_discriminant_is_square(m, r)
    print(
        f"PASS (m,r)={pair}: exact modular degree sieve; "
        f"{jordan_type} mod {jordan_prime} isolates a {cycle_prime}-cycle; "
        f"{primitive_description}; Gal=A{degree}"
    )

print("PASS: odd square-family certificates Gal=A49 and Gal=A97")
