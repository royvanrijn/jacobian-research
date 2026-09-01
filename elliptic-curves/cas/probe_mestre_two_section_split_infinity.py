#!/usr/bin/env python3
"""Finite-field split-infinity audit on the recovered two-section surface.

This is an experiment, not a characteristic-zero square certificate.  It
enumerates only the exact root surface ``M=H=0, edQ != 0`` recovered in
``verify_mestre_two_section_root_surface.sing`` and records whether its
leading invariant D is zero, a square, or a nonsquare at root-distinct points.
It also checks D exactly at the two known rational points on the surface.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations, product


DEFAULT_PRIMES = (17, 29, 31, 37)


def elementary(roots: tuple[object, ...], degree: int) -> object:
    return sum(product_value(items) for items in combinations(roots, degree))


def product_value(items: tuple[object, ...]) -> object:
    answer = 1
    for item in items:
        answer *= item
    return answer


def mestre_and_leading(roots4: tuple[object, object, object, object]) -> tuple[object, object]:
    roots = (0, 1, *roots4)
    a1 = -elementary(roots, 1)
    a2 = elementary(roots, 2)
    a3 = -elementary(roots, 3)
    a4 = elementary(roots, 4)
    a5 = -elementary(roots, 5)
    mestre = (
        a1**5 - 6 * a1**3 * a2 + 7 * a1**2 * a3 + 8 * a1 * a2**2
        - 8 * a1 * a4 - 12 * a2 * a3 + 24 * a5
    )
    leading = 5 * a1**4 - 24 * a1**2 * a2 + 32 * a1 * a3 + 16 * a2**2 - 64 * a4
    return mestre, leading


def separator(r3: int, r4: int, r5: int) -> int:
    return (
        r3**3 * r4 - 2 * r3**2 * r4**2 + r3 * r4**3 - r3**3 * r5
        + r3**2 * r4 * r5 + r3 * r4**2 * r5 - r4**3 * r5
        + r3**2 * r5**2 - 2 * r3 * r4 * r5**2 + r4**2 * r5**2
        - r3**3 + 2 * r3**2 * r4 - 2 * r3 * r4**2 + r4**3
        - r3**2 * r5 + r4**2 * r5 + 2 * r3 * r5**2 - 2 * r4 * r5**2
        + 2 * r3 * r4 - 2 * r4**2 - r3 * r5 + r4 * r5 + r5**2
        - r3 + r4 - r5
    )


def audit_prime(prime: int) -> dict[str, int]:
    squares = {value * value % prime for value in range(1, prime)}
    counts = {"regular_points": 0, "D_zero": 0, "D_square_nonzero": 0, "D_nonsquare": 0}
    for roots4 in product(range(prime), repeat=4):
        r3, r4, r5, r6 = roots4
        if len({0, 1, *roots4}) != 6:
            continue
        e = (r3 - r4 + 1) % prime
        d = (r3 - r5 - r6 + 1) % prime
        if not e or not d or separator(r3, r4, r5) % prime == 0:
            continue
        mestre, leading = mestre_and_leading(roots4)
        if mestre % prime:
            continue
        sparse = (
            r3 * pow(e, -1, prime) + r5 + r6 - r4
            - (r3 - r5 * r6) * pow(d, -1, prime)
        ) % prime
        if sparse:
            continue
        counts["regular_points"] += 1
        leading %= prime
        if leading == 0:
            counts["D_zero"] += 1
        elif leading in squares:
            counts["D_square_nonzero"] += 1
        else:
            counts["D_nonsquare"] += 1
    return counts


def rational_checks() -> list[dict[str, str]]:
    points = (
        (
            (Fraction(19, 5), Fraction(143, 25), Fraction(168, 25), Fraction(41, 5)),
            Fraction(528, 25),
        ),
        (
            (Fraction(175, 23), Fraction(93, 23), Fraction(128, 23), Fraction(133, 23)),
            Fraction(14400, 529),
        ),
    )
    result = []
    for roots4, expected_root in points:
        mestre, leading = mestre_and_leading(roots4)
        if mestre != 0 or leading != expected_root**2:
            raise AssertionError("a rational surface point or its declared square root failed")
        result.append(
            {
                "roots4": [str(value) for value in roots4],
                "D": str(leading),
                "sqrt_D": str(expected_root),
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", nargs="*", type=int, default=DEFAULT_PRIMES)
    args = parser.parse_args()
    result = {
        "status": "finite-field experiment; not a function-field square proof",
        "rational_points": rational_checks(),
        "finite_fields": {str(prime): audit_prime(prime) for prime in args.primes},
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
