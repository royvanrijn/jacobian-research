#!/usr/bin/env python3
"""Bounded exact search in the fixed integral Hasse pencil.

For the target (u,v,1), the inverse polynomial is

    E = S^5 - 3*S^4 + 8*S^3 + 4*v*S^2 - 8*S + 16*u.

Every monic integral factorization E=(cubic)*(quadratic) is parameterized by
the two leading cubic coefficients ``a,b``.  This script exhausts
|a|,|b| <= 2000, retains irreducible 3+2 factorizations whose discriminant
fields agree, and checks every bad prime for an integral source point.

The bounded search is not an infinitude theorem.  Each returned target,
however, receives an exact all-prime Hasse certificate.
"""

from __future__ import annotations

from itertools import product
import json
from math import isqrt
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "integral_hasse_pencil_search.json"
)
BOUND = 2000
S = sp.symbols("S")


def cubic_discriminant(a: int, b: int, c: int) -> int:
    return (
        a * a * b * b
        - 4 * b**3
        - 4 * a**3 * c
        - 27 * c * c
        + 18 * a * b * c
    )


def squarefree_part(value: int) -> int:
    result = -1 if value < 0 else 1
    for prime, exponent in sp.factorint(abs(value)).items():
        if exponent % 2:
            result *= int(prime)
    return result


def cubic_has_integer_root(a: int, b: int, c: int) -> bool:
    for positive in sp.divisors(abs(c)):
        for root in (int(positive), -int(positive)):
            if root**3 + a * root**2 + b * root + c == 0:
                return True
    return False


def inverse_value(u: int, v: int, value: int, modulus: int) -> int:
    return (
        value**5
        - 3 * value**4
        + 8 * value**3
        + 4 * v * value**2
        - 8 * value
        + 16 * u
    ) % modulus


def inverse_derivative(v: int, value: int, modulus: int) -> int:
    return (
        5 * value**4
        - 12 * value**3
        + 24 * value**2
        + 8 * v * value
        - 8
    ) % modulus


def first_simple_root(u: int, v: int, prime: int) -> int | None:
    for residue in range(prime):
        if (
            inverse_value(u, v, residue, prime) == 0
            and inverse_derivative(v, residue, prime) != 0
        ):
            return residue
    return None


def map_mod(
    source: tuple[int, int, int],
    modulus: int,
) -> tuple[int, int, int]:
    x, y, z = (coordinate % modulus for coordinate in source)
    t = (1 + 2 * x * y) % modulus
    q = (t**2 * z - y**2 * (1 + 3 * t)) % modulus
    return (
        (
            x * (1 - 3 * x * y)
            + 2 * x**3 * z
            - 3 * x**4 * q**4
            + 3 * x**5 * q**5
        )
        % modulus,
        (
            y
            - 6 * x * q
            + 6 * t**2 * x**2 * q**4
            - 5 * t**2 * x**3 * q**5
        )
        % modulus,
        t * q % modulus,
    )


def source_at_odd_prime(
    u: int,
    v: int,
    root: int,
    prime: int,
) -> tuple[int, int, int]:
    inverse_2 = pow(2, -1, prime)
    inverse_8 = pow(8, -1, prime)
    d = inverse_derivative(v, root, prime) * pow(-8, -1, prime) % prime
    t = pow(d, -1, prime)
    x = root * pow(2 * d, -1, prime) % prime
    beta = (
        1
        - 4 * root
        + 3 * inverse_2 * root**2
        - 5 * inverse_8 * root**3
    ) % prime
    y = (v + 1 - beta - root) % prime
    z = d**2 * (d + y**2 * (1 + 3 * t)) % prime
    return x, y, z


def source_at_two(u: int, v: int) -> tuple[int, int, int] | None:
    target = (u % 2, v % 2, 1)
    for source in product(range(2), repeat=3):
        if map_mod(source, 2) == target:
            return source
    return None


def candidate_record(
    a: int,
    b: int,
    c: int,
    d: int,
    e: int,
    u: int,
    v: int,
) -> dict[str, object] | None:
    cubic = S**3 + a * S**2 + b * S + c
    quadratic = S**2 + d * S + e
    inverse = S**5 - 3 * S**4 + 8 * S**3 + 4 * v * S**2 - 8 * S + 16 * u
    assert sp.expand(cubic * quadratic - inverse) == 0
    disc_cubic = cubic_discriminant(a, b, c)
    disc_quadratic = d * d - 4 * e
    resultant = int(sp.resultant(cubic, quadratic, S))
    if resultant == 0:
        return None

    bad_primes = sorted(
        int(prime)
        for prime in sp.factorint(
            abs(disc_cubic * disc_quadratic * resultant)
        )
    )
    two_source = source_at_two(u, v)
    if two_source is None:
        return None

    local_roots: dict[str, int] = {}
    local_sources: dict[str, list[int]] = {"2": list(two_source)}
    for prime in bad_primes:
        if prime == 2:
            continue
        root = first_simple_root(u, v, prime)
        if root is None:
            return None
        source = source_at_odd_prime(u, v, root, prime)
        assert map_mod(source, prime) == (u % prime, v % prime, 1)
        local_roots[str(prime)] = root
        local_sources[str(prime)] = list(source)

    return {
        "target": [u, v, 1],
        "cubic_coefficients": [1, a, b, c],
        "quadratic_coefficients": [1, d, e],
        "cubic_discriminant": disc_cubic,
        "quadratic_discriminant": disc_quadratic,
        "common_squarefree_part": squarefree_part(disc_cubic),
        "factor_resultant": resultant,
        "bad_primes": bad_primes,
        "simple_inverse_roots": local_roots,
        "source_points_mod_bad_primes": local_sources,
    }


def run_search() -> list[dict[str, object]]:
    hits: dict[tuple[int, int], dict[str, object]] = {}
    for a in range(-BOUND, BOUND + 1):
        d = -3 - a
        if d == 0:
            continue
        for b in range(-BOUND, BOUND + 1):
            e = a * a + 3 * a + 8 - b
            numerator = -8 - b * e
            if numerator % d:
                continue
            c = numerator // d
            if c == 0 or e == 0:
                continue

            v_numerator = c + b * d + a * e
            if v_numerator % 4 or c * e % 16:
                continue
            v = v_numerator // 4
            u = c * e // 16

            disc_cubic = cubic_discriminant(a, b, c)
            disc_quadratic = d * d - 4 * e
            if disc_cubic == 0 or disc_quadratic == 0:
                continue
            if squarefree_part(disc_cubic) != squarefree_part(
                disc_quadratic
            ):
                continue
            if cubic_has_integer_root(a, b, c):
                continue
            if (
                disc_quadratic >= 0
                and isqrt(disc_quadratic) ** 2 == disc_quadratic
            ):
                continue

            record = candidate_record(a, b, c, d, e, u, v)
            if record is not None:
                hits[(u, v)] = record

    return sorted(
        hits.values(),
        key=lambda record: (
            max(abs(record["target"][0]), abs(record["target"][1])),
            record["target"],
        ),
    )


results = {
    "search": {
        "factorization": (
            "(S^3+a*S^2+b*S+c)(S^2+d*S+e)"
        ),
        "coefficient_box": {"a": [-BOUND, BOUND], "b": [-BOUND, BOUND]},
        "scope": (
            "Exhaustive only for monic integral 3+2 factorizations in the "
            "stated a,b box; this is not an infinitude result."
        ),
    },
    "hits": run_search(),
}

expected = json.loads(RESULT_PATH.read_text())
assert results == expected
assert len(results["hits"]) == 3

for hit in results["hits"]:
    target = hit["target"]
    print(
        "PASS target",
        tuple(target),
        "bad primes",
        hit["bad_primes"],
    )
print("PASS: exhaustive bounded search agrees with", RESULT_PATH)
