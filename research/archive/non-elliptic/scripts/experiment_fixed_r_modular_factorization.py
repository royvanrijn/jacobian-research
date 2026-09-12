#!/usr/bin/env python3
"""Explore Frobenius support and block obstructions on fixed rows.

This is an experiment, not a proof of a density-one theorem.  For each
``(m, r)`` in a bounded box it scans good auxiliary primes and records

* a factorization type isolating a ramification prime ``p | mr+1``; and
* which proper block sizes remain compatible with every observed type.

The transitivity of the root action is an input.  Factorization types are
computed for the monic cancellation polynomial ``M_(m,r)``.
"""
from __future__ import annotations

import argparse
import warnings
from functools import lru_cache

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

from master_cancellation import parameter_polynomial


X = sp.symbols("x")


def parse_rows(value: str) -> tuple[int, ...]:
    rows = tuple(sorted({int(item) for item in value.split(",")}))
    if not rows or rows[0] < 1:
        raise argparse.ArgumentTypeError("rows must be positive integers")
    return rows


def factor_degrees(polynomial: sp.Poly, prime: int) -> tuple[int, ...] | None:
    """Return the squarefree factor degrees, or ``None`` at a bad prime."""
    reduced = sp.Poly(polynomial.as_expr(), X, modulus=prime)
    if reduced.degree() != polynomial.degree():
        return None
    if sp.gcd(reduced, reduced.diff()).degree() != 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SymPyDeprecationWarning)
        factors = sp.factor_list(reduced.as_expr(), modulus=prime)[1]
    degrees: list[int] = []
    for factor, multiplicity in factors:
        degrees.extend([int(sp.degree(factor, X))] * int(multiplicity))
    return tuple(sorted(degrees, reverse=True))


def isolated_prime_cycle(degrees: tuple[int, ...], prime: int) -> bool:
    """Whether a power of this cycle type is a single ``prime``-cycle."""
    return degrees.count(prime) == 1 and all(
        degree == prime or degree % prime != 0 for degree in degrees
    )


def block_size_compatible(degrees: tuple[int, ...], block_size: int) -> bool:
    """Test exact compatibility of one cycle type with a block size.

    For an orbit of ``ell`` blocks, the point cycles assigned to that orbit
    have lengths divisible by ``ell`` and total length ``ell*block_size``.
    Conversely, any partition of the point cycles into subsets with this
    property constructs the restriction of a permutation in the relevant
    imprimitive wreath product.  Thus the following exact-cover test is both
    necessary and sufficient for a *single* cycle type.
    """
    degree = sum(degrees)
    if degree % block_size:
        return False
    count = len(degrees)
    full_mask = (1 << count) - 1
    admissible: list[int] = []

    for mask in range(1, full_mask + 1):
        total = sum(degrees[index] for index in range(count) if mask >> index & 1)
        if total % block_size:
            continue
        block_orbit_length = total // block_size
        if block_orbit_length < 1 or block_orbit_length > degree // block_size:
            continue
        if all(
            degrees[index] % block_orbit_length == 0
            for index in range(count)
            if mask >> index & 1
        ):
            admissible.append(mask)

    by_first_index: list[list[int]] = [[] for _ in range(count)]
    for mask in admissible:
        first_index = (mask & -mask).bit_length() - 1
        by_first_index[first_index].append(mask)

    @lru_cache(maxsize=None)
    def exact_cover(remaining: int) -> bool:
        if remaining == 0:
            return True
        first_index = (remaining & -remaining).bit_length() - 1
        return any(
            subset & remaining == subset and exact_cover(remaining ^ subset)
            for subset in by_first_index[first_index]
        )

    return exact_cover(full_mask)


def candidate_ramification_primes(degree: int, row: int) -> tuple[int, ...]:
    """Primes satisfying the elementary Newton-edge hypotheses.

    This does not test the BFLT bad-pair condition or its exact valuation.
    """
    candidates = []
    for prime, exponent in sp.factorint(degree + 1).items():
        prime = int(prime)
        if exponent == 1 and prime > row and prime <= degree - 3:
            candidates.append(prime)
    return tuple(sorted(candidates, reverse=True))


def proper_block_sizes(degree: int) -> set[int]:
    return {
        size
        for size in range(2, degree // 2 + 1)
        if degree % size == 0
    }


def assert_block_compatibility_engine() -> None:
    """Small exact checks for the algebra and wreath-product engine."""
    for row in range(1, 5):
        for m in range(2, 6):
            degree = m * row
            total = degree + row + 1
            polynomial = parameter_polynomial(m, row, X)
            short_tail = sum(
                (-1) ** k * sp.binomial(total, k) * X**k
                for k in range(row + 1)
            )
            assert sp.expand(
                X ** (row + 1) * polynomial
                - (X - 1) ** total
                + (-1) ** total * short_tail
            ) == 0
    for degree in range(6, 15):
        for size in proper_block_sizes(degree):
            assert block_size_compatible((1,) * degree, size)
            assert not block_size_compatible((degree - 1, 1), size)
    assert block_size_compatible((5, 1, 1, 1, 1, 1, 1, 1), 6)
    assert not block_size_compatible((5, 1, 1, 1, 1, 1, 1, 1), 4)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=parse_rows, default=(1, 2, 3, 4))
    parser.add_argument("--min-m", type=int, default=2)
    parser.add_argument("--max-m", type=int, default=25)
    parser.add_argument("--prime-bound", type=int, default=100)
    parser.add_argument("--show-misses", action="store_true")
    args = parser.parse_args()
    if args.min_m < 1 or args.max_m < args.min_m:
        parser.error("require 1 <= min-m <= max-m")

    assert_block_compatibility_engine()
    auxiliary_primes = tuple(int(p) for p in sp.primerange(2, args.prime_bound + 1))
    total = with_ramification_prime = support_hits = primitive_hits = 0

    for row in args.rows:
        row_total = row_candidates = row_support = row_primitive = 0
        examples: list[str] = []
        for m in range(args.min_m, args.max_m + 1):
            degree = m * row
            if degree < 7:
                continue
            row_total += 1
            total += 1
            candidates = candidate_ramification_primes(degree, row)
            if not candidates:
                continue
            row_candidates += 1
            with_ramification_prime += 1

            polynomial = sp.Poly(parameter_polynomial(m, row, X), X, domain=sp.ZZ)
            surviving_blocks = proper_block_sizes(degree)
            support_certificate: tuple[int, int, tuple[int, ...]] | None = None
            for auxiliary_prime in auxiliary_primes:
                degrees = factor_degrees(polynomial, auxiliary_prime)
                if degrees is None:
                    continue
                if support_certificate is None:
                    for ramification_prime in candidates:
                        if isolated_prime_cycle(degrees, ramification_prime):
                            support_certificate = (
                                ramification_prime,
                                auxiliary_prime,
                                degrees,
                            )
                            break
                if surviving_blocks:
                    surviving_blocks = {
                        size
                        for size in surviving_blocks
                        if block_size_compatible(degrees, size)
                    }
                if support_certificate is not None and not surviving_blocks:
                    break

            if support_certificate is not None:
                row_support += 1
                support_hits += 1
            if not surviving_blocks:
                row_primitive += 1
                primitive_hits += 1
            if args.show_misses and (
                support_certificate is None or surviving_blocks
            ):
                print(
                    f"  MISS m={m}, d={degree}, p={candidates}, "
                    f"support={support_certificate is not None}, "
                    f"blocks={tuple(sorted(surviving_blocks))}"
                )
            if len(examples) < 4 and (support_certificate or not surviving_blocks):
                support_text = "support=-"
                if support_certificate:
                    p, q, degrees = support_certificate
                    support_text = f"support p={p} mod {q}: {degrees}"
                block_text = (
                    "blocks=none"
                    if not surviving_blocks
                    else f"blocks={tuple(sorted(surviving_blocks))}"
                )
                examples.append(f"m={m}, d={degree}, {support_text}, {block_text}")

        print(
            f"r={row}: tested={row_total}, arithmetic-candidates={row_candidates}, "
            f"support={row_support}, blocks-eliminated={row_primitive}"
        )
        for example in examples:
            print(f"  {example}")

    print(
        f"TOTAL tested={total}, arithmetic-candidates={with_ramification_prime}, "
        f"support={support_hits}, blocks-eliminated={primitive_hits}"
    )


if __name__ == "__main__":
    main()
