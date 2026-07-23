#!/usr/bin/env python3
"""Slow exact replay for the diagonal columns 301 <= m <= 499.

The ordinary parameter-irreducibility regression proves the complete block
through m=300.  Dusart's explicit prime-interval theorem handles mr>=89693
through m=499.  This script certifies every residual pair below that threshold.
"""
from __future__ import annotations

import bisect
import math
import multiprocessing as mp
import time

from flint import nmod_poly
import sympy as sp


DUSART_THRESHOLD = 89_693
LOW_COLUMN = 301
HIGH_COLUMN = 499
AUXILIARY_PRIMES = [int(prime) for prime in sp.primerange(2, 300)]


def binomial_mod_prime(total: int, index: int, prime: int) -> int:
    """Return binomial(total,index) modulo a prime via Lucas's theorem."""
    result = 1
    while total or index:
        total_digit = total % prime
        index_digit = index % prime
        if index_digit > total_digit:
            return 0
        result = result * math.comb(total_digit, index_digit) % prime
        total //= prime
        index //= prime
    return result


def modular_parameter_polynomial(
    degree: int, total: int, prime: int
) -> nmod_poly:
    """Construct M_(m,r) modulo ``prime`` without huge integer coefficients."""
    return nmod_poly(
        [
            (-1) ** (degree - index)
            * binomial_mod_prime(total, degree - index, prime)
            % prime
            for index in range(degree + 1)
        ],
        prime,
    )


def subset_sum_reaches(degrees: list[int], target: int) -> bool:
    """Decide whether a submultiset of ``degrees`` sums to ``target``."""
    sums = 1
    mask = (1 << (target + 1)) - 1
    for degree in degrees:
        if degree <= target:
            sums |= (sums << degree) & mask
    return bool((sums >> target) & 1)


def modular_small_factor_degrees(
    degree: int, total: int, prime: int, maximum: int
) -> list[int] | None:
    """Return factor degrees through ``maximum`` at a squarefree prime."""
    polynomial = modular_parameter_polynomial(degree, total, prime)
    if len(polynomial.gcd(polynomial.derivative())) > 1:
        return None

    variable = nmod_poly([0, 1], prime)
    frobenius = variable
    factor_counts = [0] * (maximum + 1)
    degrees: list[int] = []
    for current_degree in range(1, maximum + 1):
        frobenius = frobenius.pow_mod(prime, polynomial)
        divisor_degree = len((frobenius - variable).gcd(polynomial)) - 1
        known_degree = sum(
            divisor * factor_counts[divisor]
            for divisor in range(1, current_degree)
            if current_degree % divisor == 0
        )
        count = (divisor_degree - known_degree) // current_degree
        assert count >= 0
        assert known_degree + current_degree * count == divisor_degree
        factor_counts[current_degree] = count
        degrees.extend([current_degree] * count)
    return degrees


def certify_residual(
    task: tuple[int, int, int, int, int | None],
) -> tuple[int, int, bool, int, int]:
    """Return a bounded exact certificate result for one residual pair."""
    m, r, degree, total, forced_degree = task
    steps = 0

    if forced_degree is not None:
        for prime in AUXILIARY_PRIMES:
            degrees = modular_small_factor_degrees(
                degree, total, prime, forced_degree
            )
            if degrees is None:
                continue
            steps += 1
            if not subset_sum_reaches(degrees, forced_degree):
                return m, r, True, prime, steps
        return m, r, False, 0, steps

    possible_degrees = (1 << (degree + 1)) - 1
    full_mask = possible_degrees
    for prime in AUXILIARY_PRIMES:
        polynomial = modular_parameter_polynomial(degree, total, prime)
        if len(polynomial.gcd(polynomial.derivative())) > 1:
            continue
        factor_degrees: list[int] = []
        for factor, multiplicity in polynomial.factor()[1]:
            factor_degrees.extend([len(factor) - 1] * multiplicity)
        modular_mask = 1
        for factor_degree in factor_degrees:
            modular_mask |= (modular_mask << factor_degree) & full_mask
        possible_degrees &= modular_mask
        steps += 1
        if possible_degrees == 1 | (1 << degree):
            return m, r, True, prime, steps
    return m, r, False, 0, steps


def residual_tasks() -> list[tuple[int, int, int, int, int | None]]:
    """Enumerate the exact finite frontier after prime criteria are removed."""
    prime_limit = DUSART_THRESHOLD + (DUSART_THRESHOLD - 1) // 2 + 3
    primes = [int(prime) for prime in sp.primerange(2, prime_limit)]
    prime_flags = bytearray(prime_limit)
    for prime in primes:
        prime_flags[prime] = 1
    prime_prefix = [0] * prime_limit
    for value in range(1, prime_limit):
        prime_prefix[value] = prime_prefix[value - 1] + prime_flags[value]

    tasks: list[tuple[int, int, int, int, int | None]] = []
    for m in range(LOW_COLUMN, HIGH_COLUMN + 1):
        for r in range(1, (DUSART_THRESHOLD - 1) // m + 1):
            degree = m * r
            total = degree + r + 1
            interval_prime_count = (
                prime_prefix[total - 1] - prime_prefix[degree]
            )
            if prime_flags[total] or interval_prime_count >= 2:
                continue
            forced_degree: int | None = None
            if interval_prime_count == 1:
                prime_index = bisect.bisect_right(primes, degree)
                interval_prime = primes[prime_index]
                assert degree < interval_prime < total
                forced_degree = total - interval_prime
            tasks.append((m, r, degree, total, forced_degree))
    return tasks


def main() -> None:
    assert 1001**2 * HIGH_COLUMN < (HIGH_COLUMN + 1) * 1000**2
    tasks = residual_tasks()
    assert len(tasks) == 2_192
    assert sum(task[4] is None for task in tasks) == 816
    assert sum(task[4] is not None for task in tasks) == 1_376
    assert max(task[2] for task in tasks) == 31_395
    assert max(task[4] or 0 for task in tasks) == 68
    assert max(task[2] for task in tasks if task[4] is None) == 19_614

    start = time.time()
    failures: list[tuple[int, int, int]] = []
    max_prime = 0
    max_steps = 0
    context = mp.get_context("fork")
    with context.Pool(processes=4) as pool:
        results = pool.imap_unordered(certify_residual, tasks, chunksize=1)
        for index, result in enumerate(results, 1):
            m, r, certified, prime, steps = result
            if not certified:
                failures.append((m, r, steps))
            max_prime = max(max_prime, prime)
            max_steps = max(max_steps, steps)
            if index % 200 == 0:
                print(
                    f"frontier progress: {index}/{len(tasks)} "
                    f"({time.time() - start:.1f}s)",
                    flush=True,
                )

    assert not failures
    assert max_prime == 103
    assert max_steps == 13
    print(
        "PASS: all 2192 residual pairs in 301<=m<=499 are irreducible "
        f"(largest auxiliary prime {max_prime}, at most {max_steps} good steps)",
        flush=True,
    )


if __name__ == "__main__":
    main()
