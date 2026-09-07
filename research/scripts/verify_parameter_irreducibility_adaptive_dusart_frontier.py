#!/usr/bin/env python3
"""Slow exact replay for the adaptive diagonal columns 742 <= m <= 1000.

For each column, 128-bit Arb arithmetic finds the first integer r at which
two iterations of Dusart's prime-interval bound fit inside the diagonal
interval.  The bound handles that r and every larger one.  This script
certifies every residual pair at smaller r.
"""
from __future__ import annotations

import bisect
from fractions import Fraction
import multiprocessing as mp
import time

from flint import arb, ctx
import sympy as sp

try:
    from verify_parameter_irreducibility_dusart_frontier import certify_residual
except ModuleNotFoundError:
    from scripts.verify_parameter_irreducibility_dusart_frontier import (
        certify_residual,
    )


DUSART_THRESHOLD = 89_693
LOW_COLUMN = 742
HIGH_COLUMN = 1000


def certify_general_region_formulas() -> None:
    """Check the exact Dusart curve and the BHP exponent conversion."""
    theta = Fraction(21, 40)
    assert theta / (1 - theta) == Fraction(21, 19)
    assert 1 / (1 - theta) == Fraction(40, 19)

    previous_precision = ctx.prec
    try:
        ctx.prec = 128
        for degree in (DUSART_THRESHOLD, 90_524, 297_000):
            logarithm = arb(degree).log()
            reciprocal_cube = 1 / logarithm**3
            exact_delta = reciprocal_cube + (
                (1 + reciprocal_cube)
                / (logarithm + (1 + reciprocal_cube).log()) ** 3
            )
            coarse_delta = 2 * reciprocal_cube + reciprocal_cube**2
            assert exact_delta < coarse_delta
    finally:
        ctx.prec = previous_precision


def analytic_starts() -> list[tuple[int, int, int]]:
    """Certify the minimal Dusart-tail start r in every selected column."""
    previous_precision = ctx.prec
    starts: list[tuple[int, int, int]] = []
    try:
        ctx.prec = 128
        for m in range(LOW_COLUMN, HIGH_COLUMN + 1):
            r = (DUSART_THRESHOLD + m - 1) // m
            target_ratio = 1 + arb(1) / m
            while True:
                degree = m * r
                two_step_ratio = (
                    1 + 1 / arb(degree).log() ** 3
                ) ** 2
                if two_step_ratio < target_ratio:
                    break
                assert two_step_ratio > target_ratio
                r += 1
            assert degree >= DUSART_THRESHOLD
            if r > (DUSART_THRESHOLD + m - 1) // m:
                preceding_ratio = (
                    1 + 1 / arb(m * (r - 1)).log() ** 3
                ) ** 2
                assert preceding_ratio > target_ratio
            starts.append((m, r, degree))
    finally:
        ctx.prec = previous_precision
    return starts


def residual_tasks() -> list[tuple[int, int, int, int, int | None]]:
    """Enumerate the exact finite frontier after prime criteria are removed."""
    starts = analytic_starts()
    prime_limit = max(
        m * (start_r - 1) + start_r + 1 for m, start_r, _ in starts
    ) + 1
    primes = [int(prime) for prime in sp.primerange(2, prime_limit)]
    prime_flags = bytearray(prime_limit)
    for prime in primes:
        prime_flags[prime] = 1
    prime_prefix = [0] * prime_limit
    for value in range(1, prime_limit):
        prime_prefix[value] = prime_prefix[value - 1] + prime_flags[value]

    tasks: list[tuple[int, int, int, int, int | None]] = []
    for m, start_r, _ in starts:
        for r in range(1, start_r):
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
    certify_general_region_formulas()
    starts = analytic_starts()
    assert starts[0] == (742, 122, 90_524)
    assert starts[-1] == (1000, 297, 297_000)

    tasks = residual_tasks()
    assert len(tasks) == 3_335
    assert sum(task[4] is None for task in tasks) == 1_294
    assert sum(task[4] is not None for task in tasks) == 2_041
    assert max(task[2] for task in tasks) == 58_823
    assert max(task[4] or 0 for task in tasks) == 54
    assert max(task[2] for task in tasks if task[4] is None) == 44_298

    start = time.time()
    failures: list[tuple[int, int, int]] = []
    max_prime = 0
    max_steps = 0
    context = mp.get_context("fork")
    with context.Pool(processes=6) as pool:
        results = pool.imap_unordered(certify_residual, tasks, chunksize=1)
        for index, result in enumerate(results, 1):
            m, r, certified, prime, steps = result
            if not certified:
                failures.append((m, r, steps))
            max_prime = max(max_prime, prime)
            max_steps = max(max_steps, steps)
            if index % 100 == 0:
                print(
                    f"adaptive frontier progress: {index}/{len(tasks)} "
                    f"({time.time() - start:.1f}s)",
                    flush=True,
                )

    assert not failures
    assert max_prime == 107
    assert max_steps == 18
    print(
        "PASS: all 3335 residual pairs in 742<=m<=1000 are irreducible "
        f"(largest auxiliary prime {max_prime}, at most {max_steps} good steps)",
        flush=True,
    )


if __name__ == "__main__":
    main()
