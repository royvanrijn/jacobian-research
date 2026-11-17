#!/usr/bin/env python3
"""Slow exact replay for the diagonal columns 500 <= m <= 741.

An Arb interval proves that two iterations of Dusart's prime-interval bound
fit in the diagonal interval exactly through m=741.  This script certifies
every residual pair below the analytic threshold in the newly added columns.
"""
from __future__ import annotations

import bisect
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
LOW_COLUMN = 500
HIGH_COLUMN = 741


def certify_analytic_endpoint() -> None:
    """Prove the two strict endpoint inequalities with Arb intervals."""
    previous_precision = ctx.prec
    try:
        ctx.prec = 128
        logarithm = arb(DUSART_THRESHOLD).log()
        two_step_ratio = (1 + 1 / logarithm**3) ** 2
        assert two_step_ratio < 1 + arb(1) / HIGH_COLUMN
        assert two_step_ratio > 1 + arb(1) / (HIGH_COLUMN + 1)
    finally:
        ctx.prec = previous_precision


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
    certify_analytic_endpoint()
    tasks = residual_tasks()
    assert len(tasks) == 2_899
    assert sum(task[4] is None for task in tasks) == 1_126
    assert sum(task[4] is not None for task in tasks) == 1_773
    assert max(task[2] for task in tasks) == 58_806
    assert max(task[4] or 0 for task in tasks) == 57
    assert max(task[2] for task in tasks if task[4] is None) == 34_068

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
            if index % 200 == 0:
                print(
                    f"sharp frontier progress: {index}/{len(tasks)} "
                    f"({time.time() - start:.1f}s)",
                    flush=True,
                )

    assert not failures
    assert max_prime == 127
    assert max_steps == 20
    print(
        "PASS: all 2899 residual pairs in 500<=m<=741 are irreducible "
        f"(largest auxiliary prime {max_prime}, at most {max_steps} good steps)",
        flush=True,
    )


if __name__ == "__main__":
    main()
