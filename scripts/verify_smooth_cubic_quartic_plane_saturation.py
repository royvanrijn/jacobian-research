#!/usr/bin/env python3
"""Exact quartic coordinate-plane audit for the smooth cubic symbol.

For every pair in a fixed basis of the 24-dimensional order-four tensor
kernel, test the full polynomial family

    phi_h + u * psi_i + v * psi_j

over Q[u,v,x,y,z].  These 276 coordinate planes include every specialization
on every basis pair.  They do not test directions supported on three or more
basis elements.
"""

from __future__ import annotations

import itertools
import multiprocessing

import verify_cubic_symbol_double_saturation as cubic_audit


def audit_pair(pair: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, ...]]:
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    directions = cubic_audit.quartic_kernel_basis_tensors()
    first, second = pair
    result = cubic_audit.run_singular_plane(
        cubic_audit.CUBIC_STRATA["smooth"],
        directions[first],
        directions[second],
        timeout=180,
    )
    return pair, result


def main() -> None:
    pairs = list(itertools.combinations(range(24), 2))
    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        for pair, result in pool.imap_unordered(audit_pair, pairs):
            assert result == (0, 6, 0, 0, 3), (pair, result)

    print(
        "PASS: all 276 quartic coordinate planes for the smooth cubic "
        "have saturated cotangent presentation and a relative length-six "
        "Ext presentation pulled back from the parameter-plane origin"
    )


if __name__ == "__main__":
    main()
