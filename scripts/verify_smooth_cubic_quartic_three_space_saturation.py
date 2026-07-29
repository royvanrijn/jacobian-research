#!/usr/bin/env python3
"""Exact quartic coordinate-three-space audit for the smooth cubic.

For every triple in a fixed basis of the 24-dimensional order-four tensor
kernel, test the full polynomial family

    phi_h + p0 * psi_i + p1 * psi_j + p2 * psi_k

over Q[p0,p1,p2,x,y,z].  The 2,024 coordinate three-spaces include every
specialization supported on the selected basis triple.  They do not test
directions supported on four or more basis elements.
"""

from __future__ import annotations

import itertools
import multiprocessing

import verify_cubic_symbol_double_saturation as cubic_audit


def audit_triple(
    triple: tuple[int, int, int],
) -> tuple[tuple[int, int, int], tuple[int, ...]]:
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    directions = cubic_audit.quartic_kernel_basis_tensors()
    result = cubic_audit.run_singular_subspace(
        cubic_audit.CUBIC_STRATA["smooth"],
        tuple(directions[index] for index in triple),
        timeout=300,
    )
    return triple, result


def main() -> None:
    triples = list(itertools.combinations(range(24), 3))
    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        for triple, result in pool.imap_unordered(audit_triple, triples):
            assert result == (0, 6, 0, 0, 3), (triple, result)

    print(
        "PASS: all 2024 quartic coordinate three-spaces for the smooth "
        "cubic have saturated cotangent presentation and a minimal "
        "relative length-six Ext presentation pulled back from the "
        "parameter-space origin"
    )


if __name__ == "__main__":
    main()
