#!/usr/bin/env python3
"""Exact quartic coordinate-plane audit for singular squarefree cubics.

The ordinary rows test the pruned relative Ext presentation on every
coordinate plane.  Four planes have nonconstant ambient R-presentations;
their m^2-truncated parameter modules are tested instead by
Fitt_6=(1), Fitt_5=(0).
"""

from __future__ import annotations

import itertools
import multiprocessing

import verify_cubic_symbol_double_saturation as cubic_audit


SINGULAR_SQUAREFREE_STRATA = tuple(
    sorted(cubic_audit.SQUAREFREE_STRATA - {"smooth"})
)
FITTING_PAIRS = {
    ("cuspidal", (3, 8)),
    ("cuspidal", (4, 9)),
    ("line-tangent-conic", (2, 9)),
    ("nodal", (5, 11)),
}


def audit_plane(
    task: tuple[str, tuple[int, int]],
) -> tuple[tuple[str, tuple[int, int]], tuple[int, ...]]:
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    stratum, pair = task
    directions = cubic_audit.quartic_kernel_basis_tensors()
    result = cubic_audit.run_singular_plane(
        cubic_audit.CUBIC_STRATA[stratum],
        directions[pair[0]],
        directions[pair[1]],
        timeout=300,
    )
    return task, result


def main() -> None:
    pairs = list(itertools.combinations(range(24), 2))
    tasks = [
        (stratum, pair)
        for stratum in SINGULAR_SQUAREFREE_STRATA
        for pair in pairs
        if (stratum, pair) not in FITTING_PAIRS
    ]
    assert len(tasks) == 1652

    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        for task, result in pool.imap_unordered(audit_plane, tasks):
            assert result == (0, 6, 0, 0, 3), (task, result)

    directions = cubic_audit.quartic_kernel_basis_tensors()
    for stratum, pair in sorted(FITTING_PAIRS):
        fitting_result = cubic_audit.run_plane_parameter_fitting(
            cubic_audit.CUBIC_STRATA[stratum],
            directions[pair[0]],
            directions[pair[1]],
            timeout=300,
        )
        assert fitting_result == (0, 0, 0), (
            stratum,
            pair,
            fitting_result,
        )

    print(
        "PASS: 1652 ordinary singular-squarefree coordinate planes have "
        "the central pruned rank-three length-six Ext presentation"
    )
    print(
        "PASS: four exceptional-presentation planes are parameter-flat "
        "of rank six by Fitt_6=(1) and Fitt_5=(0)"
    )


if __name__ == "__main__":
    main()
