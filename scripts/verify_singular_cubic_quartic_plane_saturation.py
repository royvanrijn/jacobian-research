#!/usr/bin/env python3
"""Exact quartic coordinate-plane audit for singular squarefree cubics.

The ordinary rows test the pruned relative Ext presentation on every
coordinate plane.  One nodal plane, (5,11), has a nonconstant ambient
R-presentation; its m^2-truncated parameter module is tested instead by
Fitt_6=(1), Fitt_5=(0).  The nodal plane (21,22) is explicitly unresolved
because exact Ext elimination exceeds the bounded reproduction time.
"""

from __future__ import annotations

import itertools
import multiprocessing

import verify_cubic_symbol_double_saturation as cubic_audit


SINGULAR_SQUAREFREE_STRATA = tuple(
    sorted(cubic_audit.SQUAREFREE_STRATA - {"smooth"})
)
NODAL_FITTING_PAIR = (5, 11)
NODAL_UNRESOLVED_PAIR = (21, 22)


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
        if not (
            stratum == "nodal"
            and pair in {NODAL_FITTING_PAIR, NODAL_UNRESOLVED_PAIR}
        )
    ]
    assert len(tasks) == 1654

    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        for task, result in pool.imap_unordered(audit_plane, tasks):
            assert result == (0, 6, 0, 0, 3), (task, result)

    directions = cubic_audit.quartic_kernel_basis_tensors()
    fitting_result = cubic_audit.run_plane_parameter_fitting(
        cubic_audit.CUBIC_STRATA["nodal"],
        directions[NODAL_FITTING_PAIR[0]],
        directions[NODAL_FITTING_PAIR[1]],
        timeout=300,
    )
    assert fitting_result == (0, 0, 0)

    print(
        "PASS: 1654 ordinary singular-squarefree coordinate planes have "
        "the central pruned rank-three length-six Ext presentation"
    )
    print(
        "PASS: nodal plane (5,11) is parameter-flat of rank six by "
        "Fitt_6=(1) and Fitt_5=(0)"
    )
    print(
        "UNRESOLVED: nodal plane (21,22) exceeds the bounded exact Ext "
        "elimination time"
    )


if __name__ == "__main__":
    main()
