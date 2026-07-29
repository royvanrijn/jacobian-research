#!/usr/bin/env python3
"""Exact full-support quartic-plane audit for squarefree cubic symbols.

Let psi_0,...,psi_23 be the fixed primitive basis of the order-four
generalized-triple-cover tensor kernel.  Define

    psi_plus  = sum_i psi_i,
    psi_minus = sum_i (-1)^i psi_i.

For each squarefree ternary-cubic symbol h, this checker audits the complete
polynomial family

    phi_h + u*psi_plus + v*psi_minus

over Q[u,v,x,y,z].  On u^2-v^2 != 0 every one of the 24 basis coordinates
is nonzero, so the generic tensor on this plane has full basis support.
"""

from __future__ import annotations

import multiprocessing
from functools import cache

import sympy as sp

import verify_cubic_symbol_double_saturation as cubic_audit


@cache
def dense_directions() -> tuple[
    dict[tuple[int, int, int], sp.Expr],
    dict[tuple[int, int, int], sp.Expr],
]:
    """Return the low-height full-support sum and alternating-sum tensors."""

    basis = cubic_audit.quartic_kernel_basis_tensors()
    triples = tuple(basis[0])
    plus = {
        triple: sp.expand(sum(tensor[triple] for tensor in basis))
        for triple in triples
    }
    minus = {
        triple: sp.expand(
            sum(
                (-1) ** index * tensor[triple]
                for index, tensor in enumerate(basis)
            )
        )
        for triple in triples
    }

    # The coefficient rows have rank two and every coordinate is nonzero.
    coefficient_matrix = sp.Matrix(
        (
            (sp.Integer(1),) * len(basis),
            tuple(sp.Integer((-1) ** index) for index in range(len(basis))),
        )
    )
    assert coefficient_matrix.rank() == 2
    assert all(entry != 0 for entry in coefficient_matrix)
    return plus, minus


def audit_stratum(name: str) -> tuple[str, tuple[int, ...]]:
    """Run the exact polynomial-plane audit for one cubic-symbol stratum."""

    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    result = cubic_audit.run_singular_subspace(
        cubic_audit.CUBIC_STRATA[name],
        dense_directions(),
        timeout=600,
    )
    return name, result


def main() -> None:
    names = sorted(cubic_audit.SQUAREFREE_STRATA)
    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        for name, result in pool.imap_unordered(audit_stratum, names):
            assert result == (0, 6, 0, 0, 3), (name, result)
            print(
                f"PASS: {name}: dense quartic plane has saturated "
                "cotangent presentation and central rank-three Ext block"
            )

    print(
        "PASS: all seven squarefree symbols retain the flat relative "
        "length-six support defect on one full-support quartic plane"
    )
    print(
        "PASS: u^2-v^2 != 0 gives nonzero coordinates in all 24 fixed "
        "quartic-kernel basis directions"
    )


if __name__ == "__main__":
    main()
