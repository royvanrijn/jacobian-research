#!/usr/bin/env python3
"""Relative rank-one fiber and Hamiltonian-connection deformation complex.

In canonical fiber coordinates (S,T) over the central coordinate R, an
order-m correction is a triple (s,t,a): two fiber-symbol corrections and a
connection-Hamiltonian correction.  Its three raw defects are

    F = s_S + t_T,
    G = s_R - a_T,
    H = t_R + a_S.

They satisfy the Bianchi identity F_R-G_S-H_T=0.  This script constructs the
closed-defect module exactly, verifies that Hamiltonian gauges map to zero,
and proves the bounded polynomial complex is exact in degrees H1 and H2.
The matrices have integer entries, so the result base-changes to every
characteristic-zero seed-parameter algebra.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.deformation_complex import (  # noqa: E402
    ThreeTermComplex,
    polynomial_action_matrix,
    polynomial_coordinate_column,
    total_degree_monomials,
)


S, T, R = sp.symbols("S T R")
VARIABLES = (S, T, R)


def block_column(*columns: sp.Matrix) -> sp.Matrix:
    return sp.Matrix.vstack(*columns)


def main() -> None:
    gauge_basis = total_degree_monomials(VARIABLES, 4)
    scalar_corrections = total_degree_monomials(VARIABLES, 3)
    scalar_defects = total_degree_monomials(VARIABLES, 2)
    bianchi_targets = total_degree_monomials(VARIABLES, 1)

    gauge_columns = []
    for hamiltonian in gauge_basis:
        gauge_columns.append(
            block_column(
                polynomial_coordinate_column(
                    -sp.diff(hamiltonian, T),
                    scalar_corrections,
                    VARIABLES,
                ),
                polynomial_coordinate_column(
                    sp.diff(hamiltonian, S),
                    scalar_corrections,
                    VARIABLES,
                ),
                polynomial_coordinate_column(
                    -sp.diff(hamiltonian, R),
                    scalar_corrections,
                    VARIABLES,
                ),
            )
        )
    d0 = sp.Matrix.hstack(*gauge_columns)

    raw_defect_columns = []
    for component in range(3):
        for correction in scalar_corrections:
            if component == 0:
                defects = (
                    sp.diff(correction, S),
                    sp.diff(correction, R),
                    sp.Integer(0),
                )
            elif component == 1:
                defects = (
                    sp.diff(correction, T),
                    sp.Integer(0),
                    sp.diff(correction, R),
                )
            else:
                defects = (
                    sp.Integer(0),
                    -sp.diff(correction, T),
                    sp.diff(correction, S),
                )
            raw_defect_columns.append(
                block_column(
                    *(
                        polynomial_coordinate_column(
                            defect,
                            scalar_defects,
                            VARIABLES,
                        )
                        for defect in defects
                    )
                )
            )
    raw_d1 = sp.Matrix.hstack(*raw_defect_columns)

    # The Bianchi map sends (F,G,H) to F_R-G_S-H_T.
    bianchi_blocks = (
        polynomial_action_matrix(
            scalar_defects,
            bianchi_targets,
            VARIABLES,
            lambda polynomial: sp.diff(polynomial, R),
        ),
        polynomial_action_matrix(
            scalar_defects,
            bianchi_targets,
            VARIABLES,
            lambda polynomial: -sp.diff(polynomial, S),
        ),
        polynomial_action_matrix(
            scalar_defects,
            bianchi_targets,
            VARIABLES,
            lambda polynomial: -sp.diff(polynomial, T),
        ),
    )
    bianchi = sp.Matrix.hstack(*bianchi_blocks)
    assert bianchi * raw_d1 == sp.zeros(bianchi.rows, raw_d1.cols)

    # Use a basis of closed defects as C2.  Since the basis matrix has full
    # column rank, the coordinates of every raw d1 column are unique.
    closed_columns = bianchi.nullspace()
    closed_basis = sp.Matrix.hstack(*closed_columns)
    d1_closed, parameters = closed_basis.gauss_jordan_solve(raw_d1)
    assert not list(parameters)
    assert closed_basis * d1_closed == raw_d1

    complex_ = ThreeTermComplex(
        d0,
        d1_closed,
        "relative fiber-plus-connection complex",
    )
    assert complex_.dimensions == (35, 60, 26)
    assert complex_.ranks == (34, 26)
    assert complex_.cohomology_dimensions == (1, 0, 0)
    assert complex_.prime_rank_profile((31991, 32003, 65521)) == {
        prime: complex_.ranks for prime in (31991, 32003, 65521)
    }

    print("PASS: raw fiber/connection defects satisfy the Bianchi identity")
    print("PASS: closed defect space has dimension 26")
    print("PASS: relative complex dimensions are 35 -> 60 -> 26")
    print("PASS: H0,H1,H2 dimensions are (1,0,0)")
    print("PASS: ranks agree over Q and three good primes")
    print("SCOPE: bounded canonical/formal chart; no boundary-lattice descent")


if __name__ == "__main__":
    main()
