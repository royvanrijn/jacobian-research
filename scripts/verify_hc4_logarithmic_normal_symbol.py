#!/usr/bin/env python3
"""Regression for the uniform HC4 caustic normal-symbol matrix.

At pure-normal tensor order n, coefficient extraction from

    T((alpha, 4*(2*W+9*b**2)/3), ..., (alpha, 4*(2*W+9*b**2)/3))/2

gives the matrix M_n in Section 9.7 of
HC4_PC2_GRAPH_POLARIZATION_AUDIT.md.  The proof that its determinant is a
nonzero scalar is the triangular calculation in that note.  This script
checks the formula exactly through the stated finite regression range in
both X-caustic charts.
"""

from __future__ import annotations

import sympy as sp


W, b = sp.symbols("W b")


def symbol_matrix(n: int, alpha: int) -> sp.Matrix:
    """Return rows [W^j] and columns T(c^(n-i),d^i)."""
    columns = []
    tau = 2 * W + 9 * b**2
    for i in range(n + 1):
        polynomial = (
            sp.binomial(n, i)
            * alpha ** (n - i)
            * (sp.Rational(4, 3) * tau) ** i
            / 2
        )
        columns.append([sp.expand(polynomial).coeff(W, j) for j in range(n + 1)])
    return sp.Matrix(n + 1, n + 1, lambda j, i: columns[i][j])


def determinant_formula(n: int, alpha: int) -> sp.Expr:
    return (
        sp.Rational(1, 2) ** (n + 1)
        * sp.prod(sp.binomial(n, i) for i in range(n + 1))
        * alpha ** (n * (n + 1) // 2)
        * sp.Rational(8, 3) ** (n * (n + 1) // 2)
    )


def main() -> None:
    for alpha in (1, 2):
        for n in range(2, 11):
            matrix = symbol_matrix(n, alpha)
            assert sp.factor(matrix.det() - determinant_formula(n, alpha)) == 0
            assert not matrix.det().has(b)
    print("PASS: HC4 normal-symbol matrices are unimodular through order 10")
    print("PASS: both X-caustic chart determinants match the uniform formula")
    print("SCOPE: the all-order proof is the triangular calculation in Section 9.7")


if __name__ == "__main__":
    main()
