#!/usr/bin/env python3
"""Exact finite audit for the first rank-ten parallel-class obstructions.

The note's general arguments use one-relation rigidity, Witt index, Gale-dual
rank, the classical ternary Hesse theorem, and tensor independence.  This
checker independently verifies the finite polynomial identities underneath
those arguments:

* diagonal and square-pair rigidity on a four-term full-support hyperplane;
* factorisation and nonvanishing of a 3+3 block Hessian determinant;
* independence of the nine cross square-pair tensors for two triple classes.

It does not replace the matroid or Hesse-theorem steps in the written proof.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


def coefficient_rank(polys: list[sp.Expr], variables: tuple[sp.Symbol, ...]) -> int:
    monomials = sorted(
        set().union(*(sp.Poly(poly, *variables).monoms() for poly in polys))
    )
    matrix = sp.Matrix(
        [
            [sp.Poly(poly, *variables).coeff_monomial(monomial) for poly in polys]
            for monomial in monomials
        ]
    )
    return matrix.rank()


def four_term_hyperplane_rigidity() -> None:
    x, y, z = sp.symbols("x y z")
    forms = (x, y, z, -x - y - z)
    diagonal = [form**2 for form in forms]
    square_pairs = [forms[i] ** 2 * forms[j] ** 2 for i, j in combinations(range(4), 2)]
    dr = coefficient_rank(diagonal, (x, y, z))
    pr = coefficient_rank(square_pairs, (x, y, z))
    if dr != 4 or pr != 6:
        raise AssertionError((dr, pr))
    print("QUARTIC_HN_RANK10_CLASS4_DIAGONAL_RANK", dr)
    print("QUARTIC_HN_RANK10_CLASS4_SQUARE_PAIR_RANK", pr)


def block_hessian_factorisation() -> None:
    # Two essential ternary Waring blocks: four terms with one relation and
    # six terms with three relations.  The determinant factorisation is
    # coordinate-free; this exact instance audits the Cauchy--Binet algebra.
    A = sp.Matrix(
        [
            [1, 0, 0, 1],
            [0, 1, 0, 1],
            [0, 0, 1, 1],
        ]
    )
    B = sp.Matrix(
        [
            [1, 0, 0, 1, 1, 2],
            [0, 1, 0, 1, 2, 3],
            [0, 0, 1, 1, 3, 5],
        ]
    )
    a = sp.symbols("a0:4")
    b = sp.symbols("b0:6")
    HA = A * sp.diag(*(value**2 for value in a)) * A.T
    HB = B * sp.diag(*(value**2 for value in b)) * B.T
    H = sp.diag(1, 1, 1, 1, 1, 1)
    H[:3, :3] = HA
    H[3:, 3:] = HB
    actual = sp.factor(H.det())
    expected = sp.factor(HA.det() * HB.det())
    if sp.expand(actual - expected) != 0:
        raise AssertionError((actual, expected))

    # A full-support point on the one-relation A-value hyperplane, and one
    # ordinary point for the B block.
    aval = {a[0]: 1, a[1]: 2, a[2]: -4, a[3]: 1}
    bval = {b[i]: i + 1 for i in range(6)}
    left = sp.factor(HA.det().subs(aval))
    right = sp.factor(HB.det().subs(bval))
    if not left or not right:
        raise AssertionError((left, right))
    print("QUARTIC_HN_RANK10_CLASS4_BLOCK_FACTOR_PASS", left, right)


def two_triple_tensor_independence() -> None:
    x, y, u, v = sp.symbols("x y u v")
    left = (x, y, -x - y)
    right = (u, v, -u - v)
    tensors = [a**2 * b**2 for a in left for b in right]
    rank = coefficient_rank(tensors, (x, y, u, v))
    if rank != 9:
        raise AssertionError(rank)
    print("QUARTIC_HN_RANK10_TWO_TRIPLE_TENSOR_RANK", rank)


def main() -> None:
    four_term_hyperplane_rigidity()
    block_hessian_factorisation()
    two_triple_tensor_independence()
    print("QUARTIC_HN_WARING_RANK10_PARALLEL_GATE_PASS")


if __name__ == "__main__":
    main()
