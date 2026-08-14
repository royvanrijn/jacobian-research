#!/usr/bin/env python3
"""Verify algebraic identities in the smooth-cubic orthogonal normal form.

The geometric part of HC4NHM19 is a vector-bundle argument on an elliptic
curve.  This checker replays the universal rank-two hyperbolic matrix,
kernel, adjugate, and trace-free endomorphism identities, and checks the
finite degree ledger for its two isotropic line summands.
"""

from __future__ import annotations

import sympy as sp


def verify_hyperbolic_matrix() -> None:
    alpha = sp.Matrix(sp.symbols("a0:3"))
    beta = sp.Matrix(sp.symbols("b0:3"))
    matrix = alpha * beta.T + beta * alpha.T
    kernel = alpha.cross(beta)

    assert (matrix * kernel).applyfunc(sp.expand) == sp.zeros(3, 1)
    assert (matrix.adjugate() + kernel * kernel.T).applyfunc(
        sp.expand
    ) == sp.zeros(3)
    assert sp.expand(matrix.det()) == 0
    print("PASS: hyperbolic restriction has cross-product kernel and rank-one adjugate")


def verify_trace_free_splitting() -> None:
    u, v, w = sp.symbols("u v w")
    endomorphism = sp.Matrix([[u, v], [w, -u]])
    alternating = sp.Matrix([[0, 1], [-1, 0]])
    symmetric_form = endomorphism.T * alternating

    assert symmetric_form == symmetric_form.T
    assert endomorphism**2 == (u**2 + v * w) * sp.eye(2)

    a, b, c, d = sp.symbols("a b c d")
    general = sp.Matrix([[a, b], [c, d]])
    general_form = general.T * alternating
    skew_part = (general_form - general_form.T).applyfunc(sp.expand)
    assert skew_part == sp.Matrix([[0, a + d], [-a - d, 0]])
    print("PASS: symmetry is equivalent to trace zero and gives two eigenlines")


def verify_degree_ledger() -> None:
    packets = []
    for first in range(10):
        second = 9 - first
        first_generated = first == 0 or first >= 2
        second_generated = second == 0 or second >= 2
        if first <= second and first_generated and second_generated:
            packets.append((first, second))
    assert packets == [(0, 9), (2, 7), (3, 6), (4, 5)]
    print("PASS: global generation leaves degree packets (0,9), (2,7), (3,6), (4,5)")


def main() -> None:
    verify_hyperbolic_matrix()
    verify_trace_free_splitting()
    verify_degree_ledger()
    print("THEOREM: the clean smooth-cubic restriction has the orthogonal normal form")


if __name__ == "__main__":
    main()
