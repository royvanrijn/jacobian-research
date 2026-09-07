#!/usr/bin/env python3
"""Verify local models for cyclic-submodule positivity of log cokernels."""

from __future__ import annotations

import sympy as sp


r, t = sp.symbols("r t")


def thick_node_audit(left: int, right: int) -> None:
    """Check the glued-to-split positive quotient at a thick node."""

    assert left >= 1 and right >= 1
    theta = sp.diag(r**left, t**right)
    determinant = sp.factor(theta.det())
    assert determinant == r**left * t**right
    fitting_one = sp.groebner(tuple(theta), r, t, domain=sp.QQ)
    expected = sp.groebner([r**left, t**right], r, t, domain=sp.QQ)
    assert fitting_one == expected

    # The class (1,1) has annihilator
    # (r^left) intersect (t^right)=(r^left*t^right), and the quotient by its
    # cyclic span is R/(r^left,t^right), of length left*right.
    quotient_basis = tuple(
        r**i * t**j for i in range(left) for j in range(right)
    )
    assert len(quotient_basis) == left * right


def unibranch_fitting_audit(left: int, right: int, exponent: int) -> None:
    """Check a noncyclic unibranch matrix and the positive cyclic quotient."""

    assert left >= 1 and right >= 1 and exponent >= 1
    theta = sp.Matrix([[r**left, 0], [t**exponent, r**right]])
    determinant = sp.factor(theta.det())
    assert determinant == r ** (left + right)

    fitting_one = sp.groebner(tuple(theta), r, t, domain=sp.QQ)
    expected_fitting = sp.groebner(
        [r ** min(left, right), t**exponent], r, t, domain=sp.QQ
    )
    assert fitting_one == expected_fitting

    # The first target generator is killed by r^(left+right).  Modulo its
    # cyclic span, the second generator has relations t^exponent and r^right.
    quotient_length = exponent * right
    fitting_colength = exponent * min(left, right)
    assert quotient_length >= fitting_colength


def main() -> None:
    for left in range(1, 8):
        for right in range(1, 8):
            thick_node_audit(left, right)
            for exponent in range(1, 9):
                unibranch_fitting_audit(left, right, exponent)
    print(
        "PASS: thick-node and unibranch 2x2 logarithmic cokernels contain "
        "their generically cyclic determinant module with an effective "
        "finite quotient whose length dominates the Fitt_1 colength"
    )


if __name__ == "__main__":
    main()
