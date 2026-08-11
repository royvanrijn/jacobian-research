#!/usr/bin/env python3
"""Verify tangential-coordinate trivialization of the F2 root kernel."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.log_node_profiles import (  # noqa: E402
    tangential_kernel_trivialization_profile,
)


def logarithmic_divisibility_audit() -> None:
    # A unit is kept nonconstant to check full logarithmic divisibility, not
    # just the leading monomial.  For z=s*b, both logarithmic derivatives are
    # multiples of the full nonreduced Cartier equation s=W^3*U^18.
    W, U, a, c, d = sp.symbols("W U a c d")
    boundary_equation = W**3 * U**18
    unit = 1 + a * W + c * U + d * W * U
    tangential_pullback = sp.expand(boundary_equation * unit)
    log_w = sp.expand(W * sp.diff(tangential_pullback, W))
    log_u = sp.expand(U * sp.diff(tangential_pullback, U))
    assert sp.rem(log_w, boundary_equation, W) == 0
    assert sp.cancel(log_w / boundary_equation).is_polynomial(W, U)
    assert sp.cancel(log_u / boundary_equation).is_polynomial(W, U)

    # With T=W*U^5, the log matrix has first column (1,5) and its dz column
    # is divisible by s.  The residual determinant is the certified unit 3.
    residual = sp.cancel((log_u - 5 * log_w) / boundary_equation)
    assert sp.expand(residual.subs({W: 0, U: 0})) == 3


def f2_root_audit() -> None:
    profile = tangential_kernel_trivialization_profile(
        54,
        (3, 18),
        (3, 18),
    )
    assert profile.tangential_excess_orders == (0, 0)
    assert profile.kernel_degree == 0
    assert profile.gauss_degree == 0
    assert profile.cokernel_ch2 == 27

    # Extra tangential vanishing still gives the same fixed kernel direction.
    overvanishing = tangential_kernel_trivialization_profile(
        54,
        (3, 18),
        (4, 21),
    )
    assert overvanishing.tangential_excess_orders == (1, 3)
    assert overvanishing.gauss_degree == 0

    try:
        tangential_kernel_trivialization_profile(54, (3, 18), (2, 18))
    except ValueError:
        pass
    else:
        raise AssertionError("an insufficient tangential divisor was accepted")


def main() -> None:
    logarithmic_divisibility_audit()
    f2_root_audit()
    print(
        "PASS: the fixed target covector dz trivializes the full thickened "
        "F2 root kernel; e_root=0 and ch2(root)=27"
    )


if __name__ == "__main__":
    main()
