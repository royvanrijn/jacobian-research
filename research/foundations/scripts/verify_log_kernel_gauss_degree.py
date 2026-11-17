#!/usr/bin/env python3
"""Verify the Gauss-degree formula for a contracted cyclic log cokernel."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.log_node_profiles import (  # noqa: E402
    contracted_cyclic_cokernel_profile,
)


def basepoint_free_pencil_audit() -> None:
    # The column (s^d,t^d) embeds O(-d) in O^2 and has no common zero on P1.
    s, t = sp.symbols("s t")
    for degree in range(7):
        left = s**degree if degree else sp.Integer(1)
        right = t**degree if degree else sp.Integer(1)
        assert sp.gcd(left, right) == 1
        assert sp.Poly(left, s, t).total_degree() == degree
        assert sp.Poly(right, s, t).total_degree() == degree


def f2_root_audit() -> None:
    constant_direction = contracted_cyclic_cokernel_profile(54, 0)
    assert constant_direction.kernel_degree == 0
    assert constant_direction.cokernel_line_degree == 54
    assert constant_direction.cokernel_ch2 == 27

    exact_cancellation = contracted_cyclic_cokernel_profile(54, 27)
    assert exact_cancellation.kernel_degree == -27
    assert exact_cancellation.cokernel_line_degree == 27
    assert exact_cancellation.cokernel_ch2 == 0

    negative_contribution = contracted_cyclic_cokernel_profile(54, 28)
    assert negative_contribution.cokernel_ch2 == -1

    try:
        contracted_cyclic_cokernel_profile(54, -1)
    except ValueError:
        pass
    else:
        raise AssertionError("negative Gauss degree was accepted")


def main() -> None:
    basepoint_free_pencil_audit()
    f2_root_audit()
    print(
        "PASS: a contracted cyclic kernel is the tautological line of a "
        "Gauss map; the F2 root has ch2=27-deg(gamma_root)"
    )


if __name__ == "__main__":
    main()
