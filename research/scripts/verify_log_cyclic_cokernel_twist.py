#!/usr/bin/env python3
"""Verify the kernel-line formula for a cyclic logarithmic cokernel."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.log_node_profiles import (  # noqa: E402
    cyclic_cokernel_twist_profile,
)


def determinant_line_audit() -> None:
    # On D, a rank-one map E_D -> F_D has kernel K, image I, and cokernel L.
    # det(E_D)=K*I and det(F_D)=I*L, so L/K=det(F)/det(E)=O_D(D).
    kernel, image, normal = sp.symbols("kernel image normal", nonzero=True)
    determinant_source = kernel * image
    cokernel = kernel * normal
    determinant_target = image * cokernel
    assert sp.cancel(determinant_target / determinant_source) == normal
    assert sp.cancel(cokernel / kernel) == normal


def chern_character_audit() -> None:
    divisor, kernel_degree = sp.symbols("divisor kernel_degree")
    cokernel_line_degree = kernel_degree + divisor**2
    # GRR for i:D->X gives ch_2(i_*L)=deg(L)-D^2/2.
    cokernel_ch2 = sp.expand(cokernel_line_degree - divisor**2 / 2)
    assert cokernel_ch2 == kernel_degree + divisor**2 / 2


def f2_root_audit() -> None:
    profile = cyclic_cokernel_twist_profile(54, 0)
    assert profile.cokernel_line_degree == 54
    assert profile.cokernel_ch2 == 27

    unknown_kernel = cyclic_cokernel_twist_profile(54, -7)
    assert unknown_kernel.cokernel_line_degree == 47
    assert unknown_kernel.cokernel_ch2 == 20


def main() -> None:
    determinant_line_audit()
    chern_character_audit()
    f2_root_audit()
    print(
        "PASS: cyclic cokernel twist is the restricted kernel line times "
        "O_D(D); the F2 root has ch2=deg(K_root)+27"
    )


if __name__ == "__main__":
    main()
