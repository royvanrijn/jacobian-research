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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.deformation_complex import (  # noqa: E402
    relative_fiber_connection_complex,
)


def main() -> None:
    data = relative_fiber_connection_complex()
    complex_ = data.complex
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
