#!/usr/bin/env sage
"""Canonical exact Q80 final-q6 resolved Riemann--Roch certificate.

Input is the exact horizontal recovered by
`recover_q80_final_q6_via_basis_sections.sage`.  This entry point verifies the
transported modular gauge only as a regression, derives the exact whole-A4
and connected-A5 quotient rows, and certifies ambient dimension 4, condition
rank 2, kernel dimension 2 and h0(D)=2 over QQ(sqrt(-3)).

The old direct Hensel/resultant section lift is superseded.  The preserved
working implementation lives in the adjacent
`certify_q80_final_q6_char0_rr_from_basis_impl.sage`; this wrapper is the
canonical status-bearing entry point.
"""
from pathlib import Path
HERE = Path(__file__).resolve().parent
load(str(HERE / "certify_q80_final_q6_char0_rr_from_basis_impl.sage"))
