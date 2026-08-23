#!/usr/bin/env sage
"""Canonical recovery of the exact Q80 terminal horizontal.

This entry point implements the successful characteristic-zero construction:
recover easier high-incidence P.O=0 sections on the selected candidate1
parent, identify the correct pair by reduction to the transported historical
GF(73) P2-P3 marking, and take their exact elliptic-curve difference over
QQ(sqrt(-3))(W).

The direct three-node Hensel/resultant lift is superseded and is not the
construction used here.  The preserved implementation lives in the adjacent
`recover_q80_final_q6_via_basis_sections_impl.sage`; this wrapper is the
canonical status-bearing entry point.
"""
from pathlib import Path
HERE = Path(__file__).resolve().parent
load(str(HERE / "recover_q80_final_q6_via_basis_sections_impl.sage"))
