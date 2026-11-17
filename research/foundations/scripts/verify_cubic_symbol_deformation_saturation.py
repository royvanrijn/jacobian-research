#!/usr/bin/env python3
"""Exact one-parameter deformation audit for squarefree cubic symbols."""

from verify_cubic_symbol_double_saturation import (
    CUBIC_STRATA,
    SQUAREFREE_STRATA,
    run_singular_family,
)


for stratum in sorted(SQUAREFREE_STRATA):
    result = run_singular_family(CUBIC_STRATA[stratum])
    assert result == (0, 6, 0, 0, 0)
    print(
        f"PASS: {stratum}: uniform cotangent saturation=0; "
        "relative Ext2 has a parameter-independent presentation on the "
        "collision axis and multiplicity 6"
    )

print(
    "PASS: every squarefree stratum retains a flat relative length-six "
    "support defect along phi_h+t*psi_4"
)
