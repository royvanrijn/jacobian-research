#!/usr/bin/env python3
"""Exact quartic tangent-direction audit for squarefree cubic symbols.

For every one of the seven squarefree ternary-cubic representatives and
every element of a primitive integral basis of the 24-dimensional
order-four tensor kernel, test the family phi_h+t*psi over Q[t,x,y,z].

The four invariant assertions are uniform cotangent saturation, relative
Ext^2 multiplicity six, no t-torsion, and radical support equal to the
collision axis.  Literal equality with the central Groebner presentation
is recorded separately because it depends on presentation coordinates.
"""

import verify_cubic_symbol_double_saturation as cubic_audit


# Sparse tangent tensors are faster to serialize expanded than to refactor
# separately in each of the 168 exact family constructions.
cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
CUBIC_STRATA = cubic_audit.CUBIC_STRATA
SQUAREFREE_STRATA = cubic_audit.SQUAREFREE_STRATA
quartic_kernel_basis_tensors = cubic_audit.quartic_kernel_basis_tensors
run_singular_family = cubic_audit.run_singular_family


EXPECTED_PRESENTATION_CHANGES = {
    "concurrent-lines": [(13, 4), (17, 4)],
    "cuspidal": [],
    "line-tangent-conic": [(10, 2)],
    "line-transverse-conic": [],
    "nodal": [],
    "smooth": [],
    "triangle": [(17, 2)],
}


directions = quartic_kernel_basis_tensors()
assert len(directions) == 24

for stratum in sorted(SQUAREFREE_STRATA):
    presentation_changes = []
    for direction_index, direction in enumerate(directions):
        result = run_singular_family(CUBIC_STRATA[stratum], direction)
        assert result[:4] == (0, 6, 0, 0)
        if result[4]:
            presentation_changes.append((direction_index, result[4]))
    assert presentation_changes == EXPECTED_PRESENTATION_CHANGES[stratum]
    print(
        f"PASS: {stratum}: all 24 quartic basis directions preserve "
        "cotangent saturation and a flat relative length-six support defect; "
        f"presentation changes={presentation_changes}"
    )

print(
    "PASS: all 168 squarefree symbol/quartic tangent families preserve "
    "the invariant double-saturation defect"
)
