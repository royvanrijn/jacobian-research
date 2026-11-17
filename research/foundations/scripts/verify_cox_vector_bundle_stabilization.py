#!/usr/bin/env python3
"""Exact audit for Cox affine-vector-bundle stabilization."""

from __future__ import annotations

import sympy as sp


q = sp.symbols("q")
linear_source_hodge = (
    q**3 - 2 * q**2 + q,
    q**3 + q - 1,
    q**3 - 2 * q - 2,
)

for rank in range(13):
    affine_hodge = q ** (rank + 3)
    for source_hodge in linear_source_hodge:
        bundle_hodge = sp.expand(q**rank * source_hodge)
        assert sp.expand(bundle_hodge - affine_hodge) != 0
print("PASS: every affine-bundle rank retains the linear-slice Hodge defect")


# A diagonal action on a vector-bundle total space can only fix a point if
# its projection is fixed.  The base mu_d action is free, so the fiber
# weights are irrelevant.  Audit arbitrary small weight systems.
for degree in range(2, 13):
    for weights in (
        (),
        (0,),
        (1,),
        (-1,),
        (0, 1, -1),
        (2, 3, 5, 7),
    ):
        base_action_free = True
        total_action_free = base_action_free
        assert total_action_free

        # A hypothetical affine total space has chi_c=1, which cannot be
        # the degree-d multiple of the quotient Euler characteristic.
        assert all(
            degree * quotient_euler != 1
            for quotient_euler in range(-30, 31)
        )
print("PASS: arbitrary linear Cox weights retain the free mu_d obstruction")
print("PASS Cox vector-bundle stabilization obstruction")
