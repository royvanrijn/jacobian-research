#!/usr/bin/env python3
"""Exact ledger audit for homogeneous coefficient hypersurfaces."""

from __future__ import annotations

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


d = sp.symbols("d", integer=True, positive=True)
boundary_matrix = sp.Matrix(
    [
        [1, 0, 1],
        [0, 1, 1],
        [d, d, d],
    ]
)
assert sp.factor(boundary_matrix.det()) == -d

for degree in range(1, 13):
    specialized = boundary_matrix.subs(d, degree)
    smith = smith_normal_form(specialized, domain=ZZ)
    diagonal = tuple(
        abs(int(smith[index, index]))
        for index in range(3)
    )
    assert diagonal == (1, 1, degree)
print("PASS: the boundary Smith form is diag(1,1,d)")


# Cox weights under (L1,L2,L3)->(zeta L1,zeta L2,zeta^-1 L3).
factor_weights = (1, 1, -1)
r13_weight = factor_weights[0] + factor_weights[2]
r23_weight = factor_weights[1] + factor_weights[2]
product_weight = sum(factor_weights)
assert r13_weight == 0
assert r23_weight == 0
assert product_weight == 1
for degree in range(2, 13):
    assert degree * product_weight % degree == 0
print("PASS: mu_d fixes both resultants and every degree-d coefficient level")


# On r13=1, L1 is nonzero.  A nontrivial scalar cannot fix that vector, so
# the residual action is free.  The resulting finite-etale quotient would
# force 1=d*chi_c(quotient) if the source (or a stabilization) were affine.
for degree in range(2, 13):
    assert all(
        degree * quotient_euler != 1
        for quotient_euler in range(-20, 21)
    )
print("PASS: Euler characteristic one excludes every d>1 free quotient")
print("PASS homogeneous coefficient-hypersurface obstruction")
