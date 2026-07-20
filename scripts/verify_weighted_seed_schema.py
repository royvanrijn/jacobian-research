#!/usr/bin/env python3
"""Exact regression checks for the general weighted-seed model."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, canonical_seed, w


A, B, C, s, t = sp.symbols("A B C s t")

for degree in range(2, 6):
    model = WeightedSeedModel(canonical_seed(degree))
    assert sp.factor(model.primitive - w**degree*(1-w)) == 0
    assert model.seed_degree == degree
    assert model.fiber_degree == degree + 1
    m0, m1, extra = model.zero_profile()
    assert (m0, m1, extra) == (degree, 1, -1)

    inverse = model.inverse_polynomial(A, B, C)
    assert sp.factor(sp.diff(inverse, w) - (model.seed-B*C)) == 0
    branch_s, branch_t = model.branch_parameterization()
    pencil = model.primitive - s*w + t
    assert sp.factor(pencil.subs({s: branch_s, t: branch_t})) == 0
    assert sp.factor(sp.diff(pencil, w).subs(s, branch_s)) == 0

    boundary = model.boundary_map()
    expected_y_coefficient = -model.c/(model.kappa+2)
    expected_z_coefficient = model.b*(model.kappa+2)
    assert sp.factor(boundary[1] - expected_y_coefficient*sp.Symbol("y")) == 0
    assert sp.factor(sp.diff(boundary[0], sp.Symbol("z")) - expected_z_coefficient) == 0
    assert boundary[2] == 0

# A noncanonical quartic-sheet seed has an additional primitive zero, which is
# precisely the C=0 boundary complication detected by the scan.
alternative = WeightedSeedModel(w-2*w**3)
m0, m1, extra = alternative.zero_profile()
assert (m0, m1) == (2, 1)
assert sp.factor(extra + (w+1)/2) == 0

print("PASS: canonical H_d=w^d(1-w) models through seed degree five")
print("PASS: inverse pencil and discriminant parameterization identities")
print("PASS: every canonical x=0 boundary map is triangular and invertible")
print("PASS: additional primitive zeros are exposed by the seed profile")
