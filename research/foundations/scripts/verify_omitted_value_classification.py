#!/usr/bin/env python3
"""Exact regression checks for the omitted-value classifier."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.omitted import classify_omitted_values, multiplicity_partitions  # noqa: E402
from jcsearch.weighted import (  # noqa: E402
    WeightedSeedModel,
    canonical_seed,
    deformation_basis,
    w,
)


assert list(multiplicity_partitions(3)) == [(3,)]
assert list(multiplicity_partitions(4)) == [(4,), (2, 2)]
assert list(multiplicity_partitions(6)) == [(6,), (4, 2), (3, 3), (2, 2, 2)]

# Recover every previously proved low-degree exceptional value.
canonical_cubic = classify_omitted_values(w**2 * (1 - w), w)
assert [(item.slope, item.intercept) for item in canonical_cubic] == [
    (sp.Rational(1, 3), sp.Rational(1, 27))
]

canonical_quartic = classify_omitted_values(w**3 * (1 - w), w)
assert [(item.slope, item.intercept) for item in canonical_quartic] == [
    (sp.Rational(1, 8), -sp.Rational(1, 64))
]
assert classify_omitted_values(w**4 * (1 - w), w) == ()

exceptional_quartic = classify_omitted_values(w**2 * (1 - w**2), w)
assert [(item.slope, item.intercept) for item in exceptional_quartic] == [
    (0, -sp.Rational(1, 4))
]

# Reproduce the complete one-extra-root formula at four rational parameters.
for rho in (-2, -1, sp.Rational(1, 2), 3):
    H = sp.factor(w**2 * (1 - w) * (w - rho) / (1 - rho))
    omitted = classify_omitted_values(H, w)
    assert len(omitted) == 1
    assert sp.factor(omitted[0].slope - (1 - rho**2) / 8) == 0
    assert sp.factor(omitted[0].intercept + (1 - rho) ** 3 / 64) == 0

# Generic seeds with two and three simple extra roots are surjective.
generic_two = WeightedSeedModel(canonical_seed(2) + deformation_basis(2))
generic_three = WeightedSeedModel(canonical_seed(2) + deformation_basis(3))
assert generic_two.zero_profile()[2] != 1
assert generic_three.zero_profile()[2] != 1
assert classify_omitted_values(generic_two.primitive, w) == ()
assert classify_omitted_values(generic_three.primitive, w) == ()

# A boundary-clean two-extra-root seed on the exceptional coefficient locus.
# Its only omitted inverse polynomial has multiplicities (3,2).
special_H = sp.Rational(1, 27) * w**2 * (1 - w) * (w**2 + 8 * w + 18)
special_model = WeightedSeedModel(sp.diff(special_H, w))
special = classify_omitted_values(special_model.primitive, w)
assert len(special) == 1
assert (special[0].slope, special[0].intercept) == (-1, -1)
assert sp.expand(
    special_model.primitive + w - 1
    + sp.Rational(1, 27) * (w + 3) ** 3 * (w - 1) ** 2
) == 0

# Derive the exact coefficient-space condition for two extra roots.  Write
# their quadratic as W^2-RW+P and an omitted polynomial as
# -K(W-a)^3(W-b)^2.  The resultant is a nonzero scalar times F(R,P).
a, b, R, P = sp.symbols("a b R P")
b_from_top = (R + 1 - 3 * a) / 2
equation_2 = sp.factor(3 * a**2 + 6 * a * b_from_top + b_from_top**2 - P - R)
equation_3 = sp.factor(a**3 + 6 * a**2 * b_from_top + 3 * a * b_from_top**2 - P)
condition = sp.factor(sp.resultant(equation_2, equation_3, a))
F = (
    20 * P**3
    + 219 * P**2 * R**2
    - 312 * P**2 * R
    + 84 * P**2
    - 120 * P * R**4
    + 174 * P * R**3
    + 54 * P * R**2
    - 204 * P * R
    + 96 * P
    + 16 * R**6
    - 24 * R**5
    - 21 * R**4
    + 58 * R**3
    - 21 * R**2
    - 24 * R
    + 16
)
assert sp.factor(condition / F + sp.Rational(5, 64)) == 0
assert F.subs({R: -8, P: 18}) == 0

# Repeated extra roots are accepted by the same classifier; this admissible
# representative has no omitted values.
repeated_H = sp.Rational(1, 4) * w**2 * (1 - w) * (w + 1) ** 2
repeated_model = WeightedSeedModel(sp.diff(repeated_H, w))
assert classify_omitted_values(repeated_model.primitive, w) == ()

print("PASS: multiplicity partitions give a complete exact classifier")
print("PASS: all previously proved cubic and quartic omitted values are recovered")
print("PASS: generic two- and three-extra-root representatives are surjective")
print("PASS: F(R,P)=0 exactly detects the two-extra-root exceptional subfamily")
