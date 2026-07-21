#!/usr/bin/env python3
"""Exact regression checks for the general weighted-seed model."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, canonical_seed, w


A, B, C, s, t, alpha = sp.symbols("A B C s t alpha")

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

# The split degreewise seed H_N has an exact C^2 discriminant factor.  These
# bounded symbolic cases regress the all-degree root-valuation proof in the
# canonical multiplicity paper.
for fiber_degree in range(4, 9):
    split_primitive = w**2 * (1-w) * (1+w**(fiber_degree-3))
    split_inverse = split_primitive - B*C*w + 2*A*C**2

    # The paper's explicit target is (A,B,C)=(pi,0,1).  Algebraically, a
    # transcendental alpha cannot be a critical value of the rational seed:
    # H_N+2*alpha and H_N' are coprime over Q(alpha).  Since C=1, these are
    # exactly the simultaneous simple-root and reconstruction-denominator
    # conditions.
    explicit_fiber = split_primitive + 2*alpha
    assert sp.gcd(
        sp.Poly(explicit_fiber, w, domain=sp.QQ.frac_field(alpha)),
        sp.Poly(sp.diff(split_primitive, w), w, domain=sp.QQ.frac_field(alpha)),
    ).degree() == 0

    raw_discriminant = sp.factor(sp.discriminant(split_inverse, w))
    c_terms = sp.Poly(raw_discriminant, C).terms()
    assert min(monomial[0] for monomial, coefficient in c_terms if coefficient) == 2
    residual_discriminant = sp.cancel(raw_discriminant/C**2)
    residual_at_c0 = sp.factor(residual_discriminant.subs(C, 0))
    assert sp.rem(residual_at_c0, B**2-8*A, A) == 0
    assert sp.factor(residual_at_c0/(B**2-8*A)).is_number
    assert sp.Poly(residual_discriminant, A, B, C).is_irreducible

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
print("PASS: split H_N discriminants have exact C^2 saturation and prime residuals through degree eight")
print("PASS: the explicit target (alpha,0,1) has only simple, regularly reconstructing roots")
