#!/usr/bin/env python3
"""Exact audit of boundary saturation with repeated extra primitive roots."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.boundary import boundary_saturation_profile  # noqa: E402
from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


A, B, C = sp.symbols("A B C")


def normalized_primitive(zero_multiplicity, extras):
    """Return H with H'(1)=-1 for ``extras=[(root,multiplicity),...]``."""
    normalization = sp.prod((1 - root) ** multiplicity for root, multiplicity in extras)
    return sp.factor(
        w**zero_multiplicity
        * (1 - w)
        * sp.prod((w - root) ** multiplicity for root, multiplicity in extras)
        / normalization
    )


cases = (
    (2, ((-1, 2),)),
    (2, ((-1, 3),)),
    (3, ((-1, 2),)),
    (4, ((-1, 3),)),
    (2, ((-1, 2), (-2, 3))),
)

for zero_multiplicity, extras in cases:
    H = normalized_primitive(zero_multiplicity, extras)
    profile = boundary_saturation_profile(H, w)
    repetition_excess = sum(multiplicity - 1 for _, multiplicity in extras)
    assert profile.zero_multiplicity == zero_multiplicity
    assert profile.one_multiplicity == 1
    assert profile.extra_repetition_excess == repetition_excess
    assert profile.saturation_exponent == zero_multiplicity + repetition_excess

    # Directly compute the pulled-back discriminant and its first nonzero
    # C-coefficient.  It must agree with the predicted trace up to a unit.
    inverse = H - B * C * w + A * C**2
    discriminant = sp.Poly(sp.expand(sp.discriminant(inverse, w)), C)
    actual_exponent = min(monomial[0] for monomial, _ in discriminant.terms())
    assert actual_exponent == profile.saturation_exponent
    boundary_coefficient = sp.factor(
        sum(
            coefficient
            for (power,), coefficient in discriminant.terms()
            if power == actual_exponent
        )
    )
    predicted_trace = profile.boundary_trace(A, B)
    ratio = sp.factor(boundary_coefficient / predicted_trace)
    assert ratio != 0 and not ratio.has(A, B, C)

# Repeated nonzero roots do not create new finite C=0 charts.  Audit the
# actual polynomial map once in each zero-multiplicity regime.
for zero_multiplicity, extras in (
    (2, ((-1, 2),)),
    (3, ((-1, 2),)),
):
    H = normalized_primitive(zero_multiplicity, extras)
    model = WeightedSeedModel(sp.diff(H, w))
    gamma_zero_z = -(1 + model.a * x * y) / x**2
    gamma_image = tuple(sp.factor(component.subs(z, gamma_zero_z)) for component in model.mapping())
    u = 1 + x * y
    if zero_multiplicity == 2:
        h2 = sp.Poly(H, w).coeff_monomial(w**2)
        expected = (
            sp.factor((u + h2 * u**2) / x**2),
            sp.factor((1 + 2 * h2 * u) / x),
            0,
        )
    else:
        expected = (sp.factor(u / x**2), 1 / x, 0)
    assert gamma_image == expected

# Irreducible extra factors are counted over the algebraic closure.  A square
# of a quadratic contributes two roots times one excess multiplicity.
irreducible_H = sp.Rational(1, 25) * w**2 * (1 - w) * (w**2 + 4) ** 2
irreducible_profile = boundary_saturation_profile(irreducible_H, w)
assert irreducible_profile.extra_degree == 4
assert irreducible_profile.extra_repetition_excess == 2
assert irreducible_profile.saturation_exponent == 4

print("PASS: e=m_0+sum(mu_i-1) is the exact C-saturation exponent")
print("PASS: saturated boundary traces retain every conic and B-power stratum")
print("PASS: repeated extra roots add escaping branches but no finite C=0 chart")
print("PASS: nonsplit repeated factors are counted over the algebraic closure")
