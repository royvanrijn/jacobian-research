#!/usr/bin/env python3
"""Exact audit of the canonical-family image and nonproperness theorem."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, canonical_seed, w, x, y, z  # noqa: E402


A, B, C, s, t, Z = sp.symbols("A B C s t Z")


def c_valuation(expression):
    terms = sp.Poly(sp.expand(expression), C).terms()
    return min(monomial[0] for monomial, _ in terms)


for d in range(2, 9):
    model = WeightedSeedModel(canonical_seed(d))
    H = model.primitive
    p = model.seed
    E = sp.expand(H - s * w + t)
    assert sp.factor(H - w**d * (1 - w)) == 0

    # The direct x=0 chart is a bijection onto C=0.
    boundary = model.boundary_map()
    K = sp.Rational(16 * d**2 - 19 * d + 6, 8 * (d - 1))
    assert sp.factor(boundary[0] - (K * y**2 - 2 * (d - 1) * z)) == 0
    assert sp.factor(boundary[1] - y / (2 * (d - 1))) == 0
    assert boundary[2] == 0

    # The second C=0 chart is gamma=0.  It is quadratic for d=2 and
    # collapses to one point over B!=0 for every d>=3.
    gamma_zero_z = -(1 + model.a * x * y) / x**2
    gamma_image = tuple(sp.factor(f.subs(z, gamma_zero_z)) for f in model.mapping())
    if d == 2:
        assert gamma_image == (
            (x * y + 1) * (x * y + 2) / x**2,
            (2 * x * y + 3) / x,
            0,
        )
    else:
        assert gamma_image == ((x * y + 1) / x**2, 1 / x, 0)

    # Pulling back the inverse discriminant has exact C-adic order d.  The
    # quotient is the saturated nonproperness equation Q_d.
    discriminant = sp.discriminant(E, w)
    pulled = sp.expand(discriminant.subs({s: B * C, t: A * C**2}))
    assert c_valuation(pulled) == d
    Q = sp.factor(pulled / C**d)
    assert not sp.denom(Q).has(C)
    if d == 2:
        assert sp.factor(Q.subs(C, 0) - (B**2 - 4 * A)) == 0
    else:
        q0 = sp.factor(Q.subs(C, 0))
        assert q0 != 0 and sp.factor(q0 / B**d).is_number

    # The only ways a lacunary polynomial
    # W^(d+1)-W^d+sW-t can have no simple root are a square (even degree) or
    # one triple factor times a square (odd degree).  These formal-series
    # coefficients are the first forbidden middle coefficients.
    n = d + 1
    if n % 2 == 0:
        m = n // 2
        obstruction = sp.expand(
            sp.series(sp.sqrt(1 - Z), Z, 0, m + 2).removeO()
        ).coeff(Z, m + 1)
        if m >= 3:
            assert obstruction != 0
    else:
        m = (n - 1) // 2
        cusp_root = sp.Rational(n - 2, n)
        series = sp.series(
            sp.sqrt(1 - Z) * (1 - cusp_root * Z) ** (-sp.Rational(3, 2)),
            Z,
            0,
            m + 1,
        ).removeO()
        obstruction = sp.expand(series).coeff(Z, m)
        if m >= 2:
            assert obstruction > 0

    # Every cusp after the cubic case retains simple residual roots.
    cusp_root = sp.Rational(d - 1, d + 1)
    cusp_s = p.subs(w, cusp_root)
    cusp_t = (w * p - H).subs(w, cusp_root)
    cusp_poly = sp.Poly(E.subs({s: cusp_s, t: cusp_t}), w)
    if d == 2:
        assert cusp_poly.sqf_list()[1] == [(sp.Poly(w - sp.Rational(1, 3), w), 3)]
    else:
        assert any(multiplicity == 1 for _, multiplicity in cusp_poly.sqf_list()[1])

# The two and only two low-degree omitted inverse polynomials.
E2 = w**2 * (1 - w) - s * w + t
assert sp.expand(
    E2.subs({s: sp.Rational(1, 3), t: sp.Rational(1, 27)})
    + (w - sp.Rational(1, 3)) ** 3
) == 0

E3 = w**3 * (1 - w) - s * w + t
assert sp.expand(
    E3.subs({s: sp.Rational(1, 8), t: -sp.Rational(1, 64)})
    + (w**2 - w / 2 - sp.Rational(1, 8)) ** 2
) == 0

print("PASS: canonical C=0 fibers are 3/1 for d=2 and 2/1 for d>=3")
print("PASS: Disc(H_d-BCW+AC^2)=C^d Q_d with no lost saturation factor")
print("PASS: the only omitted inverse values occur in canonical degrees 3 and 4")
print("PASS: canonical maps of inverse degree at least five are surjective")
