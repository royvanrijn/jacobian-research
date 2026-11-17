#!/usr/bin/env python3
"""Exact regression for the easy degree-five Hessian-root witness pair."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


H_red = sp.expand(sp.Rational(1, 2) * w**2 * (w - 1) * (w**2 + 1))
H_dbl = sp.expand(w**2 * (w - 1) * (2 * w**2 - 2 * w + 1))

models = []
for H, expected_a in (
    (H_red, -sp.Rational(5, 4)),
    (H_dbl, -sp.Rational(7, 6)),
):
    p = sp.diff(H, w)
    model = WeightedSeedModel(p, c=-1)
    assert sp.factor(model.primitive - H) == 0
    assert model.a == expected_a
    assert model.fiber_degree == 5

    F = model.mapping()
    assert all(sp.denom(component) == 1 for component in F)
    jacobian = sp.factor(sp.Matrix(F).jacobian((x, y, z)).det())
    assert jacobian == -1
    models.append((model, F))

red_model, F_red = models[0]
dbl_model, F_dbl = models[1]

# The double-Hessian witness has an exact rational collision over (-20,20,1).
target = (-20, 20, 1)
P1 = (-sp.Rational(1, 19), 20, -sp.Rational(22990, 3))
P2 = (sp.Rational(1, 44), -42, 81092)

for point in (P1, P2):
    image = tuple(
        sp.factor(component.subs(dict(zip((x, y, z), point))))
        for component in F_dbl
    )
    assert image == target

E = sp.expand(H_dbl - 20 * w + 20)
assert E.subs(w, 1) == 0
assert E.subs(w, 2) == 0
assert sp.diff(E, w).subs(w, 1) == -19
assert sp.diff(E, w).subs(w, 2) == 44

# Both split primitives have two distinct, boundary-clean additional roots.
for H, residual, expected_resultant in (
    (H_red, w**2 + 1, 1),
    (H_dbl, 2 * w**2 - 2 * w + 1, 16),
):
    assert sp.discriminant(residual, w) == -4
    resultant = sp.factor(sp.resultant(residual, -w + sp.diff(H, w), w))
    assert resultant == expected_resultant

# Stable obstruction: three Hessian-root components versus two.
K_red = sp.factor(sp.diff(H_red, w, 2))
K_dbl = sp.factor(sp.diff(H_dbl, w, 2))

assert K_red == 10 * w**3 - 6 * w**2 + 3 * w - 1
assert K_dbl == 2 * (2 * w - 1) ** 2 * (5 * w - 1)
assert sp.discriminant(K_red, w) == -1080
assert sp.discriminant(K_dbl, w) == 0
assert sp.gcd(K_red, sp.diff(K_red, w)) == 1
assert sp.factor(sp.gcd(K_dbl, sp.diff(K_dbl, w))) == 2 * w - 1
assert sp.degree(sp.sqf_part(K_red), w) == 3
assert sp.degree(sp.sqf_part(K_dbl), w) == 2

print("PASS: both rational weighted maps are polynomial Keller maps of degree five")
print("PASS: the double-Hessian map has the stated rational collision")
print("PASS: both additional-root pairs are boundary-clean")
print("PASS: three versus two Hessian-root components separates the stable classes")
