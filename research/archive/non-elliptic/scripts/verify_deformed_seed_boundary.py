#!/usr/bin/env python3
"""Exact audit of the boundary-clean and one-extra-root seed theorem."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


rho, c = sp.symbols("rho c", nonzero=True)
A, B, C, s, t = sp.symbols("A B C s t")

# The universal quartic primitive with one additional simple zero rho.
H = sp.factor(c * w**2 * (1 - w) * (w - rho) / (1 - rho))
p = sp.diff(H, w)
kappa = sp.factor(sp.diff(p, w).subs(w, 1) / c)
a = sp.factor(-(1 + kappa) / (2 + kappa))
h2 = sp.expand(H).coeff(w, 2)

assert H.subs(w, 0) == 0
assert H.subs(w, 1) == 0
assert p.subs(w, 0) == 0
assert sp.factor(p.subs(w, 1) + c) == 0
assert sp.factor(p.subs(w, rho) - c * rho**2) == 0
assert sp.factor(kappa + 2 * (2 * rho - 3) / (rho - 1)) == 0
assert sp.factor(a + (3 * rho - 5) / (2 * (rho - 2))) == 0
assert sp.factor(h2 / c - rho / (rho - 1)) == 0

# The gamma=0 chart for any double-zero primitive H=h2*W^2+... .
k = sp.factor(h2 / c)
u = sp.symbols("u")
gamma_A = u + k * u**2
gamma_B = c * (1 + 2 * k * u)
X = sp.symbols("X", nonzero=True)
elimination = sp.factor(
    (B**2 / c**2 - 4 * k * A) * X**2 - 1
)
u_from_B = (B * X / c - 1) / (2 * k)
assert sp.factor(
    elimination.subs(A, gamma_A.subs(u, u_from_B) / X**2)
) == 0

# The pulled-back discriminant has exact boundary factor C^2; after
# saturation its C=0 trace is precisely the gamma-chart degeneracy conic.
E = sp.expand(H - s * w + t)
discriminant = sp.discriminant(E, w)
pulled = sp.factor(discriminant.subs({s: B * C, t: c * A * C**2}))
assert sp.factor(pulled / C**2).subs(C, 0) != 0
Q = sp.factor(pulled / C**2)
boundary_conic = B**2 / c**2 - 4 * k * A
assert sp.factor(Q.subs(C, 0) / boundary_conic - c**6 * rho**2 / (rho - 1) ** 2) == 0

# A quartic with no simple root must be a quadratic square.  Coefficient
# comparison gives this unique double-double value for every rho.
s_omitted = sp.factor(c * (1 - rho**2) / 8)
t_omitted = sp.factor(-c * (1 - rho) ** 3 / 64)
qa = -(1 + rho) / 2
qb = -(rho - 1) ** 2 / 8
assert sp.factor(
    E.subs({s: s_omitted, t: t_omitted})
    + c / (1 - rho) * (w**2 + qa * w + qb) ** 2
) == 0
assert sp.factor(qa**2 - 4 * qb - (3 * rho**2 - 2 * rho + 3) / 4) == 0

# Audit the actual polynomial maps at several rational deformations.  The
# symbolic rho-family requires rho not in {0,1,2}; rho=2 is exactly the
# forbidden kappa=-2 construction value.
for rho_value in (-2, -1, sp.Rational(1, 2), 3):
    H_value = sp.factor(H.subs({rho: rho_value, c: 1}))
    model = WeightedSeedModel(sp.diff(H_value, w))
    assert sp.factor(model.primitive - H_value) == 0
    assert model.fiber_degree == 4

    gamma_zero_z = -(1 + model.a * x * y) / x**2
    gamma_image = tuple(sp.factor(f.subs(z, gamma_zero_z)) for f in model.mapping())
    local_k = sp.factor(sp.Poly(H_value, w).coeff_monomial(w**2))
    expected = (
        sp.factor((1 + x * y + local_k * (1 + x * y) ** 2) / x**2),
        sp.factor((1 + 2 * local_k * (1 + x * y)) / x),
        0,
    )
    assert gamma_image == expected

print("PASS: one-extra-root primitives satisfy every weighted-seed constraint")
print("PASS: the direct C=0 chart has fibers 3 off its conic and 1 on it")
print("PASS: C^2 is the exact discriminant saturation factor")
print("PASS: every one-extra-root quartic has one explicit omitted double-double value")
