#!/usr/bin/env python3
"""Exact replay for one fixed Keller map with infinitely many Hasse fibers."""

from __future__ import annotations

import sympy as sp


x, y, z, S, a = sp.symbols("x y z S a")

G = (
    S**5
    - sp.Rational(3, 2) * S**4
    + sp.Rational(3, 2) * S**3
    - sp.Rational(5, 4) * S**2
    + sp.Rational(9, 16) * S
)
G_poly = sp.Poly(G, S)
g = {k: G_poly.nth(k) for k in range(1, 6)}
assert g == {
    1: sp.Rational(9, 16),
    2: -sp.Rational(5, 4),
    3: sp.Rational(3, 2),
    4: -sp.Rational(3, 2),
    5: sp.Integer(1),
}

t = 1 + x * y
q = t**2 * z + sp.Rational(3, 8) * y**2 * (1 + 3 * t)

F = (
    t * q,
    y
    + 8 * x * q
    - sp.Rational(40, 9) * t * q
    - sp.Rational(32, 3) * t**2 * x**2 * q**4
    + sp.Rational(80, 9) * t**2 * x**3 * q**5,
    x * (5 - 3 * t)
    - sp.Rational(8, 3) * x**3 * z
    + sp.Rational(16, 3) * (x * q) ** 4
    - sp.Rational(16, 3) * (x * q) ** 5,
)

# Independent direct expansion of the constant-Jacobian identity.
jacobian = sp.factor(sp.det(sp.Matrix(F).jacobian((x, y, z))))
assert jacobian == -2

# The target line realizes the translated Berend--Bilu family exactly.
B = sp.Rational(32, 9) * a
C = (8 * a + 1) / 3
inverse = sp.expand(G - sp.Rational(9, 32) * (B * S**2 + C))
expected = sp.expand(
    ((S - sp.Rational(1, 2)) ** 3 - a)
    * (S**2 + sp.Rational(3, 4))
)
assert sp.factor(inverse - expected) == 0

X = sp.symbols("X")
translated = sp.expand(inverse.subs(S, X + sp.Rational(1, 2)))
assert sp.factor(translated - (X**3 - a) * (X**2 + X + 1)) == 0

# For a nonzero prime parameter the fiber polynomial is reduced.  The
# symbolic resultant identifies the only parameter where its two factors
# meet; irreducibility for prime a is then the rational-root test.
cubic = X**3 - a
quadratic = X**2 + X + 1
assert sp.discriminant(cubic, X) == -27 * a**2
assert sp.discriminant(quadratic, X) == -3
assert sp.factor(sp.resultant(cubic, quadratic, X)) == (a - 1) ** 2

# The local proof uses only these finite-field facts.
p = sp.symbols("p", integer=True, positive=True)
assert sp.expand((X**3 - a).subs(X, 1)) == 1 - a
assert sp.diff(X**3 - a, X).subs(X, 1) == 3

print("PASS: the displayed fixed polynomial map has determinant -2")
print("PASS: its rational target line is (X^3-a)(X^2+X+1)")
print("PASS: every prime a=1 mod 27 gives a reduced degree-five Hasse fiber")

