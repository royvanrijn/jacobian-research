#!/usr/bin/env python3
"""Exact F_2 audit of the Huq--Kuruvilla marked-root normalization.

The map, collision, inverse cubic, and generic rational reconstruction are
credited to Irit Huq-Kuruvilla, arXiv:2607.20968.  This script additionally
checks the discriminant, projective collision, normalization charts,
boundary factorization, reconstruction pole, and determinant ledger used in
verified/HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md.
"""

from __future__ import annotations

import sympy as sp


def mod2(expr: sp.Expr, *generators: sp.Symbol) -> sp.Expr:
    """Return the canonical polynomial representative over F_2."""

    return sp.Poly(sp.expand(expr), *generators, modulus=2).as_expr()


def assert_poly_zero(expr: sp.Expr, *generators: sp.Symbol) -> None:
    assert mod2(expr, *generators) == 0


def assert_rat_zero(expr: sp.Expr, *generators: sp.Symbol) -> None:
    numerator, _ = sp.cancel(expr).as_numer_denom()
    assert_poly_zero(numerator, *generators)


x, y, z = sp.symbols("x y z")
P = x + x**2 * y
Q = y + x * z + x**2 * y * z
R = z + x**2 * z**2

jacobian = sp.det(sp.Matrix((P, Q, R)).jacobian((x, y, z)))
assert_poly_zero(jacobian - 1, x, y, z)
print("PASS: Huq--Kuruvilla map has determinant one over F_2")

points = ((0, 1, 0), (1, 1, 0), (1, 1, 1))
images = {
    tuple(
        int(mod2(component.subs({x: px, y: py, z: pz}), x, y, z))
        for component in (P, Q, R)
    )
    for px, py, pz in points
}
assert images == {(0, 1, 0)}
print("PASS: the credited three-point collision is exact")

# Huq--Kuruvilla's triangular target coordinates and primitive root.
V_source = mod2(Q + P * R, x, y, z)
t_source = (1 + x**2 * z) / x
a_source = P * t_source**2 + V_source
assert_rat_zero(x * a_source - 1, x, y, z)

U, V, W, t = sp.symbols("U V W t")
D = U * t**3 + t**2 + V * t + W
discriminant = sp.discriminant(D, t)
assert_poly_zero(discriminant - (V + U * W) ** 2, U, V, W)
print("PASS: inverse cubic discriminant is (V+UW)^2=Q^2")

# At (U,V,W)=(0,1,0), D_h=T*S*(T+S), so the roots are 0,1,infinity.
T, S = sp.symbols("T S")
D_h = U * T**3 + T**2 * S + V * T * S**2 + W * S**3
collision_binary = mod2(D_h.subs({U: 0, V: 1, W: 0}), T, S)
assert_poly_zero(collision_binary - T * S * (T + S), T, S)
print("PASS: projective inverse fiber displays roots 0, 1, and infinity")

# The finite normalization charts.
W_t = U * t**3 + t**2 + V * t
a = U * t**2 + V
q_t = V + U * W_t
assert_poly_zero(q_t - a * (1 + U * t), U, V, t)

s = sp.symbols("s")
U_s = s + V * s**2 + W * s**3
delta = 1 + W * s**2
eta = V + W * s
q_s = V + U_s * W
assert_poly_zero(q_s - delta * eta, V, W, s)
assert_poly_zero(U_s.subs(V, W * s) - s, W, s)
print("PASS: Q pulls back as E+A in both normalization charts")

# Reconstruction on the finite-root chart.
q = q_t
x_rec = 1 / a
z_rec = a * (t + a)
y_rec = q + U * z_rec

P_rec = x_rec + x_rec**2 * y_rec
Q_rec = y_rec + x_rec * z_rec + x_rec**2 * y_rec * z_rec
R_rec = z_rec + x_rec**2 * z_rec**2
for recovered, target in ((P_rec, U), (Q_rec, q), (R_rec, W_t)):
    assert_rat_zero(recovered - target, U, V, t)
print("PASS: rational reconstruction inverts the map on a!=0")

# The infinity chart includes the simple root at infinity and identifies the
# same boundary: on t=1/s, a=delta/s.
x_s = s / delta
z_s = W * delta
y_s = q_s + U_s * z_s
P_s = x_s + x_s**2 * y_s
Q_s = y_s + x_s * z_s + x_s**2 * y_s * z_s
R_s = z_s + x_s**2 * z_s**2
for recovered, target in ((P_s, U_s), (Q_s, q_s), (R_s, W)):
    assert_rat_zero(recovered - target, V, W, s)
assert_rat_zero(
    (U_s / s**2 + V) - delta / s,
    V,
    W,
    s,
)
print("PASS: infinity-chart reconstruction has the same sole pole delta=0")

# Boundary parameterization: a=0 gives V=Ut^2, W=t^2 and Q=0.
boundary_substitution = {V: U * t**2}
assert_poly_zero(
    (W_t - t**2).subs(boundary_substitution),
    U,
    t,
)
assert_poly_zero(q_t.subs(boundary_substitution), U, t)
print("PASS: E=(a=0) maps by (U,t)->(U,W=t^2) onto Q=0")

# Source-to-core chart and zero-pole determinant cancellation.
alpha = (P, V_source, t_source)
jacobian_alpha = sp.det(sp.Matrix(alpha).jacobian((x, y, z)))
assert_rat_zero(jacobian_alpha - x, x, y, z)
assert_rat_zero(a_source * jacobian_alpha - 1, x, y, z)
print("PASS: marked-root determinant ledger is a * J_alpha = 1")

# The identical integer formulas are not a characteristic-zero Keller lift.
integer_jacobian = sp.factor(jacobian)
expected_integer_jacobian = (
    1
    + 2 * x * y
    + 2 * x**2 * z
    + 4 * x**3 * y * z
    + 2 * x**4 * z**2
    + 2 * x**5 * y * z**2
)
assert sp.expand(integer_jacobian - expected_integer_jacobian) == 0
print("PASS: the naive integral lift has nonconstant Jacobian")

