#!/usr/bin/env python3
"""Exact F_2 replay of the Huq--Kuruvilla map and Mondello's plane theorem.

The map, collision, inverse cubic, and generic rational reconstruction are
credited to Irit Huq-Kuruvilla, arXiv:2607.20968.  This script additionally
checks the discriminant, projective collision, normalization charts,
boundary factorization, reconstruction pole, and determinant ledger used in
verified/HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md.

The plane theorem and proof architecture are due to Romy Mondello,
arXiv:2608.02634v1.  The final block internally replays the preserved-coordinate
reduction, plane Jacobian and collision, hidden cubic, recovery identities,
and separability witness.  Irreducibility of the hidden cubic is the written
degree-one-in-the-actual-target-parameter argument in the canonical note; the
checker verifies its polynomial coprimality certificate.  Because these are
polynomial identities over F_2, the replay also supports the repository's
separately proved arbitrary-characteristic-two base-field corollary.
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

# Mondello's coordinate-permuted skew-product form.  Use fresh symbols to
# keep this calculation independent of the normalization-chart variables.
xp, yp, zp = sp.symbols("xp yp zp")
aa, bb, cc = sp.symbols("aa bb cc")

phi_x = xp + xp**2 * zp
phi_y = yp + xp**2 * yp**2
phi_z = zp + xp * yp + xp**2 * yp * zp

aa_source = xp
bb_source = zp + yp**2
cc_source = xp + yp + xp**2 * bb_source

xp_inverse = aa
yp_inverse = cc + aa + aa**2 * bb
zp_inverse = bb + yp_inverse**2
inverse_substitution = {xp: xp_inverse, yp: yp_inverse, zp: zp_inverse}
source_substitution = {aa: aa_source, bb: bb_source, cc: cc_source}

for recovered, target in (
    (aa_source.subs(inverse_substitution), aa),
    (bb_source.subs(inverse_substitution), bb),
    (cc_source.subs(inverse_substitution), cc),
    (xp_inverse.subs(source_substitution), xp),
    (yp_inverse.subs(source_substitution), yp),
    (zp_inverse.subs(source_substitution), zp),
):
    assert_poly_zero(recovered - target, xp, yp, zp, aa, bb, cc)

A_skew = mod2(phi_x.subs(inverse_substitution), aa, bb, cc)
B_skew = mod2(phi_z.subs(inverse_substitution), aa, bb, cc)
C_skew = mod2((phi_x + phi_y).subs(inverse_substitution), aa, bb, cc)
assert_poly_zero(C_skew - cc, aa, bb, cc)

expected_plane_a = aa + aa**2 * bb + aa**4 + aa**6 * bb**2
expected_plane_b = bb + aa**5 + aa**6 * bb + aa**7 * bb**2 + aa**8 * bb**3
assert_poly_zero(A_skew.subs(cc, 0) - expected_plane_a, aa, bb)
assert_poly_zero(B_skew.subs(cc, 0) - expected_plane_b, aa, bb)
print("PASS: the coordinate-permuted threefold map is a skew product with the stated plane fiber")

# Direct plane Keller and collision certificates.
plane_x, plane_y = sp.symbols("plane_x plane_y")
P2 = plane_x + plane_x**2 * plane_y + plane_x**4 + plane_x**6 * plane_y**2
Q2 = (
    plane_y
    + plane_x**5
    + plane_x**6 * plane_y
    + plane_x**7 * plane_y**2
    + plane_x**8 * plane_y**3
)
plane_jacobian = sp.det(sp.Matrix((P2, Q2)).jacobian((plane_x, plane_y)))
assert_poly_zero(plane_jacobian - 1, plane_x, plane_y)

plane_points = ((0, 1), (1, 0), (1, 1))
plane_images = {
    tuple(
        int(mod2(component.subs({plane_x: px, plane_y: py}), plane_x, plane_y))
        for component in (P2, Q2)
    )
    for px, py in plane_points
}
assert plane_images == {(0, 1)}
print("PASS: the plane map has determinant one and the three-point collision")

# Hidden cubic, source-field recovery, and separability witness.
plane_r = 1 + plane_x * plane_y
plane_u = 1 + plane_x**3 * plane_r
plane_w = plane_r * plane_u**2
assert_poly_zero(P2 - plane_x * plane_r * plane_u, plane_x, plane_y)
assert_poly_zero(Q2 - (plane_y + plane_x**5 * plane_r**3), plane_x, plane_y)
assert_poly_zero(plane_x * Q2 - (1 + plane_w), plane_x, plane_y)
assert_poly_zero(P2**2 - plane_x**2 * plane_r * plane_w, plane_x, plane_y)
assert_poly_zero(P2 * Q2 + P2**3 - (plane_r * plane_u + plane_w**2), plane_x, plane_y)

hidden_cubic_at_w = (
    plane_w**3
    + plane_w**2
    + (P2 * Q2 + P2**3) * plane_w
    + P2**3
)
assert_poly_zero(hidden_cubic_at_w, plane_x, plane_y)
hidden_derivative_at_w = plane_w**2 + P2 * Q2 + P2**3
assert_poly_zero(hidden_derivative_at_w - plane_r * plane_u, plane_x, plane_y)

# In K=k(P_target), a Q_target-independent factor would divide these two
# coefficient polynomials.  Their gcd is one already over F_2[P_target,T].
P_target, T_plane = sp.symbols("P_target T_plane")
q_coefficient = P_target * T_plane
q_constant = T_plane**3 + T_plane**2 + P_target**3 * T_plane + P_target**3
coprimality = sp.gcd(
    sp.Poly(q_coefficient, P_target, T_plane, modulus=2),
    sp.Poly(q_constant, P_target, T_plane, modulus=2),
)
assert coprimality.as_expr() == 1

elimination = Q2 * (
    Q2**2 * plane_x**3
    + (P2**3 + P2 * Q2 + 1) * plane_x
    + P2
)
assert_poly_zero(elimination, plane_x, plane_y)
print("PASS: the plane hidden cubic, recovery identities, irreducibility certificate, and separability witness hold")

# The same plane formulas do not give a characteristic-zero Keller map.
expected_plane_integer_jacobian = (
    1
    + 2 * plane_x * plane_y
    + 4 * plane_x**3
    - 4 * plane_x**6
    - 2 * plane_x**7 * plane_y
    + 4 * plane_x**9
    - 2 * plane_x**9 * plane_y**3
    - 2 * plane_x**10 * plane_y
    + 6 * plane_x**5 * plane_y**2
    + 6 * plane_x**11 * plane_y**2
    - 2 * plane_x**12 * plane_y**3
    + 2 * plane_x**13 * plane_y**4
)
assert sp.expand(plane_jacobian - expected_plane_integer_jacobian) == 0
print("PASS: the naive integral plane lift has nonconstant Jacobian")
