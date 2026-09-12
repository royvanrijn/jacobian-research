#!/usr/bin/env python3
"""Exact checks for the A4 determinant-ledger reduction and rigidity."""

import sympy as sp


U, V, W = sp.symbols("U V W")

H = (
    8 * U**3
    - 6 * U * V**2
    - 18 * U * V
    - 54 * U
    - 2 * V**3
    - 9 * V**2
    - 27 * V
    - 27
)
K = 4 * U**2 + 4 * U * V + 6 * U + V**2 + 3 * V + 9
M = U**2 + 2 * V**2 + 6 * V + 18
L = (
    U**3
    - 3 * U * V**2
    - 9 * U * V
    - 27 * U
    + 2 * V**3
    + 9 * V**2
    + 27 * V
    + 27
)
N1 = sp.expand(M * K)
N2 = (
    8 * U**3 * V
    + 12 * U**2 * V**2
    + 36 * U**2 * V
    + 108 * U**2
    + 6 * U * V**3
    + 36 * U * V**2
    + 108 * U * V
    + 162 * U
    + V**4
    + 9 * V**3
    + 27 * V**2
    + 54 * V
)


# ---------------------------------------------------------------------------
# 1. The JLY B-polynomial absorbs K^3 and one copy of L
# ---------------------------------------------------------------------------

a = N1 / H
b = N2 / H
B_rational = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
assert sp.factor(B_rational - K**3 * L**2 / H**3) == 0

P, Q, R = sp.symbols("P Q R")
B_homogeneous = (
    P**3
    - 3 * P * Q**2
    + 2 * Q**3
    - 9 * P * Q * R
    + 9 * Q**2 * R
    - 27 * P * R**2
    + 27 * Q * R**2
    + 27 * R**3
)
cone = {P: W * N1, Q: W * N2, R: W * H}
B_cone_pullback = sp.factor(B_homogeneous.subs(cone))
assert B_cone_pullback == W**3 * K**3 * L**2

cone_map = sp.Matrix([W * N1, W * N2, W * H])
cone_jacobian = sp.factor(cone_map.jacobian((U, V, W)).det())
assert cone_jacobian == 4 * W**2 * K**3 * L
assert sp.factor(
    sp.cancel(B_cone_pullback / cone_jacobian) - W * L / 4
) == 0


# ---------------------------------------------------------------------------
# 2. Elementary geometry of K and L
# ---------------------------------------------------------------------------

assert sp.expand(K - ((2 * U + V) ** 2 + 3 * (2 * U + V) + 9)) == 0

t = sp.symbols("t")
V_param = (-t**3 + 27 * t - 27) / (3 * t * (t - 3))
U_param = (2 * t - 3) * (t**2 - 3 * t + 9) / (3 * t * (t - 3))
assert sp.factor(L.subs({U: U_param, V: V_param})) == 0
assert sp.factor(U_param - V_param - t) == 0
assert sp.factor(L.subs(U, V)) == 27
assert sp.factor(L.subs(U, V + 3)) == -27


# ---------------------------------------------------------------------------
# 3. The pole boundary H has a nonconstant projective image
# ---------------------------------------------------------------------------

assert sp.Poly(H, U, V, domain=sp.QQ).is_irreducible
assert sp.rem(N1, H, U) != 0
assert sp.rem(N2, H, U) != 0

# Tangent derivation to H=0: D = H_V d/dU - H_U d/dV.
tangent_N1 = sp.diff(H, V) * sp.diff(N1, U) - sp.diff(H, U) * sp.diff(N1, V)
tangent_N2 = sp.diff(H, V) * sp.diff(N2, U) - sp.diff(H, U) * sp.diff(N2, V)
ratio_derivative_numerator = sp.expand(tangent_N1 * N2 - N1 * tangent_N2)
assert sp.rem(ratio_derivative_numerator, H, U) != 0


# ---------------------------------------------------------------------------
# 4. Every defect-multiple ambient correction has a singular derivative
# ---------------------------------------------------------------------------

x, y, d = sp.symbols("x y d")
delta = x**2 * y**2 - 4 * x**3 - 4 * y**3 + 18 * x * y - 27
defect = d**2 - delta
oriented_map = sp.Matrix([
    x**2 - 2 * y,
    y**2 - 2 * x,
    d * (x * y - 1),
])
point = {x: -1, y: -1, d: 0}

base_derivative = oriented_map.jacobian((x, y, d)).subs(point)
defect_gradient = sp.Matrix([[
    sp.diff(defect, variable).subs(point)
    for variable in (x, y, d)
]])
assert defect.subs(point) == 0
assert base_derivative == sp.Matrix([
    [-2, -2, 0],
    [-2, -2, 0],
    [0, 0, 0],
])
assert defect_gradient == sp.Matrix([[32, 32, 0]])

f_value, g_value, h_value = sp.symbols("f_value g_value h_value")
arbitrary_rank_one_correction = (
    sp.Matrix([f_value, g_value, h_value]) * defect_gradient
)
assert sp.factor(
    (base_derivative + arbitrary_rank_one_correction).det()
) == 0


print("PASS: the target B-divisor pulls back as W^3*K^3*L^2")
print("PASS: the residual determinant ledger is W*L/4")
print("PASS: K splits geometrically into two parallel lines")
print("PASS: the normalization of L is P1 minus {0,3,infinity}")
print("PASS: the H-boundary ratio [N1:N2] is nonconstant")
print("PASS: every defect-multiple correction has zero Jacobian at (-1,-1,0)")
