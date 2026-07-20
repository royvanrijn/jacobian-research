#!/usr/bin/env python3
"""Exact SymPy certificate for the marked-projective-root model."""

import sympy as sp


x, y, z = sp.symbols("x y z")
a, b, c = sp.symbols("a b c")
U, V, T = sp.symbols("U V T")
u, v, s = sp.symbols("u v s")

w = 1 + x * y
F = (
    w**3 * z + y**2 * w * (4 + 3 * x * y),
    y + 3 * x * w**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)

Q = c * U**3 - 2 * U**2 * V + b * U * V**2 - 2 * a * V**3
marked_identity = sp.expand(
    Q.subs({a: F[0], b: F[1], c: F[2], U: w, V: x})
)
assert marked_identity == 0
source_s = x / w
source_d = sp.factor(
    (1 - F[1] * source_s + 3 * F[0] * source_s**2) - 1 / w
)
assert source_d == 0
print("PASS: source marking is a projective root")

# Affine root chart: substitute x=1/v, y=u-v and solve c=F_3 for z.
x_aff = 1 / v
y_aff = u - v
z_aff = 5 * v**2 - 3 * u * v - c * v**3
F_aff = tuple(sp.factor(expr.subs({x: x_aff, y: y_aff, z: z_aff})) for expr in F)
a_aff = u**2 + u * v - c * u**3
b_aff = 4 * u + 2 * v - 3 * c * u**2
assert sp.factor(F_aff[0] - a_aff) == 0
assert sp.factor(F_aff[1] - b_aff) == 0
assert sp.factor(F_aff[2] - c) == 0
P = c * T**3 - 2 * T**2 + b * T - 2 * a
assert sp.expand(Q.subs({U: T, V: 1}) - P) == 0
print("PASS: affine chart formulas for a, b, and P")

P_aff = P.subs({a: a_aff, b: b_aff})
assert sp.factor(P_aff.subs(T, u)) == 0
print("PASS: P(u) = 0")

Pprime_u = sp.factor(sp.diff(P, T).subs({a: a_aff, b: b_aff, T: u}))
assert Pprime_u == 2 * v
print("PASS: P'(u) = 2v")

v_rec = sp.diff(P, T).subs(T, u) / 2
reconstruction = (
    1 / v_rec,
    u - v_rec,
    5 * v_rec**2 - 3 * u * v_rec - c * v_rec**3,
)
reconstruction_aff = tuple(
    sp.factor(expr.subs({a: a_aff, b: b_aff}) - expected)
    for expr, expected in zip(reconstruction, (x_aff, y_aff, z_aff))
)
assert reconstruction_aff == (0, 0, 0)
print("PASS: affine simple-root reconstruction recovers x, y, z")

J_target = sp.Matrix([a_aff, b_aff, c]).jacobian((u, v, c)).det()
c_source = F[2]
J_chart = sp.Matrix([y + 1 / x, 1 / x, c_source]).jacobian((x, y, z)).det()
assert sp.factor(J_target - 2 * v) == 0
assert sp.factor(J_chart.subs(x, 1 / v) + 1 / v) == 0
assert sp.factor(J_target * J_chart.subs(x, 1 / v) + 2) == 0
print("PASS: determinant factors as (2v)(-1/v) = -2")

# U != 0 chart.  The root equation is used in solved form for c.
root_s = sp.expand(Q.subs({U: 1, V: s}))
assert root_s == c - 2 * s + b * s**2 - 2 * a * s**3
c_s = 2 * s - b * s**2 + 2 * a * s**3
d = 1 - b * s + 3 * a * s**2
assert sp.factor(sp.diff(root_s, s).subs(c, c_s) + 2 * d) == 0
print("PASS: projective s-chart root and simple-root equations")

x_s = s / d
y_s = b - 3 * a * s
z_factor = (
    a
    - 4 * b**2
    + (b**3 + 22 * a * b) * s
    - (30 * a**2 + 8 * a * b**2) * s**2
    + 21 * a**2 * b * s**3
    - 18 * a**3 * s**4
)
z_s = sp.expand(d * z_factor)
raw_z = sp.cancel((2 * x_s - 3 * x_s**2 * y_s - c_s) / x_s**3)
assert sp.factor(raw_z - z_s) == 0
for got, want in zip(F, (a, b, c_s)):
    assert sp.factor(got.subs({x: x_s, y: y_s, z: z_s}) - want) == 0
print("PASS: s-chart reconstruction formulas recover the target")

assert sp.denom(z_s) == 1
assert (d.subs(s, 0), x_s.subs(s, 0), y_s.subs(s, 0), z_s.subs(s, 0)) == (
    1,
    0,
    b,
    a - 4 * b**2,
)
print("PASS: s-chart reconstruction specializes regularly at s = 0")

infinity_source = (0, b, a - 4 * b**2)
assert tuple(sp.factor(expr.subs({x: infinity_source[0], y: infinity_source[1], z: infinity_source[2]})) for expr in F) == (
    a,
    b,
    0,
)
assert root_s.subs({s: 0, c: 0}) == 0
print("PASS: [1:0] reconstructs (0, b, a - 4b^2)")

Qdisc = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
raw_discriminant = sp.factor(sp.discriminant(P, T))
assert sp.factor(raw_discriminant / 4 + Qdisc) == 0
print("PASS: normalized cubic discriminant is -Qdisc (raw Disc(P) = -4Qdisc)")

triple_a = sp.Rational(4, 27) / c**2
triple_b = sp.Rational(4, 3) / c
triple_root = sp.Rational(2, 3) / c
assert sp.factor(3 * triple_b * c - 4) == 0
assert sp.factor(27 * triple_a * c**2 - 4) == 0
assert sp.solve((3 * b * c - 4, 27 * a * c**2 - 4), (a, b), dict=True) == [
    {a: triple_a, b: triple_b}
]
assert sp.expand(P.subs({a: triple_a, b: triple_b}) - c * (T - triple_root) ** 3) == 0
print("PASS: triple roots are exactly 3bc = 4 and 27ac^2 = 4")

collision_Q = sp.factor(Q.subs({a: -sp.Rational(1, 4), b: 0, c: 0}))
assert collision_Q == -V * (2 * U - V) * (2 * U + V) / 2
collision_roots = ((1, 0), (sp.Rational(1, 2), 1), (-sp.Rational(1, 2), 1))
assert all(collision_Q.subs({U: root_u, V: root_v}) == 0 for root_u, root_v in collision_roots)
collision_sources = (
    (0, 0, -sp.Rational(1, 4)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
    (1, -sp.Rational(3, 2), sp.Rational(13, 2)),
)
source_markings = tuple((1 + px * py, px) for px, py, _ in collision_sources)
assert all(
    marked_u * root_v == marked_v * root_u
    for (marked_u, marked_v), (root_u, root_v) in zip(source_markings, collision_roots)
)
print("PASS: collision target has the three stated projective roots")
