#!/usr/bin/env python3
"""Exact invariant and Jacobian checks for the affine A4 Keller frontier."""

import sympy as sp


# ---------------------------------------------------------------------------
# 1. Polynomial squaring on the oriented cubic quotient
# ---------------------------------------------------------------------------

x, y, d, T = sp.symbols("x y d T")

delta = x**2 * y**2 - 4 * x**3 - 4 * y**3 + 18 * x * y - 27
X = x**2 - 2 * y
Y = y**2 - 2 * x
D = d * (x * y - 1)
delta_target = delta.subs({x: X, y: Y}, simultaneous=True)

assert sp.factor(delta_target - delta * (x * y - 1) ** 2) == 0

plane_jacobian = sp.factor(
    sp.Matrix([X, Y]).jacobian((x, y)).det()
)
assert sp.expand(plane_jacobian - 4 * (x * y - 1)) == 0
assert sp.cancel(plane_jacobian * d / D) == 4

# Given (X,Y,D) on the oriented target, the source x-coordinate satisfies
# this quartic. Its discriminant is (64D)^2 on D^2=delta(X,Y).
inverse_quartic = T**4 - 2 * x * T**2 - 8 * T + x**2 - 4 * y
assert sp.factor(
    sp.discriminant(inverse_quartic, T) - 4096 * delta
) == 0

example_quartic = sp.Poly(
    inverse_quartic.subs({x: -6, y: 3}),
    T,
    domain=sp.QQ,
)
example_resolvent = sp.Poly(
    T**3 - 12 * T**2 - 96 * T + 1088,
    T,
    domain=sp.QQ,
)
assert example_quartic.as_expr() == T**4 + 12 * T**2 - 8 * T + 24
assert example_quartic.is_irreducible
assert example_resolvent.is_irreducible
assert sp.discriminant(example_quartic.as_expr(), T) == 1728**2
assert delta.subs({x: -6, y: 3}) == 27**2


# ---------------------------------------------------------------------------
# 2. Jensen--Ledet--Yui's affine two-parameter generic A4 polynomial
# ---------------------------------------------------------------------------

a, b, Z = sp.symbols("a b Z")

Acoef = a**3 - b**3 - 9 * b**2 - 27 * b - 54
Bcoef = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
Ccoef = a**3 - b**3 + 27

# This is B^4 F(a,b,Z/B), so all coefficients are polynomial in (a,b).
generic_a4 = (
    Z**4
    - 6 * Acoef * Bcoef * Z**2
    - 8 * Bcoef**3 * Z
    + Bcoef**2 * (9 * Acoef**2 - 12 * Ccoef * Bcoef)
)
square_factor = (
    2 * a**3 * b
    + 3 * a**3
    - 3 * a**2 * b**2
    - 9 * a**2 * b
    - 27 * a**2
    + b**4
    + 6 * b**3
    + 27 * b**2
    + 54 * b
    + 81
)
expected_discriminant_root = (
    1728 * (b**2 + 3 * b + 9) * Bcoef**4 * square_factor
)
assert sp.factor(
    sp.discriminant(generic_a4, Z) - expected_discriminant_root**2
) == 0


# ---------------------------------------------------------------------------
# 3. Rational C3-invariant source coordinates and polynomial homogenization
# ---------------------------------------------------------------------------

s, t, U, V, W = sp.symbols("s t U V W")


def cyclic_invariant_1(first, second):
    denominator = second * (first - 1) * (second - 1) * (first * second - 1)
    numerator = (
        first**3 * second**3
        - 3 * first * second**2
        + second**3
        + 1
    )
    return sp.cancel(numerator / denominator)


def cyclic_invariant_2(first, second):
    denominator = second * (first - 1) * (second - 1) * (first * second - 1)
    numerator = (
        first**3 * second**3
        - 3 * first**2 * second**3
        + 6 * first * second**2
        - 3 * first * second
        + second**3
        - 3 * second**2
        + 1
    )
    return sp.cancel(numerator / denominator)


source_u = cyclic_invariant_1(s, t)
source_v = cyclic_invariant_2(s, t)
target_a = cyclic_invariant_1(s**2, t**2)
target_b = cyclic_invariant_2(s**2, t**2)

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

source_substitution = {U: source_u, V: source_v}
assert sp.factor(
    target_a - (N1 / H).subs(source_substitution)
) == 0
assert sp.factor(
    target_b - (N2 / H).subs(source_substitution)
) == 0

rational_map = sp.Matrix([N1 / H, N2 / H])
rational_jacobian = sp.factor(
    rational_map.jacobian((U, V)).det()
)
assert sp.factor(
    rational_jacobian - 4 * K**3 * L / H**3
) == 0

# Homogenizing the common denominator gives a literal polynomial A3 -> A3
# with the same generic four-sheet A4 extension. Its determinant divisor is
# now completely explicit.
cone_map = sp.Matrix([W * N1, W * N2, W * H])
cone_jacobian = sp.factor(
    cone_map.jacobian((U, V, W)).det()
)
assert cone_jacobian == 4 * W**2 * K**3 * L


# ---------------------------------------------------------------------------
# 4. First ambient-lift no-go around the oriented surface core
# ---------------------------------------------------------------------------

defect = d**2 - delta
unknowns = sp.symbols("f0:4 g0:4 h0:4")
kappa = sp.symbols("kappa")
f = unknowns[0] + unknowns[1] * x + unknowns[2] * y + unknowns[3] * d
g = unknowns[4] + unknowns[5] * x + unknowns[6] * y + unknowns[7] * d
h = unknowns[8] + unknowns[9] * x + unknowns[10] * y + unknowns[11] * d

ambient_lift = sp.Matrix([
    X + defect * f,
    Y + defect * g,
    D + defect * h,
])
ambient_determinant = sp.Poly(
    sp.expand(
        ambient_lift.jacobian((x, y, d)).det() - kappa
    ),
    x,
    y,
    d,
)
coefficient_ideal = [
    coefficient for _, coefficient in ambient_determinant.terms()
]
linear_lift_groebner = sp.groebner(
    coefficient_ideal,
    *unknowns,
    kappa,
    order="grevlex",
)
assert len(linear_lift_groebner.polys) == 1
assert linear_lift_groebner.polys[0].as_expr() == 1

print("PASS: oriented cubic squaring is polynomial with residue Jacobian 4")
print("PASS: its inverse quartic has discriminant (64D)^2 and A4 monodromy")
print("PASS: the scaled Jensen--Ledet--Yui generic A4 polynomial is polynomial")
print("PASS: the rational A4 quotient map has Jacobian 4*K^3*L/H^3")
print("PASS: its polynomial cone has determinant 4*W^2*K^3*L")
print("PASS: every defect-preserving affine-linear ambient correction is empty")
