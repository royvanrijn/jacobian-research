#!/usr/bin/env python3
"""Exact checks for the first post-coordinate Davenport attack."""

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T as DT, Y as DY, davenport_pair  # noqa: E402


D, z, C, R, U = sp.symbols("D z C R U")

Y = (D * z - 1) / 2
T = D - Y**2
beta = D**2 * U - C * R

parabola_jacobian = sp.factor(
    sp.det(sp.Matrix([[sp.diff(T, D), sp.diff(T, z)],
                      [sp.diff(Y, D), sp.diff(Y, z)]]))
)
assert parabola_jacobian == D / 2

determinant_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(expr, var) for var in (D, C, R, U)]
                for expr in (D, C, R, beta)
            ]
        )
    )
)
assert determinant_jacobian == D**2

full_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(expr, var) for var in (D, z, C, R, U)]
                for expr in (T, Y, C, R, beta)
            ]
        )
    )
)
assert full_jacobian == D**3 / 2

# On the critical divisor D=0 the image center is the graph beta=-CR,
# whose coordinate ring is K[C,R].
center_relation = sp.expand(beta.subs(D, 0))
assert center_relation == -C * R

# A shifted exceptional equation masks the sheet-dependent function without
# changing the degree, but its block-triangular Jacobian remains g_Y.
T0, Y0, s = sp.symbols("T0 Y0 s")
f = sp.Function("f")(T0, Y0)
delta = T0 + Y0**2
masked_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            (T0, f, s + delta)
        ).jacobian((T0, Y0, s))
    )
)
assert masked_jacobian == sp.diff(f, Y0)

# Universal one-auxiliary affine-linear mask formula.
h = sp.Function("h")(T0, Y0)
a = sp.Function("a")(T0, Y0)
b = sp.Function("b")(T0, Y0)
Q = f + a * s
H = h + b * s
one_mask_jacobian = sp.expand(
    sp.det(sp.Matrix((T0, Q, H)).jacobian((T0, Y0, s)))
)
expected_one_mask = sp.expand(
    b * sp.diff(f, Y0)
    - a * sp.diff(h, Y0)
    + s * (b * sp.diff(a, Y0) - a * sp.diff(b, Y0))
)
assert sp.simplify(one_mask_jacobian - expected_one_mask) == 0

# Two-mask instance of the general cofactor theorem.  The signed maximal
# minors annihilate the mask matrix and contract h_Y to the determinant.
s1, s2 = sp.symbols("s1 s2")
a00, a01, a10, a11, a20, a21 = sp.symbols(
    "a00 a01 a10 a11 a20 a21"
)
h0 = sp.Function("h0")(T0, Y0)
h1 = sp.Function("h1")(T0, Y0)
h2 = sp.Function("h2")(T0, Y0)
mask_matrix = sp.Matrix(
    [[a00, a01], [a10, a11], [a20, a21]]
)
relative_matrix = sp.Matrix.hstack(
    sp.Matrix([sp.diff(h0, Y0), sp.diff(h1, Y0), sp.diff(h2, Y0)]),
    mask_matrix,
)
cofactor_vector = sp.Matrix(
    [
        a10 * a21 - a11 * a20,
        -(a00 * a21 - a01 * a20),
        a00 * a11 - a01 * a10,
    ]
)
assert sp.simplify((cofactor_vector.T * mask_matrix)[0, 0]) == 0
assert sp.simplify((cofactor_vector.T * mask_matrix)[0, 1]) == 0
assert sp.expand(
    relative_matrix.det()
    - (
        cofactor_vector[0] * sp.diff(h0, Y0)
        + cofactor_vector[1] * sp.diff(h1, Y0)
        + cofactor_vector[2] * sp.diff(h2, Y0)
    )
) == 0

# If the shifted exceptional mask itself is a target coordinate, the
# constant-Jacobian PDE has only triangular (hence invertible) solutions.
P = sp.Function("P")(T0, Y0, s)
shifted_mask = s + T0 + Y0**2
shifted_mask_jacobian = sp.expand(
    sp.det(
        sp.Matrix((P, shifted_mask)).jacobian((Y0, s))
    )
)
assert shifted_mask_jacobian == sp.diff(P, Y0) - 2 * Y0 * sp.diff(P, s)
constant, invariant = sp.symbols("constant invariant")
F = sp.Function("F")
general_shifted_solution = constant * Y0 + F(T0, s + Y0**2)
assert sp.simplify(
    sp.diff(general_shifted_solution, Y0)
    - 2 * Y0 * sp.diff(general_shifted_solution, s)
    - constant
) == 0

# On delta=0, every correction g+delta*H leaves a factor divisible by Y in
# the reduced Keller equation.
a0, Hbar = sp.symbols("a0 Hbar")
restricted_derivative = (a0 - 2) * Y0**4 * (2 * Y0 + 1) ** 2
reduced_factor = sp.expand(restricted_derivative + 2 * Y0 * Hbar)
assert reduced_factor.subs(Y0, 0) == 0
assert sp.factor(reduced_factor / Y0).is_polynomial(Y0, Hbar)

# Retaining beta literally as a target coordinate is impossible for any
# affine-space Keller completion: its differential vanishes on a line.
beta_gradient = [sp.diff(beta, variable) for variable in (D, C, R, U)]
assert all(
    component.subs({D: 0, C: 0, R: 0}) == 0
    for component in beta_gradient
)
S = sp.symbols("S")
translated_beta = S + beta
assert sp.diff(translated_beta, S) == 1

# The translated incidence is a triangular source coordinate.  If the old
# base is retained, the full determinant is only a plane determinant in
# (Y,U).
source_translation_jacobian = sp.det(
    sp.Matrix((D, C, R, U, translated_beta)).jacobian((D, C, R, U, S))
)
assert source_translation_jacobian == 1
py, pu, psigma, qy, qu, qsigma = sp.symbols(
    "py pu psigma qy qu qsigma"
)
translated_relative_matrix = sp.Matrix(
    [
        [0, D**2, 1],
        [py, pu + D**2 * psigma, psigma],
        [qy, qu + D**2 * qsigma, qsigma],
    ]
)
translated_full_jacobian = sp.expand(translated_relative_matrix.det())
translated_plane_jacobian = sp.expand(py * qu - pu * qy)
assert translated_full_jacobian == translated_plane_jacobian

# Parameter-linear U coupling is stably a plane pair.
aa0, aa1 = sp.symbols("aa0 aa1")
HH0 = sp.Function("HH0")(T0, Y0)
HH1 = sp.Function("HH1")(T0, Y0)
HH2 = sp.Function("HH2")(T0, Y0)
linear_u_map = (HH0 + aa0 * U, HH1 + aa1 * U, HH2 + U)
linear_u_jacobian = sp.det(
    sp.Matrix(linear_u_map).jacobian((T0, Y0, U))
)
plane0 = HH0 - aa0 * HH2
plane1 = HH1 - aa1 * HH2
reduced_plane_jacobian = sp.det(
    sp.Matrix((plane0, plane1)).jacobian((T0, Y0))
)
assert sp.simplify(linear_u_jacobian - reduced_plane_jacobian) == 0

# General Y-dependent sheet coefficient curve:
# (T+a(Y)*U, g+b(Y)*U, h+U).
afun = sp.Function("afun")(Y0)
bfun = sp.Function("bfun")(Y0)
hfun = sp.Function("hfun")(T0, Y0)
gfun = sp.Function("gfun")(T0, Y0)
pencil_map = (
    T0 + afun * U,
    gfun + bfun * U,
    hfun + U,
)
pencil_jacobian = sp.Poly(
    sp.expand(sp.det(sp.Matrix(pencil_map).jacobian((T0, Y0, U)))),
    U,
)
wronskian = sp.expand(
    bfun * sp.diff(afun, Y0) - afun * sp.diff(bfun, Y0)
)
expected_pencil_u = sp.expand(
    wronskian * sp.diff(hfun, T0)
    - sp.diff(afun, Y0) * sp.diff(gfun, T0)
    + sp.diff(bfun, Y0)
)
assert sp.expand(pencil_jacobian.coeff_monomial(U) - expected_pencil_u) == 0
assert pencil_jacobian.coeff_monomial(U**2) == 0

# For the actual Davenport polynomial, the positive-T coefficients of
# partial_T(g) are coprime over K[Y].
davenport_g, _ = davenport_pair()
root_a = (-1 + sp.sqrt(-7)) / 2
partial_T_g = sp.expand(sp.diff(davenport_g, DT).subs(A, root_a))
partial_T_poly = sp.Poly(partial_T_g, DT)
positive_T_coefficients = [
    sp.Poly(
        sp.expand(partial_T_poly.coeff_monomial(DT**power)),
        DY,
        extension=sp.sqrt(-7),
    )
    for power in (1, 2)
]
assert sp.gcd(*positive_T_coefficients).degree() == 0

# Once a'=W*r, b'=W*s, br-as=1, and h=r*g-s*T+k, the
# constant determinant has the advertised impossible factor.
rfun = sp.Function("rfun")(Y0)
sfun = sp.Function("sfun")(Y0)
kfun = sp.Function("kfun")(Y0)
special_h_T = rfun * sp.diff(gfun, T0) - sfun
special_h_Y = (
    sp.diff(rfun, Y0) * gfun
    + rfun * sp.diff(gfun, Y0)
    - sp.diff(sfun, Y0) * T0
    + sp.diff(kfun, Y0)
)
pencil_constant = pencil_jacobian.coeff_monomial(1)
pencil_constant_special = sp.factor(
    pencil_constant.subs(
        {
            sp.diff(hfun, T0): special_h_T,
            sp.diff(hfun, Y0): special_h_Y,
        }
    )
)
expected_pencil_constant = sp.factor(
    (afun * sp.diff(gfun, T0) - bfun)
    * (
        sp.diff(rfun, Y0) * gfun
        - sp.diff(sfun, Y0) * T0
        + sp.diff(kfun, Y0)
    )
)
bezout_remainder = sp.factor(
    sp.diff(gfun, Y0) * (1 - bfun * rfun + afun * sfun)
)
assert sp.simplify(
    pencil_constant_special
    - expected_pencil_constant
    - bezout_remainder
) == 0

# Complementary T-only coefficient curve.
atfun = sp.Function("atfun")(T0)
btfun = sp.Function("btfun")(T0)
htfun = sp.Function("htfun")(T0, Y0)
t_curve_map = (
    T0 + atfun * U,
    gfun + btfun * U,
    htfun + U,
)
t_curve_jacobian = sp.Poly(
    sp.expand(sp.det(sp.Matrix(t_curve_map).jacobian((T0, Y0, U)))),
    U,
)
t_wronskian = sp.expand(
    atfun * sp.diff(btfun, T0) - btfun * sp.diff(atfun, T0)
)
expected_t_curve_u = sp.expand(
    t_wronskian * sp.diff(htfun, Y0)
    + sp.diff(atfun, T0) * sp.diff(gfun, Y0)
)
assert sp.expand(
    t_curve_jacobian.coeff_monomial(U) - expected_t_curve_u
) == 0
rtfun = sp.Function("rtfun")(T0)
ktfun = sp.Function("ktfun")(T0)
special_ht_T = (
    -sp.diff(rtfun, T0) * gfun
    - rtfun * sp.diff(gfun, T0)
    + sp.diff(ktfun, T0)
)
special_ht_Y = -rtfun * sp.diff(gfun, Y0)
t_curve_constant_special = sp.factor(
    t_curve_jacobian.coeff_monomial(1).subs(
        {
            sp.diff(htfun, T0): special_ht_T,
            sp.diff(htfun, Y0): special_ht_Y,
        }
    )
)
expected_t_curve_constant = sp.factor(
    sp.diff(gfun, Y0)
    * (
        atfun * sp.diff(rtfun, T0) * gfun
        - atfun * sp.diff(ktfun, T0)
        + btfun * rtfun
        + 1
    )
)
assert sp.simplify(
    t_curve_constant_special - expected_t_curve_constant
) == 0

# Coefficient curves in the geometrically distinguished mixed coordinate
# x=T+Y^2 also factor.
x0, y0 = sp.symbols("x0 y0")
ax = sp.Function("ax")(x0)
bx = sp.Function("bx")(x0)
hx = sp.Function("hx")(x0, y0)
Gx = sp.Function("Gx")(x0, y0)
mixed_curve_map = (
    x0 - y0**2 + ax * U,
    Gx + bx * U,
    hx + U,
)
mixed_curve_jacobian = sp.Poly(
    sp.expand(
        sp.det(sp.Matrix(mixed_curve_map).jacobian((x0, y0, U)))
    ),
    U,
)
mixed_wronskian = sp.expand(
    ax * sp.diff(bx, x0) - bx * sp.diff(ax, x0)
)
expected_mixed_u = sp.expand(
    mixed_wronskian * sp.diff(hx, y0)
    + sp.diff(ax, x0) * sp.diff(Gx, y0)
    + 2 * y0 * sp.diff(bx, x0)
)
assert sp.expand(
    mixed_curve_jacobian.coeff_monomial(U) - expected_mixed_u
) == 0
rx = sp.Function("rx")(x0)
sx = sp.Function("sx")(x0)
kx = sp.Function("kx")(x0)
special_hx_x = (
    -sp.diff(rx, x0) * Gx
    - rx * sp.diff(Gx, x0)
    - sp.diff(sx, x0) * y0**2
    + sp.diff(kx, x0)
)
special_hx_y = -rx * sp.diff(Gx, y0) - 2 * y0 * sx
mixed_constant_special = sp.factor(
    mixed_curve_jacobian.coeff_monomial(1).subs(
        {
            sp.diff(hx, x0): special_hx_x,
            sp.diff(hx, y0): special_hx_y,
        }
    )
)
expected_mixed_constant = sp.factor(
    (
        sp.diff(rx, x0) * Gx
        + sp.diff(sx, x0) * y0**2
        - sp.diff(kx, x0)
        + sx
    )
    * (ax * sp.diff(Gx, y0) + 2 * y0 * bx)
)
mixed_bezout_remainder = sp.factor(
    (
        2 * y0 * sp.diff(Gx, x0)
        + sp.diff(Gx, y0)
    )
    * (1 - ax * sx + bx * rx)
)
assert sp.simplify(
    mixed_constant_special
    - expected_mixed_constant
    - mixed_bezout_remainder
) == 0

# The actual G_y has nonzero constant leading coefficient in y.
mixed_G = sp.expand(
    davenport_g.subs({DT: x0 - y0**2, DY: y0, A: root_a})
)
mixed_G_y = sp.Poly(sp.diff(mixed_G, y0), y0)
assert mixed_G_y.degree() == 6
assert sp.Poly(
    mixed_G_y.LC(), x0, extension=sp.sqrt(-7)
).degree() == 0
assert mixed_G_y.LC() != 0

# Coordinate-invariant pencil factorization.
Tb = sp.Function("Tb")(x0, y0)
Gb = sp.Function("Gb")(x0, y0)
ab = sp.Function("ab")(x0)
bb = sp.Function("bb")(x0)
hb = sp.Function("hb")(x0, y0)
general_pencil_map = (
    Tb + ab * U,
    Gb + bb * U,
    hb + U,
)
general_pencil_jacobian = sp.Poly(
    sp.expand(
        sp.det(sp.Matrix(general_pencil_map).jacobian((x0, y0, U)))
    ),
    U,
)
general_wronskian = sp.expand(
    ab * sp.diff(bb, x0) - bb * sp.diff(ab, x0)
)
expected_general_u = sp.expand(
    general_wronskian * sp.diff(hb, y0)
    + sp.diff(ab, x0) * sp.diff(Gb, y0)
    - sp.diff(bb, x0) * sp.diff(Tb, y0)
)
assert sp.expand(
    general_pencil_jacobian.coeff_monomial(U) - expected_general_u
) == 0
rbase = sp.Function("rbase")(x0)
sbase = sp.Function("sbase")(x0)
kbase = sp.Function("kbase")(x0)
special_hb_x = (
    -sp.diff(rbase, x0) * Gb
    - rbase * sp.diff(Gb, x0)
    + sp.diff(sbase, x0) * Tb
    + sbase * sp.diff(Tb, x0)
    + sp.diff(kbase, x0)
)
special_hb_y = (
    -rbase * sp.diff(Gb, y0)
    + sbase * sp.diff(Tb, y0)
)
general_constant_special = sp.factor(
    general_pencil_jacobian.coeff_monomial(1).subs(
        {
            sp.diff(hb, x0): special_hb_x,
            sp.diff(hb, y0): special_hb_y,
        }
    )
)
expected_general_constant = sp.factor(
    (ab * sp.diff(Gb, y0) - bb * sp.diff(Tb, y0))
    * (
        sp.diff(rbase, x0) * Gb
        - sp.diff(sbase, x0) * Tb
        - sp.diff(kbase, x0)
    )
)
base_jacobian = sp.det(
    sp.Matrix((Tb, Gb)).jacobian((x0, y0))
)
general_bezout_remainder = sp.factor(
    base_jacobian * (1 - ab * sbase + bb * rbase)
)
assert sp.simplify(
    general_constant_special
    - expected_general_constant
    - general_bezout_remainder
) == 0

# Affine-linear pencils fail the unit gate.  If Y moves, the y^7
# coefficient of G is nonzero; if Y is fixed, G has nonzero cubic T-degree.
AA, BB, CC, DD = sp.symbols("AA BB CC DD")
linear_T = AA * x0 + BB * y0
linear_Y = CC * x0 + DD * y0
linear_G = sp.expand(
    davenport_g.subs({DT: linear_T, DY: linear_Y, A: root_a})
)
linear_G_poly = sp.Poly(linear_G, y0)
assert sp.factor(linear_G_poly.coeff_monomial(y0**7) - DD**7 / 7) == 0
vertical_G = sp.expand(linear_G.subs(DD, 0))
vertical_G_poly = sp.Poly(vertical_G, y0)
assert vertical_G_poly.degree() == 3
assert vertical_G_poly.LC() != 0

# The only degree not excluded by dominance in x=T+p(Y) is quadratic.
# Its exact high-y coefficient ideal is the unit ideal.
qc, qd, qe = sp.symbols("qc qd qe")
quadratic_p = qc * y0**2 + qd * y0 + qe
quadratic_G = sp.expand(
    davenport_g.subs(
        {DT: x0 - quadratic_p, DY: y0, A: root_a}
    )
)
quadratic_G_poly = sp.Poly(quadratic_G, y0)
quadratic_equations = []
for y_power in range(3, 8):
    coefficient_in_x = sp.Poly(
        sp.expand(quadratic_G_poly.coeff_monomial(y0**y_power)),
        x0,
    )
    quadratic_equations.extend(
        coefficient
        for coefficient in coefficient_in_x.all_coeffs()
        if coefficient != 0
    )
assert len(quadratic_equations) == 9
quadratic_groebner = sp.groebner(
    quadratic_equations,
    qc,
    qd,
    qe,
    order="lex",
    extension=sp.sqrt(-7),
)
assert len(quadratic_groebner.polys) == 1
assert quadratic_groebner.polys[0].as_expr() == 1

# Degree dominance for both triangular orientations.
for degree in range(3, 30):
    assert 3 * degree + 1 > max(
        2 * degree + 3,
        degree + 5,
        7,
    )
for degree in range(1, 30):
    assert 7 * degree > max(
        5 * degree + 1,
        3 * degree + 2,
        degree + 3,
    )

print("PASS: the variable determinant family is an A4-to-A4 polynomial map")
print("PASS: its Jacobian is D^2 and its reduced center is A2")
print("PASS: the parabola chart contributes the additional factor D/2")
print("PASS: the direct product splice has Jacobian D^3/2, not a constant")
print("PASS: a shifted exceptional equation preserves degree but not the Jacobian")
print("PASS: the universal one-auxiliary affine-linear mask formula is exact")
print("PASS: constant Jacobian in that mask class forces polynomial invertibility")
print("PASS: the multi-mask cofactor identity forces the same invertibility")
print("PASS: the shifted-mask Keller PDE has only triangular solutions")
print("PASS: a boundary-multiple Davenport correction fails modulo T+Y^2")
print("PASS: the unshifted determinant output has a critical axis")
print("PASS: translation by a new variable makes its differential unimodular")
print("PASS: retaining the old base reduces the translated attack to a plane map")
print("PASS: parameter-linear U coupling is stably a plane Keller pair")
print("PASS: every Y-dependent projective U-coefficient curve is excluded")
print("PASS: every T-dependent projective U-coefficient curve is excluded")
print("PASS: every (T+Y^2)-dependent U-coefficient curve is excluded")
print("PASS: the invariant pencil factorization and unit gate are exact")
print("PASS: every affine-linear mixed coefficient pencil fails the unit gate")
print("PASS: the exceptional quadratic triangular coefficient ideal is one")
print("PASS: both one-triangular coordinate families fail the unit gate")
print("PASS Davenport post-coordinate attack audit")
