#!/usr/bin/env python3
"""Verify global and cubic higher-degree HC4 pencil obstructions.

The checker covers eleven exact packets.

* A constant-Hessian pencil is a polynomial nilpotent, Hessian-integrable
  deformation after division by its unit base Hessian.
* A pencil direction of generic Hessian rank one reduces in every degree to
  the ternary singular-Hessian gate already used for quadratic directions.
* The genuinely moving, leading-rank-three cubic direction has a unique
  normal form.  Its last three determinant faces give a unit contradiction.
* In the residual constant-kernel cubic packet, the two tangent second
  fundamental forms synchronize to one rank-one direction; a fixed ruling is
  an HC2/JC2 endpoint.
* Homogeneous rulings in every degree and arbitrary rulings through degree
  three are fixed cylinders, so those residual subpackets are closed.
* The finite binary Schur faces and the degenerate fourth-power charts close
  every quartic border coefficient as well.
* The simple-root square and four repeated-root Schur charts close every
  quintic border coefficient whose leading binary form is not a fifth power.
* In every degree, a squarefree leading binary form has too many simple
  roots for any transverse coefficient, so the border is a fixed cylinder.
* The same holds on the generic discriminant stratum with one double root:
  an exact valuation gap forces vanishing at the double root as well.
* An all-degree root-valuation formula isolates every possible repeated-root
  resonance in the first transverse Schur face.
* In degree six, the valuation sieve and the full weighted Schur face close
  every non-pure binary leading form; only the pure sixth-power chart remains.

The small-rank Hessian normal forms and the homogeneous four-variable Hesse
theorem are external structural inputs; every determinant and moving-chart
identity used after those reductions is checked here exactly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_higher_degree_pencil_obstructions.json"
)


# A four-by-four constant determinant pencil is exactly a nilpotent relative
# endomorphism.  The determinant below records all four characteristic
# coefficients independently.
s = sp.symbols("s")
n = sp.symbols("n0:16")
N = sp.Matrix(4, 4, n)
relative_characteristic = sp.Poly((sp.eye(4) + s * N).det(), s)
assert relative_characteristic.degree() == 4
assert relative_characteristic.coeff_monomial(1) == 1
assert sp.factor(relative_characteristic.coeff_monomial(s) - sp.trace(N)) == 0
assert sp.factor(relative_characteristic.coeff_monomial(s**4) - N.det()) == 0


# Rank-one directions in every degree.  If Hess(A)=h''*e0*e0^T, the only
# pencil coefficient is h'' times the passive ternary Hessian determinant.
h2 = sp.symbols("h2")
q00 = sp.symbols("q00")
d0, d1, d2 = sp.symbols("d0 d1 d2")
e00, e01, e02, e11, e12, e22 = sp.symbols(
    "e00 e01 e02 e11 e12 e22"
)
E3 = sp.Matrix(
    [[e00, e01, e02], [e01, e11, e12], [e02, e12, e22]]
)
rank_one_base = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[q00]]), sp.Matrix([[d0, d1, d2]])),
    sp.Matrix.hstack(sp.Matrix([d0, d1, d2]), E3),
)
rank_one_direction = sp.diag(h2, 0, 0, 0)
rank_one_polynomial = sp.Poly(
    (rank_one_base + s * rank_one_direction).det(method="domain-ge"), s
)
assert sp.factor(
    rank_one_polynomial.coeff_monomial(s) - h2 * E3.det(method="domain-ge")
) == 0


# The unique genuinely moving leading-rank-three cubic direction is
# A=w*z+y*b(z)+G(x,z), with q=b'(z), m=b'' != 0.  At the second-jet level its
# Hessian is the matrix T below.  The coefficient called g is G_xx, c is
# G_xz, and d=G_zz+m*y.
g, c, d, q = sp.symbols("g c d q")
T = sp.Matrix(
    [
        [g, 0, c, 0],
        [0, 0, q, 0],
        [c, q, d, 1],
        [0, 0, 1, 0],
    ]
)
v = sp.Matrix([0, 1, 0, -q])
assert T.det(method="domain-ge") == 0
assert T * v == sp.zeros(4, 1)
assert T[:3, :3].det(method="domain-ge") == -g * q**2

xx, xy, xz, xw, yy, yz, yw, zz, zw, ww = sp.symbols(
    "xx xy xz xw yy yz yw zz zw ww"
)
S = sp.Matrix(
    [
        [xx, xy, xz, xw],
        [xy, yy, yz, yw],
        [xz, yz, zz, zw],
        [xw, yw, zw, ww],
    ]
)
cubic_polynomial = sp.Poly((S + s * T).det(method="domain-ge"), s)
null_form = sp.expand((v.T * S * v)[0])
assert sp.factor(
    cubic_polynomial.coeff_monomial(s**3) + g * null_form
) == 0

# Work modulo the leading null equation and use the invariant components of
# S*v.  They are X=(S*v)_x, U=(S*v)_w, and Z=(S*v)_z.
X, U, Z = sp.symbols("X U Z")
invariant_substitution = {
    yw: U + q * ww,
    xy: X + q * xw,
    yz: Z + q * zw,
    yy: 2 * q * (U + q * ww) - q**2 * ww,
}
second_face = sp.factor(
    cubic_polynomial.coeff_monomial(s**2).subs(invariant_substitution)
)
expected_second_face = sp.expand(
    (X - c * U) ** 2 + g * U * (2 * Z - d * U)
)
assert sp.factor(second_face - expected_second_face) == 0

first_face = sp.factor(
    cubic_polynomial.coeff_monomial(s).subs(invariant_substitution)
)
constant_face = sp.factor(
    cubic_polynomial.coeff_monomial(1).subs(invariant_substitution)
)
assert sp.factor(first_face.subs({U: 0, X: 0}) + g * Z**2 * ww) == 0
assert sp.factor(
    constant_face.subs({U: 0, X: 0})
    - Z**2 * (xw**2 - ww * xx)
) == 0


# Check the moving-chart integration identities on a degree-unbounded
# symbolic representative psi=y*C(x,z,r)+D(x,z,r), r=w+q(z)*y, with q'=m.
# Independent jet symbols suffice because only the displayed derivatives are
# used in the determinant faces.
m, y = sp.symbols("m y")
Cr, Cx, Cz, Dr = sp.symbols("Cr Cx Cz Dr")
integrated_U = Cr
integrated_X = Cx
integrated_Z = Cz + 2 * m * y * Cr + m * Dr
integrated_d = sp.symbols("Gzz") + m * y
integrated_second_face = sp.expand(
    expected_second_face.subs(
        {
            U: integrated_U,
            X: integrated_X,
            Z: integrated_Z,
            d: integrated_d,
        }
    )
)
assert sp.factor(sp.diff(integrated_second_face, y) - 3 * m * g * Cr**2) == 0


# Exact normal-form calibration with b=z^2/2 and G=x^2*z/2.  The Hessian
# direction has generic rank three and moving kernel (0,1,0,-z).
x, z, w = sp.symbols("x z w")
A_moving = w * z + y * z**2 / 2 + x**2 * z / 2
moving_hessian = sp.hessian(A_moving, (x, y, z, w))
assert moving_hessian.det(method="domain-ge") == 0
assert moving_hessian * sp.Matrix([0, 1, 0, -z]) == sp.zeros(4, 1)
assert moving_hessian[:3, :3].det(method="domain-ge") == -z**3


# The residual constant-kernel packet is
#   A=a(x,y,z), psi=w*c(x,y,z)+D(x,y,z).
# On the two-dimensional tangent plane p^perp, its determinant is a binary
# pencil det(E+s*F+w*G).  The three top coefficients put F and G in the
# null cone of Sym_2 and make them orthogonal for its polarization.  The
# squares of all 2-by-2 minors below lie in that coefficient ideal, proving
# that F and G are proportional over the fraction field.
fa, fb, fc, ga, gb, gc = sp.symbols("fa fb fc ga gb gc")
F2 = sp.Matrix([[fa, fb], [fb, fc]])
G2 = sp.Matrix([[ga, gb], [gb, gc]])
binary_mixed = sp.expand(fa * gc + fc * ga - 2 * fb * gb)
null_pair_ideal = [F2.det(), G2.det(), binary_mixed]
null_pair_groebner = sp.groebner(
    null_pair_ideal, fa, fb, fc, ga, gb, gc, order="grevlex"
)
proportional_minors = [
    fa * gb - fb * ga,
    fa * gc - fc * ga,
    fb * gc - fc * gb,
]
for minor in proportional_minors:
    assert null_pair_groebner.reduce(sp.expand(minor**2))[1] == 0

# The coefficient linear in a rank-one tangent form is the value of the base
# metric on its kernel.  This is the intrinsic null-ruling statement.
er, es, et, ur, ut = sp.symbols("er es et ur ut")
E2 = sp.Matrix([[er, es], [es, et]])
rank_one_form = sp.Matrix([[ur**2, ur * ut], [ur * ut, ut**2]])
rank_one_kernel = sp.Matrix([-ut, ur])
binary_rank_one_pencil = sp.Poly((E2 + s * rank_one_form).det(), s)
assert sp.factor(
    binary_rank_one_pencil.coeff_monomial(s)
    - (rank_one_kernel.T * E2 * rank_one_kernel)[0]
) == 0


# Fixed-cylinder closure.  Put c=c(x,y), p=(c_x,c_y,0), and let z be the
# fixed cylinder direction.  The coefficient linear in Hess(c) in the
# bordered determinant is D_zz times the curvature of the plane polynomial
# c.  If the curvature is nonzero, D is affine in z and the full potential is
# the cotangent form w*c+z*L+M, whose Hessian determinant is Jac(c,L)^2.
px, py = sp.symbols("px py")
cxx, cxy, cyy = sp.symbols("cxx cxy cyy")
dxx, dxy, dxz, dyy, dyz, dzz = sp.symbols(
    "dxx dxy dxz dyy dyz dzz"
)
p3 = sp.Matrix([px, py, 0])
Hc_cylinder = sp.Matrix(
    [[cxx, cxy, 0], [cxy, cyy, 0], [0, 0, 0]]
)
HD_general = sp.Matrix(
    [[dxx, dxy, dxz], [dxy, dyy, dyz], [dxz, dyz, dzz]]
)
cylinder_border = sp.Poly(
    (p3.T * (HD_general + s * Hc_cylinder).adjugate() * p3)[0], s
)
tangent_vector = sp.Matrix([-py, px, 0])
cylinder_curvature = sp.expand(
    (tangent_vector.T * Hc_cylinder * tangent_vector)[0]
)
assert sp.factor(
    cylinder_border.coeff_monomial(s) - dzz * cylinder_curvature
) == 0

lx, ly = sp.symbols("lx ly")
cxx0, cxy0, cyy0 = sp.symbols("cxx0 cxy0 cyy0")
lxx, lxy, lyy = sp.symbols("lxx lxy lyy")
mxx, mxy, myy = sp.symbols("mxx mxy myy")
cotangent_hessian = sp.Matrix(
    [
        [w * cxx0 + z * lxx + mxx, w * cxy0 + z * lxy + mxy, lx, px],
        [w * cxy0 + z * lxy + mxy, w * cyy0 + z * lyy + myy, ly, py],
        [lx, ly, 0, 0],
        [px, py, 0, 0],
    ]
)
plane_jacobian = sp.expand(lx * py - ly * px)
assert sp.factor(cotangent_hessian.det() - plane_jacobian**2) == 0


# The universal-field (bordered-Hessian) equation for c is exactly the
# singular-Hessian equation for tau*c.  For homogeneous c, the four-variable
# Hesse theorem supplies a constant kernel and therefore a fixed cylinder.
tau = sp.symbols("tau")
p0, p1, p2 = sp.symbols("p0 p1 p2")
h00, h01, h02, h11, h12, h22 = sp.symbols(
    "h00 h01 h02 h11 h12 h22"
)
p_generic = sp.Matrix([p0, p1, p2])
H_generic = sp.Matrix(
    [[h00, h01, h02], [h01, h11, h12], [h02, h12, h22]]
)
hessian_tau_c = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(1, 1), p_generic.T),
    sp.Matrix.hstack(p_generic, tau * H_generic),
)
generic_border = sp.expand((p_generic.T * H_generic.adjugate() * p_generic)[0])
assert sp.factor(hessian_tau_c.det() + tau**2 * generic_border) == 0

# Euler's identities give the leading homogeneous reduction
# p^T adj(H) p = d/(d-1) * c * det(H).
degree = sp.symbols("degree", nonzero=True)
u0, u1, u2 = sp.symbols("u0 u1 u2")
u_generic = sp.Matrix([u0, u1, u2])
p_euler = H_generic * u_generic / (degree - 1)
c_euler = (u_generic.T * H_generic * u_generic)[0] / (
    degree * (degree - 1)
)
assert sp.factor(
    (p_euler.T * H_generic.adjugate() * p_euler)[0]
    - degree * c_euler * H_generic.det() / (degree - 1)
) == 0


# Exact degree-three classification after the leading ternary Hesse theorem.
# A binary cubic has three root-partition charts over an algebraic closure.
# For the distinct and double-root charts, the top two bordered faces kill
# every quadratic z term; the remaining linear z coefficient is then killed
# by the nonzero binary Hessian.  The cube chart has two passive-rank charts.
A0, B0, C0, D0, E0, F0, L0, M0, N0 = sp.symbols(
    "A0 B0 C0 D0 E0 F0 L0 M0 N0"
)
quadratic_tail = (
    A0 * x**2
    + 2 * B0 * x * y
    + 2 * C0 * x * z
    + D0 * y**2
    + 2 * E0 * y * z
    + F0 * z**2
) / 2
linear_tail = L0 * x + M0 * y + N0 * z


def bordered_polynomial(poly: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(poly, variable) for variable in (x, y, z)])
    hessian = sp.hessian(poly, (x, y, z))
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


for cubic_chart in (x * y * (x - y), x**2 * y):
    chart_border = sp.Poly(
        bordered_polynomial(cubic_chart + quadratic_tail + linear_tail),
        x,
        y,
        z,
    )
    degree_five = sp.factor(
        sum(
            coefficient * x ** monomial[0] * y ** monomial[1] * z ** monomial[2]
            for monomial, coefficient in chart_border.terms()
            if sum(monomial) == 5
        )
    )
    if cubic_chart == x * y * (x - y):
        expected_degree_five = -6 * F0 * x * y * (x - y) * (
            x**2 - x * y + y**2
        )
    else:
        expected_degree_five = -6 * F0 * x**4 * y
    assert sp.factor(degree_five - expected_degree_five) == 0
    degree_four_after_f0 = sp.factor(
        sum(
            coefficient.subs(F0, 0)
            * x ** monomial[0]
            * y ** monomial[1]
            * z ** monomial[2]
            for monomial, coefficient in chart_border.terms()
            if sum(monomial) == 4
        )
    )
    if cubic_chart == x * y * (x - y):
        expected_square = -(
            C0 * x**2 - 2 * C0 * x * y - 2 * E0 * x * y + E0 * y**2
        ) ** 2
    else:
        expected_square = -x**2 * (C0 * x - 2 * E0 * y) ** 2
    assert sp.factor(degree_four_after_f0 - expected_square) == 0
    cylinder_chart = sp.expand(
        (cubic_chart + quadratic_tail + linear_tail).subs(
            {F0: 0, C0: 0, E0: 0}
        )
    )
    binary_hessian = sp.hessian(cylinder_chart, (x, y)).det()
    assert sp.factor(
        bordered_polynomial(cylinder_chart) - N0**2 * binary_hessian
    ) == 0

# Triple-root chart, passive quadratic rank zero.
cube_rank_zero = x**3 / 3 + (
    A0 * x**2 + 2 * B0 * x * y + 2 * C0 * x * z
) / 2 + linear_tail
assert sp.factor(
    bordered_polynomial(cube_rank_zero) + (B0 * N0 - C0 * M0) ** 2
) == 0

# Triple-root chart, passive quadratic rank one, normalized to y^2/2.
cube_rank_one = cube_rank_zero + y**2 / 2
cube_rank_one_border = sp.Poly(bordered_polynomial(cube_rank_one), x, y, z)
assert sp.factor(cube_rank_one_border.coeff_monomial(y**2) + C0**2) == 0
assert sp.factor(
    cube_rank_one_border.coeff_monomial(x).subs(C0, 0) - 2 * N0**2
) == 0


# Degree four.  When the binary quartic top is not a fourth power, the first
# two transverse faces reduce to a finite binary Schur equation
#   q*det(B_f) = b_g^T*adj(B_f)*b_g,
# where g is binary quadratic and q=c_{2,zz}.  Its four root-partition charts
# are exact and the sole apparent exception is f=x^2*y^2, g=b*x*y,
# q=b^2/4.
qa, qb, qc, qq, qlam = sp.symbols("qa qb qc qq qlam")
quartic_g = qa * x**2 + qb * x * y + qc * y**2


def binary_bordered_matrix(poly: sp.Expr) -> sp.Matrix:
    gradient = sp.Matrix([sp.diff(poly, x), sp.diff(poly, y)])
    hessian = sp.hessian(poly, (x, y))
    return sp.Matrix.vstack(
        sp.Matrix.hstack(sp.zeros(1, 1), gradient.T),
        sp.Matrix.hstack(gradient, hessian),
    )


quartic_root_charts = {
    "1111": x * y * (x - y) * (x - qlam * y),
    "211": x**2 * y * (x - y),
    "22": x**2 * y**2,
    "31": x**3 * y,
}
quartic_schur_groebners: dict[str, sp.GroebnerBasis] = {}
for chart_name, quartic_top in quartic_root_charts.items():
    bordered_top = binary_bordered_matrix(quartic_top)
    border_column = sp.Matrix(
        [quartic_g, sp.diff(quartic_g, x), sp.diff(quartic_g, y)]
    )
    schur_face = sp.expand(
        qq * bordered_top.det()
        - (border_column.T * bordered_top.adjugate() * border_column)[0]
    )
    schur_equations = sp.Poly(schur_face, x, y).coeffs()
    quartic_schur_groebners[chart_name] = sp.groebner(
        schur_equations, qa, qb, qc, qq, order="grevlex"
    )

for chart_name in ("1111", "211", "31"):
    chart_basis = quartic_schur_groebners[chart_name]
    assert chart_basis.reduce(qa**2)[1] == 0
    assert chart_basis.reduce(qc**2)[1] == 0
    assert chart_basis.reduce(qq**2)[1] == 0
    assert chart_basis.reduce(qb**4)[1] == 0

double_double_basis = quartic_schur_groebners["22"]
assert double_double_basis.reduce(qa**2)[1] == 0
assert double_double_basis.reduce(qc**2)[1] == 0
assert double_double_basis.reduce((4 * qq - qb**2) ** 2)[1] == 0

# Once g=q=0, the same binary bordered matrices kill the linear and constant
# z tails.  The displayed chart factors make this a radical statement over
# an algebraic closure (with qlam*(qlam-1) nonzero in the 1111 chart).
qu, qv, qn = sp.symbols("qu qv qn")
for chart_name, quartic_top in quartic_root_charts.items():
    bordered_top = binary_bordered_matrix(quartic_top)
    linear_tail = qu * x + qv * y
    linear_column = sp.Matrix([linear_tail, qu, qv])
    linear_face = sp.Poly(
        sp.expand((linear_column.T * bordered_top.adjugate() * linear_column)[0]),
        x,
        y,
    )
    if chart_name == "1111":
        assert linear_face.coeff_monomial(x**6) == -4 * qu**2
        assert sp.factor(
            linear_face.coeff_monomial(y**6) + 4 * qlam**2 * qv**2
        ) == 0
    elif chart_name == "211":
        assert linear_face.coeff_monomial(x**6) == -4 * qu**2
        assert linear_face.coeff_monomial(x**2 * y**4) == -8 * qv**2
    elif chart_name == "22":
        assert linear_face.coeff_monomial(x**4 * y**2) == -8 * qu**2
        assert linear_face.coeff_monomial(x**2 * y**4) == -8 * qv**2
    else:
        assert linear_face.coeff_monomial(x**6) == -4 * qu**2
        assert linear_face.coeff_monomial(x**4 * y**2) == -12 * qv**2
    constant_column = sp.Matrix([qn, 0, 0])
    constant_face = sp.factor(
        (constant_column.T * bordered_top.adjugate() * constant_column)[0]
    )
    assert constant_face != 0 and constant_face.has(qn**2)

# The apparent 22 exception dies on the very next face.  Include arbitrary
# binary cubic/quadratic tails and arbitrary lower z-linear tails to verify
# that the decisive coefficient cannot be canceled.
eb = sp.symbols("eb")
er0, er1, er2, er3 = sp.symbols("er0:4")
eA, eB, eC, eU, eV, eL, eM, eN = sp.symbols(
    "eA eB eC eU eV eL eM eN"
)
double_double_candidate = (
    x**2 * y**2
    + er0 * x**3
    + er1 * x**2 * y
    + er2 * x * y**2
    + er3 * y**3
    + eb * x * y * z
    + (eA * x**2 + 2 * eB * x * y + eC * y**2) / 2
    + z * (eU * x + eV * y)
    + eb**2 * z**2 / 8
    + eL * x
    + eM * y
    + eN * z
)
double_double_border = sp.Poly(
    bordered_polynomial(double_double_candidate), x, y, z
)
assert sp.factor(
    double_double_border.coeff_monomial(x**2 * y**2 * z**2)
    + 3 * eb**4 / 4
) == 0


# Pure fourth-power top.  The first face says that the passive cubic P3 is a
# cube (or zero), and synchronizes the x*Q2 term with it.  The six equations
# below are the complete coefficient list.
pa0, pa1, pa2, pa3, pa4, pa5, pa6 = sp.symbols("pa0:7")
passive_cubic = (
    pa0 * z**3 + pa1 * y * z**2 + pa2 * y**2 * z + pa3 * y**3
)
x_passive_quadratic = x * (pa4 * z**2 + pa5 * y * z + pa6 * y**2)
pure_top_test = x**4 / 12 + passive_cubic + x_passive_quadratic
pure_top_border = sp.Poly(bordered_polynomial(pure_top_test), x, y, z)
pure_top_degree_eight = sp.factor(
    sum(
        coefficient * x ** monomial[0] * y ** monomial[1] * z ** monomial[2]
        for monomial, coefficient in pure_top_border.terms()
        if sum(monomial) == 8
    )
)
pure_top_quadratic = (
    12 * pa0 * pa2 * z**2
    + 36 * pa0 * pa3 * y * z
    + 12 * pa0 * pa6 * x * z
    - 4 * pa1**2 * z**2
    - 4 * pa1 * pa2 * y * z
    + 12 * pa1 * pa3 * y**2
    - 4 * pa1 * pa5 * x * z
    + 4 * pa1 * pa6 * x * y
    - 4 * pa2**2 * y**2
    + 4 * pa2 * pa4 * x * z
    - 4 * pa2 * pa5 * x * y
    + 12 * pa3 * pa4 * x * y
    + 4 * pa4 * pa6 * x**2
    - pa5**2 * x**2
)
assert sp.factor(pure_top_degree_eight - x**6 * pure_top_quadratic / 9) == 0

# Pure top, nonzero passive cube P3=y^3/3.  The next two faces force the
# x^2*z direction and the yz/z^2 quadratic coefficients to vanish.
ra, ru, rv, rk = sp.symbols("ra ru rv rk")
rA, rB, rC, rD, rE, rF, rL, rM, rN = sp.symbols(
    "rA rB rC rD rE rF rL rM rN"
)
passive_cube_candidate = (
    x**4 / 12
    + y**3 / 3
    + ra * x * y**2
    + x**2 * (ru * y + rv * z)
    + rk * x**3
    + (
        rA * x**2
        + 2 * rB * x * y
        + 2 * rC * x * z
        + rD * y**2
        + 2 * rE * y * z
        + rF * z**2
    )
    / 2
    + rL * x
    + rM * y
    + rN * z
)
passive_cube_border = sp.Poly(
    bordered_polynomial(passive_cube_candidate), x, y, z
)
passive_cube_degree_seven = sp.factor(
    sum(
        coefficient * x ** monomial[0] * y ** monomial[1] * z ** monomial[2]
        for monomial, coefficient in passive_cube_border.terms()
        if sum(monomial) == 7
    )
)
assert sp.factor(
    passive_cube_degree_seven
    + 2 * x**6 * (ra * x + y) * (-rF + 3 * rv**2) / 9
) == 0
passive_cube_degree_six = sp.factor(
    sum(
        coefficient.subs(rF, 3 * rv**2)
        * x ** monomial[0]
        * y ** monomial[1]
        * z ** monomial[2]
        for monomial, coefficient in passive_cube_border.terms()
        if sum(monomial) == 6
    )
)
passive_cube_square = -x**2 * (
    -rE * x**2
    + 6 * ra * rv * x * y
    + 3 * ru * rv * x**2
    + 3 * rv * y**2
) ** 2 / 9
assert sp.factor(passive_cube_degree_six - passive_cube_square) == 0
passive_cube_reduced = sp.Poly(
    bordered_polynomial(
        passive_cube_candidate.subs({rv: 0, rE: 0, rF: 0})
    ),
    x,
    y,
    z,
)
assert sp.factor(passive_cube_reduced.coeff_monomial(y**4) + rC**2) == 0

# Pure top, zero P3 and nonzero square Q2=ra*y^2.  The analogous faces force
# the same fixed cylinder.
square_q_candidate = passive_cube_candidate - y**3 / 3
square_q_border = sp.Poly(bordered_polynomial(square_q_candidate), x, y, z)
square_q_degree_seven = sp.factor(
    sum(
        coefficient * x ** monomial[0] * y ** monomial[1] * z ** monomial[2]
        for monomial, coefficient in square_q_border.terms()
        if sum(monomial) == 7
    )
)
assert sp.factor(
    square_q_degree_seven + 2 * ra * x**7 * (-rF + 3 * rv**2) / 9
) == 0
square_q_degree_six = sp.factor(
    sum(
        coefficient.subs(rF, 3 * rv**2)
        * x ** monomial[0]
        * y ** monomial[1]
        * z ** monomial[2]
        for monomial, coefficient in square_q_border.terms()
        if sum(monomial) == 6
    )
)
square_q_expected = -x**4 * (
    -rE * x + 6 * ra * rv * y + 3 * ru * rv * x
) ** 2 / 9
assert sp.factor(square_q_degree_six - square_q_expected) == 0
square_q_reduced = sp.Poly(
    bordered_polynomial(square_q_candidate.subs({rv: 0, rE: 0, rF: 0})),
    x,
    y,
    z,
)
assert sp.factor(
    square_q_reduced.coeff_monomial(x**5) - 2 * rC**2 * ra / 3
) == 0

# Pure top, Q2=0 and nonzero x^2*y.  The passive quadratic relation has one
# apparent full-rank chart, but its next face contains 9*rF*x^2*y^2.
line_cubic_candidate = (
    x**4 / 12
    + x**2 * y
    + rk * x**3
    + (
        rA * x**2
        + 2 * rB * x * y
        + 2 * rC * x * z
        + rD * y**2
        + 2 * rE * y * z
        + rF * z**2
    )
    / 2
    + rL * x
    + rM * y
    + rN * z
)
line_cubic_border = sp.Poly(
    bordered_polynomial(line_cubic_candidate), x, y, z
)
assert sp.factor(
    line_cubic_border.coeff_monomial(x**6)
    - (rD * rF - rE**2 - 3 * rF) / 9
) == 0
line_full_rank = sp.Poly(
    bordered_polynomial(line_cubic_candidate.subs({rD: 3, rE: 0})),
    x,
    y,
    z,
)
assert line_full_rank.coeff_monomial(x**2 * y**2) == 9 * rF
line_rank_one = sp.Poly(
    bordered_polynomial(line_cubic_candidate.subs({rE: 0, rF: 0})),
    x,
    y,
    z,
)
assert sp.factor(
    line_rank_one.coeff_monomial(x**4) - rC**2 * (rD - 3) / 3
) == 0
line_resonant = sp.Poly(
    bordered_polynomial(
        line_cubic_candidate.subs({rD: 3, rE: 0, rF: 0})
    ),
    x,
    y,
    z,
)
assert line_resonant.coeff_monomial(y**2) == -9 * rC**2
assert sp.factor(line_resonant.coeff_monomial(y).subs(rC, 0) - 6 * rN**2) == 0

# Finally, if the cubic correction is absent, the passive quadratic has rank
# at most one.  Its rank-zero and rank-one normalizations are the same fixed
# cylinder alternatives already seen in degree three.
pure_quadratic_rank_zero = (
    x**4 / 12
    + (rA * x**2 + 2 * rB * x * y + 2 * rC * x * z) / 2
    + rL * x
    + rM * y
    + rN * z
)
assert sp.factor(
    bordered_polynomial(pure_quadratic_rank_zero)
    + (rB * rN - rC * rM) ** 2
) == 0
pure_quadratic_rank_one = pure_quadratic_rank_zero + y**2 / 2
pure_quadratic_rank_one_border = sp.Poly(
    bordered_polynomial(pure_quadratic_rank_one), x, y, z
)
assert sp.factor(
    pure_quadratic_rank_one_border.coeff_monomial(x**4) - rC**2 / 3
) == 0
assert sp.factor(
    pure_quadratic_rank_one_border.coeff_monomial(x**2).subs(rC, 0)
    - rN**2
) == 0


# Degree five with a non-pure binary top.  At a simple root x=0, write
# f=x*h.  The Schur numerator restricts to the square below.  For a
# transverse cubic g, homogeneity turns it into a nonzero scalar times g^2,
# so every simple root of f is also a root of g.
sh, sha, shb, sg, sgx, sgy = sp.symbols("sh sha shb sg sgx sgy")
simple_root_border = sp.Matrix(
    [[0, sh, 0], [sh, sha, shb], [0, shb, 0]]
)
simple_root_column = sp.Matrix([sg, sgx, sgy])
assert sp.factor(
    (simple_root_column.T * simple_root_border.adjugate() * simple_root_column)[0]
    + (shb * sg - sh * sgy) ** 2
) == 0

# At a double root, use the affine chart y=1 and Euler to reconstruct the
# homogeneous binary derivatives from F(x)=f(x,1) and G(x)=g(x,1).  The
# bordered determinant begins in order four while the Schur numerator begins
# in order two.  Its scalar is positive for integer d>=3 and e<=d-2:
# writing k=d-1-e>=1 gives
#
#   E=2*(k-d/2)^2+d*(d-2)/2 != 0.
ld, le = sp.symbols("ld le")
lh0, lh1, lh2, lg0, lg1, lg2 = sp.symbols(
    "lh0 lh1 lh2 lg0 lg1 lg2"
)
local_f = x**2 * (lh0 + lh1 * x + lh2 * x**2)
local_g = lg0 + lg1 * x + lg2 * x**2
local_fx = sp.diff(local_f, x)
local_fy = ld * local_f - x * local_fx
local_fxx = sp.diff(local_f, x, 2)
local_fxy = (ld - 1) * local_fx - x * local_fxx
local_fyy = (
    ld * (ld - 1) * local_f
    - 2 * (ld - 1) * x * local_fx
    + x**2 * local_fxx
)
local_gx = sp.diff(local_g, x)
local_gy = le * local_g - x * local_gx
double_root_border = sp.Matrix(
    [
        [0, local_fx, local_fy],
        [local_fx, local_fxx, local_fxy],
        [local_fy, local_fxy, local_fyy],
    ]
)
double_root_column = sp.Matrix([local_g, local_gx, local_gy])
double_root_det = sp.Poly(sp.expand(double_root_border.det()), x)
double_root_numerator = sp.Poly(
    sp.expand(
        (double_root_column.T * double_root_border.adjugate() * double_root_column)[
            0
        ]
    ),
    x,
)
double_root_scalar = (
    ld**2 - 2 * ld * le - 3 * ld + 2 * le**2 + 4 * le + 2
)
assert sp.factor(
    double_root_det.coeff_monomial(x**4)
    - 2 * ld * (ld - 2) * lh0**3
) == 0
assert sp.factor(
    double_root_numerator.coeff_monomial(x**2)
    + 2 * lg0**2 * lh0**2 * double_root_scalar
) == 0

# The remaining non-pure root partitions reduce to four exact Schur charts.
vA, vB, vC, vD, vQ, vR = sp.symbols("vA vB vC vD vQ vR")
quintic_g = vA * x**3 + vB * x**2 * y + vC * x * y**2 + vD * y**3
quintic_q = vQ * x + vR * y
quintic_root_charts = {
    "221": x**2 * y**2 * (x - y),
    "311": x**3 * y * (x - y),
    "32": x**3 * y**2,
    "41": x**4 * y,
}
quintic_schur_groebners: dict[str, sp.GroebnerBasis] = {}
for chart_name, quintic_top in quintic_root_charts.items():
    bordered_top = binary_bordered_matrix(quintic_top)
    border_column = sp.Matrix(
        [quintic_g, sp.diff(quintic_g, x), sp.diff(quintic_g, y)]
    )
    schur_face = sp.expand(
        quintic_q * bordered_top.det()
        - (border_column.T * bordered_top.adjugate() * border_column)[0]
    )
    quintic_schur_groebners[chart_name] = sp.groebner(
        sp.Poly(schur_face, x, y).coeffs(),
        vA,
        vB,
        vC,
        vD,
        vQ,
        vR,
        order="grevlex",
    )


def assert_radical_member(
    basis: sp.GroebnerBasis, polynomial: sp.Expr, maximum_power: int = 6
) -> None:
    assert any(
        basis.reduce(sp.expand(polynomial**power))[1] == 0
        for power in range(1, maximum_power + 1)
    )


for chart_name in ("221", "311"):
    chart_basis = quintic_schur_groebners[chart_name]
    for coefficient in (vA, vB, vC, vD, vQ, vR):
        assert_radical_member(chart_basis, coefficient)

triple_double_basis = quintic_schur_groebners["32"]
for coefficient in (vA, vC, vD, vR):
    assert_radical_member(triple_double_basis, coefficient)
assert_radical_member(triple_double_basis, 30 * vQ - 11 * vB**2)

quadruple_simple_basis = quintic_schur_groebners["41"]
for coefficient in (vA, vC, vD, vQ):
    assert_radical_member(quadruple_simple_basis, coefficient)
assert_radical_member(quadruple_simple_basis, 5 * vR - vB**2)

# Both apparent quintic exceptions die one face later, even after arbitrary
# binary quartic corrections, arbitrary next z-linear quadratics, and the
# available scalar z^2 correction are included.
vb = sp.symbols("vb")
vr = sp.symbols("vr0:5")
vh = sp.symbols("vh0:3")
vq0 = sp.symbols("vq0")
binary_quartic_tail = sum(vr[i] * x ** (4 - i) * y**i for i in range(5))
binary_quadratic_tail = vh[0] * x**2 + vh[1] * x * y + vh[2] * y**2
quintic_exception_data = {
    "32": (
        x**3 * y**2,
        vb * x**2 * y,
        sp.Rational(11, 30) * vb**2 * x,
        x**7 * y**3 * z,
        -sp.Rational(8, 15) * vb**3,
    ),
    "41": (
        x**4 * y,
        vb * x**2 * y,
        sp.Rational(1, 5) * vb**2 * y,
        x**8 * y**2 * z,
        sp.Rational(6, 5) * vb**3,
    ),
}
for quintic_top, transverse_cubic, forced_q, monomial, expected in (
    quintic_exception_data.values()
):
    exception_candidate = (
        quintic_top
        + binary_quartic_tail
        + z * transverse_cubic
        + z * binary_quadratic_tail
        + z**2 * forced_q / 2
        + vq0 * z**2 / 2
    )
    exception_border = sp.Poly(
        bordered_polynomial(exception_candidate), x, y, z
    )
    assert sp.factor(exception_border.coeff_monomial(monomial) - expected) == 0

# After those cubic exceptions vanish, a possible next transverse quadratic
# must be a multiple of x*y in the 32 and 41 charts.  Its leading face is a
# pure square, so it vanishes as well.
va = sp.symbols("va")
vlow4 = sp.symbols("vlow4_0:5")
vlow3 = sp.symbols("vlow3_0:4")
v2a, v2b, v2c, v1u, v1v, v0x, v0y, v0z = sp.symbols(
    "v2a v2b v2c v1u v1v v0x v0y v0z"
)
lower_binary_four = sum(
    vlow4[i] * x ** (4 - i) * y**i for i in range(5)
)
lower_binary_three = sum(
    vlow3[i] * x ** (3 - i) * y**i for i in range(4)
)
lower_binary_two = (v2a * x**2 + 2 * v2b * x * y + v2c * y**2) / 2
for chart_name, quintic_top, monomial, expected in (
    ("32", x**3 * y**2, x**6 * y**4, -va**2),
    ("41", x**4 * y, x**8 * y**2, -9 * va**2),
):
    secondary_candidate = (
        quintic_top
        + lower_binary_four
        + lower_binary_three
        + va * x * y * z
        + lower_binary_two
        + z * (v1u * x + v1v * y)
        + v0x * x
        + v0y * y
        + v0z * z
    )
    secondary_border = sp.Poly(
        bordered_polynomial(secondary_candidate), x, y, z
    )
    assert sp.factor(secondary_border.coeff_monomial(monomial) - expected) == 0


# An arbitrary root of multiplicity m gives an all-degree valuation sieve.
# After dehomogenizing at the root, the leading coefficients of the bordered
# determinant and Schur numerator are obtained from these two coefficient
# matrices.  If g has root order n<m/2 and C(d,m,e,n) is nonzero, then the
# numerator order 2*n+2*m-2 is smaller than the determinant order 3*m-2,
# contradicting polynomial divisibility.
vd, vm, ve, vn = sp.symbols("vd vm ve vn")
root_border_coefficient = sp.Matrix(
    [
        [0, vm, vd - vm],
        [vm, vm * (vm - 1), vm * (vd - vm)],
        [
            vd - vm,
            vm * (vd - vm),
            (vd - vm) * (vd - vm - 1),
        ],
    ]
)
root_column_coefficient = sp.Matrix([1, vn, ve - vn])
root_resonance = (
    vd**2 * vm
    + vd**2 * vn**2
    - 2 * vd * ve * vm * vn
    - 2 * vd * ve * vm
    - vd * vm**2
    - vd * vm
    + ve**2 * vm**2
    + 2 * ve * vm**2
    + vm**2
)
assert sp.factor(
    root_border_coefficient.det() - vd * vm * (vd - vm)
) == 0
assert sp.factor(
    (root_column_coefficient.T
     * root_border_coefficient.adjugate()
     * root_column_coefficient)[0]
    + root_resonance
) == 0


def root_weight(degree: int, multiplicity: int, transverse_degree: int) -> int:
    """Return the first root order not excluded by the valuation face."""

    half_ceiling = (multiplicity + 1) // 2
    for order in range(half_ceiling):
        scalar = root_resonance.subs(
            {
                vd: degree,
                vm: multiplicity,
                ve: transverse_degree,
                vn: order,
            }
        )
        if scalar == 0:
            return order
    return half_ceiling


# For the first layer e=d-2 the scalar has no resonance when 1<=m<d and
# d>=4.  On n<m/2 it is decreasing in n.  Its last possible value is the
# displayed square for odd m; for even m it is a sum of positive terms after
# writing d=m+a.  The finite assertions below calibrate the formulas used in
# the proof.
odd_boundary = sp.factor(
    4
    * root_resonance.subs(
        {ve: vd - 2, vn: (vm - 1) / 2}
    )
)
assert odd_boundary == (vd * vm - vd - 2 * vm) ** 2
va_positive = sp.symbols("va_positive")
even_boundary = sp.factor(
    4
    * root_resonance.subs(
        {ve: vd - 2, vn: vm / 2 - 1}
    ).subs(vd, vm + va_positive)
)
assert sp.expand(even_boundary) == sp.expand(
    va_positive**2 * (vm**2 + 4)
    + 2 * va_positive * vm * (vm**2 - 2 * vm + 2)
    + vm**2 * (vm - 2) ** 2
)


def binary_schur_face(
    top: sp.Expr, transverse: sp.Expr, passive_quadratic: sp.Expr
) -> sp.Expr:
    border = binary_bordered_matrix(top)
    column = sp.Matrix(
        [transverse, sp.diff(transverse, x), sp.diff(transverse, y)]
    )
    return sp.expand(
        passive_quadratic * border.det()
        - (column.T * border.adjugate() * column)[0]
    )


def coefficient_ideal(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> sp.GroebnerBasis:
    coefficients = sp.Poly(sp.expand(polynomial), x, y, z).coeffs()
    return sp.groebner(coefficients, *variables, order="lex")


# Degree six.  RSD25--RSD26 already remove 111111 and 21111, while the
# valuation weights remove 3111.  The seven remaining non-pure partitions
# have the following weight sums for e=0,...,4.  Consequently the explicit
# candidates checked below exhaust every first nonbinary weighted face.
sextic_partitions = {
    "2211": (2, 2, 1, 1),
    "411": (4, 1, 1),
    "321": (3, 2, 1),
    "51": (5, 1),
    "222": (2, 2, 2),
    "42": (4, 2),
    "33": (3, 3),
}
expected_weight_sums = {
    "2211": (4, 4, 4, 4, 4),
    "411": (4, 4, 4, 4, 4),
    "321": (4, 4, 4, 3, 4),
    "51": (4, 4, 4, 4, 4),
    "222": (3, 3, 3, 3, 3),
    "42": (3, 3, 3, 3, 3),
    "33": (4, 4, 4, 2, 4),
}
for partition_name, multiplicities in sextic_partitions.items():
    assert tuple(
        sum(root_weight(6, multiplicity, degree) for multiplicity in multiplicities)
        for degree in range(5)
    ) == expected_weight_sums[partition_name]

# The 3111 partition is already globally impossible: its weight is five at
# e=4, four at the sole e=3 resonance, and five at every lower degree.
assert tuple(
    sum(root_weight(6, multiplicity, degree) for multiplicity in (3, 1, 1, 1))
    for degree in range(5)
) == (5, 5, 5, 4, 5)

s6a, s6u, s6v, s6A, s6B, s6C, s6b, s6Q, s6r = sp.symbols(
    "s6a s6u s6v s6A s6B s6C s6b s6Q s6r"
)
s6_quadratic = s6A * x**2 + s6B * x * y + s6C * y**2

# In 2211 the local sieve leaves the squarefree radical times a scalar.
# Keeping the cross-ratio s6a symbolic, even the constant Schur equation has
# unit coefficient ideal.  The same happens in the 411 and 321 top layers.
sextic_no_schur_solution = {
    "2211": (
        x**2 * y**2 * (x - y) * (x - s6a * y),
        x * y * (x - y) * (x - s6a * y),
        (s6A, s6B, s6C, s6a),
    ),
    "411": (x**4 * y * (x - y), x**2 * y * (x - y), (s6A, s6B, s6C)),
    "321": (x**3 * y**2 * (x - y), x**2 * y * (x - y), (s6A, s6B, s6C)),
}
for sextic_top, sextic_g, parameters in sextic_no_schur_solution.values():
    face = binary_schur_face(sextic_top, sextic_g, s6_quadratic)
    basis = sp.groebner(sp.Poly(face, x, y).coeffs(), *parameters, order="lex")
    assert len(basis.polys) == 1 and basis.polys[0].as_expr() == 1

# The 222 top layer also dies in the constant Schur equation.  Its locally
# admissible quartic is radical*(u*x+v*y); the radical of the face ideal is
# (u,v), independently of the passive quadratic q.
top_222 = x**2 * y**2 * (x - y) ** 2
g_222_e4 = x * y * (x - y) * (s6u * x + s6v * y)
basis_222_e4 = sp.groebner(
    sp.Poly(binary_schur_face(top_222, g_222_e4, s6_quadratic), x, y).coeffs(),
    s6A,
    s6B,
    s6C,
    s6u,
    s6v,
    order="lex",
)
assert_radical_member(basis_222_e4, s6u)
assert_radical_member(basis_222_e4, s6v)

# Four apparent Schur solutions remain.  They are weighted-homogeneous
# potentials f+z*g+z^2*q/2 plus every higher z-power of the same weight, so
# the whole bordered polynomial, not only its z^0 coefficient, is an initial
# face.  At e=4 the scalar z^3 coefficient is retained.  Its coefficient
# ideal has radical equal to all transverse parameters in every case.
weighted_sextic_faces = {
    "51_e4": (
        x**5 * y,
        s6b * x**3 * y,
        sp.Rational(11, 30) * s6b**2 * x * y,
        s6r * z**3,
        (s6b, s6r),
    ),
    "33_e4": (
        x**3 * y**3,
        s6b * x**2 * y**2,
        sp.Rational(1, 2) * s6b**2 * x * y,
        s6r * z**3,
        (s6b, s6r),
    ),
    "222_e3": (
        top_222,
        s6b * x * y * (x - y),
        sp.Rational(1, 6) * s6b**2,
        0,
        (s6b,),
    ),
    "42_e3": (
        x**4 * y**2,
        s6b * x**2 * y,
        sp.Rational(1, 6) * s6b**2,
        0,
        (s6b,),
    ),
    "33_e3": (
        x**3 * y**3,
        x * y * (s6u * x + s6v * y),
        sp.Rational(2, 3) * s6u * s6v,
        0,
        (s6u, s6v),
    ),
}
for sextic_top, sextic_g, sextic_q, higher_tail, parameters in weighted_sextic_faces.values():
    initial_potential = sextic_top + z * sextic_g + z**2 * sextic_q / 2 + higher_tail
    basis = coefficient_ideal(bordered_polynomial(initial_potential), parameters)
    for parameter in parameters:
        assert_radical_member(basis, parameter)

# In the 42 top layer every locally admissible quartic solves the constant
# Schur equation with the displayed q.  The positive-z coefficients of the
# full weighted face have radical (u,v), closing this last two-parameter
# family.
top_42 = x**4 * y**2
g_42_e4 = x**2 * y * (s6u * x + s6v * y)
q_42_e4 = (
    sp.Rational(5, 12) * s6u**2 * x**2
    + sp.Rational(4, 3) * s6u * s6v * x * y
    + sp.Rational(1, 6) * s6v**2 * y**2
)
assert sp.expand(binary_schur_face(top_42, g_42_e4, q_42_e4)) == 0
basis_42_e4 = coefficient_ideal(
    bordered_polynomial(
        top_42 + z * g_42_e4 + z**2 * q_42_e4 / 2 + s6r * z**3
    ),
    (s6u, s6v, s6r),
)
assert_radical_member(basis_42_e4, s6u)
assert_radical_member(basis_42_e4, s6v)
assert_radical_member(basis_42_e4, s6r)

# The resonant 321 layer does not even solve its scalar Schur equation.
face_321_e3 = binary_schur_face(
    x**3 * y**2 * (x - y), s6b * x * y * (x - y), s6Q
)
basis_321_e3 = sp.groebner(
    sp.Poly(face_321_e3, x, y).coeffs(), s6Q, s6b, order="lex"
)
assert_radical_member(basis_321_e3, s6b)


payload = {
    "format": "hc4-higher-degree-pencil-obstructions-v8",
    "status": [
        {
            "id": "HC4RSD17",
            "kind": "exact theorem",
            "scope": "arbitrary four-variable constant-Hessian pencils",
            "result": (
                "the relative Hessian endomorphism is polynomial, nilpotent, "
                "self-adjoint for the base Hessian metric, and Hessian-integrable"
            ),
        },
        {
            "id": "HC4RSD18",
            "kind": "hybrid theorem",
            "scope": "pencil directions of generic Hessian rank one in every degree",
            "result": "complete reduction to HC2 or the exact JC2 cotangent packet",
        },
        {
            "id": "HC4RSD19",
            "kind": "hybrid theorem",
            "scope": (
                "cubic directions whose leading homogeneous Hessian has rank three"
            ),
            "result": (
                "the moving-kernel normal form is impossible; the sole residual "
                "is the constant-kernel ternary Hessian-eikonal packet"
            ),
        },
        {
            "id": "HC4RSD20",
            "kind": "exact theorem",
            "scope": "residual constant-kernel leading-rank-three cubic packet",
            "result": (
                "the tangent Hessians of the cubic direction and border "
                "coefficient synchronize to one ruling; every fixed ruling "
                "reduces to HC2 or the exact JC2 cotangent packet"
            ),
        },
        {
            "id": "HC4RSD21",
            "kind": "hybrid theorem",
            "scope": (
                "homogeneous border coefficients in every degree and arbitrary "
                "border coefficients of degree at most three"
            ),
            "result": (
                "the universal-field equation forces a fixed cylinder, so the "
                "corresponding cubic eikonal packets reduce to HC2 or JC2"
            ),
        },
        {
            "id": "HC4RSD22",
            "kind": "hybrid theorem",
            "scope": "arbitrary border coefficients c of degree at most four",
            "result": (
                "the binary quartic Schur equation and all fourth-power-top "
                "correction charts force a fixed cylinder, hence HC2 or JC2"
            ),
        },
        {
            "id": "HC4RSD23",
            "kind": "hybrid theorem",
            "scope": (
                "degree-five border coefficients whose leading binary quintic "
                "is not a fifth power"
            ),
            "result": (
                "simple-root vanishing, the repeated-root Schur ideals, and "
                "two immutable next-face coefficients force a fixed cylinder"
            ),
        },
        {
            "id": "HC4RSD25",
            "kind": "hybrid theorem",
            "scope": (
                "arbitrary-degree border coefficients with squarefree leading "
                "binary form"
            ),
            "result": (
                "the simple-root square and degree deficit kill every "
                "transverse layer, forcing a fixed cylinder"
            ),
        },
        {
            "id": "HC4RSD26",
            "kind": "hybrid theorem",
            "scope": (
                "arbitrary-degree border coefficients whose leading binary "
                "form has exactly one double root and all other roots simple"
            ),
            "result": (
                "the double-root valuation gap adds the repeated root to the "
                "simple-root vanishing set, forcing every transverse layer to zero"
            ),
        },
        {
            "id": "HC4RSD27",
            "kind": "exact theorem",
            "scope": (
                "the first nonbinary weighted face over an arbitrary root of "
                "multiplicity m in a degree-d binary leading form"
            ),
            "result": (
                "an explicit resonance polynomial C(d,m,e,n) gives a root-order "
                "sieve for every transverse degree e and closes the sextic 3111 stratum"
            ),
        },
        {
            "id": "HC4RSD28",
            "kind": "hybrid theorem",
            "scope": (
                "degree-six border coefficients whose leading binary sextic is "
                "not a sixth power"
            ),
            "result": (
                "the valuation sieve and seven exact weighted Schur charts force "
                "a fixed cylinder; only the pure-sixth top remains"
            ),
        },
    ],
    "global_pencil": {
        "base": "S=Hess(psi), det(S)=delta in K^*",
        "direction": "T=Hess(A)",
        "relative_endomorphism": "N=S^{-1}*T=adj(S)*T/delta",
        "equivalence": "det(S+s*T)=delta iff det(I+s*N)=1 iff N is nilpotent",
        "extra_structure": ["N^T*S=S*N", "S and S*N are Hessians"],
    },
    "all_degree_rank_one": {
        "normal_form": "A=h(x)+affine",
        "pencil_face": "h''(x)*det Hess_(y,z,w)(psi)=0",
        "conclusion": (
            "HC4RSD15--HC4RSD16 apply without a degree bound; rank one is HC2/JC2"
        ),
    },
    "cubic_leading_rank_three": {
        "normal_forms": [
            "A=a(x,y,z)+affine(w) (constant kernel)",
            "A=w*z+y*b(z)+G(x,z), deg(b)<=2 (exceptional)",
        ],
        "moving_kernel": "v=(0,1,0,-b'(z))",
        "faces": {
            "s3": "v^T*S*v=0",
            "s2": "(X-G_xz*U)^2+G_xx*U*(2*Z-A_zz*U)=0",
            "s1_after_U_X_zero": "-G_xx*Z^2*psi_ww=0",
            "s0_after_U_X_zero": "Z^2*(psi_xw^2-psi_ww*psi_xx)=delta",
        },
        "moving_chart": {
            "coordinate": "r=w+b'(z)*y",
            "integrated_base": "psi=y*C(x,z,r)+D(x,z,r)",
            "s2_y_coefficient": "3*b''*G_xx*C_r^2",
            "unit_chain": [
                "C_r=0",
                "C_x=0",
                "D_rr=0",
                "D_xr and C'(z)+b''*D_r are units",
                "d_x(C'(z)+b''*D_r)=b''*D_xr != 0, contradiction",
            ],
            "result": "no genuinely moving leading-rank-three cubic direction survives",
        },
        "residual": (
            "A=a(x,y,z) and psi=w*b(x,y,z)+C(x,y,z), with "
            "grad(b)^T*adj(Hess(a))*grad(b)=0"
        ),
    },
    "constant_kernel_cubic": {
        "form": "A=a(x,y,z), psi=w*c(x,y,z)+D(x,y,z)",
        "tangent_pencil": "det(E+s*F+w*G)=unit on T=ker(dc)",
        "synchronization": (
            "F=Hess(a)|T and G=Hess(c)|T are proportional rank-one forms "
            "over the fraction field and have one common ruling"
        ),
        "fixed_ruling": (
            "c is a cylinder; affine c is triangular over HC2, while nonaffine "
            "c forces D=z*L+M and det Hess(psi)=Jac(c,L)^2"
        ),
        "universal_field_identity": (
            "det Hess_(tau,x,y,z)(tau*c)=-tau^2*grad(c)^T*adj(Hess(c))*grad(c)"
        ),
        "closed_border_classes": [
            "homogeneous c of every degree",
            "arbitrary c of degree at most four",
            "degree-five c with non-pure leading quintic",
            "arbitrary-degree c with squarefree leading binary form",
            "arbitrary-degree c with one double leading root and all others simple",
            "degree-six c with non-pure leading sextic",
        ],
        "quartic_schur_charts": {
            "1111": "g=0, q=0",
            "211": "g=0, q=0",
            "31": "g=0, q=0",
            "22": "g=b*x*y, q=b^2/4, then b=0 on the next face",
            "4": "passive cube/square/line charts all force a fixed cylinder",
        },
        "quintic_schur_charts": {
            "11111/2111": "simple-root count forces g=0 and q=0",
            "221": "radical(g,q)",
            "311": "radical(g,q)",
            "32": "g=b*x^2*y, q=11*b^2*x/30, then b=0",
            "41": "g=b*x^2*y, q=b^2*y/5, then b=0",
        },
        "root_valuation_sieve": {
            "determinant_order": "3*m-2",
            "numerator_order": "2*n+2*m-2",
            "resonance": (
                "C=d^2*m+d^2*n^2-2*d*e*m*n-2*d*e*m-d*m^2-d*m"
                "+e^2*m^2+2*e*m^2+m^2"
            ),
            "consequence": (
                "if n<m/2 and C is nonzero, polynomial Schur divisibility fails"
            ),
        },
        "sextic_schur_charts": {
            "3111": "root-weight sum exceeds e in every transverse degree",
            "2211/411": "unit first-face ideal",
            "321": "unit at e=4 and radical(g) at the resonant e=3 face",
            "51": "unique e=4 Schur solution and scalar z^3 tail die in the full weighted face",
            "222": "e=4 radical(g); unique e=3 solution dies in the weighted face",
            "42": "e=4 two-parameter solution with scalar z^3 and the e=3 scalar solution die in weighted faces",
            "33": "e=4 scalar solution with scalar z^3 and the e=3 two-parameter solution die in weighted faces",
        },
        "companion_closure": (
            "HC4BL5 closes the pure-fifth degree-five chart in the separate "
            "quintic bordered-lemma checker"
        ),
        "residual": (
            "genuinely moving nonhomogeneous ruling of degree at least seven; "
            "the pure-sixth branch closes in companion HC4RSD29--HC4RSD32 checkers"
        ),
    },
    "open_frontier": [
        (
            "constant-kernel cubic packet with a repeated-root leading border "
            "form of degree at least seven"
        ),
        "cubic directions of leading Hessian rank at most two",
        "directions of generic Hessian ranks two and three in higher degree",
        "matrix pivots and direct degree-five HC4",
    ],
    "calibrations": {
        "rank_one_pencil_coefficient": "h''*det(E3)",
        "moving_cubic_kernel": [0, 1, 0, "-z"],
        "moving_cubic_active_determinant": "-z^3",
        "moving_chart_s2_y_coefficient": "3*m*G_xx*C_r^2",
        "binary_null_pair_minor_power": 2,
        "fixed_cylinder_determinant": "Jac(c,L)^2",
        "cubic_border_root_charts": ["111", "21", "3"],
        "quartic_border_root_charts": ["1111", "211", "22", "31", "4"],
        "double_double_next_face": "-3*b^4*x^2*y^2*z^2/4",
        "quintic_border_root_charts": [
            "11111",
            "2111",
            "221",
            "311",
            "32",
            "41",
            "5",
        ],
        "triple_double_next_face": "-8*b^3*x^7*y^3*z/15",
        "quadruple_simple_next_face": "6*b^3*x^8*y^2*z/5",
        "sextic_root_partitions_closed": [
            "111111",
            "21111",
            "3111",
            "2211",
            "411",
            "321",
            "51",
            "222",
            "42",
            "33",
        ],
    },
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: constant-Hessian pencils are polynomial nilpotent deformations")
print("PASS: every all-degree rank-one direction reduces to HC2 or JC2")
print("PASS: classified the leading-rank-three cubic direction normal forms")
print("PASS: the moving cubic kernel forces the global null equation")
print("PASS: the last determinant faces give a unit contradiction")
print("THEOREM: no moving leading-rank-three cubic direction survives")
print("PASS: residual tangent Hessians synchronize to one ruling")
print("PASS: every fixed ruling reduces to HC2 or JC2")
print("THEOREM: homogeneous and degree-at-most-three border coefficients close")
print("PASS: classified every binary quartic transverse Schur face")
print("THEOREM: every degree-at-most-four border coefficient closes")
print("PASS: classified every non-pure binary quintic transverse Schur face")
print("THEOREM: every non-pure degree-five border coefficient closes")
print("THEOREM: every all-degree squarefree leading border coefficient closes")
print("THEOREM: every all-degree one-double-root border coefficient closes")
print("PASS: the root-valuation sieve isolates every repeated-root resonance")
print("THEOREM: every non-pure degree-six border coefficient closes")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
