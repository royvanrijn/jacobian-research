#!/usr/bin/env python3
"""Verify the algebra in the dense common-kernel sextic reduction.

The conceptual inputs are the characteristic-zero Gordan--Noether theorem
for homogeneous forms in at most four variables, HC_3, HC_2, and the
degree-at-most-100 plane Jacobian theorem.  This script checks the exact
determinant identities and the small binary-cubic elimination used between
those inputs.  The rank-three sextic theorem is the automatic common-kernel
corollary of the first determinant-layer identity.
"""

from __future__ import annotations

import sympy as sp


z, t = sp.symbols("z t")


# 1. If the sextic Hessian has a constant kernel direction, its rank-three
# block is a generic symmetric 3x3 matrix C.  The z^7 coefficient of
# det(H0 + z*A + z^2*B) is A_tt*det(C).
a00 = sp.symbols("a00")
c11, c12, c13, c22, c23, c33 = sp.symbols(
    "c11 c12 c13 c22 c23 c33"
)
C3 = sp.Matrix(
    [
        [c11, c12, c13],
        [c12, c22, c23],
        [c13, c23, c33],
    ]
)
h_entries = sp.symbols("h0:10")
a_entries = sp.symbols("a0:9")
H0 = sp.Matrix(
    [
        [h_entries[0], h_entries[1], h_entries[2], h_entries[3]],
        [h_entries[1], h_entries[4], h_entries[5], h_entries[6]],
        [h_entries[2], h_entries[5], h_entries[7], h_entries[8]],
        [h_entries[3], h_entries[6], h_entries[8], h_entries[9]],
    ]
)
A = sp.Matrix(
    [
        [a00, a_entries[0], a_entries[1], a_entries[2]],
        [a_entries[0], a_entries[3], a_entries[4], a_entries[5]],
        [a_entries[1], a_entries[4], a_entries[6], a_entries[7]],
        [a_entries[2], a_entries[5], a_entries[7], a_entries[8]],
    ]
)
B = sp.diag(0, 1, 1, 1)
B[1:4, 1:4] = C3
scaled_determinant = sp.Poly(
    (H0 + z * A + z**2 * B).det(method="berkowitz"), z
)
assert sp.expand(scaled_determinant.coeff_monomial(z**7) - a00 * C3.det()) == 0


# The two lower-rank leading coefficients identify the residual envelope.
# Rank two leaves the binary Hessian of h4 on the two-dimensional sextic
# kernel; rank one leaves the ternary Hessian on its three-dimensional
# kernel.
C2_leading = sp.Matrix([[c11, c12], [c12, c22]])
B_rank_two = sp.zeros(4)
B_rank_two[2:4, 2:4] = C2_leading
rank_two_scaled = sp.Poly(
    (H0 + z * A + z**2 * B_rank_two).det(method="berkowitz"),
    z,
)
P2_leading = A[0:2, 0:2]
assert sp.expand(
    rank_two_scaled.coeff_monomial(z**6)
    - C2_leading.det() * P2_leading.det()
) == 0

B_rank_one = sp.diag(0, 0, 0, c11)
rank_one_scaled = sp.Poly(
    (H0 + z * A + z**2 * B_rank_one).det(method="berkowitz"),
    z,
)
P3_leading = A[0:3, 0:3]
assert sp.expand(
    rank_one_scaled.coeff_monomial(z**5)
    - c11 * P3_leading.det()
) == 0

# The rank-two leading condition genuinely permits a variable null
# direction.  This near-miss has kernel (-y,x) in the (t,m)-plane and no
# nonzero constant kernel vector.
m_dual, x_base, y_base = sp.symbols("m_dual x_base y_base")
variable_kernel_quartic = (x_base * t + y_base * m_dual) ** 2
variable_kernel_hessian = sp.hessian(
    variable_kernel_quartic, (t, m_dual)
)
assert sp.expand(variable_kernel_hessian.det()) == 0
constant_first, constant_second = sp.symbols(
    "constant_first constant_second"
)
constant_kernel_product = (
    variable_kernel_hessian
    * sp.Matrix([constant_first, constant_second])
)
constant_kernel_equations = [
    coefficient
    for entry in constant_kernel_product
    for coefficient in sp.Poly(entry, x_base, y_base).coeffs()
]
constant_kernel_ideal = sp.groebner(
    constant_kernel_equations,
    constant_first,
    constant_second,
)
assert constant_kernel_ideal == sp.groebner(
    [constant_first, constant_second],
    constant_first,
    constant_second,
)


# 2. The isotropic case has psi=t*s(u)+phi(u).  Its t^2 coefficient is the
# bordered Hessian q^T adj(C) q.  Once a cubic kernel m is selected, the
# t coefficient factors as phi_mm times the binary bordered invariant R.
q1, q2, q3 = sp.symbols("q1 q2 q3")
q3_vector = sp.Matrix([q1, q2, q3])
p11, p12, p13, p22, p23, p33 = sp.symbols(
    "p11 p12 p13 p22 p23 p33"
)
P3 = sp.Matrix(
    [
        [p11, p12, p13],
        [p12, p22, p23],
        [p13, p23, p33],
    ]
)
bordered = sp.Matrix.vstack(
    sp.Matrix([[0, q1, q2, q3]]),
    sp.Matrix.hstack(q3_vector, P3 + t * C3),
)
bordered_polynomial = sp.Poly(bordered.det(method="berkowitz"), t)
assert sp.expand(
    bordered_polynomial.coeff_monomial(t**2)
    + (q3_vector.T * C3.adjugate() * q3_vector)[0]
) == 0

c11b, c12b, c22b = sp.symbols("c11b c12b c22b")
C2 = sp.Matrix([[c11b, c12b], [c12b, c22b]])
q2_vector = sp.Matrix([q1, q2])
C3_binary = sp.diag(1, 1, 0)
C3_binary[0:2, 0:2] = C2
q3_binary = sp.Matrix([q1, q2, 0])
bordered_binary = sp.Matrix.vstack(
    sp.Matrix([[0, q1, q2, 0]]),
    sp.Matrix.hstack(q3_binary, P3 + t * C3_binary),
)
binary_polynomial = sp.Poly(
    bordered_binary.det(method="berkowitz"), t
)
R_generic = (q2_vector.T * C2.adjugate() * q2_vector)[0]
assert sp.expand(
    binary_polynomial.coeff_monomial(t) + p33 * R_generic
) == 0


# 3. For a homogeneous ternary cubic a, Euler's identity converts the
# degree-six bordered term into (3/2)*a*det(Hess(a)).
x, y, m = sp.symbols("x y m")
ternary_coefficients = sp.symbols("u0:10")
ternary_monomials = (
    x**3,
    x**2 * y,
    x**2 * m,
    x * y**2,
    x * y * m,
    x * m**2,
    y**3,
    y**2 * m,
    y * m**2,
    m**3,
)
ternary_cubic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        ternary_coefficients, ternary_monomials, strict=True
    )
)
ternary_gradient = sp.Matrix(
    [sp.diff(ternary_cubic, variable) for variable in (x, y, m)]
)
ternary_hessian = sp.hessian(ternary_cubic, (x, y, m))
assert sp.expand(
    (ternary_gradient.T * ternary_hessian.adjugate() * ternary_gradient)[0]
    - sp.Rational(3, 2) * ternary_cubic * ternary_hessian.det()
) == 0


# 4. After Gordan--Noether supplies a kernel m orthogonal to the linear
# covector, write s=x+a(x,y).  The binary bordered invariant R vanishes only
# for a=alpha*x^3.  Three low coefficients give a triangular certificate.
alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta")
binary_cubic = (
    alpha * x**3
    + beta * x**2 * y
    + gamma * x * y**2
    + delta * y**3
)
binary_gradient = sp.Matrix(
    [1 + sp.diff(binary_cubic, x), sp.diff(binary_cubic, y)]
)
binary_hessian = sp.hessian(binary_cubic, (x, y))
binary_R = sp.Poly(
    sp.expand(
        (
            binary_gradient.T
            * binary_hessian.adjugate()
            * binary_gradient
        )[0]
    ),
    x,
    y,
)
assert binary_R.coeff_monomial(x) == 2 * gamma
assert binary_R.coeff_monomial(y) == 6 * delta
assert sp.expand(
    binary_R.coeff_monomial(x**3).subs({gamma: 0, delta: 0})
    + 4 * beta**2
) == 0
assert sp.expand(binary_R.as_expr().subs({beta: 0, gamma: 0, delta: 0})) == 0


# 5. If phi_mm=0, then phi=m*g(x,y)+h(x,y), and the four-variable
# determinant is the square of the plane Jacobian of (s,g), independently
# of h.  The displayed second derivatives are generic placeholders.
sx, sy, gx, gy = sp.symbols("sx sy gx gy")
uxx, uxy, uyy = sp.symbols("uxx uxy uyy")
cotangent_hessian = sp.Matrix(
    [
        [0, 0, sx, sy],
        [0, 0, gx, gy],
        [sx, gx, uxx, uxy],
        [sy, gy, uxy, uyy],
    ]
)
plane_jacobian = sx * gy - sy * gx
assert sp.expand(cotangent_hessian.det() - plane_jacobian**2) == 0


# 6. If s=x, expansion along the t,x hyperbolic block leaves the negative
# two-variable Hessian determinant.
fxx, fxy, fxm, fyy, fym, fmm = sp.symbols(
    "fxx fxy fxm fyy fym fmm"
)
linear_hessian = sp.Matrix(
    [
        [0, 1, 0, 0],
        [1, fxx, fxy, fxm],
        [0, fxy, fyy, fym],
        [0, fxm, fym, fmm],
    ]
)
assert sp.expand(
    linear_hessian.det() + (fyy * fmm - fym**2)
) == 0


print("PASS: rank-three leading sextic forces a quartic common-kernel direction")
print("PASS: rank-two/rank-one layers reduce to binary/ternary kernel Hessians")
print("PASS: (x*t+y*m)^2 certifies the residual variable-kernel phenomenon")
print("PASS: the isotropic t^2 and t bordered-Hessian coefficients factor exactly")
print("PASS: the ternary-cubic Euler/Gordan--Noether reduction is exact")
print("PASS: R=0 forces the binary cubic to be alpha*x^3")
print("PASS: the two terminal determinants are plane-Jacobian and HC_2 blocks")
print("SCOPE: dense common-kernel theorem and automatic rank-three corollary")
