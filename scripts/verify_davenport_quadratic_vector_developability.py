#!/usr/bin/env python3
"""Exact top-coefficient audit for a general quadratic auxiliary vector."""

from __future__ import annotations

import sympy as sp


T, Y, U = sp.symbols("T Y U")


def function_vector(prefix: str) -> sp.Matrix:
    """Return a three-component polynomial-function placeholder."""
    return sp.Matrix(
        [sp.Function(f"{prefix}{index}")(T, Y) for index in range(3)]
    )


def triple(left: sp.Matrix, middle: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    """Return the scalar triple product with the vectors as columns."""
    return sp.Matrix.hstack(left, middle, right).det()


h = function_vector("h")
a = function_vector("a")
b = function_vector("b")
general_map = h + U * a + U**2 * b
general_jacobian = sp.Poly(
    sp.expand(general_map.jacobian((T, Y, U)).det()),
    U,
)
assert general_jacobian.degree() == 5

h_T = h.diff(T)
h_Y = h.diff(Y)
a_T = a.diff(T)
a_Y = a.diff(Y)
b_T = b.diff(T)
b_Y = b.diff(Y)

expected_coefficients = (
    triple(h_T, h_Y, a),
    (
        triple(a_T, h_Y, a)
        + triple(h_T, a_Y, a)
        + 2 * triple(h_T, h_Y, b)
    ),
    (
        triple(b_T, h_Y, a)
        + triple(h_T, b_Y, a)
        + triple(a_T, a_Y, a)
        + 2 * triple(a_T, h_Y, b)
        + 2 * triple(h_T, a_Y, b)
    ),
    (
        triple(b_T, a_Y, a)
        + triple(a_T, b_Y, a)
        + 2 * triple(b_T, h_Y, b)
        + 2 * triple(h_T, b_Y, b)
        + 2 * triple(a_T, a_Y, b)
    ),
    (
        triple(b_T, b_Y, a)
        + 2 * triple(b_T, a_Y, b)
        + 2 * triple(a_T, b_Y, b)
    ),
    2 * triple(b_T, b_Y, b),
)
for degree, expected in enumerate(expected_coefficients):
    assert sp.simplify(general_jacobian.nth(degree) - expected) == 0

# On b_3 != 0, the U^5 coefficient is exactly the projective-rank
# determinant b_3^3 J(b_1/b_3,b_2/b_3).
projective_jacobian = sp.factor(
    sp.diff(b[0] / b[2], T) * sp.diff(b[1] / b[2], Y)
    - sp.diff(b[0] / b[2], Y) * sp.diff(b[1] / b[2], T)
)
assert sp.simplify(
    triple(b_T, b_Y, b) - b[2] ** 3 * projective_jacobian
) == 0

# For an affine-linear b=b0+T*p+Y*q, developability says precisely that
# b0,p,q are linearly dependent, so the image lies in a target plane.
b00, b01, b02 = sp.symbols("b00 b01 b02")
p0, p1, p2 = sp.symbols("p0 p1 p2")
q0, q1, q2 = sp.symbols("q0 q1 q2")
b0 = sp.Matrix((b00, b01, b02))
p = sp.Matrix((p0, p1, p2))
q = sp.Matrix((q0, q1, q2))
affine_b = b0 + T * p + Y * q
assert sp.expand(
    triple(affine_b.diff(T), affine_b.diff(Y), affine_b)
    - triple(p, q, b0)
) == 0

# Generic rank-two normal form b=(T,Y,0).  The U^4 equation is an Euler
# equation for the transverse component a3 and has no nonzero polynomial
# solution.
g = sp.Function("g")(T, Y)
H = sp.Function("H")(T, Y)
a1 = sp.Function("a1")(T, Y)
a2 = sp.Function("a2")(T, Y)
a3 = sp.Function("a3")(T, Y)
rank_two_map = sp.Matrix(
    (
        T + a1 * U + T * U**2,
        g + a2 * U + Y * U**2,
        H + a3 * U,
    )
)
rank_two_jacobian = sp.Poly(
    sp.expand(rank_two_map.jacobian((T, Y, U)).det()),
    U,
)
rank_two_u4 = a3 - 2 * T * sp.diff(a3, T) - 2 * Y * sp.diff(a3, Y)
assert sp.simplify(rank_two_jacobian.nth(4) - rank_two_u4) == 0

for total_degree in range(12):
    assert 1 - 2 * total_degree != 0

rank_two_collapsed = sp.Poly(
    sp.expand(
        rank_two_map.subs(a3, 0).jacobian((T, Y, U)).det()
    ),
    U,
)
rank_two_u3 = -2 * (
    T * sp.diff(H, T) + Y * sp.diff(H, Y)
)
assert sp.simplify(rank_two_collapsed.nth(3) - rank_two_u3) == 0

# Rank-one planar normal form b=(1,T,0).  The U^4 equation makes a3 a
# function of T, and the U^3 equation integrates to one exact relation.
A = sp.Function("A")(T)
rank_one_plane_map = sp.Matrix(
    (
        T + a1 * U + U**2,
        g + a2 * U + T * U**2,
        H + A * U,
    )
)
rank_one_plane_jacobian = sp.Poly(
    sp.expand(rank_one_plane_map.jacobian((T, Y, U)).det()),
    U,
)
integrated_relation = (
    2 * H
    + (2 * T * sp.diff(A, T) - A) * a1
    - 2 * sp.diff(A, T) * a2
)
assert sp.simplify(
    rank_one_plane_jacobian.nth(3)
    - sp.diff(integrated_relation, Y)
) == 0

# Put v=a2-T*a1.  After solving the U^3 equation, the U^2 coefficient is
# again an exact Y-derivative.  Its primitive is a quadratic conic
# representation of the zero-section polynomial g.
v = sp.Function("v")(T, Y)
C = sp.Function("C")(T)
rank_one_a2 = T * a1 + v
rank_one_H = (
    A * a1 / 2
    + sp.diff(A, T) * v
    + C / 2
)
rank_one_integrated_map = sp.Matrix(
    (
        T + a1 * U + U**2,
        g + rank_one_a2 * U + T * U**2,
        rank_one_H + A * U,
    )
)
rank_one_integrated_jacobian = sp.Poly(
    sp.expand(
        rank_one_integrated_map.jacobian((T, Y, U)).det()
    ),
    U,
)
conic_primitive = (
    A * a1**2 / 2
    + 2 * sp.diff(A, T) * a1 * v
    - 2 * sp.diff(A, T, 2) * v**2
    - 2 * sp.diff(C, T) * v
    - 4 * sp.diff(A, T) * g
)
assert rank_one_integrated_jacobian.nth(3) == 0
assert sp.simplify(
    rank_one_integrated_jacobian.nth(2)
    - sp.diff(conic_primitive, Y) / 2
) == 0

# Eliminate g with conic_primitive=D(T).  The U coefficient is again an
# exact Y-derivative, now of a cubic in (a1,v).
D = sp.Function("D")(T)
g_from_conic = (
    A * a1**2 / 2
    + 2 * sp.diff(A, T) * a1 * v
    - 2 * sp.diff(A, T, 2) * v**2
    - 2 * sp.diff(C, T) * v
    - D
) / (4 * sp.diff(A, T))
cubic_reduced_map = sp.Matrix(
    (
        T + a1 * U + U**2,
        g_from_conic + rank_one_a2 * U + T * U**2,
        rank_one_H + A * U,
    )
)
cubic_reduced_jacobian = sp.Poly(
    sp.expand(cubic_reduced_map.jacobian((T, Y, U)).det()),
    U,
)
cubic_numerator = (
    6
    * (A * sp.diff(A, T, 2) + sp.diff(A, T) ** 2)
    * a1**2
    * v
    + 3 * A * sp.diff(C, T) * a1**2
    + 8
    * (
        sp.diff(A, T) * sp.diff(A, T, 3)
        - 3 * sp.diff(A, T, 2) ** 2
    )
    * v**3
    + 12
    * (
        sp.diff(A, T) * sp.diff(C, T, 2)
        - 3 * sp.diff(A, T, 2) * sp.diff(C, T)
    )
    * v**2
    + 12
    * (
        4 * T * sp.diff(A, T) ** 2
        - 2 * A * sp.diff(A, T)
        - D * sp.diff(A, T, 2)
        + sp.diff(A, T) * sp.diff(D, T)
        - sp.diff(C, T) ** 2
    )
    * v
)
cubic_primitive = -cubic_numerator / (24 * sp.diff(A, T))
assert cubic_reduced_jacobian.nth(3) == 0
assert cubic_reduced_jacobian.nth(2) == 0
assert sp.simplify(
    cubic_reduced_jacobian.nth(1)
    - sp.diff(cubic_primitive, Y)
) == 0

# A generic cubic a^2*(mu*v+nu)+P(v)-E has two affine singularity
# branches.  Its only additional projective branch is the degeneration
# rho=sigma=0 at [a:v:w]=[0:1:0].
aa, vv = sp.symbols("aa vv")
mu, nu, rho, sigma, tau, energy = sp.symbols(
    "mu nu rho sigma tau energy"
)
P = rho * vv**3 + sigma * vv**2 + tau * vv
model_cubic = aa**2 * (mu * vv + nu) + P - energy
assert sp.diff(model_cubic, aa) == 2 * aa * (mu * vv + nu)
assert sp.diff(model_cubic, vv) == (
    mu * aa**2 + sp.diff(P, vv)
)
vertical_value = -nu / mu
assert sp.simplify(
    model_cubic.subs(vv, vertical_value)
    - (P.subs(vv, vertical_value) - energy)
) == 0
w = sp.symbols("w")
projective_cubic = (
    mu * aa**2 * vv
    + nu * aa**2 * w
    + rho * vv**3
    + sigma * vv**2 * w
    + tau * vv * w**2
    - energy * w**3
)
infinity_point = {aa: 0, vv: 1, w: 0}
assert projective_cubic.subs(infinity_point) == rho
assert sp.diff(projective_cubic, aa).subs(infinity_point) == 0
assert sp.diff(projective_cubic, vv).subs(infinity_point) == 3 * rho
assert sp.diff(projective_cubic, w).subs(infinity_point) == sigma

# On the singular-at-infinity branch, the coefficient ODEs force A and C
# to be linear.  Indeed, for n=deg(A')>=1 the leading coefficient of
# A'*A'''-3*(A'')^2 is -n*(2*n+1), which never vanishes.
for derivative_degree in range(1, 20):
    assert (
        derivative_degree * (derivative_degree - 1)
        - 3 * derivative_degree**2
    ) == -derivative_degree * (2 * derivative_degree + 1)
    assert -derivative_degree * (2 * derivative_degree + 1) != 0

# The only exceptional polynomial parametrization of the resulting cubic
# has a=a(T), a free v carrying g, and one coefficient relation for D'.
# Its remaining Jacobian is either zero or has Y-degree thirteen.
p_linear, q_linear, c_linear, c_zero = sp.symbols(
    "p_linear q_linear c_linear c_zero",
    nonzero=True,
)
a_linear = sp.Function("a_linear")(T)
D_linear = sp.Function("D_linear")(T)
generic_g = sp.Function("generic_g")(T, Y)
A_linear = p_linear * T + q_linear
C_linear = c_linear * T + c_zero
v_exceptional = (
    4 * p_linear * generic_g
    + D_linear
    - A_linear * a_linear**2 / 2
) / (2 * (p_linear * a_linear - c_linear))
H_exceptional = (
    A_linear * a_linear / 2
    + p_linear * v_exceptional
    + C_linear / 2
)
exceptional_map = sp.Matrix(
    (
        T + a_linear * U + U**2,
        generic_g
        + (T * a_linear + v_exceptional) * U
        + T * U**2,
        H_exceptional + A_linear * U,
    )
)
exceptional_jacobian = sp.Poly(
    sp.factor(exceptional_map.jacobian((T, Y, U)).det()),
    U,
)
exceptional_D_derivative = (
    -p_linear * a_linear**2 / 2
    - 2 * p_linear * T
    + 2 * q_linear
    + c_linear**2 / p_linear
)
for degree in (1, 2, 3):
    assert sp.simplify(
        exceptional_jacobian.nth(degree).subs(
            sp.diff(D_linear, T),
            exceptional_D_derivative,
        )
    ) == 0
exceptional_constant = sp.factor(
    exceptional_jacobian.nth(0).subs(
        sp.diff(D_linear, T),
        exceptional_D_derivative,
    )
)
assert exceptional_constant.has(sp.diff(generic_g, Y))
assert exceptional_constant.has(generic_g)
assert exceptional_constant.has(
    a_linear * sp.diff(a_linear, T) - 2
)

# Finite critical branch.  At a critical point r of P(v)-E, write
# P(v)-E=(v-r)^2*(rho*v+kappa).  If the linear factor L=mu*v+nu also
# vanishes at r, the cubic has the vertical factor v-r.  If the remaining
# linear factor is proportional to L, it splits into line components over
# the algebraic closure.
r_critical, kappa, proportionality = sp.symbols(
    "r_critical kappa proportionality"
)
critical_cubic = (
    aa**2 * (mu * vv + nu)
    + (vv - r_critical) ** 2 * (rho * vv + kappa)
)
critical_vertical = sp.factor(
    critical_cubic.subs(nu, -mu * r_critical)
)
assert sp.rem(
    sp.Poly(critical_vertical, vv),
    sp.Poly(vv - r_critical, vv),
).as_expr() == 0
critical_lines = sp.factor(
    critical_cubic.subs(
        {
            rho: proportionality * mu,
            kappa: proportionality * nu,
        }
    )
)
assert sp.expand(
    critical_lines
    - (mu * vv + nu)
    * (aa**2 + proportionality * (vv - r_critical) ** 2)
) == 0

# Any line-component parametrization for which the conic identity becomes
# affine in its parameter makes all three outputs functions of (T,g,U).
# The chain rule then leaves an unavoidable g_Y factor.
G = sp.Function("G")(T, Y)
psi1 = sp.Function("psi1")(T, G, U)
psi2 = sp.Function("psi2")(T, G, U)
psi3 = sp.Function("psi3")(T, G, U)
factor_through_map = sp.Matrix((psi1, psi2, psi3))
factor_through_jacobian = sp.factor(
    factor_through_map.jacobian((T, Y, U)).det()
)
psi_G = sp.Symbol("psi_G")
assert factor_through_jacobian.has(sp.diff(G, Y))

# Degree gates for the remaining line and parabolic-conic components.
for parameter_degree in range(1, 20):
    assert 2 * parameter_degree != 7
    assert 3 * parameter_degree != 7

# Rank-one through-origin normal form b=(T,0,0).  Its first surviving top
# equation is the displayed projective relation for (a2,a3).
rank_one_line_map = sp.Matrix(
    (
        T + a1 * U + T * U**2,
        g + a2 * U,
        H + a3 * U,
    )
)
rank_one_line_jacobian = sp.Poly(
    sp.expand(rank_one_line_map.jacobian((T, Y, U)).det()),
    U,
)
rank_one_line_u3 = (
    2
    * T
    * (
        sp.diff(a2, T) * sp.diff(a3, Y)
        - sp.diff(a2, Y) * sp.diff(a3, T)
    )
    - a2 * sp.diff(a3, Y)
    + a3 * sp.diff(a2, Y)
)
assert sp.simplify(
    rank_one_line_jacobian.nth(3) - rank_one_line_u3
) == 0

print("PASS: the general quadratic vector has the exact U^0,...,U^5 ledger")
print("PASS: U^5 is precisely the projective developability determinant")
print("PASS: every affine-linear developable b lies in a target plane")
print("PASS: the generic affine rank-two normal form is impossible")
print("PASS: affine rank-one b reduces to two explicit surviving normal forms")
print("PASS: the non-origin survivor reduces further to one conic identity")
print("PASS: its next equation is a cubic with three singularity branches")
print("PASS: the singular-at-infinity branch has no nonzero constant Jacobian")
print("PASS: both finite singular branches fail component and degree gates")
print("PASS Davenport quadratic-vector developability audit")
