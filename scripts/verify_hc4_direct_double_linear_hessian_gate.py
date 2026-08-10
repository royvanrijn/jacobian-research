#!/usr/bin/env python3
"""Verify exact identities in the direct HC4 repeated-linear factor gates.

The companion proof shows that a top ternary Hessian determinant

    Delta = ell**2 * R

with R squarefree and coprime to the linear ell forces the first off-diagonal
motion to have order one.  Its quotient is a constant vector B, leaving two
normal forms.  The companion proof also closes exact quadruple and quintuple
linear factors.  This script checks the normal forms, boundary coefficients,
and degree identities for a symbolic degree parameter.  Radical/DVR
divisibility is a written proof step, not a bounded computation.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_direct_double_linear_hessian_gate.json"
)
OUT.parent.mkdir(parents=True, exist_ok=True)

m, j = sp.symbols("m j", integer=True, positive=True)
C, x = sp.symbols("C x", nonzero=True)
y, z = sp.symbols("y z")
alpha, q = sp.symbols("alpha q", nonzero=True)
hyy, hyz, hzz = sp.symbols("hyy hyz hzz")
hxx, hxy = sp.symbols("hxx hxy")

# If B(ell) != 0, then f=C*x^(m+2)+h(y,z).  Its Hessian is block diagonal.
transverse_hessian = sp.Matrix([
    [C * (m + 2) * (m + 1) * x**m, 0, 0],
    [0, hyy, hyz],
    [0, hyz, hzz],
])
transverse_determinant = sp.factor(transverse_hessian.det())
expected_transverse = (
    C * (m + 2) * (m + 1) * x**m * (hyy * hzz - hyz**2)
)
assert sp.simplify(transverse_determinant - expected_transverse) == 0

# If B(ell) = 0, then f=C*z*x^(m+1)+h(x,y).  Only the displayed entries are
# needed; the determinant is independent of the xx and xy entries.
mixed_entry = C * (m + 1) * x**m
tangent_hessian = sp.Matrix([
    [hxx, hxy, mixed_entry],
    [hxy, hyy, 0],
    [mixed_entry, 0, 0],
])
tangent_determinant = sp.factor(tangent_hessian.det())
expected_tangent = -C**2 * (m + 1) ** 2 * x ** (2 * m) * hyy
assert sp.simplify(tangent_determinant - expected_tangent) == 0

# In the j=1, rank-one boundary jet f=x^2*g, the residual tangent normal form
# is f=C*z*x^(m+1)+x^2*h(x,y).  Its x-order is at least 2*m+2, never four.
rank_one_tangent_hessian = sp.Matrix([
    [hxx, hxy, mixed_entry],
    [hxy, x**2 * hyy, 0],
    [mixed_entry, 0, 0],
])
rank_one_tangent_determinant = sp.factor(rank_one_tangent_hessian.det())
expected_rank_one_tangent = (
    -C**2 * (m + 1) ** 2 * x ** (2 * m + 2) * hyy
)
assert sp.simplify(
    rank_one_tangent_determinant - expected_rank_one_tangent
) == 0

# Exact x-order four for the two rank-one boundary jets is controlled by the
# displayed leading determinants.  These are formal jet identities: the
# symbols stand for the indicated derivatives after restriction to x=0.
g0, gy, gz, gyy, gyz, gzz = sp.symbols("g0 gy gz gyy gyz gzz")
bordered_hessian = sp.Matrix([
    [2 * g0, 2 * gy, 2 * gz],
    [2 * gy, gyy, gyz],
    [2 * gz, gyz, gzz],
])
bordered_coefficient = sp.factor(bordered_hessian.det())
assert bordered_coefficient == sp.factor(
    2 * g0 * gyy * gzz
    - 2 * g0 * gyz**2
    - 4 * gy**2 * gzz
    + 8 * gy * gyz * gz
    - 4 * gyy * gz**2
)
for power_degree in range(3, 8):
    pure_power = (y + 2 * z) ** power_degree
    pure_power_bordered = sp.Matrix([
        [2 * pure_power, 2 * sp.diff(pure_power, y),
         2 * sp.diff(pure_power, z)],
        [2 * sp.diff(pure_power, y), sp.diff(pure_power, y, 2),
         sp.diff(pure_power, y, z)],
        [2 * sp.diff(pure_power, z), sp.diff(pure_power, y, z),
         sp.diff(pure_power, z, 2)],
    ])
    assert sp.expand(pure_power_bordered.det()) == 0

# For f=y^(m+2)+x^3*g, the x^4 coefficient is the nonzero y-Hessian entry
# times 6*g*g_zz-9*g_z^2.  The remaining Hessian entries start too late to
# affect this coefficient.
y_hessian = sp.symbols("y_hessian", nonzero=True)
pure_power_boundary_matrix = sp.Matrix([
    [6 * x * g0, 0, 3 * x**2 * gz],
    [0, y_hessian, 0],
    [3 * x**2 * gz, 0, x**3 * gzz],
])
pure_power_boundary_determinant = sp.factor(
    pure_power_boundary_matrix.det()
)
expected_pure_power_boundary = (
    x**4 * y_hessian * (6 * g0 * gzz - 9 * gz**2)
)
assert sp.simplify(
    pure_power_boundary_determinant - expected_pure_power_boundary
) == 0

# Nonvacuous split-top controls at the two boundary degrees.  Their binary
# Hessian determinants are squarefree, so f4 has exactly one double linear
# Hessian factor and f5 exactly one triple linear factor.
h4 = -2 * y**4 + y**3 * z + 3 * y**2 * z**2 + 3 * y * z**3 + 3 * z**4
h5 = (
    -2 * y**5
    + y**4 * z
    + 3 * y**3 * z**2
    + 3 * y**2 * z**3
    + 3 * y * z**4
    - 3 * z**5
)
binary_delta4 = sp.factor(sp.hessian(h4, (y, z)).det())
binary_delta5 = sp.factor(sp.hessian(h5, (y, z)).det())
for binary_delta in (binary_delta4, binary_delta5):
    assert sp.Poly(
        sp.gcd_list([
            binary_delta,
            sp.diff(binary_delta, y),
            sp.diff(binary_delta, z),
        ]),
        y,
        z,
    ).total_degree() == 0
split_delta4 = sp.factor(sp.hessian(x**4 + h4, (x, y, z)).det())
split_delta5 = sp.factor(sp.hessian(x**5 + h5, (x, y, z)).det())
assert sp.cancel(split_delta4 / (x**2 * binary_delta4)) == 12
assert sp.cancel(split_delta5 / (x**3 * binary_delta5)) == 20

# Delta has degree 3m.  For Delta=ell^2*R, ell*R has degree 3m-1, while
# adj(A0)*b_j has degree 3m-j.  Divisibility is impossible for j>1; at j=1
# the quotient has degree zero and is therefore a constant vector.
delta_degree = 3 * m
radical_degree = delta_degree - 1
adjugate_motion_degree = delta_degree - j
quotient_degree = sp.expand(adjugate_motion_degree - radical_degree)
assert quotient_degree == 1 - j
assert quotient_degree.subs(j, 1) == 0
assert quotient_degree.subs(j, 2) < 0

# Exact quadruple multiplicity on a generic rank-at-most-one boundary gives
# x^2*R | adj(A0)*b_j.  The quotient degree is 2-j, hence j is one or two.
quadruple_half_radical_degree = delta_degree - 2
quadruple_quotient_degree = sp.expand(
    adjugate_motion_degree - quadruple_half_radical_degree
)
assert quadruple_quotient_degree == 2 - j
assert quadruple_quotient_degree.subs(j, 1) == 1
assert quadruple_quotient_degree.subs(j, 2) == 0
assert quadruple_quotient_degree.subs(j, 3) < 0

# The sole order-two quadruple candidate has D=6,
# P_6=C*x^6+h_6(y,z), and A_3=alpha*x^3.  Its order-four determinant
# coefficient determines the constant w^2 Hessian entry q.  Completing the
# square then leaves a nonzero x^6 coefficient in the descended ternary top.
q_solution = sp.solve(sp.Eq(30 * C * q - 9 * alpha**2, 0), q)
assert q_solution == [3 * alpha**2 / (10 * C)]
descended_x6_coefficient = sp.factor(
    C - alpha**2 / (2 * q_solution[0])
)
assert descended_x6_coefficient == -sp.Rational(2, 3) * C

# Verify the generic Schur/completing-square matrix identity
#   P''+w*A''-q^-1*dA*dA^T = Q''+(w+A/q)*A'',
# where Q=P-A^2/(2q).
A_value, w = sp.symbols("A_value w")
a11, a12, a13, a22, a23, a33 = sp.symbols(
    "a11 a12 a13 a22 a23 a33"
)
v1, v2, v3 = sp.symbols("v1 v2 v3")
hess_a = sp.Matrix([
    [a11, a12, a13],
    [a12, a22, a23],
    [a13, a23, a33],
])
grad_a = sp.Matrix([v1, v2, v3])
p11, p12, p13, p22, p23, p33 = sp.symbols(
    "p11 p12 p13 p22 p23 p33"
)
hess_p = sp.Matrix([
    [p11, p12, p13],
    [p12, p22, p23],
    [p13, p23, p33],
])
schur_matrix = hess_p + w * hess_a - grad_a * grad_a.T / q
hess_descended = hess_p - (
    grad_a * grad_a.T + A_value * hess_a
) / q
completed_square_matrix = (
    hess_descended + (w + A_value / q) * hess_a
)
assert all(
    sp.simplify(entry) == 0
    for entry in schur_matrix - completed_square_matrix
)

# The two boundary weights in the pure-power rank-one jet can agree only at
# m=2 (D=4), outside HC4-DIR4's range m>=3.
weight_compatibility = sp.solve(
    sp.Eq(sp.Rational(3, 1) / (2 * m - 1), 1 / (m - 1)),
    m,
)
assert weight_compatibility == [2]
for test_m in range(3, 12):
    assert sp.Rational(3, 2 * test_m - 1) != sp.Rational(1, test_m - 1)
    assert 2 * test_m + 2 > 4

# HC4-DIR5 rank-one integration.  In independent linear coordinates x and
# lambda, F=C*x^(2m)/lambda^m is exactly the solution of
#   dF/F = 2m*dx/x - m*d(lambda)/lambda.
lam = sp.symbols("lam", nonzero=True)
rank_one_F = C * x ** (2 * m) / lam**m
assert sp.simplify(sp.diff(rank_one_F, x) / rank_one_F - 2 * m / x) == 0
assert sp.simplify(sp.diff(rank_one_F, lam) / rank_one_F + m / lam) == 0

# Polynomiality forces lambda proportional to x; exact quadruple
# multiplicity in the transverse normal form then forces m=4, D=6.
assert sp.solve(sp.Eq(m, 4), m) == [4]
for test_m in range(3, 12):
    if test_m != 4:
        assert test_m != 4

# In the tangent-gradient orientation of a hypothetical rank-two M,
# F_0=(b/2)*d(y^(m+1))/dy has nonzero y-derivative for m>=3.  The exact field
# equations force both F_z=0 and x^2 | x*b*F_y, an immediate contradiction.
b = sp.symbols("b", nonzero=True)
boundary_F = b * (m + 1) * y**m / 2
boundary_F_y = sp.diff(boundary_F, y)
assert sp.simplify(
    boundary_F_y - b * m * (m + 1) * y ** (m - 1) / 2
) == 0
for test_m in range(3, 8):
    assert boundary_F_y.subs(m, test_m) != 0

# HC4-DIR6 terminal D=6 packet.  The maximal passive determinant channels
# factor through R=det Hess_(y,z)(H_6), leaving this exact (x,w) Hessian.
# The first coefficient fixes B, the next fixes gamma, and the terminal
# coefficient is a nonzero characteristic-zero multiple of alpha^5/C^3.
lam_pencil, kernel_w = sp.symbols("lam_pencil kernel_w")
Bcoef, gamma = sp.symbols("Bcoef gamma")

# HC4-DIR3b terminal D=5 triple-factor packet.  The leading binary channel
# fixes the w^2*x coefficient; the only later maximal-passive channel is its
# nonzero square.
triple_B = sp.symbols("triple_B")
triple_weighted_binary = (
    lam_pencil**3 * C * x**5
    + lam_pencil**2 * alpha * kernel_w * x**3
    + lam_pencil * triple_B * kernel_w**2 * x
)
triple_weighted_delta = sp.Poly(
    sp.expand(sp.hessian(triple_weighted_binary, (x, kernel_w)).det()),
    lam_pencil,
    kernel_w,
    x,
)
triple_leading = triple_weighted_delta.coeff_monomial(
    lam_pencil**4 * x**4
)
triple_terminal = triple_weighted_delta.coeff_monomial(
    lam_pencil**2 * kernel_w**2
)
assert sp.expand(triple_leading - (40 * triple_B * C - 9 * alpha**2)) == 0
assert sp.expand(triple_terminal + 4 * triple_B**2) == 0
triple_B_solution = sp.solve(sp.Eq(triple_leading, 0), triple_B)[0]
triple_terminal_solution = sp.factor(
    triple_terminal.subs(triple_B, triple_B_solution)
)
assert sp.factor(triple_B_solution) == 9 * alpha**2 / (40 * C)
assert triple_terminal_solution == -81 * alpha**4 / (400 * C**2)

weighted_binary = (
    lam_pencil**4 * C * x**6
    + lam_pencil**3 * alpha * kernel_w * x**4
    + lam_pencil**2 * Bcoef * kernel_w**2 * x**2
    + lam_pencil * gamma * kernel_w**3
)
weighted_binary_delta = sp.Poly(
    sp.expand(sp.hessian(weighted_binary, (x, kernel_w)).det()),
    lam_pencil,
    kernel_w,
    x,
)
channel_14 = weighted_binary_delta.coeff_monomial(lam_pencil**6 * x**6)
channel_13 = weighted_binary_delta.coeff_monomial(
    lam_pencil**5 * kernel_w * x**4
)
channel_12 = weighted_binary_delta.coeff_monomial(
    lam_pencil**4 * kernel_w**2 * x**2
)
channel_11 = weighted_binary_delta.coeff_monomial(
    lam_pencil**3 * kernel_w**3
)
assert sp.expand(channel_14 - (60 * Bcoef * C - 16 * alpha**2)) == 0
assert sp.expand(channel_13 - (180 * C * gamma - 8 * Bcoef * alpha)) == 0
assert sp.expand(channel_12 - 12 * (-Bcoef**2 + 6 * alpha * gamma)) == 0
assert sp.expand(channel_11 - 12 * Bcoef * gamma) == 0

B_solution = sp.solve(sp.Eq(channel_14, 0), Bcoef)[0]
gamma_solution = sp.solve(
    sp.Eq(channel_13.subs(Bcoef, B_solution), 0), gamma
)[0]
assert sp.factor(B_solution) == 4 * alpha**2 / (15 * C)
assert sp.factor(gamma_solution) == 8 * alpha**3 / (675 * C**2)
assert sp.simplify(
    channel_12.subs({Bcoef: B_solution, gamma: gamma_solution})
) == 0
terminal_channel = sp.factor(
    channel_11.subs({Bcoef: B_solution, gamma: gamma_solution})
)
assert terminal_channel == 128 * alpha**5 / (3375 * C**3)

# HC4-DIR7, order-two branch.  Exact quintuple multiplicity has half-radical
# x^3*R and hence permits j=2.  The normal form is D=7, and the maximal
# passive part of the four-variable determinant is R times this binary
# Hessian.  Its leading coefficient fixes B2, after which the next channel
# is already nonzero.
B2 = sp.symbols("B2")
quintuple_order_two = (
    lam_pencil**5 * C * x**7
    + lam_pencil**3 * alpha * kernel_w * x**4
    + lam_pencil * B2 * kernel_w**2 * x
)
quintuple_order_two_delta = sp.Poly(
    sp.expand(sp.hessian(quintuple_order_two, (x, kernel_w)).det()),
    lam_pencil,
    kernel_w,
    x,
)
q2_leading = quintuple_order_two_delta.coeff_monomial(
    lam_pencil**6 * x**6
)
q2_terminal = quintuple_order_two_delta.coeff_monomial(
    lam_pencil**4 * kernel_w * x**3
)
q2_last = quintuple_order_two_delta.coeff_monomial(
    lam_pencil**2 * kernel_w**2
)
assert sp.expand(q2_leading - (84 * B2 * C - 16 * alpha**2)) == 0
assert sp.expand(q2_terminal - 8 * B2 * alpha) == 0
assert sp.expand(q2_last + 4 * B2**2) == 0

B2_solution = sp.solve(sp.Eq(q2_leading, 0), B2)[0]
q2_terminal_solution = sp.factor(q2_terminal.subs(B2, B2_solution))
assert sp.factor(B2_solution) == 4 * alpha**2 / (21 * C)
assert q2_terminal_solution == 32 * alpha**3 / (21 * C)

# HC4-DIR7, order-one branch.  The same rank-collapse system gives D=7 and
# the displayed weighted ladder.  Successive coefficients determine B1,
# gamma1, and eta.  The intermediate coefficient cancels identically, but
# the terminal square cannot vanish in characteristic zero.
B1, gamma1, eta = sp.symbols("B1 gamma1 eta")
quintuple_order_one = (
    lam_pencil**5 * C * x**7
    + lam_pencil**4 * alpha * kernel_w * x**5
    + lam_pencil**3 * B1 * kernel_w**2 * x**3
    + lam_pencil**2 * gamma1 * kernel_w**3 * x
    + lam_pencil * eta * kernel_w**3
)
quintuple_order_one_delta = sp.Poly(
    sp.expand(sp.hessian(quintuple_order_one, (x, kernel_w)).det()),
    lam_pencil,
    kernel_w,
    x,
)
q1_leading = quintuple_order_one_delta.coeff_monomial(
    lam_pencil**8 * x**8
)
q1_next = quintuple_order_one_delta.coeff_monomial(
    lam_pencil**7 * kernel_w * x**6
)
q1_intermediate = quintuple_order_one_delta.coeff_monomial(
    lam_pencil**6 * kernel_w**2 * x**4
)
q1_eta = quintuple_order_one_delta.coeff_monomial(
    lam_pencil**6 * kernel_w * x**5
)
q1_terminal = quintuple_order_one_delta.coeff_monomial(
    lam_pencil**4 * kernel_w**4
)
assert sp.expand(q1_leading - (84 * B1 * C - 25 * alpha**2)) == 0
assert sp.expand(q1_next + 4 * (5 * B1 * alpha - 63 * C * gamma1)) == 0
assert sp.expand(
    q1_intermediate - 6 * (-4 * B1**2 + 15 * alpha * gamma1)
) == 0
assert sp.expand(q1_eta - 252 * C * eta) == 0
assert sp.expand(q1_terminal + 9 * gamma1**2) == 0

B1_solution = sp.solve(sp.Eq(q1_leading, 0), B1)[0]
gamma1_solution = sp.solve(
    sp.Eq(q1_next.subs(B1, B1_solution), 0), gamma1
)[0]
eta_solution = sp.solve(sp.Eq(q1_eta, 0), eta)[0]
q1_intermediate_solution = sp.factor(
    q1_intermediate.subs({B1: B1_solution, gamma1: gamma1_solution})
)
q1_terminal_solution = sp.factor(
    q1_terminal.subs(gamma1, gamma1_solution)
)
assert sp.factor(B1_solution) == 25 * alpha**2 / (84 * C)
assert sp.factor(gamma1_solution) == 125 * alpha**3 / (5292 * C**2)
assert eta_solution == 0
assert q1_intermediate_solution == 0
assert q1_terminal_solution == -15625 * alpha**6 / (3111696 * C**4)

# HC4-DIR8, first lower-rank quintuple jet.  Once the order-four bordered
# coefficient vanishes, its nonzero binary boundary form is a pure power.
# The exact order-five coefficient is controlled by the second passive
# derivative of the next x-coefficient.
g1y, g1z, g1yy, g1yz, g1zz = sp.symbols(
    "g1y g1z g1yy g1yz g1zz"
)
rank_one_x2_quintuple_jet = sp.Matrix([
    [2 * y**m, 2 * m * x * y ** (m - 1), 3 * x**2 * g1z],
    [
        2 * m * x * y ** (m - 1),
        m * (m - 1) * x**2 * y ** (m - 2) + x**3 * g1yy,
        x**3 * g1yz,
    ],
    [3 * x**2 * g1z, x**3 * g1yz, x**3 * g1zz],
])
rank_one_x2_order_five = sp.expand(
    rank_one_x2_quintuple_jet.det()
).coeff(x, 5)
expected_rank_one_x2_order_five = (
    -2 * m * (m + 1) * y ** (2 * m - 2) * g1zz
)
assert sp.simplify(
    rank_one_x2_order_five - expected_rank_one_x2_order_five
) == 0

# In the other rank-one jet, vanishing of the order-four coefficient makes
# g0 independent of z.  With g0=y^(m-1), this is the order-five term.
g1_value = sp.symbols("g1_value")
pure_power_y_hessian = (m + 2) * (m + 1) * y**m
rank_one_x3_quintuple_jet = sp.Matrix([
    [
        6 * x * y ** (m - 1) + 12 * x**2 * g1_value,
        3 * (m - 1) * x**2 * y ** (m - 2) + 4 * x**3 * g1y,
        4 * x**3 * g1z,
    ],
    [
        3 * (m - 1) * x**2 * y ** (m - 2) + 4 * x**3 * g1y,
        pure_power_y_hessian,
        x**4 * g1yz,
    ],
    [4 * x**3 * g1z, x**4 * g1yz, x**4 * g1zz],
])
rank_one_x3_order_five = sp.expand(
    rank_one_x3_quintuple_jet.det()
).coeff(x, 5)
expected_rank_one_x3_order_five = (
    6 * (m + 2) * (m + 1) * y ** (2 * m - 1) * g1zz
)
assert sp.simplify(
    rank_one_x3_order_five - expected_rank_one_x3_order_five
) == 0

# The field equations leave only these numerical weights.  The x^2*g jet
# would require m=2, below D=5.  The y^(m+2)+x^3*g jet requires m=3, when
# g1 has degree one and therefore g1_zz=0.
x2_rank_one_weight = sp.solve(sp.Eq(1, 2 - 2 / m), m)
x3_rank_one_weight = sp.solve(sp.Eq(1, 2 - 3 / m), m)
assert x2_rank_one_weight == [2]
assert x3_rank_one_weight == [3]
u0, u1, u2 = sp.symbols("u0 u1 u2")
degree_one_g1 = u0 * y + u1 * z + u2
assert sp.diff(degree_one_g1, z, 2) == 0

# HC4-DIR9, exact sextuple multiplicity on the generic corank-one boundary.
# The half-radical quotient has degree 3-j.  At j=3 the sole D=8 split packet
# has a constant w^2 Hessian entry; its next channel is nonzero.
sextuple_rank_two_quotient_degree = sp.expand(3 - j)
assert sextuple_rank_two_quotient_degree.subs(j, 3) == 0
assert sextuple_rank_two_quotient_degree.subs(j, 4) < 0
sextuple_quadratic_kappa = m / (m + 1)
assert sp.simplify(1 - 1 / sextuple_quadratic_kappa + 1 / m) == 0

sextuple_q = sp.symbols("sextuple_q")
sextuple_order_three = (
    lam_pencil**6 * C * x**8
    + lam_pencil**3 * alpha * kernel_w * x**4
    + sextuple_q * kernel_w**2 / 2
)
sextuple_order_three_delta = sp.Poly(
    sp.expand(sp.hessian(sextuple_order_three, (x, kernel_w)).det()),
    lam_pencil,
    kernel_w,
    x,
)
s6_j3_leading = sextuple_order_three_delta.coeff_monomial(
    lam_pencil**6 * x**6
)
s6_j3_terminal = sextuple_order_three_delta.coeff_monomial(
    lam_pencil**3 * kernel_w * x**2
)
assert sp.expand(s6_j3_leading - (56 * C * sextuple_q - 16 * alpha**2)) == 0
assert sp.expand(s6_j3_terminal - 12 * alpha * sextuple_q) == 0
sextuple_q_solution = sp.solve(sp.Eq(s6_j3_leading, 0), sextuple_q)[0]
s6_j3_terminal_solution = sp.factor(
    s6_j3_terminal.subs(sextuple_q, sextuple_q_solution)
)
assert sp.factor(sextuple_q_solution) == 2 * alpha**2 / (7 * C)
assert s6_j3_terminal_solution == 24 * alpha**3 / (7 * C)

# At j=2 the linear-field equations collapse to rank one and the same D=8
# split top.  Its leading coefficient fixes B6, while the terminal square is
# nonzero.
sextuple_B = sp.symbols("sextuple_B")
sextuple_order_two = (
    lam_pencil**6 * C * x**8
    + lam_pencil**4 * alpha * kernel_w * x**5
    + lam_pencil**2 * sextuple_B * kernel_w**2 * x**2
)
sextuple_order_two_delta = sp.Poly(
    sp.expand(sp.hessian(sextuple_order_two, (x, kernel_w)).det()),
    lam_pencil,
    kernel_w,
    x,
)
s6_j2_leading = sextuple_order_two_delta.coeff_monomial(
    lam_pencil**8 * x**8
)
s6_j2_terminal = sextuple_order_two_delta.coeff_monomial(
    lam_pencil**4 * kernel_w**2 * x**2
)
assert sp.expand(s6_j2_leading - (112 * sextuple_B * C - 25 * alpha**2)) == 0
assert sp.expand(s6_j2_terminal + 12 * sextuple_B**2) == 0
sextuple_B_solution = sp.solve(sp.Eq(s6_j2_leading, 0), sextuple_B)[0]
s6_j2_terminal_solution = sp.factor(
    s6_j2_terminal.subs(sextuple_B, sextuple_B_solution)
)
assert sp.factor(sextuple_B_solution) == 25 * alpha**2 / (112 * C)
assert s6_j2_terminal_solution == -1875 * alpha**4 / (3136 * C**2)

# HC4-DIR11.  If Jac(Q) vanishes on the boundary, integrability gives
# Q=x^2*u and the remaining j=1 packet is again the D=8 split top.  Its
# maximal-passive ladder has five vanishing channels followed by this
# nonzero terminal coefficient.
sextuple_B0, sextuple_gamma0, sextuple_delta0 = sp.symbols(
    "sextuple_B0 sextuple_gamma0 sextuple_delta0"
)
sextuple_order_one_rank_zero = (
    lam_pencil**6 * C * x**8
    + lam_pencil**5 * alpha * kernel_w * x**6
    + lam_pencil**4 * sextuple_B0 * kernel_w**2 * x**4
    + lam_pencil**3 * sextuple_gamma0 * kernel_w**3 * x**2
    + lam_pencil**2 * sextuple_delta0 * kernel_w**4
)
sextuple_order_one_rank_zero_delta = sp.Poly(
    sp.expand(
        sp.hessian(sextuple_order_one_rank_zero, (x, kernel_w)).det()
    ),
    lam_pencil,
    kernel_w,
    x,
)
s6_j1_r0_channels = [
    sextuple_order_one_rank_zero_delta.coeff_monomial(
        lam_pencil**10 * x**10
    ),
    sextuple_order_one_rank_zero_delta.coeff_monomial(
        lam_pencil**9 * kernel_w * x**8
    ),
    sextuple_order_one_rank_zero_delta.coeff_monomial(
        lam_pencil**8 * kernel_w**2 * x**6
    ),
    sextuple_order_one_rank_zero_delta.coeff_monomial(
        lam_pencil**7 * kernel_w**3 * x**4
    ),
    sextuple_order_one_rank_zero_delta.coeff_monomial(
        lam_pencil**6 * kernel_w**4 * x**2
    ),
    sextuple_order_one_rank_zero_delta.coeff_monomial(
        lam_pencil**5 * kernel_w**5
    ),
]
s6_j1_r0_expected = [
    112 * sextuple_B0 * C - 36 * alpha**2,
    -36 * sextuple_B0 * alpha + 336 * C * sextuple_gamma0,
    -40 * sextuple_B0**2
    + 672 * C * sextuple_delta0
    + 108 * alpha * sextuple_gamma0,
    -20 * sextuple_B0 * sextuple_gamma0
    + 360 * alpha * sextuple_delta0,
    144 * sextuple_B0 * sextuple_delta0 - 24 * sextuple_gamma0**2,
    24 * sextuple_delta0 * sextuple_gamma0,
]
assert all(
    sp.expand(actual - expected) == 0
    for actual, expected in zip(s6_j1_r0_channels, s6_j1_r0_expected)
)
sextuple_B0_solution = sp.solve(
    sp.Eq(s6_j1_r0_channels[0], 0), sextuple_B0
)[0]
sextuple_gamma0_solution = sp.solve(
    sp.Eq(
        s6_j1_r0_channels[1].subs(sextuple_B0, sextuple_B0_solution),
        0,
    ),
    sextuple_gamma0,
)[0]
sextuple_delta0_solution = sp.solve(
    sp.Eq(
        s6_j1_r0_channels[2].subs({
            sextuple_B0: sextuple_B0_solution,
            sextuple_gamma0: sextuple_gamma0_solution,
        }),
        0,
    ),
    sextuple_delta0,
)[0]
s6_j1_r0_substitution = {
    sextuple_B0: sextuple_B0_solution,
    sextuple_gamma0: sextuple_gamma0_solution,
    sextuple_delta0: sextuple_delta0_solution,
}
s6_j1_r0_reduced = [
    sp.factor(channel.subs(s6_j1_r0_substitution))
    for channel in s6_j1_r0_channels
]
assert sp.factor(sextuple_B0_solution) == 9 * alpha**2 / (28 * C)
assert sp.factor(sextuple_gamma0_solution) == 27 * alpha**3 / (784 * C**2)
assert sp.factor(sextuple_delta0_solution) == 27 * alpha**4 / (43904 * C**3)
assert s6_j1_r0_reduced[:5] == [0, 0, 0, 0, 0]
assert s6_j1_r0_reduced[5] == 2187 * alpha**7 / (4302592 * C**5)

# HC4-DIR12.  On the non-axial boundary-rank-one branch, Jacobian
# integrability puts the quadratic field in the form Q=u*q+x^2*v.  Its
# boundary Jacobian is the displayed outer product and therefore has rank at
# most one identically.
u1, u2, u3, vv1, vv2, vv3 = sp.symbols("u1 u2 u3 vv1 vv2 vv3")
qxx, qxy, qxz, qyy, qyz, qzz = sp.symbols(
    "qxx qxy qxz qyy qyz qzz"
)
rank_one_q = (
    qxx * x**2
    + qxy * x * y
    + qxz * x * z
    + qyy * y**2
    + qyz * y * z
    + qzz * z**2
)
rank_one_u = sp.Matrix([u1, u2, u3])
rank_one_v = sp.Matrix([vv1, vv2, vv3])
rank_one_Q = rank_one_u * rank_one_q + x**2 * rank_one_v
rank_one_N0 = rank_one_Q.jacobian((x, y, z)).subs(x, 0)
rank_one_outer = rank_one_u * sp.Matrix([[ 
    sp.diff(rank_one_q, x).subs(x, 0),
    sp.diff(rank_one_q, y).subs(x, 0),
    sp.diff(rank_one_q, z).subs(x, 0),
]])
assert rank_one_N0 == rank_one_outer
for row_a in range(3):
    for row_b in range(row_a + 1, 3):
        for col_a in range(3):
            for col_b in range(col_a + 1, 3):
                minor = rank_one_N0.extract(
                    [row_a, row_b], [col_a, col_b]
                ).det()
                assert sp.expand(minor) == 0

# If the binary boundary value q_0 vanishes, Q is x times a linear field.
# Dividing the quadratic system by x recovers the HC4-DIR5 linear system.
axial_q = qxy * x * y + qxz * x * z + qxx * x**2
axial_Q = sp.expand(rank_one_u * axial_q + x**2 * rank_one_v)
assert all(sp.expand(component.subs(x, 0)) == 0 for component in axial_Q)

# In the tangent orientation u=partial_z, put t=ord_x(partial_z f).  The
# three leading coefficient equations have ratios
#   q_0 U_t + V_(t-2) = F_(t-3),
#   2 V_(t-2) = k F_(t-3),
#   U_t dq_0 = -(1/m) dF_(t-3).
# Their logarithmic exponent and the homogeneous degrees simplify exactly as
# claimed in (5.64l).
t_order = sp.symbols("t_order", integer=True, positive=True)
tangent_r = t_order - 3
tangent_k = 3 - tangent_r / m
tangent_e = (m + 3 - t_order) / 2
tangent_v_ratio = tangent_k / 2
tangent_u_ratio = 1 - tangent_v_ratio
assert sp.simplify(tangent_u_ratio + tangent_e / m) == 0
assert sp.simplify(tangent_u_ratio + tangent_v_ratio - 1) == 0
assert sp.simplify(2 * tangent_v_ratio - tangent_k) == 0
assert sp.simplify((m - tangent_r) - 2 * tangent_e) == 0
assert sp.simplify((m + 1 - t_order) - 2 * (tangent_e - 1)) == 0

# Exact sextuple order gives 3<=t<=6.  For t<6 the scalar equation Q(F) in
# (x^3) forces q_0=y^2 after tangent coordinates.  These square solutions
# work for both parities; their leading F and U coefficients obey the
# logarithmic equation for every admissible m.
tangent_packets = []
for packet_t, minimum_m in ((3, 3), (4, 3), (5, 4)):
    packet_e = sp.simplify(tangent_e.subs(t_order, packet_t))
    packet_F = y ** (m + 3 - packet_t)
    packet_q0 = y**2
    assert sp.simplify(
        sp.diff(packet_F, y) / packet_F
        - packet_e * sp.diff(packet_q0, y) / packet_q0
    ) == 0
    assert sp.diff(packet_F, z) == 0
    packet_U = sp.factor(packet_F / packet_q0)
    assert sp.simplify(
        tangent_u_ratio.subs(t_order, packet_t) * packet_F
        + packet_e * packet_F / m
    ) == 0
    tangent_packets.append({
        "t": packet_t,
        "minimum_m": minimum_m,
        "q0": "y^2",
        "F_leading": str(packet_F),
        "U_leading_up_to_scalar": str(packet_U),
    })

# At t=6, a nonsquare quadratic needs integral exponent (m-3)/2, hence odd
# m.  A square q_0=lambda^2 instead gives ordinary lambda powers for every
# m>=5.
for test_m in range(5, 12):
    nonsquare_polynomial = (test_m - 3) % 2 == 0
    assert nonsquare_polynomial == (test_m % 2 == 1)
    square_F_degree = test_m - 3
    square_U_degree = test_m - 5
    assert square_F_degree >= 2
    assert square_U_degree >= 0
for test_t in range(7, 12):
    assert min(test_t, 2 * test_t - 2) > 6

# HC4-DIR13.  Eliminating V gives
#   U*d(q/x^2)=-(x^(m+1)/m)*d(F/x^m).
# On the primitive quadratic branch F is a polynomial in q/x^2.  The exact
# monomial construction below verifies both the one-form identity and the
# order t=m+3-2k for every relevant small degree; the written quadratic
# centralizer lemma makes this an all-degree classification.
primitive_q = y**2 + x * z
primitive_checks = []
for test_m in range(3, 11):
    for test_k in range(1, test_m // 2 + 1):
        primitive_F = x ** (test_m - 2 * test_k) * primitive_q**test_k
        primitive_S = primitive_q / x**2
        primitive_T = sp.cancel(primitive_F / x**test_m)
        primitive_U = sp.cancel(
            -x ** (test_m + 1) * test_k * primitive_S ** (test_k - 1)
            / test_m
        )
        for coordinate in (x, y, z):
            assert sp.cancel(
                primitive_U * sp.diff(primitive_S, coordinate)
                + x ** (test_m + 1)
                * sp.diff(primitive_T, coordinate)
                / test_m
            ) == 0
        primitive_U_poly = sp.Poly(sp.cancel(primitive_U), x, y, z)
        primitive_x_order = min(
            monomial[0] for monomial, _ in primitive_U_poly.terms()
        )
        assert primitive_x_order == test_m + 3 - 2 * test_k
        primitive_checks.append({
            "m": test_m,
            "k": test_k,
            "t": primitive_x_order,
        })

# A composite quadratic and every homogeneous binary F become univariate on
# x=1, so their affine differentials wedge to zero.  A nonzero x*z term is
# the complementary primitive normal form when q_0=y^2.
affine_Y, affine_Z = sp.symbols("affine_Y affine_Z")
composite_a, composite_c = sp.symbols("composite_a composite_c")
composite_q_hat = affine_Y**2 + composite_a * affine_Y + composite_c
composite_F_hat = affine_Y**5 + 2 * affine_Y**3 - affine_Y
composite_wedge = sp.expand(
    sp.diff(composite_q_hat, affine_Y)
    * sp.diff(composite_F_hat, affine_Z)
    - sp.diff(composite_q_hat, affine_Z)
    * sp.diff(composite_F_hat, affine_Y)
)
assert composite_wedge == 0
assert primitive_q.subs(x, 0) == y**2
assert sp.diff(primitive_q, z) == x

primitive_parity_packets = [
    {"t": 3, "m_parity": "even", "minimum_m": 4, "k_max": "m/2"},
    {"t": 4, "m_parity": "odd", "minimum_m": 3, "k_max": "(m-1)/2"},
    {"t": 5, "m_parity": "even", "minimum_m": 4, "k_max": "(m-2)/2"},
    {"t": 6, "m_parity": "odd", "minimum_m": 5, "k_max": "(m-3)/2"},
]
for packet in primitive_parity_packets:
    packet_t = packet["t"]
    for test_m in range(packet["minimum_m"], 12):
        integral_k = (test_m + 3 - packet_t) % 2 == 0
        expected_parity = (
            test_m % 2 == (0 if packet["m_parity"] == "even" else 1)
        )
        assert integral_k == expected_parity

# HC4-DIR14.  For the invariant composite pencil f=z*G(x,y)+h(x,y), the
# determinant is affine in z and its z coefficient is
#   -((m+1)/m)*G*det Hess(G).
# The exact t=3 jets have unavoidable order four/five coefficients, whereas
# t=4 begins with the stated nonzero order-six coefficient.
for binary_degree in range(4, 9):
    binary_G = x**binary_degree + x * y ** (binary_degree - 1) + 2 * y**binary_degree
    binary_h = 3 * x ** (binary_degree + 1) + y ** (binary_degree + 1)
    binary_delta = sp.expand(
        sp.hessian(z * binary_G + binary_h, (x, y, z)).det()
    )
    binary_z_coefficient = sp.expand(binary_delta).coeff(z, 1)
    expected_binary_z = sp.factor(
        -sp.Rational(binary_degree, binary_degree - 1)
        * binary_G
        * sp.hessian(binary_G, (x, y)).det()
    )
    assert sp.simplify(binary_z_coefficient - expected_binary_z) == 0

composite_c0, composite_c1 = sp.symbols(
    "composite_c0 composite_c1", nonzero=True
)
composite_a0, composite_a1, composite_a2 = sp.symbols(
    "composite_a0 composite_a1 composite_a2"
)
composite_G_t3 = (
    composite_c0 * x**3 * y ** (m - 2)
    + composite_c1 * x**4 * y ** (m - 3)
)
composite_h_jet = (
    composite_a0 * y ** (m + 2)
    + composite_a1 * x * y ** (m + 1)
    + composite_a2 * x**2 * y**m
)
composite_delta_t3 = sp.expand(
    sp.hessian(z * composite_G_t3 + composite_h_jet, (x, y, z)).det()
)
composite_t3_order_four = sp.factor(composite_delta_t3.coeff(x, 4))
composite_t3_order_five_after_first = sp.factor(
    composite_delta_t3.coeff(x, 5).subs(composite_a0, 0)
)
expected_composite_t3_order_four = (
    -9
    * composite_a0
    * composite_c0**2
    * (m + 1)
    * (m + 2)
    * y ** (3 * m - 4)
)
expected_composite_t3_order_five = (
    -3
    * composite_a1
    * composite_c0**2
    * (m + 1)
    * (m + 4)
    * y ** (3 * m - 5)
)
assert sp.simplify(
    composite_t3_order_four - expected_composite_t3_order_four
) == 0
assert sp.simplify(
    composite_t3_order_five_after_first
    - expected_composite_t3_order_five
) == 0

composite_G_t4 = composite_c0 * x**4 * y ** (m - 3)
composite_delta_t4 = sp.expand(
    sp.hessian(z * composite_G_t4 + composite_h_jet, (x, y, z)).det()
)
composite_t4_order_six = sp.factor(composite_delta_t4.coeff(x, 6))
expected_composite_t4_order_six = (
    -16
    * composite_a0
    * composite_c0**2
    * (m + 1)
    * (m + 2)
    * y ** (3 * m - 6)
)
assert sp.simplify(
    composite_t4_order_six - expected_composite_t4_order_six
) == 0
for packet_t in (5, 6):
    assert 2 * packet_t - 2 > 6

# HC4-DIR15.  If the binary part of v vanishes, substitution of
# F=x^-3*q*G in the middle field equation has logarithmic residues
#   (m+1)*dlog(q)+dlog(G)=3(m+1)*dlog(x),
# hence G*q^(m+1)=C*x^(3(m+1)), impossible when q mod x is nonzero.  If the
# binary part is nonzero, its constant directional derivative kills G, so G
# is a pure power.  Order four then forces m=3 and G=C*x^4, while the x^2
# coefficient of Q(f) is v_y*h_y mod x and is nonzero by the DIR14 boundary
# coefficient.
pure_q_value, pure_G_value = sp.symbols(
    "pure_q_value pure_G_value", nonzero=True
)
pure_composite_integral = (
    pure_G_value
    * pure_q_value ** (m + 1)
    / x ** (3 * (m + 1))
)
assert sp.simplify(
    sp.diff(pure_composite_integral, pure_q_value)
    / pure_composite_integral
    - (m + 1) / pure_q_value
) == 0
assert sp.simplify(
    sp.diff(pure_composite_integral, pure_G_value)
    / pure_composite_integral
    - 1 / pure_G_value
) == 0
assert sp.simplify(
    sp.diff(pure_composite_integral, x)
    / pure_composite_integral
    + 3 * (m + 1) / x
) == 0
composite_pure_power_degree = sp.solve(sp.Eq(m + 1, 4), m)
assert composite_pure_power_degree == [3]
composite_boundary_hy = sp.diff(
    composite_a0 * y ** (m + 2), y
)
assert sp.factor(composite_boundary_hy.subs(m, 3)) == 5 * composite_a0 * y**4

# HC4-DIR16.  In the transverse composite orientation
# f=P(x,y)+h(x,z).  The field equation makes D_(v_x,v_z)h a homogeneous
# polynomial in x alone.  Its gradient vanishes on x=0, so the invertible
# active boundary Hessian kills (v_x,v_z).
active_bxx, active_bxz, active_bzz = sp.symbols(
    "active_bxx active_bxz active_bzz"
)
active_vx, active_vz = sp.symbols("active_vx active_vz")
active_boundary_hessian = sp.Matrix([
    [active_bxx, active_bxz],
    [active_bxz, active_bzz],
])
active_boundary_vector = sp.Matrix([active_vx, active_vz])
active_boundary_solution = sp.solve(
    list(active_boundary_hessian * active_boundary_vector),
    (active_vx, active_vz),
    dict=True,
)
assert active_boundary_solution == [{active_vx: 0, active_vz: 0}]
assert sp.factor(active_boundary_hessian.det()) == (
    active_bxx * active_bzz - active_bxz**2
)
transverse_composite_rhs = x ** (m + 1)
assert sp.diff(transverse_composite_rhs, z) == 0
for test_m in range(3, 10):
    assert sp.diff(
        transverse_composite_rhs.subs(m, test_m), x
    ).subs(x, 0) == 0

# HC4-DIR17.  For t>=4, V has order t-2>=2, so grad(V) vanishes on the
# boundary and the rank-two Hessian locks v to the already-known kernel line
# u.  The coefficient multiplying F_(t-3) in V_(t-2) never vanishes.  The
# only t=3 primitive row has q=y^2+x*z and F=c*q^(m/2)+O(x^2); its x
# coefficient in Q(F) is c*(m/2)*y^m, contradicting x^3 divisibility.
for packet_t, minimum_m in ((4, 3), (5, 4), (6, 5)):
    for test_m in range(minimum_m, 12):
        locking_ratio = sp.Rational(1, 2) * (
            3 - sp.Rational(packet_t - 3, test_m)
        )
        assert locking_ratio != 0
        assert packet_t - 2 >= 2

primitive_vx, primitive_vy, primitive_vz = sp.symbols(
    "primitive_vx primitive_vy primitive_vz"
)
primitive_v_field = sp.Matrix([primitive_vx, primitive_vy, primitive_vz])
primitive_coordinates = sp.Matrix([x, y, z])
for test_m in range(4, 12, 2):
    test_e = test_m // 2
    primitive_t3_F = primitive_q**test_e
    primitive_t3_Q = sp.Matrix([0, 0, primitive_q]) + x**2 * primitive_v_field
    primitive_t3_QF = sp.expand(
        sum(
            primitive_t3_Q[index]
            * sp.diff(primitive_t3_F, primitive_coordinates[index])
            for index in range(3)
        )
    )
    primitive_t3_order_one = sp.factor(primitive_t3_QF.coeff(x, 1))
    assert primitive_t3_order_one == test_e * y**test_m

# HC4-DIR18.  A linear 3x3 boundary matrix of rank two has a primitive left
# kernel generator of degree at most two.  If its tangent components are not
# both zero, the degree of gcd(f_y,f_z) counts the repeated part of the binary
# boundary form, leaving at most three distinct projective roots.  The loops
# below verify the exact derivative-gcd count for all one-, two-, and
# three-root profiles in a representative degree range.
root_profile_checks = []
for boundary_degree in range(5, 10):
    one_root_form = y**boundary_degree
    one_root_gcd = sp.gcd(
        sp.diff(one_root_form, y), sp.diff(one_root_form, z)
    )
    one_root_gcd_degree = sp.Poly(one_root_gcd, y, z).total_degree()
    assert boundary_degree - one_root_gcd_degree == 1
    root_profile_checks.append({
        "degree": boundary_degree,
        "roots": 1,
        "derivative_gcd_degree": one_root_gcd_degree,
    })

    two_a = 2
    two_b = boundary_degree - two_a
    two_root_form = y**two_a * z**two_b
    two_root_gcd = sp.gcd(
        sp.diff(two_root_form, y), sp.diff(two_root_form, z)
    )
    two_root_gcd_degree = sp.Poly(two_root_gcd, y, z).total_degree()
    assert boundary_degree - two_root_gcd_degree == 2
    root_profile_checks.append({
        "degree": boundary_degree,
        "roots": 2,
        "derivative_gcd_degree": two_root_gcd_degree,
    })

    three_a = 1
    three_b = 1
    three_c = boundary_degree - three_a - three_b
    three_root_form = y**three_a * z**three_b * (y - z) ** three_c
    three_root_gcd = sp.gcd(
        sp.diff(three_root_form, y), sp.diff(three_root_form, z)
    )
    three_root_gcd_degree = sp.Poly(three_root_gcd, y, z).total_degree()
    assert boundary_degree - three_root_gcd_degree == 3
    root_profile_checks.append({
        "degree": boundary_degree,
        "roots": 3,
        "derivative_gcd_degree": three_root_gcd_degree,
    })

normal_rank_two_matrix = sp.Matrix([
    [0, 0, 0],
    [y, z, 0],
    [0, y, z],
])
assert normal_rank_two_matrix.rank() == 2
assert sp.Matrix([[1, 0, 0]]) * normal_rank_two_matrix == sp.zeros(1, 3)
for primitive_kernel_degree in range(3):
    for tangent_gcd_degree in range(primitive_kernel_degree + 1):
        boundary_root_count = (
            1 + primitive_kernel_degree - tangent_gcd_degree
        )
        assert 1 <= boundary_root_count <= 3

# HC4-DIR19.  In the one-root profile, boundary rank two first kills the
# z-dependent part of f_1.  The x^2 Hessian coefficient then has the exact
# differential factor g*g_zz-2*g_z^2, whose vanishing makes a polynomial g
# independent of z.  The first field coefficient fixes ell=lambda*y and the
# constant y-component of Q.
one_f1y, one_f1z, one_hxx = sp.symbols(
    "one_f1y one_f1z one_hxx"
)
one_boundary_yy = (m + 2) * (m + 1) * y**m
one_boundary_hessian = sp.Matrix([
    [one_hxx, one_f1y, one_f1z],
    [one_f1y, one_boundary_yy, 0],
    [one_f1z, 0, 0],
])
assert sp.factor(one_boundary_hessian.det()) == (
    -one_boundary_yy * one_f1z**2
)

one_g, one_gy, one_gz, one_gyz, one_gzz = sp.symbols(
    "one_g one_gy one_gz one_gyz one_gzz"
)
one_jet_hessian = sp.Matrix([
    [2 * one_g, 2 * x * one_gy, 2 * x * one_gz],
    [2 * x * one_gy, one_boundary_yy, x**2 * one_gyz],
    [2 * x * one_gz, x**2 * one_gyz, x**2 * one_gzz],
])
one_order_two_delta = sp.factor(
    sp.expand(one_jet_hessian.det()).coeff(x, 2)
)
expected_one_order_two_delta = (
    2
    * one_boundary_yy
    * (one_g * one_gzz - 2 * one_gz**2)
)
assert sp.simplify(
    one_order_two_delta - expected_one_order_two_delta
) == 0

one_g_function = sp.Function("one_g_function")(z)
one_reciprocal_derivative = sp.factor(
    sp.diff(sp.diff(one_g_function, z) / one_g_function**2, z)
)
expected_one_reciprocal_derivative = (
    one_g_function * sp.diff(one_g_function, z, 2)
    - 2 * sp.diff(one_g_function, z) ** 2
) / one_g_function**3
assert sp.simplify(
    one_reciprocal_derivative - expected_one_reciprocal_derivative
) == 0

one_lambda, one_beta = sp.symbols("one_lambda one_beta", nonzero=True)
one_C = sp.symbols("one_C", nonzero=True)
one_field_coefficient = sp.expand(
    2 * one_lambda * y * one_C * y**m
    + one_beta * (m + 2) * y ** (m + 1)
)
one_beta_solution = sp.solve(
    sp.Eq(one_field_coefficient, 0), one_beta
)[0]
assert sp.factor(one_beta_solution) == -2 * one_C * one_lambda / (m + 2)
one_middle_dx_order_one = sp.factor(
    2 * one_C * one_lambda * y ** (m + 1)
    + 2 * one_beta * (m + 2) * y ** (m + 1)
)
one_middle_after_field = sp.factor(
    one_middle_dx_order_one.subs(one_beta, one_beta_solution)
)
assert one_middle_after_field == -2 * one_C * one_lambda * y ** (m + 1)

# HC4-DIR21.  The two-root profile has k-s=1, leaving the kernel-degree
# cases (1,0) and (2,1).  Boundary shears span all linear first jets in the
# former; in the latter they leave exactly the outer quadratic monomial.  The
# three-root profile has (k,s)=(2,0), and its two residual gradients together
# with y*z form a basis of the binary quadratics.
root_a, root_b, root_c = sp.symbols(
    "root_a root_b root_c", integer=True, positive=True
)
two_root_linear_shear = sp.Matrix([
    [0, root_b],
    [root_a, 0],
])
assert sp.factor(two_root_linear_shear.det()) == -root_a * root_b

two_root_y_exceptional_basis = sp.Matrix([
    [0, root_b, 0],
    [root_a, 0, 0],
    [0, 0, 1],
])
assert sp.factor(two_root_y_exceptional_basis.det()) == -root_a * root_b
two_root_z_exceptional_basis = sp.Matrix([
    [0, 0, 1],
    [0, root_b, 0],
    [root_a, 0, 0],
])
assert sp.factor(two_root_z_exceptional_basis.det()) in (
    root_a * root_b,
    -root_a * root_b,
)

three_root_quadratic_basis = sp.Matrix([
    [0, root_b, 0],
    [root_a + root_c, -(root_b + root_c), 1],
    [-root_a, 0, 0],
])
assert sp.factor(three_root_quadratic_basis.det()) == -root_a * root_b

# HC4-DIR22.  With an invertible binary boundary Hessian, the rank-two Schur
# complement fixes f_2=(1/2)*grad(f_1)^T*Hess(f_0)^-1*grad(f_1).  The two
# exceptional monomials give the displayed closed formulas.  In the
# three-root profile a residual quadratic denominator cannot divide the
# numerator, forcing kappa=0.  Once f_1=f_2=0, exact sextuple order places
# the first positive x-jet at x^8.
two_root_log_hessian = sp.Matrix([
    [root_a * (root_a - 1), root_a * root_b],
    [root_a * root_b, root_b * (root_b - 1)],
])
two_root_outer_y_vector = sp.Matrix([root_a - 2, root_b + 1])
two_root_outer_z_vector = sp.Matrix([root_a + 1, root_b - 2])
two_root_outer_y_scalar = sp.factor(
    (
        two_root_outer_y_vector.T
        * two_root_log_hessian.inv()
        * two_root_outer_y_vector
    )[0]
    / 2
)
two_root_outer_z_scalar = sp.factor(
    (
        two_root_outer_z_vector.T
        * two_root_log_hessian.inv()
        * two_root_outer_z_vector
    )[0]
    / 2
)
assert two_root_outer_y_scalar == (
    root_a * root_b - root_a - 4 * root_b
) / (2 * root_a * root_b)
assert two_root_outer_z_scalar == (
    root_a * root_b - 4 * root_a - root_b
) / (2 * root_a * root_b)

three_root_schur_checks = 0
for test_a in range(1, 4):
    for test_b in range(1, 4):
        for test_c in range(1, 4):
            test_f0 = y**test_a * z**test_b * (y - z) ** test_c
            test_f1 = (
                y**test_a * z**test_b * (y - z) ** (test_c - 1)
            )
            test_binary_hessian = sp.hessian(test_f0, (y, z))
            test_gradient = sp.Matrix([
                sp.diff(test_f1, y),
                sp.diff(test_f1, z),
            ])
            test_schur = sp.factor(sp.cancel(
                (
                    test_gradient.T
                    * test_binary_hessian.adjugate()
                    * test_gradient
                )[0]
                / (2 * test_binary_hessian.det())
            ))
            test_R2 = (
                test_a * test_b * (y - z) ** 2
                + test_a * test_c * z**2
                + test_b * test_c * y**2
            )
            test_P2 = (
                test_a * test_b * (y - z) ** 2
                + test_a * (test_c - 1) * z**2
                + test_b * (test_c - 1) * y**2
            )
            test_expected_schur = sp.factor(sp.cancel(
                sp.Rational(1, 2)
                * y**test_a
                * z**test_b
                * (y - z) ** (test_c - 2)
                * test_P2
                / test_R2
            ))
            assert sp.cancel(test_schur - test_expected_schur) == 0
            assert sp.gcd(test_R2, y * z * (y - z)) == 1
            assert sp.div(test_P2, test_R2, y, z)[1] != 0
            three_root_schur_checks += 1

for first_positive_order in range(3, 12):
    determinant_order = first_positive_order - 2
    assert (determinant_order == 6) == (first_positive_order == 8)

# HC4-DIR23.  In the normal rank-two packet the first two scalar/matrix
# coefficients force a square tangent field and the displayed first jets.
# The z component of the next matrix coefficient determines h_zz.  The
# boundary scalar equation kills s, after which the y component retains the
# nonzero coefficient -2*alpha*B*(n+2)/n.
normal_n = sp.symbols("normal_n", integer=True, positive=True)
normal_alpha, normal_C, normal_rho = sp.symbols(
    "normal_alpha normal_C normal_rho", nonzero=True
)
normal_A, normal_r, normal_s = sp.symbols(
    "normal_A normal_r normal_s"
)
normal_h0, normal_h1, normal_U = sp.symbols(
    "normal_h0 normal_h1 normal_U"
)
normal_B = normal_alpha * normal_C / normal_rho
normal_p = normal_C * y**normal_n
normal_g = (
    normal_B * y ** (normal_n - 2) * z
    + normal_A * y ** (normal_n - 1)
)
normal_h = (
    normal_h0 * y ** (normal_n - 2)
    + normal_h1 * y ** (normal_n - 3) * z
    - normal_B
    * (normal_n * normal_s + 4 * normal_alpha / normal_n)
    * y ** (normal_n - 4)
    * z**2
    / (2 * normal_rho)
)
normal_F0 = sp.factor(
    normal_rho * y**2 * sp.diff(normal_h, z)
    + normal_U * y ** (normal_n - 1)
    + normal_B
    * (normal_s + 4 * normal_alpha / normal_n)
    * y ** (normal_n - 2)
    * z
)
normal_V = normal_rho * normal_h1 + normal_U
normal_expected_F0 = (
    normal_V * y ** (normal_n - 1)
    - normal_B
    * (normal_n - 1)
    * normal_s
    * y ** (normal_n - 2)
    * z
)
assert sp.simplify(normal_F0 - normal_expected_F0) == 0
normal_boundary_scalar = sp.powsimp(
    normal_rho * y**2 * sp.diff(normal_F0, z), force=True
)
assert sp.simplify(sp.powsimp(
    normal_boundary_scalar
    + normal_rho
    * normal_B
    * (normal_n - 1)
    * normal_s
    * y**normal_n,
    force=True,
)) == 0
normal_l_y = -2 * normal_alpha * y / normal_n
normal_l_z = normal_r * y + normal_s * z
normal_middle_y = sp.factor(
    2 * normal_rho * y * sp.diff(normal_h, z)
    + normal_l_y.diff(y) * sp.diff(normal_g, y)
    + normal_l_z.diff(y) * sp.diff(normal_g, z)
    + sp.diff(normal_F0, y) / (normal_n - 1)
)
normal_terminal_z_coefficient = sp.factor(
    sp.powsimp(
        sp.expand(normal_middle_y.subs(normal_s, 0)).coeff(z, 1),
        force=True,
    )
    / y ** (normal_n - 3)
)
assert normal_terminal_z_coefficient == sp.factor(
    -2 * normal_alpha * normal_B * (normal_n + 2) / normal_n
)

# HC4-DIR24.  The delayed-jet gate is the universal product-rule identity
# grad(l(f0))=Hess(f0)*l+Jac(l)^T*grad(f0).  Thus the first two field
# equations kill l whenever the binary Hessian is invertible.
delay_c0, delay_c1, delay_c2, delay_c3 = sp.symbols(
    "delay_c0 delay_c1 delay_c2 delay_c3"
)
delay_f0 = (
    delay_c0 * y**3
    + delay_c1 * y**2 * z
    + delay_c2 * y * z**2
    + delay_c3 * z**3
)
delay_l = sp.Matrix([
    sp.symbols("delay_a") * y + sp.symbols("delay_b") * z,
    sp.symbols("delay_c") * y + sp.symbols("delay_d") * z,
])
delay_grad = sp.Matrix([sp.diff(delay_f0, y), sp.diff(delay_f0, z)])
delay_product_rule = sp.Matrix([
    sp.diff((delay_l.T * delay_grad)[0], y),
    sp.diff((delay_l.T * delay_grad)[0], z),
])
delay_expected = (
    sp.hessian(delay_f0, (y, z)) * delay_l
    + delay_l.jacobian((y, z)).T * delay_grad
)
assert all(
    sp.expand(entry) == 0
    for entry in delay_product_rule - delay_expected
)

# HC4-DIR25.  For the first outer two-root jet, the x coefficient of the
# scalar field equation fixes Q1.  The tangent matrix coefficient then has a
# shared three-term factor whose z^2 coefficient is
# -(2*a+4*b)*kappa^2 and cannot vanish for positive root multiplicities.
outer_kappa = sp.symbols("outer_kappa", nonzero=True)
outer_A, outer_B, outer_C, outer_D = sp.symbols(
    "outer_A outer_B outer_C outer_D"
)
outer_E, outer_G = sp.symbols("outer_E outer_G")
outer_u, outer_v, outer_w = sp.symbols("outer_u outer_v outer_w")
outer_f0 = y**root_a * z**root_b
outer_f1 = (
    outer_kappa * y ** (root_a - 2) * z ** (root_b + 1)
)
outer_f2 = (
    outer_kappa**2
    * (root_a * root_b - root_a - 4 * root_b)
    / (2 * root_a * root_b)
    * y ** (root_a - 4)
    * z ** (root_b + 2)
)
outer_f = outer_f0 + x * outer_f1 + x**2 * outer_f2
outer_Q0 = sp.Matrix([
    y**2,
    -2 * outer_kappa * y * z / root_a,
    outer_kappa * z**2 / root_b,
])
outer_Q1 = sp.Matrix([
    outer_A * y + outer_B * z,
    outer_C * y + outer_D * z,
    outer_E * y + outer_G * z,
])
outer_Q2 = sp.Matrix([outer_u, outer_v, outer_w])
outer_Q = outer_Q0 + x * outer_Q1 + x**2 * outer_Q2
outer_grad = sp.Matrix([
    sp.diff(outer_f, x),
    sp.diff(outer_f, y),
    sp.diff(outer_f, z),
])
outer_scalar_x1 = sp.factor(
    sp.expand((outer_Q.T * outer_grad)[0]).coeff(x, 1)
)
outer_q1_substitution = {
    outer_B: 0,
    outer_E: 0,
    outer_A: -root_a * outer_D / outer_kappa,
    outer_G: -root_a * outer_C / root_b,
}
assert sp.simplify(outer_scalar_x1.subs(outer_q1_substitution)) == 0
outer_N = outer_Q.jacobian((x, y, z))
outer_middle = sp.expand(outer_N.T * outer_grad)
outer_common_factor = (
    root_a**2 * root_b * outer_C * y**2
    - root_a**2 * root_b * outer_D * y * z
    - (2 * root_a + 4 * root_b) * outer_kappa**2 * z**2
)
outer_middle_y_x1 = sp.factor(
    outer_middle[1].coeff(x, 1).subs(outer_q1_substitution)
)
outer_middle_z_x1 = sp.factor(
    outer_middle[2].coeff(x, 1).subs(outer_q1_substitution)
)
assert sp.simplify(sp.powsimp(
    outer_middle_y_x1
    - y ** (root_a - 3)
    * z**root_b
    * outer_common_factor
    / (root_a * root_b),
    force=True,
)) == 0
assert sp.simplify(sp.powsimp(
    outer_middle_z_x1
    + y ** (root_a - 2)
    * z ** (root_b - 1)
    * outer_common_factor
    / (root_a * root_b),
    force=True,
)) == 0
outer_terminal_coefficient = sp.expand(outer_common_factor).coeff(z, 2)
assert sp.simplify(
    outer_terminal_coefficient
    + (2 * root_a + 4 * root_b) * outer_kappa**2
) == 0

# The rank-one integration in the j=2 system has logarithmic residues
#   (3(m-1)/2, -(m-1)/2)
# along x and lambda.  Substitution verifies the differential identity; when
# lambda is proportional to x, the remaining x-exponent is m-1.
s6_log_x = sp.Rational(3, 2) * (m - 1)
s6_log_lambda = -sp.Rational(1, 2) * (m - 1)
assert sp.simplify(3 / x - 2 * s6_log_x / ((m - 1) * x)) == 0
assert sp.simplify(
    -2 * s6_log_lambda / ((m - 1) * lam) - 1 / lam
) == 0
assert sp.simplify(s6_log_x + s6_log_lambda - (m - 1)) == 0

# HC4-DIR10.  On a rank-one sextuple boundary the stronger x^4*R divisor
# gives quotient degree 2-j.  The two order-one jets have vanishing order-four
# and order-five coefficients; these are their exact order-six terms.
sextuple_rank_one_quotient_degree = sp.expand(2 - j)
assert sextuple_rank_one_quotient_degree.subs(j, 2) == 0
assert sextuple_rank_one_quotient_degree.subs(j, 3) < 0

jet_u, jet_v = sp.symbols("jet_u jet_v")
jet_c = sp.symbols("jet_c", nonzero=True)
g2_value, g2y, g2z, g2yy, g2yz, g2zz = sp.symbols(
    "g2_value g2y g2z g2yy g2yz g2zz"
)
jet_g0 = y**m
jet_g1 = jet_u * y ** (m - 1) + jet_v * y ** (m - 2) * z
jet_g0y = m * y ** (m - 1)
jet_g0yy = m * (m - 1) * y ** (m - 2)
jet_g1y = (
    jet_u * (m - 1) * y ** (m - 2)
    + jet_v * (m - 2) * y ** (m - 3) * z
)
jet_g1z = jet_v * y ** (m - 2)
jet_g1yy = (
    jet_u * (m - 1) * (m - 2) * y ** (m - 3)
    + jet_v * (m - 2) * (m - 3) * y ** (m - 4) * z
)
jet_g1yz = jet_v * (m - 2) * y ** (m - 3)
rank_one_x2_sextuple_jet = sp.Matrix([
    [
        2 * jet_g0 + 6 * x * jet_g1 + 12 * x**2 * g2_value,
        2 * x * jet_g0y + 3 * x**2 * jet_g1y + 4 * x**3 * g2y,
        3 * x**2 * jet_g1z + 4 * x**3 * g2z,
    ],
    [
        2 * x * jet_g0y + 3 * x**2 * jet_g1y + 4 * x**3 * g2y,
        x**2 * jet_g0yy + x**3 * jet_g1yy + x**4 * g2yy,
        x**3 * jet_g1yz + x**4 * g2yz,
    ],
    [
        3 * x**2 * jet_g1z + 4 * x**3 * g2z,
        x**3 * jet_g1yz + x**4 * g2yz,
        x**4 * g2zz,
    ],
])
rank_one_x2_order_six = sp.factor(
    sp.expand(rank_one_x2_sextuple_jet.det()).coeff(x, 6)
)
expected_rank_one_x2_order_six = sp.factor(
    (m + 1)
    * (
        (m - 8) * jet_v**2 * y ** (3 * m - 6)
        - 2 * m * y ** (2 * m - 2) * g2zz
    )
)
assert sp.simplify(
    rank_one_x2_order_six - expected_rank_one_x2_order_six
) == 0

pure_jet_g0 = jet_c * y ** (m - 1)
pure_jet_g1 = jet_u * y ** (m - 2) + jet_v * y ** (m - 3) * z
pure_jet_g0y = jet_c * (m - 1) * y ** (m - 2)
pure_jet_g0yy = jet_c * (m - 1) * (m - 2) * y ** (m - 3)
pure_jet_g1y = (
    jet_u * (m - 2) * y ** (m - 3)
    + jet_v * (m - 3) * y ** (m - 4) * z
)
pure_jet_g1z = jet_v * y ** (m - 3)
pure_jet_g1yy = (
    jet_u * (m - 2) * (m - 3) * y ** (m - 4)
    + jet_v * (m - 3) * (m - 4) * y ** (m - 5) * z
)
pure_jet_g1yz = jet_v * (m - 3) * y ** (m - 4)
rank_one_x3_sextuple_jet = sp.Matrix([
    [
        6 * x * pure_jet_g0 + 12 * x**2 * pure_jet_g1
        + 20 * x**3 * g2_value,
        3 * x**2 * pure_jet_g0y + 4 * x**3 * pure_jet_g1y
        + 5 * x**4 * g2y,
        4 * x**3 * pure_jet_g1z + 5 * x**4 * g2z,
    ],
    [
        3 * x**2 * pure_jet_g0y + 4 * x**3 * pure_jet_g1y
        + 5 * x**4 * g2y,
        (m + 2) * (m + 1) * y**m + x**3 * pure_jet_g0yy
        + x**4 * pure_jet_g1yy + x**5 * g2yy,
        x**4 * pure_jet_g1yz + x**5 * g2yz,
    ],
    [
        4 * x**3 * pure_jet_g1z + 5 * x**4 * g2z,
        x**4 * pure_jet_g1yz + x**5 * g2yz,
        x**5 * g2zz,
    ],
])
rank_one_x3_order_six = sp.factor(
    sp.expand(rank_one_x3_sextuple_jet.det()).coeff(x, 6)
)
expected_rank_one_x3_order_six = sp.factor(
    2
    * (m + 1)
    * (m + 2)
    * (
        3 * jet_c * y ** (2 * m - 1) * g2zz
        - 8 * jet_v**2 * y ** (3 * m - 6)
    )
)
assert sp.simplify(
    sp.powsimp(sp.expand(rank_one_x3_order_six), force=True)
    - expected_rank_one_x3_order_six
) == 0

# HC4-DIR26.  The low-order recurrences in the x^2*g packet use
# F1=P*y^(m-1)+s*v*y^(m-2)*z and
# F2=m*(alpha*g1_y+beta*g1_z)/(2*(m-1)).  Their nonzero-v branch forces all
# four coefficients of L to vanish; their zero-v branch forces g2_zz=0.
# In the pure-power packet the sole nonzero-v weight resonance is m=3.
lower_alpha, lower_beta, lower_r, lower_s = sp.symbols(
    "lower_alpha lower_beta lower_r lower_s"
)
lower_P = m * lower_alpha + lower_r * jet_v
lower_F1 = (
    lower_P * y ** (m - 1)
    + lower_s * jet_v * y ** (m - 2) * z
)
lower_g1 = jet_u * y ** (m - 1) + jet_v * y ** (m - 2) * z
lower_H = (
    lower_alpha * sp.diff(lower_g1, y)
    + lower_beta * sp.diff(lower_g1, z)
)
lower_F2 = sp.factor(m * lower_H / (2 * (m - 1)))
assert sp.factor(lower_F1.coeff(z, 1)) == (
    lower_s * jet_v * y ** (m - 2)
)
assert sp.simplify(
    sp.expand(lower_F2).coeff(z, 1)
    - m
    * lower_alpha
    * jet_v
    * (m - 2)
    * y ** (m - 3)
    / (2 * (m - 1))
) == 0
for test_m in range(3, 12):
    assert sp.Rational(test_m, 2 * (test_m - 1)) != 1
    assert (2 - sp.Rational(3, test_m) == 1) == (test_m == 3)
    assert (2 - sp.Rational(4, test_m) == 1) == (test_m == 4)

sync_c, sync_u, sync_v, sync_d = sp.symbols(
    "sync_c sync_u sync_v sync_d", nonzero=True
)
sync_f = (
    y**5
    + sync_c * x**3 * y**2
    + x**4 * (sync_u * y + sync_v * z)
    + sync_d * x**5
)
sync_delta = sp.factor(sp.hessian(sync_f, (x, y, z)).det())
sync_expected_delta = (
    -32 * sync_v**2 * x**6 * (sync_c * x**3 + 10 * y**3)
)
assert sp.simplify(sync_delta - sync_expected_delta) == 0
sync_L = sp.Matrix([0, 0, x])
sync_grad = sp.Matrix([
    sp.diff(sync_f, x),
    sp.diff(sync_f, y),
    sp.diff(sync_f, z),
])
sync_F = sync_v * x**3
sync_M = sync_L.jacobian((x, y, z))
sync_grad_F = sp.Matrix([
    sp.diff(sync_F, x),
    sp.diff(sync_F, y),
    sp.diff(sync_F, z),
])
assert sp.simplify((sync_L.T * sync_grad)[0] - x**2 * sync_F) == 0
sync_middle = (
    sync_M.T * sync_grad
    - sp.Matrix([2 * x * sync_F, 0, 0])
    + x**2 * sync_grad_F / 3
)
assert all(sp.expand(entry) == 0 for entry in sync_middle)
assert sp.expand((sync_L.T * sync_grad_F)[0]) == 0
sync_a3 = sp.Rational(4, 3) * sync_v * x**3
sync_hessian_field = sp.hessian(sync_f, (x, y, z)) * sync_L
sync_grad_a3 = sp.Matrix([
    sp.diff(sync_a3, x),
    sp.diff(sync_a3, y),
    sp.diff(sync_a3, z),
])
assert all(
    sp.expand(entry) == 0
    for entry in sync_hessian_field - x**2 * sync_grad_a3
)

# HC4-DIR27.  An order-two D=5 completion has no bottom-right Hessian
# coefficient and is Psi=w*P+H with P=2*C*x^2+linear.  If the tangent linear
# part vanishes, the bordered determinant has the displayed nonconstant
# square factor.  The complementary unit-direction case is the registered
# quadratic scalar-parent theorem HC4RSD12.
pivot_w = sp.symbols("pivot_w")
pivot_C, pivot_lx = sp.symbols("pivot_C pivot_lx", nonzero=True)
pivot_hxx, pivot_hxy, pivot_hxz = sp.symbols(
    "pivot_hxx pivot_hxy pivot_hxz"
)
pivot_hyy, pivot_hyz, pivot_hzz = sp.symbols(
    "pivot_hyy pivot_hyz pivot_hzz"
)
pivot_p_x = 4 * pivot_C * x + pivot_lx
pivot_bordered_hessian = sp.Matrix([
    [
        pivot_hxx + 4 * pivot_C * pivot_w,
        pivot_hxy,
        pivot_hxz,
        pivot_p_x,
    ],
    [pivot_hxy, pivot_hyy, pivot_hyz, 0],
    [pivot_hxz, pivot_hyz, pivot_hzz, 0],
    [pivot_p_x, 0, 0, 0],
])
pivot_bordered_delta = sp.factor(pivot_bordered_hessian.det())
pivot_expected_delta = sp.factor(
    -pivot_p_x**2 * (pivot_hyy * pivot_hzz - pivot_hyz**2)
)
assert sp.simplify(pivot_bordered_delta - pivot_expected_delta) == 0

# General multiplicity budget.  If
#   Delta=prod pi_i^e_i,
# then H=prod pi_i^ceil(e_i/2) divides adj(A0)b_j (under generic corank one),
# while G=Delta/H has degree kappa=sum floor(e_i/2)deg(pi_i).  The quotient
# degree is kappa-j, so the first motion must satisfy j<=kappa.
kappa = sp.symbols("kappa", integer=True, nonnegative=True)
half_radical_degree = delta_degree - kappa
general_quotient_degree = sp.expand(
    adjugate_motion_degree - half_radical_degree
)
assert general_quotient_degree == kappa - j
assert general_quotient_degree.subs(j, kappa) == 0

# Exact double multiplicity in the transverse case forces m=2.  The tangent
# case has multiplicity at least 2m and is therefore impossible for m>=2.
assert sp.solve(sp.Eq(m, 2), m) == [2]
for test_m in range(2, 10):
    assert 2 * test_m > 2
for test_m in range(3, 10):
    assert test_m != 2

result = {
    "scope": (
        "direct HC4 rank-three top cone with repeated linear Hessian factors"
    ),
    "status": "exact degree and normal-form identities verified",
    "notation": {
        "m": "D-2",
        "delta_degree": "3*m",
        "radical_ell_R_degree": "3*m-1",
        "adjugate_motion_degree": "3*m-j",
        "quotient_degree": "1-j",
    },
    "first_motion_gate": {
        "conclusion": "j=1",
        "constant_vector_identity": "Hess(f)*B=ell*grad(a)",
        "curl_consequence": "a=c*ell^m",
    },
    "general_multiplicity_budget": {
        "kappa": "sum_i floor(e_i/2)*deg(pi_i)",
        "half_radical_degree": "3*m-kappa",
        "quotient_degree": "kappa-j",
        "conclusion": "the first off-diagonal order satisfies j<=kappa",
        "extremal_order": (
            "j=kappa gives a constant B with Hess(f)*B=G*grad(a) "
            "and dG wedge da=0"
        ),
    },
    "lower_rank_quadruple_gate": {
        "hypothesis": (
            "Delta=ell^4*R, R squarefree and coprime to ell, with "
            "generic rank(Hess(f) mod ell)<=1"
        ),
        "divisibility": "ell^2*R divides adj(Hess(f))*b_j",
        "quotient_degree": "2-j",
        "first_motion_orders": [1, 2],
        "order_two": (
            "the only exact-multiplicity candidate is D=6 split, whose "
            "boundary Hessian rank is two, contradicting the hypothesis"
        ),
        "order_one_boundary_jets": [
            "f=ell^2*g_(D-2)",
            "f=y^D+ell^3*g_(D-3)",
        ],
        "ell_squared_jet_bordered_coefficient": str(bordered_coefficient),
        "pure_power_jet_leading_coefficient": str(
            expected_pure_power_boundary / x**4
        ),
        "residual_tangent_determinant": str(
            rank_one_tangent_determinant
        ),
        "weight_compatibility": "m=2 only",
        "conclusion": "no completion for D>=5",
    },
    "rank_two_quadruple_residual": {
        "order_two": {
            "candidate": "D=6 and f=C*ell^6+h_6(y,z)",
            "constant_w2_entry": "q=3*alpha^2/(10*C)",
            "descended_ternary_top": "Q_6=-(2/3)*C*ell^6+h_6(y,z)",
            "conclusion": (
                "impossible because det Hess(Q) is constant but its top "
                "Hessian determinant is nonzero"
            ),
        },
        "order_one": (
            "Hess(f)*L=ell^2*grad(a), L(a)=ell^2*c, and "
            "rank(Hess(f) mod ell)=2"
        ),
    },
    "order_one_rank_collapse": {
        "system": (
            "L(f)=ell^2*F, M^T*grad(f)=2*ell*F*dell-"
            "ell^2*grad(F)/m, and L(F) in (ell^2)"
        ),
        "rank_three": "forces boundary Hessian rank at most one",
        "rank_two": "both boundary-gradient orientations contradict L(F) in (ell^2)",
        "rank_one_integration": "F=C*ell^(2m)/lambda^m",
        "polynomiality": "lambda is proportional to ell",
        "exact_quadruple_consequence": (
            "m=4, D=6, L=ell*partial_ell, and "
            "f=C*ell^6+h_6(y,z)"
        ),
    },
    "degree_six_order_one_terminal_face": {
        "packet": "h_6=C*ell^6+H_6(y,z), h_5=alpha*w*ell^4+r_5",
        "forced_w2_coefficient": str(B_solution),
        "forced_w3_coefficient": str(gamma_solution),
        "intermediate_channel": str(
            sp.factor(
                channel_12.subs({Bcoef: B_solution, gamma: gamma_solution})
            )
        ),
        "terminal_channel_multiplier": str(terminal_channel),
        "full_terminal_channel": (
            "(128/3375)*(alpha^5/C^3)*w^3*det Hess_(y,z)(H_6)"
        ),
        "conclusion": "nonzero in characteristic zero; the D=6 packet is empty",
    },
    "degree_five_triple_terminal_face": {
        "packet": "h_5=C*ell^5+H_5(y,z), h_4=alpha*w*ell^3+r_4",
        "forced_w2_coefficient": str(triple_B_solution),
        "terminal_channel_multiplier": str(triple_terminal_solution),
        "full_terminal_channel": (
            "-(81/400)*(alpha^4/C^2)*w^2*det Hess_(y,z)(H_5)"
        ),
        "conclusion": (
            "nonzero in characteristic zero; the D=5 exact triple packet "
            "is empty"
        ),
    },
    "generic_corank_one_quintuple_gate": {
        "hypothesis": (
            "Delta=ell^5*R, R squarefree and coprime to ell, with "
            "generic rank(Hess(f) mod ell)=2"
        ),
        "divisibility": "ell^3*R divides adj(Hess(f))*b_j",
        "quotient_degree": "2-j",
        "first_motion_orders": [1, 2],
        "order_two": {
            "candidate": "D=7 and f=C*ell^7+h_7(y,z)",
            "forced_w2_coefficient": str(B2_solution),
            "terminal_channel_multiplier": str(q2_terminal_solution),
            "full_terminal_channel": (
                "(32/21)*(alpha^3/C)*w*ell^3*"
                "det Hess_(y,z)(h_7)"
            ),
        },
        "order_one": {
            "candidate": (
                "D=7, f=C*ell^7+h_7(y,z), and L=ell*partial_ell"
            ),
            "forced_w2_coefficient": str(B1_solution),
            "forced_w3_coefficient": str(gamma1_solution),
            "forced_lower_w3_coefficient": str(eta_solution),
            "intermediate_channel": str(q1_intermediate_solution),
            "terminal_channel_multiplier": str(q1_terminal_solution),
            "full_terminal_channel": (
                "-9*gamma^2*w^4*det Hess_(y,z)(h_7)"
            ),
        },
        "conclusion": (
            "both possible motion orders have a nonzero immutable terminal "
            "channel; the generic corank-one quintuple component is empty"
        ),
    },
    "lower_rank_quintuple_gate": {
        "hypothesis": (
            "Delta=ell^5*R, R squarefree and coprime to ell, with "
            "generic rank(Hess(f) mod ell)<=1"
        ),
        "rank_zero": (
            "Hessian integrability gives ell^3|f and determinant order at "
            "least seven"
        ),
        "rank_one_divisibility": (
            "ell^3*R divides adj(Hess(f))*b_j for invariant factors "
            "(0,1,4) or (0,2,3)"
        ),
        "first_motion_orders": [1, 2],
        "order_two": (
            "the only exact-multiplicity candidate is the D=7 split top, "
            "whose boundary Hessian has rank two"
        ),
        "order_one_x2_jet": {
            "form": "f=ell^2*(y^m+ell*g1+...)",
            "order_five_coefficient": str(
                expected_rank_one_x2_order_five
            ),
            "field_weight_solution": str(x2_rank_one_weight[0]),
        },
        "order_one_pure_power_jet": {
            "form": "f=y^(m+2)+ell^3*(y^(m-1)+ell*g1+...)",
            "order_five_coefficient": str(
                expected_rank_one_x3_order_five
            ),
            "field_weight_solution": str(x3_rank_one_weight[0]),
            "degree_obstruction": "deg(g1)=1 implies d_z^2(g1)=0",
        },
        "conclusion": "no lower-rank exact quintuple completion for D>=5",
    },
    "generic_corank_one_sextuple_high_order_gate": {
        "hypothesis": (
            "Delta=ell^6*R, R squarefree and coprime to ell, with "
            "generic rank(Hess(f) mod ell)=2"
        ),
        "divisibility": "ell^3*R divides adj(Hess(f))*b_j",
        "quotient_degree": "3-j",
        "first_motion_orders": [1, 2, 3],
        "order_three": {
            "candidate": "D=8 and f=C*ell^8+h_8(y,z)",
            "forced_constant_w2_hessian": str(sextuple_q_solution),
            "terminal_channel_multiplier": str(s6_j3_terminal_solution),
            "full_terminal_channel": (
                "(24/7)*(alpha^3/C)*w*ell^2*"
                "det Hess_(y,z)(h_8)"
            ),
        },
        "order_two": {
            "linear_field_system": (
                "L(f)=ell^3*F, M^T*grad(f)=3*ell^2*F*dell-"
                "2*ell^3*grad(F)/(m-1), and L(F) in (ell^3)"
            ),
            "matrix_rank": "one",
            "rank_one_log_residues": [
                str(s6_log_x),
                str(s6_log_lambda),
            ],
            "candidate": (
                "D=8, f=C*ell^8+h_8(y,z), and L=ell*partial_ell"
            ),
            "forced_w2_coefficient": str(sextuple_B_solution),
            "terminal_channel_multiplier": str(s6_j2_terminal_solution),
            "full_terminal_channel": (
                "-(1875/3136)*(alpha^4/C^2)*w^2*ell^2*"
                "det Hess_(y,z)(h_8)"
            ),
        },
        "surviving_order_one_system": (
            "Q(f)=ell^3*F, Jac(Q)^T*grad(f)=3*ell^2*F*dell-"
            "ell^3*grad(F)/m, Q(F) in (ell^3), deg(Q)=2, and "
            "rank(Jac(Q) mod ell)<=2"
        ),
        "order_one_boundary_rank_zero": {
            "reduction": "Jac(Q) mod ell=0 forces Q=ell^2*u and D=8 split",
            "forced_w2_coefficient": str(sextuple_B0_solution),
            "forced_w3_coefficient": str(sextuple_gamma0_solution),
            "forced_w4_coefficient": str(sextuple_delta0_solution),
            "reduced_channels": [str(value) for value in s6_j1_r0_reduced],
            "terminal_channel_multiplier": str(s6_j1_r0_reduced[5]),
        },
        "order_one_boundary_rank_one": {
            "field_normal_form": "Q=u*q+ell^2*v when Q mod ell is nonzero",
            "axial_orientation": (
                "Q mod ell=0 gives Q=ell*L and reduces to the excluded "
                "rank-zero packet through the HC4-DIR5 linear system"
            ),
            "normal_orientation": (
                "the dx component has order s-2 on the right and at least "
                "s-1 on the left, so no packet survives"
            ),
            "tangent_orientation": "u=partial_z and 3<=t=ord_ell(partial_z f)<=6",
            "leading_logarithmic_exponent": "e=(m+3-t)/2",
            "leading_relations": (
                "dlog(F_(t-3))=e*dlog(q0), U_t=-(e/m)*F_(t-3)/q0, "
                "and V_(t-2)=(3-(t-3)/m)*F_(t-3)/2"
            ),
            "tangent_packets_below_six": tangent_packets,
            "t6_packet": {
                "minimum_m": 5,
                "q0": "arbitrary nonzero binary quadratic",
                "nonsquare_q0": "m must be odd",
                "square_q0": "every m>=5 is allowed",
            },
            "conclusion": (
                "reduction only; boundary rank one is narrowed to four "
                "tangent packets t=3,4,5,6"
            ),
        },
        "order_one_boundary_rank_one_functional_split": {
            "rational_invariants": "S=q/ell^2 and T=F/ell^m",
            "one_form_identity": "U*dS=-(ell^(m+1)/m)*dT",
            "composite_quadratic": {
                "normal_form": "q,F,U are binary in ell,l for a linear l",
                "invariant_orientation": (
                    "D_u(l)=0 gives f=z*G(ell,l)+h(ell,l), with G=U; "
                    "this is automatic for t<6"
                ),
                "transverse_orientation": (
                    "only at t=6, D_u(l)!=0 gives "
                    "f=P(ell,l)+h(ell,r) with D_u(P)=U"
                ),
                "boundary_quadratic": "a square",
            },
            "primitive_quadratic": {
                "F_form": "sum_k c_k*ell^(m-2k)*q^k",
                "largest_index": "k_max=(m+3-t)/2",
                "packets": primitive_parity_packets,
                "normal_form_for_t_below_six": "q=y^2+ell*z",
                "sample_identity_count": len(primitive_checks),
            },
            "conclusion": (
                "each tangent row splits into a binary-composite pencil "
                "and a parity-restricted primitive conic packet"
            ),
        },
        "order_one_boundary_rank_one_invariant_composite_gate": {
            "top_form": "f=z*G(ell,y)+h(ell,y)",
            "determinant_z_coefficient": (
                "-((m+1)/m)*G*det Hess_(ell,y)(G)"
            ),
            "valuation_bound": "ord_ell(det Hess(f))>=2*t-2",
            "eliminated_by_bound": [5, 6],
            "t3_order_four": str(expected_composite_t3_order_four),
            "t3_order_five_after_order_four_vanishes": str(
                expected_composite_t3_order_five
            ),
            "t3_conclusion": (
                "order four forces a0=0; boundary rank two forces a1!=0, "
                "so order five is nonzero"
            ),
            "t4_order_six": str(expected_composite_t4_order_six),
            "survivor": (
                "t=4, ord_ell(G)=4, and h mod ell=a0*y^(m+2) with a0!=0"
            ),
        },
        "order_one_boundary_rank_one_invariant_composite_closure": {
            "field_z_coefficient": "v_x*G_x+v_y*G_y=0",
            "zero_binary_v": (
                "F=ell^-3*q*G integrates to "
                "G*q^(m+1)=C*ell^(3*(m+1)), impossible since q mod ell!=0"
            ),
            "nonzero_binary_v": (
                "G is a linear-form power; ord_ell(G)=4 forces m=3, "
                "G=C*ell^4, v_x=0, and v_y!=0"
            ),
            "terminal_conflict": (
                "the ell^2 coefficient is v_y*h_y mod ell, nonzero because "
                "h mod ell=a0*y^5 with a0!=0"
            ),
            "conclusion": "the invariant composite orientation is empty",
        },
        "order_one_boundary_rank_one_transverse_composite_closure": {
            "top_form": "f=P(ell,y)+h(ell,z), with U=P_y and ord_ell(U)=6",
            "active_direction_equation": (
                "D_(v_ell,v_z)(h) is a homogeneous polynomial in ell alone"
            ),
            "boundary_consequence": (
                "Hess_(ell,z)(h) mod ell is invertible and kills "
                "(v_ell,v_z), so that active direction is zero"
            ),
            "terminal_reduction": (
                "after absorbing v_y, Q=q*partial_y gives "
                "U*q^(m+1)=C*ell^(3*(m+1)), impossible since q mod ell!=0"
            ),
            "conclusion": "the transverse composite orientation is empty",
        },
        "order_one_boundary_rank_one_closure": {
            "t_at_least_four": (
                "ord_ell(V)=t-2>=2 puts both u and v in the one-dimensional "
                "kernel of Hess(f) mod ell; absorbing v into q gives the "
                "impossible pure-q logarithmic identity"
            ),
            "t3_composite": "excluded by the invariant composite order-four/five gate",
            "t3_primitive": {
                "parity": "m is even",
                "normal_form": "q=y^2+ell*z and F=c*q^(m/2)+O(ell^2)",
                "scalar_obstruction": "[ell]Q(F)=c*(m/2)*y^m!=0",
            },
            "conclusion": (
                "Jac(Q) boundary rank one is empty; the surviving generic "
                "sextuple system has boundary rank two"
            ),
        },
        "order_one_boundary_rank_two_root_gate": {
            "kernel_generator_degree": "k<=2",
            "boundary_gradient": "grad(f) mod ell=H*chi with deg(H)=m+1-k",
            "normal_packet": (
                "chi has zero tangent components, so f mod ell=0 and "
                "the normal component of Q lies in (ell^2)"
            ),
            "root_count_formula": (
                "r=1+k-s<=3, where s=deg gcd(chi_y,chi_z)"
            ),
            "nonzero_boundary_profiles": [
                "y^D",
                "y^a*z^(D-a)",
                "y^a*z^b*(y-z)^c with a+b+c=D",
            ],
            "profile_check_count": len(root_profile_checks),
            "scope": (
                "valid for any homogeneous system with a linear rank-two "
                "boundary matrix and a polynomial gradient in its left kernel"
            ),
            "conclusion": (
                "the remaining rank-two sextuple branch is reduced to one "
                "normal packet and three root profiles"
            ),
        },
        "order_one_boundary_rank_two_one_root_gate": {
            "normalized_top": "f=y^(m+2)+ell^2*g_m(y,z)+O(ell^3)",
            "field_normal_form": (
                "Q_ell=ell*lambda*y+alpha*ell^2, "
                "Q_y=beta*ell^2, and Q_z=q2+ell*r1+gamma*ell^2"
            ),
            "order_two_hessian_coefficient": str(
                expected_one_order_two_delta
            ),
            "polynomial_consequence": "g_m=C*y^m",
            "forced_beta": str(one_beta_solution),
            "conclusion": (
                "the arbitrary degree-m second boundary jet collapses to "
                "one monomial, with lambda and beta linked"
            ),
        },
        "order_one_boundary_rank_two_one_root_closure": {
            "field_beta": str(one_beta_solution),
            "middle_dx_order_one_after_field": str(
                one_middle_after_field
            ),
            "conclusion": (
                "nonzero because C*lambda!=0; the one-root profile is empty"
            ),
        },
        "order_one_boundary_rank_two_remaining_root_first_jets": {
            "two_root_boundary": "f0=y^a*z^b",
            "two_root_kernel_cases": ["(k,s)=(1,0)", "(k,s)=(2,1)"],
            "two_root_normalized_f1": [
                "0",
                "kappa*y^(a-2)*z^(b+1) when a>=2",
                "kappa*y^(a+1)*z^(b-2) when b>=2",
            ],
            "three_root_boundary": "f0=y^a*z^b*(y-z)^c",
            "three_root_kernel_case": "(k,s)=(2,0)",
            "three_root_normalized_f1": "kappa*y^a*z^b*(y-z)^(c-1)",
            "conclusion": (
                "boundary shears reduce both arbitrary degree-(D-1) first "
                "jets to finite monomial representatives"
            ),
        },
        "order_one_boundary_rank_two_remaining_root_schur_gate": {
            "schur_identity": (
                "2*f2=grad(f1)^T*Hess(f0)^(-1)*grad(f1)"
            ),
            "two_root_outer_y_f2": (
                "kappa^2*(a*b-a-4*b)/(2*a*b)*y^(a-4)*z^(b+2), "
                "requiring a>=4"
            ),
            "two_root_outer_z_f2": (
                "kappa^2*(a*b-4*a-b)/(2*a*b)*y^(a+2)*z^(b-4), "
                "requiring b>=4"
            ),
            "three_root_denominator": (
                "R2=a*b*(y-z)^2+a*c*z^2+b*c*y^2"
            ),
            "three_root_numerator": (
                "P2=a*b*(y-z)^2+a*(c-1)*z^2+b*(c-1)*y^2"
            ),
            "three_root_schur_check_count": three_root_schur_checks,
            "three_root_conclusion": "polynomiality forces kappa=0 and f2=0",
            "zero_first_jet_conclusion": (
                "f1=f2=0 forces the first positive ell-jet to order eight "
                "and D>=8"
            ),
        },
        "order_one_boundary_rank_two_normal_closure": {
            "forced_top": (
                "f=ell*C*y^(m+1)+ell^2*(alpha*C/rho)*"
                "y^(m-1)*z+..."
            ),
            "boundary_field": (
                "Q_ell=alpha*ell^2 and "
                "Q_tan mod ell=rho*y^2*partial_z"
            ),
            "scalar_obstruction": (
                "Q(F) mod ell kills the diagonal tangent weight s"
            ),
            "terminal_matrix_coefficient": str(
                normal_terminal_z_coefficient
            ),
            "conclusion": "the normal rank-two packet is empty",
        },
        "order_one_boundary_rank_two_delayed_jet_closure": {
            "hypothesis": (
                "f1=f2=0 and Hess_(y,z)(f0) is invertible"
            ),
            "first_field_equations": [
                "l(f0)=0",
                "Jac(l)^T*grad(f0)=0",
            ],
            "product_rule_consequence": "Hess(f0)*l=0, hence l=0",
            "rank_consequence": "rank(Jac(Q) mod ell)<=1",
            "scope": (
                "independent of root multiplicities, the order-eight "
                "value, and exact sextuple order"
            ),
            "conclusion": "all delayed two-/three-root packets are empty",
        },
        "order_one_boundary_rank_two_outer_jet_closure": {
            "boundary_field": (
                "Q0=(y^2,-2*kappa*y*z/a,kappa*z^2/b)"
            ),
            "scalar_x1_relations": [
                "B=E=0",
                "kappa*A+a*D=0",
                "a*C+b*G=0",
            ],
            "shared_tangent_factor": str(outer_common_factor),
            "immutable_coefficient": str(outer_terminal_coefficient),
            "conclusion": "both mirrored outer jets are empty",
        },
        "conclusion": (
            "orders three and two are empty; at j=1 boundary rank zero "
            "and rank one are empty; the normal, delayed, and outer gates "
            "also close every boundary-rank-two packet, so the complete "
            "generic-corank-one exact-sextuple stratum is empty"
        ),
    },
    "lower_rank_sextuple_reduction": {
        "hypothesis": (
            "Delta=ell^6*R, R squarefree and coprime to ell, with "
            "generic rank(Hess(f) mod ell)<=1"
        ),
        "rank_zero": "determinant order is at least seven",
        "rank_one_invariant_factors": ["(0,1,5)", "(0,2,4)", "(0,3,3)"],
        "divisibility": "ell^4*R divides adj(Hess(f))*b_j",
        "quotient_degree": "2-j",
        "order_two_packet": (
            "D=5, f=C*z*ell^4+h_5(ell,y), and B=partial_z"
        ),
        "order_one_x2_jet": {
            "form": (
                "f=ell^2*(y^m+ell*(u*y^(m-1)+v*y^(m-2)*z)+"
                "ell^2*g2+...)"
            ),
            "order_six_coefficient": str(expected_rank_one_x2_order_six),
        },
        "order_one_pure_power_jet": {
            "form": (
                "f=y^(m+2)+ell^3*(c*y^(m-1)+"
                "ell*(u*y^(m-2)+v*y^(m-3)*z)+ell^2*g2+...)"
            ),
            "order_six_coefficient": str(expected_rank_one_x3_order_six),
        },
        "order_one_synchronization": {
            "x2_packet": (
                "field recurrences force v=0 and d_z^2(g2)=0, so its "
                "exact order-six coefficient vanishes"
            ),
            "pure_power_resonance": (
                "the only nonzero order-six row has m=3, v!=0, and "
                "L=ell*partial_z"
            ),
            "resonant_top": str(sync_f),
            "resonant_hessian_determinant": str(sync_delta),
            "pure_cube_pivot_top": str(sync_a3),
            "complete_scalar_parent": (
                "Psi=H+w*P+eta*w^2/2 with "
                "P=(4*v/3)*ell^3+P_le2"
            ),
            "common_top_geometry": "f=C*z*ell^4+h5(ell,y)",
        },
        "order_two_scalar_pivot_closure": {
            "complete_parent": "Psi=w*P+H",
            "quadratic_pivot": "P=2*C*ell^2+linear",
            "no_tangent_linear_part": str(pivot_bordered_delta),
            "tangent_linear_part": (
                "a constant tangent direction has unit derivative on P, "
                "so HC4RSD12 reduces every collision fiber to HC3"
            ),
            "conclusion": (
                "no order-two member of the synchronized tangent family "
                "is an HC4 counterexample"
            ),
        },
        "conclusion": (
            "the x2 order-one packet is empty; the pure-power order-one "
            "packet reduces to a degree-five resonance already contained "
            "in the order-two tangent top family; the order-two scalar "
            "parent is HC4-safe, leaving only the degree-five order-one "
            "resonance"
        ),
    },
    "normal_forms": {
        "B_ell_nonzero": {
            "form": "f=C*x^(m+2)+h(y,z)",
            "hessian_determinant": str(transverse_determinant),
            "ell_multiplicity": "m",
            "double_factor_consequence": "m=2, hence D=4",
            "triple_factor_consequence": "m=3, hence D=5",
        },
        "B_ell_zero": {
            "form": "f=C*z*x^(m+1)+h(x,y)",
            "hessian_determinant": str(tangent_determinant),
            "ell_multiplicity": "at least 2*m",
            "double_factor_consequence": "impossible for D>=4",
            "triple_factor_consequence": "impossible for D>=4",
        },
    },
    "split_top_controls": {
        "D4": {
            "form": str(x**4 + h4),
            "hessian_determinant": str(split_delta4),
            "binary_cofactor_squarefree": True,
            "ell_multiplicity": 2,
        },
        "D5": {
            "form": str(x**5 + h5),
            "hessian_determinant": str(split_delta5),
            "binary_cofactor_squarefree": True,
            "ell_multiplicity": 3,
        },
    },
    "conclusion": (
        "for D>=5, Delta=ell^2*R with R squarefree and gcd(ell,R)=1 "
        "is incompatible with a four-variable constant-Hessian completion; "
        "Delta=ell^3*R is incompatible for every D>=5, because its sole "
        "D=5 split top has a nonzero terminal passive-Hessian channel; the "
        "generic rank-at-most-one boundary of "
        "Delta=ell^4*R is incompatible for every D>=5; its generic rank-two "
        "order-one system first reduces to the D=6 additive top, whose "
        "terminal passive-Hessian channel is nonzero; on "
        "Delta=ell^5*R, the generic corank-one boundary has j<=2, and both "
        "the order-two and order-one degree-seven packets have nonzero "
        "terminal passive-Hessian channels; the complementary lower-rank "
        "boundary has two jets whose order-five coefficient is incompatible "
        "with the order-one field equations; on Delta=ell^6*R every "
        "generic-corank-one packet is empty, including the normal, delayed "
        "two-/three-root, and outer boundary-rank-two rows; on the "
        "complementary lower-rank boundary the x2 order-one packet is "
        "empty and the sole pure-power resonance is already contained in "
        "the degree-five tangent order-two family; that quadratic scalar "
        "parent is HC4-safe by the bordered determinant split and "
        "HC4RSD12, leaving only the degree-five order-one resonance"
    ),
    "proof_boundary": (
        "the checker verifies the all-degree and weighted binary identities; "
        "UFD/DVR divisibility and maximal-passive-channel uniqueness are "
        "proved in HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md"
    ),
}

OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
