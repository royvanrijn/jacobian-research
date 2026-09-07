#!/usr/bin/env python3
"""Exact symbolic regression for the transport-curvature frontier.

The script verifies:
  * the common elementary frame used by JC(3), Long/GMC(3), and GVC(3);
  * the block identities behind the universal rank-one-profile Hessian formula;
  * projective curvature formulas for the foundational and two Cohn normal forms;
  * the all-degree highest-Hessian coefficients for the two Cohn profiles;
  * the one-latitude composition formula, including arbitrary amplitude;
  * the minimal second-contact Hopf collapse for the second Cohn profile.

The no-negative-eigenfunction and geodesic impossibility arguments are short
symbolic proofs recorded in the accompanying markdown note; they are not
finite searches.
"""
from __future__ import annotations

import sympy as sp

# ---------------------------------------------------------------------------
# 1. Common elementary transport frame.
# ---------------------------------------------------------------------------
xi, eta = sp.symbols("xi eta")
M = sp.Matrix([[1 + xi * eta, eta], [xi, 1]])
E12 = sp.Matrix([[1, eta], [0, 1]])
E21 = sp.Matrix([[1, 0], [xi, 1]])
assert M.det() == 1
assert M == E12 * E21

# ---------------------------------------------------------------------------
# 2. Universal block algebra for L = c(X).U.
# ---------------------------------------------------------------------------
a, b, t, m, d = sp.symbols("a b t m d")
j11, j12, j21, j22 = sp.symbols("j11 j12 j21 j22")
s11, s12, s22 = sp.symbols("s11 s12 s22")
J = sp.Matrix([[j11, j12], [j21, j22]])
S = sp.Matrix([[s11, s12], [s12, s22]])
c = sp.Matrix([a, b])
U = sp.Matrix([t, m])
L = (c.T * U)[0]
v = sp.Matrix.vstack(J.T * U, c)
K = sp.Matrix.vstack(
    sp.Matrix.hstack(S, J.T),
    sp.Matrix.hstack(J, sp.zeros(2)),
)
Kinv = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(2), J.inv()),
    sp.Matrix.hstack(J.T.inv(), -J.T.inv() * S * J.inv()),
)
assert sp.simplify(K * Kinv - sp.eye(4)) == sp.zeros(4)
assert sp.factor(K.det() - J.det() ** 2) == 0
w = sp.simplify(J.inv() * c)
q = sp.simplify((w.T * S * w)[0])
assert sp.factor((v.T * Kinv * v)[0] - (2 * L - q)) == 0

# The matrix determinant lemma now gives, for Phi(z)=lambda*z**d,
# det Hess Phi(L) = (lambda*d)^4 * L^(4d-5) * det(J)^2
#                   * ((2d-1)L - (d-1)q).

# ---------------------------------------------------------------------------
# 3. Projective curvature kappa = det(w, D_w w).
# ---------------------------------------------------------------------------
x, y, gamma, k = sp.symbols("x y gamma k", nonzero=True)

def radial_vector_and_curvature(c1: sp.Expr, c2: sp.Expr):
    cvec = sp.Matrix([c1, c2])
    jac = cvec.jacobian([x, y])
    radial = sp.simplify(jac.inv() * cvec)
    wx, wy = radial
    acceleration = sp.Matrix([
        sp.diff(wx, x) * wx + sp.diff(wx, y) * wy,
        sp.diff(wy, x) * wx + sp.diff(wy, y) * wy,
    ])
    curvature = sp.factor(wx * acceleration[1] - wy * acceleration[0])
    return radial, acceleration, curvature

# Foundational unit-affine row (1+xy, x).
w_found, acc_found, kap_found = radial_vector_and_curvature(1 + x * y, x)
assert w_found == sp.Matrix([x, 1 / x])
assert kap_found == -2

# Cohn/hyperbolic family (gamma+xy, x^k).
w_b, acc_b, kap_b = radial_vector_and_curvature(gamma + x * y, x**k)
kap_b_expected = (
    gamma * k * (k - 3) + (k - 1) * (k - 2) * x * y
) / k**3
assert sp.factor(kap_b - kap_b_expected) == 0

# Cohn/quadratic family (x^2, delta+x*f(y)).
delta = sp.symbols("delta", nonzero=True)
f = sp.Function("f")(y)
w_a, acc_a, kap_a = radial_vector_and_curvature(x**2, delta + x * f)
kap_a_expected = -(
    (2 * delta + x * f) ** 2 * sp.diff(f, y, 2)
    + 2 * delta * x * sp.diff(f, y) ** 2
) / (8 * x * sp.diff(f, y) ** 3)
assert sp.factor(kap_a - kap_a_expected) == 0

# ---------------------------------------------------------------------------
# 4. Exact highest-Hessian coefficients for constant profile amplitude.
# ---------------------------------------------------------------------------
lam = sp.symbols("lam", nonzero=True)

# Generic slice t-axis formula, with abstract first/second derivatives.
ax, ay, bx, by = sp.symbols("ax ay bx by")
axx, axy, ayy, bxx, bxy, byy = sp.symbols(
    "axx axy ayy bxx bxy byy"
)
grad_t = sp.Matrix([ax * t, ay * t, a, b])
HL_t = sp.Matrix([
    [axx * t, axy * t, ax, bx],
    [axy * t, ayy * t, ay, by],
    [ax, ay, 0, 0],
    [bx, by, 0, 0],
])
core_t = (d - 1) * (grad_t * grad_t.T) + a * t * HL_t

# Substitute the hyperbolic Cohn family on y=m=0.
subs_b = {
    a: gamma,
    ax: 0,
    ay: x,
    axx: 0,
    axy: 1,
    ayy: 0,
    b: x**k,
    bx: k * x ** (k - 1),
    by: 0,
    bxx: k * (k - 1) * x ** (k - 2),
    bxy: 0,
    byy: 0,
}
core_b_det = sp.factor(core_t.det().subs(subs_b))
core_b_expected = (
    gamma**4
    * k
    * (2 * d * (k - 1) - k + 2)
    * x ** (2 * k)
    * t**4
)
assert sp.factor(core_b_det - core_b_expected) == 0
# Multiplying by (lam*d*(gamma*t)^(d-2))^4 gives the theorem coefficient.

# Generic m-axis formula for the quadratic Cohn family.
grad_m = sp.Matrix([bx * m, by * m, a, b])
HL_m = sp.Matrix([
    [bxx * m, bxy * m, ax, bx],
    [bxy * m, byy * m, ay, by],
    [ax, ay, 0, 0],
    [bx, by, 0, 0],
])
core_m = (d - 1) * (grad_m * grad_m.T) + b * m * HL_m
fp = sp.diff(f, y)
fpp = sp.diff(f, y, 2)
subs_a = {
    a: x**2,
    ax: 2 * x,
    ay: 0,
    b: delta + x * f,
    bx: f,
    by: x * fp,
    bxx: 0,
    bxy: fp,
    byy: x * fpp,
}
core_a_det = sp.factor(core_m.det().subs(subs_a))
bracket_a = (
    (d - 1) * (2 * delta + x * f) ** 2 * fpp
    - 2 * x * (2 * d * delta + (3 * d - 1) * x * f) * fp**2
)
core_a_expected = (
    -m**4 * x**3 * (delta + x * f) ** 3 * bracket_a
)
assert sp.factor(core_a_det - core_a_expected) == 0
# Multiplying by (lam*d*((delta+x*f)*m)^(d-2))^4 gives the theorem coefficient.

# ---------------------------------------------------------------------------
# 5. Composition through one latitude h(x,y).
# ---------------------------------------------------------------------------
gx, gy = sp.symbols("gx gy")
hxx, hxy, hyy = sp.symbols("hxx hxy hyy")
ph, phh, pht, phm = sp.symbols("ph phh pht phm")
ptt, ptm, pmm = sp.symbols("ptt ptm pmm")
grad_h = sp.Matrix([gx, gy])
hess_h = sp.Matrix([[hxx, hxy], [hxy, hyy]])
r_h = sp.Matrix([pht, phm])
hess_u = sp.Matrix([[ptt, ptm], [ptm, pmm]])
source_block = phh * grad_h * grad_h.T + ph * hess_h
mixed_block = grad_h * r_h.T
composition_hessian = sp.Matrix.vstack(
    sp.Matrix.hstack(source_block, mixed_block),
    sp.Matrix.hstack(mixed_block.T, hess_u),
)
R_h = sp.expand((grad_h.T * hess_h.adjugate() * grad_h)[0])
hess3 = sp.Matrix([
    [phh, pht, phm],
    [pht, ptt, ptm],
    [phm, ptm, pmm],
])
composition_expected = (
    ph**2 * hess_h.det() * hess_u.det()
    + ph * R_h * hess3.det()
)
assert sp.factor(composition_hessian.det() - composition_expected) == 0

# For Phi=(A(h)t+B(h)m)^d, the 3-variable Hessian determinant is the
# Wronskian square.  Verify the symbolic core using abstract jets.
A0, B0, A1, B1, A2, B2, hh = sp.symbols(
    "A0 B0 A1 B1 A2 B2 hh"
)
LL = A0 * t + B0 * m
LLh = A1 * t + B1 * m
LLhh = A2 * t + B2 * m
grad_LL = sp.Matrix([LLh, A0, B0])
hess_LL = sp.Matrix([
    [LLhh, A1, B1],
    [A1, 0, 0],
    [B1, 0, 0],
])
core3 = (d - 1) * grad_LL * grad_LL.T + LL * hess_LL
W = A0 * B1 - B0 * A1
assert sp.factor(core3.det() + (d - 1) * LL**2 * W**2) == 0
# Multiplying by (d*LL^(d-2))^3 gives
# det Hess_(h,t,m)(LL^d) = -d^3(d-1) LL^(3d-4) W^2.

# The same determinant is independent of all second h-jets for an arbitrary
# amplitude mu(h):
mu0, mu1, mu2 = sp.symbols("mu0 mu1 mu2")
phi_h = mu1 * LL**d + mu0 * d * LL**(d - 1) * LLh
phi_hh = (
    mu2 * LL**d
    + 2 * mu1 * d * LL**(d - 1) * LLh
    + mu0 * d * (d - 1) * LL**(d - 2) * LLh**2
    + mu0 * d * LL**(d - 1) * LLhh
)
phi_ht = (
    mu1 * d * LL**(d - 1) * A0
    + mu0 * d * ((d - 1) * LL**(d - 2) * LLh * A0 + LL**(d - 1) * A1)
)
phi_hm = (
    mu1 * d * LL**(d - 1) * B0
    + mu0 * d * ((d - 1) * LL**(d - 2) * LLh * B0 + LL**(d - 1) * B1)
)
phi_tt = mu0 * d * (d - 1) * LL**(d - 2) * A0**2
phi_tm = mu0 * d * (d - 1) * LL**(d - 2) * A0 * B0
phi_mm = mu0 * d * (d - 1) * LL**(d - 2) * B0**2
hess_mu = sp.Matrix([
    [phi_hh, phi_ht, phi_hm],
    [phi_ht, phi_tt, phi_tm],
    [phi_hm, phi_tm, phi_mm],
])
assert sp.simplify(sp.powsimp(
    hess_mu.det()
    + d**3 * (d - 1) * mu0**3 * W**2 * LL**(3 * d - 4),
    force=True,
)) == 0

# Bordered-Hessian calibration and the easy direction of the affine-latitude
# theorem.
h_fun = sp.Function("h")
alpha, beta, z = sp.symbols("alpha beta z")
f_lat = sp.Function("F")(alpha * x + beta * y)
R_lat = (
    sp.diff(f_lat, x) ** 2 * sp.diff(f_lat, y, 2)
    - 2 * sp.diff(f_lat, x) * sp.diff(f_lat, y) * sp.diff(f_lat, x, y)
    + sp.diff(f_lat, y) ** 2 * sp.diff(f_lat, x, 2)
)
assert sp.simplify(R_lat) == 0
R_xy = sp.factor(
    sp.diff(x * y, x) ** 2 * sp.diff(x * y, y, 2)
    - 2 * sp.diff(x * y, x) * sp.diff(x * y, y) * sp.diff(x * y, x, y)
    + sp.diff(x * y, y) ** 2 * sp.diff(x * y, x, 2)
)
assert R_xy == -2 * x * y

# Exact quadratic regression for Theorem 5.3.  For
# h=a2*x^2+b2*x*y+c2*y^2+d2*x+e2*y+f2, R(h)=0 is exactly the rank-one
# quadratic condition together with alignment of the linear part (the
# zero-quadratic case gives an affine linear h).
a2, b2, c2, d2, e2, f2 = sp.symbols("a2 b2 c2 d2 e2 f2")
h_quad = a2 * x**2 + b2 * x * y + c2 * y**2 + d2 * x + e2 * y + f2
R_quad = sp.Poly(
    sp.diff(h_quad, x) ** 2 * sp.diff(h_quad, y, 2)
    - 2 * sp.diff(h_quad, x) * sp.diff(h_quad, y) * sp.diff(h_quad, x, y)
    + sp.diff(h_quad, y) ** 2 * sp.diff(h_quad, x, 2),
    x,
    y,
)
Delta_quad = 4 * a2 * c2 - b2**2
assert sp.expand(R_quad.coeff_monomial(x**2) - (2 * a2 * Delta_quad)) == 0
assert sp.expand(R_quad.coeff_monomial(x*y) - (2 * b2 * Delta_quad)) == 0
assert sp.expand(R_quad.coeff_monomial(y**2) - (2 * c2 * Delta_quad)) == 0
assert sp.expand(R_quad.coeff_monomial(x) - (2 * d2 * Delta_quad)) == 0
assert sp.expand(R_quad.coeff_monomial(y) - (2 * e2 * Delta_quad)) == 0
assert sp.expand(R_quad.coeff_monomial(1) - (2 * (a2 * e2**2 - b2 * d2 * e2 + c2 * d2**2))) == 0

# ---------------------------------------------------------------------------
# 6. Calibration: polynomial rank collapse is necessary but not sufficient.
# ---------------------------------------------------------------------------
# c=(y, x*y^2), w=(-x,y), lambda=x^2, d=2.  Then w(lambda)=-2 lambda
# and every top coefficient is a polynomial in q=x*y.  Nevertheless the
# full Hessian determinant is nonzero: nonlinear latitude h=x*y retains a
# second-fundamental-form obstruction.
c_cal = sp.Matrix([y, x * y**2])
J_cal = c_cal.jacobian([x, y])
w_cal = sp.simplify(J_cal.inv() * c_cal)
assert w_cal == sp.Matrix([-x, y])
lambda_cal = x**2
assert sp.expand(
    w_cal[0] * sp.diff(lambda_cal, x)
    + w_cal[1] * sp.diff(lambda_cal, y)
    + 2 * lambda_cal
) == 0
G_cal = sp.expand(lambda_cal * (y * t + x * y**2 * m) ** 2)
q_cal = x * y
assert sp.expand(G_cal - (q_cal * t + q_cal**2 * m) ** 2) == 0
det_cal = sp.factor(sp.hessian(G_cal, (x, y, t, m)).det())
assert det_cal == 32 * x**8 * y**8 * (m * x * y + t) ** 3 * (2 * m * x * y + t)

# ---------------------------------------------------------------------------
# 7. Minimal second-contact Hopf collapse for the quadratic Cohn family.
# ---------------------------------------------------------------------------
rho = sp.symbols("rho", nonzero=True)
psi = sp.Function("psi")(y)
# Quadric relation T^2 + rho*x*y = 1 and V = delta + x*psi(y).
# D = 1 - delta^-2*T^2*V^2 after eliminating T^2.
D = sp.expand(
    1 - delta ** -2 * (1 - rho * x * y) * (delta + x * psi) ** 2
)
linear_x = sp.expand(sp.diff(D, x).subs(x, 0))
assert sp.factor(linear_x - (rho * y - 2 * psi / delta)) == 0
# Hence ord_x(D)>=2 iff psi=(rho*delta/2)y.  On that branch:
u = sp.symbols("u")
D_special = sp.factor(
    1 - (1 - rho * u) * (1 + rho * u / 2) ** 2
)
assert sp.factor(D_special - rho**2 * u**2 * (3 + rho * u) / 4) == 0

print("PASS common frame M=E12(eta)E21(xi)")
print("PASS universal block inverse and rank-one-profile determinant formula")
print("PASS foundational projective curvature kappa=-2")
print("PASS hyperbolic Cohn curvature formula")
print("PASS quadratic Cohn curvature formula")
print("PASS all-degree highest-Hessian coefficient: hyperbolic Cohn profile")
print("PASS all-degree highest-Hessian coefficient: quadratic Cohn profile")
print("PASS one-latitude Hessian composition, arbitrary-amplitude, and Wronskian formulas")
print("PASS affine-latitude bordered-Hessian theorem: exact quadratic regression")
print("PASS rank collapse alone is insufficient: nonlinear-latitude Hessian remains")
print("PASS minimal second-contact Hopf branch collapses to a torus-weight profile")
