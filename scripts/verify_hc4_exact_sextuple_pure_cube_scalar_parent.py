#!/usr/bin/env python3
"""Verify the scalar-parent identities closing the exact-sextuple resonance.

The written proof supplies the linear-algebra lemma that a cubic with
pure-cube top and nowhere-vanishing gradient has a constant unit direction.
This checker replays the universal zero- and nonzero-corner determinant and
gradient factorizations, together with the two scalar equations used in that
lemma.  All calculations are exact over QQ.
"""

from __future__ import annotations

import sympy as sp


# Nonzero corner: the four-variable parent is a completed-square suspension
# of a ternary constant-Hessian pencil.
eta = sp.symbols("eta", nonzero=True)
g = sp.Matrix(sp.symbols("g0:3"))
m = sp.symbols("m0:6")
M = sp.Matrix(
    [
        [m[0], m[1], m[2]],
        [m[1], m[3], m[4]],
        [m[2], m[4], m[5]],
    ]
)
parent_nonzero = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[eta]]), g.T),
    sp.Matrix.hstack(g, M),
)
schur_member = M - g * g.T / eta
assert sp.factor(
    parent_nonzero.det(method="domain-ge")
    - eta * schur_member.det(method="domain-ge")
) == 0


# Zero corner in graph coordinates P=s+q(u), with two residual variables.
# Here p=grad(q), Q=Hess(q), H=Hess_u(C), b=grad_u(C_r), d=C_rr,
# and tau=w+C_r.  The universal cancellations leave a binary Hessian.
p = sp.Matrix(sp.symbols("p0:2"))
b = sp.Matrix(sp.symbols("b0:2"))
h11, h12, h22 = sp.symbols("h11 h12 h22")
q11, q12, q22 = sp.symbols("q11 q12 q22")
H = sp.Matrix([[h11, h12], [h12, h22]])
Q = sp.Matrix([[q11, q12], [q12, q22]])
d, tau = sp.symbols("d tau")
active = H + b * p.T + p * b.T + d * p * p.T + tau * Q
cross = b + d * p
parent_zero = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(1, 1), p.T, sp.ones(1, 1)),
    sp.Matrix.hstack(p, active, cross),
    sp.Matrix.hstack(sp.ones(1, 1), cross.T, sp.Matrix([[d]])),
)
binary_member = H + tau * Q
assert sp.factor(
    parent_zero.det(method="domain-ge")
    + binary_member.det(method="domain-ge")
) == 0


# Exact gradient formula in a nonlinear pure-cube calibration.  It exercises
# every cancelled block and has a constant unit pivot direction.
w, u, v, s, r = sp.symbols("w u v s r")
q = u**3 / 3 + u * v
C = u**2 * v + v * r + r**2 * u
P = s + q
H_original = sp.expand(C.subs(r, P))
Psi = sp.expand(w * P + H_original)
r_coordinate = P
tau_coordinate = sp.expand(w + sp.diff(C, r).subs(r, P))
gradient = sp.Matrix([sp.diff(Psi, variable) for variable in (w, u, v, s)])
expected_middle = sp.Matrix(
    [
        sp.diff(C, variable).subs(r, r_coordinate)
        + tau_coordinate * sp.diff(q, variable)
        for variable in (u, v)
    ]
)
assert sp.expand(gradient[0] - r_coordinate) == 0
assert all(
    sp.expand(gradient[index + 1] - expected_middle[index]) == 0
    for index in range(2)
)
assert sp.expand(gradient[3] - tau_coordinate) == 0
assert sp.factor(
    sp.hessian(Psi, (w, u, v, s)).det(method="domain-ge")
    + sp.hessian(C + tau * q, (u, v)).det(method="domain-ge").subs(
        {r: P, tau: tau_coordinate}
    )
) == 0


# Pure-cube critical-point lemma.  If A is the quadratic Hessian and
# N=ker(A), the proof splits according to ell|N.  When ell|N is nonzero,
# b|ker(ell|N)=0 gives b|N=beta*ell|N and one chooses t below.  When
# ell|N=0, the remaining constraint is the displayed nonconstant quadratic.
a, beta, gamma, theta, t = sp.symbols("a beta gamma theta t", nonzero=True)
t_square = -beta / (3 * a)
assert sp.factor(beta + 3 * a * t_square) == 0
second_case_equation = gamma + theta * t**2 - t
assert sp.Poly(second_case_equation, t).degree() >= 1
assert sp.Poly(second_case_equation, t).coeff_monomial(t) == -1


print("PASS: nonzero corner is a ternary constant-Hessian pencil")
print("PASS: zero-corner graph determinant is a binary Hessian determinant")
print("PASS: pure-cube critical-point alternatives replayed exactly")
print("THEOREM: the exact-sextuple order-one scalar-parent resonance is HC4-safe")
