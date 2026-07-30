#!/usr/bin/env python3
"""Exact identities for the shared JC(2)--HC(4) isotropic boundary package.

This checker has three independent blocks.

1.  For the cotangent potential

        Psi=t*P(x,y)+m*Q(x,y)+H(x,y),

    its four-variable Hessian determinant is Jac(P,Q)^2, independently of
    the second derivatives and of t,m.

2.  For the more general isotropic form Psi=t*P(x,y)+Phi(x,y,m), the t^2
    Schur coefficient vanishes and the t coefficient is

        -Phi_mm * grad(P)^T adj(Hess(P)) grad(P).

    Thus the cotangent branch Phi_mm=0 is one exact factor of the first
    nontrivial Schur remainder.

3.  In an adapted index-two boundary chart g=r^2*ell(r,z), the reduced
    conormal residue r^(-1)*partial_r(g) mod r is 2*ell(0,z).  Hence the
    quartic packet's odd-square multiplier ell is half of the normal
    residue, and its square is the leading cotangent-Hessian coefficient.

The script proves only these algebraic identities.  It does not prove that
an arbitrary HC(4) counterexample admits the required polynomial isotropic
flag, nor that the endpoint residue-pairing class is nonzero.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# 1. Universal cotangent determinant
# ---------------------------------------------------------------------------

P_x, P_y, Q_x, Q_y = sp.symbols("P_x P_y Q_x Q_y")
A_xx, A_xy, A_yy = sp.symbols("A_xx A_xy A_yy")

cotangent_hessian = sp.Matrix(
    [
        [A_xx, A_xy, P_x, Q_x],
        [A_xy, A_yy, P_y, Q_y],
        [P_x, P_y, 0, 0],
        [Q_x, Q_y, 0, 0],
    ]
)
jacobian = P_x * Q_y - P_y * Q_x

assert sp.factor(cotangent_hessian.det() - jacobian**2) == 0


# ---------------------------------------------------------------------------
# 2. First isotropic Schur remainder
# ---------------------------------------------------------------------------

t = sp.symbols("t")
P_xx, P_xy, P_yy = sp.symbols("P_xx P_xy P_yy")
Phi_xx, Phi_xy, Phi_xm = sp.symbols("Phi_xx Phi_xy Phi_xm")
Phi_yy, Phi_ym, Phi_mm = sp.symbols("Phi_yy Phi_ym Phi_mm")

isotropic_hessian = sp.Matrix(
    [
        [0, P_x, P_y, 0],
        [P_x, t * P_xx + Phi_xx, t * P_xy + Phi_xy, Phi_xm],
        [P_y, t * P_xy + Phi_xy, t * P_yy + Phi_yy, Phi_ym],
        [0, Phi_xm, Phi_ym, Phi_mm],
    ]
)
determinant_in_t = sp.Poly(sp.expand(isotropic_hessian.det()), t)
binary_bordered_remainder = (
    P_x**2 * P_yy - 2 * P_x * P_y * P_xy + P_y**2 * P_xx
)

assert determinant_in_t.degree() <= 1
assert determinant_in_t.coeff_monomial(t**2) == 0
assert sp.factor(
    determinant_in_t.coeff_monomial(t)
    + Phi_mm * binary_bordered_remainder
) == 0


# ---------------------------------------------------------------------------
# 3. Index-two reduced conormal residue and its square
# ---------------------------------------------------------------------------

r, z = sp.symbols("r z")
ell_0 = sp.Function("ell_0")(z)
ell_1 = sp.Function("ell_1")(z)
ell = ell_0 + r * ell_1
g_pullback = sp.expand(r**2 * ell)

reduced_normal_residue = sp.expand(sp.diff(g_pullback, r) / r).subs(r, 0)
assert sp.simplify(reduced_normal_residue - 2 * ell_0) == 0
assert sp.expand(reduced_normal_residue**2 - 4 * ell_0**2) == 0

T = sp.symbols("T")
cusp_ell = 4 * r - 9 * T**2
cusp_pullback = r**2 * cusp_ell
cusp_residue = sp.expand(sp.diff(cusp_pullback, r) / r).subs(r, 0)
assert cusp_residue == -18 * T**2


print("PASS: cotangent determinant = Jac(P,Q)^2")
print("PASS: first isotropic Schur remainder = -Phi_mm*R(P)")
print("PASS: index-two reduced conormal residue = 2*ell")
print("PASS: clean cusp initial conormal residue = -18*T^2")
