#!/usr/bin/env python3
"""Verify the exact algebra in the dense cubic--quartic HC(4) reduction.

The conceptual proof uses the Gordan--Noether theorem in dimensions at
most four, HC(3), HC(2), and the degree-at-most-100 plane Jacobian theorem.
This checker verifies the determinant identities between those inputs.

For

    psi = q_2 + h_3 + h_4

write H0=Hess(q_2), A=Hess(h_3), and B=Hess(h_4).  Homogeneous scaling gives

    det(H0 + z*A + z^2*B) = det(H0).

The shared dense common-kernel checker verifies the leading identities for
rank(B)=3,2,1.  They force the Hessian determinant of h_3 restricted to the
constant kernel of B to vanish.  The binary restriction has a constant
kernel by unique factorization.  In the ternary restriction, a second
linear-pencil calculation below handles the possible variable kernel.
Together with Gordan--Noether this supplies a direction v with
D_v h_4=0 and D_v^2 h_3=0.

If v is nonisotropic for q_2, polynomial Schur descent gives HC(3).  In the
isotropic case the potential is

    psi = t*s(x,y) + phi(x,y,m),

after a constant kernel direction m is selected.  The shared checker also
verifies the bordered-Hessian, cotangent-lift, and final HC(2) block
identities.  The new calculation here proves that the binary polynomial

    s = constant + x + alpha*x^2 + beta*x*y + gamma*y^2

has vanishing bordered invariant only when beta=gamma=0; nonzero alpha
makes grad(s) vanish somewhere over the algebraic closure, contradicting a
nonzero constant Hessian determinant.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import runpy

import sympy as sp


SHARED_CHECKER = Path(__file__).with_name(
    "verify_hc4_meng_dense_rank_three_sextic_reduction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    shared = runpy.run_path(str(SHARED_CHECKER))


# The rank-zero leading case, not needed in the rank-three sextic theorem:
# if B=0, the z^4 coefficient is det(A).
z = sp.symbols("z")
H0 = shared["H0"]
A = shared["A"]
rank_zero_polynomial = sp.Poly(
    (H0 + z * A).det(method="berkowitz"),
    z,
)
assert sp.expand(
    rank_zero_polynomial.coeff_monomial(z**4) - A.det()
) == 0


# The residual rank-one quartic case leaves a ternary pencil
#
#     Hess(f_3)(u) + z*Hess(f_2).
#
# Gordan--Noether makes Hess(f_3) constant-kernel.  Its ranks two, one,
# and zero force the displayed restrictions of Hess(f_2).
d11, d12, d22 = sp.symbols("d11 d12 d22")
D2 = sp.Matrix([[d11, d12], [d12, d22]])
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

P_rank_two = sp.zeros(3)
P_rank_two[0:2, 0:2] = D2
ternary_rank_two = sp.Poly(
    (P_rank_two + z * C3).det(method="berkowitz"),
    z,
)
assert sp.expand(
    ternary_rank_two.coeff_monomial(z)
    - D2.det() * c33
) == 0

p11 = sp.symbols("p11")
P_rank_one = sp.diag(p11, 0, 0)
ternary_rank_one = sp.Poly(
    (P_rank_one + z * C3).det(method="berkowitz"),
    z,
)
assert sp.expand(
    ternary_rank_one.coeff_monomial(z**2)
    - p11 * sp.Matrix([[c22, c23], [c23, c33]]).det()
) == 0

ternary_rank_zero = sp.Poly((z * C3).det(), z)
assert sp.expand(
    ternary_rank_zero.coeff_monomial(z**3) - C3.det()
) == 0


# Binary quadratic bordered invariant in the isotropic case.
x, y = sp.symbols("x y")
alpha, beta, gamma, constant = sp.symbols(
    "alpha beta gamma constant"
)
binary_quadratic = alpha * x**2 + beta * x * y + gamma * y**2
binary_s = constant + x + binary_quadratic
binary_gradient = sp.Matrix(
    [sp.diff(binary_s, x), sp.diff(binary_s, y)]
)
binary_hessian = sp.hessian(binary_quadratic, (x, y))
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

# R=0 first forces gamma=0 and then beta=0.
assert binary_R.coeff_monomial(1) == 2 * gamma
assert sp.expand(
    binary_R.coeff_monomial(x).subs(gamma, 0) + 2 * beta**2
) == 0
assert sp.expand(
    binary_R.as_expr().subs({beta: 0, gamma: 0})
) == 0

# The surviving s=constant+x+alpha*x^2 has a critical point when alpha is
# nonzero.  At such a point the first row and column of the bordered Hessian
# vanish, so its determinant cannot be a nonzero constant.
critical_x = -sp.Rational(1, 2) / alpha
surviving_gradient = sp.Matrix(
    [
        sp.diff(
            binary_s.subs({beta: 0, gamma: 0}),
            variable,
        )
        for variable in (x, y)
    ]
)
assert sp.simplify(surviving_gradient.subs(x, critical_x)) == sp.zeros(2, 1)

q1, q2 = sp.symbols("q1 q2")
p11, p12, p22 = sp.symbols("p11 p12 p22")
bordered_hessian = sp.Matrix(
    [
        [0, q1, q2],
        [q1, p11, p12],
        [q2, p12, p22],
    ]
)
assert sp.expand(
    bordered_hessian.det()
    + (
        sp.Matrix([q1, q2]).T
        * sp.Matrix([[p11, p12], [p12, p22]]).adjugate()
        * sp.Matrix([q1, q2])
    )[0]
) == 0
assert bordered_hessian.det().subs({q1: 0, q2: 0}) == 0


print("PASS: the shared rank-three/rank-two/rank-one leading identities hold")
print("PASS: the rank-zero leading coefficient is det Hess(h_3)")
print("PASS: all three residual ternary-pencil rank identities hold")
print("PASS: the isotropic binary quadratic invariant forces beta=gamma=0")
print("PASS: nonzero alpha forces a zero bordered-Hessian row and determinant")
print("PASS: the shared cotangent-lift and terminal HC(2) blocks hold")
print(
    "SCOPE: with Gordan--Noether, HC(3), HC(2), and the plane degree "
    "bound, this excludes the full dense cubic--quartic chart"
)
