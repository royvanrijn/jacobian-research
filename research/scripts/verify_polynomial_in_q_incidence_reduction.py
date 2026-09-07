#!/usr/bin/env python3
"""Exact reduction of arbitrary polynomial Q-dependence.

On the fixed quadratic reciprocal chart, let B=Q+beta(P,S) and allow the
second target output to depend polynomially on Q to any degree.  Rewriting
it as C_hat(P,S,B) turns the Keller determinant equation into one ordinary
derivative at fixed (P,B).  Its general polynomial solution differs from
the standard S^2 marked line only by a triangular polynomial target shear.

The quadratic coefficient calculation is retained as a coefficientwise
regression of the general identity.
"""

from __future__ import annotations

import sympy as sp


P, S, Q, B, C_symbol, lam = sp.symbols(
    "P S Q B C lambda", nonzero=True
)
beta = sp.Function("beta")(P, S)
H = sp.Function("H")(P, B)

D = 1 - S * Q + P * S**2
B_expr = Q + beta
D_in_B = sp.expand(D.subs(Q, B - beta))
assert D_in_B == 1 - S * B + S * beta + P * S**2

# Passing from (P,S,Q) to (P,S,B) has determinant one.  In the latter
# coordinates a completely arbitrary polynomial incidence C_hat(P,S,B)
# satisfies det(P,B,C_hat)=-partial_S C_hat.
C_hat = sp.Function("C_hat")(P, S, B)
incidence_jacobian = sp.factor(
    sp.det(sp.Matrix([P, B, C_hat]).jacobian((P, S, B)))
)
assert incidence_jacobian == -sp.diff(C_hat, S)

# Integrate the required identity partial_S C_hat=lambda*D_in_B.  The
# integration "constant" is an arbitrary H(P,B), hence precisely a target
# shear.  SymPy leaves the beta integral unevaluated, which is the desired
# general formula.
Y0 = sp.integrate(lam * (1 + S * beta + P * S**2), S)
C_general = Y0 - lam * B * S**2 / 2 + H
assert sp.factor(sp.diff(C_general, S).doit() - lam * D_in_B) == 0

# Removing H is a triangular polynomial target automorphism whenever the
# original incidence is polynomial.  The remaining B coefficient is -X
# with X=lambda*S^2/2.
target_shear = sp.Matrix([P, B, C_symbol - H])
assert sp.factor(
    sp.det(target_shear.jacobian((P, B, C_symbol)))
) == 1
X_reduced = lam * S**2 / 2
assert sp.diff(X_reduced, S) == lam * S


# Coefficientwise quadratic regression.
c0 = sp.Function("c0")(P, S)
c1 = sp.Function("c1")(P, S)
c2 = sp.Function("c2")(P, S)

C_expr = c0 + c1 * Q + c2 * Q**2

jacobian = sp.expand(
    sp.det(sp.Matrix([P, B_expr, C_expr]).jacobian((P, S, Q)))
)
expected_jacobian = sp.expand(-lam * D)
defect = sp.Poly(sp.expand(jacobian - expected_jacobian), Q)

# The determinant identity is equivalent to these three coefficient
# equations.  In particular, its top coefficient forces c2_S=0.
q2_equation = sp.factor(defect.coeff_monomial(Q**2))
q1_equation = sp.factor(defect.coeff_monomial(Q))
q0_equation = sp.factor(defect.coeff_monomial(1))
assert q2_equation == -sp.diff(c2, S)
assert q1_equation == (
    -lam * S
    + 2 * sp.diff(beta, S) * c2
    - sp.diff(c1, S)
)
assert sp.factor(
    q0_equation
    - (
        lam * (P * S**2 + 1)
        + sp.diff(beta, S) * c1
        - sp.diff(c0, S)
    )
) == 0

# Work modulo q2_equation=0, so c2=g(P).  Integrating q1_equation=0 gives
# c1-2*g*beta=-lambda*S^2/2+f(P).
f = sp.Function("f")(P)
g = sp.Function("g")(P)
c1_solution = -lam * S**2 / 2 + 2 * g * beta + f
assert sp.factor(
    q1_equation.subs({c2: g, c1: c1_solution}).doit()
) == 0

# Rewrite C in the target coordinate B=Q+beta.  The B^2 coefficient is
# g(P), and after subtracting g(P)B^2 the B coefficient is
# c1-2*g*beta.  A further f(P)B shear leaves -lambda*S^2/2.
Q_in_B = B - beta
C_in_B = sp.expand(
    C_expr.subs({Q: Q_in_B, c2: g, c1: c1_solution})
)
C_sheared = sp.expand(C_in_B - g * B**2 - f * B)
C_sheared_poly = sp.Poly(C_sheared, B)
assert C_sheared_poly.degree() == 1
assert sp.factor(
    C_sheared_poly.coeff_monomial(B) + lam * S**2 / 2
) == 0

# Both shears are triangular polynomial target automorphisms when f and g
# are polynomials in P.  Their Jacobian is one.
quadratic_target_shear = sp.Matrix(
    [P, B, C_symbol - g * B**2 - f * B]
)
target_shear_jacobian = sp.factor(
    sp.det(quadratic_target_shear.jacobian((P, B, C_symbol)))
)
assert target_shear_jacobian == 1

print("PASS: arbitrary polynomial Q-dependence reduces at fixed (P,B)")
print("PASS: the integration freedom is exactly a polynomial target shear")
print("PASS: the all-degree reduced horizontal coordinate is lambda*S^2/2")
print("PASS: the quadratic coefficient equations recover the general theorem")
print("PASS: constant Jacobian forces the Q^2 coefficient to depend only on P")
print("PASS: coefficientwise target shears remove the apparent Q^2 dependence")
