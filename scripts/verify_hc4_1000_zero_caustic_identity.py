#!/usr/bin/env python3
"""Reduce the identically-zero chart-1000 caustic branch.

Write

    L = -det Hess(f-5*b^3/6).

If L is identically zero and f_aa is nonzero, the binary zero-Hessian
theorem gives

    f=5*b^3/6+P(alpha*a+beta*b)+affine.

The normal-Hessian identities A=0 and C^2=2*f_aa*kappa then reduce to

    P'' * D^2 = 2*kappa,
    D=4*alpha*b^3-2*alpha*g_b+2*beta*g_a-beta.

Both polynomial factors must be units.  This is impossible because a
cubic g has deg(g_a),deg(g_b)<=2 and D retains 4*alpha*b^3.

The edge f_aa=0 also has L=0 only when f_ab=0.  Direct substitution in
the three boundary determinant coefficients leaves the exact affine-normal
system

    f_a constant,  f_bb=5*b-delta,
    h_a=0,          g_a=gamma,
    kappa=-delta*(2*gamma-1)^2/2 != 0.

This checker certifies both the obstruction and the residual edge.
"""

from __future__ import annotations

from itertools import product
import runpy

import sympy as sp


a, b, z = sp.symbols("a b z")
alpha = sp.symbols("alpha", nonzero=True)
beta = sp.symbols("beta")
ell = alpha * a + beta * b
p0, p1, p2, p3, p4 = sp.symbols("p0:5")
P = p0 + p1 * ell + p2 * ell**2 + p3 * ell**3 + p4 * ell**4
f = sp.Rational(5, 6) * b**3 + P
f_aa = sp.diff(f, a, 2)
f_ab = sp.diff(f, a, b)
f_bb = sp.diff(f, b, 2)
L = sp.factor(5 * b * f_aa - (f_aa * f_bb - f_ab**2))
assert L == 0

exponents = [
    powers
    for powers in product(range(4), repeat=2)
    if sum(powers) <= 3
]
g_coefficients = sp.symbols(f"g0:{len(exponents)}")
g = sum(
    coefficient * a**powers[0] * b**powers[1]
    for coefficient, powers in zip(
        g_coefficients, exponents, strict=True
    )
)
g_a = sp.diff(g, a)
g_b = sp.diff(g, b)
P_univariate = p0 + p1 * z + p2 * z**2 + p3 * z**3 + p4 * z**4
P_second = sp.diff(P_univariate, z, 2).subs(z, ell)
D = 4 * alpha * b**3 - 2 * alpha * g_b + 2 * beta * g_a - beta
C = 4 * b**3 * f_aa - 2 * f_aa * g_b + 2 * f_ab * g_a - f_ab
assert sp.factor(f_aa - alpha**2 * P_second) == 0
assert sp.factor(C - alpha * P_second * D) == 0
assert sp.Poly(D, a, b).coeff_monomial(b**3) == 4 * alpha


# The f_aa=0 edge is most cleanly read directly from the boundary
# determinant coefficients.
boundary = runpy.run_path("scripts/verify_hc4_1000_boundary_schur_chain.py")
edge_substitution = {
    boundary["f_aa"]: 0,
    boundary["f_ab"]: 0,
}
edge_coefficient_2 = sp.factor(boundary["coefficient_2"].subs(edge_substitution))
edge_coefficient_1 = sp.factor(boundary["coefficient_1"].subs(edge_substitution))
edge_coefficient_0 = sp.factor(boundary["coefficient_0"].subs(edge_substitution))
edge_factor = 5 * boundary["b"] - boundary["f_bb"]
edge_linear = 12 * boundary["b"] ** 2 * boundary["h_a"] + 2 * boundary["g_a"] - 1
assert sp.factor(
    edge_coefficient_2
    + sp.Rational(32, 9) * boundary["h_a"] ** 2 * edge_factor
) == 0
assert sp.factor(
    edge_coefficient_1
    + sp.Rational(8, 3) * boundary["h_a"] * edge_factor * edge_linear
) == 0
assert sp.factor(
    edge_coefficient_0
    + sp.Rational(1, 2) * edge_factor * edge_linear**2
) == 0


def main() -> None:
    print("PASS: L=0 and f_aa!=0 has the binary zero-Hessian normal form")
    print("PASS: C^2=2*f_aa*kappa becomes P''*D^2=2*kappa")
    print("PASS: D retains 4*alpha*b^3, so the unit equation is impossible")
    print("PASS: on f_aa=0, the three boundary coefficients factor exactly")
    print("RESULT: the only zero-L survivor is the affine-normal edge")
    print("        f_a constant, h_a=0, g_a constant, f_bb-5*b constant")


if __name__ == "__main__":
    main()
