#!/usr/bin/env python3
"""Exclude the chart-1000 unit-caustic branch for quartic potentials.

Put phi=f-5*b^3/6.  If L is a nonzero constant, then

    det Hess(phi)=-L in k^*.

The binary homogeneous equations through degree four put phi, after a
linear shear preserving b, in one of the triangular forms

    phi=q*a*b+P(a),  or  phi=q*a*b+P(b),       q != 0.

The checker verifies the quartic leading-form calculation and then uses
only degree bounds on the already-forced normal Hessian V_cc.

In the P(a) form, writing

    P=p4*a^4+p3*a^3+p2*a^2+affine,

the coefficients [a^2*b^6], [a*b^6], [b^6] of V_cc are

    -48*p4/q^2, -24*p3/q^2, -8*p2/q^2.

Thus P is affine and the branch moves to the P(b) form.

Write P''=r2*b^2+r1*b+r0 and let g be an arbitrary cubic.  If r2!=0,
the b^6 and b^5 coefficients give g_ab^2=0 and then 99/80=0.  Hence
r2=0.  The remaining r1!=0 and r1=0 cases each have a seven-coefficient
triangular contradiction, ending in 99*q^2=0.  Therefore no unit-L
quartic boundary data survive.
"""

from __future__ import annotations

from itertools import product
import runpy

import sympy as sp


a, b = sp.symbols("a b")

# Exact quartic triangular classification in the chart where the top
# fourth power has nonzero a coefficient.
lam = sp.symbols("lam", nonzero=True)
c30, c21, c12, c03 = sp.symbols("c30 c21 c12 c03")
q20, q11, q02 = sp.symbols("q20 q11 q02")
phi = (
    lam * a**4
    + c30 * a**3
    + c21 * a**2 * b
    + c12 * a * b**2
    + c03 * b**3
    + q20 * a**2
    + q11 * a * b
    + q02 * b**2
)
binary_determinant = sp.Poly(
    sp.expand(sp.det(sp.hessian(phi, (a, b)))), a, b
)
assert binary_determinant.coeff_monomial(a**3) == 24 * c12 * lam
assert binary_determinant.coeff_monomial(a**2 * b) == 72 * c03 * lam
quartic_reduction = {c12: 0, c03: 0}
assert sp.factor(
    binary_determinant.coeff_monomial(a**2).subs(quartic_reduction)
    - 4 * (-c21**2 + 6 * lam * q02)
) == 0
assert sp.factor(
    binary_determinant.coeff_monomial(b).subs(quartic_reduction)
    - 4 * c21 * q02
) == 0

# If the quartic part vanishes, the identical cubic calculation gives the
# same triangular form.
lam3 = sp.symbols("lam3", nonzero=True)
phi_cubic = lam3 * a**3 + q20 * a**2 + q11 * a * b + q02 * b**2
cubic_determinant = sp.Poly(
    sp.expand(sp.det(sp.hessian(phi_cubic, (a, b)))), a, b
)
assert cubic_determinant.coeff_monomial(a) == 12 * lam3 * q02
assert cubic_determinant.coeff_monomial(1) == 4 * q20 * q02 - q11**2


boundary = runpy.run_path("scripts/verify_hc4_1000_boundary_schur_chain.py")
q, kappa = sp.symbols("q kappa", nonzero=True)
p4, p3, p2 = sp.symbols("p4 p3 p2")

g_exponents = [
    powers
    for powers in product(range(4), repeat=2)
    if sum(powers) <= 3
]
g_coefficients = sp.symbols(f"g0:{len(g_exponents)}")
g = sum(
    coefficient * a**powers[0] * b**powers[1]
    for coefficient, powers in zip(
        g_coefficients, g_exponents, strict=True
    )
)
g_a = sp.diff(g, a)
g_b = sp.diff(g, b)


def forced_V_cc(f_aa: sp.Expr, f_bb: sp.Expr) -> sp.Poly:
    expression = sp.cancel(
        boundary["forced_V_cc"].subs(
            {
                boundary["f_aa"]: f_aa,
                boundary["f_ab"]: q,
                boundary["f_bb"]: f_bb,
                boundary["g_a"]: g_a,
                boundary["g_b"]: g_b,
                boundary["kappa"]: kappa,
            }
        )
    )
    return sp.Poly(sp.expand(expression), a, b)


# P(a) orientation.
f_aa_a = 12 * p4 * a**2 + 6 * p3 * a + 2 * p2
V_cc_a = forced_V_cc(f_aa_a, 5 * b)
assert V_cc_a.coeff_monomial(a**2 * b**6) == -48 * p4 / q**2
assert V_cc_a.coeff_monomial(a * b**6) == -24 * p3 / q**2
assert V_cc_a.coeff_monomial(b**6) == -8 * p2 / q**2


# P(b) orientation.
r2, r1, r0 = sp.symbols("r2 r1 r0")
V_cc_b = forced_V_cc(0, 5 * b + r2 * b**2 + r1 * b + r0)
g6, g7, g8, g9 = (
    g_coefficients[index] for index in (6, 7, 8, 9)
)
assert V_cc_b.coeff_monomial(b**6) == -g6**2 * r2 / q**2
assert sp.factor(
    V_cc_b.coeff_monomial(b**5)
    - (
        -160 * g_coefficients[5] * g6 * r2
        - 80 * g6**2 * r1
        - 320 * g6 * q
        + 99 * q**2
    )
    / (80 * q**2)
) == 0

# With r2=0, record the seven coefficients used in both r1 cases.
linear_b_substitution = {r2: 0}
coeff_41 = sp.factor(
    V_cc_b.coeff_monomial(a**4 * b).subs(linear_b_substitution)
)
coeff_31 = sp.factor(
    V_cc_b.coeff_monomial(a**3 * b).subs(linear_b_substitution)
)
coeff_23 = sp.factor(
    V_cc_b.coeff_monomial(a**2 * b**3).subs(linear_b_substitution)
)
coeff_21 = sp.factor(
    V_cc_b.coeff_monomial(a**2 * b).subs(linear_b_substitution)
)
coeff_13 = sp.factor(
    V_cc_b.coeff_monomial(a * b**3).subs(linear_b_substitution)
)
coeff_05 = sp.factor(
    V_cc_b.coeff_monomial(b**5).subs(linear_b_substitution)
)
assert coeff_41 == -9 * g9**2 * r1 / q**2
assert sp.factor(coeff_31.subs(g9, 0) - 4 * g8**2 / q) == 0
assert sp.factor(coeff_23.subs(r1, 0) + 12 * g9 / q) == 0
assert sp.factor(
    coeff_21.subs({g9: 0, g8: 0})
    - 4 * g7 * (2 * g6 * q - g7 * r1) / q**2
) == 0
assert sp.factor(
    coeff_13.subs({g9: 0, g8: 0})
    - 4 * (g6**2 * q - g6 * g7 * r1 - 2 * g7 * q) / q**2
) == 0
assert sp.factor(
    coeff_05
    - (-80 * g6**2 * r1 - 320 * g6 * q + 99 * q**2)
    / (80 * q**2)
) == 0


def main() -> None:
    print("PASS: constant binary Hessian data have two triangular orientations")
    print("PASS: the P(a) orientation forces P''=0 from three V_cc terms")
    print("PASS: in the P(b) orientation, r2!=0 gives 99/80=0")
    print("PASS: the r1!=0 coefficient chain ends in 99*q^2=0")
    print("PASS: the r1=0 coefficient chain also ends in 99*q^2=0")
    print("RESULT: no quartic unit-L boundary data survive")


if __name__ == "__main__":
    main()
