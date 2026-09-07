#!/usr/bin/env python3
"""Close the asymptotic quartic characteristic conics in chart 1000.

These are the smooth conics with p_aa=0, omitted by the transverse-conic
parity argument.

Parabolic case.  After a shear preserving b,

    p=a+q(b),  deg(q)=2.

Divisor parity in k[b] forces the characteristic multiplier in
f_a-c=p*ell to be a nonzero constant z.  Caustic divisibility gives

    f=c*a+z*p^2/2+5*b^3/6+affine.

The equation C^2=2*f_aa*kappa forces the factor D=C/z to be constant.
For an arbitrary cubic g, this determines five coefficients and in
particular makes g_aa constant.  The h_a != 0 branch of the next Laurent
chain would require

    R=21*b*z^2-8*z*g_aa=0,

which is impossible.  On the h_a=0 branch, the already-forced normal
Hessian has V_cd=5*b^3/16.  A quartic potential has deg(V_cd)<=2, so this
branch is impossible as well.

Hyperbolic case.  After a shear and scaling,

    p=u*v-1,  d_a=d_u,  d_b=d_v-rho*d_u.

In k[v,v^-1], the only affine-linear multipliers for which p_a/ell can be
a square are ell=x*u and ell=y*v.  The first leaves the nonzero residue
x^2/v^4 in L|p after all available phi coefficients are fixed.  In the
second case the normal equation requires D to be constant, but every
cubic g leaves the unavoidable Laurent term 4*v^4 in D|p.

Thus no asymptotic quartic characteristic conic survives.
"""

from __future__ import annotations

from itertools import product
import runpy

import sympy as sp


def cubic(
    first: sp.Symbol, second: sp.Symbol, prefix: str
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    exponents = [
        powers
        for powers in product(range(4), repeat=2)
        if sum(powers) <= 3
    ]
    coefficients = sp.symbols(f"{prefix}0:{len(exponents)}")
    expression = sum(
        coefficient * first**powers[0] * second**powers[1]
        for coefficient, powers in zip(
            coefficients, exponents, strict=True
        )
    )
    return expression, coefficients


# Parabolic case.
a, b = sp.symbols("a b")
q0, q1, q2 = sp.symbols("q0 q1 q2", nonzero=True)
z, c = sp.symbols("z c", nonzero=True)
q = q2 * b**2 + q1 * b + q0
p_parabolic = a + q
f_parabolic = c * a + z * p_parabolic**2 / 2 + sp.Rational(5, 6) * b**3

f_aa = sp.diff(f_parabolic, a, 2)
f_ab = sp.diff(f_parabolic, a, b)
f_bb = sp.diff(f_parabolic, b, 2)
L_parabolic = 5 * b * f_aa - (f_aa * f_bb - f_ab**2)
assert sp.factor(L_parabolic + 2 * z**2 * q2 * p_parabolic) == 0

g_parabolic, gp = cubic(a, b, "gp")
g0, g1, _, _, g4, g5, _, g7, _, _ = gp
g2_solution = (
    g4 * q2
    + g5 * q1
    - g7 * (q0 * q2 + q1**2)
    - (q0 + q2**2) / (2 * q2)
)
g3_solution = g5 * q2 - g7 * q1 * q2 - q1 / (2 * q2)
g6_solution = g7 * q2 - 1 / (2 * q2)
g_parabolic = g_parabolic.subs(
    {
        gp[2]: g2_solution,
        gp[3]: g3_solution,
        gp[6]: g6_solution,
        gp[8]: 0,
        gp[9]: 0,
    }
)
g_a = sp.diff(g_parabolic, a)
g_b = sp.diff(g_parabolic, b)
p_b = sp.diff(p_parabolic, b)
D_parabolic = 4 * b**3 - 2 * g_b + p_b * (2 * g_a - 1)
D_constant = sp.factor(D_parabolic.subs(a, -q))
assert not D_constant.has(a, b)
assert sp.factor(sp.diff(g_parabolic, a, 2) - 2 * g7) == 0

# If h_a is nonzero on the conic, h_a*H2=0 makes it constant there.
# K2 and K3 force both derivatives of C^2/f_aa to vanish, and K4 then
# forces R=0.  Its b coefficient is nonzero.
R_parabolic = 21 * b * z**2 - 8 * z * sp.diff(g_parabolic, a, 2)
assert sp.Poly(R_parabolic, b).coeff_monomial(b) == 21 * z**2

# On the h_a=0 branch, h=b^2/4 modulo an irrelevant constant.  Verify the
# exact forced normal-Hessian entry using the canonical Schur certificate.
boundary = runpy.run_path("scripts/verify_hc4_1000_boundary_schur_chain.py")
C_constant = z * D_constant
kappa_parabolic = C_constant**2 / (2 * z)
h_parabolic = b**2 / 4
boundary_substitution = {
    boundary["f_aa"]: f_aa,
    boundary["f_ab"]: f_ab,
    boundary["f_bb"]: f_bb,
    boundary["g_a"]: g_a,
    boundary["g_b"]: g_b,
    boundary["h_a"]: sp.diff(h_parabolic, a),
    boundary["h_b"]: sp.diff(h_parabolic, b),
    boundary["kappa"]: kappa_parabolic,
}
forced_V_cd = sp.factor(
    sp.cancel(boundary["forced_V_cd"].subs(boundary_substitution))
)
assert sp.factor(forced_V_cd - sp.Rational(5, 16) * b**3) == 0


# Hyperbolic case.
u, v, rho = sp.symbols("u v rho")
x = sp.symbols("x", nonzero=True)
p_hyperbolic = u * v - 1


def derivative_a(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, u)


def derivative_b(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, v) - rho * sp.diff(expression, u)


def restrict_hyperbola(expression: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(expression.subs(u, 1 / v)))


phi_affine = sp.symbols("phi0") + sp.symbols("phi1") * v

# ell=x*u.  Caustic coefficient comparison fixes the quadratic-through-
# quartic part of phi, but a nonzero x^2/v^4 residue remains.
f_u_multiplier = (
    c * u
    + x * u**3 * v / 3
    - x * u**2 / 2
    + sp.Rational(5, 6) * v**3
    + phi_affine
)
F_u = derivative_a(derivative_a(f_u_multiplier))
P_u = derivative_b(derivative_a(f_u_multiplier))
fbb_u = derivative_b(derivative_b(f_u_multiplier))
L_u = 5 * v * F_u - (F_u * fbb_u - P_u**2)
assert sp.factor(restrict_hyperbola(L_u) - x**2 / v**4) == 0

# ell=x*v.  This branch satisfies the caustic equation, but D cannot be
# constant because its v^4 coefficient is always 4.
f_v_multiplier = (
    c * u
    + x * u**2 * v**2 / 2
    - x * u * v
    + sp.Rational(5, 6) * v**3
    + phi_affine
)
F_v = derivative_a(derivative_a(f_v_multiplier))
P_v = derivative_b(derivative_a(f_v_multiplier))
fbb_v = derivative_b(derivative_b(f_v_multiplier))
L_v = 5 * v * F_v - (F_v * fbb_v - P_v**2)
assert restrict_hyperbola(L_v) == 0

g_hyperbolic, _ = cubic(u, v, "gh")
p_a_hyperbolic = derivative_a(p_hyperbolic)
p_b_hyperbolic = derivative_b(p_hyperbolic)
D_hyperbolic = (
    p_a_hyperbolic * (4 * v**3 - 2 * derivative_b(g_hyperbolic))
    + p_b_hyperbolic * (2 * derivative_a(g_hyperbolic) - 1)
)
D_hyperbolic_restriction = sp.expand(
    restrict_hyperbola(D_hyperbolic)
)
assert D_hyperbolic_restriction.coeff(v, 4) == 4


def main() -> None:
    print("PASS: parabolic characteristicity gives f=c*a+z*p^2/2+5*b^3/6")
    print("PASS: D constant makes g_aa constant")
    print("PASS: the h_a!=0 Laurent branch would require an impossible R=0")
    print("PASS: the h_a=0 branch forces V_cd=5*b^3/16")
    print("PASS: the hyperbolic ell=x*u branch leaves L|p=x^2/v^4")
    print("PASS: the hyperbolic ell=x*v branch leaves 4*v^4 in D|p")
    print("RESULT: no asymptotic quartic characteristic conic survives")


if __name__ == "__main__":
    main()
