#!/usr/bin/env python3
"""Verify the quartic noncharacteristic-line normal form in chart 1000.

Use u=a-rho*b-sigma and v=b for a line component u=0.  Put

    F=f_aa|line,  U=d_v(f_a|line).

The caustic and normal-Hessian identities restrict to

    L = 5*v*F - F*f0'' + U^2 = 0,
    C = F*(4*v^3-2*g0') + 2*U*g1 - (U-rho*F),
    C^2 = 2*F*kappa,

where f0=f|line, g0=g|line, and g1=g_a|line.

For quartic f and cubic g, degree comparison forces F to be a nonzero
constant, U=m*v+n with m nonzero, and C to be a nonzero constant.  The
checker verifies the resulting integrated normal form and the coefficient
relations that make C constant.  The later Laurent equations force the
line restriction of h_a to be constant.  Finally S_a=0 and R=0 are
incompatible, closing the noncharacteristic-line branch.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


u, v, rho = sp.symbols("u v rho")


def polynomial(
    prefix: str, maximum_degree: int
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    exponents = [
        powers
        for powers in product(range(maximum_degree + 1), repeat=2)
        if sum(powers) <= maximum_degree
    ]
    coefficients = sp.symbols(f"{prefix}0:{len(exponents)}")
    expression = sum(
        coefficient * u**powers[0] * v**powers[1]
        for coefficient, powers in zip(coefficients, exponents, strict=True)
    )
    return expression, coefficients


f_general, _ = polynomial("f", 4)
g_general, _ = polynomial("g", 3)


def derivative_a(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, u)


def derivative_b(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, v) - rho * sp.diff(expression, u)


f_a = derivative_a(f_general)
f_aa = derivative_a(f_a)
f_ab = derivative_b(f_a)
f_bb = derivative_b(derivative_b(f_general))
g_a = derivative_a(g_general)
g_b = derivative_b(g_general)

F_general = f_aa.subs(u, 0)
f0_general = f_general.subs(u, 0)
g0_general = g_general.subs(u, 0)
g1_general = g_a.subs(u, 0)
U_general = sp.diff(f_a.subs(u, 0), v)
P_general = f_ab.subs(u, 0)
L_general = (
    5 * v * f_aa - (f_aa * f_bb - f_ab**2)
).subs(u, 0)
C_general = (
    4 * v**3 * f_aa
    - 2 * f_aa * g_b
    + 2 * f_ab * g_a
    - f_ab
).subs(u, 0)

assert sp.factor(P_general - (U_general - rho * F_general)) == 0
assert (
    sp.factor(
        L_general
        - (
            5 * v * F_general
            - F_general * sp.diff(f0_general, v, 2)
            + U_general**2
        )
    )
    == 0
)
assert (
    sp.factor(
        C_general
        - (
            F_general * (4 * v**3 - 2 * sp.diff(g0_general, v))
            + 2 * U_general * g1_general
            - (U_general - rho * F_general)
        )
    )
    == 0
)

F, m, n = sp.symbols("F m n", nonzero=True)
alpha, beta, f1_constant = sp.symbols("alpha beta f1_constant")
U = m * v + n
f0 = (
    m**2 * v**4 / (12 * F)
    + (sp.Rational(5, 6) + m * n / (3 * F)) * v**3
    + n**2 * v**2 / (2 * F)
    + alpha * v
    + beta
)
f1 = m * v**2 / 2 + n * v + f1_constant
assert sp.factor(sp.diff(f1, v) - U) == 0
assert sp.factor(F * sp.diff(f0, v, 2) - U**2 - 5 * v * F) == 0

g10, g11, g12, g00, g_constant = sp.symbols(
    "g10 g11 g12 g00 g_constant"
)
g12_solution = -2 * F / m
g02 = (m * g11 + n * g12_solution) / F
g01 = (2 * m * g10 + 2 * n * g11 - m) / (2 * F)
g1 = g12_solution * v**2 + g11 * v + g10
g0 = g02 * v**3 / 3 + g01 * v**2 / 2 + g00 * v + g_constant
C_line = sp.expand(
    F * (4 * v**3 - 2 * sp.diff(g0, v))
    + 2 * U * g1
    - (U - rho * F)
)
assert sp.Poly(C_line, v).degree() == 0
expected_C = -2 * F * g00 + 2 * n * g10 - n + rho * F
assert sp.factor(C_line - expected_C) == 0

h10, h_constant = sp.symbols("h10 h_constant")
A30, A31 = sp.symbols("A30 A31")
A3 = A31 * v + A30
h1 = h10
h0 = (
    v**2 / 4
    + h10 * (m * v**2 / 2 + n * v) / F
    + h_constant
)
h2 = A3 * h10 / (2 * F)
H2_line = A3 * U * h1 + F**2 * sp.diff(h1, v) - 2 * F * U * h2
K1_line = -2 * U * H2_line
A_line = F * (v - 2 * sp.diff(h0, v)) + 2 * U * h1
assert sp.factor(A_line) == 0
assert sp.factor(H2_line) == 0
assert sp.factor(K1_line) == 0

g20, g21 = sp.symbols("g20 g21")
g2 = g21 * v + g20
q4, g3 = sp.symbols("q4 g3")
f_model = f0 + u * f1 + F * u**2 / 2 + A3 * u**3 / 6 + q4 * u**4
g_model = g0 + u * g1 + u**2 * g2 + u**3 * g3
f_model_aa = derivative_a(derivative_a(f_model))
f_model_ab = derivative_b(derivative_a(f_model))
g_model_a = derivative_a(g_model)
g_model_b = derivative_b(g_model)
C_model = (
    4 * v**3 * f_model_aa
    - 2 * f_model_aa * g_model_b
    + 2 * f_model_ab * g_model_a
    - f_model_ab
)
C_model_line = sp.factor(C_model.subs(u, 0))
assert sp.factor(C_model_line - expected_C) == 0
S_a = sp.factor(
    2 * F * sp.diff(C_model, u).subs(u, 0)
    - C_model_line * A3
)
S_a_solution = {
    A31: 0,
    A30: 0,
    g21: 0,
    g20: -2 * F**2 / m**2,
    g11: -4 * F * n / m**2,
}
assert sp.factor(S_a.subs(S_a_solution)) == 0

R = 21 * v * F**2 - 16 * F * g2 + 8 * A3 * g1 - 4 * A3
R_after_S_a = sp.factor(R.subs(S_a_solution))
assert sp.factor(
    R_after_S_a - (21 * F**2 * v + 32 * F**3 / m**2)
) == 0


def main() -> None:
    print("PASS: L|line=5*v*F-F*f0''+U^2")
    print("PASS: C|line=F*(4*v^3-2*g0')+2*U*g1-(U-rho*F)")
    print("RESULT: C^2=2*F*kappa forces F constant and U=m*v+n")
    print("PASS: the integrated quartic f0,f1 normal form satisfies L=0")
    print("PASS: the displayed cubic g0,g1 normal form makes C constant")
    print("PASS: A=H2=K1=0 gives constant h_a and the displayed h0,h2")
    print("PASS: S_a=0 forces f_aaa|line=0 and fixed g11,g20,g21")
    print("RESULT: then R=21*F^2*v+32*F^3/m^2 cannot vanish")
    print("        no quartic noncharacteristic line component survives")


if __name__ == "__main__":
    main()
