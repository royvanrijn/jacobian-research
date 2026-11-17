#!/usr/bin/env python3
"""Verify the characteristic-line normal form in chart 1000.

Let p=a-rho*b-sigma and use u=p, v=b.  Along the reduced line u=0,
assume f_aa is generically nonzero and the line is characteristic for the
Hamiltonian field of f_a.  Equivalently, f_a is constant on the line.
Then

    f = c*u + phi(v) + u^2*q(u,v).

The chart-1000 caustic polynomial

    L = 5*b*f_aa - (f_aa*f_bb-f_ab^2)

restricts to

    L|_{u=0} = f_uu(0,v) * (5*v-phi''(v)).

Thus a characteristic line component of L=0 forces

    phi(v) = 5*v^3/6 + alpha*v + beta.

The next divisor-local Laurent conditions also collapse.  If

    h0(v)=h(0,v),  h1(v)=h_u(0,v),  F(v)=f_uu(0,v),

then on the line

    A  = F*(v-2*h0'),
    H1 = F^3*(2*h0''-1),
    H2 = F*(F*h1'-F'*h1).

The normal-Hessian gate A=0 therefore gives h0=v^2/4+delta.  The later
condition h1*H2=0 then gives h1=0 or h1/F constant.

For a quartic potential, g0(v)=g(0,v) has degree at most three and

    C|line = F*(4*v^3-2*g0'(v)+rho).

The parenthesis has degree exactly three.  The other normal-Hessian gate,
C^2=2*F*kappa, would make F times its square a nonzero constant.  This is
impossible.  Thus quartic characteristic line components do not occur.

For a quartic boundary potential, q has degree at most two.  This is a
complete obstruction to the characteristic-line branch.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


u, v, rho = sp.symbols("u v rho")
c = sp.symbols("c")
phi_coefficients = sp.symbols("phi0:5")
phi = sum(
    coefficient * v**degree
    for degree, coefficient in enumerate(phi_coefficients)
)

q_coefficients = sp.symbols("q0:6")
q_exponents = [
    exponents
    for exponents in product(range(3), repeat=2)
    if sum(exponents) <= 2
]
q = sum(
    coefficient * u**exponents[0] * v**exponents[1]
    for coefficient, exponents in zip(
        q_coefficients, q_exponents, strict=True
    )
)
f = c * u + phi + u**2 * q


def derivative_a(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, u)


def derivative_b(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, v) - rho * sp.diff(expression, u)


f_a = derivative_a(f)
f_aa = derivative_a(f_a)
f_ab = derivative_b(f_a)
f_bb = derivative_b(derivative_b(f))
L = 5 * v * f_aa - (f_aa * f_bb - f_ab**2)

assert sp.factor(f_a.subs(u, 0) - c) == 0
assert sp.factor(f_ab.subs(u, 0) + rho * f_aa.subs(u, 0)) == 0
assert (
    sp.factor(
        L.subs(u, 0)
        - f_aa.subs(u, 0) * (5 * v - sp.diff(phi, v, 2))
    )
    == 0
)

caustic_equations = sp.Poly(
    sp.expand(5 * v - sp.diff(phi, v, 2)), v
).all_coeffs()
phi_solution = {
    phi_coefficients[4]: 0,
    phi_coefficients[3]: sp.Rational(5, 6),
    phi_coefficients[2]: 0,
}
assert all(sp.factor(equation.subs(phi_solution)) == 0 for equation in caustic_equations)
assert sp.factor(
    phi.subs(phi_solution)
    - (
        sp.Rational(5, 6) * v**3
        + phi_coefficients[1] * v
        + phi_coefficients[0]
    )
) == 0

h_coefficients = sp.symbols("h0:10")
h_exponents = [
    exponents
    for exponents in product(range(4), repeat=2)
    if sum(exponents) <= 3
]
h = sum(
    coefficient * u**exponents[0] * v**exponents[1]
    for coefficient, exponents in zip(
        h_coefficients, h_exponents, strict=True
    )
)
h_a = derivative_a(h)
h_b = derivative_b(h)
h_aa = derivative_a(h_a)
h_ab = derivative_b(h_a)
h_bb = derivative_b(h_b)
f_aaa = derivative_a(f_aa)
f_aab = derivative_b(f_aa)
f_abb = derivative_b(f_ab)

A = f_aa * v - 2 * f_aa * h_b + 2 * f_ab * h_a
H1 = (
    2 * f_aa**3 * h_bb
    - f_aa**3
    - 4 * f_aa**2 * f_ab * h_ab
    - 4 * f_aa**2 * f_abb * h_a
    + 8 * f_aa * f_aab * f_ab * h_a
    + 2 * f_aa * f_ab**2 * h_aa
    - 4 * f_aaa * f_ab**2 * h_a
)
H2 = (
    f_aa**2 * h_ab
    - f_aa * f_aab * h_a
    - f_aa * f_ab * h_aa
    + f_aaa * f_ab * h_a
)
h0 = h.subs(u, 0)
h1 = h_a.subs(u, 0)
F = f_aa.subs(u, 0)
assert sp.factor(A.subs(u, 0) - F * (v - 2 * sp.diff(h0, v))) == 0
assert (
    sp.factor(H1.subs(u, 0) - F**3 * (2 * sp.diff(h0, v, 2) - 1))
    == 0
)
assert (
    sp.factor(
        H2.subs(u, 0)
        - F * (F * sp.diff(h1, v) - sp.diff(F, v) * h1)
    )
    == 0
)

g_coefficients = sp.symbols("g0:10")
g = sum(
    coefficient * u**exponents[0] * v**exponents[1]
    for coefficient, exponents in zip(
        g_coefficients, h_exponents, strict=True
    )
)
g_a = derivative_a(g)
g_b = derivative_b(g)
C = (
    4 * v**3 * f_aa
    - 2 * f_aa * g_b
    + 2 * f_ab * g_a
    - f_ab
)
g0 = g.subs(u, 0)
D = 4 * v**3 - 2 * sp.diff(g0, v) + rho
assert sp.factor(C.subs(u, 0) - F * D) == 0
assert sp.Poly(D, v).degree() == 3
assert sp.Poly(D, v).LC() == 4


def main() -> None:
    print("PASS: characteristicity gives f_a|line=c and f_ab=-rho*f_aa")
    print("PASS: L|line=f_aa*(5*v-phi'')")
    print("RESULT: a characteristic line forces")
    print("        f|line=5*v^3/6+alpha*v+beta")
    print("PASS: A, H1, H2 reduce to three univariate line expressions")
    print("RESULT: A=0 forces h|line=v^2/4+delta")
    print("        and h_a*H2=0 forces h_a|line=lambda*f_aa|line")
    print("PASS: C|line=F*(4*v^3-2*g0'+rho), with cubic second factor")
    print("RESULT: C^2=2*F*kappa is impossible on a quartic line")
    print("        no quartic characteristic line component survives")


if __name__ == "__main__":
    main()
