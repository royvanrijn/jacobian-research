#!/usr/bin/env python3
"""Exact J1 fiber certificate for the constant-normal R=1, v=Y branch."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, davenport_polynomial  # noqa: E402


z, U = sp.symbols("z U")
v = sp.Function("v")(T, z)
W = sp.Function("W")(T, z)
S = sp.Function("S")
L = sp.Function("L")

# Work in the function-field coordinates (T,z), where z=g_T(Y).  On the
# R=1 branch, use W=x+4v^2.  The quartic relation is W^2=D.
S_prime = sp.diff(S(z), z)
B = S_prime - 2 * (T + 1) * z
x = W - 4 * v**2
H = -T * z**2 - v**2 + W / 2 + S(z)
a2 = z * x + v
a3 = z**2 * x + 2 * z * v + 1

constant_normal_map = sp.Matrix(
    (
        T + x * U + U**2,
        z + a2 * U + z * U**2,
        H + a3 * U + z**2 * U**2,
    )
)
constant_normal_jacobian = sp.Poly(
    sp.expand(
        constant_normal_map.jacobian((T, z, U)).det()
    ),
    U,
)

discriminant = (
    16 * v**4
    + 8 * B * v
    + 4 * T
    - 2 * L(z)
)
square_relation = W**2 - discriminant
W_T = sp.solve(
    sp.diff(square_relation, T),
    sp.diff(W, T),
    dict=True,
)[0][sp.diff(W, T)]
W_z = sp.solve(
    sp.diff(square_relation, z),
    sp.diff(W, z),
    dict=True,
)[0][sp.diff(W, z)]

# The differentiated square relation is exactly the already integrated
# U^2 equation.
assert sp.simplify(
    constant_normal_jacobian.nth(2).subs(
        sp.diff(W, T),
        W_T,
    )
) == 0

# Any global solution must survive the ordinary fiber (T,z)=(1,0).
# In E=K[Y]/(g_1(Y)), implicit differentiation gives
# v_T=-g_T/g_Y and v_z=1/g_Y.
g = davenport_polynomial(Y)
g_T = sp.diff(g, T)
g_Y = sp.diff(g, Y)
fiber_polynomial = sp.Poly((7 * g).subs(T, 1), Y)
fiber_derivative = sp.Poly(g_Y.subs(T, 1), Y)
assert sp.gcd(fiber_polynomial, fiber_derivative).degree() == 0

s1, s2, ell0, ell1, w_symbol = sp.symbols(
    "s1 s2 ell0 ell1 w_symbol"
)
fiber_substitutions = {
    T: 1,
    z: 0,
    v: Y,
    W: w_symbol,
    sp.diff(v, T): -(g_T / g_Y).subs(T, 1),
    sp.diff(v, z): 1 / g_Y.subs(T, 1),
    S_prime: s1,
    sp.diff(S(z), (z, 2)): s2,
    L(z): ell0,
    sp.diff(L(z), z): ell1,
}

J1_with_square_derivatives = constant_normal_jacobian.nth(1).subs(
    {
        sp.diff(W, T): W_T,
        sp.diff(W, z): W_z,
    }
)
J1_fiber = sp.together(
    J1_with_square_derivatives.xreplace(
        fiber_substitutions
    )
)
J1_numerator = sp.fraction(J1_fiber)[0]

fiber_discriminant = (
    16 * Y**4
    + 8 * s1 * Y
    + 4
    - 2 * ell0
)
J1_linear = sp.rem(
    sp.Poly(J1_numerator, w_symbol, domain="EX"),
    sp.Poly(
        w_symbol**2 - fiber_discriminant,
        w_symbol,
        domain="EX",
    ),
).as_expr()
assert sp.Poly(J1_linear, w_symbol).degree() == 1

linear_coefficient = sp.diff(J1_linear, w_symbol)
linear_constant = J1_linear.subs(w_symbol, 0)
linear_coefficient = sp.Poly(
    linear_coefficient,
    Y,
    domain="EX",
).rem(sp.Poly(fiber_polynomial.as_expr(), Y, domain="EX")).as_expr()
linear_constant = sp.Poly(
    linear_constant,
    Y,
    domain="EX",
).rem(sp.Poly(fiber_polynomial.as_expr(), Y, domain="EX")).as_expr()

# Eliminate W from A*W+C=0 and W^2=D.  The resulting class
# C^2-A^2*D must vanish in the seven-dimensional Davenport algebra.
eliminated_square = sp.Poly(
    sp.expand(
        linear_constant**2
        - linear_coefficient**2 * fiber_discriminant
    ),
    Y,
    domain="EX",
).rem(sp.Poly(fiber_polynomial.as_expr(), Y, domain="EX"))
fiber_equations = tuple(
    eliminated_square.nth(degree)
    for degree in range(7)
)
assert len(fiber_equations) == 7

# Exact Gröbner certificate over Q(a), a^2+a+2=0.  The unit ideal proves
# that no four jets (S'(0),S''(0),L(0),L'(0)) satisfy this necessary fiber.
field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(
        field_variable**2 + field_variable + 2,
        field_variable,
    ),
    alias="root_a",
)
root_a = number_field.ext
jet_variables = (s1, s2, ell0, ell1)
number_field_equations = tuple(
    sp.Poly(
        equation.subs(A, root_a),
        *jet_variables,
        domain=number_field,
    ).as_expr()
    for equation in fiber_equations
)
groebner_certificate = sp.groebner(
    number_field_equations,
    *jet_variables,
    domain=number_field,
    order="grevlex",
    method="f5b",
)
assert len(groebner_certificate.polys) == 1
assert groebner_certificate.polys[0].as_expr() == 1

print("PASS: the R=1, v=Y square relation integrates J2 exactly")
print("PASS: the necessary (T,z)=(1,0) fiber is separable")
print("PASS: J1 is linear in the quartic square root on that fiber")
print("PASS: eliminating the square root gives seven jet equations")
print("PASS: their exact Groebner basis over Q(a) is the unit ideal")
print("PASS: the complete R=1, v=Y branch is impossible")
print("PASS Davenport constant-normal J1 fiber certificate")
