#!/usr/bin/env python3
"""Exact Davenport-basis screen for the constant-normal conic branch."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, davenport_polynomial, reduce_a


z = sp.symbols("z")
c, B0, L0 = sp.symbols("c B0 L0", nonzero=True)

# The monic degree-seven equation for the Davenport algebra
# K(T,z)[Y]/(g_T(Y)-z).
f_general = sp.Poly(
    7 * (davenport_polynomial(Y) - z),
    Y,
)
assert f_general.degree() == 7
assert f_general.LC() == 1
assert f_general.nth(6) == 0

# For R=c != 0 and the natural primitive v=Y, the quartic-square gate is
#
#   W^2 = 16Y^4 + 8c B Y + C
#
# in the degree-seven Davenport algebra.  A general reduced W gives seven
# exact basis equations.
w_coefficients = sp.symbols("w0:7")
W_full = sum(
    coefficient * Y**degree
    for degree, coefficient in enumerate(w_coefficients)
)
quartic_gate = 16 * Y**4 + 8 * c * B0 * Y + L0
full_remainder = sp.Poly(
    W_full**2 - quartic_gate,
    Y,
    domain="EX",
).rem(sp.Poly(f_general.as_expr(), Y, domain="EX"))
full_basis_equations = tuple(
    sp.expand(full_remainder.nth(degree))
    for degree in range(7)
)
assert len(full_basis_equations) == 7
assert all(
    not equation.has(Y)
    for equation in full_basis_equations
)

# Specialize z=0.  Then B=B(0) is constant and
# C=4c^2*T-2c*L(0).  A reduced square root of degree at most three does
# not invoke the degree-seven relation.  Coefficient comparison forces
# W=+/-4Y^2, B(0)=0, and C=0.  The last identity is impossible because
# c and L(0) are constants while dC/dT=4c^2 != 0.
u0, u1, u2, u3 = sp.symbols("u0 u1 u2 u3")
W_low = u0 + u1 * Y + u2 * Y**2 + u3 * Y**3
C_zero_fiber = 4 * c**2 * T - 2 * c * L0
low_difference = sp.Poly(
    W_low**2
    - (
        16 * Y**4
        + 8 * c * B0 * Y
        + C_zero_fiber
    ),
    Y,
)
assert low_difference.nth(6) == u3**2
assert low_difference.nth(4).subs(u3, 0) == u2**2 - 16
assert low_difference.nth(3).subs(u3, 0) == 2 * u1 * u2
assert low_difference.nth(2).subs(
    {u3: 0, u1: 0}
) == 2 * u0 * u2
assert low_difference.nth(1).subs(
    {u3: 0, u1: 0, u0: 0}
) == -8 * B0 * c
assert sp.diff(C_zero_fiber, T) == 4 * c**2

# Degree four.  Write p=w4 != 0 and t=w3/w4.  The coefficients of
# Y^8,...,Y^4 determine the quotient h1*Y+h0 and w2,w1,w0.  The Y^3 and
# Y^2 equations then have a resultant in k=p^2.  Its only nonconstant
# factor in t is the irreducible quintic displayed in the note.
f_zero = sp.Poly(7 * davenport_polynomial(Y), Y)
f5, f4, f3, f2, f1 = (
    f_zero.nth(degree)
    for degree in (5, 4, 3, 2, 1)
)

p, t, k = sp.symbols("p t k", nonzero=True)
w4 = p
w3 = p * t
h1 = p**2
h0 = 2 * p**2 * t
w2 = sp.factor((h1 * f5 - w3**2) / (2 * w4))
w1 = sp.factor(
    (
        h1 * f4
        + h0 * f5
        - 2 * w3 * w2
    )
    / (2 * w4)
)
w0 = sp.factor(
    (
        h1 * f3
        + h0 * f4
        - 2 * w3 * w1
        - w2**2
        + 16
    )
    / (2 * w4)
)

equation_y3 = sp.factor(
    2 * w3 * w0
    + 2 * w2 * w1
    - h1 * f2
    - h0 * f3
).subs(p**2, k)
equation_y2 = sp.factor(
    2 * w2 * w0
    + w1**2
    - h1 * f1
    - h0 * f2
).subs(p**2, k)

resultant_k = sp.resultant(
    sp.together(equation_y3).as_numer_denom()[0],
    sp.together(equation_y2).as_numer_denom()[0],
    k,
)
resultant_k = sp.factor(reduce_a(resultant_k, T, t))

quintic_gate = (
    t**5
    + sp.Rational(1, 3) * t**4
    + T * (2 + sp.Rational(4, 3) * A) * t**3
    + T * (-4 - sp.Rational(2, 3) * A) * t**2
    + (
        T**2 * (-sp.Rational(23, 3) + A)
        + T * (1 + A / 3)
    )
    * t
    + T**2 * (7 + sp.Rational(7, 3) * A)
)
assert sp.expand(
    reduce_a(
        resultant_k
        + 2688 * T * (A + 1) * quintic_gate,
        T,
        t,
    )
) == 0

# Verify irreducibility over Q(a)(T), a^2+a+2=0.
field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(
        field_variable**2 + field_variable + 2,
        field_variable,
    ),
    alias="root_a",
)
root_a = number_field.ext
function_field = number_field.frac_field(T)
quintic_over_function_field = sp.Poly(
    quintic_gate.subs(A, root_a),
    t,
    domain=function_field,
)
unit, irreducible_factors = sp.factor_list(
    quintic_over_function_field
)
assert unit != 0
assert len(irreducible_factors) == 1
assert irreducible_factors[0][0].degree() == 5
assert irreducible_factors[0][1] == 1

print("PASS: the constant-normal v=Y gate has seven Davenport-basis equations")
print("PASS: reduced square roots of basis degree at most three are impossible")
print("PASS: basis degree four forces the recorded quintic equation")
print("PASS: that quintic is irreducible over Q(a)(T)")
print("PASS: any surviving v=Y square root uses Y^5 or Y^6")
print("PASS Davenport constant-normal basis screen")
