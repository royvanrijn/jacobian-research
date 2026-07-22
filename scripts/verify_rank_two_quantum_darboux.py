#!/usr/bin/env python3
"""Exact localized Darboux splitting behind the rank-two quotient.

The operators Z and E are first-order differential operators in (x,q).
Consequently their commutators are exactly hbar times the Lie brackets of
their coefficient vector fields: there are no hidden ordering corrections.
"""

import sympy as sp


x, q = sp.symbols("x q")
variables = (x, q)


def apply(vector_field, function):
    return sp.expand(
        sum(c * sp.diff(function, v) for c, v in zip(vector_field, variables))
    )


def lie(left, right):
    return tuple(
        sp.expand(apply(left, coefficient) - apply(right, left[index]))
        for index, coefficient in enumerate(right)
    )


# Z=3*x^2*p+(2-6*x*q)*zeta and
# E=(1+3*x*q)*p/2-3*q^2*zeta in left differential-operator order.
Z = (3 * x**2, 2 - 6 * x * q)
E = ((1 + 3 * x * q) / 2, -3 * q**2)

R = 2 * x - 3 * x**2 * q
v = 1 / x
P = tuple(-coefficient / 3 for coefficient in Z)
A = sp.factor((3 - R * v) * v**2 / 2)
U = tuple(sp.factor(E[index] + A * P[index]) for index in range(2))

assert apply(P, v) == 1
assert apply(P, R) == 0
assert apply(U, R) == 1
assert apply(U, v) == 0
assert lie(U, P) == (0, 0)

# The inverse position change on x!=0.
rho, inv_x = sp.symbols("rho inv_x")
q_inverse = inv_x * (2 - rho * inv_x) / 3
assert sp.factor(R.subs({x: 1 / inv_x, q: q_inverse}) - rho) == 0

# In the original adapted variables U=E-(1+3*X*Q)Z/(6*X^2), with the
# coefficient placed to the left of Z.
assert sp.factor(A - (1 + 3 * x * q) / (2 * x**2)) == 0

print("PASS: P=-Z/3 and v=X^-1 form an exact localized Weyl pair")
print("PASS: U=E+(3-Rv)v^2*P/2 is exactly conjugate to R")
print("PASS: the two localized Weyl pairs have all mixed commutators zero")
