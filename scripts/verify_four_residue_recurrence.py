#!/usr/bin/env python3
"""Exact recurrence for pole-order-four quotient Hamiltonian residues."""

import sympy as sp


X, Q, C = sp.symbols("X Q C")

# The normalized top-pole solution.  For a Z-independent Laurent principal
# part, polynomiality of the third Hamiltonian component is D(f_-) in K[X,Q],
# where D=-3X^2*d_X+(6XQ-2)*d_Q.
theta = (
    X**-4
    + 6 * Q * X**-3
    + sp.Rational(45, 2) * Q**2 * X**-2
    + sp.Rational(135, 2) * Q**3 * X**-1
)


def third_component(f):
    return sp.expand(
        -3 * X**2 * sp.diff(f, X) + (6 * X * Q - 2) * sp.diff(f, Q)
    )


assert sp.expand(third_component(C * theta)) == sp.Rational(2835, 2) * C * Q**3

# Reconstruct the coefficients recursively with zero integration constants.
g4 = C
g3 = sp.integrate(12 * g4 / 2, Q)
g2 = sp.integrate((9 * g3 + 6 * Q * sp.diff(g3, Q)) / 2, Q)
g1 = sp.integrate((6 * g2 + 6 * Q * sp.diff(g2, Q)) / 2, Q)
assert g3 == 6 * C * Q
assert g2 == sp.Rational(45, 2) * C * Q**2
assert g1 == sp.Rational(135, 2) * C * Q**3
assert sp.expand(C * theta - (g4 / X**4 + g3 / X**3 + g2 / X**2 + g1 / X)) == 0

print("PASS: pole-order-four normalized residues form the universal theta line")
print("theta = X^-4 + 6Q X^-3 + 45Q^2 X^-2/2 + 135Q^3 X^-1/2")

# The kappa=-1 replacement chart has weights wt(X,Q,Z)=(1,-2,-1).
# Its weight-minus-four, Z-independent principal part can only use X^-4 and
# Q*X^-2.  Polynomiality for the replacement third Hamiltonian component
# forces their ratio to be -4.
theta_exceptional = X**-4 - 4 * Q * X**-2


def exceptional_third_component(f):
    return sp.expand(
        2 * X**3 * sp.diff(f, X)
        - (2 + 6 * X**2 * Q) * sp.diff(f, Q)
    )


assert exceptional_third_component(C * theta_exceptional) == 40 * C * Q
print("PASS: kappa=-1 normalized residues form the exceptional theta line")
print("theta_exceptional = X^-4 - 4Q X^-2")
