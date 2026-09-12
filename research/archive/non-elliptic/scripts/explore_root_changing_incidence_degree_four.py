#!/usr/bin/env python3
"""Exact test of the smallest root-changing incidence rechart.

An arbitrary birational change of one marked root is Mobius.  This script
lifts that change to the conjugate chart coordinate with fibre determinant
one, chooses the unique shear that keeps the controlled divisor polynomial,
and checks the resulting source-boundary pole.
"""

from __future__ import annotations

import sympy as sp


P, S, Q = sp.symbols("P S Q")
a, b, c, d = sp.symbols("a b c d", nonzero=True)
Delta = a * d - b * c

# The original quadratic reciprocal chart and a general Mobius root change.
D = 1 - S * Q + P * S**2
T = (a * S + b) / (c * S + d)
phi_prime = sp.factor(sp.diff(T, S))
assert sp.factor(phi_prime - Delta / (c * S + d) ** 2) == 0

# Invert the root change and express its derivative in the new root.
S_of_T = sp.factor((d * T - b) / (a - c * T))
phi_prime_of_T = sp.factor((a - c * T) ** 2 / Delta)
assert sp.factor(S_of_T - S) == 0
assert sp.factor(phi_prime_of_T - phi_prime) == 0

# The determinant-one cotangent lift is Q_tilde=Q/phi'+H(P,T).
# The unique shear removing the P*S^2 term from D is H=-P*S/phi'.
H = sp.factor(-P * S / phi_prime)
Q_tilde = sp.factor(Q / phi_prime + H)
fibre_jacobian = sp.factor(
    sp.det(sp.Matrix([T, Q_tilde]).jacobian((S, Q)))
)
assert fibre_jacobian == 1
assert sp.factor(Q_tilde - (Q - P * S) / phi_prime) == 0

kappa = sp.factor(S_of_T * phi_prime_of_T)
transformed_D = sp.factor(
    D.subs(Q, sp.factor(phi_prime * (Q_tilde - H)))
)
assert sp.factor(transformed_D - (1 - kappa * Q_tilde)) == 0

# Kappa is generically quadratic, so its primitive is a cubic horizontal
# coordinate.  It drops below degree two precisely on special Mobius strata.
T_symbol = sp.symbols("T")
kappa_polynomial = sp.factor(
    ((d * T_symbol - b) * (a - c * T_symbol)) / Delta
)
assert sp.Poly(kappa_polynomial, T_symbol).degree() == 2
assert sp.factor(
    sp.Poly(kappa_polynomial, T_symbol).coeff_monomial(T_symbol**2)
    + c * d / Delta
) == 0
X_cubic = sp.integrate(kappa_polynomial, T_symbol)
assert sp.Poly(X_cubic, T_symbol).degree() == 3

# Pull back to the original affine source chart.  Only the identities
# S=x/t, Q-P*S=y, and P=t*q are needed.
x, y, t, q = sp.symbols("x y t q")
T_source = sp.factor((a * x + b * t) / (c * x + d * t))
Q_tilde_source = sp.factor(y * (c * x + d * t) ** 2 / (Delta * t**2))
P_source = t * q

# At the generic point of t=0 (where x,c,Delta are units), T and P are
# regular, whereas Q_tilde has an exact pole of order two.
assert sp.factor(T_source.subs(t, 0) - a / c) == 0
assert P_source.subs(t, 0) == 0
pole_residue = sp.factor((t**2 * Q_tilde_source).subs(t, 0))
assert pole_residue == c**2 * x**2 * y / Delta

# A polynomial incidence shear beta(P,T) is regular at this divisor and
# therefore cannot cancel the nonzero Q_tilde residue.
u, v = sp.symbols("u v")
beta = (
    2
    + 3 * P_source
    + 5 * T_source
    + 7 * P_source**2
    + 11 * P_source * T_source
    + 13 * T_source**2
)
assert sp.factor((t**2 * beta).subs(t, 0)) == 0
combined_residue = sp.factor(
    sp.limit(t**2 * (Q_tilde_source + beta), t, 0)
)
assert combined_residue == pole_residue

print("PASS: the general Mobius root change has a determinant-one fibre lift")
print("PASS: the unique controlled-divisor shear gives D=1-kappa(T)*Q_tilde")
print("PASS: generic kappa is quadratic and yields a cubic horizontal coordinate")
print("PASS: the pulled-back incidence B-coordinate has an unavoidable t^-2 pole")
