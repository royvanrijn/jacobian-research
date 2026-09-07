#!/usr/bin/env python3
"""Exact identities for the primitive-root fourth-suspension reduction."""

from __future__ import annotations

import sympy as sp


P, S, Q, T, Q_tilde = sp.symbols("P S Q T Q_tilde")
a, c, lam = sp.symbols("a c lam", nonzero=True)
b = sp.Function("b")(P)
H = sp.Function("H")(P, T)

# A 1/P translation leaves a nonpolynomial coefficient of Q_tilde.
S_affine = (T - b - c / P) / a
Q_affine = a * (Q_tilde - H)
D_affine = sp.expand(1 - S_affine * Q_affine + P * S_affine**2)
assert sp.factor(sp.diff(D_affine, Q_tilde) + T - b - c / P) == 0

# The reciprocal scaling has a polynomial controlled divisor.
D_scaling = sp.expand(
    (1 - S * Q + P * S**2).subs(
        {S: P * T, Q: (Q_tilde - H) / P},
        simultaneous=True,
    )
)
assert sp.factor(D_scaling - (1 - T * Q_tilde + T * H + P**3 * T**2)) == 0

# Its forced quadratic incidence has the unavoidable Q*S^2/P coefficient.
f = sp.Function("f")(P)
X = lam * T**2 / 2 + f
q_coefficient_source = sp.factor((-P * X).subs(T, S / P))
assert sp.factor(q_coefficient_source + lam * S**2 / (2 * P) + P * f) == 0

# General unimodular Mobius transformation.
A, B, C, D = sp.symbols("A B C D")
Delta = A * D - B * C
L = D * T - B
K = A - C * T
S_inverse = L / K
phi_prime_in_T = K**2 / Delta
Q_inverse = phi_prime_in_T * (Q_tilde - H)
D_mobius = sp.factor(1 - S_inverse * Q_inverse + P * S_inverse**2)
expected = sp.factor(
    1
    - L * K * Q_tilde / Delta
    + L * K * H / Delta
    + P * L**2 / K**2
)
assert sp.factor(D_mobius - expected) == 0

# The unique principal part needed to remove P*L^2/K^2 is cubic in K.
H_polar = -Delta * P * L / K**3
assert sp.factor(L * K * H_polar / Delta + P * L**2 / K**2) == 0

print("PASS: a 1/P translation makes the controlled coefficient nonpolynomial")
print("PASS: T=S/P polynomializes the controlled divisor")
print("PASS: its forced quadratic incidence has an exact Q*S^2/P pole")
print("PASS: a nonmonomial unimodular denominator requires a K^-3 shear")
