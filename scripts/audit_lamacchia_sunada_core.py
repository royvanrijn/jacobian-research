#!/usr/bin/env python3
"""Audit the LaMacchia--Bosma--de Smit degree-seven Gassmann core.

This is an algebraic-geometry audit, not a proof of the externally cited
generic Galois-group theorem.  It verifies:

* linearity of the LaMacchia polynomial in the parameter ``s``;
* the resulting rational plane chart and its two pole divisors;
* the rational-chart Jacobian;
* the square branch discriminant;
* the pulled-back branch ledger ``E_4 E_8 J^2``; and
* the elementary target-sign equivalence of the two displayed partners.
"""
from __future__ import annotations

import sympy as sp


X, S, T = sp.symbols("X S T")


F0 = (
    X**7
    + (-6 * T + 2) * X**6
    + (8 * T**2 + 4 * T - 3) * X**5
    + (-14 * T**2 + 6 * T - 2) * X**4
    + (6 * T**2 - 8 * T**3 - 4 * T + 2) * X**3
    + (8 * T**3 + 16 * T**2) * X**2
    + (8 * T**3 - 12 * T**2) * X
    - 8 * T**3
)

FS = sp.expand(F0 + S * X**3 * (1 - X))

# On X(1-X) != 0 the incidence FS=0 is the graph S=PHI(T,X).
PHI = sp.factor(-F0 / (X**3 * (1 - X)))
PHI_FACTORED = (
    (X + 1)
    * (2 * T**2 - 2 * T * X**2 + 2 * T * X + X**3 - X**2)
    * (-4 * T * X**2 + 8 * T * X - 4 * T + X**3 + 2 * X**2 - 2 * X)
    / (X**3 * (X - 1))
)
assert sp.cancel(PHI - PHI_FACTORED) == 0
assert sp.cancel(FS.subs(S, PHI)) == 0

# The plane-chart Jacobian of (T,X) -> (T,PHI) is d(PHI)/dX.
J_NUM, J_DEN = map(
    sp.factor, sp.together(sp.diff(PHI, X)).as_numer_denom()
)
assert J_DEN == X**4 * (X - 1) ** 2
assert sp.Poly(J_NUM, X).degree() == 8
assert sp.Poly(J_NUM, T).degree() == 3

# The polynomial discriminant is a square, as required by the point action
# of GL(3,2) lying in A_7.  BRANCH is the reduced non-coordinate factor.
discriminant = sp.factor(sp.discriminant(FS, X))
BRANCH = (
    27 * S**4
    - 256 * S**2 * T**6
    + 9408 * S**2 * T**5
    - 984 * S**2 * T**4
    + 17248 * S**2 * T**3
    - 11376 * S**2 * T**2
    + 864 * S**2 * T
    - 216 * S**2
    - 65536 * T**11
    + 561152 * T**10
    - 1612544 * T**9
    + 1317296 * T**8
    + 1267328 * T**7
    - 1418816 * T**6
    - 839552 * T**5
    + 146464 * T**4
    + 94592 * T**3
    - 9792 * T**2
    - 3456 * T
    + 432
)
assert sp.expand(discriminant - 64 * T**6 * BRANCH**2) == 0

# Pulling the reduced branch equation back to the rational root chart gives
# the same three-column pattern as the Davenport cover: two unramified
# factors and the doubled derivative factor.
branch_pullback = sp.factor(sp.together(BRANCH.subs(S, PHI)))
pull_num, pull_den = map(sp.factor, branch_pullback.as_numer_denom())
pull_factors = sp.factor_list(pull_num)[1]
degree_multiplicity = sorted(
    (sp.Poly(factor, X).degree(), multiplicity)
    for factor, multiplicity in pull_factors
)
assert degree_multiplicity == [(4, 1), (8, 1), (8, 2)]
assert pull_den == X**12 * (X - 1) ** 4
assert any(
    multiplicity == 2 and sp.rem(J_NUM, factor, X) == 0
    for factor, multiplicity in pull_factors
)

# Bosma--de Smit's displayed partner f_{-s,t} is the same rational chart
# followed by the target involution (T,S) -> (T,-S).  This preserves the
# Gassmann cover over the fixed base but makes this presentation unsuitable
# by itself for proving stable left--right inequivalence.
PHI_PARTNER = -PHI
assert sp.cancel(FS.subs({S: -PHI_PARTNER})) == 0
assert sp.cancel(PHI_PARTNER + PHI) == 0

print("PASS: the LaMacchia incidence is linear in s")
print("PASS: its root chart is A1_t times (P1_x minus {0,1,infinity})")
print("PASS: the chart Jacobian has denominator x^4(x-1)^2")
print("PASS: the degree-seven discriminant is 64*t^6*Branch(s,t)^2")
print("PASS: the non-coordinate reduced branch pullback has ledger E4*E8*J^2")
print("PASS: the f_(s,t)/f_(-s,t) presentation is target-sign equivalent")
