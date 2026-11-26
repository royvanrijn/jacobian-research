#!/usr/bin/env python3
"""Exact finite controls for the critical-section threefold and Arala multiplicities."""

from fractions import Fraction
from itertools import combinations
from math import comb

import sympy as sp


def check_gauge():
    w, T, S, V, Q, Y, Z = sp.symbols("w T S V Q Y Z")
    P = T * (T - 1)**3 + w * (T + 1)
    derivative = sp.diff(P, T)
    tau0 = sp.Rational(1, 4)
    slope = sp.diff(derivative, T).subs({w: 0, T: tau0})
    gamma = sp.diff(derivative, T, 2).subs({w: 0, T: tau0}) / 6
    assert derivative.subs({w: 0, T: tau0}) == 0 and slope * gamma != 0

    tau = tau0
    for j in range(1, 6):
        target = 1 / gamma if j == 5 else 0
        current = sp.expand(derivative.subs(T, tau)).coeff(w, j)
        tau += (target - current) / slope * w**j
    coefficients = sp.Poly(P.subs(T, tau + S), S)
    c = [coefficients.nth(j) for j in range(5)]
    d = sp.cancel(c[1] / w**5)
    assert sp.Poly(d, w, domain=sp.QQ).nth(0) == 1 / gamma
    scale_correction = 1 - sp.diff(c[3] * d, w).subs(w, 0) * w / 2
    a = w**3 * d * scale_correction
    profiles = {
        3: sp.cancel((w * c[3] * d * scale_correction**2 - w) / w**3),
        4: sp.expand(c[4] * d**2 * scale_correction**3),
    }
    for profile in profiles.values():
        assert sp.Poly(profile, w, domain=sp.QQ).degree() >= 0

    # A symbolic Jacobian computation applies to any independent cubic
    # and quartic profiles. Their V derivatives disappear from the minor.
    h3, h4 = sp.symbols("h3 h4")
    G = S + V * S**3 + h3 * V**3 * S**3 + h4 * V**4 * S**4
    B = Q + sp.cancel((sp.diff(G, S) - 1 - V*S**2) / S)
    C = 2 * G - B * S**2
    D = 1 - S * Q + V * S**2
    assert sp.expand(sp.Matrix([V, -B/2, C]).jacobian((V, S, Q)).det() - D) == 0

    x, y, q = sp.symbols("x y q")
    t = 1 + x*y
    chart = {V: t*q, S: x/t, Q: y+x*q}
    assert sp.cancel(sp.Matrix(list(chart.values())).jacobian((x, y, q)).det() - 1/t) == 0
    assert sp.cancel(D.subs(chart) - 1/t) == 0
    for j, profile in ((3, h3), (4, h4)):
        assert sp.cancel(sp.diff(B, profile).subs(chart) - j*t**2*x**(j-2)*q**j) == 0
        assert sp.cancel(sp.diff(C, profile).subs(chart) + (j-2)*(x*q)**j) == 0

    w0 = 1
    actual_tau, actual_a, actual_c1 = [f.subs(w, w0) for f in (tau, a, c[1])]
    assert actual_a * actual_c1 != 0
    h = {symbol: profiles[j].subs(w, w0) for j, symbol in ((3, h3), (4, h4))}
    Y0, Z0 = (c[2]*a/c[1]).subs(w, w0), (-2*c[0]/(c[1]*a)).subs(w, w0)
    modulus = sp.Poly((G + Y*S**2 - Z/2).subs({V: w0, Y: Y0, Z: Z0, **h}), S)
    normalized = P.subs({w: w0,T: actual_tau + actual_a*S}) / (actual_a * actual_c1)
    assert modulus == sp.Poly(normalized, S)
    assert modulus.degree() == 4 and sp.gcd(modulus, modulus.diff()).degree() == 0


def check_arala_counts():
    for m in (3, 5, 7):
        n = m*m + 1
        b = n//2
        sizes = range(1, b, 2)
        A = sum(comb(n-1, k) for k in sizes) + comb(n, b)//2
        B = sum(comb(n-1, k-1) for k in sizes)
        assert A + B == 2**(n-2)
        assert A > B > 0
        t = Fraction(A, A + B)
        assert A/t + B/(t-1) == 0
        assert -A/t**2 - B/(t-1)**2 != 0
        assert 2*A/t**3 + 2*B/(t-1)**3 != 0

    # Enumerate the rank-256 resolvent expressions on the actual multiset.
    values = [int(9 in subset) for k in (1, 3) for subset in combinations(range(10), k)]
    for subset in combinations(range(10), 5):
        if 0 in subset:
            values.append((10 * int(9 in subset) - 5)**2 - 25)
    assert len(values) == 256 and values.count(0) == 219 and values.count(1) == 37


if __name__ == "__main__":
    check_gauge()
    check_arala_counts()
    print("PASS: critical-section normalization, cubic gauge, quartic fibre identity, Arala counts")
