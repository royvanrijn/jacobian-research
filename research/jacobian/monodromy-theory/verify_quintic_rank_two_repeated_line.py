#!/usr/bin/env python3
"""Check the rational coefficients for the pure-quartic formal slice."""

from sympy import Rational, expand, series, sqrt, symbols


a, b, ell, k = symbols("a b ell k")

# Binary quadratic in R_3, with Z^4+a*Z=1.
z = 1 - a / 4 - a**2 / 32
assert series(z**4 + a*z - 1, a, 0, 4).removeO() == 0
z_to_order_three = series(z * (1 + a*z / 4), a, 0, 4).removeO()
assert z_to_order_three.coeff(a, 2) == -Rational(5, 32)
assert z_to_order_three.coeff(a, 3) == 0
assert expand(z_to_order_three + Rational(5, 32)*a**2*z).coeff(a, 3) == -Rational(5, 128)

# Binary linear in R_3 and binary quadratic in R_2.
w = (-a + sqrt(a**2 + 4*(1-b))) / 2
z = sqrt(w)
pi = z * (1 + a*w / 4 + b / 4 + 5*a**2 / 32)
cubic = series(pi.subs({a: a*ell, b: a**2*k}, simultaneous=True), a, 0, 4)
assert expand(cubic.removeO().coeff(a, 3) - 5*ell*(ell**2 - 8*k) / 128) == 0
z = z.subs(b, a**2 / 8)
pi = series(pi.subs(b, a**2 / 8), a, 0, 6).removeO()
assert pi.coeff(a, 4) == -Rational(5, 1024)
assert pi.coeff(a, 5) == -Rational(1, 1024)
corrected = series(pi + Rational(5, 1024)*a**4*z, a, 0, 6).removeO()
assert corrected.coeff(a, 5) == -Rational(9, 4096)

# Binary quadratic in R_2 with R_3 removed by translation.
z = (1-b)**Rational(1, 4)
pi = series(z * (1 + b/4), b, 0, 4).removeO()
assert pi.coeff(b, 2) == -Rational(5, 32)
assert pi.coeff(b, 3) == -Rational(5, 64)
corrected = series(pi + Rational(5, 32)*b**2*z, b, 0, 4).removeO()
assert corrected.coeff(b, 3) == -Rational(15, 128)

print("Verified the pure-quartic slice coefficients and corrections")
