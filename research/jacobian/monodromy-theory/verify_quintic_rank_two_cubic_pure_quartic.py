"""Check the leading binary series for the rank-two cubic primitive pencil."""

from sympy import Rational, cancel, expand, gcd, series, sqrt, symbols


a, k, u, v = symbols("a k u v")
w = (-a + sqrt(a**2 + 4*(1-k*a**2))) / 2
z = sqrt(w)
unadjusted = z * (1 + a*w/4 + k*a**2/4 + 5*a**2/32)

# Nonzero order-one coefficient: Z-W-A/4 at excess two.
leading_one = series(z - w - a/4, a, 0, 5).removeO()
assert cancel(leading_one.coeff(a, 2) - (k/4 - Rational(3, 32))) == 0
specialized = leading_one.subs(k, Rational(3, 8))
assert specialized.coeff(a, 3) == -Rational(1, 64)
assert specialized.coeff(a, 4) == -Rational(5, 1024)
w_at_three_eighths = w.subs(k, Rational(3, 8))
corrected = series(specialized + a**3*(w_at_three_eighths + a/4)**2/64, a, 0, 5)
assert corrected.removeO().coeff(a, 4) == -Rational(13, 1024)

# Nonzero order-two coefficient: removing Q adds Z*(u*A**3+v*A**4).
excess_one = series(1 - w - a/2, a, 0, 3).removeO()
assert cancel(excess_one.coeff(a, 2) - (k/2 - Rational(1, 8))) == 0
z_at_one_quarter = sqrt(1 - a/2)
pi = series(
    unadjusted.subs(k, Rational(1, 4)) + z_at_one_quarter*(u*a**3 + v*a**4),
    a,
    0,
    8,
).removeO()
assert expand(pi.coeff(a, 4) + (512*u - 2048*v + 15)/2048) == 0
assert expand(pi.coeff(a, 5) + (128*u + 1024*v + 9)/4096) == 0
pi = pi.subs({u: -Rational(11, 256), v: -Rational(7, 2048)})
assert pi.coeff(a, 4) == pi.coeff(a, 5) == 0
assert pi.coeff(a, 7) == -Rational(81, 524288)

# Zero first two coefficients: removals at orders three and six.
pi = series(unadjusted, a, 0, 8).removeO()
assert expand(pi.coeff(a, 3) - 5*(1-8*k)/128) == 0
assert pi.coeff(a, 4).subs(k, Rational(1, 8)) == -Rational(5, 1024)
q = series(w, a, 0, 5).removeO() + (pi.coeff(a, 4)/pi.coeff(a, 3) + Rational(1, 2))*a
adjusted = expand(pi - pi.coeff(a, 3)*a**3*q)
assert cancel(adjusted.coeff(a, 4)) == 0
fifth = 960*k**2 - 320*k + 29
assert cancel(adjusted.coeff(a, 5) + fifth/4096) == 0
adjusted = expand(adjusted - adjusted.coeff(a, 6)*a**6*q**2)
seventh = (4*k-1)*(512*k**3 - 704*k**2 + 200*k - 23)
assert cancel(adjusted.coeff(a, 7) - 5*seventh/131072) == 0
assert gcd(fifth, seventh) == 1

print("Verified the primitive cubic pure-quartic slice coefficients and gcd")
