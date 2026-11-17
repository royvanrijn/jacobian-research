#!/usr/bin/env python3
"""Exact symbolic audit of the degree-five stable-moduli family."""

import sympy as sp

w, lam = sp.symbols("w lambda")
A, B, C = sp.symbols("A B C")
x, y, z = sp.symbols("x y z")
s, t = sp.symbols("s t")

H = sp.factor(
    w**2 * (w - 1) * (3 * w**2 - (5 * lam + 1) * w + 3 * lam) / 60
)
p = sp.diff(H, w)
c = sp.factor(-p.subs(w, 1))
kappa = sp.factor(sp.diff(p, w).subs(w, 1) / c)
R = 3 * w**2 - (5 * lam + 1) * w + 3 * lam

# Weighted admissibility and the deliberately prescribed Hessian roots.
assert H.subs(w, 0) == 0
assert p.subs(w, 0) == 0
assert H.subs(w, 1) == 0
assert sp.factor(c - (lam - 1) / 30) == 0
assert kappa == -9
assert sp.factor(
    sp.diff(H, w, 2) - (w - lam) * (10 * w**2 - 8 * w + 1) / 10
) == 0

# The additional primitive roots are distinct, nonzero, and different from one
# away from the first three factors of the exceptional polynomial.
assert sp.factor(sp.discriminant(R, w) - (25 * lam - 1) * (lam - 1)) == 0
assert R.subs(w, 0) == 3 * lam
assert sp.factor(R.subs(w, 1) - 2 * (1 - lam)) == 0

# At an additional root rho, y has a C^{-1} pole unless c*rho+H'(rho)=0.
# The last exceptional factor is exactly the resultant obstruction.
escape_numerator = sp.together(c * w + p).as_numer_denom()[0]
escape_resultant = sp.factor(sp.resultant(R, escape_numerator, w))
assert sp.factor(
    escape_resultant
    / (36 * lam * (lam - 1) ** 2 * (100 * lam**2 - 29 * lam + 10))
) == 1

# Exact C-adic discriminant saturation and its reduced boundary trace.
E = H - s * w + t
discriminant = sp.discriminant(E, w)
pulled = sp.factor(discriminant.subs({s: B * C, t: c * A * C**2}))
poly_c = sp.Poly(sp.expand(pulled), C)
minimum_c_power = min(
    monomial[0] for monomial, coefficient in poly_c.terms() if coefficient != 0
)
assert minimum_c_power == 2

trace = sp.factor((pulled / C**2).subs(C, 0))
expected_trace = (
    lam**2
    * (lam - 1) ** 3
    * (25 * lam - 1)
    * (A * lam * (lam - 1) + 150 * B**2)
    / sp.Integer(194400000000)
)
assert sp.factor(trace - expected_trace) == 0

# Construct the actual weighted map.  The apparent x and gamma denominators
# cancel; only a parameter unit lambda-1 remains in the A-coordinate.
u = 1 + x * y
gamma = 1 - sp.Rational(8, 7) * x * y + x**2 * z
W = sp.expand(u * gamma)
q = sp.factor((w * p - H) / c)

C_map = sp.expand(x * gamma)
B_map = sp.cancel((c + p.subs(w, W) / gamma) / x)
A_map = sp.cancel((u + q.subs(w, W) / gamma**2) / x**2)

assert sp.fraction(B_map)[1] == 1
assert sp.factor(sp.fraction(A_map)[1] / (686 * (lam - 1))) == 1
assert sp.factor(H.subs(w, W) - B_map * C_map * W + c * A_map * C_map**2) == 0

jacobian = sp.factor(
    sp.Matrix([A_map, B_map, C_map]).jacobian([x, y, z]).det()
)
assert sp.factor(jacobian - c) == 0

# The anharmonic invariant of the unordered Hessian-root triple is rational in
# lambda, despite the two fixed quadratic roots.
sqrt6 = sp.sqrt(6)
p0 = (4 - sqrt6) / 10
q0 = (4 + sqrt6) / 10
chi = sp.factor((lam - p0) / (q0 - p0))

J = sp.factor(
    sp.Rational(64, 75)
    * (50 * lam**2 - 40 * lam + 17) ** 3
    / (10 * lam**2 - 8 * lam + 1) ** 2
)
J_from_chi = sp.factor(256 * (1 - chi + chi**2) ** 3 / (chi**2 * (1 - chi) ** 2))
assert sp.simplify(J - J_from_chi) == 0

anharmonic_orbit = (
    chi,
    1 - chi,
    1 / chi,
    1 / (1 - chi),
    chi / (chi - 1),
    (chi - 1) / chi,
)

def j_invariant(value):
    return sp.factor(256 * (1 - value + value**2) ** 3 / (value**2 * (1 - value) ** 2))

for orbit_value in anharmonic_orbit:
    assert sp.simplify(j_invariant(orbit_value) - J_from_chi) == 0

print("PASS: H_lambda is an admissible degree-five weighted seed on Lambda")
print("PASS: the compact formulas define a polynomial Keller map with det=(lambda-1)/30")
print("PASS: the marked inverse equation has degree five and exact C^2 saturation")
print("PASS: all split-root and boundary exceptions are exactly the stated factors")
print("PASS: J(lambda) is the six-orbit invariant of the Hessian-root triple")
