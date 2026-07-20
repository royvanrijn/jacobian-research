#!/usr/bin/env python3
"""Exact case audit for every degeneration of the inverse cubic.

This complements ``image_nonproperness.py`` by making the exceptional case
split executable: C != 0 versus C = 0, simple/double/triple finite roots, the
projective root at infinity, and every denominator in the reconstruction.
"""

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, C, T, t, S = sp.symbols("A B C T t S")
u = 1 + x*y
F = (
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
)
P = C*T**3 - 2*T**2 + B*T - 2*A
r = sp.diff(P, T)
Q = 27*A**2*C**2 - 18*A*B*C + 16*A + B**3*C - B**2


# The x != 0 chart has no hidden reconstruction denominator: a finite root is
# reconstructible exactly when it is simple (r != 0).
xr = 2/r
yr = T - r/2
zr = sp.Rational(5, 4)*r**2 - sp.Rational(3, 2)*T*r - C*r**3/8
assert sp.denom(xr) == r
assert sp.denom(yr) == 1
assert sp.denom(zr) == 1
for got, want in zip(F, (A, B, C)):
    numerator = sp.together(got.subs({x: xr, y: yr, z: zr}) - want).as_numer_denom()[0]
    assert sp.rem(sp.Poly(numerator, T), sp.Poly(P, T)).as_expr() == 0


# For C != 0, parameterize every repeated finite root t by solving P(t)=P'(t)=0.
# This avoids relying only on the discriminant equation.
A_repeated = t**2 - C*t**3
B_repeated = 4*t - 3*C*t**2
P_repeated = sp.factor(P.subs({A: A_repeated, B: B_repeated}))
expected_factorization = sp.expand(C*(T - t)**2*(T - (sp.Rational(2)/C - 2*t)))
assert sp.expand(P_repeated - expected_factorization) == 0
assert sp.factor(r.subs({A: A_repeated, B: B_repeated, T: t})) == 0
assert sp.factor(Q.subs({A: A_repeated, B: B_repeated})) == 0

# The remaining root is simple unless all three roots coincide.  Triple roots
# occur exactly at t=2/(3C), giving the curve Gamma.
third_root = sp.Rational(2)/C - 2*t
assert sp.solve(sp.together(third_root - t), t) == [sp.Rational(2)/(3*C)]
triple_t = sp.Rational(2)/(3*C)
triple_A = sp.factor(A_repeated.subs(t, triple_t))
triple_B = sp.factor(B_repeated.subs(t, triple_t))
assert triple_A == sp.Rational(4)/(27*C**2)
assert triple_B == sp.Rational(4)/(3*C)
assert sp.expand(P.subs({A: triple_A, B: triple_B}) - C*(T - triple_t)**3) == 0
assert sp.factor(3*triple_B*C - 4) == 0
assert sp.factor(12*triple_A - triple_B**2) == 0

# Conversely, the two Gamma equations force the same triple-root coefficients
# when C != 0.
gamma_substitution = {B: sp.Rational(4)/(3*C), A: sp.Rational(4)/(27*C**2)}
assert sp.factor(P.subs(gamma_substitution) - C*(T - triple_t)**3) == 0


# For C=0 the finite polynomial is always quadratic: its T^2 coefficient never
# vanishes.  Its discriminant is -Q, and its missing cubic root is [T:S]=[1:0].
P_zero = sp.Poly(P.subs(C, 0), T)
assert P_zero.degree() == 2
assert sp.factor(sp.discriminant(P_zero.as_expr(), T) + Q.subs(C, 0)) == 0
assert sp.factor(P.subs({C: 0, A: B**2/16}) + 2*(T - B/4)**2) == 0
assert sp.factor(Q.subs({C: 0, A: B**2/16})) == 0
assert sp.factor((3*B*C - 4).subs(C, 0)) == -4  # Gamma never meets C=0.

# Homogenizing P to projective degree three shows that the infinity root for
# C=0 is simple: d/dS is -2 at [T:S]=[1:0].
P_homogeneous = C*T**3 - 2*T**2*S + B*T*S**2 - 2*A*S**3
assert P_homogeneous.subs({C: 0, T: 1, S: 0}) == 0
assert sp.diff(P_homogeneous, S).subs({C: 0, T: 1, S: 0}) == -2


# Resolve that infinity root with S=1/T.  Its chart is regular at S=C=0 and
# maps exactly to the finite x=0 source point, rather than to an escaping pole.
C_at_infinity = 2*S - B*S**2 + 2*A*S**3
D = 1 - B*S + 3*A*S**2
x_at_infinity = S/D
y_at_infinity = B - 3*A*S
z_at_infinity = sp.cancel(
    (2*x_at_infinity - 3*x_at_infinity**2*y_at_infinity - C_at_infinity)
    / x_at_infinity**3
)
assert D.subs(S, 0) == 1
assert sp.diff(C_at_infinity, S).subs(S, 0) == 2
assert sp.limit(x_at_infinity, S, 0) == 0
assert sp.limit(y_at_infinity, S, 0) == B
assert sp.limit(z_at_infinity, S, 0) == A - 4*B**2
for got, want in zip(F, (A, B, C_at_infinity)):
    assert sp.factor(got.subs({
        x: x_at_infinity,
        y: y_at_infinity,
        z: z_at_infinity,
    }) - want) == 0

# In the S chart, D=0 away from S=0 is exactly r=0: it is the same repeated-root
# pole, not an additional exceptional denominator.
r_infinity_chart = sp.factor(r.subs({T: 1/S, C: C_at_infinity}))
assert sp.factor(r_infinity_chart - 2*D/S) == 0


# Representatives exercise all five fibre types in the case table.
def finite_root_profile(a, b, c):
    polynomial = sp.Poly(P.subs({A: a, B: b, C: c}), T)
    roots = sp.roots(polynomial.as_expr(), T)
    simple = sum(1 for multiplicity in roots.values() if multiplicity == 1)
    repeated = sorted(roots.values())
    infinity = 1 if c == 0 else 0
    return simple, repeated, infinity, simple + infinity


representatives = {
    "C!=0, Q!=0": ((1, 0, 1), (3, [1, 1, 1], 0, 3)),
    "C!=0, double root": ((0, 0, 1), (1, [1, 2], 0, 1)),
    "C!=0, triple root/Gamma": (
        (sp.Rational(4, 27), sp.Rational(4, 3), 1),
        (0, [3], 0, 0),
    ),
    "C=0, two finite simple roots": (
        (-sp.Rational(1, 4), 0, 0),
        (2, [1, 1], 1, 3),
    ),
    "C=0, finite double root": ((0, 0, 0), (0, [2], 1, 1)),
}
for name, (target, expected) in representatives.items():
    actual = finite_root_profile(*target)
    assert actual == expected, (name, actual, expected)


print("PASS: C!=0 repeated roots split into double+simple or Gamma triple")
print("PASS: C=0 always has a quadratic finite part and one simple infinity root")
print("PASS: the infinity root extends exactly to the unique x=0 source point")
print("PASS: reconstruction poles occur exactly at repeated finite roots")
print("PASS: representatives cover fibre cardinalities 3, 1, and 0")
