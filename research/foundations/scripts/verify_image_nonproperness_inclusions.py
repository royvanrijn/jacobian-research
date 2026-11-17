#!/usr/bin/env python3
"""Denominator-safe audit of the image and nonproperness inclusions.

The proof covers four charts/strata explicitly:

* x=0 in the affine source;
* finite roots T with P'(T) != 0;
* finite repeated roots with P'(T) = 0 (the escaping boundary); and
* the projective root T=infinity, using S=1/T.

It also checks the generic discriminant, C=0, and Gamma strata separately so
that Groebner elimination cannot silently lose a leading-coefficient stratum.
"""

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, C, T, R, S = sp.symbols("A B C T R S")
u = 1 + x*y
F = (
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
)
P = C*T**3 - 2*T**2 + B*T - 2*A
r = sp.diff(P, T)
Q = 27*A**2*C**2 - 18*A*B*C + 16*A + B**3*C - B**2


# AFFINE CHART x=0.  This is an entire source stratum, not a denominator-cleared
# remnant. It maps bijectively to the target plane C=0.
assert tuple(sp.expand(component.subs(x, 0)) for component in F) == (
    z + 4*y**2, y, 0,
)
x0_source = {x: 0, y: B, z: A - 4*B**2}
assert tuple(sp.expand(component.subs(x0_source)) for component in F) == (A, B, 0)


# AFFINE CHART x!=0.  T=y+1/x is a root and R=P'(T)=2/x. Conversely every
# finite simple root reconstructs exactly one source point with no other
# denominator.
t_source = y + 1/x
assert sp.factor(P.subs({A: F[0], B: F[1], C: F[2], T: t_source})) == 0
assert sp.factor(r.subs({A: F[0], B: F[1], C: F[2], T: t_source}) - 2/x) == 0

xr = 2/R
yr = T - R/2
zr = sp.Rational(5, 4)*R**2 - sp.Rational(3, 2)*T*R - C*R**3/8
A_trc = T**2 - C*T**3 + R*T/2
B_trc = 4*T + R - 3*C*T**2
for got, want in zip(F, (A_trc, B_trc, C)):
    assert sp.factor(got.subs({x: xr, y: yr, z: zr}) - want) == 0
assert sp.factor(P.subs({A: A_trc, B: B_trc}) ) == 0
assert sp.factor(r.subs({A: A_trc, B: B_trc}) - R) == 0


# FINITE BOUNDARY R=0.  Its target parameterization lies in V(Q), and exact
# elimination gives no equation besides Q.
A_boundary = T**2 - C*T**3
B_boundary = 4*T - 3*C*T**2
assert sp.factor(Q.subs({A: A_boundary, B: B_boundary})) == 0
boundary_ideal = [2*A - 2*A_trc, B - B_trc, R]
boundary_gb = sp.groebner(boundary_ideal, T, R, A, B, C, order="lex")
boundary_elimination = [
    sp.factor(polynomial.as_expr())
    for polynomial in boundary_gb.polys
    if not polynomial.as_expr().has(T, R)
]
assert boundary_elimination == [Q]


# NO LOST DISCRIMINANT STRATA.  Off Gamma, the rational inverse to the boundary
# parameterization works modulo Q. Its failed denominator is handled below,
# rather than removed by saturation.
D_gamma = 4 - 3*B*C
N_gamma = B - 9*A*C
T_generic = N_gamma/D_gamma
generic_A_numerator = sp.factor(sp.together(
    A - (T_generic**2 - C*T_generic**3)
).as_numer_denom()[0])
generic_B_numerator = sp.factor(sp.together(
    B - (4*T_generic - 3*C*T_generic**2)
).as_numer_denom()[0])
assert sp.factor(generic_A_numerator/Q) == 4 - 27*A*C**2
assert sp.factor(generic_B_numerator/Q) == 9*C

# The leading-coefficient slice C=0 is present: Q=0 gives A=B^2/16 and is
# reached by T=B/4 on the finite repeated-root boundary.
T_c0 = B/4
assert sp.factor(A_boundary.subs({C: 0, T: T_c0}) - B**2/16) == 0
assert sp.factor(B_boundary.subs({C: 0, T: T_c0}) - B) == 0
assert sp.factor(Q.subs({C: 0, A: B**2/16})) == 0

# The failed generic inverse D_gamma=0 is exactly Gamma, and Gamma is reached
# by T=2/(3C); it is not lost when the denominator is cleared.
T_gamma = sp.Rational(2)/(3*C)
A_gamma = sp.Rational(4)/(27*C**2)
B_gamma = sp.Rational(4)/(3*C)
assert sp.factor(A_boundary.subs(T, T_gamma) - A_gamma) == 0
assert sp.factor(B_boundary.subs(T, T_gamma) - B_gamma) == 0
assert sp.factor(Q.subs({A: A_gamma, B: B_gamma})) == 0
assert sp.factor(D_gamma.subs(B, B_gamma)) == 0
assert sp.factor(Q.subs(B, B_gamma)) == (27*A*C**2 - 4)**2/(27*C**2)


# PROJECTIVE ROOT AT INFINITY. Homogenization shows it exists only for C=0 and
# is simple there. The S=1/T chart is regular and returns to the finite x=0
# source stratum; it is not an additional escaping boundary component.
P_homogeneous = C*T**3 - 2*T**2*S + B*T*S**2 - 2*A*S**3
assert P_homogeneous.subs({T: 1, S: 0}) == C
assert sp.diff(P_homogeneous, S).subs({C: 0, T: 1, S: 0}) == -2

C_s = 2*S - B*S**2 + 2*A*S**3
D_s = 1 - B*S + 3*A*S**2
x_s = S/D_s
y_s = B - 3*A*S
z_s = sp.cancel((2*x_s - 3*x_s**2*y_s - C_s)/x_s**3)
assert D_s.subs(S, 0) == 1
assert sp.diff(C_s, S).subs(S, 0) == 2
for got, want in zip(F, (A, B, C_s)):
    assert sp.factor(got.subs({x: x_s, y: y_s, z: z_s}) - want) == 0
assert (
    sp.limit(x_s, S, 0),
    sp.limit(y_s, S, 0),
    sp.limit(z_s, S, 0),
) == (0, B, A - 4*B**2)


# IMAGE INCLUSIONS.  The preceding chart bijections reduce image membership to
# simple finite roots plus the x=0 point. The only C!=0 cubic with no simple
# root is a cube, whose coefficient comparison is exactly Gamma. Gamma never
# meets C=0, and every C=0 target has its x=0 preimage.
assert sp.factor((3*B_gamma*C - 4)) == 0
assert sp.factor(12*A_gamma - B_gamma**2) == 0
assert sp.factor(P.subs({A: A_gamma, B: B_gamma}) - C*(T - T_gamma)**3) == 0
assert (3*B*C - 4).subs(C, 0) == -4


# NONPROPERNESS V(Q) SUBSET S_F. Every boundary point has a repeated-root
# parameter T0. Keeping B,C fixed and varying T gives exact nearby targets and
# source points; R=P'(T) tends to zero, so x=2/R escapes. R is not identically
# zero, hence the punctured path can always avoid its finitely many zeros.
A_path = (C*T**3 - 2*T**2 + B*T)/2
R_path = 3*C*T**2 - 4*T + B
assert sp.factor(P.subs(A, A_path)) == 0
assert sp.factor(r.subs(A, A_path) - R_path) == 0
assert sp.Poly(R_path, T).degree() >= 1
path_source = {
    x: 2/R_path,
    y: T - R_path/2,
    z: sp.Rational(5, 4)*R_path**2
       - sp.Rational(3, 2)*T*R_path - C*R_path**3/8,
}
for got, want in zip(F, (A_path, B, C)):
    assert sp.factor(got.subs(path_source) - want) == 0


# NONPROPERNESS S_F SUBSET V(Q). Outside V(Q), every finite projective root is
# simple, so R is a unit and the finite reconstruction is regular. The only
# possible unbounded-root chart is S=0,C=0, where D_s is a unit and the source
# extends to the bounded x=0 point checked above. These charts exhaust P^1.
assert sp.factor(sp.discriminant(P, T) + 4*Q) == 0
assert sp.factor(r.subs({A: A_trc, B: B_trc}) - R) == 0
assert D_s.subs(S, 0) == 1


print("PASS: affine source is exhausted by x=0 and simple finite-root charts")
print("PASS: finite boundary R=0 maps onto every stratum of V(Q)")
print("PASS: C=0 infinity root is regular and returns to the x=0 chart")
print("PASS: image is C^3 minus Gamma with no denominator-lost stratum")
print("PASS: V(Q) is contained in the nonproperness set by exact escape paths")
print("PASS: outside V(Q), finite and infinity root charts remain bounded")
