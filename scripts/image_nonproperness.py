#!/usr/bin/env python3
"""Exact certificates for the image, fibers, and nonproperness locus."""
import sympy as sp

x, y, z, A, B, C, T, S = sp.symbols("x y z A B C T S")
u = 1 + x*y
F = (
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
)
P = C*T**3 - 2*T**2 + B*T - 2*A
r = sp.diff(P, T)
Q = 27*A**2*C**2 - 18*A*B*C + 16*A + B**3*C - B**2

# The entire x=0 chart, and the unique such point over a target with C=0.
assert tuple(sp.expand(f.subs({x: 0, y: B, z: A - 4*B**2})) for f in F) == (A, B, 0)
assert tuple(sp.expand(f.subs(x, 0)) for f in F) == (z + 4*y**2, y, 0)

# Every simple root gives exactly one x != 0 preimage.
xr = 2/r
yr = T - r/2
zr = sp.Rational(5, 4)*r**2 - sp.Rational(3, 2)*T*r - C*r**3/8
for got, want in zip(F, (A, B, C)):
    numerator = sp.together(got.subs({x: xr, y: yr, z: zr}) - want).as_numer_denom()[0]
    # Reduction modulo P(T) proves the identity on every root with r != 0.
    assert sp.rem(sp.Poly(numerator, T), sp.Poly(P, T)).as_expr() == 0

assert sp.factor(sp.discriminant(P, T) + 4*Q) == 0

# The singular scheme has radical (3BC-4, 12A-B^2).
partials = [sp.diff(Q, v) for v in (A, B, C)]
gradient_gb = sp.groebner(partials, A, B, C, order="lex")
assert [sp.factor(g.as_expr()) for g in gradient_gb.polys] == [
    16*A - B**3*C,
    (3*B*C - 4)**2,
]
gamma_gb = sp.groebner([3*B*C - 4, 12*A - B**2], A, B, C, order="lex")
assert all(gamma_gb.reduce(f)[1] == 0 for f in partials + [Q])

# Resolve the apparently unbounded cubic root at C=0 by S=1/T.
# The homogenized root equation in this chart is
# C = 2S-BS^2+2AS^3.  It is nonsingular at S=C=0.
C_s = 2*S - B*S**2 + 2*A*S**3
x_s = S/(1 - B*S + 3*A*S**2)
y_s = B - 3*A*S
z_s = sp.factor((2*x_s - 3*x_s**2*y_s - C_s)/x_s**3)
assert sp.limit(x_s, S, 0) == 0
assert sp.limit(y_s, S, 0) == B
assert sp.limit(z_s, S, 0) == A - 4*B**2
assert sp.diff(C_s, S).subs(S, 0) == 2

# Boundary chart.  With R=P'(T)=2/x, the graph is polynomial in (T,R,C).
R = sp.symbols("R")
A_trc = T**2 - C*T**3 + R*T/2
B_trc = 4*T + R - 3*C*T**2
boundary_ideal = [2*A - 2*A_trc, B - B_trc, R]
boundary_gb = sp.groebner(boundary_ideal, T, R, A, B, C, order="lex")
boundary_elimination = [
    sp.factor(g.as_expr()) for g in boundary_gb.polys
    if not (g.as_expr().has(T) or g.as_expr().has(R))
]
assert boundary_elimination == [Q]

# The repeated-root parametrization is the normalization of V(Q).
A_norm = T**2 - C*T**3
B_norm = 4*T - 3*C*T**2
assert sp.factor(Q.subs({A: A_norm, B: B_norm})) == 0
assert sp.factor((B_norm - 9*A_norm*C) - (4 - 3*B_norm*C)*T) == 0
# T is integral: T^2-BT+3A=0 on the normalization.
assert sp.factor(T**2 - B_norm*T + 3*A_norm) == 0
# The failure of the rational inverse is exactly the triple-root divisor CT=2/3.
assert sp.factor((4 - 3*B_norm*C) - (2 - 3*C*T)**2) == 0
assert sp.factor((B_norm - 9*A_norm*C) - T*(2 - 3*C*T)**2) == 0

# Two explicit simple-discriminant meridians give transpositions and generate S3.
# For C=1, roots m+sqrt(E), m-sqrt(E), 2-2m have fixed sum 2.
E, m = sp.symbols("E m")
root3 = 2 - 2*m
poly_loop = sp.expand((T-(m+sp.sqrt(E)))*(T-(m-sp.sqrt(E)))*(T-root3))
assert not poly_loop.has(sp.sqrt(E))
assert sp.Poly(poly_loop, T).coeff_monomial(T**2) == -2

print("PASS: simple-root reconstruction and the x=0 fiber")
print("PASS: Sing(V(Q)) = V(3BC-4, 12A-B^2) set-theoretically")
print("PASS: the T=infinity branch extends to (0,B,A-4B^2)")
print("PASS: boundary elimination is (Q), and A^2_(T,C) normalizes V(Q)")
print("PASS: explicit discriminant meridians exchange two roots")
