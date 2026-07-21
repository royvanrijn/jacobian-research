#!/usr/bin/env python3
"""Positive covariance regression for hostile presentations of foundational Keller map."""

import sympy as sp


x, y, z = sp.symbols("x y z")
X, Y, Z = sp.symbols("X Y Z")
a, b, c = sp.symbols("a b c")
A, B, C = sp.symbols("A B C")
T, Theta, V, W = sp.symbols("T Theta V W")


def jacobian(expressions, variables):
    return sp.factor(sp.det(sp.Matrix(expressions).jacobian(variables)))


# Original foundational Keller map map.
u = 1 + x * y
F = (
    u**3 * z + y**2 * u * (4 + 3 * x * y),
    y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)
assert jacobian(F, (x, y, z)) == -2


# Nonlinear triangular source automorphism and its inverse.
alpha = (X + Y**2 + Z, Y + Z**2, Z)
alpha_inverse = (x - (y - z**2) ** 2 - z, y - z**2, z)
assert jacobian(alpha, (X, Y, Z)) == 1
assert tuple(sp.expand(expr.subs(dict(zip((x, y, z), alpha))))
             for expr in alpha_inverse) == (X, Y, Z)


# Nonlinear triangular target automorphism and inverse.
beta = (a + b**2 + c, b + c**2, c)
beta_inverse = (A - (B - C**2) ** 2 - C, B - C**2, C)
assert jacobian(beta, (a, b, c)) == 1
assert tuple(sp.expand(expr.subs(dict(zip((a, b, c), beta_inverse))))
             for expr in beta) == (A, B, C)


# G=beta F alpha and the transported marked root.
F_alpha = tuple(sp.expand(expr.subs(dict(zip((x, y, z), alpha)))) for expr in F)
G = tuple(sp.expand(expr.subs(dict(zip((a, b, c), F_alpha)))) for expr in beta)
assert jacobian(G, (X, Y, Z)) == -2

t_alpha = sp.cancel(alpha[1] + 1 / alpha[0])
P_old = c * T**3 - 2 * T**2 + b * T - 2 * a
marked_identity = P_old.subs({a: F_alpha[0], b: F_alpha[1], c: F_alpha[2], T: t_alpha})
assert sp.factor(marked_identity) == 0


# Target transport: transformed resolvent and discriminant are pullbacks.
P_new = sp.expand(P_old.subs({a: beta_inverse[0], b: beta_inverse[1], c: C}))
disc_old = sp.factor(sp.discriminant(P_old, T))
disc_transport = sp.factor(disc_old.subs({a: beta_inverse[0], b: beta_inverse[1], c: C}))
assert sp.factor(sp.discriminant(P_new, T) - disc_transport) == 0


# Quadratic primitive element Theta=T+T^2.
M = sp.factor(sp.resultant(P_old, Theta - T - T**2, T))
assert sp.Poly(M, Theta).degree() == 3
generator_collision = 2 * a * c - b * c - 2 * b - c**2 - 4 * c - 4
assert sp.factor(
    sp.discriminant(M, Theta) - disc_old * generator_collision**2
) == 0
t_from_theta = ((c + 2) * Theta + 2 * a) / (c * (Theta + 1) + b + 2)
linear_recovery = sp.expand(
    (c * (Theta + 1) + b + 2) * T - ((c + 2) * Theta + 2 * a)
)
reduced_recovery = sp.rem(
    sp.Poly(linear_recovery.subs(Theta, T + T**2), T),
    sp.Poly(P_old, T),
).as_expr()
assert sp.factor(reduced_recovery) == 0
assert sp.factor(M.subs(Theta, T + T**2) / P_old).is_polynomial(T)
assert sp.factor(t_from_theta.subs(Theta, T + T**2) - T).subs(
    a, (c * T**3 - 2 * T**2 + b * T) / 2
) == 0


# Möbius root coordinate V=(T+1)/(T+2), determinant one.
t_from_v = (2 * V - 1) / (1 - V)
P_mobius = sp.cancel((1 - V) ** 3 * P_old.subs(T, t_from_v))
assert sp.Poly(P_mobius, V).degree() == 3
assert sp.factor(sp.discriminant(P_mobius, V) - disc_old) == 0
v_from_t = (T + 1) / (T + 2)
assert sp.factor(t_from_v.subs(V, v_from_t) - T) == 0


# Stabilization: block Jacobian and extended target boundary ideal.
G_stable = G + (W,)
assert jacobian(G_stable, (X, Y, Z, W)) == -2
assert W not in disc_transport.free_symbols
assert sp.factor(disc_transport.subs(W, W + P_new) - disc_transport) == 0


print("PASS: nonlinear source/target conjugation transports the foundational Keller map cover")
print("PASS: quadratic primitive element defines the same field despite an extra raw discriminant square")
print("PASS: Möbius root chart preserves the projective incidence/discriminant")
print("PASS: stabilization extends the canonical boundary object cylindrically")
