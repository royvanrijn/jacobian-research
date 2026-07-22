#!/usr/bin/env python3
"""Exact second-order contact profile for the degree-five moduli family."""

import sympy as sp

x, y, z, t, lam, w = sp.symbols("x y z t lambda w")
variables = (x, y, z)


def source_degree(vector):
    return max(sp.Poly(sp.cancel(entry), *variables).total_degree() for entry in vector)


# The degree-five family from verify_degree_five_stable_moduli.py.
H = sp.factor(
    w**2 * (w - 1) * (3 * w**2 - (5 * lam + 1) * w + 3 * lam) / 60
)
p = sp.diff(H, w)
c = sp.factor(-p.subs(w, 1))
q = sp.factor((w * p - H) / c)
u = 1 + x * y
gamma = 1 - sp.Rational(8, 7) * x * y + x**2 * z
W = sp.expand(u * gamma)

C_map = sp.expand(x * gamma)
B_map = sp.cancel((c + p.subs(w, W) / gamma) / x)
A_map = sp.cancel((u + q.subs(w, W) / gamma**2) / x**2)
F_lam = sp.Matrix((A_map, B_map, C_map))

# Center at lambda=2.  Since det(DF_lambda)=(lambda-1)/30, rescaling the
# first target coordinate by 1/(1+t) makes G_t fixed-Jacobian and based at F_2.
F0 = F_lam.subs(lam, 2).applyfunc(sp.cancel)
G_t = F_lam.subs(lam, 2 + t)
G_t[0] = G_t[0] / (1 + t)
H1 = G_t.diff(t).subs(t, 0).applyfunc(sp.cancel)
H2 = (G_t.diff(t, 2).subs(t, 0) / 2).applyfunc(sp.cancel)

DF0 = F0.jacobian(variables)
assert sp.factor(DF0.det()) == sp.Rational(1, 30)
DF0_inverse = (30 * DF0.adjugate()).applyfunc(sp.cancel)

# Canonical coefficients in G_t=F_2(id+t V1+t^2 V2+...).
V1 = (DF0_inverse * H1).applyfunc(sp.factor)
hessian_term = sp.Matrix(tuple(
    sp.expand((V1.T * sp.hessian(component, variables) * V1)[0] / 2)
    for component in F0
))
V2 = (DF0_inverse * (H2 - hessian_term)).applyfunc(sp.factor)

# The inverse jet is id-t V1+t^2((DV1)V1-V2).
inverse_1 = -V1
inverse_2 = (V1.jacobian(variables) * V1 - V2).applyfunc(sp.factor)

forward_degrees = (1, source_degree(V1), source_degree(V2))
inverse_degrees = (1, source_degree(inverse_1), source_degree(inverse_2))

assert forward_degrees == (1, 35, 69)
assert inverse_degrees == (1, 35, 69)

# Fixed determinant forces the two special-jet identities.  The second is the
# t^2 coefficient of det(I+t DV1+t^2 DV2)=1.
DV1 = V1.jacobian(variables)
DV2 = V2.jacobian(variables)
divergence_1 = sp.expand(sp.trace(DV1))
determinant_coefficient_2 = sp.factor(
    sp.trace(DV2)
    + (sp.trace(DV1) ** 2 - sp.trace(DV1 * DV1)) / 2
)
assert divergence_1 == 0
assert determinant_coefficient_2 == 0

print("PASS: target normalization makes the degree-five arc fixed-Jacobian")
print("PASS: canonical source and inverse coefficient degrees are 1,35,69")
print("PASS: the first- and second-order special-Jacobian identities hold")
