#!/usr/bin/env python3
"""Exact bounded contact-profile audit for the degree-five moduli family."""

import sympy as sp

x, y, z, t, lam = sp.symbols("x y z t lambda")
variables = (x, y, z)
order = 3


def truncate(expression, cutoff=order):
    expression = sp.cancel(expression)
    series = sp.series(expression, t, 0, cutoff).removeO()
    return sp.cancel(sp.expand(series))


def truncate_vector(vector, cutoff=order):
    return sp.Matrix(tuple(truncate(entry, cutoff) for entry in vector))


def compose(map_vector, substitution_vector, cutoff=order):
    substitution = dict(zip(variables, substitution_vector))
    return truncate_vector(sp.Matrix(tuple(
        entry.subs(substitution, simultaneous=True) for entry in map_vector
    )), cutoff)


def coefficient(vector, exponent):
    return sp.Matrix(tuple(
        sp.expand(entry).coeff(t, exponent) for entry in vector
    ))


def source_degree(vector):
    return max(sp.Poly(sp.cancel(entry), *variables).total_degree() for entry in vector)


# The degree-five family from verify_degree_five_stable_moduli.py.
w = sp.symbols("w")
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

# Center at lambda=2.  The determinant varies as (lambda-1)/30, so rescale
# the first target coordinate by c(2)/c(2+t)=1/(1+t).  This is a based target
# automorphism over the completed local parameter ring and makes the family
# fixed-Jacobian before applying the canonical source trivialization.
F0 = F_lam.subs(lam, 2).applyfunc(sp.cancel)
family = F_lam.subs(lam, 2 + t)
family[0] = family[0] / (1 + t)
family = truncate_vector(family)

DF0 = F0.jacobian(variables)
assert sp.factor(DF0.det()) == sp.Rational(1, 30)
DF0_inverse = (30 * DF0.adjugate()).applyfunc(sp.cancel)

# Recover the unique canonical source jet alpha with family=F0(alpha).
identity = sp.Matrix(variables)
alpha = identity
alpha_coefficients = [identity]
for exponent in range(1, order):
    known = compose(F0, alpha, exponent + 1)
    error = coefficient(family - known, exponent)
    correction = (DF0_inverse * error).applyfunc(sp.cancel)
    alpha_coefficients.append(correction)
    alpha = truncate_vector(alpha + t**exponent * correction)

assert all(entry == 0 for entry in compose(F0, alpha) - family)
assert truncate(alpha.jacobian(variables).det() - 1) == 0

# Recover the inverse jet beta from alpha(beta)=id.  Its coefficient at each
# new order enters with the identity matrix because alpha_0=id.
beta = identity
beta_coefficients = [identity]
for exponent in range(1, order):
    known = compose(alpha, beta, exponent + 1)
    correction = -coefficient(known - identity, exponent)
    beta_coefficients.append(correction)
    beta = truncate_vector(beta + t**exponent * correction)

assert all(entry == 0 for entry in compose(alpha, beta) - identity)
assert all(entry == 0 for entry in compose(beta, alpha) - identity)

forward_degrees = tuple(source_degree(value) for value in alpha_coefficients)
inverse_degrees = tuple(source_degree(value) for value in beta_coefficients)
cumulative_lower_bounds = tuple(
    max(forward_degrees[: exponent + 1] + inverse_degrees[: exponent + 1])
    for exponent in range(order)
)

assert forward_degrees == (1, 14, 27)
assert inverse_degrees == (1, 14, 27)
assert cumulative_lower_bounds == (1, 14, 27)

print("PASS: target normalization makes the degree-five arc fixed-Jacobian")
print("PASS: recovered the unique source trivializer and its inverse through order two")
print("PASS: forward coefficient degrees are 1,14,27")
print("PASS: inverse coefficient degrees are 1,14,27")
