#!/usr/bin/env python3
"""Bounded exact regression for formal source-orbit triviality.

The theorem is formal and does not depend on this example.  This script uses
the foundational noninjective Keller map, composes it with a determinant-one
polynomial curve, and recovers the curve uniquely from the deformed map
through order three using only (DF)^(-1).
"""

import sympy as sp

x, y, z, t = sp.symbols("x y z t")
variables = (x, y, z)


def truncate(expression, order):
    """Expand and discard powers t^order and above."""
    polynomial = sp.Poly(sp.expand(expression), t)
    return sp.expand(sum(
        polynomial.coeff_monomial(t**degree) * t**degree
        for degree in range(order)
    ))


def truncate_vector(vector, order):
    return sp.Matrix(tuple(truncate(entry, order) for entry in vector))


one_plus_xy = 1 + x * y
F = sp.Matrix((
    one_plus_xy**3 * z + y**2 * one_plus_xy * (4 + 3 * x * y),
    y + 3 * x * one_plus_xy**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
))
DF = F.jacobian(variables)
assert DF.det() == -2
DF_inverse = DF.adjugate() / -2


def compose(map_vector, substitution_vector):
    substitution = dict(zip(variables, substitution_vector))
    return sp.Matrix(tuple(
        sp.expand(entry.subs(substitution, simultaneous=True))
        for entry in map_vector
    ))


# A reduced curve made from elementary determinant-one shears.  Its first
# velocity is a sum of LNDs; the t^2 and t^3 factors model higher corrections.
alpha = sp.Matrix((x, y, z))
alpha = sp.Matrix((alpha[0] + t * alpha[1] ** 2, alpha[1], alpha[2]))
alpha = sp.Matrix((alpha[0], alpha[1] + t * alpha[2], alpha[2]))
alpha = sp.Matrix((alpha[0], alpha[1], alpha[2] + t**2 * alpha[0]))
alpha = sp.Matrix((alpha[0], alpha[1] + t**3 * alpha[0] ** 2, alpha[2]))
assert sp.expand(alpha.jacobian(variables).det()) == 1

order = 4
alpha_jet = truncate_vector(alpha, order)
family_jet = truncate_vector(compose(F, alpha), order)

# Recover alpha coefficient by coefficient.  At order r, the new unknown
# enters F(alpha) linearly as DF*V_r.
recovered = sp.Matrix((x, y, z))
for r in range(1, order):
    known_family = truncate_vector(compose(F, recovered), r + 1)
    error_coefficient = sp.Matrix(tuple(
        sp.expand(family_jet[index] - known_family[index]).coeff(t, r)
        for index in range(3)
    ))
    correction = (DF_inverse * error_coefficient).applyfunc(sp.expand)
    recovered = truncate_vector(recovered + t**r * correction, order)

assert all(
    sp.expand(left - right) == 0
    for left, right in zip(recovered, alpha_jet)
)
assert all(
    sp.expand(entry) == 0
    for entry in truncate_vector(compose(F, recovered) - family_jet, order)
)
assert truncate(recovered.jacobian(variables).det() - 1, order) == 0

print("PASS: the foundational Keller Jacobian has polynomial inverse matrix")
print("PASS: the canonical source trivializer is recovered uniquely through order three")
print("PASS: the recovered fixed-Jacobian jet has determinant one")
