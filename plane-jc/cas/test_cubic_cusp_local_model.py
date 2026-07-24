#!/usr/bin/env python3
"""Exact regression for the cubic cusp finite-flat local model."""

import sympy as sp


t, u, v, T = sp.symbols("t u v T")

# The source is k[t,u], finite free of rank three over k[u,v] because t is
# integral with monic equation T^3 + uT - v.
minimal_equation = T**3 + u * T - v
assert sp.Poly(minimal_equation, T).monic()
assert sp.Poly(minimal_equation, T).degree() == 3

source_v = t**3 + u * t
jacobian = sp.diff(source_v, t)
assert sp.expand(jacobian) == u + 3 * t**2

# The discriminant cusp and its full pullback.
discriminant = 4 * u**3 + 27 * v**2
pulled_discriminant = sp.factor(discriminant.subs(v, source_v))
assert pulled_discriminant == (u + 3 * t**2) ** 2 * (4 * u + 3 * t**2)

# The affine companion meets the ramification boundary with contact two.
boundary_equation = u + 3 * t**2
companion_equation = 4 * u + 3 * t**2
companion_on_boundary = sp.expand(companion_equation.subs(u, -3 * t**2))
assert companion_on_boundary == -9 * t**2
assert sp.Poly(companion_on_boundary, t).degree() == 2

# On the ramification curve u=-3t^2, the image is a cusp and the derivative
# of its normalization parametrization vanishes at the origin.
cusp_u = -3 * t**2
cusp_v = sp.expand(source_v.subs(u, cusp_u))
assert cusp_v == -2 * t**3
assert sp.expand(discriminant.subs({u: cusp_u, v: cusp_v})) == 0
assert sp.diff(cusp_u, t).subs(t, 0) == 0
assert sp.diff(cusp_v, t).subs(t, 0) == 0

# The special fiber over the cusp has coordinate algebra k[t]/(t^3), hence
# length three: it is flat, curvilinear, and not a packet-length defect.
special_fiber = sp.Poly(minimal_equation.subs({u: 0, v: 0}), T)
assert special_fiber.as_expr() == T**3
assert special_fiber.degree() == 3

print("cubic cusp local model: exact checks passed")
