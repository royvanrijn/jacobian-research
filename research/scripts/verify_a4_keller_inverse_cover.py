#!/usr/bin/env python3
"""Exact checks for the oriented quartic A4 Cox--Keller inverse cover."""

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_factor


T, p, q, r, D, z, w = sp.symbols("T p q r D z w")

P = T**4 + p * T**2 + q * T + r
dP = sp.diff(P, T)
Delta = (
    256 * r**3
    - 128 * p**2 * r**2
    + 144 * p * q**2 * r
    - 27 * q**4
    + 16 * p**4 * r
    - 4 * p**3 * q**2
)

assert sp.factor(sp.discriminant(P, T) - Delta) == 0
assert sp.factor(Delta) == Delta

# The derivative is a unit on the oriented discriminant open. Compute its
# literal inverse in Q(p,q,r)[T]/(P); its denominator is Delta=D^2.
K = sp.QQ.frac_field(p, q, r)
inverse_dP = sp.cancel(
    sp.invert(
        sp.Poly(dP, T, domain=K),
        sp.Poly(P, T, domain=K),
    ).as_expr()
)
inverse_num, inverse_den = sp.fraction(inverse_dP)
assert sp.factor(inverse_den - Delta) == 0
assert sp.rem(sp.expand(dP * inverse_num - Delta), P, T) == 0

# On the root-incidence chart, solve P=0 for r. The derivative-unit
# suspension (p,q,T,z) -> (p,q,r,z/P') has constant local Jacobian -1.
r_on_source = -T**4 - p * T**2 - q * T
Z_on_source = z / dP
local_map = sp.Matrix([p, q, r_on_source, Z_on_source])
local_jacobian = sp.factor(
    local_map.jacobian((p, q, T, z)).det()
)
assert local_jacobian == -1
assert sp.factor(
    inverse_dP.subs(r, r_on_source) - 1 / dP
) == 0

# The obvious polynomial pole clearing z=Delta*w replaces Z=z/P' by the
# polynomial w*inverse_num. Its Jacobian acquires the nonconstant factor
# -Delta. Thus one-coordinate clearing cannot yield an absolute Keller map.
cleared_Z = w * inverse_num.subs(r, r_on_source)
cleared_map = sp.Matrix([p, q, r_on_source, cleared_Z])
cleared_jacobian = sp.factor(
    cleared_map.jacobian((p, q, T, w)).det()
)
Delta_on_source = sp.factor(Delta.subs(r, r_on_source))
assert sp.factor(cleared_jacobian + Delta_on_source) == 0


def factor_degrees_mod_prime(coefficients, prime):
    """Return the squarefree factor degrees over F_prime."""

    _, factors = gf_factor(coefficients, prime, ZZ)
    return tuple(sorted(
        (len(factor) - 1 for factor, exponent in factors
         for _ in range(exponent)),
        reverse=True,
    ))


# A specialization of the generic depressed quartic has S4 group:
# x^4-x-1 is irreducible mod 2 (a 4-cycle), while mod 7 it has type (3,1).
# A transitive quartic group containing a 4-cycle and a 3-cycle is S4.
s4_coefficients = [1, 0, 0, -1, -1]
assert sp.discriminant(T**4 - T - 1, T) == -283
assert factor_degrees_mod_prime(s4_coefficients, 2) == (4,)
assert factor_degrees_mod_prime(s4_coefficients, 7) == (3, 1)

# The oriented base adjoins sqrt(Delta), the sign fixed field of the generic
# S4 closure. Hence its generic group is A4 in the natural four-point action.
# This concrete A4 fiber has square discriminant and irreducible cubic
# resolvent; its good reductions exhibit all three A4 cycle partitions.
p0, q0, r0 = -7, -3, 1
P0 = sp.Poly(P.subs({p: p0, q: q0, r: r0}), T, domain=sp.QQ)
Delta0 = int(Delta.subs({p: p0, q: q0, r: r0}))
resolvent0 = sp.Poly(
    T**3 - p0 * T**2 - 4 * r0 * T + (4 * p0 * r0 - q0**2),
    T,
    domain=sp.QQ,
)

assert P0.as_expr() == T**4 - 7 * T**2 - 3 * T + 1
assert Delta0 == 183**2
assert P0.is_irreducible
assert P0.count_roots(-sp.oo, sp.oo) == 4
assert resolvent0.as_expr() == T**3 + 7 * T**2 - 4 * T - 37
assert resolvent0.is_irreducible

fiber_coefficients = [1, 0, -7, -3, 1]
assert factor_degrees_mod_prime(fiber_coefficients, 2) == (3, 1)
assert factor_degrees_mod_prime(fiber_coefficients, 11) == (2, 2)
assert factor_degrees_mod_prime(fiber_coefficients, 233) == (1, 1, 1, 1)

print("PASS: Disc(T^4+pT^2+qT+r) equals the oriented equation D^2")
print("PASS: P'(T) has an exact quotient-ring inverse with denominator D^2")
print("PASS: the derivative-unit suspension has constant local Jacobian -1")
print("PASS: direct polynomial pole clearing has Jacobian -D^2, not a unit")
print("PASS: generic oriented quartic monodromy is A4 in degree four")
print("PASS: T^4-7T^2-3T+1 is a totally real A4 fiber with local types")
print("      (3,1) at 2, (2,2) at 11, and split at 233")
