#!/usr/bin/env python3
"""Exact degree-six exceptional and Ritt geometry in Gaussian moments."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import contact_partition_incidence


def primitive_coefficients(data):
    polynomial = sp.Poly(data.primitive, data.variable)
    return {
        degree: sp.cancel(polynomial.coeff_monomial(data.variable**degree))
        for degree in range(2, data.degree + 1)
    }


def normalized_moment(H, variable, order):
    """Return mu_m=M_m/m!=(1/m)[z^(m-1)](1+H)^m."""
    coefficient = sp.Poly((1 + H) ** order, variable).coeff_monomial(
        variable ** (order - 1)
    )
    return sp.cancel(coefficient / order)


u, v, m = sp.symbols("u v m")  # mu_3, mu_4, mu_5

# Displayed equation of the all-double component in normalized moments.
moment_surface = sp.expand(
    -48 * m**2 * u**2
    - 32 * m**2 * u * v
    + 16 * m**2
    + 192 * m * u**4
    + 128 * m * u**3 * v
    - 320 * m * u**3
    - 384 * m * u**2 * v
    + 128 * m * u**2
    - 104 * m * u * v**2
    + 112 * m * u * v
    + 24 * m * u
    + 8 * m * v**3
    + 40 * m * v
    - 16 * m
    - 192 * u**6
    - 128 * u**5 * v
    + 640 * u**5
    + 768 * u**4 * v
    - 768 * u**4
    + 208 * u**3 * v**2
    - 992 * u**3 * v
    + 400 * u**3
    - 16 * u**2 * v**3
    - 384 * u**2 * v**2
    + 400 * u**2 * v
    - 64 * u**2
    - 16 * u * v**3
    + 80 * u * v**2
    + 16 * u * v
    - 16 * u
    + 21 * v**4
    - 28 * v**3
    + 46 * v**2
    - 28 * v
    + 5
)

# Re-derive the known quartic seed-coordinate equation for C_(3,0), then
# transport it through the triangular moment coordinates.
all_double = contact_partition_incidence(6, (2, 2, 2))
W = all_double.variable
x, y, z = all_double.coefficient_parameters  # h_3,h_4,h_5
h6 = sp.Poly(all_double.normalized_universal_primitive, W).coeff_monomial(W**6)
e1, e2, e3 = all_double.quotient_coordinates
top3, top4, top5 = x / h6, y / h6, z / h6
inverse_e1 = -top5 / 2
inverse_e2 = (top4 - inverse_e1**2) / 2
inverse_e3 = -(top3 + 2 * inverse_e1 * inverse_e2) / 2
seed_surface = sp.factor(
    sp.cancel(
        all_double.phi.subs(
            {e1: inverse_e1, e2: inverse_e2, e3: inverse_e3}
        )
    ).as_numer_denom()[0]
)

# On the normalized sextic slice, u=c_2, v=c_3, and m=c_4+2c_2^2.
seed_from_moments = {
    x: v,
    y: m - 2 * u**2,
    z: 1 - 3 * v - 2 * m + 4 * u**2 - 4 * u,
}
transported = sp.expand(seed_surface.subs(seed_from_moments))
content, transported_primitive = sp.primitive(transported, u, v, m)
assert content == 4
assert sp.expand(transported_primitive - moment_surface) == 0
assert sp.Poly(moment_surface, u, v, m).total_degree() == 6
assert len(sp.factor_list(moment_surface)[1]) == 1

# Transport the two vertical Ritt hypersurfaces through the same triangular
# inverse.  The 2 o 3 equation coincides exactly with the all-double
# exceptional surface; the 3 o 2 equation gives the second displayed sextic.
h3_mom = v
h4_mom = m - 2 * u**2
h5_mom = 1 - 4 * u - 3 * v - 2 * m + 4 * u**2
h6_mom = m - 2 * u**2 + 3 * u + 2 * v - 1

ritt_32_seed = (
    27 * h3_mom * h6_mom**2
    - 18 * h4_mom * h5_mom * h6_mom
    + 5 * h5_mom**3
)
ritt_23_seed = (
    32 * h3_mom * h5_mom * h6_mom**2
    + 64 * h3_mom * h6_mom**3
    + 16 * h4_mom**2 * h6_mom**2
    - 24 * h4_mom * h5_mom**2 * h6_mom
    + 64 * h4_mom * h6_mom**3
    + 5 * h5_mom**4
    + 64 * h5_mom * h6_mom**3
    + 64 * h6_mom**4
)
assert sp.expand(ritt_23_seed - moment_surface) == 0

ritt_32_moment = sp.expand(
    -4 * m**3
    + 24 * m**2 * u**2
    - 60 * m**2 * u
    - 27 * m**2 * v
    + 6 * m**2
    - 48 * m * u**4
    + 240 * m * u**3
    + 108 * m * u**2 * v
    - 288 * m * u**2
    - 252 * m * u * v
    + 114 * m * u
    - 54 * m * v**2
    + 36 * m * v
    - 12 * m
    + 32 * u**6
    - 240 * u**5
    - 108 * u**4 * v
    + 552 * u**4
    + 504 * u**3 * v
    - 548 * u**3
    + 108 * u**2 * v**2
    - 549 * u**2 * v
    + 264 * u**2
    - 216 * u * v**2
    + 198 * u * v
    - 60 * u
    - 27 * v**3
    + 27 * v**2
    - 18 * v
    + 5
)
assert sp.expand(ritt_32_seed - ritt_32_moment) == 0
assert sp.Poly(ritt_32_moment, u, v, m).total_degree() == 6

clean_ritt_witness = {u: -sp.Rational(1, 4), v: 4, m: -sp.Rational(99, 8)}
assert moment_surface.subs(clean_ritt_witness) == 0
assert ritt_32_moment.subs(clean_ritt_witness) == 0

# Parametrize the all-triple component C_(0,2), compute its first three
# normalized Gaussian moments directly, and compare with the closed formulas.
all_triple = contact_partition_incidence(6, (3, 3))
q1, q2 = all_triple.quotient_coordinates
q = sp.symbols("q")
root_parameterization = {
    q1: 1 - 3 * q**2 / (3 * q - 1),
    q2: -3 * q**3 / (3 * q - 1),
}
assert sp.factor(all_triple.phi.subs(root_parameterization)) == 0

H_triple = sp.cancel(all_triple.primitive.subs(root_parameterization))
actual_moments = tuple(
    sp.factor(normalized_moment(H_triple, all_triple.variable, order))
    for order in (3, 4, 5)
)
A = 15 * q**3 - 15 * q**2 + 6 * q - 1
B = 5 * q**2 - 5 * q + 1
expected_moments = (
    A / (3 * q * B),
    -(
        (3 * q**2 - 3 * q + 1)
        * (45 * q**4 - 15 * q**2 + 6 * q - 1)
    )
    / (27 * q**4 * B),
    A
    * (30 * q**5 - 30 * q**4 - 3 * q**3 + 18 * q**2 - 8 * q + 1)
    / (9 * q**4 * B**2),
)
assert all(
    sp.factor(actual - expected) == 0
    for actual, expected in zip(actual_moments, expected_moments)
)

# The surface meets the curve exactly at the all-six collision divisor, and
# the square proves local intersection length two at every one of its four
# simple points on this normalization chart.
collision = 45 * q**4 - 30 * q**3 + 15 * q**2 - 6 * q + 1
pullback = sp.factor(
    moment_surface.subs(dict(zip((u, v, m), expected_moments)))
)
expected_pullback = sp.factor(
    -(3 * q - 1) ** 8 * collision**2 / (3**11 * q**16 * B**4)
)
assert sp.factor(pullback - expected_pullback) == 0
assert sp.gcd(collision, sp.diff(collision, q)) == 1
assert sp.gcd(collision, q * (3 * q - 1) * B) == 1
assert sp.degree(collision, q) == 4

# At a collision the quadratic root is r=q1/2.  The all-six incidence equation
# contains exactly the collision factor, and its own weighted admissibility
# remains a unit at all four roots.
all_six = contact_partition_incidence(6, (6,))
r = all_six.quotient_coordinates[0]
collision_root = sp.factor(root_parameterization[q1] / 2)
all_six_phi_pullback = sp.factor(all_six.phi.subs(r, collision_root))
assert sp.rem(
    sp.Poly(sp.fraction(all_six_phi_pullback)[0], q), sp.Poly(collision, q)
) == 0
admissibility_numerator = sp.Poly(
    sp.fraction(
        sp.factor(all_six.weighted_admissibility_factor.subs(r, collision_root))
    )[0],
    q,
)
assert sp.gcd(admissibility_numerator, sp.Poly(collision, q)).degree() == 0

print("PASS degree-six moments: (mu_3,mu_4,mu_5) invert the normalized seed")
print("PASS degree-six moments: C_(3,0) is the displayed irreducible sextic surface")
print("PASS degree-six moments: C_(0,2) has the displayed rational parameterization")
print("PASS degree-six moments: the four all-six points have intersection length two")
print("PASS degree-six moments: both vertical Ritt loci have explicit equations")
print("PASS degree-six moments: the 2o3 Ritt locus equals C_(3,0)")
