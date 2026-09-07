#!/usr/bin/env sage
"""Test the section-projection construction of the fifth q=4 pencil.

On the compact fourth-child genus-one equation C_T(V,X)=0, the explicit
degree-one section R has new MW coordinates (1,-1).  Its inverse is the
horizontal target (-1,1) of the fifth neighbor.  The line pencil through R,

    u = (X-X_R(T))/(V-V_R(T)),

therefore has residual divisor H-R = O+(-R) whenever the compact cubic
embedding has H=3O.  This checker performs that substitution exactly over
GF(73)(T), removes the automatic V-V_R factor, and records the residual
equation and discriminant.  It deliberately avoids a global Jacobian or a
Riemann--Roch calculation on the old Weierstrass surface.
"""

import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
MOVING_ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fourth-q12-cm24-moving-cubic-gf73.json"
)
LOCAL_ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fifth-q4-local-module-ansatz-gf73.json"
)
moving_bytes = MOVING_ARTIFACT.read_bytes()
moving_hash = hashlib.sha256(moving_bytes).hexdigest()
assert moving_hash == (
    "c6560b3db2d1232866e9996fc727924090aa46293c2482885cf9f9dbf4c21c89"
)
moving_data = json.loads(moving_bytes)
local_bytes = LOCAL_ARTIFACT.read_bytes()
local_hash = hashlib.sha256(local_bytes).hexdigest()
assert local_hash == (
    "7cadf12e4035dc9325f3249158f906762aba85cebdb2e92cb4de93efc2140d15"
)
local_data = json.loads(local_bytes)

finite = GF(73, impl="modn")
parameter_field = FunctionField(finite, "t")
t = parameter_field.gen()
projection_ring = PolynomialRing(parameter_field, names=("v", "u"))
v, u = projection_ring.gens()
# Rebuild in a temporary three-variable ring so substitution is unambiguous.
compact_ring = PolynomialRing(parameter_field, names=("vv", "xx"))
vv, xx = compact_ring.gens()
compact = compact_ring(sum(
    parameter_field(coefficient)
    * t**t_degree * vv**v_degree * xx**x_degree
    for t_degree, v_degree, x_degree, coefficient
    in moving_data["moving_terms_T_v_x_coefficient"]
))

section_v = (34*t+30)/(t+4)
x_coefficients = tuple(local_data["explicit_section"]["x_coefficients"])
section_x = sum(
    finite(coefficient)*section_v**degree
    for degree, coefficient in enumerate(x_coefficients)
)
assert compact(vv=section_v, xx=section_x) == 0

substituted = projection_ring(
    compact(vv=v, xx=section_x+u*(v-section_v))
)
quotient, remainder = substituted.quo_rem(v-section_v)
assert remainder == 0
assert quotient != 0

# The compact equation is the birational image of a cubic embedding, not a
# literal total-degree-three plane equation.  Projection therefore retains
# fixed exceptional factors.  Remove the polynomial content common to all
# coefficients in u; this is the exact algebraic counterpart of stripping
# those fixed vertical/base-point components.
univariate_v = PolynomialRing(parameter_field, "z")
z = univariate_v.gen()
cubic_ring = PolynomialRing(univariate_v, "w")
w = cubic_ring.gen()
cubic = cubic_ring(compact(vv=z, xx=w))
cubic_discriminant_factors = tuple(cubic.discriminant().factor())
u_coefficients = []
for u_degree in range(quotient.degree(u)+1):
    coefficient = sum(
        value*z**v_degree
        for (v_degree, candidate_u_degree), value in quotient.dict().items()
        if candidate_u_degree == u_degree
    )
    u_coefficients.append(univariate_v(coefficient))
fixed_content = u_coefficients[0]
for coefficient in u_coefficients[1:]:
    fixed_content = fixed_content.gcd(coefficient)
fixed_content = fixed_content.monic()
fixed_content_bivariate = projection_ring(fixed_content(z=v))
moving_quotient, content_remainder = quotient.quo_rem(fixed_content_bivariate)
assert content_remainder == 0

degree_v = moving_quotient.degree(v)
degree_u = moving_quotient.degree(u)
quadratic_discriminant = None
if degree_v == 2:
    # Extract coefficients in v without relying on multivariate monomial
    # indexing conventions.
    a = sum(
        value*u**u_degree
        for (v_degree, u_degree), value in moving_quotient.dict().items()
        if v_degree == 2
    )
    b = sum(
        value*u**u_degree
        for (v_degree, u_degree), value in moving_quotient.dict().items()
        if v_degree == 1
    )
    c = sum(
        value*u**u_degree
        for (v_degree, u_degree), value in moving_quotient.dict().items()
        if v_degree == 0
    )
    quadratic_discriminant = b**2-4*a*c

print(
    "Q80FIFTHQ4PROJECTION|"
    f"compact_total_degree={compact.total_degree()}|"
    f"degree_v={degree_v}|degree_u={degree_u}|"
    f"raw_degree_v={quotient.degree(v)}|"
    f"fixed_content_degree_v={fixed_content.degree()}|"
    f"degree_v={degree_v}|degree_u={degree_u}|"
    f"term_count={len(moving_quotient.dict())}|"
    f"cubic_discriminant_factor_degrees_exponents={tuple((factor.degree(), exponent) for factor, exponent in cubic_discriminant_factors)}|"
    f"quadratic={int(degree_v == 2)}|"
    "status=PASS_SECTION_FACTOR",
    flush=True,
)
if quadratic_discriminant is not None:
    numerator = quadratic_discriminant.numerator()
    denominator = quadratic_discriminant.denominator()
    print(
        "Q80FIFTHQ4PROJECTION|"
        f"disc_num_degree_t={numerator.degree()}|"
        f"disc_den_degree_t={denominator.degree()}|"
        f"disc_num_factorization={numerator.factor()}|"
        f"disc_den_factorization={denominator.factor()}|"
        "status=PASS_QUADRATIC_RESIDUAL",
        flush=True,
    )
