#!/usr/bin/env sage -python
"""Certify the non-Cartier E7 module needed by the first H92 q=6 pencil.

The source divisor is ``O + (-P1) - F``.  At its III* fibre, the scalar
complete ideal used by the older chord experiment has boundary ``e1+2e5``.
That is the Cartier double of the required minuscule E7 class.  The local
module for the actual divisor is instead

    Z * J_-P1^dual = Z * J_P1/(U-p_U),

    J_P1  = (U-p_U(Z), Y-p_Y(Z)),
    J_-P1 = (U-p_U(Z), Y+p_Y(Z)).

Here ``(Z,U,Y)`` is the formal E7 chart and the prefactor ``Z`` is the
``-F`` term.  The two signs have the same tangent E7 prime, but differ by
their order-three branch jets.  In the standard E7 ring
``Y^2=U^3+U Z^3``, the common tangent module is represented by the height-one
prime ``J=(U,Y)`` and

    J^dual = R + R*(Y/U).

The exact chord numerator vanishes on ``P1`` and its denominator vanishes on
both signs.  The identity

    (Y-p_Y)(Y+p_Y) = (U-p_U) H

therefore identifies ``(Y-p_Y)/(U-p_U)`` as a generator of
``J_-P1^dual``.  Thus a cleared chord numerator N must be tested by
module membership

    N = Z * ((U-p_U) A + (Y-p_Y) B),

not merely by exceptional valuations.  This script certifies the exact P1
jet which makes that translation unavoidable.  It deliberately does not
claim to have solved the resulting global Riemann--Roch system.
"""

import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, PowerSeriesRing, QQ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def series_at_infinity(value, source_ring, target_ring):
    """Expand a rational function of u at Z=1/u."""
    z = target_ring.gen()
    numerator = source_ring(value.numerator())
    denominator = source_ring(value.denominator())
    reversed_numerator = sum(
        target_ring(numerator[index]) * z ** (numerator.degree() - index)
        for index in range(numerator.degree() + 1)
    )
    reversed_denominator = sum(
        target_ring(denominator[index]) * z ** (denominator.degree() - index)
        for index in range(denominator.degree() + 1)
    )
    return z ** (denominator.degree() - numerator.degree()) * (
        reversed_numerator / reversed_denominator
    )


assert digest(SECTION) == SECTION_SHA256
section = json.loads(SECTION.read_text())
anchor = SourceFileLoader("h92_q6_p1_module_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
x_p1 = u_field(polynomial(u_ring, section["x_entrance_base"]["numerator_coefficients"])) / u_field(
    polynomial(u_ring, section["x_entrance_base"]["denominator_coefficients"])
)
y_p1 = u_field(polynomial(u_ring, section["y_entrance_base"]["numerator_coefficients"])) / u_field(
    polynomial(u_ring, section["y_entrance_base"]["denominator_coefficients"])
)
a_u = A1 / u**3 + A / u**4
b_u = B1 / u**5 + B / u**6 + B2 / u**7
assert y_p1**2 == x_p1**3 + a_u * x_p1 + b_u
weierstrass_ring = PolynomialRing(u_field, names=("X", "Y"))
X, Y = weierstrass_ring.gens()
weierstrass_relation = Y**2 - X**3 - a_u * X - b_u
branch_factorization = (
    (Y - y_p1) * (Y + y_p1)
    - (X - x_p1) * (X**2 + X * x_p1 + x_p1**2 + a_u)
)
assert branch_factorization == weierstrass_relation

series_ring = PowerSeriesRing(QQ, "Z", default_prec=10)
Z = series_ring.gen()
x_series = series_at_infinity(x_p1, u_ring, series_ring)
y_series = series_at_infinity(y_p1, u_ring, series_ring)

# The formal y=0 branch is the E7 chart origin used by the existing local
# flag.  Solve its cubic recursively; division by A1 is valid at H92.
x_zero = series_ring.zero()
for degree in range(2, 8):
    residual = (
        x_zero**3
        + (A1 * Z**3 + A * Z**4) * x_zero
        + B1 * Z**5
        + B * Z**6
        + B2 * Z**7
    )
    x_zero += -residual[degree + 3] / A1 * Z**degree

p_u = x_series - x_zero
p_y = y_series
assert x_series.valuation() == 2
assert y_series.valuation() == 3
assert p_u.valuation() == 3
assert p_y.valuation() == 3

# Algebraic model of the completed E7 singularity.  J is the nontrivial
# order-two class: U has valuation two along J, whereas Y has valuation one.
local_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
z_local, u_local, y_local = local_ring.gens()
relation = y_local**2 - u_local**3 - u_local * z_local**3
branch_ideal = local_ring.ideal((relation, u_local, y_local))
assert branch_ideal.is_prime()
local_fraction_field = local_ring.fraction_field()
eta = local_fraction_field(y_local) / local_fraction_field(u_local)
assert eta * u_local == y_local
dual_relation = local_fraction_field(eta * y_local - u_local**2 - z_local**3)
assert dual_relation.numerator() % relation == 0

print(
    "H92Q6P1MODULE|normal_form=Y2-U3-UZ3|"
    "branch_class=nontrivial_E7_discriminant_Z2|"
    "module=Z*Jminus_dual=Z*Jplus/(U-pU)|"
    "cleared_ideal=Z*(U-pU,Y-pY)",
    flush=True,
)
print(
    "H92Q6P1MODULE|"
    f"valuation_x={x_series.valuation()}|valuation_y={y_series.valuation()}|"
    f"valuation_pU={p_u.valuation()}|valuation_pY={p_y.valuation()}|"
    "status=PASS_EXACT_P1_BRANCH_JET",
    flush=True,
)
