#!/usr/bin/env sage -python
"""Screen the q-frame principal-part normalization modulo a good prime.

The exact q-pole profile has pole divisor Nx/gcd(Nx,Ny).  At a prime where
this is Nx and all displayed factors remain coprime, the base Laurent
principal part is removed by

    R = Ny*(h*Dy)^(-1) mod Nx,   q_regular=q-R/Nx.

This fast reduction records the degree of R and hence the infinity order of
the normalized frame.  It is a modular shape check for a later CRT/exact
normalization, not that normalization itself.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("prime divides an input denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, coefficients):
    return ring([coefficient(field, value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, finite, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, finite, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, finite, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, finite, section["y_denominator_coefficients_low_to_high"])
h = monic_power_root(dx, 2)
delta = polynomial(ring, finite, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"])
assert nx.degree() == 96 and ny.degree() == 144 and h.degree() == 46
assert nx.gcd(ny).degree() == 0
assert nx.gcd(h).degree() == 0
assert nx.gcd(delta).degree() == 0

normalizer_numerator = (ny * (h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer_numerator*h*dy-ny) % nx == 0
normalizer_order = nx.degree()-normalizer_numerator.degree()
q_order = 44
assert normalizer_order > 0
print(
    "H92Q6CHILDQPOLEMODP|prime={}|Nx_degree={}|R_degree={}|"
    "R_infinity_order={}|q_regular_infinity_order={}|"
    "status=PASS_MODULAR_NORMALIZATION_SHAPE".format(
        args.prime, nx.degree(), normalizer_numerator.degree(), normalizer_order,
        min(q_order, normalizer_order),
    ),
    flush=True,
)
