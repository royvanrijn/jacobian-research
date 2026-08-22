#!/usr/bin/env sage -python
"""Derive the exact generic vertical-pole divisor of the q8 q-frame.

The finite q-frame ``q=(m-p)/h`` is regular at the smooth O.S collision
divisor ``h`` and is especially convenient at II* and IV*.  That local fact
does not make it a global coefficient generator.  If

    x(S)=Nx/Dx,  y(S)=Ny/Dy,

then clearing the chord denominator gives

    q = (Nx*Dy*y + Ny*Dx*x) / (h*Nx*Dy*(Dx*x-Nx)).

At a generic point over a factor of ``Nx`` away from ``h``, its valuation is
``min(ord(Nx),ord(Ny))-ord(Nx)``.  Hence its vertical pole divisor is exactly
``Nx/gcd(Nx,Ny)``.  This script records that divisor without factoring its
large rational polynomial and also records the order 44 zero of q at the
smooth fibre at infinity in the minimal Weierstrass scaling.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-profile.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

ring = PolynomialRing(QQ, "T")
T = ring.gen()
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = field(polynomial(ring, section["x_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["x_denominator_coefficients_low_to_high"])
)
sy = field(polynomial(ring, section["y_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["y_denominator_coefficients_low_to_high"])
)
nx, dx = ring(sx.numerator()), ring(sx.denominator())
ny, dy = ring(sy.numerator()), ring(sy.denominator())
h = monic_power_root(dx, 2)
assert h.degree() == 46
assert dx // h**2 in QQ and dy // h**3 in QQ
assert nx.gcd(h) in QQ and ny.gcd(h) in QQ

# q=(m-p)/h with m=(y+sy)/(x-sx) and p=-sy/sx.  The cleared expression is
# linear in the generic fibre coordinates x,y, so its coefficient valuation
# computes its valuation at every generic base divisor.
common = nx.gcd(ny)
pole_divisor, remainder = nx.quo_rem(common)
assert not remainder
assert pole_divisor.gcd(h) in QQ

delta = polynomial(ring, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"])
iv = polynomial(ring, PolynomialRing(QQ, "T")(
    next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")["factor"]
).list())
assert iv.degree() == 1 and pole_divisor.gcd(iv) in QQ

# With s=1/T, x=s^-4 X and y=s^-6 Y.  The marked section has orders
# x(S)=s^-4*unit and y(S)=s^-6*unit, while h=s^-46*unit.  Consequently
# m-p=s^-2*unit generically and q=s^44*unit.  The degree data below certifies
# those orders without choosing a point on the smooth fibre at infinity.
assert (nx.degree()-dx.degree(), ny.degree()-dy.degree()) == (4, 6)
q_infinity_order = h.degree()-2
assert q_infinity_order == 44

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-q-pole-profile.v1",
    "status": "PASS_EXACT_Q_FRAME_POLE_PROFILE",
    "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
    "q_frame": "q=(m-p)/h, m=(y+y(S))/(x-x(S)), p=-y(S)/x(S)",
    "cleared_expression": (
        "q=(Nx*Dy*y+Ny*Dx*x)/(h*Nx*Dy*(Dx*x-Nx))"
    ),
    "section_coordinate_degrees": {
        "Nx": int(nx.degree()), "Dx": int(dx.degree()),
        "Ny": int(ny.degree()), "Dy": int(dy.degree()), "h": int(h.degree()),
    },
    "generic_vertical_poles": {
        "divisor": "Nx/gcd(Nx,Ny)",
        "degree": int(pole_divisor.degree()),
        "coprime_to_h": True,
        "coprime_to_IV_star_factor": True,
        "gcd_degree_with_weierstrass_discriminant": int(pole_divisor.gcd(delta).degree()),
        "valuation_rule": "ord_f(q)=min(ord_f(Nx),ord_f(Ny))-ord_f(Nx)",
    },
    "infinity": {
        "minimal_weierstrass_scaling": "x=s^-4*X, y=s^-6*Y, s=1/T",
        "q_order": int(q_infinity_order),
        "derivation": "m-p has order -2 and h has order -46",
    },
    "regular_base_coefficient_consequence": {
        "finite_B_generator": "f_IV*",
        "uncancelled_q_pole_clearer": "f_IV*(Nx/gcd(Nx,Ny))",
        "uncancelled_q_pole_clearer_degree": int(iv.degree()+pole_divisor.degree()),
        "qualification": (
            "If C is regular at this divisor, B must contain this factor. "
            "A saturated global module may instead cancel q's base Laurent "
            "principal parts through matching poles in C."
        ),
    },
    "boundary": (
        "This is an exact generic-base pole profile. It proves that a global "
        "not derive the saturated infinity gluing, a two-dimensional q8 pencil, "
        "a rootless equation, bisections, extension collisions, or a rank claim."
        "q-frame coefficient pair must account for the displayed divisor, either "
        "by vanishing in B or by a matching base principal part in C. It does not "
        "derive the saturated infinity gluing, a two-dimensional q8 pencil, a "
        "rootless equation, bisections, extension collisions, or a rank claim."
        "not derive the saturated infinity gluing, a two-dimensional q8 pencil, "
        "a rootless equation, bisections, extension collisions, or a rank claim."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6CHILDQ8QPOLES|pole_degree={}|q_infinity_order={}|"
    "uncancelled_B_degree={}|disc_gcd_degree={}|status=PASS_EXACT_Q_FRAME_POLE_PROFILE".format(
        pole_divisor.degree(), q_infinity_order, iv.degree()+pole_divisor.degree(),
        pole_divisor.gcd(delta).degree(),
    ),
    flush=True,
)
