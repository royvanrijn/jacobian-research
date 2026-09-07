#!/usr/bin/env sage -python
"""Derive the exact smooth P1.O principal-part profile for q=8.

At the four smooth collision fibres let ``h`` be the squarefree P1 collision
polynomial and put

    p=y(P1)/x(P1)=rho/h,     q=(m-p)/h,       X=h^2*x.

The q=6 smooth module certifies that ``r`` is a unit modulo ``h`` and uses
the saturated frame ``<1,q>``.  Every endpoint q=8 ambient element is

    u^i/h^k * x^a*m^b.

Substitution ``m=r/h+h*q`` gives the exact finite expansion

    u^i*X^a * sum_j binomial(b,j) rho^(b-j)
        h^(2*j-b-k-2*a) q^j.

This script records all exponent profiles.  In particular it proves that
only jets through h^-15 can occur in the 54-dimensional endpoint ambient.
It does not decide which principal parts are allowed by the q=8 divisor;
that requires the q8 smooth line-bundle frame and yields the future condition
matrix.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, binomial, gcd


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SMOOTH = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-smooth-po-module.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-principal-parts.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--smooth", type=Path, default=SMOOTH)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
smooth = json.loads(args.smooth.read_text())
ambient = json.loads(args.ambient.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert smooth["status"] == "PASS_EXACT_SMOOTH_PO_CHORD_MODULE"
assert ambient["status"] == "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT"
assert smooth["module"]["saturated_frame"] == "<1,(m-y(P1)/x(P1))/h>"
assert ambient["ambient_dimension"] == 54

ring = PolynomialRing(QQ, "u")
field = ring.fraction_field()
h = polynomial(ring, p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and gcd(h, h.derivative()) == 1
x_p = field(polynomial(ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= field(polynomial(ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = field(polynomial(ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= field(polynomial(ring, p1["y_entrance_base"]["denominator_coefficients"]))
p = y_p/x_p
p_numerator = ring(p.numerator())
p_denominator = ring(p.denominator())
denominator_unit, remainder = p_denominator.quo_rem(h)
assert not remainder
assert gcd(h, denominator_unit) == 1 and gcd(h, p_numerator) == 1
rho = field(p_numerator)/field(denominator_unit)
assert p == rho/field(h)

profiles = []
max_pole = 0
for entry in ambient["ambient_basis"]:
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    k = int(entry["h_power"])
    terms = []
    for j in range(b+1):
        h_exponent = 2*j-b-k-2*a
        terms.append({
            "q_power": j,
            "h_exponent": h_exponent,
            "coefficient": "binomial({},{} )*u^i*rho^{}".format(b, j, b-j),
            "coefficient_nonzero_mod_h": True,
        })
        assert binomial(b, j)
        max_pole = max(max_pole, -h_exponent)
    assert terms[0]["h_exponent"] == -b-k-2*a
    assert terms[-1]["h_exponent"] == b-k-2*a
    profiles.append({
        "ambient_basis": entry,
        "local_form": "u^i*X^{}*sum_j binomial({},j)*rho^({}-j)*h^(2j-{}-{}-2*{})*q^j".format(a, b, b, b, k, a),
        "q6_saturated_coordinate": "q=(m-y(P1)/x(P1))/h",
        "h_principal_part_order": b+k+2*a,
        "terms": terms,
    })

assert len(profiles) == 54
assert max_pole == 15
assert all(profile["h_principal_part_order"] <= 15 for profile in profiles)
worst_families = {
    (
        int(profile["ambient_basis"]["x_power"]),
        int(profile["ambient_basis"]["m_power"]),
        int(profile["ambient_basis"]["h_power"]),
    )
    for profile in profiles
    if profile["h_principal_part_order"] == max_pole
}
# The ordinary m^9 family is the unique order-15 extremum.  The x*m^7
# family has order 14 only after the X=h^2*x correction; this separate check
# prevents an accidental profile in the old x-coordinate, which would omit
# the necessary 2*x_power shift.
assert worst_families == {(0, 9, 6)}
assert {
    profile["h_principal_part_order"]
    for profile in profiles
    if int(profile["ambient_basis"]["x_power"]) == 1
    and int(profile["ambient_basis"]["m_power"]) == 7
} == {14}

payload = {
    "schema": "elkies-k3.h92-q8-smooth-collision-principal-parts.v1",
    "status": "PASS_EXACT_Q8_SMOOTH_COLLISION_PRINCIPAL_PARTS",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "q6_smooth_module": {"path": str(args.smooth.relative_to(ROOT)), "sha256": digest(args.smooth)},
        "q8_endpoint_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
    },
    "collision_frame": {
        "h": str(h),
        "degree": int(h.degree()),
        "squarefree": True,
        "p": "y(P1)/x(P1)=rho/h",
        "rho_unit_mod_h": True,
        "q": "(m-p)/h",
    },
    "substitution": "u^i/h^k*x^a*m^b = u^i*X^a*sum_j binomial(b,j)*rho^(b-j)*h^(2j-b-k-2a)*q^j",
    "principal_part_profiles": profiles,
    "maximum_h_pole_order": max_pole,
    "extremal_families": [
        {"x_power": a, "m_power": b, "h_power": k}
        for a, b, k in sorted(worst_families)
    ],
    "finite_jet_instruction": (
        "A q8 smooth collision condition need only inspect h-adic principal "
        "parts through h^-15 in the displayed q-power/x*q-power coordinates. "
        "Its allowed target submodule must still be derived from the q8 divisor."
    ),
    "boundary": (
        "This is an exact expansion profile, not a smooth collision condition "
        "matrix. It does not claim which h-principal parts are allowed, a "
        "complete resolved cover, h0=2, or a q8 pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8SMOOTHPROFILE|ambient=54|max_h_pole=15|q_frame=actual_q6_saturated|"
    "status=PASS_EXACT_Q8_SMOOTH_COLLISION_PRINCIPAL_PARTS",
    flush=True,
)
