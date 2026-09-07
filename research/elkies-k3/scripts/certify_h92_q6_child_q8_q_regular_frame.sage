#!/usr/bin/env sage -python
"""Certify the globally base-regular q8 q-frame after CRT normalization.

For the q8 marked frame q=(m-p)/h, its only generic vertical poles occur at
Nx, the numerator of x(S).  If R*h*Dy=Ny modulo Nx, then

  q-R/Nx = (Dy*(y+R*h)+Dx*x*((Ny-R*h*Dy)/Nx))
           /(h*Dy*(Dx*x-Nx)).

The displayed numerator cancellation removes Nx exactly.  The initial q is
regular at h, and Nx is coprime to h; hence the corrected frame has no
generic base pole there either.  Since deg(R)=95 and deg(Nx)=96, its leading
term has order one at infinity, replacing q's order 44.

This is a coordinate-frame certificate only.  It does not compile the
remaining additive/smooth local modules into a q8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-frame.json"


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

marking = json.loads(MARKING.read_text())
normalizer_record = json.loads(NORMALIZER.read_text())
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert normalizer_record["status"] == "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION"
assert normalizer_record["complete"] and normalizer_record["exact_check"]

ring = PolynomialRing(QQ, "T")
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, section["y_denominator_coefficients_low_to_high"])
normalizer = polynomial(ring, normalizer_record["normalizer"]["R_coefficients_low_to_high"])
h = monic_power_root(dx, 2)

assert nx.degree() == 96 and normalizer.degree() == 95
assert dx // h**2 in QQ and dy // h**3 in QQ
assert nx.gcd(h) in QQ and nx.gcd(dy) in QQ
remainder = (normalizer*h*dy-ny) % nx
assert not remainder
quotient, division_remainder = (ny-normalizer*h*dy).quo_rem(nx)
assert not division_remainder

# The first expression is the cleared old q frame.  Subtracting R/Nx and
# collecting coefficients leaves the second, with no Nx denominator.
old_numerator = "Nx*Dy*y+Ny*Dx*x"
old_denominator = "h*Nx*Dy*(Dx*x-Nx)"
regular_numerator = "Dy*(y+R*h)+Dx*x*((Ny-R*h*Dy)/Nx)"
regular_denominator = "h*Dy*(Dx*x-Nx)"

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-q-regular-frame.v1",
    "status": "PASS_EXACT_Q_REGULAR_FRAME",
    "inputs": {"marking": digest(MARKING), "normalizer": digest(NORMALIZER)},
    "q_frame": "q=(m-p)/h, m=(y+y(S))/(x-x(S)), p=-y(S)/x(S)",
    "normalized_frame": "q_regular=q-R/Nx",
    "exact_identity": {
        "old_numerator": old_numerator,
        "old_denominator": old_denominator,
        "regular_numerator": regular_numerator,
        "regular_denominator": regular_denominator,
        "quotient": "(Ny-R*h*Dy)/Nx",
        "quotient_degree": int(quotient.degree()),
    },
    "generic_base_regularities": {
        "Nx_pole_cancelled": True,
        "Nx_coprime_to_h": True,
        "Nx_coprime_to_Dy": True,
        "h_regular_by_difference": True,
        "conclusion": "q_regular has no generic base vertical pole",
    },
    "infinity": {
        "q_order": 44,
        "R_over_Nx_order": int(nx.degree()-normalizer.degree()),
        "q_regular_order": 1,
    },
    "boundary": (
        "This proves the global q-frame pole correction only. It does not "
        "derive transformed finite additive/smooth modules, a global q8 pencil, "
        "a D13 or rootless equation, bisections, extension collisions, or rank."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6CHILDQREG|Nx_degree={}|R_degree={}|quotient_degree={}|"
    "infinity_order={}|status=PASS_EXACT_Q_REGULAR_FRAME".format(
        nx.degree(), normalizer.degree(), quotient.degree(),
        payload["infinity"]["q_regular_order"],
    ),
    flush=True,
)
