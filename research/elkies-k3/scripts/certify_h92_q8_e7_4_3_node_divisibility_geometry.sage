#!/usr/bin/env sage -python
"""Separate the principal node divisibility condition from a finite corner jet.

At E7_4--E7_3 the actual completed local ring is QQ[[Z,Y]], with
U=Y^2*unit and t=Z^3*U^2.  For the base q8 clearing exponent T=17,
t^T is a unit times Z^51*Y^68.  Its quotient is one-dimensional, so it is
not itself a finite coefficient quotient.  The finite rectangle
R/(Z^51,Y^68) has length 3468, but it tests a corner jet and must not be
substituted for product divisibility by t^T.

This is an exact safety certificate for the forthcoming node quotient.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
CLEARING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-clearing.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-divisibility-geometry.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--clearing", type=Path, default=CLEARING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

clearing = json.loads(args.clearing.read_text())
assert clearing["status"] == "PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_CLEARING"
assert clearing["chart"]["name"] == "E7_4--E7_3"
T = int(clearing["common_clearing"]["T"])
assert T == 17
z_order, y_order = 3*T, 4*T
assert (z_order, y_order) == (51, 68)

# The displayed H92 equation has a nonzero U derivative at the node.  Hence
# its completion is a regular two-variable ring in (Z,Y), and U=Y^2*unit.
ring = PolynomialRing(QQ, names=("Z", "U", "Y"), order="degrevlex")
Z, U, Y = ring.gens()
surface = ring(sage_eval(clearing["chart"]["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
assert surface.derivative(U)(0, 0, 0) == -1
unit_h, remainder = (Y**2-surface).quo_rem(U)
assert not remainder and unit_h(0, 0, 0) == 1

# A finite-field dimension check records the actual chart relation in the
# rectangular jet.  The monic implicit equation eliminates U, leaving the
# 51*68 standard Z,Y monomials.
finite = GF(43)
finite_ring = PolynomialRing(finite, names=("Z", "U", "Y"), order="degrevlex")
z, u, y = finite_ring.gens()
finite_surface = finite_ring(surface)
corner_ideal = finite_ring.ideal((finite_surface, z**z_order, y**y_order))
corner_length = corner_ideal.vector_space_dimension()
assert corner_length == z_order*y_order == 3468

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-3-node-divisibility-geometry.v1",
    "status": "PASS_EXACT_Q8_E7_4_3_NODE_DIVISIBILITY_GEOMETRY",
    "inputs": {"principal_clearing": {"path": str(args.clearing.relative_to(ROOT)), "sha256": digest(args.clearing)}},
    "completed_chart": {
        "parameters": ["Z", "Y"],
        "relation": "U=Y^2*unit",
        "t": "Z^3*Y^4*unit",
        "t_power": T,
        "product_divisibility": "t^T=Z^{}*Y^{}*unit".format(z_order, y_order),
    },
    "quotient_geometry": {
        "R_mod_tT_dimension": 1,
        "finite_corner_ideal": "(Z^{},Y^{})".format(z_order, y_order),
        "finite_corner_length": int(corner_length),
        "corner_dimension_certificate_prime": 43,
        "warning": "R/(t^T) is not the finite corner quotient; product divisibility cannot be replaced by rectangular-jet vanishing.",
    },
    "compiler_consequence": (
        "A finite node block must be formed only after the componentwise "
        "conditions and overlap/module data identify a finite residual corner."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E743DIVISIBILITYGEOMETRY|T={}|product=Z{}Y{}|corner_length={}|"
    "status=PASS_EXACT_Q8_E7_4_3_NODE_DIVISIBILITY_GEOMETRY".format(
        T, z_order, y_order, corner_length,
    ),
    flush=True,
)
