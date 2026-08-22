#!/usr/bin/env sage -python
"""Derive old-coordinate valuations on the actual resolved H92 E7 fibre.

This is the actual-chart counterpart of the former normal-form valuation
atlas.  Every order below is read from a certified H92 edge-chart pullback;
when a reduced strict transform has local equation ``Y=0`` we first use the
actual smooth surface equation to prove ``U`` or ``Z`` is a unit times
``Y^2``.  Thus no order is inferred from an E7/III* component label.

The resulting vectors are the raw input from which a compiler can construct
anti-nef complete ideals and finite quotient blocks for an integral vertical
correction.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import QQ, PolynomialRing, matrix, vector, ZZ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json"
RESOLUTION_SHA256 = "14378f4718d3fbe781d5b351ba4943a962a430d9827c0ce285cd1125a9e8c500"
PULLBACKS_SHA256 = "a8f19a1e205ed2250c83c312a4bf722462867d289ddd598987dd5112f1c33177"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.resolution == RESOLUTION:
    assert digest(args.resolution) == RESOLUTION_SHA256
if args.pullbacks == PULLBACKS:
    assert digest(args.pullbacks) == PULLBACKS_SHA256
resolution = json.loads(args.resolution.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"

anchor = SourceFileLoader("h92_actual_e7_valuation_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
f0 = Y**2-U**3-(A1*Z**3+A*Z**4)*U-(B1*Z**5+B*Z**6+B2*Z**7)
f1_z = strict(ring, f0, (Z, Z*U, Z*Y), Z)
f2_u = strict(ring, f1_z, (U*Z, U, U*Y), U)
f3_u = strict(ring, f2_u, (U*Z, U, U*Y), U)
f3_z = strict(ring, f2_u, (Z, Z*U, Z*Y), Z)
second = -A1/B1
third = -QQ(1)/A1
edge_equations = {
    "E7_2--E7_5": strict(ring, ring(f2_u(Z+second, U, Y)), (Z, Z*U, Z*Y), Z),
    "E7_1--E7_4": strict(ring, f3_u, (U*Z, U, U*Y), U),
    "E7_4--E7_3": strict(ring, f3_u, (Z, Z*U, Z*Y), Z),
    "E7_3--E7_7": strict(ring, f3_z, (U*Z, U, U*Y), U),
    "E7_7--E7_2": strict(ring, f3_z, (Z, Z*U, Z*Y), Z),
    "E7_3--E7_6": strict(ring, ring(f3_u(Z+third, U, Y)), (Z, Z*U, Z*Y), Z),
}
assert {name: str(value) for name, value in edge_equations.items()} == resolution["edge_charts"]

# The first four strict-transform entries require the displayed quadratic
# relation.  At their edge origins the indicated derivative is nonzero, so
# the implicit function theorem makes the listed coordinate a unit times Y^2.
assert edge_equations["E7_1--E7_4"].derivative(Z)(0, 0, 0)
assert edge_equations["E7_4--E7_3"].derivative(U)(0, 0, 0)
assert edge_equations["E7_2--E7_5"].derivative(U)(0, 0, 0)

# Each triple is the exact valuation of (t,x,y) at the generic point of the
# named reduced exceptional component.  The chart and parameter record the
# calculation, including the three cases where Y is the reduced parameter.
entries = (
    ("E7_1", "E7_1--E7_4", "Y; Z is a unit times Y^2", (2, 2, 3)),
    ("E7_2", "E7_2--E7_5", "Y; U is a unit times Y^2", (2, 4, 5)),
    ("E7_3", "E7_4--E7_3", "Y; U is a unit times Y^2", (4, 6, 9)),
    ("E7_4", "E7_1--E7_4", "U", (3, 4, 6)),
    ("E7_5", "E7_2--E7_5", "Z", (1, 2, 3)),
    ("E7_6", "E7_3--E7_6", "Z", (2, 3, 5)),
    ("E7_7", "E7_3--E7_7", "U", (3, 5, 7)),
)

base_orders = vector(ZZ, [entry[3][0] for entry in entries])
cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
assert cartan*base_orders == vector(ZZ, (1, 0, 0, 0, 0, 0, 0))
assert tuple(base_orders) == (2, 2, 4, 3, 1, 2, 3)

payload = {
    "schema": "elkies-k3.h92-q6-actual-e7-valuation-atlas.v1",
    "status": "PASS_EXACT_H92_E7_VALUATION_ATLAS",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "actual_resolution": {"path": str(args.resolution.relative_to(ROOT)), "sha256": digest(args.resolution)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
    },
    "entries": [
        {"component": name, "edge_chart": chart, "reduced_parameter": parameter,
         "old_coordinate_orders": {key: value for key, value in zip(("t", "x", "y"), orders)}}
        for name, chart, parameter, orders in entries
    ],
    "old_base_fibre_multiplicities": [int(value) for value in base_orders],
    "cartan_check": "The actual edge-chart orders satisfy E7_Cartan*m=(1,0,0,0,0,0,0), identifying the affine attachment without a Kodaira-label inference.",
    "compiler_instruction": "Use these actual valuations, together with marked smooth-point frames, to derive anti-nef complete ideals and finite quotient blocks for a chosen vertical correction.",
    "boundary": "This is a valuation atlas, not yet the all-edge transition module or a completed vertical-condition matrix.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92ACTUALE7ATLAS|t_orders=2,2,4,3,1,2,3|"
    "status=PASS_EXACT_H92_E7_VALUATION_ATLAS",
    flush=True,
)
