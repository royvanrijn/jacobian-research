#!/usr/bin/env sage -python
"""Certify the q8 E7_3--E7_6 Artinian-corner good-reduction obstruction.

The expensive local quotient R/(t^17) is not replaced here.  Instead the
Artinian quotient (surface,Z^34,U^34) is a *one-way* target: it contains
(t^17) and is supported at the actual chart origin.  A full-rank ambient map
into this quotient therefore excludes any true local solution in that fixed
ambient after good reduction.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
CLEARINGS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-node-principal-clearings.json"
CORNER = ROOT / "artifacts/local/elkies-k3-h92-q8-e7-3-6-corner-obstruction-mod-43.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-3-6-corner-obstruction-good-reduction.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def common_monomial_exponents(value):
    terms = list(value.dict())
    assert terms
    return tuple(min(exponent[index] for exponent in terms) for index in range(3))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--clearings", type=Path, default=CLEARINGS)
parser.add_argument("--corner", type=Path, default=CORNER)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
args.pullbacks = args.pullbacks.resolve()
args.clearings = args.clearings.resolve()
args.corner = args.corner.resolve()
args.output = args.output.resolve()

pullbacks = json.loads(args.pullbacks.read_text())
clearings = json.loads(args.clearings.read_text())
corner = json.loads(args.corner.read_text())
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert clearings["status"] == "PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS"
assert corner["status"] == "EXPERIMENTAL_MODULAR_Q8_E7_NODE_FINITE_CORNER_OBSTRUCTION"
assert corner["local_ring"]["chart"] == "E7_3--E7_6"
assert corner["prime"] == 43
assert corner["inputs"]["node_clearings"]["sha256"] == digest(args.clearings)

T = int(clearings["common_parameters"]["T"])
assert T == 17
chart = next(item for item in pullbacks["charts"] if item["name"] == "E7_3--E7_6")
ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
t_value = ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Z, "U": U, "Y": Y}))
t_exponents = common_monomial_exponents(t_value)
assert t_exponents == (2, 2, 0)
t_monomial = ring.monomial(*t_exponents)
t_unit = t_value // t_monomial
assert t_value == t_monomial*t_unit and t_unit(0, 0, 0) != 0

z_power, u_power = T*t_exponents[0], T*t_exponents[1]
assert (z_power, u_power) == (34, 34)
assert corner["local_ring"]["principal_ideal"] == "(surface,Z^34,U^34)"
# At Z=U=0 the actual surface is Y^2.  Hence in the displayed corner,
# Z and U are nilpotent and so is Y; the quotient is supported only at the
# chart origin and all local units map to units.
assert surface(0, 0, Y) == Y**2
# t^T=Z^34 U^34 times a chart unit, so (t^T) is contained in the corner.
t_power_quotient, t_power_remainder = (t_value**T).quo_rem(Z**z_power*U**u_power)
assert not t_power_remainder and t_power_quotient(0, 0, 0) != 0

image = corner["finite_ambient_image"]
ambient_dimension = int(image["ambient_dimension"])
assert ambient_dimension == 54
assert int(image["rank"]) == ambient_dimension
assert int(image["kernel_dimension"]) == 0
for residue in corner["good_reduction"]["common_clearing_unit_residues"].values():
    assert int(residue) % int(corner["prime"]) != 0

payload = {
    "schema": "elkies-k3.h92-q8-e7-3-6-corner-obstruction-good-reduction.v1",
    "status": "PASS_EXACT_Q8_E7_3_6_CORNER_OBSTRUCTION_INJECTIVITY",
    "inputs": {
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "node_clearings": {"path": str(args.clearings.relative_to(ROOT)), "sha256": digest(args.clearings)},
        "corner_image": {"path": str(args.corner.relative_to(ROOT)), "sha256": digest(args.corner)},
    },
    "actual_local_geometry": {
        "chart": "E7_3--E7_6",
        "t": str(t_value),
        "t_monomial": "Z^2*U^2",
        "t_unit_at_origin_nonzero": True,
        "T": T,
        "corner_ideal": "(surface,Z^34,U^34)",
        "containment": "(t^17) is contained in (Z^34,U^34) after localizing at the chart origin",
        "support": "surface(0,0,Y)=Y^2, so the corner quotient is Artinian and supported at (Z,U,Y)=0",
    },
    "finite_ambient_image": {
        "prime": int(corner["prime"]),
        "ambient_dimension": ambient_dimension,
        "rows": int(image["rows"]),
        "rank": int(image["rank"]),
        "kernel_dimension": int(image["kernel_dimension"]),
    },
    "good_reduction_argument": (
        "A true characteristic-zero local relation has common-cleared "
        "numerator in (t^17), hence maps to zero in the Artinian corner. "
        "After primitive normalization it reduces nontrivially modulo 43; "
        "the displayed common-clearing factors remain units there. Full "
        "column rank of the corner image therefore excludes a nonzero "
        "characteristic-zero true local relation in this fixed 54-column ambient."
    ),
    "boundary": (
        "This proves a bounded-ambient one-way local obstruction only. It "
        "does not identify R/(t^17) with the corner quotient, give a "
        "characteristic-zero coordinate matrix, impose overlaps, or construct "
        "a q8 pencil or child model."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E736CORNERGOODREDUCTION|prime=43|ambient={}|rank={}|kernel={}|"
    "status=PASS_EXACT_Q8_E7_3_6_CORNER_OBSTRUCTION_INJECTIVITY".format(
        ambient_dimension, image["rank"], image["kernel_dimension"]
    ),
    flush=True,
)
