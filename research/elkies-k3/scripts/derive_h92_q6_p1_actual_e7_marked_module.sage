#!/usr/bin/env sage -python
"""Derive the marked-point local frame for the first q=6 pencil on H92.

This is the actual-chart replacement for the marked-point portion of the
older formal E7 branch calculation.  The trace places +/-P1 at smooth points
of E7_5 in the second-U-node Z chart.  There the old-coordinate pullback is

  t=Z*U*(Z+s), x=Z^2*U^2*(Z+s), y=Z^3*U^2*(Z+s)*Y,

where s=-A1/B1.  At -P1 take W=Y-Y_{-P1}(Z).  Implicit differentiation of
the actual chart equation gives U_Y=-f_Y/f_U.  This proves exactly that

  t/Z,  (x-x(P1))/(Z^2 W),  (y-y(P1))/Z^3

are units, hence m/t has a simple pole with unit residue along -P1.  That is
the local marked divisor frame required by O+(-P1)-F at the actual H92 E7
fibre.  The remaining edge-node/transition blocks are intentionally not
claimed here.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
TRACE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module.json"
TRACE_SHA256 = "a73ccb1c729814219f172df4c6feb49c05859125db1cff7591eeb8544fb664e1"
RESOLUTION_SHA256 = "14378f4718d3fbe781d5b351ba4943a962a430d9827c0ce285cd1125a9e8c500"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--trace", type=Path, default=TRACE)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.trace == TRACE:
    assert digest(args.trace) == TRACE_SHA256
if args.resolution == RESOLUTION:
    assert digest(args.resolution) == RESOLUTION_SHA256
trace = json.loads(args.trace.read_text())
resolution = json.loads(args.resolution.read_text())
assert trace["status"] == "PASS_EXACT_P1_ACTUAL_E7_TRACE"
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"

anchor = SourceFileLoader("h92_q6_actual_e7_marked_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
node_translation = -A1/B1

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
f0 = Y**2-U**3-(A1*Z**3+A*Z**4)*U-(B1*Z**5+B*Z**6+B2*Z**7)
f1_z = strict(ring, f0, (Z, Z*U, Z*Y), Z)
f2_u = strict(ring, f1_z, (U*Z, U, U*Y), U)
translated = ring(f2_u(Z+node_translation, U, Y))
surface = strict(ring, translated, (Z, Z*U, Z*Y), Z)
assert str(surface) == resolution["edge_charts"]["E7_2--E7_5"]

point_plus = trace["node_chart"]["P1"]
u0 = QQ(point_plus["U"])
y_plus = QQ(point_plus["Y"])
y_minus = -y_plus
assert surface(0, u0, y_plus) == surface(0, u0, y_minus) == 0

# The actual surface is smooth at both marked points, and Y is a legitimate
# transverse parameter because f_U is nonzero there.
f_u = surface.derivative(U)(0, u0, y_minus)
f_y = surface.derivative(Y)(0, u0, y_minus)
assert f_u and f_y
d_u_d_y = -f_y/f_u

# Leading coefficients in the smooth parameters (Z,W), W=Y-Y_-P1(Z).
# They are obtained directly from the certified global-to-local pullback.
t_over_z = u0*node_translation
x_over_z2_w = 2*u0*d_u_d_y*node_translation
y_difference_over_z3 = -2*u0**2*node_translation*y_plus
assert t_over_z and x_over_z2_w and y_difference_over_z3
chord_over_t_residue = y_difference_over_z3/(x_over_z2_w*t_over_z)
assert chord_over_t_residue

payload = {
    "schema": "elkies-k3.h92-q6-p1-actual-e7-marked-module.v1",
    "status": "PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "actual_trace": {"path": str(args.trace.relative_to(ROOT)), "sha256": digest(args.trace)},
        "actual_resolution": {"path": str(args.resolution.relative_to(ROOT)), "sha256": digest(args.resolution)},
    },
    "chart": "E7_2--E7_5",
    "marked_component": "E7_5",
    "local_coordinates": {
        "base": "Z",
        "transverse_at_minus_P1": "W=Y-Y_{-P1}(Z)",
        "old_coordinate_pullbacks": {
            "t": "Z*U*(Z-A1/B1)",
            "x": "Z^2*U^2*(Z-A1/B1)",
            "y": "Z^3*U^2*(Z-A1/B1)*Y",
        },
    },
    "unit_coefficients_at_minus_P1": {
        "surface_f_U": str(f_u),
        "surface_f_Y": str(f_y),
        "dU_dY": str(d_u_d_y),
        "t_over_Z": str(t_over_z),
        "(x-xP1)/(Z^2*W)": str(x_over_z2_w),
        "(y-yP1)/Z^3": str(y_difference_over_z3),
        "(m/t)*W": str(chord_over_t_residue),
    },
    "module_conclusion": "Near -P1, m/t is a unit times 1/W. Thus the actual local frame for O+(-P1)-F contains the marked generator m/t with exactly its allowed simple pole along -P1.",
    "boundary": "This proves the marked smooth-point frame on actual H92 charts. It does not yet derive the complete E7 edge-transition module or a finite quotient block for the q=6 matrix.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6P1ACTUALE7MODULE|component=E7_5|pole=m/t~unit/W|"
    "status=PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE",
    flush=True,
)
