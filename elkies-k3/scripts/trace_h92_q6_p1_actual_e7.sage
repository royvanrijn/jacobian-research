#!/usr/bin/env sage -python
"""Trace +/-P1 through the *actual* H92 E7 blow-up tree.

The first q=6 pencil uses the marked divisor O+(-P1)-F_infinity.  Its
earlier E7 branch-module note used a formal normal form, which is not an H92
chart transport.  This script closes the first indispensable gap: it follows
the certified exact P1 rational functions through the actual first, second,
and second-U-nonzero blow-ups from
``derive_h92_q6_actual_e7_resolution_full.sage``.

Both signs reach the second-U nonzero ordinary node and, after its Z-chart
blow-up, meet the actual exceptional curve E7_5 at two distinct smooth,
non-nodal points.  Thus the minuscule branch datum is attached to an actual
exceptional curve, rather than inferred from the III* label or a normal form.

This is deliberately a geometric-location certificate.  Constructing the
complete local module/trivialization on all H92 charts remains a separate
condition-matrix step.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, PowerSeriesRing, QQ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json"
P1_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
RESOLUTION_SHA256 = "14378f4718d3fbe781d5b351ba4943a962a430d9827c0ce285cd1125a9e8c500"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def at_old_base_zero(value, source_ring, target_ring):
    """Expand a rational function in u at t=1/u=0 exactly."""
    t = target_ring.gen()
    numerator = source_ring(value.numerator())
    denominator = source_ring(value.denominator())
    reversed_numerator = sum(
        target_ring(numerator[index]) * t ** (numerator.degree() - index)
        for index in range(numerator.degree() + 1)
    )
    reversed_denominator = sum(
        target_ring(denominator[index]) * t ** (denominator.degree() - index)
        for index in range(denominator.degree() + 1)
    )
    return t ** (denominator.degree() - numerator.degree()) * (
        reversed_numerator / reversed_denominator
    )


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.p1 == P1:
    assert digest(args.p1) == P1_SHA256
if args.resolution == RESOLUTION:
    assert digest(args.resolution) == RESOLUTION_SHA256
p1 = json.loads(args.p1.read_text())
resolution = json.loads(args.resolution.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"

anchor = SourceFileLoader("h92_p1_actual_e7_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(
    u_ring, p1["x_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, p1["x_entrance_base"]["denominator_coefficients"]
))
y_p = u_field(polynomial(
    u_ring, p1["y_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, p1["y_entrance_base"]["denominator_coefficients"]
))
assert y_p**2 == x_p**3 + (A1/u**3 + A/u**4)*x_p + B1/u**5 + B/u**6 + B2/u**7

series_ring = PowerSeriesRing(QQ, "t", default_prec=12)
t = series_ring.gen()
x = at_old_base_zero(x_p, u_ring, series_ring)
y = at_old_base_zero(y_p, u_ring, series_ring)
assert x.valuation() == 2 and y.valuation() == 3

# After the first Z chart and its second U chart the actual coordinates are
# (Z2,U2,Y2)=(t^2/x, x/t, y/x).  P1 hits the nonzero second-U node.
u2 = x/t
z2 = t**2/x
y2 = y/x
node_z = -A1/B1
assert z2[0] == node_z
delta = z2-node_z
assert delta.valuation() == 1

# The node's Z-chart is the E7_2--E7_5 edge chart.  Its map is
# (delta,U2,Y2)=(Z,Z*U,Z*Y), so P1 lies at Z=0 with the following nonzero
# coordinates on E7_5.  Negation only changes the Y coordinate.
z5 = delta
u5 = u2/delta
y5 = y2/delta
assert u5.valuation() == y5.valuation() == 0
u5_0, y5_0 = QQ(u5[0]), QQ(y5[0])
assert u5_0 and y5_0

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
f0 = Y**2-U**3-(A1*Z**3+A*Z**4)*U-(B1*Z**5+B*Z**6+B2*Z**7)
f1_z = strict(ring, f0, (Z, Z*U, Z*Y), Z)
f2_u = strict(ring, f1_z, (U*Z, U, U*Y), U)
second_translated = ring(f2_u(Z+node_z, U, Y))
n2_z = strict(ring, second_translated, (Z, Z*U, Z*Y), Z)
assert str(n2_z) == resolution["edge_charts"]["E7_2--E7_5"]

for sign, point in (("P1", (u5_0, y5_0)), ("-P1", (u5_0, -y5_0))):
    point_u, point_y = point
    assert n2_z(0, point_u, point_y) == 0
    gradient = (n2_z.derivative(Z)(0, point_u, point_y),
                n2_z.derivative(U)(0, point_u, point_y),
                n2_z.derivative(Y)(0, point_u, point_y))
    assert any(gradient)
    # The edge node is (Z,U,Y)=(0,0,0).  Both displayed coordinates are
    # nonzero, so the point is away from it and meets E7_5 alone.
    assert point_u != 0 and point_y != 0

payload = {
    "schema": "elkies-k3.h92-q6-p1-actual-e7-trace.v1",
    "status": "PASS_EXACT_P1_ACTUAL_E7_TRACE",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_resolution": {"path": str(args.resolution.relative_to(ROOT)), "sha256": digest(args.resolution)},
    },
    "blowup_path": [
        "first_Z_origin", "second_U", "second_U_nonzero_node",
        "second_U_nonzero_node_Z_chart",
    ],
    "second_U_node_coordinate": str(node_z),
    "node_chart": {
        "name": "E7_2--E7_5",
        "equation": str(n2_z),
        "coordinates": "(Z,U,Y)=(Z2-(-A1/B1), U2/(Z2+A1/B1), Y2/(Z2+A1/B1))",
        "P1": {"Z": "0", "U": str(u5_0), "Y": str(y5_0)},
        "-P1": {"Z": "0", "U": str(u5_0), "Y": str(-y5_0)},
    },
    "resolved_incidence": {
        "P1": "meets E7_5 at a smooth point away from E7_2",
        "-P1": "meets E7_5 at the distinct smooth point with opposite Y",
        "component": "E7_5",
    },
    "boundary": "This transports the marked P1 branch to the actual H92 resolution. It does not yet derive the full E7 line-bundle transition module or stack its quotient block into the q=6 global matrix.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6P1ACTUALE7|component=E7_5|edge=E7_2--E7_5|"
    "signs=distinct_smooth|status=PASS_EXACT_P1_ACTUAL_E7_TRACE",
    flush=True,
)
