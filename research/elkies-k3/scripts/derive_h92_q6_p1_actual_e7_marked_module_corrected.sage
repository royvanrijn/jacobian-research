#!/usr/bin/env sage -python
"""Derive the corrected marked q=6 frame on the actual E7_2--E7_5 chart.

For the raw chord ``m=(y-y(P1))/(x-x(P1))``, the actual pullback has
``ord_Z(m/t)=-1`` along the generic E7_5 component.  The q=6 E7 line-bundle
trivialization therefore uses

    n = Z*m/t,

not ``m/t``.  At the marked point ``-P1``, with
``W=Y-Y_{-P1}(Z)``, the exact local calculation gives ``n=unit/W``.  Since
``t/Z`` is a unit, the raw generic frame ``<1,m>`` is still a valid local
frame after this corrected trivialization.
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
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
TRACE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


def invert_base(rational_u):
    """Rewrite a QQ(u) function as a QQ(t) function, t=1/u."""
    numerator = rational_u.numerator()
    denominator = rational_u.denominator()
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    t_field = t_ring.fraction_field()
    return t_field(
        t**(denominator.degree()-numerator.degree())
        * t_ring(list(reversed(numerator.list())))
        / t_ring(list(reversed(denominator.list())))
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--trace", type=Path, default=TRACE)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
trace = json.loads(args.trace.read_text())
resolution = json.loads(args.resolution.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert trace["status"] == "PASS_EXACT_P1_ACTUAL_E7_TRACE"
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"

anchor = SourceFileLoader("h92_q6_corrected_marked_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
second = -A1/B1

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p_u = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p_u /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p_u = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p_u /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p_u)
y_p_t = invert_base(y_p_u)
t_ring = x_p_t.parent()
t = t_ring.gen()
x_series = x_p_t/t**2
y_series = y_p_t/t**3
c2 = QQ(x_series(0))
c3 = QQ(x_series.derivative()(0))
d3 = QQ(y_series(0))
assert c2 == -B1/A1 and c2*second == 1

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
f0 = Y**2-U**3-(A1*Z**3+A*Z**4)*U-(B1*Z**5+B*Z**6+B2*Z**7)
f1_z = strict(ring, f0, (Z, Z*U, Z*Y), Z)
f2_u = strict(ring, f1_z, (U*Z, U, U*Y), U)
translated = ring(f2_u(Z+second, U, Y))
surface = strict(ring, translated, (Z, Z*U, Z*Y), Z)
assert str(surface) == resolution["edge_charts"]["E7_2--E7_5"]

point_plus = trace["node_chart"]["P1"]
u0 = QQ(point_plus["U"])
y_plus = QQ(point_plus["Y"])
y_minus = -y_plus
assert surface(0, u0, y_minus) == 0
f_u = surface.derivative(U)(0, u0, y_minus)
f_y = surface.derivative(Y)(0, u0, y_minus)
assert f_u and f_y
d_u_d_y = -f_y/f_u

# t/Z is a unit.  Only the first two x(P1) Laurent coefficients can affect
# the Z^3 term; higher terms contain Z^4.  The leading x difference vanishes
# at -P1 and has a nonzero W derivative, while the leading y difference is a
# unit there.
t_over_z = u0*second
x_z3_coefficient = ring(
    -c2*U**2*second - c3*U**3*second**3
)
assert x_z3_coefficient.subs({U: u0}) == 0
x_over_z3_w = x_z3_coefficient.derivative(U).subs({U: u0})*d_u_d_y
y_z3_coefficient = ring(U**2*second*Y-d3*U**3*second**3)
y_over_z3 = y_z3_coefficient.subs({U: u0, Y: y_minus})
assert t_over_z and x_over_z3_w and y_over_z3
corrected_residue = y_over_z3/(x_over_z3_w*t_over_z)
assert corrected_residue

payload = {
    "schema": "elkies-k3.h92-q6-p1-actual-e7-marked-module-corrected.v1",
    "status": "PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE_CORRECTED",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
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
    "corrected_frame": {
        "raw_chord": "m=(y-y(P1))/(x-x(P1))",
        "generator": "Z*m/t",
        "generic_E7_5_order_of_m_over_t": -1,
        "unit_relation": "m=(t/Z)*(Z*m/t), with t/Z a unit",
    },
    "unit_coefficients_at_minus_P1": {
        "surface_f_U": str(f_u),
        "surface_f_Y": str(f_y),
        "dU_dY": str(d_u_d_y),
        "t_over_Z": str(t_over_z),
        "(x-xP1)/(Z^3*W)": str(x_over_z3_w),
        "(y-yP1)/Z^3": str(y_over_z3),
        "(Z*m/t)*W": str(corrected_residue),
    },
    "module_conclusion": (
        "Near -P1, Z*m/t is a unit times 1/W. Since t/Z is a unit, the raw "
        "chord m is a unit times the corrected generator; hence <1,m> is the "
        "valid local q6 frame on this marked chart."
    ),
    "boundary": (
        "This corrects the marked E7_5 trivialization only. The all-edge E7 "
        "transition module and any q8 ninth-power quotient remain to be derived."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6P1ACTUALE7MODULECORRECTED|component=E7_5|generator=Z*m/t~unit/W|"
    "status=PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE_CORRECTED",
    flush=True,
)
