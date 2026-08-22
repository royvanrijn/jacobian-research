#!/usr/bin/env sage -python
"""Certify transport of the exact nef finite q8 module to the component-nef divisor.

The actual component-nef fibre has exceptional cycles equal to the negatives
of the previously derived abstract-nef valuation cycles.  Hence, relative to
the horizontal divisor, its local line bundle is O(-Z), where Z is represented
by the complete ideal (u^2,X,Y) at both II* and IV*.

The physical horizontal marking differs from the displayed Weierstrass marking
by translation by the transported old zero P0.  P0 specializes to the smooth
locus of both additive Weierstrass fibres, hence to the identity component.
Translation by P0 therefore preserves every component and the complete ideal
O(-Z).  Pulling back the q_regular finite module gives the same polynomial
module on the component-nef chord.

This uses the standard Neron-model fact that translation by an identity-
component section preserves the component-labelled vertical line bundle.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
ZERO = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json"
DEFAULT_COMPONENT = LOCAL / "q8-target-component-nef.json"
DEFAULT_NEF = LOCAL / "q8-target-nef.json"
DEFAULT_FINITE = LOCAL / "q8-qregular-finite-component-nef.json"
DEFAULT_OUTPUT = LOCAL / "q8-component-nef-finite-module-certificate.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational(field, ring, data, numerator, denominator):
    return field(polynomial(ring, data[numerator])) / field(polynomial(ring, data[denominator]))


def order_at(value, factor):
    return int(value.numerator().valuation(factor) - value.denominator().valuation(factor))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--component-target", type=Path, default=DEFAULT_COMPONENT)
parser.add_argument("--nef-target", type=Path, default=DEFAULT_NEF)
parser.add_argument("--finite", type=Path, default=DEFAULT_FINITE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("component_target", "nef_target", "finite", "output"):
    setattr(args, name, getattr(args, name).resolve())

component = json.loads(args.component_target.read_text())
nef = json.loads(args.nef_target.read_text())
finite = json.loads(args.finite.read_text())
child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())

assert component["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert nef["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert component["normalization"]["representative"] == "component-nef"
assert nef["normalization"]["representative"] == "nef"
assert finite["status"] == "PASS_EXACT_Q_REGULAR_FINITE_MODULE"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"

component_selected = component["selected_q8"]
nef_selected = nef["selected_q8"]
for root_name in ("E6", "E8"):
    c_cycle = list(map(int, component_selected[root_name]["vertical_cycle"]))
    n_cycle = list(map(int, nef_selected[root_name]["vertical_cycle"]))
    c_degrees = list(map(int, component_selected[root_name]["component_degrees"]))
    n_degrees = list(map(int, nef_selected[root_name]["component_degrees"]))
    assert c_cycle == [-value for value in n_cycle]
    assert c_degrees == [-value for value in n_degrees]

assert finite["module"]["basis"] == [
    ["1", "lift(R/Nx)"],
    ["0", "f_II*^2*f_IV*^2"],
]
assert finite["module"]["smith_degrees"] == [0, 4]
assert finite["module"]["finite_codimension"] == 4

ring = PolynomialRing(QQ, "T")
field = ring.fraction_field()
A = field(polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]))
B = field(polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]))
curve = EllipticCurve(field, [0, 0, 0, A, B])
p0_data = zero["section"]
p0 = curve(
    rational(field, ring, p0_data,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(field, ring, p0_data,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)
px, py = p0.xy()

additive = {}
for fibre in child["finite_fibres"]:
    if fibre["kodaira"] not in ("II*", "IV*"):
        continue
    factor = ring(fibre["factor"])
    assert factor.degree() == 1
    point = -factor[0] / factor[1]
    x_order = order_at(px, factor)
    y_order = order_at(py, factor)
    assert x_order >= 0 and y_order >= 0
    x0, y0 = QQ(px(point)), QQ(py(point))
    a0, b0 = QQ(A(point)), QQ(B(point))
    assert y0**2 == x0**3 + a0*x0 + b0
    gradient = (-3*x0**2-a0, 2*y0)
    assert gradient != (0, 0)
    additive[fibre["kodaira"]] = {
        "factor": str(factor),
        "P0_orders": [x_order, y_order],
        "P0_specialization": [str(x0), str(y0)],
        "gradient": [str(value) for value in gradient],
        "identity_component": True,
    }
assert set(additive) == {"II*", "IV*"}

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-component-nef-finite-module.v1",
    "status": "PASS_EXACT_COMPONENT_NEF_FINITE_MODULE_TRANSPORT",
    "inputs": {
        "component_target": digest(args.component_target),
        "nef_target": digest(args.nef_target),
        "q_regular_finite_module": digest(args.finite),
        "child": digest(CHILD),
        "old_zero": digest(ZERO),
    },
    "sign_duality": {
        "E6": "component-nef cycle = -(abstract-nef valuation cycle)",
        "E8": "component-nef cycle = -(abstract-nef valuation cycle)",
        "line_bundle": "O(component-nef vertical correction)=O(-Z)",
        "complete_ideal_at_II_star": "(u^2,X,Y)",
        "complete_ideal_at_IV_star": "(u^2,X,Y)",
    },
    "translation": {
        "centre": "transported old zero P0",
        "additive_fibres": additive,
        "component_group_action": "identity",
        "conclusion": (
            "tau_-P0 preserves the component-labelled complete ideals and "
            "transports the q_regular finite module unchanged"
        ),
    },
    "module": finite["module"],
    "frame_change": finite["frame_change"],
    "additive_CRT": finite["additive_CRT"],
    "boundary": (
        "This certifies the finite additive module for the component-nef "
        "degree-two divisor.  It does not certify the infinity lattice, a "
        "two-dimensional global intersection, a q8 equation, or rank."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q8COMPNEFFINITE|sign_dual_nef=1|P0_identity_II=1|P0_identity_IV=1|"
    "smith=0,4|codim=4|status=PASS_EXACT_COMPONENT_NEF_FINITE_MODULE_TRANSPORT",
    flush=True,
)
