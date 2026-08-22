#!/usr/bin/env sage -python
"""Record additive-fibre jets of the old-zero translation centre P0.

The component-nef chord is the pullback of the standard marked chord through
the group translation tau_-P0.  This script checks, in the exact minimal
Weierstrass coordinates, whether P0 specializes to the smooth locus at II*
and IV*.  It also records the translated marked point Q=P0+S.  These are the
inputs needed to pull resolved chart modules through that translation; no
chart pullback or global pencil is asserted here.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
ZERO = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json"
CHORD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-chord.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-translation-additive-jets.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational(field, ring, data, numerator, denominator):
    return field(polynomial(ring, data[numerator])) / field(polynomial(ring, data[denominator]))


def order_at(value, factor):
    return int(value.numerator().valuation(factor) - value.denominator().valuation(factor))


def affine_point_data(point, ring, factor):
    x_value, y_value = point.xy()
    return {
        "x_order": order_at(x_value, factor),
        "y_order": order_at(y_value, factor),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--zero", type=Path, default=ZERO)
parser.add_argument("--chord", type=Path, default=CHORD)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "zero", "chord", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
zero = json.loads(args.zero.read_text())
chord = json.loads(args.chord.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert chord["status"] == "PASS_EXACT_COMPONENT_NEF_OLD_ZERO_CHORD"

ring = PolynomialRing(QQ, "T")
field = ring.fraction_field()
A = field(polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]))
B = field(polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]))
curve = EllipticCurve(field, [0, 0, 0, A, B])
p0_data = chord["points_in_standard_weierstrass_group"]["transported_old_zero_P0"]
q_data = chord["points_in_standard_weierstrass_group"]["physical_second_horizontal_point_Q_equals_P0_plus_S"]
p0 = curve(
    rational(field, ring, p0_data, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
    rational(field, ring, p0_data, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
)
q = curve(
    rational(field, ring, q_data, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
    rational(field, ring, q_data, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
)

results = {}
for fibre in child["finite_fibres"]:
    if fibre["kodaira"] not in {"II*", "IV*"}:
        continue
    factor = ring(fibre["factor"])
    assert factor.degree() == 1
    base_point = -factor[0] / factor[1]
    p0_local = affine_point_data(p0, ring, factor)
    q_local = affine_point_data(q, ring, factor)
    assert p0_local["x_order"] >= 0 and p0_local["y_order"] >= 0
    px, py = p0.xy()
    x0, y0 = QQ(px(base_point)), QQ(py(base_point))
    a0, b0 = QQ(A(base_point)), QQ(B(base_point))
    assert y0**2 == x0**3 + a0*x0 + b0
    # The affine fibre is singular only at (0,0) for both additive short
    # models.  A nonzero gradient is the exact smooth-locus criterion.
    gradient = (-3*x0**2-a0, 2*y0)
    smooth = gradient != (0, 0)
    assert smooth
    results[fibre["kodaira"]] = {
        "base_factor": str(factor),
        "base_point": str(base_point),
        "P0": {
            **p0_local,
            "specialization": [str(x0), str(y0)],
            "weierstrass_gradient": [str(value) for value in gradient],
            "smooth_locus": True,
        },
        "Q": q_local,
        "translation": "tau_-P0 extends over the Neron smooth locus",
    }
assert set(results) == {"II*", "IV*"}

payload = {
    "schema": "elkies-k3.h92-q6-child-component-nef-translation-additive-jets.v1",
    "status": "PASS_EXACT_COMPONENT_NEF_TRANSLATION_ADDITIVE_JETS",
    "inputs": {
        "child": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "old_zero": {"path": str(args.zero.relative_to(ROOT)), "sha256": digest(args.zero)},
        "component_nef_chord": {"path": str(args.chord.relative_to(ROOT)), "sha256": digest(args.chord)},
    },
    "additive_fibres": results,
    "conclusion": "The translation centre P0 is smooth at both additive Weierstrass fibres.",
    "boundary": (
        "This does not provide resolved-chart pullbacks of tau_-P0, translated "
        "II*/IV* quotient modules, an infinity module, a global pencil, branch "
        "divisor, extension collision, or rank claim."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDCOMPNEFJETS|II*=P0_smooth|IV*=P0_smooth|"
    "status=PASS_EXACT_COMPONENT_NEF_TRANSLATION_ADDITIVE_JETS",
    flush=True,
)
