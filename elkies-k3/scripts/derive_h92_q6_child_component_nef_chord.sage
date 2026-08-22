#!/usr/bin/env sage -python
"""Transport the marked q8 chord to the component-nef old-zero divisor.

The displayed short Weierstrass child uses its infinity section as origin,
whereas the physical component-nef divisor is marked relative to the finite
transported old zero ``P0``.  If ``S`` is the existing standard-coordinate
marked point, its old-zero divisor is ``P0+Q`` with ``Q=P0+S``.  Translation
by ``-P0`` sends that divisor to ``O+S``.  Pulling the usual chord back gives
the exact generic degree-two function for the physical class:

  lambda=(y+y(P0))/(x-x(P0)),
  x'=lambda^2-x-x(P0),  y'=lambda*(x-x')-y,
  m_old=(y'+y(S))/(x'-x(S)).

This is a generic-fibre marking only.  It deliberately does not infer the
required resolved II*/IV* chart trivializations after this translation.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
ZERO = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-physical-root-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-chord.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational(field, ring, data, numerator, denominator):
    return field(polynomial(ring, data[numerator])) / field(polynomial(ring, data[denominator]))


def point_data(point, ring):
    x_value, y_value = point.xy()
    return {
        "x_numerator_coefficients_low_to_high": [str(value) for value in ring(x_value.numerator()).list()],
        "x_denominator_coefficients_low_to_high": [str(value) for value in ring(x_value.denominator()).list()],
        "y_numerator_coefficients_low_to_high": [str(value) for value in ring(y_value.numerator()).list()],
        "y_denominator_coefficients_low_to_high": [str(value) for value in ring(y_value.denominator()).list()],
        "coordinate_degrees": {
            "x": [int(ring(x_value.numerator()).degree()), int(ring(x_value.denominator()).degree())],
            "y": [int(ring(y_value.numerator()).degree()), int(ring(y_value.denominator()).degree())],
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--zero", type=Path, default=ZERO)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "zero", "marking", "target", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
zero = json.loads(args.zero.read_text())
marking = json.loads(args.marking.read_text())
target = json.loads(args.target.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

ring = PolynomialRing(QQ, "T")
field = ring.fraction_field()
A = field(polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]))
B = field(polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]))
curve = EllipticCurve(field, [0, 0, 0, A, B])
p0_data = zero["section"]
p0 = curve(
    rational(field, ring, p0_data, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
    rational(field, ring, p0_data, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
)
s_data = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
s = curve(
    rational(field, ring, s_data, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
    rational(field, ring, s_data, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
)
assert not p0.is_zero() and not s.is_zero()
q = p0 + s
assert q - p0 == s and q != p0

# On the generic Weierstrass curve, these formulas are the group-law map
# tau_{-P0}; the formula for m_old is m_S after that translation.  Keep the
# variables formal to avoid conflating this generic function with its later
# resolved chart pullbacks.
x0, y0 = p0.xy()
sx, sy = s.xy()
formula = {
    "translation": "tau_-P0",
    "lambda": "(y+y0)/(x-x0)",
    "x_prime": "lambda^2-x-x0",
    "y_prime": "lambda*(x-x_prime)-y",
    "m_component_nef": "(y_prime+yS)/(x_prime-xS)",
}

payload = {
    "schema": "elkies-k3.h92-q6-child-component-nef-chord.v1",
    "status": "PASS_EXACT_COMPONENT_NEF_OLD_ZERO_CHORD",
    "inputs": {
        "child": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "old_zero": {"path": str(args.zero.relative_to(ROOT)), "sha256": digest(args.zero)},
        "standard_chord_marking": {"path": str(args.marking.relative_to(ROOT)), "sha256": digest(args.marking)},
        "component_nef_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
    },
    "points_in_standard_weierstrass_group": {
        "transported_old_zero_P0": point_data(p0, ring),
        "translated_marked_point_S": point_data(s, ring),
        "physical_second_horizontal_point_Q_equals_P0_plus_S": point_data(q, ring),
        "checks": ["Q-P0=S", "P0 is nonzero", "S is nonzero"],
    },
    "generic_divisor": {
        "physical_old_zero_marking": "P0+Q",
        "translation_image": "O_standard+S",
        "translation": "tau_-P0",
        "generic_rr_basis": ["1", "m_component_nef"],
        "justification": "tau_-P0 pulls L(O_standard+S) back isomorphically to L(P0+Q)",
    },
    "generic_chord_formula": formula,
    "boundary": (
        "This identifies the exact generic chord for the component-nef lattice "
        "pencil. It does not transport the II*/IV* resolved chart modules through "
        "tau_-P0, impose the infinity module, construct a global pencil or branch "
        "divisor, or supply an extension collision or rank result."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDCOMPNEFCHORD|translation=tau_-P0|generic_basis=2|"
    "status=PASS_EXACT_COMPONENT_NEF_OLD_ZERO_CHORD",
    flush=True,
)
