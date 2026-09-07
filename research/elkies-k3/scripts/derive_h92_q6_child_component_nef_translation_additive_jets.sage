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

from sage.all import EllipticCurve, PolynomialRing, QQ, matrix


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


def translated_germ_data(p0, ring, factor, base_point, coefficient_a, coefficient_b, orders):
    """Certify the local singular-germ action of ``tau_-P0``.

    This is deliberately one stage before a resolved-chart pullback.  At the
    additive cusp, ``x(P0)`` is a unit, so the group-law denominator
    ``X-x(P0)`` is a unit in the local surface germ.  The returned tangent map
    identifies the direction that a subsequent actual blow-up-chart compiler
    must transport; it does not claim to have already performed that lift.
    """

    u_ring = PolynomialRing(QQ, "u")
    u = u_ring.gen()
    u_field = u_ring.fraction_field()
    px, py = p0.xy()

    def localize(value):
        return u_field(u_ring(value.numerator()(base_point + u))) / u_field(
            u_ring(value.denominator()(base_point + u))
        )

    x0_u, y0_u = localize(px), localize(py)
    assert x0_u.valuation() == y0_u.valuation() == 0
    assert x0_u.numerator()(0) and x0_u.denominator()(0)
    assert y0_u.numerator()(0) and y0_u.denominator()(0)
    a_u = u_ring(coefficient_a(base_point + u))
    b_u = u_ring(coefficient_b(base_point + u))
    assert a_u.valuation() == orders[0] and b_u.valuation() == orders[1]
    assert y0_u**2 == x0_u**3 + u_field(a_u)*x0_u + u_field(b_u)

    x0, y0 = QQ(x0_u(0)), QQ(y0_u(0))
    assert x0 and y0 and y0**2 == x0**3
    lambda0 = -y0/x0
    # The image of the singular cusp under the rational group-law formula.
    # The Weierstrass equation of P0 forces both coordinates to vanish to
    # order at least three, so there is no linear u-term in the tangent map.
    cusp_x_image = y0_u**2/x0_u**2-x0_u
    cusp_y_image = (y0_u/x0_u)*cusp_x_image
    assert cusp_x_image.valuation() >= min(orders)
    assert cusp_y_image.valuation() >= min(orders)
    x_linear_y = 2*y0/x0**2
    assert lambda0**2 == x0
    # For lambda=(Y+y0)/(X-x0), x'=lambda^2-X-x0 and
    # y'=lambda*(X-x')-Y, these are the exact tangent coefficients at the
    # cusp. The displayed matrix has determinant one.
    tangent = matrix(QQ, [
        [1, 0, 0],
        [0, 1, x_linear_y],
        [0, 0, 1],
    ])
    assert tangent.det() == 1
    return {
        "formula": {
            "lambda": "(Y+y_P0(u))/(X-x_P0(u))",
            "X_image": "lambda^2-X-x_P0(u)",
            "Y_image": "lambda*(X-X_image)-Y",
        },
        "cusp_denominator": "X-x_P0(u)",
        "denominator_unit_at_cusp": True,
        "cusp_image_orders": {
            "X": int(cusp_x_image.valuation()),
            "Y": int(cusp_y_image.valuation()),
        },
        "tangent_pullback_mod_maximal_ideal_squared": {
            "u": "u",
            "X": "X+({})*Y".format(x_linear_y),
            "Y": "Y",
            "matrix_in_(u,X,Y)": [[str(value) for value in row] for row in tangent.rows()],
            "determinant": int(1),
        },
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
    germ = translated_germ_data(
        p0, ring, factor, base_point, A, B,
        (4, 5) if fibre["kodaira"] == "II*" else (3, 4),
    )
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
        "singular_germ_translation": germ,
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
    "conclusion": (
        "The translation centre P0 is smooth at both additive Weierstrass fibres, "
        "and its group-law formula is regular at each additive cusp with the "
        "recorded invertible tangent action."
    ),
    "boundary": (
        "This is a singular-germ and tangent-direction prerequisite only. It does "
        "not provide the required lifts through the actual II*/IV* blow-up charts, "
        "translated quotient modules, an infinity module, a global pencil, branch "
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
