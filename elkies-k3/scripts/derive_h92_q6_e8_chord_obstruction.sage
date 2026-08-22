#!/usr/bin/env sage -python
"""Certify the E8 transverse obstruction for the first H92 q=6 chord.

At the II* place put ``u=1/t``, ``X=u^4*x`` and ``Y=u^6*y``.  If ``P`` is
the marked section, then its integral coordinates ``X_P,Y_P`` are units at
the E8 singular point and the chord has the local form

    q = u^-2 * (Y-Y_P)/(X-X_P).

It is natural, but wrong, to try to correct this by only two base Laurent
terms ``a/u^2+b/u``.  The terminal node at infinity in the explicit E8
resolution has exceptional valuation

    ord_E(u,X,Y) = (5,8,12).

After the unique choices of ``a`` and ``b`` that cancel the constant and
base-linear terms, the coefficient of ``X/u^2`` remains nonzero and has
order ``8-2*5=-2`` on ``E``.  Thus base-only compensation is not even
regular on that resolved chart.  This is a local obstruction only: it does
not construct the required non-Cartier E8 module or rule out the neighbour.
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
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e8-chord-obstruction.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(SECTION) == SECTION_SHA256
section = json.loads(SECTION.read_text())
anchor = SourceFileLoader("h92_e8_chord_obstruction_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)
assert B2

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, section["x_entrance_base"]["numerator_coefficients"])) / u_field(
    polynomial(u_ring, section["x_entrance_base"]["denominator_coefficients"])
)
y_p = u_field(polynomial(u_ring, section["y_entrance_base"]["numerator_coefficients"])) / u_field(
    polynomial(u_ring, section["y_entrance_base"]["denominator_coefficients"])
)

# The marked section meets the smooth identity component of the II* fibre,
# rather than the singular point which is resolved by the E8 chart tree.
x_p_integral = u**4 * x_p
y_p_integral = u**6 * y_p
assert x_p_integral.valuation() == y_p_integral.valuation() == 0
x0 = QQ(x_p_integral(0))
y0 = QQ(y_p_integral(0))
assert x0 and y0 and y0**2 == x0**3

# Restrict the unit factor of q to the singular branch X=Y=0.  Its first two
# coefficients are the only possible base-only Laurent corrections.
series_ring = PowerSeriesRing(QQ, "u", default_prec=3)
u_series = series_ring.gen()
y_series = series_ring(y_p_integral.numerator()(u_series)) / series_ring(
    y_p_integral.denominator()(u_series)
)
x_series = series_ring(x_p_integral.numerator()(u_series)) / series_ring(
    x_p_integral.denominator()(u_series)
)
ratio_series = y_series / x_series
a = QQ(ratio_series[0])
b = QQ(ratio_series[1])
assert a == y0 / x0

# In G=(Y-Y_P)/(X-X_P), the X derivative at the singular point is y0/x0^2.
# It is nonzero, so its pullback is the first transverse surviving term.
x_coefficient = y0 / x0**2
y_coefficient = -1 / x0
assert x_coefficient and y_coefficient

# This is the final ordinary-node chart from derive_h92_q6_e8_resolution.sage:
# (u,X,Y)=(x^3*v^2,x^5*v^3,x^7*v^4*w).  In any ordinary blow-up chart at
# its origin, the new exceptional parameter divides x, v, and w once.  The
# exceptional valuation is therefore the sum of the displayed exponents.
pre_node_exponents = {
    "u": (3, 2, 0),
    "X": (5, 3, 0),
    "Y": (7, 4, 1),
}
terminal_infinity_valuation = tuple(
    sum(pre_node_exponents[name]) for name in ("u", "X", "Y")
)
assert terminal_infinity_valuation == (5, 8, 12)
u_order, x_order, y_order = terminal_infinity_valuation
assert x_order - 2 * u_order == -2
assert y_order - 2 * u_order == 2

payload = {
    "schema": "elkies-k3.h92-q6-e8-chord-obstruction.v1",
    "status": "PASS_BASE_ONLY_COMPENSATION_OBSTRUCTED",
    "inputs": {
        "h92_source": {
            "path": str(H92.relative_to(ROOT)),
            "sha256": digest(H92),
        },
        "marked_section": {
            "path": str(SECTION.relative_to(ROOT)),
            "sha256": SECTION_SHA256,
        },
        "resolution_charts": "elkies-k3/scripts/derive_h92_q6_e8_resolution.sage",
    },
    "integral_coordinates": "u=1/t, X=u^4*x, Y=u^6*y",
    "chord": "q=u^-2*(Y-Y_P)/(X-X_P)",
    "base_laurent_correction": {
        "a_over_u2": str(a),
        "b_over_u": str(b),
        "uniqueness": "a and b are forced by the X=Y=0 u-series through order one",
    },
    "terminal_infinity_exceptional_valuation": {
        "u_X_Y": list(terminal_infinity_valuation),
        "X_over_u2": x_order - 2 * u_order,
        "Y_over_u2": y_order - 2 * u_order,
        "X_coefficient": str(x_coefficient),
    },
    "conclusion": (
        "q-a/u^2-b/u has a nonzero X/u^2 term of exceptional order -2; "
        "therefore no correction by base Laurent terms alone is regular on this E8 chart."
    ),
    "not_claimed": [
        "the complete E8 non-Cartier module",
        "the full vertical-condition matrix",
        "a q6 neighbour pencil or child equation",
    ],
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "H92Q6E8CHORD|terminal_infinity_valuation=(5,8,12)|"
    "X_over_u2=-2|Y_over_u2=2|"
    "status=PASS_BASE_ONLY_COMPENSATION_OBSTRUCTED",
    flush=True,
)
