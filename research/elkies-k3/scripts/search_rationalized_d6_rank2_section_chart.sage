#!/usr/bin/env sage-python
"""Replay a rational D6 marked-section chart and its exact two-section gate.

status: ACTIVE_SEARCH
claim: the Kimura D6 chart is rational over QQ, one polynomial section has a
       two-parameter chart, and that chart contains no nontrivial section pair
inputs: --height-bound (default 30)
outputs: elkies-k3-rationalized-d6-section-chart-search-v1.json

The exact obstruction concerns only pairs obtained from this particular
polynomial marked-section chart.  It is not a nonexistence theorem for a
second section in a larger D6 chart or for a rank-sum-four family.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import gcd
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-rationalized-d6-section-chart-search-v1.json"
)


R = PolynomialRing(QQ, names=("z", "q", "u"))
z, q, u = R.gens()
p = q + 2
y3 = -(64 * z**4 + q**4) / (64 * z**3 * q)
a = (-64 * z**4 + 3 * q**4 + 24 * q**3) / (32 * z**2 * q)
b = -(-64 * z**4 * q + 3 * q**5 + 256 * z**4 + 12 * q**4) / (256 * z**4)
g = y3**2

A = u**2 * (-u**2 + a * u - 3)
B = u**3 * (g * u**3 + b * u**2 + a * u - 2)
x_section = z**2 + p * u
y_section = (
    z**3
    + 3 * z * p * u / 2
    + 3 * (p**2 - 4) * u**2 / (8 * z)
    + y3 * u**3
)
if y_section**2 != x_section**3 + A * x_section + B:
    raise ArithmeticError("rationalized D6 section identity failed")

# At u=0, the reduced discriminant has order exactly two beyond the common
# u^6 factor, hence the fibre is generically I_2^*.  The remaining factor has
# degree four.  At infinity the leading discriminant is -4+27*g^2, so it is
# generically smooth; Kimura's sqrt(3) occurs only when this leading term is
# normalized to vanish.
F = -u**2 + a * u - 3
G = g * u**3 + b * u**2 + a * u - 2
reduced = 4 * F**3 + 27 * G**2
if (
    reduced(u=0) != 0
    or reduced.derivative(u)(u=0) != 0
    or reduced.derivative(u, 2)(u=0) == 0
):
    raise ArithmeticError("generic I2* valuation at u=0 changed")
residual = reduced / u**2
if residual.numerator().degree(u) - residual.denominator().degree(u) != 4:
    raise ArithmeticError("D6 residual discriminant degree changed")

# Put h=q/z.  The leading y coefficient is -(h^4+64)/(64*h).  Two sections
# in this chart can lie on the same D6 equation only if their leading
# coefficients agree up to sign.  Apart from the duplicate j=h and the
# sign duplicate j=-h, the two necessary correspondences are
#
#   h*j*(h^2+h*j+j^2)=64,
#   h*j*(h^2-h*j+j^2)=-64.
H_RING = PolynomialRing(QQ, names=("h", "j"))
h, j = H_RING.gens()
same_numerator = (h**4 + 64) * j - (j**4 + 64) * h
opposite_numerator = (h**4 + 64) * j + (j**4 + 64) * h
same_expected = -(j - h) * (h**3 * j + h**2 * j**2 + h * j**3 - 64)
opposite_expected = (h + j) * (h**3 * j - h**2 * j**2 + h * j**3 + 64)
if same_numerator != same_expected or opposite_numerator != opposite_expected:
    raise ArithmeticError("D6 leading-coefficient correspondence changed")

# Both nontrivial correspondence curves are birational to
#
#     E: Y^2=X^3+X^2+X.
#
# For the same-sign curve put k=j/h, x=h*j and Y=8*k/x, X=k.  For the
# opposite-sign curve use the same k,x,Y and X=-k.  E(Q) has rank zero and
# torsion Z/2={O,(0,0)}.  The affine torsion point forces k=0, while O lies at
# the omitted boundary, so neither curve has a nondegenerate rational point.
correspondence_jacobian = EllipticCurve(QQ, [0, 1, 0, 1, 0])
correspondence_rank = correspondence_jacobian.rank(proof=True)
correspondence_torsion = correspondence_jacobian.torsion_points()
if correspondence_rank != 0:
    raise ArithmeticError("D6 correspondence Jacobian rank changed")
if len(correspondence_torsion) != 2:
    raise ArithmeticError("D6 correspondence Jacobian torsion changed")
affine_torsion = [point for point in correspondence_torsion if point != correspondence_jacobian(0)]
if len(affine_torsion) != 1 or affine_torsion[0][0] != 0 or affine_torsion[0][1] != 0:
    raise ArithmeticError("D6 correspondence torsion support changed")


def rationals_of_height(bound):
    values = set()
    for denominator in range(1, bound + 1):
        for numerator in range(-bound, bound + 1):
            if numerator and gcd(abs(numerator), denominator) == 1:
                values.add(Fraction(numerator, denominator))
    return sorted(values)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--height-bound", type=int, default=30)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.height_bound < 1:
    raise SystemExit("--height-bound must be positive")

values = rationals_of_height(arguments.height_bound)
same_hits = []
opposite_hits = []
for h_value in values:
    for j_value in values:
        if j_value != h_value and h_value * j_value * (
            h_value**2 + h_value * j_value + j_value**2
        ) == 64:
            same_hits.append((h_value, j_value))
        if j_value != -h_value and h_value * j_value * (
            h_value**2 - h_value * j_value + j_value**2
        ) == -64:
            opposite_hits.append((h_value, j_value))

payload = {
    "schema": "elkies-k3.rationalized-d6-section-chart-search.v1",
    "status": "PASS_EXACT_CHART_TWO_SECTION_OBSTRUCTION",
    "rational_d6_model": {
        "equation": (
            "y^2=x^3+u^2*(-u^2+a*u-3)*x+"
            "u^3*(g*u^3+b*u^2+a*u-2)"
        ),
        "fibre_profile_generic": "I2*+4I1",
        "root_lattice": "D6",
        "infinity_smooth_gate": "-4+27*g^2 != 0",
        "coefficients": {"a": str(a), "b": str(b), "g": str(g)},
        "section": {"x": str(x_section), "y": str(y_section)},
    },
    "two_section_leading_gate": {
        "height_bound": arguments.height_bound,
        "rational_values_tested": len(values),
        "ordered_pairs_tested": len(values) ** 2,
        "same_sign_equation": "h*j*(h^2+h*j+j^2)=64",
        "opposite_sign_equation": "h*j*(h^2-h*j+j^2)=-64",
        "nontrivial_same_sign_hits": [
            [str(left), str(right)] for left, right in same_hits
        ],
        "nontrivial_opposite_sign_hits": [
            [str(left), str(right)] for left, right in opposite_hits
        ],
        "exact_correspondence_obstruction": {
            "elliptic_curve": "Y^2=X^3+X^2+X",
            "same_sign_map": "X=j/h, x=h*j, Y=8*(j/h)/x",
            "opposite_sign_map": "X=-j/h, x=h*j, Y=8*(j/h)/x",
            "rank_over_QQ": int(correspondence_rank),
            "torsion_points": [
                "O" if point == correspondence_jacobian(0)
                else [str(point[0]), str(point[1])]
                for point in correspondence_torsion
            ],
            "nondegenerate_rational_points": 0,
        },
    },
    "proof_boundary": {
        "proved": (
            "The QQ D6 equation, generic fibre valuation, marked-section "
            "identity, correspondence factorization, exact rank-zero elliptic-curve "
            "obstruction for both nontrivial correspondence curves, and exhaustive "
            "stated height-box regression."
        ),
        "not_proved": (
            "A second rational section is impossible in a larger D6 section "
            "chart; or rank sum four is impossible on a D6 surface."
        ),
    },
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output_path = arguments.output.resolve()
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "D6QQ|surface=I2star+4I1|marked_sections=1|"
    f"height_box={arguments.height_bound}|same_hits={len(same_hits)}|"
    f"opposite_hits={len(opposite_hits)}|correspondence_rank=0|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
