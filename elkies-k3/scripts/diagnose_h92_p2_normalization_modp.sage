#!/usr/bin/env sage -python
"""Compare the plane-cubic hyperplane class with the certified H92 P1.

The degree-10 residual construction gives ``D+R ~ 10H`` for the marked
degree-29 divisor ``D`` and plane-hyperplane class ``H``.  With ``R`` as
origin, its tangent class is ``C=H-3R``.  Consequently the intrinsic
degree-zero class is exactly

    3D-29H = C.

This diagnostic evaluates ``C`` and the independently certified H92 ``P1``
on common good finite fibers, looking only for a small multiple relation.  It
does *not* identify ``C`` with ``P2``; a negative result is evidence that the
residual-origin normalization cannot simply be relabelled as the known
height-21/2 direction.
"""

from sage.all import *

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
CANDIDATE = ROOT / "artifacts/generated-results/h92-p2-candidate-mod-100003.json"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"


def evaluate(coefficients, value, field):
    return sum(field(QQ(coefficient)) * value**index for index, coefficient in enumerate(coefficients))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--candidate", type=Path, default=CANDIDATE)
parser.add_argument("--multiple-bound", type=int, default=12)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
candidate = json.loads(args.candidate.read_text())
prime = ZZ(candidate["prime"])
field = GF(prime)
anchor = SourceFileLoader("h92_p2_normalization_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (field(QQ(value(r, s))) for value in formulas)


def p1_at(target_t):
    u = 1 / target_t
    x = evaluate(p1["x_entrance_base"]["numerator_coefficients"], u, field)
    x /= evaluate(p1["x_entrance_base"]["denominator_coefficients"], u, field)
    y = evaluate(p1["y_entrance_base"]["numerator_coefficients"], u, field)
    y /= evaluate(p1["y_entrance_base"]["denominator_coefficients"], u, field)
    curve = EllipticCurve([
        0, 0, 0,
        A1 * target_t**3 + A * target_t**4,
        B1 * target_t**5 + B * target_t**6 + B2 * target_t**7,
    ])
    return curve(x, y)


common_multiples = None
tested = 0
for integer, x_value, y_value in candidate["values"]:
    target_t = field(integer)
    try:
        point_p1 = p1_at(target_t)
    except (ArithmeticError, ZeroDivisionError, ValueError):
        continue
    point_candidate = point_p1.curve()(field(x_value), field(y_value))
    matches = {
        multiple
        for multiple in range(-args.multiple_bound, args.multiple_bound + 1)
        if multiple * point_p1 == point_candidate
    }
    common_multiples = matches if common_multiples is None else common_multiples & matches
    tested += 1
    if tested >= 12:
        break

assert tested >= 3
# From D+R ~ 10H, D-29R = 10(H-3R), so
# 3(D-29R)-29(H-3R) = H-3R.  Thus the sampled tangent class is the intrinsic
# degree-zero combination 3D-29H, rather than the nonzero-degree D-7H.
print(
    "H92P2NORMALIZATION|"
    f"prime={prime}|fibers={tested}|multiple_bound={args.multiple_bound}|"
    f"canonical_class=3D-29H|common_small_multiples={sorted(common_multiples)}"
)
