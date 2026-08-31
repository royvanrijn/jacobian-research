#!/usr/bin/env sage
"""Evaluate the explicit second Q80 q=4 pencil at one exact curve point.

This is an exact rational witness probe for the second-neighbor discriminant.
It is not a generic-family certificate until the predicted open conditions
are promoted to QQ(u), as done for the first q=4 step.
"""

import argparse
import json
from pathlib import Path

from sage.all import *
from sage.misc.persist import load


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--parameter",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-PDQE-parameter.json",
)
parser.add_argument("--u", type=int, default=0)
args = parser.parse_args()

load(str(ROOT / "elkies-k3/scripts/derive_q80_second_q4_pencil.sage"))

payload = json.loads(args.parameter.read_text())
assert payload["schema"] == "q80-cm24-qq-PDQE-parameter-v1"
parameter_ring = PolynomialRing(QQ, "u")
u = parameter_ring.gen()
parameter_field = parameter_ring.fraction_field()


def parameter_function(name):
    record = payload["original_functions"][name]
    numerator = parameter_ring(record["numerator"])
    denominator = parameter_ring(record["denominator"])
    assert numerator.gcd(denominator) == 1
    return parameter_field(numerator / denominator)


point = QQ(args.u)
values = tuple(
    QQ(parameter_function(name)(point)) for name in ("d", "p", "q", "e")
)
map_to_qq = parameters.hom(values, QQ)


def specialize(value):
    value = K(value)
    denominator = map_to_qq(parameters(value.denominator()))
    if denominator == 0:
        raise ZeroDivisionError("chosen parameter point hits a Q80 chart denominator")
    return QQ(map_to_qq(parameters(value.numerator())) / denominator)


W_ring = PolynomialRing(QQ, "W")
W = W_ring.gen()
v_ring = PolynomialRing(W_ring, "v")
v = v_ring.gen()
curve = v_ring([
    W_ring([specialize(value) for value in KW(coefficient).list()])
    for coefficient in second_curve.list()
])
assert curve.degree(v) == 3
q0, q1, q2, q3, q4 = [curve[index] for index in range(5)]
invariant_i = 12*q4*q0-3*q3*q1+q2**2
invariant_j = (
    72*q4*q2*q0+9*q3*q2*q1-27*q4*q1**2
    -27*q3**2*q0-2*q2**3
)
child_a = -27*invariant_i
child_b = -27*invariant_j
child_delta = 4*child_a**3+27*child_b**2
factor_data = tuple(
    (int(factor.degree()), int(exponent),
     int(child_a.valuation(factor)), int(child_b.valuation(factor)))
    for factor, exponent in child_delta.factor()
)
infinity = (
    int(8-child_a.degree()),
    int(12-child_b.degree()),
    int(24-child_delta.degree()),
)
assert sorted(factor_data) == [(1, 7, 2, 3), (8, 1, 0, 0)]
assert infinity == (2, 3, 9)
print(
    "Q80UNMARKEDSECONDQ4PROBE|u={}|curve_degree={}|"
    "finite=(I1*,8I1)|infinity=I3*|status=PASS_EXACT_WITNESS".format(
        point, curve.degree(v)
    ),
    flush=True,
)
