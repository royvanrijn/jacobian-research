#!/usr/bin/env sage -python
"""Derive the complete E8 module for the source-nef H3 q=8 pencil.

The q=8 E8 target cycle, in the explicit chart-component order, is
``(-2,-4,-6,-10,-4,-7,-5,-8)``.  The q=6 marked module is ``u*<1,Q>`` and
``Q`` is a unit at the E8 singularity.  This script proves that the required
integral exceptional twist is the complete ideal ``(u^2,X,Y)`` in the actual
integral chart.  Hence the q=8 local module is exactly

    u^9 * (u^2, X, Y).

Its quotient has basis ``1,u``.  This supplies the E8 finite quotient target;
assembling it on a bounded global coefficient ambient remains separate.
"""

import argparse
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-local-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-complete-module.json"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

target = json.loads(args.target.read_text())
assert target["status"] == "PASS_EXACT_Q8_E8_SOURCE_TARGET"
assert tuple(target["q8"]["chart_exceptional_cycle"]) == (
    "-2", "-4", "-6", "-10", "-4", "-7", "-5", "-8",
)

# These are the exceptional valuations of u,X,Y in the chart order
# (B1,B2,B3,B4,N3,N40,N4B,N4inf).  Each is C^{-1}e_j for the strict transform
# met by that coordinate: u meets B1, X meets B2, and Y meets N4B.
u_value = vector(ZZ, (2, 2, 4, 6, 3, 4, 3, 5))
x_value = vector(ZZ, (2, 4, 6, 10, 4, 7, 5, 8))
y_value = vector(ZZ, (3, 5, 9, 15, 6, 10, 8, 12))
twist_cycle = -vector(ZZ, [ZZ(value) for value in target["q8"]["chart_exceptional_cycle"]])
assert twist_cycle == x_value

chart_cartan = matrix(ZZ, [
    [2, 0, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, 0, -1, 0, 0],
    [0, 0, 2, 0, -1, 0, 0, -1],
    [0, 0, 0, 2, 0, -1, -1, -1],
    [-1, 0, -1, 0, 2, 0, 0, 0],
    [0, -1, 0, -1, 0, 2, 0, 0],
    [0, 0, 0, -1, 0, 0, 2, 0],
    [0, 0, -1, -1, 0, 0, 0, 2],
])
assert chart_cartan * u_value == vector(ZZ, (1, 0, 0, 0, 0, 0, 0, 0))
assert chart_cartan * x_value == vector(ZZ, (0, 1, 0, 0, 0, 0, 0, 0))
assert chart_cartan * y_value == vector(ZZ, (0, 0, 0, 0, 0, 0, 1, 0))

# The E8 Weierstrass relation respects these valuations: each component has
# at least two terms of minimal order.  This makes the normal monomials
# u^a X^b Y^e (e=0,1) a valid valuation test set.
for index in range(8):
    orders = (
        2 * y_value[index], 3 * x_value[index],
        4 * u_value[index] + x_value[index], 5 * u_value[index],
        6 * u_value[index], 7 * u_value[index],
    )
    assert orders.count(min(orders)) >= 2

minimal_generators = []
for y_exponent in range(2):
    for u_exponent in range(6):
        for x_exponent in range(6):
            valuation = u_exponent * u_value + x_exponent * x_value + y_exponent * y_value
            if min(valuation - twist_cycle) < 0:
                continue
            candidate = (u_exponent, x_exponent, y_exponent)
            if any(
                generator[0] <= u_exponent
                and generator[1] <= x_exponent
                and generator[2] <= y_exponent
                for generator in minimal_generators
            ):
                continue
            minimal_generators.append(candidate)
assert minimal_generators == [(0, 1, 0), (2, 0, 0), (0, 0, 1)]

# The complete-ideal colength equals C.Cart an.C/2=2.  The displayed ideal
# has exactly the same quotient length, so the valuation-selected ideal has
# no hidden generator.
ring = PolynomialRing(QQ, names=("u", "X", "Y"), order="degrevlex")
u, X, Y = ring.gens()
relation = Y**2 - X**3 - u**4 * X - u**5
ideal = ring.ideal((relation, u**2, X, Y))
assert ideal.vector_space_dimension() == 2
assert ZZ(twist_cycle * chart_cartan * twist_cycle) // 2 == 2

payload = {
    "schema": "elkies-k3.h92-q8-e8-complete-module.v1",
    "status": "PASS_EXACT_Q8_E8_COMPLETE_MODULE",
    "inputs": {"e8_target": str(args.target.relative_to(ROOT))},
    "chart_component_order": target["chart_component_map"]["chart_component_order"],
    "coordinate_valuations": {
        "u": list(map(int, u_value)),
        "X": list(map(int, x_value)),
        "Y": list(map(int, y_value)),
    },
    "twist_cycle": list(map(int, twist_cycle)),
    "complete_ideal": {
        "generators": ["Y", "X", "u^2"],
        "quotient_basis": ["1", "u"],
        "colength": 2,
        "cycle_colength": 2,
    },
    "module": {
        "q6_ninth_power_local_unit_factor": "u^9",
        "q8_E8_module": "u^9*(u^2,X,Y)",
        "finite_quotient": "(u^9*R)/(u^9*(u^2,X,Y)) has basis 1,u",
    },
    "boundary": "This is the E8 local module only; no bounded global ambient or q8 kernel is asserted.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q8E8MODULE|ideal=(u2,X,Y)|colength=2|status=PASS_EXACT_Q8_E8_COMPLETE_MODULE", flush=True)
