#!/usr/bin/env sage -python
"""Attach the q=8 E7 integral twist to actual H92 resolved edge charts.

The exact q=8 E7 target satisfies

    c8 = 9*c6 + w,     w=(2,5,6,4,6,3,5).

Here ``c6`` is the non-Cartier q=6 marked cycle and ``w`` is an *integral*
exceptional divisor.  It is not anti-nef, so it cannot be replaced by a
complete ideal in the singular E7 germ.  On every actual H92 edge chart, if
``s=0`` and ``r=0`` cut the two listed exceptional components, the Cartier
equation for ``w`` is ``g=s^w_s*r^w_r``.  A q=8 local representative ``f``
is therefore compared with the ninth q=6 module by

    g*f in (q6_marked_module)^9.

This is the local line-bundle gluing direction for O(w): its local generator
is ``g^-1``.  The output supplies actual equations and factors for a later
resolved-chart evaluator; it deliberately does not assert a finite q=8
quotient or global kernel.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, vector


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-local-target.json"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

target = json.loads(args.target.read_text())
resolution = json.loads(args.resolution.read_text())
assert target["status"] == "PASS_EXACT_Q8_E7_LOCAL_TARGET"
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"
twist = vector(QQ, [QQ(value) for value in target["tensor_comparison"]["integral_exceptional_twist"]])
assert twist == vector(QQ, (2, 5, 6, 4, 6, 3, 5))
assert target["tensor_comparison"]["identity"] == "c8=9*c6+(2,5,6,4,6,3,5)"

# These component equations are transported from the actual, not formal, H92
# blow-up tree.  They are the same resolved edge coordinates used for the
# q=6 third-divisor Cartier calculation, but the sign below is opposite: w
# is added to the q6^9 cycle rather than subtracted from a divisor.
edge_data = (
    ("E7_2--E7_5", "E7_5", "Z", "E7_2", "Y"),
    ("E7_1--E7_4", "E7_4", "U", "E7_1", "Y"),
    ("E7_4--E7_3", "E7_4", "Z", "E7_3", "Y"),
    ("E7_3--E7_7", "E7_7", "U", "E7_3", "Y"),
    ("E7_7--E7_2", "E7_7", "Z", "E7_2", "Y"),
    ("E7_3--E7_6", "E7_6", "Z", "E7_3", "Y"),
)
assert set(resolution["edge_charts"]) == {entry[0] for entry in edge_data}

charts = []
for name, first_component, first_variable, second_component, second_variable in edge_data:
    first_order = int(twist[int(first_component[-1])-1])
    second_order = int(twist[int(second_component[-1])-1])
    cartier_equation = "{}^{}*{}^{}".format(
        first_variable, first_order, second_variable, second_order
    )
    charts.append({
        "name": name,
        "actual_h92_surface_equation": resolution["edge_charts"][name],
        "components": [
            {"name": first_component, "equation": first_variable, "w_coefficient": first_order},
            {"name": second_component, "equation": second_variable, "w_coefficient": second_order},
        ],
        "w_cartier_equation": cartier_equation,
        "O_w_local_generator": "({})^-1".format(cartier_equation),
        "comparison_condition": "({})*f belongs to (q6_marked_module)^9".format(cartier_equation),
    })

payload = {
    "schema": "elkies-k3.h92-q8-actual-e7-gluing.v1",
    "status": "PASS_EXACT_Q8_ACTUAL_E7_GLUING",
    "inputs": {
        "q8_e7_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
        "actual_h92_resolution": {"path": str(args.resolution.relative_to(ROOT)), "sha256": digest(args.resolution)},
    },
    "cycle_comparison": {
        "identity": "c8=9*c6+w",
        "integral_twist_w": [int(value) for value in twist],
        "orientation": (
            "w is added to the q6^9 exceptional cycle. Thus O(w) has local "
            "generator g^-1 and a q8 representative f is compared to the q6^9 "
            "module by g*f."
        ),
        "non_antinef_warning": target["tensor_comparison"]["non_antinef_obstruction"],
    },
    "actual_edge_chart_gluing": charts,
    "compiler_instruction": (
        "On each listed actual chart, evaluate a q8 ambient element f after "
        "multiplication by the displayed w_cartier_equation, then impose "
        "membership in the ninth q6 marked module with compatible chart gluing. "
        "Do not replace these six factors by one complete ideal downstairs."
    ),
    "boundary": (
        "This is an actual resolved-chart gluing specification for the q8 E7 "
        "integral twist. It does not evaluate the degree-18 basis, derive its "
        "finite chart quotients, prove a complete cover, or compute a pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8ACTUALE7GLUING|charts=6|twist=2,5,6,4,6,3,5|"
    "status=PASS_EXACT_Q8_ACTUAL_E7_GLUING",
    flush=True,
)
