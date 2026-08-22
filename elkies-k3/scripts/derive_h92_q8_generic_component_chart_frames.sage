#!/usr/bin/env sage -python
"""Pin actual H92 generic E7 chart frames for the q=8 residue compiler.

The q=8 all-component valuation template says which ambient terms have a
negative generic order.  To evaluate a non-singleton group, one must also
know an *actual blow-up chart*, its reduced component equation, and a normal
weight vector reproducing the transported orders of ``t,x,y``.  This script
extracts and verifies precisely that input from the six H92 chart pullbacks.

For a component cut by Z or U the normal weight is the corresponding unit
vector.  For a component cut by Y, the resolved surface makes the adjacent
coordinate quadratic in Y; the displayed weights record that actual normal
branch.  They are not inferred from an E7 Kodaira diagram.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-component-chart-frames.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def weighted_order(value, weights):
    """Return the order of a nonzero polynomial for a positive weight vector."""
    monomials = value.dict()
    assert monomials
    return min(sum(exponent * weight for exponent, weight in zip(monomial, weights))
               for monomial in monomials)


# Each selection is transported from an actual edge chart.  ``Y`` rows carry
# the local relation that makes one neighbouring coordinate quadratic in Y.
SELECTIONS = {
    "E7_1": {
        "chart": "E7_1--E7_4", "equation": "Y", "weights": (2, 0, 1),
        "normal_branch": "Z is a unit times Y^2; U is a generic unit coordinate",
    },
    "E7_2": {
        "chart": "E7_2--E7_5", "equation": "Y", "weights": (0, 2, 1),
        "normal_branch": "U is a unit times Y^2; Z is a generic unit coordinate",
    },
    "E7_3": {
        "chart": "E7_4--E7_3", "equation": "Y", "weights": (0, 2, 1),
        "normal_branch": "U is a unit times Y^2; Z is a generic unit coordinate",
    },
    "E7_4": {
        "chart": "E7_4--E7_3", "equation": "Z", "weights": (1, 0, 0),
        "normal_branch": "Z is the reduced normal parameter",
    },
    "E7_5": {
        "chart": "E7_2--E7_5", "equation": "Z", "weights": (1, 0, 0),
        "normal_branch": "Z is the reduced normal parameter",
    },
    "E7_6": {
        "chart": "E7_3--E7_6", "equation": "Z", "weights": (1, 0, 0),
        "normal_branch": "Z is the reduced normal parameter",
    },
    "E7_7": {
        "chart": "E7_3--E7_7", "equation": "U", "weights": (0, 1, 0),
        "normal_branch": "U is the reduced normal parameter",
    },
}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--atlas", type=Path, default=ATLAS)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

pullbacks = json.loads(args.pullbacks.read_text())
atlas = json.loads(args.atlas.read_text())
gluing = json.loads(args.gluing.read_text())
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert atlas["status"] == "PASS_EXACT_H92_E7_VALUATION_ATLAS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
chart_data = {entry["name"]: entry for entry in pullbacks["charts"]}
edge_data = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
atlas_data = {entry["component"]: entry for entry in atlas["entries"]}
assert set(SELECTIONS) == set(atlas_data)

records = []
for component in tuple("E7_{}".format(index) for index in range(1, 8)):
    selection = SELECTIONS[component]
    chart = chart_data[selection["chart"]]
    edge = edge_data[selection["chart"]]
    assert any(
        item["name"] == component and item["equation"] == selection["equation"]
        for item in edge["components"]
    )
    pullback = {
        name: ring(sage_eval(expression, locals={"Z": Z, "U": U, "Y": Y}))
        for name, expression in chart["old_coordinate_pullback"].items()
    }
    weights = tuple(selection["weights"])
    orders = {name: int(weighted_order(value, weights)) for name, value in pullback.items()}
    expected = {name: int(atlas_data[component]["old_coordinate_orders"][name])
                for name in ("t", "x", "y")}
    assert orders == expected
    surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
    records.append({
        "component": component,
        "actual_edge_chart": selection["chart"],
        "component_equation": selection["equation"],
        "normal_weights_Z_U_Y": list(weights),
        "normal_branch": selection["normal_branch"],
        "surface_equation": chart["surface_equation"],
        "surface_initial_weight": int(weighted_order(surface, weights)),
        "old_coordinate_pullback": chart["old_coordinate_pullback"],
        "old_coordinate_weight_orders": orders,
    })

payload = {
    "schema": "elkies-k3.h92-q8-generic-component-chart-frames.v1",
    "status": "PASS_EXACT_Q8_GENERIC_COMPONENT_CHART_FRAMES",
    "inputs": {
        "actual_chart_pullbacks": {"path": path_label(args.pullbacks), "sha256": digest(args.pullbacks)},
        "actual_valuation_atlas": {"path": path_label(args.atlas), "sha256": digest(args.atlas)},
        "q8_actual_gluing": {"path": path_label(args.gluing), "sha256": digest(args.gluing)},
    },
    "component_frames": records,
    "compiler_instruction": (
        "For each non-singleton negative-order group, substitute its old-model "
        "term in the selected chart, divide by the recorded normal order, then "
        "reduce its leading residue in the component function field using the "
        "displayed actual surface equation. Add the resulting exact relation "
        "before any edge-node or overlap condition."
    ),
    "boundary": (
        "This pins and verifies actual generic-component chart frames only. It "
        "does not yet evaluate marked-chord leading residues, form their linear "
        "relations, test nodes or overlaps, or certify a q8 global kernel."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8COMPONENTFRAMES|components=7|charts={}|status="
    "PASS_EXACT_Q8_GENERIC_COMPONENT_CHART_FRAMES".format(
        len(set(record["actual_edge_chart"] for record in records))
    ),
    flush=True,
)
