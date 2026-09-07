#!/usr/bin/env sage -python
"""Derive actual sibling-chart transitions in the resolved H92 E7 atlas.

Two pairs of final E7 edge charts are the U- and Z-charts of the same ordinary
blow-up: E7_1--E7_4 with E7_4--E7_3, and E7_3--E7_7 with E7_7--E7_2.
Their overlap maps are part of the resolved geometry needed for a later
Čech/residual compiler block.  They are derived here from the actual H92
blow-up equations and old-coordinate pullbacks, not from Dynkin labels.

The emitted Cartier-factor ratios are transition functions, not assertions
that the two factors are units on the whole exceptional component.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-sibling-chart-transitions.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_expression(value, Z, U, Y):
    return sage_eval(str(value), locals={"Z": Z, "U": U, "Y": Y})


def normalized_fraction(field, value):
    """Cancel polynomial content for readable, deterministic overlap ratios."""
    ring = field.ring()
    numerator = ring(value.numerator())
    denominator = ring(value.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    scale = denominator.leading_coefficient()
    assert scale
    return field(numerator/scale) / field(denominator/scale)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

pullbacks = json.loads(args.pullbacks.read_text())
gluing = json.loads(args.gluing.read_text())
exec(compile(CORE.read_text(), str(CORE), "exec"))
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
charts = {item["name"]: item for item in pullbacks["charts"]}
twists = {item["name"]: item for item in gluing["actual_edge_chart_gluing"]}

# In source U-chart coordinates (Z,U,Y), the target Z-chart coordinates are
# (U*Z, 1/Z, Y/Z); its inverse is (1/U, Z*U, Y/U).  Both formulae are valid
# on the displayed principal open set and are derived from the shared blow-up,
# not guessed from the graph edge names.
pairs = (
    ("E7_1--E7_4", "E7_4--E7_3", "E7_4"),
    ("E7_3--E7_7", "E7_7--E7_2", "E7_7"),
)

records = []
for source_name, target_name, shared_component in pairs:
    source = charts[source_name]
    target = charts[target_name]
    source_twist = twists[source_name]
    target_twist = twists[target_name]
    ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
    Z, U, Y = ring.gens()
    field = ring.fraction_field()
    z, u, y = (field(value) for value in (Z, U, Y))
    forward = {"Z": u*z, "U": 1/z, "Y": y/z}
    inverse = {"Z": 1/u, "U": z*u, "Y": y/u}

    source_surface = field(rational_expression(source["surface_equation"], Z, U, Y))
    target_surface = field(rational_expression(target["surface_equation"], Z, U, Y))
    source_pullbacks = {}
    target_pullbacks = {}
    for coordinate in ("t", "x", "y"):
        source_pullbacks[coordinate] = field(rational_expression(
            source["old_coordinate_pullback"][coordinate], Z, U, Y
        ))
        target_pullbacks[coordinate] = field(rational_expression(
            target["old_coordinate_pullback"][coordinate], Z, U, Y
        ))

    source_g = field(rational_expression(
        source_twist["w_cartier_equation"], Z, U, Y
    ))
    target_g = field(rational_expression(
        target_twist["w_cartier_equation"], Z, U, Y
    ))
    transition = verify_resolved_chart_transition(
        "{} to {}".format(source_name, target_name),
        ring,
        source_surface,
        target_surface,
        (forward["Z"], forward["U"], forward["Y"]),
        source_pullbacks,
        target_pullbacks,
        source_g,
        target_g,
        "actual sibling blow-up charts and q8 Cartier factors",
    )
    surface_ratio = normalized_fraction(field, transition["surface_ratio"])
    cartier_ratio = normalized_fraction(field, transition["frame_ratio"])

    records.append({
        "source_chart": source_name,
        "target_chart": target_name,
        "shared_component": shared_component,
        "overlap": "source Z is nonzero; equivalently target U is nonzero",
        "target_coordinates_in_source": {key: str(value) for key, value in forward.items()},
        "source_coordinates_in_target": {key: str(value) for key, value in inverse.items()},
        "strict_surface_ratio_target_over_source": str(surface_ratio),
        "old_coordinate_pullbacks_in_source": {
            key: str(value) for key, value in transition["transported_pullbacks"].items()
        },
        "target_w_cartier_over_source_w_cartier": str(cartier_ratio),
        "transition_instruction": (
            "For a q8 local representative, transport the target chart "
            "trivialization through this rational map and multiply by the "
            "displayed Cartier-factor ratio before comparing residues."
        ),
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-sibling-chart-transitions.v1",
    "status": "PASS_EXACT_H92_Q8_E7_SIBLING_CHART_TRANSITIONS",
    "inputs": {
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "q8_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
    },
    "transitions": records,
    "boundary": (
        "These are the two sibling-chart overlaps among the final edge-node "
        "charts. Other E7 component overlap covers, finite residual quotients, "
        "and the complete Čech compatibility matrix remain to be derived."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E7SIBLINGTRANSITIONS|overlaps=2|status="
    "PASS_EXACT_H92_Q8_E7_SIBLING_CHART_TRANSITIONS",
    flush=True,
)
