#!/usr/bin/env sage -python
"""Resolve the formal III* E7 model and identify its component sections.

The formal E7 chart is ``Y^2=U^3+U*Z^3``. This script follows its ordinary
blow-up charts rather than using an E7 label as a substitute for a resolution.
It is a normal-form regression model only: the companion
``derive_h92_q6_actual_e7_resolution.sage`` proves that this coordinate chart
has not yet been transported from the exact H92 germ. The resulting graph
still identifies the formal E7_7 curve and formal affine component.
"""

import argparse
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e7-resolution.json"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()


def strict(value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


def hessian(value, point):
    return matrix(
        QQ,
        [[value.derivative(left, right)(*point) for right in (Z, U, Y)]
         for left in (Z, U, Y)],
    ).det()


f0 = Y**2 - U**3 - U * Z**3

# First blow-up.  The Z-chart carries the singular continuation.  The
# U-chart is smooth and contains the strict transform of the old fibre Z=0.
f1_z = strict(f0, (Z, Z * U, Z * Y), Z)
f1_u = strict(f0, (U * Z, U, U * Y), U)
assert f1_z == Y**2 - Z * U**3 - Z**2 * U
assert f1_u == Y**2 - U - U**2 * Z**3
assert f1_u.subs({Z: 0, U: Y**2}) == 0

# Second blow-up at the origin of the Z-chart.  Its Z-chart is an ordinary
# node, while its U-chart continues to the third blow-up.
f2_z = strict(f1_z, (Z, Z * U, Z * Y), Z)
f2_u = strict(f1_z, (U * Z, U, U * Y), U)
assert f2_z == Y**2 - Z**2 * U**3 - Z * U
assert f2_u == Y**2 - U**2 * Z - U * Z**2

# Third blow-up at the origin of f2_u.  There are three distinct ordinary
# nodes after this step: two endpoint nodes and one generic E3 node.  The
# apparent fourth node in the other chart is the overlap of the generic one.
f3_u = strict(f2_u, (U * Z, U, U * Y), U)
f3_z = strict(f2_u, (Z, Z * U, Z * Y), Z)
assert f3_u == Y**2 - U * Z - U * Z**2
assert f3_z == Y**2 - Z * U**2 - Z * U

nodes = (
    ("second_Z_node", f2_z, (QQ(0), QQ(0), QQ(0))),
    ("third_U_endpoint", f3_u, (QQ(0), QQ(0), QQ(0))),
    ("third_Z_endpoint", f3_z, (QQ(0), QQ(0), QQ(0))),
    ("third_generic", f3_u, (-QQ(1), QQ(0), QQ(0))),
)
for unused_name, equation, point in nodes:
    assert equation(*point) == 0
    assert all(equation.derivative(variable)(*point) == 0 for variable in (Z, U, Y))
    assert hessian(equation, point)

# The alternate chart sees the generic node at U=-1; it is the overlap of
# (Z,U)=(-1,0) in f3_u, not another exceptional curve.
assert f3_z(-QQ(0), -QQ(1), QQ(0)) == 0
assert hessian(f3_z, (QQ(0), -QQ(1), QQ(0)))

# Incidence from the displayed chart centers.  E1,E2,E3 are the first three
# exceptional curves.  Blowing the four nodes adds N2, N3u, N3z, N3g.  The
# map to the repository's E7 Cartan numbering is explicit below.
edges = (
    ("E1", "N3u"), ("N3u", "E3"), ("E3", "N3z"),
    ("N3z", "E2"), ("E2", "N2"), ("E3", "N3g"),
)
label_map = {
    "E1": "E7_1", "E2": "E7_2", "E3": "E7_3",
    "N3u": "E7_4", "N2": "E7_5", "N3g": "E7_6", "N3z": "E7_7",
}
expected_edges = {
    tuple(sorted(pair)) for pair in ((1, 4), (2, 5), (2, 7), (3, 4), (3, 6), (3, 7))
}
actual_edges = {
    tuple(sorted((int(label_map[left][-1]), int(label_map[right][-1]))) )
    for left, right in edges
}
assert actual_edges == expected_edges

payload = {
    "schema": "elkies-k3.h92-q6-e7-resolution.v1",
    "status": "PASS_FORMAL_E7_BLOWUP_CHARTS",
    "provenance": "Formal E7 normal form only; no explicit H92 coordinate transport is asserted.",
    "equation": "Y^2=U^3+U*Z^3",
    "charts": {
        "blow1_Z": str(f1_z), "blow1_U": str(f1_u),
        "blow2_Z": str(f2_z), "blow2_U": str(f2_u),
        "blow3_U": str(f3_u), "blow3_Z": str(f3_z),
    },
    "nodes": [name for name, unused_equation, unused_point in nodes],
    "exceptional_graph_edges": [list(edge) for edge in edges],
    "cartan_label_map": label_map,
    "component_sections": {
        "affine_E7": "strict transform Z=0 in blow1_U: U=Y^2",
        "E7_7": "N3z, the ordinary blow-up at the origin of blow3_Z",
        "reason": "These are distinct resolved curves; their D-intersections are certified separately by the q=6 component-section lattice artifact.",
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q6E7CHART|blowups=7|nodes=4|E7_7=N3z|affine=blow1_U|status=PASS_FORMAL_E7_BLOWUP_CHARTS", flush=True)
