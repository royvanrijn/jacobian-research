#!/usr/bin/env sage -python
"""Resolve the exact H92 E7 fibre through its actual ordinary blow-up tree.

The standard E7 normal form has a node in its second Z-chart.  For the exact
H92 germ that chart is smooth; its fourth node instead occurs in the second
U-chart at the nonzero coordinate ``Z=-A1/B1``.  This script follows that
actual tree:

* blow the origin, then its Z-chart origin, then the second-U origin;
* blow the separate second-U node at ``Z=-A1/B1``;
* blow the three ordinary nodes after the third blow-up.

It retains every exact H92 coefficient and records the six final edge charts.
Those are the coordinate charts required for a genuine H92 E7 divisor
trivialization and replace the untransported normal-form model for future
matrix work.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


def hessian(value, variables, point):
    return matrix(
        QQ,
        [[value.derivative(left, right)(*point) for right in variables]
         for left in variables],
    ).det()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

anchor = SourceFileLoader("h92_actual_e7_full_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
assert A1 and B1

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
variables = (Z, U, Y)
f0 = Y**2 - U**3 - (A1 * Z**3 + A * Z**4) * U - (
    B1 * Z**5 + B * Z**6 + B2 * Z**7
)
f1_z = strict(ring, f0, (Z, Z * U, Z * Y), Z)
f2_z = strict(ring, f1_z, (Z, Z * U, Z * Y), Z)
f2_u = strict(ring, f1_z, (U * Z, U, U * Y), U)
f3_u = strict(ring, f2_u, (U * Z, U, U * Y), U)
f3_z = strict(ring, f2_u, (Z, Z * U, Z * Y), Z)

assert f2_z.derivative(Z)(0, 0, 0) == -B1
assert f2_u == Y**2 - U**2 * Z - A1 * U * Z**2 - A * U**2 * Z**3 - B1 * U * Z**3 - B * U**2 * Z**4 - B2 * U**3 * Z**5
assert f3_u == Y**2 - U * Z - A1 * U * Z**2 - A * U**3 * Z**3 - B1 * U**2 * Z**3 - B * U**4 * Z**4 - B2 * U**6 * Z**5
assert f3_z == Y**2 - Z * U**2 - A1 * Z * U - A * Z**3 * U**2 - B1 * Z**2 * U - B * Z**4 * U**2 - B2 * Z**6 * U**3

second_u_node = -A1 / B1
third_generic_node = -QQ(1) / A1
nodes = (
    ("second_U_nonzero", f2_u, (second_u_node, QQ(0), QQ(0))),
    ("third_U_endpoint", f3_u, (QQ(0), QQ(0), QQ(0))),
    ("third_Z_endpoint", f3_z, (QQ(0), QQ(0), QQ(0))),
    ("third_U_generic", f3_u, (third_generic_node, QQ(0), QQ(0))),
)
for unused_name, equation, point in nodes:
    assert equation(*point) == 0
    assert all(equation.derivative(variable)(*point) == 0 for variable in variables)
    assert hessian(equation, variables, point)

# The f2_U nonzero node is translated to the origin before its node blow-up.
# Its Z-chart contains the E2--N2 edge because the strict transform of U=0
# is U=0 after substituting U=Z*U.
second_translated = ring(f2_u(Z + second_u_node, U, Y))
n2_z = strict(ring, second_translated, (Z, Z * U, Z * Y), Z)

# The two endpoint nodes after blow-up three give E1--N3u/E3--N3u and
# E3--N3z/E2--N3z respectively.
n3u_u = strict(ring, f3_u, (U * Z, U, U * Y), U)
n3u_z = strict(ring, f3_u, (Z, Z * U, Z * Y), Z)
n3z_u = strict(ring, f3_z, (U * Z, U, U * Y), U)
n3z_z = strict(ring, f3_z, (Z, Z * U, Z * Y), Z)

# Translate the third generic node before taking its Z-chart; it contains the
# E3--N3g edge.
third_translated = ring(f3_u(Z + third_generic_node, U, Y))
n3g_z = strict(ring, third_translated, (Z, Z * U, Z * Y), Z)

edge_charts = {
    "E7_2--E7_5": n2_z,
    "E7_1--E7_4": n3u_u,
    "E7_4--E7_3": n3u_z,
    "E7_3--E7_7": n3z_u,
    "E7_7--E7_2": n3z_z,
    "E7_3--E7_6": n3g_z,
}
for name, equation in edge_charts.items():
    assert equation(0, 0, 0) == 0
    assert any(equation.derivative(variable)(0, 0, 0) for variable in variables), name

payload = {
    "schema": "elkies-k3.h92-q6-actual-e7-resolution-full.v1",
    "status": "PASS_EXACT_H92_E7_BLOWUP_TREE",
    "inputs": {"h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)}},
    "equation": str(f0),
    "coefficients": {name: str(value) for name, value in zip(("A1", "A", "B1", "B", "B2"), (A1, A, B1, B, B2))},
    "intermediate_charts": {"blow1_Z": str(f1_z), "blow2_Z": str(f2_z), "blow2_U": str(f2_u), "blow3_U": str(f3_u), "blow3_Z": str(f3_z)},
    "smooth_standard_blow2_Z_derivative": str(f2_z.derivative(Z)(0, 0, 0)),
    "nodes": [{"name": name, "point": [str(value) for value in point]} for name, unused_equation, point in nodes],
    "edge_charts": {name: str(equation) for name, equation in edge_charts.items()},
    "component_origin": {
        "E7_1": "first exceptional", "E7_2": "second exceptional", "E7_3": "third exceptional",
        "E7_4": "third_U_endpoint node", "E7_5": "second_U_nonzero node",
        "E7_6": "third_U_generic node", "E7_7": "third_Z_endpoint node",
    },
    "boundary": "These are actual H92 blow-up charts and graph incidences. Pulling the high-degree marked chord and the integral vertical divisor into their local rings remains the next condition-matrix step.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92ACTUALE7FULL|nodes=4|blowups=7|second_U_Z={}|third_U_Z={}|status=PASS_EXACT_H92_E7_BLOWUP_TREE".format(
        second_u_node, third_generic_node
    ),
    flush=True,
)
