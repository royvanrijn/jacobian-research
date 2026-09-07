#!/usr/bin/env sage -python
"""Derive formal-chart Cartier factors for the third q=6 E7 correction.

For the third marked class the resolved E7 vertical divisor is

    V_E7 = -22 E1 -44 E2 -66 E3 -44 E4 -33 E5 -33 E6 -55 E7.

This script does not turn that coefficient vector into a component-label
heuristic. It follows the ordinary blow-up charts of the *formal* model
``Y^2=U^3+U*Z^3`` through every E7 intersection and writes the local Cartier
factor in that model. For example, at E1--E4 the
surface chart is ``Y^2-Z-U*Z^2=0``; E4 is ``U=0`` and E1 is ``Y=0``, hence a
section of O(V_E7) must lie in ``(U^44*Y^22)`` in that local ring.

The resulting factors are a normal-form regression input only. They cannot
be used as H92 quotient blocks until an explicit transport from the exact H92
germ is exhibited. They also deliberately exclude the horizontal marked point
and smooth P2 condition.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-cartier-charts.json"
TARGET_SHA256 = "a0699e4ec75930cc93a9706ddf96f4ffc744954809e53a24029bd8c6668843f7"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.target == TARGET:
    assert digest(args.target) == TARGET_SHA256
target = json.loads(args.target.read_text())
assert target["status"] == "PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET"
orders = -vector(ZZ, target["resolved_exceptional_coefficients"])
assert orders == vector(ZZ, (22, 44, 66, 44, 33, 33, 55))

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()


def strict(value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


f0 = Y**2 - U**3 - U * Z**3
f1_z = strict(f0, (Z, Z * U, Z * Y), Z)
f2_z = strict(f1_z, (Z, Z * U, Z * Y), Z)
f2_u = strict(f1_z, (U * Z, U, U * Y), U)
f3_u = strict(f2_u, (U * Z, U, U * Y), U)
f3_z = strict(f2_u, (Z, Z * U, Z * Y), Z)
assert f2_z == Y**2 - Z**2 * U**3 - Z * U
assert f3_u == Y**2 - U * Z - U * Z**2
assert f3_z == Y**2 - Z * U**2 - Z * U

# Each record is an honest smooth chart after the indicated node blow-up.
# ``exceptional_variable`` cuts the newly introduced component and
# ``old_component_variable`` cuts its adjacent strict transform.  The surface
# equation has a nonzero coordinate derivative at the origin, so the displayed
# component equations are regular local parameters after eliminating the
# remaining coordinate (which differs from chart to chart).
records = []


def record(name, equation, exceptional_component, exceptional_variable,
           old_component, old_component_variable):
    assert equation(0, 0, 0) == 0
    assert any(equation.derivative(variable)(0, 0, 0) in (QQ(1), -QQ(1))
               for variable in (Z, U, Y))
    exceptional_index = int(exceptional_component[-1]) - 1
    old_index = int(old_component[-1]) - 1
    exceptional_order = int(orders[exceptional_index])
    old_order = int(orders[old_index])
    records.append({
        "name": name,
        "surface_equation": str(equation),
        "components": {
            "exceptional": {"name": exceptional_component, "equation": exceptional_variable,
                            "required_vanishing_order": exceptional_order},
            "strict_transform": {"name": old_component, "equation": old_component_variable,
                                 "required_vanishing_order": old_order},
        },
        "cartier_factor": "{}^{}*{}^{}".format(
            exceptional_variable, exceptional_order, old_component_variable, old_order
        ),
        "membership_condition": "regular representative belongs to ({})".format(
            "{}^{}*{}^{}".format(
                exceptional_variable, exceptional_order, old_component_variable, old_order
            )
        ),
    })


# Blow up the f3_U endpoint node.  Its U-chart contains E1--E4 and its
# Z-chart contains E4--E3.
n3u_u = strict(f3_u, (U * Z, U, U * Y), U)
n3u_z = strict(f3_u, (Z, Z * U, Z * Y), Z)
assert n3u_u == Y**2 - Z - U * Z**2
assert n3u_z == Y**2 - U - Z * U
record("E7_1--E7_4", n3u_u, "E7_4", "U", "E7_1", "Y")
record("E7_4--E7_3", n3u_z, "E7_4", "Z", "E7_3", "Y")

# Blow up the f3_Z endpoint node.  Its U-chart contains E3--E7 and its
# Z-chart contains E7--E2.
n3z_u = strict(f3_z, (U * Z, U, U * Y), U)
n3z_z = strict(f3_z, (Z, Z * U, Z * Y), Z)
assert n3z_u == Y**2 - Z - U * Z
assert n3z_z == Y**2 - U - Z * U**2
record("E7_3--E7_7", n3z_u, "E7_7", "U", "E7_3", "Y")
record("E7_7--E7_2", n3z_z, "E7_7", "Z", "E7_2", "Y")

# Blow up the f2_Z endpoint node; its U-chart contains E2--E5.
n2_u = strict(f2_z, (U * Z, U, U * Y), U)
assert n2_u == Y**2 - Z - U**3 * Z**2
record("E7_2--E7_5", n2_u, "E7_5", "U", "E7_2", "Y")

# At the generic f3_U node put W=Z+1.  The W-chart of its blow-up contains
# E3--E6.  Substituting W=Z in the translated germ yields the displayed
# actual chart equation.
generic = ring(f3_u(Z - 1, U, Y))
assert generic == Y**2 + U * Z - U * Z**2
n3g_w = strict(generic, (Z, Z * U, Z * Y), Z)
assert n3g_w == Y**2 + U - Z * U
record("E7_3--E7_6", n3g_w, "E7_6", "Z", "E7_3", "Y")

assert {record_value["name"] for record_value in records} == {
    "E7_1--E7_4", "E7_4--E7_3", "E7_3--E7_7",
    "E7_7--E7_2", "E7_2--E7_5", "E7_3--E7_6",
}

payload = {
    "schema": "elkies-k3.h92-q6-third-e7-cartier-charts.v1",
    "status": "PASS_FORMAL_Q6_THIRD_E7_CARTIER_CHARTS",
    "equation": "Y^2=U^3+U*Z^3",
    "resolved_vanishing_orders": [int(value) for value in orders],
    "charts": records,
    "compiler_instruction": (
        "For each displayed formal chart, first apply the stated line-bundle "
        "trivialization (the inverse of the Cartier factor), then reduce the "
        "ambient basis in the corresponding finite local quotient. These are "
        "only the integral E7 vertical conditions; add smooth marked-point, E8, "
        "and finite-base blocks before asserting h0."
    ),
    "boundary": (
        "This is a formal-model Cartier atlas for V_E7, not a transported H92 "
        "atlas. It does not evaluate the degree-44 marked chord or build a finite "
        "matrix, so it does not claim a kernel dimension."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDE7CARTIER|charts=6|orders=22,44,66,44,33,33,55|"
    "status=PASS_FORMAL_Q6_THIRD_E7_CARTIER_CHARTS",
    flush=True,
)
