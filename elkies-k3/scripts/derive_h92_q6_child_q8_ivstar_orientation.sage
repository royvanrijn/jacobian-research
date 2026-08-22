#!/usr/bin/env sage -python
"""Attach the physical E6 q8 target to its IV* resolved chart arm.

The IV* vertical-ideal-pair checker leaves only the E6 arm involution open.
The transported old ``E7_7`` component removes it: its source NS class meets
physical E6 root five, while its exact standard-Jacobian jet has
``(X/u^2,Y/u^2)`` finite with ``Y/u^2=c``.  In the ordinary IV* resolution
this is the ``plus_outer`` branch.  Therefore physical root five is the
plus outer component, selecting the second candidate from the pair checker.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, EllipticCurve, PolynomialRing, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
COMPONENTS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-e7-infinity-sections.json"
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-physical-root-target.json"
PAIR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-vertical-ideal.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-orientation.json"
REFLECTIONS = (1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational_function(field, ring, entry, numerator_key, denominator_key):
    return field(polynomial(ring, entry[numerator_key])) / field(polynomial(ring, entry[denominator_key]))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--components", type=Path, default=COMPONENTS)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--pair", type=Path, default=PAIR)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "components", "target", "pair", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
components = json.loads(args.components.read_text())
target = json.loads(args.target.read_text())
pair = json.loads(args.pair.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert pair["status"] == "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_VERTICAL_IDEAL_PAIR"

# First identify the old E7_7 curve in the physical E6 lattice.  This is a
# source-side intersection computation, independent of its Weierstrass jet.
ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -load_gram(FRAME))
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)


def reflect(class_value, nodes):
    result = vector(ZZ, class_value)
    for node in nodes:
        curve = simple[node - 1]
        result += (result * ns * curve) * curve
    return result


raw_e7_7 = reflect(simple[6], tuple(reversed(REFLECTIONS)))
e6_roots = matrix(ZZ, target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"])
physical_intersections = raw_e7_7 * ns * e6_roots.transpose()
assert physical_intersections == vector(ZZ, (0, 0, 0, 0, 1, 0))

# Then compute its exact IV* entrance branch from the independently exported
# child section coordinate.
T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
base_field = T_ring.fraction_field()
A = polynomial(T_ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(T_ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
curve = EllipticCurve(base_field, [0, 0, 0, A, B])
e7_7_data = next(
    entry for entry in components["sections"]
    if entry["sign"] == components["source"]["E7_7_sign"]
)
e7_7 = curve(
    rational_function(base_field, T_ring, e7_7_data,
                      "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
    rational_function(base_field, T_ring, e7_7_data,
                      "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
)
iv_star = next(fibre for fibre in child["finite_fibres"] if fibre["kodaira"] == "IV*")
factor = T_ring(iv_star["factor"])
assert factor.degree() == 1 and tuple(iv_star["minimal_orders"]) == (3, 4, 8)
base_point = -factor[0] / factor[1]
u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
x, y = e7_7.xy()


def translate(value):
    return u_field(u_ring(value.numerator()(base_point + u))) / u_field(
        u_ring(value.denominator()(base_point + u))
    )


x_u, y_u = translate(x), translate(y)
assert x_u.valuation() == y_u.valuation() == 2
A_u = u_ring(A(base_point + u))
B_u = u_ring(B(base_point + u))
a, remainder = A_u.quo_rem(u**3)
assert not remainder and a(0)
b, remainder = B_u.quo_rem(u**4)
assert not remainder and b(0)
c = QQ(b(0)).sqrt()
assert c and c**2 == b(0)
X_entrance, Y_entrance = x_u / u**2, y_u / u**2
assert X_entrance.denominator()(0) and Y_entrance.denominator()(0)
assert Y_entrance(0) == c

# The pair checker's component order defines chart component two to be this
# plus branch.  Since root five is the physical component met above, its map
# must be the second of the two E6 diagram isomorphisms.
selected = next(
    candidate for candidate in pair["orientation_candidates"]
    if candidate["physical_E6_i_to_chart_component"] == [3, 5, 6, 4, 2, 1]
)
assert selected["generators"] == ["Y+c*u^2", "u*X", "X^2", "u^3"]
assert selected["chart_component_degrees"] == [0, 0, -1, 0, -1, 0]
assert selected["chart_cycle"] == [3, 2, 3, 4, 5, 6]

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-ivstar-orientation.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_ORIENTATION",
    "inputs": {
        "source_frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": digest(FRAME)},
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "child_e7_components": {"path": str(args.components.relative_to(ROOT)), "sha256": digest(args.components)},
        "physical_root_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
        "ivstar_ideal_pair": {"path": str(args.pair.relative_to(ROOT)), "sha256": digest(args.pair)},
    },
    "orientation_witness": {
        "section": "transported old E7_7",
        "physical_E6_component_intersections": list(map(int, physical_intersections)),
        "resolved_branch": "plus_outer (Y/u^2=c)",
        "entrance_coordinates": {"X_over_u2_at_u0": str(X_entrance(0)), "Y_over_u2_at_u0": str(Y_entrance(0))},
    },
    "selected_q8_ivstar_module": {
        "physical_E6_i_to_chart_component": selected["physical_E6_i_to_chart_component"],
        "chart_cycle": selected["chart_cycle"],
        "chart_component_degrees": selected["chart_component_degrees"],
        "complete_ideal_generators": selected["generators"],
        "quotient_basis": pair["common_quotient_basis"],
        "colength": 4,
    },
    "boundary": (
        "This selects the IV* vertical ideal only. The generic-chord trivialization, "
        "combination with II* and smooth conditions, q8 pencil, rootless equation, "
        "bisection covers, extension collisions, and generic rank 18 or 19 remain unproved."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8IVORIENT|physical_root=5|branch=plus|ideal=(Y+c*u2,uX,X2,u3)|"
    "status=PASS_EXACT_Q6_CHILD_Q8_IVSTAR_ORIENTATION",
    flush=True,
)
