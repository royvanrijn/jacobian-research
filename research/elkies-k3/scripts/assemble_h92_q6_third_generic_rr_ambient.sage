#!/usr/bin/env sage -python
"""Construct the exact symbolic generic ambient for the third q=6 divisor.

The Weyl certificate gives the old generic restriction

    43*O + P,     P = 22*(-P1)-P2.

There is no mathematical need to expand this high-height point into a single
reduced pair of rational functions before imposing local conditions. Such an
expansion makes Sage repeatedly take enormous gcds in ``QQ(t)`` and hides the
marked-divisor structure needed on a resolution. The compiler instead stores
the point as an exact group-law expression DAG over already certified P1/P2
coordinates. A chart evaluator must evaluate that DAG in *its own* resolved
local ring, where cancellation is controlled by the chosen trivialization.

On the old generic fibre, the 43 standard Weierstrass monomials span
``L(43*O)``. If ``Q=-P``, the exact marked chord
``m_Q=(y-y(Q))/(x-x(Q))`` has one pole at O and one at P; adjoining it gives
the 44-dimensional space ``L(43*O+P)``. This is an ambient certificate, not
a completed resolved-condition matrix or a child equation.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
P2 = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json"
THIRD_E7 = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-generic-rr-ambient.json"
P1_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
P2_SHA256 = "e02e2803387d3a7f53907f548b275bb592d366f653f630f6ba8c9ef2611f3e37"
THIRD_E7_SHA256 = "a0699e4ec75930cc93a9706ddf96f4ffc744954809e53a24029bd8c6668843f7"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--p2", type=Path, default=P2)
parser.add_argument("--third-e7", type=Path, default=THIRD_E7)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.p1 == P1:
    assert digest(args.p1) == P1_SHA256
if args.p2 == P2:
    assert digest(args.p2) == P2_SHA256
if args.third_e7 == THIRD_E7:
    assert digest(args.third_e7) == THIRD_E7_SHA256
p1 = json.loads(args.p1.read_text())
p2 = json.loads(args.p2.read_text())
third_e7 = json.loads(args.third_e7.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1" and p1["exact_weierstrass_square"]
assert p2["schema"] == "elkies-k3.h92-p2-hensel-lift.v1" and p2["complete"]
assert third_e7["status"] == "PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET"
assert third_e7["generic_restriction"] == "43*O + (22*(-P1)-P2)"

# Every reference is explicit about base coordinates and signs. P1 uses the
# reciprocal entrance coordinate u=1/t; P2 is reconstructed in t. A resolved
# evaluator is required to make that substitution before applying this tree.
point_dag = {
    "operation": "add",
    "left": {"operation": "scalar", "scalar": 22,
             "point": {"operation": "negate", "point": "P1"}},
    "right": "reconstructed_-P2",
}
assert p1["x_entrance_base"]["degrees"] == [10, 12]
assert p1["y_entrance_base"]["degrees"] == [15, 18]
assert len(p2["X"]) == 47 and len(p2["Y"]) == 70 and len(p2["Z"]) == 22
# In the original (P1,P2) frame the marked point has vector (-22,-1).
# Its positive height certifies it is neither zero nor 2-torsion, so the
# marked chord has the asserted simple distinct pole at P.
source_height = matrix(QQ, [[QQ(21) / 2, 3], [3, 46]])
point_height = vector(QQ, (-22, -1)) * source_height * vector(QQ, (-22, -1))
assert point_height == 5260

basis = []
for b in (0, 1):
    for a in range((43 - 3 * b) // 2 + 1):
        basis.append({"kind": "monomial", "x_power": a, "y_power": b,
                      "pole_order_at_O": 2 * a + 3 * b})
assert len(basis) == 43
basis.append({
    "kind": "marked_chord",
    "point": {"operation": "negate", "point": point_dag},
    "formula": "(y-y(-P))/(x-x(-P))",
    "pole_order_at_O": 1,
    "pole_order_at_P": 1,
})
assert len(basis) == 44

payload = {
    "schema": "elkies-k3.h92-q6-third-generic-rr-ambient.v2",
    "status": "PASS_EXACT_Q6_THIRD_SYMBOLIC_GENERIC_RR_AMBIENT",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "p2": {"path": str(args.p2.relative_to(ROOT)), "sha256": digest(args.p2)},
        "third_e7_target": {"path": str(args.third_e7.relative_to(ROOT)), "sha256": digest(args.third_e7)},
    },
    "generic_fibre_divisor": "43*O + P, P=22*(-P1)-P2",
    "generic_fibre_degree": 44,
    "coordinate_orientation": "reconstructed_p2=-P2 in the H3 frame",
    "point_height_in_source_frame": str(point_height),
    "point_expression_dag": point_dag,
    "point_evaluation_instruction": (
        "First change the P1 reciprocal entrance coordinate u to the local "
        "base coordinate. Then evaluate this group-law DAG in the selected "
        "resolved chart or its finite quotient; do not require a globally "
        "reduced QQ(t) coordinate as an intermediate representation."
    ),
    "marked_chord": "m_-P=(y-y(-P))/(x-x(-P))",
    "basis": basis,
    "dimension": 44,
    "riemann_roch": "genus_one_degree_44_implies_h0=44",
    "next_required_step": (
        "Apply the exact E7 exceptional target, the smooth P2 chart, and the "
        "remaining E8/finite conditions as resolved quotient blocks."
    ),
    "boundary": (
        "This fixes the exact generic ambient and its local-evaluation interface. "
        "It does not claim that the complete vertical-condition matrix, its kernel, "
        "a new parameter, or transported child coordinates have been computed."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDAMBIENT|degree=44|basis=44|point=DAG|"
    "status=PASS_EXACT_Q6_THIRD_SYMBOLIC_GENERIC_RR_AMBIENT",
    flush=True,
)
