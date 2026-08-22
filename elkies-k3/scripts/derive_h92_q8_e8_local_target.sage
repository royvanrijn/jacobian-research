#!/usr/bin/env sage -python
"""Derive the source E8 target for the source-nef H3 q=8 pencil.

At the H92 II* fibre the q=6 marked module is ``u*<1,Q>`` and has degree zero
on every finite E8 component.  This script computes the q=8 restriction
degrees directly from the source-nef class.  The unique exceptional cycle
with those degrees is integral because E8 is unimodular:

    (-4,-5,-7,-10,-8,-6,-4,-2).

Consequently the q=8 E8 module is the ninth q=6 module tensor power followed
by this exact exceptional twist.  The script also identifies every source E8
component with the explicit ordinary-blow-up components.  Deriving the finite
quotient map remains a later gate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
Q8_AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-local-target.json"


def load_frame(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=Q8_AMBIENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
q8 = vector(ZZ, ambient["source_q8_lattice_class"])
support = ambient["generic_fibre_support_certificate"]
assert support["old_fiber_degree"] == 18
assert support["fiber_twist"] == -11
assert [(entry["name"], entry["multiplicity"]) for entry in support["horizontal_support"]] == [
    ("O", 9), ("-P1", 9),
]
assert [entry["coefficient"] for entry in support["vertical_support"][7:]] == [
    -4, -5, -7, -10, -8, -6, -4, -2,
]

source_ns = block_diagonal_matrix(
    matrix(ZZ, ((0, 1), (1, 0))), -load_frame(FRAME)
)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
e8_components = simple[7:15]
e8_cartan = -matrix(ZZ, [
    [left * source_ns * right for right in e8_components]
    for left in e8_components
])
assert e8_cartan == matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])
assert e8_cartan.det() == 1

q8_degrees = vector(ZZ, [q8 * source_ns * component for component in e8_components])
assert q8_degrees == vector(ZZ, (1, 0, 0, 0, 0, 0, 0, 0))
q6_degrees = vector(ZZ, (0,) * 8)
q8_cycle = -e8_cartan.inverse() * q8_degrees
assert q8_cycle == vector(QQ, (-4, -5, -7, -10, -8, -6, -4, -2))
assert -e8_cartan * q8_cycle == q8_degrees

# The E8 coordinates of the source-nef divisor are exactly this cycle: the
# old horizontal divisors have no E8 finite-component coefficient.
assert vector(QQ, q8[9:17]) == q8_cycle

# The E8 blow-up sequence has components B1,B2,B3,B4, the ordinary node N3
# in the third x-chart, and the three terminal-node components N40,N4B,N4inf.
# Its edges follow directly from the displayed centres: B3 separates B1/B2,
# N3 separates B1/B3, B4 separates B2/B3, and N40/N4inf separate its B2/B3
# intersections while N4B is its remaining smooth-point branch.
chart_component_order = ("B1", "B2", "B3", "B4", "N3", "N40", "N4B", "N4inf")
chart_edges = (
    ("B1", "N3"), ("N3", "B3"), ("B3", "N4inf"), ("N4inf", "B4"),
    ("B4", "N40"), ("N40", "B2"), ("B4", "N4B"),
)
chart_cartan = matrix(ZZ, 8, 8, lambda left, right:
                      ZZ(2 if left == right else 0))
for left, right in chart_edges:
    left_index = chart_component_order.index(left)
    right_index = chart_component_order.index(right)
    chart_cartan[left_index, right_index] = chart_cartan[right_index, left_index] = -1
assert chart_cartan.det() == 1

# There is a unique E8 graph identification compatible with the old affine
# fibre: source node eight meets the affine component and maps to B1, the
# first exceptional component in the actual infinity chart.
source_to_chart = (2, 7, 6, 4, 8, 3, 5, 1)
permutation = matrix(ZZ, 8, 8, lambda source, chart:
                     ZZ(source_to_chart[source] - 1 == chart))
assert permutation.transpose() * e8_cartan * permutation == chart_cartan
source_highest = vector(ZZ, (2, 3, 4, 6, 5, 4, 3, 2))
chart_u_cycle = permutation.transpose() * source_highest
assert chart_u_cycle == vector(ZZ, (2, 2, 4, 6, 3, 4, 3, 5))
assert e8_cartan * source_highest == vector(ZZ, (0, 0, 0, 0, 0, 0, 0, 1))
assert chart_cartan * chart_u_cycle == vector(ZZ, (1, 0, 0, 0, 0, 0, 0, 0))
q8_chart_degrees = permutation.transpose() * q8_degrees
q8_chart_cycle = permutation.transpose() * q8_cycle
assert q8_chart_degrees == vector(ZZ, (0, 1, 0, 0, 0, 0, 0, 0))
assert q8_chart_cycle == vector(QQ, (-2, -4, -6, -10, -4, -7, -5, -8))
assert -chart_cartan * q8_chart_cycle == q8_chart_degrees

payload = {
    "schema": "elkies-k3.h92-q8-e8-local-target.v1",
    "status": "PASS_EXACT_Q8_E8_SOURCE_TARGET",
    "inputs": {
        "q8_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "source_frame": str(FRAME.relative_to(ROOT)),
        "q6_marked_module": "u*<1,Q>",
        "q6_marked_module_source": "elkies-k3/scripts/derive_h92_q6_e8_p1_branch_module.sage",
        "resolved_chart_source": "elkies-k3/scripts/derive_h92_q6_e8_resolution.sage",
    },
    "source_component_order": "pinned E8 simple-root order in kumar_e7e8_mw2_frame_3.txt",
    "chart_component_map": {
        "chart_component_order": list(chart_component_order),
        "chart_edges": [list(edge) for edge in chart_edges],
        "source_E8_i_to_chart_component": list(source_to_chart),
        "geometric_normalization": "source affine attachment node eight maps to B1",
        "chart_u_exceptional_cycle": [str(value) for value in chart_u_cycle],
    },
    "q8": {
        "component_degrees": list(map(int, q8_degrees)),
        "exceptional_cycle": [str(value) for value in q8_cycle],
        "source_nef_E8_coordinates": [int(value) for value in q8[9:17]],
        "chart_component_degrees": list(map(int, q8_chart_degrees)),
        "chart_exceptional_cycle": [str(value) for value in q8_chart_cycle],
        "support_certificate_E8_vertical_coefficients": [
            entry["coefficient"] for entry in support["vertical_support"][7:]
        ],
    },
    "q6": {"component_degrees": list(map(int, q6_degrees))},
    "tensor_comparison": {
        "identity": "c8=9*0+(-4,-5,-7,-10,-8,-6,-4,-2)",
        "integral_exceptional_twist": [str(value) for value in q8_cycle],
        "compiler_instruction": (
            "Represent the q8 E8 module as the ninth tensor power of "
            "u*<1,Q>, followed by the displayed source-E8 exceptional twist."
        ),
    },
    "boundary": (
        "This maps the target to the actual E8 blow-up components, but does "
        "not derive the finite q8 quotient map or global kernel."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8E8TARGET|degrees=1,0,0,0,0,0,0,0|"
    "twist=-4,-5,-7,-10,-8,-6,-4,-2|status=PASS_EXACT_Q8_E8_SOURCE_TARGET",
    flush=True,
)
