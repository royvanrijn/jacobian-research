#!/usr/bin/env sage -python
"""Derive the resolved E7 target for the source-nef H3 q=8 pencil.

The q=8 class is first stripped of its old E7+E8 fixed components by
``derive_h92_q8_generic_rr_ambient.sage``.  At the H92 III* point, this
script identifies the source E7 component numbering with the explicit
blow-up numbering used by ``derive_h92_q6_e7_resolution.sage`` and computes
the precise local line-bundle class.

In the resolved numbering, the q=6 module has exceptional cycle c6 and the
q=8 target has c8.  The key exact relation is

    c8 - 9*c6 = (2,5,6,4,6,3,5),

an integral exceptional divisor.  Thus the future q=8 quotient compiler must
start from the ninth tensor power of the non-Cartier q=6 marked module and
apply this displayed integral twist.  This is local target data only: it does
not construct its finite quotient matrix or a q=8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
Q8_AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-local-target.json"


def load_frame(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=Q8_AMBIENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
q8 = vector(ZZ, ambient["source_q8_lattice_class"])
assert len(q8) == 19
support = ambient["generic_fibre_support_certificate"]
assert support["old_fiber_degree"] == 18
assert support["fiber_twist"] == -11
assert [(entry["name"], entry["multiplicity"]) for entry in support["horizontal_support"]] == [
    ("O", 9), ("-P1", 9),
]
assert [entry["coefficient"] for entry in support["vertical_support"][:7]] == [
    2, 3, 4, 6, 5, 5, 6,
]

source_ns = block_diagonal_matrix(
    matrix(ZZ, ((0, 1), (1, 0))), -load_frame(FRAME)
)
source_simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
source_cartan = -matrix(ZZ, [
    [left * source_ns * right for right in source_simple[:7]]
    for left in source_simple[:7]
])
assert source_cartan == matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0],
    [0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, -1, 2],
])

# This is the Cartan order arising from the actual chart resolution.  The
# unique marked graph isomorphism preserves the affine attachment: source
# node 1 and the resolved divisor of Z are both its attachment node.
resolved_cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
source_to_resolved = (1, 6, 4, 3, 7, 2, 5)
permutation = matrix(ZZ, 7, 7, lambda source, resolved:
                     ZZ(source_to_resolved[source] - 1 == resolved))
assert permutation.transpose() * source_cartan * permutation == resolved_cartan

# The multiplicity vector of the old affine E7 component in the source
# numbering maps to div(Z) in the resolved germ.  This pins the graph
# isomorphism geometrically rather than merely up to a Dynkin automorphism.
source_highest = vector(ZZ, (2, 2, 3, 4, 3, 2, 1))
resolved_vz = vector(ZZ, (2, 2, 4, 3, 1, 2, 3))
assert source_cartan * source_highest == vector(ZZ, (1, 0, 0, 0, 0, 0, 0))
assert resolved_cartan * resolved_vz == vector(ZZ, (1, 0, 0, 0, 0, 0, 0))
assert permutation.transpose() * source_highest == resolved_vz

q8_source_degrees = vector(ZZ, [
    q8 * source_ns * component for component in source_simple[:7]
])
assert q8_source_degrees == vector(ZZ, (0, 0, 0, 0, 1, 1, 2))
q8_resolved_degrees = permutation.transpose() * q8_source_degrees
assert q8_resolved_degrees == vector(ZZ, (0, 1, 0, 0, 2, 0, 1))

# The certified q=6 source divisor has one degree at source E7_7.  Its
# non-Cartier marked module is Z*J_-P1^dual; this establishes the q=6 base
# class for the tensor comparison below.
q6_source_degrees = vector(ZZ, (0, 0, 0, 0, 0, 0, 1))
q6_resolved_degrees = permutation.transpose() * q6_source_degrees
assert q6_resolved_degrees == vector(ZZ, (0, 0, 0, 0, 1, 0, 0))

# For an exceptional coefficient vector c, O(sum c_i E_i) has degree
# -Cartan*c on the exceptional curves.  The rational cycles here encode the
# unavoidable non-Cartier part; their difference is integral.
q6_cycle = -resolved_cartan.inverse() * q6_resolved_degrees
q8_cycle = -resolved_cartan.inverse() * q8_resolved_degrees
integral_twist = q8_cycle - 9 * q6_cycle
assert q6_cycle == vector(QQ, (-1, -2, -3, -2, QQ(-3) / 2, QQ(-3) / 2, QQ(-5) / 2))
assert q8_cycle == vector(QQ, (-7, -13, -21, -14, QQ(-15) / 2, QQ(-21) / 2, QQ(-35) / 2))
assert integral_twist == vector(QQ, (2, 5, 6, 4, 6, 3, 5))
assert -resolved_cartan * q6_cycle == q6_resolved_degrees
assert -resolved_cartan * q8_cycle == q8_resolved_degrees

# This correction is an integral *line-bundle* twist on the resolution, but
# it is not an anti-nef exceptional cycle.  Consequently it cannot be
# replaced by the complete ideal of one exceptional cycle in the singular
# E7 local ring.  This is an exact obstruction to a tempting one-chart
# quotient implementation: the compiler must retain the resolved-chart
# trivializations and their gluing.
twist_intersections = -resolved_cartan * integral_twist
assert twist_intersections == vector(ZZ, (0, 1, 0, 0, -7, 0, 1))
assert any(value > 0 for value in twist_intersections)
assert any(value < 0 for value in twist_intersections)

# The cycle relation is a statement about the exceptional class only.  It
# does not lower the horizontal marked-branch pole order: at the generic
# point of the -P1 branch, J_-P1^dual has a simple pole, hence its ninth
# tensor power has pole order nine.  Retaining this datum prevents an invalid
# replacement of (J_-P1^dual)^9 by one class-group-equivalent branch factor.
branch_pole_order = 9
assert branch_pole_order * 1 == 9

payload = {
    "schema": "elkies-k3.h92-q8-e7-local-target.v1",
    "status": "PASS_EXACT_Q8_E7_LOCAL_TARGET",
    "inputs": {
        "q8_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "source_frame": str(FRAME.relative_to(ROOT)),
        "q6_marked_module": "Z*J_-P1^dual",
        "q6_marked_module_source": "elkies-k3/scripts/derive_h92_q6_e7_p1_branch_module.sage",
        "resolved_chart_source": "elkies-k3/scripts/derive_h92_q6_e7_resolution.sage",
    },
    "source_to_resolved_component_map": {
        "source_E7_i_to_resolved_E7_j": list(source_to_resolved),
        "geometric_normalization": "source affine multiplicity vector maps to div(Z)",
    },
    "q8": {
        "source_component_degrees": list(map(int, q8_source_degrees)),
        "resolved_component_degrees": list(map(int, q8_resolved_degrees)),
        "exceptional_cycle": [str(value) for value in q8_cycle],
        "support_certificate_E7_vertical_coefficients": [
            entry["coefficient"] for entry in support["vertical_support"][:7]
        ],
    },
    "q6": {
        "source_component_degrees": list(map(int, q6_source_degrees)),
        "resolved_component_degrees": list(map(int, q6_resolved_degrees)),
        "exceptional_cycle": [str(value) for value in q6_cycle],
    },
    "tensor_comparison": {
        "identity": "c8=9*c6+(2,5,6,4,6,3,5)",
        "integral_exceptional_twist": [str(value) for value in integral_twist],
        "intersection_degrees": list(map(int, twist_intersections)),
        "non_antinef_obstruction": (
            "The integral correction has exceptional intersection degrees "
            "(0,1,0,0,-7,0,1), so it is not anti-nef and cannot be "
            "represented by one complete ideal downstairs in the singular "
            "E7 local ring. A q8 compiler must use resolved-chart "
            "trivializations and gluing."
        ),
        "compiler_instruction": (
            "Represent the q8 E7 module as the ninth tensor power of "
            "J_-P1^dual, with marked-branch pole order nine, followed by "
            "the displayed integral exceptional twist."
        ),
        "marked_branch_pole_order": branch_pole_order,
        "forbidden_simplification": "Do not replace (J_-P1^dual)^9 by J_-P1^dual; their branch pole orders are 9 and 1.",
    },
    "boundary": (
        "This certifies the resolved E7 target only. It does not derive the "
        "finite q8 quotient map, the E8 target, a global kernel, or a pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8E7TARGET|degrees=0,1,0,0,2,0,1|twist=2,5,6,4,6,3,5|"
    "status=PASS_EXACT_Q8_E7_LOCAL_TARGET",
    flush=True,
)
