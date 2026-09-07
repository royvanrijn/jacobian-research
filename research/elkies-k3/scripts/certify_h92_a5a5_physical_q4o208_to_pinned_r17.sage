#!/usr/bin/env sage -python
"""Certify the physical q4/orbit208 2A5 -> 3A3 splice to pinned R17."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
MARKING = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
NEIGHBORS = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-neighbors.json"
FRONTIER = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-frontier.json"
COST = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-equation-cost.json"
GENERIC = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-3a3-lattice-certificate.json"
CURRENT = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
P1229_QQ = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


marking = json.loads(MARKING.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
frontier = json.loads(FRONTIER.read_text())
cost = json.loads(COST.read_text())
generic = json.loads(GENERIC.read_text())
current = json.loads(CURRENT.read_text())
p1229_qq = json.loads(P1229_QQ.read_text())
assert marking["status"] == "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING"
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS_TARGET_FILTERED"
assert frontier["status"] == "PASS_EXACT_MARKED_ROOT_ADAPTED_FRONTIER_RANKING"
assert cost["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
assert generic["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert current["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"
assert p1229_qq["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"

source_frame_path = ROOT / marking["frame_output"]
source_frame = load_matrix(source_frame_path)
source_gram = block_diagonal_matrix(U2, -source_frame)
current_frame_path = ROOT / current["frame_output"]
current_frame = load_matrix(current_frame_path)
current_gram = block_diagonal_matrix(U2, -current_frame)
pinned_frame = load_matrix(PINNED_FRAME)
pinned_gram = block_diagonal_matrix(U2, -pinned_frame)

candidate_id = {"q": 4, "old_fibre_degree": 2, "orbit_index": 208}
raw = next(
    item for item in neighbors["neighbors"]
    if (int(item["q"]), int(item["old_fiber_degree"]), int(item["orbit_index"])) == (4, 2, 208)
)
scored = next(item for item in cost["retained_candidates"] if item["candidate_id"] == candidate_id)
ranked = next(item for item in frontier["ranked_candidates"] if item["candidate_id"] == candidate_id)
fibre = vector(ZZ, raw["fiber"])
assert fibre == vector(ZZ, [2, 2] + [1] * 10 + [0, 0, -1, 0, 0, 0, 1])
assert fibre * source_gram * fibre == 0
assert raw["child_root_data"] == [9, 36, 64]
assert raw["child_ade"] == "A3+A3+A3" and int(raw["child_mw_rank"]) == 8
assert generic["candidate_id"]["q"] == 4
assert generic["first_edge_nef_audit"]["nef_in_selected_component_chamber"]
assert generic["first_edge_exact_horizontal_nef_gate"]
assert generic["first_edge_exact_negative_horizontal_walls"] == []
assert ranked["minimum_section_intersection"] == "0"

curves = {
    name: vector(ZZ, value)
    for name, value in marking["equation_explicit_curves_in_child"].items()
}
degrees = {name: int(curve * source_gram * fibre) for name, curve in curves.items()}
assert all(value >= 0 for value in degrees.values())
degree_one_names = sorted(name for name, value in degrees.items() if value == 1)
assert degree_one_names == [
    "first_I6_affine_component", "old_A11_component_5",
    "old_A11_component_7", "second_I6_affine_component",
]

# The special I4 member is literally four already-explicit curves.  This is
# stronger than a quotient-MW assertion and fixes the RR signs geometrically.
old_zero = curves["old_zero"]
p1229 = curves["P1229"]
c10 = curves["old_A11_component_10"]
c8 = curves["old_A11_component_8"]
special_names = ["old_zero", "P1229", "old_A11_component_10", "old_A11_component_8"]
special_curves = [old_zero, p1229, c10, c8]
assert fibre == sum(special_curves, vector(ZZ, [0] * 19))
special_intersection = matrix(ZZ, [
    [left * source_gram * right for right in special_curves]
    for left in special_curves
])
assert special_intersection == matrix(ZZ, (
    (-2, 0, 1, 1), (0, -2, 1, 1),
    (1, 1, -2, 0), (1, 1, 0, -2),
))
assert p1229 == vector(ZZ, scored["horizontal"]["section"])
assert scored["horizontal"]["P_dot_O"] == 0
assert scored["horizontal"]["vertical_layers"] == 2
assert scored["expected_RR_ambient"] == 4
assert p1229_qq["P1229"]["NS_coordinates"]
assert p1229_qq["P1229"]["degrees_X_Y_Z"] == [4, 6, 0]


def zero_frame(section):
    mate = fibre + section
    kernel = matrix(ZZ, [list(fibre * source_gram), list(mate * source_gram)]).right_kernel_matrix()
    basis = matrix(ZZ, [list(fibre), list(mate)] + [list(row) for row in kernel.rows()])
    assert abs(basis.det()) == 1
    child_frame = -(kernel * source_gram * kernel.transpose())
    child_gram = block_diagonal_matrix(U2, -child_frame)
    assert basis * source_gram * basis.transpose() == child_gram
    roots = matrix(ZZ, pari(child_frame).qfminim(2)[2]).transpose()
    root_module = roots.row_module()
    root_basis = root_module.basis_matrix()
    root_gram = root_basis * child_frame * root_basis.transpose()
    root_data = [int(root_module.rank()), int(2 * roots.nrows()), int(abs(root_gram.det()))]
    assert root_data == [9, 36, 64]

    # PARI returns Q with Q^t * current_frame * Q = child_frame.
    qiso = matrix(ZZ, pari(child_frame).qfisom(pari(current_frame)))
    assert qiso and qiso.transpose() * current_frame * qiso == child_frame
    current_tail_in_child = qiso.transpose().inverse().change_ring(ZZ)
    current_in_child = block_diagonal_matrix(identity_matrix(ZZ, 2), current_tail_in_child)
    assert current_in_child * child_gram * current_in_child.transpose() == current_gram
    current_in_source = current_in_child * basis
    return {
        "basis": basis, "inverse": basis.inverse().change_ring(ZZ),
        "child_frame": child_frame, "child_gram": child_gram,
        "root_basis": root_basis, "root_gram": root_gram, "root_data": root_data,
        "current_in_child": current_in_child,
        "child_in_current": current_in_child.inverse().change_ring(ZZ),
        "current_in_source": current_in_source,
        "metrics": {
            "child_frame_max": int(max(abs(item) for item in child_frame.list())),
            "current_tail_transport_max": int(max(abs(item) for item in current_tail_in_child.list())),
            "current_basis_in_source_max": int(max(abs(item) for item in current_in_source.list())),
        },
    }


zero_candidates = []
for name in degree_one_names:
    data = zero_frame(curves[name])
    zero_candidates.append({"name": name, "section": curves[name], **data})
zero_candidates.sort(key=lambda item: (
    item["metrics"]["current_basis_in_source_max"],
    item["metrics"]["child_frame_max"], item["name"],
))
selected = zero_candidates[0]
assert selected["name"] == "old_A11_component_5"
assert [(item["name"], item["metrics"]) for item in zero_candidates] == [
    ("old_A11_component_5", {
        "child_frame_max": 44, "current_tail_transport_max": 5,
        "current_basis_in_source_max": 5,
    }),
    ("first_I6_affine_component", {
        "child_frame_max": 40, "current_tail_transport_max": 5,
        "current_basis_in_source_max": 6,
    }),
    ("second_I6_affine_component", {
        "child_frame_max": 30, "current_tail_transport_max": 7,
        "current_basis_in_source_max": 7,
    }),
    ("old_A11_component_7", {
        "child_frame_max": 84, "current_tail_transport_max": 9,
        "current_basis_in_source_max": 9,
    }),
]

# Compose the canonical suffix and the complete pinned endpoint in both
# directions.  No ADE/MW-only identification is used.
pinned_in_current = matrix(ZZ, current["pinned_R17_basis_in_source"])
pinned_in_source = pinned_in_current * selected["current_in_source"]
source_in_pinned = pinned_in_source.inverse().change_ring(ZZ)
assert abs(pinned_in_source.det()) == 1
assert pinned_in_source * source_gram * pinned_in_source.transpose() == pinned_gram
assert source_in_pinned * pinned_gram * source_in_pinned.transpose() == source_gram
physical_in_equation = matrix(ZZ, marking["physical_component_chamber_basis_in_equation_A11"])
pinned_in_equation = pinned_in_source * physical_in_equation
equation_in_pinned = pinned_in_equation.inverse().change_ring(ZZ)
assert abs(pinned_in_equation.det()) == 1

inputs = (
    MARKING, NEIGHBORS, FRONTIER, COST, GENERIC, CURRENT,
    current_frame_path, P1229_QQ, PINNED_FRAME, source_frame_path,
)
payload = {
    "schema": "elkies-k3.h3-a5a5-physical-q4o208-to-pinned-r17.v1",
    "status": "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17",
    "candidate_id": candidate_id,
    "parent": "equation-effective component9-zero physical 2A5/MW7",
    "fibre": {
        "class_in_parent": entries(fibre), "square": 0, "old_fibre_degree": 2,
        "q": 4, "primitive": True, "nef": True,
        "all_section_minimum_intersection": ranked["minimum_section_intersection"],
        "finite_negative_horizontal_walls": [],
        "explicit_curve_degrees": degrees,
    },
    "compiler_profile": {
        "P_dot_O": 0, "expected_RR_ambient": 4, "vertical_connected_layers": 2,
        "horizontal_section": "P1229", "horizontal_section_degrees_X_Y_Z": [4, 6, 0],
        "literal_special_I4": {
            "identity": "F_q4 = old_zero + P1229 + old_A11_component_10 + old_A11_component_8",
            "components": special_names,
            "intersection_gram": rows(special_intersection),
        },
        "explicit_degree_zero_curves": scored["explicit_degree_zero_curves"],
        "explicit_degree_one_curves": scored["explicit_degree_one_curves"],
        "operational_equation_cost_terms": scored["equation_cost_terms"],
        "operational_equation_cost_score": int(scored["equation_cost_score"]),
    },
    "effective_zero_candidates": [
        {"zero": item["name"], "section": entries(item["section"]), **item["metrics"]}
        for item in zero_candidates
    ],
    "selection": {
        "zero": selected["name"],
        "reason": (
            "minimum exact canonical-current-3A3 basis growth among the four already-explicit "
            "physical degree-one curves; C5 also has an existing equation-level component slice"
        ),
        "marked_U": {
            "fibre_in_parent": entries(fibre),
            "mate_in_parent": entries(selected["basis"].row(1)),
            "zero_in_parent": entries(selected["section"]),
        },
        "child_root_data": selected["root_data"], "child_ADE": "3A3", "child_MW_rank": 8,
        "child_root_basis": rows(selected["root_basis"]),
        "child_root_gram": rows(selected["root_gram"]),
        "child_frame": rows(selected["child_frame"]),
    },
    "transport": {
        "parent_to_effective_zero_child_basis": rows(selected["basis"]),
        "effective_zero_child_to_parent_basis": rows(selected["inverse"]),
        "current_3A3_basis_in_effective_zero_child": rows(selected["current_in_child"]),
        "effective_zero_child_basis_in_current_3A3": rows(selected["child_in_current"]),
        "current_3A3_basis_in_parent": rows(selected["current_in_source"]),
        "all_determinants_absolute_one": True, "all_Gram_transports_exact": True,
    },
    "endpoint": {
        "name": "pinned_R17", "root_data": [0, 0, 1],
        "pinned_R17_basis_in_parent": rows(pinned_in_source),
        "parent_basis_in_pinned_R17": rows(source_in_pinned),
        "pinned_R17_basis_in_equation_A11": rows(pinned_in_equation),
        "equation_A11_basis_in_pinned_R17": rows(equation_in_pinned),
        "forward_determinant": int(pinned_in_source.det()),
        "inverse_determinant": int(source_in_pinned.det()),
        "Gram_transport_exact": True,
    },
    "route": {
        "q_sequence_from_A11": [8, 4, 4, 4, 4, 4, 4, 4, 6],
        "splice": "A11 --q8/orbit12--> 2A5 --physical q4/orbit208--> canonical 3A3",
        "suffix": "unchanged canonical 3A3-to-pinned-R17 certified suffix",
    },
    "proof_boundary": (
        "Exact physical-component/all-section/finite-horizontal-wall nef certificate, literal "
        "four-curve I4 divisor, exact marked U with effective C5 zero, complete 3A3 root data, "
        "bidirectional determinant-one NS transports, and full pinned-R17 Gram identification. "
        "The q4 H0 basis, quartic, Jacobian, remaining fibres, and equation marking are not yet compiled."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
assert payload["compiler_profile"]["operational_equation_cost_score"] == -1412
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5Q4O208PINNED|q=4|PO=0|RR=4|special=O+P1229+C10+C8|zero=C5|"
    "child=3A3/MW8|score=-1412|landing_det={}|endpoint_det={}|status={}|output={}".format(
        selected["current_in_source"].det(), pinned_in_source.det(), payload["status"], OUTPUT,
    ), flush=True,
)
