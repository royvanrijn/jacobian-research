#!/usr/bin/env sage -python
"""Compose the exact q25/MW7 -> degree-3 intermediate -> pinned-R17 detour."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
FIRST = GENERATED / "elkies-k3-h3-pinned-r17-q12d3-mw3-cvp-lattice-certificate.json"
SECOND = GENERATED / "elkies-k3-h3-pinned-r17-q12d3-q28d2-q25mw7-lattice-certificate.json"
SEMISTABLE = GENERATED / "elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json"
OUTPUT = GENERATED / "elkies-k3-h3-q25mw7-pinned-r17-degree2-degree3-detour.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def root_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return [0, 0, 1]
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = half + [-root for root in half]
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    return [int(basis.rank()), int(count), int(abs((basis * gram * basis.transpose()).det()))]


first = json.loads(FIRST.read_text())
second = json.loads(SECOND.read_text())
semistable = json.loads(SEMISTABLE.read_text())
assert first["status"] == "PASS_EXACT_PINNED_R17_TARGETED_CANDIDATE_CERTIFICATE"
assert second["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert semistable["status"] == "PASS_EXACT_SEMISTABLE_MW2_TO_PINNED_R17_NEF_SUFFIX"

pinned_frame = load_matrix(ROOT / "elkies-k3/data/lattice/rank17_gram.txt")
intermediate_frame = load_matrix(ROOT / first["frame_output"])
q25_candidate_frame = load_matrix(ROOT / second["frame_output"])
g_pinned = block_diagonal_matrix(U2, -pinned_frame)
g_intermediate = block_diagonal_matrix(U2, -intermediate_frame)
g_candidate = block_diagonal_matrix(U2, -q25_candidate_frame)

# Rows are always the new basis written in the preceding basis.
intermediate_in_pinned = matrix(ZZ, first["source_to_child_basis"])
candidate_in_intermediate = matrix(ZZ, second["source_to_child_basis"])
candidate_in_pinned = candidate_in_intermediate * intermediate_in_pinned
assert intermediate_in_pinned * g_pinned * intermediate_in_pinned.transpose() == g_intermediate
assert candidate_in_intermediate * g_intermediate * candidate_in_intermediate.transpose() == g_candidate
assert candidate_in_pinned * g_pinned * candidate_in_pinned.transpose() == g_candidate
assert all(abs(value.det()) == 1 for value in (
    intermediate_in_pinned, candidate_in_intermediate, candidate_in_pinned
))

q25_fibre_in_intermediate = vector(ZZ, first["target_fibres_in_child"]["q25_mw7"])
assert q25_fibre_in_intermediate == candidate_in_intermediate.row(0)
q25_fibre_in_candidate = q25_fibre_in_intermediate * candidate_in_intermediate.inverse()
assert q25_fibre_in_candidate == vector(ZZ, [1, 0] + [0] * 17)
assert second["target_profiles"]["q25_mw7"]["old_fibre_degree"] == 0
assert second["target_profiles"]["q25_mw7"]["nef_audit"]["same_fibre_ray"]

# This is the physical q25 basis already chamber-certified in the canonical
# reverse suffix, transported first to the intermediate and then to the new
# candidate frame.
q25_in_intermediate = matrix(
    ZZ, first["target_transports"]["q25_mw7"]["target_basis_in_child"]
)
q25_in_candidate = q25_in_intermediate * candidate_in_intermediate.inverse()
q25_in_pinned = q25_in_intermediate * intermediate_in_pinned
q25_gram = q25_in_candidate * g_candidate * q25_in_candidate.transpose()
canonical_q25_in_pinned = matrix(ZZ, semistable["steps"][3]["inverse_transport"])
assert q25_in_pinned.row(0) == canonical_q25_in_pinned.row(0)
canonical_q25_gram = canonical_q25_in_pinned * g_pinned * canonical_q25_in_pinned.transpose()
q25_physical_change = q25_in_pinned * canonical_q25_in_pinned.inverse()
assert q25_physical_change in MatrixSpace(ZZ, 19) and abs(q25_physical_change.det()) == 1
assert q25_physical_change * canonical_q25_gram * q25_physical_change.transpose() == q25_gram
assert q25_in_candidate.row(0) == vector(ZZ, [1, 0] + [0] * 17)
q25_frame = -q25_gram[2:, 2:]
assert q25_gram[:2, :2] == U2 and q25_gram[:2, 2:] == 0
assert root_data(q25_frame) == [10, 26, 512]
assert root_data(q25_candidate_frame) == [10, 26, 512]

q25_reverse = first["target_transports"]["q25_mw7"]["reverse_edge_profile"]
pinned_reverse = first["target_profiles"]["pinned_R17"]
assert q25_reverse["old_fibre_degree"] == 2
assert q25_reverse["P_dot_O"] == 1
assert q25_reverse["nef_audit"]["nef_in_selected_component_chamber"]
assert pinned_reverse["old_fibre_degree"] == 3
assert pinned_reverse["P_dot_O"] == 10
assert pinned_reverse["nef_audit"]["nef_in_selected_component_chamber"]
assert first["first_edge_nef_audit"]["nef"]
assert second["first_edge_nef_audit"]["nef_in_selected_component_chamber"]

inputs = (FIRST, SECOND, SEMISTABLE)
payload = {
    "schema": "elkies-k3.h3-q25mw7-pinned-r17-degree2-degree3-detour.v1",
    "status": "PASS_EXACT_Q25_MW7_TO_PINNED_R17_DEGREE2_DEGREE3_DETOUR",
    "route": [
        {
            "source": "q25_mw7",
            "target": "q12d3_intermediate_6roots_mw11",
            "q": q25_reverse["q"],
            "old_fibre_degree": q25_reverse["old_fibre_degree"],
            "P_dot_O": q25_reverse["P_dot_O"],
            "source_root_data": [10, 26, 512],
            "target_root_data": first["child"]["root_data"],
            "nef_audit": q25_reverse["nef_audit"],
            "target_basis_in_source": rows(q25_in_intermediate.inverse().change_ring(ZZ)),
            "source_basis_in_target": rows(q25_in_intermediate),
        },
        {
            "source": "q12d3_intermediate_6roots_mw11",
            "target": "pinned_R17",
            "q": pinned_reverse["q"],
            "old_fibre_degree": pinned_reverse["old_fibre_degree"],
            "P_dot_O": pinned_reverse["P_dot_O"],
            "source_root_data": first["child"]["root_data"],
            "target_root_data": [0, 0, 1],
            "nef_audit": pinned_reverse["nef_audit"],
            "target_basis_in_source": rows(intermediate_in_pinned.inverse().change_ring(ZZ)),
            "source_basis_in_target": rows(intermediate_in_pinned),
        },
    ],
    "endpoint_identification": {
        "q25_fibre_in_candidate_frame": list(map(int, q25_fibre_in_candidate)),
        "q25_physical_basis_in_candidate": rows(q25_in_candidate),
        "candidate_basis_in_q25_physical": rows(q25_in_candidate.inverse().change_ring(ZZ)),
        "q25_physical_basis_in_pinned_R17": rows(q25_in_pinned),
        "same_fibre_as_canonical_q25_basis": True,
        "canonical_q25_to_selected_physical_basis": rows(q25_physical_change),
        "selected_physical_to_canonical_q25_basis": rows(
            q25_physical_change.inverse().change_ring(ZZ)
        ),
        "full_unimodular_canonical_q25_identification": True,
        "q25_root_data": [10, 26, 512],
        "pinned_R17_root_data": [0, 0, 1],
        "candidate_basis_in_pinned_R17": rows(candidate_in_pinned),
        "pinned_R17_basis_in_candidate": rows(candidate_in_pinned.inverse().change_ring(ZZ)),
    },
    "comparison": {
        "known_direct_q25_to_pinned": {"q": 40, "old_fibre_degree": 5, "P_dot_O": 3},
        "new_detour_edge_profiles": [
            {"q": q25_reverse["q"], "old_fibre_degree": 2, "P_dot_O": 1},
            {"q": pinned_reverse["q"], "old_fibre_degree": 3, "P_dot_O": 10},
        ],
        "promotion_decision": "NOT_PROMOTED_NO_CERTIFIED_A11_TO_Q25_PREFIX_AND_PO10_TRADEOFF",
    },
    "proof_boundary": (
        "Exact two-edge lattice route with primitive nef fibres, component and "
        "all-horizontal-curve gates, complete determinant-one transports both ways, "
        "and exact equality with the canonical physical q25 basis and pinned R17. "
        "It is not a full replacement route from equation A11."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q25D3DETOUR|edges=q{}d{}po{},q{}d{}po{}|q25_exact=1|pinned_exact=1|"
    "det={}|status={}|output={}".format(
        q25_reverse["q"], q25_reverse["old_fibre_degree"], q25_reverse["P_dot_O"],
        pinned_reverse["q"], pinned_reverse["old_fibre_degree"], pinned_reverse["P_dot_O"],
        candidate_in_pinned.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
