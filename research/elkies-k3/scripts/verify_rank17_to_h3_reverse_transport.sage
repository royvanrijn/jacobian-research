#!/usr/bin/env sage -python
"""Certify and export the lossless reverse transport from pinned R17 to H3.

The existing forward corridor checker proves a sequence of primitive
U-neighbours but stops at an unnamed rootless determinant-948 frame.  This
checker performs the missing endpoint comparison with the recovered pinned
``rank17_gram.txt``, prepends the exact H3 q6 and q8 transports, and inverts
the complete 19-dimensional Neron--Severi transport.

The generated ledger deliberately retains every stage Gram, every incoming
transition and its inverse, and both the H3 and pinned-R17 coordinates of
every stage basis.  It is a lattice/marking certificate, not an execution of
the eleven characteristic-zero pencils after D13.
"""

from sage.all import *

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"
load(str(ENGINE))

U = matrix(ZZ, ((0, 1), (1, 0)))
H3_FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
D13_FRAME = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
PINNED_R17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
ENDPOINT_ISOMETRY = (
    ROOT
    / "elkies-k3/data/fibrations/h3_rootless_mw17_to_pinned_rank17_isometry.txt"
)
H3_Q6_Q8 = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
FORWARD_SUFFIX = (
    ROOT / "artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rank17-to-h3-reverse-transport.json"
)

EXPECTED_HASHES = {
    H3_Q6_Q8: "289b4d7e3e2556e4680f8bf523cbf75e28e7034225a5535a39470467261c8ff6",
    FORWARD_SUFFIX: "f6eac2339c86de84b79a0ddfec3229df9b9c1617110bdd9c474443e7e39fd484",
    H3_FRAME: "ba09ec834a7229e11e4ca687d187f663b6368c3e2fac9b5133bb1570e7031599",
    PINNED_R17: "ea1330ce903c55624de1705ebedf36f6f5964cdd6ce911927001db82626ff237",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display(path):
    return str(path.resolve().relative_to(ROOT))


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def ns(frame):
    return block_diagonal_matrix(U, -frame)


def frame_from_transport(transport, parent):
    child_ns = transport * ns(parent) * transport.transpose()
    assert child_ns[:2, :2] == U
    assert not any(child_ns[:2, 2:].list())
    assert not any(child_ns[2:, :2].list())
    return -child_ns[2:, 2:]


def stage_record(
    stage_id,
    ade,
    mw_rank,
    frame,
    cumulative,
    pinned_to_h3,
    incoming=None,
):
    in_h3 = cumulative
    h3_in_stage = cumulative.inverse()
    in_pinned = cumulative * pinned_to_h3
    pinned_in_stage = in_pinned.inverse()
    for value in (h3_in_stage, in_pinned, pinned_in_stage):
        assert value.change_ring(ZZ) == value
    in_h3 = in_h3.change_ring(ZZ)
    h3_in_stage = h3_in_stage.change_ring(ZZ)
    in_pinned = in_pinned.change_ring(ZZ)
    pinned_in_stage = pinned_in_stage.change_ring(ZZ)
    assert in_h3 * ns(load_matrix(H3_FRAME)) * in_h3.transpose() == ns(frame)
    assert in_pinned * ns(load_matrix(PINNED_R17)) * in_pinned.transpose() == ns(frame)
    record = {
        "stage_id": stage_id,
        "ade": ade,
        "mw_rank": mw_rank,
        "positive_frame": rows(frame),
        "stage_basis_in_h3_ns": rows(in_h3),
        "h3_basis_in_stage_ns": rows(h3_in_stage),
        "stage_basis_in_pinned_rank17_ns": rows(in_pinned),
        "pinned_rank17_basis_in_stage_ns": rows(pinned_in_stage),
    }
    if incoming is not None:
        record["incoming_neighbor"] = incoming
    return record


for path, expected in EXPECTED_HASHES.items():
    assert digest(path) == expected, (path, digest(path), expected)

h3 = load_matrix(H3_FRAME)
d13 = load_matrix(D13_FRAME)
pinned = load_matrix(PINNED_R17)
endpoint_isometry = load_matrix(ENDPOINT_ISOMETRY)
assert h3.dimensions() == d13.dimensions() == pinned.dimensions() == (17, 17)
assert h3.det() == d13.det() == pinned.det() == 948

entrance = json.loads(H3_Q6_Q8.read_text())
assert entrance["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"

# H3 -> raw q6, followed by the stored root/MW and deterministic E6+E8
# simple-root adaptations used by the pinned lattice-corridor q8 search.
t6_raw = matrix(ZZ, entrance["q6"]["neighbor_basis_in_source_ns"])
root_mw = matrix(ZZ, entrance["q6"]["root_mw_basis_in_child"])
simple_root_change = matrix(
    ZZ, entrance["q8"]["simple_root_change_in_root_block"]
)
simple_change = block_diagonal_matrix(
    simple_root_change, identity_matrix(ZZ, 3)
)
a6 = block_diagonal_matrix(
    identity_matrix(ZZ, 2), simple_change * root_mw
)
t6 = a6 * t6_raw
q6 = matrix(ZZ, entrance["q8"]["simple_frame_gram"])
assert frame_from_transport(t6, h3) == q6

# The pinned D13 frame used by the eleven-step lattice corridor is the first
# stored dominant q8 hit.  This is intentionally not silently replaced by
# the separate component-nef equation representative.
q8 = entrance["q8"]["d13_mw4_hits"][0]
t8_raw = matrix(ZZ, q8["neighbor_basis_in_q6_ns"])
a8 = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, q8["d13_root_adapted_basis_in_child"])
)
t8 = a8 * t8_raw
assert t8.change_ring(ZZ) == t8
t8 = t8.change_ring(ZZ)
t_d13 = t8 * t6
assert frame_from_transport(t8, q6) == d13
assert frame_from_transport(t_d13, h3) == d13

# Keep the exact component-nef equation marking rather than collapsing it to
# the pinned dominant lattice marking.  Both are U-embeddings in NS(H3), and
# their integral bridge generally changes the standard fibre/zero U.
q8_nef = entrance["q8"]["nef_representative"]
t8_nef_raw = matrix(ZZ, q8_nef["neighbor_basis_in_q6_ns"])
a8_nef = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, q8_nef["d13_root_adapted_basis_in_child"]),
)
t_nef_d13 = a8_nef * t8_nef_raw * t6_raw
nef_d13 = matrix(ZZ, q8_nef["d13_root_adapted_gram"])
assert frame_from_transport(t_nef_d13, h3) == nef_d13
dominant_to_nef = t_nef_d13 * t_d13.inverse()
assert dominant_to_nef.change_ring(ZZ) == dominant_to_nef
dominant_to_nef = dominant_to_nef.change_ring(ZZ)
assert dominant_to_nef * ns(d13) * dominant_to_nef.transpose() == ns(nef_d13)
standard_u_rows = identity_matrix(ZZ, 19)[:2, :]
preserves_standard_u = dominant_to_nef[:2, :] == standard_u_rows

stage_data = [
    {
        "stage_id": "H3-SOURCE",
        "ade": "E7+E8",
        "mw_rank": 2,
        "frame": h3,
        "cumulative": identity_matrix(ZZ, 19),
        "incoming": None,
    },
    {
        "stage_id": "H3-01-E8E6",
        "ade": "E8+E6",
        "mw_rank": 3,
        "frame": q6,
        "cumulative": t6,
        "incoming": {
            "q": 6,
            "factor_order": [3, 2],
            "old_fiber_degree": 2,
            "witness_in_parent_frame": entrance["q6"]["witness"],
            "transition_basis_in_parent_ns": rows(t6),
            "parent_basis_in_child_ns": rows(t6.inverse().change_ring(ZZ)),
        },
    },
    {
        "stage_id": "H3-02-D13",
        "ade": "D13",
        "mw_rank": 4,
        "frame": d13,
        "cumulative": t_d13,
        "incoming": {
            "q": 8,
            "factor_order": [4, 2],
            "old_fiber_degree": 2,
            "representative": "first stored dominant D13 hit (pinned lattice frame)",
            "witness_in_simple_q6_frame": q8["witness_simple_frame"],
            "witness_in_raw_q6_frame": q8["witness_q6_child"],
            "transition_basis_in_root_adapted_parent_ns": rows(t8),
            "parent_basis_in_child_ns": rows(t8.inverse().change_ring(ZZ)),
        },
    },
]

suffix = json.loads(FORWARD_SUFFIX.read_text())
assert suffix["status"] == "PASS_H3_D13_TO_MW17_LATTICE_PATH"
current = d13
cumulative = t_d13
suffix_composite = identity_matrix(ZZ, 19)

stage_ids = (
    "H3-03-D12",
    "H3-04-A11",
    "H3-05-2A5",
    "H3-06-3A3",
    "H3-07-A3-2A2",
    "H3-08-5A1",
    "H3-09-4A1",
    "H3-10-3A1",
    "H3-11-2A1",
    "H3-12-A1",
    "H3-13-ROOTLESS",
)

for stage_id, step in zip(stage_ids, suffix["steps"]):
    artifact = ROOT / step["artifact"]
    assert digest(artifact) == step["artifact_sha256"]
    search = json.loads(artifact.read_text())
    assert load_matrix(ROOT / search["frame"]) == current
    matches = [
        record
        for record in search["neighbors"]
        if record["orbit_index"] == step["orbit_index"]
    ]
    assert len(matches) == 1
    selected = matches[0]
    assert selected["q"] == step["q"]
    assert selected["old_fiber_degree"] == 2
    result = degree_two_neighbor(
        ns(current), matrix(ZZ, [selected["fiber"]]).row(0),
        matrix(ZZ, [[1, 0] + [0] * 17]).row(0), curves=()
    )
    raw_child = result["child_frame"]
    adaptation = matrix(ZZ, selected["child_root_adapted_basis"])
    transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * result["transport"]
    child = adaptation * raw_child * adaptation.transpose()
    assert child == matrix(ZZ, selected["child_root_adapted_frame"])
    assert transition.det() in (-1, 1)
    assert frame_from_transport(transition, current) == child
    cumulative = transition * cumulative
    suffix_composite = transition * suffix_composite
    current = child
    stage_data.append(
        {
            "stage_id": stage_id,
            "ade": step["ade"],
            "mw_rank": step["mw_rank"],
            "frame": child,
            "cumulative": cumulative,
            "incoming": {
                "q": step["q"],
                "factor_order": step["factor_order"],
                "old_fiber_degree": 2,
                "orbit_index": step["orbit_index"],
                "witness_in_parent_frame": step["witness"],
                "fiber_in_parent_ns": selected["fiber"],
                "transition_basis_in_parent_ns": rows(transition),
                "parent_basis_in_child_ns": rows(
                    transition.inverse().change_ring(ZZ)
                ),
                "source_artifact": step["artifact"],
                "source_artifact_sha256": step["artifact_sha256"],
            },
        }
    )

assert rows(suffix_composite) == suffix["composite_transport"]
assert current == matrix(ZZ, suffix["final_frame"])
endpoint = current

# PARI qfisom(endpoint, pinned) produced C with C^t P C = endpoint.  Pin C
# as data and verify it directly, so replay does not depend on a fresh
# isometry search or on a particular PARI enumeration order.
c = endpoint_isometry
assert c.det() == 1
assert c.transpose() * pinned * c == endpoint
endpoint_basis_in_pinned = block_diagonal_matrix(
    identity_matrix(ZZ, 2), c.transpose()
)
assert (
    endpoint_basis_in_pinned
    * ns(pinned)
    * endpoint_basis_in_pinned.transpose()
    == ns(endpoint)
)

total = cumulative
pinned_to_h3 = total.inverse() * endpoint_basis_in_pinned
assert pinned_to_h3.change_ring(ZZ) == pinned_to_h3
pinned_to_h3 = pinned_to_h3.change_ring(ZZ)
assert pinned_to_h3.det() in (-1, 1)
assert pinned_to_h3 * ns(pinned) * pinned_to_h3.transpose() == ns(h3)

stages = [
    stage_record(
        item["stage_id"],
        item["ade"],
        item["mw_rank"],
        item["frame"],
        item["cumulative"],
        pinned_to_h3,
        item["incoming"],
    )
    for item in stage_data
]

payload = {
    "schema": "elkies-k3.rank17-to-h3-reverse-transport.v1",
    "status": "PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT",
    "proof_boundary": (
        "This exactly identifies the selected H3 corridor endpoint with the "
        "pinned recovered rank17 positive frame and gives a lossless integral "
        "NS transport in both directions. It does not execute the eleven "
        "characteristic-zero equation pencils after D13 or transport explicit "
        "Weierstrass section functions through them."
    ),
    "source_files": {
        display(path): "sha256:" + digest(path)
        for path in (*EXPECTED_HASHES.keys(), ENDPOINT_ISOMETRY)
    },
    "endpoint_identification": {
        "relation": "C^t * pinned_rank17_gram * C = h3_corridor_endpoint_frame",
        "positive_frame_isometry_C": rows(c),
        "determinant": int(c.det()),
        "endpoint_basis_in_pinned_rank17_ns": rows(endpoint_basis_in_pinned),
    },
    "q8_marking_distinction": {
        "pinned_lattice_corridor_representative": "first stored dominant D13 hit",
        "component_nef_equation_representative": "q8.nef_representative",
        "positive_frames_equal": bool(d13 == nef_d13),
        "component_nef_positive_frame": rows(nef_d13),
        "dominant_to_component_nef_ns_basis": rows(dominant_to_nef),
        "component_nef_to_dominant_ns_basis": rows(
            dominant_to_nef.inverse().change_ring(ZZ)
        ),
        "preserves_standard_fiber_and_u_mate": bool(preserves_standard_u),
        "interpretation": (
            "The bridge is an exact NS isometry but changes the embedded U; "
            "the two D13 records must not be substituted as the same marked "
            "elliptic fibration without transporting fibre, zero, components, "
            "and sections through this matrix."
        ),
    },
    "complete_reverse_transport": {
        "relation": "R * NS(pinned_rank17) * R^t = NS(H3)",
        "pinned_rank17_to_h3_basis_R": rows(pinned_to_h3),
        "h3_to_pinned_rank17_basis": rows(
            pinned_to_h3.inverse().change_ring(ZZ)
        ),
        "determinant": int(pinned_to_h3.det()),
    },
    "stages": stages,
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "R17H3|endpoint_positive_isometry=PASS|det={}|stages={}|".format(
        c.det(), len(stages)
    )
    + "reverse_ns_transport=PASS|artifact={}|status={}".format(
        display(OUTPUT), payload["status"]
    ),
    flush=True,
)
