#!/usr/bin/env sage -python
"""Price direct Gate-B moves to the alternate A1 and rootless frames.

All fibrations are transported into the pinned R17 Neron--Severi lattice.
The reported intersection is therefore the literal old-fibre degree of the
selected neighbour-frame graph ray, not an ADE-label comparison or an
unrelated fresh frame isometry.  Equation-chamber reductions are recorded
separately and are not inferred from these raw graph coordinates.

This is the first (direct-edge) layer of the requested target-specific
meet-in-the-middle audit.  It does not enumerate new intermediate neighbours.
"""

import csv
import hashlib
import json
from pathlib import Path

from sage.all import (
    ZZ,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    matrix,
    vector,
    xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
DATA = ROOT / "elkies-k3/data/fibrations"
OUTPUT = GENERATED / "elkies-k3-other-r17-gate-b-direct-costs.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


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


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coefficients = [-value for value in coefficients]
    return vector(ZZ, coefficients)


def neighbour(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U2, -parent)
    fibre = vector(ZZ, [a, b] + list(coordinates))
    assert a * b == qnorm
    assert coordinates * parent * coordinates == 2 * qnorm
    assert fibre * ns * fibre == 0
    assert gcd([abs(ZZ(value)) for value in ns * fibre]) == 1
    mate = bezout_vector(list(ns * fibre))
    assert fibre * ns * mate == 1
    mate -= ZZ(mate * ns * mate) // 2 * fibre
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(ZZ, [list(fibre), list(mate)] + complement.rows())
    assert abs(transport.det()) == 1
    assert transport * ns * transport.transpose() == block_diagonal_matrix(
        U2, -child
    )
    return child, transport


PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SUFFIX = GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
MANIFEST = LOCAL / "q24-equation-d13-to-pinned-r17.json"
OTHER = GENERATED / "elkies-k3-other-rank17-candidate.json"
ALTERNATE = GENERATED / "q80-alternate-fifth-q6-rootless-transport.json"
ALTERNATE_Q4_SEARCH = (
    ROOT / "artifacts/local/q80-fifth-q4-low-degree-neighbor-search-gf73-v1.json"
)
PINNED_TO_Q80 = DATA / "kumar_q80_rootless_target_to_q80_ns_transport.txt"
Q80_FRAME = DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt"
Q80_PATH = DATA / "kumar_q80_to_rootless_path.tsv"
PHYSICAL_CERTIFICATES = (
    GENERATED / "elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json",
    GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json",
    GENERATED / "elkies-k3-h3-q4o164-c8-q8o376-4a1-certificate.json",
    GENERATED / "elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json",
)
PHYSICAL_ROUTE = (
    GENERATED
    / "elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-"
    "q12o5867-pinned-r17-route-certificate.json"
)
INPUTS = (
    PINNED_FRAME,
    SUFFIX,
    MANIFEST,
    OTHER,
    ALTERNATE,
    ALTERNATE_Q4_SEARCH,
    PINNED_TO_Q80,
    Q80_FRAME,
    Q80_PATH,
    PHYSICAL_ROUTE,
    *PHYSICAL_CERTIFICATES,
)


pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
suffix = json.loads(SUFFIX.read_text())
manifest = json.loads(MANIFEST.read_text())
other = json.loads(OTHER.read_text())
assert suffix["status"] == "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING"
assert manifest["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert other["status"] == "PASS_EXACT_DISTINCT_ROOTLESS_U_EMBEDDING"

# Every row of a stored basis matrix is a local basis vector in pinned NS
# coordinates.  The first row is its nef fibre and the second is O+F.
nodes = {}

equation_d13_in_pinned = matrix(
    ZZ, manifest["equation_d13_to_pinned_r17_transition"]
).inverse().change_ring(ZZ)
nodes["equation_D13"] = {
    "basis": equation_d13_in_pinned,
    "equation_status": "exact QQ",
    "source": str(MANIFEST.relative_to(ROOT)),
}

for name, data in suffix["current_suffix_stages"].items():
    # These are the current-route frames.  Only stages with an actual QQ
    # equation are admitted below; lattice-only middle stages are excluded.
    if name not in {
        "current_D12",
        "current_A11",
        "current_A5A5",
        "current_4A1",
        "current_3A1",
        "current_2A1",
        "current_A1",
        "current_rootless",
    }:
        continue
    nodes[name] = {
        "basis": matrix(ZZ, data["basis_in_pinned_R17"]),
        "equation_status": "exact QQ",
        "source": str(SUFFIX.relative_to(ROOT)),
    }

# The physical q323-free branch has different markings from the fixed suffix.
# Pin its equation-A11 basis through its own exact rootless endpoint rather
# than conflating the fixed-route A11 basis with the equation A11 basis.
physical_route = json.loads(PHYSICAL_ROUTE.read_text())
assert physical_route["status"] == (
    "PASS_EXACT_Q323_FREE_Q4O1584_Q4O164_Q8O376_Q12O5867_PINNED_R17_ROUTE"
)
rootless_child_to_pinned_frame = matrix(
    ZZ, physical_route["rootless_child_to_pinned_r17_isometry"]
)
pinned_basis_in_rootless_child = block_diagonal_matrix(
    identity_matrix(ZZ, 2), rootless_child_to_pinned_frame
)
rootless_child_basis_in_pinned = (
    pinned_basis_in_rootless_child.inverse().change_ring(ZZ)
)
equation_a11_in_pinned = (
    matrix(ZZ, physical_route["rootless_child_to_equation_A11_basis"])
    * rootless_child_basis_in_pinned
)
assert (
    equation_a11_in_pinned * g_pinned * equation_a11_in_pinned.transpose()
    == matrix(ZZ, suffix["root_adapted_hub_to_equation_A11_basis"])
    * g_pinned
    * matrix(ZZ, suffix["root_adapted_hub_to_equation_A11_basis"]).transpose()
)
physical = [json.loads(path.read_text()) for path in PHYSICAL_CERTIFICATES]
first_physical = physical[0]
source_to_child_stored = matrix(ZZ, first_physical["source_to_child_basis"])
child_in_equation_a11 = matrix(
    ZZ, first_physical["equation_A11_to_child_basis"]
)
nodes["physical_q4o208_3A3"] = {
    "basis": source_to_child_stored.inverse().change_ring(ZZ)
    * child_in_equation_a11
    * equation_a11_in_pinned,
    "equation_status": "exact QQ",
    "source": str(PHYSICAL_CERTIFICATES[0].relative_to(ROOT)),
}
physical_names = (
    "physical_q4o1584_D4_A3_3A1",
    "physical_q4o164_2A3_2A1",
    "physical_q8o376_4A1",
    "physical_q12o5867_rootless",
)
for name, certificate, path in zip(
    physical_names, physical, PHYSICAL_CERTIFICATES
):
    nodes[name] = {
        "basis": matrix(ZZ, certificate["equation_A11_to_child_basis"])
        * equation_a11_in_pinned,
        "equation_status": "exact QQ",
        "source": str(path.relative_to(ROOT)),
    }

for name, record in nodes.items():
    basis = record["basis"]
    assert abs(basis.det()) == 1, name
    local_gram = basis * g_pinned * basis.transpose()
    if not (local_gram[:2, :2] == U2 and local_gram[:2, 2:] == 0):
        raise AssertionError(
            f"{name}: leading Gram={local_gram[:2, :2]}, "
            f"cross={local_gram[:2, 2:]}"
        )

# Reconstruct the alternate A1 parent using the exact q80 witnesses.  This
# avoids inferring its fibre from the rootless frame or from an ADE label.
# The transition's first row is the raw q4 divisor, so it must not itself be
# priced as the physical target fibre.  The certified one-reflection chamber
# reduction below is the actual nef A1 fibre ray.
q80 = load_matrix(Q80_FRAME)
with Q80_PATH.open() as handle:
    q80_steps = list(csv.DictReader(handle, delimiter="\t"))
assert len(q80_steps) == 6
q80_to_fourth = identity_matrix(ZZ, 19)
frame = q80
for row in q80_steps[:4]:
    frame, transition = neighbour(
        frame,
        ZZ(row["q"]),
        ZZ(row["a"]),
        ZZ(row["b"]),
        vector(ZZ, [ZZ(value) for value in row["v"].split(",")]),
    )
    q80_to_fourth = transition * q80_to_fourth
alternate_q4_v = vector(
    ZZ, (-9, 8, -11, 10, -4, 0, 5, 1, -6, 6, 1, -2, -1, -1, 1, 2, 0)
)
alternate_a1, alternate_q4_transport = neighbour(
    frame, ZZ(4), ZZ(2), ZZ(2), alternate_q4_v
)
alternate_a1_to_q80 = alternate_q4_transport * q80_to_fourth
pinned_to_q80 = load_matrix(PINNED_TO_Q80)
alternate_a1_to_pinned = (
    alternate_a1_to_q80 * pinned_to_q80.inverse()
).change_ring(ZZ)
assert abs(alternate_a1_to_pinned.det()) == 1
assert (
    alternate_a1_to_pinned * g_pinned * alternate_a1_to_pinned.transpose()
    == block_diagonal_matrix(U2, -alternate_a1)
)

physical_rootless_in_a1 = vector(
    ZZ,
    (
        3, 2, -1, -2, 4, 2, -1, 2, 1, -1,
        1, 0, 1, -1, 1, 0, 0, 0, 0,
    ),
)
targets = {
    # Graph matching uses the exact neighbour-frame fibre rays.  Equation
    # pricing must separately use the displayed chamber-reduced divisors.
    "alternate_A1_MW16": vector(ZZ, alternate_a1_to_pinned.row(0)),
    "alternate_rootless_MW17": vector(
        ZZ,
        matrix(
            ZZ,
            other["pinned_ns_coordinates"]["alternate_basis_to_pinned_ns"],
        ).row(0),
    ),
}
assert all(target * g_pinned * target == 0 for target in targets.values())

records = []
for node_name, node in sorted(nodes.items()):
    basis = node["basis"]
    pinned_in_node = basis.inverse().change_ring(ZZ)
    for target_name, target_pinned in targets.items():
        target = vector(ZZ, target_pinned * pinned_in_node)
        assert target * (basis * g_pinned * basis.transpose()) * target == 0
        old_degree = ZZ(target[1])
        zero_pairing = ZZ(target[0] - target[1])
        records.append(
            {
                "equation_node": node_name,
                "target": target_name,
                "target_fibre_in_node": list(map(int, target)),
                "old_fibre_degree": int(old_degree),
                "target_zero_pairing": int(zero_pairing),
                "presentation_q": int(target[0] * target[1]),
                "coordinate_growth_max": int(max(abs(value) for value in target)),
                "equation_status": node["equation_status"],
                "source": node["source"],
            }
        )

rankings = {
    target: sorted(
        (record for record in records if record["target"] == target),
        key=lambda record: (
            record["old_fibre_degree"],
            record["presentation_q"],
            record["coordinate_growth_max"],
            record["equation_node"],
        ),
    )
    for target in targets
}

payload = {
    "schema": "elkies-k3.other-r17-gate-b-direct-costs.v1",
    "status": "PASS_EXACT_OTHER_R17_GATE_B_DIRECT_COST_AUDIT",
    "search_scope": {
        "type": "direct-edge layer of target-specific meet-in-the-middle audit",
        "equation_explicit_node_count": len(nodes),
        "corridor_boundary": (
            "All equation-explicit nodes at and after the equation-D13 common "
            "marking are included. The earlier H3 and E8+E6 equations are "
            "represented by their first common pinned node equation_D13; this "
            "layer does not silently identify the separately marked first-q8 "
            "D13 presentation with that equation frame."
        ),
        "targets": list(targets),
        "neighbour_enumeration": "not performed in this layer",
        "physical_nef_interpretation": (
            "The graph targets are the exact fibre rays of the selected "
            "neighbour frames in the repository's fixed integral marking. "
            "They are appropriate for literal meet-in-the-middle matching. "
            "Physical equation pricing is separate: the alternate q4 has a "
            "degree-47 CM24 specialization, and the final q6 has a certified "
            "two-reflection parent-chamber reduction. Source-local RR scores "
            "are not inferred from raw target.O here."
        ),
    },
    "targets_in_pinned_ns": {
        name: list(map(int, value)) for name, value in targets.items()
    },
    "target_graph_bases_in_pinned_ns": {
        "alternate_A1_MW16": rows(alternate_a1_to_pinned),
        "alternate_rootless_MW17": other["pinned_ns_coordinates"]
        ["alternate_basis_to_pinned_ns"],
    },
    "parent_chamber_reduced_targets": {
        "alternate_rootless_MW17_in_alternate_A1_frame": list(
            map(int, physical_rootless_in_a1)
        ),
        "alternate_A1_MW16_equation_boundary": {
            "generic_graph_ray_in_fourth_frame": [
                2, 2, *map(int, alternate_q4_v)
            ],
            "CM24_specialization_old_fibre_degree": 47,
            "CM24_specialization_old_zero_pairing": 15,
            "source": str(ALTERNATE_Q4_SEARCH.relative_to(ROOT)),
        },
        "boundary": (
            "The alternate q4 graph ray has an exact generic A1/MW16 child, "
            "but its physical equation divisor is presently known through "
            "the degree-47 CM24 specialization rather than a rank-19 Weyl "
            "word. The rootless vector is the certified alternate final-q6 "
            "parent-chamber reduction."
        ),
    },
    "equation_nodes": {
        name: {
            "basis_in_pinned_ns": rows(record["basis"]),
            "equation_status": record["equation_status"],
            "source": record["source"],
        }
        for name, record in sorted(nodes.items())
    },
    "alternate_A1_basis_in_pinned_ns": rows(alternate_a1_to_pinned),
    "records": records,
    "rankings": rankings,
    "conclusion": (
        "This artifact prices only literal direct moves. A negative result does "
        "not exclude a two-edge low-degree bridge; the next Gate-B layer must "
        "enumerate bounded physically-nef shells from both sides."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): file_sha256(path) for path in INPUTS
        },
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/audit_other_r17_gate_b_direct_costs.sage"
    ),
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target_name, ranking in rankings.items():
    best = ranking[0]
    print(
        "OTHERR17GATEB|target={}|nodes={}|best_node={}|degree={}|q={}|status=PASS".format(
            target_name,
            len(nodes),
            best["equation_node"],
            best["old_fibre_degree"],
            best["presentation_q"],
        ),
        flush=True,
    )
print(f"OTHERR17GATEB|output={OUTPUT}|status=PASS_ARTIFACT_WRITE", flush=True)
