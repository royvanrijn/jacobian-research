#!/usr/bin/env sage -python
"""Certify a lattice-foundry route from a rootful source to a rootless frame.

status: ACTIVE_PROOF
claim: exact primitive-nef neighbour edges and lossless integral NS transport
inputs: foundry database, exact source certificate, ordered neighbour searches
output: caller-selected exact JSON route ledger

The neighbour searches deliberately separate orbit discovery from this proof
step.  Here every selected edge is replayed, its physical component and
all-section gates are checked, Proposition C2's finite horizontal-wall list is
exhausted, and the full determinant-one NS markings are composed.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def display_path(path):
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
        component = []
        while todo:
            node = todo.pop()
            component.append(node)
            for other in tuple(unseen):
                if cartan[node, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def highest_roots(cartan):
    if cartan.nrows() == 0:
        return ()
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-root for root in half)
    answer = []
    for component in connected_components(cartan):
        candidates = [
            root
            for root in roots
            if all(value >= 0 for value in root)
            and all(
                index in component or root[index] == 0
                for index in range(cartan.nrows())
            )
        ]
        answer.append(max(candidates, key=lambda root: sum(root)))
    return tuple(answer)


def physical_nef_profile(fibre, frame, root_rank):
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    labels = w * frame[:, :root_rank]
    cartan = frame[:root_rank, :root_rank]
    affine = [degree - top * labels for top in highest_roots(cartan)]
    center = vector(QQ, w) / degree
    closest = vector(ZZ, next(IntegralLattice(frame).enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = degree * (distance - 2) / 2
    return {
        "component_pairings": list(map(int, labels)),
        "affine_pairings": list(map(int, affine)),
        "minimum_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "passes": bool(
            min(tuple(labels) + tuple(affine) + (ZZ(0),)) >= 0
            and minimum_section >= 0
        ),
    }


def negative_horizontal_walls(fibre, frame):
    """Enumerate all negative old-horizontal (-2)-walls (Proposition C2)."""
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    walls = []
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * frame * w.column()
        augmented = block_matrix(
            ZZ,
            [
                [degree**2 * frame, cross],
                [
                    cross.transpose(),
                    matrix(ZZ, [[m**2 * (w * frame * w) + 1]]),
                ],
            ],
        )
        result = pari(augmented).qfminim(2 * degree**2 - 1)
        normalized = set()
        for raw in matrix(ZZ, result[2]).transpose().rows():
            if abs(raw[-1]) != 1:
                continue
            value = raw if raw[-1] == 1 else -raw
            normalized.add(tuple(value))
        for value in normalized:
            x = vector(ZZ, value[:-1])
            x_norm = ZZ(x * frame * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ(
                (w * frame * w // (2 * degree)) * m
                + degree * k
                - w * frame * x
            )
            if intersection < 0:
                walls.append(
                    {
                        "old_fibre_degree": int(m),
                        "curve": [int(k), int(m)] + list(map(int, x)),
                        "intersection": int(intersection),
                    }
                )
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


def integral_isometry(left, right):
    raw = pari(left).qfisom(pari(right))
    if raw == 0:
        return None
    candidate = matrix(ZZ, raw)
    possibilities = (candidate, candidate.transpose())
    for value in possibilities:
        if value * left * value.transpose() == right:
            return value
        if value * right * value.transpose() == left:
            return value.inverse().change_ring(ZZ)
    raise AssertionError("PARI qfisom returned an unrecognized orientation")


def rootless_intrinsics(frame, norm_bound=4):
    result = pari(frame).qfminim(norm_bound)
    columns = tuple(vector(ZZ, column) for column in matrix(ZZ, result[2]).columns())
    by_norm = Counter(int(value * frame * value) for value in columns)
    theta = {"0": 1}
    theta.update({str(norm): 2 * count for norm, count in sorted(by_norm.items())})
    minimum_norm = min(by_norm) if by_norm else None
    norm_four_half = tuple(value for value in columns if value * frame * value == 4)
    pairing_histogram = Counter()
    for left_index, left in enumerate(norm_four_half):
        for right in norm_four_half[left_index + 1 :]:
            pairing_histogram[abs(int(left * frame * right))] += 1
    short_cosets = {
        tuple(int(entry % 2) for entry in value) for value in columns
    }
    return {
        "minimum_squared_norm": minimum_norm,
        "theta_coefficients_by_squared_norm_through_bound": theta,
        "theta_squared_norm_bound": norm_bound,
        "norm_four_vectors": 2 * len(norm_four_half),
        "norm_four_unoriented_pairs": len(norm_four_half),
        "automorphism_group_order": int(pari(frame).qfauto()[0]),
        "short_cosets_mod_2_hit_through_bound": len(short_cosets),
        "short_coset_squared_norm_bound": norm_bound,
        "low_height_degree_two_q2_oriented_candidates": 2 * len(norm_four_half),
        "minimal_vector_absolute_pairing_histogram_unoriented_representatives": {
            str(pairing): count for pairing, count in sorted(pairing_histogram.items())
        },
        "hermite_invariant": float(
            QQ(minimum_norm) / (QQ(frame.det()) ** (QQ(1) / 17))
        ),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", type=Path)
parser.add_argument("--database", type=Path)
parser.add_argument("--source-hunt", type=Path)
parser.add_argument("--step", action="append", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

if args.manifest is not None:
    assert not any((args.database, args.source_hunt, args.step, args.output))
    specification = json.loads(args.manifest.resolve().read_text())
    assert specification["schema"] == "elkies-k3.lattice-foundry-route-specification.v1"
    DATABASE = (ROOT / specification["database"]).resolve()
    SOURCE_HUNT = (ROOT / specification["source_hunt"]).resolve()
    STEPS = tuple((ROOT / path).resolve() for path in specification["steps"])
    EXHAUSTED_SHELLS = tuple(
        (ROOT / path).resolve()
        for path in specification.get("exhausted_cheaper_shells", ())
    )
    OUTPUT = (ROOT / specification["output"]).resolve()
    MANIFEST = args.manifest.resolve()
else:
    assert args.database is not None and args.source_hunt is not None
    assert args.step and args.output is not None
    DATABASE = args.database.resolve()
    SOURCE_HUNT = args.source_hunt.resolve()
    STEPS = tuple(path.resolve() for path in args.step)
    EXHAUSTED_SHELLS = ()
    OUTPUT = args.output.resolve()
    MANIFEST = None

database = json.loads(DATABASE.read_text())
source_hunt = json.loads(SOURCE_HUNT.read_text())
assert database["status"] == "PASS_EXACT_DECLARED_SHELL_NEW_K3_TARGETS_ROUTE_GATE_OPEN"
assert source_hunt["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_MW5_SOURCE_AND_NIEMEIER_CERTIFICATE"
assert source_hunt["niemeier_certificate"]["primitive_auxiliary_embedding"]
assert source_hunt["niemeier_certificate"]["saturated_orthogonal_complement"]
assert source_hunt["niemeier_certificate"]["complement_integrally_isometric_to_source"]

ns_id = source_hunt["target"]["ns_id"]
ns_record = next(item for item in database["ns_classes"] if item["ns_id"] == ns_id)
source_frame = matrix(ZZ, source_hunt["source"]["root_adapted_gram"])
assert abs(source_frame.det()) == ns_record["determinant"]
assert source_hunt["source"]["root_rank"] > 0

current_frame = source_frame
current_root_rank = int(source_hunt["source"]["root_rank"])
current_to_source = identity_matrix(ZZ, 19)
edge_records = []

for edge_index, path in enumerate(STEPS, start=1):
    search = json.loads(path.read_text())
    assert search["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    assert len(search["neighbors"]) == 1
    record = search["neighbors"][0]
    input_frame = load_matrix((ROOT / search["frame"]).resolve())
    assert input_frame == current_frame
    assert int(search["input_root_data"][0]) == current_root_rank

    ns = block_diagonal_matrix(U2, -current_frame)
    fibre = vector(ZZ, record["fiber"])
    assert fibre * ns * fibre == 0
    assert gcd(tuple(ns * fibre)) == 1
    assert int(record["old_fiber_degree"]) == int(fibre[1])
    assert int(record["q"]) == int(fibre[0] * fibre[1])

    nef = physical_nef_profile(fibre, current_frame, current_root_rank)
    horizontal_walls = negative_horizontal_walls(fibre, current_frame)
    assert nef["passes"] and not horizontal_walls

    raw_child = matrix(ZZ, record["child_frame"])
    adapted_child = matrix(ZZ, record["child_root_adapted_frame"])
    neighbor_basis = matrix(ZZ, record["neighbor_basis"])
    adaptation = matrix(ZZ, record["child_root_adapted_basis"])
    edge_transport = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * neighbor_basis
    child_ns = block_diagonal_matrix(U2, -adapted_child)
    assert abs(edge_transport.det()) == 1
    assert edge_transport * ns * edge_transport.transpose() == child_ns
    assert adaptation * raw_child * adaptation.transpose() == adapted_child

    child_root_rank = int(record["child_root_data"][0])
    assert int(record["child_mw_rank"]) == 17 - child_root_rank
    new_to_source = edge_transport * current_to_source
    assert new_to_source * block_diagonal_matrix(U2, -source_frame) * new_to_source.transpose() == child_ns

    edge_records.append(
        {
            "edge_index": edge_index,
            "search_artifact": display_path(path),
            "q": int(record["q"]),
            "orbit_index": int(record["orbit_index"]),
            "old_fibre_degree": int(record["old_fiber_degree"]),
            "physical_weyl_reflection_count": 0,
            "predicted_horizontal_P_dot_O": int(fibre[0] - fibre[1]),
            "source_root_type": search["input_ade"],
            "source_root_rank": current_root_rank,
            "source_mw_rank": 17 - current_root_rank,
            "target_root_type": record["child_ade"],
            "target_root_rank": child_root_rank,
            "target_mw_rank": int(record["child_mw_rank"]),
            "root_rank_mutation": child_root_rank - current_root_rank,
            "determinant_before": int(abs(current_frame.det())),
            "determinant_after": int(abs(adapted_child.det())),
            "determinant_mutation": int(abs(adapted_child.det()) - abs(current_frame.det())),
            "expected_resolved_RR_dimension": 2,
            "expected_resolved_RR_dimension_status": "equation-planning heuristic",
            "zero_changing_loop": False,
            "physical_nef_profile": nef,
            "negative_horizontal_walls": horizontal_walls,
            "primitive_isotropic": True,
            "edge_transport_child_to_parent": rows(edge_transport),
            "composed_transport_child_to_source": rows(new_to_source),
        }
    )
    current_frame = adapted_child
    current_root_rank = child_root_rank
    current_to_source = new_to_source

assert current_root_rank == 0
rootless_matches = []
for target in database["rootless_targets"]:
    if target["ns_id"] != ns_id:
        continue
    frame_record = next(
        item for item in ns_record["frames"] if item["frame_id"] == target["frame_id"]
    )
    target_gram = matrix(ZZ, frame_record["gram"])
    isometry = integral_isometry(current_frame, target_gram)
    if isometry is not None:
        rootless_matches.append(
            {
                "frame_id": target["frame_id"],
                "terminal_to_catalogue_basis": rows(isometry),
                "catalogue_invariants": target["invariants"],
            }
        )
terminal_gram_sha256 = hashlib.sha256(
    ("\n".join(" ".join(map(str, row)) for row in current_frame.rows()) + "\n").encode()
).hexdigest()
terminal_intrinsics = rootless_intrinsics(current_frame)

exhausted_shell_records = []
for path in EXHAUSTED_SHELLS:
    shell = json.loads(path.read_text())
    assert shell["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    assert not shell["neighbors"] and len(shell["summaries"]) == 1
    summary = shell["summaries"][0]
    assert summary["mw_enumeration_complete"]
    assert summary["mw_vector_cap"] is None
    assert not summary["search_stopped_early"]
    assert not summary["stream_limit_reached"]
    assert summary["stream_skip"] == 0
    assert summary["stream_tested"] == summary["screened_orbits"]
    assert summary["stream_tested"] == summary["dominant_orbits"]
    exhausted_shell_records.append(
        {
            "artifact": display_path(path),
            "input_root_type": shell["input_ade"],
            "input_mw_rank": int(shell["input_mw_rank"]),
            "q": int(summary["q"]),
            "old_fibre_degree": int(summary["factor_order"][1]),
            "tested_dominant_classes": int(summary["stream_tested"]),
            "rank_growth_found": False,
            "exact_stream_exhaustion": True,
        }
    )

q_values = [edge["q"] for edge in edge_records]
degrees = [edge["old_fibre_degree"] for edge in edge_records]
cost_vector = {
    "maximum_old_fibre_degree": max(degrees),
    "maximum_q": max(q_values),
    "sum_q": sum(q_values),
    "edge_count": len(edge_records),
    "physical_weyl_reflections": sum(
        edge["physical_weyl_reflection_count"] for edge in edge_records
    ),
    "root_rank_area": sum(edge["source_root_rank"] for edge in edge_records),
}

inputs = (
    ((MANIFEST,) if MANIFEST is not None else ())
    + (DATABASE, SOURCE_HUNT)
    + STEPS
    + EXHAUSTED_SHELLS
)
payload = {
    "schema": "elkies-k3.lattice-foundry-nef-route.v1",
    "status": "PASS_EXACT_NEW_K3_ROOTFUL_TO_ROOTLESS_NEF_ROUTE",
    "ns_id": ns_id,
    "determinant": int(ns_record["determinant"]),
    "source": {
        "root_type": source_hunt["source"]["root_type"],
        "root_rank": int(source_hunt["source"]["root_rank"]),
        "mw_rank": int(source_hunt["source"]["mw_rank_for_rho_19"]),
        "niemeier_ambient": source_hunt["niemeier_certificate"]["ambient_label"],
    },
    "terminal": {
        "root_type": "0",
        "mw_rank": 17,
        "catalogue_isometry_matches": rootless_matches,
        "novel_relative_to_declared_foundry_shell": not bool(rootless_matches),
        "frame_gram_sha256": terminal_gram_sha256,
        "intrinsics": terminal_intrinsics,
        "frame_gram": rows(current_frame),
        "composed_transport_terminal_to_source": rows(current_to_source),
    },
    "edges": edge_records,
    "exhausted_cheaper_shells": exhausted_shell_records,
    "equation_cost_vector": cost_vector,
    "search_policy": {
        "objective": "lexicographic equation cost, not graph distance",
        "key_order": [
            "maximum_old_fibre_degree",
            "maximum_q",
            "physical_weyl_reflections",
            "root_rank_area",
            "edge_count",
            "sum_q",
        ],
        "status": "selected route certified; global optimality not claimed",
    },
    "proof_boundary": {
        "proved": (
            "Every listed edge is primitive, isotropic, physically nef in the complete "
            "component/all-section/finite-horizontal-wall sense, has a unimodular integral "
            "NS transport, and composes to a rootless catalogue frame on the same exact "
            "K3 NS class."
        ),
        "not_proved": (
            "No equation, effective equation-side zero, resolved Riemann--Roch pencil, "
            "or global optimality among all neighbour routes is asserted.  The recorded "
            "resolved-RR dimensions are planning estimates only."
        ),
    },
    "inputs": {
        "paths": [display_path(path) for path in inputs],
        "sha256": {
            display_path(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    assert OUTPUT.read_text() == encoded
else:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(encoded)

print(
    "FOUNDRYROUTE|ns={}|source={}/MW{}|edges={}|max_q={}|degree={}|"
    "target=MW17|matches={}|status=PASS".format(
        ns_id,
        source_hunt["source"]["root_type"],
        source_hunt["source"]["mw_rank_for_rho_19"],
        len(edge_records),
        max(q_values),
        max(degrees),
        (
            ",".join(item["frame_id"] for item in rootless_matches)
            if rootless_matches else "NEW_FRAME_CLASS"
        ),
    ),
    flush=True,
)
