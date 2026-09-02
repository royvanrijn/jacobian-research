#!/usr/bin/env sage-python
"""Search same-NS fibration graphs by actual low-pole compiler cost.

This is a bounded route-discovery driver around
``search_root_adapted_weyl_neighbors_targeted.sage``.  An expanded shell is
exact and exhaustive for its declared ``(q, old degree)`` pair unless the
optional Mordell--Weil-vector cap is active.  Candidate edges are retained
only when the physical finite/affine component gate, the
all-section gate, and Proposition C2's finite horizontal-wall gate pass.

The path cost is

    (maximum horizontal P.O., maximum q, maximum old-fibre degree,
     edge count, total physical Weyl repairs).

Horizontal P.O. is not inferred from q.  It is computed from the exact
minimum frame norm in the edge's Mordell--Weil quotient class.  The current
search deliberately retains only already-nef presentations, so every stored
edge has zero physical Weyl repairs.  This makes the last cost coordinate
exact while leaving repaired presentations to a later widening.

The beam pruning makes graph coverage bounded; a miss is not an obstruction.
Every emitted route nevertheless retains exact unimodular transports and an
integral isometry from its terminal frame to the requested target of the same
root rank.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from sage.all import QQ, ZZ, block_diagonal_matrix, gcd, matrix, pari, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))
NEIGHBOR_SCRIPT = ROOT / "elkies-k3/scripts/search_root_adapted_weyl_neighbors_targeted.sage"
SAGE_PYTHON = Path("/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python")
FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
FOUNDRY_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
GOLAY_TARGET = ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
GOLAY_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
K304B_TARGET = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-source-search-target-v1.json"
)
K304B_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json"
)
K36CE_TARGET = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-6ce16abb9de3c7c5-source-search-target-partner1-lattice-only-v1.json"
)
K36CE_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
)
K314AD_TARGET = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-14ad03cd7c1848b2-source-search-target-partner1-lattice-only-v1.json"
)
K314AD_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-14ad03cd7c1848b2-semistable-mw0-2-sources-large-a-partner1-v1.json"
)
K3CF7F_TARGET = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-cf7f6c91a3a40d32-source-search-target-partner2-lattice-only-v1.json"
)
K3CF7F_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-cf7f6c91a3a40d32-semistable-mw0-2-sources-large-a-partner2-v1.json"
)


CASE_DATA = {
    "golay720": {
        "source_artifact": GOLAY_SOURCES,
        "source_id": "G720-S0128",
        "target_ids": ("G720-F001",),
        "default_output": ROOT
        / "artifacts/generated-results/elkies-k3-golay-det720-same-ns-compiler-routes-v1.json",
    },
    "ns0031": {
        "source_artifact": FOUNDRY_SOURCES,
        "source_id": "NS0031-S001",
        "target_ids": (
            "NS0031-F002",
            "NS0031-F003",
            "NS0031-F004",
            "NS0031-F006",
            "NS0031-F007",
            "NS0031-F008",
            "NS0031-F011",
            "NS0031-F015",
            "NS0031-F017",
            "NS0031-F018",
            "NS0031-F019",
            "NS0031-F026",
            "NS0031-F028",
        ),
        "preferred_target_ids": ("NS0031-F017", "NS0031-F018"),
        "default_output": ROOT
        / "artifacts/generated-results/elkies-k3-ns0031-same-ns-compiler-routes-v1.json",
    },
    "k304b": {
        "source_artifact": K304B_SOURCES,
        "source_id": "K3-04b86146cc6b284b-S0160",
        "target_ids": ("K3-04b86146cc6b284b-F001",),
        "target_artifact": K304B_TARGET,
        "default_output": ROOT
        / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-same-ns-compiler-routes-v1.json",
    },
    "k36ce": {
        "source_artifact": K36CE_SOURCES,
        "source_id": "K3-6ce16abb9de3c7c5-S0008",
        "target_ids": ("K3-6ce16abb9de3c7c5-F001",),
        "target_artifact": K36CE_TARGET,
        "default_output": ROOT
        / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-same-ns-compiler-routes-v1.json",
    },
    "k314ad": {
        "source_artifact": K314AD_SOURCES,
        "source_id": "K3-14ad03cd7c1848b2-S0050",
        "target_ids": ("K3-14ad03cd7c1848b2-F001",),
        "target_artifact": K314AD_TARGET,
        "default_output": ROOT
        / "artifacts/generated-results/elkies-k3-k3-14ad03cd7c1848b2-same-ns-compiler-routes-v1.json",
    },
    "k3cf7f": {
        "source_artifact": K3CF7F_SOURCES,
        "source_id": "K3-cf7f6c91a3a40d32-S0223",
        "target_ids": ("K3-cf7f6c91a3a40d32-F001",),
        "target_artifact": K3CF7F_TARGET,
        "default_output": ROOT
        / "artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-same-ns-compiler-routes-v1.json",
    },
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rows(value) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in value.rows()]


def matrix_key(value) -> str:
    raw = json.dumps(rows(value), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
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
    result = []
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
        result.append(max(candidates, key=lambda root: sum(root)))
    return tuple(result)


def edge_context(frame, root_rank):
    root = matrix(QQ, frame[:root_rank, :root_rank])
    return {
        "frame": frame,
        "root_rank": root_rank,
        "highest_roots": highest_roots(frame[:root_rank, :root_rank]),
        "frame_lattice": IntegralLattice(frame),
        "root_inverse": root.inverse(),
        "coupling": matrix(QQ, frame[:root_rank, root_rank:]),
        "root_lattice": IntegralLattice(root.change_ring(ZZ)),
    }


def physical_nef_profile(fibre, context):
    frame = context["frame"]
    root_rank = context["root_rank"]
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    labels = w * frame[:, :root_rank]
    affine = [degree - top * labels for top in context["highest_roots"]]
    center = vector(QQ, w) / degree
    closest = vector(ZZ, next(context["frame_lattice"].enumerate_close_vectors(center)))
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
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    walls = []
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * frame * w.column()
        augmented = block_diagonal_matrix(degree**2 * frame, matrix(ZZ, [[1]]))
        augmented[:17, 17] = cross
        augmented[17, :17] = cross.transpose()
        augmented[17, 17] += m**2 * (w * frame * w)
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


def section_pole_order(context, tail_values):
    if not tail_values or not any(tail_values):
        raise ValueError("neighbor horizontal has zero MW quotient class")
    tail = vector(QQ, tail_values)
    target = -(context["root_inverse"] * context["coupling"] * tail)
    closest = vector(ZZ, next(context["root_lattice"].enumerate_close_vectors(target)))
    representative = vector(ZZ, list(closest) + list(map(ZZ, tail_values)))
    norm = ZZ(representative * context["frame"] * representative)
    if norm < 4 or norm % 2:
        raise ArithmeticError("invalid minimum section-frame norm")
    return int((norm - 4) // 2), int(norm), list(map(int, closest))


def integral_isometry(left, right):
    raw = pari(left).qfisom(pari(right))
    if raw == 0:
        return None
    candidate = matrix(ZZ, raw)
    for value in (candidate, candidate.transpose()):
        if value * left * value.transpose() == right:
            return value
        if value * right * value.transpose() == left:
            return value.inverse().change_ring(ZZ)
    raise ArithmeticError("PARI returned an unrecognized isometry orientation")


def load_case(case_name):
    config = CASE_DATA[case_name]
    source_payload = json.loads(config["source_artifact"].read_text())
    source_row = next(
        row for row in source_payload["sources"] if row["source_id"] == config["source_id"]
    )
    source = source_row["source"]
    source_frame = matrix(ZZ, source["root_adapted_gram"])
    source_root_rank = int(source["root_rank"])
    if config.get("target_artifact") is not None:
        target_payload = json.loads(config["target_artifact"].read_text())
        target_row = target_payload["frame"]
        targets = {
            target_row["frame_id"]: {
                "frame": matrix(ZZ, target_row["gram"]),
                "root_rank": int(target_row["root_rank"]),
            }
        }
        target_inputs = {
            relative(config["target_artifact"]): digest(config["target_artifact"])
        }
    elif case_name == "golay720":
        target_payload = json.loads(GOLAY_TARGET.read_text())
        targets = {
            "G720-F001": {
                "frame": matrix(ZZ, target_payload["frame"]["gram"]),
                "root_rank": int(target_payload["frame"].get("root_rank", 0)),
            }
        }
        target_inputs = {relative(GOLAY_TARGET): digest(GOLAY_TARGET)}
    else:
        foundry = json.loads(FOUNDRY.read_text())
        ns = next(row for row in foundry["ns_classes"] if row["ns_id"] == "NS0031")
        targets = {
            row["frame_id"]: {
                "frame": matrix(ZZ, row["gram"]),
                "root_rank": int(row["root_rank"]),
            }
            for row in ns["frames"]
            if row["frame_id"] in config["target_ids"]
        }
        target_inputs = {relative(FOUNDRY): digest(FOUNDRY)}
    if set(targets) != set(config["target_ids"]):
        raise ValueError("requested target frame is absent")
    return config, source, source_frame, source_root_rank, targets, target_inputs


def write_frame(path: Path, frame):
    path.write_text("\n".join(" ".join(map(str, row)) for row in frame.rows()) + "\n")


def expand_state(
    state, q_values, degree, workdir, expansion_index, mw_vector_cap=None
):
    frame_path = workdir / f"frame-{expansion_index:05d}.txt"
    output_path = workdir / f"neighbors-{expansion_index:05d}.json"
    write_frame(frame_path, state["frame"])
    command = [
        str(SAGE_PYTHON),
        str(NEIGHBOR_SCRIPT),
        "--frame",
        str(frame_path),
        "--root-rank",
        str(state["root_rank"]),
        "--degree",
        str(degree),
        "--adapt-mw-at-least",
        "1",
        "--rank-growth-only",
        "--output",
        str(output_path),
    ]
    for q in q_values:
        if q % degree == 0:
            command.extend(("--q", str(q)))
    if mw_vector_cap is not None:
        command.extend(("--mw-vector-cap", str(mw_vector_cap)))
    if not any(command[index] == "--q" for index in range(len(command))):
        return [], []
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stdout + "\n" + completed.stderr)
    payload = json.loads(output_path.read_text())
    return payload["neighbors"], payload["summaries"]


def route_cost(edges):
    return (
        max((edge["horizontal_P_dot_O"] for edge in edges), default=-1),
        max((edge["q"] for edge in edges), default=0),
        max((edge["old_fibre_degree"] for edge in edges), default=0),
        len(edges),
        sum(edge["physical_weyl_repairs"] for edge in edges),
    )


def search_case(
    case_name,
    q_values,
    degree,
    max_pole,
    beam_width,
    max_depth,
    mw_vector_cap,
    cap_from_mw_rank,
    rank_first,
    accepted_per_state_cap,
    retain_frontier_witnesses,
    resume_frontier_artifact,
    resume_root_rank,
    resume_frontier_depth,
):
    config, source, source_frame, source_root_rank, targets, target_inputs = load_case(case_name)
    source_ns = block_diagonal_matrix(U2, -source_frame)
    resume_inputs = {}
    if resume_frontier_artifact is None:
        initial_states = [
            {
                "frame": source_frame,
                "root_rank": source_root_rank,
                "to_source": matrix(ZZ, 19, 19, 1),
                "edges": [],
            }
        ]
    else:
        resume_payload = json.loads(resume_frontier_artifact.read_text())
        resume_results = [
            row for row in resume_payload["results"] if row["case"] == case_name
        ]
        if len(resume_results) != 1:
            raise ValueError("resume artifact must contain exactly one matching case")
        if resume_frontier_depth is None:
            witness = resume_results[0]["frontier_witnesses_by_root_rank"].get(
                str(resume_root_rank)
            )
            if witness is None:
                raise ValueError(
                    "requested root-rank witness is absent from resume artifact"
                )
            resume_rows = [witness]
        else:
            depth_rows = [
                row
                for row in resume_results[0]["retained_frontiers_by_depth"]
                if int(row["depth"]) == resume_frontier_depth
            ]
            if len(depth_rows) != 1:
                raise ValueError("requested retained frontier depth is absent or ambiguous")
            resume_rows = depth_rows[0]["states"]
            if resume_root_rank is not None:
                resume_rows = [
                    row
                    for row in resume_rows
                    if int(row["root_rank"]) == resume_root_rank
                ]
            if not resume_rows:
                raise ValueError("no retained frontier states match the resume filter")
        initial_states = []
        for witness in resume_rows:
            state = {
                "frame": matrix(ZZ, witness["frame"]),
                "root_rank": int(witness["root_rank"]),
                "to_source": matrix(ZZ, witness["composed_transport_to_source"]),
                "edges": witness["edges"],
            }
            resumed_ns = block_diagonal_matrix(U2, -state["frame"])
            if state["to_source"] * source_ns * state["to_source"].transpose() != resumed_ns:
                raise ArithmeticError("resume witness transport does not recover the source NS")
            initial_states.append(state)
        resume_inputs[relative(resume_frontier_artifact)] = digest(
            resume_frontier_artifact
        )
    starting_depths = {len(state["edges"]) for state in initial_states}
    if len(starting_depths) != 1:
        raise ValueError("all resumed states must have the same route depth")
    starting_depth = starting_depths.pop()
    if max_depth <= starting_depth:
        raise ValueError("maximum depth must exceed the resumed route depth")
    frontier = initial_states
    accounting = []
    hits = []
    frontier_witnesses = {}
    retained_frontiers = []
    expansion_index = 0
    seen = {
        matrix_key(state["frame"]): route_cost(state["edges"])
        for state in initial_states
    }
    with tempfile.TemporaryDirectory(prefix=f"{case_name}-compiler-route-") as temporary:
        workdir = Path(temporary)
        for depth in range(starting_depth + 1, max_depth + 1):
            candidates = []
            raw_edges = accepted_edges = 0
            shell_summaries = []
            for state_index, state in enumerate(frontier, start=1):
                print(
                    "SAME_NS_ROUTE_EXPAND|case={}|depth={}|state={}/{}|root_rank={}|status=RUNNING".format(
                        case_name,
                        depth,
                        state_index,
                        len(frontier),
                        state["root_rank"],
                    ),
                    flush=True,
                )
                expansion_index += 1
                active_cap = (
                    mw_vector_cap
                    if mw_vector_cap is not None
                    and 17 - state["root_rank"] >= cap_from_mw_rank
                    else None
                )
                records, summaries = expand_state(
                    state,
                    q_values,
                    degree,
                    workdir,
                    expansion_index,
                    active_cap,
                )
                context = edge_context(state["frame"], state["root_rank"])
                shell_summaries.extend(summaries)
                raw_edges += len(records)
                if rank_first:
                    records.sort(
                        key=lambda record: (
                            int(record.get("child_root_data", [10**9])[0]),
                            int(record.get("q", 10**9)),
                            int(record.get("orbit_index", 10**9)),
                        )
                    )
                state_accepted = 0
                for record in records:
                    if "child_root_adapted_frame" not in record:
                        continue
                    fibre = vector(ZZ, record["fiber"])
                    if gcd(tuple(block_diagonal_matrix(U2, -state["frame"]) * fibre)) != 1:
                        continue
                    nef = physical_nef_profile(fibre, context)
                    if not nef["passes"]:
                        continue
                    pole, section_norm, closest_root = section_pole_order(
                        context, record["mw_projection"]
                    )
                    if pole > max_pole:
                        continue
                    walls = negative_horizontal_walls(fibre, state["frame"])
                    if walls:
                        continue
                    child = matrix(ZZ, record["child_root_adapted_frame"])
                    adaptation = matrix(ZZ, record["child_root_adapted_basis"])
                    neighbor_basis = matrix(ZZ, record["neighbor_basis"])
                    transport = block_diagonal_matrix(
                        matrix(ZZ, 2, 2, 1), adaptation
                    ) * neighbor_basis
                    child_ns = block_diagonal_matrix(U2, -child)
                    if abs(transport.det()) != 1:
                        raise ArithmeticError("edge transport is not unimodular")
                    if transport * block_diagonal_matrix(U2, -state["frame"]) * transport.transpose() != child_ns:
                        raise ArithmeticError("edge transport Gram identity failed")
                    to_source = transport * state["to_source"]
                    if to_source * source_ns * to_source.transpose() != child_ns:
                        raise ArithmeticError("composed source transport failed")
                    edge = {
                        "edge_index": depth,
                        "q": int(record["q"]),
                        "factor_order": list(map(int, record["factor_order"])),
                        "old_fibre_degree": int(record["old_fiber_degree"]),
                        "orbit_index": int(record["orbit_index"]),
                        "horizontal_P_dot_O": pole,
                        "horizontal_minimum_frame_norm": section_norm,
                        "horizontal_mw_quotient": list(map(int, record["mw_projection"])),
                        "horizontal_closest_root_coordinates": closest_root,
                        "physical_weyl_repairs": 0,
                        "source_root_type": None if not state["edges"] else state["edges"][-1]["target_root_type"],
                        "source_root_rank": int(state["root_rank"]),
                        "source_mw_rank": 17 - int(state["root_rank"]),
                        "target_root_type": record["child_ade"],
                        "target_root_rank": int(record["child_root_data"][0]),
                        "target_mw_rank": int(record["child_mw_rank"]),
                        "fibre": list(map(int, fibre)),
                        "physical_nef_profile": nef,
                        "negative_horizontal_walls": walls,
                        "edge_transport_child_to_parent": rows(transport),
                    }
                    edges = state["edges"] + [edge]
                    accepted_edges += 1
                    state_accepted += 1
                    child_state = {
                        "frame": child,
                        "root_rank": int(record["child_root_data"][0]),
                        "to_source": to_source,
                        "edges": edges,
                    }
                    if retain_frontier_witnesses:
                        witness = {
                            "root_rank": child_state["root_rank"],
                            "mw_rank": 17 - child_state["root_rank"],
                            "cost": list(route_cost(edges)),
                            "frame": rows(child),
                            "composed_transport_to_source": rows(to_source),
                            "edges": edges,
                        }
                        previous = frontier_witnesses.get(str(child_state["root_rank"]))
                        if previous is None or tuple(witness["cost"]) < tuple(previous["cost"]):
                            frontier_witnesses[str(child_state["root_rank"])] = witness
                    matching_targets = {
                        target_id: target
                        for target_id, target in targets.items()
                        if target["root_rank"] == child_state["root_rank"]
                    }
                    matched_target = False
                    if matching_targets:
                        for target_id, target in matching_targets.items():
                            isometry = integral_isometry(child, target["frame"])
                            if isometry is None:
                                continue
                            matched_target = True
                            hits.append(
                                {
                                    "target_frame_id": target_id,
                                    "cost": list(route_cost(edges)),
                                    "edges": edges,
                                    "terminal_frame": rows(child),
                                    "terminal_to_target_frame_isometry": rows(isometry),
                                    "composed_transport_terminal_to_source": rows(to_source),
                                }
                            )
                    if matched_target:
                        continue
                    key = matrix_key(child)
                    cost = route_cost(edges)
                    if key in seen and seen[key] <= cost:
                        continue
                    seen[key] = cost
                    candidates.append(child_state)
                    if (
                        accepted_per_state_cap is not None
                        and state_accepted >= accepted_per_state_cap
                    ):
                        break
            candidates.sort(
                key=lambda state: (
                    (
                        state["root_rank"],
                        route_cost(state["edges"]),
                        matrix_key(state["frame"]),
                    )
                    if rank_first
                    else (
                        route_cost(state["edges"]),
                        state["root_rank"],
                        matrix_key(state["frame"]),
                    )
                )
            )
            # Preserve some route-shape diversity instead of allowing one ADE
            # profile to consume the entire beam.
            eligible_by_max_q = defaultdict(list)
            per_profile = defaultdict(int)
            for state in candidates:
                cost = route_cost(state["edges"])
                profile = (
                    state["root_rank"],
                    state["edges"][-1]["target_root_type"],
                    cost[:3],
                )
                if per_profile[profile] >= 3:
                    continue
                per_profile[profile] += 1
                eligible_by_max_q[cost[1]].append(state)
            next_frontier = []
            bucket_index = 0
            ordered_buckets = [
                eligible_by_max_q[value] for value in sorted(eligible_by_max_q)
            ]
            while len(next_frontier) < beam_width:
                added = False
                for bucket in ordered_buckets:
                    if bucket_index < len(bucket):
                        next_frontier.append(bucket[bucket_index])
                        added = True
                        if len(next_frontier) == beam_width:
                            break
                if not added:
                    break
                bucket_index += 1
            if retain_frontier_witnesses:
                retained_frontiers.append(
                    {
                        "depth": depth,
                        "states": [
                            {
                                "root_rank": state["root_rank"],
                                "mw_rank": 17 - state["root_rank"],
                                "cost": list(route_cost(state["edges"])),
                                "frame": rows(state["frame"]),
                                "composed_transport_to_source": rows(state["to_source"]),
                                "edges": state["edges"],
                            }
                            for state in next_frontier
                        ],
                    }
                )
            accounting.append(
                {
                    "depth": depth,
                    "frontier_in": len(frontier),
                    "raw_neighbor_edges": raw_edges,
                    "accepted_zero_repair_edges": accepted_edges,
                    "distinct_candidates_before_beam": len(candidates),
                    "frontier_out": len(next_frontier),
                    "best_root_rank": min(
                        (state["root_rank"] for state in next_frontier), default=None
                    ),
                    "hits_so_far": len(hits),
                    "shell_summaries": shell_summaries,
                }
            )
            print(
                "SAME_NS_ROUTE|case={}|depth={}|frontier={}|raw={}|accepted={}|"
                "next={}|best_root_rank={}|hits={}|status=RUNNING".format(
                    case_name,
                    depth,
                    len(frontier),
                    raw_edges,
                    accepted_edges,
                    len(next_frontier),
                    accounting[-1]["best_root_rank"],
                    len(hits),
                ),
                flush=True,
            )
            if hits or not next_frontier:
                break
            frontier = next_frontier
    hits.sort(key=lambda hit: (tuple(hit["cost"]), hit["target_frame_id"]))
    best_by_target = {}
    for hit in hits:
        best_by_target.setdefault(hit["target_frame_id"], hit)
    return {
        "case": case_name,
        "source": {
            "source_id": config["source_id"],
            "root_type": source["root_type"],
            "root_rank": int(source["root_rank"]),
            "mw_rank": int(source["mw_rank_for_rho_19"]),
            "frame_gram_sha256": source["gram_sha256"],
        },
        "requested_targets": list(config["target_ids"]),
        "search": {
            "q_values": list(q_values),
            "old_fibre_degree": degree,
            "maximum_horizontal_P_dot_O": max_pole,
            "beam_width": beam_width,
            "maximum_depth": max_depth,
            "starting_depth": starting_depth,
            "resumed_root_rank": (
                resume_root_rank if resume_frontier_artifact is not None else None
            ),
            "resumed_frontier_depth": resume_frontier_depth,
            "resumed_state_count": (
                len(initial_states) if resume_frontier_artifact is not None else None
            ),
            "mw_vector_cap": mw_vector_cap,
            "mw_vector_cap_from_mw_rank": (
                cap_from_mw_rank if mw_vector_cap is not None else None
            ),
            "accepted_edge_cap_per_expanded_state": accepted_per_state_cap,
            "zero_repair_presentations_only": True,
            "beam_priority": (
                "root_rank_then_compiler_cost"
                if rank_first
                else "compiler_cost_then_root_rank"
            ),
            "cost_key_order": [
                "maximum_horizontal_P_dot_O",
                "maximum_q",
                "maximum_old_fibre_degree",
                "edge_count",
                "total_physical_weyl_repairs",
            ],
        },
        "accounting": accounting,
        "best_routes_by_target": best_by_target,
        **(
            {
                "frontier_witnesses_by_root_rank": frontier_witnesses,
                "retained_frontiers_by_depth": retained_frontiers,
            }
            if retain_frontier_witnesses
            else {}
        ),
        "inputs": {
            relative(config["source_artifact"]): digest(config["source_artifact"]),
            relative(NEIGHBOR_SCRIPT): digest(NEIGHBOR_SCRIPT),
            **resume_inputs,
            **target_inputs,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=tuple(CASE_DATA), action="append")
    parser.add_argument("--q", type=int, action="append", default=[])
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--max-pole", type=int, default=0)
    parser.add_argument("--beam-width", type=int, default=12)
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--mw-vector-cap", type=int)
    parser.add_argument("--cap-from-mw-rank", type=int, default=7)
    parser.add_argument(
        "--accepted-per-state-cap",
        type=int,
        help=(
            "after rank-first ordering, stop a state expansion after this many "
            "fully physical accepted edges; this is an explicit additional beam cap"
        ),
    )
    parser.add_argument(
        "--rank-first",
        action="store_true",
        help=(
            "prioritize root-rank descent inside the declared pole/q bounds; "
            "the default prioritizes compiler cost before root rank"
        ),
    )
    parser.add_argument(
        "--retain-frontier-witnesses",
        action="store_true",
        help="retain the cheapest full marked state encountered at every root rank",
    )
    parser.add_argument(
        "--resume-frontier-artifact",
        type=Path,
        help="resume from a retained full marked state in an earlier route artifact",
    )
    parser.add_argument(
        "--resume-root-rank",
        type=int,
        help="root-rank witness or retained-frontier filter to resume",
    )
    parser.add_argument(
        "--resume-frontier-depth",
        type=int,
        help="resume all retained states at this depth, optionally filtered by root rank",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    cases = tuple(args.case or CASE_DATA)
    q_values = tuple(sorted(set(args.q or (4,))))
    if args.degree <= 0 or args.max_pole < 0 or args.beam_width <= 0 or args.max_depth <= 0:
        parser.error("degree, beam width, and depth must be positive; max pole nonnegative")
    if args.mw_vector_cap is not None and args.mw_vector_cap <= 0:
        parser.error("MW vector cap must be positive")
    if args.accepted_per_state_cap is not None and args.accepted_per_state_cap <= 0:
        parser.error("accepted-per-state cap must be positive")
    if not 1 <= args.cap_from_mw_rank <= 16:
        parser.error("cap-from MW rank must lie between 1 and 16")
    if any(q <= 0 or q % args.degree for q in q_values):
        parser.error("every q must be a positive multiple of the old-fibre degree")
    if len(cases) > 1 and args.output is not None:
        parser.error("--output is available only with one --case")
    if args.resume_frontier_artifact is None and (
        args.resume_root_rank is not None or args.resume_frontier_depth is not None
    ):
        parser.error("resume filters require --resume-frontier-artifact")
    if args.resume_frontier_artifact is not None and (
        args.resume_root_rank is None and args.resume_frontier_depth is None
    ):
        parser.error("resume requires --resume-root-rank or --resume-frontier-depth")
    if args.resume_frontier_artifact is not None and len(cases) != 1:
        parser.error("frontier resume is available only with one --case")
    results = [
        search_case(
            case_name,
            q_values,
            args.degree,
            args.max_pole,
            args.beam_width,
            args.max_depth,
            args.mw_vector_cap,
            args.cap_from_mw_rank,
            args.rank_first,
            args.accepted_per_state_cap,
            args.retain_frontier_witnesses,
            args.resume_frontier_artifact.resolve()
            if args.resume_frontier_artifact is not None
            else None,
            args.resume_root_rank,
            args.resume_frontier_depth,
        )
        for case_name in cases
    ]
    payload = {
        "schema": "elkies-k3.same-ns-compiler-route-search.v1",
        "status": (
            "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_HIT"
            if any(result["best_routes_by_target"] for result in results)
            else "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_EMPTY"
        ),
        "results": results,
        "proof_boundary": {
            "proved": (
                "Every retained edge is exact and has primitive fibre, complete physical "
                "component/all-section/finite-horizontal-wall gates, "
                "zero physical Weyl repairs, and bidirectional unimodular NS transport. Every "
                "reported terminal is integrally isometric to the named target frame."
            ),
            "not_proved": (
                "Beam pruning is not a complete graph search. Repaired presentations, other q "
                "values, other old-fibre degrees, equation-level pencils, and global route "
                "optimality remain open. A lattice route does not prove a QQ equation lift."
            ),
        },
        "reproduce": " ".join(
            [
                str(SAGE_PYTHON),
                relative(Path(__file__)),
                *(sum((["--case", case_name] for case_name in cases), [])),
                *(sum((["--q", str(q)] for q in q_values), [])),
                "--degree",
                str(args.degree),
                "--max-pole",
                str(args.max_pole),
                "--beam-width",
                str(args.beam_width),
                "--max-depth",
                str(args.max_depth),
                *(
                    [
                        "--mw-vector-cap",
                        str(args.mw_vector_cap),
                        "--cap-from-mw-rank",
                        str(args.cap_from_mw_rank),
                    ]
                    if args.mw_vector_cap is not None
                    else []
                ),
                *(["--rank-first"] if args.rank_first else []),
                *(
                    ["--retain-frontier-witnesses"]
                    if args.retain_frontier_witnesses
                    else []
                ),
                *(
                    [
                        "--resume-frontier-artifact",
                        relative(args.resume_frontier_artifact),
                        *(
                            ["--resume-root-rank", str(args.resume_root_rank)]
                            if args.resume_root_rank is not None
                            else []
                        ),
                        *(
                            [
                                "--resume-frontier-depth",
                                str(args.resume_frontier_depth),
                            ]
                            if args.resume_frontier_depth is not None
                            else []
                        ),
                    ]
                    if args.resume_frontier_artifact is not None
                    else []
                ),
                *(
                    ["--accepted-per-state-cap", str(args.accepted_per_state_cap)]
                    if args.accepted_per_state_cap is not None
                    else []
                ),
            ]
        ),
    }
    if args.output is not None:
        output = args.output.resolve()
    elif len(cases) == 1:
        output = CASE_DATA[cases[0]]["default_output"]
    else:
        output = (
            ROOT
            / "artifacts/generated-results/elkies-k3-golay720-ns0031-same-ns-compiler-routes-v1.json"
        )
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit("same-NS compiler-route artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "SAME_NS_ROUTE|cases={}|targets_hit={}|status={}".format(
            len(cases),
            sum(len(result["best_routes_by_target"]) for result in results),
            payload["status"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
