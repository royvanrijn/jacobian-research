#!/usr/bin/env sage -python
"""Probe the bounded two-edge Gate-B intersection shell by exact CVP checks.

For every equation-explicit H3 node and each alternate graph endpoint, seek a
primitive isotropic ray M which is simultaneously a low-degree/low-q neighbour
on both sides.  Candidate discovery is target-directed and bounded; every
retained match is then subjected to exact component and horizontal-wall nef
tests on both sides.
"""

import hashlib
import json
from itertools import islice
from pathlib import Path

from sage.all import (
    QQ, ZZ, RealField, block_diagonal_matrix, block_matrix, gcd,
    matrix, pari, vector,
)
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-other-r17-gate-b-direct-costs.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-other-r17-gate-b-mitm.json"
)
U2 = matrix(ZZ, ((0, 1), (1, 0)))
DEGREE_MAX = ZZ(4)
Q_VALUES = (ZZ(4), ZZ(6), ZZ(8), ZZ(12))
SCALE_VALUES = tuple(range(20, 181, 4))
CLOSE_COUNT = 64
REAL = RealField(200)
RATIONAL_DENOMINATOR = ZZ(10) ** 40


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
        component = []
        while todo:
            index = todo.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other]]
            for other in adjacent:
                unseen.remove(other)
                todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def rational_approximation(value):
    return QQ(ZZ((value * RATIONAL_DENOMINATOR).round())) / RATIONAL_DENOMINATOR


def root_rank_and_tops(frame):
    result = pari(frame).qfminim(2)
    root_rank = matrix(ZZ, result[2]).rank() if int(result[0]) else 0
    if root_rank == 0:
        return 0, ()
    cartan = frame[:root_rank, :root_rank]
    assert set(cartan.diagonal()) == {2}
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(root_rank)
        for column in range(root_rank)
        if row != column
    )
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    tops = []
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
        tops.append(max(candidates, key=lambda root: sum(root)))
    return root_rank, tuple(tops)


def negative_horizontal_walls(frame, tail, degree):
    walls = []
    tail_norm = ZZ(tail * frame * tail)
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * frame * tail.column()
        augmented = block_matrix(
            ZZ,
            [
                [degree**2 * frame, cross],
                [
                    cross.transpose(),
                    matrix(ZZ, [[m**2 * tail_norm + 1]]),
                ],
            ],
        )
        result = pari(augmented).qfminim(2 * degree**2 - 1)
        normalized = set()
        for candidate in matrix(ZZ, result[2]).transpose().rows():
            if abs(candidate[-1]) != 1:
                continue
            value = candidate if candidate[-1] == 1 else -candidate
            normalized.add(tuple(value))
        for value in normalized:
            section_tail = vector(ZZ, value[:-1])
            section_norm = ZZ(section_tail * frame * section_tail)
            if (section_norm - 2) % (2 * m):
                continue
            k = ZZ((section_norm - 2) // (2 * m))
            intersection = ZZ(
                (tail_norm // (2 * degree)) * m
                + degree * k
                - tail * frame * section_tail
            )
            if intersection < 0:
                walls.append((int(m), int(intersection), list(map(int, section_tail))))
    return walls


def nef_gate(frame, root_rank, tops, divisor):
    a, degree = map(ZZ, divisor[:2])
    tail = vector(ZZ, divisor[2:])
    if a - degree < 0:
        return False, {"reason": "negative_zero_pairing"}
    labels = tail * frame[:, :root_rank]
    affine = [ZZ(degree - top * labels) for top in tops]
    if min(tuple(labels) + tuple(affine) + (ZZ(0),)) < 0:
        return False, {"reason": "negative_component_pairing"}
    lattice = IntegralLattice(frame)
    center = vector(QQ, tail) / degree
    closest = vector(ZZ, next(lattice.enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = degree * (distance - 2) / 2
    if minimum_section < 0:
        return False, {"reason": "negative_closest_section"}
    walls = negative_horizontal_walls(frame, tail, degree)
    if walls:
        return False, {"reason": "negative_horizontal_wall", "walls": walls}
    return True, {
        "zero_pairing": int(a - degree),
        "component_pairings": list(map(int, labels)),
        "affine_component_pairings": list(map(int, affine)),
        "minimum_section_intersection": str(minimum_section),
        "negative_horizontal_walls": [],
    }


data = json.loads(INPUT.read_text())
assert data["status"] == "PASS_EXACT_OTHER_R17_GATE_B_DIRECT_COST_AUDIT"
pinned = matrix(
    ZZ,
    [
        [ZZ(value) for value in line.split()]
        for line in (
            ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
        ).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ],
)
g_pinned = block_diagonal_matrix(U2, -pinned)

targets = {
    name: vector(ZZ, values)
    for name, values in data["targets_in_pinned_ns"].items()
}
target_bases = {
    name: matrix(ZZ, values)
    for name, values in data["target_graph_bases_in_pinned_ns"].items()
}
for name in targets:
    assert target_bases[name].row(0) == targets[name]
    assert abs(target_bases[name].det()) == 1

candidate_records = []
node_summaries = []
target_geometry = {}
for target_name, target_basis in target_bases.items():
    target_gram = target_basis * g_pinned * target_basis.transpose()
    target_frame = -target_gram[2:, 2:]
    target_root_rank, target_tops = root_rank_and_tops(target_frame)
    target_geometry[target_name] = (
        target_gram,
        target_frame,
        target_root_rank,
        target_tops,
    )

for node_name, node in sorted(data["equation_nodes"].items()):
    basis = matrix(ZZ, node["basis_in_pinned_ns"])
    local_gram = basis * g_pinned * basis.transpose()
    assert local_gram[:2, :2] == U2
    frame = -local_gram[2:, 2:]
    assert frame.is_positive_definite()
    root_rank, tops = root_rank_and_tops(frame)
    lattice = IntegralLattice(frame)
    pinned_in_node = basis.inverse().change_ring(ZZ)

    before = len(candidate_records)
    close_vectors = 0
    exact_norm_vectors = set()
    retained_keys = set()
    two_sided_lattice_matches = 0
    for target_name, target_basis in target_bases.items():
        target_local = vector(ZZ, targets[target_name] * pinned_in_node)
        target_tail = vector(ZZ, target_local[2:])
        target_tail_norm = ZZ(target_tail * frame * target_tail)
        assert target_tail_norm == 2 * target_local[0] * target_local[1]
        if target_tail_norm <= 0:
            continue
        target_inverse = target_basis.inverse().change_ring(ZZ)
        target_gram, target_frame, target_root_rank, target_tops = (
            target_geometry[target_name]
        )
        for source_degree in range(1, int(DEGREE_MAX) + 1):
            b = ZZ(source_degree)
            for q in Q_VALUES:
                if q % b:
                    continue
                a = q // b
                optimum = (REAL(2 * q) / REAL(target_tail_norm)).sqrt()
                seen = set()
                for scale_percent in SCALE_VALUES:
                    scalar = rational_approximation(
                        optimum * REAL(scale_percent) / 100
                    )
                    center = scalar * target_tail
                    for close in islice(
                        lattice.enumerate_close_vectors(center), CLOSE_COUNT
                    ):
                        close_vectors += 1
                        tail = vector(ZZ, close)
                        key = tuple(tail)
                        if key in seen:
                            continue
                        seen.add(key)
                        if tail * frame * tail != 2 * q:
                            continue
                        exact_norm_vectors.add((target_name, int(q), int(b), key))
                        middle_local = vector(ZZ, [a, b] + list(tail))
                        if gcd(
                            [abs(ZZ(value)) for value in local_gram * middle_local]
                        ) != 1:
                            continue
                        middle_pinned = vector(ZZ, middle_local * basis)
                        middle_target = vector(
                            ZZ, middle_pinned * target_inverse
                        )
                        assert middle_target * target_gram * middle_target == 0
                        reverse_degree = ZZ(middle_target[1])
                        reverse_q = ZZ(middle_target[0] * middle_target[1])
                        if not (
                            1 <= reverse_degree <= DEGREE_MAX
                            and reverse_q in Q_VALUES
                        ):
                            continue
                        retained_key = (
                            target_name,
                            tuple(middle_pinned),
                            int(q),
                            int(reverse_q),
                        )
                        if retained_key in retained_keys:
                            continue
                        retained_keys.add(retained_key)
                        two_sided_lattice_matches += 1
                        source_nef, source_nef_data = nef_gate(
                            frame, root_rank, tops, middle_local
                        )
                        target_nef, target_nef_data = nef_gate(
                            target_frame,
                            target_root_rank,
                            target_tops,
                            middle_target,
                        )
                        candidate_records.append(
                            {
                                "equation_node": node_name,
                                "target": target_name,
                                "middle_in_source": list(map(int, middle_local)),
                                "middle_in_target": list(map(int, middle_target)),
                                "middle_in_pinned_ns": list(map(int, middle_pinned)),
                                "source_degree": int(source_degree),
                                "source_q": int(q),
                                "target_degree": int(reverse_degree),
                                "target_q": int(reverse_q),
                                "discovery_scale_percent": scale_percent,
                                "source_nef": source_nef,
                                "source_nef_data": source_nef_data,
                                "target_nef": target_nef,
                                "target_nef_data": target_nef_data,
                                "physically_nef_both_sides": (
                                    source_nef and target_nef
                                ),
                            }
                        )
    node_summaries.append(
        {
            "equation_node": node_name,
            "root_rank": int(root_rank),
            "close_vectors_with_repetitions": close_vectors,
            "distinct_exact_norm_presentations": len(exact_norm_vectors),
            "two_sided_lattice_matches_with_possible_scale_repetition": (
                two_sided_lattice_matches
            ),
            "retained_matches": len(candidate_records) - before,
        }
    )
    print(
        "OTHERR17MITM|node={}|close={}|exact_norm={}|lattice_matches={}|"
        "retained={}|status=PASS_NODE".format(
            node_name,
            close_vectors,
            len(exact_norm_vectors),
            two_sided_lattice_matches,
            len(candidate_records) - before,
        ),
        flush=True,
    )

physical_candidates = [
    record for record in candidate_records
    if record["physically_nef_both_sides"]
]
payload = {
    "schema": "elkies-k3.other-r17-gate-b-mitm-cvp.v2",
    "status": "PASS_BOUNDED_OTHER_R17_GATE_B_MITM_CVP",
    "bounds": {
        "source_old_fibre_degree_max": int(DEGREE_MAX),
        "target_old_fibre_degree_max": int(DEGREE_MAX),
        "source_presentation_q_values": list(map(int, Q_VALUES)),
        "target_presentation_q_values": list(map(int, Q_VALUES)),
        "scale_percent_values": list(SCALE_VALUES),
        "close_vectors_per_scale": CLOSE_COUNT,
        "rational_approximation_denominator": str(RATIONAL_DENOMINATOR),
    },
    "method": (
        "Target-directed closest-vector sampling around the exact Lagrange ray "
        "for each equation node, endpoint, degree, and q. Every discovered "
        "vector is checked exactly for norm, primitive isotropy, and reverse "
        "degree/q. Retained matches then pass exact simple/affine component, "
        "closest-section, and augmented-lattice negative-horizontal-wall gates "
        "on both sides. Discovery is bounded and is not shell-exhaustive."
    ),
    "equation_node_count": len(data["equation_nodes"]),
    "targets": sorted(targets),
    "node_summaries": node_summaries,
    "candidates": candidate_records,
    "candidate_count": len(candidate_records),
    "physically_nef_candidate_count": len(physical_candidates),
    "physically_nef_candidates": physical_candidates,
    "conclusion": (
        "The declared target-directed CVP probe found no physically nef "
        "two-edge bridge; this is a bounded route-search result, not an "
        "exhaustive nonexistence theorem. Resume the retained Q80 equation route."
        if not physical_candidates
        else "Physically nef two-edge matches exist and require exact RR scoring."
    ),
    "input": str(INPUT.relative_to(ROOT)),
    "input_sha256": digest(INPUT),
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/search_other_r17_gate_b_mitm.sage"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "OTHERR17MITM|nodes={}|targets={}|matches={}|output={}|status=PASS".format(
        len(node_summaries), len(targets), len(candidate_records), OUTPUT
    ),
    flush=True,
)
