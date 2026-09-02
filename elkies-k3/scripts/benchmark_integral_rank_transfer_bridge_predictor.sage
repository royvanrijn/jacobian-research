#!/usr/bin/env sage-python
"""Benchmark two rank-two bridge predictors on normalized corridor data.

For a proposed new hyperbolic plane U', every root in

    K = U_old^perp intersect U'^perp

survives in the new frame.  Hence rank(Phi(K)) and |Phi(K)| are exact lower
bounds for the child root rank and signed root count.  This script replays
the deterministic Weyl-quotient order used by four historical H3 suffix
searches and measures how many full child classifications remain after those
two bridge-core bounds are imposed.  It never enumerates roots of a candidate
child, so the benchmark is independent of the outcome calculation it is
intended to screen.

The first benchmark replays historical H3 first-hit streams.  The second fixes
each observed rootless terminal core and determinant, exhausts its admissible
positive even binary bridge classes and oriented cyclic graph labels, and
tests the exploratory rule "maximize the bridge minimum."  Root enumeration
is used only to label the outcomes, not by either predictor.

Both tests are retrospective.  In particular, the terminal cores and their
determinants come from successful edges, so the second test measures fixed-core
reglue enrichment rather than prospective q-neighbor construction.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from sage.all import Genus, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix, lcm, matrix, pari, vector, xgcd


ROOT = Path(__file__).resolve().parents[2]
U = matrix(ZZ, ((0, 1), (1, 0)))
H3 = ROOT / "artifacts/generated-results/elkies-k3-rank17-to-h3-reverse-transport.json"
BRIDGES = ROOT / "artifacts/generated-results/elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-integral-rank-transfer-bridge-predictor-benchmark-v1.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def signed_roots(gram):
    result = pari(gram).qfminim(2)
    positive = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2].sage()).columns()
    ]
    answer = positive + [-root for root in positive]
    assert len(answer) == int(result[0])
    return answer


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
    return tuple(sorted(result, key=lambda component: (len(component), component)))


def dominant_weights(cartan, component, bound):
    block = cartan.matrix_from_rows_and_columns(component, component)
    inverse = block.inverse()
    assert all(value >= 0 for value in inverse.list())
    weights = []

    def recurse(prefix, norm):
        index = len(prefix)
        if index == len(component):
            weights.append((tuple(prefix), norm))
            return
        value = 0
        while True:
            added = inverse[index, index] * value**2
            added += 2 * value * sum(
                inverse[index, previous] * prefix[previous]
                for previous in range(index)
            )
            new_norm = norm + added
            if new_norm > bound:
                break
            recurse(prefix + [value], new_norm)
            value += 1

    recurse([], QQ(0))
    return tuple(weights)


def bezout_mate(ns, fibre):
    current = ZZ(0)
    coefficients = [ZZ(0)] * ns.nrows()
    for index, pairing in enumerate(ns * fibre):
        if not pairing:
            continue
        common, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = common
    if abs(current) != 1:
        return None
    if current == -1:
        coefficients = [-value for value in coefficients]
    mate = vector(ZZ, coefficients)
    mate -= ZZ(mate * ns * mate) // 2 * fibre
    assert fibre * ns * mate == 1 and mate * ns * mate == 0
    return mate


def dominant_orbits(frame, root_rank, q, degree):
    """Reproduce search_root_adapted_weyl_neighbors.sage's sorted orbit order."""
    target = ZZ(2 * q)
    cartan = frame[:root_rank, :root_rank]
    coupling = frame[:root_rank, root_rank:]
    tail = frame[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    height_scale = lcm(entry.denominator() for entry in height.list())
    scaled_height = (height_scale * height).change_ring(ZZ)
    mw_result = pari(scaled_height).qfminim(height_scale * target)
    mw_vector_map = {}
    for column in matrix(ZZ, mw_result[2].sage()).columns():
        for sign in (1, -1):
            value = sign * vector(ZZ, column)
            if value == 0 or value * height * value > target:
                continue
            canonical = min(tuple(value), tuple(-value))
            mw_vector_map[canonical] = vector(ZZ, canonical)
    mw_vectors = tuple(
        sorted(
            mw_vector_map.values(),
            key=lambda value: (value * height * value, tuple(value)),
        )
    )

    components = connected_components(cartan)
    component_weights = tuple(
        dominant_weights(cartan, component, QQ(target))
        for component in components
    )
    combined_by_norm = {}

    def combine(index, choices, norm):
        if index == len(component_weights):
            combined_by_norm.setdefault(norm, []).append(tuple(choices))
            return
        for values, weight_norm in component_weights[index]:
            new_norm = norm + weight_norm
            if new_norm <= target:
                combine(index + 1, choices + [(values, weight_norm)], new_norm)

    combine(0, [], QQ(0))
    cartan_inverse = cartan.inverse()
    orbit_map = {}
    for mw in mw_vectors:
        mw_norm = mw * height * mw
        for choices in combined_by_norm.get(target - mw_norm, ()):
            labels = vector(ZZ, [0] * root_rank)
            for component, (values, _) in zip(components, choices):
                for index, value in zip(component, values):
                    labels[index] = value
            root_coordinates = cartan_inverse * (labels - coupling * mw)
            if not all(value in ZZ for value in root_coordinates):
                continue
            witness = vector(ZZ, list(root_coordinates) + list(mw))
            assert witness * frame * witness == target
            orbit_map[tuple(witness)] = {
                "witness": witness,
                "mw_coordinates": mw,
                "dominant_labels": labels,
            }
    return (
        [orbit_map[key] for key in sorted(orbit_map)],
        {
            "mw_projection_representatives": len(mw_vectors),
            "pari_signed_mw_vector_count": int(mw_result[0]),
            "dominant_orbit_count": len(orbit_map),
            "height_gram": rational_rows(height),
        },
    )


def dominant_stream_window(
    frame, root_rank, q, degree, mw_vector_cap, first_index, last_index
):
    """Reproduce the capped deterministic stream without materializing it."""
    target = ZZ(2 * q)
    cartan = frame[:root_rank, :root_rank]
    coupling = frame[:root_rank, root_rank:]
    tail = frame[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    height_scale = lcm(entry.denominator() for entry in height.list())
    scaled_height = (height_scale * height).change_ring(ZZ)
    mw_result = pari(scaled_height).qfminim(
        height_scale * target, mw_vector_cap
    )
    mw_vector_map = {}
    for column in matrix(ZZ, mw_result[2].sage()).columns():
        for sign in (1, -1):
            value = sign * vector(ZZ, column)
            if value == 0 or value * height * value > target:
                continue
            canonical = min(tuple(value), tuple(-value))
            mw_vector_map[canonical] = vector(ZZ, canonical)
    mw_vectors = tuple(
        sorted(
            mw_vector_map.values(),
            key=lambda value: (value * height * value, tuple(value)),
        )
    )
    components = connected_components(cartan)
    component_weights = tuple(
        dominant_weights(cartan, component, QQ(target))
        for component in components
    )
    combined_by_norm = {}

    def combine(index, choices, norm):
        if index == len(component_weights):
            combined_by_norm.setdefault(norm, []).append(tuple(choices))
            return
        for values, weight_norm in component_weights[index]:
            new_norm = norm + weight_norm
            if new_norm <= target:
                combine(index + 1, choices + [(values, weight_norm)], new_norm)

    combine(0, [], QQ(0))
    cartan_inverse = cartan.inverse()
    seen = set()
    selected = []
    for mw in mw_vectors:
        mw_norm = mw * height * mw
        for choices in combined_by_norm.get(target - mw_norm, ()):
            labels = vector(ZZ, [0] * root_rank)
            for component, (values, _) in zip(components, choices):
                for index, value in zip(component, values):
                    labels[index] = value
            root_coordinates = cartan_inverse * (labels - coupling * mw)
            if not all(value in ZZ for value in root_coordinates):
                continue
            witness = vector(ZZ, list(root_coordinates) + list(mw))
            key = tuple(witness)
            if key in seen:
                continue
            seen.add(key)
            index = len(seen)
            if index >= first_index:
                selected.append(
                    {
                        "global_stream_index": index,
                        "witness": witness,
                        "mw_coordinates": mw,
                        "dominant_labels": labels,
                    }
                )
            if index == last_index:
                return (
                    selected,
                    {
                        "mw_projection_representatives": len(mw_vectors),
                        "pari_signed_mw_vector_count": int(mw_result[0]),
                        "mw_vector_cap": mw_vector_cap,
                        "stream_first_global_index": first_index,
                        "stream_last_global_index": last_index,
                        "height_gram": rational_rows(height),
                    },
                )
    raise ArithmeticError(f"stream ended at {len(seen)} before {last_index}")


def core_root_lower_bound(frame, roots, q, degree, witness):
    ns = block_diagonal_matrix(U, -frame)
    fibre = vector(ZZ, [q // degree, degree] + list(witness))
    assert fibre * ns * fibre == 0
    mate = bezout_mate(ns, fibre)
    if mate is None:
        return None
    survivors = []
    for root in roots:
        ambient = vector(ZZ, [0, 0] + list(root))
        if ambient * ns * fibre == 0 and ambient * ns * mate == 0:
            survivors.append(root)
    survivor_rank = (
        0
        if not survivors
        else matrix(ZZ, [list(root) for root in survivors]).rank()
    )
    return {
        "root_rank_lower_bound": int(survivor_rank),
        "signed_root_count_lower_bound": len(survivors),
    }


def bridge_split_lower_bound(frame, roots, q, degree, witness):
    """Return the exact root lower bound from the orthogonal split K+C."""
    ns = block_diagonal_matrix(U, -frame)
    fibre = vector(ZZ, [q // degree, degree] + list(witness))
    mate = bezout_mate(ns, fibre)
    if mate is None:
        return None
    child_basis = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(child_basis * ns * child_basis.transpose())
    assert child.is_positive_definite() and child.det() == frame.det()
    old_basis = identity_matrix(ZZ, 19)[2:, :]
    core_basis = old_basis.row_module(ZZ).intersection(
        child_basis.row_module(ZZ)
    ).basis_matrix()
    assert core_basis.nrows() == 15
    core_coordinates = child_basis.solve_left(core_basis).change_ring(ZZ)
    bridge_coordinates = (core_coordinates * child).right_kernel_matrix()
    assert bridge_coordinates.nrows() == 2
    bridge_gram = bridge_coordinates * child * bridge_coordinates.transpose()
    bridge_roots = signed_roots(bridge_gram)
    bridge_root_rank = (
        0
        if not bridge_roots
        else matrix(ZZ, [list(root) for root in bridge_roots]).rank()
    )
    core = core_root_lower_bound(frame, roots, q, degree, witness)
    assert core is not None
    split_coordinates = core_coordinates.stack(bridge_coordinates)
    split_index = abs(int(split_coordinates.det()))
    lll_basis = matrix(ZZ, pari(bridge_gram).qflllgram()).transpose()
    reduced_bridge_gram = lll_basis * bridge_gram * lll_basis.transpose()
    return {
        "core_root_rank": core["root_rank_lower_bound"],
        "core_signed_root_count": core["signed_root_count_lower_bound"],
        "bridge_root_rank": int(bridge_root_rank),
        "bridge_signed_root_count": len(bridge_roots),
        "split_root_rank_lower_bound": int(
            core["root_rank_lower_bound"] + bridge_root_rank
        ),
        "split_signed_root_count_lower_bound": int(
            core["signed_root_count_lower_bound"] + len(bridge_roots)
        ),
        "cyclic_glue_order_if_maximal": split_index,
        "bridge_determinant_absolute": abs(int(bridge_gram.det())),
        "reduced_bridge_gram": rows(reduced_bridge_gram),
    }


def discriminant_form_key(gram):
    normal = Genus(gram).discriminant_form().normal_form()
    return (
        tuple(map(int, normal.invariants())),
        tuple(map(str, normal.gram_matrix_quadratic().list())),
        str(normal.value_module_qf()),
    )


def reduced_even_binary_forms(determinant):
    """Enumerate Minkowski-reduced positive even binary forms of odd determinant."""
    result = []
    for half_left in range(1, determinant + 1):
        for half_right in range(half_left, determinant + 1):
            for off_diagonal in range(0, half_left + 1):
                if (
                    4 * half_left * half_right - off_diagonal**2
                    == determinant
                ):
                    result.append(
                        matrix(
                            ZZ,
                            [
                                [2 * half_left, off_diagonal],
                                [off_diagonal, 2 * half_right],
                            ],
                        )
                    )
    return result


def shortest_norm(gram):
    result = pari(gram).qfminim(2 * abs(int(gram.det())))
    vectors = matrix(ZZ, result[2].sage()).columns()
    return min(int(vector(ZZ, value) * gram * vector(ZZ, value)) for value in vectors)


def glued_frame(core, bridge, glue_vector):
    split_gram = block_diagonal_matrix(core, bridge)
    denominator = lcm(value.denominator() for value in glue_vector)
    generators = (
        denominator * identity_matrix(QQ, split_gram.nrows())
    ).stack(matrix(QQ, [denominator * glue_vector])).change_ring(ZZ)
    basis = (
        generators.row_module(ZZ).basis_matrix().change_ring(QQ)
        / denominator
    )
    candidate = basis * split_gram * basis.transpose()
    if not all(value in ZZ for value in candidate.list()):
        return None
    candidate = candidate.change_ring(ZZ)
    if any(value % 2 for value in candidate.diagonal()):
        return None
    return candidate


def terminal_binary_bridge_census(bridges):
    records = []
    for edge in bridges["edges"]:
        if int(edge["target_root_rank"]) != 0:
            continue
        order = ZZ(edge["bridge_replacement"]["common_cyclic_glue_order"])
        assert order.is_prime()
        core = matrix(ZZ, edge["core"]["gram"])
        assert not signed_roots(core)
        selected_bridge = matrix(ZZ, edge["new_frame"]["bridge_gram"])
        generators = edge["new_frame"]["glue_generators"]
        assert len(generators) == 1
        selected_glue = vector(
            QQ,
            [QQ(value) for value in generators[0]["K_plus_C_dual_coordinates"]],
        )
        core_glue = selected_glue[:-2]
        selected_frame = glued_frame(core, selected_bridge, selected_glue)
        assert selected_frame is not None
        target_form = discriminant_form_key(selected_frame)

        classes = []
        for binary in reduced_even_binary_forms(int(order)):
            bridge_generator = vector(ZZ, [1, 0]) * binary.inverse()
            if bridge_generator in ZZ**2:
                bridge_generator = vector(ZZ, [0, 1]) * binary.inverse()
            assert order * bridge_generator in ZZ**2
            oriented = []
            for multiplier in range(1, int(order)):
                if gcd(multiplier, order) != 1:
                    continue
                glue = vector(
                    QQ,
                    list(multiplier * core_glue) + list(bridge_generator),
                )
                candidate = glued_frame(core, binary, glue)
                if candidate is None or abs(candidate.det()) != abs(selected_frame.det()):
                    continue
                if discriminant_form_key(candidate) != target_form:
                    continue
                oriented.append(
                    {
                        "core_glue_multiplier": multiplier,
                        "signed_root_count": int(pari(candidate).qfminim(2)[0]),
                    }
                )
            assert len(oriented) == 2
            assert oriented[0]["signed_root_count"] == oriented[1]["signed_root_count"]
            classes.append(
                {
                    "bridge_gram": rows(binary),
                    "bridge_minimum": shortest_norm(binary),
                    "admissible_oriented_graph_labels": oriented,
                    "signed_root_count": oriented[0]["signed_root_count"],
                    "rootless": oriented[0]["signed_root_count"] == 0,
                }
            )
        maximum_minimum = max(row["bridge_minimum"] for row in classes)
        retained = [row for row in classes if row["bridge_minimum"] == maximum_minimum]
        records.append(
            {
                "corridor": edge["corridor"],
                "edge_index": int(edge["edge_index"]),
                "cyclic_glue_order": int(order),
                "binary_bridge_class_count": len(classes),
                "rootless_binary_bridge_classes": sum(row["rootless"] for row in classes),
                "maximum_bridge_minimum": maximum_minimum,
                "maximum_minimum_retained_classes": len(retained),
                "maximum_minimum_rootless_classes": sum(row["rootless"] for row in retained),
                "classes": classes,
            }
        )
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    h3 = json.loads(H3.read_text())
    bridges = json.loads(BRIDGES.read_text())
    assert h3["status"] == "PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT"
    assert bridges["status"] == "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES"
    assert bridges["aggregate"]["edge_count"] == 42
    bridge_by_edge = {
        (row["corridor"], int(row["edge_index"])): row
        for row in bridges["edges"]
    }

    # These are the four deterministic --stop-after-first-growth searches
    # whose complete MW quotient shells were used (no vector cap or stream
    # skip).  Edge 8 was a selected exhaustive-shell orbit rather than a first
    # hit, and edge 13 used a separately capped/streamed shell.
    edge_indices = (9, 10, 11, 12)
    records = []
    total_raw = 0
    total_rank_pass = 0
    total_rank_count_pass = 0
    total_split_rank_pass = 0
    total_split_rank_count_pass = 0
    for edge_index in edge_indices:
        old_stage = h3["stages"][edge_index - 1]
        new_stage = h3["stages"][edge_index]
        incoming = new_stage["incoming_neighbor"]
        frame = matrix(ZZ, old_stage["positive_frame"])
        root_rank = 17 - int(old_stage["mw_rank"])
        roots = signed_roots(frame)
        assert matrix(ZZ, frame[:root_rank, :root_rank]) == frame[:root_rank, :root_rank]
        assert matrix(ZZ, [list(root) for root in roots]).rank() == root_rank

        q = int(incoming["q"])
        degree = int(incoming["old_fiber_degree"])
        selected_orbit = int(incoming["orbit_index"])
        selected_witness = tuple(map(int, incoming["witness_in_parent_frame"]))
        orbit_rows, enumeration = dominant_orbits(frame, root_rank, q, degree)
        assert selected_orbit <= len(orbit_rows)
        assert tuple(map(int, orbit_rows[selected_orbit - 1]["witness"])) == selected_witness

        bridge = bridge_by_edge[("H3", edge_index)]
        target_root_rank = int(bridge["target_root_rank"])
        target_root_count = int(bridge["new_frame"]["root_count_signed"])
        screened = []
        for orbit_index, orbit in enumerate(orbit_rows[:selected_orbit], 1):
            lower = bridge_split_lower_bound(
                frame, roots, q, degree, orbit["witness"]
            )
            if lower is None:
                lower = {
                    "core_root_rank": root_rank + 1,
                    "core_signed_root_count": len(roots) + 1,
                    "bridge_root_rank": 0,
                    "bridge_signed_root_count": 0,
                    "split_root_rank_lower_bound": root_rank + 1,
                    "split_signed_root_count_lower_bound": len(roots) + 1,
                }
            core_rank_pass = lower["core_root_rank"] <= target_root_rank
            core_rank_count_pass = (
                core_rank_pass
                and lower["core_signed_root_count"] <= target_root_count
            )
            split_rank_pass = lower["split_root_rank_lower_bound"] <= target_root_rank
            split_rank_count_pass = (
                split_rank_pass
                and lower["split_signed_root_count_lower_bound"] <= target_root_count
            )
            screened.append(
                (
                    orbit_index,
                    lower,
                    core_rank_pass,
                    core_rank_count_pass,
                    split_rank_pass,
                    split_rank_count_pass,
                )
            )

        selected = screened[-1]
        assert selected[0] == selected_orbit
        assert all(selected[2:])
        rank_pass_count = sum(row[2] for row in screened)
        rank_count_pass_count = sum(row[3] for row in screened)
        split_rank_pass_count = sum(row[4] for row in screened)
        split_rank_count_pass_count = sum(row[5] for row in screened)
        total_raw += selected_orbit
        total_rank_pass += rank_pass_count
        total_rank_count_pass += rank_count_pass_count
        total_split_rank_pass += split_rank_pass_count
        total_split_rank_count_pass += split_rank_count_pass_count
        lower_histogram = Counter(
            (
                row[1]["core_root_rank"],
                row[1]["core_signed_root_count"],
                row[1]["bridge_root_rank"],
                row[1]["bridge_signed_root_count"],
            )
            for row in screened
        )
        records.append(
            {
                "corridor": "H3",
                "edge_index": edge_index,
                "source": old_stage["stage_id"],
                "target": new_stage["stage_id"],
                "q": q,
                "old_fibre_degree": degree,
                "source_root_rank": root_rank,
                "source_signed_root_count": len(roots),
                "requested_target_root_rank_at_most": target_root_rank,
                "requested_target_signed_root_count_at_most": target_root_count,
                "historical_first_hit_orbit": selected_orbit,
                "raw_full_child_classifications_through_first_hit": selected_orbit,
                "bridge_core_rank_bound_passes_through_first_hit": rank_pass_count,
                "bridge_core_rank_and_count_bound_passes_through_first_hit": rank_count_pass_count,
                "bridge_split_rank_bound_passes_through_first_hit": split_rank_pass_count,
                "bridge_split_rank_and_count_bound_passes_through_first_hit": split_rank_count_pass_count,
                "classification_reduction_fraction": str(
                    QQ(selected_orbit - split_rank_count_pass_count) / selected_orbit
                ),
                "classification_speedup_factor": str(
                    QQ(selected_orbit) / split_rank_count_pass_count
                ),
                "selected_candidate_bridge_split": selected[1],
                "lower_bound_histogram_through_first_hit": [
                    {
                        "core_root_rank": core_rank,
                        "core_signed_root_count": core_count,
                        "bridge_root_rank": bridge_rank,
                        "bridge_signed_root_count": bridge_count,
                        "candidate_count": multiplicity,
                    }
                    for (
                        core_rank,
                        core_count,
                        bridge_rank,
                        bridge_count,
                    ), multiplicity in sorted(lower_histogram.items())
                ],
                "enumeration": enumeration,
            }
        )

    # The terminal q6 search used a capped MW-vector sample and streamed only
    # global dominant indices 1001 onward.  Its pinned first hit is index 2247,
    # so this separately benchmarks the exactly tested 1,247-candidate window.
    edge_index = 13
    old_stage = h3["stages"][edge_index - 1]
    new_stage = h3["stages"][edge_index]
    incoming = new_stage["incoming_neighbor"]
    frame = matrix(ZZ, old_stage["positive_frame"])
    root_rank = 17 - int(old_stage["mw_rank"])
    roots = signed_roots(frame)
    selected_global_index = int(incoming["orbit_index"])
    stream_rows, enumeration = dominant_stream_window(
        frame,
        root_rank,
        int(incoming["q"]),
        int(incoming["old_fiber_degree"]),
        10_000,
        1_001,
        selected_global_index,
    )
    assert tuple(map(int, stream_rows[-1]["witness"])) == tuple(
        map(int, incoming["witness_in_parent_frame"])
    )
    bridge = bridge_by_edge[("H3", edge_index)]
    assert int(bridge["target_root_rank"]) == 0
    assert int(bridge["new_frame"]["root_count_signed"]) == 0
    screened = []
    for orbit in stream_rows:
        lower = bridge_split_lower_bound(
            frame,
            roots,
            int(incoming["q"]),
            int(incoming["old_fiber_degree"]),
            orbit["witness"],
        )
        if lower is None:
            lower = {
                "core_root_rank": root_rank + 1,
                "core_signed_root_count": len(roots) + 1,
                "bridge_root_rank": 0,
                "bridge_signed_root_count": 0,
                "split_root_rank_lower_bound": root_rank + 1,
                "split_signed_root_count_lower_bound": len(roots) + 1,
            }
        core_pass = (
            lower["core_root_rank"] == 0
            and lower["core_signed_root_count"] == 0
        )
        split_pass = (
            lower["split_root_rank_lower_bound"] == 0
            and lower["split_signed_root_count_lower_bound"] == 0
        )
        screened.append((orbit["global_stream_index"], lower, core_pass, split_pass))
    assert screened[-1][0] == selected_global_index and screened[-1][2] and screened[-1][3]
    raw_count = len(screened)
    core_pass_count = sum(row[2] for row in screened)
    split_pass_count = sum(row[3] for row in screened)
    total_raw += raw_count
    total_rank_pass += core_pass_count
    total_rank_count_pass += core_pass_count
    total_split_rank_pass += split_pass_count
    total_split_rank_count_pass += split_pass_count
    lower_histogram = Counter(
        (
            row[1]["core_root_rank"],
            row[1]["core_signed_root_count"],
            row[1]["bridge_root_rank"],
            row[1]["bridge_signed_root_count"],
        )
        for row in screened
    )
    records.append(
        {
            "corridor": "H3",
            "edge_index": edge_index,
            "source": old_stage["stage_id"],
            "target": new_stage["stage_id"],
            "q": int(incoming["q"]),
            "old_fibre_degree": int(incoming["old_fiber_degree"]),
            "source_root_rank": root_rank,
            "source_signed_root_count": len(roots),
            "requested_target_root_rank_at_most": 0,
            "requested_target_signed_root_count_at_most": 0,
            "historical_first_hit_global_stream_index": selected_global_index,
            "tested_stream_window": [1_001, selected_global_index],
            "raw_full_child_classifications_through_first_hit": raw_count,
            "bridge_core_rank_bound_passes_through_first_hit": core_pass_count,
            "bridge_core_rank_and_count_bound_passes_through_first_hit": core_pass_count,
            "bridge_split_rank_bound_passes_through_first_hit": split_pass_count,
            "bridge_split_rank_and_count_bound_passes_through_first_hit": split_pass_count,
            "classification_reduction_fraction": str(QQ(raw_count - split_pass_count) / raw_count),
            "classification_speedup_factor": str(QQ(raw_count) / split_pass_count),
            "selected_candidate_bridge_split": screened[-1][1],
            "lower_bound_histogram_through_first_hit": [
                {
                    "core_root_rank": core_rank,
                    "core_signed_root_count": core_count,
                    "bridge_root_rank": bridge_rank,
                    "bridge_signed_root_count": bridge_count,
                    "candidate_count": multiplicity,
                }
                for (
                    core_rank,
                    core_count,
                    bridge_rank,
                    bridge_count,
                ), multiplicity in sorted(lower_histogram.items())
            ],
            "enumeration": enumeration,
        }
    )

    terminal_records = terminal_binary_bridge_census(bridges)
    terminal_class_count = sum(
        row["binary_bridge_class_count"] for row in terminal_records
    )
    terminal_rootless_count = sum(
        row["rootless_binary_bridge_classes"] for row in terminal_records
    )
    terminal_retained_count = sum(
        row["maximum_minimum_retained_classes"] for row in terminal_records
    )
    terminal_retained_rootless_count = sum(
        row["maximum_minimum_rootless_classes"] for row in terminal_records
    )

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-bridge-predictor-benchmark.v1",
        "status": "PASS_EXACT_RETROSPECTIVE_BRIDGE_PREDICTOR_BENCHMARKS",
        "predictor": {
            "name": "mandatory orthogonal bridge-split root budget",
            "input": "the old marked frame, proposed U-prime, and requested child root-rank/root-count budget",
            "rule": (
                "Reject when the root rank or signed root count of the orthogonal "
                "subsystem Phi(K)+Phi(C_new) exceeds the requested child budget."
            ),
            "safety": (
                "Exact with no false negatives for the declared budget because "
                "K+C_new is a sublattice of every child frame completing this U-prime."
            ),
            "excluded_information": (
                "No candidate child frame or candidate child roots are constructed; "
                "new roots introduced by bridge/glue cosets are deliberately unknown."
            ),
        },
        "inputs": {
            str(H3.relative_to(ROOT)): digest(H3),
            str(BRIDGES.relative_to(ROOT)): digest(BRIDGES),
        },
        "benchmark_design": {
            "edges": list(edge_indices) + [13],
            "selection": (
                "Four consecutive historical H3 suffix searches with complete "
                "MW quotient shells and deterministic first-hit stopping, plus "
                "the pinned capped/streamed terminal q6 first-hit window."
            ),
            "negative_labels": (
                "Every earlier orbit was rejected by the historical exact full-child "
                "classification; the selected orbit was the first target hit."
            ),
            "limitations": (
                "Retrospective, one corridor, five edges, and target budgets taken from "
                "the observed route. It tests screening efficiency, not out-of-sample "
                "route discovery or completeness."
            ),
            "q80_exclusion": (
                "The retained Q80 score tables contain only already rank-growing "
                "children, so they do not supply an unbiased negative candidate set."
            ),
        },
        "edges": records,
        "aggregate": {
            "edge_count": len(records),
            "raw_full_child_classifications_through_first_hits": total_raw,
            "bridge_core_rank_bound_passes_through_first_hits": total_rank_pass,
            "bridge_core_rank_and_count_bound_passes_through_first_hits": total_rank_count_pass,
            "bridge_split_rank_bound_passes_through_first_hits": total_split_rank_pass,
            "bridge_split_rank_and_count_bound_passes_through_first_hits": total_split_rank_count_pass,
            "classification_reduction_fraction": str(
                QQ(total_raw - total_split_rank_count_pass) / total_raw
            ),
            "classification_speedup_factor": str(
                QQ(total_raw) / total_split_rank_count_pass
            ),
            "selected_hits_preserved": len(records),
        },
        "terminal_binary_bridge_census": {
            "design": (
                "For each of the four observed rootless terminal cores, enumerate "
                "all Minkowski-reduced positive even binary bridge classes of the "
                "observed prime determinant and every compatible oriented cyclic "
                "graph label in the frame discriminant genus."
            ),
            "predictor": (
                "Retain the binary bridge class or classes having maximum lattice minimum."
            ),
            "selection_warning": (
                "The core and bridge determinant come from the observed successful "
                "edge. This is a complete fixed-core reglue census but is not an "
                "out-of-sample method for choosing a core from the old frame."
            ),
            "aggregate": {
                "corridor_count": len(terminal_records),
                "binary_bridge_class_count": terminal_class_count,
                "rootless_binary_bridge_classes": terminal_rootless_count,
                "unfiltered_rootless_precision": str(
                    QQ(terminal_rootless_count) / terminal_class_count
                ),
                "maximum_minimum_retained_classes": terminal_retained_count,
                "maximum_minimum_rootless_classes": terminal_retained_rootless_count,
            "maximum_minimum_rootless_precision": str(
                QQ(terminal_retained_rootless_count) / terminal_retained_count
            ),
            "maximum_minimum_rootless_recall": str(
                QQ(terminal_retained_rootless_count) / terminal_rootless_count
            ),
            "candidate_reduction_fraction": str(
                1 - QQ(terminal_retained_count) / terminal_class_count
            ),
            "projected_classification_speedup_factor": str(
                QQ(terminal_class_count) / terminal_retained_count
            ),
            "precision_enrichment_factor": str(
                (QQ(terminal_retained_rootless_count) / terminal_retained_count)
                / (QQ(terminal_rootless_count) / terminal_class_count)
                ),
            },
            "corridors": terminal_records,
        },
        "conclusion": {
            "result": "mixed_retrospective_evidence",
            "observed_answer": (
                "The bridge core alone rejects no candidate. Adding the roots of "
                "the rank-two bridge rejects only 178 of 2892 classifications; "
                "this is not a substantial improvement over raw enumeration."
            ),
            "fixed_core_reglue_answer": (
                "Across the four complete terminal fixed-core binary-bridge "
                "censuses, 5 of 14 bridge classes are rootless. The exploratory "
                "maximum-minimum rule retains 5 classes, 4 rootless. This is "
                "strong enrichment inside a selected successful core, but not a "
                "prospective q-neighbor construction test."
            ),
            "new_construction_algorithm_gate": (
                "Not passed. A glue-coset-minimum score must next be predeclared and "
                "tested on an untouched prospective shell, with its runtime compared "
                "against direct child-root enumeration."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "benchmark_integral_rank_transfer_bridge_predictor.sage --check"
        ),
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS integral rank-transfer bridge predictor benchmark")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
    try:
        print(output.relative_to(ROOT))
    except ValueError:
        print(output)


if __name__ == "__main__":
    main()
