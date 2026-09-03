#!/usr/bin/env sage-python
"""Search bounded relative-U markings for the NS0024 completed-frame chain.

The input frames are the four exact positive-definite completions certified by
``certify_ns0024_new_rootless_source_route.sage``.  For each consecutive target
this script enumerates the declared intersection box in lexicographic order,
represents ``G_A=A^t*J*A-J`` in the current positive frame, constructs the
literal primitive target U, and tests its orthogonal frame against the declared
target by roots and exact integral isometry.

This is a bounded lattice search.  It does not certify nefness, an effective
irreducible zero, equations, or rational maps.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import runpy

import numpy as np

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
OUTPUT = GENERATED / "elkies-k3-ns0024-relative-u-bridge-lifts-v1.json"
CERT_SCRIPT = ROOT / "elkies-k3/scripts/certify_ns0024_new_rootless_source_route.sage"
KNOWN_ROUTE = GENERATED / "elkies-k3-lattice-foundry-ns0024-r13-nef-route.json"
KNOWN_SOURCE = GENERATED / "elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
J = matrix(ZZ, [[0, 1], [1, 0]])


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def relative(path):
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text())


def completed_ns0024_frames():
    """Rebuild the four pinned completions without copying their Gram data."""
    certificate = runpy.run_path(str(CERT_SCRIPT))
    base = runpy.run_path(str(certificate["BASE_SCRIPT"]))
    search = runpy.run_path(str(certificate["SEARCH_SCRIPT"]))
    core = runpy.run_path(str(certificate["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(certificate["REVERSE_SCRIPT"]))
    bridges = load_json(certificate["BRIDGES"])
    theta = load_json(certificate["THETA"])
    masked = load_json(certificate["MASKED"])

    prepared = search["prepare_corridor"](
        "NS0024", bridges, theta, base, core, reverse
    )
    search["configure_order"](base, prepared["order"])
    masked_row = next(
        row for row in masked["corridors"] if row["corridor"] == "NS0024"
    )
    bridge_index = masked_row["completion"]["bridge_class_index"]
    terminal_multiplier = masked_row["completion"]["glue_multiplier"]
    bridge = next(
        row
        for row in prepared["viable_bridges"]
        if row["bridge_class_index"] == bridge_index
    )

    quadratic = base["quadratic_form"](prepared["seed"])
    cores = [base["lll_reduce"](quadratic.Hessian_matrix())]
    for prime, raw_witness in certificate["PATH"]:
        witness = vector(ZZ, raw_witness)
        quadratic = quadratic.find_p_neighbor_from_vec(prime, witness)
        cores.append(base["lll_reduce"](quadratic.Hessian_matrix()))

    expected = ((13, 280, 4), (5, 12, 24), (5, 12, 24), (0, 0, 1))
    frames = []
    multipliers = []
    for index, core_gram in enumerate(cores):
        choices = certificate["completed_frames"](
            core_gram, bridge, prepared["order"], base, core
        )
        matching = [
            row
            for row in choices
            if certificate["root_data"](row[1]) == expected[index]
        ]
        assert matching
        if index == len(cores) - 1:
            selected = next(row for row in matching if row[0] == terminal_multiplier)
        else:
            selected = min(matching, key=lambda row: row[0])
        multipliers.append(int(selected[0]))
        frames.append(selected[1])
    return frames, multipliers, certificate, expected


def known_ns0024_route_frames():
    source = load_json(KNOWN_SOURCE)
    route = load_json(KNOWN_ROUTE)
    current = matrix(ZZ, source["source"]["root_adapted_gram"])
    current_basis = identity_matrix(ZZ, 19)
    frames = [current]
    bases = [current_basis]
    for edge in route["edges"]:
        transport = matrix(ZZ, edge["edge_transport_child_to_parent"])
        ns = block_diagonal_matrix(J, -current)
        child_ns = transport * ns * transport.transpose()
        assert child_ns[:2, :2] == J
        assert child_ns[:2, 2:] == 0 and child_ns[2:, :2] == 0
        current = -child_ns[2:, 2:]
        current_basis = transport * current_basis
        frames.append(current)
        bases.append(current_basis)
    return frames, bases


def signed_vectors_of_norm(gram, norm, cache):
    key = (tuple(map(tuple, rows(gram))), int(norm))
    if key in cache:
        return cache[key]
    result = pari(gram).qfminim(int(norm))
    positive = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2].sage()).columns()
    ]
    answer = []
    for value in positive:
        for signed in (value, -value):
            if int(signed * gram * signed) == norm:
                answer.append(signed)
    cache[key] = answer
    return answer


def root_data(gram):
    result = pari(gram).qfminim(2)
    count = int(result[0])
    if count == 0:
        return 0, 0, 1
    roots = matrix(ZZ, result[2].sage()).transpose()
    root_basis = roots.row_module(ZZ).basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return int(root_basis.rank()), count, int(abs(root_gram.det()))


def deterministic_simple_roots(roots):
    positive = [
        root
        for root in roots
        if next(value for value in root if value != 0) > 0
    ]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root
        for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    ]
    return simple


def weyl_dominant(value, gram, simple_roots):
    result = vector(ZZ, value)
    while True:
        changed = False
        for root in simple_roots:
            pairing = int(result * gram * root)
            if pairing < 0:
                result -= pairing * root
                changed = True
                break
        if not changed:
            return result


def root_adaptation(gram):
    roots = signed_vectors_of_norm(gram, 2, {})
    simple = matrix(ZZ, [list(value) for value in deterministic_simple_roots(roots)])
    rank = simple.nrows()
    root_module = matrix(ZZ, [list(value) for value in roots]).row_module(ZZ)
    smith, _, smith_right = root_module.basis_matrix().smith_form()
    diagonal = [abs(int(smith[index, index])) for index in range(rank)]
    if diagonal != [1] * rank:
        return {"primitive_root_lattice": False, "smith_diagonal": diagonal}
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[rank:])
    adapted = adapted_basis * gram * adapted_basis.transpose()
    cartan = adapted[:rank, :rank]
    coupling = adapted[:rank, rank:]
    tail = adapted[rank:, rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    minimum = int(pari((scale * height).change_ring(ZZ)).qfminim()[1]) / scale
    return {
        "primitive_root_lattice": True,
        "root_rank": rank,
        "adapted_basis": rows(adapted_basis),
        "adapted_gram": rows(adapted),
        "mw_height_gram": [[str(entry) for entry in row] for row in height.rows()],
        "mw_height_minimum": str(minimum),
    }


def verified_qfisom(left, right):
    raw = pari(left).qfisom(pari(right))
    if raw == 0:
        return None
    candidate = matrix(ZZ, raw.sage())
    for value in (candidate, candidate.transpose()):
        if value.transpose() * right * value == left:
            return value
        if value * right * value.transpose() == left:
            return value.transpose()
    raise ArithmeticError("PARI qfisom returned an unrecognized orientation")


def exact_int64_pairings(left_rows, gram, right_rows):
    """Return ``left_rows * gram * right_rows^t`` after an overflow proof.

    The relative-U search spends most of its time rejecting representations by
    their surviving source roots.  The matrices in this search are small enough
    for exact signed 64-bit arithmetic, but make that a checked hypothesis rather
    than silently relying on NumPy overflow behavior.
    """
    left = np.asarray(rows(left_rows), dtype=np.int64)
    middle = np.asarray(rows(gram), dtype=np.int64)
    right = np.asarray(rows(right_rows), dtype=np.int64)
    if left.ndim == 1:
        left = left.reshape((1, -1))
    if right.ndim == 1:
        right = right.reshape((1, -1))
    factors = (
        left.shape[1] ** 2,
        int(np.max(np.abs(left), initial=0)),
        int(np.max(np.abs(middle), initial=0)),
        int(np.max(np.abs(right), initial=0)),
    )
    bound = 1
    for factor in factors:
        bound *= factor
    if bound >= 2**62:
        raise OverflowError(
            "relative-U pairing block does not have a certified int64 bound"
        )
    return left @ middle @ right.transpose()


def construct_relative_u(frame, A, w_1, w_2):
    ns = block_diagonal_matrix(J, -frame)
    w_rows = matrix(ZZ, [list(w_1), list(w_2)])
    u_part = A.transpose() * J
    u_prime = u_part.augment(w_rows)
    assert u_prime * ns * u_prime.transpose() == J
    complement = (u_prime * ns).right_kernel_matrix()
    transport = u_prime.stack(complement)
    assert abs(int(transport.det())) == 1
    child_ns = transport * ns * transport.transpose()
    assert child_ns[:2, :2] == J
    assert child_ns[:2, 2:] == 0 and child_ns[2:, :2] == 0
    child = -child_ns[2:, 2:]
    return child, transport, u_prime, complement


def search_edge(
    source,
    target,
    expected_root_data,
    degrees,
    s_values,
    t_values,
    z_values,
    shell_cache,
    pairing_chunk_size=20000,
):
    counters = Counter()
    degree_summaries = []
    source_roots = signed_vectors_of_norm(source, 2, shell_cache)
    source_root_matrix = matrix(ZZ, [list(value) for value in source_roots])
    simple_roots = deterministic_simple_roots(source_roots)
    for degree in degrees:
        degree_start = dict(counters)
        for s in s_values:
            for t in t_values:
                for z in z_values:
                    counters["matrices_tested"] += 1
                    A = matrix(
                        ZZ,
                        [
                            [degree, degree + s],
                            [degree + t, degree + s + t + z],
                        ],
                    )
                    G_A = A.transpose() * J * A - J
                    if not G_A.is_positive_definite():
                        counters["non_positive_grams"] += 1
                        continue
                    reduction = matrix(ZZ, pari(G_A).qflllgram()).transpose()
                    assert abs(int(reduction.det())) == 1
                    reduced_gram = reduction * G_A * reduction.transpose()
                    shell_1 = signed_vectors_of_norm(
                        source, int(reduced_gram[0, 0]), shell_cache
                    )
                    shell_2 = signed_vectors_of_norm(
                        source, int(reduced_gram[1, 1]), shell_cache
                    )
                    shell_1_orbits = {}
                    for value in shell_1:
                        dominant = weyl_dominant(value, source, simple_roots)
                        shell_1_orbits[tuple(dominant)] = dominant
                    shell_1 = [shell_1_orbits[key] for key in sorted(shell_1_orbits)]
                    counters["positive_grams"] += 1
                    counters["first_shell_vectors"] += len(shell_1)
                    counters["second_shell_vectors"] += len(shell_2)
                    inverse_reduction = reduction.inverse().change_ring(ZZ)
                    shell_2_matrix = matrix(ZZ, [list(value) for value in shell_2])
                    for reduced_w_1 in shell_1:
                        root_pairings_1 = source_root_matrix * source * reduced_w_1
                        roots_orthogonal_to_first = matrix(
                            ZZ,
                            [
                                list(source_root_matrix.row(index))
                                for index, pairing in enumerate(root_pairings_1)
                                if pairing == 0
                            ],
                        )
                        shell_2_pairings = exact_int64_pairings(
                            matrix(ZZ, [list(reduced_w_1)]),
                            source,
                            shell_2_matrix,
                        )[0]
                        matching_second_indices = [
                            index
                            for index, pairing in enumerate(shell_2_pairings)
                            if pairing == reduced_gram[0, 1]
                        ]
                        for chunk_start in range(
                            0, len(matching_second_indices), pairing_chunk_size
                        ):
                            chunk_indices = matching_second_indices[
                                chunk_start : chunk_start + pairing_chunk_size
                            ]
                            chunk = matrix(
                                ZZ,
                                [list(shell_2_matrix.row(index)) for index in chunk_indices],
                            )
                            root_pairing_block = exact_int64_pairings(
                                roots_orthogonal_to_first, source, chunk
                            )
                            zero_counts = np.count_nonzero(
                                root_pairing_block == 0, axis=0
                            )
                            for local_index, second_index in enumerate(chunk_indices):
                                reduced_w_2 = shell_2_matrix.row(second_index)
                                counters["representations"] += 1
                                core_root_count = int(zero_counts[local_index])
                                if core_root_count > expected_root_data[1]:
                                    counters["rejected_by_core_root_count"] += 1
                                    continue
                                zero_indices = np.flatnonzero(
                                    root_pairing_block[:, local_index] == 0
                                )
                                core_roots = matrix(
                                    ZZ,
                                    [
                                        list(roots_orthogonal_to_first.row(int(root_index)))
                                        for root_index in zero_indices
                                    ],
                                )
                                core_root_rank = (
                                    0 if not core_root_count else core_roots.rank()
                                )
                                if core_root_rank > expected_root_data[0]:
                                    counters["rejected_by_core_root_rank"] += 1
                                    continue
                                w_rows = inverse_reduction * matrix(
                                    ZZ, [list(reduced_w_1), list(reduced_w_2)]
                                )
                                w_1 = w_rows.row(0)
                                w_2 = w_rows.row(1)
                                assert w_rows * source * w_rows.transpose() == G_A
                                counters["core_root_bound_survivors"] += 1
                                child, transport, u_prime, complement = construct_relative_u(
                                    source, A, w_1, w_2
                                )
                                counters["primitive_u_markings"] += 1
                                if root_data(child) != expected_root_data:
                                    continue
                                counters["root_profile_matches"] += 1
                                isometry = verified_qfisom(child, target)
                                if isometry is None:
                                    continue
                                counters["integral_isometry_matches"] += 1
                                return {
                                    "status": "HIT",
                                    "degree": degree,
                                    "intersection_coordinates": {
                                        "F_dot_F_prime": degree,
                                        "F_dot_O_prime": s,
                                        "O_dot_F_prime": t,
                                        "O_dot_O_prime": z,
                                    },
                                    "cross_pairing_A": rows(A),
                                    "positive_projection_gram_G_A": rows(G_A),
                                    "binary_reduction": {
                                        "basis_change": rows(reduction),
                                        "reduced_gram": rows(reduced_gram),
                                    },
                                    "projected_vectors_in_source_frame": [
                                        list(map(int, w_1)),
                                        list(map(int, w_2)),
                                    ],
                                    "primitive_u_basis_in_source_ns": rows(u_prime),
                                    "child_frame_basis_in_source_ns": rows(complement),
                                    "child_frame_gram": rows(child),
                                    "target_isometry": rows(isometry),
                                    "counters_through_hit": dict(counters),
                                    "degree_summaries_before_hit": degree_summaries,
                                    "search_order": "degree,s,t,z,source-shell vectors",
                                }, transport
        degree_summaries.append(
            {
                "degree": degree,
                "complete_no_hit": True,
                "counter_delta": {
                    key: counters[key] - degree_start.get(key, 0)
                    for key in sorted(counters)
                },
            }
        )
    return {
        "status": "NO_HIT_IN_DECLARED_BOX",
        "counters": dict(counters),
        "degree_summaries": degree_summaries,
    }, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degrees", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument("--max-s", type=int, default=4)
    parser.add_argument("--max-t", type=int, default=4)
    parser.add_argument("--max-z", type=int, default=4)
    parser.add_argument("--s-values", type=int, nargs="+")
    parser.add_argument("--t-values", type=int, nargs="+")
    parser.add_argument("--z-values", type=int, nargs="+")
    parser.add_argument("--compare-known-only", action="store_true")
    parser.add_argument("--export-adapted-source", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    assert all(value >= 0 for value in (arguments.max_s, arguments.max_t, arguments.max_z))
    s_values = sorted(set(arguments.s_values or range(arguments.max_s + 1)))
    t_values = sorted(set(arguments.t_values or range(arguments.max_t + 1)))
    z_values = sorted(set(arguments.z_values or range(arguments.max_z + 1)))
    assert all(value >= 0 for value in s_values + t_values + z_values)

    frames, multipliers, certificate, expected = completed_ns0024_frames()
    known_frames, known_bases = known_ns0024_route_frames()
    known_route_matches = []
    for target_index, target in enumerate(frames):
        matches = []
        for known_index, known in enumerate(known_frames):
            if root_data(known) != expected[target_index]:
                continue
            isometry = verified_qfisom(target, known)
            if isometry is not None:
                matches.append(
                    {
                        "known_route_stage": known_index,
                        "target_to_known_frame_isometry": rows(isometry),
                        "known_u_basis_in_known_source_ns": rows(
                            known_bases[known_index][:2, :]
                        ),
                    }
                )
        known_route_matches.append(
            {"target_stage": target_index, "matches": matches}
        )

    if arguments.compare_known_only:
        source_adaptation = root_adaptation(frames[0])
        if arguments.export_adapted_source is not None:
            adapted_path = (
                arguments.export_adapted_source
                if arguments.export_adapted_source.is_absolute()
                else ROOT / arguments.export_adapted_source
            )
            adapted_path.parent.mkdir(parents=True, exist_ok=True)
            adapted_path.write_text(
                "\n".join(
                    " ".join(map(str, row))
                    for row in source_adaptation["adapted_gram"]
                )
                + "\n"
            )
        output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(
                {
                    "schema": "elkies-k3.ns0024-relative-u-known-route-comparison.v1",
                    "status": "PASS_EXACT_KNOWN_ROUTE_FRAME_COMPARISON",
                    "matches": known_route_matches,
                    "completed_source_root_adaptation": source_adaptation,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(relative(output))
        return
    current_frame = frames[0]
    current_basis_in_source_ns = identity_matrix(ZZ, 19)
    shell_cache = {}
    edges = []
    for index, target in enumerate(frames[1:], 1):
        result, local_transport = search_edge(
            current_frame,
            target,
            expected[index],
            arguments.degrees,
            s_values,
            t_values,
            z_values,
            shell_cache,
        )
        result["edge_index"] = index
        result["source_root_data"] = list(root_data(current_frame))
        result["target_root_data"] = list(expected[index])
        if local_transport is None:
            edges.append(result)
            break
        composed = local_transport * current_basis_in_source_ns
        result["primitive_u_basis_in_initial_source_ns"] = rows(composed[:2, :])
        result["child_frame_basis_in_initial_source_ns"] = rows(composed[2:, :])
        edges.append(result)
        current_frame = matrix(ZZ, result["child_frame_gram"])
        current_basis_in_source_ns = composed

    all_hit = len(edges) == 3 and all(row["status"] == "HIT" for row in edges)
    payload = {
        "schema": "elkies-k3.ns0024-relative-u-bridge-lifts.v1",
        "status": (
            "PASS_EXACT_THREE_EDGE_RELATIVE_U_LIFT"
            if all_hit
            else "PASS_EXACT_BOUNDED_RELATIVE_U_SEARCH_NO_COMPLETE_CHAIN"
        ),
        "source": {
            "ns": "NS0024",
            "frame_sequence": ["D5+E8", "3A1+A2", "3A1+A2", "W_new_950"],
            "completion_glue_multipliers": multipliers,
        },
        "known_route_frame_comparison": known_route_matches,
        "search_box": {
            "degrees": arguments.degrees,
            "F_dot_O_prime_values": s_values,
            "O_dot_F_prime_values": t_values,
            "O_dot_O_prime_values": z_values,
            "nonnegative_intersections_only": True,
        },
        "edges": edges,
        "input_hashes": {
            relative(CERT_SCRIPT): digest(CERT_SCRIPT),
            relative(certificate["MASKED"]): digest(certificate["MASKED"]),
            relative(certificate["BRIDGES"]): digest(certificate["BRIDGES"]),
            relative(certificate["THETA"]): digest(certificate["THETA"]),
            relative(KNOWN_ROUTE): digest(KNOWN_ROUTE),
            relative(KNOWN_SOURCE): digest(KNOWN_SOURCE),
        },
        "proof_boundary": {
            "proved": (
                "Every tested pair is an exact representation of G_A in the current "
                "positive frame. Each reported hit gives a literal primitive U basis, "
                "an integral orthogonal frame, a unimodular full-NS transport, the exact "
                "declared root profile, and an integral isometry to the target frame."
            ),
            "not_proved": (
                "The search is complete only for degrees and intersection intervals in "
                "search_box. A hit is a lattice-level marking; nefness, irreducibility "
                "and effectivity of the displayed zero, horizontal walls, equations, "
                "and rational maps are not certified."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/search_ns0024_relative_u_bridge_lifts.sage "
            + "--degrees "
            + " ".join(map(str, arguments.degrees))
            + " --s-values "
            + " ".join(map(str, s_values))
            + " --t-values "
            + " ".join(map(str, t_values))
            + " --z-values "
            + " ".join(map(str, z_values))
            + " --output "
            + relative(
                arguments.output
                if arguments.output.is_absolute()
                else ROOT / arguments.output
            )
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(relative(output))
    print(payload["status"])


if __name__ == "__main__":
    main()
