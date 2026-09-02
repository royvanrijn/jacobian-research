#!/usr/bin/env sage-python
"""Blind bridge/glue benchmark on the held-out determinant-78 E6 NS.

The four corridors used to formulate the bridge predictor are not inputs.
Starting only from the explicit ``2E6+A1/MW4`` frame, this script generates
the complete zero-neutral old-degree-two shell modulo the source root Weyl
group.  For every primitive candidate it computes the common core ``K``, the
new rank-two bridge ``C``, and the exact minimum in every nonzero glue coset
of ``K+C``.  Candidates are ranked by the weakest nonzero-coset minimum.

Only after the ranking is fixed does the evaluation phase enumerate child
roots and identify the child in the independently mass-closed 1,549-class J2
catalogue.  Thus neither a successful corridor nor a target frame is used to
generate cores or scores.  The experiment is a prospective negative control:
the truth catalogue contains no rootless class.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from itertools import product
from pathlib import Path
from time import perf_counter

from sage.all import (
    CartanMatrix,
    IntegralLattice,
    QQ,
    ZZ,
    floor,
    gcd,
    identity_matrix,
    matrix,
    pari,
    vector,
    zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
FRAME_PATH = ROOT / "elkies-k3/data/lattice/e6_rank4_det78_frame.txt"
TRUTH = GEN / "elkies-k3-e6-rank4-det78-niemeier-frames-v1.json"
OUTPUT = GEN / "elkies-k3-e6-rank4-det78-prospective-bridge-predictor-v1.json"

_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        pending = [min(unseen)]
        unseen.remove(pending[0])
        component = []
        while pending:
            index = pending.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other] != 0]
            for other in adjacent:
                unseen.remove(other)
                pending.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def component_label_vectors(cartan, norm_bound):
    inverse = cartan.inverse()
    bounds = [
        int(floor((QQ(norm_bound) / inverse[index, index]).sqrt()))
        for index in range(inverse.nrows())
    ]
    result = []
    for values in product(*[range(bound + 1) for bound in bounds]):
        label = vector(ZZ, values)
        norm = label * inverse * label
        if norm <= norm_bound:
            result.append((label, norm))
    return tuple(result)


def enumerate_dominant_norm_vectors(frame, simple, target_norm):
    """Complete source norm shell modulo its root Weyl group."""
    cartan = simple * frame * simple.transpose()
    component_labels = []
    for component in connected_components(cartan):
        block = cartan.matrix_from_rows_and_columns(component, component)
        component_labels.append(
            (component, component_label_vectors(block, target_norm))
        )
    combined = []
    for choices in product(*[records for unused, records in component_labels]):
        if sum((norm for unused, norm in choices), QQ(0)) > target_norm:
            continue
        label = vector(ZZ, simple.nrows())
        for (component, unused), (part, unused_norm) in zip(
            component_labels, choices
        ):
            for index, value in zip(component, part):
                label[index] = value
        combined.append(label)

    pairing = simple * frame
    smith, left, right = pairing.smith_form()
    assert smith == left * pairing * right
    rank = pairing.rank()
    diagonal = tuple(abs(ZZ(smith[index, index])) for index in range(rank))
    kernel = right[:, rank:]
    kernel_gram = kernel.transpose() * frame * kernel
    lattice = IntegralLattice(kernel_gram)
    answers = []
    for label in combined:
        rhs = left * label
        if any(rhs[index] % diagonal[index] for index in range(rank)):
            continue
        coordinates = vector(ZZ, frame.nrows())
        for index in range(rank):
            coordinates[index] = rhs[index] // diagonal[index]
        particular = right * coordinates
        centre = -kernel_gram.inverse() * kernel.transpose() * frame * particular
        for close in lattice.enumerate_close_vectors(centre):
            candidate = particular + kernel * vector(ZZ, close)
            norm = candidate * frame * candidate
            if norm > target_norm:
                break
            if norm == target_norm:
                assert pairing * candidate == label
                answers.append(tuple(map(int, candidate)))
    return tuple(vector(ZZ, item) for item in sorted(set(answers)))


def integral_ns():
    """Rebuild the source O,F,2E6,A1,P,Q,R1,R2 marking."""
    ns = zero_matrix(ZZ, 19)
    ns[0, 0] = -2
    ns[0, 1] = ns[1, 0] = 1
    ns[2:8, 2:8] = -CartanMatrix(["E", 6])
    ns[8:14, 8:14] = -CartanMatrix(["E", 6])
    ns[14, 14] = -2
    for section_index in range(15, 19):
        ns[section_index, section_index] = -2
        ns[1, section_index] = ns[section_index, 1] = 1
    for section_index in (15, 16):
        for component_index in (2, 8):
            ns[section_index, component_index] = 1
            ns[component_index, section_index] = 1
    for section_index in (17, 18):
        for component_index in (2, 14):
            ns[section_index, component_index] = 1
            ns[component_index, section_index] = 1
    for left, right in ((15, 18), (16, 17)):
        ns[left, right] = ns[right, left] = 1
    return ns


def minimum_norm(gram):
    """Exact lattice minimum, using PARI only to enumerate the first shell."""
    data = pari(gram).qfminim(2)
    if int(data[0]):
        return 2
    bound = 4
    while True:
        data = pari(gram).qfminim(bound)
        if int(data[0]):
            candidates = matrix(ZZ, data[2].sage()).columns()
            return min(int(vector(ZZ, item) * gram * vector(ZZ, item)) for item in candidates)
        bound += 2


def exact_closest_coset_norm(gram, target):
    """Minimize (z-target) G (z-target) exactly after PARI CVP proposes z."""
    closest = pari(gram).qfcvp(pari(target), flag=0)
    candidates = matrix(ZZ, closest[2].sage()).columns()
    if not candidates:
        raise ArithmeticError("PARI qfcvp returned no closest vector")
    exact = [
        (vector(QQ, item) - target) * gram * (vector(QQ, item) - target)
        for item in candidates
    ]
    result = min(exact)
    # PARI's small-integral-Gram CVP uses floating Fincke--Pohst internally;
    # the returned distance is always re-evaluated in QQ here.
    if result < 0 or result.denominator() != 1:
        raise ArithmeticError("glue-coset minimum is not a nonnegative integer")
    return int(result)


def prospective_bridge_score(old_basis, child_basis, ns):
    """Compute K, C_new and every nonzero graph-glue coset minimum."""
    old_module = old_basis.row_module(ZZ)
    child_module = child_basis.row_module(ZZ)
    core_basis = old_module.intersection(child_module).basis_matrix()
    if core_basis.nrows() != 15:
        raise ArithmeticError("prospective click does not have rank-15 core")

    child = -(child_basis * ns * child_basis.transpose())
    core_coordinates = child_basis.solve_left(core_basis).change_ring(ZZ)
    bridge_coordinates = (core_coordinates * child).right_kernel_matrix()
    if bridge_coordinates.nrows() != 2:
        raise ArithmeticError("prospective click does not have rank-two bridge")
    split_coordinates = core_coordinates.stack(bridge_coordinates)
    split_gram = split_coordinates * child * split_coordinates.transpose()
    index = abs(int(split_coordinates.det()))
    diagonal_matrix, smith_left, smith_right = split_coordinates.smith_form()
    assert smith_left * split_coordinates * smith_right == diagonal_matrix
    diagonal = [abs(int(value)) for value in diagonal_matrix.diagonal()]
    nontrivial = [(i, order) for i, order in enumerate(diagonal) if order > 1]
    if len(nontrivial) != 1:
        return {
            "admissible": False,
            "reason": "noncyclic K_plus_C quotient",
            "glue_group_invariants": [order for unused, order in nontrivial],
        }
    coordinate_index, order = nontrivial[0]
    generator = vector(ZZ, identity_matrix(ZZ, 17).row(coordinate_index))
    frame_generator = generator * smith_right.inverse().change_ring(ZZ)
    split_generator = frame_generator * split_coordinates.inverse()
    assert order * split_generator in ZZ**17

    lll = matrix(ZZ, pari(split_gram).qflllgram()).transpose()
    reduced_gram = lll * split_gram * lll.transpose()
    inverse_lll = lll.inverse().change_ring(ZZ)
    coset_minima = []
    for multiplier in range(1, order // 2 + 1):
        fractional = multiplier * split_generator
        fractional -= vector(ZZ, [value.floor() for value in fractional])
        target = -fractional * inverse_lll
        coset_minima.append(exact_closest_coset_norm(reduced_gram, target))
    if not coset_minima:
        raise ArithmeticError("nontrivial glue quotient has no nonzero coset")

    core_gram = -(core_basis * ns * core_basis.transpose())
    bridge_gram = bridge_coordinates * child * bridge_coordinates.transpose()
    return {
        "admissible": True,
        "core_minimum": minimum_norm(core_gram),
        "bridge_minimum": minimum_norm(bridge_gram),
        "glue_group_invariants": [order],
        "glue_coset_minimum": min(coset_minima),
        "glue_coset_minima_histogram": {
            str(key): value for key, value in sorted(Counter(coset_minima).items())
        },
    }


def truth_class(child, root_rank, signed_root_count, truth_by_signature):
    candidates = truth_by_signature[(root_rank, signed_root_count)]
    matches = [
        row
        for row in candidates
        if pari(child).qfisom(pari(matrix(ZZ, row["gram"]))) != 0
    ]
    if len(matches) != 1:
        raise ArithmeticError(
            f"expected one J2 truth match, found {len(matches)} for "
            f"signature {(root_rank, signed_root_count)}"
        )
    return matches[0]


def top_k_summary(ranked, k):
    selected = ranked[: min(k, len(ranked))]
    return {
        "k": k,
        "candidate_count": len(selected),
        "distinct_j2_classes": len({row["j2_class_index"] for row in selected}),
        "root_rank_distribution": {
            str(key): value
            for key, value in sorted(Counter(row["child_root_rank"] for row in selected).items())
        },
        "minimum_child_root_rank": min(row["child_root_rank"] for row in selected),
        "rootless_count": sum(row["child_root_rank"] == 0 for row in selected),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--limit",
        type=int,
        help="development-only prefix limit; cannot write or check a certificate",
    )
    arguments = parser.parse_args()
    if arguments.limit is not None and arguments.limit <= 0:
        parser.error("--limit must be positive")
    if arguments.limit is not None and arguments.check:
        parser.error("--limit cannot be combined with --check")

    source = json.loads(SOURCE.read_text())
    assert source.get("status", "").startswith("PASS_EXACT")

    ns = integral_ns()
    old_fibre = vector(ZZ, [0, 1] + [0] * 17)
    old_zero = vector(ZZ, [1] + [0] * 18)
    constraints = matrix(ZZ, [ns * old_fibre, ns * (old_zero + old_fibre)])
    old_basis = constraints.right_kernel().basis_matrix()
    old_transport = matrix(ZZ, [old_fibre, old_zero + old_fibre, *old_basis.rows()])
    assert abs(old_transport.det()) == 1
    frame = -(old_basis * ns * old_basis.transpose())
    assert frame == load_gram(FRAME_PATH) and frame.det() == 78
    simple, unused_positive, unused_cartan = deterministic_simple_roots(frame)
    assert roots_and_data(frame)[2] == (13, 146, 18)

    horizontals = enumerate_dominant_norm_vectors(frame, simple, ZZ(8))
    assert len(horizontals) == 280
    if arguments.limit is not None:
        horizontals = horizontals[: arguments.limit]

    # Phase A is blind: generate and rank without opening the J2 truth file or
    # enumerating a child root.
    scored = []
    nonprimitive = 0
    score_seconds = 0.0
    for orbit, horizontal in enumerate(horizontals):
        divisor_split = vector(ZZ, [2, 2] + [-entry for entry in horizontal])
        divisor = divisor_split * old_transport
        if gcd(list(divisor)) != 1:
            nonprimitive += 1
            continue
        neighbor = primitive_hyperbolic_split(ns, divisor)
        child = matrix(ZZ, neighbor["child_frame"])
        child_basis = matrix(ZZ, neighbor["transport"])[2:, :]

        started = perf_counter()
        score = prospective_bridge_score(old_basis, child_basis, ns)
        score_seconds += perf_counter() - started
        if not score["admissible"]:
            raise ArithmeticError(score["reason"])

        scored.append(
            {
                "source_weyl_orbit": orbit,
                **score,
                "_child": child,
            }
        )
        if len(scored) % 25 == 0:
            print(
                f"DET78PROSPECTIVE|scored={len(scored)}|"
                f"glue_min={score['glue_coset_minimum']}|status=BLIND",
                flush=True,
            )

    if arguments.limit is None:
        assert nonprimitive == 3 and len(scored) == 277
    ranked_internal = sorted(
        scored,
        key=lambda row: (-row["glue_coset_minimum"], row["source_weyl_orbit"]),
    )

    # Phase B opens the independently mass-closed truth set only after the
    # complete prospective order has been fixed.
    truth = json.loads(TRUTH.read_text())
    assert truth["status"] == "PASS_COMPLETE_DET78_J2_FRAME_CLASSIFICATION"
    assert truth["accounting"]["genus_mass_closed"]
    assert truth["accounting"]["rootless_frame_classes"] == 0
    truth_by_signature = defaultdict(list)
    for row in truth["frames"]:
        truth_by_signature[(int(row["root_rank"]), int(row["signed_root_count"]))].append(row)

    ranked = []
    root_classification_seconds = 0.0
    j2_matching_seconds = 0.0
    for rank_index, scored_row in enumerate(ranked_internal, 1):
        child = scored_row.pop("_child")
        started = perf_counter()
        root_result = pari(child).qfminim(2)
        signed_root_count = int(root_result[0])
        root_vectors = matrix(ZZ, root_result[2].sage()).columns()
        signed_roots = [vector(ZZ, item) for item in root_vectors]
        signed_roots += [-item for item in signed_roots]
        root_rank = matrix(ZZ, [list(item) for item in signed_roots]).rank()
        root_classification_seconds += perf_counter() - started

        started = perf_counter()
        truth_row = truth_class(child, root_rank, signed_root_count, truth_by_signature)
        j2_matching_seconds += perf_counter() - started
        scored_row.update(
            {
                "prospective_rank": rank_index,
                "child_root_rank": int(root_rank),
                "child_signed_root_count": signed_root_count,
                "j2_class_index": int(truth_row["class_index"]),
                "j2_root_type": truth_row["root_type"],
            }
        )
        ranked.append(scored_row)
        if rank_index % 25 == 0:
            print(
                f"DET78PROSPECTIVE|evaluated={rank_index}|"
                f"root_rank={root_rank}|status=TRUTH",
                flush=True,
            )
    score_distribution = Counter(row["glue_coset_minimum"] for row in ranked)
    best_score = max(score_distribution)
    best = [row for row in ranked if row["glue_coset_minimum"] == best_score]
    best_root_rank = min(row["child_root_rank"] for row in ranked)
    best_truth = [row for row in ranked if row["child_root_rank"] == best_root_rank]
    root_rank_distribution = Counter(row["child_root_rank"] for row in ranked)
    distinct_j2_classes = {row["j2_class_index"] for row in ranked}
    if arguments.limit is None:
        assert score_distribution == Counter({2: 277})
        assert root_rank_distribution == Counter({12: 18, 13: 168, 14: 71, 15: 20})
        assert len(distinct_j2_classes) == 31

    payload = {
        "schema": "elkies-k3.e6-det78-prospective-bridge-predictor.v1",
        "status": "PASS_BLIND_DET78_PROSPECTIVE_BRIDGE_PREDICTOR_NEGATIVE_CONTROL",
        "inputs": {
            relative(SOURCE): digest(SOURCE),
            relative(FRAME_PATH): digest(FRAME_PATH),
            relative(TRUTH): digest(TRUTH),
        },
        "holdout_design": {
            "held_out_ns": "determinant-78 E6 rank-four Neron--Severi lattice",
            "training_excluded": ["H3", "Q80", "NS0024", "Golay720"],
            "core_generation": (
                "Complete zero-neutral old-degree-two norm-eight shell modulo the "
                "source 2E6+A1 root Weyl group; no corridor or target is supplied."
            ),
            "predeclared_score": (
                "Descending exact minimum over the nonzero graph-glue cosets of "
                "K+C_new; source Weyl-orbit order breaks ties."
            ),
            "label_separation": (
                "Child root enumeration and exact isometry matching to the mass-closed "
                "J2 catalogue are performed only after each score is computed."
            ),
        },
        "truth_set": {
            "j2_frame_classes": int(truth["accounting"]["frame_isometry_classes"]),
            "rootless_frame_classes": int(truth["accounting"]["rootless_frame_classes"]),
            "genus_mass": truth["accounting"]["target_genus_mass"],
            "mass_closed": bool(truth["accounting"]["genus_mass_closed"]),
        },
        "shell": {
            "dominant_orbits": len(horizontals),
            "nonprimitive_orbits": nonprimitive,
            "primitive_candidates": len(ranked),
            "distinct_j2_classes_reached": len(distinct_j2_classes),
            "glue_coset_minimum_distribution": {
                str(key): value for key, value in sorted(score_distribution.items())
            },
            "child_root_rank_distribution": {
                str(key): value
                for key, value in sorted(root_rank_distribution.items())
            },
            "best_glue_coset_minimum": best_score,
            "best_score_candidate_count": len(best),
            "best_score_root_rank_distribution": {
                str(key): value
                for key, value in sorted(Counter(row["child_root_rank"] for row in best).items())
            },
            "minimum_child_root_rank": best_root_rank,
            "minimum_root_rank_candidate_count": len(best_truth),
            "minimum_root_rank_candidates_at_best_score": sum(
                row["glue_coset_minimum"] == best_score for row in best_truth
            ),
            "top_k": [top_k_summary(ranked, k) for k in (1, 5, 10, 25, 50, 100)],
        },
        "timing_seconds": {
            "prospective_core_bridge_and_glue_coset_scoring": score_seconds,
            "direct_child_root_classification": root_classification_seconds,
            "independent_j2_truth_matching": j2_matching_seconds,
            "scoring_over_direct_root_classification_ratio": float(
                score_seconds / root_classification_seconds
            ),
            "boundary": (
                "Wall-clock timings are informative workstation measurements, not "
                "mathematical certificate fields."
            ),
        },
        "ranked_candidates": ranked,
        "conclusion": {
            "construction_algorithm_gate": "not_testable_on_rootless_positive_cases",
            "observed_answer": (
                "This genuinely held-out, mass-complete NS is a negative control: its "
                "J2 truth set has no rootless frame, so no unseen rootless endpoint can "
                "be recovered. The score is assessed only for enrichment toward the "
                "lowest root rank present in the declared source shell and for runtime."
            ),
            "next_positive_gate": (
                "Repeat the identical predeclared protocol on a held-out NS with a "
                "mass-complete J2 catalogue containing at least one rootless class."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "benchmark_e6_det78_prospective_bridge_predictor.sage --check"
        ),
    }
    if arguments.limit is not None:
        print(json.dumps(payload["shell"], indent=2, sort_keys=True))
        return

    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        stored = json.loads(output.read_text())
        # Timings are deliberately non-certificate metadata.
        payload["timing_seconds"] = stored["timing_seconds"]
        encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS blind determinant-78 prospective bridge predictor")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
