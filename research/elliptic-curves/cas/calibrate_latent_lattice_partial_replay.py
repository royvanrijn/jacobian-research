#!/usr/bin/env python3
"""Calibrate exact proper-subspace replay on the rank-25 R17 control.

No wgxli record is loaded.  The published R17 embedding is used only to build
one supervised rank-16 relation path.  The resulting partial map is saturated,
lifted, Smith-tested, replayed on the full finite clouds, and decorated with
source-free good-reduction code invariants.  A separately bounded oracle-center
beam then tests whether the same invariant rescues the edgewise proposal
generator without supplying the path itself.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from latent_lattice import (  # noqa: E402
    EllipticCurve,
    FiniteQuotientBlock,
    bounded_metric_relation_search,
    build_relation_complex,
    canonical_unoriented,
    exact_partial_relation_replay,
    height_gram,
    partial_replay_finite_signature,
    rational_rank,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
CONSENSUS = ARTIFACTS / "latent_lattice_relation_consensus_v1.json"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
FINITE = ARTIFACTS / "latent_lattice_finite_calibration_v1.json"
BASELINE = ARTIFACTS / "latent_lattice_metric_relation_search_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_partial_replay_v1.json"
LABEL = "rank_at_least_25"
FORCED_GROWTH_STEPS = 287
SEARCH_BOUNDS = {
    "source_center_indices": [116],
    "target_center_indices": [15],
    "seed_edges_per_center": 8,
    "initial_states_per_center_pair": 128,
    "minimum_states_per_center_pair": 1,
    "beam_width": 128,
    "maximum_steps": 350,
    "maximum_exact_lift_attempts": 500,
    "maximum_embeddings": 1,
    "minimum_global_replay_rays": 100,
    "reseed_after_skips": 20,
    "reseed_source_limit": 16,
    "reseed_target_limit": 4,
    "reseed_state_limit": 32,
    "norm_log_tolerance": 0.45,
    "angle_tolerance": 0.16,
    "angle_hard_tolerance": 0.36,
    "maximum_partial_replay_attempts": 400,
    "minimum_partial_replay_rank": 10,
    "partial_replay_candidate_limit": 128,
    "partial_replay_rank_stride": 2,
    "partial_replay_vertex_stride": 12,
    "preserve_initial_state_lineages": True,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def finite_block(record: dict[str, object]) -> FiniteQuotientBlock:
    return FiniteQuotientBlock(
        reduction_prime=int(record["reduction_prime"]),
        relation_prime=int(record["relation_prime"]),
        group_order=int(record["group_order"]),
        multiple_subgroup_order=int(record["multiple_subgroup_order"]),
        quotient_dimension=int(record["quotient_dimension"]),
        rows=tuple(tuple(map(int, row)) for row in record["rows"]),
    )


def gram_inputs(truth_document):
    forms = []
    for record in truth_document["positive_controls"]:
        if record["label"] == LABEL:
            continue
        form = np.asarray(record["canonical_height_gram"], dtype=float)
        forms.append(form / np.trace(form))
    inferred = sum(forms) / len(forms)
    module = importlib.import_module("elkies_rank25")
    target = height_gram(
        EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS)),
        tuple(module.POINTS),
        digits=60,
    )
    return inferred, target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    consensus = json.loads(CONSENSUS.read_text())
    truth_document = json.loads(TRUTH.read_text())
    finite_document = json.loads(FINITE.read_text())
    held = next(
        record
        for record in consensus["leave_one_fibre_out"]
        if record["held_out_label"] == LABEL
    )
    control = next(
        record for record in consensus["controls"] if record["label"] == LABEL
    )
    truth = next(
        record
        for record in truth_document["positive_controls"]
        if record["label"] == LABEL
    )
    finite = next(
        record for record in finite_document["controls"] if record["label"] == LABEL
    )
    source = build_relation_complex(held["training_two_of_three_vectors"])
    target = build_relation_complex(
        tuple(
            tuple(map(int, row.split()))
            for row in control["ambient_short_vector_coordinate_rows"]
        )
    )
    matrix = tuple(
        tuple(map(int, row)) for row in truth["embedding_matrix_columns"]
    )
    target_index = {vector: index for index, vector in enumerate(target.vertices)}
    truth_map = []
    for vector in source.vertices:
        image = canonical_unoriented(
            tuple(
                sum(int(vector[row]) * matrix[row][column] for row in range(17))
                for column in range(len(matrix[0]))
            )
        )
        truth_map.append(target_index.get(image))

    source_by_vertex = [[] for _ in source.vertices]
    for edge_index, edge in enumerate(source.ternary_relations):
        for vertex in edge:
            source_by_vertex[vertex].append(edge_index)
    target_edges = set(target.ternary_relations)
    source_center = SEARCH_BOUNDS["source_center_indices"][0]
    seed_candidates = sorted(
        source_by_vertex[source_center],
        key=lambda edge: (
            -sum(
                source.additive_degrees[index]
                for index in source.ternary_relations[edge]
            ),
            edge,
        ),
    )[: SEARCH_BOUNDS["seed_edges_per_center"]]
    seed_edge = next(
        edge
        for edge in seed_candidates
        if all(truth_map[index] is not None for index in source.ternary_relations[edge])
        and tuple(sorted(truth_map[index] for index in source.ternary_relations[edge]))
        in target_edges
    )
    mapped = {
        index: truth_map[index] for index in source.ternary_relations[seed_edge]
    }
    processed = {seed_edge}
    for _step in range(FORCED_GROWTH_STEPS):
        frontier = set()
        for vertex in mapped:
            frontier.update(source_by_vertex[vertex])
        frontier.difference_update(processed)
        mapped_set = set(mapped)
        edge_index = max(
            frontier,
            key=lambda edge: (
                len(mapped_set & set(source.ternary_relations[edge])),
                sum(
                    source.additive_degrees[index]
                    for index in source.ternary_relations[edge]
                ),
                -edge,
            ),
        )
        processed.add(edge_index)
        edge = source.ternary_relations[edge_index]
        if (
            all(truth_map[index] is not None for index in edge)
            and tuple(sorted(truth_map[index] for index in edge)) in target_edges
        ):
            mapped.update({index: truth_map[index] for index in edge})
    partial_map = [-1] * len(source.vertices)
    for source_index, target_vertex in mapped.items():
        partial_map[source_index] = target_vertex
    replays = exact_partial_relation_replay(source, target, partial_map)
    if not replays:
        raise ArithmeticError("supervised rank-16 component no longer lifts")
    replay = max(
        replays,
        key=lambda item: (
            len(item.replayed_source_vertex_indices),
            item.replayed_relation_count,
        ),
    )
    blocks = tuple(
        finite_block(record)
        for record in finite["development_blocks"] + finite["held_out_blocks"]
    )
    finite_signature = partial_replay_finite_signature(
        replay, target, finite_blocks=blocks
    )

    source_gram, target_gram = gram_inputs(truth_document)
    ledger = bounded_metric_relation_search(
        source,
        target,
        source_gram,
        target_gram,
        **SEARCH_BOUNDS,
    )
    exact_validator_passed = (
        len(mapped) == 103
        and rational_rank(tuple(source.vertices[index] for index in mapped)) == 16
        and replay.source_rank == 16
        and replay.source_subspace_ray_count == 362
        and len(replay.replayed_source_vertex_indices) == 194
        and replay.replayed_relation_count == 318
        and replay.primitive_target_image
    )
    selector_passed = bool(ledger.embeddings)
    status = (
        "PASS_EXACT_PARTIAL_REPLAY_AND_SELECTOR"
        if exact_validator_passed and selector_passed
        else "PASS_EXACT_PARTIAL_REPLAY_SELECTOR_FAIL"
        if exact_validator_passed
        else "FAIL_EXACT_PARTIAL_REPLAY_VALIDATOR"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-partial-replay.v1",
        "status": status,
        "scope": "Phase-0 rank-25 R17 control only; no wgxli target is loaded",
        "supervised_exact_validator": {
            "source_center_index": source_center,
            "target_center_index": SEARCH_BOUNDS["target_center_indices"][0],
            "seed_edge_index": seed_edge,
            "forced_growth_steps": FORCED_GROWTH_STEPS,
            "mapped_ray_count": len(mapped),
            "mapped_source_rank": rational_rank(
                tuple(source.vertices[index] for index in mapped)
            ),
            "intrinsic_relation_count": replay.intrinsic_relation_count,
            "source_subspace_ray_count": replay.source_subspace_ray_count,
            "exact_partial_lift_count": len(replays),
            "replayed_ray_count": len(replay.replayed_source_vertex_indices),
            "replayed_source_rank": replay.replayed_source_rank,
            "replayed_relation_count": replay.replayed_relation_count,
            "target_smith_invariant_factors": list(
                replay.target_smith_invariant_factors
            ),
            "primitive_target_image": replay.primitive_target_image,
            "finite_signature": finite_signature.to_record(),
        },
        "oracle_center_bounded_selector": {
            "bounds": SEARCH_BOUNDS,
            "initial_state_count": ledger.initial_state_count,
            "expanded_state_count": ledger.expanded_state_count,
            "partial_replay_attempt_count": ledger.partial_replay_attempt_count,
            "maximum_partial_replay_ray_count": (
                ledger.maximum_partial_replay_ray_count
            ),
            "maximum_partial_replay_relation_count": (
                ledger.maximum_partial_replay_relation_count
            ),
            "maximum_global_replay_ray_count": ledger.maximum_global_replay_ray_count,
            "accepted_embedding_count": len(ledger.embeddings),
        },
        "baseline_blind_search_artifact": str(BASELINE.relative_to(ROOT)),
        "gate_decision": (
            "CLOSED. Exact partial replay certifies a 194-ray primitive rank-16 "
            "truth component, but even the declared oracle-center, seed-lineage "
            "beam loses that component and accepts no full embedding. Replace the "
            "edgewise continuation invariant before any target use."
        ),
        "proof_boundary": (
            "Primitive source closure, intrinsic integer coordinates, ternary "
            "lifting, Smith factors, global replay, finite quotient restrictions, "
            "and all counts in the supervised validator are exact within the "
            "supplied finite clouds. The path and center pair use withheld R17 "
            "truth and are validation only. Height pruning and finite-profile "
            "comparison are heuristic. The failed bounded selector disproves only "
            "this edgewise proposal box, not existence of the R17 subgroup."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CONSENSUS,
                TRUTH,
                FINITE,
                BASELINE,
                ELLIPTIC / "cas/elkies_rank25.py",
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice partial-replay artifact is stale")
        print(
            "LATENTPARTIAL|check=PASS|"
            f"sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTPARTIAL|status={status}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
