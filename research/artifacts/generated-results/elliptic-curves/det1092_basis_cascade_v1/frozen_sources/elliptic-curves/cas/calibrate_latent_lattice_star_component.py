#!/usr/bin/env python3
"""Freeze the finite-aware center-star R17 rank-25 calibration.

No wgxli target is loaded.  A bounded center-star assignment uses height only
for proposal order, enforces exact mod-2/mod-3 rank compatibility, and audits
proper subspaces by exact integral replay.  The published embedding is revealed
afterward to measure recall and to certify the best possible visible truth star.
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
    bounded_metric_star_component_search,
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
OUTPUT = ARTIFACTS / "latent_lattice_star_component_v1.json"
LABEL = "rank_at_least_25"
SEARCH_BOUNDS = {
    "source_center_indices": [116],
    "target_center_indices": [15],
    "source_edge_limit": 20,
    "target_edges_per_source": 64,
    "beam_width": 512,
    "minimum_states_per_match_rank_group": 8,
    "maximum_partial_replay_attempts": 500,
    "minimum_partial_replay_rank": 10,
    "partial_replay_candidate_limit": 64,
    "minimum_partial_replay_rays": 0,
    "maximum_candidates": 500,
    "norm_log_tolerance": 0.45,
    "angle_tolerance": 0.16,
    "angle_hard_tolerance": 0.36,
    "finite_matroid_primes": [2, 3],
    "finite_matroid_subset_size": 1,
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
    target_edges = set(target.ternary_relations)
    source_center = SEARCH_BOUNDS["source_center_indices"][0]
    incident_edges = tuple(
        edge_index
        for edge_index, edge in enumerate(source.ternary_relations)
        if source_center in edge
    )
    visible_truth_edges = tuple(
        edge_index
        for edge_index in incident_edges
        if all(
            truth_map[index] is not None
            for index in source.ternary_relations[edge_index]
        )
        and tuple(
            sorted(
                truth_map[index]
                for index in source.ternary_relations[edge_index]
            )
        )
        in target_edges
    )
    truth_star_indices = tuple(
        sorted(
            {
                index
                for edge_index in visible_truth_edges
                for index in source.ternary_relations[edge_index]
            }
        )
    )
    truth_star_map = [-1] * len(source.vertices)
    for index in truth_star_indices:
        truth_star_map[index] = truth_map[index]
    truth_star_replays = exact_partial_relation_replay(
        source, target, truth_star_map
    )
    if not truth_star_replays:
        raise ArithmeticError("visible supervised truth star no longer lifts")
    truth_star = max(
        truth_star_replays,
        key=lambda replay: (
            len(replay.replayed_source_vertex_indices),
            replay.replayed_relation_count,
        ),
    )
    blocks = tuple(
        finite_block(record)
        for record in finite["development_blocks"] + finite["held_out_blocks"]
    )
    truth_star_finite = partial_replay_finite_signature(
        truth_star, target, finite_blocks=blocks
    )

    training_forms = []
    for record in truth_document["positive_controls"]:
        if record["label"] == LABEL:
            continue
        form = np.asarray(record["canonical_height_gram"], dtype=float)
        training_forms.append(form / np.trace(form))
    source_gram = sum(training_forms) / len(training_forms)
    module = importlib.import_module("elkies_rank25")
    target_gram = height_gram(
        EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS)),
        tuple(module.POINTS),
        digits=60,
    )
    ledger = bounded_metric_star_component_search(
        source,
        target,
        source_gram,
        target_gram,
        **SEARCH_BOUNDS,
    )
    postselection = []
    for source_rank, candidate in enumerate(ledger.candidates):
        candidate_rank = candidate.replay.source_rank
        overlap = candidate_rank + 17 - rational_rank(
            tuple(candidate.replay.integral_matrix) + matrix
        )
        postselection.append(
            {
                "source_rank": source_rank,
                "candidate_rank": candidate_rank,
                "withheld_truth_overlap": overlap,
                "replayed_ray_count": len(
                    candidate.replay.replayed_source_vertex_indices
                ),
                "replayed_relation_count": candidate.replay.replayed_relation_count,
            }
        )
    postselection.sort(
        key=lambda item: (
            item["withheld_truth_overlap"],
            item["replayed_ray_count"],
            item["replayed_relation_count"],
        ),
        reverse=True,
    )
    maximum_truth_overlap = max(
        (item["withheld_truth_overlap"] for item in postselection), default=0
    )
    truth_recalled = maximum_truth_overlap >= truth_star.source_rank
    status = (
        "PASS_STAR_COMPONENT_RECALL"
        if truth_recalled
        else "FAIL_STAR_COMPONENT_RECALL_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-star-component.v1",
        "status": status,
        "scope": "Phase-0 rank-25 R17 control only; no wgxli target is loaded",
        "supervised_truth_star_ceiling": {
            "source_center_index": source_center,
            "target_center_index": SEARCH_BOUNDS["target_center_indices"][0],
            "visible_incident_relation_count": len(visible_truth_edges),
            "mapped_ray_count": len(truth_star_indices),
            "mapped_rank": rational_rank(
                tuple(source.vertices[index] for index in truth_star_indices)
            ),
            "exact_lift_count": len(truth_star_replays),
            "source_subspace_ray_count": truth_star.source_subspace_ray_count,
            "replayed_ray_count": len(truth_star.replayed_source_vertex_indices),
            "replayed_relation_count": truth_star.replayed_relation_count,
            "primitive_target_image": truth_star.primitive_target_image,
            "finite_signature": truth_star_finite.to_record(),
        },
        "bounded_search": {
            "bounds": SEARCH_BOUNDS,
            "processed_star_layer_count": ledger.processed_star_layer_count,
            "maximum_beam_population": ledger.maximum_beam_population,
            "expanded_state_count": ledger.expanded_state_count,
            "finite_rank_rejection_count": ledger.finite_rank_rejection_count,
            "partial_replay_attempt_count": ledger.partial_replay_attempt_count,
            "maximum_source_rank_reached": ledger.maximum_source_rank_reached,
            "maximum_mapped_vertex_count": ledger.maximum_mapped_vertex_count,
            "maximum_partial_replay_ray_count": (
                ledger.maximum_partial_replay_ray_count
            ),
            "maximum_partial_replay_relation_count": (
                ledger.maximum_partial_replay_relation_count
            ),
            "candidate_count": len(ledger.candidates),
            "maximum_withheld_truth_overlap": maximum_truth_overlap,
            "truth_star_recalled": truth_recalled,
            "top_ten_withheld_overlap_diagnostic": postselection[:10],
        },
        "rank_three_matroid_diagnostic": {
            "status": "REJECTED_AS_DEFAULT",
            "declared_change": "finite_matroid_subset_size=3 with all other bounds fixed",
            "elapsed_seconds": "150.6614167690277",
            "finite_rank_rejection_count": 4152,
            "maximum_withheld_truth_overlap": 9,
            "conclusion": (
                "Pair/triple ranks rejected no branch beyond whole-state ranks "
                "in this box and did not improve recall."
            ),
        },
        "gate_decision": (
            "CLOSED. The exact truth star has rank 11, but the bounded, "
            "finite-aware star ledger reaches withheld overlap only "
            f"{maximum_truth_overlap}. Whole-star assignment is not yet a "
            "calibrated proposal generator."
        ),
        "proof_boundary": (
            "Ternary relations, mod-2/mod-3 ranks, primitive closures, integral "
            "lifts, Smith factors, replay, finite signatures, and withheld overlap "
            "dimensions are exact within the supplied clouds and declared search "
            "box. Height ordering is numerical. The center pair and postselection "
            "overlap use withheld R17 truth for calibration only. This bounded "
            "failure is not a nonexistence result for R17 or any target lattice."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CONSENSUS,
                TRUTH,
                FINITE,
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
            raise SystemExit("latent-lattice star-component artifact is stale")
        print(
            "LATENTSTAR|check=PASS|"
            f"sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTSTAR|status={status}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
