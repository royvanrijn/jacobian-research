#!/usr/bin/env python3
"""Calibrate joint rank-16 relation components before rank-17 completion.

No wgxli record is loaded.  Candidate rank-16 spaces are proposed inside a
fixed, source-free pool of rank-17 R17-control proposals, globally replayed,
and deduplicated by their exact rational spaces.  Finite ranks and a generous
height-angle bound are rejection filters.  Compact relation/height/code
profiles are then matched jointly across the four fibres.  Only after the
joint selection and one-rank completion are fixed are published R17
coordinates revealed for audit.
"""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC))

from latent_lattice import (  # noqa: E402
    CandidateRelationFingerprint,
    FiniteQuotientBlock,
    bounded_dense_hyperplanes,
    build_relation_complex,
    candidate_finite_signature,
    candidate_finite_signatures,
    candidate_relation_fingerprint,
    height_angle_profile,
    height_angle_profile_distance,
    joint_component_bundle_ledger,
    joint_nearest_candidate_scores,
    modular_rank,
    primitive_hermite_signatures,
    rational_rank,
    replay_and_deduplicate_components,
    row_basis_coordinates,
    row_embedding_smith_invariant_factors,
)
from latent_lattice.subspace import (  # noqa: E402
    exact_span_mask,
    independent_row_basis,
    primitive_span_basis,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
CONSENSUS = ARTIFACTS / "latent_lattice_relation_consensus_v1.json"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
GRAPH = ARTIFACTS / "latent_lattice_graph_walk_calibration_v1.json"
FINITE = ARTIFACTS / "latent_lattice_finite_calibration_v1.json"
PROPOSALS = ARTIFACTS / "latent_lattice_joint_fingerprint_ledger_v1.json.gz"
OUTPUT = ARTIFACTS / "latent_lattice_joint_component_calibration_v1.json"
LABELS = tuple(f"rank_at_least_{rank}" for rank in range(25, 29))
SOURCE_LABEL = "rank_at_least_25"
SOURCE_SAMPLE_COUNT = 2_000
TARGET_SAMPLE_COUNT = 400
TARGET_COMPONENTS_PER_PARENT = 8
PARENT_POOL_SIZE = 16
RANDOM_SEED = 20_260_901
DEVELOPMENT_NUMERATOR = 3
DEVELOPMENT_DENOMINATOR = 4
MAXIMUM_SOURCE_ANGLE_DISTANCE = 0.03
FINGERPRINT_QUANTILES = 16
FINGERPRINT_PROJECTIVE_MULTIPLICITIES = 16
STRUCTURAL_BUNDLE_QUANTILE = 0.10


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


def ambient_rows(component_rows, parent_basis):
    raw = tuple(
        tuple(
            sum(int(vector[row]) * int(parent_basis[row][column]) for row in range(len(parent_basis)))
            for column in range(len(parent_basis[0]))
        )
        for vector in component_rows
    )
    return primitive_span_basis(independent_row_basis(raw))


def embedded_rows(left, right):
    return tuple(
        tuple(
            sum(int(left[row][middle]) * int(right[middle][column]) for middle in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def percentile_ranks(values, *, larger_is_better=True):
    ordered = sorted(set(values), reverse=larger_is_better)
    denominator = max(1, len(ordered) - 1)
    return tuple(1.0 - ordered.index(value) / denominator for value in values)


def extended_fingerprint(fingerprint, angle_profile):
    return CandidateRelationFingerprint(
        dimension=fingerprint.dimension,
        ray_count=fingerprint.ray_count,
        ternary_relation_count=fingerprint.ternary_relation_count,
        scaled_relation_count=fingerprint.scaled_relation_count,
        integral_ray_count=fingerprint.integral_ray_count,
        feature_names=fingerprint.feature_names
        + tuple(f"absolute_height_angle_q{index}" for index in range(len(angle_profile))),
        feature_values=fingerprint.feature_values
        + tuple(f"{float(value):.17g}" for value in angle_profile),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    consensus = json.loads(CONSENSUS.read_text())
    truth_document = json.loads(TRUTH.read_text())
    graph_document = json.loads(GRAPH.read_text())
    finite_document = json.loads(FINITE.read_text())
    with gzip.open(PROPOSALS, "rt") as stream:
        proposal_document = json.load(stream)

    held = next(
        record
        for record in consensus["leave_one_fibre_out"]
        if record["held_out_label"] == SOURCE_LABEL
    )
    source_complex = build_relation_complex(held["training_two_of_three_vectors"])
    source_components = bounded_dense_hyperplanes(
        source_complex,
        sample_count=SOURCE_SAMPLE_COUNT,
        random_seed=RANDOM_SEED,
        maximum_components=8,
    )
    source_component = source_components.components[0]
    source_basis = primitive_span_basis(
        independent_row_basis(
            tuple(source_complex.vertices[index] for index in source_component.vertex_indices)
        )
    )
    training_forms = [
        np.asarray(record["canonical_height_gram"], dtype=float)
        for record in truth_document["positive_controls"]
        if record["label"] != SOURCE_LABEL
    ]
    source_gram = sum(form / np.trace(form) for form in training_forms) / len(training_forms)
    source_rows = tuple(source_complex.vertices[index] for index in source_component.vertex_indices)
    source_angle_profile = height_angle_profile(source_rows, source_gram)

    truth_by_label = {
        record["label"]: record for record in truth_document["positive_controls"]
    }
    control_by_label = {
        record["label"]: record for record in consensus["controls"] if record["label"] in LABELS
    }
    graph_by_label = {
        record["label"]: record for record in graph_document["default_controls"] if record["label"] in LABELS
    }
    finite_by_label = {
        record["label"]: record for record in finite_document["controls"] if record["label"] in LABELS
    }
    proposal_by_label = {
        record["label"]: record for record in proposal_document["fibres"] if record["label"] in LABELS
    }

    states = []
    for fibre_index, label in enumerate(LABELS):
        control = control_by_label[label]
        proposal_fibre = proposal_by_label[label]
        candidates_by_source = {
            int(record["source_index"]): record for record in proposal_fibre["candidates"]
        }
        graph_control = graph_by_label[label]
        parent_score_records = graph_control["default_ledger"]["candidates"][:PARENT_POOL_SIZE]
        rows = tuple(tuple(map(int, row.split())) for row in control["ambient_short_vector_coordinate_rows"])
        heights = tuple(map(float, control["ambient_short_vector_heights"]))
        full_complex = build_relation_complex(rows)
        order = tuple(sorted(range(len(rows)), key=lambda index: (heights[index], index)))
        cutoff = DEVELOPMENT_NUMERATOR * len(order) // DEVELOPMENT_DENOMINATOR
        development = tuple(order[:cutoff])
        held_out = tuple(order[cutoff:])
        proposed_bases = []
        origins = []
        parent_bases = {}
        parent_masks = {}
        for parent_position, parent_score in enumerate(parent_score_records):
            source_index = int(parent_score["source_index"])
            parent_basis = tuple(
                tuple(map(int, row)) for row in candidates_by_source[source_index]["basis_rows"]
            )
            parent_bases[source_index] = parent_basis
            parent_mask = exact_span_mask(rows, parent_basis)
            parent_masks[source_index] = parent_mask
            # Proposal cutoffs are adaptive inside each parent so comparable
            # numbers of parent rays participate.  Acceptance/replay below is
            # still measured against the one global fibre split.
            parent_order = tuple(
                index for index in order if parent_mask[index]
            )
            parent_cutoff = DEVELOPMENT_NUMERATOR * len(parent_order) // DEVELOPMENT_DENOMINATOR
            parent_development = parent_order[:parent_cutoff]
            coordinates = row_basis_coordinates(
                tuple(rows[index] for index in parent_development), parent_basis
            )
            parent_complex = build_relation_complex(coordinates)
            sampled = bounded_dense_hyperplanes(
                parent_complex,
                sample_count=TARGET_SAMPLE_COUNT,
                random_seed=RANDOM_SEED + source_index,
                maximum_components=TARGET_COMPONENTS_PER_PARENT,
            )
            for component_index, component in enumerate(sampled.components):
                basis = ambient_rows(
                    tuple(parent_complex.vertices[index] for index in component.vertex_indices),
                    parent_basis,
                )
                if len(basis) != 16:
                    continue
                proposed_bases.append(basis)
                origins.append(
                    {
                        "parent_source_index": source_index,
                        "parent_position": parent_position,
                        "component_index": component_index,
                        "sampled_development_ray_count": len(component.vertex_indices),
                        "sampled_development_relation_count": len(component.relation_indices),
                    }
                )
        deduplicated = replay_and_deduplicate_components(
            rows,
            full_complex,
            proposed_bases,
            development_indices=development,
            held_out_indices=held_out,
        )
        development_blocks = tuple(
            finite_block(record) for record in finite_by_label[label]["development_blocks"]
        )
        held_out_blocks = tuple(
            finite_block(record) for record in finite_by_label[label]["held_out_blocks"]
        )
        fingerprints = []
        development_finite_signatures = []
        candidate_records = []
        for candidate_index, candidate in enumerate(deduplicated):
            replayed = candidate.full_replayed_ray_indices
            angle_profile = height_angle_profile(
                tuple(rows[index] for index in replayed), proposal_fibre["ambient_height_gram"]
            )
            angle_distance = height_angle_profile_distance(source_angle_profile, angle_profile)
            relation_fingerprint = candidate_relation_fingerprint(
                rows,
                heights,
                tuple({} for _ in rows),
                replayed,
                full_complex,
                dimension=16,
                quantiles=FINGERPRINT_QUANTILES,
                projective_multiplicities=FINGERPRINT_PROJECTIVE_MULTIPLICITIES,
                finite_primes=(2, 3),
            )
            fingerprints.append(extended_fingerprint(relation_fingerprint, angle_profile))
            containing_parents = []
            for parent_score in parent_score_records:
                source_index = int(parent_score["source_index"])
                parent_basis = parent_bases[source_index]
                if rational_rank((*candidate.basis_rows, *parent_basis)) != 17:
                    continue
                parent_mask = parent_masks[source_index]
                parent_held = int(sum(parent_mask[index] for index in held_out))
                containing_parents.append(
                    {
                        "source_index": source_index,
                        "added_rank": 1,
                        "held_out_replayed_ray_count": parent_held,
                        "held_out_rays_gained": parent_held - len(candidate.held_out_replayed_ray_indices),
                        "held_out_rays_gained_per_added_rank": parent_held - len(candidate.held_out_replayed_ray_indices),
                        "full_replayed_ray_count": int(sum(parent_mask)),
                        "combined_graph_score": parent_score["combined_score"],
                        "shape_percentile": parent_score["shape_percentile"],
                    }
                )
            containing_parents.sort(
                key=lambda record: (
                    -record["held_out_rays_gained_per_added_rank"],
                    -float(record["combined_graph_score"]),
                    -float(record["shape_percentile"]),
                    record["source_index"],
                )
            )
            candidate_records.append(
                {
                    "candidate_index": candidate_index,
                    "primitive_embedding_matrix_rows": [list(row) for row in candidate.basis_rows],
                    "smith_invariant_factors": list(row_embedding_smith_invariant_factors(candidate.basis_rows)),
                    "origin_records": [origins[index] for index in candidate.origin_indices],
                    "origin_multiplicity": len(candidate.origin_indices),
                    "development_replayed_ray_count": len(candidate.development_replayed_ray_indices),
                    "held_out_replayed_ray_count": len(candidate.held_out_replayed_ray_indices),
                    "held_out_rays_per_rank": f"{candidate.held_out_rays_per_rank:.17g}",
                    "full_replayed_ray_count": len(replayed),
                    "full_replayed_relation_count": len(candidate.full_replayed_relation_indices),
                    "modular_ranks": {str(prime): rank for prime, rank in candidate.modular_ranks},
                    "finite_rank_compatible": all(rank == 16 for _prime, rank in candidate.modular_ranks),
                    "source_height_angle_profile_distance": f"{angle_distance:.17g}",
                    "passes_height_angle_filter": angle_distance <= MAXIMUM_SOURCE_ANGLE_DISTANCE,
                    "relation_height_fingerprint": fingerprints[-1].to_record(),
                    "rank17_completion_candidates": containing_parents,
                }
            )
        basis_matrices = tuple(
            tuple(tuple(map(int, row)) for row in record["primitive_embedding_matrix_rows"])
            for record in candidate_records
        )
        development_finite_signatures = list(
            candidate_finite_signatures(
                basis_matrices, full_complex, finite_blocks=development_blocks
            )
        )
        held_out_finite_signatures = candidate_finite_signatures(
            basis_matrices, full_complex, finite_blocks=held_out_blocks
        )
        hermite_signatures = primitive_hermite_signatures(
            proposal_fibre["ambient_height_gram"],
            basis_matrices,
            digits=60,
            batch_size=128,
        )
        hermite_values = []
        for record, signature, finite_development, finite_held_out in zip(
            candidate_records,
            hermite_signatures,
            development_finite_signatures,
            held_out_finite_signatures,
        ):
            record["primitive_hermite_signature"] = signature.to_record()
            record["development_finite_signature"] = finite_development.to_record()
            record["held_out_finite_signature"] = finite_held_out.to_record()
            hermite_values.append(float(signature.hermite.log_hermite_invariant))
        states.append(
            {
                "label": label,
                "rows": rows,
                "heights": heights,
                "complex": full_complex,
                "development": development,
                "held_out": held_out,
                "fingerprints": tuple(fingerprints),
                "development_finite_signatures": tuple(development_finite_signatures),
                "hermite_values": tuple(hermite_values),
                "candidates": candidate_records,
                "parent_bases": parent_bases,
                "graph_control": graph_control,
                "proposed_count": len(proposed_bases),
            }
        )
        print(
            f"LATENTJOINTCOMPONENTPROGRESS|label={label}|proposed={len(proposed_bases)}|"
            f"deduplicated={len(deduplicated)}",
            flush=True,
        )

    bundle_ledger = joint_component_bundle_ledger(
        tuple(state["fingerprints"] for state in states),
        tuple(state["hermite_values"] for state in states),
        tuple(state["development_finite_signatures"] for state in states),
        tuple(
            tuple(record["held_out_replayed_ray_count"] for record in state["candidates"])
            for state in states
        ),
        (16, 16, 16, 16),
        eligible_families=tuple(
            tuple(
                record["finite_rank_compatible"]
                and record["passes_height_angle_filter"]
                and all(value == 1 for value in record["smith_invariant_factors"])
                for record in state["candidates"]
            )
            for state in states
        ),
        structural_quantile=STRUCTURAL_BUNDLE_QUANTILE,
    )
    selected_bundle = bundle_ledger.selected

    fibre_results = []
    exact_component_count = 0
    exact_completion_count = 0
    for fibre_index, (state, candidate_index) in enumerate(
        zip(states, selected_bundle.candidate_indices)
    ):
        candidate = state["candidates"][candidate_index]
        completions = candidate["rank17_completion_candidates"]
        if not completions:
            raise ArithmeticError("selected rank-16 component has no one-rank completion")
        completion = completions[0]
        truth = truth_by_label[state["label"]]
        published = tuple(tuple(map(int, row)) for row in truth["embedding_matrix_columns"])
        source_image = embedded_rows(source_basis, published)
        selected_basis = tuple(tuple(map(int, row)) for row in candidate["primitive_embedding_matrix_rows"])
        completion_basis = state["parent_bases"][completion["source_index"]]
        component_exact = rational_rank((*selected_basis, *source_image)) == 16
        completion_exact = rational_rank((*completion_basis, *published)) == 17
        exact_component_count += int(component_exact)
        exact_completion_count += int(completion_exact)
        parent_change = row_basis_coordinates(completion_basis, published) if completion_exact else ()
        fibre_results.append(
            {
                "label": state["label"],
                "proposed_component_occurrence_count": state["proposed_count"],
                "deduplicated_component_count": len(state["candidates"]),
                "selected_candidate_index": candidate_index,
                "selected_rank16_embedding_matrix_rows": candidate["primitive_embedding_matrix_rows"],
                "selected_rank16_held_out_replayed_ray_count": candidate["held_out_replayed_ray_count"],
                "selected_rank16_full_replayed_ray_count": candidate["full_replayed_ray_count"],
                "selected_rank17_completion_source_index": completion["source_index"],
                "selected_rank17_completion_embedding_matrix_rows": [list(row) for row in completion_basis],
                "completion_held_out_rays_gained_per_added_rank": completion["held_out_rays_gained_per_added_rank"],
                "completion_basis_change_to_published_rows": [list(map(int, row)) for row in parent_change],
                "postselection_rank16_equals_source_component_image": component_exact,
                "postselection_rank17_equals_published_r17_space": completion_exact,
                "candidate_components": state["candidates"],
            }
        )

    passed = exact_component_count == len(LABELS) and exact_completion_count == len(LABELS)
    status = (
        "PASS_JOINT_CENTER_FREE_R17_RANK16_THEN_RANK17"
        if passed
        else "FAIL_JOINT_COMPONENT_CALIBRATION_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-joint-component-calibration.v1",
        "status": status,
        "scope": "Phase-0 R17 controls only; no wgxli target is loaded",
        "bounds": {
            "source_sample_count": SOURCE_SAMPLE_COUNT,
            "target_sample_count_per_parent": TARGET_SAMPLE_COUNT,
            "target_components_per_parent": TARGET_COMPONENTS_PER_PARENT,
            "parent_pool_size": PARENT_POOL_SIZE,
            "random_seed": RANDOM_SEED,
            "global_development_fraction": f"{DEVELOPMENT_NUMERATOR}/{DEVELOPMENT_DENOMINATOR}",
            "adaptive_parent_proposal_fraction": f"{DEVELOPMENT_NUMERATOR}/{DEVELOPMENT_DENOMINATOR}",
            "finite_rank_primes": [2, 3],
            "maximum_source_height_angle_profile_distance": f"{MAXIMUM_SOURCE_ANGLE_DISTANCE:.17g}",
            "structural_bundle_quantile": f"{STRUCTURAL_BUNDLE_QUANTILE:.17g}",
            "structural_channel_weights": {
                "relation_height_component": "1/3",
                "primitive_hermite": "1/3",
                "development_finite_code": "1/3",
            },
        },
        "center_free_source_component": {
            "rank": source_component.rational_rank,
            "ray_count": len(source_component.vertex_indices),
            "relation_count": len(source_component.relation_indices),
            "primitive_basis_rows": [list(row) for row in source_basis],
            "oracle_center_used": False,
        },
        "joint_selection": {
            "candidate_indices": list(selected_bundle.candidate_indices),
            "mean_combined_distance_percentile": f"{selected_bundle.mean_combined_distance_percentile:.17g}",
            "maximum_combined_distance_percentile": f"{selected_bundle.maximum_combined_distance_percentile:.17g}",
            "relation_distance_percentile": f"{selected_bundle.relation_distance_percentile:.17g}",
            "hermite_distance_percentile": f"{selected_bundle.hermite_distance_percentile:.17g}",
            "finite_distance_percentile": f"{selected_bundle.finite_distance_percentile:.17g}",
            "held_out_replayed_ray_count": selected_bundle.held_out_replayed_ray_count,
            "held_out_rays_per_rank": f"{selected_bundle.held_out_rays_per_rank:.17g}",
            "structurally_retained_bundle_count": bundle_ledger.structurally_retained_bundle_count,
            "eligible_bundle_count": bundle_ledger.eligible_bundle_count,
            "all_bundle_count": bundle_ledger.generated_bundle_count,
        },
        "controls": fibre_results,
        "postselection_truth_audit": {
            "exact_rank16_component_count": exact_component_count,
            "exact_rank17_completion_count": exact_completion_count,
        },
        "gate_decision": (
            "OPEN FOR FERMIGIER COMPONENT CALIBRATION; keep wgxli closed."
            if passed
            else "CLOSED. The fixed joint component selector did not recover every R17 control."
        ),
        "proof_boundary": (
            "Rational-space keys, global replay memberships, induced relations, F_2/F_3 ranks, "
            "finite quotient restrictions, Smith factors, and postselection space comparisons are exact. "
            "Dense-hyperplane enumeration is bounded to the declared deterministic samples. Height profiles, "
            "joint fingerprint distances, the parent shortlist, and all selection thresholds are heuristic. "
            "The common source complex uses supervised R17 coefficient alignments, but no vertex center or "
            "published target embedding is supplied during proposal or selection."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CONSENSUS,
                TRUTH,
                GRAPH,
                FINITE,
                PROPOSALS,
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice joint-component artifact is stale")
        print(f"LATENTJOINTCOMPONENT|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTJOINTCOMPONENT|status={status}|rank16={exact_component_count}/4|"
        f"rank17={exact_completion_count}/4|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
