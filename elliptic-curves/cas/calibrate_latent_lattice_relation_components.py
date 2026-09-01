#!/usr/bin/env python3
"""Calibrate center-free dense-component selection on the rank-25 R17 fibre.

No wgxli record is loaded.  A dense rank-16 source component is proposed by
deterministic exact hyperplane sampling.  The same proposal is made inside 16
height-prefiltered rank-17 parent candidates in the held-out public subgroup.
Finite ranks and height-angle/height-shape consistency prune candidates before
the selector maximizes exact replayed held-out rays per added rank.  Published
R17 coordinates are read only for the final postselection audit.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC))

from latent_lattice import (  # noqa: E402
    biconnected_components,
    bounded_dense_hyperplanes,
    build_relation_complex,
    height_angle_profile,
    height_angle_profile_distance,
    modular_rank,
    primitive_hermite_signatures,
    rational_rank,
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
LEDGER = ARTIFACTS / "latent_lattice_joint_fingerprint_ledger_v1.json.gz"
OUTPUT = ARTIFACTS / "latent_lattice_relation_component_calibration_v1.json"
LABEL = "rank_at_least_25"
SOURCE_SAMPLE_COUNT = 2_000
TARGET_SAMPLE_COUNT = 400
TARGET_COMPONENTS_PER_PARENT = 4
RANDOM_SEED = 20_260_901
PARENT_POOL_SIZE = 16
DEVELOPMENT_FRACTION_NUMERATOR = 3
DEVELOPMENT_FRACTION_DENOMINATOR = 4
ANGLE_DISTANCE_SLACK = 0.0021
HERMITE_DISTANCE_SLACK = 0.007


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ambient_rows(component, parent_basis):
    raw = tuple(
        tuple(
            sum(
                int(component_vertex[row]) * int(parent_basis[row][column])
                for row in range(len(parent_basis))
            )
            for column in range(len(parent_basis[0]))
        )
        for component_vertex in component
    )
    return primitive_span_basis(independent_row_basis(raw))


def select_candidate(records, angle_slack, hermite_slack):
    finite = [record for record in records if record["finite_rank_compatible"]]
    if not finite:
        return None
    best_angle = min(float(record["height_angle_profile_distance"]) for record in finite)
    best_hermite = min(float(record["hermite_distance"]) for record in finite)
    survivors = [
        record
        for record in finite
        if float(record["height_angle_profile_distance"]) <= best_angle + angle_slack
        and float(record["hermite_distance"]) <= best_hermite + hermite_slack
    ]
    if not survivors:
        return None
    return max(
        survivors,
        key=lambda record: (
            float(record["held_out_rays_per_added_rank"]),
            record["full_replay_ray_count"],
            -float(record["hermite_distance"]),
            -float(record["height_angle_profile_distance"]),
            -record["parent_source_index"],
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    consensus = json.loads(CONSENSUS.read_text())
    truth_document = json.loads(TRUTH.read_text())
    graph = json.loads(GRAPH.read_text())
    with gzip.open(LEDGER, "rt") as stream:
        proposal_ledger = json.load(stream)

    held = next(
        record
        for record in consensus["leave_one_fibre_out"]
        if record["held_out_label"] == LABEL
    )
    control = next(
        record for record in consensus["controls"] if record["label"] == LABEL
    )
    source = build_relation_complex(held["training_two_of_three_vectors"])
    source_ledger = bounded_dense_hyperplanes(
        source,
        sample_count=SOURCE_SAMPLE_COUNT,
        random_seed=RANDOM_SEED,
        maximum_components=8,
    )
    if not source_ledger.components:
        raise ArithmeticError("source component sampler returned no hyperplane")
    source_component = source_ledger.components[0]
    source_basis = primitive_span_basis(
        independent_row_basis(
            tuple(source.vertices[index] for index in source_component.vertex_indices)
        )
    )
    source_rows = tuple(
        source.vertices[index] for index in source_component.vertex_indices
    )
    source_biconnected = biconnected_components(source)

    training_forms = []
    for record in truth_document["positive_controls"]:
        if record["label"] == LABEL:
            continue
        form = np.asarray(record["canonical_height_gram"], dtype=float)
        training_forms.append(form / np.trace(form))
    source_gram = sum(training_forms) / len(training_forms)
    source_angle_profile = height_angle_profile(source_rows, source_gram)
    source_hermite = primitive_hermite_signatures(
        source_gram, (source_basis,), digits=60
    )[0]

    fibre = next(
        record
        for record in proposal_ledger["fibres"]
        if record["label"] == LABEL
    )
    candidates_by_source = {
        int(record["source_index"]): record for record in fibre["candidates"]
    }
    graph_control = next(
        record for record in graph["default_controls"] if record["label"] == LABEL
    )
    parent_indices = tuple(
        int(record["source_index"])
        for record in graph_control["default_ledger"]["candidates"][:PARENT_POOL_SIZE]
    )
    rows = tuple(
        tuple(map(int, row.split()))
        for row in control["ambient_short_vector_coordinate_rows"]
    )
    heights = tuple(map(float, control["ambient_short_vector_heights"]))
    target_gram = fibre["ambient_height_gram"]

    records = []
    target_bases = []
    for parent_source_index in parent_indices:
        parent = candidates_by_source[parent_source_index]
        parent_basis = tuple(
            tuple(map(int, row)) for row in parent["basis_rows"]
        )
        parent_mask = exact_span_mask(rows, parent_basis)
        parent_ray_indices = tuple(
            index for index, retained in enumerate(parent_mask) if retained
        )
        height_order = sorted(
            parent_ray_indices, key=lambda index: (heights[index], index)
        )
        cutoff = (
            DEVELOPMENT_FRACTION_NUMERATOR * len(height_order)
        ) // DEVELOPMENT_FRACTION_DENOMINATOR
        development_indices = tuple(height_order[:cutoff])
        held_out_indices = frozenset(height_order[cutoff:])
        development_rows = tuple(rows[index] for index in development_indices)
        development_coordinates = row_basis_coordinates(
            development_rows, parent_basis
        )
        development_complex = build_relation_complex(development_coordinates)
        component_ledger = bounded_dense_hyperplanes(
            development_complex,
            sample_count=TARGET_SAMPLE_COUNT,
            random_seed=RANDOM_SEED + parent_source_index,
            maximum_components=TARGET_COMPONENTS_PER_PARENT,
        )
        if not component_ledger.components:
            continue
        for component_index, component in enumerate(component_ledger.components):
            development_component_rows = tuple(
                development_complex.vertices[index]
                for index in component.vertex_indices
            )
            target_basis = ambient_rows(development_component_rows, parent_basis)
            replay_mask = exact_span_mask(rows, target_basis)
            held_out_count = int(
                sum(replay_mask[index] for index in held_out_indices)
            )
            full_count = int(sum(replay_mask))
            outside_parent_count = int(
                sum(
                    retained and not parent_retained
                    for retained, parent_retained in zip(replay_mask, parent_mask)
                )
            )
            target_mod2 = modular_rank(target_basis, 2)
            target_mod3 = modular_rank(target_basis, 3)
            angle_profile = height_angle_profile(
                tuple(
                    tuple(
                        sum(
                            int(vector[row]) * int(parent_basis[row][column])
                            for row in range(len(parent_basis))
                        )
                        for column in range(len(parent_basis[0]))
                    )
                    for vector in development_component_rows
                ),
                target_gram,
            )
            record = {
                "parent_source_index": parent_source_index,
                "component_index_within_parent": component_index,
                "parent_primitive_embedding_matrix_rows": [
                    list(row) for row in parent_basis
                ],
                "parent_ray_count": len(parent_ray_indices),
                "development_ray_count": len(development_indices),
                "held_out_ray_count": len(held_out_indices),
                "component_development_ray_count": len(component.vertex_indices),
                "component_development_relation_count": len(component.relation_indices),
                "component_embedding_matrix_rows": [list(row) for row in target_basis],
                "component_smith_invariant_factors": list(
                    row_embedding_smith_invariant_factors(target_basis)
                ),
                "component_mod2_rank": target_mod2,
                "component_mod3_rank": target_mod3,
                "finite_rank_compatible": (
                    target_mod2 == source_component.mod2_rank
                    and target_mod3 == source_component.mod3_rank
                ),
                "height_angle_profile_distance": f"{height_angle_profile_distance(source_angle_profile, angle_profile):.17g}",
                "held_out_replayed_ray_count": held_out_count,
                "held_out_rays_per_added_rank": f"{held_out_count / 16:.17g}",
                "full_replay_ray_count": full_count,
                "outside_parent_replay_ray_count": outside_parent_count,
                "held_out_replayed_ray_indices": [
                    index for index in sorted(held_out_indices) if replay_mask[index]
                ],
                "component_sampler": {
                    "sample_count": component_ledger.sample_count,
                    "independent_sample_count": component_ledger.independent_sample_count,
                    "distinct_hyperplane_count": component_ledger.distinct_hyperplane_count,
                },
            }
            records.append(record)
            target_bases.append(target_basis)

    target_hermites = primitive_hermite_signatures(
        target_gram,
        target_bases,
        digits=60,
        batch_size=PARENT_POOL_SIZE * TARGET_COMPONENTS_PER_PARENT,
    )
    source_hermite_value = float(source_hermite.hermite.log_hermite_invariant)
    for record, signature in zip(records, target_hermites):
        value = float(signature.hermite.log_hermite_invariant)
        record["primitive_hermite_signature"] = signature.to_record()
        record["hermite_distance"] = f"{abs(value - source_hermite_value):.17g}"

    selected = select_candidate(
        records, ANGLE_DISTANCE_SLACK, HERMITE_DISTANCE_SLACK
    )
    if selected is None:
        raise ArithmeticError("component filters returned no candidate")
    best_angle = min(float(record["height_angle_profile_distance"]) for record in records if record["finite_rank_compatible"])
    best_hermite = min(float(record["hermite_distance"]) for record in records if record["finite_rank_compatible"])
    for record in records:
        record["passes_angle_filter"] = (
            record["finite_rank_compatible"]
            and float(record["height_angle_profile_distance"])
            <= best_angle + ANGLE_DISTANCE_SLACK
        )
        record["passes_hermite_filter"] = (
            record["finite_rank_compatible"]
            and float(record["hermite_distance"])
            <= best_hermite + HERMITE_DISTANCE_SLACK
        )
        record["eligible_for_replay_optimization"] = (
            record["passes_angle_filter"] and record["passes_hermite_filter"]
        )

    truth = next(
        record
        for record in truth_document["positive_controls"]
        if record["label"] == LABEL
    )
    published_embedding = tuple(
        tuple(map(int, row)) for row in truth["embedding_matrix_columns"]
    )
    selected_component_basis = tuple(
        tuple(map(int, row)) for row in selected["component_embedding_matrix_rows"]
    )
    selected_parent_basis = tuple(
        tuple(map(int, row))
        for row in selected["parent_primitive_embedding_matrix_rows"]
    )
    source_component_image = tuple(
        tuple(
            sum(
                int(source_basis[row][middle])
                * int(published_embedding[middle][column])
                for middle in range(17)
            )
            for column in range(len(published_embedding[0]))
        )
        for row in range(16)
    )
    component_matches_truth = (
        rational_rank(selected_component_basis + source_component_image) == 16
    )
    parent_matches_truth = (
        rational_rank(selected_parent_basis + published_embedding) == 17
    )
    parent_change = (
        row_basis_coordinates(selected_parent_basis, published_embedding)
        if parent_matches_truth
        else ()
    )
    parent_change_smith = (
        row_embedding_smith_invariant_factors(parent_change)
        if parent_change
        else ()
    )
    selected_component_primitive = all(
        value == 1
        for value in row_embedding_smith_invariant_factors(selected_component_basis)
    )
    completion_primitive = all(
        value == 1
        for value in row_embedding_smith_invariant_factors(selected_parent_basis)
    )

    stability = []
    for angle_slack in (0.0021, 0.0023, 0.0025):
        for hermite_slack in (0.0065, 0.007, 0.008):
            choice = select_candidate(records, angle_slack, hermite_slack)
            stability.append(
                {
                    "angle_distance_slack": f"{angle_slack:.17g}",
                    "hermite_distance_slack": f"{hermite_slack:.17g}",
                    "selected_parent_source_index": (
                        choice["parent_source_index"] if choice else None
                    ),
                }
            )
    stability_failure_count = sum(
        record["selected_parent_source_index"]
        != graph_control["truth_source_index"]
        for record in stability
    )
    passed = (
        source_component.rational_rank == 16
        and len(source_component.vertex_indices) == 362
        and component_matches_truth
        and parent_matches_truth
        and selected_component_primitive
        and completion_primitive
        and selected["parent_source_index"] == graph_control["truth_source_index"]
    )
    status = (
        "PASS_CENTER_FREE_R17_RANK16_COMPONENT_AND_RANK17_COMPLETION"
        if passed
        else "FAIL_RELATION_COMPONENT_CALIBRATION_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-relation-component-calibration.v1",
        "status": status,
        "scope": "Phase-0 rank-25 R17 control only; no wgxli target is loaded",
        "bounds": {
            "source_sample_count": SOURCE_SAMPLE_COUNT,
            "target_sample_count_per_parent": TARGET_SAMPLE_COUNT,
            "target_components_per_parent": TARGET_COMPONENTS_PER_PARENT,
            "random_seed": RANDOM_SEED,
            "height_prefiltered_parent_count": PARENT_POOL_SIZE,
            "development_fraction": f"{DEVELOPMENT_FRACTION_NUMERATOR}/{DEVELOPMENT_FRACTION_DENOMINATOR}",
            "angle_distance_slack": f"{ANGLE_DISTANCE_SLACK:.17g}",
            "hermite_distance_slack": f"{HERMITE_DISTANCE_SLACK:.17g}",
            "finite_rank_primes": [2, 3],
        },
        "center_free_source_component": {
            "rank": source_component.rational_rank,
            "ray_count": len(source_component.vertex_indices),
            "relation_count": len(source_component.relation_indices),
            "mod2_rank": source_component.mod2_rank,
            "mod3_rank": source_component.mod3_rank,
            "primitive_basis_rows": [list(row) for row in source_basis],
            "leading_support_counts": [
                len(component.vertex_indices)
                for component in source_ledger.components
            ],
            "sample_ledger": {
                "independent_sample_count": source_ledger.independent_sample_count,
                "distinct_hyperplane_count": source_ledger.distinct_hyperplane_count,
            },
            "biconnected_component_sizes": [
                len(component.vertex_indices) for component in source_biconnected[:16]
            ],
            "primitive_hermite_signature": source_hermite.to_record(),
            "oracle_center_used": False,
        },
        "candidate_components": records,
        "selection": {
            "selected_parent_source_index": selected["parent_source_index"],
            "selected_rank16_embedding_matrix_rows": selected[
                "component_embedding_matrix_rows"
            ],
            "selected_rank16_smith_invariant_factors": selected[
                "component_smith_invariant_factors"
            ],
            "selected_rank16_held_out_replayed_ray_count": selected[
                "held_out_replayed_ray_count"
            ],
            "selected_rank16_full_replay_ray_count": selected[
                "full_replay_ray_count"
            ],
            "selected_rank17_completion_embedding_matrix_rows": selected[
                "parent_primitive_embedding_matrix_rows"
            ],
            "rank17_completion_smith_invariant_factors": list(
                row_embedding_smith_invariant_factors(selected_parent_basis)
            ),
            "rank17_completion_basis_change_to_published_rows": [
                [int(value) for value in row] for row in parent_change
            ],
            "rank17_completion_basis_change_smith_invariant_factors": list(
                parent_change_smith
            ),
        },
        "postselection_truth_audit": {
            "withheld_truth_parent_source_index": graph_control["truth_source_index"],
            "selected_rank16_equals_image_of_blind_source_component": component_matches_truth,
            "selected_rank17_equals_published_r17_rational_space": parent_matches_truth,
            "selected_rank16_primitive_in_public_subgroup": selected_component_primitive,
            "selected_rank17_primitive_in_public_subgroup": completion_primitive,
        },
        "stability_box": {
            "configurations": stability,
            "failure_count": stability_failure_count,
        },
        "gate_decision": (
            "OPEN FOR THE NEXT CONTROL. A center-free dense relation component "
            "selects the exact primitive R17 rank-16 image by finite/height "
            "pruning followed by held-out replay per rank, and its parent gives "
            "the exact primitive rank-17 completion. Do not apply to wgxli until "
            "the same component selector is checked on the other R17 fibres."
            if passed
            else "CLOSED. The declared relation-component box did not recover the R17 control."
        ),
        "proof_boundary": (
            "Ranks over Q, F_2, and F_3; sampled hyperplane supports; induced "
            "relations; rational-space replay; Smith factors; held-out ray "
            "membership; and the final basis changes are exact. Hyperplane "
            "exhaustion is only over the declared deterministic samples. Canonical "
            "heights, angle profiles, Hermite filters, parent prefiltering, and "
            "thresholds are numerical or heuristic. The common source coordinate "
            "complex was built from supervised R17 control alignments. Published "
            "coordinates and the truth source index are used only after selection."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CONSENSUS,
                TRUTH,
                GRAPH,
                LEDGER,
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
            raise SystemExit("latent-lattice relation-component artifact is stale")
        print(
            "LATENTCOMPONENT|check=PASS|"
            f"sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTCOMPONENT|status={status}|selected_parent="
        f"{selected['parent_source_index']}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
