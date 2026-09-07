#!/usr/bin/env python3
"""Create the content-addressed no-tuning freeze for the first wgxli run.

This program never imports or reads an ICARM 351/356/376/377/385 record.  Its
output is the one-way methodological boundary: target observations may not be
used to change any field or hashed implementation under this tag.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
OUTPUT = ARTIFACTS / "latent_lattice_target_method_freeze_v1.json"
ALGORITHM_TAG = "LATENT-LATTICE-WGXLI-FROZEN-2026-09-01-v1"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    pinned_sources = (
        ELLIPTIC / "cas/calibrate_latent_lattice_method.py",
        ELLIPTIC / "cas/calibrate_latent_lattice_dimension_persistence.py",
        ELLIPTIC / "cas/calibrate_latent_lattice_graph_walk_consensus.py",
        ELLIPTIC / "cas/calibrate_latent_lattice_joint_components.py",
        ARTIFACTS / "latent_lattice_calibration_v2.json",
        ARTIFACTS / "latent_lattice_dimension_persistence_v1.json",
        ARTIFACTS / "latent_lattice_graph_walk_calibration_v1.json",
        ARTIFACTS / "latent_lattice_joint_component_calibration_v1.json",
        ARTIFACTS / "latent_lattice_finite_calibration_v1.json",
        *library_sources,
    )
    payload = {
        "schema": "elliptic-curves.latent-lattice-target-method-freeze.v1",
        "algorithm_tag": ALGORITHM_TAG,
        "status": "FROZEN_EXPERIMENTAL_TARGET_METHOD_NO_TUNING",
        "created_date": "2026-09-01",
        "target_records": [351, 356, 376, 377, 385],
        "target_data_read_before_freeze": False,
        "no_tuning_rule": (
            "After this manifest is written, no target observation may change a bound, seed, "
            "threshold, channel weight, score ordering, dimension rule, parent/component pool, "
            "finite-prime policy, or hold-out rule under this algorithm tag. A changed method "
            "requires a new tag and remains invalid as confirmation on these five records."
        ),
        "calibration_boundary": {
            "r17_rank17_graph_selector": "passes 4/4 positive-control fibres",
            "fermigier_dimension_persistence": "selects dimension 12 at both declared cutoffs",
            "center_free_rank16_component_proposal_recall": "exact withheld component occurs uniquely in 4/4 ledgers",
            "equal_weight_rank16_component_selection": "recovers the particular withheld component in 2/4 controls",
            "one_rank_completion_after_component_selection": "recovers R17 in 3/4 controls",
            "interpretation": (
                "The target run is exploratory and fail-closed. Component-selection calibration is "
                "not a theorem and is not strong enough for a provenance claim by itself."
            ),
        },
        "cloud_protocol": {
            "population": "full displayed independent subgroup; primitive unoriented short rays",
            "height_candidate_bounds": list(range(20, 82, 2)),
            "adaptive_bound_rule": (
                "choose the first bound with at least 1800 rays; reject the fibre if none occurs "
                "before bound 80 or if enumeration reaches 100000 rays"
            ),
            "persistence_cutoffs_by_height_order": ["3/4", "7/8", "1"],
            "canonical_height_digits": 80,
        },
        "dimension_protocol": {
            "candidate_dimensions": list(range(10, 21)),
            "beam_subspace_scan": {
                "pool": 300,
                "beam_width": 8,
                "branch_width": 80,
                "seed_rule": "10000 + ICARM record number",
            },
            "score": (
                "for k, take the minimum integrality-LLR over k-1,k,k+1 at every "
                "persistence cutoff; maximize this cross-cutoff bottleneck, then the k score, "
                "then prefer smaller k"
            ),
            "rank17_assumed": False,
        },
        "finite_protocol": {
            "relation_primes": [2, 3],
            "good_reduction_prime_bound": 251,
            "block_rule": (
                "for each relation prime, first three usable one-dimensional quotient blocks "
                "are development and next three are validation; reduction-prime identities and "
                "quotient bases are forgotten in cross-fibre signatures"
            ),
            "validation_blocks_used_in_selection": False,
            "primitive_embedding_required": True,
            "coordinate_modular_rank_prunes": [2, 3],
        },
        "component_protocol": {
            "parent_pool_size": 16,
            "adaptive_parent_proposal_fraction": "3/4",
            "global_replay_development_fraction": "3/4",
            "dense_hyperplane_samples_per_parent": 400,
            "components_retained_per_parent": 8,
            "random_seed": 20260901,
            "exact_rational_space_deduplication": True,
            "exact_saturation_after_proposal_or_merge": True,
            "maximum_source_height_angle_profile_distance": "0.029999999999999999",
            "structural_channels": [
                "abstract relation/height component fingerprint",
                "primitive Hermite signature after exact saturation",
                "development finite-code signature",
            ],
            "channel_weights": ["1/3", "1/3", "1/3"],
            "distance_normalization": "empirical lower-is-better percentile separately for every fibre pair/channel",
            "structural_pruning_quantile": "0.10000000000000001",
            "post_prune_objective": "maximize exact held-out replayed rays per total candidate rank",
        },
        "cross_validation_protocol": {
            "folds": "all five choices of one held-out fibre",
            "training_fibre_count": 4,
            "dimension_acceptance": "one k must recur in at least four independent fibre scans",
            "training_acceptance": "a compatible primitive component bundle must occur in all four training fibres",
            "held_out_prediction": (
                "select the held-out component using training structural signatures and development "
                "finite blocks only; evaluate global held-out rays and validation finite blocks afterward"
            ),
            "strong_signal": (
                "the same k and compatible component embeddings pass at least four folds, including "
                "successful prediction in each corresponding untouched fibre"
            ),
        },
        "forbidden_under_tag": [
            "displayed-label sign/permutation search",
            "unrestricted GL(k,Z)",
            "forcing k=17",
            "R17/Q80 provenance inference from shells",
            "equation or first-jet interpolation",
            "changing a parameter after inspecting any target outcome",
        ],
        "frozen_inputs": {
            str(path.relative_to(ROOT)): digest(path) for path in pinned_sources
        },
        "software": {"python": platform.python_version()},
        "proof_boundary": (
            "The manifest freezes a bounded experimental algorithm, not a lattice-existence theorem. "
            "Exact arithmetic subroutines remain exact within their declared clouds; height selection, "
            "dimension persistence, component sampling, structural distances, and cutoffs are heuristic."
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice target-method freeze is stale")
        print(f"LATENTFREEZE|check=PASS|tag={ALGORITHM_TAG}|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"LATENTFREEZE|status=FROZEN|tag={ALGORITHM_TAG}|output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
