#!/usr/bin/env sage -python
"""Freeze the half-lattice replacement detector for the existing CRT cohort.

This program reads no outcome from the 2,560 prospective fibres.  It fixes the
native-074d9 coordinates of the 43 deepest generic R17 parity classes, checks
the already-frozen +12 positive-control ablation, and commits the two-stage
escalation rule and every finite search/certification limit.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
from typing import Any

from sage.all import Matrix, ZZ


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
LINEAGE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
OLD_SENSITIVITY = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-search-sensitivity-v1.json"
ABLATION_BLIND = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_search_ablation_rank29_holdout_blind_v1.json"
ABLATION_VERIFIED = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_search_ablation_rank29_holdout_verification_v1.json"
ABLATION_SUMMARY = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_search_ablation_summary_v1.json"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
RUNNER = ROOT / "elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage"
ANALYZER = ROOT / "elkies-k3/scripts/analyze_r17_prospective_crt_half_lattice_experiment.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-protocol-v3.json"

EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
EXPECTED_MANIFEST_SHA256 = "7e8c43a6f67eac96dd9dede333f94e0cce139fa685b421f83ad7e4d69c1a75d4"
EXPECTED_OLD_SENSITIVITY_SHA256 = "9787d6010c8384b7ce7f13915345b03cff30c87bdc7fea64b3c32861036a7a01"
EXPECTED_ABLATION_BLIND_SHA256 = "1ee832ce6ecebc0550c008f8a10ccc2d75e727dfe9d5625802624c160e7969e6"
EXPECTED_ABLATION_VERIFIED_SHA256 = "aae62d0a582aa55f0c16d7bd0fbd728e6a3f274fa64ce6b826e0873578b3f599"
EXPECTED_ABLATION_SUMMARY_SHA256 = "fbdfa24b14bc86ee33a576f5e3c3e894dd91dd5e0d1fbfb47bf208e167a7282a"

DIMENSION = 17
ARM_SIZE = 43
SPECIALIZED_SCALE = 1_000_000
SPECIALIZED_AUDIT_SCALE = 100_000
HEIGHT_BOUND = 100_000
PER_COVER_TIMEOUT_SECONDS = 15
GP_STACK_BYTES = 1_000_000_000
FIBRE_WORKER_TIMEOUT_SECONDS = 1_800
FIBRE_WORKER_ADDRESS_SPACE_BYTES = None
CERTIFICATE_PRIME_BOUND = 1_000
RETRIES = 0


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def build() -> dict[str, Any]:
    pinned_hashes = {
        MANIFEST: EXPECTED_MANIFEST_SHA256,
        OLD_SENSITIVITY: EXPECTED_OLD_SENSITIVITY_SHA256,
        ABLATION_BLIND: EXPECTED_ABLATION_BLIND_SHA256,
        ABLATION_VERIFIED: EXPECTED_ABLATION_VERIFIED_SHA256,
        ABLATION_SUMMARY: EXPECTED_ABLATION_SUMMARY_SHA256,
    }
    for path, expected in pinned_hashes.items():
        if digest(path) != expected:
            raise ArithmeticError(f"pinned protocol input changed: {relative(path)}")

    manifest = json.loads(MANIFEST.read_text())
    lineage = json.loads(LINEAGE.read_text())
    old_sensitivity = json.loads(OLD_SENSITIVITY.read_text())
    blind = json.loads(ABLATION_BLIND.read_text())
    verified = json.loads(ABLATION_VERIFIED.read_text())
    summary = json.loads(ABLATION_SUMMARY.read_text())
    if manifest.get("status") != "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS":
        raise ArithmeticError("the original CRT cohort is not frozen and unopened")
    if manifest["commitment"]["candidate_list_sha256"] != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the original CRT candidate commitment changed")
    if len(manifest["rows"]) != 2_560 or any(
        row["outcome_status"] != "NOT_OPENED" for row in manifest["rows"]
    ):
        raise ArithmeticError("the original manifest was modified or opened")
    if old_sensitivity.get("status") != "FAILED_TO_REDETECT_BOTH_KNOWN_PLUS12_POSITIVE_CONTROLS":
        raise ArithmeticError("the superseded detector sensitivity record changed")
    if blind.get("status") != "PASS_BLIND_ABLATION_SEARCH":
        raise ArithmeticError("the half-lattice holdout search is not complete")
    if verified.get("status") != "PASS_EXACT_PUBLIC_QUOTIENT_ABLATION":
        raise ArithmeticError("the half-lattice holdout verifier is not exact and passing")
    if summary.get("status") != "PASS_EQUAL_BUDGET_BLIND_ABLATION":
        raise ArithmeticError("the equal-budget ablation summary is not passing")
    if verified["phase_boundary"]["blind_artifact_sha256_before_fixture_import"] != digest(
        ABLATION_BLIND
    ):
        raise ArithmeticError("the +12 verifier did not hash the pinned blind bytes")

    positive = {}
    for case in verified["results"]:
        if case["label"] not in {"curve356-rank29", "curve385-rank29"}:
            continue
        arms = {arm["id"]: arm for arm in case["arms"]}
        positive[case["label"]] = {
            "stage_a_exact_quotient_rank": arms["generic-deepest43"][
                "exact_quotient_rank_over_Q"
            ],
            "stage_a_exact_quotient_rank_mod2": arms["generic-deepest43"][
                "exact_quotient_rank_mod2"
            ],
            "stage_a_class_count": arms["generic-deepest43"]["class_count"],
            "stage_b_union_exact_quotient_rank": arms["deep-union"][
                "exact_quotient_rank_over_Q"
            ],
            "stage_b_union_exact_quotient_rank_mod2": arms["deep-union"][
                "exact_quotient_rank_mod2"
            ],
            "stage_b_union_class_count": arms["deep-union"]["class_count"],
        }
    if positive != {
        "curve356-rank29": {
            "stage_a_exact_quotient_rank": 12,
            "stage_a_exact_quotient_rank_mod2": 12,
            "stage_a_class_count": 43,
            "stage_b_union_exact_quotient_rank": 12,
            "stage_b_union_exact_quotient_rank_mod2": 12,
            "stage_b_union_class_count": 75,
        },
        "curve385-rank29": {
            "stage_a_exact_quotient_rank": 3,
            "stage_a_exact_quotient_rank_mod2": 3,
            "stage_a_class_count": 43,
            "stage_b_union_exact_quotient_rank": 4,
            "stage_b_union_exact_quotient_rank_mod2": 4,
            "stage_b_union_class_count": 75,
        },
    }:
        raise ArithmeticError("the frozen +12 positive-control acceptance result changed")

    if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ArithmeticError("the exact 074d9 lineage input is not passing")
    gram = Matrix(ZZ, lineage["generic_basis"]["height_gram"])
    if gram.dimensions() != (DIMENSION, DIMENSION) or gram.det() != 948:
        raise ArithmeticError("the native 074d9 R17 height form changed")
    engine = SourceFileLoader("r17_crt_half_lattice_protocol_engine", str(ENGINE)).load_module()
    oracle = engine.CosetOracle(tuple(tuple(int(value) for value in row) for row in gram.rows()))
    rows = []
    histogram: Counter[int] = Counter()
    for mask in range(1 << DIMENSION):
        norm, representative, error = oracle.solve(mask)
        if error > 1.0e-6:
            raise ArithmeticError("generic CVP error exceeded its exact-recompute tolerance")
        rows.append((norm, mask, representative))
        histogram[norm] += 1
    rows.sort(key=lambda row: (-row[0], row[1]))
    maximum_norm = rows[0][0]
    deepest = rows[:ARM_SIZE]
    if maximum_norm != 12 or histogram[maximum_norm] != ARM_SIZE:
        raise ArithmeticError("the native R17 deepest stratum stopped having norm 12 and size 43")
    deepest_masks = [mask for unused_norm, mask, unused_rep in deepest]
    deepest_representatives = {
        str(mask): list(map(int, representative))
        for unused_norm, mask, representative in deepest
    }

    protocol = {
        "schema": "elkies-k3.r17-prospective-crt-half-lattice-protocol.v3",
        "status": "FROZEN_AFTER_POSITIVE_CONTROLS_BEFORE_NEW_COHORT_OUTCOMES",
        "protocol_id": "r17-prospective-crt-half-lattice-two-stage-v3",
        "candidate_list_sha256": EXPECTED_CANDIDATE_HASH,
        "cohort_selection": {
            "scheduled_candidate_count": 2_560,
            "original_manifest_reused_without_extension_or_rebalancing": True,
            "parameters_cohorts_match_sets_and_primary_contrasts_unchanged": True,
            "new_detector_outcomes_used_for_selection": False,
        },
        "freeze_boundary": {
            "new_2560_fibre_half_lattice_searches_completed_before_freeze": 0,
            "new_2560_fibre_half_lattice_points_seen_before_freeze": 0,
            "positive_controls_are_external_predeclared_rank29_fibres": True,
            "positive_control_public_points_loaded_only_after_blind_search_froze": True,
        },
        "superseded_detector": {
            "protocol": "direct completed-square x-search at height 10000",
            "sensitivity_status": old_sensitivity["status"],
            "curve356_certified_escape_count": 0,
            "curve385_certified_escape_count": 0,
            "interpretation": "DETECTOR_LIMITED_NOT_EVIDENCE_AGAINST_CRT_LOCAL_CONDITIONS",
        },
        "positive_control_acceptance": positive,
        "native_generic_lattice": {
            "basis": "ordered native polynomial sections 1,...,17 on norm12-orbit-074d9",
            "height_gram": [[int(value) for value in row] for row in gram.rows()],
            "height_gram_determinant": int(gram.det()),
            "height_gram_sha256": canonical_hash(
                [[int(value) for value in row] for row in gram.rows()]
            ),
            "exact_minimum_norm_histogram": {
                str(norm): count for norm, count in sorted(histogram.items())
            },
            "deepest_minimum_norm": maximum_norm,
            "deepest_half_lattice_depth": "3",
            "deepest_class_count": ARM_SIZE,
            "deepest_masks_in_norm_then_mask_order": deepest_masks,
            "deepest_masks_sha256": canonical_hash(deepest_masks),
            "generic_shortest_representatives_audit_only": deepest_representatives,
        },
        "specialized_representative_policy": {
            "canonical_height_engine": "PARI ellheightmatrix on exact short specialization",
            "real_precision_decimal_digits": 110,
            "operative_integer_rounding_scale": SPECIALIZED_SCALE,
            "audit_integer_rounding_scale": SPECIALIZED_AUDIT_SCALE,
            "coset_cvp_degree": 2,
            "representative_rule": (
                "for each mask, exact parity residue plus twice the closest lattice vector "
                "returned by the fplll dd oracle; recompute the rounded quadratic norm exactly"
            ),
        },
        "stage_a": {
            "selection": "the fixed 43 exact generic-deepest masks",
            "representative": "specialized shortest representative at operative scale 10^6",
            "covers_per_fibre": ARM_SIZE,
            "all_scheduled_fibres_receive_stage_a": True,
        },
        "stage_b": {
            "gate": "at least one Stage-A point has a full exact finite-reduction independence certificate beyond specialized MW17",
            "gate_is_frozen_before_any_new_cohort_outcome": True,
            "full_specialized_ranking_masks": 1 << DIMENSION,
            "specialized_top_count": ARM_SIZE,
            "operative_order": "descending actual decimal canonical depth of the scale-10^6 CVP representative, then mask",
            "stability_audit": "also compute scale-10^5 top-43 set and record equality without changing the operative set",
            "search_set": "generic top-43 union specialized top-43",
            "incremental_execution": "do not rerun Stage-A masks; search only specialized top-43 masks outside the generic top-43",
            "no_other_escalation_rule": True,
        },
        "cover_pipeline": {
            "specialization_normalization": (
                "exact Sage local_data(2).minimal_model(), first exact isomorphism, "
                "then the canonical integral short model [-27*c4,-54*c6]"
            ),
            "model": "w^2=m^4-6*x_P*m^2-8*y_P*m-3*x_P^2-4*A",
            "base_point": "exact specialized linear combination given by the selected representative",
            "denominator_clearing": "exact square clearing",
            "minimization": "PARI hyperellminimalmodel",
            "reduction": "PARI hyperellred",
            "search": "PARI hyperellratpoints on the reduced model",
            "retries": RETRIES,
            "height_bound_each_cover": HEIGHT_BOUND,
            "wall_timeout_seconds_each_cover_including_minimize_reduce_search": PER_COVER_TIMEOUT_SECONDS,
            "gp_stack_bytes_each_cover": GP_STACK_BYTES,
            "same_pipeline_and_limits_for_every_cover": True,
            "every_returned_point_mapped_back_by_exact_inverse_model_changes": True,
        },
        "fibre_worker_envelope": {
            "wall_timeout_seconds": FIBRE_WORKER_TIMEOUT_SECONDS,
            "address_space_bytes": FIBRE_WORKER_ADDRESS_SPACE_BYTES,
            "memory_control": "per-cover GP stack bound; no RLIMIT_AS because loaded numerical libraries reserve large virtual arenas",
            "retries": RETRIES,
            "worker_timeout_or_failure_is_censored_not_a_bounded_miss": True,
        },
        "point_acceptance": {
            "exact_original_short_curve_equation_required": True,
            "finite_reduction_prime_bound": CERTIFICATE_PRIME_BOUND,
            "combined_mod2_rank_must_equal_17_plus_all_counted_directions": True,
            "full_rank_certifies_nonmembership_and_Q_linear_independence": True,
            "uncertified_returned_points_do_not_count_as_escapes": True,
        },
        "primary_outcome": {
            "unit": "frozen fibre",
            "event": "at least one exactly certified Stage-A quotient direction",
            "primary_comparison": "pooled A_356_full plus B_385_full versus C_matched_ordinary",
            "primary_estimand": "certified Stage-A detector yield per scheduled fibre (intention to search)",
            "censoring_rule": (
                "an exact event still counts on a partially failed row; otherwise every scheduled "
                "row remains in the primary denominator, with complete-case rates reported as sensitivity"
            ),
            "secondary_fixed_cohort_order": [
                "F_random_equal_codimension",
                "D_two_only",
                "E_odd_only",
                "A_356_full+B_385_full",
            ],
            "stage_b_is_conditional_recovery_depth_not_an_unconditional_event_comparison": True,
        },
        "claim_boundary": [
            "The generic CVP census and every accepted point/reduction certificate are exact.",
            "Specialized canonical-height rankings are numerical and retain the declared scale-stability audit.",
            "All point searches are bounded; a miss is not rank 17 or a Selmer upper bound.",
            "The pointed quartics are search charts, not nontrivial locally soluble 2-covering torsors.",
            "The frozen experiment tests enrichment in detector-visible escapes, not a theorem of p-adic cylinder constancy.",
            "Stage B is missing-by-design on Stage-A-negative fibres and cannot be analyzed as an unconditional union-search response.",
        ],
    }
    return {
        **protocol,
        "protocol_definition_sha256": canonical_hash(protocol),
        "inputs": {
            relative(path): digest(path)
            for path in (
                MANIFEST,
                LINEAGE,
                OLD_SENSITIVITY,
                ABLATION_BLIND,
                ABLATION_VERIFIED,
                ABLATION_SUMMARY,
                ENGINE,
                RUNNER,
                ANALYZER,
            )
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_r17_prospective_crt_half_lattice_protocol.sage"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored half-lattice protocol differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17CRTHALFPROTOCOL"
        f"|hash={document['protocol_definition_sha256']}"
        "|status=FROZEN_BEFORE_NEW_COHORT_OUTCOMES",
        flush=True,
    )


if __name__ == "__main__":
    main()
