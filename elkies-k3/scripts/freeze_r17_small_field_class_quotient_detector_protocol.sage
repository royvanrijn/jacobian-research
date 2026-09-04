#!/usr/bin/env sage-python
"""Freeze the uniform detector only after every pre-search Q_t is certified."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
from typing import Any

from sage.all import Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
COHORT = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-cohort-v1.json"
FEATURES = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-features-v1.json"
TARGET = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
PINNED_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
BASE_RUNNER = ROOT / "elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage"
RUNNER = ROOT / "elkies-k3/scripts/run_r17_small_field_class_quotient_detector.sage"
ANALYZER = ROOT / "elkies-k3/scripts/analyze_r17_small_field_class_quotient_experiment.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-detector-protocol-v1.json"

FEATURE_STATUS = "FROZEN_COMPLETE_UNCONDITIONAL_PRE_SEARCH_FEATURES"
STATUS = "FROZEN_AFTER_ALL_Q_BEFORE_ANY_POINT_SEARCH"
DIMENSION = 17
TOP_COUNT = 43


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def public_height_gram(target):
    pinned = Matrix(ZZ, [
        [int(value) for value in line.split()]
        for line in PINNED_GRAM.read_text().splitlines()
        if line.strip()
    ])
    change = Matrix(ZZ, target["pinned_identification"]["basis_change_matrix"])
    if target["pinned_identification"]["gram_identity_orientation"] != "M^T*Gpub*M=Gpinned":
        raise ArithmeticError("the published-basis Gram orientation changed")
    public = change.transpose().inverse() * pinned * change.inverse()
    if public not in Matrix(ZZ, DIMENSION, DIMENSION).parent():
        raise ArithmeticError("the recovered published height Gram is not integral")
    public = Matrix(ZZ, public)
    if public.det() != 948 or not public.is_symmetric():
        raise ArithmeticError("the published height Gram lost the R17 invariants")
    return public


def build(features_path: Path):
    cohort = json.loads(COHORT.read_text())
    features = json.loads(features_path.read_text())
    target = json.loads(TARGET.read_text())
    if cohort.get("status") != "FROZEN_RANK_BLIND_PRE_CLASS_GROUP_COHORT":
        raise ArithmeticError("the cohort is not at the sealed pre-feature boundary")
    if any(row.get("outcome_status") != "SEALED_UNTIL_ALL_FEATURES_FREEZE" for row in cohort["rows"]):
        raise ArithmeticError("a point-search outcome was opened before protocol freeze")
    if features.get("status") != FEATURE_STATUS or not features.get("point_search_unlocked"):
        raise ArithmeticError("all unconditional Q_t features must freeze before detector protocol creation")
    if features.get("candidate_list_sha256") != cohort["commitment"]["candidate_list_sha256"]:
        raise ArithmeticError("the feature ledger names another cohort")
    if canonical_hash(features["records"]) != features.get("feature_commitment_sha256"):
        raise ArithmeticError("the feature commitment hash does not replay")
    if len(features["records"]) != len(cohort["rows"]):
        raise ArithmeticError("the feature ledger does not cover the cohort")
    for feature in features["records"]:
        if (
            feature.get("status") != "PASS_COMPLETE_UNCONDITIONAL_CLASS_QUOTIENT"
            or feature.get("dim_Q") is None
            or not feature.get("class_group", {}).get("bnfcertify_passed")
        ):
            raise ArithmeticError("an incomplete class quotient entered the detector freeze")

    gram = public_height_gram(target)
    engine = SourceFileLoader("small_field_detector_protocol_engine", str(ENGINE)).load_module()
    oracle = engine.CosetOracle(tuple(tuple(int(value) for value in row) for row in gram.rows()))
    census = []
    histogram = Counter()
    for mask in range(1 << DIMENSION):
        norm, representative, error = oracle.solve(mask)
        if error > 1.0e-6:
            raise ArithmeticError("generic R17 CVP error exceeded the exact-recompute tolerance")
        census.append((norm, mask, representative))
        histogram[norm] += 1
    census.sort(key=lambda row: (-row[0], row[1]))
    deepest_norm = census[0][0]
    deepest = [row for row in census if row[0] == deepest_norm]
    if deepest_norm != 12 or len(deepest) != TOP_COUNT:
        raise ArithmeticError("the published R17 basis lost its 43 deepest classes")
    deepest_masks = [mask for _norm, mask, _representative in deepest]

    protocol = {
        "schema": "elkies-k3.r17-small-field-class-quotient-detector-protocol.v1",
        "status": STATUS,
        "protocol_id": "r17-small-field-class-quotient-half-lattice-v1",
        "candidate_list_sha256": cohort["commitment"]["candidate_list_sha256"],
        "phase_boundary": {
            "feature_artifact": relative(features_path),
            "feature_artifact_sha256": digest(features_path),
            "feature_commitment_sha256": features["feature_commitment_sha256"],
            "all_dim_Q_values_frozen_before_protocol": True,
            "point_search_outcomes_seen_before_protocol": 0,
            "detector_loads_feature_values": False,
            "detector_may_hash_feature_file_without_parsing_it": True,
        },
        "detector_manifest": [
            {
                "sample_id": row["sample_id"],
                "manifest_index": index,
                "family": row["family"],
                "parameter": row["parameter"],
            }
            for index, row in enumerate(cohort["rows"])
        ],
        "native_generic_lattice": {
            "basis": "published ordered sections P1,...,P17",
            "height_gram": [[int(value) for value in row] for row in gram.rows()],
            "height_gram_determinant": int(gram.det()),
            "exact_minimum_norm_histogram": {
                str(norm): count for norm, count in sorted(histogram.items())
            },
            "deepest_minimum_norm": deepest_norm,
            "deepest_class_count": len(deepest),
            "deepest_masks_in_norm_then_mask_order": deepest_masks,
            "deepest_masks_sha256": canonical_hash(deepest_masks),
        },
        "specialized_representative_policy": {
            "canonical_height_engine": "PARI ellheightmatrix on exact specialized points",
            "real_precision_decimal_digits": 110,
            "operative_integer_rounding_scale": 1_000_000,
            "audit_integer_rounding_scale": 100_000,
            "representative_rule": (
                "exact parity residue plus twice the closest vector at the declared scale"
            ),
        },
        "stage_a": {
            "selection": "all 43 published-basis generic deepest masks",
            "representative": "specialized shortest representative at scale 10^6",
            "covers_per_fibre": TOP_COUNT,
            "all_scheduled_fibres_receive_stage_a": True,
        },
        "stage_b": {
            "gate": "at least one Stage-A point has a full exact finite-reduction independence certificate beyond specialized MW17",
            "full_specialized_ranking_masks": 1 << DIMENSION,
            "specialized_top_count": TOP_COUNT,
            "search_set": "generic top-43 union specialized top-43",
            "no_other_escalation_rule": True,
        },
        "cover_pipeline": {
            "specialization_normalization": (
                "exact Sage local_data(2).minimal_model(), first exact isomorphism, "
                "then canonical integral short model [-27*c4,-54*c6]"
            ),
            "model": "w^2=m^4-6*x_P*m^2-8*y_P*m-3*x_P^2-4*A",
            "minimization": "PARI hyperellminimalmodel then hyperellred",
            "search": "PARI hyperellratpoints on the reduced model",
            "height_bound_each_cover": 100_000,
            "wall_timeout_seconds_each_cover_including_minimize_reduce_search": 15,
            "gp_stack_bytes_each_cover": 512_000_000,
            "same_pipeline_and_limits_for_every_cover": True,
        },
        "fibre_worker_envelope": {
            "wall_timeout_seconds": 1800,
            "address_space_bytes": None,
            "worker_timeout_or_failure_is_censored_not_a_bounded_miss": True,
        },
        "point_acceptance": {
            "exact_original_short_curve_equation_required": True,
            "finite_reduction_prime_bound": 1000,
            "combined_mod2_rank_must_equal_17_plus_all_counted_directions": True,
            "uncertified_returned_points_do_not_count_as_escapes": True,
        },
        "predeclared_analysis": {
            "predictor": "dim_Q",
            "negative_control_predictor": "dim Cl(K)[2] = dim Cl(K)/2Cl(K)",
            "primary_outcome": "integer Stage-A exact certified quotient gain",
            "primary_direction": "larger dim_Q predicts larger gain",
            "primary_statistic": "Kendall tau-b",
            "primary_randomization_test": (
                "100000 deterministic permutations within cubic-signature by "
                "field-discriminant-quartile blocks"
            ),
            "secondary_outcomes": [
                "any exact Stage-A certified escape",
                "conditional Stage-B incremental exact gain",
            ],
            "confirmatory_censoring_rule": (
                "if any scheduled Stage-A row is censored, the confirmatory statistic "
                "and p-value are null; complete-case summaries are exploratory only"
            ),
            "no_threshold_training": True,
            "compare_quotient_against_total_class_group_without_model_selection": True,
        },
        "claim_boundary": [
            "Every counted direction receives an exact equation and finite-reduction independence certificate.",
            "The detector is bounded; a miss is not absence of Mordell-Weil escape or a rank upper bound.",
            "The primary experiment tests detector-visible certified gain, not the full unknown Mordell-Weil quotient.",
            "Stage B is conditional and is not an unconditional primary endpoint.",
            "Association does not by itself prove that Q_t is causal.",
        ],
    }
    return {
        **protocol,
        "protocol_definition_sha256": canonical_hash(protocol),
        "inputs": {
            relative(path): digest(path)
            for path in (
                COHORT,
                features_path,
                TARGET,
                PINNED_GRAM,
                ENGINE,
                BASE_RUNNER,
                RUNNER,
                ANALYZER,
            )
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                "sage -python elkies-k3/scripts/"
                "freeze_r17_small_field_class_quotient_detector_protocol.sage"
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=FEATURES)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build(args.features.resolve())
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored detector protocol differs from exact replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        f"R17SMALLFIELDPROTOCOL|rows={len(document['detector_manifest'])}"
        f"|hash={document['protocol_definition_sha256']}|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
