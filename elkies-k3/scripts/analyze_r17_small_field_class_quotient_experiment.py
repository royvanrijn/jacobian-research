#!/usr/bin/env python3
"""Join the frozen Q_t features to blind detector outcomes exactly once."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import statistics
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
FEATURES = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-features-v1.json"
PROTOCOL = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-detector-protocol-v1.json"
OUTCOMES = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-detector-ledger-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-analysis-v1.json"

PERMUTATIONS = 100_000
SEED = 20_260_904


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def kendall_counts(left: Sequence[int], right: Sequence[int]):
    concordant = discordant = left_only_ties = right_only_ties = 0
    for first in range(len(left)):
        for second in range(first + 1, len(left)):
            delta_left = left[first] - left[second]
            delta_right = right[first] - right[second]
            if delta_left == 0 and delta_right == 0:
                continue
            if delta_left == 0:
                left_only_ties += 1
            elif delta_right == 0:
                right_only_ties += 1
            elif delta_left * delta_right > 0:
                concordant += 1
            else:
                discordant += 1
    return concordant, discordant, left_only_ties, right_only_ties


def kendall_tau_b(left: Sequence[int], right: Sequence[int]):
    concordant, discordant, left_ties, right_ties = kendall_counts(left, right)
    denominator = math.sqrt(
        (concordant + discordant + left_ties)
        * (concordant + discordant + right_ties)
    )
    return (
        (concordant - discordant) / denominator if denominator else None,
        {
            "concordant_pairs": concordant,
            "discordant_pairs": discordant,
            "dim_Q_tie_only_pairs": left_ties,
            "gain_tie_only_pairs": right_ties,
        },
    )


def quartile_blocks(rows):
    by_signature = defaultdict(list)
    for index, row in enumerate(rows):
        signature = tuple(row["field_signature"])
        by_signature[signature].append((abs(int(row["field_discriminant"])), index))
    blocks = [None] * len(rows)
    for signature, members in by_signature.items():
        members.sort()
        count = len(members)
        for rank, (_discriminant, index) in enumerate(members):
            quartile = min(3, (4 * rank) // count)
            blocks[index] = f"signature-{signature[0]}-{signature[1]}-quartile-{quartile + 1}"
    return blocks


def stratified_randomization(left, right, blocks, permutations=PERMUTATIONS, seed=SEED):
    observed_counts = kendall_counts(left, right)
    observed_score = observed_counts[0] - observed_counts[1]
    block_indices = defaultdict(list)
    for index, block in enumerate(blocks):
        block_indices[block].append(index)
    rng = random.Random(seed)
    at_least = 0
    permuted = list(right)
    for _iteration in range(permutations):
        for indices in block_indices.values():
            values = [right[index] for index in indices]
            rng.shuffle(values)
            for index, value in zip(indices, values):
                permuted[index] = value
        counts = kendall_counts(left, permuted)
        if counts[0] - counts[1] >= observed_score:
            at_least += 1
    return {
        "alternative": "positive association",
        "permutations": permutations,
        "seed": seed,
        "observed_concordant_minus_discordant": observed_score,
        "permuted_scores_at_least_observed": at_least,
        "p_value_plus_one": (at_least + 1) / (permutations + 1),
        "blocks": dict(sorted(Counter(blocks).items())),
    }


def grouped_summary(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["dim_Q"]].append(row["stage_a_certified_gain"])
    return {
        str(dimension): {
            "row_count": len(values),
            "escape_count": sum(value > 0 for value in values),
            "escape_rate": sum(value > 0 for value in values) / len(values),
            "certified_gain_sum": sum(values),
            "certified_gain_mean": statistics.fmean(values),
            "certified_gain_median": statistics.median(values),
        }
        for dimension, values in sorted(groups.items())
    }


def load_and_join():
    features = json.loads(FEATURES.read_text())
    protocol = json.loads(PROTOCOL.read_text())
    outcomes = json.loads(OUTCOMES.read_text())
    if features.get("status") != "FROZEN_COMPLETE_UNCONDITIONAL_PRE_SEARCH_FEATURES":
        raise ArithmeticError("the pre-search feature ledger is not complete and frozen")
    if protocol.get("status") != "FROZEN_AFTER_ALL_Q_BEFORE_ANY_POINT_SEARCH":
        raise ArithmeticError("the detector protocol is not frozen at the phase boundary")
    if outcomes.get("status") != "COMPLETE_FROZEN_SMALL_FIELD_DETECTOR_LEDGER":
        raise ArithmeticError("the detector ledger is incomplete")
    if digest(FEATURES) != protocol["phase_boundary"]["feature_artifact_sha256"]:
        raise ArithmeticError("the analyzed feature bytes differ from the protocol commitment")
    if protocol["protocol_definition_sha256"] != outcomes["protocol_definition_sha256"]:
        raise ArithmeticError("the outcomes were produced under another protocol")
    if protocol.get("inputs", {}).get(relative(Path(__file__))) != digest(Path(__file__)):
        raise ArithmeticError("this analyzer differs from the source frozen before point search")
    features_by_id = {row["sample_id"]: row for row in features["records"]}
    outcomes_by_id = {row["sample_id"]: row for row in outcomes["records"]}
    expected_ids = [row["sample_id"] for row in protocol["detector_manifest"]]
    if set(features_by_id) != set(expected_ids) or set(outcomes_by_id) != set(expected_ids):
        raise ArithmeticError("the feature/outcome join does not cover the detector manifest")
    rows = []
    for sample_id in expected_ids:
        feature = features_by_id[sample_id]
        outcome = outcomes_by_id[sample_id]
        complete = bool(outcome.get("analysis_eligible_complete_stage_a", False))
        rows.append(
            {
                "sample_id": sample_id,
                "family": feature["family"],
                "parameter": feature["parameter"],
                "dim_Q": int(feature["dim_Q"]),
                "dim_Cl_mod_2Cl": int(feature["class_group"]["dim_Cl_mod_2Cl"]),
                "field_signature": feature["cubic_field"]["signature"],
                "field_discriminant": feature["cubic_field"]["field_discriminant"],
                "detector_status": outcome["status"],
                "complete_stage_a": complete,
                "stage_a_certified_gain": (
                    int(outcome["stage_a"]["certified_quotient_gain"])
                    if complete
                    else None
                ),
                "stage_b_incremental_certified_gain": (
                    int(outcome.get("stage_b", {}).get("incremental_certified_quotient_gain", 0))
                    if complete
                    else None
                ),
            }
        )
    return features, protocol, outcomes, rows


def analysis_for_rows(rows, *, confirmatory):
    dimensions = [row["dim_Q"] for row in rows]
    total_dimensions = [row["dim_Cl_mod_2Cl"] for row in rows]
    gains = [row["stage_a_certified_gain"] for row in rows]
    tau, pair_counts = kendall_tau_b(dimensions, gains)
    total_tau, total_pair_counts = kendall_tau_b(total_dimensions, gains)
    blocks = quartile_blocks(rows)
    return {
        "role": "confirmatory" if confirmatory else "exploratory_complete_case",
        "row_count": len(rows),
        "kendall_tau_b": tau,
        "kendall_pair_counts": pair_counts,
        "stratified_randomization": stratified_randomization(dimensions, gains, blocks),
        "predeclared_total_class_group_negative_control": {
            "predictor": "dim Cl(K)[2] (= dim Cl(K)/2Cl(K))",
            "kendall_tau_b": total_tau,
            "kendall_pair_counts": total_pair_counts,
            "stratified_randomization": stratified_randomization(
                total_dimensions, gains, blocks
            ),
        },
        "by_dim_Q": grouped_summary(rows),
        "distinct_dim_Q_values": sorted(set(dimensions)),
        "any_escape_count": sum(gain > 0 for gain in gains),
        "certified_gain_sum": sum(gains),
    }


def build():
    features, protocol, outcomes, rows = load_and_join()
    complete_rows = [row for row in rows if row["complete_stage_a"]]
    all_complete = len(complete_rows) == len(rows)
    confirmatory = analysis_for_rows(rows, confirmatory=True) if all_complete else {
        "role": "confirmatory",
        "row_count": len(rows),
        "kendall_tau_b": None,
        "stratified_randomization": None,
        "reason_null": "AT_LEAST_ONE_SCHEDULED_STAGE_A_ROW_WAS_CENSORED",
    }
    exploratory = (
        None
        if all_complete or not complete_rows
        else analysis_for_rows(complete_rows, confirmatory=False)
    )
    result = {
        "schema": "elkies-k3.r17-small-field-class-quotient-analysis.v1",
        "status": (
            "COMPLETE_CONFIRMATORY_CLASS_QUOTIENT_PREDICTION_TEST"
            if all_complete
            else "CENSORED_CONFIRMATORY_TEST_NULL"
        ),
        "question": protocol["predeclared_analysis"]["predictor"] + " predicts future certified MW escape?",
        "primary": confirmatory,
        "exploratory_complete_case": exploratory,
        "censoring": {
            "scheduled_rows": len(rows),
            "complete_stage_a_rows": len(complete_rows),
            "censored_rows": len(rows) - len(complete_rows),
            "status_counts": dict(sorted(Counter(row["detector_status"] for row in rows).items())),
            "censored_rows_counted_as_zero_gain": False,
        },
        "joined_rows": rows,
        "inputs": {
            relative(path): digest(path) for path in (FEATURES, PROTOCOL, OUTCOMES)
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                ".venv/bin/python elkies-k3/scripts/"
                "analyze_r17_small_field_class_quotient_experiment.py"
            ),
        },
        "claim_boundary": protocol["claim_boundary"],
    }
    result["analysis_sha256"] = canonical_hash(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored small-field analysis differs from exact replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        f"R17SMALLFIELDANALYSIS|status={document['status']}"
        f"|rows={document['censoring']['scheduled_rows']}|output={relative(output)}"
    )


if __name__ == "__main__":
    main()
