#!/usr/bin/env python3
"""Replay the frozen v2 ladder tests without requiring public containment.

The frozen score is the exact blind rank gain in E(Q)/MW17.  It was never
defined as the rank of the blind points inside the later-opened displayed
subgroup.  The original verifier imposed that stronger condition and found
two informative failures.  This analyzer preserves the predeclared Kendall
and tail formulas while treating integral public containment descriptively.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import comb, factorial, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v2.json"
BLIND = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v2.json"
VERIFICATION = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_verification_v2.json"
FROZEN_ANALYZER = ROOT / "elliptic-curves/cas/analyze_r17_refresh_jump_ladder.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_analysis_v2.json"

FRAME_BY_J_CLASS = {
    "norm12-orbit-074d9": "published-R17",
    "norm12-orbit-07ca9": "published-R17",
    "norm12-orbit-08234": "published-R17",
    "norm12-orbit-0e80b": "published-R17",
    "norm12-orbit-08f72": "alternate-Q80",
    "norm12-orbit-11952": "alternate-Q80",
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def exact_kendall_distribution(frozen, scores, truths):
    """Return the exact tied-margin distribution used by the frozen test."""

    truth_values = sorted(set(truths))
    truth_counts = tuple(truths.count(value) for value in truth_values)
    score_counts = sorted(Counter(scores).items())
    states = {((0,) * len(truth_values), 0): 1}
    for _score, block_size in score_counts:
        updated = defaultdict(int)
        for (used, statistic), ways in states.items():
            remaining = tuple(
                truth_counts[index] - used[index]
                for index in range(len(truth_values))
            )
            for allocation in frozen.allocations(block_size, remaining):
                contribution = 0
                for current_index, amount in enumerate(allocation):
                    lower = sum(used[:current_index])
                    higher = sum(used[current_index + 1 :])
                    contribution += amount * (lower - higher)
                within_ways = factorial(block_size)
                for amount in allocation:
                    within_ways //= factorial(amount)
                new_used = tuple(
                    used[index] + allocation[index]
                    for index in range(len(truth_values))
                )
                updated[(new_used, statistic + contribution)] += ways * within_ways
        states = updated
    distribution = Counter()
    for (used, statistic), ways in states.items():
        if used != truth_counts:
            raise ArithmeticError("a stratified Kendall block ended early")
        distribution[statistic] += ways
    expected_total = factorial(len(truths))
    for count in truth_counts:
        expected_total //= factorial(count)
    if sum(distribution.values()) != expected_total:
        raise ArithmeticError("a stratified Kendall block has the wrong mass")
    return distribution


def tail_summary(frozen, scores, truths, truth_threshold, score_threshold=10):
    tails = [truth >= truth_threshold for truth in truths]
    positives = [score >= score_threshold for score in scores]
    a = sum(tail and positive for tail, positive in zip(tails, positives))
    c = sum(tail and not positive for tail, positive in zip(tails, positives))
    b = sum(not tail and positive for tail, positive in zip(tails, positives))
    d = sum(not tail and not positive for tail, positive in zip(tails, positives))
    table = {
        "true_tail_detector_positive": a,
        "true_tail_detector_negative": c,
        "non_tail_detector_positive": b,
        "non_tail_detector_negative": d,
    }
    if a + c == 0 or b + d == 0:
        return {
            "true_tail_threshold": truth_threshold,
            "detector_positive_threshold": score_threshold,
            "table": table,
            "estimable": False,
            "reason_null": "STRATUM_HAS_NO_TAIL_OR_NO_NONTAIL_CASES",
            "detector_positive_rate_in_true_tail": None,
            "detector_positive_rate_outside_true_tail": None,
            "risk_difference": None,
            "odds_ratio": None,
            "fisher_exact_one_sided_p": None,
            "exact_p_numerator": None,
            "exact_p_denominator": None,
        }
    return {
        "true_tail_threshold": truth_threshold,
        "detector_positive_threshold": score_threshold,
        "estimable": True,
        **frozen.fisher_greater(a, b, c, d),
    }


def stratum_summary(frozen, key, rows):
    scores = [
        int(row["exact_quotient_rank_recovered_before_public_complement"])
        for row in rows
    ]
    truths = [int(row["true_displayed_jump_opened_after_blind_freeze"]) for row in rows]
    return {
        "stratum": key,
        "case_count": len(rows),
        "curve_ids": [int(row["curve_id"]) for row in rows],
        "S_q_pairs_by_curve": [
            {
                "curve_id": int(row["curve_id"]),
                "S_exact_blind_recovered_quotient_rank": score,
                "q_exact_displayed_subgroup_quotient_rank": truth,
            }
            for row, score, truth in zip(rows, scores, truths)
        ],
        "ordinal_association": frozen.exact_kendall_permutation(scores, truths),
        "high_S_enrichment": {
            "q_at_least_10": tail_summary(frozen, scores, truths, 10),
            "q_at_least_11": tail_summary(frozen, scores, truths, 11),
        },
    }


def stratified_ordinal(frozen, groups):
    """Exact randomization of q within each fixed fibration or j-class."""

    distribution = Counter({0: 1})
    observed = 0
    score_comparable_pairs = 0
    truth_comparable_pairs = 0
    for rows in groups:
        scores = [
            int(row["exact_quotient_rank_recovered_before_public_complement"])
            for row in rows
        ]
        truths = [
            int(row["true_displayed_jump_opened_after_blind_freeze"])
            for row in rows
        ]
        observed += frozen.kendall_s(scores, truths)
        n0 = comb(len(rows), 2)
        score_comparable_pairs += n0 - sum(
            comb(count, 2) for count in Counter(scores).values()
        )
        truth_comparable_pairs += n0 - sum(
            comb(count, 2) for count in Counter(truths).values()
        )
        block = exact_kendall_distribution(frozen, scores, truths)
        updated = Counter()
        for left, left_ways in distribution.items():
            for right, right_ways in block.items():
                updated[left + right] += left_ways * right_ways
        distribution = updated
    total = sum(distribution.values())
    extreme = sum(
        ways for statistic, ways in distribution.items() if statistic >= observed
    )
    denominator = sqrt(score_comparable_pairs * truth_comparable_pairs)
    return {
        "statistic": "sum of within-stratum concordance-minus-discordance",
        "cross_stratum_pairs_excluded": True,
        "observed_statistic": observed,
        "block_restricted_kendall_tau_b": observed / denominator if denominator else None,
        "exact_one_sided_permutation_p": float(Fraction(extreme, total)),
        "exact_p_numerator": extreme,
        "exact_p_denominator": total,
        "null_assignment_count": total,
    }


def stratified_tail(groups, truth_threshold, score_threshold=10):
    """Exact tail randomization conditional on each stratum's margins."""

    distribution = Counter({0: 1})
    observed = 0
    strata = []
    for key, rows in groups:
        scores = [
            int(row["exact_quotient_rank_recovered_before_public_complement"])
            for row in rows
        ]
        truths = [
            int(row["true_displayed_jump_opened_after_blind_freeze"])
            for row in rows
        ]
        tail_count = sum(truth >= truth_threshold for truth in truths)
        positive_count = sum(score >= score_threshold for score in scores)
        observed_hits = sum(
            truth >= truth_threshold and score >= score_threshold
            for score, truth in zip(scores, truths)
        )
        observed += observed_hits
        block = Counter()
        lower = max(0, positive_count - (len(rows) - tail_count))
        upper = min(tail_count, positive_count)
        for hits in range(lower, upper + 1):
            block[hits] = comb(tail_count, hits) * comb(
                len(rows) - tail_count, positive_count - hits
            )
        if sum(block.values()) != comb(len(rows), positive_count):
            raise ArithmeticError("a stratified tail block has the wrong mass")
        updated = Counter()
        for left, left_ways in distribution.items():
            for right, right_ways in block.items():
                updated[left + right] += left_ways * right_ways
        distribution = updated
        strata.append(
            {
                "stratum": key,
                "case_count": len(rows),
                "true_tail_count": tail_count,
                "detector_positive_count": positive_count,
                "observed_tail_detector_positive": observed_hits,
            }
        )
    total = sum(distribution.values())
    extreme = sum(ways for statistic, ways in distribution.items() if statistic >= observed)
    return {
        "true_tail_threshold": truth_threshold,
        "detector_positive_threshold": score_threshold,
        "statistic": "total high-S cases in the true tail after conditioning on stratum margins",
        "observed_statistic": observed,
        "strata": strata,
        "exact_one_sided_permutation_p": float(Fraction(extreme, total)),
        "exact_p_numerator": extreme,
        "exact_p_denominator": total,
        "null_assignment_count": total,
    }


def build():
    protocol = json.loads(PROTOCOL.read_text())
    blind = json.loads(BLIND.read_text())
    verified = json.loads(VERIFICATION.read_text())
    if protocol.get("freeze_boundary") != "AFTER_CURVE478_BEFORE_OTHER_RECOVERY_OUTCOMES":
        raise ArithmeticError("the disclosed v2 freeze boundary changed")
    if blind.get("status") != "PASS_COMPLETE_BLIND_RECOVERY_BEFORE_PUBLIC_COMPLEMENT":
        raise ArithmeticError("the blind recovery artifact is incomplete")
    if verified["phase_boundary"]["blind_artifact_sha256_before_public_complement_import"] != digest(BLIND):
        raise ArithmeticError("the public complement was opened against another blind artifact")
    if not verified["phase_boundary"][
        "public_complement_opened_only_after_blind_status_and_hash_were_fixed"
    ]:
        raise ArithmeticError("the public-complement phase boundary failed")

    frozen = SourceFileLoader("r17_jump_ladder_frozen_analysis", str(FROZEN_ANALYZER)).load_module()
    blind_by_id = {int(row["curve_id"]): row for row in blind["results"]}
    rows = sorted(verified["results"], key=lambda row: int(row["curve_id"]))
    if [int(row["curve_id"]) for row in rows] != protocol["eligible_curve_ids"]:
        raise ArithmeticError("the analyzed cases differ from the frozen panel")
    observed_j_classes = {row["representative_class"] for row in rows}
    if observed_j_classes != set(FRAME_BY_J_CLASS):
        raise ArithmeticError("the six-class fibration split changed")
    for row in rows:
        row["fibration_class"] = FRAME_BY_J_CLASS[row["representative_class"]]
    scores = [
        int(blind_by_id[int(row["curve_id"])][
            "exact_quotient_rank_recovered_before_public_complement"
        ])
        for row in rows
    ]
    truths = [int(row["true_displayed_jump_opened_after_blind_freeze"]) for row in rows]
    ordinal = frozen.exact_kendall_permutation(scores, truths)
    ordinal_pass = (
        ordinal["kendall_tau_b"] is not None
        and ordinal["kendall_tau_b"] >= 0.35
        and ordinal["exact_one_sided_permutation_p"] <= 0.05
    )

    tails = [truth >= 10 for truth in truths]
    positives = [score >= 10 for score in scores]
    a = sum(tail and positive for tail, positive in zip(tails, positives))
    c = sum(tail and not positive for tail, positive in zip(tails, positives))
    b = sum(not tail and positive for tail, positive in zip(tails, positives))
    d = sum(not tail and not positive for tail, positive in zip(tails, positives))
    tail = frozen.fisher_greater(a, b, c, d)
    tail_pass = (
        tail["risk_difference"] >= 0.25
        and tail["fisher_exact_one_sided_p"] <= 0.05
    )
    joint_pass = ordinal_pass and tail_pass

    by_fibration_groups = defaultdict(list)
    by_j_class_groups = defaultdict(list)
    for row in rows:
        by_fibration_groups[row["fibration_class"]].append(row)
        by_j_class_groups[row["representative_class"]].append(row)
    by_fibration = [
        stratum_summary(frozen, key, stratum_rows)
        for key, stratum_rows in sorted(by_fibration_groups.items())
    ]
    by_j_class = [
        stratum_summary(frozen, key, stratum_rows)
        for key, stratum_rows in sorted(by_j_class_groups.items())
    ]
    fibration_group_items = sorted(by_fibration_groups.items())
    j_class_group_items = sorted(by_j_class_groups.items())

    by_jump = []
    for jump in sorted(set(truths)):
        values = [score for score, truth in zip(scores, truths) if truth == jump]
        by_jump.append(
            {
                "true_displayed_jump": jump,
                "case_count": len(values),
                "exact_blind_recovered_ranks": sorted(values),
                "mean_exact_blind_recovered_rank": float(
                    Fraction(sum(values), len(values))
                ),
            }
        )

    return {
        "schema": "elliptic-curves.r17-refresh-jump-ladder-analysis.v2",
        "status": (
            "PASS_USABLE_EXTREME_JUMP_DETECTOR"
            if joint_pass
            else "FAIL_STOP_SERIOUS_RANK32_HALF_LATTICE_BUDGET"
        ),
        "response": [
            {
                "curve_id": int(row["curve_id"]),
                "fibration_class": row["fibration_class"],
                "j_map_class": row["representative_class"],
                "exact_quotient_rank_recovered_before_public_complement": score,
                "exact_displayed_subgroup_quotient_rank_opened_after_blind_freeze": truth,
            }
            for row, score, truth in zip(rows, scores, truths)
        ],
        "confirmatory": {
            "ordinal_association": {
                **ordinal,
                "predeclared_acceptance": "tau_b >= 0.35 and exact one-sided p <= 0.05",
                "passed": ordinal_pass,
            },
            "upper_tail_enrichment": {
                "true_tail_threshold": 10,
                "detector_positive_threshold": 10,
                **tail,
                "predeclared_acceptance": "risk difference >= 0.25 and one-sided Fisher exact p <= 0.05",
                "passed": tail_pass,
            },
            "joint_decision": {
                "both_endpoints_required": True,
                "passed": joint_pass,
                "action": (
                    "retain half-lattice recovery as a usable extreme-jump detector, subject to the separate residual 2-Selmer gate"
                    if joint_pass
                    else "stop serious rank-32 point-search spending on the half-lattice extreme-jump hypothesis"
                ),
            },
        },
        "descriptive": {
            "by_displayed_jump": by_jump,
            "post_freeze_fibration_and_j_class_sensitivity": {
                "role": "descriptive_post_freeze_not_confirmatory",
                "score_definition": (
                    "S = exact blind rank of the discovered subgroup modulo "
                    "specialized MW17"
                ),
                "truth_definition": (
                    "q = exact free rank of the displayed public subgroup modulo "
                    "specialized MW17"
                ),
                "pooled_q_at_least_11": tail_summary(frozen, scores, truths, 11),
                "by_fibration": by_fibration,
                "by_j_map_class": by_j_class,
                "conditional_on_fibration": {
                    "ordinal_association": stratified_ordinal(
                        frozen, list(by_fibration_groups.values())
                    ),
                    "high_S_enrichment": {
                        "q_at_least_10": stratified_tail(
                            fibration_group_items, 10
                        ),
                        "q_at_least_11": stratified_tail(
                            fibration_group_items, 11
                        ),
                    },
                },
                "conditional_on_j_map_class": {
                    "ordinal_association": stratified_ordinal(
                        frozen, list(by_j_class_groups.values())
                    ),
                    "high_S_enrichment": {
                        "q_at_least_10": stratified_tail(j_class_group_items, 10),
                        "q_at_least_11": stratified_tail(j_class_group_items, 11),
                    },
                },
                "operational_read": (
                    "High S may promote scheduling inside the calibrated norm-twelve "
                    "R17 setting, but low S may not veto a candidate. Transfer to "
                    "alternate-Q80 or an untested j-class is unvalidated because the "
                    "alternate-Q80 rows contain no q>=10 case."
                ),
            },
            "total_exact_blind_recovered_rank": sum(scores),
            "total_displayed_jump": sum(truths),
            "overall_recovery_fraction": f"{sum(scores)}/{sum(truths)}",
            "full_344_chart_cases": sum(
                int(blind_by_id[int(row["curve_id"])]["attempted_chart_count"])
                == 344
                for row in rows
            ),
            "structural_zero_43_chart_cases": [
                int(row["curve_id"])
                for row in rows
                if int(blind_by_id[int(row["curve_id"])]["attempted_chart_count"])
                == 43
            ],
            "total_timeouts": sum(
                int(blind_by_id[int(row["curve_id"])].get("timeout_chart_count") or 0)
                for row in rows
            ),
            "total_pari_failures": sum(
                int(blind_by_id[int(row["curve_id"])].get("pari_failure_chart_count") or 0)
                for row in rows
            ),
            "all_blind_basis_points_integral_in_opened_displayed_subgroup_count": sum(
                bool(row["all_final_blind_basis_points_in_opened_public_subgroup"])
                for row in rows
            ),
            "integral_public_containment_failure_curve_ids": [
                int(row["curve_id"])
                for row in rows
                if not row["all_final_blind_basis_points_in_opened_public_subgroup"]
            ],
        },
        "analysis_boundary_correction": {
            "frozen_score_definition_changed": False,
            "frozen_statistics_or_thresholds_changed": False,
            "removed_extra_unfrozen_requirement": (
                "the v1 verifier required every blind point to be an integral "
                "combination of the displayed public basis, although the protocol "
                "defined the score only as exact rank gain in E(Q)/MW17"
            ),
            "effect": (
                "curve 478 certifies rank at least 23, two above its displayed rank-21 "
                "subgroup; curve 539 has equal blind and displayed ranks but a "
                "nonintegral-containment issue. Both exact blind rank scores remain valid."
            ),
        },
        "inputs": {
            relative(PROTOCOL): digest(PROTOCOL),
            relative(BLIND): digest(BLIND),
            relative(VERIFICATION): digest(VERIFICATION),
            relative(FROZEN_ANALYZER): digest(FROZEN_ANALYZER),
            relative(Path(__file__).resolve()): digest(Path(__file__).resolve()),
        },
        "claim_boundary": [
            "The exact permutation test conditions on this fixed atlas-refresh panel and tied margins.",
            "The displayed jumps are certified subgroup quotients, not full Mordell-Weil ranks.",
            "Passing does not replace the separately required residual 2-Selmer gate for rank-32 promotion.",
            "Protocol v2 is transparently post-curve-478 because v1 used a false cross-class deepest-count assertion.",
            "The fibration/j-class splits and the q>=11 endpoint were computed after unblinding and are sensitivity analyses, not frozen confirmatory endpoints.",
            "Only the published-R17 frame, and only j-class 08234 within a mixed tail/body class, directly tests high-S enrichment; alternate-Q80 has no q>=10 row in this panel.",
            "A low detector score cannot reject a candidate because every search is bounded; curve 544 is the explicit q=11, S=0 false negative.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != encoded:
            raise ArithmeticError("stored v2 jump-ladder analysis changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(
        f"R17JUMPLADDERANALYSISV2|status={payload['status']}|"
        f"tau={payload['confirmatory']['ordinal_association']['kendall_tau_b']}|"
        f"ordinal_p={payload['confirmatory']['ordinal_association']['exact_one_sided_permutation_p']}|"
        f"tail_p={payload['confirmatory']['upper_tail_enrichment']['fisher_exact_one_sided_p']}"
    )


if __name__ == "__main__":
    main()
