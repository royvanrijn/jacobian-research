#!/usr/bin/env python3
"""Run the predeclared ordinal and upper-tail jump-ladder tests."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb, factorial, sqrt
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v1.json"
VERIFICATION = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_verification_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_analysis_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_hash(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def sign(value: int) -> int:
    return (value > 0) - (value < 0)


def kendall_s(scores, truths) -> int:
    total = 0
    for left in range(len(scores)):
        for right in range(left + 1, len(scores)):
            total += sign(scores[right] - scores[left]) * sign(
                truths[right] - truths[left]
            )
    return total


def allocations(total: int, upper_bounds):
    for values in product(*(range(bound + 1) for bound in upper_bounds)):
        if sum(values) == total:
            yield values


def exact_kendall_permutation(scores, truths):
    """Exact one-sided null distribution, retaining ties in both margins.

    Score-tied cases form blocks.  The dynamic program assigns the fixed
    multiset of truth labels to each block, multiplying by the within-block
    multinomial count.  Contributions from pairs inside a score-tied block
    are zero by construction.
    """

    truth_values = sorted(set(truths))
    truth_counts = tuple(truths.count(value) for value in truth_values)
    score_counts = sorted(Counter(scores).items())
    # states map (used truth counts, concordance-minus-discordance) to the
    # number of distinct assignments to the named fixed cases.
    states = {((0,) * len(truth_values), 0): 1}
    for unused_score, block_size in score_counts:
        updated = defaultdict(int)
        for (used, statistic), ways in states.items():
            remaining = tuple(
                truth_counts[index] - used[index]
                for index in range(len(truth_values))
            )
            for allocation in allocations(block_size, remaining):
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

    observed = kendall_s(scores, truths)
    distribution = Counter()
    for (used, statistic), ways in states.items():
        if used != truth_counts:
            raise ArithmeticError("the exact Kendall dynamic program ended early")
        distribution[statistic] += ways
    expected_total = factorial(len(truths))
    for count in truth_counts:
        expected_total //= factorial(count)
    observed_total = sum(distribution.values())
    if observed_total != expected_total:
        raise ArithmeticError(
            f"exact Kendall assignment count {observed_total} != {expected_total}"
        )
    extreme = sum(ways for statistic, ways in distribution.items() if statistic >= observed)
    n0 = comb(len(scores), 2)
    score_ties = sum(comb(count, 2) for count in Counter(scores).values())
    truth_ties = sum(comb(count, 2) for count in Counter(truths).values())
    denominator = sqrt((n0 - score_ties) * (n0 - truth_ties))
    tau_b = observed / denominator if denominator else None
    return {
        "concordance_minus_discordance": observed,
        "kendall_tau_b": tau_b,
        "exact_one_sided_permutation_p": float(Fraction(extreme, observed_total)),
        "exact_p_numerator": extreme,
        "exact_p_denominator": observed_total,
        "null_assignment_count": observed_total,
        "score_tied_pair_count": score_ties,
        "truth_tied_pair_count": truth_ties,
    }


def fisher_greater(a: int, b: int, c: int, d: int):
    # Rows: true tail/body. Columns: detector positive/negative.
    tail_total = a + c
    body_total = b + d
    positives = a + b
    total = tail_total + body_total
    denominator = comb(total, positives)

    def probability(value):
        return Fraction(
            comb(tail_total, value) * comb(body_total, positives - value),
            denominator,
        )

    upper = min(tail_total, positives)
    p_value = sum((probability(value) for value in range(a, upper + 1)), Fraction(0))
    tail_rate = Fraction(a, tail_total)
    body_rate = Fraction(b, body_total)
    odds_ratio = (
        None if b * c == 0 and a * d == 0 else (
            "Infinity" if b * c == 0 else float(Fraction(a * d, b * c))
        )
    )
    return {
        "table": {
            "true_tail_detector_positive": a,
            "true_tail_detector_negative": c,
            "non_tail_detector_positive": b,
            "non_tail_detector_negative": d,
        },
        "detector_positive_rate_in_true_tail": float(tail_rate),
        "detector_positive_rate_outside_true_tail": float(body_rate),
        "risk_difference": float(tail_rate - body_rate),
        "odds_ratio": odds_ratio,
        "fisher_exact_one_sided_p": float(p_value),
        "exact_p_numerator": p_value.numerator,
        "exact_p_denominator": p_value.denominator,
    }


def build():
    protocol = json.loads(PROTOCOL.read_text())
    verified = json.loads(VERIFICATION.read_text())
    if protocol.get("status") != "FROZEN_BEFORE_BLIND_RECOVERY":
        raise ArithmeticError("the protocol was not frozen")
    if verified.get("status") != "PASS_ALL_BLIND_RECOVERY_RANKS_EXACT_IN_OPENED_DISPLAYED_QUOTIENTS":
        raise ArithmeticError("the opened public-complement verification did not pass")
    rows = sorted(verified["results"], key=lambda row: int(row["curve_id"]))
    if [row["curve_id"] for row in rows] != protocol["eligible_curve_ids"]:
        raise ArithmeticError("the analyzed cases differ from the frozen panel")
    if len(rows) != protocol["confirmatory_analysis"]["case_count"]:
        raise ArithmeticError("the analyzed panel size changed")

    scores = [
        int(row["exact_quotient_rank_recovered_before_public_complement"])
        for row in rows
    ]
    truths = [int(row["true_displayed_jump_opened_after_blind_freeze"]) for row in rows]
    ordinal = exact_kendall_permutation(scores, truths)
    ordinal_pass = (
        ordinal["kendall_tau_b"] is not None
        and ordinal["kendall_tau_b"] >= 0.35
        and ordinal["exact_one_sided_permutation_p"] <= 0.05
    )

    tail = [truth >= 10 for truth in truths]
    positive = [score >= 10 for score in scores]
    a = sum(is_tail and is_positive for is_tail, is_positive in zip(tail, positive))
    c = sum(is_tail and not is_positive for is_tail, is_positive in zip(tail, positive))
    b = sum(not is_tail and is_positive for is_tail, is_positive in zip(tail, positive))
    d = sum(not is_tail and not is_positive for is_tail, is_positive in zip(tail, positive))
    upper_tail = fisher_greater(a, b, c, d)
    tail_pass = (
        upper_tail["risk_difference"] >= 0.25
        and upper_tail["fisher_exact_one_sided_p"] <= 0.05
    )
    joint_pass = ordinal_pass and tail_pass

    by_jump = []
    for jump in sorted(set(truths)):
        stratum_scores = [
            score for score, truth in zip(scores, truths) if truth == jump
        ]
        by_jump.append(
            {
                "displayed_jump": jump,
                "case_count": len(stratum_scores),
                "recovered_ranks": sorted(stratum_scores),
                "mean_recovered_rank": float(
                    Fraction(sum(stratum_scores), len(stratum_scores))
                ),
                "minimum_recovered_rank": min(stratum_scores),
                "maximum_recovered_rank": max(stratum_scores),
            }
        )

    body = {
        "schema": "elliptic-curves.r17-refresh-jump-ladder-analysis.v1",
        "status": (
            "PASS_USABLE_EXTREME_JUMP_DETECTOR"
            if joint_pass
            else "FAIL_STOP_SERIOUS_RANK32_HALF_LATTICE_BUDGET"
        ),
        "response": [
            {
                "curve_id": row["curve_id"],
                "exact_quotient_rank_recovered_before_public_complement": row[
                    "exact_quotient_rank_recovered_before_public_complement"
                ],
            }
            for row in rows
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
                **upper_tail,
                "predeclared_acceptance": "risk difference >= 0.25 and one-sided Fisher exact p <= 0.05",
                "passed": tail_pass,
            },
            "joint_decision": {
                "both_endpoints_required": True,
                "passed": joint_pass,
                "action": (
                    "half-lattice recovery may be retained as a usable extreme-jump detector, subject to the separate residual-Selmer gate"
                    if joint_pass
                    else "stop serious rank-32 point-search spending on the half-lattice extreme-jump hypothesis"
                ),
            },
        },
        "descriptive": {
            "by_displayed_jump": by_jump,
            "total_recovered_rank": sum(scores),
            "total_displayed_jump": sum(truths),
            "overall_recovery_fraction": f"{sum(scores)}/{sum(truths)}",
            "all_cases_used_full_344_chart_budget": all(
                int(row["attempted_chart_count"]) == 344 for row in rows
            ),
            "total_timeouts": sum(int(row["timeout_chart_count"]) for row in rows),
            "total_pari_failures": sum(
                int(row["pari_failure_chart_count"]) for row in rows
            ),
        },
        "claim_boundary": [
            "The p-values condition on this fixed sixteen-fibre panel and its tied margins.",
            "Association does not prove that the displayed subgroups are full Mordell-Weil groups.",
            "A passing detector endpoint would not replace the residual 2-Selmer rank-32 promotion gate.",
        ],
    }
    return {
        **body,
        "analysis_definition_sha256": canonical_hash(body),
        "inputs": {
            relative(PROTOCOL): digest(PROTOCOL),
            relative(VERIFICATION): digest(VERIFICATION),
            relative(Path(__file__).resolve()): digest(Path(__file__).resolve()),
        },
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
            raise ArithmeticError("the stored jump-ladder analysis changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(
        f"R17JUMPLADDERANALYSIS|status={payload['status']}|"
        f"tau={payload['confirmatory']['ordinal_association']['kendall_tau_b']}|"
        f"tail_p={payload['confirmatory']['upper_tail_enrichment']['fisher_exact_one_sided_p']}"
    )


if __name__ == "__main__":
    main()
