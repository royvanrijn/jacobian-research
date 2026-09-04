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
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v2.json"
BLIND = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v2.json"
VERIFICATION = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_verification_v2.json"
FROZEN_ANALYZER = ROOT / "elliptic-curves/cas/analyze_r17_refresh_jump_ladder.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_analysis_v2.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


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
                "exact_quotient_rank_recovered_before_public_complement": score,
            }
            for row, score in zip(rows, scores)
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
