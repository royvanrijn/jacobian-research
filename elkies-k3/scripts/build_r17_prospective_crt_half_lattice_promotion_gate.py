#!/usr/bin/env python3
"""Fail closed on rank-32 promotion from the binary half-lattice endpoint."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-protocol-v3.json"
)
ABLATION = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/half_lattice_search_ablation_summary_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-promotion-gate-v1.json"
)
POLICY_SOURCE = ROOT / "elliptic-curves/cas/half_lattice_chart_policy.py"
ANALYZER_SOURCE = (
    ROOT
    / "elkies-k3/scripts/"
    "analyze_r17_prospective_crt_half_lattice_censoring_gated.py"
)

sys.path.insert(0, str(POLICY_SOURCE.parent))
from half_lattice_chart_policy import (  # noqa: E402
    bind_ordering,
    policy_document,
    validate_ordering,
)

EXPECTED_PROTOCOL_SHA256 = (
    "a402b1a286dd72ad579c753315a55309a92d03886f60ac8fee84e434119da626"
)
EXPECTED_ABLATION_SHA256 = (
    "fbdfa24b14bc86ee33a576f5e3c3e894dd91dd5e0d1fbfb47bf208e167a7282a"
)
EXPECTED_PROTOCOL_DEFINITION_SHA256 = (
    "9584174de7625031e5f95ce73d0117a9caf8341d91063061ea672f2e4e36e521"
)
GENERIC_RANK = 17
TARGET_RANK = 32
REQUIRED_JUMP = TARGET_RANK - GENERIC_RANK


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def selected_score(case: dict[str, Any]) -> int:
    arms = {arm["id"]: arm for arm in case["arms"]}
    return int(arms["generic-deepest43"]["exact_quotient_rank_over_Q"])


def calibration_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "label": case["label"],
            "displayed_quotient_dimension": int(case["target_quotient_dimension"]),
            "stage_a_exact_certified_quotient_gain": selected_score(case),
        }
        for case in cases
    ]


def strict_order_reversals(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    reversals = []
    for left_index, left in enumerate(rows):
        for right in rows[left_index + 1 :]:
            low, high = sorted(
                (left, right), key=lambda row: row["displayed_quotient_dimension"]
            )
            if (
                low["displayed_quotient_dimension"]
                < high["displayed_quotient_dimension"]
                and low["stage_a_exact_certified_quotient_gain"]
                > high["stage_a_exact_certified_quotient_gain"]
            ):
                reversals.append(
                    {
                        "lower_jump_label": low["label"],
                        "lower_jump": low["displayed_quotient_dimension"],
                        "lower_jump_score": low[
                            "stage_a_exact_certified_quotient_gain"
                        ],
                        "higher_jump_label": high["label"],
                        "higher_jump": high["displayed_quotient_dimension"],
                        "higher_jump_score": high[
                            "stage_a_exact_certified_quotient_gain"
                        ],
                    }
                )
    return reversals


def executable_policy_audit() -> dict[str, Any]:
    fixture = {
        "basis_records": [{"x": "1", "y": "2"}, {"x": "3", "y": "4"}],
        "height_gram_rows": [["2", "1/2"], ["1/2", "3"]],
        "generic_coordinate_rows": [[1, 0]],
        "quotient_coordinate_rows": [[0, 1]],
        "chart_universe_id": "promotion-gate-self-audit",
        "ordered_chart_ids": ["chart-3", "chart-1", "chart-2"],
        "heuristics": [
            "legacy_half_lattice_depth",
            "old_deep_43",
            "quotient_hamming_weight",
        ],
    }
    certificate = bind_ordering(**fixture)
    validation_args = {
        key: value for key, value in fixture.items() if key != "heuristics"
    }
    validate_ordering(certificate, **validation_args)
    rejected_changes = []
    for label, field, changed_value in (
        (
            "lattice_enlargement",
            "basis_records",
            fixture["basis_records"] + [{"x": "5", "y": "6"}],
        ),
        ("basis_change", "basis_records", list(reversed(fixture["basis_records"]))),
        ("height_gram_change", "height_gram_rows", [["2", "0"], ["0", "3"]]),
        ("quotient_basis_change", "quotient_coordinate_rows", [[1, 1]]),
        ("chart_order_change", "ordered_chart_ids", ["chart-1", "chart-2", "chart-3"]),
    ):
        changed = {**validation_args, field: changed_value}
        try:
            validate_ordering(certificate, **changed)
        except ValueError:
            rejected_changes.append(label)
        else:
            raise AssertionError(f"chart policy accepted stale state: {label}")
    return {
        "unchanged_exact_state_accepted": True,
        "rejected_changes": rejected_changes,
        "all_declared_test_changes_rejected": len(rejected_changes) == 5,
        "certificate": certificate,
    }


def build() -> dict[str, Any]:
    pinned = {
        PROTOCOL: EXPECTED_PROTOCOL_SHA256,
        ABLATION: EXPECTED_ABLATION_SHA256,
    }
    for path, expected in pinned.items():
        if digest(path) != expected:
            raise ArithmeticError(f"pinned input changed: {relative(path)}")

    protocol = json.loads(PROTOCOL.read_text())
    ablation = json.loads(ABLATION.read_text())
    if protocol.get("protocol_definition_sha256") != EXPECTED_PROTOCOL_DEFINITION_SHA256:
        raise ArithmeticError("the assessed half-lattice protocol definition changed")
    if protocol.get("status") != "FROZEN_AFTER_POSITIVE_CONTROLS_BEFORE_NEW_COHORT_OUTCOMES":
        raise ArithmeticError("the assessed half-lattice protocol is not frozen")
    if ablation.get("status") != "PASS_EQUAL_BUDGET_BLIND_ABLATION":
        raise ArithmeticError("the half-lattice calibration artifact is not passing")

    primary = protocol["primary_outcome"]
    if primary.get("event") != "at least one exactly certified Stage-A quotient direction":
        raise ArithmeticError("the assessed endpoint is no longer the frozen binary event")

    development = calibration_rows(ablation["development"]["cases"])
    holdout = calibration_rows(ablation["holdout"]["cases"])
    all_rows = development + holdout
    holdout_strata = sorted(
        {row["displayed_quotient_dimension"] for row in holdout}
    )
    maximum_calibrated_jump = max(
        row["displayed_quotient_dimension"] for row in all_rows
    )
    development_reversals = strict_order_reversals(development)

    current_checks = {
        "assessed_v3_explicitly_forbids_binary_endpoint_for_promotion": False,
        "analyzer_requires_balanced_censoring_for_any_prospective_contrast": True,
        "score_was_validated_on_an_independent_panel_with_multiple_jump_strata": (
            len(holdout_strata) >= 2
        ),
        "score_has_predeclared_magnitude_tracking_acceptance_rule": False,
        "score_has_predeclared_upper_tail_acceptance_rule": False,
        "prospective_rank_or_exact_quotient_tail_outcomes_are_available": False,
    }
    if all(current_checks.values()):
        raise AssertionError("v3 unexpectedly passed the fail-closed promotion audit")
    promotion_authorized = all(current_checks.values())

    body = {
        "schema": "elkies-k3.r17-prospective-crt-half-lattice-promotion-gate.v1",
        "status": "FAIL_CLOSED_NO_RANK32_CANDIDATE_PROMOTION",
        "assessed_protocol": {
            "protocol_id": protocol["protocol_id"],
            "protocol_definition_sha256": protocol["protocol_definition_sha256"],
            "candidate_list_sha256": protocol["candidate_list_sha256"],
            "primary_event": primary["event"],
            "primary_estimand": primary["primary_estimand"],
            "permitted_inference": (
                "a complete-case comparison of bounded Stage-A detector-visible "
                "escape yield only when every censor-status proportion is exactly "
                "balanced between the compared frozen cohorts"
            ),
        },
        "amendment_boundary": {
            "kind": "POST_FREEZE_RESTRICTIVE_INTERPRETATION_GUARD",
            "does_not_change_candidates_searches_or_v3_estimand": True,
            "applies_regardless_of_any_partial_or_final_v3_outcome": True,
            "cannot_convert_v3_into_a_confirmatory_magnitude_or_tail_study": True,
            "may_only_narrow_the_inference_from_v3": True,
        },
        "rank_32_target": {
            "generic_rank": GENERIC_RANK,
            "target_rank": TARGET_RANK,
            "required_quotient_jump": REQUIRED_JUMP,
            "direct_certificate_rule": (
                "Fifteen exactly certified independent directions beyond specialized "
                "MW17 directly prove rank at least 32; that is a lower-bound certificate, "
                "not score-based candidate promotion."
            ),
        },
        "detector_score": {
            "definition_for_any_successor": (
                "exact finite-reduction-certified Stage-A quotient gain beyond "
                "specialized MW17 per scheduled fibre"
            ),
            "raw_point_count_forbidden": True,
            "binary_any_escape_forbidden": True,
            "conditional_stage_b_gain_cannot_be_an_unconditional_cohort_score": True,
        },
        "censoring_policy": {
            "scheduled_rows_may_not_be_counted_as_completed_non_events": True,
            "primary_effect_denominator": "rows with complete Stage-A outcomes",
            "balance_rule": (
                "each distinct censored Stage-A status must have exactly the same "
                "scheduled-row proportion in both compared arms"
            ),
            "failure_action": (
                "set risk difference, risk ratio, odds ratio, and Fisher exact p to "
                "null and authorize no prospective cohort conclusion"
            ),
            "observed_events_on_censored_rows": (
                "retain as exact descriptive lower-bound events only; do not place "
                "the row in an inferential event-rate denominator"
            ),
        },
        "chart_order_policy": policy_document(),
        "executable_chart_order_policy_audit": executable_policy_audit(),
        "legacy_field_interpretation": {
            "depth_fields": (
                "retained only as fixed-presentation CVP chart-priority scores; "
                "they are not arithmetic depths or filtrations"
            ),
            "old_deep_43": (
                "a calibrated search-order list for its recorded generic basis, "
                "not a basis-invariant arithmetic subset"
            ),
            "quotient_hamming_weight": (
                "an enumeration order in one recorded quotient basis, not a "
                "Selmer or Mordell--Weil filtration"
            ),
            "historical_exact_point_certificates_remain_valid": True,
            "historical_misses_gain_no_new_meaning": True,
        },
        "available_calibration": {
            "development_rows_used_to_choose_the_detector": development,
            "sealed_holdout_rows": holdout,
            "sealed_holdout_distinct_displayed_jump_strata": holdout_strata,
            "sealed_holdout_score_range_at_jump_12": [
                min(row["stage_a_exact_certified_quotient_gain"] for row in holdout),
                max(row["stage_a_exact_certified_quotient_gain"] for row in holdout),
            ],
            "development_strict_order_reversals": development_reversals,
            "maximum_calibrated_displayed_jump": maximum_calibrated_jump,
            "gap_from_rank32_jump": REQUIRED_JUMP - maximum_calibrated_jump,
            "interpretation": (
                "The sealed controls establish bounded-search sensitivity and class-set "
                "enrichment, but they all have displayed jump 12 and therefore cannot "
                "identify score-versus-jump magnitude or upper-tail calibration."
            ),
        },
        "promotion_requirements": {
            "independent_magnitude_validation": (
                "Freeze the score and an outcome-independent validation panel spanning "
                "multiple exact displayed-quotient or certified rank-jump strata; "
                "predeclare and pass a directional association rule with uncertainty."
            ),
            "tail_validation": (
                "Before opening outcomes, define the upper-tail rank/jump threshold and "
                "minimum effect, then show that the score enriches that tail rather than "
                "only the event jump>=1."
            ),
            "candidate_rule": (
                "Freeze a score threshold justified by those magnitude and tail checks. "
                "Neither one detected direction nor a significant binary cohort contrast "
                "may promote a fibre."
            ),
            "existing_arithmetic_gate": (
                "Any expensive rank-32 follow-up still requires the separately declared "
                "completed residual 2-Selmer quotient gate on the same minimal curve."
            ),
            "state_bound_chart_ordering": (
                "After every independent-point addition, finite-index enlargement, "
                "generator or basis change, height-form change, or quotient-basis "
                "change, discard the old order. Recompute chart identities, shortest "
                "representatives, scores, and the intended order, then bind them to "
                "the new lattice-state fingerprint before searching. Efficiency or "
                "enrichment calibration does not transfer to the changed state."
            ),
            "balanced_censoring": (
                "Before reporting any prospective cohort contrast, pass the exact "
                "censor-status balance rule and use only completed Stage-A rows in "
                "the event-rate denominator."
            ),
        },
        "current_gate_checks": current_checks,
        "rank_32_candidate_promotion_authorized": promotion_authorized,
        "claim_boundary": [
            "Protocol v3 remains an immutable bounded-search record, but its scheduled-denominator contrast is not authorized for inference.",
            "No v3 binary success, effect size, or p-value can authorize a rank-32 candidate.",
            "No prospective cohort effect is estimable unless censoring is exactly balanced by failure status between its arms.",
            "The available controls do not validate magnitude tracking or extreme-tail incidence.",
            "Legacy depth, old-deep-43, and quotient-Hamming fields order birational search charts only; they encode no covering, Selmer, or rank structure.",
            "A chart miss cannot prove rational-point absence, local insolubility, Selmer structure, Mordell--Weil saturation, or a rank upper bound.",
            (
                "A score g exactly certifies rank at least 17+g on that fibre, but g<15 "
                "cannot promote it on the hope of enough unseen directions; a bounded "
                "miss is not a rank upper bound."
            ),
        ],
    }
    return {
        **body,
        "gate_definition_sha256": canonical_hash(body),
        "inputs": {
            **{relative(path): digest(path) for path in pinned},
            relative(POLICY_SOURCE): digest(POLICY_SOURCE),
            relative(ANALYZER_SOURCE): digest(ANALYZER_SOURCE),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                "python3 elkies-k3/scripts/"
                "build_r17_prospective_crt_half_lattice_promotion_gate.py"
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
            raise ArithmeticError("stored half-lattice promotion gate differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17CRTHALFPROMOTION"
        f"|status={document['status']}"
        f"|authorized={document['rank_32_candidate_promotion_authorized']}"
        f"|hash={document['gate_definition_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
