#!/usr/bin/env python3
"""Analyze the frozen R17 CRT half-lattice detector without retuning it."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from math import comb, sqrt
import json
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
PROTOCOL = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-protocol-v3.json"
LEDGER = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-ledger-v3.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-analysis-v3.json"

EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
EXPECTED_PROTOCOL_STATUS = "FROZEN_AFTER_POSITIVE_CONTROLS_BEFORE_NEW_COHORT_OUTCOMES"
Z_975 = 1.959963984540054


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def ratio_record(numerator: int, denominator: int):
    if denominator == 0:
        return None
    value = Fraction(numerator, denominator)
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def wilson_interval(successes: int, total: int):
    if total == 0:
        return None
    estimate = successes / total
    z2 = Z_975**2
    denominator = 1 + z2 / total
    center = (estimate + z2 / (2 * total)) / denominator
    half_width = (
        Z_975
        * sqrt(estimate * (1 - estimate) / total + z2 / (4 * total**2))
        / denominator
    )
    return [max(0.0, center - half_width), min(1.0, center + half_width)]


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    row1 = a + b
    row2 = c + d
    successes = a + c
    total = row1 + row2
    denominator = comb(total, successes)

    def probability(value: int) -> Fraction:
        return Fraction(comb(row1, value) * comb(row2, successes - value), denominator)

    lower = max(0, successes - row2)
    upper = min(successes, row1)
    observed = probability(a)
    return float(
        sum(
            (
                probability(value)
                for value in range(lower, upper + 1)
                if probability(value) <= observed
            ),
            Fraction(0),
        )
    )


def outcome_summary(rows: Iterable[dict[str, Any]]):
    rows = list(rows)
    eligible = [row for row in rows if row["analysis_eligible_complete_stage_a"]]
    events = [row for row in rows if row.get("stage_a", {}).get("certified_quotient_gain", 0) > 0]
    eligible_events = [
        row
        for row in eligible
        if row.get("stage_a", {}).get("certified_quotient_gain", 0) > 0
    ]
    stage_b_rows = [
        row
        for row in rows
        if row.get("stage_b", {}).get("full_specialized_ranking_computed", False)
    ]
    return {
        "scheduled": len(rows),
        "complete_stage_a": len(eligible),
        "censored_stage_a": len(rows) - len(eligible),
        "status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
        "stage_a_certified_escape_rows": len(events),
        "stage_a_certified_directions": sum(
            row.get("stage_a", {}).get("certified_quotient_gain", 0) for row in rows
        ),
        "intention_to_search_detector_yield": ratio_record(len(events), len(rows)),
        "intention_to_search_wilson_95_percent_interval": wilson_interval(
            len(events), len(rows)
        ),
        "complete_case_escape_rows": len(eligible_events),
        "complete_case_escape_rate": ratio_record(len(eligible_events), len(eligible)),
        "complete_case_wilson_95_percent_interval": wilson_interval(
            len(eligible_events), len(eligible)
        ),
        "stage_b_gate_rows": len(stage_b_rows),
        "stage_b_incremental_certified_directions": sum(
            row.get("stage_b", {}).get("incremental_certified_quotient_gain", 0)
            for row in stage_b_rows
        ),
        "largest_certified_rank_lower_bound": max(
            (row.get("largest_certified_rank_lower_bound", 17) for row in rows),
            default=None,
        ),
    }


def comparison_record(label: str, exposed: list[dict], control: list[dict]):
    exposed_summary = outcome_summary(exposed)
    control_summary = outcome_summary(control)
    exposed_events = exposed_summary["stage_a_certified_escape_rows"]
    control_events = control_summary["stage_a_certified_escape_rows"]
    exposed_total = exposed_summary["scheduled"]
    control_total = control_summary["scheduled"]
    exposed_rate = exposed_events / exposed_total
    control_rate = control_events / control_total
    risk_ratio = exposed_rate / control_rate if control_rate else None
    odds_ratio = (
        exposed_events * (control_total - control_events)
        / ((exposed_total - exposed_events) * control_events)
        if exposed_total > exposed_events
        and control_total > control_events
        and control_events
        else None
    )
    return {
        "contrast": label,
        "estimand": "INTENTION_TO_SEARCH_STAGE_A_DETECTOR_YIELD_PER_SCHEDULED_FIBRE",
        "exposed": exposed_summary,
        "control": control_summary,
        "risk_difference": exposed_rate - control_rate,
        "risk_ratio": risk_ratio,
        "odds_ratio": odds_ratio,
        "fisher_exact_two_sided_p": fisher_two_sided(
            exposed_events,
            exposed_total - exposed_events,
            control_events,
            control_total - control_events,
        ),
        "complete_case_sensitivity": {
            "exposed_events": exposed_summary["complete_case_escape_rows"],
            "exposed_total": exposed_summary["complete_stage_a"],
            "control_events": control_summary["complete_case_escape_rows"],
            "control_total": control_summary["complete_stage_a"],
        },
    }


def build():
    manifest = json.loads(MANIFEST.read_text())
    protocol = json.loads(PROTOCOL.read_text())
    ledger = json.loads(LEDGER.read_text())
    if manifest["commitment"]["candidate_list_sha256"] != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the original candidate commitment changed")
    if protocol.get("status") != EXPECTED_PROTOCOL_STATUS:
        raise ArithmeticError("the half-lattice protocol is not frozen")
    if protocol.get("candidate_list_sha256") != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the protocol names another candidate list")
    if ledger.get("status") != "COMPLETE_FROZEN_HALF_LATTICE_DETECTOR_LEDGER":
        raise ArithmeticError("the half-lattice ledger is incomplete")
    if ledger.get("candidate_list_sha256") != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the ledger names another candidate list")
    if ledger.get("protocol_definition_sha256") != protocol["protocol_definition_sha256"]:
        raise ArithmeticError("the outcomes used another detector protocol")

    manifest_by_id = {
        row["sample_id"]: (index, row) for index, row in enumerate(manifest["rows"])
    }
    outcome_by_id = {row["sample_id"]: row for row in ledger["records"]}
    if len(manifest_by_id) != 2_560 or set(manifest_by_id) != set(outcome_by_id):
        raise ArithmeticError("manifest and ledger do not cover the same 2,560 fibres")

    by_cohort = defaultdict(list)
    by_cohort_anchor = defaultdict(list)
    for sample_id, (manifest_index, manifest_row) in manifest_by_id.items():
        outcome = outcome_by_id[sample_id]
        if (
            outcome["manifest_index"] != manifest_index
            or outcome["parameter"] != manifest_row["parameter"]
            or outcome["cohort"] != manifest_row["cohort"]
            or outcome["anchor_curve_id"] != manifest_row["anchor_curve_id"]
        ):
            raise ArithmeticError("a joined outcome changed its frozen row identity")
        by_cohort[manifest_row["cohort"]].append(outcome)
        by_cohort_anchor[
            (manifest_row["cohort"], manifest_row["anchor_curve_id"])
        ].append(outcome)

    pooled_full = by_cohort["A_356_full"] + by_cohort["B_385_full"]
    primary = comparison_record(
        "pooled full 356/385 CRT fingerprint versus anchor-matched ordinary controls",
        pooled_full,
        by_cohort["C_matched_ordinary"],
    )
    body = {
        "schema": "elkies-k3.r17-prospective-crt-half-lattice-analysis.v3",
        "status": "COMPLETE_FROZEN_HALF_LATTICE_DETECTOR_ANALYSIS",
        "candidate_list_sha256": EXPECTED_CANDIDATE_HASH,
        "protocol_definition_sha256": protocol["protocol_definition_sha256"],
        "analysis": {
            "primary_comparison": primary,
            "anchor_specific_comparisons": [
                comparison_record(
                    "356 full fingerprint versus 356-anchor ordinary controls",
                    by_cohort["A_356_full"],
                    by_cohort_anchor[("C_matched_ordinary", 356)],
                ),
                comparison_record(
                    "385 full fingerprint versus 385-anchor ordinary controls",
                    by_cohort["B_385_full"],
                    by_cohort_anchor[("C_matched_ordinary", 385)],
                ),
            ],
            "cohort_outcomes": {
                cohort: outcome_summary(rows) for cohort, rows in sorted(by_cohort.items())
            },
            "fixed_ablation_order": [
                {
                    "order": order,
                    "cohort": cohort,
                    **outcome_summary(
                        pooled_full if cohort == "FULL_POOLED_A_PLUS_B" else by_cohort[cohort]
                    ),
                }
                for order, cohort in enumerate(
                    (
                        "F_random_equal_codimension",
                        "D_two_only",
                        "E_odd_only",
                        "FULL_POOLED_A_PLUS_B",
                    ),
                    start=1,
                )
            ],
            "stage_b_conditional_recovery": {
                "gate": "one or more exactly certified Stage-A quotient directions",
                "executed_rows": sum(
                    row.get("stage_b", {}).get("full_specialized_ranking_computed", False)
                    for row in ledger["records"]
                ),
                "incremental_certified_directions": sum(
                    row.get("stage_b", {}).get("incremental_certified_quotient_gain", 0)
                    for row in ledger["records"]
                ),
                "not_an_unconditional_cohort_response": True,
            },
            "censoring_audit": {
                "censored_rows": sum(
                    not row["analysis_eligible_complete_stage_a"] for row in ledger["records"]
                ),
                "primary_denominator_keeps_all_scheduled_rows": True,
                "exact_events_on_partially_failed_rows_still_count": True,
                "complete_case_results_are_sensitivity_only": True,
            },
        },
        "claim_boundary": [
            "The primary response is detector-visible Stage-A escape, not the true rank-jump incidence.",
            "Every event is an exact point with finite-reduction nonmembership and independence certification.",
            "A bounded miss is not rank 17 and gives no Selmer upper bound.",
            "Stage B is outcome-gated and is reported only as conditional incremental recovery.",
            "No cohorts, masks, bounds, gates, or contrasts are retuned after outcomes are opened.",
        ],
    }
    return {
        **body,
        "analysis_definition_sha256": canonical_hash(body),
        "inputs": {
            relative(MANIFEST): digest(MANIFEST),
            relative(PROTOCOL): digest(PROTOCOL),
            relative(LEDGER): digest(LEDGER),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": "python3 elkies-k3/scripts/analyze_r17_prospective_crt_half_lattice_experiment.py",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored half-lattice analysis differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    primary = document["analysis"]["primary_comparison"]
    print(
        "R17CRTHALFANALYSIS"
        f"|full={primary['exposed']['stage_a_certified_escape_rows']}/{primary['exposed']['scheduled']}"
        f"|ordinary={primary['control']['stage_a_certified_escape_rows']}/{primary['control']['scheduled']}"
        f"|hash={document['analysis_definition_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
