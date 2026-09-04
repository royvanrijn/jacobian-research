#!/usr/bin/env python3
"""Analyze the frozen v3 detector with a fail-closed censoring gate."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
LEGACY = ROOT / "elkies-k3/scripts/analyze_r17_prospective_crt_half_lattice_experiment.py"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-prospective-crt-half-lattice-analysis-v4.json"
)


def load_legacy():
    spec = importlib.util.spec_from_file_location("r17_half_lattice_v3_analyzer", LEGACY)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy = load_legacy()


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def ratio_record(numerator: int, denominator: int):
    if denominator == 0:
        return None
    value = Fraction(numerator, denominator)
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def outcome_summary(rows: Iterable[dict[str, Any]]):
    rows = list(rows)
    eligible = [row for row in rows if row["analysis_eligible_complete_stage_a"]]
    events = [
        row
        for row in rows
        if row.get("stage_a", {}).get("certified_quotient_gain", 0) > 0
    ]
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
        "observed_event_lower_bound_per_scheduled_fibre": ratio_record(
            len(events), len(rows)
        ),
        "scheduled_denominator_is_not_an_inferential_event_rate": True,
        "complete_case_escape_rows": len(eligible_events),
        "complete_case_escape_rate": ratio_record(len(eligible_events), len(eligible)),
        "complete_case_wilson_95_percent_interval": legacy.wilson_interval(
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


def censoring_balance(exposed: list[dict], control: list[dict]) -> dict[str, Any]:
    """Require exact equality of every censor-status proportion."""

    exposed_scheduled = len(exposed)
    control_scheduled = len(control)
    exposed_censored = Counter(
        row["status"] for row in exposed if not row["analysis_eligible_complete_stage_a"]
    )
    control_censored = Counter(
        row["status"] for row in control if not row["analysis_eligible_complete_stage_a"]
    )
    statuses = sorted(set(exposed_censored) | set(control_censored))
    status_rows = []
    balanced = exposed_scheduled > 0 and control_scheduled > 0
    for status in statuses:
        exposed_count = exposed_censored[status]
        control_count = control_censored[status]
        status_balanced = (
            exposed_count * control_scheduled
            == control_count * exposed_scheduled
        )
        balanced = balanced and status_balanced
        status_rows.append(
            {
                "status": status,
                "exposed_count": exposed_count,
                "control_count": control_count,
                "exact_proportions_equal": status_balanced,
            }
        )
    return {
        "rule": (
            "each distinct censored Stage-A status must have exactly the same "
            "scheduled-row proportion in both arms"
        ),
        "exposed_scheduled": exposed_scheduled,
        "control_scheduled": control_scheduled,
        "exposed_censored": sum(exposed_censored.values()),
        "control_censored": sum(control_censored.values()),
        "censored_statuses": status_rows,
        "balanced": balanced,
    }


def comparison_record(label: str, exposed: list[dict], control: list[dict]):
    exposed_summary = outcome_summary(exposed)
    control_summary = outcome_summary(control)
    exposed_events = exposed_summary["complete_case_escape_rows"]
    control_events = control_summary["complete_case_escape_rows"]
    exposed_total = exposed_summary["complete_stage_a"]
    control_total = control_summary["complete_stage_a"]
    balance = censoring_balance(exposed, control)
    authorized = balance["balanced"] and exposed_total > 0 and control_total > 0
    exposed_rate = exposed_events / exposed_total if authorized else None
    control_rate = control_events / control_total if authorized else None
    risk_ratio = (
        exposed_rate / control_rate
        if authorized and control_rate not in (None, 0)
        else None
    )
    odds_ratio = (
        exposed_events * (control_total - control_events)
        / ((exposed_total - exposed_events) * control_events)
        if authorized
        and exposed_total > exposed_events
        and control_total > control_events
        and control_events
        else None
    )
    return {
        "contrast": label,
        "estimand": (
            "COMPLETE_CASE_STAGE_A_DETECTOR_YIELD_ONLY_AFTER_EXACT_CENSORING_BALANCE"
        ),
        "censoring_gate": balance,
        "prospective_conclusion_authorized": authorized,
        "exposed": exposed_summary,
        "control": control_summary,
        "risk_difference": exposed_rate - control_rate if authorized else None,
        "risk_ratio": risk_ratio,
        "odds_ratio": odds_ratio,
        "fisher_exact_two_sided_p": (
            legacy.fisher_two_sided(
                exposed_events,
                exposed_total - exposed_events,
                control_events,
                control_total - control_events,
            )
            if authorized
            else None
        ),
        "gated_complete_case_counts": {
            "exposed_events": exposed_events,
            "exposed_total": exposed_total,
            "control_events": control_events,
            "control_total": control_total,
        },
        "interpretation": (
            "ESTIMABLE_COMPLETE_CASE_CONTRAST_AFTER_BALANCED_CENSORING"
            if authorized
            else "GATED_NO_PROSPECTIVE_CONCLUSION_UNBALANCED_OR_EMPTY_CENSORING"
        ),
    }


def build() -> dict[str, Any]:
    legacy.outcome_summary = outcome_summary
    legacy.comparison_record = comparison_record
    old = legacy.build()
    body = {
        key: value
        for key, value in old.items()
        if key not in {"analysis_definition_sha256", "inputs", "generation"}
    }
    body["schema"] = "elkies-k3.r17-prospective-crt-half-lattice-analysis.v4"
    body["status"] = "COMPLETE_FROZEN_HALF_LATTICE_DETECTOR_ANALYSIS_CENSORING_GATED"
    primary = body["analysis"]["primary_comparison"]
    anchor_specific = body["analysis"]["anchor_specific_comparisons"]
    body["analysis"]["prospective_conclusions_gate"] = {
        "primary_balanced_censoring": primary["prospective_conclusion_authorized"],
        "both_anchor_specific_balanced_censoring": all(
            row["prospective_conclusion_authorized"] for row in anchor_specific
        ),
        "authorized": primary["prospective_conclusion_authorized"]
        and all(row["prospective_conclusion_authorized"] for row in anchor_specific),
    }
    body["analysis"]["censoring_audit"] = {
        "censored_rows": sum(
            cohort["censored_stage_a"]
            for cohort in body["analysis"]["cohort_outcomes"].values()
        ),
        "scheduled_rows_are_never_counted_as_completed_non_events": True,
        "exact_events_on_partially_failed_rows_still_count": True,
        "events_on_censored_rows_are_descriptive_lower_bounds_only": True,
        "complete_case_contrasts_require_exact_censor_status_balance": True,
    }
    body["claim_boundary"].extend(
        [
            "No prospective contrast is reported unless every censor-status proportion is exactly balanced between its arms.",
            "Scheduled-denominator observed-event yields are descriptive lower bounds, never event-rate estimates.",
        ]
    )
    return {
        **body,
        "analysis_definition_sha256": canonical_hash(body),
        "inputs": {
            **old["inputs"],
            relative(LEGACY): digest(LEGACY),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                "python3 elkies-k3/scripts/"
                "analyze_r17_prospective_crt_half_lattice_censoring_gated.py"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    encoded = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != encoded:
            raise SystemExit(f"stale or missing censoring-gated analysis: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    primary = document["analysis"]["primary_comparison"]
    counts = primary["gated_complete_case_counts"]
    print(
        "R17CRTHALFCENSOR"
        f"|full={counts['exposed_events']}/{counts['exposed_total']}"
        f"|ordinary={counts['control_events']}/{counts['control_total']}"
        f"|gate={str(primary['prospective_conclusion_authorized']).lower()}"
        f"|hash={document['analysis_definition_sha256']}"
    )


if __name__ == "__main__":
    main()
