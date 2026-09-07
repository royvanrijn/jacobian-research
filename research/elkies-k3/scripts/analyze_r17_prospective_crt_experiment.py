#!/usr/bin/env python3
"""Analyze the frozen prospective R17 CRT experiment without retuning it."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
from hashlib import sha256
from math import comb, sqrt
import json
from pathlib import Path
from statistics import median
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
FEATURES = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz"
PROTOCOL = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-search-protocol-v2.json"
POINTS = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-point-search-ledger-v2.json"
SENSITIVITY = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-search-sensitivity-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-statistical-analysis-v1.json"

SCHEMA = "elkies-k3.r17-prospective-crt-statistical-analysis.v1"
EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
EXPECTED_PROTOCOL_HASH = "63d6b9e83f52bc7208b9057298e05941dfcedc85d53f5681186c953498947d4b"
Z_975 = 1.959963984540054


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text())


def load_gzip_json(path: Path):
    with gzip.open(path, "rt") as source:
        return json.load(source)


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


def clopper_pearson_zero_event_upper(total: int, alpha: float = 0.05):
    if total == 0:
        return None
    return 1 - (alpha / 2) ** (1 / total)


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
    return float(sum((probability(value) for value in range(lower, upper + 1)
                      if probability(value) <= observed), Fraction(0)))


def outcome_summary(rows: Iterable[dict]):
    rows = list(rows)
    status_counts = Counter(row["status"] for row in rows)
    positives = sum(row["status"] == "CERTIFIED_MW17_ESCAPE" for row in rows)
    completed = sum(not row["status"].startswith("CENSORED_") for row in rows)
    censored = len(rows) - completed
    extra_directions = sum(row.get("certified_independent_extra_directions", 0) for row in rows)
    return {
        "scheduled": len(rows),
        "completed": completed,
        "censored": censored,
        "status_counts": dict(sorted(status_counts.items())),
        "certified_escape_rows": positives,
        "certified_extra_directions": extra_directions,
        "intention_to_search_escape_rate": ratio_record(positives, len(rows)),
        "complete_case_escape_rate": ratio_record(positives, completed),
        "wilson_95_percent_interval_complete_case": wilson_interval(positives, completed),
        "two_sided_95_percent_clopper_pearson_upper_if_zero_events": (
            clopper_pearson_zero_event_upper(completed) if positives == 0 else None
        ),
        "largest_certified_rank_lower_bound": max(
            (row.get("largest_certified_rank_lower_bound", 17) for row in rows),
            default=None,
        ),
    }


def comparison_record(label: str, exposed: list[dict], control: list[dict]):
    exposed_summary = outcome_summary(exposed)
    control_summary = outcome_summary(control)
    exposed_events = exposed_summary["certified_escape_rows"]
    control_events = control_summary["certified_escape_rows"]
    exposed_total = exposed_summary["completed"]
    control_total = control_summary["completed"]
    exposed_rate = exposed_events / exposed_total if exposed_total else None
    control_rate = control_events / control_total if control_total else None
    risk_difference = (
        exposed_rate - control_rate
        if exposed_rate is not None and control_rate is not None
        else None
    )
    risk_ratio = (
        exposed_rate / control_rate
        if exposed_rate is not None and control_rate not in (None, 0)
        else None
    )
    odds_ratio = (
        exposed_events * (control_total - control_events)
        / ((exposed_total - exposed_events) * control_events)
        if exposed_total > exposed_events and control_total > control_events and control_events
        else None
    )
    fisher = (
        fisher_two_sided(
            exposed_events,
            exposed_total - exposed_events,
            control_events,
            control_total - control_events,
        )
        if exposed_total and control_total
        else None
    )
    exposed_upper = (
        clopper_pearson_zero_event_upper(exposed_total) if exposed_events == 0 else None
    )
    control_upper = (
        clopper_pearson_zero_event_upper(control_total) if control_events == 0 else None
    )
    return {
        "contrast": label,
        "exposed": exposed_summary,
        "control": control_summary,
        "risk_difference": risk_difference,
        "risk_ratio": risk_ratio,
        "odds_ratio": odds_ratio,
        "fisher_exact_two_sided_p": fisher,
        "zero_event_cartesian_95_percent_risk_difference_bounds": (
            [-control_upper, exposed_upper]
            if exposed_upper is not None and control_upper is not None
            else None
        ),
        "interpretation": (
            "NO_EVENTS_EITHER_ARM_EFFECT_RATIOS_UNDEFINED"
            if exposed_events == control_events == 0
            else "ESTIMABLE_FROM_COMPLETED_ROWS"
        ),
    }


def distribution(values: Iterable[int | str]):
    return dict(sorted(Counter(str(value) for value in values).items()))


def exact_mean(values: list[int]):
    if not values:
        return None
    return ratio_record(sum(values), len(values))


def build():
    manifest = load_json(MANIFEST)
    features = load_gzip_json(FEATURES)
    protocol = load_json(PROTOCOL)
    points = load_json(POINTS)
    sensitivity = load_json(SENSITIVITY)
    for document in (manifest, features, points):
        candidate_hash = document.get("candidate_list_sha256") or document.get("commitment", {}).get(
            "candidate_list_sha256"
        )
        if candidate_hash != EXPECTED_CANDIDATE_HASH:
            raise ArithmeticError("an analysis input names another candidate list")
    if protocol.get("protocol_definition_sha256") != EXPECTED_PROTOCOL_HASH:
        raise ArithmeticError("the reviewed amended search protocol changed")
    if points.get("search_protocol_sha256") != EXPECTED_PROTOCOL_HASH:
        raise ArithmeticError("the point ledger used another search protocol")
    if sensitivity.get("search_protocol_sha256") != EXPECTED_PROTOCOL_HASH:
        raise ArithmeticError("the sensitivity audit used another search protocol")
    if points.get("status") != "COMPLETE_FROZEN_BOUNDED_POINT_SEARCH_LEDGER":
        raise ArithmeticError("the point ledger is incomplete")

    manifest_by_id = {row["sample_id"]: row for row in manifest["rows"]}
    feature_by_id = {row["sample_id"]: row for row in features["rows"]}
    point_by_id = {row["sample_id"]: row for row in points["records"]}
    if not (set(manifest_by_id) == set(feature_by_id) == set(point_by_id)):
        raise ArithmeticError("the frozen manifest, features, and outcomes do not cover the same rows")
    if len(manifest_by_id) != len(manifest["rows"]):
        raise ArithmeticError("duplicate sample id in frozen manifest")

    by_cohort = defaultdict(list)
    by_cohort_anchor = defaultdict(list)
    feature_by_cohort = defaultdict(list)
    for sample_id in manifest_by_id:
        manifest_row = manifest_by_id[sample_id]
        point_row = point_by_id[sample_id]
        feature_row = feature_by_id[sample_id]
        if (
            manifest_row["cohort"] != point_row["cohort"]
            or manifest_row["cohort"] != feature_row["cohort"]
            or manifest_row["parameter"] != point_row["parameter"]
            or manifest_row["parameter"] != feature_row["parameter"]
        ):
            raise ArithmeticError("a joined row changed cohort or parameter")
        by_cohort[manifest_row["cohort"]].append(point_row)
        by_cohort_anchor[(manifest_row["cohort"], manifest_row["anchor_curve_id"])].append(point_row)
        feature_by_cohort[manifest_row["cohort"]].append(feature_row)

    full = by_cohort["A_356_full"] + by_cohort["B_385_full"]
    primary = comparison_record(
        "pooled full 356/385 fingerprint versus anchor-matched ordinary controls",
        full,
        by_cohort["C_matched_ordinary"],
    )
    anchor_specific = [
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
    ]
    ablation_order = [
        "F_random_equal_codimension",
        "D_two_only",
        "E_odd_only",
        "FULL_POOLED_A_PLUS_B",
    ]
    ablation_rows = {
        "F_random_equal_codimension": by_cohort["F_random_equal_codimension"],
        "D_two_only": by_cohort["D_two_only"],
        "E_odd_only": by_cohort["E_odd_only"],
        "FULL_POOLED_A_PLUS_B": full,
    }
    ablation = [
        {"order": index, "cohort": cohort, **outcome_summary(ablation_rows[cohort])}
        for index, cohort in enumerate(ablation_order, start=1)
    ]

    presearch = {}
    for cohort, rows in sorted(feature_by_cohort.items()):
        stacked_ranks = [
            row["localization_intersections"]["full_stacked_localization_rank"]
            for row in rows
        ]
        simultaneous_kernels = [
            row["localization_intersections"]["full_simultaneous_source_kernel_dimension"]
            for row in rows
        ]
        nagao_scores = [row["nagao_comparison"]["total_score_units_1e12"] for row in rows]
        fingerprint_configuration_hashes = [
            canonical_hash(
                [
                    [local["rational_prime"], local["comparison_sha256"]]
                    for local in row["local_fingerprints"]
                ]
            )
            for row in rows
        ]
        presearch[cohort] = {
            "rows": len(rows),
            "matched_intended_place_count_distribution": distribution(
                row["anchor_fingerprint_survival"]["matched_intended_place_count"]
                for row in rows
            ),
            "full_stacked_known_MW17_localization_rank_distribution": distribution(stacked_ranks),
            "full_simultaneous_known_MW17_source_kernel_dimension_distribution": distribution(
                simultaneous_kernels
            ),
            "distinct_full_local_fingerprint_configurations": len(set(fingerprint_configuration_hashes)),
            "nagao_total_score_units_1e12": {
                "minimum": min(nagao_scores),
                "median": median(nagao_scores),
                "maximum": max(nagao_scores),
                "exact_mean": exact_mean(nagao_scores),
                "status": "HEURISTIC_COMPARISON_FEATURE_ONLY_NOT_USED_FOR_SELECTION",
            },
        }

    all_rows = list(point_by_id.values())
    analysis = {
        "primary_comparison": primary,
        "anchor_specific_comparisons": anchor_specific,
        "cohort_outcomes": {
            cohort: outcome_summary(rows) for cohort, rows in sorted(by_cohort.items())
        },
        "ablation_order": ablation,
        "presearch_arithmetic_feature_summary": presearch,
        "zero_event_diagnostics": {
            "all_candidates": outcome_summary(all_rows),
            "post_experiment_positive_control_sensitivity": {
                "status": sensitivity["status"],
                "controls": [
                    {
                        "curve_id": row["curve_id"],
                        "historically_certified_residual_MW_dimension": row[
                            "historically_certified_residual_MW_dimension"
                        ],
                        "certified_escape_count_at_frozen_bound": row[
                            "certified_escape_count"
                        ],
                    }
                    for row in sensitivity["controls"]
                ],
                "changes_frozen_experiment": False,
            },
            "meaning": (
                "The frozen height-10000 protocol produced no event variance and subsequently "
                "failed to rediscover an escape on either known +12 positive-control fibre. It "
                "therefore cannot estimate relative enrichment or support fitting an outcome predictor."
            ),
        },
        "predictor_fit": {
            "status": "NOT_FIT_ZERO_OUTCOME_VARIANCE",
            "features_inspected_for_retuning": False,
            "reason": "All 2560 frozen outcomes are completed bounded misses.",
        },
        "chain_assessment": {
            "parameter_congruences_to_stable_local_conditions": (
                "SUPPORTED_ON_THE_FROZEN_FINITE_SAMPLE_AFTER_PREOUTCOME_EXPONENT_REFINEMENT"
            ),
            "stable_local_conditions_to_unusually_large_residual_global_space": (
                "NOT_TESTED_NO_COMPLETE_SELMER_GROUP_OR_FINITE_PROVED_RESIDUAL_UPPER_BOUND"
            ),
            "residual_global_space_to_increased_new_rational_directions": (
                "NOT_RESOLVED_ZERO_EVENTS_AND_FAILED_POSTHOC_POSITIVE_CONTROL_SENSITIVITY"
            ),
            "prospective_enrichment_claim": "NO_EVIDENCE_DETECTOR_LIMITED_AT_THE_FROZEN_SEARCH_BOUND",
        },
    }
    body = {
        "schema": SCHEMA,
        "status": "COMPLETE_ZERO_EVENT_DETECTOR_LIMITED_FROZEN_EXPERIMENT_ANALYSIS",
        "candidate_list_sha256": EXPECTED_CANDIDATE_HASH,
        "search_protocol_sha256": EXPECTED_PROTOCOL_HASH,
        "analysis": analysis,
        "claim_boundary": [
            "Zero bounded-search events provide no evidence of relative enrichment at this search bound.",
            "A completed bounded miss is not rank 17 and gives no Selmer upper bound.",
            "The refined local cylinders are empirically stable on the sampled rows, not proved constant cylinders.",
            "Known-MW17 localization intersections are not residual-Selmer dimensions.",
            "No predictor is fit because the frozen outcome has zero variance.",
            "The post-experiment positive-control sensitivity failure does not alter the frozen protocol or rescue an enrichment claim.",
            "No cohort is extended or rebalanced after outcomes were opened.",
        ],
    }
    return {
        **body,
        "analysis_definition_sha256": canonical_hash(body),
        "inputs": {
            relative(MANIFEST): digest(MANIFEST),
            relative(FEATURES): digest(FEATURES),
            relative(PROTOCOL): digest(PROTOCOL),
            relative(POINTS): digest(POINTS),
            relative(SENSITIVITY): digest(SENSITIVITY),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": "python3 elkies-k3/scripts/analyze_r17_prospective_crt_experiment.py",
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
            raise ArithmeticError("stored prospective CRT analysis differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    primary = document["analysis"]["primary_comparison"]
    print(
        "R17CRTANALYSIS"
        f"|full={primary['exposed']['certified_escape_rows']}/{primary['exposed']['completed']}"
        f"|ordinary={primary['control']['certified_escape_rows']}/{primary['control']['completed']}"
        f"|hash={document['analysis_definition_sha256']}"
        f"|status={document['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
