#!/usr/bin/env python3
"""Freeze the prospective zero-gain rescue arm for MW17-jump-v2.

The source campaign rows and their order are reused byte-for-byte.  Treatment
assignment is a deterministic hash bucket and reads no detector outcome.  No
parameter census, specialization, descent, or point search occurs here.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CAMPAIGN = GENERATED / "elkies-k3-mw17-jump-v2-campaign-v1.json"
KNOWN_LEDGER = GENERATED / "elkies-k3-mw17-jump-v2-ledger-v1.json"
RUNNER = ROOT / "elliptic-curves/cas/run_mw17_jump_v2_zero_gain_rescue.sage"
BASE_RUNNER = ROOT / "elliptic-curves/cas/run_mw17_jump_v2.sage"
PRODUCTION_GATES = ROOT / "elliptic-curves/cas/production_search_gates.py"
OUTPUT = GENERATED / "elkies-k3-mw17-jump-v2-zero-gain-rescue-arm-v1.json"

SCHEMA = "elkies-k3.mw17-jump-v2-zero-gain-rescue-arm.v1"
STATUS = "FROZEN_PROSPECTIVE_ZERO_GAIN_RESCUE_BEFORE_RESCUE_SEARCH"
CAMPAIGN_SCHEMA = "elkies-k3.mw17-jump-v2-campaign.v1"
CAMPAIGN_STATUS = "FROZEN_SELECTED_IMMUTABLE_POPULATIONS_BEFORE_V2_EVALUATION"
ASSIGNMENT_NAMESPACE = "mw17-jump-v2-zero-gain-rescue-v1"
ASSIGNMENT_DENOMINATOR = 8
ASSIGNMENT_BUCKET = 0


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def assignment_hash(sample_id: str) -> str:
    return sha256(f"{ASSIGNMENT_NAMESPACE}\0{sample_id}".encode()).hexdigest()


def assignment_bucket(sample_id: str) -> int:
    return int(assignment_hash(sample_id), 16) % ASSIGNMENT_DENOMINATOR


def build() -> dict[str, Any]:
    campaign = json.loads(CAMPAIGN.read_text())
    if (
        campaign.get("schema") != CAMPAIGN_SCHEMA
        or campaign.get("status") != CAMPAIGN_STATUS
        or campaign.get("candidate_count") != 2_239
    ):
        raise ArithmeticError("the immutable MW17-jump-v2 campaign changed")
    if campaign.get("candidate_list_sha256") != canonical_hash(campaign["rows"]):
        raise ArithmeticError("the source candidate-list commitment does not replay")

    assignments = []
    for row in campaign["rows"]:
        bucket = assignment_bucket(row["sample_id"])
        assignments.append(
            {
                "campaign_index": row["campaign_index"],
                "sample_id": row["sample_id"],
                "source_population": row["source_population"],
                "family": row["family"],
                "priority_tranche": row["priority_tranche"],
                "assignment_sha256": assignment_hash(row["sample_id"]),
                "assignment_bucket_mod_8": bucket,
                "assigned_to_rescue_arm": bucket == ASSIGNMENT_BUCKET,
            }
        )
    if [row["sample_id"] for row in assignments] != [
        row["sample_id"] for row in campaign["rows"]
    ]:
        raise ArithmeticError("the rescue assignment changed the source population")

    treated = [row for row in assignments if row["assigned_to_rescue_arm"]]
    known_context = None
    if KNOWN_LEDGER.exists():
        ledger = json.loads(KNOWN_LEDGER.read_text())
        known_context = {
            "ledger": relative(KNOWN_LEDGER),
            "ledger_sha256_at_protocol_freeze": digest(KNOWN_LEDGER),
            "completed_worker_count_at_protocol_freeze": ledger.get(
                "completed_worker_count"
            ),
            "measured_gain_histogram_at_protocol_freeze": dict(
                sorted(
                    Counter(
                        str(row["actual_certified_quotient_rank_gain"])
                        for row in ledger.get(
                            "measurements_ranked_only_by_actual_certified_quotient_rank_gain",
                            [],
                        )
                    ).items()
                )
            ),
            "used_for_assignment": False,
        }

    document: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "campaign_name": "MW17-jump-v2 zero-gain rescue arm",
        "source_campaign": {
            "path": relative(CAMPAIGN),
            "sha256": digest(CAMPAIGN),
            "candidate_count": campaign["candidate_count"],
            "candidate_list_sha256": campaign["candidate_list_sha256"],
            "rows_and_order_reused_without_change": True,
        },
        "assignment": {
            "namespace": ASSIGNMENT_NAMESPACE,
            "rule": "SHA256(namespace || NUL || sample_id) interpreted as an integer modulo 8",
            "rescue_bucket": ASSIGNMENT_BUCKET,
            "denominator": ASSIGNMENT_DENOMINATOR,
            "prospective_treatment_probability": "1/8",
            "uses_detector_outcomes": False,
            "uses_rank_selmer_point_or_public_control_data": False,
            "assigned_candidate_count": len(treated),
            "assigned_counts_by_source_population": dict(
                sorted(Counter(row["source_population"] for row in treated).items())
            ),
            "assigned_counts_by_family": dict(
                sorted(Counter(row["family"] for row in treated).items())
            ),
        },
        "assignments": assignments,
        "eligibility_after_assignment": {
            "required_base_status": "PASS_EXACT_CERTIFIED_QUOTIENT_GAIN",
            "required_initial_43_chart_certified_gain": 0,
            "required_base_timeout_count": 0,
            "required_base_backend_failure_count": 0,
            "required_base_attempted_chart_count": 43,
            "unmeasured_or_censored_base_rows": "not eligible until a clean initial result exists",
            "eligibility_does_not_change_assignment": True,
        },
        "rescue_detector": {
            "generic_rank": 17,
            "base_initial_chart_count": 43,
            "additional_budget_chart_count": 301,
            "maximum_total_chart_count_including_base": 344,
            "generic_rescue_class_ranks_one_based": [44, 344],
            "generic_rescue_batch_count": 7,
            "generic_rescue_charts_per_batch": 43,
            "generic_rescue_policy": (
                "search the next 301 classes in the already-computed exact generic "
                "half-class order, reranked within that fixed set by specialized "
                "canonical depth; no discovered quotient coordinate is required"
            ),
            "switch_rule": (
                "after each complete 43-chart rescue batch, classify all returned points "
                "exactly; on the first certified quotient escape, assign every unused slot "
                "from the 301-chart rescue allowance to the existing adaptive quotient policy"
            ),
            "height_bound_each_chart": 100_000,
            "wall_timeout_seconds_each_chart": 15.0,
            "gp_stack_bytes_each_chart": 1_000_000_000,
            "relation_chunk_size": 64,
            "relation_timeout_seconds": 180.0,
            "ranking_field": "actual_certified_quotient_rank_gain",
            "rank_lower_bound_formula": "17 + actual_certified_quotient_rank_gain",
            "global_termination": "write stop sentinel immediately after any certified gain at least 15",
        },
        "proof_and_budget_policy": {
            "production_policy_schema": "elliptic-curves.production-rank-search-gates.v1",
            "certified_rank_or_selmer_upper_bound_below_32": "mathematical exclusion",
            "incomplete_or_conditional_descent": "scheduling information only",
            "certified_32_independent_points": "unconditional rank-at-least-32 success",
            "complete_descent_required_to_run_this_bounded_rescue": False,
            "experiment_specific_seals_inherited": False,
            "small_field_100_row_laboratory_changed": False,
        },
        "known_context_at_freeze": known_context,
        "implementation_hashes": {
            relative(path): digest(path)
            for path in (Path(__file__), RUNNER, BASE_RUNNER, PRODUCTION_GATES)
        },
        "resource_gate": {
            "checkpoint_unit": "one rescued fibre per chunk rewrite",
            "worker_wall_timeout_seconds": 7_200,
            "worker_address_space_bytes": 24_000_000_000,
            "forbidden_parameter_work": (
                "do not rerun the 121,589,944-parameter census or rerank the source populations"
            ),
        },
        "claim_boundary": [
            "Treatment assignment is frozen independently of zero-gain outcomes; only eligibility depends on a later clean exact zero.",
            "A rescued nonzero gain is an exact rank lower bound after equation and finite-reduction independence checks.",
            "A zero, timeout, or backend failure is bounded search evidence and never a rank upper bound.",
            "The rescue arm changes detector exposure, not the immutable candidate population or exact-gain leaderboard.",
        ],
        "generation": {
            "script": relative(Path(__file__)),
            "command": "python3 elkies-k3/scripts/build_mw17_jump_v2_zero_gain_rescue.py",
        },
    }
    definition = {
        key: value
        for key, value in document.items()
        if key != "protocol_definition_sha256"
    }
    document["protocol_definition_sha256"] = canonical_hash(definition)
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = build()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != text:
            raise SystemExit("MW17-jump-v2 zero-gain rescue artifact is stale")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(
        "MW17JUMPV2RESCUE|status=FROZEN"
        f"|assigned={payload['assignment']['assigned_candidate_count']}"
        f"|sha256={digest(args.output)}|output={relative(args.output)}"
    )


if __name__ == "__main__":
    main()
