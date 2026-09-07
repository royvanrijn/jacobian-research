#!/usr/bin/env python3
"""Freeze the R17 multi-stratum blind jump-ladder protocol."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
RUNNER = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
POLICY = ROOT / "elliptic-curves/cas/half_lattice_chart_policy.py"
VERIFIER = ROOT / "elliptic-curves/cas/verify_r17_refresh_jump_ladder.sage"
ANALYZER = ROOT / "elliptic-curves/cas/analyze_r17_refresh_jump_ladder.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_hash(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build():
    blind = json.loads(INPUT.read_text())
    if blind.get("status") != "FROZEN_MW17_ONLY_NO_PUBLIC_COMPLEMENT":
        raise ArithmeticError("the MW17-only input is not frozen")
    if blind["case_count"] != 16:
        raise ArithmeticError("the quotient-eligible ladder size changed")
    if blind["pre_search_exclusion"]["curve_id"] != 499:
        raise ArithmeticError("the pre-search quotient-ineligible case changed")
    implementation_paths = [RUNNER, LEGACY, ENGINE, POLICY, VERIFIER, ANALYZER]
    for path in implementation_paths:
        if not path.exists():
            raise FileNotFoundError(path)
    body = {
        "schema": "elliptic-curves.r17-refresh-jump-ladder-protocol.v1",
        "status": "FROZEN_BEFORE_BLIND_RECOVERY",
        "objective": (
            "Test whether exact quotient rank recovered by one fixed adaptive "
            "half-lattice policy is ordinally associated with exact displayed "
            "jump and enriches the predeclared +10/+11/+12 tail."
        ),
        "eligible_curve_ids": [row["curve_id"] for row in blind["cases"]],
        "pre_search_exclusion": blind["pre_search_exclusion"],
        "blind_input_sha256": digest(INPUT),
        "information_boundary": {
            "available_to_runner": [
                "curve id and atlas class",
                "short Weierstrass equation",
                "seventeen specialized generic MW17 points",
                "exact generic MW17 height Gram",
            ],
            "forbidden_until_blind_artifact_is_frozen": [
                "every public point P18 or beyond",
                "displayed subgroup rank",
                "displayed quotient rank or jump label",
                "coordinates of blind discoveries in the displayed public basis",
            ],
            "blind_response": (
                "exact_quotient_rank_recovered_before_public_complement = "
                "rank(<MW17, all returned points>)-17"
            ),
        },
        "search_policy": {
            "policy_id": "r17-refresh-adaptive-half-lattice-43-plus-301-v1",
            "generic_rank": 17,
            "initial_chart_count": 43,
            "initial_universe": (
                "all 43 maximum-norm parity classes in the complete exact "
                "generic MW17 half-lattice census"
            ),
            "initial_representative_and_order": (
                "shortest representatives at rounded specialized canonical-height "
                "scale 10^6; decreasing exact-decimal depth, then mask"
            ),
            "adaptive_trigger": "at least one exactly certified initial quotient direction",
            "adaptive_chart_count": 301,
            "adaptive_pool": (
                "the first 301 generic parity masks by decreasing exact generic "
                "depth then mask, paired cyclically with nonzero quotient words "
                "ordered by Hamming weight then integer value"
            ),
            "adaptive_representative_and_order": (
                "recompute the current canonical-height form, generic coordinates, "
                "deterministic mod-2 quotient complement, every selected shortest "
                "representative, and the complete 301-chart order for that exact "
                "state; order by decreasing exact-decimal depth then mask and word"
            ),
            "state_binding": (
                "bind both orders to the full basis, height Gram, generic-coordinate "
                "rows, quotient-complement rows, universe id, and ordered chart ids"
            ),
            "total_chart_cap_per_fibre": 344,
            "height_bound_each_quartic": 100_000,
            "wall_timeout_seconds_each_quartic": 15.0,
            "gp_stack_bytes_each_quartic": 1_000_000_000,
            "relation_chunk_size": 64,
            "relation_timeout_seconds_each_chunk": 180.0,
            "retries": 0,
            "audit_rounding_scale": 100_000,
            "operative_rounding_scale": 1_000_000,
            "structural_zero_rule": (
                "if the initial exact recovered quotient rank is zero, the adaptive "
                "quotient space is empty and the fixed policy stops after 43 charts"
            ),
        },
        "confirmatory_analysis": {
            "case_count": 16,
            "score": "exact_quotient_rank_recovered_before_public_complement",
            "truth": "exact displayed-subgroup free quotient rank by specialized MW17",
            "ordinal_endpoint": {
                "statistic": "Kendall tau-b between blind score and displayed jump",
                "null": "truth labels are exchangeable across the sixteen fixed cases",
                "p_value": (
                    "exact one-sided permutation probability of concordance-minus-"
                    "discordance at least as large as observed, preserving all ties"
                ),
                "direction": "positive",
                "acceptance": "tau_b >= 0.35 and exact one-sided p <= 0.05",
            },
            "upper_tail_endpoint": {
                "true_tail": "displayed jump >= 10, namely +10/+11/+12",
                "detector_positive": "blind recovered rank >= 10",
                "statistic": (
                    "detector-positive risk difference, tail minus non-tail, "
                    "with one-sided Fisher exact p"
                ),
                "acceptance": (
                    "risk difference >= 0.25 and one-sided Fisher exact p <= 0.05"
                ),
            },
            "joint_decision": {
                "usable_extreme_jump_detector": "both confirmatory endpoints pass",
                "stop_rule": (
                    "if either endpoint fails, do not spend serious rank-32 point-"
                    "search budget on the half-lattice extreme-jump hypothesis"
                ),
            },
            "descriptive_only": [
                "per-stratum score summaries",
                "recovery fractions relative to displayed jumps",
                "initial versus adaptive incremental gains",
                "timeouts and PARI failures under the equal maximum budget",
            ],
        },
        "implementation_hashes": {
            relative(path): digest(path) for path in implementation_paths
        },
        "claim_boundary": [
            "This freeze occurred before the new blind runner was executed.",
            "The score is an exact rank lower bound for the discovered quotient, not a Mordell-Weil rank upper bound.",
            "The public displayed jump is a certified displayed-subgroup quotient, not an assertion of full E(Q).",
            "The sixteen cases are a fixed atlas-refresh panel rather than a random population sample.",
        ],
    }
    return {
        **body,
        "protocol_definition_sha256": canonical_hash(body),
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
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
            raise ArithmeticError("the stored jump-ladder protocol changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(
        f"R17JUMPLADDERPROTOCOL|status=FROZEN|cases={len(payload['eligible_curve_ids'])}|"
        f"hash={payload['protocol_definition_sha256']}"
    )


if __name__ == "__main__":
    main()
