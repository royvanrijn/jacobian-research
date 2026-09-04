#!/usr/bin/env python3
"""Freeze the cross-class correction after curve 478 and before curve 498."""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
V1_PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v1.json"
V1_PARTIAL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v1.json"
V1_RUNNER = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage"
V2_RUNNER = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind_v2.sage"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
POLICY = ROOT / "elliptic-curves/cas/half_lattice_chart_policy.py"
VERIFIER = ROOT / "elliptic-curves/cas/verify_r17_refresh_jump_ladder.sage"
ANALYZER = ROOT / "elliptic-curves/cas/analyze_r17_refresh_jump_ladder.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v2.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_hash(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build():
    base = json.loads(V1_PROTOCOL.read_text())
    partial = json.loads(V1_PARTIAL.read_text())
    if base.get("status") != "FROZEN_BEFORE_BLIND_RECOVERY":
        raise ArithmeticError("the v1 protocol was not frozen")
    if len(partial.get("results", [])) != 1:
        raise ArithmeticError("the v1 stop did not occur after exactly one case")
    first = partial["results"][0]
    if (
        int(first["curve_id"]) != 478
        or first.get("status") != "PASS_COMPLETE_EQUAL_BUDGET_BLIND_RECOVERY"
        or int(first["exact_quotient_rank_recovered_before_public_complement"]) != 6
    ):
        raise ArithmeticError("the disclosed v1 curve-478 outcome changed")
    search_policy = deepcopy(base["search_policy"])
    search_policy["policy_id"] = "r17-refresh-adaptive-half-lattice-top43-plus-301-v2"
    search_policy["initial_universe"] = (
        "the first 43 parity classes in the complete exact generic MW17 "
        "half-lattice order: decreasing exact generic norm, then mask"
    )
    paths = [V1_RUNNER, V2_RUNNER, LEGACY, ENGINE, POLICY, VERIFIER, ANALYZER]
    body = {
        **{
            key: deepcopy(value)
            for key, value in base.items()
            if key
            not in {
                "schema",
                "status",
                "search_policy",
                "implementation_hashes",
                "protocol_definition_sha256",
                "generation",
            }
        },
        "schema": "elliptic-curves.r17-refresh-jump-ladder-protocol.v2",
        "status": "FROZEN_BEFORE_BLIND_RECOVERY",
        "freeze_boundary": "AFTER_CURVE478_BEFORE_OTHER_RECOVERY_OUTCOMES",
        "search_policy": search_policy,
        "cross_class_amendment": {
            "trigger": (
                "the complete generic census for curve 498 contradicted the v1 "
                "cross-class assertion before any curve-498 chart was searched"
            ),
            "known_blind_outcomes_at_amendment": [
                {
                    "curve_id": 478,
                    "exact_quotient_rank_recovered_before_public_complement": 6,
                }
            ],
            "change": (
                "replace 'all 43 maximum-norm classes' by 'first 43 classes in "
                "decreasing exact generic norm then mask'"
            ),
            "curve478_selected_set_and_order_unchanged": True,
            "height_time_and_chart_budgets_unchanged": True,
            "confirmatory_endpoints_and_acceptance_rules_unchanged": True,
            "v2_reruns_curve478_instead_of_importing_its_outcome": True,
            "remaining_fifteen_recovery_outcomes_known_at_amendment": False,
        },
        "v1_freeze_and_stop_hashes": {
            relative(V1_PROTOCOL): digest(V1_PROTOCOL),
            relative(V1_PARTIAL): digest(V1_PARTIAL),
        },
        "implementation_hashes": {relative(path): digest(path) for path in paths},
        "claim_boundary": base["claim_boundary"]
        + [
            "The initial-set correction is outcome-independent for the remaining fifteen cases but is transparently post-curve-478.",
            "Curve 478 is rerun from the redacted input under v2; its known response is not imported into the runner.",
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
            raise ArithmeticError("stored v2 jump-ladder protocol changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(
        f"R17JUMPLADDERPROTOCOLV2|status={payload['status']}|"
        f"hash={payload['protocol_definition_sha256']}"
    )


if __name__ == "__main__":
    main()
