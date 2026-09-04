#!/usr/bin/env python3
"""Freeze the independent-restart-budget amendment for curve 385."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
PROTOCOL = ART / "curve385_sparse_quotient_rank32_protocol_v1.json"
PRIMARY = ART / "curve385_sparse_quotient_rank32_primary_v1.json"
POLICY = CAS / "curve385_sparse_restart_policy.py"
RUNNER = CAS / "run_curve385_sparse_quotient_rank32_search_v2.sage"
OUTPUT = ART / "curve385_sparse_restart_budget_v2.json"

EXPECTED_PROTOCOL_SHA256 = (
    "2c9150f50f305b8aa3763590cd5e81c4d7e121f9373177827780789ce472834f"
)
EXPECTED_PRIMARY_SHA256 = (
    "ff1b0a2e8dd29b9a34a3b81cf8db3bed5350d12ee0dcfd9dd3936f226d255d61"
)

sys.path.insert(0, str(CAS))
from curve385_sparse_restart_policy import (  # noqa: E402
    RANK_CHANGING,
    SATURATION_ONLY,
    policy_document,
    simulate_unit_rank_path,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def build() -> dict[str, Any]:
    if digest(PROTOCOL) != EXPECTED_PROTOCOL_SHA256:
        raise ArithmeticError("the frozen v1 sparse protocol changed")
    if digest(PRIMARY) != EXPECTED_PRIMARY_SHA256:
        raise ArithmeticError("the completed v1 primary manifest changed")
    protocol = json.loads(PROTOCOL.read_text())
    primary = json.loads(PRIMARY.read_text())
    old_cap = protocol["protocol_definition"]["point_search_budget"][
        "maximum_lattice_states"
    ]
    if old_cap != 4:
        raise ArithmeticError("the audited combined v1 state cap is no longer four")
    if primary.get("status") != "PASS_COMPLETE_PRIMARY_SPARSE_CAMPAIGN_BOUNDED_NO_GROWTH":
        raise ArithmeticError("the historical v1 primary campaign is not frozen")

    adverse_path = [
        SATURATION_ONLY,
        SATURATION_ONLY,
        RANK_CHANGING,
        RANK_CHANGING,
        RANK_CHANGING,
    ]
    repaired = simulate_unit_rank_path(adverse_path)
    if repaired["status"] != "TARGET_REACHED" or repaired["rank"] != 32:
        raise AssertionError("the independent restart budgets do not reach rank 32")
    old_combined_cap_stops_before_rank32 = len(adverse_path) > old_cap

    definition = {
        "amends_protocol": relative(PROTOCOL),
        "amends_protocol_sha256": EXPECTED_PROTOCOL_SHA256,
        "scope": (
            "future curve-385 sparse rank-32 searches and resumes of v2 "
            "checkpoints; completed v1 evidence remains immutable"
        ),
        "superseded_rule": {
            "field": "point_search_budget.maximum_lattice_states",
            "value": old_cap,
            "reason": (
                "one combined counter charges rank changes and finite-index-only "
                "saturations against the same nearly minimal state allowance"
            ),
        },
        "replacement": policy_document(),
        "regression": {
            "event_sequence": adverse_path,
            "meaning": "two saturation-only changes followed by unit rank gains 29->30->31->32",
            "old_combined_cap_stops_before_final_rank_gain": (
                old_combined_cap_stops_before_rank32
            ),
            "repaired_result": repaired,
        },
        "historical_boundary": {
            "completed_v1_primary_manifest": relative(PRIMARY),
            "completed_v1_primary_manifest_sha256": EXPECTED_PRIMARY_SHA256,
            "v1_no_growth_result_is_unchanged": True,
            "no_new_curve385_search_outcomes_opened": True,
        },
    }
    body = {
        "schema": "elliptic-curves.curve385-sparse-restart-budget.v2",
        "status": "FROZEN_INDEPENDENT_RESTART_BUDGETS_BEFORE_FUTURE_SEARCH",
        "definition": definition,
        "definition_sha256": canonical_hash(definition),
        "inputs": {
            relative(PROTOCOL): digest(PROTOCOL),
            relative(PRIMARY): digest(PRIMARY),
            relative(POLICY): digest(POLICY),
            relative(RUNNER): digest(RUNNER),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "command": (
                "python3 elliptic-curves/cas/"
                "build_curve385_sparse_restart_budget.py"
            ),
        },
    }
    return body


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    encoded = canonical_bytes(build())
    if args.check:
        if not args.output.exists() or args.output.read_bytes() != encoded:
            raise SystemExit(f"stale or missing restart-budget artifact: {args.output}")
        print(f"C385RESTART|status=PASS|sha256={sha256(encoded).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(encoded)
    print(f"C385RESTART|status=WROTE|sha256={sha256(encoded).hexdigest()}")


if __name__ == "__main__":
    main()
