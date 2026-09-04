#!/usr/bin/env python3
"""Freeze the pre-outcome backend-feasibility amendment for the CRT search."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
FEATURES = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-search-protocol-v2.json"
EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_hash(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build():
    manifest = json.loads(MANIFEST.read_text())
    if manifest["commitment"]["candidate_list_sha256"] != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the reviewed candidate list changed")
    if any(row["outcome_status"] != "NOT_OPENED" for row in manifest["rows"]):
        raise ArithmeticError("a Mordell--Weil outcome was opened before the amendment")
    protocol = {
        "schema": "elkies-k3.r17-prospective-crt-search-protocol.v2",
        "candidate_list_sha256": EXPECTED_CANDIDATE_HASH,
        "frozen_before_any_point_search_call_completed_or_returned_points": True,
        "selection_or_rebalancing_changed": False,
        "reason_for_amendment": (
            "The frozen v1 Sage/eclib canary spent the complete 300-second envelope "
            "inside mwrank_EllipticCurve initialization and never reached the search call."
        ),
        "infeasible_v1_canary": {
            "sample_id": manifest["rows"][0]["sample_id"],
            "cohort": manifest["rows"][0]["cohort"],
            "declared_wall_seconds": 300,
            "last_completed_phase": "exact specialization and construction of the integral p=2-minimal model",
            "last_entered_phase": "eclib initialization",
            "point_search_call_reached": False,
            "points_returned_or_inspected": False,
            "mathematical_outcome": "NOT_OBSERVED_BACKEND_FEASIBILITY_ONLY",
        },
        "amended_uniform_bounded_search": {
            "protocol_id": "r17-prospective-crt-direct-rational-points-v2",
            "engine": "PARI hyperellratpoints",
            "search_model": (
                "For the deterministic exact integral p=2-minimal model "
                "y^2+a1*x*y+a3*y=x^3+a2*x^2+a4*x+a6, search the exact completed-square "
                "cubic Y^2=(2y+a1*x+a3)^2."
            ),
            "x_numerator_denominator_height": 10_000,
            "wall_clock_limit_seconds_including_setup": 30,
            "memory_limit_bytes": 8_000_000_000,
            "retries": 0,
            "adaptive_stopping": False,
            "same_algorithm_and_limits_for_every_row": True,
            "monotone_residual_selmer_gate_required": True,
            "finite_quotient_certificate_prime_bound": 1000,
            "positive_counting_rule": (
                "Count a direction only after exact point transport/equation verification and "
                "a full mod-2 finite-quotient independence certificate for generic MW17 plus "
                "all counted directions."
            ),
        },
        "outcome_labels": {
            "positive": "CERTIFIED_MW17_ESCAPE",
            "completed_miss": "BOUNDED_PROTOCOL_NO_ESCAPE_FOUND",
            "uncertified_point": "COMPLETED_UNCERTIFIED_CANDIDATES_NO_COUNTED_ESCAPE",
            "timeout": "CENSORED_TIMEOUT",
            "backend_failure": "CENSORED_BACKEND_FAILURE",
        },
        "claim_boundary": [
            "The v1 canary is an operational failure, not a bounded point-search miss.",
            "The amendment was frozen without seeing a returned point or completed search outcome.",
            "The direct height-10000 search is bounded and cannot prove rank 17 or a Selmer upper bound.",
            "Candidate lists, cohorts, matching, local fingerprints, and analysis contrasts are unchanged.",
        ],
    }
    return {
        **protocol,
        "protocol_definition_sha256": canonical_hash(protocol),
        "inputs": {relative(MANIFEST): digest(MANIFEST), relative(FEATURES): digest(FEATURES)},
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": "python3 elkies-k3/scripts/build_r17_prospective_crt_search_protocol.py",
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
            raise ArithmeticError("stored amended search protocol differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17CRTSEARCHPROTOCOL"
        f"|hash={document['protocol_definition_sha256']}"
        "|status=FROZEN_BEFORE_POINT_SEARCH",
        flush=True,
    )


if __name__ == "__main__":
    main()
