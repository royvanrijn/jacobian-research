#!/usr/bin/env sage-python
"""Run the bounded maximum-depth MW16 search on anonymous sampled fibres."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import lcm
from pathlib import Path
import platform
import shutil
import subprocess
import time


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/elliptic-curves/a1_mw16_target_free_parameter_candidates_h300_v1.json"
LADDER = ROOT / "elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
OUTPUT = ROOT / "artifacts/local/elliptic-curves/pointed-quartic-search/campaigns/run_a1_mw16_target_free_parameter_search.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def write_payload(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def direct_integral_quartic_search(engine, **kwargs):
    """Compatibility name; all active searches use PointedQuarticSearch."""
    from pointed_quartic_search import run_quartic_search
    return run_quartic_search(**kwargs)


def initial_payload(args, inputs, candidate_ids):
    return {
        "schema": "elliptic-curves.a1-mw16-target-free-parameter-search.v1",
        "status": "SEARCHING",
        "source_candidate_count": len(inputs["candidates"]),
        "candidate_ids": candidate_ids,
        "declared_budget": {
            "height_bound_each_quartic": args.height_bound,
            "timeout_seconds_each_quartic": args.timeout_seconds,
            "stack_bytes_each_quartic": args.stack_bytes,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "chart_rule": (
                "complete exact maximum-depth stratum of the specialized generic "
                "MW16 half-lattice, recomputed separately for each sampled fibre; "
                "PARI searches the exact denominator-cleared quartic directly"
            ),
            "adaptive_quotient_lifts": 0,
            "unrestricted_point_search": False,
        },
        "results": [],
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, LADDER, LEGACY, ENGINE, Path(__file__))
        },
        "software": {"python": platform.python_version()},
        "next_gate": {
            "positive_definition": "exact_quotient_rank_recovered >= 1",
            "action": (
                "compute a minimal model, transport MW16 and recovered points, "
                "then run the complete residual 2-Selmer group on that fibre"
            ),
            "expensive_continuation_authorized_before_complete_selmer": False,
        },
        "claim_boundary": [
            "The search consumes only the anonymous target-free candidate ledger.",
            "Every exact positive is checked by rational group law and finite-reduction independence certificates.",
            "Every miss is bounded and gives no rank upper bound or point-absence claim.",
            "A timed-out chart is censored rather than counted as a negative search.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/run_a1_mw16_target_free_parameter_search.sage"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--maximum-candidates", type=int, default=0)
    parser.add_argument("--candidate-start", type=int, default=0)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic-search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")
    if args.maximum_candidates < 0 or args.candidate_start < 0:
        raise SystemExit("candidate slice bounds must be nonnegative")
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")
    args.exact_generic_order_only = True

    inputs = json.loads(args.input.read_text())
    if inputs.get("status") != "PASS_TARGET_FREE_A1_MW16_PARAMETER_CANDIDATES":
        raise ArithmeticError("target-free parameter ledger is not passing")
    if inputs["next_gate"]["stage"] != "bounded_half_lattice_jump_recovery":
        raise ArithmeticError("candidate ledger does not authorize this search")
    source_candidates = inputs["candidates"]
    candidates = source_candidates[args.candidate_start :]
    if args.maximum_candidates:
        candidates = candidates[: args.maximum_candidates]
    if not candidates:
        raise SystemExit("the requested candidate slice is empty")
    candidate_ids = [row["candidate_id"] for row in candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ArithmeticError("candidate identifiers are not unique")

    fresh = initial_payload(args, inputs, candidate_ids)
    if args.output.is_file():
        payload = json.loads(args.output.read_text())
        for key in ("schema", "source_candidate_count", "candidate_ids", "declared_budget", "inputs"):
            if payload.get(key) != fresh.get(key):
                raise ArithmeticError(f"existing checkpoint differs in immutable field {key}")
    else:
        payload = fresh
        write_payload(args.output, payload)
    completed = {row["candidate_id"] for row in payload["results"]}
    if not completed.issubset(set(candidate_ids)):
        raise ArithmeticError("checkpoint contains a candidate outside the frozen slice")

    ladder = SourceFileLoader("a1_mw16_target_free_ladder", str(LADDER)).load_module()
    legacy = SourceFileLoader("a1_mw16_target_free_legacy", str(LEGACY)).load_module()
    legacy.GENERIC_DIMENSION = 16
    legacy.engine.run_quartic_search = lambda **kwargs: direct_integral_quartic_search(
        legacy.engine, **kwargs
    )
    for index, candidate in enumerate(candidates, 1):
        if candidate["candidate_id"] in completed:
            print(
                f"A1MW16SEARCH|candidate={candidate['candidate_id']}|"
                f"index={index}/{len(candidates)}|status=RESUME_SKIP",
                flush=True,
            )
            continue
        # run_parent's historical field names are internal adapter keys only;
        # the anonymous candidate ledger and emitted result contain no record IDs.
        adapted = {
            "parent_id": candidate["candidate_id"],
            "curve_id": 0,
            "priority_rank": 0,
            "target_short_model": candidate["raw_short_model"],
            "specialized_generic_points": candidate["raw_generic_points"],
            "generic_height_gram": candidate["generic_height_gram"],
        }
        try:
            result = ladder.run_parent(adapted, legacy, args)
            result["candidate_id"] = result.pop("parent_id")
            result.pop("curve_id")
            result.pop("priority_rank")
            result["presentation_id"] = candidate["presentation_id"]
            result["fibration_id"] = candidate["fibration_id"]
            result["parameter"] = candidate["parameter"]
            result["q_isomorphism_class_id"] = candidate["q_isomorphism_class_id"]
            result["nagao"] = candidate["nagao"]
        except ArithmeticError as error:
            result = {
                "candidate_id": candidate["candidate_id"],
                "presentation_id": candidate["presentation_id"],
                "fibration_id": candidate["fibration_id"],
                "parameter": candidate["parameter"],
                "q_isomorphism_class_id": candidate["q_isomorphism_class_id"],
                "nagao": candidate["nagao"],
                "status": "REJECTED_FAIL_CLOSED_BEFORE_EXACT_QUOTIENT_RANK",
                "reason": str(error),
                "exact_quotient_rank_recovered": None,
            }
        payload["results"].append(result)
        write_payload(args.output, payload)
        print(
            f"A1MW16SEARCH|candidate={candidate['candidate_id']}|"
            f"index={index}/{len(candidates)}|"
            f"quotient_rank={result['exact_quotient_rank_recovered']}|"
            f"status={result['status']}",
            flush=True,
        )

    payload["positive_candidate_ids"] = [
        row["candidate_id"]
        for row in payload["results"]
        if row["exact_quotient_rank_recovered"] is not None
        and int(row["exact_quotient_rank_recovered"]) >= 1
    ]
    payload["positive_candidate_count"] = len(payload["positive_candidate_ids"])
    payload["completed_candidate_count"] = len(payload["results"])
    payload["failed_closed_candidate_count"] = sum(
        row["exact_quotient_rank_recovered"] is None for row in payload["results"]
    )
    payload["status"] = (
        "PASS_COMPLETE_TARGET_FREE_A1_MW16_PARAMETER_SEARCH"
        if args.candidate_start == 0 and len(payload["results"]) == len(source_candidates)
        else "PASS_TARGET_FREE_A1_MW16_PARAMETER_SEARCH_SHARD"
    )
    write_payload(args.output, payload)
    print(
        f"A1MW16SEARCH|completed={len(payload['results'])}|"
        f"positive={payload['positive_candidate_count']}|"
        f"failed_closed={payload['failed_closed_candidate_count']}|"
        f"output={relative(args.output)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
