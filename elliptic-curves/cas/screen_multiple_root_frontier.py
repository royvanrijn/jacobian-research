#!/usr/bin/env python3
"""Conductor-screen the short-height frontier from a pinned CRT experiment.

The multiple-root CRT search deliberately retains the union of a score prefix
and a short-height prefix, but its canonical run asks PARI for only four
conductors.  This follow-up consumes that pinned, finite candidate set and
screens a declared height/score union without rerunning or changing the source
experiment.

Root number is recorded as a parity-triage feature only.  Neither it nor a
small conductor is promoted to a Mordell--Weil rank claim.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Sequence

from fermigier_mestre import FermigierMestreFamily
from pari_bridge import minimal_curve_data, pari_version


TARGET_LOG_CONDUCTOR = Decimal("182.72")


def candidate_key(record: dict[str, Any]) -> tuple[int, int, int, int]:
    """Stable order for a union of height and score prefixes."""

    return (
        int(record["height_rank"]),
        int(record["score_rank_within_height_pool"]),
        int(record["height"]),
        int(record["numerator"]),
    )


def select_candidates(
    records: Sequence[dict[str, Any]],
    *,
    height_count: int,
    score_count: int,
) -> list[dict[str, Any]]:
    """Return the exact union of the declared rank prefixes."""

    if height_count < 0 or score_count < 0:
        raise ValueError("selection counts must be nonnegative")
    selected = [
        record
        for record in records
        if int(record["height_rank"]) <= height_count
        or int(record["score_rank_within_height_pool"]) <= score_count
    ]
    selected.sort(key=candidate_key)
    return selected


def verify_source_artifact(source: dict[str, Any], source_path: Path) -> None:
    if source.get("schema_version") != 1:
        raise ValueError("unsupported source-artifact schema")
    if not isinstance(source.get("candidates"), list):
        raise ValueError("source artifact has no candidate list")
    if source.get("target", {}).get("hits"):
        raise ValueError("source artifact unexpectedly claims a target hit")
    for record in source["candidates"]:
        required = {
            "t",
            "numerator",
            "denominator",
            "height",
            "height_rank",
            "score_rank_within_height_pool",
            "h_valuations",
        }
        if not required <= record.keys():
            missing = sorted(required - record.keys())
            raise ValueError(f"candidate is missing required fields: {missing}")
        if int(record["denominator"]) <= 0:
            raise ValueError("candidate denominator must be positive")
    if not source_path.is_file():
        raise ValueError("source artifact does not exist")


def verify_engineered_local_data(
    record: dict[str, Any],
    pari: dict[str, Any],
    constraint_groups: Sequence[dict[str, Any]],
) -> dict[str, bool]:
    """Replay every clean split-multiplicative local certificate."""

    checks: dict[str, bool] = {}
    reductions = pari["local_reduction"]
    for group in constraint_groups:
        prime = str(group["prime"])
        scaling = int(group["presented_model_scaling"])
        expected_delta = int(record["h_valuations"][prime]) - 12 * scaling
        local = reductions[prime]
        valid = (
            int(local["conductor_exponent"]) == 1
            and int(local["minimal_c4_valuation"]) == 0
            and int(local["minimal_discriminant_valuation"]) == expected_delta
            and int(local["ellap"]) == 1
        )
        checks[prime] = valid
        if not valid:
            raise AssertionError(
                f"T={record['t']}: PARI contradicted the p={prime} certificate"
            )
    return checks


def frontier_key(record: dict[str, Any]) -> tuple[Any, ...]:
    """Put low-conductor odd-parity triage candidates first."""

    pari = record["pari"]
    below = Decimal(pari["log_conductor"]) < TARGET_LOG_CONDUCTOR
    odd_parity_heuristic = int(pari["root_number"]) == -1
    return (
        not below,
        not odd_parity_heuristic,
        Decimal(pari["log_conductor"]),
        int(record["height_rank"]),
        record["t"],
    )


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_multiple_root_crt.json"
        ),
    )
    parser.add_argument("--height-count", type=int, default=12)
    parser.add_argument("--score-count", type=int, default=0)
    parser.add_argument("--timeout", type=float, default=45.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_multiple_root_frontier.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.height_count < 0 or args.score_count < 0:
        raise SystemExit("selection counts must be nonnegative")
    if args.height_count == 0 and args.score_count == 0:
        raise SystemExit("at least one selection count must be positive")
    if args.timeout <= 0 or args.stack_bytes < 8_000_000:
        raise SystemExit("PARI timeout and stack bounds must be positive")

    source_bytes = args.input.read_bytes()
    source = json.loads(source_bytes)
    verify_source_artifact(source, args.input)
    selected = select_candidates(
        source["candidates"],
        height_count=args.height_count,
        score_count=args.score_count,
    )
    if not selected:
        raise SystemExit("the declared prefixes selected no retained candidates")

    completed: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    local_primes = tuple(
        int(group["prime"]) for group in source["constraint_groups"]
    )
    for source_record in selected:
        record = {
            key: source_record[key]
            for key in (
                "t",
                "numerator",
                "denominator",
                "height",
                "height_rank",
                "score_rank_within_height_pool",
                "h_valuations",
            )
        }
        parameter = Fraction(record["t"])
        try:
            pari = minimal_curve_data(
                FermigierMestreFamily.coefficients(parameter),
                timeout=args.timeout,
                local_primes=local_primes,
                stack_bytes=args.stack_bytes,
            )
            record["pari"] = pari
            record["engineered_local_checks"] = verify_engineered_local_data(
                record, pari, source["constraint_groups"]
            )
            record["below_strict_log_conductor_target"] = (
                Decimal(pari["log_conductor"]) < TARGET_LOG_CONDUCTOR
            )
            record["odd_rank_parity_heuristic"] = pari["root_number"] == -1
            record["rank_status"] = "not computed or inferred"
            completed.append(record)
            print(
                f"T={record['t']} logN={pari['log_conductor']} "
                f"root={pari['root_number']}",
                flush=True,
            )
        except Exception as error:
            errors.append({"t": record["t"], "error": str(error)})
            print(f"T={record['t']} error={error}", flush=True)

    frontier = sorted(completed, key=frontier_key)
    low_conductor = [
        record for record in completed if record["below_strict_log_conductor_target"]
    ]
    parity_priority = [
        record
        for record in low_conductor
        if record["odd_rank_parity_heuristic"]
    ]
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded conductor screen of a pinned finite candidate union; local "
            "and conductor data are PARI computations, root number is used only "
            "for heuristic triage, and no Mordell--Weil rank is claimed"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
            "reason": "no rank computation or independence certificate was made",
        },
        "selection": {
            "height_rank_at_most": args.height_count,
            "score_rank_at_most": args.score_count,
            "selected": len(selected),
            "completed": len(completed),
            "errors": len(errors),
            "definition": "union of the declared prefixes in the source artifact",
        },
        "summary": {
            "below_strict_log_conductor_target": len(low_conductor),
            "below_target_and_root_number_minus_one": len(parity_priority),
            "best_log_conductor": (
                str(
                    min(
                        Decimal(record["pari"]["log_conductor"])
                        for record in completed
                    )
                )
                if completed
                else None
            ),
            "frontier_order": [record["t"] for record in frontier],
        },
        "method": {
            "frontier_priority": (
                "below conductor target, then root number -1, then increasing "
                "log conductor; this is resource triage, not rank inference"
            ),
            "engineered_local_replay": (
                "all source CRT primes must remain split multiplicative with "
                "conductor exponent one and the predicted minimal discriminant "
                "valuation"
            ),
        },
        "source": {
            "path": str(args.input),
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
            "reproducing_command": source.get("reproducing_command"),
        },
        "parameters": {
            "timeout_seconds_per_candidate": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "errors": errors,
        "candidates": frontier,
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
