#!/usr/bin/env sage-python
"""Run direct-reduction MW16 half-lattice recovery on exact Nagao finalists.

The input is the exact 104-fibre raw-specialization ledger produced after the
height-300 local sieve.  Every fibre first rechecks specialization independence
and its complete generic M/2M depth spectrum, then searches only the exact
maximum-depth stratum.  Pointed quartics use exact square-content normalization
and direct Cremona--Stoll reduction, with no quartic-minimalization call.
Results are checkpointed after every candidate.

No adaptive quotient wave, unrestricted point search, global curve
minimalization, or Selmer calculation is performed.  Even a positive must
first be transported to an exact global minimal model with renewed section
checks before reaching the same-fibre residual-Selmer gate; zero cannot reject
the fibre mathematically.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import platform
import shutil


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json"
LADDER = ROOT / "elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
DIRECT = ROOT / "elliptic-curves/cas/half_lattice_direct_reduction.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_direct_reduction_h300_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def write_payload(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def initial_payload(args, inputs, candidate_ids):
    return {
        "schema": "elliptic-curves.icarm-mw16-nagao-finalist-half-lattice.v1",
        "status": "SEARCHING",
        "source_finalist_count": len(inputs["candidates"]),
        "candidate_ids": candidate_ids,
        "declared_budget": {
            "height_bound_each_quartic": args.height_bound,
            "timeout_seconds_each_quartic": args.timeout_seconds,
            "stack_bytes_each_quartic": args.stack_bytes,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "chart_rule": (
                "complete exact maximum-depth stratum of the specialized generic "
                "MW16 half-lattice, recomputed separately for each fibre"
            ),
            "adaptive_quotient_lifts": 0,
            "unrestricted_point_search": False,
            "quartic_backend": "exact_square_content_then_direct_hyperellred_v1",
            "quartic_minimalization_called": False,
        },
        "results": [],
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, LADDER, LEGACY, ENGINE, DIRECT, Path(__file__))
        },
        "software": {"python": platform.python_version()},
        "next_gate": {
            "positive_definition": "exact_quotient_rank_recovered >= 1",
            "action": (
                "compute the global minimal model, transport the specialized MW16 "
                "basis and every recovered independent point, then request the "
                "complete residual 2-Selmer group on that same curve"
            ),
            "expensive_continuation_authorized_before_complete_selmer": False,
        },
        "claim_boundary": [
            "Nagao supplied only the frozen candidate order; it contributes no rank evidence.",
            "Every positive quotient rank is supported by exact rational points, exact group law, and finite-reduction independence certificates.",
            "Every miss is bounded and gives no rank upper bound, point absence, saturation, covering, or Selmer information.",
            "Direct hyperelliptic reduction is exact but makes no global-minimal-model claim.",
            "The stage stops before adaptive quotient lifts or any unrestricted point search.",
            "No positive candidate is authorized for expensive continuation until its complete same-minimal-curve residual 2-Selmer gate finishes.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/run_icarm_mw16_nagao_finalist_half_lattice.sage"
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
    parser.add_argument(
        "--maximum-candidates",
        type=int,
        default=0,
        help="0 means the entire frozen finalist ledger; positive values are a prefix",
    )
    parser.add_argument(
        "--candidate-start",
        type=int,
        default=0,
        help="zero-based start in the frozen finalist ledger (for checkpoint shards)",
    )
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
    if inputs.get("status") != "PASS_EXACT_MW16_NAGAO_FINALIST_SPECIALIZATIONS":
        raise ArithmeticError("exact finalist specialization ledger is not passing")
    if (
        inputs["next_gate"]["stage"]
        != "bounded_half_lattice_arithmetic_size_diagnostic"
    ):
        raise ArithmeticError("finalist ledger does not authorize the raw-model diagnostic")
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
        for key in (
            "schema",
            "source_finalist_count",
            "candidate_ids",
            "declared_budget",
            "inputs",
        ):
            if payload.get(key) != fresh.get(key):
                raise ArithmeticError(
                    f"existing checkpoint differs in immutable field {key}"
                )
    else:
        payload = fresh
        write_payload(args.output, payload)
    completed = {row["candidate_id"] for row in payload["results"]}
    if not completed.issubset(set(candidate_ids)):
        raise ArithmeticError("checkpoint contains a candidate outside the frozen prefix")

    ladder = SourceFileLoader("mw16_finalist_ladder", str(LADDER)).load_module()
    legacy = SourceFileLoader("mw16_finalist_legacy", str(LEGACY)).load_module()
    direct = SourceFileLoader("mw16_finalist_direct", str(DIRECT)).load_module()
    legacy.GENERIC_DIMENSION = 16
    legacy.engine = direct
    for index, candidate in enumerate(candidates, 1):
        if candidate["candidate_id"] in completed:
            print(
                f"MW16FINALISTHL|candidate={candidate['candidate_id']}|"
                f"index={index}/{len(candidates)}|status=RESUME_SKIP",
                flush=True,
            )
            continue
        adapted = {
            "parent_id": candidate["candidate_id"],
            "curve_id": int(candidate["parent_curve_id"]),
            "priority_rank": int(
                candidate["parent_id"].rsplit("p", 1)[1]
            ),
            "target_short_model": candidate["raw_short_model"],
            "specialized_generic_points": candidate["raw_generic_points"],
            "generic_height_gram": candidate["generic_height_gram"],
        }
        try:
            result = ladder.run_parent(adapted, legacy, args)
            result["candidate_id"] = result.pop("parent_id")
            result["parent_id"] = candidate["parent_id"]
            result["parameter"] = candidate["parameter"]
            result["q_isomorphism_class_id"] = candidate[
                "q_isomorphism_class_id"
            ]
            result["nagao"] = candidate["nagao"]
        except ArithmeticError as error:
            result = {
                "candidate_id": candidate["candidate_id"],
                "parent_id": candidate["parent_id"],
                "parameter": candidate["parameter"],
                "q_isomorphism_class_id": candidate[
                    "q_isomorphism_class_id"
                ],
                "nagao": candidate["nagao"],
                "status": "REJECTED_FAIL_CLOSED_BEFORE_EXACT_QUOTIENT_RANK",
                "reason": str(error),
                "exact_quotient_rank_recovered": None,
            }
        payload["results"].append(result)
        write_payload(args.output, payload)
        print(
            f"MW16FINALISTHL|candidate={candidate['candidate_id']}|"
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
        row["exact_quotient_rank_recovered"] is None
        for row in payload["results"]
    )
    payload["status"] = (
        "PASS_COMPLETE_FROZEN_NAGAO_FINALIST_HALF_LATTICE_GATE"
        if args.candidate_start == 0
        and len(payload["results"]) == len(source_candidates)
        else "PASS_BOUNDED_PREFIX_NAGAO_FINALIST_HALF_LATTICE_GATE"
    )
    write_payload(args.output, payload)
    print(
        f"MW16FINALISTHL|completed={len(payload['results'])}|"
        f"positive={payload['positive_candidate_count']}|"
        f"failed_closed={payload['failed_closed_candidate_count']}|"
        f"output={relative(args.output)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
