#!/usr/bin/env python3
"""Leakage-free integer-parameter scan in Nagao's rank-13 base change.

Nagao's six-root tuple ``(148,116,104,57,25,0)`` has an additional displayed
point, and after

``T = (-u^2 + 23550)/(2u)``

the cited family has thirteen independent sections over ``Q(u)``.  This script
scores a finite set of positive integer ``u`` in three declared stages, then
computes conductors only for the final survivors.

Generic independence does not by itself certify any specialization.  Scores,
root numbers, and low conductors are therefore triage data, never target hits.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
from typing import Any, Sequence

from ek_k3 import primes_up_to, rational_to_string
from fermigier_mestre import FermigierMestreFamily, NORMALIZED_RECORD_PARAMETER
from mestre_root_tuples import SixRootMestreConstruction
from pari_bridge import minimal_curve_data, pari_version


NAGAO_ROOTS = (148, 116, 104, 57, 25, 0)
NAGAO_CONSTRUCTION = SixRootMestreConstruction(
    tuple(Fraction(root) for root in NAGAO_ROOTS)
)
TARGET_LOG_CONDUCTOR = Decimal("182.72")


@dataclass(frozen=True)
class Candidate:
    identifier: str
    parameter_u: int | None
    parameter_t: Fraction | None
    coefficients: tuple[Fraction, ...]
    provenance: str


def nagao_base_change(parameter_u: int) -> Fraction:
    if not isinstance(parameter_u, int) or parameter_u == 0:
        raise ValueError("Nagao's integer parameter u must be nonzero")
    return Fraction(-parameter_u**2 + 23550, 2 * parameter_u)


def nagao_curve(parameter_u: int) -> Candidate:
    parameter_t = nagao_base_change(parameter_u)
    invariant_i, invariant_j = NAGAO_CONSTRUCTION.binary_invariants(parameter_t)
    return Candidate(
        identifier=f"nagao-u-{parameter_u}",
        parameter_u=parameter_u,
        parameter_t=parameter_t,
        coefficients=(
            Fraction(0),
            Fraction(0),
            Fraction(0),
            -27 * invariant_i,
            -27 * invariant_j,
        ),
        provenance="Nagao rank-13 base-changed family",
    )


def benchmark_curves() -> tuple[Candidate, ...]:
    return (
        Candidate(
            "published-e22",
            None,
            NORMALIZED_RECORD_PARAMETER,
            FermigierMestreFamily.coefficients(NORMALIZED_RECORD_PARAMETER),
            "published Fermigier rank-at-least-22 benchmark",
        ),
        Candidate(
            "fermigier-1666-9",
            None,
            Fraction(1666, 9),
            FermigierMestreFamily.coefficients(Fraction(1666, 9)),
            "best Fermigier conductor-only candidate",
        ),
    )


def parse_stages(value: str) -> tuple[int, ...]:
    try:
        stages = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("stages must be integers") from error
    if not stages or any(stage < 5 for stage in stages):
        raise argparse.ArgumentTypeError("every score stage must be at least 5")
    if any(left >= right for left, right in zip(stages, stages[1:])):
        raise argparse.ArgumentTypeError("score stages must be strictly increasing")
    return stages


def parse_keep_counts(value: str) -> tuple[int, ...]:
    try:
        counts = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("keep counts must be integers") from error
    if not counts or any(count <= 0 for count in counts):
        raise argparse.ArgumentTypeError("keep counts must be positive")
    if any(left < right for left, right in zip(counts, counts[1:])):
        raise argparse.ArgumentTypeError("keep counts must be nonincreasing")
    return counts


def _gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def score_candidates_with_pari(
    candidates: Sequence[Candidate],
    *,
    cutoff: int,
    timeout: float,
    stack_bytes: int,
) -> dict[str, dict[str, Any]]:
    """Score every input curve in one deterministic PARI process."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if cutoff < 5 or timeout <= 0 or stack_bytes < 8_000_000:
        raise ValueError("invalid score-process bounds")
    last_prime = primes_up_to(cutoff)[-1]
    commands = ["default(realprecision,80);"]
    for index, candidate in enumerate(candidates):
        vector = ",".join(_gp_rational(value) for value in candidate.coefficients)
        commands.extend(
            [
                f"E=ellminimalmodel(ellinit([{vector}]));",
                "S=0;USED=0;BAD=0;",
                (
                    f"forprime(p=5,{cutoff},"
                    "if(valuation(E.disc,p)>0,BAD++,"
                    "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++));"
                ),
                (
                    f'print("ROW|{index}|",S,"|",USED,"|",BAD,'
                    f'"|{last_prime}");'
                ),
            ]
        )
    commands.append("quit")
    result = subprocess.run(
        [executable, "-q", "-s", str(stack_bytes)],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP failed: {result.stderr.strip()}")
    records: dict[str, dict[str, Any]] = {}
    for line in result.stdout.splitlines():
        if not line.startswith("ROW|"):
            continue
        _, index_text, score, used, bad, observed_last = line.split("|")
        index = int(index_text)
        candidate = candidates[index]
        records[candidate.identifier] = {
            "score": score,
            "good_primes_used": int(used),
            "bad_primes_skipped": int(bad),
            "last_numerical_prime": int(observed_last),
        }
    if len(records) != len(candidates):
        raise RuntimeError("PARI omitted one or more score records")
    return records


def score_sort_key(record: dict[str, Any]) -> tuple[Decimal, int, str]:
    return (
        -Decimal(record["score"]),
        int(record.get("parameter_u") or 0),
        record["candidate_id"],
    )


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--u-bound", type=int, default=200)
    parser.add_argument("--stages", type=parse_stages, default=(200, 1_000, 10_000))
    parser.add_argument(
        "--keep-counts", type=parse_keep_counts, default=(40, 12, 12)
    )
    parser.add_argument("--score-timeout", type=float, default=120.0)
    parser.add_argument("--conductor-timeout", type=float, default=45.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_integer_u.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.u_bound < 1:
        raise SystemExit("--u-bound must be positive")
    if len(args.stages) != len(args.keep_counts):
        raise SystemExit("--stages and --keep-counts must have equal length")
    if args.keep_counts[0] > args.u_bound:
        raise SystemExit("the first keep count exceeds the population")
    if args.score_timeout <= 0 or args.conductor_timeout <= 0:
        raise SystemExit("timeouts must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("PARI stack bound is too small")

    candidates = [nagao_curve(parameter_u) for parameter_u in range(1, args.u_bound + 1)]
    by_id = {candidate.identifier: candidate for candidate in candidates}
    survivors = candidates
    stage_records: list[dict[str, Any]] = []
    latest_scores: dict[str, dict[str, Any]] = {}
    for cutoff, keep_count in zip(args.stages, args.keep_counts):
        scores = score_candidates_with_pari(
            survivors,
            cutoff=cutoff,
            timeout=args.score_timeout,
            stack_bytes=args.stack_bytes,
        )
        ranked = []
        for candidate in survivors:
            record = {
                "candidate_id": candidate.identifier,
                "parameter_u": candidate.parameter_u,
                "parameter_t": rational_to_string(candidate.parameter_t),
                **scores[candidate.identifier],
            }
            ranked.append(record)
            latest_scores[candidate.identifier] = record
        ranked.sort(key=score_sort_key)
        retained = ranked[: min(keep_count, len(ranked))]
        stage_records.append(
            {
                "cutoff": cutoff,
                "population_scored": len(ranked),
                "keep_count": len(retained),
                "ranked_population": ranked,
                "retained_candidate_ids": [record["candidate_id"] for record in retained],
            }
        )
        survivors = [by_id[record["candidate_id"]] for record in retained]

    conductor_records: list[dict[str, Any]] = []
    conductor_errors: list[dict[str, str]] = []
    for candidate in survivors:
        record = dict(latest_scores[candidate.identifier])
        try:
            pari = minimal_curve_data(
                candidate.coefficients,
                timeout=args.conductor_timeout,
                stack_bytes=args.stack_bytes,
            )
            record["pari"] = pari
            record["below_strict_log_conductor_target"] = (
                Decimal(pari["log_conductor"]) < TARGET_LOG_CONDUCTOR
            )
            record["rank_status"] = (
                "generic rank-at-least-13 is cited; specialization independence "
                "was not checked in this artifact"
            )
            conductor_records.append(record)
        except Exception as error:
            conductor_errors.append(
                {"candidate_id": candidate.identifier, "error": str(error)}
            )

    benchmark_scores: list[dict[str, Any]] = []
    for benchmark in benchmark_curves():
        per_cutoff = {}
        for cutoff in args.stages:
            per_cutoff[str(cutoff)] = score_candidates_with_pari(
                (benchmark,),
                cutoff=cutoff,
                timeout=args.score_timeout,
                stack_bytes=args.stack_bytes,
            )[benchmark.identifier]
        benchmark_scores.append(
            {
                "candidate_id": benchmark.identifier,
                "parameter": rational_to_string(benchmark.parameter_t),
                "provenance": benchmark.provenance,
                "scores": per_cutoff,
            }
        )

    conductor_records.sort(
        key=lambda record: (
            not record["below_strict_log_conductor_target"],
            -Decimal(record["score"]),
            Decimal(record["pari"]["log_conductor"]),
        )
    )
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded leakage-free score/conductor experiment; Nagao's generic "
            "rank-at-least-13 theorem is cited, specialization independence and "
            "extra rank are not certified, and no target hit is claimed"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "family": {
            "root_tuple": list(NAGAO_ROOTS),
            "base_change": "T=(-u^2+23550)/(2u)",
            "published_generic_rank_lower_bound": 13,
            "status": "cited theorem; not reproved by this search",
            "source": "https://doi.org/10.24546/E0003610",
        },
        "selection_protocol": {
            "population": f"every positive integer 1<=u<={args.u_bound}",
            "stages": list(args.stages),
            "keep_counts": list(args.keep_counts),
            "leakage_control": (
                "each cutoff scores only survivors selected at the preceding "
                "cutoff; conductor work occurs only after the final cutoff"
            ),
            "score": "sum ((2-a_p)/(p+1-a_p))*log(p) over good numerical primes",
        },
        "stages": stage_records,
        "final_conductor_candidates": conductor_records,
        "conductor_errors": conductor_errors,
        "benchmarks": benchmark_scores,
        "summary": {
            "population": len(candidates),
            "finalists": len(survivors),
            "conductor_calls_completed": len(conductor_records),
            "below_strict_log_conductor_target": sum(
                record["below_strict_log_conductor_target"]
                for record in conductor_records
            ),
            "best_below_target_candidate": (
                conductor_records[0]["candidate_id"]
                if conductor_records
                and conductor_records[0]["below_strict_log_conductor_target"]
                else None
            ),
        },
        "parameters": {
            "u_bound": args.u_bound,
            "score_timeout_seconds_per_stage": args.score_timeout,
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    for record in conductor_records:
        print(
            f"u={record['parameter_u']} T={record['parameter_t']} "
            f"score={record['score']} logN={record['pari']['log_conductor']} "
            f"root={record['pari']['root_number']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
