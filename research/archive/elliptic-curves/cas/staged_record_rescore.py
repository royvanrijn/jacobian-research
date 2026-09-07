#!/usr/bin/env python3
"""Leakage-free staged rescore of Fermigier's record residue class.

Every primitive ``T=a/b`` in the declared projective-height box and the
fixed class ``a=2142*b (mod 4403)`` is scored at the first cutoff.  Only then
are the best survivors rescored at successively larger prime cutoffs.  This
avoids selecting at 500 and retrospectively presenting a larger-cutoff score
as though the full population had been searched there.

At every stage the good-prime heuristic is

    sum ((2-a_p)/(p+1-a_p))*log(p)

over numerical primes ``5 <= p <= B``, excluding parameter-denominator and
bad-reduction primes.  PARI/GP supplies exact ``a_p`` through ``ellap``; the
logarithms and sums are approximate.  Scores do not imply rank.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import json
from math import gcd
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Sequence

from compare_score_cutoffs import last_primes_for_cutoffs, parse_cutoffs
from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily, NORMALIZED_RECORD_PARAMETER
from pari_bridge import minimal_curve_data, pari_version
from search_record_residue_class import (
    CRT_MODULUS,
    CRT_RESIDUE,
    primitive_rationals_in_class,
)


TARGET_LOG_CONDUCTOR = Decimal("182.72")
BENCHMARK_ID = "published-e22"


@dataclass(frozen=True)
class Candidate:
    numerator: int
    denominator: int

    @property
    def parameter(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)

    @property
    def height(self) -> int:
        return max(abs(self.numerator), self.denominator)

    @property
    def identifier(self) -> str:
        return rational_to_string(self.parameter)


def parse_keep_counts(value: str) -> tuple[int, ...]:
    try:
        counts = tuple(int(part) for part in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("keep counts must be integers") from error
    if not counts or any(count < 1 for count in counts):
        raise argparse.ArgumentTypeError("every keep count must be positive")
    if any(left < right for left, right in zip(counts, counts[1:])):
        raise argparse.ArgumentTypeError("keep counts must be nonincreasing")
    return counts


def candidate_sort_key(record: dict[str, Any]) -> tuple[Any, ...]:
    """Order by descending score, then exact deterministic height tie-breaks."""

    return (
        -Decimal(record["score"]),
        record["height"],
        record["numerator"],
        record["denominator"],
    )


def _gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(value)})"


def score_candidates_with_pari(
    candidates: Sequence[Candidate],
    cutoff: int,
    *,
    timeout: float,
    stack_bytes: int,
) -> list[dict[str, Any]]:
    """Score a finite candidate list in one exact-model PARI/GP batch."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if cutoff < 5:
        raise ValueError("the numerical prime cutoff must be at least 5")
    if stack_bytes < 8_000_000:
        raise ValueError("the PARI stack must be at least 8,000,000 bytes")
    last_prime = last_primes_for_cutoffs((cutoff,))[0]
    commands = ["default(realprecision,80);"]
    for index, candidate in enumerate(candidates):
        coefficients = FermigierMestreFamily.coefficients(candidate.parameter)
        vector = ",".join(_gp_rational(value) for value in coefficients)
        commands.extend(
            [
                f"E=ellminimalmodel(ellinit([{vector}]));",
                "S=0;USED=0;DENOM=0;BAD=0;",
                (
                    f"forprime(p=5,{last_prime},"
                    f"if({candidate.denominator}%p==0,DENOM++,"
                    "if(valuation(E.disc,p)>0,BAD++,"
                    "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++)));"
                ),
                (
                    f'print("SCORE|{index}|",S,"|",USED,"|",'
                    'DENOM,"|",BAD);'
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

    score_rows: dict[int, tuple[str, int, int, int]] = {}
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line.startswith("SCORE|"):
            continue
        fields = line.split("|")
        if len(fields) != 6:
            raise RuntimeError(f"unexpected PARI score row: {line}")
        score_rows[int(fields[1])] = (
            fields[2],
            int(fields[3]),
            int(fields[4]),
            int(fields[5]),
        )
    if set(score_rows) != set(range(len(candidates))):
        raise RuntimeError("PARI omitted or duplicated candidate score rows")

    records: list[dict[str, Any]] = []
    for index, candidate in enumerate(candidates):
        score, used, skipped_denominator, skipped_bad = score_rows[index]
        records.append(
            {
                "t": candidate.identifier,
                "numerator": candidate.numerator,
                "denominator": candidate.denominator,
                "height": candidate.height,
                "score": score,
                "primes_used": used,
                "skipped_denominator_primes": skipped_denominator,
                "skipped_bad_primes": skipped_bad,
            }
        )
    records.sort(key=candidate_sort_key)
    return records


def assert_population(candidates: Sequence[Candidate], height: int) -> None:
    """Check exact normalization and uniqueness of an enumerated population."""

    identifiers = set()
    for candidate in candidates:
        if candidate.denominator <= 0 or candidate.height > height:
            raise AssertionError("a candidate is outside the projective-height box")
        if gcd(candidate.numerator, candidate.denominator) != 1:
            raise AssertionError("a candidate is not primitive")
        if (
            candidate.numerator - CRT_RESIDUE * candidate.denominator
        ) % CRT_MODULUS:
            raise AssertionError("a candidate is outside the record residue class")
        if candidate.identifier in identifiers:
            raise AssertionError("the candidate population contains a duplicate")
        identifiers.add(candidate.identifier)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=5000)
    parser.add_argument(
        "--stage-cutoffs",
        type=parse_cutoffs,
        default=(2000, 10_000, 100_000),
    )
    parser.add_argument(
        "--keep-counts",
        type=parse_keep_counts,
        default=(50, 10, 10),
        help="one retained count for every score stage",
    )
    parser.add_argument("--conductor-count", type=int, default=0)
    parser.add_argument(
        "--score-timeout",
        type=float,
        default=300.0,
        help="timeout for each whole score-stage PARI batch",
    )
    parser.add_argument(
        "--conductor-timeout",
        type=float,
        default=60.0,
        help="timeout for each finalist conductor/point check",
    )
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_record_rescore_h5000.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.height < 1:
        raise SystemExit("--height must be positive")
    if len(args.keep_counts) != len(args.stage_cutoffs):
        raise SystemExit("--keep-counts must have one entry per --stage-cutoffs entry")
    if args.conductor_count < 0:
        raise SystemExit("--conductor-count must be nonnegative")
    if args.score_timeout <= 0 or args.conductor_timeout <= 0:
        raise SystemExit("PARI timeouts must be positive")

    population = tuple(
        Candidate(numerator, denominator)
        for numerator, denominator in primitive_rationals_in_class(args.height)
    )
    assert_population(population, args.height)
    benchmark = Candidate(
        NORMALIZED_RECORD_PARAMETER.numerator,
        NORMALIZED_RECORD_PARAMETER.denominator,
    )
    benchmark_in_population = benchmark in population

    stage_input = population
    stages: list[dict[str, Any]] = []
    survivor_histories: dict[str, dict[str, Any]] = {}
    benchmark_scores: dict[str, dict[str, Any]] = {}
    for cutoff, keep_count in zip(args.stage_cutoffs, args.keep_counts):
        start = time.monotonic()
        records = score_candidates_with_pari(
            stage_input,
            cutoff,
            timeout=args.score_timeout,
            stack_bytes=args.pari_stack_bytes,
        )
        elapsed = time.monotonic() - start
        retained = records[: min(keep_count, len(records))]
        for position, record in enumerate(records, 1):
            if record["t"] in {item["t"] for item in retained}:
                history = survivor_histories.setdefault(
                    record["t"],
                    {
                        key: record[key]
                        for key in ("t", "numerator", "denominator", "height")
                    },
                )
                history.setdefault("stages", {})[str(cutoff)] = {
                    "score": record["score"],
                    "position_among_stage_input": position,
                    "primes_used": record["primes_used"],
                    "skipped_denominator_primes": record[
                        "skipped_denominator_primes"
                    ],
                    "skipped_bad_primes": record["skipped_bad_primes"],
                }

        benchmark_record = next(
            (record for record in records if record["t"] == benchmark.identifier),
            None,
        )
        benchmark_source = "stage population"
        if benchmark_record is None:
            benchmark_record = score_candidates_with_pari(
                (benchmark,),
                cutoff,
                timeout=args.score_timeout,
                stack_bytes=args.pari_stack_bytes,
            )[0]
            benchmark_source = "separate comparison; not in stage population"
        benchmark_scores[str(cutoff)] = {
            "score": benchmark_record["score"],
            "source": benchmark_source,
            "position_among_stage_input": (
                next(
                    (
                        index
                        for index, record in enumerate(records, 1)
                        if record["t"] == benchmark.identifier
                    ),
                    None,
                )
            ),
        }

        retained_ids = {record["t"] for record in retained}
        stage_input = tuple(
            Candidate(record["numerator"], record["denominator"])
            for record in retained
        )
        stages.append(
            {
                "numeric_prime_cutoff": cutoff,
                "input_count": len(records),
                "retained_count": len(retained),
                "requested_keep_count": keep_count,
                "elapsed_seconds": round(elapsed, 6),
                "best": retained[0] if retained else None,
                "worst_retained_score": retained[-1]["score"] if retained else None,
                "benchmark": benchmark_scores[str(cutoff)],
                "benchmark_survived": benchmark.identifier in retained_ids,
            }
        )

    final_records = [
        survivor_histories[candidate.identifier] for candidate in stage_input
    ]
    final_cutoff = args.stage_cutoffs[-1]
    final_records.sort(
        key=lambda record: (
            Decimal(record["stages"][str(final_cutoff)]["score"]) * -1,
            record["height"],
            record["numerator"],
            record["denominator"],
        )
    )
    conductor_errors: list[dict[str, str]] = []
    for record in final_records[: args.conductor_count]:
        parameter = Fraction(record["numerator"], record["denominator"])
        try:
            known_points = FermigierMestreFamily.known_jacobian_points(parameter)[1:]
            record["pari"] = minimal_curve_data(
                FermigierMestreFamily.coefficients(parameter),
                timeout=args.conductor_timeout,
                known_points=known_points,
                stack_bytes=args.pari_stack_bytes,
            )
            record["below_log_conductor_target"] = (
                Decimal(record["pari"]["log_conductor"])
                < TARGET_LOG_CONDUCTOR
            )
        except Exception as error:
            conductor_errors.append({"t": record["t"], "error": str(error)})

    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exhaustive first-stage and staged computational rank-triage "
            "experiment; scores and numerical height determinants imply no rank"
        ),
        "family": "normalized Fermigier--Mestre family",
        "residue_class": {
            "condition": "a=2142*b mod 4403 for primitive T=a/b with b>0",
            "crt_residue": CRT_RESIDUE,
            "crt_modulus": CRT_MODULUS,
        },
        "score_definition": {
            "formula": (
                "sum_{5<=p<=B, p not dividing b, good reduction} "
                "((2-a_p)/(p+1-a_p))*log(p)"
            ),
            "prime_bound_semantics": "numerical prime cutoff p<=B",
            "trace_source": "exact PARI/GP ellap on exact minimal models",
            "approximation": "real logarithms and their cumulative sums only",
            "rank_inference": "none",
        },
        "selection_protocol": (
            "The full declared height-box population is scored at the first "
            "cutoff. Each later stage sees only the explicitly retained leaders "
            "from its immediate predecessor."
        ),
        "parameters": {
            "height": args.height,
            "stage_cutoffs": list(args.stage_cutoffs),
            "keep_counts": list(args.keep_counts),
            "conductor_count": args.conductor_count,
            "score_timeout": args.score_timeout,
            "conductor_timeout": args.conductor_timeout,
            "pari_stack_bytes": args.pari_stack_bytes,
            "output": str(args.output),
        },
        "population": {
            "primitive_candidates": len(population),
            "benchmark_in_population": benchmark_in_population,
        },
        "benchmark": {
            "id": BENCHMARK_ID,
            "parameter": benchmark.identifier,
            "published_rank_lower_bound": 22,
            "rank_source": "published metadata; not recomputed here",
            "scores": benchmark_scores,
        },
        "stages": stages,
        "finalists": final_records,
        "conductor_errors": conductor_errors,
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(f"population={len(population)} finalists={len(final_records)}")
    for stage in stages:
        print(
            f"B={stage['numeric_prime_cutoff']} "
            f"input={stage['input_count']} retained={stage['retained_count']} "
            f"best={stage['best']['t']} score={stage['best']['score']}"
        )
    for record in final_records[: args.conductor_count]:
        print(
            f"finalist={record['t']} "
            f"logN={record.get('pari', {}).get('log_conductor')} "
            f"below_target={record.get('below_log_conductor_target')}"
        )


if __name__ == "__main__":
    main()
