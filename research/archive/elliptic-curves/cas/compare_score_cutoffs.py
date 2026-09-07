#!/usr/bin/env python3
"""Compare good-prime scores at growing numerical prime cutoffs.

The comparison is deliberately narrow.  It contains Fermigier's published
rank-22 specialization and the first ten parameters retained by the bounded
``search_record_residue_class.py`` run at score bound 500.  For ``T=a/b`` it
computes

    S_B(T) = sum (2-a_p)/(p+1-a_p) log(p),

where the sum is over numerical primes ``5 <= p <= B`` at which ``p`` does
not divide ``b`` and the specialized curve has good reduction.  PARI/GP's
``ellap`` supplies the exact trace ``a_p``; only the real logarithms and the
resulting sum are approximate.

These scores are rank-triage heuristics.  In particular, extending a score
to more primes does not prove a Mordell--Weil rank or independence of known
points.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
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
from pari_bridge import pari_version


DEFAULT_CUTOFFS = (500, 1_000, 2_000, 5_000, 10_000, 100_000)
PRIMARY_CANDIDATE_ID = "record-top-1"
BENCHMARK_ID = "published-e22"


@dataclass(frozen=True)
class ComparisonCurve:
    identifier: str
    parameter: Fraction
    provenance: str
    record_scan_position_at_500: int | None


# These are pinned inputs, not dynamically read from the earlier artifact.
# That makes this calculation reproducible even if a later search overwrites
# its top-K file with a larger height box or a different score cutoff.
RECORD_CLASS_CURVES = (
    ComparisonCurve("record-top-1", Fraction(1666, 9), "record scan", 1),
    ComparisonCurve("record-top-2", Fraction(1666, 4227), "record scan", 2),
    ComparisonCurve("record-top-3", Fraction(3332, 2571), "record scan", 3),
    ComparisonCurve("record-top-4", Fraction(-2975, 1493), "record scan", 4),
    ComparisonCurve("record-top-5", Fraction(4760, 1289), "record scan", 5),
    ComparisonCurve("record-top-6", Fraction(595, 471), "record scan", 6),
    ComparisonCurve("record-top-7", Fraction(4879, 1213), "record scan", 7),
    ComparisonCurve("record-top-8", Fraction(714, 617), "record scan", 8),
    ComparisonCurve("record-top-9", Fraction(-2023, 2217), "record scan", 9),
    ComparisonCurve("record-top-10", Fraction(1547, 3156), "record scan", 10),
)


def parse_cutoffs(value: str) -> tuple[int, ...]:
    """Parse a strictly increasing comma-separated list of prime cutoffs."""

    try:
        cutoffs = tuple(int(part) for part in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("cutoffs must be integers") from error
    if not cutoffs or any(cutoff < 5 for cutoff in cutoffs):
        raise argparse.ArgumentTypeError("every cutoff must be at least 5")
    if any(left >= right for left, right in zip(cutoffs, cutoffs[1:])):
        raise argparse.ArgumentTypeError("cutoffs must be strictly increasing")
    return cutoffs


def last_primes_for_cutoffs(cutoffs: Sequence[int]) -> tuple[int, ...]:
    """Return the largest prime not exceeding each numerical cutoff."""

    if not cutoffs or any(cutoff < 5 for cutoff in cutoffs):
        raise ValueError("every cutoff must be at least 5")
    if any(left >= right for left, right in zip(cutoffs, cutoffs[1:])):
        raise ValueError("cutoffs must be strictly increasing")
    prime_list = primes_up_to(max(cutoffs))
    last_primes: list[int] = []
    index = 0
    current = 0
    for cutoff in cutoffs:
        while index < len(prime_list) and prime_list[index] <= cutoff:
            current = prime_list[index]
            index += 1
        if current < 5:
            raise AssertionError("no scoring prime exists below a cutoff")
        last_primes.append(current)
    return tuple(last_primes)


def _gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(value)})"


def score_curve_with_pari(
    parameter: Fraction,
    cutoffs: Sequence[int],
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[int, dict[str, Any]]:
    """Return cumulative score records from one exact PARI/GP run."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if stack_bytes < 8_000_000:
        raise ValueError("the PARI stack must be at least 8,000,000 bytes")
    last_primes = last_primes_for_cutoffs(cutoffs)
    coefficients = FermigierMestreFamily.coefficients(parameter)
    vector = ",".join(_gp_rational(value) for value in coefficients)
    cutoff_vector = ",".join(str(value) for value in cutoffs)
    last_prime_vector = ",".join(str(value) for value in last_primes)
    commands = [
        "default(realprecision,80);",
        f"E=ellminimalmodel(ellinit([{vector}]));",
        f"C=[{cutoff_vector}];",
        f"L=[{last_prime_vector}];",
        "S=0;USED=0;DENOM=0;BAD=0;J=1;",
        (
            f"forprime(p=5,{max(cutoffs)},"
            f"if({parameter.denominator}%p==0,DENOM++,"
            "if(valuation(E.disc,p)>0,BAD++,"
            "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++));"
            "while(J<=#C && p==L[J],"
            'print("SCORE|",C[J],"|",S,"|",USED,"|",DENOM,"|",BAD);'
            "J++));"
        ),
        "quit",
    ]
    result = subprocess.run(
        [executable, "-q", "-s", str(stack_bytes)],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP failed: {result.stderr.strip()}")

    records: dict[int, dict[str, Any]] = {}
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line.startswith("SCORE|"):
            continue
        fields = line.split("|")
        if len(fields) != 6:
            raise RuntimeError(f"unexpected PARI score row: {line}")
        cutoff = int(fields[1])
        records[cutoff] = {
            "score": fields[2],
            "primes_used": int(fields[3]),
            "skipped_denominator_primes": int(fields[4]),
            "skipped_bad_primes": int(fields[5]),
        }
    missing = set(cutoffs).difference(records)
    if missing:
        raise RuntimeError(f"PARI omitted score cutoffs: {sorted(missing)}")
    return records


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cutoffs",
        type=parse_cutoffs,
        default=DEFAULT_CUTOFFS,
        help="strictly increasing numerical prime cutoffs",
    )
    parser.add_argument(
        "--record-candidate-count",
        type=int,
        default=len(RECORD_CLASS_CURVES),
        help="number of the pinned top-ten record-scan candidates to compare",
    )
    parser.add_argument("--pari-timeout", type=float, default=30.0)
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_score_cutoffs.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.record_candidate_count <= len(RECORD_CLASS_CURVES):
        raise SystemExit(
            f"--record-candidate-count must be between 1 and "
            f"{len(RECORD_CLASS_CURVES)}"
        )
    if args.pari_timeout <= 0:
        raise SystemExit("--pari-timeout must be positive")
    if args.pari_stack_bytes < 8_000_000:
        raise SystemExit("--pari-stack-bytes must be at least 8000000")

    curves = (
        ComparisonCurve(
            BENCHMARK_ID,
            NORMALIZED_RECORD_PARAMETER,
            "Fermigier's published E22 specialization",
            None,
        ),
        *RECORD_CLASS_CURVES[: args.record_candidate_count],
    )
    score_records: dict[str, dict[int, dict[str, Any]]] = {}
    for curve in curves:
        score_records[curve.identifier] = score_curve_with_pari(
            curve.parameter,
            args.cutoffs,
            timeout=args.pari_timeout,
            stack_bytes=args.pari_stack_bytes,
        )

    comparisons: list[dict[str, Any]] = []
    for cutoff in args.cutoffs:
        benchmark_score = Decimal(
            score_records[BENCHMARK_ID][cutoff]["score"]
        )
        ranked = sorted(
            curves,
            key=lambda curve: (
                -Decimal(score_records[curve.identifier][cutoff]["score"]),
                curve.identifier,
            ),
        )
        ranks = {
            curve.identifier: position for position, curve in enumerate(ranked, 1)
        }
        candidate_score = Decimal(
            score_records[PRIMARY_CANDIDATE_ID][cutoff]["score"]
        )
        comparisons.append(
            {
                "numeric_prime_cutoff": cutoff,
                "candidate_score": str(candidate_score),
                "benchmark_score": str(benchmark_score),
                "candidate_minus_benchmark": str(
                    candidate_score - benchmark_score
                ),
                "candidate_rank_among_compared": ranks[PRIMARY_CANDIDATE_ID],
                "curves_compared": len(curves),
                "ranking": [curve.identifier for curve in ranked],
            }
        )

    curve_output = []
    for curve in curves:
        curve_output.append(
            {
                "id": curve.identifier,
                "parameter": rational_to_string(curve.parameter),
                "provenance": curve.provenance,
                "record_scan_position_at_numeric_cutoff_500": (
                    curve.record_scan_position_at_500
                ),
                "scores": {
                    str(cutoff): score_records[curve.identifier][cutoff]
                    for cutoff in args.cutoffs
                },
            }
        )

    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded computational comparison of rank-triage heuristics; "
            "no score is a rank or point-independence inference"
        ),
        "selection_warning": (
            "record-top-1 through record-top-10 were selected by their scores "
            "at numerical cutoff 500 in an earlier bounded height-5000 scan, "
            "so comparisons involving that cutoff are selection-biased"
        ),
        "family": "normalized Fermigier--Mestre family",
        "score_definition": {
            "formula": (
                "sum_{5<=p<=B, p not dividing b, good reduction} "
                "((2-a_p)/(p+1-a_p))*log(p), for T=a/b in lowest terms"
            ),
            "prime_bound_semantics": "numerical prime cutoff p<=B",
            "trace_source": "exact PARI/GP ellap on an exact minimal model",
            "exclusions": (
                "p=2,3; primes dividing the parameter denominator; and all "
                "primes dividing the minimal discriminant"
            ),
            "approximation": "PARI real logarithms and their cumulative sums only",
        },
        "parameters": {
            "cutoffs": list(args.cutoffs),
            "record_candidate_count": args.record_candidate_count,
            "pari_timeout": args.pari_timeout,
            "pari_stack_bytes": args.pari_stack_bytes,
            "output": str(args.output),
        },
        "primary_candidate": {
            "id": PRIMARY_CANDIDATE_ID,
            "parameter": "1666/9",
            "comparison_to_benchmark": comparisons,
        },
        "curves": curve_output,
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    for comparison in comparisons:
        print(
            f"B={comparison['numeric_prime_cutoff']} "
            f"candidate_rank={comparison['candidate_rank_among_compared']}/"
            f"{comparison['curves_compared']} "
            f"candidate_minus_benchmark="
            f"{comparison['candidate_minus_benchmark']}"
        )


if __name__ == "__main__":
    main()
