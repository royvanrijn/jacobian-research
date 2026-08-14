#!/usr/bin/env python3
"""Bounded CRT/lattice search in Nagao's generic-rank-13 base change.

The search discovers certified split-multiplicative root balls at
``p=7,11,13,19,31``, combines every declared local-symbol choice by CRT, and
uses exact two-dimensional Gauss reduction to recover small rational ``u``.
Candidates are scored in staged PARI batches and conductors are computed only
for the final survivors.

Nagao's generic-rank result is context, not a specialization certificate.
This script makes no rank claim and never reports a target hit from score,
root number, or conductor alone.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import log
import os
from pathlib import Path
import platform
import shlex
import shutil
import signal
import subprocess
import sys
from typing import Any, Sequence

from crt_lattice import crt_pair, short_rational_representatives
from ek_k3 import primes_up_to, rational_to_string
from fermigier_mestre import FermigierMestreFamily, NORMALIZED_RECORD_PARAMETER
from nagao_1994 import (
    PRIMARY_SOURCE,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
)
from nagao_rank13_local import (
    DEFAULT_CRT_PRIMES,
    NagaoLocalBall,
    default_crt_balls,
    default_local_discoveries,
    rational_discriminant_valuation,
)
from pari_bridge import minimal_curve_data, pari_version


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank13_local_crt.py"
)


@dataclass(frozen=True)
class CRTNagaoCandidate:
    identifier: str
    numerator: int
    denominator: int
    height: int
    crt_residue: int
    crt_modulus: int
    choices: tuple[NagaoLocalBall, ...]

    @property
    def parameter_u(self) -> Fraction:
        return Q(self.numerator, self.denominator)

    @property
    def parameter_t(self) -> Fraction:
        return rank13_base_parameter(self.parameter_u)

    @property
    def coefficients(self) -> tuple[Fraction, ...]:
        return rank13_base_changed_short_jacobian_coefficients(self.parameter_u)


@dataclass(frozen=True)
class ScoreCandidate:
    identifier: str
    coefficients: tuple[Fraction, ...]


def enumerate_crt_candidates(
    *,
    coefficient_radius: int,
    representatives_per_class: int,
) -> tuple[CRTNagaoCandidate, ...]:
    """Exhaust every default local-symbol combination and reconstruct rational u."""

    if coefficient_radius < 1 or representatives_per_class < 1:
        raise ValueError("lattice enumeration bounds must be positive")
    groups = default_crt_balls()
    ordered_groups = tuple(groups[prime] for prime in DEFAULT_CRT_PRIMES)
    candidates: dict[tuple[int, int], CRTNagaoCandidate] = {}
    serial = 0
    for choices in product(*ordered_groups):
        residue, modulus = 0, 1
        for choice in choices:
            residue, modulus = crt_pair(
                residue, modulus, choice.residue, choice.modulus
            )
        representatives = short_rational_representatives(
            residue,
            modulus,
            coefficient_radius=coefficient_radius,
            limit=representatives_per_class,
        )
        for representative in representatives:
            if representative.numerator == 0:
                continue
            key = (representative.numerator, representative.denominator)
            if key in candidates:
                continue
            for choice in choices:
                if (
                    representative.numerator
                    - choice.residue * representative.denominator
                ) % choice.modulus:
                    raise AssertionError("a reconstructed parameter missed its ball")
                observed = rational_discriminant_valuation(
                    representative.numerator,
                    representative.denominator,
                    choice.prime,
                )
                if observed < choice.requested_discriminant_exponent:
                    raise AssertionError("a reconstructed parameter lost a p-adic power")
            serial += 1
            candidates[key] = CRTNagaoCandidate(
                identifier=f"nagao-crt-{serial:05d}",
                numerator=representative.numerator,
                denominator=representative.denominator,
                height=representative.height,
                crt_residue=residue,
                crt_modulus=modulus,
                choices=tuple(choices),
            )
    ordered = sorted(
        candidates.values(),
        key=lambda candidate: (
            candidate.height,
            abs(candidate.numerator),
            candidate.denominator,
            candidate.numerator < 0,
            candidate.crt_residue,
        ),
    )
    # The quartic coefficients are even in T, while u -> -u sends T -> -T.
    # Retaining both signs would score every curve twice.
    distinct_curves: dict[Fraction, CRTNagaoCandidate] = {}
    for candidate in ordered:
        distinct_curves.setdefault(abs(candidate.parameter_t), candidate)
    return tuple(distinct_curves.values())


def _gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def _run_gp_with_group_cleanup(
    executable: str,
    arguments: Sequence[str],
    program: str,
    *,
    timeout: float,
) -> tuple[int, str, str]:
    """Run one GP process group and reap it on timeout or Python interruption."""

    process = subprocess.Popen(
        [executable, *arguments],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(program, timeout=timeout)
    except BaseException:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.communicate(timeout=2.0)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.communicate()
        raise
    return process.returncode, stdout, stderr


def score_candidates_with_pari(
    candidates: Sequence[ScoreCandidate],
    *,
    cutoff: int,
    timeout: float,
    stack_bytes: int,
) -> dict[str, dict[str, Any]]:
    """Compute the declared good-prime score for a batch in one PARI process."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if not candidates or cutoff < 5 or timeout <= 0 or stack_bytes < 8_000_000:
        raise ValueError("invalid PARI score bounds")
    last_prime = primes_up_to(cutoff)[-1]
    commands = ["default(realprecision,80);"]
    for index, candidate in enumerate(candidates):
        vector = ",".join(_gp_rational(value) for value in candidate.coefficients)
        commands.extend(
            (
                f"E=ellminimalmodel(ellinit([{vector}]));",
                "S=0;USED=0;BAD=0;",
                (
                    f"forprime(p=5,{cutoff},"
                    "if(valuation(E.disc,p)>0,BAD++,"
                    "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++));"
                ),
                f'print("ROW|{index}|",S,"|",USED,"|",BAD,"|{last_prime}");',
            )
        )
    commands.append("quit")
    returncode, stdout, stderr = _run_gp_with_group_cleanup(
        executable,
        ("-q", "-s", str(stack_bytes)),
        "\n".join(commands) + "\n",
        timeout=timeout,
    )
    if returncode != 0 or "***" in stderr:
        raise RuntimeError(f"PARI/GP score failed: {stderr.strip()}")
    answer: dict[str, dict[str, Any]] = {}
    for line in stdout.splitlines():
        if not line.startswith("ROW|"):
            continue
        _, index_text, score, used, bad, observed_last = line.split("|")
        candidate = candidates[int(index_text)]
        answer[candidate.identifier] = {
            "score": score,
            "good_primes_used": int(used),
            "bad_primes_skipped": int(bad),
            "last_numerical_prime": int(observed_last),
        }
    if len(answer) != len(candidates):
        raise RuntimeError("PARI omitted a score record")
    return answer


def parse_integer_tuple(value: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not answer or any(integer < 1 for integer in answer):
        raise argparse.ArgumentTypeError("all values must be positive")
    return answer


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficient-radius", type=int, default=6)
    parser.add_argument("--representatives-per-class", type=int, default=1)
    parser.add_argument("--height-proxy-keep", type=int, default=8)
    parser.add_argument("--stages", type=parse_integer_tuple, default=(200, 1000))
    parser.add_argument("--keep-counts", type=parse_integer_tuple, default=(64, 8))
    parser.add_argument("--score-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=45.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_local_crt.json"
        ),
    )
    return parser


def _candidate_score_input(candidate: CRTNagaoCandidate) -> ScoreCandidate:
    return ScoreCandidate(candidate.identifier, candidate.coefficients)


def benchmark_candidates() -> tuple[tuple[str, ScoreCandidate, str], ...]:
    benchmarks = []
    for parameter_u in (Q(1), Q(42), Q(128)):
        identifier = f"nagao-u-{parameter_u.numerator}"
        benchmarks.append(
            (
                identifier,
                ScoreCandidate(
                    identifier,
                    rank13_base_changed_short_jacobian_coefficients(parameter_u),
                ),
                "Nagao base-changed generic-rank-13 family",
            )
        )
    benchmarks.extend(
        (
            (
                "fermigier-1666-9",
                ScoreCandidate(
                    "fermigier-1666-9",
                    FermigierMestreFamily.coefficients(Q(1666, 9)),
                ),
                "Fermigier conductor-only candidate",
            ),
            (
                "published-e22",
                ScoreCandidate(
                    "published-e22",
                    FermigierMestreFamily.coefficients(NORMALIZED_RECORD_PARAMETER),
                ),
                "published rank-at-least-22 benchmark",
            ),
        )
    )
    return tuple(benchmarks)


def main() -> None:
    args = build_parser().parse_args()
    if len(args.stages) != len(args.keep_counts):
        raise SystemExit("--stages and --keep-counts must have equal length")
    if any(left >= right for left, right in zip(args.stages, args.stages[1:])):
        raise SystemExit("--stages must be strictly increasing")
    if any(left < right for left, right in zip(args.keep_counts, args.keep_counts[1:])):
        raise SystemExit("--keep-counts must be nonincreasing")
    if args.height_proxy_keep < 0:
        raise SystemExit("--height-proxy-keep must be nonnegative")

    discoveries = default_local_discoveries()
    candidates = enumerate_crt_candidates(
        coefficient_radius=args.coefficient_radius,
        representatives_per_class=args.representatives_per_class,
    )
    if not candidates:
        raise RuntimeError("CRT/lattice reconstruction produced no candidates")
    if args.keep_counts[0] > len(candidates):
        raise SystemExit("the first keep count exceeds the candidate population")
    by_identifier = {candidate.identifier: candidate for candidate in candidates}
    height_proxy_candidates = sorted(
        candidates,
        key=lambda candidate: (
            max(abs(candidate.parameter_t.numerator), candidate.parameter_t.denominator),
            candidate.height,
            candidate.identifier,
        ),
    )[: args.height_proxy_keep]
    height_proxy_ids = {
        candidate.identifier for candidate in height_proxy_candidates
    }
    benchmarks = benchmark_candidates()
    benchmark_stage_scores: dict[str, dict[str, Any]] = {
        identifier: {} for identifier, _, _ in benchmarks
    }

    survivors = list(candidates)
    stages: list[dict[str, Any]] = []
    latest_scores: dict[str, dict[str, Any]] = {}
    for cutoff, keep_count in zip(args.stages, args.keep_counts):
        score_input = tuple(
            [_candidate_score_input(candidate) for candidate in survivors]
            + [candidate for _, candidate, _ in benchmarks]
        )
        scores = score_candidates_with_pari(
            score_input,
            cutoff=cutoff,
            timeout=args.score_timeout,
            stack_bytes=args.stack_bytes,
        )
        ranked = sorted(
            (
                {
                    "candidate_id": candidate.identifier,
                    "parameter_u": rational_to_string(candidate.parameter_u),
                    "height": candidate.height,
                    **scores[candidate.identifier],
                }
                for candidate in survivors
            ),
            key=lambda record: (
                -Decimal(record["score"]),
                record["height"],
                record["candidate_id"],
            ),
        )
        retained = ranked[: min(keep_count, len(ranked))]
        retained_ids = {record["candidate_id"] for record in retained}
        retained.extend(
            record
            for record in ranked
            if record["candidate_id"] in height_proxy_ids
            and record["candidate_id"] not in retained_ids
        )
        for identifier, _, _ in benchmarks:
            benchmark_stage_scores[identifier][str(cutoff)] = scores[identifier]
        for record in retained:
            latest_scores[record["candidate_id"]] = record
        stages.append(
            {
                "cutoff": cutoff,
                "population_scored": len(ranked),
                "keep_count": len(retained),
                "retained": retained,
            }
        )
        survivors = [by_identifier[record["candidate_id"]] for record in retained]

    conductor_records = []
    conductor_errors = []
    for candidate in survivors:
        try:
            curve = minimal_curve_data(
                candidate.coefficients,
                timeout=args.conductor_timeout,
                stack_bytes=args.stack_bytes,
                local_primes=DEFAULT_CRT_PRIMES,
            )
            local_checks = []
            for choice in candidate.choices:
                observed = rational_discriminant_valuation(
                    candidate.numerator, candidate.denominator, choice.prime
                )
                local_checks.append(
                    {
                        "prime": choice.prime,
                        "ball_exponent": choice.exponent,
                        "ball_residue": choice.residue,
                        "forced_discriminant_valuation": (
                            choice.forced_discriminant_valuation
                        ),
                        "observed_discriminant_valuation": observed,
                        "reduction_proved_on_ball": choice.reduction,
                        "conductor_exponent_proved_on_ball": choice.conductor_exponent,
                    }
                )
            score = latest_scores[candidate.identifier]
            conductor_records.append(
                {
                    **score,
                    "parameter_t": rational_to_string(candidate.parameter_t),
                    "crt_residue": candidate.crt_residue,
                    "crt_modulus": candidate.crt_modulus,
                    "local_checks": local_checks,
                    "curve": curve,
                    "below_strict_log_conductor_target": (
                        Decimal(curve["log_conductor"]) < TARGET_LOG_CONDUCTOR
                    ),
                    "rank_claim": None,
                    "generic_family_context": (
                        "Nagao cites generic rank at least 13 after base change; "
                        "specialized independence and extra rank were not checked"
                    ),
                }
            )
        except Exception as error:
            conductor_errors.append(
                {"candidate_id": candidate.identifier, "error": str(error)}
            )
    conductor_records.sort(
        key=lambda record: (
            not record["below_strict_log_conductor_target"],
            -Decimal(record["score"]),
            Decimal(record["curve"]["log_conductor"]),
        )
    )

    benchmark_records = []
    for identifier, candidate, provenance in benchmarks:
        curve = minimal_curve_data(
            candidate.coefficients,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        benchmark_records.append(
            {
                "candidate_id": identifier,
                "provenance": provenance,
                "scores": benchmark_stage_scores[identifier],
                "curve": curve,
            }
        )

    search_space_size = 1
    for balls in default_crt_balls().values():
        search_space_size *= len(balls)
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded experiment: exact local root balls and CRT/lattice replay; "
            "scores and conductors are computational triage, with no rank claim"
        ),
        "primary_source": PRIMARY_SOURCE,
        "strict_target": {"rank_at_least": 21, "log_conductor_less_than": "182.72"},
        "target_hits": [],
        "local_discoveries": [discovery.as_json() for discovery in discoveries],
        "crt_search": {
            "primes": list(DEFAULT_CRT_PRIMES),
            "certified_reduction_required": "split multiplicative",
            "ball_choice_combinations": search_space_size,
            "coefficient_radius": args.coefficient_radius,
            "representatives_per_class": args.representatives_per_class,
            "sign_pairs_identified_via_even_T_family": True,
            "distinct_curves": len(candidates),
            "height_proxy_keep": args.height_proxy_keep,
            "height_minimum": min(candidate.height for candidate in candidates),
            "height_maximum": max(candidate.height for candidate in candidates),
            "local_power_savings_proxy": sum(
                (
                    min(ball.forced_discriminant_valuation for ball in balls) - 1
                )
                * log(prime)
                for prime, balls in default_crt_balls().items()
            ),
        },
        "score_definition": (
            "sum over good numerical primes 5<=p<=cutoff of "
            "((2-a_p)/(p+1-a_p))*log(p)"
        ),
        "stages": stages,
        "final_conductor_records": conductor_records,
        "conductor_errors": conductor_errors,
        "benchmarks": benchmark_records,
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "reproducing_command": REPRODUCING_COMMAND,
        "invocation": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"CRT classes={search_space_size}, rational candidates={len(candidates)}, "
        f"conductor survivors={len(conductor_records)}"
    )
    if conductor_records:
        best = conductor_records[0]
        print(
            f"best={best['candidate_id']} u={best['parameter_u']} "
            f"score={best['score']} log(N)={best['curve']['log_conductor']}"
        )


if __name__ == "__main__":
    main()
