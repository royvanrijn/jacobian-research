#!/usr/bin/env python3
"""Scan the multiple-root residue class containing Fermigier's record.

The normalized parameter ``T=a/b`` is constrained by

``T=0 (mod 7), T=0 (mod 17), T=33 (mod 37)``,

or equivalently ``a=2142*b (mod 4403)``.  The scan is exhaustive only in
the declared projective-height box.  Its good-prime score is a heuristic;
neither a high score nor a numerical height-pairing determinant proves a
Mordell--Weil rank.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import heapq
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterator

from ek_k3 import primes_up_to, rational_to_string, valuation
from fermigier_mestre import (
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
    PUBLISHED_PARAMETER,
    ROOTS,
)
from multiple_root_lifting import (
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from pari_bridge import minimal_curve_data, pari_version


RESIDUE_CONSTRAINTS = ((7, 0), (17, 0), (37, 33))
CRT_RESIDUE = 2142
CRT_MODULUS = 4403
PUBLISHED_RANK_LOWER_BOUND = 22
TARGET_LOG_CONDUCTOR = Decimal("182.72")
PUBLISHED_BENCHMARK_LOG_CONDUCTOR = Decimal(
    "182.724910950637428796"
)


@dataclass(frozen=True)
class ScoreCell:
    """Precomputed good-reduction data for one parameter residue."""

    trace: int
    point_count: int
    term: float


@dataclass(frozen=True)
class PrimeScoreTable:
    """Score and inverse lookup tables for one prime."""

    prime: int
    inverses: tuple[int | None, ...]
    cells: tuple[ScoreCell | None, ...]


def score_term(trace: int, point_count: int, prime: int, score: str) -> float:
    """Return one good-prime score summand."""

    if score == "fermigier-good":
        return (2 - trace) / point_count * log(prime)
    if score == "nagao-log":
        return log(point_count / prime)
    raise ValueError(f"unknown score {score}")


def build_score_tables(bound: int, score: str) -> tuple[PrimeScoreTable, ...]:
    """Precompute good-reduction score cells for every ``5 <= p <= bound``."""

    if bound < 5:
        raise ValueError("the score bound must be at least 5")
    tables: list[PrimeScoreTable] = []
    for prime in primes_up_to(bound):
        if prime < 5:
            continue
        inverses: list[int | None] = [None]
        inverses.extend(pow(value, -1, prime) for value in range(1, prime))
        cells: list[ScoreCell | None] = []
        for residue in range(prime):
            local = FermigierMestreFamily.local_data(residue, prime)
            if not local.good_reduction:
                cells.append(None)
                continue
            cells.append(
                ScoreCell(
                    trace=local.trace,
                    point_count=local.point_count,
                    term=score_term(
                        local.trace,
                        local.point_count,
                        prime,
                        score,
                    ),
                )
            )
        tables.append(PrimeScoreTable(prime, tuple(inverses), tuple(cells)))
    return tuple(tables)


def score_rational(
    numerator: int,
    denominator: int,
    tables: tuple[PrimeScoreTable, ...],
    *,
    include_traces: bool = False,
) -> dict[str, Any]:
    """Score a reduced rational, excluding denominator and bad primes.

    A prime dividing ``denominator`` has no affine residue and is counted in
    ``skipped_denominator_primes``.  A residue at which the specialized
    discriminant vanishes is counted in ``skipped_bad_primes`` and contributes
    no good-reduction score term.
    """

    total = 0.0
    used = 0
    skipped_denominator = 0
    skipped_bad = 0
    traces: list[dict[str, int]] = []
    for table in tables:
        prime = table.prime
        denominator_mod_prime = denominator % prime
        inverse = table.inverses[denominator_mod_prime]
        if inverse is None:
            skipped_denominator += 1
            continue
        residue = (numerator % prime) * inverse % prime
        cell = table.cells[residue]
        if cell is None:
            skipped_bad += 1
            continue
        total += cell.term
        used += 1
        if include_traces:
            traces.append(
                {
                    "prime": prime,
                    "trace": cell.trace,
                }
            )
    answer: dict[str, Any] = {
        "value": total,
        "primes_used": used,
        "skipped_denominator_primes": skipped_denominator,
        "skipped_bad_primes": skipped_bad,
    }
    if include_traces:
        answer["traces"] = traces
    return answer


def primitive_rationals_in_class(height: int) -> Iterator[tuple[int, int]]:
    """Yield every reduced ``a/b`` in the residue class with ``H(a:b)<=height``.

    The canonical representative has ``b>0`` and
    ``H(a:b)=max(abs(a),b)``.  Enumeration is by increasing denominator and
    then increasing numerator; no rational appears twice.
    """

    if height < 1:
        raise ValueError("the projective height bound must be positive")

    def ceiling_division(numerator: int, denominator: int) -> int:
        return -((-numerator) // denominator)

    for denominator in range(1, height + 1):
        residue = CRT_RESIDUE * denominator % CRT_MODULUS
        minimum_multiplier = ceiling_division(-height - residue, CRT_MODULUS)
        maximum_multiplier = (height - residue) // CRT_MODULUS
        for multiplier in range(minimum_multiplier, maximum_multiplier + 1):
            numerator = residue + multiplier * CRT_MODULUS
            if gcd(numerator, denominator) != 1:
                continue
            if gcd(denominator, CRT_MODULUS) != 1:
                raise AssertionError(
                    "a primitive constrained denominator was not a unit"
                )
            yield numerator, denominator


def divide_by_linear_mod(
    coefficients: tuple[int, ...], residue: int, prime: int
) -> tuple[tuple[int, ...], int]:
    """Synthetic division by ``T-residue`` over ``F_prime``.

    Coefficients are in ascending order.  The returned second component is
    the remainder.
    """

    if len(coefficients) < 2:
        return (), coefficients[0] % prime
    reduced = tuple(coefficient % prime for coefficient in coefficients)
    quotient = [0] * (len(reduced) - 1)
    quotient[-1] = reduced[-1]
    for index in range(len(reduced) - 2, 0, -1):
        quotient[index - 1] = (
            reduced[index] + residue * quotient[index]
        ) % prime
    remainder = (reduced[0] + residue * quotient[0]) % prime
    return tuple(quotient), remainder


def root_multiplicity_mod(prime: int, residue: int) -> int:
    """Return the exact multiplicity of a root of ``H(T)`` modulo ``prime``."""

    coefficients = DISCRIMINANT_FACTOR_COEFFICIENTS
    multiplicity = 0
    while len(coefficients) > 1:
        quotient, remainder = divide_by_linear_mod(coefficients, residue, prime)
        if remainder:
            break
        multiplicity += 1
        coefficients = quotient
    return multiplicity


def ranking_key(
    score_data: dict[str, Any], numerator: int, denominator: int
) -> tuple[Any, ...]:
    """Return the deterministic ascending key used for final score ranking."""

    return (
        -score_data["value"],
        max(abs(numerator), denominator),
        numerator,
        denominator,
    )


def exact_candidate_record(
    numerator: int,
    denominator: int,
    tables: tuple[PrimeScoreTable, ...],
    minimum_valuations: dict[int, int],
) -> dict[str, Any]:
    """Build a retained-candidate record with exact ``H(T)`` valuations."""

    parameter = Fraction(numerator, denominator)
    h_value = FermigierMestreFamily.discriminant_factor(parameter)
    if h_value == 0:
        raise ValueError("the retained specialization is singular")
    h_valuations = {
        str(prime): valuation(h_value, prime)
        for prime, _ in RESIDUE_CONSTRAINTS
    }
    for prime, guaranteed in minimum_valuations.items():
        if h_valuations[str(prime)] < guaranteed:
            raise AssertionError(f"the residue class lost its v_{prime}(H) guarantee")
    score_data = score_rational(
        numerator,
        denominator,
        tables,
        include_traces=True,
    )
    return {
        "t": rational_to_string(parameter),
        "numerator": numerator,
        "denominator": denominator,
        "height": max(abs(numerator), denominator),
        "congruence_remainder": (
            numerator - CRT_RESIDUE * denominator
        ) % CRT_MODULUS,
        "residues": {
            str(prime): numerator * pow(denominator, -1, prime) % prime
            for prime, _ in RESIDUE_CONSTRAINTS
        },
        "discriminant_factor": rational_to_string(h_value),
        "h_valuations": h_valuations,
        "score": score_data,
        "is_published_benchmark": parameter == NORMALIZED_RECORD_PARAMETER,
    }


def select_for_pari(
    records: list[dict[str, Any]], count: int
) -> list[dict[str, Any]]:
    """Select the first ``count`` top-K records, deterministically."""

    return [record for record in records if record.get("in_top_k", False)][:count]


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=5000)
    parser.add_argument("--score-bound", type=int, default=500)
    parser.add_argument(
        "--score",
        choices=("fermigier-good", "nagao-log"),
        default="fermigier-good",
    )
    parser.add_argument("--keep", type=int, default=50)
    parser.add_argument("--pari-count", type=int, default=0)
    parser.add_argument("--pari-rank-count", type=int, default=0)
    parser.add_argument("--pari-rank-effort", type=int, default=0)
    parser.add_argument("--pari-timeout", type=float, default=60.0)
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_record_residue_class.json"
        ),
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if args.height < 1:
        raise SystemExit("--height must be positive")
    if args.score_bound < 5:
        raise SystemExit("--score-bound must be at least 5")
    if args.keep < 1:
        raise SystemExit("--keep must be positive")
    if args.pari_count < 0 or args.pari_rank_count < 0:
        raise SystemExit("PARI candidate counts must be nonnegative")
    if args.pari_rank_effort < 0:
        raise SystemExit("--pari-rank-effort must be nonnegative")


def main() -> None:
    args = build_parser().parse_args()
    validate_arguments(args)
    if CRT_MODULUS != 7 * 17 * 37:
        raise AssertionError("the pinned CRT modulus changed")
    for prime, residue in RESIDUE_CONSTRAINTS:
        if CRT_RESIDUE % prime != residue:
            raise AssertionError("the pinned CRT residue changed")

    multiplicities = {
        prime: root_multiplicity_mod(prime, residue)
        for prime, residue in RESIDUE_CONSTRAINTS
    }
    if any(multiplicity < 2 for multiplicity in multiplicities.values()):
        raise AssertionError("the pinned class no longer consists of multiple roots")
    minimum_valuations = {
        prime: fixed_divisor_valuation(
            affine_variable_coefficients(
                DISCRIMINANT_FACTOR_COEFFICIENTS,
                residue,
                prime,
            ),
            prime,
        )
        for prime, residue in RESIDUE_CONSTRAINTS
    }
    tables = build_score_tables(args.score_bound, args.score)

    benchmark_height = max(
        abs(NORMALIZED_RECORD_PARAMETER.numerator),
        NORMALIZED_RECORD_PARAMETER.denominator,
    )
    benchmark_eligible = benchmark_height <= args.height
    benchmark_score = None
    benchmark_key = None
    if benchmark_eligible:
        benchmark_score = score_rational(
            NORMALIZED_RECORD_PARAMETER.numerator,
            NORMALIZED_RECORD_PARAMETER.denominator,
            tables,
        )
        benchmark_key = ranking_key(
            benchmark_score,
            NORMALIZED_RECORD_PARAMETER.numerator,
            NORMALIZED_RECORD_PARAMETER.denominator,
        )

    heap: list[tuple[tuple[Any, ...], int, int]] = []
    scanned = 0
    scored = 0
    benchmark_seen = False
    candidates_before_benchmark = 0
    for numerator, denominator in primitive_rationals_in_class(args.height):
        scanned += 1
        score_data = score_rational(numerator, denominator, tables)
        # A good reduction at even one prime proves H(a/b) is nonzero.  This
        # also avoids ranking a singular specialization without evaluating a
        # degree-20 rational polynomial for every point in a large box.
        if score_data["primes_used"] == 0:
            continue
        scored += 1
        key = ranking_key(score_data, numerator, denominator)
        if benchmark_key is not None and key < benchmark_key:
            candidates_before_benchmark += 1
        if (
            numerator == NORMALIZED_RECORD_PARAMETER.numerator
            and denominator == NORMALIZED_RECORD_PARAMETER.denominator
        ):
            benchmark_seen = True

        # Larger quality is better, so the min-heap root is the worst retained
        # candidate.  The final ranking is recomputed from detailed records.
        quality = (
            score_data["value"],
            -max(abs(numerator), denominator),
            -numerator,
            -denominator,
        )
        item = (quality, numerator, denominator)
        if len(heap) < args.keep:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)

    if benchmark_eligible and not benchmark_seen:
        raise AssertionError("the exhaustive box failed to enumerate the benchmark")

    top_pairs = [(item[1], item[2]) for item in heap]
    records = [
        exact_candidate_record(
            numerator,
            denominator,
            tables,
            minimum_valuations,
        )
        for numerator, denominator in top_pairs
    ]
    records.sort(
        key=lambda record: ranking_key(
            record["score"], record["numerator"], record["denominator"]
        )
    )
    for position, record in enumerate(records, 1):
        record["score_position"] = position
        record["in_top_k"] = True

    benchmark_record = None
    if benchmark_eligible:
        benchmark_record = next(
            (record for record in records if record["is_published_benchmark"]),
            None,
        )
        if benchmark_record is None:
            benchmark_record = exact_candidate_record(
                NORMALIZED_RECORD_PARAMETER.numerator,
                NORMALIZED_RECORD_PARAMETER.denominator,
                tables,
                minimum_valuations,
            )
            benchmark_record["in_top_k"] = False
            records.append(benchmark_record)
        benchmark_record["score_position"] = candidates_before_benchmark + 1
        benchmark_record["published_rank_metadata"] = {
            "lower_bound": PUBLISHED_RANK_LOWER_BOUND,
            "source": "Fermigier 1996; not recomputed by this scan",
        }

    pari_errors: list[dict[str, str]] = []
    selected = select_for_pari(records, args.pari_count)
    rank_runs = 0
    for record in selected:
        parameter = Fraction(record["numerator"], record["denominator"])
        rank_effort = None
        if (
            rank_runs < args.pari_rank_count
            and not record["is_published_benchmark"]
        ):
            rank_effort = args.pari_rank_effort
            rank_runs += 1
        try:
            # PARI checks all twelve exact images and records their numerical
            # height determinant even when no descent/rank call is requested.
            # Neither operation is promoted to an independence proof.
            known_points = FermigierMestreFamily.known_jacobian_points(
                parameter
            )[1:]
            record["pari"] = minimal_curve_data(
                FermigierMestreFamily.coefficients(parameter),
                timeout=args.pari_timeout,
                rank_effort=rank_effort,
                known_points=known_points,
                stack_bytes=args.pari_stack_bytes,
            )
            record["below_log_conductor_target"] = (
                Decimal(record["pari"]["log_conductor"])
                < TARGET_LOG_CONDUCTOR
            )
        except Exception as error:
            pari_errors.append({"t": record["t"], "error": str(error)})

    target_hits: list[dict[str, Any]] = []
    for record in records:
        computed_rank = record.get("pari", {}).get("pari_ellrank")
        if computed_rank is None:
            continue
        lower_bound = computed_rank["lower_bound"]
        if (
            lower_bound >= 21
            and record.get("below_log_conductor_target", False)
        ) or lower_bound >= 30:
            target_hits.append(record)

    benchmark_summary: dict[str, Any] = {
        "normalized_parameter": rational_to_string(NORMALIZED_RECORD_PARAMETER),
        "published_parameter": rational_to_string(PUBLISHED_PARAMETER),
        "height": benchmark_height,
        "inside_height_box": benchmark_eligible,
        "inside_residue_class": (
            NORMALIZED_RECORD_PARAMETER.numerator
            - CRT_RESIDUE * NORMALIZED_RECORD_PARAMETER.denominator
        )
        % CRT_MODULUS
        == 0,
        "published_rank_lower_bound": PUBLISHED_RANK_LOWER_BOUND,
        "rank_source": "published metadata only; no rank value inferred by this scan",
        "published_log_conductor": str(PUBLISHED_BENCHMARK_LOG_CONDUCTOR),
        "strict_small_conductor_target_met": (
            PUBLISHED_BENCHMARK_LOG_CONDUCTOR < TARGET_LOG_CONDUCTOR
        ),
        "included_in_candidates": benchmark_record is not None,
        "score_position_among_scored_candidates": (
            benchmark_record["score_position"] if benchmark_record else None
        ),
        "in_top_k": (
            benchmark_record["in_top_k"] if benchmark_record else None
        ),
    }

    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exhaustive residue-class experiment; good-prime scores are "
            "heuristics, PARI rank fields are computational bounds only, and the "
            "benchmark rank is published metadata"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "family": {
            "name": "normalized Fermigier--Mestre rank-at-least-12 family",
            "root_tuple": list(ROOTS),
            "source": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        },
        "residue_class": {
            "local_constraints": [
                {"prime": prime, "residue": residue}
                for prime, residue in RESIDUE_CONSTRAINTS
            ],
            "crt_residue": CRT_RESIDUE,
            "crt_modulus": CRT_MODULUS,
            "root_multiplicities_mod_p": {
                str(prime): multiplicity
                for prime, multiplicity in multiplicities.items()
            },
            "guaranteed_h_valuations": {
                str(prime): minimum
                for prime, minimum in minimum_valuations.items()
            },
        },
        "method": {
            "enumeration": (
                "all reduced a/b with b>0, max(abs(a),b)<=height, "
                "and a=2142*b mod 4403"
            ),
            "score": args.score,
            "score_scope": (
                "good reduction only at numerical primes 5<=p<=score_bound; "
                "primes dividing b and bad-reduction residues are excluded"
            ),
            "exact_postprocessing": (
                "H(T) is evaluated exactly for retained candidates and its "
                "7-, 17-, and 37-adic valuations are recorded"
            ),
            "rank_seed": (
                "optional PARI ellrank calls receive exact Jacobian images "
                "2--13; their height determinant is numerical evidence only"
            ),
        },
        "parameters": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in vars(args).items()
        },
        "counts": {
            "primitive_candidates_scanned": scanned,
            "candidates_with_at_least_one_good_score_prime": scored,
            "candidates_without_a_good_score_prime": scanned - scored,
            "top_k_requested": args.keep,
            "records_output": len(records),
        },
        "benchmark": benchmark_summary,
        "pari_errors": pari_errors,
        "target_hits": target_hits,
        "candidates": records,
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"scanned={scanned} scored={scored} retained={len(records)} "
        f"target_hits={len(target_hits)}"
    )
    if benchmark_record is not None:
        print(
            "benchmark score position "
            f"{benchmark_record['score_position']}/{scored}; "
            f"top_k={benchmark_record['in_top_k']}"
        )
    for record in selected:
        if "pari" in record:
            print(
                f"t={record['t']} score={record['score']['value']:.6f} "
                f"logN={record['pari']['log_conductor']}"
            )


if __name__ == "__main__":
    main()
