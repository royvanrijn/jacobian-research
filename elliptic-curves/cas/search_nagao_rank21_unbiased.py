#!/usr/bin/env python3
"""Unbiased rational score-and-point search in Nagao's rank-21 family.

This program exhausts every positive primitive ``T=a/b`` in one declared
rectangle.  It imposes no CRT, conductor, or previously successful residue
condition.  Candidate selection is leakage-free:

1. exact finite-field lookup tables score the complete projective parameter
   line at every prime ``5 <= p <= 200``;
2. the retained residue-score tail is recomputed by PARI through ``p=2000``;
3. exact quartic points are searched in uniform ``H=50k,250k,1m`` stages;
4. two-precision height matrices rank only the exact point pools returned by
   the preceding stage; and
5. every stable numerical rank at least 18 triggers a saturation,
   finite-reduction independence search, and exact conductor replay.

The known parameters 531/2, 956/9, 1637/12, and 5777/32 are excluded from all
selection populations and run separately as calibration.  Numerical height
rank and bounded-search failures are never promoted to rank theorems.
Every PARI child is synchronous, joined, and protected by a finite timeout.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import heapq
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from ek_k3 import legendre_symbol, primes_up_to, rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    short_jacobian_coefficients,
)
from nagao_rank13_local import (
    polynomial_add,
    polynomial_multiply,
    polynomial_scale,
)
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank21_neighborhood import (
    ScoreCandidate,
    score_candidates_with_pari,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    stable_height_rank,
)
from triage_nagao_rank21_neighbor import (
    bounded_quartic_points,
    exact_visible_seeds,
    map_and_deduplicate,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
CALIBRATION_PARAMETERS = (Q(531, 2), Q(956, 9), Q(1637, 12), Q(5777, 32))
DEFAULT_A_MAX = 10_000
DEFAULT_B_MAX = 100
RESIDUE_CUTOFF = 200
EXACT_SCORE_CUTOFF = 2_000
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_unbiased.py"
)


# Ascending coefficients of the primitive rank-21 quartic
# e(T)+d(T)X+c(T)X^2+b(T)X^3+a(T)X^4.  The formulas are independently
# checked against the constructor in the tests below this module.
QUARTIC_E = (197_804_341_504, 0, -5_788_761_232, 0, 4_306_112, 0, -23)
QUARTIC_D = (-37_842_168_672, 0, 98_440_030, 0, -18_124)
QUARTIC_C = (6_945_772_145, 0, -3_814_698, 0, 46)
QUARTIC_B = (-35_986_906, 0, 18_124)
QUARTIC_A = (47_342, 0, -23)


def rank21_invariant_polynomials() -> tuple[tuple[int, ...], tuple[int, ...]]:
    invariant_i = polynomial_add(
        polynomial_scale(polynomial_multiply(QUARTIC_A, QUARTIC_E), 12),
        polynomial_scale(polynomial_multiply(QUARTIC_B, QUARTIC_D), -3),
        polynomial_multiply(QUARTIC_C, QUARTIC_C),
    )
    invariant_j = polynomial_add(
        polynomial_scale(
            polynomial_multiply(QUARTIC_A, QUARTIC_C, QUARTIC_E), 72
        ),
        polynomial_scale(
            polynomial_multiply(QUARTIC_B, QUARTIC_C, QUARTIC_D), 9
        ),
        polynomial_scale(
            polynomial_multiply(QUARTIC_A, QUARTIC_D, QUARTIC_D), -27
        ),
        polynomial_scale(
            polynomial_multiply(QUARTIC_B, QUARTIC_B, QUARTIC_E), -27
        ),
        polynomial_scale(
            polynomial_multiply(QUARTIC_C, QUARTIC_C, QUARTIC_C), -2
        ),
    )
    if len(invariant_i) != 9 or len(invariant_j) != 13:
        raise AssertionError("the rank-21 invariant degrees changed")
    return invariant_i, invariant_j


INVARIANT_I, INVARIANT_J = rank21_invariant_polynomials()


def polynomial_value(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def homogeneous_value_mod(
    coefficients: Sequence[int], numerator: int, denominator: int, prime: int
) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient
        * pow(numerator, power, prime)
        * pow(denominator, degree - power, prime)
        for power, coefficient in enumerate(coefficients)
    ) % prime


@dataclass(frozen=True)
class ResidueSymbol:
    prime: int
    projective_index: int
    ellap: int | None
    good_reduction: bool
    contribution: float

    @property
    def label(self) -> str:
        return "infinity" if self.projective_index == self.prime else str(self.projective_index)


def residue_table(prime: int) -> tuple[ResidueSymbol, ...]:
    """Score every point of ``P^1(F_p)``, with infinity stored at index p."""

    if prime < 5 or prime not in primes_up_to(prime):
        raise ValueError("the score modulus must be a prime at least five")
    symbols = []
    for projective_index in range(prime + 1):
        numerator, denominator = (
            (projective_index, 1)
            if projective_index < prime
            else (1, 0)
        )
        coefficient_a = (
            -27
            * homogeneous_value_mod(
                INVARIANT_I, numerator, denominator, prime
            )
            % prime
        )
        coefficient_b = (
            -27
            * homogeneous_value_mod(
                INVARIANT_J, numerator, denominator, prime
            )
            % prime
        )
        discriminant_core = (
            4 * coefficient_a**3 + 27 * coefficient_b**2
        ) % prime
        if discriminant_core == 0:
            symbols.append(
                ResidueSymbol(prime, projective_index, None, False, 0.0)
            )
            continue
        trace = -sum(
            legendre_symbol(
                (x_value**3 + coefficient_a * x_value + coefficient_b) % prime,
                prime,
            )
            for x_value in range(prime)
        )
        contribution = (2 - trace) / (prime + 1 - trace) * log(prime)
        symbols.append(
            ResidueSymbol(prime, projective_index, trace, True, contribution)
        )
    return tuple(symbols)


def build_residue_tables(cutoff: int = RESIDUE_CUTOFF) -> dict[int, tuple[ResidueSymbol, ...]]:
    if cutoff < 5:
        raise ValueError("the residue cutoff must be at least five")
    return {
        prime: residue_table(prime)
        for prime in primes_up_to(cutoff)
        if prime >= 5
    }


def projective_index(numerator: int, denominator: int, prime: int) -> int:
    if denominator % prime == 0:
        if numerator % prime == 0:
            raise ValueError("a nonprimitive pair has no projective reduction")
        return prime
    return numerator * pow(denominator, -1, prime) % prime


def residue_score(
    numerator: int,
    denominator: int,
    tables: dict[int, tuple[ResidueSymbol, ...]],
) -> tuple[float, int, int]:
    score = 0.0
    good_count = 0
    bad_count = 0
    for prime, table in tables.items():
        symbol = table[projective_index(numerator, denominator, prime)]
        if symbol.good_reduction:
            score += symbol.contribution
            good_count += 1
        else:
            bad_count += 1
    return score, good_count, bad_count


@dataclass(frozen=True)
class PrefilterCandidate:
    numerator: int
    denominator: int
    residue_score_b200: float
    residue_good_primes: int
    residue_bad_primes: int

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)

    @property
    def identifier(self) -> str:
        return f"unbiased-{self.numerator}-{self.denominator}"


def primitive_population_count(a_max: int, b_max: int) -> int:
    if a_max < 1 or b_max < 1:
        raise ValueError("population bounds must be positive")
    return sum(
        gcd(numerator, denominator) == 1
        for denominator in range(1, b_max + 1)
        for numerator in range(1, a_max + 1)
    )


def prefilter_population(
    *,
    a_max: int,
    b_max: int,
    keep_count: int,
    tables: dict[int, tuple[ResidueSymbol, ...]],
    excluded_parameters: Iterable[Fraction] = CALIBRATION_PARAMETERS,
) -> tuple[tuple[PrefilterCandidate, ...], dict[str, Any]]:
    """Exhaust the rectangle and retain its deterministic score tail."""

    if keep_count < 1:
        raise ValueError("the prefilter keep count must be positive")
    excluded = {Q(parameter) for parameter in excluded_parameters}
    heap: list[tuple[float, int, int, int, int, int, int, int]] = []
    complete_count = 0
    eligible_count = 0
    excluded_seen: list[str] = []
    digest = hashlib.sha256()
    primes = tuple(tables)
    for denominator in range(1, b_max + 1):
        # Denominator-specific projective reductions avoid repeated inversions.
        denominator_maps = {
            prime: (
                None
                if denominator % prime == 0
                else pow(denominator, -1, prime)
            )
            for prime in primes
        }
        for numerator in range(1, a_max + 1):
            if gcd(numerator, denominator) != 1:
                continue
            complete_count += 1
            parameter = Q(numerator, denominator)
            if parameter in excluded:
                excluded_seen.append(rational_to_string(parameter))
                continue
            eligible_count += 1
            score = 0.0
            good_count = 0
            bad_count = 0
            for prime in primes:
                inverse = denominator_maps[prime]
                index = (
                    prime
                    if inverse is None
                    else numerator % prime * inverse % prime
                )
                symbol = tables[prime][index]
                if symbol.good_reduction:
                    score += symbol.contribution
                    good_count += 1
                else:
                    bad_count += 1
            digest.update(
                f"{numerator}/{denominator}|{score.hex()}|{good_count}|{bad_count}\n".encode()
            )
            # Larger tuple is better: score first, then smaller height/a/b.
            key = (
                score,
                -max(numerator, denominator),
                -numerator,
                -denominator,
                numerator,
                denominator,
                good_count,
                bad_count,
            )
            if len(heap) < keep_count:
                heapq.heappush(heap, key)
            elif key > heap[0]:
                heapq.heapreplace(heap, key)
    if len(heap) != keep_count:
        raise ValueError("the prefilter keep count exceeds the eligible population")
    retained = tuple(
        PrefilterCandidate(
            numerator=item[4],
            denominator=item[5],
            residue_score_b200=item[0],
            residue_good_primes=item[6],
            residue_bad_primes=item[7],
        )
        for item in sorted(heap, reverse=True)
    )
    return retained, {
        "complete_primitive_population_count": complete_count,
        "excluded_calibration_parameters_seen": sorted(excluded_seen),
        "eligible_population_count": eligible_count,
        "eligible_population_score_stream_sha256": digest.hexdigest(),
        "retained_count": len(retained),
        "retention_boundary_score": retained[-1].residue_score_b200,
        "best_score": retained[0].residue_score_b200,
    }


@dataclass(frozen=True)
class ExactScoredCandidate:
    prefilter: PrefilterCandidate
    exact_score_b2000: str
    exact_good_primes: int
    exact_bad_primes: int
    exact_last_prime: int

    @property
    def parameter(self) -> Fraction:
        return self.prefilter.parameter

    @property
    def identifier(self) -> str:
        return self.prefilter.identifier


def exact_score_candidates(
    candidates: Sequence[PrefilterCandidate],
    *,
    cutoff: int,
    batch_size: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[ExactScoredCandidate, ...]:
    if batch_size < 1:
        raise ValueError("the exact-score batch size must be positive")
    records: dict[str, dict[str, Any]] = {}
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        inputs = tuple(
            ScoreCandidate(
                candidate.identifier,
                short_jacobian_coefficients(
                    RANK21_CONSTRUCTION, candidate.parameter
                ),
            )
            for candidate in batch
        )
        observed = score_candidates_with_pari(
            inputs,
            cutoff=cutoff,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        overlap = set(records) & set(observed)
        if overlap:
            raise AssertionError("an exact-score candidate was repeated")
        records.update(observed)
    if len(records) != len(candidates):
        raise AssertionError("PARI omitted an exact-score candidate")
    answer = tuple(
        ExactScoredCandidate(
            prefilter=candidate,
            exact_score_b2000=records[candidate.identifier]["score"],
            exact_good_primes=records[candidate.identifier]["good_primes_used"],
            exact_bad_primes=records[candidate.identifier]["bad_primes_skipped"],
            exact_last_prime=records[candidate.identifier]["last_numerical_prime"],
        )
        for candidate in candidates
    )
    return tuple(
        sorted(
            answer,
            key=lambda candidate: (
                -Decimal(candidate.exact_score_b2000),
                max(candidate.prefilter.numerator, candidate.prefilter.denominator),
                candidate.prefilter.numerator,
                candidate.prefilter.denominator,
            ),
        )
    )


@dataclass(frozen=True)
class PointPool:
    candidate: ExactScoredCandidate
    height_bound: int
    status: str
    signed_points: int
    signless_points: int
    visible_abscissas_returned: int
    new_images: tuple[tuple[Fraction, Fraction], ...]
    seed_images: tuple[tuple[Fraction, Fraction], ...]
    coefficients: tuple[Fraction, ...]
    wall_seconds: float
    pari_milliseconds: int
    error: str | None = None

    @property
    def pool(self) -> tuple[tuple[Fraction, Fraction], ...]:
        return self.seed_images + self.new_images


def search_one_pool(
    candidate: ExactScoredCandidate,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> PointPool:
    started = time.monotonic()
    try:
        seed_quartic, seed_images, coefficients = exact_visible_seeds(
            candidate.parameter
        )
        raw_points, wall_seconds, milliseconds = bounded_quartic_points(
            candidate.parameter,
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        mapped, new_images, _ = map_and_deduplicate(
            candidate.parameter, raw_points, seed_quartic, seed_images
        )
        return PointPool(
            candidate=candidate,
            height_bound=height_bound,
            status="completed",
            signed_points=len(raw_points),
            signless_points=len({point[0] for point in raw_points}),
            visible_abscissas_returned=sum(
                bool(record["visible_section_abscissa"]) for record in mapped
            ),
            new_images=new_images,
            seed_images=seed_images,
            coefficients=coefficients,
            wall_seconds=wall_seconds,
            pari_milliseconds=milliseconds,
        )
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return PointPool(
            candidate=candidate,
            height_bound=height_bound,
            status="timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            signed_points=0,
            signless_points=0,
            visible_abscissas_returned=0,
            new_images=(),
            seed_images=(),
            coefficients=(),
            wall_seconds=time.monotonic() - started,
            pari_milliseconds=0,
            error=str(error)[:500],
        )


def parallel_point_search(
    candidates: Sequence[ExactScoredCandidate],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> tuple[PointPool, ...]:
    if workers < 1:
        raise ValueError("the worker count must be positive")
    by_identifier: dict[str, PointPool] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                search_one_pool,
                candidate,
                height_bound=height_bound,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            result = future.result()
            by_identifier[result.candidate.identifier] = result
    if len(by_identifier) != len(candidates):
        raise AssertionError("the point stage lost a candidate")
    return tuple(by_identifier[candidate.identifier] for candidate in candidates)


def rank_one_pool(
    point_pool: PointPool,
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    if point_pool.status != "completed":
        return {"status": "not-run", "reason": "point search did not complete"}
    try:
        runs = height_matrix_replay(
            point_pool.coefficients,
            point_pool.pool,
            precisions=precisions,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        rank = stable_height_rank(runs)
        indices = tuple(runs[-1]["subset_indices_one_based"])
        selected = tuple(point_pool.pool[index - 1] for index in indices)
        return {
            "status": "completed",
            "stable_numerical_rank": rank,
            "precision_runs": list(runs),
            "selected_subset_indices_one_based": list(indices),
            "selected_points": selected,
            "selected_point_sha256": point_digest(selected),
        }
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return {
            "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            "error": str(error)[:500],
        }


def parallel_rank_replay(
    pools: Sequence[PointPool],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                rank_one_pool,
                point_pool,
                precisions=precisions,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): point_pool.candidate.identifier
            for point_pool in pools
        }
        for future in as_completed(futures):
            records[futures[future]] = future.result()
    if len(records) != len(pools):
        raise AssertionError("the height stage lost a candidate")
    return records


def pool_priority(
    point_pool: PointPool, rank_record: dict[str, Any]
) -> tuple[Any, ...]:
    numerical_rank = (
        int(rank_record["stable_numerical_rank"])
        if rank_record.get("status") == "completed"
        else -1
    )
    return (
        rank_record.get("status") != "completed",
        -numerical_rank,
        -len(point_pool.new_images),
        -Decimal(point_pool.candidate.exact_score_b2000),
        max(
            point_pool.candidate.prefilter.numerator,
            point_pool.candidate.prefilter.denominator,
        ),
        point_pool.candidate.identifier,
    )


def exact_score_record(candidate: ExactScoredCandidate) -> dict[str, Any]:
    return {
        "candidate_id": candidate.identifier,
        "constructor_parameter": rational_to_string(candidate.parameter),
        "height": max(candidate.prefilter.numerator, candidate.prefilter.denominator),
        "residue_table_b200": {
            "score": candidate.prefilter.residue_score_b200,
            "good_primes_used": candidate.prefilter.residue_good_primes,
            "bad_primes_skipped": candidate.prefilter.residue_bad_primes,
        },
        "exact_pari_b2000": {
            "score": candidate.exact_score_b2000,
            "good_primes_used": candidate.exact_good_primes,
            "bad_primes_skipped": candidate.exact_bad_primes,
            "last_numerical_prime": candidate.exact_last_prime,
        },
    }


def point_pool_record(
    point_pool: PointPool,
    rank_record: dict[str, Any],
    *,
    include_points: bool,
) -> dict[str, Any]:
    record = {
        **exact_score_record(point_pool.candidate),
        "point_search": {
            "height_bound": point_pool.height_bound,
            "status": point_pool.status,
            "signed_points": point_pool.signed_points,
            "distinct_quartic_abscissas": point_pool.signless_points,
            "visible_abscissas_returned": point_pool.visible_abscissas_returned,
            "new_distinct_jacobian_sign_pairs": len(point_pool.new_images),
            "new_image_sha256": point_digest(point_pool.new_images),
            "wall_seconds": point_pool.wall_seconds,
            "pari_milliseconds": point_pool.pari_milliseconds,
        },
        "height_rank": {
            key: value
            for key, value in rank_record.items()
            if key != "selected_points"
        },
    }
    if point_pool.error is not None:
        record["point_search"]["error"] = point_pool.error
    if include_points:
        record["exact_selected_points"] = [
            {
                "jacobian_x": rational_to_string(point[0]),
                "jacobian_y": rational_to_string(point[1]),
                "exact_membership_checked": True,
            }
            for point in rank_record.get("selected_points", ())
        ]
    return record


def finite_reduction_certificate(
    point_pool: PointPool,
    rank_record: dict[str, Any],
    *,
    saturation_timeout: float,
    certificate_prime_bound: int,
    stack_bytes: int,
) -> dict[str, Any]:
    points = tuple(rank_record["selected_points"])
    saturated, saturation = saturate_exact_basis(
        point_pool.coefficients,
        points,
        prime_bound=20,
        timeout=saturation_timeout,
        stack_bytes=stack_bytes,
    )
    signatures = find_mod2_reduction_certificate(
        point_pool.coefficients,
        saturated,
        prime_bound=certificate_prime_bound,
    )
    binary_rank = combined_mod2_rank(signatures, len(saturated))
    certified = binary_rank == len(saturated)
    return {
        "status": "certified" if certified else "bounded-search-rank-deficient",
        "input_numerical_subset_count": len(points),
        "saturated_point_count": len(saturated),
        "saturation": saturation,
        "saturated_point_sha256": point_digest(saturated),
        "certificate_prime_bound": certificate_prime_bound,
        "certificate_primes": [signature.prime for signature in signatures],
        "combined_exact_rank_over_F2": binary_rank,
        "two_torsion_certificate_prime": (
            find_two_torsion_certificate_prime(point_pool.coefficients)
            if certified
            else None
        ),
        "signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "certified_algebraic_rank_lower_bound": len(saturated) if certified else None,
    }


def parse_positive_ints(value: str) -> tuple[int, ...]:
    try:
        result = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not result or any(item < 1 for item in result):
        raise argparse.ArgumentTypeError("all integers must be positive")
    return result


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-max", type=int, default=DEFAULT_A_MAX)
    parser.add_argument("--b-max", type=int, default=DEFAULT_B_MAX)
    parser.add_argument("--residue-keep", type=int, default=10_000)
    parser.add_argument("--exact-keep", type=int, default=3_000)
    parser.add_argument("--exact-score-batch", type=int, default=500)
    parser.add_argument("--exact-score-timeout", type=float, default=60.0)
    parser.add_argument("--stage-heights", type=parse_positive_ints, default=(50_000, 250_000, 1_000_000))
    parser.add_argument("--stage-keeps", type=parse_positive_ints, default=(160, 16))
    parser.add_argument("--stage-timeouts", type=parse_positive_ints, default=(5, 20, 90))
    parser.add_argument("--stage-workers", type=parse_positive_ints, default=(4, 4, 2))
    parser.add_argument("--height-precisions", type=parse_positive_ints, default=(72, 120))
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts/generated-results/elliptic_nagao_rank21_unbiased.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.a_max < DEFAULT_A_MAX or args.b_max < DEFAULT_B_MAX:
        raise SystemExit("the declared unbiased population may not be reduced below 10000x100")
    if len(args.stage_heights) != 3 or len(args.stage_keeps) != 2:
        raise SystemExit("provide three stage heights and two intermediate keep counts")
    if len(args.stage_timeouts) != 3 or len(args.stage_workers) != 3:
        raise SystemExit("provide one timeout and worker count per point stage")
    if tuple(sorted(set(args.stage_heights))) != args.stage_heights:
        raise SystemExit("point-stage heights must be strictly increasing")
    if tuple(sorted(set(args.height_precisions))) != args.height_precisions:
        raise SystemExit("height precisions must be strictly increasing")
    if args.residue_keep < args.exact_keep or args.stage_keeps[0] > args.exact_keep:
        raise SystemExit("selection keep counts must be nonincreasing")
    if args.stage_keeps[1] > args.stage_keeps[0]:
        raise SystemExit("point-stage keep counts must be nonincreasing")
    if min(
        args.exact_score_timeout,
        *args.stage_timeouts,
        args.height_timeout,
        args.saturation_timeout,
        args.conductor_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if args.stack_bytes < 64_000_000 or args.certificate_prime_bound < 3:
        raise SystemExit("invalid stack or certificate-prime bound")

    tables = build_residue_tables()
    prefiltered, population_audit = prefilter_population(
        a_max=args.a_max,
        b_max=args.b_max,
        keep_count=args.residue_keep,
        tables=tables,
    )
    exactly_scored = exact_score_candidates(
        prefiltered,
        cutoff=EXACT_SCORE_CUTOFF,
        batch_size=args.exact_score_batch,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )
    exact_survivors = exactly_scored[: args.exact_keep]

    # Calibrations pass through the same scorers but never consume a survivor
    # position or influence a selection boundary.
    calibration_prefilter = tuple(
        PrefilterCandidate(
            parameter.numerator,
            parameter.denominator,
            *residue_score(parameter.numerator, parameter.denominator, tables),
        )
        for parameter in CALIBRATION_PARAMETERS
    )
    calibration_exact = exact_score_candidates(
        calibration_prefilter,
        cutoff=EXACT_SCORE_CUTOFF,
        batch_size=len(calibration_prefilter),
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )

    stage_artifacts: list[dict[str, Any]] = []
    current = tuple(exact_survivors)
    final_pools: tuple[PointPool, ...] = ()
    final_ranks: dict[str, dict[str, Any]] = {}
    threshold_pools: dict[str, tuple[PointPool, dict[str, Any]]] = {}
    for stage_index, (height_bound, search_timeout, workers) in enumerate(
        zip(args.stage_heights, args.stage_timeouts, args.stage_workers), start=1
    ):
        pools = parallel_point_search(
            current,
            height_bound=height_bound,
            timeout=float(search_timeout),
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ranks = parallel_rank_replay(
            pools,
            precisions=args.height_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ordered = tuple(
            sorted(
                pools,
                key=lambda point_pool: pool_priority(
                    point_pool, ranks[point_pool.candidate.identifier]
                ),
            )
        )
        for point_pool in ordered:
            record = ranks[point_pool.candidate.identifier]
            if (
                record.get("status") == "completed"
                and int(record["stable_numerical_rank"]) >= 18
            ):
                threshold_pools[point_pool.candidate.identifier] = (
                    point_pool,
                    record,
                )
        keep_count = (
            args.stage_keeps[stage_index - 1]
            if stage_index <= len(args.stage_keeps)
            else len(ordered)
        )
        retained = ordered[: min(keep_count, len(ordered))]
        stage_artifacts.append(
            {
                "stage": stage_index,
                "quartic_naive_height_bound": height_bound,
                "population_searched": len(pools),
                "completed_point_searches": sum(pool.status == "completed" for pool in pools),
                "point_search_timeouts": sum(pool.status == "timeout" for pool in pools),
                "point_search_errors": sum(pool.status == "error" for pool in pools),
                "completed_height_replays": sum(record.get("status") == "completed" for record in ranks.values()),
                "height_replay_timeouts": sum(record.get("status") == "timeout" for record in ranks.values()),
                "height_replay_errors": sum(record.get("status") == "error" for record in ranks.values()),
                "ranked_population": [
                    point_pool_record(
                        point_pool,
                        ranks[point_pool.candidate.identifier],
                        include_points=(stage_index == 3),
                    )
                    for point_pool in ordered
                ],
                "retained_candidate_ids": [pool.candidate.identifier for pool in retained],
            }
        )
        current = tuple(pool.candidate for pool in retained)
        final_pools, final_ranks = ordered, ranks

    # Run all calibration parameters at all point heights, separately.
    calibration_stages = []
    calibration_current = tuple(calibration_exact)
    for height_bound, search_timeout, workers in zip(
        args.stage_heights, args.stage_timeouts, args.stage_workers
    ):
        pools = parallel_point_search(
            calibration_current,
            height_bound=height_bound,
            timeout=float(search_timeout),
            stack_bytes=args.stack_bytes,
            workers=min(workers, len(calibration_current)),
        )
        ranks = parallel_rank_replay(
            pools,
            precisions=args.height_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            workers=min(workers, len(calibration_current)),
        )
        calibration_stages.append(
            {
                "quartic_naive_height_bound": height_bound,
                "records": [
                    point_pool_record(
                        pool, ranks[pool.candidate.identifier], include_points=False
                    )
                    for pool in pools
                ],
            }
        )

    certificates: dict[str, dict[str, Any]] = {}
    conductors: dict[str, dict[str, Any]] = {}
    certified_hits = []
    certificate_targets = dict(threshold_pools)
    # Final leaders also get conductors even below rank 18, so the frontier is
    # immediately interpretable against the target.
    final_leaders = final_pools[: min(16, len(final_pools))]
    replay_targets = {
        point_pool.candidate.identifier: (
            point_pool,
            final_ranks[point_pool.candidate.identifier],
        )
        for point_pool in final_leaders
    }
    # A rank>=18 pool is replayed even if an unusually crowded threshold set
    # later pushes it outside a stage keep boundary.
    replay_targets.update(certificate_targets)
    for identifier in sorted(replay_targets):
        point_pool, rank_record = replay_targets[identifier]
        try:
            conductors[identifier] = {
                "status": "completed",
                **minimal_curve_data(
                    point_pool.coefficients,
                    timeout=args.conductor_timeout,
                    stack_bytes=args.stack_bytes,
                ),
            }
            conductors[identifier]["below_strict_log_conductor_target"] = (
                Decimal(conductors[identifier]["log_conductor"])
                < TARGET_LOG_CONDUCTOR
            )
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
            conductors[identifier] = {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:500],
            }
        if identifier in certificate_targets:
            try:
                certificate = finite_reduction_certificate(
                    point_pool,
                    rank_record,
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                    stack_bytes=args.stack_bytes,
                )
            except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
                certificate = {
                    "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                    "error": str(error)[:500],
                }
            certificates[identifier] = certificate
            if (
                certificate.get("status") == "certified"
                and (
                    int(certificate["certified_algebraic_rank_lower_bound"]) >= 30
                    or (
                        int(certificate["certified_algebraic_rank_lower_bound"]) >= 21
                        and conductors[identifier].get("below_strict_log_conductor_target") is True
                    )
                )
            ):
                certified_hits.append(
                    {
                        "candidate_id": identifier,
                        "constructor_parameter": rational_to_string(point_pool.candidate.parameter),
                        "certified_rank_lower_bound": certificate["certified_algebraic_rank_lower_bound"],
                        "conductor": conductors[identifier].get("conductor"),
                        "log_conductor": conductors[identifier].get("log_conductor"),
                    }
                )

    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded unbiased score-and-point search; numerical heights are "
            "triage only, and target hits require exact finite-reduction and "
            "conductor certificates"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": certified_hits,
        },
        "population": {
            "numerator_range": [1, args.a_max],
            "denominator_range": [1, args.b_max],
            "primitive_positive_pairs_only": True,
            "no_crt_or_conductor_filter": True,
            "calibration_parameters_excluded_before_selection": [
                rational_to_string(parameter) for parameter in CALIBRATION_PARAMETERS
            ],
            **population_audit,
        },
        "residue_prefilter": {
            "cutoff": RESIDUE_CUTOFF,
            "projective_parameter_line_includes_infinity": True,
            "prime_count": len(tables),
            "lookup_tables": {
                str(prime): [
                    {
                        "projective_symbol": symbol.label,
                        "good_reduction": symbol.good_reduction,
                        "ellap": symbol.ellap,
                        "contribution": symbol.contribution,
                    }
                    for symbol in table
                ]
                for prime, table in tables.items()
            },
            "retained_count": len(prefiltered),
            "ranked_tail": [
                {
                    "constructor_parameter": rational_to_string(candidate.parameter),
                    "score": candidate.residue_score_b200,
                    "good_primes_used": candidate.residue_good_primes,
                    "bad_primes_skipped": candidate.residue_bad_primes,
                }
                for candidate in prefiltered
            ],
        },
        "exact_b2000_rescore": {
            "cutoff": EXACT_SCORE_CUTOFF,
            "population_rescored": len(exactly_scored),
            "retained_for_points": len(exact_survivors),
            "ranked_population": [exact_score_record(candidate) for candidate in exactly_scored],
            "retained_candidate_ids": [candidate.identifier for candidate in exact_survivors],
        },
        "point_stages": stage_artifacts,
        "calibration": {
            "selection_role": "excluded before every novel-candidate selection boundary",
            "exact_scores": [exact_score_record(candidate) for candidate in calibration_exact],
            "point_stages": calibration_stages,
        },
        "finite_reduction_certificates": certificates,
        "final_conductor_replays": conductors,
        "summary": {
            "complete_primitive_population": population_audit["complete_primitive_population_count"],
            "eligible_novel_population": population_audit["eligible_population_count"],
            "residue_prefilter_survivors": len(prefiltered),
            "exact_b2000_point_survivors": len(exact_survivors),
            "deep_search_population": len(final_pools),
            "maximum_stable_numerical_rank": max(
                (
                    int(record["stable_numerical_rank"])
                    for record in final_ranks.values()
                    if record.get("status") == "completed"
                ),
                default=0,
            ),
            "stable_rank_at_least_18_candidate_ids": sorted(threshold_pools),
            "certified_target_hit": bool(certified_hits),
        },
        "bounds": {
            "residue_keep_count": args.residue_keep,
            "exact_b2000_keep_count": args.exact_keep,
            "exact_score_batch_size": args.exact_score_batch,
            "exact_score_timeout_seconds_per_batch": args.exact_score_timeout,
            "quartic_naive_height_stages": list(args.stage_heights),
            "point_stage_keep_counts": list(args.stage_keeps),
            "search_timeout_seconds_per_candidate": list(args.stage_timeouts),
            "worker_counts": list(args.stage_workers),
            "height_precisions": list(args.height_precisions),
            "height_timeout_seconds_per_candidate": args.height_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "pari_stack_bytes_per_process": args.stack_bytes,
            "process_policy": (
                "bounded thread pools launch only synchronous subprocess.run "
                "calls; exiting each pool joins every worker and child"
            ),
        },
        "interpretation": (
            "The score stages select a rare-event tail but do not estimate or "
            "bound rank.  H-bounded point searches are not exhaustive over Q, "
            "and numerical height rank is not an independence proof."
        ),
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "invocation": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True, default=str) + "\n")
    print(f"wrote {args.output}")
    print(json.dumps(artifact["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
