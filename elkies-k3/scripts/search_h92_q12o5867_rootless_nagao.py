#!/usr/bin/env python3
"""Projective rational-parameter Nagao sieve for the rootless R17 family.

By default the exact short model is the compact coordinate published by
Elkies in arXiv:2608.25406v1.  The independently reconstructed q12/orbit5867
model remains accepted through ``--model`` as a regression chart.  The two
models are exactly related by the repository's pinned Mobius and Weierstrass
scaling certificate.  Searching in the published coordinate is materially
better: the four published rank-25--28 fibres have projective height at most
9529 there but roughly 200-bit coordinates in the raw q12 chart.

For every usable prime and every point of P^1(F_p), this script records the
exact point count, Frobenius trace, the repository's Nagao contribution, and
whether the reduced cubic is singular. Rational parameters are then scored
using table lookups only; the endpoint sections are deliberately not loaded.

This is a heuristic search tool.  A high score or survival through a stage is
not a rank statement.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, log
from pathlib import Path
import shlex
import sys
from time import perf_counter
from typing import Iterable, Iterator, Sequence


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL = (
    ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
)
DEFAULT_PRIME_BLOCKS = (
    (19, 41, 43, 61, 71, 73, 79, 83),
    (89, 107, 113, 127, 131, 137, 139, 151),
    (157, 163, 167, 173, 179, 181, 191, 193, 197),
)
SCORE_SCALE = 10**12


def _allow_large_exact_coefficients() -> None:
    # The pinned endpoint coefficients contain decimal integers beyond Python's
    # default 4,300-digit parsing guard.  They are trusted local exact data, not
    # user-supplied expressions.
    setter = getattr(sys, "set_int_max_str_digits", None)
    if setter is not None:
        setter(100_000)


@dataclass(frozen=True)
class FamilyModel:
    source: Path
    source_sha256: str
    a_coefficients: tuple[Fraction, ...]
    b_coefficients: tuple[Fraction, ...]
    a_degree: int
    b_degree: int
    coordinate: str
    coefficient_source_keys: tuple[str, str]


@dataclass(frozen=True)
class LocalSymbol:
    prime: int
    projective_index: int
    point_count: int | None
    trace: int | None
    contribution_units: int
    good_reduction: bool
    singular_mod_prime: bool

    @property
    def label(self) -> str:
        return "infinity" if self.projective_index == self.prime else str(self.projective_index)


@dataclass(frozen=True)
class Candidate:
    numerator: int
    denominator: int
    height: int
    block_score_units: tuple[int, ...] = ()
    good_primes: int = 0
    bad_primes: int = 0

    @property
    def parameter(self) -> str:
        if self.denominator == 0:
            return "infinity"
        return f"{self.numerator}/{self.denominator}"

    @property
    def total_score_units(self) -> int:
        return sum(self.block_score_units)


def load_family_model(path: Path = DEFAULT_MODEL) -> FamilyModel:
    """Load A and B from the published chart or exact q12 replay artifact."""

    _allow_large_exact_coefficients()
    payload = path.read_bytes()
    document = json.loads(payload)
    status = document.get("status")
    if status == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
        child = document
        coordinate = document["coordinate"]
        coefficient_keys = (
            "A_coefficients_low_to_high",
            "B_coefficients_low_to_high",
        )
    elif status == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN":
        child = document["child"]
        coordinate = "q12o5867_raw_u"
        coefficient_keys = (
            "minimal_A_coefficients_low_to_high",
            "minimal_B_coefficients_low_to_high",
        )
    else:
        raise ValueError("the input is not a certified rootless R17 model")
    degrees = child["degrees_A_B_Delta"]
    if degrees[:2] != [8, 12]:
        raise ValueError(f"expected short-model degrees (8,12), found {degrees[:2]}")
    a_coefficients = tuple(Fraction(value) for value in child[coefficient_keys[0]])
    b_coefficients = tuple(Fraction(value) for value in child[coefficient_keys[1]])
    if len(a_coefficients) != 9 or len(b_coefficients) != 13:
        raise ValueError("the exact coefficient arrays have unexpected lengths")
    return FamilyModel(
        source=path.resolve(),
        source_sha256=sha256(payload).hexdigest(),
        a_coefficients=a_coefficients,
        b_coefficients=b_coefficients,
        a_degree=8,
        b_degree=12,
        coordinate=coordinate,
        coefficient_source_keys=coefficient_keys,
    )


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def coefficient_mod(coefficient: Fraction, prime: int) -> int:
    denominator = coefficient.denominator % prime
    if denominator == 0:
        raise ZeroDivisionError(f"coefficient denominator is divisible by p={prime}")
    return coefficient.numerator % prime * pow(denominator, -1, prime) % prime


def reduced_coefficients(model: FamilyModel, prime: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if prime < 5 or not is_prime(prime):
        raise ValueError("local tables require a prime at least five")
    return (
        tuple(coefficient_mod(value, prime) for value in model.a_coefficients),
        tuple(coefficient_mod(value, prime) for value in model.b_coefficients),
    )


def homogeneous_value_mod(
    coefficients: Sequence[int], degree: int, numerator: int, denominator: int, prime: int
) -> int:
    """Evaluate the degree-``degree`` homogenization at ``(numerator:denominator)``."""

    if len(coefficients) > degree + 1:
        raise ValueError("coefficient list exceeds the declared homogeneous degree")
    return sum(
        coefficient
        * pow(numerator, power, prime)
        * pow(denominator, degree - power, prime)
        for power, coefficient in enumerate(coefficients)
    ) % prime


def quadratic_character_table(prime: int) -> tuple[int, ...]:
    characters = [-1] * prime
    characters[0] = 0
    for value in range(1, prime):
        characters[value * value % prime] = 1
    return tuple(characters)


def local_symbol(
    prime: int,
    projective_index: int,
    a_coefficients: Sequence[int],
    b_coefficients: Sequence[int],
    characters: Sequence[int] | None = None,
) -> LocalSymbol:
    if not 0 <= projective_index <= prime:
        raise ValueError("projective index must lie in [0,p]")
    numerator, denominator = (
        (projective_index, 1) if projective_index < prime else (1, 0)
    )
    coefficient_a = homogeneous_value_mod(a_coefficients, 8, numerator, denominator, prime)
    coefficient_b = homogeneous_value_mod(b_coefficients, 12, numerator, denominator, prime)
    discriminant_core = (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime
    if discriminant_core == 0:
        return LocalSymbol(prime, projective_index, None, None, 0, False, True)
    if characters is None:
        characters = quadratic_character_table(prime)
    trace = -sum(
        characters[(x_value**3 + coefficient_a * x_value + coefficient_b) % prime]
        for x_value in range(prime)
    )
    point_count = prime + 1 - trace
    contribution = (2.0 - trace) / point_count * log(float(prime))
    return LocalSymbol(
        prime,
        projective_index,
        point_count,
        trace,
        int(round(contribution * SCORE_SCALE)),
        True,
        False,
    )


def residue_table(model: FamilyModel, prime: int) -> tuple[LocalSymbol, ...]:
    """Return all p+1 local symbols, with infinity stored at index p."""

    a_coefficients, b_coefficients = reduced_coefficients(model, prime)
    characters = quadratic_character_table(prime)
    table = tuple(
        local_symbol(prime, index, a_coefficients, b_coefficients, characters)
        for index in range(prime + 1)
    )
    if len(table) != prime + 1 or table[-1].projective_index != prime:
        raise AssertionError("incomplete projective local table")
    return table


def build_residue_tables(
    model: FamilyModel, prime_blocks: Sequence[Sequence[int]]
) -> tuple[tuple[dict[int, tuple[LocalSymbol, ...]], ...], tuple[dict[str, object], ...]]:
    """Build usable tables and report primes where this rational model is not integral."""

    blocks: list[dict[int, tuple[LocalSymbol, ...]]] = []
    rejected: list[dict[str, object]] = []
    seen: set[int] = set()
    for block_number, primes in enumerate(prime_blocks, start=1):
        block: dict[int, tuple[LocalSymbol, ...]] = {}
        for prime in primes:
            if prime in seen:
                raise ValueError(f"prime p={prime} occurs in more than one block")
            seen.add(prime)
            if prime < 5 or not is_prime(prime):
                raise ValueError(f"p={prime} is not a prime at least five")
            try:
                table = residue_table(model, prime)
            except ZeroDivisionError:
                rejected.append(
                    {
                        "prime": prime,
                        "block": block_number,
                        "reason": "coefficient_denominator_not_invertible",
                    }
                )
                continue
            if not any(symbol.good_reduction for symbol in table):
                rejected.append(
                    {
                        "prime": prime,
                        "block": block_number,
                        "reason": "short_model_degenerate_at_every_projective_parameter",
                    }
                )
                continue
            block[prime] = table
        if not block:
            raise ValueError(f"prime block {block_number} has no usable primes")
        blocks.append(block)
    return tuple(blocks), tuple(rejected)


def projective_index(numerator: int, denominator: int, prime: int) -> int:
    if denominator % prime == 0:
        if numerator % prime == 0:
            raise ValueError("a nonprimitive pair has no projective reduction")
        return prime
    return numerator % prime * pow(denominator % prime, -1, prime) % prime


def primitive_parameters(numerator_bound: int, denominator_bound: int) -> Iterator[Candidate]:
    """Enumerate canonical primitive ``(a:b)``, including zero and infinity."""

    if numerator_bound < 0 or denominator_bound < 1:
        raise ValueError("require numerator_bound >= 0 and denominator_bound >= 1")
    yield Candidate(1, 0, 1)  # infinity
    for denominator in range(1, denominator_bound + 1):
        for numerator in range(-numerator_bound, numerator_bound + 1):
            if gcd(abs(numerator), denominator) == 1:
                yield Candidate(numerator, denominator, max(abs(numerator), denominator))


def score_block(
    candidate: Candidate,
    tables: dict[int, tuple[LocalSymbol, ...]],
    inverse_cache: dict[tuple[int, int], int | None],
) -> Candidate:
    score_units = 0
    good = 0
    bad = 0
    for prime, table in tables.items():
        if candidate.denominator == 0 or candidate.denominator % prime == 0:
            index = prime
        else:
            cache_key = (prime, candidate.denominator % prime)
            inverse = inverse_cache.get(cache_key)
            if inverse is None:
                inverse = pow(cache_key[1], -1, prime)
                inverse_cache[cache_key] = inverse
            index = candidate.numerator % prime * inverse % prime
        symbol = table[index]
        if symbol.good_reduction:
            score_units += symbol.contribution_units
            good += 1
        else:
            bad += 1
    return replace(
        candidate,
        block_score_units=candidate.block_score_units + (score_units,),
        good_primes=candidate.good_primes + good,
        bad_primes=candidate.bad_primes + bad,
    )


def dominates(left: Candidate, right: Candidate) -> bool:
    if len(left.block_score_units) != len(right.block_score_units):
        raise ValueError("Pareto comparison requires the same completed stages")
    weak = (
        left.total_score_units >= right.total_score_units
        and left.good_primes >= right.good_primes
        and left.bad_primes <= right.bad_primes
        and left.height <= right.height
    )
    strict = (
        left.total_score_units != right.total_score_units
        or left.good_primes != right.good_primes
        or left.bad_primes != right.bad_primes
        or left.height != right.height
    )
    return weak and strict


def candidate_sort_key(candidate: Candidate) -> tuple[object, ...]:
    return (
        -candidate.total_score_units,
        -min(candidate.block_score_units),
        -candidate.good_primes,
        candidate.bad_primes,
        candidate.height,
        candidate.denominator,
        candidate.numerator,
    )


def retain_bucketed_pareto(
    candidates: Iterable[Candidate], *, bucket_width: int, cap_per_bucket: int
) -> tuple[list[Candidate], dict[str, int]]:
    """Retain a deterministic capped Pareto frontier in each height bucket."""

    if bucket_width < 1 or cap_per_bucket < 1:
        raise ValueError("bucket width and cap must be positive")
    buckets: dict[int, list[Candidate]] = {}
    population = 0
    pareto_insertions = 0
    for candidate in candidates:
        population += 1
        bucket = (candidate.height - 1) // bucket_width
        frontier = buckets.setdefault(bucket, [])
        if any(dominates(retained, candidate) for retained in frontier):
            continue
        frontier[:] = [retained for retained in frontier if not dominates(candidate, retained)]
        frontier.append(candidate)
        pareto_insertions += 1
        if len(frontier) > cap_per_bucket:
            frontier.sort(key=candidate_sort_key)
            del frontier[cap_per_bucket:]
    retained = sorted(
        (candidate for frontier in buckets.values() for candidate in frontier),
        key=candidate_sort_key,
    )
    return retained, {
        "population_scored": population,
        "height_bucket_count": len(buckets),
        "pareto_insertions": pareto_insertions,
        "retained_count": len(retained),
    }


def run_staged_sieve(
    *,
    numerator_bound: int,
    denominator_bound: int,
    table_blocks: Sequence[dict[int, tuple[LocalSymbol, ...]]],
    keep_per_bucket: Sequence[int],
    bucket_width: int,
) -> tuple[list[Candidate], list[dict[str, object]]]:
    if len(table_blocks) != len(keep_per_bucket):
        raise ValueError("one keep-per-bucket value is required for each prime block")
    survivors: Iterable[Candidate] = primitive_parameters(numerator_bound, denominator_bound)
    stages: list[dict[str, object]] = []
    for stage_number, (tables, cap) in enumerate(zip(table_blocks, keep_per_bucket), start=1):
        inverse_cache: dict[tuple[int, int], int | None] = {}
        started = perf_counter()
        scored = (score_block(candidate, tables, inverse_cache) for candidate in survivors)
        retained, summary = retain_bucketed_pareto(
            scored, bucket_width=bucket_width, cap_per_bucket=cap
        )
        elapsed = perf_counter() - started
        stages.append(
            {
                "stage": stage_number,
                "primes": list(tables),
                "cap_per_height_bucket": cap,
                **summary,
                "runtime_seconds": elapsed,
                "parameters_per_second": summary["population_scored"] / elapsed if elapsed else None,
            }
        )
        survivors = retained
    return list(survivors), stages


def parse_prime_blocks(text: str) -> tuple[tuple[int, ...], ...]:
    """Parse semicolon-separated blocks of comma-separated primes or ranges."""

    blocks: list[tuple[int, ...]] = []
    for raw_block in text.split(";"):
        primes: list[int] = []
        for token in raw_block.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                start_text, end_text = token.split("-", 1)
                start, end = int(start_text), int(end_text)
                primes.extend(value for value in range(start, end + 1) if is_prime(value))
            else:
                primes.append(int(token))
        if not primes:
            raise ValueError("empty prime block")
        blocks.append(tuple(primes))
    if not blocks:
        raise ValueError("no prime blocks supplied")
    return tuple(blocks)


def local_symbol_record(symbol: LocalSymbol) -> dict[str, object]:
    return {
        "projective_index": symbol.projective_index,
        "parameter_mod_p": symbol.label,
        "good_reduction": symbol.good_reduction,
        "bad_reduction": not symbol.good_reduction,
        "singular_mod_prime": symbol.singular_mod_prime,
        "point_count": symbol.point_count,
        "a_p": symbol.trace,
        "nagao_contribution_units_1e12": symbol.contribution_units,
    }


def candidate_record(candidate: Candidate) -> dict[str, object]:
    return {
        "parameter": candidate.parameter,
        "projective_pair": [candidate.numerator, candidate.denominator],
        "projective_height": candidate.height,
        "block_score_units_1e12": list(candidate.block_score_units),
        "total_score_units_1e12": candidate.total_score_units,
        "total_score": candidate.total_score_units / SCORE_SCALE,
        "good_prime_count": candidate.good_primes,
        "bad_reduction_prime_count": candidate.bad_primes,
    }


def default_prime_blocks_text() -> str:
    return ";".join(",".join(map(str, block)) for block in DEFAULT_PRIME_BLOCKS)


def export_cpp_tables(
    path: Path,
    model: FamilyModel,
    table_blocks: Sequence[dict[int, tuple[LocalSymbol, ...]]],
) -> None:
    """Write a compact, complete table file for the C++ lookup scanner."""

    lines = [
        "H92_Q12O5867_PROJECTIVE_NAGAO_TABLE_V1",
        f"M {model.source_sha256} {model.a_degree} {model.b_degree} {SCORE_SCALE}",
        f"C {len(table_blocks)}",
    ]
    for block_number, block in enumerate(table_blocks, start=1):
        lines.append(f"B {block_number} {len(block)}")
        for prime, table in block.items():
            if len(table) != prime + 1:
                raise AssertionError("cannot export an incomplete projective table")
            lines.append(f"P {prime} {len(table)}")
            for symbol in table:
                point_count = symbol.point_count if symbol.point_count is not None else -1
                trace = symbol.trace if symbol.trace is not None else 0
                lines.append(
                    f"{int(symbol.good_reduction)} {int(symbol.singular_mod_prime)} "
                    f"{point_count} {trace} {symbol.contribution_units}"
                )
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Projective Nagao sieve for the exact rootless R17 family."
    )
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--export-cpp-tables",
        type=Path,
        help="write complete projective local tables for the C++ hot loop",
    )
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="stop after --export-cpp-tables (no rational-parameter scan)",
    )
    parser.add_argument("--numerator-bound", type=int, default=1000)
    parser.add_argument("--denominator-bound", type=int, default=1000)
    parser.add_argument("--prime-blocks", default=default_prime_blocks_text())
    parser.add_argument(
        "--keep-per-bucket",
        default="32,16,8",
        help="comma-separated capped Pareto sizes, one per prime block",
    )
    parser.add_argument("--height-bucket-width", type=int, default=100)
    parser.add_argument("--finalists", type=int, default=1000)
    args = parser.parse_args()

    if args.tables_only and args.export_cpp_tables is None:
        raise SystemExit("--tables-only requires --export-cpp-tables")
    if not args.tables_only and args.output is None:
        raise SystemExit("--output is required unless --tables-only is used")

    started = perf_counter()
    prime_blocks = parse_prime_blocks(args.prime_blocks)
    keep_per_bucket = tuple(int(value) for value in args.keep_per_bucket.split(","))
    if len(keep_per_bucket) != len(prime_blocks):
        raise SystemExit("--keep-per-bucket must have one entry per prime block")
    if args.finalists < 1:
        raise SystemExit("--finalists must be positive")

    model = load_family_model(args.model)
    tables_started = perf_counter()
    table_blocks, rejected_primes = build_residue_tables(model, prime_blocks)
    table_seconds = perf_counter() - tables_started
    if args.export_cpp_tables is not None:
        export_cpp_tables(args.export_cpp_tables, model, table_blocks)
    if args.tables_only:
        print(
            f"PASS tables={sum(len(block) for block in table_blocks)} "
            f"rejected={len(rejected_primes)} output={args.export_cpp_tables}"
        )
        return
    survivors, stages = run_staged_sieve(
        numerator_bound=args.numerator_bound,
        denominator_bound=args.denominator_bound,
        table_blocks=table_blocks,
        keep_per_bucket=keep_per_bucket,
        bucket_width=args.height_bucket_width,
    )
    finalists = sorted(survivors, key=candidate_sort_key)[: args.finalists]
    total_seconds = perf_counter() - started

    document = {
        "schema": "h92-q12o5867-rootless-projective-nagao-sieve-v1",
        "status": "PASS_BOUNDED_HEURISTIC_PROJECTIVE_NAGAO_SIEVE",
        "proof_boundary": (
            "Nagao scores and staged survival are search heuristics, not rank bounds. "
            "No endpoint section was evaluated and no candidate was promoted."
        ),
        "model": {
            "source": str(model.source),
            "source_sha256": model.source_sha256,
            "coordinate": model.coordinate,
            "degrees_A_B": [model.a_degree, model.b_degree],
            "coefficient_source_keys": list(model.coefficient_source_keys),
        },
        "nagao_contribution": {
            "formula": "((2-a_p)/(p+1-a_p))*log(p)",
            "integer_scale": SCORE_SCALE,
            "bad_reduction_contribution": 0,
        },
        "search": {
            "numerator_interval": [-args.numerator_bound, args.numerator_bound],
            "denominator_interval": [1, args.denominator_bound],
            "includes_zero": True,
            "includes_infinity": True,
            "primitive_pairs_only": True,
            "height": "max(abs(a),b)",
            "height_bucket_width": args.height_bucket_width,
            "requested_prime_blocks": [list(block) for block in prime_blocks],
            "usable_prime_blocks": [list(block) for block in table_blocks],
            "rejected_primes": list(rejected_primes),
            "keep_per_bucket": list(keep_per_bucket),
        },
        "local_tables": {
            str(prime): [local_symbol_record(symbol) for symbol in table]
            for block in table_blocks
            for prime, table in block.items()
        },
        "stages": stages,
        "final_survivor_count": len(survivors),
        "finalists": [candidate_record(candidate) for candidate in finalists],
        "runtime": {
            "table_seconds": table_seconds,
            "total_seconds": total_seconds,
        },
        "reproducing_command": shlex.join(sys.argv),
    }
    assert args.output is not None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS tables={sum(len(block) for block in table_blocks)} "
        f"survivors={len(survivors)} finalists={len(finalists)} "
        f"seconds={total_seconds:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
