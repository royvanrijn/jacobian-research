#!/usr/bin/env python3
"""Nagao-style census of the bisection quadratic twists.

For a good finite fibre and a nonzero twist value, the trace identity is

    a_p(E^q_t) = chi_p(q(t)) * a_p(E_t).

The expensive point counts on the published rootless R17 family are therefore
computed once per prime.  This script then scores all 39,120 individual
bisection twists and the 5,566 product twists belonging to the immediate-point
paired bases by quadratic-character dot products.

The output is a bounded heuristic ranking, not a Mordell--Weil rank statement.
In particular, a high score neither constructs a second section nor proves a
nonzero product-twist rank.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, log, lcm
from pathlib import Path
import shlex
import sys
from time import perf_counter
from typing import Sequence

from search_h92_q12o5867_rootless_nagao import (
    DEFAULT_MODEL,
    load_family_model,
    quadratic_character_table,
    reduced_coefficients,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BISECTIONS = (
    ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
)
DEFAULT_PAIRS = (
    ROOT
    / "artifacts/generated-results/elkies-2026-immediate-point-pair-catalogue-full.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-2026-quadratic-twist-rank-census.json"
)
DEFAULT_PRIME_BLOCKS = (
    (211, 223, 227, 229, 233, 239, 241, 251),
    (257, 263, 269, 271, 277, 281, 283, 293),
    (307, 311, 313, 317, 331, 337, 347, 349),
    (353, 359, 367, 373, 379, 383, 389, 397),
    (401, 409, 419, 421, 431, 433, 439, 443),
    (449, 457, 461, 463, 467, 479, 487, 491),
)
EXPECTED_BISECTION_SCHEMA = "elkies-k3.bisection-extension-input.v1"
EXPECTED_PAIR_SCHEMA = "elkies-k3.elkies-2026-immediate-point-pair-catalogue.v1"
OUTPUT_SCHEMA = "elkies-k3.elkies-2026-quadratic-twist-rank-census.v1"
RANK_NINE_PAIR_KEYS = {"42110:43109", "71804:81769"}


@dataclass(frozen=True)
class Candidate:
    kind: str
    key: str
    masks: tuple[int, ...]
    coefficients: tuple[int, ...]
    forced_twist_rank: int
    metadata: dict[str, object]


@dataclass(frozen=True)
class PrimeData:
    prime: int
    traces: tuple[int, ...]
    characters: tuple[int, ...]
    singular_fibre_count: int


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def square_equivalent_integer_polynomial(values: Sequence[str]) -> tuple[int, ...]:
    """Return integer coefficients after multiplying by a rational square.

    If D clears all denominators, D^2*q is integral and represents the same
    class in QQ(t)^*/QQ(t)^*2.  We intentionally do not primitive-normalize:
    a nonsquare rational content is part of the quadratic twist.
    """

    coefficients = tuple(Fraction(value) for value in values)
    denominator = lcm(*(value.denominator for value in coefficients))
    integers = tuple(int(value * denominator * denominator) for value in coefficients)
    if not integers or all(value == 0 for value in integers):
        raise ValueError("zero twist polynomial")
    return trim_integer_polynomial(integers)


def trim_integer_polynomial(values: Sequence[int]) -> tuple[int, ...]:
    coefficients = list(map(int, values))
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def multiply_integer_polynomials(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    product = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product[left_index + right_index] += left_value * right_value
    return trim_integer_polynomial(product)


def valuation(value: int, prime: int) -> int:
    if value == 0:
        return 10**9
    value = abs(value)
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def trim_mod(values: Sequence[int]) -> list[int]:
    result = list(values)
    while result and result[-1] == 0:
        result.pop()
    return result


def polynomial_remainder_mod(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    remainder = trim_mod([value % prime for value in left])
    divisor = trim_mod([value % prime for value in right])
    if not divisor:
        raise ZeroDivisionError("zero polynomial divisor")
    inverse = pow(divisor[-1], -1, prime)
    while len(remainder) >= len(divisor):
        scale = remainder[-1] * inverse % prime
        offset = len(remainder) - len(divisor)
        for index, value in enumerate(divisor):
            remainder[offset + index] = (remainder[offset + index] - scale * value) % prime
        remainder = trim_mod(remainder)
    return remainder


def polynomial_squarefree_mod(values: Sequence[int], declared_degree: int, prime: int) -> bool:
    polynomial = trim_mod([value % prime for value in values])
    if not polynomial:
        return False
    # The homogeneous form has an infinity root of multiplicity
    # declared_degree-degree(polynomial).  More than one is not squarefree.
    if len(polynomial) - 1 < declared_degree - 1:
        return False
    derivative = [(index * polynomial[index]) % prime for index in range(1, len(polynomial))]
    left = polynomial
    right = trim_mod(derivative)
    while right:
        left, right = right, polynomial_remainder_mod(left, right, prime)
    return len(left) == 1


def reduce_twist_polynomial(
    coefficients: Sequence[int], prime: int
) -> tuple[int, ...] | None:
    """Reduce a squareclass at p, or return None at a bad twist prime."""

    minimum = min(valuation(value, prime) for value in coefficients)
    if minimum >= 10**9 or minimum % 2:
        return None
    divisor = prime**minimum
    reduced = tuple((value // divisor) % prime for value in coefficients)
    if not polynomial_squarefree_mod(reduced, len(coefficients) - 1, prime):
        return None
    return tuple(reduced)


def evaluate_mod(coefficients: Sequence[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def build_prime_data(model, prime: int) -> PrimeData:
    sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
    from research_runtime.finite_fields import family_traces
    record = family_traces(model.a_coefficients, model.b_coefficients, prime,
                           a_degree=model.a_degree, b_degree=model.b_degree)
    affine = record["fibres"][:-1]
    return PrimeData(prime=prime,
        traces=tuple(0 if row["singular"] else row["trace"] for row in affine),
        characters=tuple(record["quadratic_characters"]),
        singular_fibre_count=sum(row["singular"] for row in affine))


def canonical_mod_square(values: Sequence[int], prime: int) -> tuple[int, tuple[int, ...]]:
    pivot = next(value for value in values if value % prime)
    multiplier = 1 if pow(pivot, (prime - 1) // 2, prime) == 1 else -1
    inverse = pow(pivot, -1, prime)
    normalized = tuple(value * inverse % prime for value in values)
    return multiplier, normalized


def local_twist_average(
    candidate: Candidate,
    data: PrimeData,
    cache: dict[tuple[int, ...], int],
) -> Fraction | None:
    reduced = reduce_twist_polynomial(candidate.coefficients, data.prime)
    if reduced is None:
        return None
    multiplier, normalized = canonical_mod_square(reduced, data.prime)
    trace_sum = cache.get(normalized)
    if trace_sum is None:
        trace_sum = sum(
            data.characters[evaluate_mod(normalized, parameter, data.prime)]
            * data.traces[parameter]
            for parameter in range(data.prime)
        )
        cache[normalized] = trace_sum
    # A_p=(1/p) sum_t a_p(E^q_t); Nagao uses -A_p.
    return Fraction(-multiplier * trace_sum, data.prime)


def parse_prime_block(text: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in text.split(",") if value.strip())
    if not values:
        raise argparse.ArgumentTypeError("prime block is empty")
    return values


def load_candidates(
    bisection_path: Path, pair_path: Path
) -> tuple[list[Candidate], list[Candidate], dict[str, object]]:
    bisections = json.loads(bisection_path.read_text())
    if bisections.get("schema") != EXPECTED_BISECTION_SCHEMA:
        raise ValueError("unexpected bisection batch schema")
    by_mask: dict[int, Candidate] = {}
    singleton_candidates = []
    for record in bisections["bisections"]:
        mask = int(record["lattice_orbit_mask"])
        candidate = Candidate(
            kind="singleton",
            key=str(mask),
            masks=(mask,),
            coefficients=square_equivalent_integer_polynomial(
                record["residual_chord"]["q_coefficients"]
            ),
            forced_twist_rank=1,
            metadata={
                "orbit_hex": f"0x{mask:05x}",
                "priority_rank": int(record["priority_rank"]),
                "equation_rank": int(record["equation_rank"]),
            },
        )
        if mask in by_mask:
            raise ValueError(f"duplicate bisection mask {mask}")
        by_mask[mask] = candidate
        singleton_candidates.append(candidate)

    pairs = json.loads(pair_path.read_text())
    if pairs.get("schema") != EXPECTED_PAIR_SCHEMA:
        raise ValueError("unexpected immediate-point pair catalogue schema")
    product_candidates = []
    for row in pairs["pairs"]:
        masks = tuple(map(int, row["orbit_masks"]))
        if len(masks) != 2 or any(mask not in by_mask for mask in masks):
            raise ValueError(f"unknown pair masks {masks}")
        key = f"{masks[0]}:{masks[1]}"
        if key != row["pair_key"]:
            raise ValueError(f"pair key mismatch: {key} != {row['pair_key']}")
        product_candidates.append(
            Candidate(
                kind="product",
                key=key,
                masks=masks,
                coefficients=multiply_integer_polynomials(
                    by_mask[masks[0]].coefficients,
                    by_mask[masks[1]].coefficients,
                ),
                forced_twist_rank=0,
                metadata={
                    "orbit_hex": list(row["orbit_hex"]),
                    "arithmetic_complexity_rank": int(row["arithmetic_complexity_rank"]),
                    "base_rank_lower_bound_nine": key in RANK_NINE_PAIR_KEYS,
                    "base_global_root_number": int(row["global_root_number"]),
                },
            )
        )

    if len(singleton_candidates) != 39120 or len(product_candidates) != 5566:
        raise ValueError(
            f"unexpected census sizes: {len(singleton_candidates)}, {len(product_candidates)}"
        )
    return singleton_candidates, product_candidates, {
        "bisection_schema": bisections["schema"],
        "pair_schema": pairs["schema"],
    }


def score_candidates(
    candidates: Sequence[Candidate],
    prime_blocks: Sequence[Sequence[int]],
    prime_data: dict[int, PrimeData],
) -> list[dict[str, object]]:
    caches = {prime: {} for prime in prime_data}
    records = []
    for candidate in candidates:
        block_scores = []
        usable_by_block = []
        local_values = {}
        for block in prime_blocks:
            weighted_sum = 0.0
            weight = 0.0
            usable = []
            for prime in block:
                value = local_twist_average(candidate, prime_data[prime], caches[prime])
                if value is None:
                    continue
                numeric = float(value)
                local_values[str(prime)] = {
                    "minus_A_p_numerator": value.numerator,
                    "minus_A_p_denominator": value.denominator,
                }
                weighted_sum += numeric * log(prime)
                weight += log(prime)
                usable.append(prime)
            score = weighted_sum / weight if weight else float("nan")
            block_scores.append(score)
            usable_by_block.append(usable)
        finite_scores = [score for score in block_scores if isfinite(score)]
        weakest = min(finite_scores) if len(finite_scores) == len(prime_blocks) else float("nan")
        mean_score = (
            sum(finite_scores) / len(finite_scores)
            if len(finite_scores) == len(prime_blocks)
            else float("nan")
        )
        records.append(
            {
                "kind": candidate.kind,
                "key": candidate.key,
                "masks": list(candidate.masks),
                "forced_twist_rank": candidate.forced_twist_rank,
                "block_scores": block_scores,
                "mean_block_score": mean_score,
                "weakest_block_score": weakest,
                "weakest_block_excess_over_forced_rank": (
                    weakest - candidate.forced_twist_rank if isfinite(weakest) else weakest
                ),
                "usable_primes_by_block": usable_by_block,
                "local_minus_A_p": local_values,
                **candidate.metadata,
            }
        )
    return records


def ranking_key(record: dict[str, object]) -> tuple[object, ...]:
    score = float(record["weakest_block_score"])
    return (
        0 if isfinite(score) else 1,
        -score if isfinite(score) else 0.0,
        -sum(float(value) for value in record["block_scores"] if isfinite(float(value))),
        record["key"],
    )


def summarize_records(records: Sequence[dict[str, object]], top_count: int) -> dict[str, object]:
    ordered = sorted(records, key=ranking_key)
    for rank, record in enumerate(ordered, start=1):
        record["census_rank"] = rank
    finite = [
        float(record["weakest_block_score"])
        for record in ordered
        if isfinite(float(record["weakest_block_score"]))
    ]
    means = sorted(
        float(record["mean_block_score"])
        for record in ordered
        if isfinite(float(record["mean_block_score"]))
    )
    finite_sorted = sorted(finite)

    def quantile(values: Sequence[float], numerator: int, denominator: int) -> float | None:
        if not values:
            return None
        position = (len(values) - 1) * numerator / denominator
        lower = int(position)
        upper = min(lower + 1, len(values) - 1)
        fraction = position - lower
        return values[lower] * (1.0 - fraction) + values[upper] * fraction

    score_ledger_sha256 = sha256(
        "\n".join(
            record["key"]
            + "|"
            + "|".join(format(float(value), ".17g") for value in record["block_scores"])
            for record in sorted(records, key=lambda item: item["key"])
        ).encode("ascii")
    ).hexdigest()
    designated = [
        record
        for record in ordered
        if record.get("base_rank_lower_bound_nine") or record["key"] in RANK_NINE_PAIR_KEYS
    ]
    return {
        "candidate_count": len(records),
        "finite_score_count": len(finite),
        "weakest_block_score_range": [min(finite), max(finite)] if finite else None,
        "complete_score_ledger_sha256": score_ledger_sha256,
        "score_distribution": {
            "weakest_block_quantiles": {
                "q00": quantile(finite_sorted, 0, 100),
                "q25": quantile(finite_sorted, 25, 100),
                "q50": quantile(finite_sorted, 50, 100),
                "q75": quantile(finite_sorted, 75, 100),
                "q90": quantile(finite_sorted, 90, 100),
                "q99": quantile(finite_sorted, 99, 100),
                "q100": quantile(finite_sorted, 100, 100),
            },
            "mean_block_quantiles": {
                "q00": quantile(means, 0, 100),
                "q25": quantile(means, 25, 100),
                "q50": quantile(means, 50, 100),
                "q75": quantile(means, 75, 100),
                "q90": quantile(means, 90, 100),
                "q99": quantile(means, 99, 100),
                "q100": quantile(means, 100, 100),
            },
            "weakest_block_threshold_counts": {
                format(threshold, ".2f"): sum(value >= threshold for value in finite)
                for threshold in (-0.5, 0.0, 0.25, 0.5, 0.75, 1.0)
            },
        },
        "top": ordered[:top_count],
        "designated_rank_nine_bases": designated,
    }


def direct_trace(prime: int, coefficient_a: int, coefficient_b: int, characters: Sequence[int]) -> int:
    return -sum(
        characters[(x_value**3 + coefficient_a * x_value + coefficient_b) % prime]
        for x_value in range(prime)
    )


def verify_trace_identity(model, candidate: Candidate, data: PrimeData) -> dict[str, object]:
    reduced = reduce_twist_polynomial(candidate.coefficients, data.prime)
    if reduced is None:
        raise ValueError("self-test candidate is bad at the selected prime")
    a_coefficients, b_coefficients = reduced_coefficients(model, data.prime)
    checked = 0
    for parameter in range(data.prime):
        q_value = evaluate_mod(reduced, parameter, data.prime)
        coefficient_a = evaluate_mod(a_coefficients, parameter, data.prime)
        coefficient_b = evaluate_mod(b_coefficients, parameter, data.prime)
        if q_value == 0 or (4 * coefficient_a**3 + 27 * coefficient_b**2) % data.prime == 0:
            continue
        twist_a = coefficient_a * q_value**2 % data.prime
        twist_b = coefficient_b * q_value**3 % data.prime
        twist_trace = direct_trace(data.prime, twist_a, twist_b, data.characters)
        expected = data.characters[q_value] * data.traces[parameter]
        if twist_trace != expected:
            raise ArithmeticError(
                f"trace identity failed at p={data.prime}, t={parameter}: {twist_trace} != {expected}"
            )
        checked += 1
        if checked == 5:
            break
    if checked < 5:
        raise ArithmeticError("too few good fibres for trace self-test")
    return {"prime": data.prime, "candidate_key": candidate.key, "fibres_checked": checked}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bisections", type=Path, default=DEFAULT_BISECTIONS)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--prime-block",
        type=parse_prime_block,
        action="append",
        help="comma-separated disjoint prime block; repeat for multiple blocks",
    )
    parser.add_argument("--top", type=int, default=100)
    parser.add_argument("--singletons-only", action="store_true")
    parser.add_argument("--products-only", action="store_true")
    parser.add_argument(
        "--singleton-key",
        action="append",
        default=[],
        help="score only this singleton mask; repeat to retain several",
    )
    parser.add_argument(
        "--product-key",
        action="append",
        default=[],
        help="score only this colon-separated pair key; repeat to retain several",
    )
    args = parser.parse_args()
    if args.top < 1 or (args.singletons_only and args.products_only):
        raise ValueError("invalid output size or mutually exclusive census flags")

    prime_blocks = tuple(args.prime_block or DEFAULT_PRIME_BLOCKS)
    flattened = [prime for block in prime_blocks for prime in block]
    if len(flattened) != len(set(flattened)):
        raise ValueError("prime blocks must be pairwise disjoint")

    started = perf_counter()
    model = load_family_model(args.model)
    all_singleton_candidates, all_product_candidates, schemas = load_candidates(
        args.bisections, args.pairs
    )
    singleton_keys = set(args.singleton_key)
    product_keys = set(args.product_key)
    singleton_candidates = [
        candidate
        for candidate in all_singleton_candidates
        if not singleton_keys or candidate.key in singleton_keys
    ]
    product_candidates = [
        candidate
        for candidate in all_product_candidates
        if not product_keys or candidate.key in product_keys
    ]
    if singleton_keys and singleton_keys != {candidate.key for candidate in singleton_candidates}:
        missing = sorted(singleton_keys - {candidate.key for candidate in singleton_candidates})
        raise ValueError(f"unknown singleton keys: {missing}")
    if product_keys and product_keys != {candidate.key for candidate in product_candidates}:
        missing = sorted(product_keys - {candidate.key for candidate in product_candidates})
        raise ValueError(f"unknown product keys: {missing}")
    prime_data = {prime: build_prime_data(model, prime) for prime in flattened}
    self_test = None
    for prime in flattened:
        for candidate in all_singleton_candidates:
            if reduce_twist_polynomial(candidate.coefficients, prime) is not None:
                self_test = verify_trace_identity(model, candidate, prime_data[prime])
                break
        if self_test is not None:
            break
    if self_test is None:
        raise ArithmeticError("no usable singleton/prime pair for trace self-test")

    sections = {}
    if not args.products_only:
        singleton_records = score_candidates(singleton_candidates, prime_blocks, prime_data)
        sections["singleton_twists"] = summarize_records(singleton_records, args.top)
    if not args.singletons_only:
        product_records = score_candidates(product_candidates, prime_blocks, prime_data)
        sections["immediate_pair_product_twists"] = summarize_records(product_records, args.top)

    result = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_BOUNDED_HEURISTIC_QUADRATIC_TWIST_RANK_CENSUS",
        "proof_boundary": (
            "These finite-prime Nagao-style scores are search heuristics, not rank bounds. "
            "A singleton candidate needs an exact second twist section to prove twist rank at "
            "least two; a product candidate needs an exact non-torsion product-twist section "
            "to prove paired-base surface rank at least twenty."
        ),
        "method": {
            "trace_identity": "a_p(E^q_t)=chi_p(q(t))*a_p(E_t)",
            "fibral_average": "A_p=(1/p)*sum_{t in F_p} a_p(E^q_t)",
            "local_score": "-A_p",
            "block_score": "sum_p(-A_p*log(p))/sum_p(log(p)) over usable primes",
            "ranking": "descending weakest block score, then descending block-score sum",
            "singular_fibres": "trace set to zero",
            "bad_twist_primes": "skipped after p-adic squareclass normalization or branch collision",
            "constant_squareclass_preserved": True,
            "self_test": self_test,
        },
        "inputs": {
            "bisections": {
                "path": display_path(args.bisections),
                "sha256": digest(args.bisections),
                "schema": schemas["bisection_schema"],
            },
            "pairs": {
                "path": display_path(args.pairs),
                "sha256": digest(args.pairs),
                "schema": schemas["pair_schema"],
            },
            "model": {
                "path": display_path(args.model),
                "sha256": model.source_sha256,
                "coordinate": model.coordinate,
            },
        },
        "prime_blocks": [list(block) for block in prime_blocks],
        "prime_tables": {
            str(prime): {
                "singular_fibre_count": data.singular_fibre_count,
                "trace_sha256": sha256(
                    ",".join(map(str, data.traces)).encode("ascii")
                ).hexdigest(),
            }
            for prime, data in prime_data.items()
        },
        **sections,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    leaders = []
    for name, section in sections.items():
        leaders.append(f"{name}={section['top'][0]['key']}")
    print(
        f"PASS primes={len(flattened)} {' '.join(leaders)} "
        f"seconds={result['runtime_seconds']:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
