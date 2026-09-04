#!/usr/bin/env python3
"""Bounded Nagao sieve for rational R17 quadratic characters.

The search family is

    q(t) = t^2 + b*t + c,  |b|,|c| <= H,  b^2-4*c != 0.

The monic leading coefficient gives two rational points above infinity on
``u^2=q(t)``, so every retained character has rational base.  Scores are
finite-prime heuristics only; they are not Mordell--Weil rank bounds.
"""

# status: ACTIVE_SEARCH
# claim: bounded heuristic ranking in one compact published-R17 coordinate
# inputs: elkies-k3/data/fibrations/elkies_2026_published_r17_model.json
# outputs: artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-nagao-h100-v1.json

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import log
from pathlib import Path
import shlex
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
COMMON_PATH = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
MODEL_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-nagao-h100-v1.json"
)
DEFAULT_PRIME_BLOCKS = (
    (211, 223, 227, 229, 233, 239, 241, 251),
    (257, 263, 269, 271, 277, 281, 283, 293),
    (307, 311, 313, 317, 331, 337, 347, 349),
    (353, 359, 367, 373, 379, 383, 389, 397),
    (401, 409, 419, 421, 431, 433, 439, 443),
    (449, 457, 461, 463, 467, 479, 487, 491),
)


@dataclass
class Candidate:
    b: int
    c: int
    block_scores: list[float]

    @property
    def height(self) -> int:
        return max(abs(self.b), abs(self.c))

    @property
    def key(self) -> tuple[int, int]:
        return self.b, self.c


def load_common():
    spec = importlib.util.spec_from_file_location("r17_quadratic_twist_common", COMMON_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {COMMON_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_prime_blocks(text: str | None) -> tuple[tuple[int, ...], ...]:
    if text is None:
        return DEFAULT_PRIME_BLOCKS
    blocks = []
    for raw_block in text.split(";"):
        values = tuple(int(value) for value in raw_block.split(",") if value.strip())
        if not values:
            raise argparse.ArgumentTypeError("empty prime block")
        blocks.append(values)
    if not blocks:
        raise argparse.ArgumentTypeError("no prime blocks")
    return tuple(blocks)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def polynomial_value(coefficients: tuple[int, ...], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def build_prime_data(common, model, prime: int) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    a_coefficients, b_coefficients = common.reduced_coefficients(model, prime)
    characters = common.quadratic_character_table(prime)
    traces = []
    singular = 0
    for parameter in range(prime):
        coefficient_a = polynomial_value(a_coefficients, parameter, prime)
        coefficient_b = polynomial_value(b_coefficients, parameter, prime)
        if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
            traces.append(0)
            singular += 1
            continue
        trace = -sum(
            characters[(x_value**3 + coefficient_a * x_value + coefficient_b) % prime]
            for x_value in range(prime)
        )
        traces.append(trace)
    return tuple(traces), characters, singular


def local_score(candidate: Candidate, prime_data, prime: int) -> Fraction | None:
    traces, characters, _ = prime_data
    b = candidate.b % prime
    c = candidate.c % prime
    # A repeated branch point is bad reduction for the twist surface.
    if (b * b - 4 * c) % prime == 0:
        return None
    trace_sum = 0
    for parameter, trace in enumerate(traces):
        q_value = (parameter * parameter + b * parameter + c) % prime
        trace_sum += characters[q_value] * trace
    return Fraction(-trace_sum, prime)


def block_score(candidate: Candidate, block, data) -> tuple[float, list[dict[str, object]]]:
    numerator = 0.0
    denominator = 0.0
    rows = []
    for prime in block:
        score = local_score(candidate, data[prime], prime)
        if score is None:
            rows.append({"prime": prime, "status": "skipped_bad_twist_reduction"})
            continue
        weight = log(prime)
        numerator += float(score) * weight
        denominator += weight
        rows.append(
            {
                "prime": prime,
                "negative_fibral_average": f"{score.numerator}/{score.denominator}",
            }
        )
    if denominator == 0:
        return float("-inf"), rows
    return numerator / denominator, rows


def candidate_sort_key(candidate: Candidate):
    weakest = min(candidate.block_scores)
    mean = sum(candidate.block_scores) / len(candidate.block_scores)
    return (-weakest, -mean, candidate.height, candidate.b, candidate.c)


def retain_per_height(candidates: list[Candidate], keep: int) -> list[Candidate]:
    buckets: dict[int, list[Candidate]] = {}
    for candidate in candidates:
        buckets.setdefault(candidate.height, []).append(candidate)
    retained = []
    for height in sorted(buckets):
        retained.extend(sorted(buckets[height], key=candidate_sort_key)[:keep])
    return retained


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficient-bound", type=int, default=100)
    parser.add_argument("--keep-per-height", type=int, default=32)
    parser.add_argument("--finalists", type=int, default=200)
    parser.add_argument("--prime-blocks")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.coefficient_bound < 1 or args.keep_per_height < 1 or args.finalists < 1:
        parser.error("bounds and retention counts must be positive")

    common = load_common()
    model = common.load_family_model(MODEL_PATH)
    prime_blocks = parse_prime_blocks(args.prime_blocks)
    primes = tuple(dict.fromkeys(prime for block in prime_blocks for prime in block))
    started = perf_counter()
    data = {prime: build_prime_data(common, model, prime) for prime in primes}

    bound = args.coefficient_bound
    candidates = [
        Candidate(b=b, c=c, block_scores=[])
        for b in range(-bound, bound + 1)
        for c in range(-bound, bound + 1)
        if b * b != 4 * c
    ]
    initial_count = len(candidates)
    stages = []
    for block_number, block in enumerate(prime_blocks, start=1):
        for candidate in candidates:
            score, _ = block_score(candidate, block, data)
            candidate.block_scores.append(score)
        before = len(candidates)
        candidates = retain_per_height(candidates, args.keep_per_height)
        stages.append(
            {
                "block_number": block_number,
                "primes": list(block),
                "before": before,
                "after": len(candidates),
                "best_weakest_block_score": max(
                    min(candidate.block_scores) for candidate in candidates
                )
                if candidates
                else None,
            }
        )

    candidates.sort(key=candidate_sort_key)
    finalist_records = []
    for rank, candidate in enumerate(candidates[: args.finalists], start=1):
        local_rows = []
        for block_number, block in enumerate(prime_blocks, start=1):
            _, rows = block_score(candidate, block, data)
            local_rows.append({"block_number": block_number, "rows": rows})
        weakest = min(candidate.block_scores)
        mean = sum(candidate.block_scores) / len(candidate.block_scores)
        finalist_records.append(
            {
                "rank": rank,
                "q_coefficients_low_to_high": [candidate.c, candidate.b, 1],
                "formula": f"t^2{candidate.b:+d}*t{candidate.c:+d}",
                "coefficient_height": candidate.height,
                "discriminant": candidate.b * candidate.b - 4 * candidate.c,
                "rational_base_witness": "two QQ-points above t=infinity with u/t=+/-1",
                "block_scores": candidate.block_scores,
                "weakest_block_score": weakest,
                "mean_block_score": mean,
                "local_scores": local_rows,
            }
        )

    payload = {
        "schema": "elkies-k3.r17-rational-quadratic-twist-nagao.v1",
        "status": "PASS_BOUNDED_HEURISTIC_RATIONAL_QUADRATIC_TWIST_SIEVE",
        "model": {
            "source": str(MODEL_PATH.relative_to(ROOT)),
            "source_sha256": digest(MODEL_PATH),
            "coordinate": model.coordinate,
        },
        "search": {
            "family": "q(t)=t^2+b*t+c",
            "coefficient_bound": bound,
            "initial_squarefree_character_count": initial_count,
            "keep_per_exact_coefficient_height_per_stage": args.keep_per_height,
            "prime_blocks": [list(block) for block in prime_blocks],
            "singular_fibre_counts": {
                str(prime): data[prime][2] for prime in primes
            },
        },
        "stages": stages,
        "retained_count": len(candidates),
        "finalists": finalist_records,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(
            argument for argument in sys.argv if argument != "--check"
        ),
        "proof_boundary": (
            "Every searched character has a rational genus-zero base and nonzero "
            "characteristic-zero discriminant. The finite-prime Nagao scores are "
            "heuristic rankings only. No section, twist-rank lower bound, MW20 "
            "surface, specialization transport, or tail-survival claim follows."
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        existing = json.loads(args.output.read_text())
        existing.pop("runtime_seconds", None)
        payload.pop("runtime_seconds", None)
        if existing != payload:
            raise SystemExit("stored artifact differs from replay")
        print(f"PASS check {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(
        "PASS bounded rational quadratic twist sieve "
        f"initial={initial_count} retained={len(candidates)} "
        f"best={finalist_records[0]['weakest_block_score']:.6f} "
        f"seconds={payload['runtime_seconds']:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
