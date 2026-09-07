#!/usr/bin/env python3
"""Search rational R17 quadratic characters in all eight lineage coordinates.

For every exact rational-PGL2 coordinate ``z`` in the published-R17 lineage,
the default bounded family is

    q_z(z) = z^2 + b*z + c,  |b|,|c| <= H,  b^2-4*c != 0.

An optional control anchor instead uses

    q_z(z) = 1 + b*(z-z0) + c*(z-z0)^2,

which supplies a rational cover point above the selected exact control.
Pullback through ``z=(a*t+b0)/(c0*t+d)`` gives a quadratic polynomial in the
native 074d9 coordinate after discarding a square denominator. Scores are
finite-prime heuristics only; they are not Mordell--Weil rank bounds.
"""

# status: ACTIVE_SEARCH
# claim: bounded heuristic ranking in the eight exact 074d9-lineage coordinates
# inputs: artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json
# outputs: artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-*-v1.json

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import shlex
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
BASE_SCRIPT = ROOT / "elkies-k3/scripts/search_r17_rational_quadratic_twist_nagao.py"
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-pgl8-nagao-h100-v1.json"
)


@dataclass
class Candidate:
    chart: str
    anchor_curve_id: int | None
    anchor_parameter: Fraction | None
    b: int
    c: int
    q_coefficients: tuple[int, int, int]
    block_scores: list[float]

    @property
    def height(self) -> int:
        return max(abs(self.b), abs(self.c))


def load_base():
    spec = importlib.util.spec_from_file_location("r17_rational_twist_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def representative_model(base, common, lineage):
    record = lineage["representative"]
    if record.get("chart") != "norm12-orbit-074d9":
        raise ValueError("the exact lineage representative changed")
    a_coefficients = tuple(Fraction(value) for value in record["A_coefficients_low_to_high"])
    b_coefficients = tuple(Fraction(value) for value in record["B_coefficients_low_to_high"])
    if len(a_coefficients) != 9 or len(b_coefficients) != 13:
        raise ValueError("the lineage representative is not a degree-(8,12) model")
    return common.FamilyModel(
        source=LINEAGE.resolve(),
        source_sha256=digest(LINEAGE),
        a_coefficients=a_coefficients,
        b_coefficients=b_coefficients,
        a_degree=8,
        b_degree=12,
        coordinate="norm12-orbit-074d9 native u",
        coefficient_source_keys=(
            "representative.A_coefficients_low_to_high",
            "representative.B_coefficients_low_to_high",
        ),
    )


def pulled_quadratic(matrix: tuple[int, int, int, int], b: int, c: int):
    a, b0, c0, d = matrix
    return (
        b0 * b0 + b * b0 * d + c * d * d,
        2 * a * b0 + b * (a * d + b0 * c0) + 2 * c * c0 * d,
        a * a + b * a * c0 + c * c0 * c0,
    )


def anchored_pulled_quadratic(
    matrix: tuple[int, int, int, int],
    anchor: Fraction,
    b: int,
    c: int,
):
    """Pull back 1+b*(z-z0)+c*(z-z0)^2, clearing a square denominator."""

    a, b0, c0, d = matrix
    numerator = anchor.numerator
    denominator = anchor.denominator
    # D=c0*t+d and L=denominator*N-numerator*D.
    d_linear = (d, c0)
    l_linear = (
        denominator * b0 - numerator * d,
        denominator * a - numerator * c0,
    )

    def product(left, right):
        return (
            left[0] * right[0],
            left[0] * right[1] + left[1] * right[0],
            left[1] * right[1],
        )

    d_squared = product(d_linear, d_linear)
    ld = product(l_linear, d_linear)
    l_squared = product(l_linear, l_linear)
    return tuple(
        denominator * denominator * d_squared[index]
        + b * denominator * ld[index]
        + c * l_squared[index]
        for index in range(3)
    )


def local_score(base, candidate: Candidate, prime_data, prime: int):
    traces, characters, _ = prime_data
    q0, q1, q2 = (value % prime for value in candidate.q_coefficients)
    if q0 == q1 == q2 == 0 or (q1 * q1 - 4 * q0 * q2) % prime == 0:
        return None
    trace_sum = 0
    for parameter, trace in enumerate(traces):
        q_value = (q2 * parameter * parameter + q1 * parameter + q0) % prime
        trace_sum += characters[q_value] * trace
    return Fraction(-trace_sum, prime)


def block_score(base, candidate: Candidate, block, data):
    from math import log

    numerator = 0.0
    denominator = 0.0
    rows = []
    for prime in block:
        score = local_score(base, candidate, data[prime], prime)
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
    anchor = -1 if candidate.anchor_curve_id is None else candidate.anchor_curve_id
    return (
        -weakest,
        -mean,
        candidate.height,
        candidate.chart,
        anchor,
        candidate.b,
        candidate.c,
    )


def retain_per_chart_height(candidates: list[Candidate], keep: int):
    buckets: dict[tuple[str, int | None, int], list[Candidate]] = {}
    for candidate in candidates:
        buckets.setdefault(
            (candidate.chart, candidate.anchor_curve_id, candidate.height), []
        ).append(candidate)
    retained = []
    for key in sorted(buckets):
        retained.extend(sorted(buckets[key], key=candidate_sort_key)[:keep])
    return retained


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficient-bound", type=int, default=100)
    parser.add_argument("--keep-per-chart-height", type=int, default=32)
    parser.add_argument("--finalists", type=int, default=200)
    parser.add_argument("--prime-blocks")
    parser.add_argument(
        "--anchor-control",
        type=int,
        action="append",
        choices=(351, 356, 376, 377, 385),
        help=(
            "search q(z)=1+b*(z-z0)+c*(z-z0)^2 above this exact control; "
            "repeat for several controls; omission searches q(z)=z^2+b*z+c"
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.coefficient_bound < 1 or args.keep_per_chart_height < 1 or args.finalists < 1:
        parser.error("bounds and retention counts must be positive")

    base = load_base()
    common = base.load_common()
    lineage = json.loads(LINEAGE.read_text())
    if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ValueError("the exact R17 lineage certificate is unavailable")
    model = representative_model(base, common, lineage)
    maps = {
        chart: tuple(int(value) for value in record["representative_to_member_base_map"]["a_b_c_d"])
        for chart, record in lineage["chart_transports"].items()
    }
    if len(maps) != 8:
        raise ValueError("expected eight exact PGL2 coordinates in the R17 lineage")
    anchor_ids = tuple(dict.fromkeys(args.anchor_control or ()))
    anchor_parameters = {
        (record["chart"], int(record["curve_id"])): Fraction(record["parameter"])
        for record in lineage["target_isomorphisms"]
        if int(record["curve_id"]) in anchor_ids
    }
    if anchor_ids and len(anchor_parameters) != 8 * len(anchor_ids):
        raise ValueError("the selected controls are not present in all lineage coordinates")

    prime_blocks = base.parse_prime_blocks(args.prime_blocks)
    primes = tuple(dict.fromkeys(prime for block in prime_blocks for prime in block))
    started = perf_counter()
    data = {prime: base.build_prime_data(common, model, prime) for prime in primes}
    bound = args.coefficient_bound
    candidates = []
    for chart, matrix in sorted(maps.items()):
        chart_anchors = anchor_ids or (None,)
        for anchor_curve_id in chart_anchors:
            anchor = (
                None
                if anchor_curve_id is None
                else anchor_parameters[(chart, anchor_curve_id)]
            )
            for b in range(-bound, bound + 1):
                for c in range(-bound, bound + 1):
                    if anchor is None:
                        if b * b == 4 * c:
                            continue
                        coefficients = pulled_quadratic(matrix, b, c)
                    else:
                        # Discriminant of 1+b*w+c*w^2 is b^2-4c.
                        if b * b == 4 * c:
                            continue
                        coefficients = anchored_pulled_quadratic(
                            matrix, anchor, b, c
                        )
                    # Keep the exact scalar. Dividing common nonsquare content would
                    # change the arithmetic quadratic character.
                    candidates.append(
                        Candidate(
                            chart,
                            anchor_curve_id,
                            anchor,
                            b,
                            c,
                            coefficients,
                            [],
                        )
                    )
    initial_count = len(candidates)

    stages = []
    for block_number, block in enumerate(prime_blocks, start=1):
        for candidate in candidates:
            score, _ = block_score(base, candidate, block, data)
            candidate.block_scores.append(score)
        before = len(candidates)
        candidates = retain_per_chart_height(candidates, args.keep_per_chart_height)
        stages.append(
            {
                "block_number": block_number,
                "primes": list(block),
                "before": before,
                "after": len(candidates),
                "best_weakest_block_score": max(
                    min(candidate.block_scores) for candidate in candidates
                ),
            }
        )

    candidates.sort(key=candidate_sort_key)
    finalist_records = []
    for rank, candidate in enumerate(candidates[: args.finalists], start=1):
        local_rows = []
        for block_number, block in enumerate(prime_blocks, start=1):
            _, rows = block_score(base, candidate, block, data)
            local_rows.append({"block_number": block_number, "rows": rows})
        if candidate.anchor_parameter is None:
            q_lineage = [candidate.c, candidate.b, 1]
            witness = "two QQ-points above infinity in the displayed lineage coordinate"
        else:
            z0 = candidate.anchor_parameter
            q_lineage = [
                Fraction(1) - candidate.b * z0 + candidate.c * z0 * z0,
                Fraction(candidate.b) - 2 * candidate.c * z0,
                Fraction(candidate.c),
            ]
            q_lineage = [f"{value.numerator}/{value.denominator}" for value in q_lineage]
            witness = (
                f"q(z0)=1 above curve {candidate.anchor_curve_id}, "
                f"z0={z0.numerator}/{z0.denominator}"
            )
        finalist_record = {
            "rank": rank,
            "lineage_coordinate": candidate.chart,
            "anchor_curve_id": candidate.anchor_curve_id,
            "anchor_parameter": (
                None
                if candidate.anchor_parameter is None
                else f"{candidate.anchor_parameter.numerator}/{candidate.anchor_parameter.denominator}"
            ),
            "q_in_lineage_coordinate_coefficients_low_to_high": q_lineage,
            "q_in_074d9_coefficients_low_to_high": list(candidate.q_coefficients),
            "coefficient_height_in_lineage_coordinate": candidate.height,
            "discriminant_in_lineage_coordinate": candidate.b * candidate.b - 4 * candidate.c,
            "rational_base_witness": witness,
            "block_scores": candidate.block_scores,
            "weakest_block_score": min(candidate.block_scores),
            "mean_block_score": sum(candidate.block_scores) / len(candidate.block_scores),
            "local_scores": local_rows,
        }
        if candidate.anchor_parameter is not None:
            finalist_record["family_formula"] = "q(z)=1+b*(z-z0)+c*(z-z0)^2"
        finalist_records.append(finalist_record)

    payload = {
        "schema": "elkies-k3.r17-rational-quadratic-twist-pgl8-nagao.v1",
        "status": "PASS_BOUNDED_HEURISTIC_RATIONAL_QUADRATIC_TWIST_PGL8_SIEVE",
        "model": {
            "source": str(LINEAGE.relative_to(ROOT)),
            "source_sha256": digest(LINEAGE),
            "coordinate": model.coordinate,
        },
        "coordinate_maps": {
            chart: list(matrix) for chart, matrix in sorted(maps.items())
        },
        "search": {
            "family": (
                "q_z(z)=z^2+b*z+c in each of eight exact R17 lineage coordinates"
                if not anchor_ids
                else "q_z(z)=1+b*(z-z0)+c*(z-z0)^2 above selected exact controls"
            ),
            "anchor_control_ids": list(anchor_ids),
            "coefficient_bound": bound,
            "initial_squarefree_character_count_with_multiplicity": initial_count,
            "keep_per_lineage_coordinate_and_exact_height_per_stage": args.keep_per_chart_height,
            "prime_blocks": [list(block) for block in prime_blocks],
            "singular_fibre_counts": {str(prime): data[prime][2] for prime in primes},
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
            "characteristic-zero discriminant. Exact PGL2 maps preserve the quadratic "
            "squareclass, including its scalar. The finite-prime Nagao scores are "
            "heuristic rankings only. No section, twist-rank lower bound, MW20 surface, "
            "specialization transport, or tail-survival claim follows."
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
        "PASS bounded PGL8 rational quadratic twist sieve "
        f"initial={initial_count} retained={len(candidates)} "
        f"best={finalist_records[0]['weakest_block_score']:.6f} "
        f"seconds={payload['runtime_seconds']:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
