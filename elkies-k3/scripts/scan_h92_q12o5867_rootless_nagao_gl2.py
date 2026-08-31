#!/usr/bin/env python3
"""Run the q12/orbit5867 Nagao sieve in an integral PGL2 chart.

For a primitive matrix ``[[alpha,beta],[gamma,delta]]`` this wrapper searches
balanced primitive chart parameters ``v=(a:b)`` and maps them to

    u = (alpha*a + beta*b : gamma*a + delta*b).

It permutes the exact P^1(F_p) lookup tables before invoking the existing C++
hot loop, then transports every complete finalist back to the original
projective u-coordinate.  The matrix determinant must be nonzero over Q and
invertible at every discovery prime.  This is a bounded heuristic scan, not
rank evidence.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from typing import Sequence


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
sys.path.insert(0, str(SCRIPTS))

from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    DEFAULT_PRIME_BLOCKS,
    build_residue_tables,
    export_cpp_tables,
    load_family_model,
)


DEFAULT_SCANNER = ROOT / "artifacts/local/elkies-k3/scan-q12o5867-rootless-nagao-skew"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_matrix(values: Sequence[int]) -> tuple[int, int, int, int]:
    if len(values) != 4:
        raise ValueError("a chart matrix requires exactly four entries")
    alpha, beta, gamma, delta = map(int, values)
    common = gcd(gcd(abs(alpha), abs(beta)), gcd(abs(gamma), abs(delta)))
    if common:
        alpha, beta, gamma, delta = (
            alpha // common,
            beta // common,
            gamma // common,
            delta // common,
        )
    for value in (alpha, beta, gamma, delta):
        if value:
            if value < 0:
                alpha, beta, gamma, delta = (
                    -alpha,
                    -beta,
                    -gamma,
                    -delta,
                )
            break
    if alpha * delta - beta * gamma == 0:
        raise ValueError("chart matrix must have nonzero determinant")
    return alpha, beta, gamma, delta


def map_projective_pair(
    numerator: int,
    denominator: int,
    matrix: Sequence[int],
) -> tuple[int, int, int]:
    alpha, beta, gamma, delta = matrix
    mapped_numerator = alpha * numerator + beta * denominator
    mapped_denominator = gamma * numerator + delta * denominator
    if mapped_denominator == 0:
        if mapped_numerator == 0:
            raise ValueError("singular chart matrix maps a point to (0:0)")
        return 1, 0, 1
    common = gcd(abs(mapped_numerator), abs(mapped_denominator))
    mapped_numerator //= common
    mapped_denominator //= common
    if mapped_denominator < 0:
        mapped_numerator, mapped_denominator = -mapped_numerator, -mapped_denominator
    return (
        mapped_numerator,
        mapped_denominator,
        max(abs(mapped_numerator), mapped_denominator),
    )


def transformed_blocks(
    blocks: Sequence[dict[int, tuple[object, ...]]],
    matrix: Sequence[int],
) -> tuple[dict[int, tuple[object, ...]], ...]:
    alpha, beta, gamma, delta = matrix
    determinant = alpha * delta - beta * gamma
    transformed = []
    for block in blocks:
        transformed_block = {}
        for prime, table in block.items():
            if determinant % prime == 0:
                raise ValueError(
                    f"chart determinant is not invertible at discovery prime p={prime}"
                )
            symbols = []
            mapped_indices = []
            for chart_index in range(prime + 1):
                a, b = (
                    (chart_index, 1)
                    if chart_index < prime
                    else (1, 0)
                )
                mapped_a = (alpha * a + beta * b) % prime
                mapped_b = (gamma * a + delta * b) % prime
                if mapped_a == mapped_b == 0:
                    raise AssertionError("invertible chart produced (0:0) modulo p")
                original_index = (
                    prime
                    if mapped_b == 0
                    else mapped_a * pow(mapped_b, -1, prime) % prime
                )
                mapped_indices.append(original_index)
                symbols.append(table[original_index])
            if len(set(mapped_indices)) != prime + 1:
                raise AssertionError("PGL2 chart did not permute P1(F_p)")
            transformed_block[prime] = tuple(symbols)
        transformed.append(transformed_block)
    return tuple(transformed)


def finalist_record(
    raw: dict[str, object], matrix: Sequence[int]
) -> dict[str, object]:
    chart_numerator, chart_denominator = map(int, raw["projective_pair"])
    numerator, denominator, height = map_projective_pair(
        chart_numerator, chart_denominator, matrix
    )
    result = dict(raw)
    result.update(
        {
            "parameter": (
                "infinity" if denominator == 0 else f"{numerator}/{denominator}"
            ),
            "projective_pair": [numerator, denominator],
            "projective_height": height,
            "chart_parameter": raw["parameter"],
            "chart_projective_pair": raw["projective_pair"],
            "chart_height": raw["projective_height"],
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bounded PGL2-chart Nagao scan for q12/orbit5867."
    )
    parser.add_argument("--matrix", nargs=4, type=int, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--scanner", type=Path, default=DEFAULT_SCANNER)
    parser.add_argument("--numerator-bound", type=int, default=10000)
    parser.add_argument("--denominator-bound", type=int, default=10000)
    parser.add_argument("--bucket-width", type=int, default=100)
    parser.add_argument("--keep-per-bucket", default="32,16,8")
    parser.add_argument("--finalists", type=int, default=10000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if min(
        args.numerator_bound,
        args.denominator_bound,
        args.bucket_width,
        args.finalists,
    ) < 1:
        raise SystemExit("bounds, bucket width, and finalists must be positive")
    matrix = normalize_matrix(args.matrix)
    determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2]
    model = load_family_model()
    base_blocks, rejected = build_residue_tables(model, DEFAULT_PRIME_BLOCKS)
    if rejected:
        raise ValueError(f"unexpected rejected discovery primes: {rejected}")
    transformed = transformed_blocks(base_blocks, matrix)

    with tempfile.TemporaryDirectory(prefix="q12o5867-gl2-") as temporary:
        temporary_path = Path(temporary)
        table_path = temporary_path / "tables.txt"
        raw_path = temporary_path / "raw.json"
        export_cpp_tables(table_path, model, transformed)
        command = [
            str(args.scanner),
            str(table_path),
            str(args.numerator_bound),
            str(args.denominator_bound),
            str(args.bucket_width),
            args.keep_per_bucket,
            str(args.finalists),
            str(raw_path),
        ]
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
        raw = json.loads(raw_path.read_text())

    if len(raw["finalists"]) != raw["final_survivor_count"]:
        raise ValueError(
            "scanner output truncated the survivor population; increase --finalists"
        )
    finalists = [finalist_record(record, matrix) for record in raw["finalists"]]
    if len({tuple(record["projective_pair"]) for record in finalists}) != len(finalists):
        raise AssertionError("an invertible chart produced duplicate rational parameters")
    population_text = "".join(f"{record['parameter']}\n" for record in finalists).encode()

    output = {
        "schema": "h92-q12o5867-rootless-projective-nagao-cpp-gl2-v1",
        "status": "PASS_BOUNDED_HEURISTIC_PROJECTIVE_NAGAO_CPP_SIEVE",
        "proof_boundary": (
            "This is a bounded table-lookup Nagao scan in one exact PGL2 chart. "
            "Scores and survival are heuristics, not rank evidence; no section "
            "or specialized point was evaluated."
        ),
        "model_sha256": model.source_sha256,
        "scanner": {
            "path": str(args.scanner.resolve()),
            "sha256": file_sha256(args.scanner),
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
        },
        "search": {
            "chart_label": args.label,
            "parameter_chart": "u=(alpha*v_num+beta*v_den)/(gamma*v_num+delta*v_den)",
            "chart_matrix_alpha_beta_gamma_delta": list(matrix),
            "chart_matrix_determinant": determinant,
            "chart_matrix_primitive": True,
            "determinant_invertible_at_all_discovery_primes": True,
            "numerator_interval": [-args.numerator_bound, args.numerator_bound],
            "denominator_interval": [1, args.denominator_bound],
            "primitive_chart_pairs_only": True,
            "chart_height": "max(abs(v_num),v_den)",
            "height_bucket_width": args.bucket_width,
            "keep_per_bucket": [int(value) for value in args.keep_per_bucket.split(",")],
            "discovery_prime_blocks": [list(block) for block in DEFAULT_PRIME_BLOCKS],
        },
        "stages": raw["stages"],
        "final_survivor_count": len(finalists),
        "complete_finalist_population_sha256": sha256(population_text).hexdigest(),
        "finalists": finalists,
        "runtime_seconds": raw["runtime_seconds"],
        "reproducing_command": shlex.join(sys.argv),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS chart={args.label} population={raw['stages'][0]['population_scored']} "
        f"survivors={len(finalists)} seconds={raw['runtime_seconds']:.3f} "
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
