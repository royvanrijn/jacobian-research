#!/usr/bin/env python3
"""Two-stage Nagao search for rational specializations of orbit 103.

The orbit-103 surface is

    y^2 = x^3 - 27 (r^2-4)^2 A4(k,r) x
                + 54 (r^2-4)^3 B6(k,r).

Its generic arithmetic Mordell--Weil group over ``QQ(k)(r)`` has the two
known directions ``Q_plus,Q_minus``.  The third geometric direction is
anti-invariant over ``QQ(sqrt(-3))`` and is deliberately absent here.

Stage one scores a rational ``k`` by aggregating exact local traces over all
``r`` in ``P1(F_p)``.  This is a finite Nagao-style heuristic for a rank jump
of the elliptic surface over ``QQ(r)``.  Stage two fixes the leading ``k``
values and performs the repository's ordinary centered/standardized local
trace scan in ``r``.  Three pairwise-disjoint prime ensembles are used at
each stage and candidates are ranked by their weakest block.

The output is a bounded search artifact, not a rank certificate.  Exact
point discovery and finite-quotient tests against ``Q_plus,Q_minus`` are the
separate next gate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt, log, sqrt
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
SOURCE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json"
)
CPP_SOURCE = SCRIPTS / "scan_h92_q12o5867_rootless_nagao.cpp"
DEFAULT_LOCAL = ROOT / "artifacts/local/elkies-k3/e6a1-orbit103-specialization-search"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-e6a1-orbit103-specialization-search-v1.json"
)
SCORE_SCALE = 10**12


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def primes_in_interval(lower: int, upper: int) -> tuple[int, ...]:
    return tuple(value for value in range(lower, upper + 1) if is_prime(value))


def prime_ensembles(lower: int, upper: int, count: int) -> tuple[tuple[int, ...], ...]:
    primes = primes_in_interval(lower, upper)
    if count < 3 or len(primes) < count:
        raise ValueError("at least three nonempty prime ensembles are required")
    return tuple(tuple(primes[index::count]) for index in range(count))


def quadratic_characters(prime: int) -> tuple[int, ...]:
    values = [-1] * prime
    values[0] = 0
    for value in range(1, prime):
        values[value * value % prime] = 1
    return tuple(values)


def projective_value(index: int, prime: int) -> tuple[int, int]:
    if not 0 <= index <= prime:
        raise ValueError("projective index is outside P1(F_p)")
    return (index, 1) if index < prime else (1, 0)


def _powers(value: int, maximum: int, prime: int) -> tuple[int, ...]:
    result = [1]
    for _ in range(maximum):
        result.append(result[-1] * value % prime)
    return tuple(result)


def orbit103_coefficients_mod(
    k_pair: tuple[int, int], r_pair: tuple[int, int], prime: int
) -> tuple[int, int]:
    """Evaluate the separate bihomogenizations of short ``a`` and ``b``."""

    kn, kd = (value % prime for value in k_pair)
    rn, rd = (value % prime for value in r_pair)
    kp = _powers(kn, 6, prime)
    kdp = _powers(kd, 6, prime)
    rp = _powers(rn, 12, prime)
    rdp = _powers(rd, 12, prime)

    # Homogeneous A4 of bidegree (4,4).
    a4 = (
        kp[4] * rp[4]
        + 16 * kn * kd * (kp[2] + 12 * kdp[2]) * rp[3] * rd
        + 8 * kp[2] * kdp[2] * (kp[2] + 152 * kdp[2]) * rp[2] * rdp[2]
        + 64 * kn * kdp[3] * (31 * kp[2] - 12 * kdp[2]) * rn * rdp[3]
        + 16 * kp[2] * kdp[2] * (61 * kp[2] - 48 * kdp[2]) * rdp[4]
    ) % prime
    r2_minus_4 = (rp[2] - 4 * rdp[2]) % prime
    coefficient_a = -27 * r2_minus_4**2 * a4 % prime

    # Homogeneous B6 of bidegree (6,6).
    b6 = (
        kp[6] * rp[6]
        + 24 * kp[3] * kd * (kp[2] + 12 * kdp[2]) * rp[5] * rd
        + 12
        * (
            kp[6]
            + 160 * kp[4] * kdp[2]
            + 192 * kp[2] * kdp[4]
            + 1152 * kdp[6]
        )
        * rp[4]
        * rdp[2]
        + kn
        * kd
        * (3072 * kp[4] + 35072 * kp[2] * kdp[2] + 27648 * kdp[4])
        * rp[3]
        * rdp[3]
        + (
            1488 * kp[6]
            + 91776 * kp[4] * kdp[2]
            + 4608 * kp[2] * kdp[4]
            - 55296 * kdp[6]
        )
        * rp[2]
        * rdp[4]
        + kn
        * kd
        * (85632 * kp[4] - 13824 * kp[2] * kdp[2] - 110592 * kdp[4])
        * rn
        * rdp[5]
        + kp[2]
        * kdp[2]
        * (26560 * kp[4] - 4608 * kp[2] * kdp[2] - 55296 * kdp[4])
        * rdp[6]
    ) % prime
    coefficient_b = 54 * r2_minus_4**3 * b6 % prime
    return coefficient_a, coefficient_b


def trace_and_count(
    coefficient_a: int,
    coefficient_b: int,
    prime: int,
    characters: Sequence[int],
) -> tuple[int, int, bool]:
    discriminant_core = (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime
    if discriminant_core == 0:
        return 0, -1, False
    trace = -sum(
        characters[(x**3 + coefficient_a * x + coefficient_b) % prime]
        for x in range(prime)
    )
    return trace, prime + 1 - trace, True


def local_contribution_units(trace: int, point_count: int, prime: int) -> int:
    return int(round(((2.0 - trace) / point_count) * log(float(prime)) * SCORE_SCALE))


def excluded_k_mod(k_pair: tuple[int, int], prime: int) -> bool:
    kn, kd = (value % prime for value in k_pair)
    if kd == 0:
        return True
    k2 = kn * kn % prime
    d2 = kd * kd % prime
    return (
        kn == 0
        or (k2 - 4 * d2) % prime == 0
        or (3 * k2 - 4 * d2) % prime == 0
        or (3 * k2 + 4 * d2) % prime == 0
    )


def stage1_surface_symbol(prime: int, k_index: int) -> tuple[bool, int, int, int, int]:
    """Return good,singular,placeholder_count,aggregate_trace,score_units."""

    k_pair = projective_value(k_index, prime)
    if excluded_k_mod(k_pair, prime):
        return False, 1, -1, 0, 0
    characters = quadratic_characters(prime)
    contribution = 0
    trace_sum = 0
    good_count = 0
    for r_index in range(prime + 1):
        coefficient_a, coefficient_b = orbit103_coefficients_mod(
            k_pair, projective_value(r_index, prime), prime
        )
        trace, point_count, good = trace_and_count(
            coefficient_a, coefficient_b, prime, characters
        )
        if not good:
            continue
        trace_sum += trace
        contribution += local_contribution_units(trace, point_count, prime)
        good_count += 1
    if good_count < prime - 8:
        return False, 1, -1, 0, 0
    # Averaging keeps the integer scale comparable with ordinary local tables.
    score = int(round(contribution / good_count))
    return True, 0, good_count, trace_sum, score


def stage2_fibre_symbol(
    prime: int, k_pair: tuple[int, int], r_index: int
) -> tuple[bool, int, int, int, int]:
    characters = quadratic_characters(prime)
    coefficient_a, coefficient_b = orbit103_coefficients_mod(
        (k_pair[0] % prime, k_pair[1] % prime),
        projective_value(r_index, prime),
        prime,
    )
    trace, point_count, good = trace_and_count(
        coefficient_a, coefficient_b, prime, characters
    )
    if not good:
        return False, 1, -1, 0, 0
    return True, 0, point_count, trace, local_contribution_units(trace, point_count, prime)


def write_table(
    path: Path,
    *,
    model_tag: str,
    ensembles: Sequence[Sequence[int]],
    symbol_builder,
) -> tuple[tuple[int, ...], ...]:
    lines = [
        "H92_Q12O5867_PROJECTIVE_NAGAO_TABLE_V1",
        f"M {model_tag} 8 12 {SCORE_SCALE}",
        f"C {len(ensembles)}",
    ]
    usable_blocks = []
    serialized_blocks = []
    for block in ensembles:
        serialized_primes = []
        for prime in block:
            symbols = [symbol_builder(prime, index) for index in range(prime + 1)]
            good_scores = [symbol[4] for symbol in symbols if symbol[0]]
            # A rational k can reduce to an exceptional member modulo p.  Such
            # a prime may leave no variance in r and cannot be standardized.
            if len(good_scores) < 2 or min(good_scores) == max(good_scores):
                continue
            prime_lines = [f"P {prime} {prime + 1}"]
            for good, singular, point_count, trace, score in symbols:
                prime_lines.append(
                    f"{int(good)} {int(singular)} {point_count} {trace} {score}"
                )
            serialized_primes.append((prime, prime_lines))
        if not serialized_primes:
            raise ArithmeticError("a requested prime ensemble has no usable tables")
        usable_blocks.append(tuple(prime for prime, _ in serialized_primes))
        serialized_blocks.append(serialized_primes)
    for block_index, serialized_primes in enumerate(serialized_blocks, start=1):
        lines.append(f"B {block_index} {len(serialized_primes)}")
        for _prime, prime_lines in serialized_primes:
            lines.extend(prime_lines)
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return tuple(usable_blocks)


def compile_scanner(local_directory: Path) -> tuple[Path, list[str]]:
    compiler = shutil.which("g++")
    if compiler is None:
        raise FileNotFoundError("g++ is required for the deep projective scan")
    binary = local_directory / "scan-orbit103-nagao"
    command = [
        compiler,
        "-O3",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-pedantic",
        str(CPP_SOURCE),
        "-o",
        str(binary),
    ]
    subprocess.run(command, check=True)
    return binary, command


def run_complete_scan(
    binary: Path,
    table: Path,
    output: Path,
    *,
    height: int,
    finalists: int,
    ensemble_count: int,
) -> tuple[list[dict[str, object]], list[str], float]:
    command = [
        str(binary),
        str(table),
        str(height),
        str(height),
        "1",
        ",".join("1" for _ in range(ensemble_count)),
        str(finalists),
        str(output),
        "1",
        "--rank-region",
        "1",
        "1",
    ]
    started = perf_counter()
    subprocess.run(command, check=True)
    elapsed = perf_counter() - started
    document = json.loads(output.read_text())
    records = [row["score"] for row in document["ranked_prefix"]]
    return records, command, elapsed


def parse_projective_pair(record: dict[str, object]) -> tuple[int, int]:
    pair = record["projective_pair"]
    if not isinstance(pair, list) or len(pair) != 2:
        raise ValueError("scanner finalist has no projective pair")
    return int(pair[0]), int(pair[1])


def weakest_score(record: dict[str, object]) -> float:
    return float(record["worst_block_signal"])


def mean_score(record: dict[str, object]) -> float:
    return float(record["mean_block_signal"])


def exact_fibre_record(k_pair: tuple[int, int], r_pair: tuple[int, int]) -> dict[str, object]:
    """Record exact integral coefficients and the two homogeneous known points."""

    ka, kb = k_pair
    ra, rb = r_pair
    if kb <= 0 or rb <= 0 or gcd(abs(ka), kb) != 1 or gcd(abs(ra), rb) != 1:
        raise ValueError("specialization pairs must be primitive with positive denominator")
    k = Fraction(ka, kb)
    r = Fraction(ra, rb)
    a4 = (
        k**4 * r**4
        + 16 * k * (k**2 + 12) * r**3
        + 8 * k**2 * (k**2 + 152) * r**2
        + 64 * k * (31 * k**2 - 12) * r
        + 16 * k**2 * (61 * k**2 - 48)
    )
    b6 = (
        k**6 * r**6
        + 24 * k**3 * (k**2 + 12) * r**5
        + 12 * (k**6 + 160 * k**4 + 192 * k**2 + 1152) * r**4
        + k * (3072 * k**4 + 35072 * k**2 + 27648) * r**3
        + (1488 * k**6 + 91776 * k**4 + 4608 * k**2 - 55296) * r**2
        + k * (85632 * k**4 - 13824 * k**2 - 110592) * r
        + k**2 * (26560 * k**4 - 4608 * k**2 - 55296)
    )
    coefficient_a = -27 * (r**2 - 4) ** 2 * a4
    coefficient_b = 54 * (r**2 - 4) ** 3 * b6
    scale = kb * rb**2
    integral_a = coefficient_a * scale**4
    integral_b = coefficient_b * scale**6
    if integral_a.denominator != 1 or integral_b.denominator != 1:
        raise ArithmeticError("bihomogeneous orbit-103 scaling stopped being integral")

    x_plus = 3 * (r + 2) * (
        k**2 * r**3
        + (-2 * k**2 - 16 * k + 96) * r**2
        + (-20 * k**2 + 320 * k + 192) * r
        + 232 * k**2
        + 192 * k
    )
    y_plus = 216 * (r + k) * (r + 2) ** 2 * (
        (k - 2) ** 2 * r**2
        + (-12 * k**2 + 32 * k + 48) * r
        + 52 * k**2
        + 80 * k
        + 16
    )
    x_minus = 3 * (r - 2) * (
        k**2 * r**3
        + (2 * k**2 - 16 * k - 96) * r**2
        + (-20 * k**2 - 320 * k + 192) * r
        - 232 * k**2
        + 192 * k
    )
    y_minus = -216 * (r + k) * (r - 2) ** 2 * (
        (k + 2) ** 2 * r**2
        + (12 * k**2 + 32 * k - 48) * r
        + 52 * k**2
        - 80 * k
        + 16
    )
    points = []
    for x_coordinate, y_coordinate in ((x_plus, y_plus), (x_minus, y_minus)):
        x_integral = x_coordinate * scale**2
        y_integral = y_coordinate * scale**3
        if x_integral.denominator != 1 or y_integral.denominator != 1:
            raise ArithmeticError("known orbit-103 point stopped being integral")
        if y_integral**2 != x_integral**3 + integral_a * x_integral + integral_b:
            raise ArithmeticError("known orbit-103 point missed the integral fibre")
        points.append([str(x_integral.numerator), str(y_integral.numerator)])
    discriminant_core = 4 * integral_a**3 + 27 * integral_b**2
    return {
        "short_integral_model": ["0", "0", "0", str(integral_a.numerator), str(integral_b.numerator)],
        "known_points_Q_plus_Q_minus": points,
        "nonsingular": bool(discriminant_core),
        "integral_scale": str(scale),
        "model_bit_height": max(
            abs(integral_a.numerator).bit_length(),
            abs(integral_b.numerator).bit_length(),
        ),
    }


def candidate_sort_key(record: dict[str, object]) -> tuple[object, ...]:
    return (
        -float(record["confirmation_weakest_block"]),
        -float(record["selection_weakest_block"]),
        -float(record["confirmation_mean_block"]),
        int(record["combined_height"]),
        str(record["k"]),
        str(record["r"]),
    )


def rescore_candidate(
    k_pair: tuple[int, int],
    r_pair: tuple[int, int],
    ensembles: Sequence[Sequence[int]],
    table_cache: dict[tuple[tuple[int, int], int], tuple[tuple[int, ...], float, float]],
) -> tuple[list[float], int, int]:
    standardized_blocks: list[float] = []
    good = 0
    bad = 0
    for block in ensembles:
        total = 0.0
        for prime in block:
            k_mod = (k_pair[0] % prime, k_pair[1] % prime)
            if k_mod[1] != 0 and excluded_k_mod(k_mod, prime):
                bad += 1
                continue
            cache_key = (k_pair, prime)
            cached = table_cache.get(cache_key)
            if cached is None:
                scores = []
                good_scores = []
                for projective_index in range(prime + 1):
                    symbol = stage2_fibre_symbol(prime, k_mod, projective_index)
                    score = symbol[4]
                    scores.append(score)
                    if symbol[0]:
                        good_scores.append(score)
                if len(good_scores) < 2:
                    bad += 1
                    continue
                mean = sum(good_scores) / len(good_scores)
                deviation = sqrt(
                    sum((score - mean) ** 2 for score in good_scores)
                    / len(good_scores)
                )
                if not deviation:
                    bad += 1
                    continue
                cached = tuple(scores), mean, deviation
                table_cache[cache_key] = cached
            scores, mean, deviation = cached
            if r_pair[1] % prime == 0:
                index = prime
            else:
                index = r_pair[0] % prime * pow(r_pair[1] % prime, -1, prime) % prime
            symbol = stage2_fibre_symbol(prime, k_mod, index)
            if not symbol[0]:
                bad += 1
                continue
            total += (scores[index] - mean) / deviation
            good += 1
        standardized_blocks.append(total / sqrt(len(block)))
    return standardized_blocks, good, bad


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-height", type=int, default=3000)
    parser.add_argument("--r-height", type=int, default=5000)
    parser.add_argument("--k-finalists", type=int, default=8)
    parser.add_argument("--r-finalists-per-k", type=int, default=100)
    parser.add_argument("--small-r-height", type=int, default=500)
    parser.add_argument("--small-r-finalists", type=int, default=1000)
    parser.add_argument("--finalists", type=int, default=100)
    parser.add_argument("--prime-lower", type=int, default=19)
    parser.add_argument("--prime-upper", type=int, default=199)
    parser.add_argument("--confirmation-prime-lower", type=int, default=211)
    parser.add_argument("--confirmation-prime-upper", type=int, default=401)
    parser.add_argument("--ensemble-count", type=int, default=3)
    parser.add_argument("--local-directory", type=Path, default=DEFAULT_LOCAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if min(
        args.k_height,
        args.r_height,
        args.k_finalists,
        args.r_finalists_per_k,
        args.small_r_height,
        args.small_r_finalists,
        args.finalists,
    ) < 1:
        raise SystemExit("all height and finalist bounds must be positive")
    source = json.loads(SOURCE.read_text())
    if source.get("status") != "PASS_EXACT_RESOLVED_RR_QUARTIC_AND_WEIERSTRASS":
        raise SystemExit("the exact orbit-103 source certificate is missing")

    args.local_directory.mkdir(parents=True, exist_ok=True)
    binary, compile_command = compile_scanner(args.local_directory)
    selection_ensembles = prime_ensembles(
        args.prime_lower, args.prime_upper, args.ensemble_count
    )
    confirmation_ensembles = prime_ensembles(
        args.confirmation_prime_lower,
        args.confirmation_prime_upper,
        args.ensemble_count,
    )
    if set(sum(selection_ensembles, ())) & set(sum(confirmation_ensembles, ())):
        raise AssertionError("selection and confirmation prime windows overlap")

    total_started = perf_counter()
    stage1_table = args.local_directory / "stage1-k-surface-tables.txt"
    table_started = perf_counter()
    stage1_usable_ensembles = write_table(
        stage1_table,
        model_tag=sha256((SOURCE.read_bytes() + b"|stage1-k")).hexdigest(),
        ensembles=selection_ensembles,
        symbol_builder=stage1_surface_symbol,
    )
    stage1_table_seconds = perf_counter() - table_started
    stage1_raw = args.local_directory / "stage1-k-complete-ranking.json"
    stage1_records_raw, stage1_command, stage1_scan_seconds = run_complete_scan(
        binary,
        stage1_table,
        stage1_raw,
        height=args.k_height,
        finalists=3 * args.k_finalists,
        ensemble_count=args.ensemble_count,
    )
    # The equation is fixed by (k,r)->(-k,-r).  Keep k>0 so that the second
    # stage does not spend half its budget on isomorphic duplicate scans.
    stage1_records = [
        record
        for record in stage1_records_raw
        if parse_projective_pair(record)[0] > 0
    ][: args.k_finalists]
    if len(stage1_records) < args.k_finalists:
        raise ArithmeticError("the symmetry-reduced k prefix is unexpectedly short")

    k_rows = []
    fibre_rows = []
    confirmation_table_cache: dict[
        tuple[tuple[int, int], int], tuple[tuple[int, ...], float, float]
    ] = {}
    stage2_commands = []
    stage2_usable_ensembles = {}
    stage2_table_seconds = 0.0
    stage2_scan_seconds = 0.0
    for k_rank, k_record in enumerate(stage1_records, start=1):
        k_pair = parse_projective_pair(k_record)
        if k_pair[1] == 0:
            continue
        k_text = f"{k_pair[0]}/{k_pair[1]}"
        safe = k_text.replace("-", "m").replace("/", "_")
        table_path = args.local_directory / f"stage2-r-k-{safe}-tables.txt"
        table_started = perf_counter()
        usable = write_table(
            table_path,
            model_tag=sha256((SOURCE.read_bytes() + f"|stage2|{k_text}".encode())).hexdigest(),
            ensembles=selection_ensembles,
            symbol_builder=lambda prime, index, pair=k_pair: stage2_fibre_symbol(
                prime, pair, index
            ),
        )
        stage2_usable_ensembles[k_text] = [list(block) for block in usable]
        stage2_table_seconds += perf_counter() - table_started
        scan_path = args.local_directory / f"stage2-r-k-{safe}-complete-ranking.json"
        r_records, command, seconds = run_complete_scan(
            binary,
            table_path,
            scan_path,
            height=args.r_height,
            finalists=args.r_finalists_per_k,
            ensemble_count=args.ensemble_count,
        )
        stage2_commands.append(command)
        stage2_scan_seconds += seconds
        tagged_r_records = [("deep_trace", record) for record in r_records]
        if k_pair == (1, 1):
            small_scan_path = args.local_directory / "stage2-r-k-1_1-small-model-ranking.json"
            small_records, small_command, small_seconds = run_complete_scan(
                binary,
                table_path,
                small_scan_path,
                height=args.small_r_height,
                finalists=args.small_r_finalists,
                ensemble_count=args.ensemble_count,
            )
            stage2_commands.append(small_command)
            stage2_scan_seconds += small_seconds
            tagged_r_records.extend(("small_model", record) for record in small_records)
        k_rows.append(
            {
                "symmetry_reduced_rank": k_rank,
                "parameter": k_text,
                "projective_pair": list(k_pair),
                "weakest_standardized_block": weakest_score(k_record),
                "mean_standardized_block": mean_score(k_record),
                "complete_scan_record": k_record,
            }
        )
        seen_r_pairs = set()
        for r_rank, (lane, r_record) in enumerate(tagged_r_records, start=1):
            r_pair = parse_projective_pair(r_record)
            if r_pair[1] == 0:
                continue
            if r_pair in seen_r_pairs:
                continue
            seen_r_pairs.add(r_pair)
            confirmation_blocks, good, bad = rescore_candidate(
                k_pair, r_pair, confirmation_ensembles, confirmation_table_cache
            )
            exact = exact_fibre_record(k_pair, r_pair)
            if not exact["nonsingular"]:
                continue
            fibre_rows.append(
                {
                    "k": k_text,
                    "k_projective_pair": list(k_pair),
                    "k_symmetry_reduced_rank": k_rank,
                    "r": f"{r_pair[0]}/{r_pair[1]}",
                    "r_projective_pair": list(r_pair),
                    "r_population_rank_within_k": r_rank,
                    "stage2_lane": lane,
                    "selection_weakest_block": weakest_score(r_record),
                    "selection_mean_block": mean_score(r_record),
                    "confirmation_standardized_blocks": confirmation_blocks,
                    "confirmation_weakest_block": min(confirmation_blocks),
                    "confirmation_mean_block": sum(confirmation_blocks) / len(confirmation_blocks),
                    "confirmation_good_primes": good,
                    "confirmation_bad_primes": bad,
                    "combined_height": max(abs(k_pair[0]), k_pair[1], abs(r_pair[0]), r_pair[1]),
                    "exact_specialization": exact,
                }
            )
    fibre_rows.sort(key=candidate_sort_key)
    finalists = fibre_rows[: args.finalists]
    small_coefficient_finalists = sorted(
        (row for row in fibre_rows if row["stage2_lane"] == "small_model"),
        key=lambda row: (
            row["exact_specialization"]["model_bit_height"],
            -row["confirmation_weakest_block"],
            -row["selection_weakest_block"],
        ),
    )[: args.small_r_finalists]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "elkies-k3.e6a1-orbit103-specialization-search.v1",
        "status": "PASS_BOUNDED_TWO_STAGE_NAGAO_SEARCH",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": sha256(SOURCE.read_bytes()).hexdigest(),
            "generic_arithmetic_rank": 2,
            "known_rational_directions": ["Q_plus", "Q_minus"],
            "excluded_geometric_direction": "anti-invariant over QQ(sqrt(-3))",
        },
        "method": {
            "stage1": (
                "For each k mod p, average the exact repository Nagao contribution "
                "over good r in P1(F_p), then center and population-standardize in k. "
                "The involution (k,r)->(-k,-r) is quotiented by retaining k>0."
            ),
            "stage2": (
                "For each retained rational k, score every rational r in the declared "
                "height box by centered/standardized exact local traces."
            ),
            "ranking": "weakest of three or more pairwise-disjoint prime ensembles",
            "bad_reduction_policy": "mean imputation in standardized scans",
            "confirmation": (
                "Independent larger primes rescore only the retained (k,r) pairs; "
                "each local score is centered and population-standardized over "
                "P1(F_p) before weakest-block ranking."
            ),
        },
        "bounds": {
            "k_projective_height": args.k_height,
            "r_projective_height": args.r_height,
            "k_finalists": args.k_finalists,
            "r_finalists_per_k": args.r_finalists_per_k,
            "small_model_k": "1/1",
            "small_model_r_projective_height": args.small_r_height,
            "small_model_finalists": args.small_r_finalists,
            "reported_fibre_finalists": args.finalists,
            "selection_prime_ensembles": [list(block) for block in selection_ensembles],
            "stage1_usable_prime_ensembles": [
                list(block) for block in stage1_usable_ensembles
            ],
            "stage2_usable_prime_ensembles_by_k": stage2_usable_ensembles,
            "confirmation_prime_ensembles": [list(block) for block in confirmation_ensembles],
        },
        "k_surface_candidates": k_rows,
        "fibre_candidates_scored": len(fibre_rows),
        "finalists": finalists,
        "small_coefficient_finalists": small_coefficient_finalists,
        "commands": {
            "compile": compile_command,
            "stage1": stage1_command,
            "stage2": stage2_commands,
            "driver": sys.argv,
        },
        "runtime_seconds": {
            "stage1_table": stage1_table_seconds,
            "stage1_scan": stage1_scan_seconds,
            "stage2_tables": stage2_table_seconds,
            "stage2_scans": stage2_scan_seconds,
            "total": perf_counter() - total_started,
        },
        "proof_boundary": (
            "All finite-field point counts, rational boxes, exact fibre coefficients, "
            "and Q_plus/Q_minus specializations are reproducible. Nagao scores are "
            "ranking heuristics only. No finalist has a certified rank jump until an "
            "additional rational point passes a combined exact finite-quotient test "
            "against Q_plus,Q_minus."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "E6A1O103SEARCH|"
        f"k_population={json.loads(stage1_raw.read_text())['population_count']}|"
        f"k_retained={len(k_rows)}|fibres={len(fibre_rows)}|"
        f"finalists={len(finalists)}|seconds={payload['runtime_seconds']['total']:.3f}|"
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
