#!/usr/bin/env python3
"""Bounded exact screen of four split-infinity six-root Mestre families.

For each family this program reconstructs the primitive Mestre quartic from
its roots.  If its leading coefficient is ``m^2*(T^2+C)``, the base change

``T=(C-u^2)/(2u)``

splits quartic infinity with square root ``m*(C+u^2)/(2u)``.  The closed
search population is every positive reduced ``u=a/b`` with ``a<=512`` and
``b<=32``.  The C++ companion scores that population at all usable primes
11..251.  Per family, 48 score leaders and 16 smallest-T-height anchors get
exact conductor calls.  Every completed curve below log(N)=182.72 gets a
raw/+T/-T ratpoints search at H=2,000,000 and denominator at most 13,000,
followed by an exact mod-3 finite-reduction certificate.

The local score and bounded point search are heuristics.  A negative result
is not a rank upper bound, and the discriminant is not used as a conductor.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, isqrt, lcm
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
from typing import Any, Iterable, Sequence

from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    capped_minimal_curve_data,
    point_digest,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
    run_capped_process,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import local_trace
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
TARGET_LOG_CONDUCTOR = Decimal("182.72")
NUMERATOR_BOUND = 512
DENOMINATOR_BOUND = 32
SCORE_KEEP = 256
SCORE_SELECTION_PER_FAMILY = 48
HEIGHT_SELECTION_PER_FAMILY = 16
POINT_HEIGHT = 2_000_000
POINT_DENOMINATOR_BOUND = 13_000
RELATION_PRIME_BOUND = 499
STACK_BYTES = 512_000_000
DEFAULT_OUTPUT_DIRECTORY = Path(
    "artifacts/local/elliptic-curves/mestre-dsquare-four-v1"
)
RATPOINTS = ROOT / "tmp/ratpoints/root/usr/bin/ratpoints"
RATPOINTS_LIBRARY = ROOT / "tmp/ratpoints/root/usr/lib/x86_64-linux-gnu"
POINT_PATTERN = re.compile(r"\((-?\d+) : (\d+)\)")


@dataclass(frozen=True)
class Family:
    index: int
    label: str
    roots: tuple[int, ...]
    base_constant: Fraction
    leading_square_multiplier: int

    @property
    def construction(self) -> SixRootMestreConstruction:
        return SixRootMestreConstruction(tuple(Q(root) for root in self.roots))


FAMILIES = (
    Family(0, "r0_7_225_232_235_265", (0, 7, 225, 232, 235, 265), Q(84_878, 3), 3),
    Family(1, "r0_9_213_247_256_291", (0, 9, 213, 247, 256, 291), Q(196_250, 3), 3),
    Family(2, "r0_25_95_143_168_205", (0, 25, 95, 143, 168, 205), Q(39_146), 1),
    Family(3, "r0_43_128_197_231_289", (0, 43, 128, 197, 231, 289), Q(55_950), 3),
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def gp_version() -> str:
    completed = subprocess.run(
        ["gp", "--version"], capture_output=True, text=True, check=True
    )
    lines = [
        line.strip()
        for line in (completed.stdout + completed.stderr).splitlines()
        if line.strip()
    ]
    return lines[0] if lines else "gp version unavailable"


def primes(lower: int, upper: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (upper + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(upper) + 1):
        if sieve[prime]:
            sieve[prime * prime : upper + 1 : prime] = b"\x00" * (
                (upper - prime * prime) // prime + 1
            )
    return tuple(index for index in range(lower, upper + 1) if sieve[index])


SCORE_PRIMES = primes(11, 251)


def base_parameter(family: Family, parameter_u: Fraction) -> Fraction:
    parameter_u = Q(parameter_u)
    if not parameter_u:
        raise ValueError("the split-infinity base change has a pole at u=0")
    return (family.base_constant - parameter_u**2) / (2 * parameter_u)


def leading_square(family: Family, parameter_u: Fraction) -> Fraction:
    parameter_u = Q(parameter_u)
    if not parameter_u:
        raise ValueError("the split-infinity base change has a pole at u=0")
    return (
        family.leading_square_multiplier
        * (family.base_constant + parameter_u**2)
        / (2 * parameter_u)
    )


def split_infinity_point(
    family: Family, parameter_u: Fraction
) -> tuple[Fraction, Fraction]:
    construction = family.construction
    parameter_t = base_parameter(family, parameter_u)
    _, d, c, b, a = construction.primitive_quartic_coefficients(parameter_t)
    square_root = leading_square(family, parameter_u)
    if square_root**2 != a:
        raise AssertionError("the claimed base change did not split infinity")
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    point = (
        36 * g0 / a,
        54 * (a * g1 - b * g0) / square_root**3,
    )
    coefficient_a = construction.primitive_jacobian_coefficients(parameter_t)[3]
    coefficient_b = construction.primitive_jacobian_coefficients(parameter_t)[4]
    if point[1] ** 2 != point[0] ** 3 + coefficient_a * point[0] + coefficient_b:
        raise AssertionError("split infinity missed the exact Jacobian")
    return point


def known_jacobian_points(
    family: Family, parameter_u: Fraction
) -> tuple[tuple[Fraction, Fraction], ...]:
    construction = family.construction
    parameter_t = base_parameter(family, parameter_u)
    affine = tuple(
        quartic_point_to_jacobian(construction, parameter_t, point)
        for point in primitive_visible_points(construction, parameter_t)
    )
    return affine + (split_infinity_point(family, parameter_u),)


def family_geometry(family: Family) -> dict[str, Any]:
    construction = family.construction
    if construction.quartic_condition != 0:
        raise AssertionError("a declared family missed the Mestre obstruction")
    for parameter_t in (Q(1), Q(2), Q(-7, 3), Q(11, 5)):
        leading = construction.primitive_quartic_coefficients(parameter_t)[4]
        expected = family.leading_square_multiplier**2 * (
            parameter_t**2 + family.base_constant
        )
        if leading != expected:
            raise AssertionError("a declared leading coefficient changed")
    if (
        len(construction.primitive_discriminant_polynomial) != 21
        or any(construction.primitive_discriminant_polynomial[1::2])
    ):
        raise AssertionError("a family lost its even degree-20 discriminant")
    parameter_u = Q(1)
    parameter_t = base_parameter(family, parameter_u)
    points = known_jacobian_points(family, parameter_u)
    certificate = mod3_independence_certificate(
        construction.primitive_jacobian_coefficients(parameter_t),
        points,
        prime_bound=251,
    )
    if certificate["combined_exact_rank_over_F3"] != 12:
        raise AssertionError("the visible-plus-infinity generic baseline changed")
    return {
        "index": family.index,
        "label": family.label,
        "roots": list(family.roots),
        "quartic_condition": "0",
        "quartic_fixed_square_content": rational_text(construction.quartic_content),
        "quartic_fixed_square_scale": rational_text(construction.quartic_square_scale),
        "leading_coefficient": (
            f"{family.leading_square_multiplier**2}*(T^2+"
            f"{rational_text(family.base_constant)})"
        ),
        "base_change": (
            f"T=(({rational_text(family.base_constant)})-u^2)/(2u)"
        ),
        "leading_square_root": (
            f"{family.leading_square_multiplier}*"
            f"(({rational_text(family.base_constant)})+u^2)/(2u)"
        ),
        "primitive_discriminant_degree_in_T": 20,
        "primitive_discriminant_even": True,
        "u1_visible_plus_infinity_point_count": len(points),
        "u1_exact_mod3_rank_lower_bound": certificate[
            "combined_exact_rank_over_F3"
        ],
        "u1_certificate_primes": certificate["certificate_primes"],
        "u1_point_sha256": certificate["point_sha256"],
    }


def compile_and_run_scorer(
    source: Path, output: Path, *, threads: int
) -> dict[str, Any]:
    compiler = shutil.which("c++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="mestre-dsquare-four-") as directory:
        binary = Path(directory) / "score"
        compile_command = (
            compiler,
            "-std=c++17",
            "-O3",
            "-DNDEBUG",
            "-fopenmp",
            str(source),
            "-o",
            str(binary),
        )
        run_capped_process(compile_command, timeout=60)
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = str(threads)
        environment["OMP_DYNAMIC"] = "FALSE"
        completed = subprocess.run(
            [str(binary), str(NUMERATOR_BOUND), str(DENOMINATOR_BOUND), str(SCORE_KEEP)],
            capture_output=True,
            text=True,
            env=environment,
            timeout=600,
            check=True,
        )
    if completed.stderr.strip():
        raise RuntimeError(f"scorer wrote stderr: {completed.stderr.strip()}")
    output.write_text(completed.stdout)
    return {
        "compiler": subprocess.run(
            [compiler, "--version"], capture_output=True, text=True, check=True
        ).stdout.splitlines()[0],
        "source_sha256": sha256_file(source),
        "output_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "wall_seconds": round(time.monotonic() - started, 6),
        "threads": threads,
    }


def parse_score_output(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_DSQUARE_FOUR_SCORE_V1" or lines[-1] != "DONE":
        raise AssertionError("malformed four-family score output")
    observed_primes = tuple(map(int, lines[1].split()[1:]))
    if observed_primes != SCORE_PRIMES:
        raise AssertionError("the score-prime band changed")
    domain = tuple(map(int, lines[2].split()[1:]))
    if domain[:2] != (NUMERATOR_BOUND, DENOMINATOR_BOUND):
        raise AssertionError("the score domain changed")
    expected_population = sum(
        gcd(numerator, denominator) == 1
        for denominator in range(1, DENOMINATOR_BOUND + 1)
        for numerator in range(1, NUMERATOR_BOUND + 1)
    )
    if domain[2] != expected_population:
        raise AssertionError("the reduced score population changed")
    records = []
    rank_by_family = Counter()
    for line in lines[3:-1]:
        fields = line.split()
        if fields[0] != "C" or len(fields) != 7:
            raise AssertionError("malformed score candidate")
        family_index, score_rank, numerator, denominator, good = (
            int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4]), int(fields[6])
        )
        rank_by_family[family_index] += 1
        if score_rank != rank_by_family[family_index]:
            raise AssertionError("score ranks ceased to be consecutive")
        records.append(
            {
                "family_index": family_index,
                "score_rank": score_rank,
                "numerator": numerator,
                "denominator": denominator,
                "u": rational_text(Q(numerator, denominator)),
                "local_score": fields[5],
                "good_prime_coverage": good,
            }
        )
    if rank_by_family != Counter({index: SCORE_KEEP for index in range(4)}):
        raise AssertionError("the scorer retained the wrong family counts")
    return records, {
        "primes": list(observed_primes),
        "numerator_bound": NUMERATOR_BOUND,
        "denominator_bound": DENOMINATOR_BOUND,
        "primitive_positive_rational_population_per_family": expected_population,
        "retained_per_family": SCORE_KEEP,
    }


def replay_score(record: dict[str, Any]) -> None:
    family = FAMILIES[int(record["family_index"])]
    parameter_u = Q(record["numerator"], record["denominator"])
    parameter_t = base_parameter(family, parameter_u)
    coefficients = family.construction.primitive_jacobian_coefficients(parameter_t)
    score = Q(0)
    good = 0
    for prime in SCORE_PRIMES:
        trace = local_trace(coefficients, prime)
        if trace is None:
            continue
        score += Q(2 - trace, prime + 1 - trace)
        good += 1
    if good != record["good_prime_coverage"]:
        raise AssertionError("Python score replay changed good-prime coverage")
    if abs(float(score) - float(record["local_score"])) > 5e-13:
        raise AssertionError("Python score replay disagreed with C++")


def t_height_record(family: Family, numerator: int, denominator: int) -> dict[str, Any]:
    parameter_u = Q(numerator, denominator)
    parameter_t = base_parameter(family, parameter_u)
    return {
        "family_index": family.index,
        "score_rank": None,
        "numerator": numerator,
        "denominator": denominator,
        "u": rational_text(parameter_u),
        "local_score": None,
        "good_prime_coverage": None,
        "base_T": rational_text(parameter_t),
        "T_height": max(abs(parameter_t.numerator), parameter_t.denominator),
    }


def select_candidates(score_records: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: dict[tuple[int, int, int], dict[str, Any]] = {}
    for family in FAMILIES:
        leaders = [
            dict(record)
            for record in score_records
            if record["family_index"] == family.index
            and record["score_rank"] <= SCORE_SELECTION_PER_FAMILY
        ]
        for record in leaders:
            parameter_t = base_parameter(
                family, Q(record["numerator"], record["denominator"])
            )
            record["base_T"] = rational_text(parameter_t)
            record["T_height"] = max(abs(parameter_t.numerator), parameter_t.denominator)
            record["selection_strata"] = ["top-local-score"]
            selected[(family.index, record["numerator"], record["denominator"])] = record
        anchors = sorted(
            (
                t_height_record(family, numerator, denominator)
                for denominator in range(1, DENOMINATOR_BOUND + 1)
                for numerator in range(1, NUMERATOR_BOUND + 1)
                if gcd(numerator, denominator) == 1
            ),
            key=lambda row: (
                row["T_height"],
                max(row["numerator"], row["denominator"]),
                row["denominator"],
                row["numerator"],
            ),
        )[:HEIGHT_SELECTION_PER_FAMILY]
        score_lookup = {
            (record["numerator"], record["denominator"]): record
            for record in score_records
            if record["family_index"] == family.index
        }
        for anchor in anchors:
            key = (family.index, anchor["numerator"], anchor["denominator"])
            if key in selected:
                selected[key]["selection_strata"].append("smallest-T-height")
                continue
            scored = score_lookup.get((anchor["numerator"], anchor["denominator"]))
            if scored is not None:
                anchor.update(
                    score_rank=scored["score_rank"],
                    local_score=scored["local_score"],
                    good_prime_coverage=scored["good_prime_coverage"],
                )
            anchor["selection_strata"] = ["smallest-T-height"]
            selected[key] = anchor
    result = sorted(
        selected.values(),
        key=lambda row: (
            row["family_index"],
            "top-local-score" not in row["selection_strata"],
            row["score_rank"] if row["score_rank"] is not None else 10**9,
            row["T_height"],
            row["numerator"],
            row["denominator"],
        ),
    )
    for index, record in enumerate(result, start=1):
        record["search_index"] = index
    return result


def candidate_identifier(candidate: dict[str, Any]) -> str:
    return (
        f"f{candidate['family_index']}_u{candidate['numerator']}_"
        f"{candidate['denominator']}"
    )


def conductor_worker(candidate: dict[str, Any], timeout: float) -> dict[str, Any]:
    family = FAMILIES[int(candidate["family_index"])]
    parameter_u = Q(candidate["numerator"], candidate["denominator"])
    parameter_t = base_parameter(family, parameter_u)
    coefficients = family.construction.primitive_jacobian_coefficients(parameter_t)
    global_data = capped_minimal_curve_data(
        coefficients, timeout=timeout, stack_bytes=STACK_BYTES
    )
    return {
        **candidate,
        "status": "completed exact PARI minimal-model/conductor computation",
        "global_curve": global_data,
        "below_strict_log_conductor_182_72": (
            Decimal(str(global_data["log_conductor"])) < TARGET_LOG_CONDUCTOR
        ),
    }


def translated_quartic(
    coefficients: Sequence[Fraction], offset: Fraction
) -> tuple[Fraction, ...]:
    e, d, c, b, a = map(Q, coefficients)
    offset = Q(offset)
    return (
        e + offset * d + offset**2 * c + offset**3 * b + offset**4 * a,
        d + 2 * offset * c + 3 * offset**2 * b + 4 * offset**3 * a,
        c + 3 * offset * b + 6 * offset**2 * a,
        b + 4 * offset * a,
        a,
    )


def integral_ratpoints_coefficients(
    coefficients: Sequence[Fraction],
) -> tuple[tuple[int, ...], int]:
    denominator = 1
    for value in coefficients:
        denominator = lcm(denominator, Q(value).denominator)
    integral = tuple(int(Q(value) * denominator**2) for value in coefficients)
    if any(
        Q(value) * denominator**2 != integer
        for value, integer in zip(coefficients, integral)
    ):
        raise AssertionError("ratpoints square scaling failed")
    return integral, denominator


def rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def run_ratpoints_chart(
    family: Family,
    parameter_t: Fraction,
    chart: str,
    raw_directory: Path,
    timeout: float,
) -> dict[str, Any]:
    offsets = {"raw": Q(0), "plus-T": parameter_t, "minus-T": -parameter_t}
    offset = offsets[chart]
    quartic = family.construction.primitive_quartic_coefficients(parameter_t)
    translated = translated_quartic(quartic, offset)
    integral, ordinate_scale = integral_ratpoints_coefficients(translated)
    path = raw_directory / f"{chart}.out"
    started = time.monotonic()
    if path.exists():
        output = path.read_text()
        cached = True
    else:
        environment = os.environ.copy()
        environment["LD_LIBRARY_PATH"] = str(RATPOINTS_LIBRARY)
        completed = subprocess.run(
            [
                str(RATPOINTS),
                " ".join(map(str, integral)),
                str(POINT_HEIGHT),
                "-du",
                str(POINT_DENOMINATOR_BOUND),
                "-q",
                "-y",
            ],
            capture_output=True,
            text=True,
            env=environment,
            timeout=timeout,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"ratpoints {chart} exited {completed.returncode}: {completed.stderr[:500]}"
            )
        output = completed.stdout
        path.write_text(output)
        cached = False
    finite_x = []
    infinity_count = 0
    for line in output.splitlines():
        match = POINT_PATTERN.fullmatch(line.strip())
        if match is None:
            raise AssertionError(f"unexpected ratpoints output line {line!r}")
        numerator, denominator = map(int, match.groups())
        if denominator == 0:
            infinity_count += 1
        else:
            finite_x.append(Q(numerator, denominator) + offset)
    return {
        "chart": chart,
        "offset": rational_text(offset),
        "integral_coefficients": list(integral),
        "ordinate_square_scale": ordinate_scale,
        "finite_abscissae": finite_x,
        "finite_abscissa_count": len(finite_x),
        "projective_infinity_count": infinity_count,
        "raw_output_sha256": hashlib.sha256(output.encode()).hexdigest(),
        "cached": cached,
        "wall_seconds": round(time.monotonic() - started, 6),
    }


def point_worker(
    candidate: dict[str, Any], raw_root_text: str, timeout: float
) -> dict[str, Any]:
    family = FAMILIES[int(candidate["family_index"])]
    parameter_u = Q(candidate["numerator"], candidate["denominator"])
    parameter_t = base_parameter(family, parameter_u)
    construction = family.construction
    coefficients = construction.primitive_jacobian_coefficients(parameter_t)
    raw_directory = Path(raw_root_text) / candidate_identifier(candidate)
    raw_directory.mkdir(parents=True, exist_ok=True)
    charts = tuple(
        run_ratpoints_chart(family, parameter_t, chart, raw_directory, timeout)
        for chart in ("raw", "plus-T", "minus-T")
    )
    quartic = construction.primitive_quartic_coefficients(parameter_t)
    searched_x = sorted(
        {x_value for chart in charts for x_value in chart["finite_abscissae"]}
    )
    known = known_jacobian_points(family, parameter_u)
    by_jacobian_x: dict[Fraction, tuple[tuple[Fraction, Fraction], str]] = {}
    for index, point in enumerate(known):
        by_jacobian_x.setdefault(point[0], (point, f"known-{index + 1:02d}"))
    known_column_count = len(by_jacobian_x)
    searched_records = []
    for index, x_value in enumerate(searched_x):
        square_root = rational_square_root(quartic_value(quartic, x_value))
        if square_root is None:
            raise AssertionError("a ratpoints abscissa failed the original quartic")
        if square_root == 0:
            continue
        jacobian_point = quartic_point_to_jacobian(
            construction, parameter_t, (x_value, square_root)
        )
        source = f"searched-{index + 1:04d}"
        novel = jacobian_point[0] not in by_jacobian_x
        by_jacobian_x.setdefault(jacobian_point[0], (jacobian_point, source))
        searched_records.append(
            {
                "quartic_x": rational_text(x_value),
                "quartic_y_positive": rational_text(square_root),
                "jacobian_x": rational_text(jacobian_point[0]),
                "novel_jacobian_x": novel,
            }
        )
    pool = tuple(item[0] for item in by_jacobian_x.values())
    sources = tuple(item[1] for item in by_jacobian_x.values())
    known_certificate = mod3_independence_certificate(
        coefficients, pool[:known_column_count], prime_bound=RELATION_PRIME_BOUND
    )
    full_certificate = mod3_independence_certificate(
        coefficients, pool, prime_bound=RELATION_PRIME_BOUND
    )
    return {
        **candidate,
        "status": "completed three-chart bounded search and exact mod-3 certification",
        "point_search": {
            "height_bound": POINT_HEIGHT,
            "denominator_bound": POINT_DENOMINATOR_BOUND,
            "charts": [
                {key: value for key, value in chart.items() if key != "finite_abscissae"}
                for chart in charts
            ],
            "finite_abscissa_union_count": len(searched_x),
            "searched_points": searched_records,
            "novel_jacobian_x_count": sum(
                record["novel_jacobian_x"] for record in searched_records
            ),
        },
        "known_column_count_modulo_inverse": known_column_count,
        "known_exact_rank_lower_bound": known_certificate[
            "combined_exact_rank_over_F3"
        ],
        "pool_point_count_modulo_inverse": len(pool),
        "pool_sources": list(sources),
        "pool_point_sha256": point_digest(pool),
        "exact_specialization_rank_lower_bound": full_certificate[
            "combined_exact_rank_over_F3"
        ],
        "finite_reduction_certificate": full_certificate,
        "bounded_search_is_not_a_rank_upper_bound": True,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--ratpoints-timeout", type=float, default=600.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 12:
        raise SystemExit("--workers must lie in [1,12]")
    output = ROOT / args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    score_path = output / "score-top256.tsv"
    scorer_source = Path(__file__).with_name("score_mestre_dsquare_four.cpp")
    if score_path.exists():
        scorer_run = {
            "source_sha256": sha256_file(scorer_source),
            "output_sha256": sha256_file(score_path),
            "cached": True,
            "threads": args.workers,
        }
    else:
        scorer_run = compile_and_run_scorer(
            scorer_source, score_path, threads=args.workers
        )
        scorer_run["cached"] = False
    score_records, score_scope = parse_score_output(score_path)
    for family in FAMILIES:
        replay_score(
            next(record for record in score_records if record["family_index"] == family.index)
        )
    candidates = select_candidates(score_records)
    input_payload = {
        "schema_version": 1,
        "scope": {
            **score_scope,
            "family_count": len(FAMILIES),
            "scored_family_parameter_pairs": (
                score_scope["primitive_positive_rational_population_per_family"]
                * len(FAMILIES)
            ),
            "score_selection_per_family": SCORE_SELECTION_PER_FAMILY,
            "smallest_T_height_selection_per_family": HEIGHT_SELECTION_PER_FAMILY,
            "selection_uses_no_conductor_or_point_data": True,
        },
        "scorer": scorer_run,
        "families": [family_geometry(family) for family in FAMILIES],
        "candidate_count": len(candidates),
        "candidates": candidates,
    }
    input_payload["result_sha256"] = canonical_digest(input_payload)
    atomic_json(output / "candidate-input.json", input_payload)
    print(
        f"SCORE population={score_scope['primitive_positive_rational_population_per_family']}x4 "
        f"selected={len(candidates)}",
        flush=True,
    )

    conductor_directory = output / "conductor-records"
    conductor_directory.mkdir(exist_ok=True)
    conductor_results = []
    pending = []
    for candidate in candidates:
        path = conductor_directory / f"{candidate_identifier(candidate)}.json"
        if path.exists():
            conductor_results.append(json.loads(path.read_text()))
        else:
            pending.append(candidate)
    print(f"CONDUCTOR cached={len(conductor_results)} pending={len(pending)}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(conductor_worker, candidate, args.conductor_timeout): candidate
            for candidate in pending
        }
        for index, future in enumerate(as_completed(futures), start=1):
            candidate = futures[future]
            try:
                record = future.result()
            except CappedProcessTimeout as error:
                record = {**candidate, "status": "timeout", "error": str(error)}
            except Exception as error:
                record = {**candidate, "status": "error", "error": repr(error)}
            record["result_sha256"] = canonical_digest(record)
            atomic_json(
                conductor_directory / f"{candidate_identifier(candidate)}.json", record
            )
            conductor_results.append(record)
            print(
                f"CONDUCTOR {index}/{len(pending)} {candidate_identifier(candidate)} "
                f"lnN={record.get('global_curve', {}).get('log_conductor')} "
                f"status={record['status']}",
                flush=True,
            )
    completed_conductors = [
        record
        for record in conductor_results
        if record["status"].startswith("completed")
    ]
    qualified = sorted(
        (
            record
            for record in completed_conductors
            if record["below_strict_log_conductor_182_72"]
        ),
        key=lambda record: (
            record["family_index"],
            -float(record["local_score"] or "-1e99"),
            Decimal(str(record["global_curve"]["log_conductor"])),
            record["T_height"],
        ),
    )
    point_input = {
        "scope": "every completed exact conductor below log(N)=182.72",
        "height_bound": POINT_HEIGHT,
        "denominator_bound": POINT_DENOMINATOR_BOUND,
        "chart_set": ["raw", "x=T+z", "x=-T+z"],
        "candidate_count": len(qualified),
        "candidates": qualified,
    }
    point_input["result_sha256"] = canonical_digest(point_input)
    atomic_json(output / "point-search-input.json", point_input)

    point_directory = output / "point-certificates"
    raw_root = output / "ratpoints-raw"
    point_directory.mkdir(exist_ok=True)
    raw_root.mkdir(exist_ok=True)
    point_results = []
    pending_points = []
    retried_point_errors = []
    repaired_legacy_point_digests = []
    for candidate in qualified:
        path = point_directory / f"{candidate_identifier(candidate)}.json"
        if path.exists():
            cached_record = json.loads(path.read_text())
            if cached_record["status"].startswith("completed"):
                declared_digest = cached_record.pop("result_sha256")
                replayed_digest = canonical_digest(cached_record)
                cached_record["result_sha256"] = replayed_digest
                if declared_digest != replayed_digest:
                    repaired_legacy_point_digests.append(
                        {
                            "candidate": candidate_identifier(candidate),
                            "prior_file_sha256": sha256_file(path),
                            "prior_declared_result_sha256": declared_digest,
                            "replayed_result_sha256": replayed_digest,
                            "reason": (
                                "legacy point worker inherited and then overwrote the "
                                "conductor record result_sha256"
                            ),
                        }
                    )
                    atomic_json(path, cached_record)
                point_results.append(cached_record)
            else:
                retried_point_errors.append(
                    {
                        "candidate": candidate_identifier(candidate),
                        "prior_status": cached_record["status"],
                        "prior_error": cached_record.get("error"),
                        "prior_file_sha256": sha256_file(path),
                    }
                )
                pending_points.append(candidate)
        else:
            pending_points.append(candidate)
    print(f"POINT cached={len(point_results)} pending={len(pending_points)}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                point_worker, candidate, str(raw_root), args.ratpoints_timeout
            ): candidate
            for candidate in pending_points
        }
        for index, future in enumerate(as_completed(futures), start=1):
            candidate = futures[future]
            try:
                record = future.result()
            except subprocess.TimeoutExpired as error:
                record = {**candidate, "status": "timeout", "error": repr(error)}
            except Exception as error:
                record = {**candidate, "status": "error", "error": repr(error)}
            record.pop("result_sha256", None)
            record["result_sha256"] = canonical_digest(record)
            atomic_json(
                point_directory / f"{candidate_identifier(candidate)}.json", record
            )
            point_results.append(record)
            rank = record.get("exact_specialization_rank_lower_bound")
            print(
                f"POINT {index}/{len(pending_points)} {candidate_identifier(candidate)} "
                f"rank={rank} status={record['status']}",
                flush=True,
            )
            if rank is not None and rank >= 17:
                print(
                    f"ALERT exact LB={rank} {candidate_identifier(candidate)} "
                    f"lnN={candidate['global_curve']['log_conductor']}",
                    flush=True,
                )

    completed_points = [
        record for record in point_results if record["status"].startswith("completed")
    ]
    completed_points.sort(
        key=lambda record: (
            -record["exact_specialization_rank_lower_bound"],
            Decimal(str(record["global_curve"]["log_conductor"])),
            record["family_index"],
            record["numerator"],
            record["denominator"],
        )
    )
    promotions = [
        record
        for record in completed_points
        if record["exact_specialization_rank_lower_bound"] >= 17
    ]

    def compact_frontier_record(record: dict[str, Any]) -> dict[str, Any]:
        certificate_path = point_directory / f"{candidate_identifier(record)}.json"
        return {
            "family_index": record["family_index"],
            "u": record["u"],
            "T": record["base_T"],
            "exact_rank_lower_bound": record["exact_specialization_rank_lower_bound"],
            "conductor": record["global_curve"]["conductor"],
            "log_conductor": record["global_curve"]["log_conductor"],
            "root_number": record["global_curve"]["root_number"],
            "point_certificate": str(certificate_path.relative_to(ROOT)),
            "point_certificate_sha256": sha256_file(certificate_path),
        }

    pareto_records = sorted(
        (
            record
            for record in completed_points
            if not any(
                other["exact_specialization_rank_lower_bound"]
                >= record["exact_specialization_rank_lower_bound"]
                and Decimal(str(other["global_curve"]["log_conductor"]))
                <= Decimal(str(record["global_curve"]["log_conductor"]))
                and (
                    other["exact_specialization_rank_lower_bound"]
                    > record["exact_specialization_rank_lower_bound"]
                    or Decimal(str(other["global_curve"]["log_conductor"]))
                    < Decimal(str(record["global_curve"]["log_conductor"]))
                )
                for other in completed_points
            )
        ),
        key=lambda record: (
            Decimal(str(record["global_curve"]["log_conductor"])),
            record["exact_specialization_rank_lower_bound"],
        ),
    )
    family_maxima = []
    for family in FAMILIES:
        family_records = [
            record for record in completed_points if record["family_index"] == family.index
        ]
        maximum_rank = max(
            record["exact_specialization_rank_lower_bound"] for record in family_records
        )
        leader = min(
            (
                record
                for record in family_records
                if record["exact_specialization_rank_lower_bound"] == maximum_rank
            ),
            key=lambda record: Decimal(str(record["global_curve"]["log_conductor"])),
        )
        family_maxima.append(compact_frontier_record(leader))
    summary = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed bounded four-family split-infinity screen",
        "claim_level": "bounded search plus exact certificates; negative result is not an upper bound",
        "candidate_input": {
            "path": str((output / "candidate-input.json").relative_to(ROOT)),
            "sha256": sha256_file(output / "candidate-input.json"),
            "result_sha256": input_payload["result_sha256"],
        },
        "point_search_input": {
            "path": str((output / "point-search-input.json").relative_to(ROOT)),
            "sha256": sha256_file(output / "point-search-input.json"),
            "result_sha256": point_input["result_sha256"],
        },
        "conductor_summary": {
            "selected": len(candidates),
            "completed": len(completed_conductors),
            "timeouts": sum(record["status"] == "timeout" for record in conductor_results),
            "errors": sum(record["status"] == "error" for record in conductor_results),
            "below_log_182_72": len(qualified),
            "below_by_family": dict(
                sorted(Counter(record["family_index"] for record in qualified).items())
            ),
        },
        "point_summary": {
            "selected": len(qualified),
            "completed": len(completed_points),
            "timeouts": sum(record["status"] == "timeout" for record in point_results),
            "errors": sum(record["status"] == "error" for record in point_results),
            "rank_lower_bound_distribution": dict(
                sorted(
                    Counter(
                        record["exact_specialization_rank_lower_bound"]
                        for record in completed_points
                    ).items()
                )
            ),
            "maximum_exact_rank_lower_bound": max(
                (
                    record["exact_specialization_rank_lower_bound"]
                    for record in completed_points
                ),
                default=None,
            ),
            "rank_at_least_17_promotions": [
                {
                    "family_index": record["family_index"],
                    "u": record["u"],
                    "T": record["base_T"],
                    "exact_rank_lower_bound": record[
                        "exact_specialization_rank_lower_bound"
                    ],
                    "conductor": record["global_curve"]["conductor"],
                    "log_conductor": record["global_curve"]["log_conductor"],
                    "root_number": record["global_curve"]["root_number"],
                    "point_certificate": str(
                        (
                            point_directory
                            / f"{candidate_identifier(record)}.json"
                        ).relative_to(ROOT)
                    ),
                }
                for record in promotions
            ],
            "retried_initial_error_records": retried_point_errors,
            "repaired_legacy_point_digests": repaired_legacy_point_digests,
        },
        "leaders": [
            {
                "family_index": record["family_index"],
                "u": record["u"],
                "T": record["base_T"],
                "exact_rank_lower_bound": record[
                    "exact_specialization_rank_lower_bound"
                ],
                "log_conductor": record["global_curve"]["log_conductor"],
                "root_number": record["global_curve"]["root_number"],
                "certificate_sha256": sha256_file(
                    point_directory / f"{candidate_identifier(record)}.json"
                ),
            }
            for record in completed_points[:20]
        ],
        "rank_conductor_pareto_frontier": {
            "scope": "61 completed exact-conductor-qualified and point-searched fibers",
            "order": "minimize exact conductor while maximizing certified rank lower bound",
            "records": [compact_frontier_record(record) for record in pareto_records],
        },
        "family_rank_maxima_with_lowest_conductor_tiebreak": family_maxima,
        "software": {
            "python": os.sys.version,
            "ratpoints": str(RATPOINTS.relative_to(ROOT)),
            "ratpoints_sha256": sha256_file(RATPOINTS),
            "pari_gp": gp_version(),
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves/cas python3 "
            "elliptic-curves/cas/search_mestre_dsquare_four.py --workers 8"
        ),
    }
    summary["result_sha256"] = canonical_digest(summary)
    atomic_json(output / "summary.json", summary)
    print(
        f"DONE completed_points={len(completed_points)} "
        f"max_rank={summary['point_summary']['maximum_exact_rank_lower_bound']} "
        f"promotions={len(promotions)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
