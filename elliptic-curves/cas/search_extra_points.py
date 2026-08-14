#!/usr/bin/env python3
"""Bounded extra-point search on selected Fermigier quartics.

PARI's ``hyperellratpoints`` enumerates affine points on the exact genus-one
quartic up to a configured naive-height bound.  One point from each ``z`` / 
``-z`` pair is mapped exactly to the short Jacobian, then PARI recomputes the
Neron--Tate height matrix at two precisions and extracts a numerically
full-rank subset.  This is an exploratory bounded computation: numerical
height rank is evidence, not an exact independence certificate.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import time
from typing import Any

from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version


DEFAULT_PARAMETERS = ("1666/9",)
ELLRANK_OVERFLOW = {
    "parameter": "1666/9",
    "effort": 0,
    "seed_count": 12,
    "stack_bytes": 1_000_000_000,
    "timeout_seconds": 300,
    "elapsed_seconds_approx": 300,
    "result": "failed before returning rank bounds",
    "error": (
        "ellrank: the PARI stack overflows; current stack size: 1000000000 "
        "(953.674 Mbytes)"
    ),
}


def parse_parameters(value: str) -> tuple[str, ...]:
    parameters = tuple(item for item in value.split(",") if item)
    if not parameters:
        raise argparse.ArgumentTypeError("at least one rational parameter is required")
    try:
        for parameter in parameters:
            Fraction(parameter)
    except ValueError as error:
        raise argparse.ArgumentTypeError("parameters must be rational numbers") from error
    return parameters


def parse_precisions(value: str) -> tuple[int, ...]:
    try:
        precisions = tuple(int(item) for item in value.split(",") if item)
    except ValueError as error:
        raise argparse.ArgumentTypeError("precisions must be integers") from error
    if len(precisions) < 2 or precisions != tuple(sorted(set(precisions))):
        raise argparse.ArgumentTypeError("provide at least two increasing precisions")
    return precisions


def gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(value)})"


def gp_vector(point: tuple[Fraction, Fraction]) -> str:
    return f"[{gp_rational(point[0])},{gp_rational(point[1])}]"


def quartic_polynomial(parameter: Fraction) -> str:
    coefficients = FermigierMestreFamily.quartic_coefficients(parameter)
    return "+".join(
        f"{gp_rational(coefficient)}*x^{4-index}"
        for index, coefficient in enumerate(coefficients)
    )


def run_gp(program: str, *, timeout: float, stack_bytes: int) -> tuple[str, float]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    started = time.monotonic()
    result = subprocess.run(
        [executable, "-q", "-s", str(stack_bytes)],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - started
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP failed: {result.stderr.strip()}")
    return result.stdout, elapsed


def parse_point_vector(output: str) -> tuple[tuple[Fraction, Fraction], ...]:
    pairs = re.findall(
        r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]", output
    )
    return tuple((Fraction(x_value), Fraction(y_value)) for x_value, y_value in pairs)


def search_quartic_points(
    parameter: Fraction, *, height_bound: int, timeout: float, stack_bytes: int
) -> tuple[tuple[tuple[Fraction, Fraction], ...], float, int]:
    program = "\n".join(
        (
            f"Q={quartic_polynomial(parameter)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS ",R);',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    match = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if match is None:
        raise AssertionError("PARI did not report the search CPU time")
    points = parse_point_vector(output.split("POINTS ", 1)[1])
    return points, wall_seconds, int(match.group(1))


def signless_quartic_points(
    points: tuple[tuple[Fraction, Fraction], ...]
) -> tuple[tuple[Fraction, Fraction], ...]:
    answers: list[tuple[Fraction, Fraction]] = []
    seen_x: set[Fraction] = set()
    for point in points:
        if point[0] in seen_x:
            continue
        answers.append(point)
        seen_x.add(point[0])
    return tuple(answers)


def parse_vecsmall(text: str) -> list[int]:
    match = re.search(r"Vecsmall\(\[(.*?)\]\)", text)
    if match is None:
        raise AssertionError("PARI did not emit the subset index vector")
    return [int(value) for value in match.group(1).split(",")]


def height_replay(
    parameter: Fraction,
    jacobian_points: tuple[tuple[Fraction, Fraction], ...],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], ...]:
    coefficients = ",".join(
        gp_rational(value) for value in FermigierMestreFamily.coefficients(parameter)
    )
    points = ",".join(gp_vector(point) for point in jacobian_points)
    commands = [f"E=ellinit([{coefficients}]);", f"P=[{points}];"]
    for precision in precisions:
        commands.extend(
            (
                f"default(realprecision,{precision});",
                "H=ellheightmatrix(E,P);",
                "IX=matindexrank(H);",
                "K=vecextract(P,IX[2]);",
                "HK=ellheightmatrix(E,K);",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(matrank(H));",
                "print(IX[2]);",
                "print(matdet(HK));",
                "EV=mateigen(HK,1)[1];print(vecmin(EV));print(vecmax(EV));",
                f'print("HEIGHT_{precision}_END");',
            )
        )
    commands.append("quit")
    output, _ = run_gp(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    records = []
    for precision in precisions:
        start = lines.index(f"HEIGHT_{precision}_BEGIN") + 1
        end = lines.index(f"HEIGHT_{precision}_END")
        values = lines[start:end]
        records.append(
            {
                "decimal_precision": precision,
                "numerical_rank": int(values[0]),
                "subset_indices_one_based": parse_vecsmall(values[1]),
                "subset_height_determinant": values[2],
                "subset_smallest_eigenvalue": values[3],
                "subset_largest_eigenvalue": values[4],
            }
        )
    return tuple(records)


def exact_point_record(
    parameter: Fraction,
    quartic_point: tuple[Fraction, Fraction],
    jacobian_point: tuple[Fraction, Fraction],
) -> dict[str, str]:
    quartic_residual = (
        quartic_point[1] ** 2
        - FermigierMestreFamily.quartic_value(parameter, quartic_point[0])
    )
    _, _, _, coefficient_a, coefficient_b = FermigierMestreFamily.coefficients(
        parameter
    )
    jacobian_residual = (
        jacobian_point[1] ** 2
        - jacobian_point[0] ** 3
        - coefficient_a * jacobian_point[0]
        - coefficient_b
    )
    if quartic_residual or jacobian_residual:
        raise AssertionError("an exact point check failed")
    return {
        "quartic_x": rational_to_string(quartic_point[0]),
        "quartic_z": rational_to_string(quartic_point[1]),
        "jacobian_x": rational_to_string(jacobian_point[0]),
        "jacobian_y": rational_to_string(jacobian_point[1]),
        "quartic_residual": rational_to_string(quartic_residual),
        "jacobian_residual": rational_to_string(jacobian_residual),
    }


def explore_parameter(
    parameter_text: str, args: argparse.Namespace
) -> dict[str, Any]:
    parameter = Fraction(parameter_text)
    raw_points, wall_seconds, pari_milliseconds = search_quartic_points(
        parameter,
        height_bound=args.height_bound,
        timeout=args.timeout,
        stack_bytes=args.stack_bytes,
    )
    signless = signless_quartic_points(raw_points)
    jacobian = tuple(
        FermigierMestreFamily.quartic_point_to_jacobian(parameter, point)
        for point in signless
    )
    known_quartic_x = {
        point[0] for point in FermigierMestreFamily.known_quartic_points(parameter)
    }
    height_runs = height_replay(
        parameter,
        jacobian,
        precisions=args.precisions,
        timeout=args.timeout,
        stack_bytes=args.stack_bytes,
    )
    ranks = {run["numerical_rank"] for run in height_runs}
    subsets = {tuple(run["subset_indices_one_based"]) for run in height_runs}
    if len(ranks) != 1 or len(subsets) != 1:
        raise AssertionError("numerical rank or selected subset changed with precision")
    subset_indices = height_runs[-1]["subset_indices_one_based"]
    subset = [
        exact_point_record(parameter, signless[index - 1], jacobian[index - 1])
        for index in subset_indices
    ]
    return {
        "t": rational_to_string(parameter),
        "quartic_height_bound": args.height_bound,
        "search_scope": (
            "PARI hyperellratpoints affine rational points of naive height at most "
            f"{args.height_bound} on the implemented quartic model"
        ),
        "search_wall_seconds": wall_seconds,
        "pari_reported_milliseconds": pari_milliseconds,
        "signed_quartic_points_found": len(raw_points),
        "distinct_quartic_x_values": len(signless),
        "visible_section_x_values_found": sum(
            point[0] in known_quartic_x for point in signless
        ),
        "new_x_values_beyond_visible_sections": sum(
            point[0] not in known_quartic_x for point in signless
        ),
        "all_mapped_jacobian_points_checked_exactly": True,
        "height_matrix_runs": list(height_runs),
        "stable_numerical_rank": next(iter(ranks)),
        "explicit_numerically_independent_subset": subset,
        "status": (
            "bounded point enumeration plus high-precision numerical independence "
            "evidence; not an exact rank certificate"
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parameters", type=parse_parameters, default=DEFAULT_PARAMETERS
    )
    parser.add_argument("--height-bound", type=int, default=1_000_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(96, 192))
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_extra_points.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.height_bound < 1 or args.timeout <= 0 or args.stack_bytes < 8_000_000:
        raise SystemExit("height, timeout, and stack bounds must be positive")
    results = [explore_parameter(parameter, args) for parameter in args.parameters]
    script_path = Path(__file__).resolve()
    maximum_rank = max(result["stable_numerical_rank"] for result in results)
    reproducing_command = " ".join(
        (
            "PYTHONPATH=elliptic-curves/cas",
            ".venv/bin/python",
            "elliptic-curves/cas/search_extra_points.py",
            "--parameters",
            shlex.quote(",".join(args.parameters)),
            "--height-bound",
            str(args.height_bound),
            "--precisions",
            ",".join(str(value) for value in args.precisions),
            "--timeout",
            str(args.timeout),
            "--stack-bytes",
            str(args.stack_bytes),
            "--output",
            shlex.quote(str(args.output)),
        )
    )
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded experiment; exact curve membership is verified, but numerical "
            "height rank is not promoted to an exact Mordell--Weil rank claim"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": "182.72",
            "hit": False,
            "reason": (
                f"largest numerical rank found is {maximum_rank}, below 21, and "
                "no exact independence certificate was generated"
            ),
        },
        "family": "normalized Fermigier--Mestre family",
        "results": results,
        "prior_ellrank_attempt": ELLRANK_OVERFLOW,
        "methods_not_used": {
            "ellsearch": "PARI curve-database lookup, not a rational-point search",
            "ellsaturation": (
                "requires independent finite-index input and would not itself "
                "establish the missing exact-rank premise"
            ),
        },
        "parameters": {
            "specializations": list(args.parameters),
            "quartic_height_bound": args.height_bound,
            "precisions": list(args.precisions),
            "timeout_seconds_per_pari_call": args.timeout,
            "stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "source": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        "reproducing_command": reproducing_command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    for result in results:
        print(
            f"T={result['t']} signed_points={result['signed_quartic_points_found']} "
            f"x_values={result['distinct_quartic_x_values']} "
            f"new_x={result['new_x_values_beyond_visible_sections']} "
            f"numerical_rank={result['stable_numerical_rank']}"
        )


if __name__ == "__main__":
    main()
