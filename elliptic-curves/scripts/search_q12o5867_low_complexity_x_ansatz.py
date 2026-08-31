#!/usr/bin/env python3
"""Search q12o5867 rank jumps with low-complexity section x-interpolants.

For ordered generic section abscissas ``x_i(u)`` and a small rational lambda,
set

    x(u) = (1-lambda)*x_i(u) + lambda*x_j(u).

Then ``x(u)^3 + A(u)*x(u) + B(u)`` has degree at most twelve.  PARI's
``hyperellratpoints`` searches a declared rational-parameter box on this
exact hyperelliptic curve.  Every hit is homogenized back to the projective
q12o5867 specialization, checked by exact substitution, and measured in
finite good-reduction quotients relative to the specialized generic rank-17
subgroup.  No eclib or Mordell--Weil initialization is used.

This is a bounded search.  Positive quotient escape is exact; no hit and
finite-quotient non-escape are not rank upper bounds.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt, lcm
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Sequence


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    integral_scaling_factor,
    integralize_specialization,
    load_q12o5867_data,
)
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    verify_finite_quotient_certificate,
)
from finite_quotient_escape import QuotientBlock, analyze_escape  # noqa: E402
from pari_bridge import pari_version  # noqa: E402


Q = Fraction
Point = tuple[Fraction, Fraction]
DEFAULT_MODEL = REPOSITORY / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
DEFAULT_SECTIONS = (
    REPOSITORY
    / "artifacts/local/elkies-k3/q12o5867-rootless-selected-basis-qq.json"
)
DEFAULT_ANCHORS = ((-7801, 1463), (601, 418), (677, 3402), (-267, 847))


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_text(value: Fraction | int) -> str:
    value = Q(value)
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def point_record(point: Point) -> list[str]:
    return [fraction_text(point[0]), fraction_text(point[1])]


def parse_fractions(text: str) -> tuple[Fraction, ...]:
    try:
        answer = tuple(Q(value) for value in text.split(",") if value)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError("expected comma-separated rationals") from error
    if not answer:
        raise argparse.ArgumentTypeError("the rational list must be nonempty")
    return answer


def parse_int_tuple(text: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not answer:
        raise argparse.ArgumentTypeError("the integer list must be nonempty")
    return answer


def parse_anchors(text: str) -> tuple[tuple[int, int], ...]:
    answer = []
    try:
        for item in text.split(","):
            numerator, denominator = item.split("/", 1)
            value = Q(int(numerator), int(denominator))
            answer.append((value.numerator, value.denominator))
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(
            "anchors must be comma-separated rational parameters"
        ) from error
    if not answer:
        raise argparse.ArgumentTypeError("at least one anchor is required")
    return tuple(answer)


def polynomial_add(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    answer = [Q(0) for _ in range(max(len(left), len(right)))]
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def polynomial_scale(values: Sequence[Fraction], scale: Fraction) -> list[Fraction]:
    return [scale * value for value in values]


def polynomial_multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> list[Fraction]:
    answer = [Q(0) for _ in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def polynomial_power(values: Sequence[Fraction], exponent: int) -> list[Fraction]:
    answer = [Q(1)]
    for _ in range(exponent):
        answer = polynomial_multiply(answer, values)
    return answer


def polynomial_value(values: Sequence[Fraction], parameter: Fraction) -> Fraction:
    answer = Q(0)
    for value in reversed(values):
        answer = answer * parameter + value
    return answer


def pad(values: Sequence[Fraction], length: int) -> tuple[Fraction, ...]:
    if len(values) > length:
        raise ValueError("polynomial exceeds its declared degree")
    return tuple(values) + (Q(0),) * (length - len(values))


def ansatz_polynomials(
    data: Any, section_i: int, section_j: int, interpolation: Fraction
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    x_i = data.section_coefficients[section_i][0]
    x_j = data.section_coefficients[section_j][0]
    x_values = pad(
        polynomial_add(
            polynomial_scale(x_i, 1 - interpolation),
            polynomial_scale(x_j, interpolation),
        ),
        5,
    )
    right_side = pad(
        polynomial_add(
            polynomial_add(
                polynomial_power(x_values, 3),
                polynomial_multiply(data.a_coefficients, x_values),
            ),
            data.b_coefficients,
        ),
        13,
    )
    return x_values, right_side


def primitive_complexity(values: Sequence[Fraction]) -> tuple[int, int, int]:
    denominator = lcm(*(value.denominator for value in values))
    integers = [int(value * denominator) for value in values]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    if common:
        integers = [value // common for value in integers]
    bits = [abs(value).bit_length() for value in integers if value]
    return max(bits, default=0), sum(bits), denominator.bit_length()


def select_directions(
    data: Any,
    interpolations: Sequence[Fraction],
    direction_count: int,
) -> tuple[dict[str, Any], ...]:
    directions = []
    for section_i in range(17):
        for section_j in range(section_i + 1, 17):
            for interpolation in interpolations:
                if interpolation in (0, 1):
                    continue
                x_values, right_side = ansatz_polynomials(
                    data, section_i, section_j, interpolation
                )
                complexity = primitive_complexity(right_side)
                directions.append(
                    {
                        "section_i": section_i,
                        "section_j": section_j,
                        "lambda": interpolation,
                        "x_coefficients": x_values,
                        "right_side_coefficients": right_side,
                        "complexity": complexity,
                    }
                )
    directions.sort(
        key=lambda row: (
            row["complexity"],
            row["section_i"],
            row["section_j"],
            row["lambda"].numerator,
            row["lambda"].denominator,
        )
    )
    return tuple(directions[:direction_count])


def rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def gp_polynomial(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"({fraction_text(value)})*x^{index}"
        for index, value in enumerate(coefficients)
    )


def run_hyperellratpoints(
    coefficients: Sequence[Fraction],
    *,
    numerator_bound: int,
    denominator_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[Point, ...], dict[str, Any]]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    height = f"[{numerator_bound},{denominator_bound}]"
    program = "\n".join(
        (
            f"Q={gp_polynomial(coefficients)};gettime();",
            f"R=hyperellratpoints(Q,{height});",
            'print("PARI_MILLISECONDS|",gettime());',
            'print("POINTS_BEGIN");print(R);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [executable, "-q", "-f", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        return (), {
            "status": "timeout",
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
            "stdout_tail": str(error.stdout or "")[-1000:],
            "stderr_tail": str(error.stderr or "")[-1000:],
            "program_sha256": sha256(program.encode()).hexdigest(),
        }
    wall_seconds = time.monotonic() - started
    fatal = [line for line in completed.stderr.splitlines() if "***" in line]
    if completed.returncode != 0 or fatal:
        return (), {
            "status": "pari_error",
            "returncode": completed.returncode,
            "wall_seconds": wall_seconds,
            "error": " ".join(fatal or completed.stderr.splitlines())[:2000],
            "program_sha256": sha256(program.encode()).hexdigest(),
        }
    marker = re.search(
        r"POINTS_BEGIN\s*(.*?)\s*POINTS_END", completed.stdout, re.DOTALL
    )
    milliseconds = re.search(
        r"^PARI_MILLISECONDS\|(\d+)$", completed.stdout, re.MULTILINE
    )
    if marker is None or milliseconds is None:
        raise RuntimeError("PARI omitted x-ansatz search markers")
    points = tuple(
        (Q(parameter), Q(ordinate))
        for parameter, ordinate in re.findall(
            r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]", marker.group(1)
        )
    )
    if any(
        ordinate**2 != polynomial_value(coefficients, parameter)
        for parameter, ordinate in points
    ):
        raise AssertionError("PARI returned a point off an exact ansatz curve")
    return points, {
        "status": "completed",
        "height_specification": height,
        "wall_seconds": wall_seconds,
        "pari_milliseconds": int(milliseconds.group(1)),
        "signed_point_count": len(points),
        "program_sha256": sha256(program.encode()).hexdigest(),
    }


def exact_quotient_profile(
    model: Sequence[Fraction],
    baseline: Sequence[Point],
    candidates: Sequence[Point],
    *,
    relation_primes: Sequence[int],
    reduction_prime_bound: int,
) -> dict[str, Any]:
    if not candidates:
        return {
            "candidate_count": 0,
            "maximum_marginal_dimension": 0,
            "independent_escape_basis_labels": [],
            "relation_prime_profiles": [],
        }
    all_points = (*baseline, *candidates)
    labels = tuple(f"candidate-{index:03d}" for index in range(len(candidates)))
    attempts = []
    maximum = 0
    basis: list[str] = []
    for relation_prime in relation_primes:
        try:
            certificate = build_finite_quotient_certificate(
                model,
                all_points,
                relation_prime=relation_prime,
                prime_bound=reduction_prime_bound,
            )
            verify_finite_quotient_certificate(model, all_points, certificate)
            blocks = tuple(
                QuotientBlock.build(
                    modulus=relation_prime,
                    rows=signature["rows"],
                    column_count=len(all_points),
                    source=f"good-reduction-p={signature['prime']}",
                )
                for signature in certificate["signatures"]
            )
            profile = analyze_escape(
                blocks,
                known_column_count=len(baseline),
                candidate_labels=labels,
            )
            record = profile.to_record()
            record["finite_quotient_certificate"] = certificate
            attempts.append(record)
            if profile.marginal_dimension > maximum:
                maximum = profile.marginal_dimension
                basis = list(profile.independent_escape_basis_labels)
        except (RuntimeError, ValueError, AssertionError) as error:
            attempts.append(
                {
                    "modulus": relation_prime,
                    "status": "bounded-certificate-error",
                    "exception": type(error).__name__,
                    "detail": str(error),
                }
            )
    return {
        "candidate_count": len(candidates),
        "candidate_labels": list(labels),
        "maximum_marginal_dimension": maximum,
        "independent_escape_basis_labels": basis,
        "relation_prime_profiles": attempts,
        "promotion_threshold": 15,
        "promotion_eligible": maximum >= 15,
        "claim_boundary": (
            "positive finite-quotient escape is exact; non-escape is not a "
            "dependence proof"
        ),
    }


def certify_parameter(
    data: Any,
    parameter: tuple[int, int],
    hits: Sequence[dict[str, Any]],
    args: argparse.Namespace,
) -> dict[str, Any]:
    a, b = parameter
    specialization = evaluate_projective_specialization(data, a, b)
    projective_model = specialization.model
    baseline = specialization.points
    extra_by_x: dict[Fraction, dict[str, Any]] = {}
    for hit in hits:
        direction = hit["direction"]
        affine_parameter = Q(a, b)
        affine_x = polynomial_value(direction["x_coefficients"], affine_parameter)
        affine_y = Q(hit["ordinate"])
        projective_point = (b**4 * affine_x, b**6 * affine_y)
        if not is_on_weierstrass_curve(projective_model, projective_point):
            raise AssertionError("an ansatz hit missed its exact projective specialization")
        if projective_point[0] in {point[0] for point in baseline}:
            continue
        canonical = min(
            (projective_point, (projective_point[0], -projective_point[1])),
            key=lambda point: (
                point[0].numerator,
                point[0].denominator,
                point[1].numerator,
                point[1].denominator,
            ),
        )
        extra_by_x.setdefault(
            canonical[0],
            {
                "point": canonical,
                "direction_id": hit["direction_id"],
                "source": hit["source"],
            },
        )
    scale, scale_record = integral_scaling_factor(
        specialization.coefficient_a,
        specialization.coefficient_b,
        timeout=min(args.certificate_timeout, 120.0),
        stack_bytes=min(args.stack_bytes, 128_000_000),
    )
    integral_model, integral_baseline = integralize_specialization(specialization, scale)
    extras = tuple(
        (record["point"][0] * scale**2, record["point"][1] * scale**3)
        for _x, record in sorted(extra_by_x.items(), key=lambda item: item[0])
    )
    if any(not is_on_weierstrass_curve(integral_model, point) for point in extras):
        raise AssertionError("an integralized ansatz point missed the curve")
    quotient = exact_quotient_profile(
        integral_model,
        integral_baseline,
        extras,
        relation_primes=args.relation_primes,
        reduction_prime_bound=args.reduction_prime_bound,
    )
    extra_records = []
    for label, (x_coordinate, source), point in zip(
        quotient.get("candidate_labels", []),
        sorted(extra_by_x.items(), key=lambda item: item[0]),
        extras,
        strict=True,
    ):
        extra_records.append(
            {
                "label": label,
                "point_on_integral_short_model": point_record(point),
                "direction_id": source["direction_id"],
                "source": source["source"],
                "exact_curve_membership_verified": True,
            }
        )
    return {
        "parameter": {
            "normalized_projective": [a, b],
            "affine_value": fraction_text(Q(a, b)),
        },
        "hit_count_before_point_deduplication": len(hits),
        "integralization": scale_record,
        "integral_short_model": [fraction_text(value) for value in integral_model],
        "exact_points_beyond_the_listed_17": extra_records,
        "finite_quotient_escape": quotient,
        "status": (
            "EXACT_QUOTIENT_ESCAPE_FOUND"
            if quotient["maximum_marginal_dimension"] > 0
            else "EXACT_HITS_WITH_NO_DETECTED_QUOTIENT_GAIN"
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--sections", type=Path, default=DEFAULT_SECTIONS)
    parser.add_argument(
        "--lambdas", type=parse_fractions, default=(Q(1, 2), Q(-1), Q(2))
    )
    parser.add_argument("--directions", type=int, default=24)
    parser.add_argument("--numerator-bound", type=int, default=20_000)
    parser.add_argument("--denominator-bound", type=int, default=200)
    parser.add_argument("--direction-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument("--anchors", type=parse_anchors, default=DEFAULT_ANCHORS)
    parser.add_argument("--relation-primes", type=parse_int_tuple, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic-curves"
            / "q12o5867-low-complexity-x-ansatz-search.json"
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    sys.set_int_max_str_digits(0)
    args = build_parser().parse_args()
    if args.directions < 1:
        raise SystemExit("--directions must be positive")
    if args.numerator_bound < 1 or args.denominator_bound < 1:
        raise SystemExit("search bounds must be positive")
    if args.denominator_bound > args.numerator_bound:
        raise SystemExit("PARI requires denominator bound <= numerator bound")
    if args.direction_timeout <= 0 or args.certificate_timeout <= 0:
        raise SystemExit("timeouts must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("--stack-bytes is too small")

    started = time.monotonic()
    data = load_q12o5867_data(args.model.resolve(), args.sections.resolve())
    directions = select_directions(data, args.lambdas, args.directions)
    hits_by_parameter: dict[tuple[int, int], list[dict[str, Any]]] = {}
    direction_records = []
    for direction_id, direction in enumerate(directions):
        right_side = direction["right_side_coefficients"]
        points, search = run_hyperellratpoints(
            right_side,
            numerator_bound=args.numerator_bound,
            denominator_bound=args.denominator_bound,
            timeout=args.direction_timeout,
            stack_bytes=args.stack_bytes,
        )
        anchor_hits = []
        for a, b in args.anchors:
            value = polynomial_value(right_side, Q(a, b))
            root = rational_square_root(value)
            if root is None:
                continue
            anchor_hits.append((a, b))
            hits_by_parameter.setdefault((a, b), []).append(
                {
                    "direction_id": direction_id,
                    "direction": direction,
                    "ordinate": root,
                    "source": "explicit-current-anchor-evaluation",
                }
            )
        parameter_hits = set()
        for parameter, ordinate in points:
            key = (parameter.numerator, parameter.denominator)
            parameter_hits.add(key)
            hits_by_parameter.setdefault(key, []).append(
                {
                    "direction_id": direction_id,
                    "direction": direction,
                    "ordinate": ordinate,
                    "source": "PARI-hyperellratpoints",
                }
            )
        direction_records.append(
            {
                "direction_id": direction_id,
                "section_i": direction["section_i"],
                "section_j": direction["section_j"],
                "lambda": fraction_text(direction["lambda"]),
                "complexity_max_sum_denominator_bits": list(direction["complexity"]),
                "x_coefficients_low_to_high": [
                    fraction_text(value) for value in direction["x_coefficients"]
                ],
                "right_side_coefficients_low_to_high": [
                    fraction_text(value) for value in right_side
                ],
                "right_side_sha256": sha256(
                    json.dumps(
                        [fraction_text(value) for value in right_side],
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest(),
                "search": search,
                "distinct_parameter_hits": len(parameter_hits),
                "current_anchor_hits": [[a, b] for a, b in anchor_hits],
            }
        )

    certifications = []
    rejected_singular = []
    for parameter, hits in sorted(hits_by_parameter.items()):
        try:
            certifications.append(certify_parameter(data, parameter, hits, args))
        except ValueError as error:
            if "singular fibre" not in str(error):
                raise
            rejected_singular.append(
                {"parameter": list(parameter), "detail": str(error)}
            )

    command = " ".join(
        (
            ".venv/bin/python",
            "elliptic-curves/scripts/search_q12o5867_low_complexity_x_ansatz.py",
            "--lambdas",
            shlex.quote(",".join(fraction_text(value) for value in args.lambdas)),
            "--directions",
            str(args.directions),
            "--numerator-bound",
            str(args.numerator_bound),
            "--denominator-bound",
            str(args.denominator_bound),
            "--direction-timeout",
            str(args.direction_timeout),
            "--certificate-timeout",
            str(args.certificate_timeout),
            "--stack-bytes",
            str(args.stack_bytes),
            "--anchors",
            shlex.quote(",".join(f"{a}/{b}" for a, b in args.anchors)),
            "--relation-primes",
            ",".join(str(value) for value in args.relation_primes),
            "--reduction-prime-bound",
            str(args.reduction_prime_bound),
            "--output",
            shlex.quote(str(args.output)),
            "--overwrite",
        )
    )
    maximum_gain = max(
        (
            row["finite_quotient_escape"]["maximum_marginal_dimension"]
            for row in certifications
        ),
        default=0,
    )
    artifact = {
        "schema": "elliptic-curves.q12o5867-low-complexity-x-ansatz.v1",
        "status": (
            "EXACT_QUOTIENT_ESCAPE_FOUND"
            if maximum_gain > 0
            else "COMPLETE_BOUNDED_SEARCH_NO_DETECTED_QUOTIENT_GAIN"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "model": {"path": str(args.model), "sha256": data.model_sha256},
            "sections": {"path": str(args.sections), "sha256": data.sections_sha256},
        },
        "method": {
            "description": "x=(1-lambda)*x_i+lambda*x_j, degree-at-most-12 parameter curves",
            "engine": "PARI/GP hyperellratpoints",
            "pari_version": pari_version(),
            "eclib_used": False,
            "lambdas": [fraction_text(value) for value in args.lambdas],
            "selected_direction_count": len(directions),
            "selection_order": "primitive coefficient max bits, sum bits, denominator bits",
            "parameter_numerator_bound": args.numerator_bound,
            "parameter_denominator_bound": args.denominator_bound,
            "direction_timeout_seconds": args.direction_timeout,
        },
        "directions": direction_records,
        "distinct_nonsingular_hit_parameters": len(certifications),
        "rejected_singular_parameters": rejected_singular,
        "certifications": certifications,
        "maximum_exact_quotient_gain": maximum_gain,
        "wall_seconds": time.monotonic() - started,
        "script": {
            "path": str(Path(__file__).resolve().relative_to(REPOSITORY)),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
        "reproducing_command": command,
        "claim_boundary": [
            "Each direction is completely searched only inside its declared box.",
            "Positive finite-quotient escape is exact.",
            "No hit or detected escape is not a dependence result or rank upper bound.",
        ],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with output.open(mode) as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={artifact['status']}")
    print(f"directions={len(directions)}")
    print(f"distinct_nonsingular_hit_parameters={len(certifications)}")
    print(f"maximum_exact_quotient_gain={maximum_gain}")
    for row in certifications:
        print(
            f"parameter={row['parameter']['normalized_projective']} "
            f"hits={row['hit_count_before_point_deduplication']} "
            f"new_points={len(row['exact_points_beyond_the_listed_17'])} "
            "quotient_gain="
            f"{row['finite_quotient_escape']['maximum_marginal_dimension']}"
        )
    print(f"output={output}")


if __name__ == "__main__":
    main()
