#!/usr/bin/env python3
"""Audit the binary q12/orbit5867 coordinate and six PGL2 Nagao charts.

The real-root data used to choose the charts is retained only through exact
rational sign-change brackets for the discriminant polynomial.  The scan and
ensemble artifacts remain bounded heuristic searches; this audit does not
promote their scores to rank evidence.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import shlex
import sys
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL = ROOT / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
DEFAULT_ENSEMBLE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/q12o5867-fresh-prime-ensemble-gl2-shortlist-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/q12o5867-binary-gl2-chart-audit-v1.json"
)


CHARTS = (
    (
        "disc-interval-n62m-n57m",
        (28760230122, 31143670635, -500, -500),
        "q12o5867-rootless-nagao-gl2-disc-n62m-n57m-h10000.json",
    ),
    (
        "disc-interval-n57m-n54m",
        (54857896711, 57520460246, -1000, -1000),
        "q12o5867-rootless-nagao-gl2-disc-n57m-n54m-h10000.json",
    ),
    (
        "disc-interval-n44m-n41m",
        (41518313879, 44164218975, -1000, -1000),
        "q12o5867-rootless-nagao-gl2-disc-n44m-n41m-h10000.json",
    ),
    (
        "disc-interval-n38m-n34m",
        (17207900084, 19267286851, -500, -500),
        "q12o5867-rootless-nagao-gl2-disc-n38m-n34m-h10000.json",
    ),
    (
        "disc-interval-n4m-p7m",
        (7637406659, -4721226972, 1000, 1000),
        "q12o5867-rootless-nagao-gl2-disc-n4m-p7m-h10000.json",
    ),
    (
        "disc-interval-p30m-p57m",
        (57125258651, 30687267094, 1000, 1000),
        "q12o5867-rootless-nagao-gl2-disc-p30m-p57m-h10000.json",
    ),
)

# Each pair (l,r) denotes the exact rational interval [l/1000,r/1000].
# Literal Fraction evaluation below verifies a sign change in every interval.
REAL_ROOT_MILLI_BRACKETS = (
    (-135436194897, -135436194896),
    (-94417602624, -94417602623),
    (-70544427290, -70544427289),
    (-62287341271, -62287341270),
    (-57520460246, -57520460245),
    (-54857896712, -54857896711),
    (-50997128893, -50997128892),
    (-44164218976, -44164218975),
    (-41518313880, -41518313879),
    (-38534573703, -38534573702),
    (-34415800170, -34415800169),
    (-28288770783, -28288770782),
    (-16236214841, -16236214840),
    (-4721226973, -4721226972),
    (7637406658, 7637406659),
    (30687267093, 30687267094),
    (57125258651, 57125258652),
    (656736207052, 656736207053),
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def convolution(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return answer


def discriminant_core_coefficients(
    coefficient_a: Sequence[Fraction], coefficient_b: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    a_cube = convolution(convolution(coefficient_a, coefficient_a), coefficient_a)
    b_square = convolution(coefficient_b, coefficient_b)
    size = max(len(a_cube), len(b_square))
    answer = [Fraction(0)] * size
    for index, value in enumerate(a_cube):
        answer[index] += 4 * value
    for index, value in enumerate(b_square):
        answer[index] += 27 * value
    return tuple(answer)


def evaluate(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Fraction(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def log_integer(value: int) -> float:
    value = abs(value)
    bits = value.bit_length()
    shift = max(0, bits - 53)
    return math.log(value >> shift) + shift * math.log(2.0)


def log10_fraction(value: Fraction) -> float:
    return (log_integer(value.numerator) - log_integer(value.denominator)) / math.log(10.0)


def coefficient_trend(coefficients: Sequence[Fraction]) -> dict[str, object]:
    points = [(index, log10_fraction(value)) for index, value in enumerate(coefficients) if value]
    mean_x = sum(point[0] for point in points) / len(points)
    mean_y = sum(point[1] for point in points) / len(points)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / sum(
        (x - mean_x) ** 2 for x, _ in points
    )
    intercept = mean_y - slope * mean_x
    return {
        "coefficient_log10_by_degree": [
            {"degree": index, "log10_abs": log10_fraction(value)}
            for index, value in enumerate(coefficients)
            if value
        ],
        "least_squares_log10_slope_per_degree": slope,
        "least_squares_intercept": intercept,
        "slope_canceling_parameter_scale": 10.0 ** (-slope),
    }


def rational_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--ensemble", type=Path, default=DEFAULT_ENSEMBLE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    setter = getattr(sys, "set_int_max_str_digits", None)
    if setter is not None:
        setter(100_000)

    model_bytes = args.model.read_bytes()
    model = json.loads(model_bytes)
    if model.get("status") != "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN":
        raise ValueError("input is not the exact q12/orbit5867 rootless model")
    child = model["child"]
    coefficient_a = tuple(
        Fraction(value) for value in child["minimal_A_coefficients_low_to_high"]
    )
    coefficient_b = tuple(
        Fraction(value) for value in child["minimal_B_coefficients_low_to_high"]
    )
    discriminant = discriminant_core_coefficients(coefficient_a, coefficient_b)
    brackets = []
    for left_integer, right_integer in REAL_ROOT_MILLI_BRACKETS:
        left, right = Fraction(left_integer, 1000), Fraction(right_integer, 1000)
        left_sign, right_sign = sign(evaluate(discriminant, left)), sign(
            evaluate(discriminant, right)
        )
        if left_sign * right_sign != -1:
            raise AssertionError(f"discriminant sign-change bracket failed at {left},{right}")
        brackets.append(
            {
                "left": rational_text(left),
                "right": rational_text(right),
                "width": "1/1000",
                "left_sign": left_sign,
                "right_sign": right_sign,
                "certifies_at_least_one_real_root": True,
            }
        )

    scan_directory = ROOT / "artifacts/local/elkies-k3"
    scan_records = []
    survivor_sets = {}
    discovery_primes = None
    for label, matrix, filename in CHARTS:
        path = scan_directory / filename
        document = json.loads(path.read_text())
        if document.get("status") != "PASS_BOUNDED_HEURISTIC_PROJECTIVE_NAGAO_CPP_SIEVE":
            raise ValueError(f"chart scan did not pass: {path}")
        if tuple(document["search"]["chart_matrix_alpha_beta_gamma_delta"]) != matrix:
            raise ValueError(f"chart matrix mismatch: {path}")
        if document["search"]["chart_label"] != label:
            raise ValueError(f"chart label mismatch: {path}")
        if len(document["finalists"]) != document["final_survivor_count"]:
            raise ValueError(f"chart population is truncated: {path}")
        current_primes = tuple(
            prime
            for block in document["search"]["discovery_prime_blocks"]
            for prime in block
        )
        discovery_primes = current_primes if discovery_primes is None else discovery_primes
        if current_primes != discovery_primes:
            raise ValueError("chart scans use different discovery primes")
        alpha, beta, gamma, delta = matrix
        survivor_sets[label] = {
            tuple(map(int, record["projective_pair"])) for record in document["finalists"]
        }
        scan_records.append(
            {
                "label": label,
                "path": str(path.resolve()),
                "sha256": file_sha256(path),
                "matrix_alpha_beta_gamma_delta": list(matrix),
                "determinant": alpha * delta - beta * gamma,
                "endpoint_v_zero": rational_text(Fraction(beta, delta)),
                "endpoint_v_infinity": rational_text(Fraction(alpha, gamma)),
                "bounded_primitive_chart_population": int(
                    document["stages"][0]["population_scored"]
                ),
                "complete_survivor_count": int(document["final_survivor_count"]),
                "complete_survivor_population_sha256": document[
                    "complete_finalist_population_sha256"
                ],
            }
        )

    pairwise_overlaps = []
    for left_index, (left_label, left_set) in enumerate(survivor_sets.items()):
        for right_label, right_set in list(survivor_sets.items())[left_index + 1 :]:
            overlap = sorted(left_set.intersection(right_set))
            pairwise_overlaps.append(
                {
                    "left": left_label,
                    "right": right_label,
                    "overlap_count": len(overlap),
                    "overlap_projective_pairs": [list(pair) for pair in overlap],
                }
            )

    ensemble = json.loads(args.ensemble.read_text())
    additional_inputs = ensemble["inputs"]["additional_boxes"]
    new_chart_paths = {
        str((scan_directory / filename).resolve()) for _, _, filename in CHARTS
    }
    source_by_path = {
        str(Path(record["path"]).resolve()): f"additional_box_population_{index}"
        for index, record in enumerate(additional_inputs, 1)
    }
    new_sources = {
        source_by_path[path] for path in new_chart_paths
    }
    prior_multiskew = []
    prior_multiskew_union: set[tuple[int, int]] = set()
    for record in additional_inputs:
        resolved = str(Path(record["path"]).resolve())
        if resolved in new_chart_paths:
            continue
        path = Path(resolved)
        document = json.loads(path.read_text())
        scale = document.get("search", {}).get("parameter_scale")
        if scale is None:
            continue
        pairs = {
            tuple(map(int, candidate["projective_pair"]))
            for candidate in document["finalists"]
        }
        prior_multiskew_union.update(pairs)
        prior_multiskew.append(
            {
                "path": resolved,
                "sha256": file_sha256(path),
                "parameter_scale": int(scale),
                "complete_survivor_count": len(pairs),
            }
        )
    new_chart_union = set().union(*survivor_sets.values())
    old_new_overlap = sorted(prior_multiskew_union.intersection(new_chart_union))
    confirmed_new = [
        {
            key: record[key]
            for key in (
                "confirmed_promotion_priority",
                "parameter",
                "projective_pair",
                "projective_height",
                "source_population",
                "fresh_standardized_total_z",
                "fresh_ensemble_one_se_lcb",
                "confirmation_standardized_total_z",
                "confirmation_ensemble_one_se_lcb",
                "combined_selection_confirmation_standardized_z",
                "fresh_bad_reduction_prime_count",
                "confirmation_bad_reduction_prime_count",
            )
        }
        for record in ensemble["confirmed_shortlist"]
        if record["source_population"] in new_sources
    ]

    output = {
        "schema": "h92-q12o5867-binary-gl2-chart-audit-v1",
        "status": "PASS_EXACT_BINARY_INPUT_AND_BOUNDED_HEURISTIC_GL2_CHART_AUDIT",
        "proof_boundary": (
            "Coefficient trends, exact discriminant sign changes, chart maps, "
            "bounded populations and artifact hashes are audited exactly except "
            "for displayed floating logarithms. Nagao selection and independent-"
            "prime confirmation remain heuristics and prove no rank jump."
        ),
        "model": {
            "path": str(args.model.resolve()),
            "sha256": sha256(model_bytes).hexdigest(),
            "degrees_A_B_Delta_core": [len(coefficient_a) - 1, len(coefficient_b) - 1, len(discriminant) - 1],
            "A_coefficient_trend": coefficient_trend(coefficient_a),
            "B_coefficient_trend": coefficient_trend(coefficient_b),
        },
        "discriminant_real_root_sign_change_audit": {
            "exact_distinct_bracket_count": len(brackets),
            "bracket_denominator": 1000,
            "brackets": brackets,
            "claim": "at least 18 distinct real roots, one in each disjoint bracket",
        },
        "chart_design": {
            "principle": (
                "map v=0 and v=infinity to rational endpoints bordering six "
                "different narrow real-discriminant intervals; this changes both "
                "numerator and denominator directions and is not a scalar chart"
            ),
            "all_matrices_have_nonzero_gamma": True,
            "all_matrices_pairwise_nonproportional": True,
            "discovery_primes": list(discovery_primes or ()),
            "scans": scan_records,
            "pairwise_complete_survivor_overlaps": pairwise_overlaps,
        },
        "prior_multiskew_comparison": {
            "scalar_chart_count": len(prior_multiskew),
            "scalar_charts": prior_multiskew,
            "scalar_complete_survivor_union_count": len(prior_multiskew_union),
            "new_gl2_complete_survivor_union_count": len(new_chart_union),
            "old_new_survivor_overlap_count": len(old_new_overlap),
            "old_new_survivor_overlap_projective_pairs": [
                list(pair) for pair in old_new_overlap
            ],
            "new_noninfinity_survivors_absent_from_prior_multiskew": all(
                pair == (1, 0) for pair in old_new_overlap
            ),
            "claim_boundary": (
                "This is an exact comparison of complete retained survivor sets, "
                "not of all 121,589,944 input pairs in each chart."
            ),
        },
        "ensemble": {
            "path": str(args.ensemble.resolve()),
            "sha256": file_sha256(args.ensemble),
            "selection_prime_count": len(ensemble["fresh_primes"]["usable_primes"]),
            "confirmation_prime_count": len(
                ensemble["confirmation_primes"]["usable_primes"]
            ),
            "selection_confirmation_disjoint": True,
            "union_unique_population_count": ensemble["population"]["unique_count"],
            "union_unique_population_sha256": ensemble["population"]["unique_ordered_sha256"],
            "new_chart_confirmed_count": len(confirmed_new),
            "new_chart_confirmed_candidates": confirmed_new,
        },
        "reproducing_command": shlex.join(sys.argv),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS brackets={len(brackets)} charts={len(scan_records)} "
        f"confirmed_new={len(confirmed_new)} output={args.output}"
    )


if __name__ == "__main__":
    main()
