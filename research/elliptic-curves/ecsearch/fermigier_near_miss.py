"""Pinned rank-at-least-20, low-conductor Fermigier near miss."""

from __future__ import annotations

import ast
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import shutil
import subprocess

from .fermigier_rank import (
    parse_ratpoints_output,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from .rank_certification import select_independent_subset


FERMIGIER_RANK20_PARAMETER = Fraction(28917, 20)
FERMIGIER_RANK20_SEARCH_HEIGHT = 2_000_000
FERMIGIER_RANK20_DENOMINATOR_BOUND = 13_000
FERMIGIER_RANK20_MAXIMUM_REDUCTION_PRIME = 2_000
FERMIGIER_RANK20_RELATION_PRIME = 5
STRICT_LOG_CONDUCTOR_MAXIMUM_INTEGER = int(
    "22609332114411420624526008180120083289443726642551045721326266745468787166869234"
)


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def canonical_ratpoints_output(abscissas: tuple[Fraction, ...]) -> str:
    return "".join(
        f"({value.numerator} : {value.denominator})\n" for value in abscissas
    )


def _gp_global_data(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction, Fraction],
) -> dict[str, object]:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP executable 'gp' is required")
    model = ",".join(
        f"({value.numerator}/{value.denominator})" for value in coefficients
    )
    program = f"""
default(realprecision,100);
E=ellinit([{model}]);v=0;M=ellminimalmodel(E,&v);N=ellglobalred(M)[1];
print("MODEL=",[M.a1,M.a2,M.a3,M.a4,M.a6]);
print("CONDUCTOR=",N);
print("DISCRIMINANT=",M.disc);
print("ROOT_NUMBER=",ellrootno(M));
print("VERSION=",version());
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout + completed.stderr:
        raise RuntimeError(completed.stdout + completed.stderr)
    values = dict(
        line.split("=", 1) for line in completed.stdout.splitlines() if "=" in line
    )
    return {
        "version": values["VERSION"],
        "minimal_model": ast.literal_eval(values["MODEL"]),
        "conductor": values["CONDUCTOR"],
        "minimal_discriminant": values["DISCRIMINANT"],
        "root_number": int(values["ROOT_NUMBER"]),
    }


def build_fermigier_rank20_manifest(
    ratpoints_output: str,
    *,
    maximum_reduction_prime: int = FERMIGIER_RANK20_MAXIMUM_REDUCTION_PRIME,
) -> dict[str, object]:
    specialization = specialize_fermigier_rank_sections(
        FERMIGIER_RANK20_PARAMETER
    )
    searched_points = parse_ratpoints_output(
        specialization.quartic_model, ratpoints_output
    )
    abscissas: list[Fraction] = []
    for point in searched_points:
        if not abscissas or point[0] != abscissas[-1]:
            abscissas.append(point[0])
    canonical_output = canonical_ratpoints_output(tuple(abscissas))
    if canonical_output != ratpoints_output:
        raise ValueError("ratpoints output is not in canonical quiet format")

    cloud = section_and_point_cloud_differences(
        specialization, searched_points
    )
    selected_indices, certificate = select_independent_subset(
        specialization.canonical_model,
        cloud,
        relation_prime=FERMIGIER_RANK20_RELATION_PRIME,
        maximum_reduction_prime=maximum_reduction_prime,
    )
    if len(selected_indices) != 20 or not set(range(12)).issubset(selected_indices):
        raise ArithmeticError("the pinned point cloud did not reproduce rank at least 20")

    common_denominator = 1_562_500
    integer_polynomial = [
        int(coefficient * common_denominator**2)
        for coefficient in specialization.quartic_model.quartic
    ]
    global_data = _gp_global_data(specialization.canonical_model)
    conductor = int(global_data["conductor"])
    with localcontext() as context:
        context.prec = 80
        log_conductor = str(Decimal(conductor).ln())
    return {
        "schema": "elliptic-curves.fermigier-rank20-near-miss.v1",
        "claim_level": "exact_rank_lower_bound_and_conductor",
        "generator": "elliptic-curves/scripts/run_fermigier_rank20_near_miss.py",
        "canonical_pinned_command": (
            "python3 elliptic-curves/scripts/run_fermigier_rank20_near_miss.py "
            "--output artifacts/generated-results/elliptic-curves/"
            "fermigier_rank20_near_miss_v1.json"
        ),
        "family": {
            "adapter_parameter": _fraction_text(FERMIGIER_RANK20_PARAMETER),
            "literal_shift": _fraction_text(2 * FERMIGIER_RANK20_PARAMETER),
            "canonical_model": [
                _fraction_text(value) for value in specialization.canonical_model
            ],
        },
        "bounded_search": {
            "engine": "ratpoints 2.1.3",
            "height_bound": FERMIGIER_RANK20_SEARCH_HEIGHT,
            "denominator_bound": FERMIGIER_RANK20_DENOMINATOR_BOUND,
            "quartic_coefficient_common_denominator": common_denominator,
            "integer_polynomial_coefficients_low_to_high": integer_polynomial,
            "options": ["-du", "13000", "-q", "-y"],
            "reported_runtime_seconds_on_search_host": "4.15",
            "abscissa_count": len(abscissas),
            "abscissas": [_fraction_text(value) for value in abscissas],
            "canonical_output_sha256": hashlib.sha256(
                canonical_output.encode()
            ).hexdigest(),
            "interpretation": "bounded search only; absence of further points is not proved",
        },
        "point_cloud": {
            "signed_quartic_point_count": len(searched_points),
            "deduplicated_difference_count": len(cloud),
            "maximum_reduction_prime": maximum_reduction_prime,
            "selected_indices": list(selected_indices),
            "selected_count": len(selected_indices),
            "all_twelve_baseline_differences_selected": True,
            "certificate": certificate.to_json_object(),
        },
        "global_curve": {
            **global_data,
            "log_conductor": log_conductor,
            "strict_log_target": "182.72",
            "maximum_integer_conductor_under_strict_target": str(
                STRICT_LOG_CONDUCTOR_MAXIMUM_INTEGER
            ),
            "below_strict_log_target": (
                conductor <= STRICT_LOG_CONDUCTOR_MAXIMUM_INTEGER
            ),
        },
        "conclusion": (
            "this specialization has at least 20 independent rational points "
            "and log conductor below 182.72; it is one point short of the "
            "rank-at-least-21 target"
        ),
        "limitations": {
            "target_status": "near miss, not a solution",
            "rank": "no twenty-first independent point and no upper bound",
            "saturation": "the selected subgroup is not claimed saturated",
            "normalization": "the Fermigier factor-two source discrepancy remains unresolved",
        },
        "randomness": "none; the bounded search and reduction scans are deterministic",
    }
