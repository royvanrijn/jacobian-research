#!/usr/bin/env python3
"""Cross-check family metadata against exact code and a PARI specialization."""

from __future__ import annotations

import ast
import json
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
    FERMIGIER_E22_RECONSTRUCTION_SHIFT,
    FERMIGIER_REPORTED_PARAMETER,
    FERMIGIER_ROOTS,
    fermigier_canonical_coefficients,
    fermigier_quartic,
    thirteenth_visible_point,
)


def evaluate(coefficients: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def verify_calibration_family() -> None:
    path = PROGRAM_ROOT / "families" / "calibration_rank_two.json"
    family = json.loads(path.read_text())
    assert family["id"] == "calibration-rank-at-least-two"
    assert family["generic_rank_lower_bound"] == 2
    assert [item["label"] for item in family["known_independent_sections"]] == [
        "Q",
        "R",
    ]
    for parameter in range(-7, 8):
        if parameter == 0:
            continue
        for x_coordinate, y_coordinate in (
            (parameter, parameter),
            (-parameter, parameter),
            (0, parameter),
            (1, -1),
        ):
            assert y_coordinate**2 == (
                x_coordinate**3
                - parameter**2 * x_coordinate
                + parameter**2
            )
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input="setrand(1);E=ellinit([0,0,0,-25,25]);print(ellrank(E,2));\n",
        text=True,
        capture_output=True,
        check=True,
    )
    assert "***" not in completed.stdout + completed.stderr
    rank = ast.literal_eval(completed.stdout.strip())
    # PARI 2.15.4 returned Q,R while 2.17.4 returns R,Q for the same seeded
    # calculation.  Generator order is not mathematical output: retain the
    # exact rank bounds and require precisely the two pinned points.
    assert rank[:3] == [2, 2, 0]
    assert {tuple(point) for point in rank[3]} == {(5, 5), (-5, 5)}
    witness = family["independence_witness"]
    assert witness["specialized_model"] == [0, 0, 0, -25, 25]
    assert {tuple(point) for point in witness["specialized_points"]} == {
        tuple(point) for point in rank[3]
    }
    assert witness["pari_gp_setrand_seed"] == 1


def verify_fermigier_family() -> None:
    path = PROGRAM_ROOT / "families" / "fermigier_mestre_rank12.json"
    family = json.loads(path.read_text())
    assert tuple(family["fixed_roots"]) == FERMIGIER_ROOTS
    normalization = family["parameter_normalization"]
    assert normalization["paper_reported_e22_parameter"] == "t=19754/39"
    assert FERMIGIER_REPORTED_PARAMETER == Fraction(19754, 39)
    assert FERMIGIER_E22_RECONSTRUCTION_SHIFT == Fraction(39508, 39)
    model = family["canonical_weierstrass_model"]
    assert tuple(model["discriminant_coefficients_low_to_high"]) == (
        FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS
    )
    coefficient_data = model["coefficients_low_to_high"]
    for parameter in range(-6, 7):
        expected = fermigier_canonical_coefficients(parameter)
        actual = (
            Fraction(evaluate(coefficient_data["a1"], parameter)),
            Fraction(evaluate(coefficient_data["a2"], parameter)),
            Fraction(evaluate(coefficient_data["a3"], parameter)),
            Fraction(evaluate(coefficient_data["a4"], parameter)),
            Fraction(evaluate(coefficient_data["a6"], parameter)),
        )
        assert actual == expected
    assert family["generic_independence_reproduced_here"] is True
    construction = family["construction"]
    assert construction["source_parameters_for_fixed_roots"].startswith("(u,v)=(3,5)")
    assert construction["thirteenth_point"].startswith("x=1256/5-(17/35)*s")
    # Both sides of the thirteenth-point square identity have degree at most
    # twelve after its linear abscissa substitution.  Thirteen exact values
    # therefore check the identity coefficientwise without a CAS dependency.
    for shift in range(-6, 8):
        if shift:
            thirteenth_visible_point(fermigier_quartic(shift))
    assert "unresolved factor-two discrepancy" in normalization["status"]
    discrepancy = family["source_discrepancy_replay"]
    assert discrepancy["literal_printed_shift"] == "19754/39"
    assert discrepancy["reconstruction_shift"] == "39508/39"
    assert discrepancy["literal_minimal_model"][3:] == [
        -1223348097402005168062873899944,
        -213263015130965060475376699543914227367884158,
    ]
    benchmark = family["benchmark_specialization"]
    assert benchmark["adapter_parameter_u"] == "19754/39"
    assert benchmark["literal_shift_s"] == "39508/39"
    assert benchmark["rank_lower_bound_reproduced_here"] == 22
    assert benchmark["generic_section_difference_rank_reproduced_here"] == 12
    near_miss = family["rank20_near_miss"]
    assert near_miss["adapter_parameter_u"] == "28917/20"
    assert near_miss["rank_lower_bound_reproduced_here"] == 20
    assert near_miss["minimal_model"][3:] == [
        "-4437412060110743641525245114305",
        "3586842216822165612930264910099076801587288127",
    ]
    assert near_miss["target_status"].endswith("not a solution")


def main() -> None:
    verify_calibration_family()
    verify_fermigier_family()
    print(
        "PASS family data: calibration sections/rank-two specialization and "
        "Fermigier equation/discriminant/normalization metadata agree with code"
    )


if __name__ == "__main__":
    main()
