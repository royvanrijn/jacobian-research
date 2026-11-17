#!/usr/bin/env python3
"""Certify an exact-rank-two two-separated-row SIC exclusion.

On the dense support

    {0} x {1,2,3,4}  union  {4} x {0,1,2,3},

the contraction-preserving two-dimensional torus normalizes c_01=c_43=1.
Saturating the six residual coefficients makes c_04*c_40 nonzero, hence
the coefficient matrix has exact rank two.  Over QQ, moments through order
seven give a zero-dimensional scheme of degree 604.  Adding mu_8 gives
the unit ideal.  Thus every exact finite component on this dense chart is
killed by the eighth moment, before any all-order recurrence is needed.

Coordinate boundaries are not classified here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from research_two_pair_sic_bidegree44_two_row_channel import (
    RESIDUALS,
    SUPPORT,
    moment,
    solve,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_two_row_channel.json"
)


def compact_record(record: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in record.items()
        if key not in {"result", "stdout_tail", "stderr_tail"}
    }


def validate_existing_artifact(path: Path) -> None:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("format") != "two-pair-sic-bidegree44-two-row-channel-v1":
        raise AssertionError("unexpected two-row channel artifact format")
    if artifact.get("coefficient_support") != [list(position) for position in SUPPORT]:
        raise AssertionError("stored dense support does not match the exact chart")

    factor_chart = artifact.get("factor_chart")
    if not isinstance(factor_chart, dict):
        raise AssertionError("stored artifact has no factor-chart metadata")
    if factor_chart.get("residual_coordinates") != [str(value) for value in RESIDUALS]:
        raise AssertionError("stored residual coordinates do not match the chart")
    if factor_chart.get("coefficient_torus_saturation") != (
        "z0*z1*z2*z3*z4*z5 != 0"
    ):
        raise AssertionError("stored chart does not retain full torus saturation")
    if factor_chart.get("exact_rank_two_minor") != (
        "det C[rows 0,4; columns 0,4]=-z2*z3"
    ):
        raise AssertionError("stored chart has the wrong exact-rank-two minor")

    moments = [moment(order) for order in range(1, 9)]
    if not moments[0].is_zero or any(value.is_zero for value in moments[1:]):
        raise AssertionError("unexpected dense-chart moment profile")
    profiles = [
        [order, value.total_degree(), len(value.terms())]
        for order, value in enumerate(moments, start=1)
    ]
    seven = artifact.get("through_seven")
    if (
        not isinstance(seven, dict)
        or seven.get("returncode") != 0
        or seven.get("through") != 7
        or seven.get("status") != "zero_dimensional"
        or seven.get("variables") != 7
        or seven.get("scheme_degree") != 604
        or seven.get("profiles") != profiles[:7]
    ):
        raise AssertionError("stored seven-moment degree-604 fibre is incomplete")
    eight = artifact.get("through_eight")
    if (
        not isinstance(eight, dict)
        or eight.get("returncode") != 0
        or eight.get("through") != 8
        or eight.get("status") != "unit_ideal"
        or eight.get("profiles") != profiles
    ):
        raise AssertionError("stored eighth-moment unit certificate is incomplete")

    print("PASS stored dense support, torus saturation, and rank minor match the chart")
    print("PASS stored degree-604 prefix and mu_8 unit outcomes match exact profiles")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help=(
            "validate the stored chart, moment profiles, and solver outcomes "
            "without rerunning msolve"
        ),
    )
    arguments = parser.parse_args()
    if arguments.audit_existing_only:
        validate_existing_artifact(arguments.output)
        return

    # The first moment vanishes identically because this support has no
    # diagonal coefficient.  The next moments are genuine quotient
    # equations, with primitive integer coefficients.
    moments = [moment(order) for order in range(1, 9)]
    assert moments[0].is_zero
    assert all(not value.is_zero for value in moments[1:])
    assert len(RESIDUALS) == 6

    seven = solve(7, 600, 4, rational_parametrization=True)
    assert seven["returncode"] == 0
    assert seven["status"] == "zero_dimensional"
    assert seven["variables"] == 7
    assert seven["scheme_degree"] == 604

    eight = solve(8, 600, 4, rational_parametrization=False)
    assert eight["returncode"] == 0
    assert eight["status"] == "unit_ideal"

    artifact = {
        "format": "two-pair-sic-bidegree44-two-row-channel-v1",
        "field": "characteristic zero",
        "coefficient_support": [list(position) for position in SUPPORT],
        "factor_chart": {
            "U": "rows e_0 and e_4",
            "B": "two coefficient rows on columns 1..4 and 0..3",
            "internal_GL2_gauge": "already quotiented by the fixed U pivot",
            "orbit_normalization": "c_01=c_43=1",
            "residual_coordinates": [str(value) for value in RESIDUALS],
            "coefficient_torus_saturation": "z0*z1*z2*z3*z4*z5 != 0",
            "exact_rank_two_minor": "det C[rows 0,4; columns 0,4]=-z2*z3",
        },
        "moment_formula": (
            "mu_m=sum_I I!*(4m-I)!*[x^I*y^I]Phi_C(x,y)^m"
        ),
        "through_seven": compact_record(seven),
        "through_eight": compact_record(eight),
        "component_conclusion": (
            "the exact QQ seven-moment coefficient-torus scheme is finite "
            "of degree 604, and mu_8 removes every component"
        ),
        "recurrence_gate": (
            "no component survives the eighth pure moment, so there is no "
            "all-order candidate on which creative telescoping is required"
        ),
        "conclusion": (
            "the dense two-separated-row exact-rank-two chart contains no "
            "all-order pure-moment point and hence no SIC counterexample"
        ),
        "scope": (
            "exact dense coefficient-torus exclusion; coordinate boundaries "
            "and other rank-two factor charts remain open"
        ),
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )

    print("PASS direct two-row chart has exact coefficient rank two")
    print("PASS moments through seven give an exact degree-604 scheme over QQ")
    print("PASS adjoining mu_8 gives the unit ideal over QQ")
    print("PASS dense two-separated-row rank-two chart is SIC-safe")


if __name__ == "__main__":
    main()
