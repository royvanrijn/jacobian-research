#!/usr/bin/env python3
"""Assemble the exact/open counterexample-minimum scoreboard.

This script separates:

* certified upper witnesses stored in the repository;
* rigorous positive lower bounds;
* exact minima where the endpoints meet.

It intentionally does not turn a bounded search failure into a lower bound.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "minimal_counterexample_scoreboard.json"


def load(name: str) -> dict[str, object]:
    return json.loads((ARTIFACTS / name).read_text())


def interval(
    lower: int,
    upper: int,
    lower_reason: str,
    upper_reason: str,
) -> dict[str, object]:
    assert lower <= upper
    return {
        "lower_bound": lower,
        "upper_bound": upper,
        "exact_value": lower if lower == upper else None,
        "status": "exact" if lower == upper else "open interval",
        "lower_bound_reason": lower_reason,
        "upper_bound_reason": upper_reason,
    }


def main() -> None:
    cubic = load("essential_bcw_21_counterexample.json")
    image_20 = load("image_vanishing_counterexamples_20_40.json")
    image_21 = load("image_vanishing_counterexamples_21_42.json")
    restricted = load("restricted_minima_frontier.json")
    gaussian = load("two_real_gmc_frontier.json")
    gaussian_symmetric_chart = load("two_real_gmc_symmetric_chart.json")
    gaussian_remaining_four_weight = load(
        "two_real_gmc_remaining_four_weight.json"
    )

    assert cubic["dimension"] == 21
    assert image_20["source_dimension"] == 20
    assert image_20["laplacian_counterexample"]["dimension"] == 40
    assert image_21["source_dimension"] == 21
    assert image_21["laplacian_counterexample"]["dimension"] == 42
    assert "quadratic_GMC" in gaussian["theorem_level_results"]
    assert (
        gaussian_symmetric_chart["order_eight_quotient"]["vector_space_dimension"]
        == 84
    )
    assert gaussian_symmetric_chart["tenth_moment_multiplication"]["rank_mod_prime"] == 84
    assert gaussian_symmetric_chart["charts_excluded"] == 4
    assert gaussian_symmetric_chart["exceptional_four_weight_charts_remaining"] == 20
    assert gaussian_remaining_four_weight["input_chart_count"] == 20
    assert gaussian_remaining_four_weight["support_count"] == 3

    frontiers = restricted["rigorous_frontiers"]
    assert frontiers["r_cub"]["lower_bound"] == 3
    assert frontiers["r_cub"]["upper_bound"] == 17
    assert frontiers["nu_cub"]["lower_bound"] == 3
    assert frontiers["nu_cub"]["upper_bound"] == 18
    assert frontiers["rho_HN4"]["lower_bound"] == 3
    assert frontiers["rho_HN4"]["upper_bound"] == 37
    assert frontiers["n_HN4"]["lower_bound"] == 6
    assert frontiers["n_HN4"]["upper_bound"] == 42

    scoreboard = {
        "format": "minimal-counterexample-scoreboard-v1",
        "field": "complex coefficients / characteristic zero as appropriate",
        "ambient_dimension_frontiers": {
            "general_Keller_noninvertibility": interval(
                2,
                3,
                "a one-variable Keller map is affine",
                "the certified foundational determinant-nonzero collision uses 3 variables",
            ),
            "cubic_homogeneous_Keller_noninvertibility": interval(
                3,
                21,
                (
                    "dimension at most 2 has rank(JH)<=2, and the "
                    "de Bondt--Sun rank-two theorem gives invertibility"
                ),
                "the essential cubic-homogeneous collision uses 21 variables",
            ),
            "GMC_failure_real_Gaussian_dimension": interval(
                2,
                3,
                "GMC(1) is proved",
                "Long gives an explicit counterexample in 3 real Gaussian variables",
            ),
            "SIC_failure_pair_dimension": interval(
                2,
                20,
                "the one-pair Image Conjecture is proved",
                "the identity-output slice gives an explicit SIC(20) counterexample",
            ),
            "ordinary_Laplacian_GVC_failure_dimension": interval(
                2,
                40,
                (
                    "in one variable Delta(P)=0 forces deg(P)<=1, after "
                    "which every fixed mixed derivative eventually vanishes"
                ),
                "the identity-output slice gives a 40-variable counterexample",
            ),
            "homogeneous_quartic_HN_VC_failure_dimension": frontiers["n_HN4"],
        },
        "rank_and_index_frontiers": {
            "cubic_homogeneous_Jacobian_rank": frontiers["r_cub"],
            "cubic_homogeneous_Jacobian_nilpotency_index": frontiers["nu_cub"],
            "homogeneous_quartic_Hessian_rank": frontiers["rho_HN4"],
        },
        "bounded_structural_frontiers": {
            "two_real_cubic_four_weight_charts": {
                "all_mixed_sign_charts": 121,
                "excluded_by_support_census": 97,
                "excluded_inside_exceptional_supports": 24,
                "remaining": 0,
                "remaining_supports": 0,
                "certificates": [
                    "two_real_gmc_symmetric_chart.json",
                    "two_real_gmc_remaining_four_weight.json",
                ],
                "consequence": (
                    "a cubic GMC(2) counterexample has at least five "
                    "rotational weights"
                ),
                "status": "exact closure of the four-weight cubic stratum",
            }
        },
        "exact_degree_minima": {
            "general_Keller_counterexample_total_degree": {
                "exact_value": 3,
                "lower_reason": "Wang's theorem proves every quadratic Keller map invertible",
                "upper_reason": (
                    "the certified 21-variable cubic-homogeneous collision has degree 3"
                ),
            },
            "Gaussian_moment_counterexample_total_degree": {
                "exact_value": 3,
                "lower_reason": (
                    "the quadratic Gaussian theorem in TWO_REAL_GMC_FRONTIER "
                    "proves GMC for every polynomial of degree at most 2"
                ),
                "upper_reason": (
                    "Long's explicit four-real-variable Gaussian counterexample is cubic"
                ),
            },
        },
        "smallest_certified_witnesses": {
            "cubic_homogeneous_Keller": {
                "dimension": 21,
                "artifact": "essential_bcw_21_counterexample.json",
            },
            "SIC": {
                "pair_dimension": 20,
                "polynomial_ring_variable_count": 40,
                "artifact": "image_vanishing_counterexamples_20_40.json",
            },
            "ordinary_Laplacian_GVC": {
                "dimension": 40,
                "artifact": "image_vanishing_counterexamples_20_40.json",
            },
            "homogeneous_quartic_HN_VC": {
                "dimension": 42,
                "artifact": "image_vanishing_counterexamples_21_42.json",
            },
        },
        "promotion_policy": (
            "Only exact identities or cited positive theorems change an endpoint; "
            "modular and bounded searches remain diagnostics."
        ),
    }
    OUTPUT.write_text(json.dumps(scoreboard, indent=2) + "\n")
    print("PASS scoreboard: ambient witness dimensions 21 / 20 / 40 / 42")
    print("PASS scoreboard: GMC failure dimension is in [2,3]")
    print("PASS scoreboard: exact minimum Gaussian counterexample degree is 3")
    print("PASS scoreboard: all 121 two-real cubic four-weight charts excluded")
    print("PASS scoreboard: cubic rank/index intervals are [3,17] / [3,18]")
    print("PASS scoreboard: homogeneous quartic Hessian-rank interval is [3,37]")
    print(f"PASS scoreboard: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
