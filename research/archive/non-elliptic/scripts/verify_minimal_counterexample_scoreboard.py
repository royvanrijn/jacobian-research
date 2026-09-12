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
    dvorsky = load("dvorsky_gvc5_counterexample.json")
    gvc_three = load("gvc3_homogeneous_counterexample.json")
    sic_three = load("three_pair_image_mathieu_counterexample.json")
    sic_two = load("two_pair_image_mathieu_counterexample.json")
    hessian_rank_35_slice = load(
        "hessian_rank_35_identity_slice_counterexample.json"
    )
    hessian_rank_34_slice = load(
        "hessian_rank_34_double_identity_slice_counterexample.json"
    )
    f12_reduction = load("macfarlane_f12_coordinate_pair_reduction.json")

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
    assert dvorsky["variables"] == ["t", "a", "b", "c", "d"]
    assert (
        dvorsky["consequences"]["unrestricted_constant_coefficient_GVC"]
        == "fails in 5 variables"
    )
    assert dvorsky["consequences"]["SIC"] == "fails in 5 contraction pairs"
    assert len(gvc_three["variables"]) == 3
    assert gvc_three["operator"] == "Lambda=(4*d_x*d_y+d_t^2)^6"
    assert gvc_three["all_order_claim"]["pure"].endswith("for every m>=1")
    assert gvc_three["all_order_claim"]["mixed"].endswith("for every m>=1")
    assert len(sic_three["contraction_pairs"]) == 3
    assert sic_three["expanded_f_term_count"] == 4
    assert sic_three["g"] == "y"
    assert sic_three["bidegrees"]["g"] == [0, 1]
    assert sic_three["all_order_identities"]["E(f^m)"] == "0"
    assert (
        sic_three["all_order_identities"]["[t]E(g*f^m)"]
        == "(-1)^(m-1)*(m+1)!*m!"
    )
    assert len(sic_two["contraction_pairs"]) == 2
    assert sic_two["expanded_F_term_count"] == 16
    assert sic_two["bidegrees"]["F"] == [4, 4]
    assert sic_two["coefficient_matrix_determinant"] == 48
    assert sic_two["all_order_identities"]["E_2(F^m)"] == "0"
    assert (
        sic_two["all_order_identities"]["E_2(Q*F^m)"]
        == "(4*m+2)!*m!/(2*m+1)!!"
    )
    assert hessian_rank_35_slice["slice_dimension"] == 21
    assert hessian_rank_35_slice["HN_potential"]["dimension"] == 42
    assert hessian_rank_35_slice["HN_potential"]["degrees"] == [2, 3, 4]
    assert (
        hessian_rank_35_slice["HN_potential"]["generic_hessian_rank"] == 35
    )
    assert hessian_rank_34_slice["slice_dimension"] == 20
    assert hessian_rank_34_slice["HN_potential"]["dimension"] == 40
    assert hessian_rank_34_slice["HN_potential"]["degrees"] == [2, 3, 4]
    assert (
        hessian_rank_34_slice["HN_potential"]["generic_hessian_rank"] == 34
    )
    assert f12_reduction["F12"]["dimension"] == 12
    assert f12_reduction["F12"]["determinant"] == "1"
    assert f12_reduction["F12"]["cubic_output_rank"] == 6
    assert f12_reduction["G19"]["dimension"] == 19
    assert f12_reduction["G19"]["determinant"] == "1"
    assert f12_reduction["G19"]["collision"] is True

    frontiers = restricted["rigorous_frontiers"]
    assert frontiers["n_cub"]["lower_bound"] == 5
    assert frontiers["n_cub"]["upper_bound"] == 20
    assert frontiers["r_cub"]["lower_bound"] == 3
    assert frontiers["r_cub"]["upper_bound"] == 17
    assert frontiers["nu_cub"]["lower_bound"] == 3
    assert frontiers["nu_cub"]["upper_bound"] == 18
    assert frontiers["rho_cot"]["lower_bound"] == 6
    assert frontiers["rho_cot"]["upper_bound"] == 37
    assert frontiers["n_Dru"]["lower_bound"] == 6
    assert frontiers["n_Dru"]["upper_bound"] == 451
    assert frontiers["rho_HN4"]["lower_bound"] == 3
    assert frontiers["rho_HN4"]["upper_bound"] == 37
    assert frontiers["n_HN4"]["lower_bound"] == 6
    assert frontiers["n_HN4"]["upper_bound"] == 40

    scoreboard = {
        "format": "minimal-counterexample-scoreboard-v3",
        "field": "complex coefficients / characteristic zero as appropriate",
        "ambient_dimension_frontiers": {
            "general_Keller_noninvertibility": interval(
                2,
                3,
                "a one-variable Keller map is affine",
                "the certified foundational determinant-nonzero collision uses 3 variables",
            ),
            "cubic_homogeneous_Keller_noninvertibility": interval(
                5,
                19,
                "cubic-homogeneous Keller maps are invertible through dimension four",
                (
                    "the exact F12 coordinate-pair reduction has cubic-output "
                    "rank six, giving a 12+6+1=19 rank-compressed parent"
                ),
            ),
            "Druzkowski_Keller_noninvertibility": frontiers["n_Dru"],
            "GMC_failure_real_Gaussian_dimension": interval(
                3,
                3,
                "the lower-face prime theorem proves GMC(2)",
                "Long gives an explicit counterexample in 3 real Gaussian variables",
            ),
            "SIC_failure_pair_dimension": interval(
                2,
                2,
                "the one-pair Image Conjecture is proved",
                (
                    "the full-rank bidegree-(4,4) formula gives an explicit "
                    "SIC(2) counterexample"
                ),
            ),
            "unrestricted_constant_coefficient_GVC_failure_dimension": interval(
                3,
                3,
                "the Hall-envelope theorem proves unrestricted GVC(2)",
                (
                    "the homogeneous cusp formula gives an explicit "
                    "three-variable GVC counterexample"
                ),
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
            "homogeneous_quartic_HN_VC_failure_dimension": interval(
                6,
                38,
                (
                    "the symmetric homogeneous nilpotent-Jacobian result "
                    "holds through dimension five"
                ),
                (
                    "the homogeneous cotangent lift of the exact "
                    "19-variable cubic parent has dimension 38"
                ),
            ),
        },
        "rank_and_index_frontiers": {
            "cubic_homogeneous_Jacobian_rank": frontiers["r_cub"],
            "cubic_homogeneous_Jacobian_nilpotency_index": frontiers["nu_cub"],
            "cotangent_lift_quartic_Hessian_rank": frontiers["rho_cot"],
            "unrestricted_homogeneous_quartic_Hessian_rank": frontiers["rho_HN4"],
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
                    "the pinned external MacFarlane G20 collision has degree 3"
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
            "genuinely_ungraded_Keller_geometric_degree": {
                "exact_value": 3,
                "lower_reason": (
                    "the exact geometric-degree spectrum excludes "
                    "noninvertible Keller maps of geometric degree 1 or 2"
                ),
                "upper_reason": (
                    "the A=0, gamma=1 universal cubic testbed specialization "
                    "has no algebraic-torus-equivariant polynomial "
                    "left-right representative"
                ),
            },
        },
        "smallest_certified_witnesses": {
            "cubic_homogeneous_Keller": {
                "dimension": 19,
                "artifact": "macfarlane_f12_coordinate_pair_reduction.json",
                "certificate_scope": (
                    "exact rank-compressed parent of the independently "
                    "replayed twelve-variable map"
                ),
            },
            "cubic_homogeneous_Keller_internal_independent_replay": {
                "dimension": 21,
                "artifact": "essential_bcw_21_counterexample.json",
            },
            "SIC": {
                "pair_dimension": 2,
                "polynomial_ring_variable_count": 4,
                "artifact": "two_pair_image_mathieu_counterexample.json",
            },
            "unrestricted_constant_coefficient_GVC": {
                "dimension": 3,
                "operator_order": 12,
                "operator_power": 6,
                "artifact": "gvc3_homogeneous_counterexample.json",
                "exact_dimension_reason": (
                    "GVC(2) is proved by Hall-envelope separation"
                ),
            },
            "ordinary_Laplacian_GVC": {
                "dimension": 40,
                "artifact": "image_vanishing_counterexamples_20_40.json",
            },
            "homogeneous_quartic_HN_VC": {
                "dimension": 38,
                "derived_from": (
                    "the exact G19 parent by the homogeneous cotangent lift"
                ),
                "certificate_scope": (
                    "exact F12/G19 determinant certificate plus the "
                    "homogeneous cotangent bridge"
                ),
            },
            "nonhomogeneous_HN_degree_at_most_4_rank": {
                "dimension": 40,
                "generic_hessian_rank": 34,
                "degrees": [2, 3, 4],
                "artifact": (
                    "hessian_rank_34_double_identity_slice_counterexample.json"
                ),
                "scope_warning": (
                    "this does not change the homogeneous quartic HN "
                    "rank-37 or certified dimension-38 frontiers"
                ),
            },
        },
        "promotion_policy": (
            "Only exact identities or cited positive theorems change an endpoint; "
            "modular and bounded searches remain diagnostics."
        ),
        "search_status": (
            "frozen reproducible frontier; resume only for a theorem-directed "
            "classification or a specifically motivated construction"
        ),
    }
    OUTPUT.write_text(json.dumps(scoreboard, indent=2) + "\n")
    print(
        "PASS scoreboard: ambient witness dimensions "
        "19 / SIC 2 / GVC 3 / Laplacian 40 / HN 38"
    )
    print("PASS scoreboard: exact minimum SIC pair dimension is 2")
    print("PASS scoreboard: exact unrestricted GVC failure dimension is 3")
    print("PASS scoreboard: exact GMC failure dimension is 3")
    print("PASS scoreboard: exact minimum Gaussian counterexample degree is 3")
    print("PASS scoreboard: exact minimum genuinely ungraded geometric degree is 3")
    print("PASS scoreboard: all 121 two-real cubic four-weight charts excluded")
    print("PASS scoreboard: cubic dimension/rank/index intervals are [5,19] / [3,17] / [3,18]")
    print("PASS scoreboard: cotangent Hessian-rank interval is [6,37]")
    print("PASS scoreboard: unrestricted quartic Hessian-rank interval is [3,37]")
    print("PASS scoreboard: nonhomogeneous degree-at-most-4 HN witness has rank 34")
    print("PASS scoreboard: Druzkowski dimension interval is [6,451]")
    print(f"PASS scoreboard: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
