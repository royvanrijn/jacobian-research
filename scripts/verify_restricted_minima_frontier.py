#!/usr/bin/env python3
"""Verify and serialize the present restricted-minima frontier.

This checker keeps three logically different statements separate:

1. exact invariants of the certified repository witnesses;
2. rigorous lower bounds for the unrestricted restricted classes;
3. narrower positive theorems (for example the power-linear index-three
   theorem) which must not be promoted to the full cubic-homogeneous class.

The output is a machine-readable research ledger, not a claim that any of the
open minima has already been determined.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "restricted_minima_frontier.json"


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
    cubic = load("low_complexity_bcw_21_counterexample.json")
    source = load("essential_bcw_21_counterexample.json")
    index_reduced = load("index_reduced_bcw_22_counterexample.json")
    rank_reduced = load("rank_reduced_bcw_24_counterexample.json")
    hessian_rank_reduced = load(
        "hessian_rank_reduced_bcw_22_counterexample.json"
    )
    xxs_search = load("restricted_bcw_circuit_search_xxs_w32_d28.json")
    yvyb_search = load(
        "restricted_bcw_circuit_search_yvyb_structural_w22.json"
    )
    shared_yb_search = load(
        "restricted_bcw_circuit_search_ayb_yvyb_w36.json"
    )
    shared_vz_structural_search = load(
        "restricted_bcw_circuit_search_xvvz_v2h_structural_w25.json"
    )
    shared_vz_mixed_search = load(
        "restricted_bcw_circuit_search_xvvz_v2h_mixed_w25.json"
    )
    exceptional_as_search = load(
        "restricted_bcw_circuit_search_aspert_m12.json"
    )
    rank37_perturbation_search = load(
        "rank37_gate_perturbation_search.json"
    )
    kernel_excess = load("cotangent_kernel_excess_frontier.json")
    index_three_model = load("index_three_inverse_degree_model.json")
    index_three_degree_counterexample = load(
        "index_three_degree_bound_counterexample.json"
    )
    index_three_trees = load("index_three_inverse_tree_obstruction.json")
    index_three_normal_form = load(
        "index_three_rank_normal_form_exclusion.json"
    )
    xxs_rank_hybrid = load(
        "restricted_bcw_circuit_search_xxs_rank_hybrid_w24.json"
    )
    xxs_v2r_hybrid = load("restricted_bcw_circuit_search_xxs_v2r_w16.json")
    xxs_y2vb_hybrid = load(
        "restricted_bcw_circuit_search_xxs_y2vb_w16.json"
    )
    all_circuits_search = load(
        "restricted_bcw_circuit_search_all_w64.json"
    )
    assert cubic["dimension"] == source["dimension"] == 21
    statistics = cubic["statistics"]
    assert statistics["generic_rank_JH_over_QQ_x"] == 18
    assert statistics["nilpotency_index_JH"] == 19
    assert statistics["exact_point_power_ranks"] == list(range(18, -1, -1))
    assert len(cubic["collision_points"]) == 3
    assert len({tuple(point) for point in cubic["collision_points"]}) == 3
    assert source["jacobian_determinant"] == "1"
    assert index_reduced["dimension"] == 22
    assert index_reduced["statistics"]["generic_rank_JH_over_QQ_x"] == 18
    assert index_reduced["statistics"]["nilpotency_index_JH"] == 18
    assert index_reduced["jacobian_determinant"] == "1"
    assert len({tuple(point) for point in index_reduced["collision_points"]}) == 3
    assert rank_reduced["dimension"] == 24
    assert rank_reduced["statistics"]["generic_rank_JH_over_QQ_x"] == 17
    assert rank_reduced["statistics"]["nilpotency_index_JH"] == 18
    assert rank_reduced["jacobian_determinant"] == "1"
    assert len({tuple(point) for point in rank_reduced["collision_points"]}) == 3
    assert hessian_rank_reduced["dimension"] == 22
    assert hessian_rank_reduced["homogeneous_quartic_HN_lift"]["dimension"] == 44
    assert hessian_rank_reduced["homogeneous_quartic_HN_lift"]["degree"] == 4
    assert hessian_rank_reduced["statistics"]["generic_rank_JH_over_QQ_x"] == 18
    assert hessian_rank_reduced["statistics"]["nilpotency_index_JH"] == 18
    assert (
        hessian_rank_reduced["statistics"][
            "generic_cotangent_hessian_rank_over_QQ_xy"
        ]
        == 37
    )
    assert (
        hessian_rank_reduced["statistics"][
            "transformed_hessian_sampled_nilpotency_index"
        ]
        == 35
    )
    assert hessian_rank_reduced["jacobian_determinant"] == "1"
    assert (
        len(
            {
                tuple(point)
                for point in hessian_rank_reduced["collision_points"]
            }
        )
        == 3
    )
    for search in (
        xxs_search,
        yvyb_search,
        shared_yb_search,
        shared_vz_structural_search,
        shared_vz_mixed_search,
        exceptional_as_search,
        xxs_rank_hybrid,
        xxs_v2r_hybrid,
        xxs_y2vb_hybrid,
        all_circuits_search,
    ):
        assert search["format"] == "restricted-bcw-circuit-search-v2"
        assert not any(
            terminal["beats_a_current_incumbent_diagnostically"]
            for terminal in search["pareto_terminals"]
        )
    assert sum(row["count"] for row in xxs_search["terminal_histogram"]) == 244
    assert xxs_search["pareto_terminals"][0]["objective"][:4] == [
        19,
        18,
        38,
        23,
    ]
    assert (
        sum(row["count"] for row in yvyb_search["terminal_histogram"])
        == 29
    )
    assert yvyb_search["pareto_terminals"][0]["objective"][:4] == [
        18,
        18,
        38,
        23,
    ]
    assert (
        sum(row["count"] for row in shared_yb_search["terminal_histogram"])
        == 46
    )
    assert [
        terminal["objective"][:4]
        for terminal in shared_yb_search["pareto_terminals"]
    ] == [[18, 18, 40, 25], [18, 18, 41, 24]]
    for search, count in (
        (shared_vz_structural_search, 47),
        (shared_vz_mixed_search, 35),
    ):
        assert (
            sum(row["count"] for row in search["terminal_histogram"])
            == count
        )
        assert search["pareto_terminals"][0]["objective"][:4] == [
            18,
            18,
            41,
            23,
        ]
        assert search["search_scope"]["required_atoms"] == ["v2h", "xvvz"]
    assert (
        sum(
            row["count"]
            for row in exceptional_as_search["terminal_histogram"]
        )
        == 33
    )
    assert exceptional_as_search["pareto_terminals"][0]["objective"][:4] == [
        18,
        18,
        41,
        23,
    ]
    assert exceptional_as_search["search_scope"]["required_atoms"] == [
        "aspert_m12",
        "qb",
        "x2s",
    ]
    assert rank37_perturbation_search["format"] == (
        "rank37-gate-perturbation-search-v1"
    )
    assert rank37_perturbation_search["coefficients"] == [
        -12,
        -6,
        -3,
        -1,
        1,
        3,
        6,
        12,
    ]
    assert (
        sum(
            row["count"]
            for row in rank37_perturbation_search["terminal_histogram"]
        )
        == 32
    )
    assert rank37_perturbation_search["pareto_terminals"][0]["objective"] == [
        18,
        18,
        37,
        22,
        35,
    ]
    assert len(rank37_perturbation_search["coefficient_family_leaders"]) == 8
    assert all(
        row["objective"] == [18, 18, 37, 22, 35]
        for row in rank37_perturbation_search["coefficient_family_leaders"]
    )
    assert kernel_excess["format"] == "cotangent-kernel-excess-frontier-v1"
    assert [
        (
            row["certified_generic_rank_JH"],
            row["sampled_kernel_excess"],
            row["sampled_cotangent_hessian_rank"],
        )
        for row in kernel_excess["profiles"]
    ] == [(18, 2, 38), (18, 2, 38), (17, 4, 38), (18, 1, 37)]
    assert kernel_excess["rank_index_coupling"]["current_saturation"] == (
        "the rank-17 witness has index 18, so it saturates index=rank+1"
    )
    assert index_three_model["format"] == (
        "cubic-index-three-inverse-degree-model-v1"
    )
    assert index_three_model["nilpotency_index_JH"] == 3
    assert index_three_model["inverse_degree"] == 9
    assert index_three_model["nonzero_inverse_correction_degrees"] == [
        3,
        5,
        7,
        9,
    ]
    assert index_three_degree_counterexample["format"] == (
        "cubic-index-three-degree-bound-counterexample-v1"
    )
    assert index_three_degree_counterexample["dimension"] == 5
    assert index_three_degree_counterexample["generic_rank_JH"] == 3
    assert (
        index_three_degree_counterexample["weak_nilpotency_index_JH"] == 3
    )
    assert index_three_degree_counterexample["inverse_degree"] == 13
    assert index_three_degree_counterexample["omega_11_evaluation"] == (
        "(-4*x4*x5^10,0,0,0,0)"
    )
    assert index_three_trees["format"] == (
        "cubic-index-three-inverse-tree-obstruction-v1"
    )
    assert index_three_trees["inverse_tree_counts_by_degree"] == {
        "3": 1,
        "5": 1,
        "7": 2,
        "9": 4,
        "11": 8,
    }
    assert index_three_trees["universal_relation_counts_by_degree"] == {
        "7": 1,
        "9": 2,
        "11": 7,
    }
    assert not index_three_trees["degree_11_lies_in_relation_span"]
    assert index_three_trees["degree_11_quotient"]["relation_rank"] == 5
    assert index_three_trees["degree_11_quotient"]["quotient_dimension"] == 3
    assert index_three_trees["degree_11_quotient"][
        "target_normal_form_formula"
    ] == {
        "C(N(H),H,H)": "-1/2",
        "B(B(H,H),H)": "-1/2",
        "N(C(H,H,H))": "-1/6",
    }
    assert index_three_trees["van_den_essen_realization"][
        "degree_11_evaluation"
    ] == "(-4*y4*y5^10,0,0,0,0)"
    assert index_three_normal_form["format"] == (
        "cubic-index-three-rank-normal-form-exclusion-v1"
    )
    assert index_three_normal_form["index_three_condition"] == "lambda=0"
    assert index_three_normal_form["index_three_locus_rank_bound"] == 2
    assert (
        index_three_normal_form["index_three_locus_inverse_degree_bound"]
        == 5
    )
    for search, count, objective, excess in (
        (xxs_rank_hybrid, 48, [22, 19, 45, 28], 1),
        (xxs_v2r_hybrid, 32, [21, 19, 44, 27], 2),
        (xxs_y2vb_hybrid, 34, [19, 18, 43, 26], 5),
    ):
        assert (
            sum(row["count"] for row in search["terminal_histogram"])
            == count
        )
        assert search["pareto_terminals"][0]["objective"][:4] == objective
        assert (
            search["pareto_terminals"][0]["profile"][
                "cotangent_kernel_excess_mod_1000003"
            ]
            == excess
        )
    assert all_circuits_search["search_scope"]["width"] == 64
    assert all_circuits_search["search_scope"]["prebeam_factor"] == 4
    assert all_circuits_search["search_scope"]["max_steps"] == 24
    assert len(all_circuits_search["search_scope"]["circuit_atoms"]) == 9
    assert (
        sum(
            row["count"]
            for row in all_circuits_search["terminal_histogram"]
        )
        == 140
    )
    assert [
        terminal["objective"][:4]
        for terminal in all_circuits_search["pareto_terminals"]
    ] == [
        [17, 18, 38, 24],
        [17, 18, 39, 23],
        [18, 18, 37, 22],
    ]

    frontier = {
        "format": "restricted-jacobian-minima-frontier-v1",
        "field": "characteristic zero",
        "definitions": {
            "r_cub": (
                "minimum generic rank(JH) over noninjective maps F=X+H "
                "with H cubic homogeneous and det JF=1"
            ),
            "nu_cub": (
                "minimum polynomial-matrix nilpotency index nu with (JH)^nu=0 "
                "over the same full cubic-homogeneous class"
            ),
            "rho_HN4": (
                "minimum generic Hessian rank among homogeneous quartic "
                "Hessian-nilpotent Vanishing counterexamples"
            ),
            "n_HN4": (
                "minimum ambient dimension among homogeneous quartic "
                "Hessian-nilpotent Vanishing counterexamples"
            ),
        },
        "rigorous_frontiers": {
            "r_cub": interval(
                3,
                17,
                (
                    "de Bondt--Sun prove every cubic-homogeneous Keller map "
                    "with rank(JH)<=2 is invertible"
                ),
                (
                    "the certified 24-variable circuit-level collision has "
                    "exact generic rank(JH)=17"
                ),
            ),
            "nu_cub": interval(
                3,
                18,
                (
                    "index at most two gives a quasi-translation and hence an "
                    "inverse; no full-class index-three theorem is recorded"
                ),
                (
                    "the certified 22-variable circuit-level collision has "
                    "exact polynomial nilpotency index 18"
                ),
            ),
            "rho_HN4": interval(
                3,
                37,
                (
                    "rank at most two reduces to the invertible "
                    "cubic-homogeneous gradient case"
                ),
                (
                    "the certified 44-variable HN quartic has exact generic "
                    "Hessian rank 37"
                ),
            ),
            "n_HN4": interval(
                6,
                42,
                (
                    "the homogeneous symmetric nilpotent-Jacobian case is "
                    "known through dimension five"
                ),
                "the certified homogeneous HN quartic counterexample uses 42 variables",
            ),
        },
        "certified_rank_upper_witness": {
            "cubic_artifact": "rank_reduced_bcw_24_counterexample.json",
            "cubic_generator": "scripts/verify_rank_reduced_bcw_24_route.py",
            "cubic_independent_replay": (
                "scripts/audit_rank_reduced_bcw_24_independent.py"
            ),
            "dimension": 24,
            "generic_rank_JH": 17,
            "nilpotency_index_JH": 18,
            "specialized_JH_jordan_type": [18, 1, 1, 1, 1, 1, 1],
        },
        "smallest_ambient_cubic_witness_and_HN_source": {
            "cubic_artifact": "low_complexity_bcw_21_counterexample.json",
            "source_artifact": "essential_bcw_21_counterexample.json",
            "dimension": 21,
            "generic_rank_JH": 18,
            "nilpotency_index_JH": 19,
            "specialized_JH_jordan_type": [19, 1, 1],
            "quartic_source_artifact": "essential_bcw_21_counterexample.json",
            "quartic_generator": "scripts/generate_image_vanishing_counterexamples.py",
            "quartic_rank_certificate": "scripts/audit_fixed_rank_hessian_witness.py",
            "quartic_dimension": 42,
            "generic_hessian_rank": 38,
            "sampled_hessian_nilpotency_index": 35,
            "sampled_index_status": (
                "diagnostic lower bound on the polynomial index, not an exact "
                "polynomial-matrix nilpotency certificate"
            ),
        },
        "certified_index_upper_witness": {
            "artifact": "index_reduced_bcw_22_counterexample.json",
            "generator": "scripts/verify_index_reduced_bcw_22_route.py",
            "independent_replay": (
                "scripts/audit_index_reduced_bcw_22_independent.py"
            ),
            "dimension": 22,
            "generic_rank_JH": 18,
            "nilpotency_index_JH": 18,
            "specialized_JH_jordan_type": [18, 2, 1, 1],
            "direct_polynomial_power_certificate": (
                "(JH)^17 has six nonzero matrix entries while (JH)^18=0"
            ),
        },
        "certified_HN_rank_upper_witness": {
            "cubic_artifact": (
                "hessian_rank_reduced_bcw_22_counterexample.json"
            ),
            "generator": (
                "scripts/verify_hessian_rank_reduced_bcw_22_route.py"
            ),
            "independent_replay": (
                "scripts/audit_hessian_rank_reduced_bcw_22_independent.py"
            ),
            "cubic_dimension": 22,
            "quartic_dimension": 44,
            "generic_rank_JH": 18,
            "nilpotency_index_JH": 18,
            "generic_hessian_rank": 37,
            "sampled_transformed_hessian_nilpotency_index": 35,
            "sampled_index_status": (
                "diagnostic specialization only; not an exact "
                "polynomial-matrix nilpotency certificate"
            ),
            "exact_hessian_certificate": (
                "12 polynomial syzygy generators contain 7 generically "
                "independent kernel columns, while an exact rational "
                "specialization has rank 37"
            ),
        },
        "class_scope_warning": {
            "full_cubic_index_three": (
                "invertibility is open in this ledger; the uniform degree-9 "
                "bound is false by van den Essen's dimension-5 example"
            ),
            "power_linear_index_three": (
                "invertible by the cited power-linear results, a strictly "
                "narrower class than arbitrary cubic-homogeneous H"
            ),
            "symmetric_jacobian_index_three": (
                "invertible by Wright's homogeneous symmetric theorem, also "
                "strictly narrower than the arbitrary tensor class"
            ),
        },
        "attack_reductions": {
            "rank_index_coupling": {
                "identity": "nilpotency index(JH) <= rank(JH)+1",
                "consequence": (
                    "invertibility through index d implies "
                    "nu_cub>=d+1 and r_cub>=d"
                ),
                "current_rank_witness_saturates": "18=17+1",
            },
            "index_three_inverse_degree": {
                "lower_calibration_artifact": (
                    "index_three_inverse_degree_model.json"
                ),
                "counterexample_artifact": (
                    "index_three_degree_bound_counterexample.json"
                ),
                "exact_result": (
                    "van den Essen's dimension-5 generic-rank-3 cubic map "
                    "has weak index 3 and inverse degree 13"
                ),
                "proposed_uniform_bound_9": "false",
                "invertibility_only_status": "open for the full class",
            },
            "index_three_degree_11_obstruction": {
                "artifact": "index_three_inverse_tree_obstruction.json",
                "tree_space_dimension": 8,
                "known_relation_rank": 5,
                "quotient_dimension": 3,
                "normal_form": (
                    "-1/2 C(NH,H,H) -1/2 B(B(H,H),H) "
                    "-1/6 N(C(H,H,H))"
                ),
                "status": (
                    "realized nontrivially by van den Essen's exact tensor, "
                    "so it is not killed by the full coefficient ideal"
                ),
            },
            "rank_three_index_three_normal_form_exclusion": {
                "artifact": "index_three_rank_normal_form_exclusion.json",
                "family": (
                    "four-variable lambda plus two arbitrary binary-cubic "
                    "normal form"
                ),
                "result": (
                    "index three forces lambda=0, rank at most 2, and an "
                    "explicit inverse of degree at most 5; rank three forces "
                    "index at least 4"
                ),
            },
            "cotangent_kernel_excess": {
                "artifact": "cotangent_kernel_excess_frontier.json",
                "identity": (
                    "rank Hess(y.H)=2 rank(JH)+rank(K^t A K), "
                    "columns(K)=ker(JH)"
                ),
                "next_strict_targets": [
                    "rank(JH)=18 and kernel excess 0",
                    "rank(JH)=17 and kernel excess at most 2",
                ],
            },
            "cotangent_dimension": (
                "beating quartic dimension 42 by a cotangent lift requires "
                "a cubic source of dimension at most 20"
            ),
        },
        "primary_sources": {
            "rank_two_cubic": "https://arxiv.org/abs/1803.05551",
            "power_linear_index_three": "https://arxiv.org/abs/1302.5864",
            "cubic_linear_index_three": "https://arxiv.org/abs/1508.02012",
            "symmetric_jacobian_index_three": (
                "https://arxiv.org/abs/math/0511214"
            ),
            "index_three_degree_bound_counterexample": (
                "https://eudml.org/doc/262620"
            ),
            "HN_vanishing_equivalence": "https://arxiv.org/abs/math/0409534",
            "symmetric_low_dimension": (
                "https://doi.org/10.1016/j.jpaa.2004.08.030"
            ),
        },
        "search_policy": {
            "script": "scripts/search_restricted_bcw_circuits.py",
            "required_change": (
                "alter the polynomial BCW gate circuit before homogenization"
            ),
            "partial_score": (
                "rank and complete sampled power-rank tuple of J(F-id)"
            ),
            "terminal_score": (
                "post-homogenization/post-kernel cubic rank/index plus "
                "cotangent-Hessian rank and transformed-Hessian power ranks"
            ),
            "promotion_rule": (
                "freeze any modular improvement, replay its collision "
                "independently, and certify generic rank and polynomial powers over QQ"
            ),
            "bounded_negative_records": [
                {
                    "artifact": (
                        "restricted_bcw_circuit_search_xxs_w32_d28.json"
                    ),
                    "scope": "linear-cubic first-component gate split",
                    "terminal_count": 244,
                    "best_objective": [19, 18, 38, 23],
                },
                {
                    "artifact": (
                        "restricted_bcw_circuit_search_yvyb_structural_w22.json"
                    ),
                    "scope": "balanced (y*v)*(y*b) structural beam",
                    "terminal_count": 29,
                    "best_objective": [18, 18, 38, 23],
                },
                {
                    "artifact": (
                        "restricted_bcw_circuit_search_ayb_yvyb_w36.json"
                    ),
                    "scope": "shared y*b gate across two target components",
                    "terminal_count": 46,
                    "best_objective": [18, 18, 40, 25],
                },
                {
                    "artifacts": [
                        (
                            "restricted_bcw_circuit_search_xvvz_v2h_"
                            "structural_w25.json"
                        ),
                        (
                            "restricted_bcw_circuit_search_xvvz_v2h_"
                            "mixed_w25.json"
                        ),
                    ],
                    "scope": (
                        "shared v*z gate across the second and third target "
                        "components, under structural and mixed beams"
                    ),
                    "terminal_counts": [47, 35],
                    "best_objective": [18, 18, 41, 23],
                },
                {
                    "artifact": "rank37_gate_perturbation_search.json",
                    "scope": (
                        "eight nonzero coefficients of the no-new-gate "
                        "lambda*a*q target shear continued from the frozen "
                        "rank-37 prefix"
                    ),
                    "terminal_count": 32,
                    "best_objective": [18, 18, 37, 22, 35],
                    "conclusion": (
                        "all eight coefficient-family leaders reproduce "
                        "the incumbent profile; lambda*q is removed exactly "
                        "by affine normalization"
                    ),
                },
                {
                    "artifact": (
                        "restricted_bcw_circuit_search_aspert_m12.json"
                    ),
                    "scope": (
                        "exceptional -12*a*s neutral shear, whose coefficient "
                        "cancels a step of the frozen cleanup"
                    ),
                    "terminal_count": 33,
                    "best_objective": [18, 18, 41, 23],
                },
                {
                    "artifacts": [
                        "restricted_bcw_circuit_search_xxs_rank_hybrid_w24.json",
                        "restricted_bcw_circuit_search_xxs_v2r_w16.json",
                        "restricted_bcw_circuit_search_xxs_y2vb_w16.json",
                    ],
                    "scope": (
                        "zero-excess xxs gate combined with the two known "
                        "rank-reducing atoms, jointly and pairwise"
                    ),
                    "terminal_counts": [48, 32, 34],
                    "best_objectives": [
                        [22, 19, 45, 28],
                        [21, 19, 44, 27],
                        [19, 18, 43, 26],
                    ],
                    "kernel_excesses": [1, 2, 5],
                    "conclusion": (
                        "the existing zero-excess and rank-reducing atoms "
                        "are not naively composable"
                    ),
                },
                {
                    "artifact": (
                        "restricted_bcw_circuit_search_all_w64.json"
                    ),
                    "scope": (
                        "all nine circuit atoms, width 64, prebeam factor 4, "
                        "and 24 cleanup depths"
                    ),
                    "terminal_count": 140,
                    "pareto_objectives": [
                        [17, 18, 38, 24],
                        [17, 18, 39, 23],
                        [18, 18, 37, 22],
                    ],
                    "conclusion": (
                        "no modular profile beats a certified endpoint; "
                        "the rank-37 qb+x2s route is rediscovered"
                    ),
                },
            ],
        },
    }
    OUTPUT.write_text(json.dumps(frontier, indent=2) + "\n")
    print("PASS restricted minima: rank witness has exact generic rank 17")
    print("PASS restricted minima: circuit witness has exact nilpotency index 18")
    print("PASS restricted minima: HN quartic upper bounds are dimension 42 and rank 37")
    print("PASS restricted minima: full cubic and power-linear index-three scopes are separated")
    print("PASS restricted minima: van den Essen disproves the uniform degree-9 bound")
    print(
        "OPEN exact intervals: "
        "3<=r_cub<=17, 3<=nu_cub<=18, "
        "3<=rho_HN4<=37, 6<=n_HN4<=42"
    )
    print(f"PASS restricted minima: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
