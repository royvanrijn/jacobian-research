#!/usr/bin/env python3
"""Probe the localized F2 nonlinear circuit ideal over a good split prime.

The deterministic default reduction is GF(31), with ``rho=14`` and ``y=3``.
It starts from a genuine descent-eight local seed: ``a=1`` for the cubic
coefficient of K_P7 at rho, ``K_Q1(rho)=0``, and an explicit Rabinowitsch
inverse imposing ``a != 0``.  Sparse full-Jacobian Newton corrections are
then computed for all displayed equations.  Because Newton iteration inside a
single finite field is not a Hensel contraction, each tangent direction is
searched exactly along its full affine GF(p)-line; only a strict reduction of
the residual support is accepted.

A point found here is an exact point of one finite-field reduction.  It does
not lift automatically to characteristic zero and is not a Keller map.
Conversely, failure of this heuristic Newton chart would not prove that the
localized ideal is empty.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path
import re

import sympy as sp

from compile_f2_75_125_nonlinear_forcing import (
    ARTIFACT as NONLINEAR_ARTIFACT,
    ROOT,
    CircuitDAG,
    RHO_FIELD,
    build_band_circuits,
    coupled_forcing_audit,
    descent_eight_incidence_audit,
    evaluate_polynomial,
    target_and_layer_zero_audit,
)
from sparse_circuit_modp import (
    apply_row_functionals,
    apply_scaled_correction,
    evaluate,
    interpolate_consecutive_values,
    left_cokernel_basis,
    polynomial_gcd_modp,
    quadratic_key_embedding,
    solve_linearization,
)


ARTIFACT = ROOT / "artifacts/generated-results/jc2_f2_75_125_modular_probe.json"

P_STAIRCASE_LAYERS = frozenset(range(-17, 12, 4))
Q_STAIRCASE_LAYERS = frozenset(range(-3, 14, 4))
STAIRCASE_PROBE_FREE_VARIABLES = (
    "P_-17_d0",
    "Q_13_d2",
    "Q_9_d0",
    "Q_5_d2",
    "Q_1_d2",
    "Q_-3_d0",
)


def roots_mod_quadratic(a: int, b: int, c: int, prime: int) -> list[int]:
    return [
        value
        for value in range(prime)
        if (a * value * value + b * value + c) % prime == 0
    ]


def initial_localized_point(prime: int, rho: int, y: int) -> dict[str, int]:
    """Set K_P7=(w-rho)^3, K_Q1(rho)=0, y, and a_inverse=1."""

    offset = (1 - rho) % prime
    return {
        "P_7_d0": offset**3 % prime,
        "P_7_d1": 3 * offset**2 % prime,
        "P_7_d2": 3 * offset % prime,
        "P_7_d3": 1,
        "Q_1_d1": pow(rho - 1, -1, prime),
        "descent8_y": y,
        "localize_a_inverse": 1,
    }


def spacing_four_staircase_variables(variable_names: list[str]) -> set[str]:
    """Return the natural P/Q step-four chart seen by the first tangent."""

    selected = {"descent8_y", "localize_a_inverse"}
    for name in variable_names:
        match = re.fullmatch(r"([PQ])_(-?\d+)_d\d+", name)
        if match is None:
            continue
        side, layer_text = match.groups()
        layer = int(layer_text)
        if (
            side == "P" and layer in P_STAIRCASE_LAYERS
        ) or (
            side == "Q" and layer in Q_STAIRCASE_LAYERS
        ):
            selected.add(name)
    return selected


def build_modular_presentation() -> tuple[
    CircuitDAG,
    list[int],
    dict[str, tuple[int, int]],
    int,
    list[str],
]:
    dag = CircuitDAG()
    bands = build_band_circuits(dag)
    coupled_audit = coupled_forcing_audit(dag, bands)
    final_audit = target_and_layer_zero_audit(dag, bands)
    incidence_audit = descent_eight_incidence_audit(dag, bands)
    coupled = coupled_audit["all_equation_nodes"]
    final = final_audit["all_equation_nodes"]
    incidence = incidence_audit["all_equation_nodes"]

    p_seven = bands["K_polynomials"]["P", 7]
    a = dag.scale(
        RHO_FIELD.convert(sp.Rational(1, 6)),
        evaluate_polynomial(dag, p_seven, RHO_FIELD.unit, 3),
    )
    inverse = dag.variable("localize_a_inverse")
    localization = dag.add(
        dag.multiply(a, inverse),
        dag.scale(RHO_FIELD.convert(-1), dag.one),
    )
    roots = coupled + final + incidence + [localization]
    groups = {
        "coupled_Laurent": (0, len(coupled)),
        "target_and_Hermite": (len(coupled), len(coupled) + len(final)),
        "incidence_and_defect": (
            len(coupled) + len(final),
            len(coupled) + len(final) + len(incidence),
        ),
        "a_localization": (len(roots) - 1, len(roots)),
    }
    equation_labels = [
        f"coupled.descent_{row['descent']}.divisibility_{coordinate}"
        for row in coupled_audit["rows"]
        for coordinate in range(row["divisibility_coordinate_count"])
    ] + [
        f"coupled.descent_{row['descent']}.quotient_cokernel_{coordinate}"
        for row in coupled_audit["rows"]
        for coordinate in range(row["quotient_cokernel_coordinate_count"])
    ]
    equation_labels += [
        *[f"target.rho_jet_{order}" for order in range(5)],
        "target.quotient_residue_0",
        "target.quotient_residue_1",
        "layer_zero.H_rho_minus_H_zero",
        *[f"layer_zero.H_rho_derivative_{order}" for order in range(1, 6)],
        "incidence.K_P7_rho",
        "incidence.K_P7_prime_rho",
        "incidence.K_P7_second_rho",
        "incidence.K_Q1_rho",
        "incidence.K_Pminus1_ratio",
        "incidence.relative_quadratic_defect",
        "localization.a_inverse_minus_one",
    ]
    if len(equation_labels) != len(roots):
        raise AssertionError("modular equation label ledger changed")
    return dag, roots, groups, a, equation_labels


def residual_record(
    values: list[int],
    roots: list[int],
    groups: dict[str, tuple[int, int]],
    *,
    equation_overrides: dict[int, int] | None = None,
) -> dict[str, object]:
    overrides = equation_overrides or {}
    residuals = [
        overrides.get(index, values[root])
        for index, root in enumerate(roots)
    ]
    return {
        "nonzero_total": sum(bool(value) for value in residuals),
        "nonzero_by_group": {
            name: sum(bool(value) for value in residuals[start:end])
            for name, (start, end) in groups.items()
        },
        "nonzero_indices": [
            index for index, value in enumerate(residuals) if value
        ],
        "residual_digest_sha256": hashlib.sha256(
            json.dumps(residuals, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def build_payload(
    *, prime: int = 31, rho: int = 14, y: int = 3, maximum_steps: int = 8
) -> dict[str, object]:
    rho_roots = roots_mod_quadratic(1, -3, 1, prime)
    y_roots = roots_mod_quadratic(27, -9, 1, prime)
    if rho not in rho_roots or y not in y_roots:
        raise ValueError("rho and y must define split quadratic embeddings")
    if prime in (3, 5) or rho in (0, 1):
        raise ValueError("the selected prime violates a carrier localization")

    dag, roots, groups, a_node, equation_labels = build_modular_presentation()
    embedding = quadratic_key_embedding(prime, rho)
    point = initial_localized_point(prime, rho, y)
    steps: list[dict[str, object]] = []
    converged = False

    for step in range(maximum_steps + 1):
        evaluated = evaluate(
            dag, point, prime, embedding, with_jacobian=True
        )
        residual = residual_record(evaluated.values, roots, groups)
        record: dict[str, object] = {
            "step": step,
            **residual,
            "a_value": evaluated.values[a_node],
            "point_support_size": len(point),
            "Jacobian_gradient_entry_count": evaluated.gradient_entry_count,
        }
        if residual["nonzero_total"] == 0:
            converged = True
            record["Jacobian_rank"] = solve_linearization(
                evaluated, roots, prime
            ).rank
            steps.append(record)
            break
        linearized = solve_linearization(evaluated, roots, prime)
        record.update(
            {
                "Jacobian_rank": linearized.rank,
                "linearized_inconsistent_rows": (
                    linearized.inconsistent_rows
                ),
                "correction_support_size": len(linearized.correction),
                "correction_support_names": sorted(linearized.correction),
            }
        )
        if step == 0:
            staircase_variables = spacing_four_staircase_variables(
                evaluated.variable_names
            )
            staircase_linearized = solve_linearization(
                evaluated,
                roots,
                prime,
                allowed_variables=staircase_variables,
            )
            record["spacing_four_staircase_tangent"] = {
                "P_layers": sorted(P_STAIRCASE_LAYERS, reverse=True),
                "Q_layers": sorted(Q_STAIRCASE_LAYERS, reverse=True),
                "variable_count": len(staircase_variables),
                "Jacobian_rank": staircase_linearized.rank,
                "affine_tangent_dimension": (
                    len(staircase_variables) - staircase_linearized.rank
                ),
                "linearized_inconsistent_rows": (
                    staircase_linearized.inconsistent_rows
                ),
                "particular_correction_support_size": len(
                    staircase_linearized.correction
                ),
                "particular_correction_support_names": sorted(
                    staircase_linearized.correction
                ),
            }
            full_cokernel = left_cokernel_basis(
                evaluated, roots, prime
            )
            if full_cokernel.rank != linearized.rank:
                raise AssertionError("Jacobian rank/cokernel rank mismatch")
            initial_projection = apply_row_functionals(
                full_cokernel.functionals,
                [evaluated.values[root] for root in roots],
                prime,
            )
            if any(initial_projection):
                raise AssertionError(
                    "consistent tangent forcing has nonzero cokernel image"
                )
            serialized_functionals = [
                sorted(functional.items())
                for functional in full_cokernel.functionals
            ]
            record["full_Jacobian_cokernel"] = {
                "dimension": len(full_cokernel.functionals),
                "Jacobian_rank": full_cokernel.rank,
                "initial_forcing_projection_is_zero": True,
                "basis_digest_sha256": hashlib.sha256(
                    json.dumps(
                        serialized_functionals, separators=(",", ":")
                    ).encode()
                ).hexdigest(),
            }
        steps.append(record)
        if linearized.inconsistent_rows or not linearized.correction:
            break

        current_score = int(residual["nonzero_total"])
        best_score = current_score
        best_scale: int | None = None
        best_point: dict[str, int] | None = None
        candidate_scores: list[int | None] = []
        geometric_line_values = [
            [evaluated.values[root] for root in roots[:-1]]
        ]
        for scale in range(1, prime):
            trial = apply_scaled_correction(
                point, linearized.correction, scale, prime
            )
            trial_evaluation = evaluate(
                dag, trial, prime, embedding, with_jacobian=False
            )
            trial_a = trial_evaluation.values[a_node]
            geometric_line_values.append(
                [trial_evaluation.values[root] for root in roots[:-1]]
            )
            if not trial_a:
                candidate_scores.append(None)
                continue

            # The inverse is an auxiliary localization coordinate.  Eliminate
            # it exactly after moving the geometric variables, rather than
            # mistaking the quadratic Rabinowitsch residual for geometry.
            trial["localize_a_inverse"] = pow(trial_a, -1, prime)
            trial_residual = residual_record(
                trial_evaluation.values,
                roots,
                groups,
                equation_overrides={len(roots) - 1: 0},
            )
            trial_score = int(trial_residual["nonzero_total"])
            candidate_scores.append(trial_score)
            if trial_score < best_score:
                best_score = trial_score
                best_scale = scale
                best_point = trial

        degree_bound = max(dag.degrees[root] for root in roots[:-1])
        line_polynomials = [
            interpolate_consecutive_values(
                [values[equation] for values in geometric_line_values],
                prime,
                degree_bound,
            )
            for equation in range(len(roots) - 1)
        ]
        line_gcd = [0]
        for polynomial in line_polynomials:
            line_gcd = polynomial_gcd_modp(line_gcd, polynomial, prime)
            if line_gcd == [1]:
                break
        if step == 0:
            projected_line_values = [
                apply_row_functionals(
                    full_cokernel.functionals,
                    values + [0],
                    prime,
                )
                for values in geometric_line_values
            ]
            projected_line_polynomials = [
                interpolate_consecutive_values(
                    [
                        values[coordinate]
                        for values in projected_line_values
                    ],
                    prime,
                    degree_bound,
                )
                for coordinate in range(len(full_cokernel.functionals))
            ]
            projected_line_gcd = [0]
            for polynomial in projected_line_polynomials:
                projected_line_gcd = polynomial_gcd_modp(
                    projected_line_gcd, polynomial, prime
                )
                if projected_line_gcd == [1]:
                    break
        record["tangent_line_search"] = {
            "candidate_nonzero_totals_for_scales_1_through_p_minus_1": (
                candidate_scores
            ),
            "accepted_scale": best_scale,
            "accepted_nonzero_total": (
                best_score if best_scale is not None else None
            ),
            "strict_improvement": best_scale is not None,
            "geometric_equation_restriction": {
                "degree_bound": degree_bound,
                "maximum_observed_degree": max(
                    len(polynomial) - 1 for polynomial in line_polynomials
                ),
                "nonzero_polynomial_count": sum(
                    polynomial != [0] for polynomial in line_polynomials
                ),
                "common_gcd_coefficients_low_to_high": line_gcd,
                "common_gcd_degree": len(line_gcd) - 1,
                "disjoint_from_zero_locus_over_algebraic_closure": (
                    line_gcd == [1]
                ),
                "verified_sample_count": prime,
                "restriction_digest_sha256": hashlib.sha256(
                    json.dumps(
                        line_polynomials, separators=(",", ":")
                    ).encode()
                ).hexdigest(),
            },
        }
        if step == 0:
            record["tangent_line_search"]["Jacobian_cokernel_restriction"] = {
                "coordinate_count": len(full_cokernel.functionals),
                "maximum_observed_degree": max(
                    len(polynomial) - 1
                    for polynomial in projected_line_polynomials
                ),
                "nonzero_polynomial_count": sum(
                    polynomial != [0]
                    for polynomial in projected_line_polynomials
                ),
                "nonzero_at_scale_one": sum(
                    bool(value) for value in projected_line_values[1]
                ),
                "common_gcd_coefficients_low_to_high": (
                    projected_line_gcd
                ),
                "nonzero_coordinate_polynomials": [
                    {
                        "coordinate": coordinate,
                        "coefficients_low_to_high": polynomial,
                    }
                    for coordinate, polynomial in enumerate(
                        projected_line_polynomials
                    )
                    if polynomial != [0]
                ],
                "restriction_digest_sha256": hashlib.sha256(
                    json.dumps(
                        projected_line_polynomials,
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest(),
            }
            coordinate_profiles = []
            for free_name in STAIRCASE_PROBE_FREE_VARIABLES:
                fiber_point = solve_linearization(
                    evaluated,
                    roots,
                    prime,
                    allowed_variables=staircase_variables,
                    free_values={free_name: 1},
                )
                sampled_values = [
                    [evaluated.values[root] for root in roots[:-1]]
                ]
                sampled_a_values = [evaluated.values[a_node]]
                for scale in range(1, degree_bound + 1):
                    trial = apply_scaled_correction(
                        point, fiber_point.correction, scale, prime
                    )
                    trial_evaluation = evaluate(
                        dag, trial, prime, embedding, with_jacobian=False
                    )
                    sampled_values.append(
                        [
                            trial_evaluation.values[root]
                            for root in roots[:-1]
                        ]
                    )
                    sampled_a_values.append(
                        trial_evaluation.values[a_node]
                    )
                fiber_polynomials = [
                    interpolate_consecutive_values(
                        [
                            values[equation]
                            for values in sampled_values
                        ],
                        prime,
                        degree_bound,
                    )
                    for equation in range(len(roots) - 1)
                ]
                fiber_gcd = [0]
                for polynomial in fiber_polynomials:
                    fiber_gcd = polynomial_gcd_modp(
                        fiber_gcd, polynomial, prime
                    )
                    if fiber_gcd == [1]:
                        break
                projected_fiber_values = [
                    apply_row_functionals(
                        full_cokernel.functionals,
                        values + [0],
                        prime,
                    )
                    for values in sampled_values
                ]
                projected_fiber_polynomials = [
                    interpolate_consecutive_values(
                        [
                            values[coordinate]
                            for values in projected_fiber_values
                        ],
                        prime,
                        degree_bound,
                    )
                    for coordinate in range(
                        len(full_cokernel.functionals)
                    )
                ]
                projected_fiber_gcd = [0]
                for polynomial in projected_fiber_polynomials:
                    projected_fiber_gcd = polynomial_gcd_modp(
                        projected_fiber_gcd, polynomial, prime
                    )
                    if projected_fiber_gcd == [1]:
                        break
                a_polynomial = interpolate_consecutive_values(
                    sampled_a_values,
                    prime,
                    1,
                )
                coordinate_profiles.append(
                    {
                        "free_variable": free_name,
                        "free_value": 1,
                        "correction_support_size": len(
                            fiber_point.correction
                        ),
                        "nonzero_at_scale_one": sum(
                            bool(value) for value in sampled_values[1]
                        ),
                        "maximum_observed_degree": max(
                            len(polynomial) - 1
                            for polynomial in fiber_polynomials
                        ),
                        "nonzero_polynomial_count": sum(
                            polynomial != [0]
                            for polynomial in fiber_polynomials
                        ),
                        "common_gcd_coefficients_low_to_high": fiber_gcd,
                        "disjoint_from_zero_locus_over_algebraic_closure": (
                            fiber_gcd == [1]
                        ),
                        "a_polynomial_coefficients_low_to_high": (
                            a_polynomial
                        ),
                        "Jacobian_cokernel_restriction": {
                            "maximum_observed_degree": max(
                                len(polynomial) - 1
                                for polynomial in projected_fiber_polynomials
                            ),
                            "nonzero_polynomial_count": sum(
                                polynomial != [0]
                                for polynomial in projected_fiber_polynomials
                            ),
                            "nonzero_at_scale_one": sum(
                                bool(value)
                                for value in projected_fiber_values[1]
                            ),
                            "common_gcd_coefficients_low_to_high": (
                                projected_fiber_gcd
                            ),
                            "nonzero_coordinate_polynomials": [
                                {
                                    "coordinate": coordinate,
                                    "coefficients_low_to_high": polynomial,
                                }
                                for coordinate, polynomial in enumerate(
                                    projected_fiber_polynomials
                                )
                                if polynomial != [0]
                            ],
                            "restriction_digest_sha256": hashlib.sha256(
                                json.dumps(
                                    projected_fiber_polynomials,
                                    separators=(",", ":"),
                                ).encode()
                            ).hexdigest(),
                        },
                        "restriction_digest_sha256": hashlib.sha256(
                            json.dumps(
                                fiber_polynomials,
                                separators=(",", ":"),
                            ).encode()
                        ).hexdigest(),
                    }
                )
            record["spacing_four_staircase_tangent"][
                "coordinate_line_profiles"
            ] = coordinate_profiles
            obstruction_support = sorted(
                {
                    item["coordinate"]
                    for profile in coordinate_profiles
                    for item in profile["Jacobian_cokernel_restriction"][
                        "nonzero_coordinate_polynomials"
                    ]
                }
            )
            record["spacing_four_staircase_tangent"][
                "sampled_first_obstruction"
            ] = {
                "nonzero_cokernel_coordinate_union": obstruction_support,
                "nonzero_cokernel_coordinate_union_count": len(
                    obstruction_support
                ),
                "detected_free_variables": [
                    profile["free_variable"]
                    for profile in coordinate_profiles
                    if profile["Jacobian_cokernel_restriction"][
                        "nonzero_polynomial_count"
                    ]
                ],
                "detected_cokernel_functionals": [
                    {
                        "coordinate": coordinate,
                        "equation_coefficients": sorted(
                            (
                                {
                                    "equation": equation,
                                    "label": equation_labels[equation],
                                    "coefficient": coefficient,
                                }
                                for equation, coefficient in full_cokernel.functionals[
                                    coordinate
                                ].items()
                            ),
                            key=lambda item: item["equation"],
                        ),
                        "support_size": len(
                            full_cokernel.functionals[coordinate]
                        ),
                        "support_by_equation_group": {
                            name: sum(
                                start <= equation < end
                                for equation in full_cokernel.functionals[
                                    coordinate
                                ]
                            )
                            for name, (start, end) in groups.items()
                        },
                    }
                    for coordinate in obstruction_support
                ],
                "claim_boundary": (
                    "six coordinate lines do not determine the full "
                    "112-parameter obstruction map"
                ),
            }
        if best_point is None:
            break
        point = best_point
        del evaluated
        gc.collect()

    final_evaluation = evaluate(
        dag, point, prime, embedding, with_jacobian=False
    )
    final_residual = residual_record(final_evaluation.values, roots, groups)
    if converged and final_residual["nonzero_total"] != 0:
        raise AssertionError("the modular point failed its independent replay")
    if final_evaluation.values[a_node] == 0:
        raise AssertionError("the modular probe left the a-localized chart")

    sorted_point = dict(sorted(point.items()))
    return {
        "schema": "plane-jc.f2-75-125-nonlinear-modular-probe.v1",
        "status": (
            "exact-finite-field-point" if converged else "exact-finite-field-probe"
        ),
        "source_presentation_artifact": str(NONLINEAR_ARTIFACT.relative_to(ROOT)),
        "field": {
            "prime": prime,
            "rho": rho,
            "rho_roots": rho_roots,
            "y": y,
            "y_roots": y_roots,
        },
        "localization": {
            "equation": "a*localize_a_inverse-1",
            "a": "K_P7'''(rho)/6",
            "final_a_value": final_evaluation.values[a_node],
        },
        "equation_count_with_localization": len(roots),
        "Newton_steps": steps,
        "iteration_method": (
            "full sparse tangent solve followed by exact GF(p) affine-line "
            "search, with exact a-inverse renormalization"
        ),
        "converged": converged,
        "final_residual": final_residual,
        "final_point_support": sorted_point,
        "final_point_digest_sha256": hashlib.sha256(
            json.dumps(sorted_point, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "claim_boundary": (
            "an exact point over one good split prime is not a characteristic-"
            "zero lift, an F2 realization, or a polynomial Keller map"
        ),
        "software": {"sympy": sp.__version__},
    }


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--maximum-steps", type=int, default=8)
    args = parser.parse_args()

    payload = build_payload(maximum_steps=args.maximum_steps)
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"WROTE {artifact.relative_to(ROOT)}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "the pinned modular-probe artifact is stale; inspect before refresh"
            )

    print("F2_MODULAR_LOCALIZED_SEED_PASS")
    print("F2_MODULAR_SPARSE_JACOBIAN_PASS")
    print(
        "F2_MODULAR_NEWTON_"
        + ("POINT_FOUND" if payload["converged"] else "PROBE_INCOMPLETE")
    )
    print(f"F2_MODULAR_PROBE_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()
