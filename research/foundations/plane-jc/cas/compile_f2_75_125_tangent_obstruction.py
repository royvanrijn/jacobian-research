#!/usr/bin/env python3
"""Compile directional slices of the localized F2 tangent obstruction map.

The source is the GF(31) localized seed and its 169-coordinate spacing-four
staircase chart.  The restricted inhomogeneous tangent fiber has dimension
112.  For every free coordinate, this compiler follows the affine line
``d0+u*k`` inside that fiber, evaluates all seven required nonzero parameter
values, and interpolates the projection of the nonlinear remainder to the
153-dimensional left cokernel of the full Jacobian.  Deterministic dense
kernel directions audit mixed terms that coordinate lines cannot see.

All interpolation is exact over the selected finite field.  Even a complete
directional census is not a characteristic-zero lift or an emptiness proof:
mixed multivariate coefficients are not reconstructed by finitely many lines.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random

import sympy as sp

from compile_f2_75_125_nonlinear_forcing import ROOT
from probe_f2_75_125_nonlinear_modular import (
    ARTIFACT as MODULAR_PROBE_ARTIFACT,
    build_modular_presentation,
    initial_localized_point,
    residual_record,
    roots_mod_quadratic,
    spacing_four_staircase_variables,
)
from sparse_circuit_modp import (
    apply_correction,
    apply_row_functionals,
    evaluate,
    interpolate_consecutive_values,
    interpolate_corrections,
    left_cokernel_basis,
    polynomial_gcd_modp,
    quadratic_key_embedding,
    solve_linearization,
)


ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_tangent_obstruction.json"
)


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def obstruction_line(
    *,
    dag: object,
    roots: list[int],
    a_node: int,
    seed: dict[str, int],
    base_correction: dict[str, int],
    unit_correction: dict[str, int],
    base_projection: list[int],
    functionals: tuple[dict[int, int], ...],
    prime: int,
    embedding: object,
    degree_bound: int,
) -> dict[str, object]:
    projected_values = [base_projection]
    a_values: list[int] = []
    raw_nonzero_counts: list[int] = []
    for parameter in range(degree_bound + 1):
        correction = interpolate_corrections(
            base_correction, unit_correction, parameter, prime
        )
        point = apply_correction(seed, correction, prime)
        if parameter == 0:
            evaluation = evaluate(
                dag, point, prime, embedding, with_jacobian=False
            )
            a_values.append(evaluation.values[a_node])
            raw_nonzero_counts.append(
                sum(bool(evaluation.values[root]) for root in roots[:-1])
            )
            continue
        evaluation = evaluate(
            dag, point, prime, embedding, with_jacobian=False
        )
        geometric = [evaluation.values[root] for root in roots[:-1]]
        projected_values.append(
            apply_row_functionals(functionals, geometric + [0], prime)
        )
        a_values.append(evaluation.values[a_node])
        raw_nonzero_counts.append(sum(bool(value) for value in geometric))

    polynomials = [
        interpolate_consecutive_values(
            [values[coordinate] for values in projected_values],
            prime,
            degree_bound,
        )
        for coordinate in range(len(functionals))
    ]
    nonzero = [
        {
            "coordinate": coordinate,
            "coefficients_low_to_high": polynomial,
        }
        for coordinate, polynomial in enumerate(polynomials)
        if polynomial != [0]
    ]
    common_gcd = [0]
    for polynomial in polynomials:
        common_gcd = polynomial_gcd_modp(common_gcd, polynomial, prime)
        if common_gcd == [1]:
            break
    a_polynomial = interpolate_consecutive_values(a_values, prime, 1)
    return {
        "maximum_obstruction_degree": max(
            (len(item["coefficients_low_to_high"]) - 1 for item in nonzero),
            default=0,
        ),
        "nonzero_obstruction_coordinate_count": len(nonzero),
        "nonzero_obstruction_coordinates": nonzero,
        "common_gcd_coefficients_low_to_high": common_gcd,
        "a_polynomial_coefficients_low_to_high": a_polynomial,
        "excluded_localization_parameters": [
            parameter
            for parameter in range(prime)
            if sum(
                coefficient * parameter**degree
                for degree, coefficient in enumerate(a_polynomial)
            )
            % prime
            == 0
        ],
        "raw_nonzero_counts_at_parameters_0_through_degree_bound": (
            raw_nonzero_counts
        ),
        "restriction_digest_sha256": digest(polynomials),
    }


def build_payload(
    *,
    prime: int = 31,
    rho: int = 14,
    y: int = 3,
    random_line_count: int = 8,
    random_seed: int = 75125112,
    coordinate_limit: int | None = None,
) -> dict[str, object]:
    if rho not in roots_mod_quadratic(1, -3, 1, prime):
        raise ValueError("rho is not a carrier-field root")
    if y not in roots_mod_quadratic(27, -9, 1, prime):
        raise ValueError("y is not a defect-field root")

    dag, roots, groups, a_node, equation_labels = build_modular_presentation()
    embedding = quadratic_key_embedding(prime, rho)
    seed = initial_localized_point(prime, rho, y)
    evaluated = evaluate(
        dag, seed, prime, embedding, with_jacobian=True
    )
    allowed = spacing_four_staircase_variables(evaluated.variable_names)
    tangent = solve_linearization(
        evaluated, roots, prime, allowed_variables=allowed
    )
    if tangent.inconsistent_rows:
        raise AssertionError("the pinned staircase tangent fiber disappeared")
    free_variables = sorted(allowed - set(tangent.pivot_variables))
    if len(free_variables) != 112:
        raise AssertionError("the staircase tangent dimension changed")

    cokernel = left_cokernel_basis(evaluated, roots, prime)
    if cokernel.rank != 214 or len(cokernel.functionals) != 153:
        raise AssertionError("the full Jacobian cokernel changed")

    base_point = apply_correction(seed, tangent.correction, prime)
    base_evaluation = evaluate(
        dag, base_point, prime, embedding, with_jacobian=False
    )
    base_geometric = [base_evaluation.values[root] for root in roots[:-1]]
    base_projection = apply_row_functionals(
        cokernel.functionals, base_geometric + [0], prime
    )
    if any(base_projection):
        raise AssertionError("the pinned particular correction is obstructed")

    degree_bound = max(dag.degrees[root] for root in roots[:-1])
    audited_free_variables = (
        free_variables
        if coordinate_limit is None
        else free_variables[:coordinate_limit]
    )
    coordinate_profiles: list[dict[str, object]] = []
    for index, free_variable in enumerate(audited_free_variables, start=1):
        unit = solve_linearization(
            evaluated,
            roots,
            prime,
            allowed_variables=allowed,
            free_values={free_variable: 1},
        )
        profile = obstruction_line(
            dag=dag,
            roots=roots,
            a_node=a_node,
            seed=seed,
            base_correction=tangent.correction,
            unit_correction=unit.correction,
            base_projection=base_projection,
            functionals=cokernel.functionals,
            prime=prime,
            embedding=embedding,
            degree_bound=degree_bound,
        )
        profile["free_variable"] = free_variable
        coordinate_profiles.append(profile)
        if index % 16 == 0 or index == len(audited_free_variables):
            print(
                f"COORDINATE_LINES {index}/{len(audited_free_variables)}",
                flush=True,
            )

    generator = random.Random(random_seed)
    random_profiles: list[dict[str, object]] = []
    for line in range(random_line_count):
        free_values = {
            name: generator.randrange(prime) for name in free_variables
        }
        if not any(free_values.values()):
            free_values[free_variables[0]] = 1
        unit = solve_linearization(
            evaluated,
            roots,
            prime,
            allowed_variables=allowed,
            free_values=free_values,
        )
        profile = obstruction_line(
            dag=dag,
            roots=roots,
            a_node=a_node,
            seed=seed,
            base_correction=tangent.correction,
            unit_correction=unit.correction,
            base_projection=base_projection,
            functionals=cokernel.functionals,
            prime=prime,
            embedding=embedding,
            degree_bound=degree_bound,
        )
        profile["line"] = line
        profile["free_assignment_digest_sha256"] = digest(free_values)
        profile["free_assignment_support_size"] = sum(
            bool(value) for value in free_values.values()
        )
        random_profiles.append(profile)
        print(f"MIXED_LINES {line + 1}/{random_line_count}", flush=True)

    coordinate_union = sorted(
        {
            item["coordinate"]
            for profile in coordinate_profiles
            for item in profile["nonzero_obstruction_coordinates"]
        }
    )
    random_union = sorted(
        {
            item["coordinate"]
            for profile in random_profiles
            for item in profile["nonzero_obstruction_coordinates"]
        }
    )
    total_union = sorted(set(coordinate_union) | set(random_union))
    functional_ledger = [
        {
            "coordinate": coordinate,
            "equation_coefficients": [
                {
                    "equation": equation,
                    "label": equation_labels[equation],
                    "coefficient": coefficient,
                }
                for equation, coefficient in sorted(
                    cokernel.functionals[coordinate].items()
                )
            ],
        }
        for coordinate in total_union
    ]
    degree_histogram = Counter(
        int(profile["maximum_obstruction_degree"])
        for profile in coordinate_profiles
    )
    random_degree_histogram = Counter(
        int(profile["maximum_obstruction_degree"])
        for profile in random_profiles
    )
    return {
        "schema": "plane-jc.f2-75-125-tangent-obstruction.v1",
        "status": "exact-finite-field-directional-obstruction-audit",
        "source_modular_probe_artifact": str(
            MODULAR_PROBE_ARTIFACT.relative_to(ROOT)
        ),
        "field": {"prime": prime, "rho": rho, "y": y},
        "presentation": {
            "equation_count_with_localization": len(roots),
            "geometric_equation_count": len(roots) - 1,
            "maximum_degree": degree_bound,
            "initial_residual": residual_record(
                evaluated.values, roots, groups
            ),
        },
        "tangent_fiber": {
            "staircase_variable_count": len(allowed),
            "restricted_Jacobian_rank": tangent.rank,
            "dimension": len(free_variables),
            "free_variables": free_variables,
            "particular_correction_support_size": len(tangent.correction),
            "particular_projected_obstruction_is_zero": True,
        },
        "full_Jacobian_cokernel": {
            "Jacobian_rank": cokernel.rank,
            "dimension": len(cokernel.functionals),
            "basis_digest_sha256": digest(
                [sorted(functional.items()) for functional in cokernel.functionals]
            ),
        },
        "coordinate_line_audit": {
            "line_count": len(coordinate_profiles),
            "complete_free_coordinate_census": (
                len(coordinate_profiles) == len(free_variables)
            ),
            "degree_histogram": {
                str(degree): count
                for degree, count in sorted(degree_histogram.items())
            },
            "nonzero_coordinate_union": coordinate_union,
            "nonzero_coordinate_union_count": len(coordinate_union),
            "profiles": coordinate_profiles,
        },
        "mixed_line_audit": {
            "line_count": len(random_profiles),
            "random_seed": random_seed,
            "degree_histogram": {
                str(degree): count
                for degree, count in sorted(random_degree_histogram.items())
            },
            "nonzero_coordinate_union": random_union,
            "nonzero_coordinate_union_count": len(random_union),
            "profiles": random_profiles,
        },
        "combined_obstruction_coordinate_union": total_union,
        "combined_obstruction_coordinate_union_count": len(total_union),
        "combined_obstruction_functional_ledger": functional_ledger,
        "claim_boundary": (
            "coordinate lines plus finitely many deterministic dense lines "
            "do not reconstruct or prove the complete 112-parameter "
            "multivariate obstruction map"
        ),
        "software": {"sympy": sp.__version__},
    }


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--random-lines", type=int, default=8)
    parser.add_argument("--coordinate-limit", type=int)
    args = parser.parse_args()

    payload = build_payload(
        random_line_count=args.random_lines,
        coordinate_limit=args.coordinate_limit,
    )
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            artifact_label = artifact.relative_to(ROOT)
        except ValueError:
            artifact_label = artifact
        print(f"WROTE {artifact_label}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {
            key: value for key, value in payload.items() if key != "software"
        }
        pinned_claim = {
            key: value for key, value in expected.items() if key != "software"
        }
        if current_claim != pinned_claim:
            raise AssertionError(
                "the pinned tangent-obstruction artifact is stale"
            )
    print("F2_TANGENT_OBSTRUCTION_DIRECTIONAL_PASS")
    print(f"F2_TANGENT_OBSTRUCTION_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()
