#!/usr/bin/env python3
"""Lift the localized F2 seed along a fixed-Jacobian formal homotopy.

Over the pinned good reduction, let F be the 367-equation localized circuit
and x0 the descent-eight seed.  The homotopy equation is

    F(x(lambda)) = (1-lambda) F(x0).

Its first coefficient is the consistent inhomogeneous tangent equation.  We
choose the zero-free-parameter point of the 112-dimensional staircase tangent
fiber, which kills the sampled first obstruction, then solve every higher
coefficient against the fixed Jacobian at x0.  A nonzero left-cokernel
projection is an exact obstruction to continuing that chosen jet.

This constructs a finite formal jet over GF(31), not a value at lambda=1, a
finite-field solution, a characteristic-zero lift, or a Keller map.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import sympy as sp

from compile_f2_75_125_nonlinear_forcing import (
    ARTIFACT as NONLINEAR_ARTIFACT,
    ROOT,
)
from probe_f2_75_125_nonlinear_modular import (
    build_modular_presentation,
    initial_localized_point,
    spacing_four_staircase_variables,
)
from sparse_circuit_modp import (
    apply_row_functionals,
    evaluate,
    evaluate_truncated_series,
    left_cokernel_basis,
    quadratic_key_embedding,
    solve_linearization,
)


ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_formal_homotopy.json"
)
REGULAR_GAUGE_VARIABLES = (
    "P_-5_d0",
    *tuple(f"P_7_d{degree}" for degree in range(2, 8)),
)


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def jacobian_action(
    evaluated: object,
    roots: list[int],
    correction: dict[str, int],
    prime: int,
) -> list[int]:
    name_by_index = evaluated.variable_names
    return [
        sum(
            coefficient * correction.get(name_by_index[index], 0)
            for index, coefficient in evaluated.gradients[root].items()
        )
        % prime
        for root in roots
    ]


def set_series_coefficient(
    variable_series: dict[str, list[int]],
    correction: dict[str, int],
    order: int,
    maximum_order: int,
) -> None:
    for name, value in correction.items():
        coefficients = variable_series.setdefault(
            name, [0] * (maximum_order + 1)
        )
        coefficients[order] = value


def build_payload(
    *,
    prime: int = 31,
    rho: int = 14,
    y: int = 3,
    maximum_order: int = 16,
    regular_gauge: bool = False,
) -> dict[str, object]:
    if maximum_order < 2:
        raise ValueError("the formal lift must include order two")
    dag, roots, groups, a_node, equation_labels = build_modular_presentation()
    bad_constant_keys = sorted(
        {
            node[1]
            for node in dag.nodes
            if node[0] in ("const", "scale")
            and (node[1][0][1] % prime == 0 or node[1][1][1] % prime == 0)
        }
    )
    if bad_constant_keys:
        raise ValueError(
            f"prime {prime} is not a good denominator reduction "
            f"({len(bad_constant_keys)} nonintegral constant keys)"
        )
    embedding = quadratic_key_embedding(prime, rho)
    seed = initial_localized_point(prime, rho, y)
    evaluated = evaluate(dag, seed, prime, embedding, with_jacobian=True)
    if evaluated.gradients is None:
        raise AssertionError("the seed Jacobian was not constructed")
    initial_values = [evaluated.values[root] for root in roots]
    allowed = spacing_four_staircase_variables(evaluated.variable_names)
    first = solve_linearization(
        evaluated, roots, prime, allowed_variables=allowed
    )
    if first.inconsistent_rows:
        raise AssertionError("the first homotopy equation is inconsistent")
    free_variables = sorted(allowed - set(first.pivot_variables))
    if len(free_variables) != 112:
        raise AssertionError("the staircase tangent dimension changed")
    if any(name in first.correction for name in ("Q_9_d0",)) or any(
        f"Q_5_d{degree}" in first.correction for degree in range(2, 19)
    ):
        raise AssertionError("the chosen first jet activated the obstruction variables")

    first_action = jacobian_action(evaluated, roots, first.correction, prime)
    if any(
        (action + value) % prime
        for action, value in zip(first_action, initial_values)
    ):
        raise AssertionError("the first homotopy coefficient failed")
    cokernel = left_cokernel_basis(evaluated, roots, prime)
    if cokernel.rank != 214:
        raise AssertionError("the full Jacobian rank changed")

    variable_series: dict[str, list[int]] = {
        name: [value] + [0] * maximum_order
        for name, value in seed.items()
    }
    set_series_coefficient(
        variable_series, first.correction, 1, maximum_order
    )
    steps: list[dict[str, object]] = [
        {
            "order": 1,
            "forcing_nonzero_count": sum(bool(value) for value in initial_values),
            "cokernel_projection_nonzero_count": 0,
            "correction_support_size": len(first.correction),
            "correction": dict(sorted(first.correction.items())),
            "equation_coefficient_after_correction_digest_sha256": digest(
                [0] * len(roots)
            ),
        }
    ]
    achieved_order = 1
    obstruction: dict[str, object] | None = None

    for order in range(2, maximum_order + 1):
        series = evaluate_truncated_series(
            dag, variable_series, prime, embedding, order
        )
        forcing = [series[root][order] for root in roots]
        projection = apply_row_functionals(
            cokernel.functionals, forcing, prime
        )
        nonzero_projection = [
            {
                "coordinate": coordinate,
                "value": value,
            }
            for coordinate, value in enumerate(projection)
            if value
        ]
        record: dict[str, object] = {
            "order": order,
            "forcing_nonzero_count": sum(bool(value) for value in forcing),
            "forcing_digest_sha256": digest(forcing),
            "cokernel_projection_nonzero_count": len(nonzero_projection),
            "cokernel_projection": nonzero_projection,
            "cokernel_projection_digest_sha256": digest(projection),
        }
        del series
        gc.collect()
        if nonzero_projection:
            obstruction = {
                "order": order,
                "coordinates": nonzero_projection,
                "functional_equation_ledgers": [
                    {
                        "coordinate": item["coordinate"],
                        "equation_coefficients": [
                            {
                                "equation": equation,
                                "label": equation_labels[equation],
                                "coefficient": coefficient,
                            }
                            for equation, coefficient in sorted(
                                cokernel.functionals[item["coordinate"]].items()
                            )
                        ],
                    }
                    for item in nonzero_projection
                ],
            }
            steps.append(record)
            break

        correction = solve_linearization(
            evaluated,
            roots,
            prime,
            right_hand_side=[-value % prime for value in forcing],
            prescribed_values=(
                {name: 0 for name in REGULAR_GAUGE_VARIABLES}
                if regular_gauge
                else None
            ),
        )
        if correction.inconsistent_rows:
            obstruction = {
                "order": order,
                "type": "gauge-incompatible fixed-Jacobian solve",
                "gauge_variables": list(REGULAR_GAUGE_VARIABLES),
                "inconsistent_rows": correction.inconsistent_rows,
            }
            record["gauge_inconsistent_rows"] = correction.inconsistent_rows
            steps.append(record)
            break
        action = jacobian_action(
            evaluated, roots, correction.correction, prime
        )
        corrected_coefficient = [
            (value + delta) % prime
            for value, delta in zip(forcing, action)
        ]
        if any(corrected_coefficient):
            raise AssertionError("a formal coefficient failed after correction")
        set_series_coefficient(
            variable_series,
            correction.correction,
            order,
            maximum_order,
        )
        record.update(
            {
                "Jacobian_rank": correction.rank,
                "correction_support_size": len(correction.correction),
                "correction": dict(sorted(correction.correction.items())),
                "equation_coefficient_after_correction_digest_sha256": digest(
                    corrected_coefficient
                ),
            }
        )
        steps.append(record)
        achieved_order = order
        print(f"FORMAL_HOMOTOPY_ORDER {order}/{maximum_order}", flush=True)

    final_series = evaluate_truncated_series(
        dag, variable_series, prime, embedding, achieved_order
    )
    root_series = [list(final_series[root]) for root in roots]
    expected_failure = False
    for equation, coefficients in enumerate(root_series):
        if coefficients[0] != initial_values[equation]:
            expected_failure = True
        if achieved_order >= 1 and coefficients[1] != -initial_values[equation] % prime:
            expected_failure = True
        if any(coefficients[2:]):
            expected_failure = True
    if expected_failure:
        raise AssertionError("the independent formal-homotopy replay failed")
    a_series = list(final_series[a_node])
    del final_series

    compact_variable_series = {
        name: coefficients[: achieved_order + 1]
        for name, coefficients in sorted(variable_series.items())
        if any(coefficients[: achieved_order + 1])
    }
    truncated_lambda_one_point = {
        name: sum(coefficients) % prime
        for name, coefficients in compact_variable_series.items()
        if sum(coefficients) % prime
    }
    lambda_one_evaluation = evaluate(
        dag,
        truncated_lambda_one_point,
        prime,
        embedding,
        with_jacobian=False,
    )
    lambda_one_residuals = [
        lambda_one_evaluation.values[root] for root in roots
    ]
    return {
        "schema": "plane-jc.f2-75-125-formal-homotopy.v1",
        "status": (
            "exact-finite-field-formal-jet"
            if obstruction is None
            else "exact-finite-field-formal-obstruction"
        ),
        "source_circuit_artifact": str(
            NONLINEAR_ARTIFACT.relative_to(ROOT)
        ),
        "field": {"prime": prime, "rho": rho, "y": y},
        "homotopy_equation": "F(x(lambda))=(1-lambda)*F(x0)",
        "equation_count_with_localization": len(roots),
        "initial_residual_nonzero_count": sum(bool(value) for value in initial_values),
        "Jacobian_rank": cokernel.rank,
        "Jacobian_cokernel_dimension": len(cokernel.functionals),
        "first_tangent_choice": {
            "chart": "spacing-four staircase",
            "free_parameter_choice": "all 112 free parameters zero",
            "active_obstruction_variables_Q9_Q5_zero": True,
        },
        "higher_order_gauge": {
            "name": (
                "seven-pole-coordinate-zero" if regular_gauge else "pivot-default"
            ),
            "variables_prescribed_zero_from_order_two": (
                list(REGULAR_GAUGE_VARIABLES) if regular_gauge else []
            ),
        },
        "requested_order": maximum_order,
        "achieved_order": achieved_order,
        "obstruction": obstruction,
        "steps": steps,
        "a_series_low_to_high": a_series,
        "a_is_formal_unit": bool(a_series[0]),
        "variable_series": compact_variable_series,
        "variable_series_digest_sha256": digest(compact_variable_series),
        "root_series_digest_sha256": digest(root_series),
        "truncated_lambda_one_evaluation": {
            "point_support_size": len(truncated_lambda_one_point),
            "point_digest_sha256": digest(truncated_lambda_one_point),
            "a_value": lambda_one_evaluation.values[a_node],
            "nonzero_total": sum(bool(value) for value in lambda_one_residuals),
            "nonzero_by_group": {
                name: sum(
                    bool(value)
                    for value in lambda_one_residuals[start:end]
                )
                for name, (start, end) in groups.items()
            },
            "residual_digest_sha256": digest(lambda_one_residuals),
            "is_exact_modular_point": not any(lambda_one_residuals),
            "claim_boundary": (
                "substituting lambda=1 into a truncated formal jet has no "
                "formal-convergence guarantee"
            ),
        },
        "claim_boundary": (
            "a finite formal jet over GF(31) does not specialize safely at "
            "lambda=1 and is not a finite-field point or characteristic-zero lift"
        ),
        "software": {"sympy": sp.__version__},
    }


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--maximum-order", type=int, default=16)
    parser.add_argument("--prime", type=int, default=31)
    parser.add_argument("--rho", type=int, default=14)
    parser.add_argument("--y", type=int, default=3)
    parser.add_argument("--regular-gauge", action="store_true")
    args = parser.parse_args()

    payload = build_payload(
        prime=args.prime,
        rho=args.rho,
        y=args.y,
        maximum_order=args.maximum_order,
        regular_gauge=args.regular_gauge,
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
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError("the pinned formal-homotopy artifact is stale")
    print(
        "F2_FORMAL_HOMOTOPY_"
        + ("JET_PASS" if payload["obstruction"] is None else "OBSTRUCTED")
    )
    print(f"F2_FORMAL_HOMOTOPY_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()
