#!/usr/bin/env python3
"""Search the smallest oblique coisotropic graphs at the exact scalar gate.

Start from the collision-centred six-variable foundational doubling used by
``search_hc4_mixed_canonical_pivots.py``.  For ordered ``i != j`` put

    K = q_i + rho*p_j,
    L = q_j + rho*p_i,
    H = tau*K*L^2.

The two mixed linear forms Poisson-commute.  The time-one Hamiltonian flow is
therefore polynomial, fixes K and L, and sends the coordinate coisotropic
``p_i=0`` to the nonlinear mixed graph ``p_i + tau*L^2=0`` (up to the
pullback/pushforward sign convention).  Nonzero rho and tau make this the
smallest reciprocal mixed-line enlargement of the coordinate graph family.

For every coordinate in which the pulled-back potential is affine,

    Phi = t*A(w) + B(w),

the script tests the weaker graph-specialized Schur gate

    det Hess_w(B + s*A) |_(s=mu+lambda*A(w)) = constant.

Unequal values modulo the stated prime are exact characteristic-zero
nonconstancy certificates because all specialized polynomials have integer
coefficients.  Agreement at the sampled points would be reported only as a
modular survivor.  No full descended determinant is formed.

This is a bounded parameter search, not a classification of arbitrary
rational rho, tau, lambda, or mu and not a theorem about all coisotropic
embeddings.
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import sympy as sp

try:
    from search_hc4_mixed_canonical_pivots import (
        PRIME,
        SCALAR_LAMBDAS,
        SCALAR_MUS,
        Jet,
        affine_coordinates,
        audit_parent_hessian,
        deterministic_points,
        first_difference,
        p,
        poisson_bracket,
        q,
        scalar_remainder_values,
        split_scalar,
        time_one_pullback,
        variables,
    )
except ModuleNotFoundError:
    from scripts.search_hc4_mixed_canonical_pivots import (
        PRIME,
        SCALAR_LAMBDAS,
        SCALAR_MUS,
        Jet,
        affine_coordinates,
        audit_parent_hessian,
        deterministic_points,
        first_difference,
        p,
        poisson_bracket,
        q,
        scalar_remainder_values,
        split_scalar,
        time_one_pullback,
        variables,
    )


SLOPE_BOX = (-2, -1, 1, 2)
TIME_BOX = (-2, -1, 1, 2)


def run_search() -> dict[str, object]:
    points = deterministic_points(5)
    parent_points = deterministic_points(6)
    chart_rows: list[dict[str, object]] = []
    summary = {
        "charts": 0,
        "affine_pivots": 0,
        "scalar_gate_trials": 0,
        "scalar_gate_modular_survivors": 0,
        "parent_hessian_nonconstant": 0,
        "parent_hessian_modular_survivors": 0,
        "full_determinants_formed": 0,
    }

    for source_index in range(3):
        for reciprocal_index in range(3):
            if source_index == reciprocal_index:
                continue
            for rho, tau in product(SLOPE_BOX, TIME_BOX):
                K = q[source_index] + rho * p[reciprocal_index]
                L = q[reciprocal_index] + rho * p[source_index]
                assert poisson_bracket(K, L) == 0
                hamiltonian = tau * K * L**2
                transformed = time_one_pullback(hamiltonian)
                transformed_polynomial = sp.Poly(
                    transformed, *variables, domain=sp.ZZ
                )
                affine = affine_coordinates(transformed_polynomial)

                summary["charts"] += 1
                summary["affine_pivots"] += len(affine)
                pivot_rows: list[dict[str, object]] = []
                for pivot_index in affine:
                    A, B0, _ = split_scalar(transformed, pivot_index)
                    jets = (Jet.from_poly(A), Jet.from_poly(B0))
                    evaluations = tuple(
                        (
                            jets[0].evaluate(point),
                            jets[1].evaluate(point),
                        )
                        for point in points
                    )
                    trial_rows: list[dict[str, object]] = []
                    for lam, mu in product(SCALAR_LAMBDAS, SCALAR_MUS):
                        summary["scalar_gate_trials"] += 1
                        values = scalar_remainder_values(
                            evaluations, lam, mu
                        )
                        witness = first_difference(values)
                        trial: dict[str, object] = {
                            "lambda": lam,
                            "mu": mu,
                        }
                        if witness is None:
                            summary[
                                "scalar_gate_modular_survivors"
                            ] += 1
                            trial["status"] = "modular_survivor"
                            trial["values_mod_p"] = sorted(set(values))
                        else:
                            trial["status"] = (
                                "exact_characteristic_zero_nonconstancy"
                            )
                            trial["witness_mod_p"] = witness
                        trial_rows.append(trial)
                    pivot_rows.append(
                        {
                            "pivot": str(variables[pivot_index]),
                            "A_term_count": len(A.terms()),
                            "B_term_count": len(B0.terms()),
                            "trials": trial_rows,
                        }
                    )

                # This is a side audit after the requested scalar gate.  A
                # nonlinear symplectic change of independent variables need
                # not preserve a constant Hessian determinant.
                parent_constant, parent_witness = audit_parent_hessian(
                    transformed,
                    parent_points,
                    quadratic_generator=False,
                )
                if parent_constant is False:
                    summary["parent_hessian_nonconstant"] += 1
                else:
                    summary["parent_hessian_modular_survivors"] += 1

                chart_rows.append(
                    {
                        "ordered_pair": [
                            source_index,
                            reciprocal_index,
                        ],
                        "rho": rho,
                        "tau": tau,
                        "hamiltonian": (
                            "tau*(q_i+rho*p_j)*(q_j+rho*p_i)^2"
                        ),
                        "term_count": len(
                            transformed_polynomial.terms()
                        ),
                        "affine_pivots": [
                            str(variables[index]) for index in affine
                        ],
                        "pivot_rows": pivot_rows,
                        "parent_hessian_status": (
                            "exact_characteristic_zero_nonconstancy"
                            if parent_constant is False
                            else "modular_survivor"
                        ),
                        "parent_hessian_witness_mod_p": parent_witness,
                    }
                )

    return {
        "status": "bounded_exact_search",
        "scope": {
            "base": "collision-centred foundational cubic Keller doubling",
            "prime": PRIME,
            "ordered_pairs": [
                [left, right]
                for left in range(3)
                for right in range(3)
                if left != right
            ],
            "rho_box": list(SLOPE_BOX),
            "tau_box": list(TIME_BOX),
            "lambda_box": list(SCALAR_LAMBDAS),
            "mu_box": list(SCALAR_MUS),
            "determinant_gate": (
                "det Hess_w(B+s*A) at s=mu+lambda*A(w)"
            ),
            "full_descended_determinants": "not formed",
        },
        "summary": summary,
        "charts": chart_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete deterministic JSON census",
    )
    arguments = parser.parse_args()

    result = run_search()
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )

    print("HC4_NONCOORDINATE_COISOTROPIC_SCALAR_GATE")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    if arguments.output is not None:
        print(f"artifact={arguments.output}")
    print(
        "SCOPE: bounded exact parameter box; arbitrary rational parameters "
        "and general coisotropic embeddings remain open"
    )


if __name__ == "__main__":
    main()
