#!/usr/bin/env python3
"""Search the unit quadratic--cubic HC(4) commutator box.

For every noncommuting pair

    H1 = (q_i + epsilon_1*p_j)^2,
    H2 = (q_k + epsilon_2*p_l)^3,

this checker forms the canonical group commutator

    T_{-H2} o T_{-H1} o T_{H2} o T_{H1}.

The 18 quadratic and 18 cubic letters give 162 noncommuting pairs.  They
are classified by their factored Poisson bracket before the commutator is
expanded.

Expanding the pulled-back foundational potential is needlessly expensive.
Instead, the transformed six-variable Hessian is evaluated exactly modulo
the good prime 1000003 using

    Hess(Phi o S)
      = DS^T Hess(Phi)(S) DS
        + sum_k partial_k(Phi)(S) Hess(S_k).

Two different determinant values modulo the prime certify that the parent
Hessian determinant is nonconstant over characteristic zero.  Agreement
on the finite deterministic point set is reported only as an unresolved
modular survivor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Sequence

import sympy as sp

import search_hc4_mixed_canonical_pivots as base
import search_hc4_mixed_quadratic_words as words


def compose_flows(
    outer: Sequence[sp.Expr],
    inner: Sequence[sp.Expr],
) -> tuple[sp.Expr, ...]:
    """Return ``outer o inner`` as a polynomial coordinate map."""

    substitution = dict(zip(base.variables, inner, strict=True))
    return tuple(
        sp.expand(expression.subs(substitution, simultaneous=True))
        for expression in outer
    )


def commutator_flow(
    h1: words.Letter,
    h2: words.Letter,
) -> tuple[sp.Expr, ...]:
    """Return ``T_-H2 o T_-H1 o T_H2 o T_H1``."""

    flow = h1.flow
    for shear in (
        h2.flow,
        words.flow_of(-h1.hamiltonian),
        words.flow_of(-h2.hamiltonian),
    ):
        flow = compose_flows(shear, flow)
    return flow


def modular_polynomial(expression: sp.Expr) -> base.ModularPolynomial:
    return base.ModularPolynomial(sp.Poly(expression, *base.variables))


BASE_GRADIENT = tuple(
    modular_polynomial(sp.diff(base.base_potential, variable))
    for variable in base.variables
)
BASE_HESSIAN = tuple(
    tuple(
        modular_polynomial(
            sp.diff(base.base_potential, left, right)
        )
        for right in base.variables
    )
    for left in base.variables
)


def modular_flow_jet(
    flow: Sequence[sp.Expr],
) -> tuple[
    tuple[base.ModularPolynomial, ...],
    tuple[tuple[base.ModularPolynomial, ...], ...],
    tuple[
        tuple[tuple[base.ModularPolynomial, ...], ...],
        ...,
    ],
]:
    values = tuple(modular_polynomial(expression) for expression in flow)
    jacobian = tuple(
        tuple(
            modular_polynomial(sp.diff(expression, variable))
            for variable in base.variables
        )
        for expression in flow
    )
    hessians = tuple(
        tuple(
            tuple(
                modular_polynomial(
                    sp.diff(expression, left, right)
                )
                for right in base.variables
            )
            for left in base.variables
        )
        for expression in flow
    )
    return values, jacobian, hessians


def transformed_hessian_determinant(
    jet: tuple[
        tuple[base.ModularPolynomial, ...],
        tuple[tuple[base.ModularPolynomial, ...], ...],
        tuple[
            tuple[tuple[base.ModularPolynomial, ...], ...],
            ...,
        ],
    ],
    point: Sequence[int],
) -> int:
    """Evaluate ``det Hess(base_potential o flow)`` modulo ``base.PRIME``."""

    values, jacobian, hessians = jet
    image = tuple(entry.evaluate(point) for entry in values)
    jacobian_value = tuple(
        tuple(entry.evaluate(point) for entry in row)
        for row in jacobian
    )
    hessian_values = tuple(
        tuple(
            tuple(entry.evaluate(point) for entry in row)
            for row in matrix
        )
        for matrix in hessians
    )
    gradient_base = tuple(
        entry.evaluate(image) for entry in BASE_GRADIENT
    )
    hessian_base = tuple(
        tuple(entry.evaluate(image) for entry in row)
        for row in BASE_HESSIAN
    )

    size = len(base.variables)
    transformed = [[0 for _ in range(size)] for _ in range(size)]
    for left in range(size):
        for right in range(size):
            congruence = sum(
                jacobian_value[image_left][left]
                * hessian_base[image_left][image_right]
                * jacobian_value[image_right][right]
                for image_left in range(size)
                for image_right in range(size)
            )
            correction = sum(
                gradient_base[image_index]
                * hessian_values[image_index][left][right]
                for image_index in range(size)
            )
            transformed[left][right] = (
                congruence + correction
            ) % base.PRIME
    return base.determinant_mod(transformed)


def preclassification() -> tuple[dict[str, object], ...]:
    """Classify the noncommuting supports before composing their flows."""

    census: Counter[tuple[str, int, int]] = Counter()
    commuting = 0
    for h1, h2 in product(
        words.letters(degree=2, tau_box=(1,)),
        words.letters(degree=3, tau_box=(1,)),
    ):
        bracket = base.poisson_bracket(
            h1.hamiltonian, h2.hamiltonian
        )
        if bracket == 0:
            commuting += 1
            continue
        signature = words.bracket_signature(h1, h2)
        census[
            (
                str(signature["incidence"]),
                int(signature["linear_pairing"]),
                int(signature["bracket_coefficient"]),
            )
        ] += 1

    assert commuting == 162
    assert sum(census.values()) == 162
    return tuple(
        {
            "degree_pair": [2, 3],
            "incidence": incidence,
            "linear_pairing": pairing,
            "bracket_coefficient": coefficient,
            "factorization": f"{coefficient}*L1*L2^2",
            "pairs": count,
        }
        for (incidence, pairing, coefficient), count in sorted(
            census.items()
        )
    )


def flow_record(flow: Sequence[sp.Expr]) -> dict[str, object]:
    polynomials = tuple(
        sp.Poly(expression, *base.variables) for expression in flow
    )
    return {
        "coordinate_degrees": [
            int(polynomial.total_degree()) for polynomial in polynomials
        ],
        "coordinate_term_counts": [
            len(polynomial.terms()) for polynomial in polynomials
        ],
        "maximum_degree": max(
            int(polynomial.total_degree()) for polynomial in polynomials
        ),
        "total_terms": sum(
            len(polynomial.terms()) for polynomial in polynomials
        ),
    }


def run(point_count: int = 7) -> dict[str, object]:
    quadratic = words.letters(degree=2, tau_box=(1,))
    cubic = words.letters(degree=3, tau_box=(1,))
    points = base.deterministic_points(6, count=point_count)

    rows: list[dict[str, object]] = []
    map_keys: set[object] = set()
    incidence_census: Counter[str] = Counter()
    parent_census: Counter[str] = Counter()
    bracket_type_parent_census: Counter[tuple[str, str]] = Counter()
    noncommuting = 0

    for h1, h2 in product(quadratic, cubic):
        if base.poisson_bracket(h1.hamiltonian, h2.hamiltonian) == 0:
            continue
        noncommuting += 1
        if noncommuting == 1 or noncommuting % 20 == 0:
            print(
                f"progress={noncommuting}/162 "
                f"pair={h1.letter_id}__{h2.letter_id}",
                flush=True,
            )

        signature = words.bracket_signature(h1, h2)
        incidence = str(signature["incidence"])
        incidence_census[incidence] += 1

        flow = commutator_flow(h1, h2)
        key = words.polynomial_map_key(flow)
        duplicate = key in map_keys
        map_keys.add(key)

        jet = modular_flow_jet(flow)
        evaluations: list[dict[str, object]] = []
        first_value: int | None = None
        witness: dict[str, object] | None = None
        for point in points:
            determinant = transformed_hessian_determinant(jet, point)
            evaluation = {
                "point": list(point),
                "determinant_mod_prime": determinant,
            }
            evaluations.append(evaluation)
            if first_value is None:
                first_value = determinant
            elif determinant != first_value:
                witness = {
                    "first": evaluations[0],
                    "second": evaluation,
                }
                break

        if witness is None:
            parent_class = "unresolved_modular_survivor"
        else:
            parent_class = "nonconstant"
        parent_census[parent_class] += 1
        bracket_type_parent_census[(incidence, parent_class)] += 1

        rows.append(
            {
                "pair_id": f"{h1.letter_id}__{h2.letter_id}",
                "quadratic_letter": h1.letter_id,
                "cubic_letter": h2.letter_id,
                "degree_pair": [2, 3],
                "composition": "T_-H2 o T_-H1 o T_H2 o T_H1",
                "poisson_bracket": signature,
                "duplicate_commutator_map": duplicate,
                "flow": flow_record(flow),
                "parent_hessian": {
                    "classification": parent_class,
                    "certificate": (
                        "unequal_determinants_mod_prime"
                        if witness is not None
                        else "finite_modular_agreement_only"
                    ),
                    "witness": witness,
                    "evaluations": evaluations,
                },
                "affine_pivot_subspace": {
                    "classification": "not_audited_after_parent_gate",
                    "reason": (
                        "parent Hessian determinant is nonconstant"
                        if witness is not None
                        else "requires exact survivor expansion"
                    ),
                },
                "reduced_hessian_corank": {
                    "classification": "not_audited_after_parent_gate",
                    "reason": (
                        "parent Hessian determinant is nonconstant"
                        if witness is not None
                        else "requires exact survivor expansion"
                    ),
                },
            }
        )

    assert noncommuting == 162
    assert len(rows) == 162
    assert len(map_keys) <= len(rows)

    return {
        "schema_version": 1,
        "claim_type": "bounded_exact_computation",
        "family": "unit_quadratic_cubic_group_commutators",
        "search_box": {
            "quadratic_letters": 18,
            "cubic_letters": 18,
            "coefficient_box": [1],
            "epsilon_box": [-1, 1],
            "raw_pairs": 324,
            "commuting_pairs_excluded": 162,
            "noncommuting_pairs": 162,
            "commutator": "T_-H2 o T_-H1 o T_H2 o T_H1",
            "point_count": point_count,
            "prime": base.PRIME,
        },
        "method": {
            "classification_before_expansion": [
                "Poisson-bracket incidence",
                "linear pairing",
                "factored bracket coefficient",
                "degree pair",
            ],
            "parent_hessian_formula": (
                "DS^T Hess(Phi)(S) DS + "
                "sum_k partial_k(Phi)(S) Hess(S_k)"
            ),
            "certificate_logic": (
                "unequal modular determinant values at two integer "
                "points imply characteristic-zero nonconstancy"
            ),
            "survivor_policy": (
                "finite modular agreement is unresolved and must be "
                "expanded for an exact identity check"
            ),
        },
        "preclassification": list(preclassification()),
        "summary": {
            "raw_pairs": 324,
            "commuting_pairs_excluded": 162,
            "noncommuting_pairs": noncommuting,
            "unique_commutator_maps": len(map_keys),
            "duplicate_commutator_maps": len(rows) - len(map_keys),
            "incidence_census": dict(sorted(incidence_census.items())),
            "parent_hessian_census": dict(sorted(parent_census.items())),
            "bracket_type_parent_census": {
                f"{incidence}:{parent_class}": count
                for (incidence, parent_class), count in sorted(
                    bracket_type_parent_census.items()
                )
            },
            "parent_survivors_requiring_expansion": parent_census[
                "unresolved_modular_survivor"
            ],
        },
        "pairs": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete commutator census as JSON",
    )
    parser.add_argument(
        "--point-count",
        type=int,
        default=7,
        help="deterministic modular points per unresolved word",
    )
    args = parser.parse_args()

    result = run(args.point_count)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        print(f"artifact={args.output}")
        print(
            "artifact_sha256="
            f"{hashlib.sha256(encoded.encode()).hexdigest()}"
        )
    print("HC4_MIXED_COMMUTATOR_WORD_SUMMARY")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
