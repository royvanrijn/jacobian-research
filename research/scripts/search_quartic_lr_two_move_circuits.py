#!/usr/bin/env python3
"""Bounded two-move circuit search for the ungraded quartic.

This exact/rational computation has two independent parts.

1. Precompose by one elementary monomial source shear with a rational
   parameter in a declared height box, then optimize one linear triangular
   target cleanup over every exact rational coefficient ratio.
2. Exhaust all three-term representatives of the essential q-boundary jet
   through a declared exponent.  The two Hermite constraints leave one
   parameter; every rational exceptional value of every expanded coefficient
   is tested exactly.  This structured family contains multi-monomial
   z-shears which can be decomposed into commuting elementary source moves.

The result is bounded search evidence, not an absolute sparsity theorem.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy.polys.polytools import ground_roots

from search_quartic_lr_sparsity import (
    base_support,
    mapping,
    parameter,
    rational_exceptional_values,
    rational_text,
    source_monomials,
    support_at,
    t,
    variables,
    x,
    y,
    z,
)


def signed_reduced_rationals(bound: int) -> tuple[sp.Rational, ...]:
    return tuple(
        sorted(
            {
                sp.Rational(numerator, denominator)
                for numerator in range(-bound, bound + 1)
                for denominator in range(1, bound + 1)
                if numerator and gcd(abs(numerator), denominator) == 1
            }
        )
    )


def polynomial_terms(expression: sp.Expr) -> dict[tuple[int, ...], sp.Expr]:
    return dict(sp.Poly(sp.expand(expression), *variables).terms())


def best_linear_target_cleanup(
    transformed_mapping: tuple[sp.Expr, ...],
) -> dict[str, object]:
    terms = tuple(polynomial_terms(component) for component in transformed_mapping)
    intermediate_support = tuple(len(component) for component in terms)
    best = {
        "total_support": sum(intermediate_support),
        "support": list(intermediate_support),
        "target_coordinate": None,
        "added_coordinate": None,
        "target_parameter": "0",
    }
    for target_coordinate in range(3):
        for added_coordinate in range(3):
            if target_coordinate == added_coordinate:
                continue
            target_terms = terms[target_coordinate]
            added_terms = terms[added_coordinate]
            common = set(target_terms) & set(added_terms)
            union = set(target_terms) | set(added_terms)
            candidate_values = {sp.Integer(0)}
            candidate_values.update(
                sp.cancel(
                    -target_terms[monomial] / added_terms[monomial]
                )
                for monomial in common
            )
            for value in candidate_values:
                new_count = sum(
                    sp.expand(
                        target_terms.get(monomial, 0)
                        + value * added_terms.get(monomial, 0)
                    )
                    != 0
                    for monomial in union
                )
                support = list(intermediate_support)
                support[target_coordinate] = new_count
                total_support = sum(support)
                if total_support < int(best["total_support"]):
                    best = {
                        "total_support": total_support,
                        "support": support,
                        "target_coordinate": target_coordinate + 1,
                        "added_coordinate": added_coordinate + 1,
                        "target_parameter": rational_text(value),
                    }
    best["intermediate_support"] = list(intermediate_support)
    best["intermediate_total_support"] = sum(intermediate_support)
    return best


def source_target_task(
    task: tuple[int, tuple[int, int], str],
) -> dict[str, object]:
    coordinate, exponents, source_parameter_text = task
    source_parameter = sp.Rational(source_parameter_text)
    other_variables = (
        (y, z),
        (x, z),
        (x, y),
    )[coordinate]
    monomial = (
        other_variables[0] ** exponents[0]
        * other_variables[1] ** exponents[1]
    )
    substitution = {
        variables[coordinate]:
        variables[coordinate] + source_parameter * monomial
    }
    transformed_mapping = tuple(
        sp.expand(component.subs(substitution, simultaneous=True))
        for component in mapping
    )
    result = best_linear_target_cleanup(transformed_mapping)
    result.update(
        {
            "source_coordinate": str(variables[coordinate]),
            "source_monomial_exponents": list(exponents),
            "source_parameter": source_parameter_text,
        }
    )
    return result


u, jet_parameter = sp.symbols("u jet_parameter")
unit = 1 + u
original_boundary = -sp.Rational(4, 7) * (4 + 3 * u)
boundary_value = -sp.Rational(4, 7)
negative_boundary_derivative = sp.Rational(12, 7)


def univariate_coefficients(expression: sp.Expr) -> list[sp.Expr]:
    return [
        coefficient
        for _, coefficient in sp.Poly(sp.expand(expression), u).terms()
        if coefficient != 0
    ]


def compressed_component_coefficients(
    boundary: sp.Expr,
) -> tuple[list[sp.Expr], list[sp.Expr], list[sp.Expr], sp.Expr]:
    """Return disjoint coefficient blocks for the three affine components."""

    source_shear_quotient = sp.cancel(
        (boundary - original_boundary) / unit**2
    )
    assert sp.denom(source_shear_quotient) == 1

    first = (
        univariate_coefficients(unit**3)
        + univariate_coefficients(unit * boundary)
    )
    second = (
        [sp.Integer(1)]
        + univariate_coefficients(unit**2)
        + univariate_coefficients(boundary)
    )
    third = (
        [sp.Integer(1)] * 3
        + univariate_coefficients(source_shear_quotient)
    )
    boundary_powers = [sp.Integer(1)]
    for _ in range(4):
        boundary_powers.append(sp.expand(boundary_powers[-1] * boundary))
    for z_degree in range(5):
        boundary_power = boundary_powers[4 - z_degree]
        second.extend(
            univariate_coefficients(
                unit ** (2 + 2 * z_degree) * boundary_power
            )
        )
        third.extend(
            univariate_coefficients(
                unit ** (2 * z_degree) * boundary_power
            )
        )
    return first, second, third, source_shear_quotient


def three_term_boundary(exponents: tuple[int, int, int]) -> sp.Expr:
    first, second, third = exponents
    signed_parameter = jet_parameter * (-1) ** third
    signed_second_coefficient = (
        negative_boundary_derivative
        - first * boundary_value
        + (first - third) * signed_parameter
    ) / (second - first)
    signed_first_coefficient = (
        boundary_value
        - signed_parameter
        - signed_second_coefficient
    )
    first_coefficient = signed_first_coefficient * (-1) ** first
    second_coefficient = signed_second_coefficient * (-1) ** second
    return sp.expand(
        first_coefficient * u**first
        + second_coefficient * u**second
        + jet_parameter * u**third
    )


def structured_jet_task(
    exponents: tuple[int, int, int],
) -> dict[str, object]:
    boundary = three_term_boundary(exponents)
    first, second, third, source_shear_quotient = (
        compressed_component_coefficients(boundary)
    )
    component_coefficients = (first, second, third)
    exceptional_values = {sp.Integer(0)}
    parameter_coefficients = {
        sp.factor(coefficient)
        for component in component_coefficients
        for coefficient in component
        if coefficient.has(jet_parameter)
    }
    for coefficient in parameter_coefficients:
        exceptional_values.update(
            ground_roots(sp.Poly(coefficient, jet_parameter))
        )

    candidates = []
    for value in exceptional_values:
        support = tuple(
            sum(
                coefficient.subs(jet_parameter, value) != 0
                for coefficient in component
            )
            for component in component_coefficients
        )
        is_identity = (
            sp.expand(
                source_shear_quotient.subs(jet_parameter, value)
            )
            == 0
        )
        candidates.append(
            {
                "parameter": rational_text(value),
                "support": list(support),
                "total_support": sum(support),
                "identity": is_identity,
            }
        )
    best = min(
        candidates,
        key=lambda item: (
            item["total_support"],
            item["identity"],
            item["parameter"],
        ),
    )
    nonidentity = [item for item in candidates if not item["identity"]]
    best_nonidentity = min(
        nonidentity,
        key=lambda item: (item["total_support"], item["parameter"]),
    )
    return {
        "exponents": list(exponents),
        "boundary": str(boundary),
        "exceptional_values": len(exceptional_values),
        "best": best,
        "best_nonidentity": best_nonidentity,
    }


def full_structured_support(
    exponents: tuple[int, int, int], value: sp.Rational
) -> tuple[int, ...]:
    """Independent full x,y,z expansion for compressed-ledger spot checks."""

    transformed_mapping = full_structured_mapping(exponents, value)
    return tuple(
        len(sp.Poly(component, *variables).terms())
        for component in transformed_mapping
    )


def full_structured_mapping(
    exponents: tuple[int, int, int], value: sp.Rational
) -> tuple[sp.Expr, ...]:
    """Return the fully expanded map for one structured source jet."""

    boundary = three_term_boundary(exponents).subs(jet_parameter, value)
    source_shear_quotient = sp.cancel(
        (boundary - original_boundary) / unit**2
    )
    source_shear = sp.expand(
        y**2 * source_shear_quotient.subs(u, x * y)
    )
    transformed_q = sp.expand(
        t**2 * z + y**2 * boundary.subs(u, x * y)
    )
    transformed_mapping = (
        -sp.Rational(1, 2) * t * transformed_q,
        y
        - sp.Rational(21, 4) * x * transformed_q
        + 3 * t**2 * x**2 * transformed_q**4,
        x * (5 - 3 * t)
        + sp.Rational(7, 4) * x**3 * (z + source_shear)
        - sp.Rational(3, 2) * (x * transformed_q) ** 4,
    )
    return tuple(sp.expand(component) for component in transformed_mapping)


def best_polynomial_target_cleanup(
    transformed_mapping: tuple[sp.Expr, ...], degree_bound: int
) -> dict[str, object]:
    """Optimize one monomial target shear through ``degree_bound``."""

    terms = tuple(polynomial_terms(component) for component in transformed_mapping)
    starting_support = tuple(len(component) for component in terms)
    nonidentity_candidates: list[dict[str, object]] = []
    for target_coordinate in range(3):
        other_coordinates = [
            coordinate
            for coordinate in range(3)
            if coordinate != target_coordinate
        ]
        for total_degree in range(1, degree_bound + 1):
            for first_degree in range(total_degree + 1):
                second_degree = total_degree - first_degree
                shear = sp.expand(
                    transformed_mapping[other_coordinates[0]] ** first_degree
                    * transformed_mapping[other_coordinates[1]] ** second_degree
                )
                shear_terms = polynomial_terms(shear)
                target_terms = terms[target_coordinate]
                common = set(target_terms) & set(shear_terms)
                union = set(target_terms) | set(shear_terms)
                candidate_values = {
                    sp.cancel(
                        -target_terms[monomial] / shear_terms[monomial]
                    )
                    for monomial in common
                }
                for value in candidate_values:
                    new_count = sum(
                        sp.expand(
                            target_terms.get(monomial, 0)
                            + value * shear_terms.get(monomial, 0)
                        )
                        != 0
                        for monomial in union
                    )
                    support = list(starting_support)
                    support[target_coordinate] = new_count
                    nonidentity_candidates.append(
                        {
                            "total_support": sum(support),
                            "support": support,
                            "target_coordinate": target_coordinate + 1,
                            "other_coordinate_exponents": [
                                first_degree,
                                second_degree,
                            ],
                            "target_parameter": rational_text(value),
                        }
                    )
    best_nonidentity = min(
        nonidentity_candidates,
        key=lambda item: (
            item["total_support"],
            item["target_coordinate"],
            item["other_coordinate_exponents"],
            item["target_parameter"],
        ),
    )
    return {
        "starting_support": list(starting_support),
        "degree_bound": degree_bound,
        "best_nonidentity": best_nonidentity,
        "improves_starting_support": (
            best_nonidentity["total_support"] < sum(starting_support)
        ),
    }


def best_second_source_cleanup(
    transformed_mapping: tuple[sp.Expr, ...],
    degree_bound: int,
    inverse_move: tuple[int, tuple[int, int], sp.Rational],
) -> dict[str, object]:
    """Test every rational exceptional second elementary source shear."""

    candidates: list[dict[str, object]] = []
    inverse_record = None
    for coordinate in range(3):
        for exponents, monomial in source_monomials(coordinate, degree_bound):
            transformed = tuple(
                sp.Poly(
                    sp.expand(
                        component.subs(
                            {
                                variables[coordinate]:
                                variables[coordinate] + parameter * monomial
                            },
                            simultaneous=True,
                        )
                    ),
                    *variables,
                )
                for component in transformed_mapping
            )
            for value in rational_exceptional_values(transformed):
                if value == 0:
                    continue
                support = support_at(transformed, value)
                record = {
                    "total_support": sum(support),
                    "support": list(support),
                    "source_coordinate": str(variables[coordinate]),
                    "source_monomial_exponents": list(exponents),
                    "source_parameter": rational_text(value),
                }
                if (coordinate, exponents, value) == inverse_move:
                    inverse_record = record
                else:
                    candidates.append(record)
    assert inverse_record is not None
    best_noninverse = min(
        candidates,
        key=lambda item: (
            item["total_support"],
            item["source_coordinate"],
            item["source_monomial_exponents"],
            item["source_parameter"],
        ),
    )
    return {
        "degree_bound": degree_bound,
        "inverse_move": inverse_record,
        "best_noninverse": best_noninverse,
    }


def structured_near_miss_cleanup(maximum_exponent: int) -> dict[str, object]:
    """Add one further structured monomial to the 97-term source jet."""

    cleanup_parameter = sp.symbols("cleanup_parameter")
    records = []
    for exponent in range(1, maximum_exponent + 1):
        source_shear_quotient = (
            sp.Rational(16, 7) + cleanup_parameter * u**exponent
        )
        boundary = sp.expand(
            original_boundary + unit**2 * source_shear_quotient
        )
        components = compressed_component_coefficients(boundary)[:3]
        exceptional_values = set()
        for coefficient in {
            sp.factor(entry)
            for component in components
            for entry in component
            if entry.has(cleanup_parameter)
        }:
            exceptional_values.update(
                ground_roots(sp.Poly(coefficient, cleanup_parameter))
            )
        candidates = []
        for value in exceptional_values:
            if value == 0:
                continue
            support = tuple(
                sum(
                    coefficient.subs(cleanup_parameter, value) != 0
                    for coefficient in component
                )
                for component in components
            )
            candidates.append(
                {
                    "parameter": rational_text(value),
                    "support": list(support),
                    "total_support": sum(support),
                }
            )
        best_nonidentity = (
            min(
                candidates,
                key=lambda item: (
                    item["total_support"],
                    item["parameter"],
                ),
            )
            if candidates
            else None
        )
        records.append(
            {
                "exponent": exponent,
                "rational_nonzero_exceptional_values": len(candidates),
                "best_nonidentity": best_nonidentity,
            }
        )
    realized = [
        record["best_nonidentity"]
        | {"exponent": record["exponent"]}
        for record in records
        if record["best_nonidentity"] is not None
    ]
    best = min(
        realized,
        key=lambda item: (
            item["total_support"],
            item["exponent"],
            item["parameter"],
        ),
    )
    return {
        "maximum_exponent": maximum_exponent,
        "best_nonidentity": best,
        "records": records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-degree", type=int, default=2)
    parser.add_argument("--source-parameter-bound", type=int, default=4)
    parser.add_argument("--jet-max-exponent", type=int, default=12)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def source_target_tasks(
    source_degree: int, source_parameter_bound: int
) -> Iterable[tuple[int, tuple[int, int], str]]:
    parameters = signed_reduced_rationals(source_parameter_bound)
    for coordinate in range(3):
        for exponents, _ in source_monomials(coordinate, source_degree):
            for source_parameter in parameters:
                yield coordinate, exponents, rational_text(source_parameter)


def main() -> None:
    args = parse_args()
    circuit_tasks = list(
        source_target_tasks(
            args.source_degree, args.source_parameter_bound
        )
    )
    jet_tasks = list(combinations(range(args.jet_max_exponent + 1), 3))
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        circuit_records = list(
            executor.map(source_target_task, circuit_tasks, chunksize=4)
        )
        jet_records = list(
            executor.map(structured_jet_task, jet_tasks, chunksize=2)
        )

    best_circuit = min(
        circuit_records,
        key=lambda item: (
            item["total_support"],
            item["intermediate_total_support"],
            item["source_coordinate"],
            item["source_monomial_exponents"],
            item["source_parameter"],
        ),
    )
    best_jet = min(
        (
            record["best_nonidentity"] | {"exponents": record["exponents"]}
            for record in jet_records
        ),
        key=lambda item: (
            item["total_support"],
            item["exponents"],
            item["parameter"],
        ),
    )

    # Verify the compressed support ledger by full expansion at the identity
    # and at the winning nonidentity structured jet.
    assert full_structured_support((0, 1, 2), sp.Integer(0)) == base_support
    winning_exponents = tuple(int(value) for value in best_jet["exponents"])
    winning_parameter = sp.Rational(best_jet["parameter"])
    assert (
        full_structured_support(winning_exponents, winning_parameter)
        == tuple(best_jet["support"])
    )
    winning_mapping = full_structured_mapping(
        winning_exponents, winning_parameter
    )
    near_miss_target_cleanup = best_polynomial_target_cleanup(
        winning_mapping, degree_bound=3
    )
    near_miss_source_cleanup = best_second_source_cleanup(
        winning_mapping,
        degree_bound=2,
        inverse_move=(2, (0, 2), -sp.Rational(16, 7)),
    )
    near_miss_structured_cleanup = structured_near_miss_cleanup(
        maximum_exponent=args.jet_max_exponent
    )
    assert not near_miss_target_cleanup["improves_starting_support"]
    assert (
        near_miss_source_cleanup["inverse_move"]["support"]
        == list(base_support)
    )
    assert (
        near_miss_source_cleanup["best_noninverse"]["total_support"]
        > sum(base_support)
    )
    assert (
        near_miss_structured_cleanup["best_nonidentity"]["total_support"]
        > sum(base_support)
    )

    result = {
        "status": "bounded experiment, not an absolute sparsity theorem",
        "base_support": list(base_support),
        "parameters": {
            "source_degree": args.source_degree,
            "source_parameter_bound": args.source_parameter_bound,
            "jet_max_exponent": args.jet_max_exponent,
            "workers": args.workers,
        },
        "source_target_circuits": {
            "task_count": len(circuit_tasks),
            "best": best_circuit,
            "records": circuit_records,
        },
        "three_term_structured_jets": {
            "task_count": len(jet_tasks),
            "best_nonidentity": best_jet,
            "records": jet_records,
        },
        "near_miss_cleanup": {
            "target_monomial_shears": near_miss_target_cleanup,
            "second_elementary_source_shears": near_miss_source_cleanup,
            "second_structured_source_shears": (
                near_miss_structured_cleanup
            ),
        },
    }
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print("PASS: bounded two-move quartic circuit search")
    print("base support:", base_support)
    print("source-to-target circuits:", len(circuit_tasks))
    print("best source-to-target circuit:", best_circuit)
    print("structured three-term jets:", len(jet_tasks))
    print("best nonidentity structured jet:", best_jet)
    print("near-miss target cleanup:", near_miss_target_cleanup)
    print("near-miss second-source cleanup:", near_miss_source_cleanup)
    print(
        "near-miss structured cleanup:",
        near_miss_structured_cleanup["best_nonidentity"],
    )
    if args.output is not None:
        print("wrote:", args.output)


if __name__ == "__main__":
    main()
