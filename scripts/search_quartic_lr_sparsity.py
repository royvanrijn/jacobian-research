#!/usr/bin/env python3
"""Bounded exact left-right sparsity search for the ungraded quartic.

This is an experiment, not an absolute minimality proof.  It searches:

* one elementary monomial source shear at a time;
* one elementary monomial target shear at a time;
* the structured two-term t-adic representatives of the essential q-jet;
* rational diagonal left-right scalings for coefficient/collision height.

Every support and height comparison is performed over ``QQ``.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Iterable

import sympy as sp


x, y, z, parameter = sp.symbols("x y z parameter")
variables = (x, y, z)
t = 1 + x * y
q = t**2 * z - sp.Rational(4, 7) * y**2 * (1 + 3 * t)
mapping = (
    -sp.Rational(1, 2) * t * q,
    y - sp.Rational(21, 4) * x * q + 3 * t**2 * x**2 * q**4,
    x * (5 - 3 * t)
    + sp.Rational(7, 4) * x**3 * z
    - sp.Rational(3, 2) * (x * q) ** 4,
)
base_polynomials = tuple(
    sp.Poly(sp.expand(component), *variables) for component in mapping
)
base_support = tuple(len(polynomial.terms()) for polynomial in base_polynomials)
base_degrees = tuple(
    polynomial.total_degree() for polynomial in base_polynomials
)
assert base_support == (7, 51, 38)
assert base_degrees == (7, 26, 24)


def rational_text(value: sp.Expr | Fraction) -> str:
    rational = sp.Rational(value)
    if rational.q == 1:
        return str(rational.p)
    return f"{rational.p}/{rational.q}"


def rational_height(value: sp.Expr | Fraction) -> int:
    rational = sp.Rational(value)
    return max(abs(int(rational.p)), int(rational.q))


def support_at(
    polynomials: tuple[sp.Poly, ...], value: sp.Rational
) -> tuple[int, ...]:
    return tuple(
        sum(
            sp.expand(coefficient.subs(parameter, value)) != 0
            for _, coefficient in polynomial.terms()
        )
        for polynomial in polynomials
    )


def rational_exceptional_values(
    polynomials: tuple[sp.Poly, ...],
) -> set[sp.Rational]:
    values = {sp.Integer(0)}
    for polynomial in polynomials:
        for _, coefficient in polynomial.terms():
            if not coefficient.has(parameter):
                continue
            numerator = sp.cancel(coefficient).as_numer_denom()[0]
            univariate = sp.Poly(numerator, parameter)
            for root in sp.roots(univariate):
                if root.is_Rational:
                    values.add(sp.Rational(root))
    return values


def source_monomials(
    coordinate: int, degree_bound: int
) -> Iterable[tuple[tuple[int, int], sp.Expr]]:
    other_variables = (
        (y, z),
        (x, z),
        (x, y),
    )[coordinate]
    for total_degree in range(1, degree_bound + 1):
        for first_degree in range(total_degree + 1):
            second_degree = total_degree - first_degree
            yield (
                (first_degree, second_degree),
                other_variables[0] ** first_degree
                * other_variables[1] ** second_degree,
            )


def search_elementary_source_shears(
    degree_bound: int,
) -> dict[str, object]:
    records: list[dict[str, object]] = []
    best = {
        "total_support": sum(base_support),
        "support": list(base_support),
        "coordinate": None,
        "monomial_exponents": None,
        "parameter": "0",
    }
    for coordinate in range(3):
        for exponents, monomial in source_monomials(coordinate, degree_bound):
            substitution = {
                variables[coordinate]:
                variables[coordinate] + parameter * monomial
            }
            transformed = tuple(
                sp.Poly(
                    sp.expand(
                        component.subs(substitution, simultaneous=True)
                    ),
                    *variables,
                )
                for component in mapping
            )
            exceptional_values = rational_exceptional_values(transformed)
            local_best = min(
                (
                    sum(support_at(transformed, value)),
                    support_at(transformed, value),
                    value,
                )
                for value in exceptional_values
            )
            record = {
                "coordinate": str(variables[coordinate]),
                "monomial_exponents": list(exponents),
                "generic_support": [
                    len(polynomial.terms()) for polynomial in transformed
                ],
                "rational_exceptional_values": len(exceptional_values),
                "best_total_support": local_best[0],
                "best_support": list(local_best[1]),
                "best_parameter": rational_text(local_best[2]),
            }
            records.append(record)
            if local_best[0] < int(best["total_support"]):
                best = {
                    "total_support": local_best[0],
                    "support": list(local_best[1]),
                    "coordinate": str(variables[coordinate]),
                    "monomial_exponents": list(exponents),
                    "parameter": rational_text(local_best[2]),
                }
    return {"best": best, "records": records}


def sparse_terms(expression: sp.Expr) -> dict[tuple[int, ...], sp.Expr]:
    return dict(sp.Poly(sp.expand(expression), *variables).terms())


def search_elementary_target_shears(
    degree_bound: int,
) -> dict[str, object]:
    records: list[dict[str, object]] = []
    best = {
        "total_support": sum(base_support),
        "support": list(base_support),
        "coordinate": None,
        "monomial_exponents": None,
        "parameter": "0",
    }
    for coordinate in range(3):
        other_coordinates = [
            index for index in range(3) if index != coordinate
        ]
        target_terms = sparse_terms(mapping[coordinate])
        for total_degree in range(1, degree_bound + 1):
            for first_degree in range(total_degree + 1):
                second_degree = total_degree - first_degree
                shear = (
                    mapping[other_coordinates[0]] ** first_degree
                    * mapping[other_coordinates[1]] ** second_degree
                )
                shear_terms = sparse_terms(shear)
                union = set(target_terms) | set(shear_terms)
                candidate_values = {sp.Integer(0)}
                for monomial in set(target_terms) & set(shear_terms):
                    candidate_values.add(
                        sp.cancel(
                            -target_terms[monomial] / shear_terms[monomial]
                        )
                    )
                local_best = None
                for value in candidate_values:
                    new_count = sum(
                        sp.expand(
                            target_terms.get(monomial, 0)
                            + value * shear_terms.get(monomial, 0)
                        )
                        != 0
                        for monomial in union
                    )
                    support = list(base_support)
                    support[coordinate] = new_count
                    candidate = (sum(support), tuple(support), value)
                    if local_best is None or candidate < local_best:
                        local_best = candidate
                assert local_best is not None
                record = {
                    "coordinate": coordinate + 1,
                    "other_coordinate_exponents": [
                        first_degree,
                        second_degree,
                    ],
                    "candidate_values": len(candidate_values),
                    "best_total_support": local_best[0],
                    "best_support": list(local_best[1]),
                    "best_parameter": rational_text(local_best[2]),
                }
                records.append(record)
                if local_best[0] < int(best["total_support"]):
                    best = {
                        "total_support": local_best[0],
                        "support": list(local_best[1]),
                        "coordinate": coordinate + 1,
                        "monomial_exponents": [
                            first_degree,
                            second_degree,
                        ],
                        "parameter": rational_text(local_best[2]),
                    }
    return {"best": best, "records": records}


def search_two_term_tadic_shears(maximum_exponent: int) -> dict[str, object]:
    auxiliary = sp.symbols("auxiliary")
    original_boundary = -sp.Rational(4, 7) * (4 + 3 * auxiliary)
    records: list[dict[str, object]] = []
    for first, second in combinations(range(maximum_exponent + 1), 2):
        coefficient_1, coefficient_2 = sp.symbols(
            "coefficient_1 coefficient_2"
        )
        two_term_boundary = (
            coefficient_1 * auxiliary**first
            + coefficient_2 * auxiliary**second
        )
        equations = (
            sp.Eq(
                two_term_boundary.subs(auxiliary, -1),
                original_boundary.subs(auxiliary, -1),
            ),
            sp.Eq(
                sp.diff(two_term_boundary, auxiliary).subs(auxiliary, -1),
                sp.diff(original_boundary, auxiliary).subs(auxiliary, -1),
            ),
        )
        solutions = sp.solve(
            equations, (coefficient_1, coefficient_2), dict=True
        )
        if not solutions:
            continue
        boundary = sp.expand(two_term_boundary.subs(solutions[0]))
        shear_quotient = sp.cancel(
            (boundary - original_boundary) / (1 + auxiliary) ** 2
        )
        if sp.denom(shear_quotient) != 1:
            continue
        source_shear = sp.expand(
            y**2 * shear_quotient.subs(auxiliary, x * y)
        )
        transformed_q = sp.expand(
            t**2 * z + y**2 * boundary.subs(auxiliary, x * y)
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
        polynomials = tuple(
            sp.Poly(sp.expand(component), *variables)
            for component in transformed_mapping
        )
        support = tuple(len(polynomial.terms()) for polynomial in polynomials)
        degrees = tuple(
            polynomial.total_degree() for polynomial in polynomials
        )
        records.append(
            {
                "exponents": [first, second],
                "boundary": str(boundary),
                "source_shear_quotient": str(shear_quotient),
                "support": list(support),
                "degrees": list(degrees),
                "total_support": sum(support),
            }
        )
    best = min(
        records,
        key=lambda record: (
            int(record["total_support"]),
            sum(record["degrees"]),
        ),
    )
    return {"best": best, "records": records}


def positive_reduced_rationals(bound: int) -> list[Fraction]:
    return sorted(
        {
            Fraction(numerator, denominator)
            for numerator in range(1, bound + 1)
            for denominator in range(1, bound + 1)
            if gcd(numerator, denominator) == 1
        }
    )


def coefficient_encoding() -> list[tuple[Fraction, int, int]]:
    c, d = sp.symbols("c d", nonzero=True)
    generic_q = t**2 * z + c * y**2 * (1 + 3 * t)
    generic_mapping = (
        -sp.Rational(1, 2) * t * generic_q,
        y + 3 * x * generic_q / c + 4 * d * t**2 * x**2 * generic_q**4,
        x * (5 - 3 * t)
        - x**3 * z / c
        - 2 * d * (x * generic_q) ** 4,
    )
    encoding: list[tuple[Fraction, int, int]] = []
    for component in generic_mapping:
        for _, coefficient in sp.Poly(
            sp.expand(component), *variables
        ).terms():
            powers = coefficient.as_powers_dict()
            c_power = int(powers.get(c, 0))
            d_power = int(powers.get(d, 0))
            constant = sp.Rational(
                coefficient / (c**c_power * d**d_power)
            )
            encoding.append(
                (Fraction(int(constant.p), int(constant.q)), c_power, d_power)
            )
    return encoding


def search_rational_scalings(bound: int) -> dict[str, object]:
    encoding = coefficient_encoding()
    collision = (
        (Fraction(0), Fraction(0), Fraction(1)),
        (Fraction(-4, 5), Fraction(9, 4), Fraction(-265, 32)),
        (Fraction(1, 2), Fraction(-3, 2), Fraction(100)),
        (
            Fraction(3, 10),
            Fraction(-29, 6),
            Fraction(-24820, 729),
        ),
    )
    candidates = []
    rationals = positive_reduced_rationals(bound)
    for alpha in rationals:
        for beta in rationals:
            c = alpha**2 * beta * Fraction(-4, 7)
            d = alpha**-3 * beta**-4 * Fraction(3, 4)
            coefficient_height = max(
                rational_height(constant * c**c_power * d**d_power)
                for constant, c_power, d_power in encoding
            )
            collision_height = max(
                rational_height(value)
                for point in collision
                for value in (
                    alpha * point[0],
                    point[1] / alpha,
                    beta * point[2],
                )
            )
            candidates.append(
                {
                    "alpha": rational_text(alpha),
                    "beta": rational_text(beta),
                    "c": rational_text(c),
                    "d": rational_text(d),
                    "coefficient_height": coefficient_height,
                    "collision_height": collision_height,
                    "maximum_height": max(
                        coefficient_height, collision_height
                    ),
                }
            )
    best_balanced = min(
        candidates,
        key=lambda item: (
            item["maximum_height"],
            item["coefficient_height"] + item["collision_height"],
        ),
    )
    best_coefficient = min(
        candidates,
        key=lambda item: (
            item["coefficient_height"],
            item["collision_height"],
        ),
    )
    return {
        "candidate_count": len(candidates),
        "best_balanced": best_balanced,
        "best_coefficient": best_coefficient,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-degree", type=int, default=2)
    parser.add_argument("--target-degree", type=int, default=2)
    parser.add_argument("--tadic-max-exponent", type=int, default=12)
    parser.add_argument("--scaling-bound", type=int, default=16)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = {
        "status": "bounded experiment, not an absolute sparsity theorem",
        "base_support": list(base_support),
        "base_degrees": list(base_degrees),
        "elementary_source_shears": search_elementary_source_shears(
            args.source_degree
        ),
        "elementary_target_shears": search_elementary_target_shears(
            args.target_degree
        ),
        "two_term_tadic_shears": search_two_term_tadic_shears(
            args.tadic_max_exponent
        ),
        "rational_scalings": search_rational_scalings(args.scaling_bound),
        "parameters": {
            "source_degree": args.source_degree,
            "target_degree": args.target_degree,
            "tadic_max_exponent": args.tadic_max_exponent,
            "scaling_bound": args.scaling_bound,
        },
    }
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print("PASS: exact bounded quartic LR sparsity search")
    print("base support:", base_support)
    print(
        "best elementary source shear:",
        result["elementary_source_shears"]["best"],
    )
    print(
        "best elementary target shear:",
        result["elementary_target_shears"]["best"],
    )
    print(
        "best structured t-adic shear:",
        result["two_term_tadic_shears"]["best"],
    )
    print(
        "best balanced scaling:",
        result["rational_scalings"]["best_balanced"],
    )
    if args.output is not None:
        print("wrote:", args.output)


if __name__ == "__main__":
    main()
