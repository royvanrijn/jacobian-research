#!/usr/bin/env python3
"""Search a direct one-monomial sparse ansatz in the second Weyl algebra.

Use PBW normal order

    x1^a x2^b d1^c d2^d,       [d_i, x_j] = delta_ij.

For each canonical generator Y_i in (x1, x2, d1, d2), the searched image is

    Y_i + c_i M_i,

where either c_i M_i is absent or M_i is one PBW monomial of Bernstein
degree 2 through ``--max-degree``.  The coefficient alphabet is finite and
declared on the command line.  Thus this is an exhaustive search only in
that finite ansatz.

Every returned tuple satisfies the six Weyl relations by direct PBW
multiplication.  A second calculation grows the span of ordered words

    P1^a P2^b Q1^c Q2^d

and records the first word degree at which all four ambient generators are
recovered.  Recovery proves that the tuple generates A_2.  Failure to recover
within the bound is only an unresolved bounded result, never a proof that the
generated subalgebra is proper.

The same calculation may be run modulo a prime.  The script separately tests
the literal least-nonnegative-residue integer lift.  Failure of that lift is
not a proof that no other characteristic-zero lift exists.  It also computes
the p-th powers, checks whether they lie in the ambient center, and performs
the same bounded generation test on the induced center tuple.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from math import comb, prod
from pathlib import Path
import platform
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "sparse_a2_weyl_endomorphism_search.json"
)

Exponent = tuple[int, int, int, int]
Polynomial = dict[Exponent, int]
ZERO_EXPONENT: Exponent = (0, 0, 0, 0)
GENERATOR_EXPONENTS: tuple[Exponent, ...] = (
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
)
GENERATOR_NAMES = ("P1", "P2", "Q1", "Q2")


def reduce_coefficient(value: int, prime: int | None) -> int:
    return value if prime is None else value % prime


def clean(poly: Polynomial, prime: int | None) -> Polynomial:
    answer = {}
    for monomial, coefficient in poly.items():
        coefficient = reduce_coefficient(coefficient, prime)
        if coefficient:
            answer[monomial] = coefficient
    return answer


def add(
    left: Polynomial,
    right: Polynomial,
    prime: int | None,
    right_scale: int = 1,
) -> Polynomial:
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, 0) + right_scale * coefficient
    return clean(answer, prime)


def scale(poly: Polynomial, coefficient: int, prime: int | None) -> Polynomial:
    return clean(
        {monomial: coefficient * value for monomial, value in poly.items()},
        prime,
    )


def falling(value: int, length: int) -> int:
    return prod(range(value - length + 1, value + 1))


def multiply(
    left: Polynomial,
    right: Polynomial,
    prime: int | None,
) -> Polynomial:
    """PBW multiplication with x variables to the left of d variables."""
    answer: Polynomial = {}
    for (a1, a2, c1, c2), left_coefficient in left.items():
        for (b1, b2, e1, e2), right_coefficient in right.items():
            for k1 in range(min(c1, b1) + 1):
                factor1 = comb(c1, k1) * falling(b1, k1)
                for k2 in range(min(c2, b2) + 1):
                    factor2 = comb(c2, k2) * falling(b2, k2)
                    monomial = (
                        a1 + b1 - k1,
                        a2 + b2 - k2,
                        c1 + e1 - k1,
                        c2 + e2 - k2,
                    )
                    coefficient = (
                        left_coefficient
                        * right_coefficient
                        * factor1
                        * factor2
                    )
                    answer[monomial] = answer.get(monomial, 0) + coefficient
    return clean(answer, prime)


def commutator(
    left: Polynomial,
    right: Polynomial,
    prime: int | None,
) -> Polynomial:
    return add(
        multiply(left, right, prime),
        multiply(right, left, prime),
        prime,
        right_scale=-1,
    )


def power(poly: Polynomial, exponent: int, prime: int | None) -> Polynomial:
    answer = {ZERO_EXPONENT: 1}
    factor = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, factor, prime)
        remaining >>= 1
        if remaining:
            factor = multiply(factor, factor, prime)
    return answer


def compositions(total: int, parts: int) -> Iterable[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def monomials_between(minimum: int, maximum: int) -> list[Exponent]:
    answer = []
    for degree in range(minimum, maximum + 1):
        answer.extend(compositions(degree, 4))
    return answer


@dataclass(frozen=True)
class Option:
    polynomial_key: tuple[tuple[Exponent, int], ...]
    correction_monomial: Exponent | None
    correction_coefficient: int

    @property
    def polynomial(self) -> Polynomial:
        return dict(self.polynomial_key)


def polynomial_key(poly: Polynomial) -> tuple[tuple[Exponent, int], ...]:
    return tuple(sorted(poly.items()))


def options_for_generator(
    generator_index: int,
    monomials: list[Exponent],
    coefficients: list[int],
    prime: int | None,
) -> list[Option]:
    base = {GENERATOR_EXPONENTS[generator_index]: 1}
    options = {
        polynomial_key(base): Option(polynomial_key(base), None, 0),
    }
    for monomial in monomials:
        for coefficient in coefficients:
            poly = add(base, {monomial: coefficient}, prime)
            key = polynomial_key(poly)
            options[key] = Option(
                key,
                monomial,
                reduce_coefficient(coefficient, prime),
            )
    return list(options.values())


def expected_commutator(
    left_index: int,
    right_index: int,
    prime: int | None,
) -> Polynomial:
    # Relations are tested in the displayed order [left, right].
    value = 0
    if left_index == 2 and right_index == 0:
        value = 1
    elif left_index == 3 and right_index == 1:
        value = 1
    elif left_index == 0 and right_index == 2:
        value = -1
    elif left_index == 1 and right_index == 3:
        value = -1
    return clean({ZERO_EXPONENT: value}, prime)


def compatible(
    left: Option,
    right: Option,
    left_index: int,
    right_index: int,
    prime: int | None,
) -> bool:
    return (
        commutator(left.polynomial, right.polynomial, prime)
        == expected_commutator(left_index, right_index, prime)
    )


def search_tuples(
    option_lists: list[list[Option]],
    prime: int | None,
) -> list[tuple[Option, Option, Option, Option]]:
    p1_options, p2_options, q1_options, q2_options = option_lists
    q1_by_p1 = {
        p1: [q1 for q1 in q1_options if compatible(q1, p1, 2, 0, prime)]
        for p1 in p1_options
    }
    q1_by_p2 = {
        p2: {q1 for q1 in q1_options if compatible(q1, p2, 2, 1, prime)}
        for p2 in p2_options
    }
    q2_by_p1 = {
        p1: [q2 for q2 in q2_options if compatible(q2, p1, 3, 0, prime)]
        for p1 in p1_options
    }
    q2_by_p2 = {
        p2: {q2 for q2 in q2_options if compatible(q2, p2, 3, 1, prime)}
        for p2 in p2_options
    }
    q2_commuting = {
        q1: {q2 for q2 in q2_options if compatible(q1, q2, 2, 3, prime)}
        for q1 in q1_options
    }

    answer = []
    for p1 in p1_options:
        for p2 in p2_options:
            if not compatible(p1, p2, 0, 1, prime):
                continue
            q1_candidates = q1_by_p2[p2]
            q2_p2_candidates = q2_by_p2[p2]
            for q1 in q1_by_p1[p1]:
                if q1 not in q1_candidates:
                    continue
                q2_candidates = q2_p2_candidates & q2_commuting[q1]
                for q2 in q2_by_p1[p1]:
                    if q2 in q2_candidates:
                        answer.append((p1, p2, q1, q2))
    return answer


class SparseSpan:
    def __init__(self, prime: int | None):
        self.prime = prime
        self.rows: dict[Exponent, dict[Exponent, Fraction | int]] = {}

    def coerce(self, value: int) -> Fraction | int:
        if self.prime is None:
            return Fraction(value)
        return value % self.prime

    def inverse(self, value: Fraction | int) -> Fraction | int:
        if self.prime is None:
            return 1 / value
        return pow(int(value), -1, self.prime)

    def normalize_value(self, value: Fraction | int) -> Fraction | int:
        if self.prime is None:
            return value
        return int(value) % self.prime

    def vector(self, poly: Polynomial) -> dict[Exponent, Fraction | int]:
        return {
            monomial: self.coerce(coefficient)
            for monomial, coefficient in poly.items()
            if self.coerce(coefficient)
        }

    def reduce(
        self,
        vector: dict[Exponent, Fraction | int],
    ) -> dict[Exponent, Fraction | int]:
        answer = dict(vector)
        while answer:
            pivot = max(answer, key=lambda item: (sum(item), item))
            if pivot not in self.rows:
                break
            factor = answer[pivot]
            for monomial, coefficient in self.rows[pivot].items():
                value = self.normalize_value(
                    answer.get(monomial, 0) - factor * coefficient
                )
                if value:
                    answer[monomial] = value
                else:
                    answer.pop(monomial, None)
        return answer

    def add(self, poly: Polynomial) -> bool:
        row = self.reduce(self.vector(poly))
        if not row:
            return False
        pivot = max(row, key=lambda item: (sum(item), item))
        inverse = self.inverse(row[pivot])
        row = {
            monomial: self.normalize_value(coefficient * inverse)
            for monomial, coefficient in row.items()
        }
        self.rows[pivot] = row
        return True

    def contains(self, poly: Polynomial) -> bool:
        return not self.reduce(self.vector(poly))


def recovery_degree(
    images: tuple[Polynomial, Polynomial, Polynomial, Polynomial],
    bound: int,
    prime: int | None,
) -> tuple[int | None, int]:
    powers = []
    for image in images:
        image_powers = [{ZERO_EXPONENT: 1}]
        for _ in range(bound):
            image_powers.append(multiply(image_powers[-1], image, prime))
        powers.append(image_powers)

    span = SparseSpan(prime)
    targets = [{exponent: 1} for exponent in GENERATOR_EXPONENTS]
    for degree in range(bound + 1):
        for exponents in compositions(degree, 4):
            word = {ZERO_EXPONENT: 1}
            for index, exponent in enumerate(exponents):
                word = multiply(word, powers[index][exponent], prime)
            span.add(word)
        if all(span.contains(target) for target in targets):
            return degree, len(span.rows)
    return None, len(span.rows)


def commutative_multiply(
    left: Polynomial,
    right: Polynomial,
    prime: int,
) -> Polynomial:
    answer: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                a + b for a, b in zip(left_monomial, right_monomial)
            )
            answer[monomial] = (
                answer.get(monomial, 0)
                + left_coefficient * right_coefficient
            )
    return clean(answer, prime)


def center_recovery_degree(
    center_images: tuple[Polynomial, Polynomial, Polynomial, Polynomial],
    bound: int,
    prime: int,
) -> tuple[int | None, int]:
    powers = []
    for image in center_images:
        image_powers = [{ZERO_EXPONENT: 1}]
        for _ in range(bound):
            image_powers.append(
                commutative_multiply(image_powers[-1], image, prime)
            )
        powers.append(image_powers)

    span = SparseSpan(prime)
    targets = [{exponent: 1} for exponent in GENERATOR_EXPONENTS]
    for degree in range(bound + 1):
        for exponents in compositions(degree, 4):
            word = {ZERO_EXPONENT: 1}
            for index, exponent in enumerate(exponents):
                word = commutative_multiply(
                    word,
                    powers[index][exponent],
                    prime,
                )
            span.add(word)
        if all(span.contains(target) for target in targets):
            return degree, len(span.rows)
    return None, len(span.rows)


def least_residue_integer_lift_satisfies_relations(
    images: tuple[Polynomial, Polynomial, Polynomial, Polynomial],
) -> bool:
    relations = (
        (0, 1),
        (2, 3),
        (2, 0),
        (2, 1),
        (3, 0),
        (3, 1),
    )
    return all(
        commutator(images[left], images[right], None)
        == expected_commutator(left, right, None)
        for left, right in relations
    )


def center_record(
    images: tuple[Polynomial, Polynomial, Polynomial, Polynomial],
    recovery_bound: int,
    prime: int | None,
) -> dict[str, object] | None:
    if prime is None:
        return None
    pth_powers = tuple(power(image, prime, prime) for image in images)
    central = all(
        all(all(exponent % prime == 0 for exponent in monomial) for monomial in poly)
        for poly in pth_powers
    )
    if not central:
        return {
            "p_power_images_are_ambient_central": False,
            "classification": "no_ambient_center_map",
        }
    center_images = tuple(
        {
            tuple(exponent // prime for exponent in monomial): coefficient
            for monomial, coefficient in poly.items()
        }
        for poly in pth_powers
    )
    recovered_at, span_dimension = center_recovery_degree(
        center_images,
        recovery_bound,
        prime,
    )
    dependency_sets = []
    univariate_degrees = []
    for poly in center_images:
        dependencies = {
            index
            for monomial in poly
            for index, exponent in enumerate(monomial)
            if exponent
        }
        dependency_sets.append(sorted(dependencies))
        if len(dependencies) == 1:
            variable = next(iter(dependencies))
            univariate_degrees.append(
                max(monomial[variable] for monomial in poly)
            )
        else:
            univariate_degrees.append(None)
    assigned_variables = [
        dependencies[0]
        for dependencies in dependency_sets
        if len(dependencies) == 1
    ]
    separable_univariate = (
        len(assigned_variables) == 4
        and sorted(assigned_variables) == list(range(4))
    )
    proper_separable_degree = (
        separable_univariate
        and any(
            degree is not None and degree > 1
            for degree in univariate_degrees
        )
    )
    if proper_separable_degree:
        classification = "proper_center_subalgebra_separable_degree"
    elif recovered_at is not None:
        classification = "generates_center_within_bound"
    else:
        classification = "unresolved_bounded_center_generation"
    return {
        "p_power_images_are_ambient_central": True,
        "classification": classification,
        "proper_subalgebra_certificate": (
            "The center map is a tensor product of four nonconstant "
            "univariate maps on distinct variables, and its fraction-field "
            "degree is the product of their degrees, which is greater than one."
            if proper_separable_degree
            else None
        ),
        "separable_univariate_degrees": (
            univariate_degrees if separable_univariate else None
        ),
        "ambient_center_generators_recovered_at_word_degree": recovered_at,
        "center_word_span_dimension_at_stop": span_dimension,
        "center_images": [
            [
                {
                    "coefficient": coefficient,
                    "exponent": list(monomial),
                }
                for monomial, coefficient in sorted(poly.items())
            ]
            for poly in center_images
        ],
    }


def option_record(option: Option) -> dict[str, object]:
    return {
        "correction_coefficient": option.correction_coefficient,
        "correction_monomial": (
            list(option.correction_monomial)
            if option.correction_monomial is not None
            else None
        ),
    }


def top_symbol_exponent(option: Option) -> Exponent:
    terms = option.polynomial
    return max(terms, key=lambda item: (sum(item), item))


def determinant(matrix: list[list[int]]) -> int:
    if len(matrix) == 1:
        return matrix[0][0]
    return sum(
        (-1) ** column
        * matrix[0][column]
        * determinant(
            [
                row[:column] + row[column + 1 :]
                for row in matrix[1:]
            ]
        )
        for column in range(len(matrix))
    )


def tuple_record(
    candidate: tuple[Option, Option, Option, Option],
    recovery_bound: int,
    prime: int | None,
) -> dict[str, object]:
    images = tuple(option.polynomial for option in candidate)
    recovered_at, span_dimension = recovery_degree(
        images,
        recovery_bound,
        prime,
    )
    exponents = [list(top_symbol_exponent(option)) for option in candidate]
    top_determinant = determinant(exponents)
    mod_p_center = center_record(images, recovery_bound, prime)
    proper_via_center = (
        mod_p_center is not None
        and mod_p_center["classification"]
        == "proper_center_subalgebra_separable_degree"
    )
    return {
        "images": {
            name: option_record(option)
            for name, option in zip(GENERATOR_NAMES, candidate)
        },
        "least_residue_integer_lift_satisfies_relations": (
            least_residue_integer_lift_satisfies_relations(images)
        ),
        "classification": (
            "proper_A2_subalgebra_via_center"
            if proper_via_center
            else (
                "generates_A2_within_bound"
                if recovered_at is not None
                else "unresolved_bounded_generation"
            )
        ),
        "ambient_generators_recovered_at_word_degree": recovered_at,
        "word_span_dimension_at_stop": span_dimension,
        "top_symbol_exponent_matrix": exponents,
        "top_symbol_jacobian_monomial_coefficient": top_determinant,
        "top_symbol_is_keller": (
            top_determinant in (-1, 1)
            and sum(sum(row) for row in exponents) == 4
        ),
        "mod_p_center": mod_p_center,
    }


def parse_coefficients(text: str, prime: int | None) -> list[int]:
    if text == "all":
        if prime is None:
            raise ValueError("'all' coefficients require a prime field")
        return list(range(1, prime))
    values = [int(value) for value in text.split(",") if value]
    if not values:
        raise ValueError("the coefficient alphabet is empty")
    if prime is not None:
        values = sorted({value % prime for value in values if value % prime})
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=3)
    parser.add_argument(
        "--coefficients",
        default="-1,1",
        help="comma-separated nonzero integers, or 'all' in prime characteristic",
    )
    parser.add_argument(
        "--prime",
        type=int,
        default=0,
        help="0 for characteristic zero; otherwise a prime modulus",
    )
    parser.add_argument("--recovery-degree", type=int, default=6)
    parser.add_argument(
        "--max-records",
        type=int,
        default=10000,
        help="refuse to truncate silently if more candidates are found",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.max_degree < 2:
        raise ValueError("--max-degree must be at least 2")
    if args.recovery_degree < 0:
        raise ValueError("--recovery-degree must be nonnegative")
    prime = args.prime or None
    coefficients = parse_coefficients(args.coefficients, prime)
    monomials = monomials_between(2, args.max_degree)
    option_lists = [
        options_for_generator(index, monomials, coefficients, prime)
        for index in range(4)
    ]
    candidates = search_tuples(option_lists, prime)
    nonidentity = [
        candidate
        for candidate in candidates
        if any(option.correction_monomial is not None for option in candidate)
    ]
    if len(nonidentity) > args.max_records:
        raise RuntimeError(
            f"found {len(nonidentity)} nonidentity tuples, exceeding "
            f"--max-records={args.max_records}; no truncated artifact was written"
        )

    records = [
        tuple_record(candidate, args.recovery_degree, prime)
        for candidate in nonidentity
    ]
    classifications: dict[str, int] = {}
    for record in records:
        label = str(record["classification"])
        classifications[label] = classifications.get(label, 0) + 1
    integer_lifts = sum(
        bool(record["least_residue_integer_lift_satisfies_relations"])
        for record in records
    )
    non_keller_top_symbols = sum(
        not bool(record["top_symbol_is_keller"]) for record in records
    )

    artifact = {
        "format": "sparse-a2-weyl-endomorphism-search-v1",
        "software": {
            "python_implementation": platform.python_implementation(),
            "python_version": sys.version.split()[0],
            "dependencies": "Python standard library only",
        },
        "scope": (
            "Each of x1,x2,d1,d2 is changed by zero or one PBW monomial; "
            "this is exhaustive only for the displayed degree and coefficient "
            "alphabet. Bounded failure to recover generators is not a "
            "proper-subalgebra certificate."
        ),
        "convention": "[d_i,x_j]=delta_ij; PBW order x1,x2,d1,d2",
        "parameters": {
            "characteristic": 0 if prime is None else prime,
            "max_correction_bernstein_degree": args.max_degree,
            "coefficient_alphabet": coefficients,
            "recovery_word_degree": args.recovery_degree,
            "correction_monomial_count": len(monomials),
            "options_per_generator": [len(options) for options in option_lists],
        },
        "summary": {
            "weyl_tuples_including_identity": len(candidates),
            "nonidentity_weyl_tuples": len(nonidentity),
            "least_residue_integer_lifts_satisfying_relations": integer_lifts,
            "modular_tuples_failing_that_literal_lift": (
                len(records) - integer_lifts
            ),
            "non_keller_top_symbol_tuples": non_keller_top_symbols,
            "classifications": classifications,
        },
        "candidates": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n")

    field = "Q" if prime is None else f"F_{prime}"
    print(
        f"PASS: exhausted {len(monomials)} correction monomials and "
        f"{len(coefficients)} coefficients over {field}"
    )
    print(
        f"PASS: found {len(nonidentity)} nonidentity exact Weyl tuples "
        f"({integer_lifts} integer lifts)"
    )
    print(f"RESULT: bounded generation classifications {classifications}")
    print(f"RESULT: {non_keller_top_symbols} tuples have non-Keller top symbols")
    try:
        displayed_output = args.output.resolve().relative_to(ROOT)
    except ValueError:
        displayed_output = args.output.resolve()
    print(f"PASS: wrote {displayed_output}")


if __name__ == "__main__":
    main()
