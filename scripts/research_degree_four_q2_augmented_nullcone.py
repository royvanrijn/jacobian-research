#!/usr/bin/env python3
"""Local synchronization test for the quartic q2-augmented moment origin.

This is a branchwise exact calculation, not a global nullcone theorem.
On the branch where the Sym^2 Casimir component is nonzero and q2=0,
SL2 normalizes that component to the highest-weight vector E.  Write the
higher components in the integral weight bases

    ad(F)^k(E^r),  r=2,3,4,  0<=k<=2r.

The expected synchronized nullcone retains k<r.  This leaves nine
allowed coordinates and twelve forbidden normal coordinates.

At one deterministic integral point of the allowed space, the script:

* constructs the Taylor jets of mu_2,...,mu_21 in the twelve normal
  coordinates modulo a good prime;
* uses the four linear pivots in mu_2,...,mu_5 to eliminate four normal
  coordinates as formal series;
* computes local standard bases for the remaining eight-variable jets.

The output is exact finite-field evidence about the completed moment
origin near a generic synchronized point.  It neither treats every point
of this branch nor the boundary Sym^2=0.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
from collections import defaultdict
from math import factorial
from pathlib import Path

import sympy as sp

from research_completed_moment_algebra import sl2_matrices


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_q2_augmented_nullcone_local.json"
)
Matrix = sp.Matrix
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, int]
PositionPolynomial = dict[tuple[int, int], Polynomial]


def add_polynomials(
    left: Polynomial,
    right: Polynomial,
    prime: int,
) -> Polynomial:
    result: defaultdict[Exponent, int] = defaultdict(int)
    result.update(left)
    for monomial, coefficient in right.items():
        result[monomial] = (
            result[monomial] + coefficient
        ) % prime
    return {
        monomial: coefficient
        for monomial, coefficient in result.items()
        if coefficient
    }


def multiply_polynomials(
    left: Polynomial,
    right: Polynomial,
    prime: int,
    max_degree: int,
) -> Polynomial:
    result: defaultdict[Exponent, int] = defaultdict(int)
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_value + right_value
                for left_value, right_value in zip(
                    left_monomial,
                    right_monomial,
                    strict=True,
                )
            )
            if sum(monomial) > max_degree:
                continue
            result[monomial] = (
                result[monomial]
                + left_coefficient * right_coefficient
            ) % prime
    return {
        monomial: coefficient
        for monomial, coefficient in result.items()
        if coefficient
    }


def multiply_numeric_positions(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
    prime: int,
) -> dict[tuple[int, int], int]:
    result: defaultdict[tuple[int, int], int] = defaultdict(int)
    for (left_i, left_j), left_coefficient in left.items():
        for (right_i, right_j), right_coefficient in right.items():
            position = (left_i + right_i, left_j + right_j)
            result[position] = (
                result[position]
                + left_coefficient * right_coefficient
            ) % prime
    return {
        position: coefficient
        for position, coefficient in result.items()
        if coefficient
    }


def primitive_matrix(matrix: Matrix) -> Matrix:
    entries = [abs(int(value)) for value in matrix if value]
    divisor = int(sp.gcd_list(entries)) if entries else 1
    return matrix.applyfunc(lambda value: sp.Rational(value, divisor))


def component_weight_bases() -> tuple[
    Matrix,
    dict[tuple[int, int], Matrix],
    list[tuple[int, int]],
    list[tuple[int, int]],
]:
    raising, lowering, _cartan = sl2_matrices(4)

    def adjoint(matrix: Matrix) -> Matrix:
        return lowering * matrix - matrix * lowering

    basis = {}
    allowed = []
    forbidden = []
    for component in range(2, 5):
        vector = raising**component
        for lowering_order in range(2 * component + 1):
            basis[component, lowering_order] = primitive_matrix(vector)
            target = allowed if lowering_order < component else forbidden
            target.append((component, lowering_order))
            vector = adjoint(vector)
    assert len(allowed) == 9
    assert len(forbidden) == 12
    return raising, basis, allowed, forbidden


def modular_rational(value: sp.Expr, prime: int) -> int:
    numerator = int(sp.numer(value)) % prime
    denominator = int(sp.denom(value)) % prime
    return numerator * pow(denominator, -1, prime) % prime


def build_moment_jets(
    prime: int,
    max_jet: int,
) -> tuple[
    dict[int, Polynomial],
    list[tuple[int, int]],
    list[tuple[int, int]],
]:
    raising, basis, allowed, forbidden = component_weight_bases()
    operator = raising.copy()
    for index, coordinate in enumerate(allowed):
        operator += (index + 2) * basis[coordinate]

    factorial_diagonal = sp.diag(24, 6, 4, 6, 24)
    coefficients = factorial_diagonal.inv() * operator.T
    normal_directions = [
        factorial_diagonal.inv() * basis[coordinate].T
        for coordinate in forbidden
    ]

    base = {
        (row, column): modular_rational(
            coefficients[row, column], prime
        )
        for row in range(5)
        for column in range(5)
        if coefficients[row, column]
    }
    normal_terms = []
    for row in range(5):
        for column in range(5):
            for variable, direction in enumerate(normal_directions):
                coefficient = modular_rational(
                    direction[row, column], prime
                )
                if coefficient:
                    normal_terms.append(
                        (row, column, variable, coefficient)
                    )

    base_powers = [{(0, 0): 1}]
    for _order in range(21):
        base_powers.append(
            multiply_numeric_positions(
                base_powers[-1], base, prime
            )
        )

    zero = (0,) * len(forbidden)
    normal_powers: list[PositionPolynomial] = [
        {(0, 0): {zero: 1}}
    ]
    for _degree in range(1, max_jet + 1):
        result: dict[
            tuple[int, int], defaultdict[Exponent, int]
        ] = {}
        for (
            (left_i, left_j),
            polynomial,
        ) in normal_powers[-1].items():
            for (
                right_i,
                right_j,
                variable,
                scalar,
            ) in normal_terms:
                target = result.setdefault(
                    (left_i + right_i, left_j + right_j),
                    defaultdict(int),
                )
                for monomial, coefficient in polynomial.items():
                    updated = list(monomial)
                    updated[variable] += 1
                    exponent = tuple(updated)
                    target[exponent] = (
                        target[exponent] + coefficient * scalar
                    ) % prime
        normal_powers.append(
            {
                position: {
                    monomial: coefficient
                    for monomial, coefficient in polynomial.items()
                    if coefficient
                }
                for position, polynomial in result.items()
            }
        )

    def moment_piece(order: int, degree: int) -> Polynomial:
        result: defaultdict[Exponent, int] = defaultdict(int)
        choose = math.comb(order, degree) % prime
        for (
            (normal_i, normal_j),
            polynomial,
        ) in normal_powers[degree].items():
            for (
                base_i,
                base_j,
            ), base_coefficient in base_powers[order - degree].items():
                final_i = normal_i + base_i
                final_j = normal_j + base_j
                if final_i != final_j:
                    continue
                scalar = (
                    choose
                    * base_coefficient
                    * factorial(final_i)
                    * factorial(4 * order - final_i)
                ) % prime
                for monomial, coefficient in polynomial.items():
                    result[monomial] = (
                        result[monomial] + scalar * coefficient
                    ) % prime
        return {
            monomial: coefficient
            for monomial, coefficient in result.items()
            if coefficient
        }

    moment_jets = {}
    for order in range(2, 22):
        moment: Polynomial = {}
        for degree in range(1, min(order, max_jet) + 1):
            moment = add_polynomials(
                moment,
                moment_piece(order, degree),
                prime,
            )
        moment_jets[order] = moment
    return moment_jets, allowed, forbidden


def invert_matrix_mod(
    matrix: list[list[int]],
    prime: int,
) -> list[list[int]]:
    size = len(matrix)
    augmented = [
        [
            *[value % prime for value in row],
            *[int(row_index == column) for column in range(size)],
        ]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row
            for row in range(column, size)
            if augmented[row][column]
        )
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        inverse = pow(augmented[column][column], -1, prime)
        augmented[column] = [
            value * inverse % prime
            for value in augmented[column]
        ]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            scalar = augmented[row][column]
            augmented[row] = [
                (left - scalar * right) % prime
                for left, right in zip(
                    augmented[row],
                    augmented[column],
                    strict=True,
                )
            ]
    return [row[size:] for row in augmented]


def formal_pivot_elimination(
    moment_jets: dict[int, Polynomial],
    prime: int,
    max_jet: int,
    compose_reduced_moments: bool,
) -> tuple[
    dict[int, Polynomial],
    list[int],
    list[int],
    list[Polynomial],
]:
    variable_count = 12
    unit_vectors = [
        tuple(int(index == variable) for index in range(variable_count))
        for variable in range(variable_count)
    ]
    linear_matrix = [
        [moment_jets[order].get(unit, 0) for unit in unit_vectors]
        for order in range(2, 6)
    ]

    work = [row[:] for row in linear_matrix]
    pivots = []
    row = 0
    for column in range(variable_count):
        pivot = next(
            (
                index
                for index in range(row, 4)
                if work[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [value * inverse % prime for value in work[row]]
        for index in range(4):
            if index == row or not work[index][column]:
                continue
            scalar = work[index][column]
            work[index] = [
                (left - scalar * right) % prime
                for left, right in zip(
                    work[index],
                    work[row],
                    strict=True,
                )
            ]
        pivots.append(column)
        row += 1
        if row == 4:
            break
    assert len(pivots) == 4
    free = [
        variable
        for variable in range(variable_count)
        if variable not in pivots
    ]

    pivot_matrix = [
        [linear_matrix[row][column] for column in pivots]
        for row in range(4)
    ]
    inverse_pivot_matrix = invert_matrix_mod(pivot_matrix, prime)

    zero = (0,) * len(free)
    substitutions: list[Polynomial] = []
    for variable in range(variable_count):
        if variable in free:
            exponent = [0] * len(free)
            exponent[free.index(variable)] = 1
            substitutions.append({tuple(exponent): 1})
        else:
            substitutions.append({})

    # Solve the linear pivot terms, including the nonpivot columns.
    for free_index, variable in enumerate(free):
        right_hand_side = [
            -linear_matrix[row][variable] % prime
            for row in range(4)
        ]
        solution = [
            sum(
                inverse_pivot_matrix[row][column]
                * right_hand_side[column]
                for column in range(4)
            )
            % prime
            for row in range(4)
        ]
        exponent = [0] * len(free)
        exponent[free_index] = 1
        monomial = tuple(exponent)
        for pivot_index, coefficient in enumerate(solution):
            if coefficient:
                substitutions[pivots[pivot_index]][monomial] = coefficient

    substitution_cache: dict[
        tuple[int, Exponent], Polynomial
    ] = {}

    def substitute_monomial(
        monomial: Exponent,
        cutoff: int,
    ) -> Polynomial:
        cache_key = (cutoff, monomial)
        cached = substitution_cache.get(cache_key)
        if cached is not None:
            return cached
        term = {zero: 1}
        for variable, exponent in enumerate(monomial):
            for _power in range(exponent):
                term = multiply_polynomials(
                    term,
                    substitutions[variable],
                    prime,
                    cutoff,
                )
                if not term:
                    break
            if not term:
                break
        substitution_cache[cache_key] = term
        return term

    def substitute(
        polynomial: Polynomial,
        cutoff: int,
    ) -> Polynomial:
        result: Polynomial = {}
        for monomial, coefficient in polynomial.items():
            term = {
                exponent: coefficient * value % prime
                for exponent, value in substitute_monomial(
                    monomial, cutoff
                ).items()
            }
            result = add_polynomials(result, term, prime)
        return result

    for degree in range(2, max_jet + 1):
        residuals = []
        for order in range(2, 6):
            residuals.append(
                {
                    monomial: coefficient
                    for monomial, coefficient in substitute(
                        moment_jets[order], degree
                    ).items()
                    if sum(monomial) == degree
                }
            )
        monomials = sorted(
            set().union(
                *(residual.keys() for residual in residuals)
            )
        )
        corrections = [defaultdict(int) for _ in range(4)]
        for monomial in monomials:
            right_hand_side = [
                -residual.get(monomial, 0) % prime
                for residual in residuals
            ]
            solution = [
                sum(
                    inverse_pivot_matrix[row][column]
                    * right_hand_side[column]
                    for column in range(4)
                )
                % prime
                for row in range(4)
            ]
            for row, coefficient in enumerate(solution):
                if coefficient:
                    corrections[row][monomial] = coefficient
        for row, pivot in enumerate(pivots):
            substitutions[pivot] = add_polynomials(
                substitutions[pivot],
                corrections[row],
                prime,
            )
        substitution_cache.clear()

        for order in range(2, 6):
            remainder = {
                monomial: coefficient
                for monomial, coefficient in substitute(
                    moment_jets[order], degree
                ).items()
                if sum(monomial) <= degree
            }
            assert not remainder

    reduced = (
        {
            order: substitute(moment_jets[order], max_jet)
            for order in range(6, 22)
        }
        if compose_reduced_moments
        else {}
    )
    return reduced, pivots, free, substitutions


def polynomial_string(
    polynomial: Polynomial,
    variables: list[str],
) -> str:
    terms = []
    for monomial, coefficient in polynomial.items():
        factors = [] if coefficient == 1 else [str(coefficient)]
        for variable, exponent in zip(
            variables, monomial, strict=True
        ):
            if exponent == 1:
                factors.append(variable)
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        terms.append("*".join(factors) if factors else "1")
    return "+".join(terms) if terms else "0"


def local_standard_basis(
    reduced: dict[int, Polynomial],
    prime: int,
    jet: int,
    timeout: int,
    ordering: str,
    basis_algorithm: str,
) -> dict[str, int | str]:
    executable = shutil.which("Singular")
    if executable is None:
        raise SystemExit("Singular is required on PATH")
    variables = [f"x{index}" for index in range(8)]
    polynomials = [
        {
            monomial: coefficient
            for monomial, coefficient in polynomial.items()
            if sum(monomial) <= jet
        }
        for polynomial in reduced.values()
    ]
    code = (
        f"ring r={prime},({','.join(variables)}),{ordering};\n"
        + "ideal I="
        + ",".join(
            polynomial_string(polynomial, variables)
            for polynomial in polynomials
        )
        + ";\n"
        + f"ideal G={basis_algorithm}(I);\n"
        + "dim(G);\n"
        + "vdim(G);\n"
        + "size(G);\n"
    )
    try:
        completed = subprocess.run(
            [executable, "-q"],
            input=code,
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "timeout_seconds": timeout}
    output = [
        line.strip()
        for line in completed.stdout.splitlines()
        if line.strip()
    ]
    if completed.stderr.strip():
        raise AssertionError(completed.stderr)
    return {
        "status": "completed",
        "dimension": int(output[0]),
        "vector_space_dimension": int(output[1]),
        "standard_basis_size": int(output[2]),
    }


def native_local_standard_basis(
    moment_jets: dict[int, Polynomial],
    substitutions: list[Polynomial],
    prime: int,
    jet: int,
    timeout: int,
    ordering: str,
    basis_algorithm: str,
) -> dict[str, int | str]:
    executable = shutil.which("Singular")
    if executable is None:
        raise SystemExit("Singular is required on PATH")
    source_variables = [f"y{index}" for index in range(12)]
    target_variables = [f"x{index}" for index in range(8)]
    source_ideal = ",".join(
        polynomial_string(
            {
                monomial: coefficient
                for monomial, coefficient in moment_jets[order].items()
                if sum(monomial) <= jet
            },
            source_variables,
        )
        for order in range(6, 22)
    )
    map_entries = ",".join(
        polynomial_string(
            {
                monomial: coefficient
                for monomial, coefficient in polynomial.items()
                if sum(monomial) <= jet
            },
            target_variables,
        )
        for polynomial in substitutions
    )
    code = (
        f"ring ra={prime},({','.join(source_variables)}),dp;\n"
        + f"ideal IA={source_ideal};\n"
        + f"ring rb={prime},({','.join(target_variables)}),"
        f"{ordering};\n"
        + f"map phi=ra,{map_entries};\n"
        + f"ideal I=jet(phi(IA),{jet});\n"
        + f"ideal G={basis_algorithm}(I);\n"
        + "dim(G);\n"
        + "vdim(G);\n"
        + "size(G);\n"
    )
    try:
        completed = subprocess.run(
            [executable, "-q"],
            input=code,
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "timeout_seconds": timeout}
    output = [
        line.strip()
        for line in completed.stdout.splitlines()
        if line.strip()
    ]
    if completed.stderr.strip():
        raise AssertionError(completed.stderr)
    return {
        "status": "completed",
        "dimension": int(output[0]),
        "vector_space_dimension": int(output[1]),
        "standard_basis_size": int(output[2]),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--max-jet", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--ordering",
        choices=("dp", "ds"),
        default="dp",
        help="global dp is faster; local ds directly studies the germ",
    )
    parser.add_argument(
        "--basis-algorithm",
        choices=("std", "slimgb"),
        default="slimgb",
    )
    parser.add_argument(
        "--composition",
        choices=("native", "python"),
        default="native",
        help="compose pivot series using Singular maps or Python",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    moment_jets, allowed, forbidden = build_moment_jets(
        arguments.prime,
        arguments.max_jet,
    )
    reduced, pivots, free, substitutions = formal_pivot_elimination(
        moment_jets,
        arguments.prime,
        arguments.max_jet,
        arguments.composition == "python",
    )
    standard_bases = {}
    for jet in range(2, arguments.max_jet + 1):
        if arguments.composition == "native":
            result = native_local_standard_basis(
                moment_jets,
                substitutions,
                arguments.prime,
                jet,
                arguments.timeout,
                arguments.ordering,
                arguments.basis_algorithm,
            )
        else:
            result = local_standard_basis(
                reduced,
                arguments.prime,
                jet,
                arguments.timeout,
                arguments.ordering,
                arguments.basis_algorithm,
            )
        standard_bases[str(jet)] = result
        print(f"jet {jet}: {result}", flush=True)

    payload = {
        "format": "degree-four-q2-augmented-nullcone-local-v1",
        "scope": (
            "one formal normal slice at a deterministic synchronized "
            "point on the nonzero Sym^2 nullcone branch"
        ),
        "prime": arguments.prime,
        "maximum_normal_jet_degree": arguments.max_jet,
        "standard_basis_ordering": arguments.ordering,
        "standard_basis_algorithm": arguments.basis_algorithm,
        "formal_series_composition": arguments.composition,
        "normalized_sym2_component": "highest-weight raising matrix E",
        "allowed_point_coefficients": list(range(2, 11)),
        "allowed_weight_coordinates": [list(value) for value in allowed],
        "forbidden_weight_coordinates": [
            list(value) for value in forbidden
        ],
        "linear_pivot_indices": pivots,
        "linear_pivot_coordinates": [
            list(forbidden[index]) for index in pivots
        ],
        "free_normal_indices": free,
        "free_normal_coordinates": [
            list(forbidden[index]) for index in free
        ],
        "pivot_series_term_counts": [
            len(substitutions[index]) for index in pivots
        ],
        "reduced_moment_term_counts": (
            {
                str(order): len(polynomial)
                for order, polynomial in reduced.items()
            }
            if reduced
            else None
        ),
        "local_standard_bases": standard_bases,
        "scope_warning": (
            "a completed finite jet can prove only the stated local "
            "normal-slice result when its Nakayama containment is also "
            "certified; no global q2-augmented nullcone equality follows"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
