#!/usr/bin/env python3
"""Completed-invariant zero-fiber probes for bidegree (3,3).

This research checker joins two previously separate calculations:

* the Hilbert/Jacobian comparison between corrected moment parameters and
  systems augmented by the quadratic Casimirs q_2 and q_4;
* modular zero-fiber probes after q_2=0, on the normalized null-quadratic
  chart F_2=X^2 and on its boundary F_2=0.

The Hilbert coefficients and nonzero modular Jacobian ranks are exact
characteristic-zero certificates.  The Singular radical probes are only
finite-field evidence unless independently reconstructed over QQ.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import sympy as sp

from research_completed_moment_algebra import (
    casimir_projectors_mod,
    component_quadratic_exact_sparse,
    component_matrices_mod,
    deterministic_point,
    hilbert_numerator,
    invariant_values_mod,
    invariant_hilbert_coefficients,
    moments_mod,
    moment_jacobian_mod,
    monomial_value,
    rank_mod,
    weighted_exponents,
)
from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    Q_POLYNOMIALS,
    WEIGHTS,
)
from verify_two_pair_sic_bidegree33_anchor_jacobians import (
    COEFFICIENT_MAP,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_casimir_fiber.json"
)
DEFAULT_PRIME = 32003
QUOTIENT_DIMENSION = 13
AMBIENT_DIMENSION = 16
MOMENT_CUTOFF = 24
HILBERT_CUTOFF = 120

# The highest-weight quadratic F_2=X^2 has phase +1.  After removing the
# unmatched biform variable, its coefficient polynomial is (1+q)^2.
NULL_QUADRATIC_WEIGHT = 1
NULL_QUADRATIC_Q = (1, 2, 1)

# On F_2=X^2, synchronization with the same destabilizing root is the
# coordinate space s3=...=s6=t2=t3=t4=0.
NULL_QUADRATIC_FORBIDDEN = (3, 4, 5, 6, 9, 10, 11)
NULL_QUADRATIC_ALLOWED = tuple(
    index
    for index in range(len(PARAMETERS))
    if index not in NULL_QUADRATIC_FORBIDDEN
)
assert NULL_QUADRATIC_ALLOWED == (0, 1, 2, 7, 8)


def convolve_mod(
    left: tuple[int, ...],
    right: tuple[int, ...],
    prime: int,
) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        if left_coefficient % prime == 0:
            continue
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] = (
                answer[left_index + right_index]
                + left_coefficient * right_coefficient
            ) % prime
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_powers_mod(
    polynomial: tuple[int, ...],
    maximum: int,
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    powers = [(1,)]
    for _ in range(maximum):
        powers.append(convolve_mod(powers[-1], polynomial, prime))
    return tuple(powers)


def restricted_moment_terms_mod(
    order: int,
    quadratic_mode: str,
    prime: int,
) -> dict[tuple[int, ...], int]:
    """Return one restricted moment as a sparse polynomial modulo prime."""

    assert quadratic_mode in {"null", "zero"}
    factorials = [
        factorial(index) % prime for index in range(3 * order + 1)
    ]
    inverse_factorials = [
        pow(factorial(index) % prime, -1, prime)
        for index in range(order + 1)
    ]
    basis_powers = tuple(
        polynomial_powers_mod(polynomial, order, prime)
        for polynomial in Q_POLYNOMIALS
    )
    quadratic_powers = polynomial_powers_mod(
        NULL_QUADRATIC_Q, order, prime
    )
    parameter_order = (0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9)
    exponents = [0] * len(PARAMETERS)
    answer: dict[tuple[int, ...], int] = defaultdict(int)
    order_factorial = factorials[order]

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(
        position: int,
        degree_left: int,
    ) -> tuple[int, int]:
        remaining = [
            WEIGHTS[index] for index in parameter_order[position:]
        ]
        if quadratic_mode == "null":
            remaining.append(NULL_QUADRATIC_WEIGHT)
        if degree_left == 0:
            return 0, 0
        if not remaining:
            return 1, 0
        return (
            degree_left * min(remaining),
            degree_left * max(remaining),
        )

    def visit(
        position: int,
        used_degree: int,
        weight: int,
        shift: int,
        inverse_denominator: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(parameter_order):
            quadratic_exponent = order - used_degree
            if quadratic_mode == "zero" and quadratic_exponent:
                return
            total_weight = weight
            product = q_polynomial
            total_shift = shift
            if quadratic_mode == "null":
                total_weight += (
                    NULL_QUADRATIC_WEIGHT * quadratic_exponent
                )
                total_shift += (
                    NULL_QUADRATIC_WEIGHT * quadratic_exponent
                )
                product = convolve_mod(
                    product,
                    quadratic_powers[quadratic_exponent],
                    prime,
                )
            if total_weight != 0:
                return

            scalar = (
                order_factorial
                * inverse_denominator
                * inverse_factorials[quadratic_exponent]
            ) % prime
            contraction = 0
            for q_degree, coefficient in enumerate(product):
                diagonal = total_shift + q_degree
                if 0 <= diagonal <= 3 * order:
                    contraction += (
                        coefficient
                        * factorials[3 * order - diagonal]
                        * factorials[diagonal]
                    )
            coefficient = scalar * contraction % prime
            if coefficient:
                exponent_tuple = tuple(exponents)
                answer[exponent_tuple] = (
                    answer[exponent_tuple] + coefficient
                ) % prime
            return

        parameter_index = parameter_order[position]
        parameter_weight = WEIGHTS[parameter_index]
        available = order - used_degree
        for exponent in range(available + 1):
            new_weight = weight + exponent * parameter_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(
                position + 1, degree_left
            )
            target_weight = -new_weight
            if not minimum <= target_weight <= maximum:
                continue
            exponents[parameter_index] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(parameter_weight, 0) * exponent,
                inverse_denominator
                * inverse_factorials[exponent]
                % prime,
                convolve_mod(
                    q_polynomial,
                    basis_powers[parameter_index][exponent],
                    prime,
                ),
            )
        exponents[parameter_index] = 0

    visit(0, 0, 0, 0, 1, (1,))
    return {
        exponent_tuple: coefficient % prime
        for exponent_tuple, coefficient in answer.items()
        if coefficient % prime
    }


def serialize_modular_polynomial(
    terms: dict[tuple[int, ...], int],
    prime: int,
) -> str:
    pieces = []
    for exponents, coefficient in sorted(terms.items()):
        signed = coefficient if coefficient <= prime // 2 else coefficient - prime
        factors = []
        for variable, exponent in zip(PARAMETERS, exponents, strict=True):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if monomial:
            if signed == 1:
                pieces.append(monomial)
            elif signed == -1:
                pieces.append(f"-{monomial}")
            else:
                pieces.append(f"{signed}*{monomial}")
        else:
            pieces.append(str(signed))
    expression = "+".join(pieces).replace("+-", "-")
    return expression or "0"


def primitive_modular_terms(
    terms: dict[tuple[int, ...], int],
    prime: int,
) -> dict[tuple[int, ...], int]:
    """Scale a nonzero modular polynomial to leading coefficient one."""

    if not terms:
        return terms
    leading = terms[min(terms)]
    inverse = pow(leading, -1, prime)
    return {
        exponents: coefficient * inverse % prime
        for exponents, coefficient in terms.items()
    }


def evaluate_sparse_mod(
    terms: dict[tuple[int, ...], int],
    values: tuple[int, ...],
    prime: int,
) -> int:
    answer = 0
    for exponents, coefficient in terms.items():
        answer += coefficient * monomial_value(
            list(values), exponents, prime
        )
    return answer % prime


def serialize_sparse_mod(
    terms: dict[tuple[int, ...], int],
    variables: tuple[str, ...],
    prime: int,
) -> str:
    pieces = []
    for exponents, coefficient in sorted(terms.items()):
        signed = (
            coefficient
            if coefficient <= prime // 2
            else coefficient - prime
        )
        factors = []
        for variable, exponent in zip(
            variables, exponents, strict=True
        ):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if monomial:
            if signed == 1:
                pieces.append(monomial)
            elif signed == -1:
                pieces.append(f"-{monomial}")
            else:
                pieces.append(f"{signed}*{monomial}")
        else:
            pieces.append(str(signed))
    return "+".join(pieces).replace("+-", "-") or "0"


def null_quadratic_normal_symbols(
    prime: int,
    maximum_order: int = 12,
) -> tuple[dict[str, object], list[list[dict[tuple[int, ...], int]]]]:
    """Construct the J_sync-linear symbols on F_2=X^2."""

    orders = tuple(range(2, maximum_order + 1))
    matrix: list[list[dict[tuple[int, ...], int]]] = []
    order_records = []
    for order in orders:
        terms = restricted_moment_terms_mod(order, "null", prime)
        normal_degrees = [
            sum(exponents[index] for index in NULL_QUADRATIC_FORBIDDEN)
            for exponents in terms
        ]
        assert normal_degrees and min(normal_degrees) >= 1
        row = []
        for forbidden_index in NULL_QUADRATIC_FORBIDDEN:
            coefficient_terms: dict[tuple[int, ...], int] = {}
            for exponents, coefficient in terms.items():
                if (
                    sum(
                        exponents[index]
                        for index in NULL_QUADRATIC_FORBIDDEN
                    )
                    != 1
                    or exponents[forbidden_index] != 1
                ):
                    continue
                allowed_exponents = tuple(
                    exponents[index]
                    for index in NULL_QUADRATIC_ALLOWED
                )
                coefficient_terms[allowed_exponents] = (
                    coefficient_terms.get(allowed_exponents, 0)
                    + coefficient
                ) % prime
            row.append({
                exponents: coefficient
                for exponents, coefficient in coefficient_terms.items()
                if coefficient
            })
        matrix.append(row)
        order_records.append({
            "order": order,
            "moment_terms": len(terms),
            "minimum_sync_normal_degree": min(normal_degrees),
            "linear_sync_normal_terms": sum(
                len(entry) for entry in row
            ),
            "nonzero_linear_normal_coordinates": [
                PARAMETERS[index]
                for index, entry in zip(
                    NULL_QUADRATIC_FORBIDDEN, row, strict=True
                )
                if entry
            ],
        })

    sample_records = []
    maximum_rank = 0
    rank_certificate = None
    target_linear_rank = 3
    for sample in range(1, 13):
        values = tuple(
            (
                sample * (17 + 2 * position)
                + (position + 1) * (position + 3)
            )
            % prime
            for position in range(len(NULL_QUADRATIC_ALLOWED))
        )
        evaluated = [
            [
                evaluate_sparse_mod(entry, values, prime)
                for entry in row
            ]
            for row in matrix
        ]
        rank = rank_mod(evaluated, prime)
        maximum_rank = max(maximum_rank, rank)
        sample_records.append({
            "sample": sample,
            "allowed_values": {
                PARAMETERS[index]: value
                for index, value in zip(
                    NULL_QUADRATIC_ALLOWED, values, strict=True
                )
            },
            "rank": rank,
        })
        if rank == target_linear_rank and rank_certificate is None:
            for selected_rows in combinations(
                range(len(orders)), target_linear_rank
            ):
                for selected_columns in combinations(
                    range(len(NULL_QUADRATIC_FORBIDDEN)),
                    target_linear_rank,
                ):
                    minor = sp.Matrix([
                        [
                            evaluated[row][column]
                            for column in selected_columns
                        ]
                        for row in selected_rows
                    ])
                    determinant = int(minor.det()) % prime
                    if determinant:
                        rank_certificate = {
                            "allowed_values": sample_records[-1][
                                "allowed_values"
                            ],
                            "moment_orders": [
                                orders[index] for index in selected_rows
                            ],
                            "normal_coordinates": [
                                PARAMETERS[
                                    NULL_QUADRATIC_FORBIDDEN[index]
                                ]
                                for index in selected_columns
                            ],
                            "determinant_mod_prime": determinant,
                        }
                        break
                if rank_certificate is not None:
                    break
        if rank_certificate is not None:
            break

    exact_linear_data, _quotient_minors = (
        exact_linear_symbol_factorization(matrix, prime)
    )
    result = {
        "status": (
            "exact characteristic-zero generic-rank certificate via "
            "J_sync-linear symbols reduced at one good prime"
        ),
        "prime": prime,
        "normalization": "F_2=X^2",
        "allowed_coordinates": [
            PARAMETERS[index] for index in NULL_QUADRATIC_ALLOWED
        ],
        "forbidden_coordinates": [
            PARAMETERS[index] for index in NULL_QUADRATIC_FORBIDDEN
        ],
        "moment_orders": list(orders),
        "orders": order_records,
        "sample_ranks": sample_records,
        "maximum_rank": maximum_rank,
        "rank_certificate": rank_certificate,
        "active_linear_symbol_matrix": [
            [
                serialize_sparse_mod(
                    entry,
                    tuple(
                        PARAMETERS[index]
                        for index in NULL_QUADRATIC_ALLOWED
                    ),
                    prime,
                )
                for entry in row
            ]
            for row in matrix
            if any(entry for entry in row)
        ],
        "exact_linear_symbol_factorization": exact_linear_data,
        "interpretation": (
            "rank three gives three generic first-order pivots on the "
            "null-quadratic chart; the other four normal directions "
            "begin in quadratic or cubic order"
        ),
    }
    return result, matrix


def exact_linear_symbol_factorization(
    matrix: list[list[dict[tuple[int, ...], int]]],
    prime: int,
) -> tuple[dict[str, object], list[str]]:
    """Reconstruct and factor the rank-three minors exactly over QQ."""

    active_matrix = [
        row for row in matrix if any(entry for entry in row)
    ]
    assert len(active_matrix) == 3
    allowed_symbols = sp.symbols(
        " ".join(PARAMETERS[index] for index in NULL_QUADRATIC_ALLOWED)
    )
    allowed_by_index = dict(
        zip(
            NULL_QUADRATIC_ALLOWED,
            allowed_symbols,
            strict=True,
        )
    )
    x, y = sp.symbols("x y")
    base = sp.Integer(0)
    normal_basis = {
        index: sp.Integer(0)
        for index in NULL_QUADRATIC_FORBIDDEN
    }
    for row, column, parameter, coefficient in COEFFICIENT_MAP:
        monomial = x**row * y**column
        if parameter in allowed_by_index:
            base += (
                coefficient
                * allowed_by_index[parameter]
                * monomial
            )
        elif parameter in normal_basis:
            normal_basis[parameter] += coefficient * monomial
    # The r0=1 highest-weight quadratic has entries
    # (1,0), 2*(2,1), (3,2).
    base += x + 2 * x**2 * y + x**3 * y**2

    exact_rows = []
    for order in range(2, 5):
        base_power = sp.Poly(sp.expand(base ** (order - 1)), x, y)
        row_entries = []
        for forbidden_index in NULL_QUADRATIC_FORBIDDEN:
            product = sp.Poly(
                sp.expand(
                    base_power.as_expr()
                    * normal_basis[forbidden_index]
                ),
                x,
                y,
            )
            entry = sp.expand(
                order
                * sum(
                    factorial(3 * order - diagonal)
                    * factorial(diagonal)
                    * product.coeff_monomial(
                        x**diagonal * y**diagonal
                    )
                    for diagonal in range(3 * order + 1)
                )
            )
            row_entries.append(entry)
        exact_rows.append(row_entries)
    symbolic_matrix = sp.Matrix(exact_rows)

    for exact_row, modular_row in zip(
        exact_rows, active_matrix, strict=True
    ):
        for exact_entry, modular_entry in zip(
            exact_row, modular_row, strict=True
        ):
            polynomial = sp.Poly(
                exact_entry, *allowed_symbols
            )
            reduced = {
                exponents: int(coefficient) % prime
                for exponents, coefficient in polynomial.terms()
                if int(coefficient) % prime
            }
            assert reduced == modular_entry

    minors = []
    minor_columns = []
    for columns in combinations(range(7), 3):
        determinant = sp.Poly(
            symbolic_matrix[:, columns].det(),
            *allowed_symbols,
            domain=sp.QQ,
        )
        if determinant.is_zero:
            continue
        minors.append(determinant)
        minor_columns.append(columns)
    assert minors
    common_factor = minors[0]
    for determinant in minors[1:]:
        common_factor = sp.gcd(common_factor, determinant)
    quotients = []
    for determinant in minors:
        quotient, remainder = sp.div(
            determinant, common_factor
        )
        assert remainder.is_zero
        quotients.append(quotient)

    quotient_strings = []
    for quotient in quotients:
        modular = sp.Poly(
            quotient.as_expr(),
            *allowed_symbols,
            modulus=prime,
        )
        quotient_strings.append(
            str(modular.as_expr()).replace("**", "^")
        )
    return {
        "active_moment_orders": [2, 3, 4],
        "exact_matrix": [
            [str(sp.factor(entry)) for entry in row]
            for row in exact_rows
        ],
        "nonzero_rank_three_minors": len(minors),
        "nonzero_minor_column_sets": [
            [
                PARAMETERS[NULL_QUADRATIC_FORBIDDEN[index]]
                for index in columns
            ]
            for columns in minor_columns
        ],
        "common_factor": str(sp.factor(common_factor.as_expr())),
        "common_factor_factorization": [
            [str(factor.as_expr()), exponent]
            for factor, exponent in sp.factor_list(
                common_factor.as_expr()
            )[1]
        ],
        "common_factor_degree": common_factor.total_degree(),
        "quotient_minor_count": len(quotients),
        "status": (
            "exact characteristic-zero biform differentiation and gcd"
        ),
    }, quotient_strings


def run_normal_rank_locus(
    singular: str,
    matrix: list[list[dict[tuple[int, ...], int]]],
    generic_rank: int,
    prime: int,
    timeout: int,
) -> dict[str, object]:
    """Compute the five-variable degeneracy locus of the linear symbols."""

    variables = tuple(
        PARAMETERS[index] for index in NULL_QUADRATIC_ALLOWED
    )
    active_matrix = [
        row for row in matrix if any(entry for entry in row)
    ]
    assert len(active_matrix) >= generic_rank
    factorization, quotient_minors = exact_linear_symbol_factorization(
        matrix, prime
    )
    entries = ",".join(
        serialize_sparse_mod(entry, variables, prime)
        for row in active_matrix
        for entry in row
    )
    code = f"""
option(redSB);
ring R={prime},({",".join(variables)}),dp;
matrix M[{len(active_matrix)}][7]={entries};
ideal D=minor(M,{generic_rank});
ideal G=std(D);
ideal Q={",".join(quotient_minors)};
ideal GQ=std(Q);
print("META "+string(size(D))+" "+string(dim(G))+" "
  +string(mult(G))+" "+string(size(G))+" "
  +string(dim(GQ))+" "+string(mult(GQ))+" "+string(size(GQ)));
exit;
"""
    with tempfile.TemporaryDirectory(
        prefix="bidegree33-null-symbols-"
    ) as temporary:
        path = Path(temporary) / "rank_locus.sing"
        path.write_text(code, encoding="utf-8")
        try:
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "status": "timeout",
                "timeout_seconds": timeout,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
    meta = re.search(
        r"(?m)^META (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    return {
        "status": (
            "completed"
            if completed.returncode == 0 and meta is not None
            else "failed"
        ),
        "returncode": completed.returncode,
        "generic_rank": generic_rank,
        "active_symbol_rows": len(active_matrix),
        "maximal_minors": int(meta.group(1)) if meta else None,
        "dimension": int(meta.group(2)) if meta else None,
        "multiplicity": int(meta.group(3)) if meta else None,
        "groebner_basis_size": int(meta.group(4)) if meta else None,
        "quotient_minor_locus_dimension": (
            int(meta.group(5)) if meta else None
        ),
        "quotient_minor_locus_multiplicity": (
            int(meta.group(6)) if meta else None
        ),
        "quotient_minor_groebner_basis_size": (
            int(meta.group(7)) if meta else None
        ),
        "minor_factorization": factorization,
        "stdout_tail": completed.stdout[-3000:],
        "stderr_tail": completed.stderr[-2000:],
        "interpretation": (
            "finite-field rank-drop locus of the exact linear normal "
            "symbol matrix; positive dimension requires higher-symbol "
            "analysis on its components"
        ),
    }


def exceptional_normal_strata() -> dict[str, object]:
    """Describe the exact support of the quotient-minor rank locus."""

    s1, s2, t0, t1 = sp.symbols("s1 s2 t0 t1")
    p_factor = (
        27 * s2**3
        - 468 * s2**2
        - 156 * s2 * t1**2
        + 429 * s2
        - 572 * t1**2
        - 429
    )
    first_quadratic = 3 * s2**2 - 3 * s2 + 7
    second_quadratic = 4 * t1**2 + 25
    w_entry = 3 * s2**2 - 3 * s2 + 4 * t1**2 + 7
    quotient_generators = [
        s2 * t1 * (s2 + 2),
        s2 * w_entry,
        t1**2 * (s2 + 2),
        3 * s1 * w_entry - 14 * t0 * t1 * (s2 + 2),
        t1 * w_entry,
    ]
    component_ideals = [
        [t1, first_quadratic],
        [s2 + 2, second_quadratic],
        [s1, s2, t1],
    ]
    for component in component_ideals:
        basis = sp.groebner(
            component, s1, s2, t0, t1, domain=sp.QQ
        )
        for generator in quotient_generators:
            assert basis.reduce(generator)[1] == 0
    first_remainder = sp.rem(
        p_factor.subs(t1, 0), first_quadratic, s2
    )
    second_remainder = sp.rem(
        p_factor.subs(s2, -2), second_quadratic, t1
    )
    lower_value = p_factor.subs({s1: 0, s2: 0, t1: 0})
    assert sp.expand(first_remainder + 75 * (s2 - 8)) == 0
    assert first_quadratic.subs(s2, 8) == 175
    assert second_remainder == -1750
    assert lower_value == -429

    p_section = sp.Poly(
        p_factor.subs(t1, 0) / 3, s2, domain=sp.QQ
    )
    assert p_section.is_irreducible
    assert sp.Poly(first_quadratic, s2, domain=sp.QQ).is_irreducible
    assert sp.Poly(second_quadratic, t1, domain=sp.QQ).is_irreducible

    return {
        "status": (
            "exact characteristic-zero set-theoretic decomposition "
            "and disjointness certificates"
        ),
        "common_cubic_divisor": str(p_factor),
        "quotient_minor_support_generators_up_to_units": [
            str(sp.factor(generator))
            for generator in quotient_generators
        ],
        "set_theoretic_branch_certificate": (
            "If W is nonzero, s2*W=t1*W=0 and "
            "3*s1*W-14*t0*t1*(s2+2)=0 force "
            "(s1,s2,t1)=0. If W=0, "
            "t1^2*(s2+2)=0 gives t1=0 or s2=-2, producing "
            "the first or second quadratic component."
        ),
        "quotient_minor_support": [
            {
                "label": "residual_component_A",
                "ideal": [
                    "t1",
                    str(first_quadratic),
                ],
                "dimension_in_allowed_base": 3,
                "field_of_sample": (
                    "Q(alpha), 3*alpha^2-3*alpha+7=0"
                ),
                "P_remainder": str(sp.factor(first_remainder)),
                "disjointness_certificate": (
                    "P=-75*(s2-8) modulo the component and "
                    "(3*s2^2-3*s2+7)|_(s2=8)=175"
                ),
            },
            {
                "label": "residual_component_B",
                "ideal": [
                    "s2+2",
                    str(second_quadratic),
                ],
                "dimension_in_allowed_base": 3,
                "field_of_sample": "Q(i), t1=5*i/2",
                "P_remainder": str(second_remainder),
                "disjointness_certificate": (
                    "P=-1750 modulo the component"
                ),
            },
            {
                "label": "residual_lower_locus",
                "ideal": ["s1", "s2", "t1"],
                "dimension_in_allowed_base": 2,
                "field_of_sample": "Q",
                "P_remainder": str(lower_value),
                "disjointness_certificate": (
                    "P=-429 on the lower locus"
                ),
            },
        ],
        "P_section_sample": {
            "ideal": [
                "t1",
                "9*s2^3-156*s2^2+143*s2-143",
            ],
            "field_of_sample": (
                "Q(beta), 9*beta^3-156*beta^2+143*beta-143=0"
            ),
            "minimal_polynomial_irreducible": True,
        },
        "P_intersection_quotient_minor_support": "empty",
        "interpretation": (
            "the linear-rank exceptional set is the disjoint union of "
            "the cubic divisor P=0 and the displayed quotient-minor "
            "support; each component can therefore be tested separately"
        ),
    }


def rref_mod(
    matrix: list[list[int]],
    prime: int,
) -> tuple[list[list[int]], list[int]]:
    reduced = [
        [entry % prime for entry in row]
        for row in matrix
    ]
    row = 0
    pivots = []
    for column in range(len(reduced[0])):
        pivot = next(
            (
                candidate
                for candidate in range(row, len(reduced))
                if reduced[candidate][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[row], reduced[pivot] = reduced[pivot], reduced[row]
        inverse = pow(reduced[row][column], -1, prime)
        reduced[row] = [
            entry * inverse % prime for entry in reduced[row]
        ]
        for other in range(len(reduced)):
            if other == row or not reduced[other][column]:
                continue
            scale = reduced[other][column]
            reduced[other] = [
                (left - scale * right) % prime
                for left, right in zip(
                    reduced[other], reduced[row], strict=True
                )
            ]
        pivots.append(column)
        row += 1
        if row == len(reduced):
            break
    return reduced, pivots


def run_generic_residual_normal_probe(
    singular: str,
    prime: int,
    maximum_order: int,
    timeout: int,
    power_bound: int,
) -> dict[str, object]:
    """Test the four residual normal directions on one exact open chart."""

    allowed_values = (20, 27, 36, 47, 60)
    normal_symbols = sp.symbols(
        " ".join(
            PARAMETERS[index]
            for index in NULL_QUADRATIC_FORBIDDEN
        )
    )
    evaluated_forms = {}
    minimum_degrees = {}
    for order in range(2, maximum_order + 1):
        terms = restricted_moment_terms_mod(order, "null", prime)
        minimum_degree = min(
            sum(
                exponents[index]
                for index in NULL_QUADRATIC_FORBIDDEN
            )
            for exponents in terms
        )
        minimum_degrees[order] = minimum_degree
        normal_terms: dict[tuple[int, ...], int] = defaultdict(int)
        for exponents, coefficient in terms.items():
            normal_exponents = tuple(
                exponents[index]
                for index in NULL_QUADRATIC_FORBIDDEN
            )
            if sum(normal_exponents) != minimum_degree:
                continue
            allowed_exponents = tuple(
                exponents[index]
                for index in NULL_QUADRATIC_ALLOWED
            )
            scalar = coefficient * monomial_value(
                list(allowed_values), allowed_exponents, prime
            )
            normal_terms[normal_exponents] = (
                normal_terms[normal_exponents] + scalar
            ) % prime
        evaluated_forms[order] = {
            exponents: coefficient
            for exponents, coefficient in normal_terms.items()
            if coefficient
        }

    linear_matrix = []
    for order in (2, 3, 4):
        row = []
        for column in range(len(normal_symbols)):
            exponent = [0] * len(normal_symbols)
            exponent[column] = 1
            row.append(
                evaluated_forms[order].get(tuple(exponent), 0)
            )
        linear_matrix.append(row)
    reduced, pivot_columns = rref_mod(linear_matrix, prime)
    assert pivot_columns == [1, 2, 3]
    residual_columns = [
        column
        for column in range(len(normal_symbols))
        if column not in pivot_columns
    ]
    residual_symbols = tuple(
        normal_symbols[column] for column in residual_columns
    )
    substitutions = {}
    for row, pivot_column in enumerate(pivot_columns):
        substitutions[normal_symbols[pivot_column]] = -sum(
            reduced[row][column] * normal_symbols[column]
            for column in residual_columns
        )

    residual_forms = {}
    for order in range(5, maximum_order + 1):
        expression = sum(
            coefficient
            * sp.prod(
                symbol**exponent
                for symbol, exponent in zip(
                    normal_symbols, exponents, strict=True
                )
            )
            for exponents, coefficient in evaluated_forms[order].items()
        )
        residual_forms[order] = sp.Poly(
            sp.expand(expression.subs(substitutions)),
            *residual_symbols,
            modulus=prime,
        )
    nonzero_forms = {
        order: polynomial
        for order, polynomial in residual_forms.items()
        if not polynomial.is_zero
    }
    assert nonzero_forms
    serialized = [
        str(polynomial.as_expr()).replace("**", "^")
        for polynomial in nonzero_forms.values()
    ]
    residual_names = tuple(map(str, residual_symbols))
    power_code = ""
    for index, variable in enumerate(residual_names):
        power_code += f"""
poly h{index}={variable};
int e{index}=1;
while ((e{index}<={power_bound}) && (reduce(h{index},G)!=0))
{{
  h{index}=h{index}*{variable};
  e{index}=e{index}+1;
}}
if (reduce(h{index},G)==0)
{{
  print("POWER {variable} "+string(e{index}));
}}
else
{{
  print("POWER {variable} 0");
}}
"""
    code = f"""
option(redSB);
ring R={prime},({",".join(residual_names)}),dp;
ideal I={",".join(serialized)};
ideal G=std(I);
int dimension=dim(G);
int length=-1;
if (dimension==0) {{ length=vdim(G); }}
print("META "+string(dimension)+" "+string(length)+" "+string(size(G)));
{power_code}
exit;
"""
    with tempfile.TemporaryDirectory(
        prefix="bidegree33-null-residual-"
    ) as temporary:
        path = Path(temporary) / "residual.sing"
        path.write_text(code, encoding="utf-8")
        try:
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "status": "timeout",
                "timeout_seconds": timeout,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
    meta = re.search(
        r"(?m)^META (-?\d+) (-?\d+) (\d+)$",
        completed.stdout,
    )
    powers = {
        match.group(1): int(match.group(2))
        for match in re.finditer(
            r"(?m)^POWER (\w+) (\d+)$", completed.stdout
        )
    }
    return {
        "status": (
            "completed"
            if completed.returncode == 0 and meta is not None
            else "failed"
        ),
        "returncode": completed.returncode,
        "prime": prime,
        "allowed_values": {
            PARAMETERS[index]: value
            for index, value in zip(
                NULL_QUADRATIC_ALLOWED,
                allowed_values,
                strict=True,
            )
        },
        "linear_pivot_coordinates": [
            str(normal_symbols[column]) for column in pivot_columns
        ],
        "residual_normal_coordinates": list(residual_names),
        "minimum_normal_degrees": {
            str(order): degree
            for order, degree in minimum_degrees.items()
        },
        "nonzero_residual_symbol_orders": list(nonzero_forms),
        "residual_symbol_term_counts": {
            str(order): len(polynomial.terms())
            for order, polynomial in nonzero_forms.items()
        },
        "dimension": int(meta.group(1)) if meta else None,
        "quotient_length": int(meta.group(2)) if meta else None,
        "groebner_basis_size": int(meta.group(3)) if meta else None,
        "coordinate_power_memberships": powers,
        "power_search_bound": power_bound,
        "stdout_tail": completed.stdout[-3000:],
        "stderr_tail": completed.stderr[-2000:],
        "interpretation": (
            "a zero-dimensional homogeneous residual symbol ideal "
            "proves transverse isolation on a nonempty characteristic-"
            "zero open subset after the three linear pivots; it does "
            "not cover the pivot or rank-drop divisors"
        ),
    }


def run_complete_normal_fiber_at_base(
    singular: str,
    prime: int,
    maximum_order: int,
    timeout: int,
    power_bound: int,
    allowed_values: tuple[int, int, int, int, int],
    stratum: str,
    characteristic_zero_lift: str,
) -> dict[str, object]:
    """Compute the complete seven-normal fiber at a specified base point."""

    allowed_values = tuple(value % prime for value in allowed_values)
    normal_names = tuple(
        PARAMETERS[index] for index in NULL_QUADRATIC_FORBIDDEN
    )
    evaluated_forms = {}
    for order in range(2, maximum_order + 1):
        terms = restricted_moment_terms_mod(order, "null", prime)
        normal_terms: dict[tuple[int, ...], int] = defaultdict(int)
        for exponents, coefficient in terms.items():
            normal_exponents = tuple(
                exponents[index]
                for index in NULL_QUADRATIC_FORBIDDEN
            )
            allowed_exponents = tuple(
                exponents[index]
                for index in NULL_QUADRATIC_ALLOWED
            )
            scalar = coefficient * monomial_value(
                list(allowed_values), allowed_exponents, prime
            )
            normal_terms[normal_exponents] = (
                normal_terms[normal_exponents] + scalar
            ) % prime
        evaluated_forms[order] = {
            exponents: coefficient
            for exponents, coefficient in normal_terms.items()
            if coefficient
        }
    serialized = [
        serialize_sparse_mod(
            evaluated_forms[order], normal_names, prime
        )
        for order in range(2, maximum_order + 1)
    ]
    power_code = ""
    for index, variable in enumerate(normal_names):
        power_code += f"""
poly h{index}={variable};
int e{index}=1;
while ((e{index}<={power_bound}) && (reduce(h{index},G)!=0))
{{
  h{index}=h{index}*{variable};
  e{index}=e{index}+1;
}}
if (reduce(h{index},G)==0)
{{
  print("POWER {variable} "+string(e{index}));
}}
else
{{
  print("POWER {variable} 0");
}}
"""
    code = f"""
option(redSB);
ring R={prime},({",".join(normal_names)}),dp;
ideal I={",".join(serialized)};
ideal G=std(I);
int dimension=dim(G);
int length=-1;
if (dimension==0) {{ length=vdim(G); }}
print("META "+string(dimension)+" "+string(length)+" "+string(size(G)));
{power_code}
exit;
"""
    with tempfile.TemporaryDirectory(
        prefix="bidegree33-null-complete-normal-"
    ) as temporary:
        path = Path(temporary) / "normal_fiber.sing"
        path.write_text(code, encoding="utf-8")
        try:
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "status": "timeout",
                "timeout_seconds": timeout,
                "prime": prime,
                "stratum": stratum,
                "characteristic_zero_lift": characteristic_zero_lift,
                "allowed_values": {
                    PARAMETERS[index]: value
                    for index, value in zip(
                        NULL_QUADRATIC_ALLOWED,
                        allowed_values,
                        strict=True,
                    )
                },
                "moment_term_counts": {
                    str(order): len(terms)
                    for order, terms in evaluated_forms.items()
                },
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
    meta = re.search(
        r"(?m)^META (-?\d+) (-?\d+) (\d+)$",
        completed.stdout,
    )
    powers = {
        match.group(1): int(match.group(2))
        for match in re.finditer(
            r"(?m)^POWER (\w+) (\d+)$", completed.stdout
        )
    }
    return {
        "status": (
            "completed"
            if completed.returncode == 0 and meta is not None
            else "failed"
        ),
        "returncode": completed.returncode,
        "prime": prime,
        "stratum": stratum,
        "characteristic_zero_lift": characteristic_zero_lift,
        "allowed_values": {
            PARAMETERS[index]: value
            for index, value in zip(
                NULL_QUADRATIC_ALLOWED,
                allowed_values,
                strict=True,
            )
        },
        "normal_coordinates": list(normal_names),
        "moment_orders": list(range(2, maximum_order + 1)),
        "moment_term_counts": {
            str(order): len(terms)
            for order, terms in evaluated_forms.items()
        },
        "dimension": int(meta.group(1)) if meta else None,
        "quotient_length": int(meta.group(2)) if meta else None,
        "groebner_basis_size": int(meta.group(3)) if meta else None,
        "coordinate_power_memberships": powers,
        "power_search_bound": power_bound,
        "stdout_tail": completed.stdout[-3000:],
        "stderr_tail": completed.stderr[-2000:],
        "interpretation": (
            "dimension zero at this good reduction proves transverse "
            "isolation at the displayed characteristic-zero algebraic "
            "lift and hence on a nonempty open subset of its stratum; "
            "coordinate power memberships are finite-field evidence"
        ),
    }


def run_generic_complete_normal_fiber(
    singular: str,
    prime: int,
    maximum_order: int,
    timeout: int,
    power_bound: int,
) -> dict[str, object]:
    """Compute the complete seven-normal fiber at one generic base point."""

    return run_complete_normal_fiber_at_base(
        singular,
        prime,
        maximum_order,
        timeout,
        power_bound,
        (20, 27, 36, 47, 60),
        "generic_allowed_base",
        "the same rational point over Q",
    )


def run_exceptional_complete_normal_fibers(
    singular: str,
    maximum_order: int,
    timeout: int,
    power_bound: int,
) -> dict[str, object]:
    """Test one good reduction on every exact exceptional rank stratum."""

    samples = [
        {
            "label": "P_divisor",
            "prime": 32003,
            "allowed_values": (20, 27, 12804, 47, 0),
            "characteristic_zero_lift": (
                "s0=20,s1=27,t0=47,t1=0 and s2=beta, where "
                "9*beta^3-156*beta^2+143*beta-143=0"
            ),
            "equations_mod_prime": [
                "t1=0",
                "9*s2^3-156*s2^2+143*s2-143=0",
            ],
        },
        {
            "label": "residual_component_A",
            "prime": 30013,
            "allowed_values": (20, 27, 16605, 47, 0),
            "characteristic_zero_lift": (
                "s0=20,s1=27,t0=47,t1=0 and s2=alpha, where "
                "3*alpha^2-3*alpha+7=0"
            ),
            "equations_mod_prime": [
                "t1=0",
                "3*s2^2-3*s2+7=0",
            ],
        },
        {
            "label": "residual_component_B",
            "prime": 30013,
            "allowed_values": (20, 27, -2, 47, 15619),
            "characteristic_zero_lift": (
                "s0=20,s1=27,t0=47,s2=-2 and t1=5*i/2"
            ),
            "equations_mod_prime": [
                "s2=-2",
                "4*t1^2+25=0",
            ],
        },
        {
            "label": "residual_lower_locus",
            "prime": 32003,
            "allowed_values": (20, 0, 0, 47, 0),
            "characteristic_zero_lift": (
                "the same rational point over Q"
            ),
            "equations_mod_prime": ["s1=0", "s2=0", "t1=0"],
        },
    ]
    results = []
    for sample in samples:
        prime = int(sample["prime"])
        values = tuple(sample["allowed_values"])
        assert len(values) == len(NULL_QUADRATIC_ALLOWED)
        s0, s1, s2, t0, t1 = (
            int(value) % prime for value in values
        )
        p_value = (
            27 * s2**3
            - 468 * s2**2
            - 156 * s2 * t1**2
            + 429 * s2
            - 572 * t1**2
            - 429
        ) % prime
        if sample["label"] == "P_divisor":
            assert p_value == 0
            root_derivative = (
                27 * s2**2 - 312 * s2 + 143
            ) % prime
            assert root_derivative != 0
        else:
            assert p_value != 0
        if sample["label"] == "residual_component_A":
            assert t1 == 0
            assert (3 * s2**2 - 3 * s2 + 7) % prime == 0
            root_derivative = (6 * s2 - 3) % prime
            assert root_derivative != 0
        elif sample["label"] == "residual_component_B":
            assert (s2 + 2) % prime == 0
            assert (4 * t1**2 + 25) % prime == 0
            root_derivative = (8 * t1) % prime
            assert root_derivative != 0
        elif sample["label"] == "residual_lower_locus":
            assert s1 == s2 == t1 == 0
            root_derivative = None
        result = run_complete_normal_fiber_at_base(
            singular,
            prime,
            maximum_order,
            timeout,
            power_bound,
            tuple(int(value) for value in values),
            str(sample["label"]),
            str(sample["characteristic_zero_lift"]),
        )
        result["stratum_equations_mod_prime"] = sample[
            "equations_mod_prime"
        ]
        result["P_value_mod_prime"] = p_value
        result["simple_algebraic_root_mod_prime"] = (
            root_derivative is None or root_derivative != 0
        )
        result["algebraic_root_derivative_mod_prime"] = root_derivative
        results.append(result)
    return {
        "status": (
            "completed"
            if all(result["status"] == "completed" for result in results)
            else "incomplete"
        ),
        "samples": results,
        "interpretation": (
            "a zero-dimensional result proves a nonempty "
            "characteristic-zero transverse-isolation open subset on "
            "the corresponding exact irreducible stratum"
        ),
    }


def quadratic_gradient_row_mod(
    point: list[list[int]],
    component: int,
    prime: int,
) -> list[int]:
    """Differentiate q_(2*component)=tr(A_component^2) by evaluation."""

    projectors = casimir_projectors_mod(3, prime)
    size = 4
    values = []
    for row in range(size):
        for column in range(size):
            plus = [entries[:] for entries in point]
            minus = [entries[:] for entries in point]
            plus[row][column] = (plus[row][column] + 1) % prime
            minus[row][column] = (minus[row][column] - 1) % prime
            plus_component = component_matrices_mod(
                plus, 3, projectors, prime
            )[component]
            minus_component = component_matrices_mod(
                minus, 3, projectors, prime
            )[component]
            plus_value = sum(
                plus_component[i][j] * plus_component[j][i]
                for i in range(size)
                for j in range(size)
            ) % prime
            minus_value = sum(
                minus_component[i][j] * minus_component[j][i]
                for i in range(size)
                for j in range(size)
            ) % prime
            values.append((plus_value - minus_value) * pow(2, -1, prime) % prime)
    return values


def degree_system_record(
    label: str,
    moment_orders: tuple[int, ...],
    casimir_components: tuple[int, ...],
    hilbert: list[int],
    moment_rows: list[list[int]],
    quadratic_rows: dict[int, list[int]],
    prime: int,
) -> dict[str, object]:
    degrees = (
        tuple(moment_orders)
        + (2,) * len(casimir_components)
    )
    assert len(degrees) == QUOTIENT_DIMENSION
    predicted_top = sum(degrees) - AMBIENT_DIMENSION
    numerator = hilbert_numerator(hilbert, degrees)
    first_negative = next(
        (
            [index, coefficient]
            for index, coefficient in enumerate(numerator)
            if coefficient < 0
        ),
        None,
    )
    first_tail = next(
        (
            [index, coefficient]
            for index, coefficient in enumerate(
                numerator[predicted_top + 1 :],
                start=predicted_top + 1,
            )
            if coefficient
        ),
        None,
    )
    rows = [moment_rows[order - 1] for order in moment_orders]
    rows.extend(quadratic_rows[component] for component in casimir_components)
    return {
        "label": label,
        "moment_orders": list(moment_orders),
        "casimirs": [f"q_{2 * component}" for component in casimir_components],
        "degrees": list(degrees),
        "total_invariant_degree": sum(degrees),
        "degree_product": (
            product(degrees)
        ),
        "predicted_hilbert_numerator_top_degree": predicted_top,
        "first_negative_through_cutoff": first_negative,
        "first_nonzero_after_predicted_top_through_cutoff": first_tail,
        "hilbert_checked_through_degree": len(hilbert) - 1,
        "numerator_coefficient_sum_through_predicted_top": sum(
            numerator[: predicted_top + 1]
        ),
        "jacobian_rank_mod_prime": rank_mod(rows, prime),
    }


def product(values: tuple[int, ...]) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def degree_comparison(prime: int) -> dict[str, object]:
    hilbert = invariant_hilbert_coefficients(3, HILBERT_CUTOFF)
    point = deterministic_point(3, 0, prime)
    moment_rows = moment_jacobian_mod(
        point, 3, MOMENT_CUTOFF, prime
    )
    quadratic_rows = {
        component: quadratic_gradient_row_mod(point, component, prime)
        for component in (1, 2)
    }

    systems = [
        (
            "corrected_moments",
            tuple(range(1, 13)) + (14,),
            (),
        ),
        (
            "low_degree_q2",
            tuple(range(1, 13)),
            (1,),
        ),
    ]
    records = [
        degree_system_record(
            label,
            moment_orders,
            casimirs,
            hilbert,
            moment_rows,
            quadratic_rows,
            prime,
        )
        for label, moment_orders, casimirs in systems
    ]

    matched_searches = []
    for casimir_components in ((1,), (1, 2)):
        extra_count = QUOTIENT_DIMENSION - len(casimir_components) - 2
        extra_sum = (
            92 - 1 - 2 - 2 * len(casimir_components)
        )
        candidates = []
        degree_sum_matches = 0
        hilbert_matches = 0
        for extra_orders in combinations(
            range(3, MOMENT_CUTOFF + 1), extra_count
        ):
            if sum(extra_orders) != extra_sum:
                continue
            degree_sum_matches += 1
            moment_orders = (1, 2) + extra_orders
            record = degree_system_record(
                "candidate",
                moment_orders,
                casimir_components,
                hilbert,
                moment_rows,
                quadratic_rows,
                prime,
            )
            if (
                record["first_negative_through_cutoff"] is not None
                or record[
                    "first_nonzero_after_predicted_top_through_cutoff"
                ]
                is not None
            ):
                continue
            hilbert_matches += 1
            if record["jacobian_rank_mod_prime"] != QUOTIENT_DIMENSION:
                continue
            candidates.append(record)
        candidates.sort(
            key=lambda record: (
                max(record["moment_orders"]),
                record["degree_product"],
                record["moment_orders"],
            )
        )
        assert candidates
        selected = candidates[0]
        selected["label"] = (
            "matched_total_q2"
            if len(casimir_components) == 1
            else "matched_total_q2_q4"
        )
        records.append(selected)
        matched_searches.append({
            "casimirs": selected["casimirs"],
            "target_total_invariant_degree": 92,
            "moment_search_range": [3, MOMENT_CUTOFF],
            "degree_sum_matches": degree_sum_matches,
            "hilbert_compatible_matches": hilbert_matches,
            "hilbert_compatible_full_rank_matches": len(candidates),
            "selection_rule": (
                "minimize largest moment order, then degree product, "
                "then lexicographic moment-order tuple"
            ),
            "selected_label": selected["label"],
        })
    assert all(
        record["jacobian_rank_mod_prime"] == QUOTIENT_DIMENSION
        for record in records
    )
    return {
        "status": (
            "exact Hilbert coefficient arithmetic and characteristic-zero "
            "Jacobian rank certificates via one good prime"
        ),
        "prime": prime,
        "quadratic_anchor_identity": quadratic_anchor_identity(),
        "systems": records,
        "matched_total_degree_searches": matched_searches,
        "mu14_casimir_class_test": mu14_casimir_class_test(prime),
    }


def quadratic_anchor_identity() -> dict[str, str]:
    """Verify the exact normalization q_2=80*(r1^2-r0*r2)."""

    r0, r1, r2 = sp.symbols("r0 r1 r2")
    coefficients = [[sp.Integer(0)] * 4 for _ in range(4)]
    for row, column, coefficient in (
        (1, 0, r0),
        (2, 1, 2 * r0),
        (3, 2, r0),
        (0, 0, -r1),
        (1, 1, -r1),
        (2, 2, r1),
        (3, 3, r1),
        (0, 1, -r2),
        (1, 2, -2 * r2),
        (2, 3, -r2),
    ):
        coefficients[row][column] += coefficient
    q2 = component_quadratic_exact_sparse(
        coefficients, 3, 1
    )
    discriminant = r1**2 - r0 * r2
    assert sp.expand(q2 - 80 * discriminant) == 0
    return {
        "q_2": str(q2),
        "Delta_2": str(discriminant),
        "identity": "q_2=80*Delta_2",
        "status": "exact symbolic",
    }


def mu14_casimir_class_test(
    prime: int,
    extra_samples: int = 3,
) -> dict[str, object]:
    """Compare mu_14 and q_2^7 modulo the lower-moment monomial ideal.

    Every degree-14 polynomial in mu_1,...,mu_12,q_2 that vanishes after
    mu_1=...=mu_12=0 is included, not a bounded subset of monomials.
    Full column rank modulo one prime is therefore an exact nonrelation
    certificate over characteristic zero for this generated subalgebra.
    """

    projectors = casimir_projectors_mod(3, prime)

    def one_test(
        components: tuple[int, ...],
        sample_offset: int,
    ) -> dict[str, object]:
        weights = tuple(range(1, 13)) + (2,) * len(components)
        exponents = weighted_exponents(weights, 14)
        lower_ideal_exponents = [
            exponent_tuple
            for exponent_tuple in exponents
            if any(exponent_tuple[:12])
        ]
        pure_casimir_exponents = [
            exponent_tuple
            for exponent_tuple in exponents
            if not any(exponent_tuple[:12])
        ]
        assert (
            len(lower_ideal_exponents) + len(pure_casimir_exponents)
            == len(exponents)
        )
        sample_count = len(exponents) + 1 + extra_samples
        lower_matrix = []
        pure_matrix = []
        mu14_column = []
        for sample_index in range(1, sample_count + 1):
            point = deterministic_point(
                3, sample_offset + sample_index, prime
            )
            moments = moments_mod(point, 3, 14, prime)
            quadratics, _odd = invariant_values_mod(
                point, 3, projectors, prime
            )
            base_values = moments[:12] + [
                quadratics[component] for component in components
            ]
            lower_matrix.append([
                monomial_value(base_values, exponent_tuple, prime)
                for exponent_tuple in lower_ideal_exponents
            ])
            pure_matrix.append([
                monomial_value(base_values, exponent_tuple, prime)
                for exponent_tuple in pure_casimir_exponents
            ])
            mu14_column.append(moments[13])

        rank_lower = rank_mod(lower_matrix, prime)
        lower_pure_matrix = [
            lower + pure
            for lower, pure in zip(
                lower_matrix, pure_matrix, strict=True
            )
        ]
        rank_lower_pure = rank_mod(lower_pure_matrix, prime)
        rank_lower_mu14 = rank_mod(
            [
                row + [mu14]
                for row, mu14 in zip(
                    lower_matrix, mu14_column, strict=True
                )
            ],
            prime,
        )
        rank_all = rank_mod(
            [
                row + [mu14]
                for row, mu14 in zip(
                    lower_pure_matrix, mu14_column, strict=True
                )
            ],
            prime,
        )
        return {
            "casimirs": [
                f"q_{2 * component}" for component in components
            ],
            "degree_14_monomials_total": len(exponents),
            "lower_moment_ideal_columns": len(lower_ideal_exponents),
            "pure_casimir_columns": len(pure_casimir_exponents),
            "samples": sample_count,
            "ranks": {
                "lower_moment_ideal": rank_lower,
                "lower_plus_pure_casimirs": rank_lower_pure,
                "lower_plus_mu14": rank_lower_mu14,
                "lower_plus_pure_casimirs_plus_mu14": rank_all,
            },
            "pure_casimir_classes_independent_mod_lower_moment_span": (
                rank_lower_pure
                == rank_lower + len(pure_casimir_exponents)
            ),
            "mu14_in_pure_casimir_span_mod_lower_moments": (
                rank_all == rank_lower_pure
            ),
            "mu14_independent_from_pure_casimir_span_mod_lower_moments": (
                rank_all == rank_lower_pure + 1
            ),
        }

    tests = [
        one_test((1,), 10_000),
        one_test((1, 2), 20_000),
    ]
    return {
        "prime": prime,
        "weight": 14,
        "lower_moment_generators": [
            f"mu_{order}" for order in range(1, 13)
        ],
        "tests": tests,
        "interpretation": (
            "If mu_14 is independent from the pure Casimir span, there "
            "is no degree-14 identity expressing mu_14 as a polynomial "
            "in the selected quadratic Casimirs modulo the degree-14 "
            "part of (mu_1,...,mu_12) inside the generated algebra. "
            "This does not decide equality of their zero divisors in "
            "the full invariant quotient."
        ),
    }


def run_singular_probe(
    singular: str,
    quadratic_mode: str,
    prime: int,
    maximum_order: int,
    power_bound: int,
    timeout: int,
) -> dict[str, object]:
    moment_orders = tuple(range(2, maximum_order + 1))
    moments = {
        order: primitive_modular_terms(
            restricted_moment_terms_mod(order, quadratic_mode, prime),
            prime,
        )
        for order in moment_orders
    }
    serialized = {
        order: serialize_modular_polynomial(terms, prime)
        for order, terms in moments.items()
    }
    ideal = ",".join(serialized.values())
    variables = ",".join(PARAMETERS)
    forbidden = (
        NULL_QUADRATIC_FORBIDDEN
        if quadratic_mode == "null"
        else ()
    )
    power_code = ""
    for index in forbidden:
        variable = PARAMETERS[index]
        suffix = str(index)
        power_code += f"""
poly h{suffix}={variable};
int exponent{suffix}=1;
while ((exponent{suffix}<={power_bound}) && (reduce(h{suffix},G)!=0))
{{
  h{suffix}=h{suffix}*{variable};
  exponent{suffix}=exponent{suffix}+1;
}}
if (reduce(h{suffix},G)==0)
{{
  print("POWER {variable} "+string(exponent{suffix}));
}}
else
{{
  print("POWER {variable} 0");
}}
"""

    if quadratic_mode == "null":
        synchronized = ",".join(PARAMETERS[index] for index in forbidden)
        containment_code = f"""
ideal J={synchronized};
ideal GJ=std(J);
print("SYNC_REMAINDER_SIZE "+string(size(reduce(I,GJ))));
"""
    else:
        containment_code = ""

    code = f"""
option(redSB);
option(prot);
ring R={prime},({variables}),dp;
ideal I={ideal};
ideal G=std(I);
print("META "+string(dim(G))+" "+string(mult(G))+" "+string(size(G)));
{containment_code}
{power_code}
exit;
"""
    with tempfile.TemporaryDirectory(
        prefix="bidegree33-casimir-"
    ) as temporary:
        path = Path(temporary) / f"{quadratic_mode}.sing"
        path.write_text(code, encoding="utf-8")
        try:
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "quadratic_mode": quadratic_mode,
                "prime": prime,
                "moment_orders": list(moment_orders),
                "moment_term_counts": {
                    str(order): len(terms)
                    for order, terms in moments.items()
                },
                "status": "timeout",
                "timeout_seconds": timeout,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
    meta = re.search(r"(?m)^META (\d+) (\d+) (\d+)$", completed.stdout)
    powers = {
        match.group(1): int(match.group(2))
        for match in re.finditer(
            r"(?m)^POWER (\w+) (\d+)$", completed.stdout
        )
    }
    synchronization = re.search(
        r"(?m)^SYNC_REMAINDER_SIZE (\d+)$", completed.stdout
    )
    return {
        "quadratic_mode": quadratic_mode,
        "prime": prime,
        "moment_orders": list(moment_orders),
        "moment_term_counts": {
            str(order): len(terms)
            for order, terms in moments.items()
        },
        "status": (
            "completed" if completed.returncode == 0 and meta else "failed"
        ),
        "returncode": completed.returncode,
        "dimension": int(meta.group(1)) if meta else None,
        "multiplicity": int(meta.group(2)) if meta else None,
        "groebner_basis_size": int(meta.group(3)) if meta else None,
        "synchronized_linear_ideal_contains_moments": (
            synchronization is not None
            and int(synchronization.group(1)) == 0
        ),
        "forbidden_coordinate_power_memberships": powers,
        "power_search_bound": power_bound,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-2000:],
        "interpretation": (
            "finite-field evidence only; a completed radical equality on "
            "the null chart does not by itself reconstruct an identity "
            "over QQ"
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    parser.add_argument("--maximum-order", type=int, default=12)
    parser.add_argument("--power-bound", type=int, default=40)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument(
        "--run-singular-probes",
        action="store_true",
        help=(
            "also run the older direct standard-basis probes on F2=X^2 "
            "and F2=0; off by default"
        ),
    )
    parser.add_argument(
        "--run-normal-rank-locus",
        action="store_true",
        help=(
            "compute only the new five-variable rank-drop locus of the "
            "J_sync-linear normal symbols"
        ),
    )
    parser.add_argument(
        "--run-residual-normal-probe",
        action="store_true",
        help=(
            "after three linear pivots, test the four residual normal "
            "directions at one exact allowed base point"
        ),
    )
    parser.add_argument(
        "--run-complete-normal-fiber",
        action="store_true",
        help=(
            "compute the full seven-normal moment fiber at one exact "
            "allowed nullcone base point"
        ),
    )
    parser.add_argument(
        "--run-exceptional-normal-fibers",
        action="store_true",
        help=(
            "compute one full seven-normal fiber on P=0 and on every "
            "exact component of the quotient-minor rank locus"
        ),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    normal_symbols, normal_matrix = null_quadratic_normal_symbols(
        arguments.prime, arguments.maximum_order
    )
    payload: dict[str, object] = {
        "status": (
            "exact parameter comparison and characteristic-zero "
            "nonrelation certificates plus the null-quadratic "
            "synchronization normal symbols; direct zero-fiber probes "
            "are optional and were not run"
        ),
        "degree_comparison": degree_comparison(arguments.prime),
        "null_quadratic_normal_symbols": normal_symbols,
        "null_quadratic_exceptional_rank_strata": (
            exceptional_normal_strata()
        ),
    }
    if arguments.run_normal_rank_locus:
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "Singular is required with --run-normal-rank-locus"
            )
        payload["null_quadratic_normal_rank_locus"] = (
            run_normal_rank_locus(
                singular,
                normal_matrix,
                int(normal_symbols["maximum_rank"]),
                arguments.prime,
                arguments.timeout,
            )
        )
    if arguments.run_residual_normal_probe:
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "Singular is required with --run-residual-normal-probe"
            )
        payload["null_quadratic_generic_residual_normal_probe"] = (
            run_generic_residual_normal_probe(
                singular,
                arguments.prime,
                arguments.maximum_order,
                arguments.timeout,
                arguments.power_bound,
            )
        )
    if arguments.run_complete_normal_fiber:
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "Singular is required with --run-complete-normal-fiber"
            )
        payload["null_quadratic_generic_complete_normal_fiber"] = (
            run_generic_complete_normal_fiber(
                singular,
                arguments.prime,
                arguments.maximum_order,
                arguments.timeout,
                arguments.power_bound,
            )
        )
    if arguments.run_exceptional_normal_fibers:
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "Singular is required with "
                "--run-exceptional-normal-fibers"
            )
        payload["null_quadratic_exceptional_complete_normal_fibers"] = (
            run_exceptional_complete_normal_fibers(
                singular,
                arguments.maximum_order,
                arguments.timeout,
                arguments.power_bound,
            )
        )
    if arguments.run_singular_probes:
        payload["status"] = (
            "exact parameter comparison and characteristic-zero "
            "nonrelation certificates plus bounded finite-field "
            "zero-fiber probes"
        )
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "Singular is required with --run-singular-probes"
            )
        payload["zero_fiber_probes"] = [
            run_singular_probe(
                singular,
                mode,
                arguments.prime,
                arguments.maximum_order,
                arguments.power_bound,
                arguments.timeout,
            )
            for mode in ("null", "zero")
        ]
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
