#!/usr/bin/env python3
"""
Exact symbolic audit for the quartic-Hessian-nilpotent Waring-rank gates.

The mathematical proof is in QUARTIC_HN_WARING_RIGIDITY.md.  This checker
verifies the only finite computer-assisted layer used there:

* all codimension-two Gale multiplicity patterns for eight Waring terms,
  with no relation supported on at most two terms;
* injectivity of the square-pair restriction map in every pattern with at
  most four nonzero Gale directions and at most three zero Gale columns;
* the unique exceptional pattern: four zero Gale columns plus four distinct
  nonzero directions, where the kernel has dimension five;
* the exact linear normal forms in the first live rank-nine stratum;
* the K_(3,3) channel split on the triple-root cubic profile.

For four distinct Gale directions we normalize them to
    0, 1, infinity, lambda
and certify a nonzero 28 x 28 minor over QQ(lambda).  Every displayed
determinant is a unit times lambda^4 or lambda^6, so it is nonzero for the
required lambda != 0,1.

Requires: Python 3.11+, SymPy.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable

import sympy as sp

LAM = sp.symbols("lambda")
X = sp.symbols("x0:6")
PAIRS = list(combinations(range(8), 2))


@dataclass(frozen=True)
class Pattern:
    zero_columns: int
    multiplicities: tuple[int, ...]

    @property
    def nonzero_count(self) -> int:
        return sum(self.multiplicities)

    @property
    def directions(self) -> int:
        return len(self.multiplicities)

    @property
    def minimum_relation_support(self) -> int:
        # A nonzero row combination vanishes precisely on one projective
        # Gale class and on every zero column.
        return self.nonzero_count - max(self.multiplicities)


def partitions_nonincreasing(n: int, maximum: int | None = None) -> Iterable[tuple[int, ...]]:
    if n == 0:
        yield ()
        return
    maximum = n if maximum is None else min(maximum, n)
    for first in range(maximum, 0, -1):
        for tail in partitions_nonincreasing(n - first, first):
            yield (first,) + tail


def enumerate_small_direction_patterns() -> list[Pattern]:
    patterns: list[Pattern] = []
    for zero_columns in range(4):
        nonzero = 8 - zero_columns
        for multiplicities in partitions_nonincreasing(nonzero):
            if not 2 <= len(multiplicities) <= 4:
                continue
            pattern = Pattern(zero_columns, multiplicities)
            if pattern.minimum_relation_support >= 3:
                patterns.append(pattern)
    return patterns


def gale_matrix(pattern: Pattern) -> sp.Matrix:
    """Return a canonical 2 x 8 Gale matrix."""
    columns: list[tuple[sp.Expr, sp.Expr]] = [(0, 0)] * pattern.zero_columns

    if pattern.directions == 2:
        directions = [(1, 0), (0, 1)]
    elif pattern.directions == 3:
        directions = [(1, 0), (1, 1), (0, 1)]
    elif pattern.directions == 4:
        directions = [(1, 0), (1, 1), (0, 1), (1, LAM)]
    else:
        raise ValueError(pattern)

    for multiplicity, direction in zip(pattern.multiplicities, directions, strict=True):
        columns.extend([direction] * multiplicity)

    if len(columns) != 8:
        raise AssertionError((pattern, columns))

    return sp.Matrix(
        [
            [column[0] for column in columns],
            [column[1] for column in columns],
        ]
    )


def square_pair_restriction_matrix(pattern: Pattern) -> sp.Matrix:
    """
    Matrix of
        (A_ij) -> sum_{i<j} A_ij l_i^2 l_j^2
    after parameterizing the six-plane ker(K).
    """
    kernel = gale_matrix(pattern).nullspace()
    if len(kernel) != 6:
        raise AssertionError((pattern, len(kernel)))

    parameterization = sp.Matrix.hstack(*kernel)  # 8 x 6
    linear_forms = [
        sp.expand(sum(parameterization[i, j] * X[j] for j in range(6)))
        for i in range(8)
    ]

    columns: list[dict[tuple[int, ...], sp.Expr]] = []
    monomial_index: dict[tuple[int, ...], int] = {}

    for i, j in PAIRS:
        polynomial = sp.Poly(sp.expand(linear_forms[i] ** 2 * linear_forms[j] ** 2), *X)
        terms = dict(polynomial.terms())
        columns.append(terms)
        for monomial in terms:
            monomial_index.setdefault(monomial, len(monomial_index))

    matrix = sp.zeros(len(monomial_index), len(PAIRS))
    for column, terms in enumerate(columns):
        for monomial, coefficient in terms.items():
            matrix[monomial_index[monomial], column] = coefficient

    return matrix


def certified_minor(matrix: sp.Matrix, parameterized: bool) -> tuple[tuple[int, ...], sp.Expr]:
    """
    Select independent rows at lambda=2 and compute the corresponding exact
    symbolic determinant.
    """
    specialization = matrix.subs(LAM, 2) if parameterized else matrix
    _, pivot_rows = specialization.T.rref()
    rows = tuple(int(row) for row in pivot_rows[:28])
    if len(rows) != 28:
        raise AssertionError(("rank below 28", len(rows)))

    determinant = sp.factor(matrix[list(rows), :].det())
    return rows, determinant


EXPECTED_CONSTANT_CASES = {
    (0, (5, 3)): -4096,
    (0, (5, 2, 1)): 4096,
    (0, (4, 4)): 4096,
    (0, (4, 3, 1)): -4096,
    (0, (4, 2, 2)): -8192,
    (0, (3, 3, 2)): -8192,
    (1, (4, 3)): -4096,
    (1, (4, 2, 1)): -4096,
    (1, (3, 3, 1)): -4096,
    (1, (3, 2, 2)): -8192,
    (2, (3, 3)): -4096,
    (2, (3, 2, 1)): -4096,
    (2, (2, 2, 2)): -8192,
    (3, (2, 2, 1)): 4096,
}

EXPECTED_FOUR_DIRECTION_CASES = {
    (0, (5, 1, 1, 1)): 4096 * LAM**6,
    (0, (4, 2, 1, 1)): -8192 * LAM**4,
    (0, (3, 3, 1, 1)): -8192 * LAM**4,
    (0, (3, 2, 2, 1)): 8192,
    (0, (2, 2, 2, 2)): -8192,
    (1, (4, 1, 1, 1)): -4096 * LAM**6,
    (1, (3, 2, 1, 1)): -8192 * LAM**4,
    (1, (2, 2, 2, 1)): -8192,
    (2, (3, 1, 1, 1)): -4096 * LAM**6,
    (2, (2, 2, 1, 1)): -8192 * LAM**4,
    (3, (2, 1, 1, 1)): 4096 * LAM**6,
}



def homogeneous_exponents(total: int, variables: int) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []

    def recurse(prefix: tuple[int, ...], remaining: int, slots: int) -> None:
        if slots == 1:
            result.append(prefix + (remaining,))
            return
        for entry in range(remaining + 1):
            recurse(prefix + (entry,), remaining - entry, slots - 1)

    recurse((), total, variables)
    return result


def coefficient_linear_rows(
    expression: sp.Expr,
    polynomial_variables: tuple[sp.Symbol, ...],
    coefficient_variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    polynomial = sp.Poly(sp.expand(expression), *polynomial_variables)
    rows = [
        [sp.expand(coefficient).coeff(variable) for variable in coefficient_variables]
        for _, coefficient in polynomial.terms()
    ]
    return sp.Matrix(rows) if rows else sp.zeros(0, len(coefficient_variables))


def rank_nine_linear_normal_forms() -> None:
    x1, x2, y1, y2, tau = sp.symbols("x1 x2 y1 y2 tau")
    variables = (x1, x2, y1, y2)
    exponents = homogeneous_exponents(4, 4)
    coefficients = sp.symbols(f"e0:{len(exponents)}")
    generic_e = sum(
        coefficient
        * x1**exponent[0]
        * x2**exponent[1]
        * y1**exponent[2]
        * y2**exponent[3]
        for coefficient, exponent in zip(coefficients, exponents, strict=True)
    )

    split_q = sp.Matrix(
        [
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    generic_m = split_q * sp.hessian(generic_e, variables)

    cubic_profiles = {
        "triple": x1**3,
        "double": x1**2 * x2,
        "squarefree": x1 * x2 * (x1 - x2),
    }
    expected_dimensions = {"triple": 16, "double": 10, "squarefree": 10}

    reduced_data: dict[str, tuple[tuple[sp.Symbol, ...], sp.Expr, sp.Matrix]] = {}

    for name, cubic in cubic_profiles.items():
        cubic_n = split_q * sp.hessian(cubic, variables)
        gradient = sp.Matrix([sp.diff(cubic, variable) for variable in variables])
        dual_gradient = split_q * gradient
        expressions = (
            sp.trace(generic_m),
            sp.trace(generic_m * cubic_n),
            (gradient.T * generic_m * dual_gradient)[0],
        )
        rows = sp.Matrix.vstack(
            *(
                coefficient_linear_rows(expression, variables, coefficients)
                for expression in expressions
            )
        )
        nullspace = rows.nullspace()
        expected_dimension = expected_dimensions[name]
        if len(nullspace) != expected_dimension:
            raise AssertionError((name, "unexpected linear dimension", len(nullspace)))

        parameters = sp.symbols(f"{name[0]}p0:{expected_dimension}")
        coefficient_vector = sp.zeros(len(coefficients), 1)
        for parameter, basis_vector in zip(parameters, nullspace, strict=True):
            coefficient_vector += parameter * basis_vector
        reduced_e = sp.expand(
            sum(
                coefficient_vector[index]
                * x1**exponent[0]
                * x2**exponent[1]
                * y1**exponent[2]
                * y2**exponent[3]
                for index, exponent in enumerate(exponents)
            )
        )

        if name in {"double", "squarefree"}:
            for left in (y1, y2):
                for right in (y1, y2):
                    if sp.expand(sp.diff(reduced_e, left, right)) != 0:
                        raise AssertionError((name, "quadratic y-term survived"))
            h1 = sp.expand(sp.diff(reduced_e, y1))
            h2 = sp.expand(sp.diff(reduced_e, y2))
            base = sp.expand(reduced_e - y1 * h1 - y2 * h2)
            if y1 in base.free_symbols or y2 in base.free_symbols:
                raise AssertionError((name, "base quartic still uses y"))
            if sp.expand(sp.diff(h1, x1) + sp.diff(h2, x2)) != 0:
                raise AssertionError((name, "divergence identity failed"))
        else:
            if sp.expand(sp.diff(reduced_e, y1, y1)) != 0:
                raise AssertionError((name, "not linear in y1"))
            a3 = sp.expand(sp.diff(reduced_e, y1))
            b4 = sp.expand(reduced_e - y1 * a3)
            if y1 in a3.free_symbols or y1 in b4.free_symbols:
                raise AssertionError((name, "y1 elimination failed"))
            if sp.expand(sp.diff(a3, x1) + sp.diff(b4, x2, y2)) != 0:
                raise AssertionError((name, "triple harmonic identity failed"))

        reduced_data[name] = (parameters, reduced_e, cubic_n)
        print(
            "QUARTIC_HN_WARING_RANK9_LINEAR_PASS",
            f"profile={name}",
            f"rank={rows.rank()}",
            f"dimension={len(nullspace)}",
        )

    parameters, reduced_e, cubic_n = reduced_data["triple"]
    reduced_m = split_q * sp.hessian(reduced_e, variables)
    trace_three_t = sp.Poly(
        sp.expand(sp.trace((reduced_m + tau * cubic_n) ** 3)),
        tau,
    ).coeff_monomial(tau)
    coefficient_polynomial = sp.Poly(sp.expand(trace_three_t), *variables)

    product_supports: set[tuple[int, int]] = set()
    for _, coefficient in coefficient_polynomial.terms():
        polynomial = sp.Poly(sp.expand(coefficient), *parameters)
        terms = polynomial.terms()
        if len(terms) != 1:
            raise AssertionError(("triple", "nonmonomial channel equation", coefficient))
        monomial, scalar = terms[0]
        if scalar == 0 or sum(monomial) != 2:
            raise AssertionError(("triple", "unexpected channel degree", coefficient))
        indices = [index for index, exponent in enumerate(monomial) for _ in range(exponent)]
        if len(indices) != 2 or indices[0] == indices[1]:
            raise AssertionError(("triple", "unexpected square channel", coefficient))
        product_supports.add(tuple(sorted(indices)))

    expected_products = {
        tuple(sorted((left, right)))
        for left in (1, 5, 9)
        for right in (2, 6, 10)
    }
    if product_supports != expected_products:
        raise AssertionError(("triple", product_supports, expected_products))

    print("QUARTIC_HN_WARING_RANK9_TRIPLE_K33_PASS")


def main() -> None:
    patterns = enumerate_small_direction_patterns()
    if len(patterns) != 25:
        raise AssertionError(("unexpected pattern count", len(patterns)))

    certificates: list[tuple[Pattern, int, sp.Expr]] = []

    for pattern in patterns:
        matrix = square_pair_restriction_matrix(pattern)
        parameterized = pattern.directions == 4
        _, determinant = certified_minor(matrix, parameterized)

        key = (pattern.zero_columns, pattern.multiplicities)
        expected = (
            EXPECTED_FOUR_DIRECTION_CASES[key]
            if parameterized
            else EXPECTED_CONSTANT_CASES[key]
        )
        if sp.expand(determinant - expected) != 0:
            raise AssertionError((pattern, determinant, expected))

        if matrix.rank() != 28:
            raise AssertionError((pattern, "restriction map is not injective"))

        certificates.append((pattern, matrix.rows, determinant))

    exceptional = Pattern(4, (1, 1, 1, 1))
    exceptional_matrix = square_pair_restriction_matrix(exceptional)
    exceptional_kernel_dimension = 28 - exceptional_matrix.rank()
    if exceptional_kernel_dimension != 5:
        raise AssertionError(("unexpected exceptional kernel", exceptional_kernel_dimension))

    print("QUARTIC_HN_WARING_PATTERN_COUNT=25")
    for pattern, rows, determinant in certificates:
        print(
            "QUARTIC_HN_WARING_MINOR_PASS",
            f"zero={pattern.zero_columns}",
            f"mult={pattern.multiplicities}",
            f"rows={rows}",
            f"det={sp.factor(determinant)}",
        )
    print("QUARTIC_HN_WARING_EXCEPTIONAL_PATTERN=(4;(1,1,1,1))")
    print("QUARTIC_HN_WARING_EXCEPTIONAL_KERNEL_DIMENSION=5")
    print("QUARTIC_HN_WARING_RANK8_GATE_PASS")
    rank_nine_linear_normal_forms()
    print("QUARTIC_HN_WARING_RANK9_LINEAR_GATE_PASS")


if __name__ == "__main__":
    main()
