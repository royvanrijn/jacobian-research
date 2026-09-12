#!/usr/bin/env python3
"""Close the residual characteristic-zero rank-ten Gale-loop types.

This checker starts from the frozen loop census and the 103 types closed by
``verify_quartic_hn_rank10_loop_survivors.py``.  It verifies three final
mechanisms:

* one disconnected active Gale matroid splits into summands of at most four
  variables;
* on three triple-class types the unique possible self-square support avoids
  the triple, so the first two traces and the six-dimensional Witt bound give
  a contradiction; and
* four residual realization-plus-self-square saturations have unit ideal over
  QQ in every required projective chart.

The last item requires Singular.  The script reconstructs all ideals from the
frozen ranklines; no Singular output is trusted as a checked-in certificate.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from itertools import combinations
from pathlib import Path

import sympy as sp

from verify_quartic_hn_rank10_loop_survivors import (
    CENSUS_PATH,
    CENSUS_SHA256,
    EMPTY_Q_KEYS,
    SELF_SQUARE_SATURATION_CLOSED,
    file_sha256,
    has_matroidal_self_square_support,
)
from verify_quartic_hn_rank10_matroid_survivors import (
    bases_from_rankline,
    normalized_realization,
    rank_table,
)


DISCONNECTED_CLOSED = {(7, 70, "1111113")}
LOOP_TRIPLE_WITT_CLOSED = {
    (7, 0, "1111113"),
    (7, 1, "1111113"),
    (7, 2, "3111111"),
}
SINGULAR_SATURATION_CLOSED = {
    (7, 3, "1111113"),
    (8, 581, "11121111"),
    (9, 188841, "111111111"),
    (9, 188846, "111111111"),
}


def direction_mask(copies: list[int], copy_mask: int) -> int:
    result = 0
    for copy, element in enumerate(copies):
        if copy_mask & (1 << copy):
            result |= 1 << element
    return result


def possible_self_square_supports(
    n: int, bases: list[int], weights: tuple[int, ...]
) -> list[int]:
    """Return every support passing the cyclic-support/Witt necessary test."""

    ranks = rank_table(n, bases)
    copies = [
        element
        for element, multiplicity in enumerate(weights)
        for _ in range(multiplicity)
    ]
    all_copies = (1 << len(copies)) - 1
    result: list[int] = []
    for support in range(1, 1 << len(copies)):
        support_rank = ranks[direction_mask(copies, support)]
        cyclic = all(
            ranks[direction_mask(copies, support ^ (1 << copy))]
            == support_rank
            for copy in range(len(copies))
            if support & (1 << copy)
        )
        if not cyclic:
            continue
        projected_waring_rank = (
            support.bit_count()
            - 4
            + ranks[direction_mask(copies, all_copies ^ support)]
        )
        if projected_waring_rank <= support.bit_count() // 2:
            result.append(support)
    return result


def matroid_components(n: int, bases: list[int]) -> list[set[int]]:
    """Compute components by joining elements that occur in one circuit."""

    ranks = rank_table(n, bases)
    parents = list(range(n))

    def root(element: int) -> int:
        while parents[element] != element:
            parents[element] = parents[parents[element]]
            element = parents[element]
        return element

    def join(left: int, right: int) -> None:
        left_root, right_root = root(left), root(right)
        if left_root != right_root:
            parents[right_root] = left_root

    for mask in range(1, 1 << n):
        if ranks[mask] != mask.bit_count() - 1:
            continue
        if any(
            ranks[mask ^ (1 << element)] < ranks[mask]
            for element in range(n)
            if mask & (1 << element)
        ):
            continue
        elements = [element for element in range(n) if mask & (1 << element)]
        for element in elements[1:]:
            join(elements[0], element)

    components: dict[int, set[int]] = {}
    for element in range(n):
        components.setdefault(root(element), set()).add(element)
    return list(components.values())


def polynomial_rank(polynomials: list[sp.Expr], variables: tuple[sp.Symbol, ...]) -> int:
    monomials = sorted(
        {
            monomial
            for polynomial in polynomials
            for monomial, _ in sp.Poly(polynomial, *variables).terms()
        }
    )
    matrix = sp.zeros(len(monomials), len(polynomials))
    for column, polynomial in enumerate(polynomials):
        terms = dict(sp.Poly(polynomial, *variables).terms())
        for row, monomial in enumerate(monomials):
            matrix[row, column] = terms.get(monomial, 0)
    return matrix.rank()


def squarefree_open_product(
    open_minors: list[sp.Expr], variables: tuple[sp.Symbol, ...]
) -> sp.Expr:
    factors: list[sp.Expr] = []
    for minor in open_minors:
        _, factorization = sp.factor_list(minor, *variables)
        for factor, _ in factorization:
            factor = sp.Poly(factor, *variables).monic().as_expr()
            if factor not in factors:
                factors.append(factor)
    product = sp.prod(factors)
    if product == 0:
        raise AssertionError("zero realization-open product")
    return product


def self_square_data(
    gale: sp.Matrix,
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...], list[sp.Expr]]:
    kernel = sp.Matrix.hstack(*gale.nullspace())
    if kernel.shape != (9, 5) or gale * kernel != sp.zeros(4, 5):
        raise AssertionError("bad symbolic active Gale kernel")
    parameters = sp.symbols("qhnw10_t0:5")
    test_variables = sp.symbols("qhnw10_u0:5")
    active_row = kernel * sp.Matrix(parameters)
    kernel_forms = kernel * sp.Matrix(test_variables)
    polynomial = sp.Poly(
        sp.expand(
            sum(
                active_row[index] ** 2 * kernel_forms[index] ** 2
                for index in range(9)
            )
        ),
        *test_variables,
    )
    equations: list[sp.Expr] = []
    for coefficient in polynomial.coeffs():
        coefficient = sp.factor(coefficient)
        if coefficient != 0 and coefficient not in equations:
            equations.append(coefficient)
    return active_row, parameters, equations


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def assert_unit_chart(
    singular: str,
    equations: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> None:
    source = (
        f"ring R=0,({','.join(map(str, variables))}),dp;\n"
        "option(redSB);\n"
        f"ideal I={','.join(singular_expression(value) for value in equations)};\n"
        "ideal J=slimgb(I);\n"
        "print(reduce(1,J));\n"
        "print(dim(J));\n"
        "print(size(J));\n"
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=source,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(completed.stdout + completed.stderr)
    if completed.stdout.split() != ["0", "-1", "1"]:
        raise AssertionError(completed.stdout)


def main() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the final loop saturations")
    census_path = CENSUS_PATH.resolve()
    if file_sha256(census_path) != CENSUS_SHA256:
        raise AssertionError("stale loop census artifact")
    census = json.loads(census_path.read_text())

    records: dict[tuple[int, int, int], dict[str, object]] = {
        (
            record["simple_points"],
            record["catalogue_index"],
            record["loops"],
        ): record
        for record in census["records"]
    }
    one_loop_types: set[tuple[int, int, str]] = set()
    support_closed: set[tuple[int, int, str]] = set()
    schemes: dict[
        tuple[int, int],
        tuple[
            sp.Matrix,
            tuple[sp.Symbol, ...],
            list[sp.Expr],
            list[sp.Expr],
            list[int],
        ],
    ] = {}
    for record in census["records"]:
        if record["loops"] != 1:
            continue
        key = (record["simple_points"], record["catalogue_index"])
        if key in EMPTY_Q_KEYS:
            continue
        bases = bases_from_rankline(record["simple_points"], record["rankline"])
        if key not in schemes:
            matrix, variables, equations, open_minors = normalized_realization(
                record["simple_points"], bases
            )
            schemes[key] = (
                matrix,
                variables,
                equations,
                open_minors,
                bases,
            )
        for encoded in record["weight_representatives"]:
            survivor = (*key, encoded)
            one_loop_types.add(survivor)
            weights = tuple(int(value) for value in encoded)
            if not has_matroidal_self_square_support(
                record["simple_points"], bases, weights
            ):
                support_closed.add(survivor)

    if len(one_loop_types) != 45 or len(support_closed) != 16:
        raise AssertionError((len(one_loop_types), len(support_closed)))
    closure_partition = (
        support_closed
        | SELF_SQUARE_SATURATION_CLOSED
        | DISCONNECTED_CLOSED
        | LOOP_TRIPLE_WITT_CLOSED
        | SINGULAR_SATURATION_CLOSED
    )
    closure_parts = [
        support_closed,
        SELF_SQUARE_SATURATION_CLOSED,
        DISCONNECTED_CLOSED,
        LOOP_TRIPLE_WITT_CLOSED,
        SINGULAR_SATURATION_CLOSED,
    ]
    if any(
        left & right
        for left, right in combinations(closure_parts, 2)
    ):
        raise AssertionError("overlapping one-loop closure mechanisms")
    if closure_partition != one_loop_types:
        raise AssertionError(one_loop_types - closure_partition)

    # The only residual disconnected type has active component Waring
    # dimensions three and two; the Gale loop adds a one-variable summand.
    for n, catalogue_index, encoded in DISCONNECTED_CLOSED:
        matrix, _, _, _, bases = schemes[(n, catalogue_index)]
        del matrix
        components = matroid_components(n, bases)
        weights = tuple(int(value) for value in encoded)
        ranks = rank_table(n, bases)
        waring_dimensions = sorted(
            sum(weights[element] for element in component)
            - ranks[sum(1 << element for element in component)]
            for component in components
        )
        if waring_dimensions != [2, 3]:
            raise AssertionError(waring_dimensions)

    # The trace restriction to a triple class uses these independent
    # quadratic and quartic forms.  It first kills the three Gram diagonals,
    # then all three off-diagonal Gram entries.
    y = sp.symbols("y0:3")
    latitude = sum(y)
    quadratics = [latitude**2, *(value**2 for value in y)]
    quartics = [
        latitude**4,
        *(latitude**2 * value**2 for value in y),
        *(y[left] ** 2 * y[right] ** 2 for left, right in combinations(range(3), 2)),
    ]
    if polynomial_rank(quadratics, y) != 4:
        raise AssertionError("dependent triple-class trace quadratics")
    if polynomial_rank(quartics, y) != 7:
        raise AssertionError("dependent triple-class trace quartics")

    for n, catalogue_index, encoded in LOOP_TRIPLE_WITT_CLOSED:
        _, _, _, _, bases = schemes[(n, catalogue_index)]
        weights = tuple(int(value) for value in encoded)
        ranks = rank_table(n, bases)
        supports = possible_self_square_supports(n, bases, weights)
        copies = [
            element
            for element, multiplicity in enumerate(weights)
            for _ in range(multiplicity)
        ]
        triple_element = weights.index(3)
        complement = ((1 << n) - 1) ^ (1 << triple_element)
        triple_waring_rank = 3 - 4 + ranks[complement]
        if triple_waring_rank != 3:
            raise AssertionError(((n, catalogue_index, encoded), triple_waring_rank))
        triple_mask = sum(
            1 << copy
            for copy, element in enumerate(copies)
            if element == triple_element
        )
        if len(supports) != 1 or supports[0] & triple_mask:
            raise AssertionError(((n, catalogue_index, encoded), supports))

    for n, catalogue_index, encoded in sorted(SINGULAR_SATURATION_CLOSED):
        matrix, realization_variables, realization_equations, open_minors, bases = schemes[
            (n, catalogue_index)
        ]
        del bases
        weights = tuple(int(value) for value in encoded)
        gale = sp.Matrix.hstack(
            *[
                matrix[:, element]
                for element, multiplicity in enumerate(weights)
                for _ in range(multiplicity)
            ]
        )
        active_row, parameters, self_square_equations = self_square_data(gale)
        inverse = sp.Symbol("qhnw10_open_inverse")
        open_product = squarefree_open_product(open_minors, realization_variables)
        common = realization_equations + self_square_equations

        if (n, catalogue_index, encoded) == (9, 188841, "111111111"):
            supports = possible_self_square_supports(n, schemes[(n, catalogue_index)][4], weights)
            if len(supports) != 1:
                raise AssertionError(supports)
            support = supports[0]
            outside_equations = [
                active_row[index]
                for index in range(9)
                if not support & (1 << index)
            ]
            pivot = next(index for index in range(9) if support & (1 << index))
            charts = [outside_equations + [active_row[pivot] - 1]]
        else:
            charts = [[parameter - 1] for parameter in parameters]

        for chart in charts:
            assert_unit_chart(
                singular,
                common
                + chart
                + [inverse * open_product - 1],
                (*parameters, inverse, *realization_variables),
            )

    print("QHNW10_ONE_LOOP_DISCONNECTED_CLOSED=1")
    print("QHNW10_ONE_LOOP_TRIPLE_WITT_CLOSED=3")
    print("QHNW10_ONE_LOOP_SINGULAR_SATURATION_CLOSED=4")
    print("QHNW10_ONE_LOOP_QBAR_TYPES_UNIVERSALLY_CLOSED=45")
    print("QHNW10_LOOP_QBAR_TYPES_UNIVERSALLY_CLOSED=111")
    print("QHNW10_LOOP_QBAR_TYPES_REMAINING=0")
    print("QUARTIC_HN_RANK10_LOOP_CLOSURE_PASS")


if __name__ == "__main__":
    main()
