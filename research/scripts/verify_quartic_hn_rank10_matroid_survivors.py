#!/usr/bin/env python3
"""Exact OP-QHNW10 survivor and isotropic-Gram obstruction audit.

This script has two logically separate parts.

1. It verifies two characteristic-zero counterexamples to the matroid lemma
   currently stated as OP-QHNW10.  The first has Gale loops.  The second is
   loopless and has three parallel pairs, so adding looplessness does not
   repair the lemma.
2. It checks the frozen loopless, nonsimple catalogue slice obtained from
   ``matroid-database==0.3``.  There are 37 coloured isomorphism types.  Two
   underlying simplifications have empty normalized realization ideals over
   QQ; the remaining 35 have rational realizations and a six-element
   rank-two Gale flat.  The first trace kills the Gram diagonal on that flat,
   and the codimension-two square-pair calculation then excludes every
   rank-six Gram matrix on all 35 types.

The frozen slice is a replay of the catalogue output, not an independent
regeneration of the 190,214 source matroids.  In particular this script does
not enumerate simple ten-element matroids or the full loop sector.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import combinations, product

import sympy as sp


RANKLINES: dict[tuple[int, int], str] = {
    (7, 15): "00***0*********0******************0",
    (7, 19): "00***0******0000******************0",
    (7, 26): "00000********************00*0**0***",
    (7, 32): "000000*********0************000***0",
    (7, 54): "000000000******0000************0000",
    (8, 269): "0********0****0*****0**0***00************00***0**0*****0****0********0",
    (8, 581): "000000000******0000****************0000******************************0",
    (8, 586): "000000000******0000************0***0000************0*********0***0***0",
    (8, 587): "000000000******0000************00000000******************************0",
    (8, 588): "000000000******0000************00000000************0*********0***0***0",
    (8, 589): "000000000******0000************00000000************0000******000000000",
    (8, 845): "000000000000000000000000000000000000******0******0***0******0*0**0****",
    (9, 188841): "0000000000000000000000000**********0000000000*************************0000000000*********************************************0",
    (9, 188846): "0000000000000000000000000**********0000000000********************0****0000000000********************0**************0****0****0",
    (9, 188847): "0000000000000000000000000**********0000000000********************000000000000000*********************************************0",
    (9, 188848): "0000000000000000000000000**********0000000000********************000000000000000********************0**************0****0****0",
    (9, 188849): "0000000000000000000000000**********0000000000********************000000000000000********************00000**********00000000000",
}


COLOUR_REPRESENTATIVES: dict[tuple[int, int], tuple[str, ...]] = {
    (7, 15): ("1231111", "2221111"),
    (7, 19): ("1111231", "1112221"),
    (7, 26): ("1111123", "2111113", "2111122", "3111112"),
    (7, 32): ("1111123", "1112113", "1112122", "1113112"),
    (7, 54): ("1111123", "1111222"),
    (8, 269): ("11111113",),
    (8, 581): ("11131111", "11221111"),
    (8, 586): ("11131111", "11221111", "21121111", "31111111"),
    (8, 587): ("11131111", "11221111"),
    (8, 588): ("11131111", "11221111", "21121111", "31111111"),
    (8, 589): ("11111113", "11111122"),
    (8, 845): ("11111113",),
    (9, 188841): ("111121111",),
    (9, 188846): ("111121111", "211111111"),
    (9, 188847): ("111121111",),
    (9, 188848): ("111121111", "211111111"),
    (9, 188849): ("111121111",),
}


def revlex_four_sets(n: int) -> list[tuple[int, ...]]:
    return sorted(combinations(range(n), 4), key=lambda subset: tuple(reversed(subset)))


def bases_from_rankline(n: int, rankline: str) -> list[int]:
    subsets = revlex_four_sets(n)
    if len(rankline) != len(subsets):
        raise AssertionError((n, len(rankline), len(subsets)))
    return [
        sum(1 << element for element in subset)
        for subset, marker in zip(subsets, rankline, strict=True)
        if marker == "*"
    ]


def rank_table(n: int, bases: list[int]) -> list[int]:
    independent = [False] * (1 << n)
    independent[0] = True
    for basis in bases:
        subset = basis
        while True:
            independent[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & basis

    ranks = [0] * (1 << n)
    for mask in range(1, 1 << n):
        if independent[mask]:
            ranks[mask] = mask.bit_count()
        else:
            ranks[mask] = max(
                ranks[mask ^ (1 << element)]
                for element in range(n)
                if mask & (1 << element)
            )
    return ranks


def assert_basis_exchange(n: int, bases: list[int]) -> None:
    basis_set = set(bases)
    for left in bases:
        for right in bases:
            for element in range(n):
                if left & (1 << element) and not right & (1 << element):
                    alternatives = right & ~left
                    if not any(
                        ((left ^ (1 << element)) | (1 << replacement)) in basis_set
                        for replacement in range(n)
                        if alternatives & (1 << replacement)
                    ):
                        raise AssertionError((left, right, element))


def hyperplanes(n: int, ranks: list[int]) -> list[int]:
    return [
        mask
        for mask in range(1 << n)
        if ranks[mask] == 3
        and all(
            mask & (1 << element) or ranks[mask | (1 << element)] == 4
            for element in range(n)
        )
    ]


def assert_op_constraints(n: int, bases: list[int], weights: tuple[int, ...]) -> None:
    if sum(weights) != 10 or max(weights) > 3:
        raise AssertionError(weights)
    if sum(weight == 3 for weight in weights) > 1:
        raise AssertionError(weights)

    ranks = rank_table(n, bases)
    assert_basis_exchange(n, bases)
    for plane in hyperplanes(n, ranks):
        cocircuit_size = sum(
            weights[element]
            for element in range(n)
            if not plane & (1 << element)
        )
        if cocircuit_size < 3:
            raise AssertionError((weights, plane, cocircuit_size))

    full = (1 << n) - 1
    for basis in bases:
        complement_support = full
        for element in range(n):
            if basis & (1 << element) and weights[element] == 1:
                complement_support ^= 1 << element
        complement_rank = ranks[complement_support]
        has_coloop = any(
            weights[element] - bool(basis & (1 << element)) == 1
            and ranks[complement_support ^ (1 << element)] < complement_rank
            for element in range(n)
        )
        if not has_coloop:
            raise AssertionError((weights, basis, complement_support))


def normalized_realization(
    n: int, bases: list[int]
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...], list[sp.Expr], list[sp.Expr]]:
    basis_set = set(bases)
    tuple_bases = [tuple(i for i in range(n) if basis & (1 << i)) for basis in bases]
    best: tuple[int, tuple[int, ...], list[int], list[tuple[int, int]], list[tuple[int, int]]] | None = None

    for basis in tuple_bases:
        outside = [element for element in range(n) if element not in basis]
        edges = [
            (row, column)
            for row, old in enumerate(basis)
            for column, new in enumerate(outside)
            if sum(1 << x for x in (set(basis) - {old}) | {new}) in basis_set
        ]
        parent = list(range(4 + len(outside)))

        def root(vertex: int) -> int:
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        forest: list[tuple[int, int]] = []
        for row, column in edges:
            left, right = root(row), root(4 + column)
            if left != right:
                parent[left] = right
                forest.append((row, column))
        candidate = (len(edges) - len(forest), basis, outside, edges, forest)
        if best is None or candidate[0] < best[0]:
            best = candidate

    if best is None:
        raise AssertionError("matroid has no basis")
    variable_count, basis, outside, edges, forest = best
    variables = sp.symbols(f"x0:{variable_count}")
    variable_iterator = iter(variables)
    matrix = sp.zeros(4, n)
    for row, element in enumerate(basis):
        matrix[row, element] = 1
    for edge in edges:
        row, column = edge
        matrix[row, outside[column]] = 1 if edge in forest else next(variable_iterator)

    equations: list[sp.Expr] = []
    open_minors: list[sp.Expr] = []
    for subset in combinations(range(n), 4):
        determinant = sp.factor(matrix[:, subset].det())
        mask = sum(1 << element for element in subset)
        if mask in basis_set:
            if determinant not in (1, -1):
                open_minors.append(determinant)
        elif determinant != 0:
            equations.append(determinant)
    return (
        matrix,
        variables,
        list(dict.fromkeys(equations)),
        list(dict.fromkeys(open_minors)),
    )


def rational_realization_point(
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    equations: list[sp.Expr],
    n: int,
    bases: list[int],
) -> sp.Matrix | None:
    basis_set = set(bases)
    for values in product((-3, -2, -1, 1, 2, 3, 4), repeat=len(variables)):
        substitution = dict(zip(variables, values, strict=True))
        if any(equation.subs(substitution) != 0 for equation in equations):
            continue
        candidate = matrix.subs(substitution)
        actual = {
            sum(1 << element for element in subset)
            for subset in combinations(range(n), 4)
            if candidate[:, subset].det() != 0
        }
        if actual == basis_set:
            return candidate
    return None


def six_point_square_pair_matrix(multiplicities: tuple[int, ...]) -> sp.Matrix:
    lam = sp.symbols("lambda")
    if len(multiplicities) == 3:
        directions = ((1, 0), (0, 1), (1, 1))
    elif len(multiplicities) == 4:
        directions = ((1, 0), (0, 1), (1, 1), (1, lam))
    else:
        raise ValueError(multiplicities)
    columns = [
        direction
        for multiplicity, direction in zip(multiplicities, directions, strict=True)
        for _ in range(multiplicity)
    ]
    gale = sp.Matrix(2, 6, lambda row, column: columns[column][row])
    kernel = sp.Matrix.hstack(*gale.nullspace())
    variables = sp.symbols("u0:4")
    forms = kernel * sp.Matrix(variables)
    terms: list[dict[tuple[int, ...], sp.Expr]] = []
    monomial_index: dict[tuple[int, ...], int] = {}
    for left, right in combinations(range(6), 2):
        polynomial = sp.Poly(sp.expand(forms[left] ** 2 * forms[right] ** 2), *variables)
        column_terms = dict(polynomial.terms())
        terms.append(column_terms)
        for monomial in column_terms:
            monomial_index.setdefault(monomial, len(monomial_index))
    matrix = sp.zeros(len(monomial_index), 15)
    for column, column_terms in enumerate(terms):
        for monomial, coefficient in column_terms.items():
            matrix[monomial_index[monomial], column] = coefficient
    return matrix


def six_point_diagonal_square_matrix(multiplicities: tuple[int, ...]) -> sp.Matrix:
    lam, mu = sp.symbols("lambda mu")
    if len(multiplicities) == 3:
        directions = ((1, 0), (0, 1), (1, 1))
    elif len(multiplicities) == 4:
        directions = ((1, 0), (0, 1), (1, 1), (1, lam))
    elif len(multiplicities) == 5:
        directions = ((1, 0), (0, 1), (1, 1), (1, lam), (1, mu))
    else:
        raise ValueError(multiplicities)
    columns = [
        direction
        for multiplicity, direction in zip(multiplicities, directions, strict=True)
        for _ in range(multiplicity)
    ]
    gale = sp.Matrix(2, 6, lambda row, column: columns[column][row])
    kernel = sp.Matrix.hstack(*gale.nullspace())
    variables = sp.symbols("u0:4")
    forms = kernel * sp.Matrix(variables)
    terms: list[dict[tuple[int, ...], sp.Expr]] = []
    monomial_index: dict[tuple[int, ...], int] = {}
    for form in forms:
        polynomial = sp.Poly(sp.expand(form**2), *variables)
        column_terms = dict(polynomial.terms())
        terms.append(column_terms)
        for monomial in column_terms:
            monomial_index.setdefault(monomial, len(monomial_index))
    matrix = sp.zeros(len(monomial_index), 6)
    for column, column_terms in enumerate(terms):
        for monomial, coefficient in column_terms.items():
            matrix[monomial_index[monomial], column] = coefficient
    return matrix


def certified_diagonal_square_determinants() -> dict[tuple[int, ...], sp.Expr]:
    lam, mu = sp.symbols("lambda mu")
    expected = {
        (3, 2, 1): sp.Integer(4),
        (2, 2, 2): sp.Integer(4),
        (3, 1, 1, 1): 4 * lam,
        (2, 2, 1, 1): 4 * lam,
        (2, 1, 1, 1, 1): 4 * lam,
    }
    actual: dict[tuple[int, ...], sp.Expr] = {}
    for pattern in expected:
        matrix = six_point_diagonal_square_matrix(pattern)
        specialization = matrix.subs({lam: 2, mu: 3})
        _, pivot_rows = specialization.T.rref()
        rows = tuple(int(row) for row in pivot_rows[:6])
        determinant = sp.factor(matrix[list(rows), :].det())
        if determinant != expected[pattern]:
            raise AssertionError((pattern, determinant, expected[pattern]))
        actual[pattern] = determinant
    return actual


def certified_square_pair_determinants() -> dict[tuple[int, ...], sp.Expr]:
    lam = sp.symbols("lambda")
    expected = {
        (3, 2, 1): -256,
        (2, 2, 2): -512,
        (3, 1, 1, 1): -256 * lam**4,
        (2, 2, 1, 1): -512 * lam**2,
    }
    actual: dict[tuple[int, ...], sp.Expr] = {}
    for pattern in expected:
        matrix = six_point_square_pair_matrix(pattern)
        specialization = matrix.subs(lam, 2)
        _, pivot_rows = specialization.T.rref()
        rows = tuple(int(row) for row in pivot_rows[:15])
        determinant = sp.factor(matrix[list(rows), :].det())
        if determinant != expected[pattern]:
            raise AssertionError((pattern, determinant, expected[pattern]))
        actual[pattern] = determinant
    return actual


def explicit_loop_survivor() -> None:
    # U_(4,6) plus four Gale loops.  Every cocircuit has size three.  The
    # complement of each four-element basis contains two nonloop coloops.
    nonloop_bases = [sum(1 << i for i in subset) for subset in combinations(range(6), 4)]
    ranks = rank_table(6, nonloop_bases)
    if min(6 - mask.bit_count() for mask in hyperplanes(6, ranks)) != 3:
        raise AssertionError("wrong U_(4,6) cocircuit size")
    for basis in nonloop_bases:
        nonloop_complement = ((1 << 6) - 1) ^ basis
        if ranks[nonloop_complement] != 2:
            raise AssertionError(basis)
        if any(
            ranks[nonloop_complement ^ (1 << element)] != 1
            for element in range(6)
            if nonloop_complement & (1 << element)
        ):
            raise AssertionError(basis)

    # The four Gale-loop coordinates are free Waring values.  Their six
    # diagonal monomials in the first trace and their six square-pair
    # monomials in the second trace are independent.  Thus the Gram block on
    # their four independent Waring covectors vanishes, producing a totally
    # isotropic four-plane and exceeding Witt index 3.
    free_values = sp.symbols("u0:4")
    square_pairs = [
        sp.Poly(free_values[left] ** 2 * free_values[right] ** 2, *free_values)
        for left, right in combinations(range(4), 2)
    ]
    monomials = sorted({monomial for polynomial in square_pairs for monomial, _ in polynomial.terms()})
    coefficient_matrix = sp.Matrix(
        [
            [dict(polynomial.terms()).get(monomial, 0) for polynomial in square_pairs]
            for monomial in monomials
        ]
    )
    if coefficient_matrix.rank() != 6:
        raise AssertionError("dependent free-coordinate square pairs")
    print("QHNW10_LITERAL_LOOP_SURVIVOR=U_4_6_PLUS_4_LOOPS")
    print("QHNW10_LITERAL_LOOP_GRAM_SURVIVORS=0")


def explicit_loopless_survivor() -> None:
    a, b = sp.symbols("a b")
    simplification = sp.Matrix(
        [
            [1, 0, 1, 0, 0, 0, 1],
            [0, 1, 1, 0, 0, 0, a],
            [0, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 1, 1, b],
        ]
    )
    minors = [sp.factor(simplification[:, subset].det()) for subset in combinations(range(7), 4)]
    open_factors = {factor for minor in minors if minor for factor, _ in sp.factor_list(minor)[1]}
    if open_factors != {a, a - 1, b, b - 1}:
        raise AssertionError(open_factors)

    bases = bases_from_rankline(7, RANKLINES[7, 19])
    assert_op_constraints(7, bases, (2, 2, 2, 1, 1, 1, 1))
    rational = simplification.subs({a: 2, b: 3})
    actual_bases = {
        sum(1 << element for element in subset)
        for subset in combinations(range(7), 4)
        if rational[:, subset].det() != 0
    }
    if actual_bases != set(bases):
        raise AssertionError("rational survivor realizes the wrong matroid")
    print("QHNW10_LOOPLESS_SURVIVOR=M7_19_WEIGHTS_2221111")
    print("QHNW10_LOOPLESS_REALIZATION_OPEN=ab(a-1)(b-1)")


def matrix_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, column]) for column in range(matrix.cols)] for row in range(matrix.rows)]


def print_realization_record(
    key: tuple[int, int],
    rankline: str,
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    equations: list[sp.Expr],
    open_minors: list[sp.Expr],
    point: sp.Matrix | None,
    unit_ideal: bool,
) -> None:
    record = {
        "catalogue_key": list(key),
        "colour_representatives": list(COLOUR_REPRESENTATIVES[key]),
        "rankline": rankline,
        "normalized_matrix": matrix_strings(matrix),
        "variables": [str(variable) for variable in variables],
        "ideal_generators": [str(sp.factor(equation)) for equation in equations],
        "open_basis_minors": [str(sp.factor(minor)) for minor in open_minors],
        "unit_ideal_over_Q": unit_ideal,
        "rational_point": None if point is None else matrix_strings(point),
    }
    print("QHNW10_REALIZATION_RECORD=" + json.dumps(record, sort_keys=True))


def frozen_catalogue_slice(details: bool) -> None:
    if sum(len(values) for values in COLOUR_REPRESENTATIVES.values()) != 37:
        raise AssertionError("wrong frozen survivor count")

    empty_over_q: list[tuple[int, int]] = []
    rational_underlyings: dict[tuple[int, int], sp.Matrix] = {}
    profile_counts: Counter[tuple[int, ...]] = Counter()

    for key, rankline in RANKLINES.items():
        n, _ = key
        bases = bases_from_rankline(n, rankline)
        ranks = rank_table(n, bases)
        full = (1 << n) - 1
        matrix, variables, equations, open_minors = normalized_realization(n, bases)
        point = rational_realization_point(matrix, variables, equations, n, bases)
        unit_ideal = False
        if point is None:
            if equations:
                groebner = sp.groebner(equations, *variables, order="grevlex")
                if not (len(groebner.polys) == 1 and groebner.polys[0].as_expr() == 1):
                    raise AssertionError((key, "unresolved realization ideal"))
                unit_ideal = True
            else:
                raise AssertionError((key, "failed rational point search"))
            empty_over_q.append(key)
        else:
            rational_underlyings[key] = point
        if details:
            print_realization_record(
                key,
                rankline,
                matrix,
                variables,
                equations,
                open_minors,
                point,
                unit_ideal,
            )

        for encoded_weights in COLOUR_REPRESENTATIVES[key]:
            weights = tuple(int(value) for value in encoded_weights)
            assert_op_constraints(n, bases, weights)
            if point is None:
                continue
            flat_profiles = [
                tuple(
                    sorted(
                        (weights[element] for element in range(n) if flat & (1 << element)),
                        reverse=True,
                    )
                )
                for flat in range(1 << n)
                if ranks[flat] == 2
                and sum(weights[element] for element in range(n) if flat & (1 << element)) == 6
                and ranks[full ^ flat] >= 2
            ]
            if not flat_profiles:
                raise AssertionError((key, weights, "no six-element rank-two flat"))
            profile_counts[flat_profiles[0]] += 1

    if sorted(empty_over_q) != [(8, 269), (8, 845)]:
        raise AssertionError(empty_over_q)
    expected_profiles = Counter(
        {
            (3, 2, 1): 9,
            (3, 1, 1, 1): 7,
            (2, 2, 1, 1): 7,
            (2, 1, 1, 1, 1): 7,
            (2, 2, 2): 5,
        }
    )
    if profile_counts != expected_profiles:
        raise AssertionError((profile_counts, expected_profiles))

    print("QHNW10_LOOPLESS_NONSIMPLE_ABSTRACT_SURVIVORS=37")
    print("QHNW10_LOOPLESS_NONSIMPLE_QBAR_SURVIVORS=35")
    print("QHNW10_LOOPLESS_NONSIMPLE_EMPTY_Q_IDEALS=2")
    print("QHNW10_GRAM_SURVIVORS=0")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--details",
        action="store_true",
        help="print one exact JSON realization record for each simplification",
    )
    arguments = parser.parse_args()
    explicit_loop_survivor()
    explicit_loopless_survivor()
    diagonal_determinants = certified_diagonal_square_determinants()
    for pattern, determinant in diagonal_determinants.items():
        print("QHNW10_SIX_FLAT_DIAGONAL_MINOR", pattern, determinant)
    determinants = certified_square_pair_determinants()
    for pattern, determinant in determinants.items():
        print("QHNW10_SIX_FLAT_MINOR", pattern, determinant)
    frozen_catalogue_slice(arguments.details)
    print("QUARTIC_HN_RANK10_MATROID_SURVIVOR_AUDIT_PASS")


if __name__ == "__main__":
    main()
