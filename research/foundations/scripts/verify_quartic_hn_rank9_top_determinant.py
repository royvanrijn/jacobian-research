#!/usr/bin/env python3
"""Exact audit for the rank-nine top-determinant obstruction.

The mathematical argument is coordinate-free and uses Cauchy--Binet plus a
rank-three matroid packing lemma.  This checker verifies the finite algebraic
identities and the only exceptional multiplicity-four partition combinatorics:

* complementary maximal minors of a Gale pair;
* the one-surviving-term determinant at a six-support kernel vector;
* every partition of five outside terms into direction classes of size at
  most two admits a 3+2 split with distinct directions on each side;
* the multiplicity-five restriction has exactly one relation on its five
  repeated terms and therefore triggers the one-relation square-pair lemma.

It is an exact regression, not a replacement for the matroid base-packing
argument in the note.
"""

from __future__ import annotations

from itertools import combinations, product

import sympy as sp


def complementary_minor_identity() -> None:
    # A generic rational 3 x 9 Gale matrix, and a canonical 6 x 9 matrix with
    # the same kernel.  Complementary minors must be proportional.
    K = sp.Matrix(
        [
            [1, 0, 0, 1, 1, 1, 2, 3, 5],
            [0, 1, 0, 1, 2, 3, 5, 7, 11],
            [0, 0, 1, 1, 3, 5, 7, 11, 13],
        ]
    )
    kernel = K.nullspace()
    V = sp.Matrix.hstack(*kernel).T  # 6 x 9, with ker(V)=row(K)^T
    if V.rank() != 6 or K.rank() != 3 or V * K.T != sp.zeros(6, 3):
        raise AssertionError("invalid Gale pair")

    ratios = []
    for T in combinations(range(9), 3):
        S = tuple(i for i in range(9) if i not in T)
        kt = K[:, T].det()
        vs = V[:, S].det()
        if kt:
            ratios.append(sp.factor((vs / kt) ** 2))
        elif vs != 0:
            raise AssertionError((T, kt, vs))
    if not ratios or any(sp.simplify(ratio - ratios[0]) != 0 for ratio in ratios):
        raise AssertionError(ratios)
    print("QUARTIC_HN_RANK9_COMPLEMENTARY_MINOR_PASS", ratios[0])


def one_surviving_term_identity() -> None:
    # Choose T={0,1,2}; the complement S has a full-support relation.
    K = sp.Matrix(
        [
            [1, 0, 0, 1, 1, 1, 2, 3, 5],
            [0, 1, 0, 1, 2, 3, 5, 7, 11],
            [0, 0, 1, 1, 3, 5, 7, 11, 13],
        ]
    )
    V = sp.Matrix.hstack(*K.nullspace()).T
    T = (0, 1, 2)
    S = tuple(i for i in range(9) if i not in T)
    KS = K[:, S]
    relations = KS.nullspace()
    trial = None
    # A finite exact search for a full-support vector in the 3-dimensional
    # relation space.  Over an infinite field its existence follows from the
    # no-coloop condition.
    for coefficients in product(range(-3, 4), repeat=len(relations)):
        if coefficients == (0,) * len(relations):
            continue
        candidate = sum(
            (coefficient * relation for coefficient, relation in zip(coefficients, relations, strict=True)),
            sp.zeros(6, 1),
        )
        if all(entry != 0 for entry in candidate):
            trial = candidate
            break
    if trial is None:
        raise AssertionError("no full-support relation found")

    l = sp.zeros(9, 1)
    for index, value in zip(S, trial, strict=True):
        l[index] = value
    if K * l != sp.zeros(3, 1):
        raise AssertionError("not a Gale relation")

    D = sp.diag(*(entry**2 for entry in l))
    hessian_core = V * D * V.T
    expected = V[:, S].det() ** 2 * sp.prod(l[index] ** 2 for index in S)
    actual = sp.factor(hessian_core.det())
    if sp.expand(actual - expected) != 0 or actual == 0:
        raise AssertionError((actual, expected))
    print("QUARTIC_HN_RANK9_ONE_TERM_CAUCHY_BINET_PASS", actual)


def partitions_of_five(maximum: int = 2) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []

    def rec(remaining: int, upper: int, prefix: tuple[int, ...]) -> None:
        if remaining == 0:
            result.append(prefix)
            return
        for part in range(min(upper, remaining, maximum), 0, -1):
            rec(remaining - part, part, prefix + (part,))

    rec(5, maximum, ())
    return result


def has_distinct_direction_split(parts: tuple[int, ...]) -> bool:
    labels = []
    for direction, multiplicity in enumerate(parts):
        labels.extend([direction] * multiplicity)
    for selected in combinations(range(5), 3):
        selected_set = set(selected)
        left = [labels[i] for i in selected]
        right = [labels[i] for i in range(5) if i not in selected_set]
        if len(set(left)) == 3 and len(set(right)) == 2:
            return True
    return False


def multiplicity_four_partition_gate() -> None:
    patterns = partitions_of_five()
    expected = {(2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)}
    if set(patterns) != expected:
        raise AssertionError(patterns)
    for pattern in patterns:
        if not has_distinct_direction_split(pattern):
            raise AssertionError(pattern)
        print("QUARTIC_HN_RANK9_MULT4_DIRECTION_SPLIT_PASS", pattern)


def multiplicity_five_one_relation_gate() -> None:
    # Five proportional Gale columns and four columns in distinct directions
    # through the repeated point.  Restricting away the four complementary
    # Waring covectors leaves five values with one full-support relation.
    K = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 1, 0, 1, 2],
            [0, 0, 0, 0, 0, 0, 1, 1, 3],
        ]
    )
    if K.rank() != 3:
        raise AssertionError("rank")
    repeated_restriction = K[:, :5]
    if repeated_restriction.rank() != 1 or len(repeated_restriction.nullspace()) != 4:
        raise AssertionError("repeated restriction")
    # The coefficient vectors of relations, restricted to the five repeated
    # entries, span one dimension and have full support.
    row_restrictions = [sp.Matrix(K[row, :5]).T for row in range(3)]
    span = sp.Matrix.hstack(*row_restrictions)
    if span.rank() != 1:
        raise AssertionError(span)
    generator = next(vector for vector in row_restrictions if vector != sp.zeros(5, 1))
    if any(entry == 0 for entry in generator):
        raise AssertionError(generator)
    print("QUARTIC_HN_RANK9_MULT5_ONE_RELATION_PASS", tuple(generator))


def main() -> None:
    complementary_minor_identity()
    one_surviving_term_identity()
    multiplicity_four_partition_gate()
    multiplicity_five_one_relation_gate()
    print("QUARTIC_HN_WARING_RANK9_TOP_DETERMINANT_GATE_PASS")


if __name__ == "__main__":
    main()
