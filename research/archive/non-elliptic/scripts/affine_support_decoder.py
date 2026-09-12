#!/usr/bin/env python3
"""Reusable exact decoder for sparse affine systems over F_2.

Rows are pairs ``(left, right)``.  ``left`` is a Python integer whose set
bits encode the variables occurring in the equation, and ``right`` is zero
or one.  The module provides exact rank/consistency, one solution, incidence
components, and componentwise minimum-Hamming-weight decoding.
"""

from __future__ import annotations

from collections import defaultdict
from functools import reduce

import z3


AffineRow = tuple[int, int]
AffineComponent = tuple[str, list[int], list[AffineRow]]


def _xor(terms: list[z3.BoolRef]) -> z3.BoolRef:
    if not terms:
        return z3.BoolVal(False)
    return reduce(z3.Xor, terms)


def row_echelon(rows: list[AffineRow]) -> dict[int, AffineRow] | None:
    """Return highest-pivot echelon rows, or ``None`` if inconsistent."""

    pivots: dict[int, AffineRow] = {}
    for left, right in rows:
        assert right in (0, 1)
        while left:
            pivot = left.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (left, right)
                break
            pivot_left, pivot_right = pivots[pivot]
            left ^= pivot_left
            right ^= pivot_right
        else:
            if right:
                return None
    return pivots


def rank_affine_rows(rows: list[AffineRow]) -> int:
    """Return the GF(2) rank and assert consistency of augmented rows."""

    pivots = row_echelon(rows)
    assert pivots is not None, "affine system is inconsistent"
    return len(pivots)


def solve_affine_rows(rows: list[AffineRow]) -> set[int] | None:
    """Return one solution with all free variables zero, or ``None``."""

    pivots = row_echelon(rows)
    if pivots is None:
        return None
    solution_mask = 0
    for pivot in sorted(pivots):
        left, right = pivots[pivot]
        lower = left ^ (1 << pivot)
        value = right ^ ((lower & solution_mask).bit_count() % 2)
        if value:
            solution_mask |= 1 << pivot
    return {
        index
        for index in range(solution_mask.bit_length())
        if solution_mask & (1 << index)
    }


def affine_inconsistency_certificate(rows: list[AffineRow]) -> set[int] | None:
    """Return row indices whose XOR is ``0 = 1``, or ``None`` if consistent.

    The returned set is a directly replayable dual certificate: XORing its
    left masks gives zero and XORing its right sides gives one.
    """

    pivots: dict[int, tuple[int, int, int]] = {}
    for row_index, (left, right) in enumerate(rows):
        combination = 1 << row_index
        while left:
            pivot = left.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (left, right, combination)
                break
            pivot_left, pivot_right, pivot_combination = pivots[pivot]
            left ^= pivot_left
            right ^= pivot_right
            combination ^= pivot_combination
        else:
            if right:
                certificate = {
                    index
                    for index in range(combination.bit_length())
                    if combination & (1 << index)
                }
                assert certificate
                assert not _xor_integer_rows(rows, certificate)[0]
                assert _xor_integer_rows(rows, certificate)[1] == 1
                return certificate
    return None


def _xor_integer_rows(
    rows: list[AffineRow], indices: set[int]
) -> AffineRow:
    left = 0
    right = 0
    for index in indices:
        row_left, row_right = rows[index]
        left ^= row_left
        right ^= row_right
    return left, right


def affine_components(
    rows: list[AffineRow],
    variable_layers: list[str],
) -> list[AffineComponent]:
    """Split rows by the connected variable-equation incidence graph."""

    parent = list(range(len(variable_layers)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    active: set[int] = set()
    for left, right in rows:
        if not left:
            assert right == 0, "affine system is inconsistent"
            continue
        support: list[int] = []
        remaining = left
        while remaining:
            bit = remaining & -remaining
            support.append(bit.bit_length() - 1)
            remaining ^= bit
        assert support
        active.update(support)
        for index in support[1:]:
            union(support[0], index)

    variables_by_root: dict[int, list[int]] = defaultdict(list)
    for index in active:
        variables_by_root[find(index)].append(index)
    rows_by_root: dict[int, list[AffineRow]] = defaultdict(list)
    for row in rows:
        left, _ = row
        if not left:
            continue
        first_index = (left & -left).bit_length() - 1
        rows_by_root[find(first_index)].append(row)

    components: list[AffineComponent] = []
    for root, indices in variables_by_root.items():
        layers = {variable_layers[index] for index in indices}
        assert len(layers) == 1
        components.append((layers.pop(), sorted(indices), rows_by_root[root]))
    for layer in set(variable_layers):
        components.extend(
            (layer, [index], [])
            for index, variable_layer in enumerate(variable_layers)
            if variable_layer == layer and index not in active
        )
    components.sort(key=lambda component: (component[0], component[1][0]))
    return components


def affine_component_statistics(
    components: list[AffineComponent],
) -> dict[str, list[tuple[int, int, int, int]]]:
    """Return (variables,equations,rank,nullity) for each component."""

    statistics: dict[str, list[tuple[int, int, int, int]]] = defaultdict(list)
    for layer, indices, rows in components:
        rank = rank_affine_rows(rows)
        statistics[layer].append(
            (len(indices), len(rows), rank, len(indices) - rank)
        )
    for layer in statistics:
        statistics[layer].sort(reverse=True)
    return dict(statistics)


def minimize_affine_components(
    components: list[AffineComponent],
    layer: str,
    known_true_indices: set[int],
    timeout_ms: int,
) -> tuple[set[int], int, int, bool, list[tuple[int, int, int, int]]]:
    """Minimize support componentwise, exactly when all SAT calls finish."""

    selected_indices: set[int] = set()
    total_lower = 0
    total_upper = 0
    all_exact = True
    records: list[tuple[int, int, int, int]] = []

    for component_number, (component_layer, indices, rows) in enumerate(
        components
    ):
        if component_layer != layer:
            continue
        known_component = known_true_indices.intersection(indices)
        known_mask = sum(1 << index for index in known_component)
        known_is_feasible = all(
            (left & known_mask).bit_count() % 2 == right
            for left, right in rows
        )
        if not rows:
            records.append((len(indices), 0, 0, 0))
            continue

        local_variables = {
            index: z3.Bool(f"component_{component_number}_{index}")
            for index in indices
        }
        solver = z3.Then(
            z3.Tactic("simplify"),
            z3.Tactic("bit-blast"),
            z3.Tactic("sat"),
        ).solver()
        solver.set(timeout=timeout_ms)
        for left, right in rows:
            terms = [
                variable
                for index, variable in local_variables.items()
                if left & (1 << index)
            ]
            solver.add(_xor(terms) == z3.BoolVal(bool(right)))

        lower = 0
        if known_is_feasible:
            best = set(known_component)
        else:
            result = solver.check()
            assert result == z3.sat, result
            model = solver.model()
            best = {
                index
                for index, variable in local_variables.items()
                if z3.is_true(model.eval(variable, model_completion=True))
            }
        upper = len(best)
        while lower < upper:
            bound = (lower + upper) // 2
            solver.push()
            solver.add(
                z3.PbLe(
                    [(variable, 1) for variable in local_variables.values()],
                    bound,
                )
            )
            result = solver.check()
            if result == z3.sat:
                model = solver.model()
                best = {
                    index
                    for index, variable in local_variables.items()
                    if z3.is_true(model.eval(variable, model_completion=True))
                }
                upper = len(best)
            elif result == z3.unsat:
                lower = bound + 1
            else:
                all_exact = False
                solver.pop()
                break
            solver.pop()

        selected_indices.update(best)
        total_lower += lower
        total_upper += upper
        records.append((len(indices), len(rows), lower, upper))

    return selected_indices, total_lower, total_upper, all_exact, records
