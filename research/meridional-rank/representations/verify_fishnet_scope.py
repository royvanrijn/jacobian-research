#!/usr/bin/env python3
"""Exact determinant obstruction for the 24 seven-sheet sharp-four knots.

This is a finite replay of Fox matrices and integer determinants, plus small
weighted-incidence controls.  The knot/Coxeter/fishnet implications are the
written statements of FISHNET_SCOPE_AUDIT.tex.  No other-family or global
novelty conclusion is produced.
"""

import argparse
from collections import Counter
from copy import deepcopy
from functools import reduce
import hashlib
from itertools import combinations
import json
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


def determinant(matrix):
    """Fraction-free integer determinant, with explicit exact divisions."""
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    previous, sign = 1, 1
    for column in range(size - 1):
        pivot_row = next((row for row in range(column, size) if work[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                value = work[row][target] * pivot - work[row][column] * work[column][target]
                assert value % previous == 0
                work[row][target] = value // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def maximal_minor_gcd(rows):
    rank = len(rows[0])
    assert len(rows) >= rank
    minors = [determinant([rows[i] for i in selected]) for selected in combinations(range(len(rows)), rank)]
    return reduce(gcd, map(abs, minors)), minors


def factorization(value):
    assert value > 0
    result, divisor = {}, 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result[divisor] = result.get(divisor, 0) + 1
            value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        result[value] = result.get(value, 0) + 1
    return result


def fox_matrix(pd):
    count = len(pd)
    labels = [value for crossing in pd for value in crossing]
    assert all(len(crossing) == 4 for crossing in pd)
    assert Counter(labels) == Counter({value: 2 for value in range(1, 2 * count + 1)}), "invalid strand incidence"
    parent = {value: value for value in labels}

    def root(value):
        while parent[value] != value:
            value = parent[value]
        return value

    for _, over1, _, over2 in pd:
        parent[root(over2)] = root(over1)
    components = {}
    for value in sorted(parent):
        components.setdefault(root(value), []).append(value)
    arcs = sorted(components.values())
    assert len(arcs) == count
    arc = {value: i for i, component in enumerate(arcs) for value in component}
    rows = []
    for incoming, over, outgoing, _ in pd:
        row = [0] * count
        row[arc[over]] += 2
        row[arc[incoming]] -= 1
        row[arc[outgoing]] -= 1
        assert sum(row) == 0
        rows.append(row[1:])  # Fix one normalized core-group arc to zero.
    return rows


def weighted_incidence(vertices, edges):
    rows = []
    for start, end, label in edges:
        assert 0 <= start < vertices and 0 <= end < vertices and start != end
        assert isinstance(label, int) and label >= 2
        row = [0] * vertices
        row[start] -= label
        row[end] += label
        rows.append(row[1:])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    source = HERE / "seven_sheet_cover_inputs.json"
    inputs = json.loads(source.read_text())
    positive = [case for case in inputs["cases"] if any(action["expected_H1_F3_dimension"] for action in case["actions"])]
    assert len(positive) == len({case["name"] for case in positive}) == 24
    assert determinant([[0, 2], [3, 4]]) == -6
    assert determinant([[2, 4], [1, 2]]) == 0
    tree_order, _ = maximal_minor_gcd(weighted_incidence(4, [(0, 1, 3), (1, 2, 5), (2, 3, 7)]))
    assert tree_order == 105
    cycle_order, cycle_minors = maximal_minor_gcd(weighted_incidence(4, [(0, 1, 3), (1, 2, 5), (2, 3, 7), (3, 0, 11)]))
    assert cycle_order == 1
    weight_control = [(0, 1, 3), (1, 2, 3), (2, 3, 3),
                      (0, 2, 2), (0, 3, 2), (1, 3, 2)]
    weight_control_order, _ = maximal_minor_gcd(weighted_incidence(4, weight_control))
    assert weight_control_order == 1 and sum(label for _, _, label in weight_control) == 15
    # The classical sum of three trefoils has determinant 27; its prime
    # multiplicity is three, so this obstruction must not exclude it.
    assert sum(factorization(27).values()) == 3
    results = []
    for case in positive:
        matrix = fox_matrix(case["pd"])
        value, minors = maximal_minor_gcd(matrix)
        assert len(minors) == 13 and value > 0 and value % 2 == 1
        assert value % 3 != 0
        assert set(map(abs, minors)) == {value}
        factors = factorization(value)
        assert reduce(lambda a, item: a * item[0] ** item[1], factors.items(), 1) == value
        omega = sum(factors.values())
        assert omega <= 2
        results.append({
            "name": case["name"],
            "determinant": value,
            "prime_factorization": {str(prime): exponent for prime, exponent in factors.items()},
            "prime_factors_with_multiplicity": omega,
            "maximal_minors": minors,
            "tree_coxeter_rank_upper_bound": omega + 1,
            "excluded_from_BKP_proposition_1_4_at_bridge_four": True,
        })
    bad = deepcopy(positive[0]["pd"])
    bad[0][0] = bad[0][1]
    try:
        fox_matrix(bad)
    except AssertionError as error:
        assert "strand incidence" in str(error)
    else:
        raise AssertionError("corrupt diagram incidence accepted")
    paths = [Path(__file__).resolve(), source, HERE / "FISHNET_SCOPE_AUDIT.tex"]
    report = {
        "status": "PASS",
        "positive_knots": 24,
        "BKP_proposition_1_4_exclusions": 24,
        "Lustig_Moriah_nonunit_plat_gcd_exclusions": 24,
        "maximum_prime_factor_multiplicity": max(result["prime_factors_with_multiplicity"] for result in results),
        "results": results,
        "controls": {
            "zero_pivot_integer_determinant": -6,
            "path_3_5_7_orientation_abelianization_order": tree_order,
            "cycle_3_5_7_11_maximal_minors": cycle_minors,
            "cycle_orientation_abelianization_order": cycle_order,
            "sharp_weight_sum_control": 15,
            "weight_sum_control_orientation_abelianization_order": weight_control_order,
            "three_trefoils_not_excluded_by_multiplicity": True,
            "corrupt_diagram_incidence_rejected": True,
        },
        "scope": "Exact Fox determinant and prime-factor counts on the retained 24 positive seven-sheet diagrams. The written arguments exclude every presentation satisfying BKP2026 Proposition 1.4 or the Lustig-Moriah nonunit plat gcd criterion as reported in Moriah2007 Theorem 4.20. Arbitrary fishnets, cyclic twisted graphs and global prior coverage remain outside this conclusion.",
        "input_sha256": {path.relative_to(HERE).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "results"}, indent=2))


if __name__ == "__main__":
    main()
