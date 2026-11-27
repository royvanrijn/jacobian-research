"""Check the high-weight parallel-extension obstruction and grade-104 scan bounds."""

from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from math import comb, isqrt
from pathlib import Path
import json

from verify_cubic_index_four_syzygy58 import jacobian, multiply, validate
from verify_subtree_shuffle import ShuffleQuotient


PRIME = 1_000_000_009
LITERAL_INPUT = Path(__file__).with_name("cubic_index_four_syzygy58.json")
OLD_DEGREE = 207


def high_weight_monomials(weights):
    groups = defaultdict(list)
    for monomial in combinations_with_replacement(range(len(weights)), 3):
        weight = sum(weights[index] for index in monomial)
        if weight > OLD_DEGREE:
            groups[weight].append(monomial)
    return groups


def reduce_cube_to_field(jacobian_cube):
    rows = defaultdict(list)
    for (row, column), polynomial in jacobian_cube.items():
        residues = {}
        for monomial, rational in polynomial.items():
            value = rational.numerator * pow(rational.denominator, -1, PRIME) % PRIME
            if value:
                residues[monomial] = value
        if residues:
            rows[row].append((column, residues))
    return rows


def gradient_row(monomial, cube_rows):
    """Coefficient column indexed by output column and degree-eight monomial."""
    coefficients = {}
    for index, multiplicity in Counter(monomial).items():
        derivative = list(monomial)
        derivative.remove(index)
        for column, polynomial in cube_rows.get(index, ()):
            for other, coefficient in polynomial.items():
                key = (column, tuple(sorted((*derivative, *other))))
                value = (coefficients.get(key, 0) + multiplicity * coefficient) % PRIME
                if value:
                    coefficients[key] = value
                else:
                    coefficients.pop(key, None)
    return coefficients


def check_full_column_rank(monomials, cube_rows):
    pivots = {}
    for monomial in monomials:
        coefficients = gradient_row(monomial, cube_rows)
        while coefficients:
            pivot = max(coefficients)
            leading = coefficients[pivot]
            if pivot not in pivots:
                inverse = pow(leading, -1, PRIME)
                pivots[pivot] = {
                    key: value * inverse % PRIME for key, value in coefficients.items()
                }
                break
            for key, coefficient in pivots[pivot].items():
                value = (coefficients.get(key, 0) - leading * coefficient) % PRIME
                if value:
                    coefficients[key] = value
                else:
                    coefficients.pop(key, None)
        else:
            raise AssertionError(f"Dependent column at {monomial}")
    return len(pivots)


def check_parallel_extension():
    assert all(PRIME % divisor for divisor in range(2, isqrt(PRIME) + 1))
    terms, weights = validate(json.loads(LITERAL_INPUT.read_text()))
    derivative = jacobian(terms)
    cube = multiply(multiply(derivative, derivative), derivative)
    assert cube and not multiply(cube, derivative)
    rows = reduce_cube_to_field(cube)
    groups = high_weight_monomials(weights)
    assert len(groups) == 128 and sum(map(len, groups.values())) == 3_788
    assert min(groups) == 209 and max(groups) == 621 and len(groups[209]) == 129
    for weight, monomials in sorted(groups.items()):
        assert check_full_column_rank(monomials, rows) == len(monomials), weight
    return len(groups), sum(map(len, groups.values()))


def count_recursive_scan():
    quotient = ShuffleQuotient(3, 4)
    known = [1] + [quotient.step(grade)["quotient_dimension"] for grade in range(1, 13)]
    assert known == [1, 1, 1, 2, 3, 6, 12, 25, 52, 114, 253, 571, 1304]
    # The all-grade Engel-four example gives at least one class in later grades.
    dimensions = known + [1] * (104 - len(known))
    counts = {(0, 0): 1}
    for grade, dimension in enumerate(dimensions):
        updated = defaultdict(int)
        for (size, total_grade), count in counts.items():
            for multiplicity in range(9 - size):
                new_grade = total_grade + grade * multiplicity
                if new_grade > 103:
                    break
                choices = comb(dimension + multiplicity - 1, multiplicity)
                updated[(size + multiplicity, new_grade)] += count * choices
        counts = updated
    columns = counts[(3, 103)]
    substitutions = sum(dimensions[grade] * counts.get((8, 100 - grade), 0)
                        for grade in range(101))
    assert columns == 2_829_774
    assert substitutions == 5_211_060_703_987_899_795_082_169
    return columns, substitutions


if __name__ == "__main__":
    groups, total = check_parallel_extension()
    columns, substitutions = count_recursive_scan()
    print(f"Full column rank for {total} high-weight cubics in {groups} weight spaces.")
    print(f"Grade 104 lower bounds: {columns} columns and {substitutions} root substitutions.")
