#!/usr/bin/env python3
"""Exclude new degree-at-most-five invariant slices of the stored 21D map.

Let V=X+H be the cubic-homogeneous map in the essential BCW artifact.  A
polynomial identity-output coordinate phi must satisfy phi(V)=phi.  Good-prime
rank certificates prove over Q that the fixed space through degree five is
exactly

    span(1, X_20, ..., X_20^5).

The low-degree Lie kernels contain the near-invariant
Q=X_18*X_20-X_6*X_8 and its products with X_20 and Q, but the full pullback
defect is a single nonzero monomial. Thus no invariant through degree five
exposes a second identity slice independent of X_20.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "essential_bcw_21_counterexample.json"
)
PRIME = 1_000_003

Exponent = tuple[int, ...]
SparsePoly = dict[Exponent, int]
RationalSparsePoly = dict[Exponent, Fraction]


def rational_mod(value: str, prime: int) -> int:
    number = Fraction(value)
    return number.numerator * pow(number.denominator, -1, prime) % prime


def add_term(poly: SparsePoly, monomial: Exponent, coefficient: int) -> None:
    coefficient %= PRIME
    if not coefficient:
        return
    updated = (poly.get(monomial, 0) + coefficient) % PRIME
    if updated:
        poly[monomial] = updated
    else:
        poly.pop(monomial, None)


def multiply(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    answer: SparsePoly = {}
    for alpha, coefficient_alpha in left.items():
        for beta, coefficient_beta in right.items():
            monomial = tuple(a + b for a, b in zip(alpha, beta))
            add_term(answer, monomial, coefficient_alpha * coefficient_beta)
    return answer


def pullback_columns(source: dict[str, object]) -> tuple[list[Exponent], list[SparsePoly]]:
    dimension = source["dimension"]
    zero = (0,) * dimension
    coordinates: list[SparsePoly] = []
    for index, component in enumerate(source["H"]):
        coordinate: SparsePoly = {}
        linear = [0] * dimension
        linear[index] = 1
        coordinate[tuple(linear)] = 1
        for term in component:
            exponent = [0] * dimension
            for variable, power in term["monomial"]:
                exponent[variable] = power
            add_term(coordinate, tuple(exponent), rational_mod(term["coefficient"], PRIME))
        coordinates.append(coordinate)

    basis: list[Exponent] = [zero]
    for index in range(dimension):
        exponent = [0] * dimension
        exponent[index] = 1
        basis.append(tuple(exponent))
    for left in range(dimension):
        for right in range(left, dimension):
            exponent = [0] * dimension
            exponent[left] += 1
            exponent[right] += 1
            basis.append(tuple(exponent))

    columns: list[SparsePoly] = []
    for exponent in basis:
        support = [index for index, power in enumerate(exponent) for _ in range(power)]
        if not support:
            pullback = {zero: 1}
        elif len(support) == 1:
            pullback = dict(coordinates[support[0]])
        else:
            assert len(support) == 2
            pullback = multiply(coordinates[support[0]], coordinates[support[1]])
        add_term(pullback, exponent, -1)
        columns.append(pullback)
    return basis, columns


def sparse_column_rank(columns: list[SparsePoly]) -> tuple[int, list[int]]:
    """Return modular column rank and the indices of independent columns."""
    pivots: dict[Exponent, SparsePoly] = {}
    independent: list[int] = []
    for column_index, original in enumerate(columns):
        column = dict(original)
        while column:
            pivot = min(column)
            if pivot not in pivots:
                inverse = pow(column[pivot], -1, PRIME)
                column = {
                    monomial: coefficient * inverse % PRIME
                    for monomial, coefficient in column.items()
                }
                pivots[pivot] = column
                independent.append(column_index)
                break
            factor = column[pivot]
            for monomial, coefficient in pivots[pivot].items():
                add_term(column, monomial, -factor * coefficient)
    return len(pivots), independent


def homogeneous_basis(dimension: int, degree: int) -> list[Exponent]:
    basis = []
    for indices in combinations_with_replacement(range(dimension), degree):
        exponent = [0] * dimension
        for index in indices:
            exponent[index] += 1
        basis.append(tuple(exponent))
    return basis


def decode_h_mod(source: dict[str, object]) -> list[SparsePoly]:
    dimension = source["dimension"]
    components = []
    for component in source["H"]:
        decoded: SparsePoly = {}
        for term in component:
            exponent = [0] * dimension
            for variable, power in term["monomial"]:
                exponent[variable] = power
            add_term(decoded, tuple(exponent), rational_mod(term["coefficient"], PRIME))
        components.append(decoded)
    return components


def lie_columns(
    source: dict[str, object], degree: int
) -> tuple[list[Exponent], list[SparsePoly]]:
    """Columns of P -> grad(P).H from Sym^degree to Sym^(degree+2)."""
    dimension = source["dimension"]
    components = decode_h_mod(source)
    basis = homogeneous_basis(dimension, degree)
    columns = []
    for exponent in basis:
        column: SparsePoly = {}
        for variable, power in enumerate(exponent):
            if not power:
                continue
            derivative = list(exponent)
            derivative[variable] -= 1
            for h_exponent, coefficient in components[variable].items():
                monomial = tuple(
                    derivative[index] + h_exponent[index]
                    for index in range(dimension)
                )
                add_term(column, monomial, power * coefficient)
        columns.append(column)
    return basis, columns


def near_invariant_defect(source: dict[str, object]) -> SparsePoly:
    """Return Q(V)-Q for Q=X_18*X_20-X_6*X_8."""
    dimension = source["dimension"]
    components = decode_h_mod(source)
    coordinates = []
    for index, component in enumerate(components):
        coordinate = dict(component)
        exponent = [0] * dimension
        exponent[index] = 1
        add_term(coordinate, tuple(exponent), 1)
        coordinates.append(coordinate)

    first = multiply(coordinates[18], coordinates[20])
    second = multiply(coordinates[6], coordinates[8])
    defect = dict(first)
    for monomial, coefficient in second.items():
        add_term(defect, monomial, -coefficient)

    first_original = [0] * dimension
    first_original[18] = 1
    first_original[20] = 1
    add_term(defect, tuple(first_original), -1)
    second_original = [0] * dimension
    second_original[6] = 1
    second_original[8] = 1
    add_term(defect, tuple(second_original), 1)
    return defect


def near_invariant_defect_exact(
    source: dict[str, object],
) -> RationalSparsePoly:
    """Compute the same defect over Q, independently of modular reduction."""
    dimension = source["dimension"]

    def add_rational(
        poly: RationalSparsePoly, monomial: Exponent, coefficient: Fraction
    ) -> None:
        updated = poly.get(monomial, Fraction(0)) + coefficient
        if updated:
            poly[monomial] = updated
        else:
            poly.pop(monomial, None)

    def multiply_rational(
        left: RationalSparsePoly, right: RationalSparsePoly
    ) -> RationalSparsePoly:
        answer: RationalSparsePoly = {}
        for alpha, coefficient_alpha in left.items():
            for beta, coefficient_beta in right.items():
                monomial = tuple(a + b for a, b in zip(alpha, beta))
                add_rational(answer, monomial, coefficient_alpha * coefficient_beta)
        return answer

    coordinates = []
    for index, component in enumerate(source["H"]):
        coordinate: RationalSparsePoly = {}
        linear = [0] * dimension
        linear[index] = 1
        coordinate[tuple(linear)] = Fraction(1)
        for term in component:
            exponent = [0] * dimension
            for variable, power in term["monomial"]:
                exponent[variable] = power
            add_rational(
                coordinate, tuple(exponent), Fraction(term["coefficient"])
            )
        coordinates.append(coordinate)

    first = multiply_rational(coordinates[18], coordinates[20])
    second = multiply_rational(coordinates[6], coordinates[8])
    defect = dict(first)
    for monomial, coefficient in second.items():
        add_rational(defect, monomial, -coefficient)

    first_original = [0] * dimension
    first_original[18] = 1
    first_original[20] = 1
    add_rational(defect, tuple(first_original), Fraction(-1))
    second_original = [0] * dimension
    second_original[6] = 1
    second_original[8] = 1
    add_rational(defect, tuple(second_original), Fraction(1))
    return defect


def main() -> None:
    source = json.loads(SOURCE.read_text())
    assert source["dimension"] == 21
    assert source["jacobian_determinant"] == "1"
    assert len(source["H"]) == 21
    assert source["H"][20] == []
    assert all(point[20] == "1" for point in source["collision_points"])

    basis, columns = pullback_columns(source)
    assert len(basis) == 253

    zero = (0,) * 21
    x20 = list(zero)
    x20[20] = 1
    x20_squared = list(zero)
    x20_squared[20] = 2
    known_invariants = {zero, tuple(x20), tuple(x20_squared)}
    known_indices = {index for index, monomial in enumerate(basis) if monomial in known_invariants}
    assert len(known_indices) == 3
    assert all(not columns[index] for index in known_indices)

    rank, independent = sparse_column_rank(columns)
    assert rank == 250
    assert len(independent) == 250
    assert not (known_indices & set(independent))

    # The three displayed invariants give nullity at least three over Q. The
    # rank-250 minor that is nonzero modulo PRIME is nonzero over Q, giving
    # rank at least 250 there.  Since there are 253 columns, equality follows.
    print("PASS quadratic invariant audit: 253 monomials of degree <= 2")
    print(f"PASS quadratic invariant audit: pullback-minus-identity rank {rank} mod {PRIME}")
    print("PASS quadratic invariant audit: Q-kernel is exactly span(1, X_20, X_20^2)")

    # Put Q=X_18*X_20-X_6*X_8.  Exact pullback gives
    # Q(V)-Q=-M*X_20^2 for M=X_0^2*X_1*X_2, so Q is a first
    # integral of the cubic vector field H but not an invariant of V.
    defect = near_invariant_defect(source)
    defect_exponent = [0] * 21
    defect_exponent[0] = 2
    defect_exponent[1] = 1
    defect_exponent[2] = 1
    defect_exponent[20] = 2
    assert defect == {tuple(defect_exponent): PRIME - 1}
    assert near_invariant_defect_exact(source) == {
        tuple(defect_exponent): Fraction(-1)
    }

    expected = {
        2: (231, 229),
        3: (1771, 1769),
        4: (10626, 10623),
        5: (53130, 53127),
    }
    computed: dict[int, tuple[list[Exponent], list[SparsePoly], int]] = {}
    for degree, (column_count, expected_rank) in expected.items():
        homogeneous, derivative_columns = lie_columns(source, degree)
        derivative_rank, _independent = sparse_column_rank(derivative_columns)
        assert len(homogeneous) == column_count
        assert derivative_rank == expected_rank
        computed[degree] = (homogeneous, derivative_columns, derivative_rank)
        print(
            f"PASS degree-{degree} Lie audit: rank {derivative_rank} "
            f"on {column_count} monomials"
        )

    # In degrees 2,...,5 the displayed nullities are attained by
    #   Q[s,Q]_d.
    # Hence those are the exact Lie kernels over Q.  To mix the defect of Q
    # with a quartic correction would require M*s^2 in im(L_H:Sym^4->Sym^6);
    # the odd analogue requires M*s^3 in im(L_H:Sym^5->Sym^7).
    # Appending either target raises modular rank, excluding both lifts over Q.
    for degree, s_power in ((4, 2), (5, 3)):
        _basis_d, derivative_columns, derivative_rank = computed[degree]
        target = [0] * 21
        target[0] = 2
        target[1] = 1
        target[2] = 1
        target[20] = s_power
        augmented_rank, _augmented_independent = sparse_column_rank(
            derivative_columns + [{tuple(target): 1}]
        )
        assert augmented_rank == derivative_rank + 1
        print(
            f"PASS degree-{degree} lift obstruction: "
            f"X_0^2*X_1*X_2*X_20^{s_power} is not in the Lie image"
        )

    print("PASS near-invariant audit: Q=X_18*X_20-X_6*X_8")
    print("  Q(V)-Q=-X_0^2*X_1*X_2*X_20^2")
    print("PASS low-degree invariant audit: fixed space through degree 5 is Q[X_20]_{<=5}")
    print("PASS low-degree invariant audit: no independent invariant identity slice through degree 5")


if __name__ == "__main__":
    main()
