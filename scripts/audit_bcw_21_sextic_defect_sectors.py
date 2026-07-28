#!/usr/bin/env python3
"""Exclude sextic corrections of the two quartic near-invariant directions.

For the stored cubic vector field H, two exact torus gradings split

    L_H : Sym^6 -> Sym^8

into small blocks.  The degree-eight defects of X_20^2*Q and Q^2, where
Q=X_18*X_20-X_6*X_8, lie in weight sectors (1,3) and (2,2).  This audit proves
modulo the good prime 1000003 that the corresponding Lie blocks have ranks
103 and 1604 (full column rank) and that the two defects are outside their
images.  The nonzero modular minors give characteristic-zero certificates.

This closes sextic correction of lower-degree near-invariants.  It does not
classify genuinely new homogeneous sextic invariants in other weight sectors.
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
DIMENSION = 21
DEGREE = 6

WEIGHT_1 = (-1, 1, 2, 1, 1, 0, 1, 2, 0, -2, 0, 1, 1, 2, -1, 1, 2, 1, 1, 0, 0)
WEIGHT_2 = (-4, 2, 5, 2, 2, -1, 2, 5, -1, -7, -1, 2, 2, 5, -6, 0, 3, 0, 0, -3, 1)

PackedPoly = dict[int, int]
DecodedComponent = list[tuple[int, int]]


def rational_mod(value: str) -> int:
    number = Fraction(value)
    return number.numerator * pow(number.denominator, -1, PRIME) % PRIME


def pack(exponent: list[int] | tuple[int, ...]) -> int:
    return sum(power << (4 * variable) for variable, power in enumerate(exponent))


def add_term(poly: PackedPoly, monomial: int, coefficient: int) -> None:
    updated = (poly.get(monomial, 0) + coefficient) % PRIME
    if updated:
        poly[monomial] = updated
    else:
        poly.pop(monomial, None)


def bidegree(exponent: list[int] | tuple[int, ...]) -> tuple[int, int]:
    return (
        sum(weight * power for weight, power in zip(WEIGHT_1, exponent)),
        sum(weight * power for weight, power in zip(WEIGHT_2, exponent)),
    )


def decode_h(source: dict[str, object]) -> list[DecodedComponent]:
    components = []
    for output, component in enumerate(source["H"]):
        decoded = []
        for term in component:
            exponent = [0] * DIMENSION
            for variable, power in term["monomial"]:
                exponent[variable] = power
            # Both gradings make every H_output term have the weight of X_output.
            assert bidegree(exponent) == (
                WEIGHT_1[output],
                WEIGHT_2[output],
            )
            decoded.append((pack(exponent), rational_mod(term["coefficient"])))
        components.append(decoded)
    return components


def reduce_column(
    column: PackedPoly, pivots: dict[int, PackedPoly]
) -> tuple[bool, PackedPoly]:
    while column:
        pivot = min(column)
        coefficient = column[pivot]
        if pivot not in pivots:
            return False, column
        for monomial, value in pivots[pivot].items():
            add_term(column, monomial, -coefficient * value)
    return True, column


def sector_columns(
    components: list[DecodedComponent], sector: tuple[int, int]
) -> list[PackedPoly]:
    columns = []
    for indices in combinations_with_replacement(range(DIMENSION), DEGREE):
        if (
            sum(WEIGHT_1[index] for index in indices),
            sum(WEIGHT_2[index] for index in indices),
        ) != sector:
            continue
        exponent = [0] * DIMENSION
        packed_exponent = 0
        for index in indices:
            exponent[index] += 1
            packed_exponent += 1 << (4 * index)
        column: PackedPoly = {}
        for variable, power in enumerate(exponent):
            if not power:
                continue
            base = packed_exponent - (1 << (4 * variable))
            for h_monomial, coefficient in components[variable]:
                add_term(column, base + h_monomial, power * coefficient)
        columns.append(column)
    return columns


def rank_and_target(
    columns: list[PackedPoly], target: PackedPoly
) -> tuple[int, bool, int]:
    pivots: dict[int, PackedPoly] = {}
    rank = 0
    for original in columns:
        inside, residual = reduce_column(dict(original), pivots)
        if inside:
            continue
        pivot = min(residual)
        inverse = pow(residual[pivot], -1, PRIME)
        pivots[pivot] = {
            monomial: coefficient * inverse % PRIME
            for monomial, coefficient in residual.items()
        }
        rank += 1
    inside, residual = reduce_column(dict(target), pivots)
    return rank, inside, len(residual)


def main() -> None:
    source = json.loads(SOURCE.read_text())
    assert source["dimension"] == DIMENSION
    assert source["H"][20] == []
    components = decode_h(source)

    # M=X_0^2*X_1*X_2 and s=X_20.
    first_target_exponent = [0] * DIMENSION
    first_target_exponent[0] = 2
    first_target_exponent[1] = 1
    first_target_exponent[2] = 1
    first_target_exponent[20] = 4
    first_target = {pack(first_target_exponent): 1}  # M*s^4
    assert bidegree(first_target_exponent) == (1, 3)

    # M*s^2*Q = M*X_18*s^3 - M*X_6*X_8*s^2.
    second_left = [0] * DIMENSION
    second_left[0] = 2
    second_left[1] = 1
    second_left[2] = 1
    second_left[18] = 1
    second_left[20] = 3
    second_right = [0] * DIMENSION
    second_right[0] = 2
    second_right[1] = 1
    second_right[2] = 1
    second_right[6] = 1
    second_right[8] = 1
    second_right[20] = 2
    assert bidegree(second_left) == bidegree(second_right) == (2, 2)
    second_target = {
        pack(second_left): 1,
        pack(second_right): PRIME - 1,
    }

    first_columns = sector_columns(components, (1, 3))
    first_rank, first_inside, first_residual = rank_and_target(
        first_columns, first_target
    )
    assert len(first_columns) == first_rank == 103
    assert not first_inside
    assert first_residual

    second_columns = sector_columns(components, (2, 2))
    second_rank, second_inside, second_residual = rank_and_target(
        second_columns, second_target
    )
    assert len(second_columns) == second_rank == 1604
    assert not second_inside
    assert second_residual

    print("PASS sextic sector audit: both stored torus gradings are exact")
    print(f"PASS sextic sector (1,3): rank {first_rank}/103 mod {PRIME}")
    print("PASS sextic sector (1,3): M*X_20^4 is outside the Lie image")
    print(f"PASS sextic sector (2,2): rank {second_rank}/1604 mod {PRIME}")
    print("PASS sextic sector (2,2): M*X_20^2*Q is outside the Lie image")
    print("THEOREM: no sextic correction rescues the lower-degree Q directions")
    print("OPEN: genuinely new homogeneous sextic invariants in other weight sectors")


if __name__ == "__main__":
    main()
