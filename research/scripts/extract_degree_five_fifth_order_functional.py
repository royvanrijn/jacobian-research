#!/usr/bin/env python3
"""Extract an exact bounded hbar^5 dual functional at a rational seed.

The complete fifth-order span test proves nonexistence when the constant
defect lies outside the span of the correction image and every linear or
quadratic lower-lift coefficient.  This script converts that rank jump into
an explicit functional Lambda:

    Lambda(span columns) = 0,    Lambda(constant defect) = 1.

It solves only for one point of this affine dual space and therefore avoids
constructing its usually very large homogeneous kernel.
"""

from __future__ import annotations

import argparse
from fractions import Fraction

from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref

from explore_degree_five_quantum_residue import (
    GENERIC_PROFILE,
    degree_five_family,
    fifth_order_coefficients,
    laurent_monomials,
    third_order_family,
)


def parse_rational(text: str):
    value = Fraction(text)
    return QQ(value.numerator, value.denominator)


def affine_particular(columns, rhs, field):
    """Solve columns*x=rhs without constructing the homogeneous kernel."""

    output_coordinates = sorted(
        set(rhs).union(*(set(column) for column in columns))
    )
    output_index = {
        coordinate: index
        for index, coordinate in enumerate(output_coordinates)
    }
    rhs_column = len(columns)
    rows = {}
    for column_index, column in enumerate(columns):
        for coordinate, coefficient in column.items():
            rows.setdefault(output_index[coordinate], {})[
                column_index
            ] = coefficient
    for coordinate, coefficient in rhs.items():
        rows.setdefault(output_index[coordinate], {})[
            rhs_column
        ] = -coefficient
    reduced, pivots, _ = sdm_irref(rows)
    if rhs_column in pivots:
        raise ValueError("dual functional equation is inconsistent")
    particular = {}
    for reduced_row, pivot in enumerate(pivots):
        rhs_value = reduced.get(reduced_row, {}).get(
            rhs_column,
            field.zero,
        )
        if rhs_value:
            particular[pivot] = -rhs_value
    return particular, len(pivots)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kappa", default="0")
    parser.add_argument("--tau", default="-3")
    args = parser.parse_args()

    field = QQ
    kappa = parse_rational(args.kappa)
    tau = parse_rational(args.tau)
    if kappa in (-field(2), -field.one):
        parser.error("this extractor uses the generic chart")
    a = -(field.one + kappa) / (field(2) + kappa)
    S, T = degree_five_family(field, a, tau)
    third = third_order_family(S, T, field)
    constant, nonconstant = fifth_order_coefficients(
        S,
        T,
        third,
        field,
    )
    s4_monomials = laurent_monomials(
        GENERIC_PROFILE.s4_degree,
        1,
        0,
        GENERIC_PROFILE.z_weight,
    )
    t4_monomials = laurent_monomials(
        GENERIC_PROFILE.t4_degree,
        0,
        0,
        GENERIC_PROFILE.z_weight,
    )
    corrections = [
        GENERIC_PROFILE.poisson({monomial: field.one}, T)
        for monomial in s4_monomials
    ]
    corrections += [
        GENERIC_PROFILE.poisson(S, {monomial: field.one})
        for monomial in t4_monomials
    ]
    span = corrections + nonconstant
    output_monomials = sorted(
        set(constant).union(*(set(column) for column in span))
    )

    # Each output monomial is a possible coordinate of Lambda.  Its column
    # records evaluation on every span vector and on the constant defect.
    dual_columns = []
    constant_equation = len(span)
    for monomial in output_monomials:
        column = {
            span_index: polynomial.get(monomial, field.zero)
            for span_index, polynomial in enumerate(span)
            if polynomial.get(monomial, field.zero)
        }
        constant_value = constant.get(monomial, field.zero)
        if constant_value:
            column[constant_equation] = constant_value
        dual_columns.append(column)
    functional_vector, dual_rank = affine_particular(
        dual_columns,
        {constant_equation: field.one},
        field,
    )
    functional = {
        output_monomials[index]: value
        for index, value in functional_vector.items()
        if value
    }
    for polynomial in span:
        value = sum(
            functional.get(monomial, field.zero) * coefficient
            for monomial, coefficient in polynomial.items()
        )
        if value:
            raise AssertionError("functional does not annihilate the span")
    period = sum(
        functional.get(monomial, field.zero) * coefficient
        for monomial, coefficient in constant.items()
    )
    if period != field.one:
        raise AssertionError(period)

    print(
        f"SEED_KAPPA={kappa},TAU={tau},A={a}",
    )
    print(
        f"H3_KERNEL={len(third.kernel)},"
        f"SPAN_COLUMNS={len(span)},"
        f"DUAL_RANK={dual_rank}",
    )
    print(f"FUNCTIONAL_SUPPORT={len(functional)}")
    print(f"FUNCTIONAL_PERIOD={period}")
    for monomial, coefficient in functional.items():
        print(f"  {monomial}: {coefficient}")


if __name__ == "__main__":
    main()
