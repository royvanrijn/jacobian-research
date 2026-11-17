#!/usr/bin/env python3
"""Prove function-field polynomial solvability of the degree-five hbar^3 equation.

The full bounded correction space has 1569 monomials.  Exact computations
at ordinary fibers reveal a much smaller stable support: 27 monomials for
``S_2`` and 15 for ``T_2``.  This checker works over the characteristic-zero
function field ``QQ(a,tau)`` and proves that those 42 columns already solve

    {S_2,T} + {S,T_2} = -Pi^3(S,T)/24

for the complete generic-chart family.  It clears common denominators from
``S`` and ``T`` before differentiating, performs one 42-variable function
field row reduction, and verifies the result against every output equation.

The parameter ``a`` is related to the normalized seed parameter by

    kappa = -(1 + 2*a)/(1 + a).

Thus the certificate is valid on the generic chart ``a != 0,-1`` (equivalently
``kappa != -1,-2``).  The same checker also proves the result over the
one-parameter replacement chart ``QQ(tau)`` at ``kappa=-1``.
"""

from __future__ import annotations

from time import perf_counter

from sympy.polys.domains import QQ

from explore_degree_five_a2_subprincipal import add
from explore_degree_five_quantum_residue import (
    degree_five_family,
    degree_five_exceptional_family,
    exceptional_pi_power,
    exceptional_poisson,
    pi_power,
    poisson,
    scale,
    solve_affine,
)


S_SUPPORT = [
    (4, 0, 0),
    (5, 1, 0),
    (6, 0, 1),
    (6, 2, 0),
    (7, 1, 1),
    (7, 3, 0),
    (8, 0, 2),
    (8, 2, 1),
    (9, 1, 2),
    (9, 3, 1),
    (9, 5, 0),
    (10, 0, 3),
    (10, 2, 2),
    (10, 4, 1),
    (10, 6, 0),
    (11, 1, 3),
    (11, 3, 2),
    (11, 5, 1),
    (11, 7, 0),
    (12, 2, 3),
    (12, 4, 2),
    (12, 6, 1),
    (12, 8, 0),
    (13, 3, 3),
    (13, 5, 2),
    (13, 7, 1),
    (13, 9, 0),
]

T_SUPPORT = [
    (5, 0, 0),
    (6, 1, 0),
    (7, 0, 1),
    (7, 2, 0),
    (8, 1, 1),
    (8, 3, 0),
    (9, 0, 2),
    (9, 2, 1),
    (10, 1, 2),
    (10, 3, 1),
    (11, 2, 2),
    (11, 4, 1),
    (11, 6, 0),
    (12, 3, 2),
    (12, 7, 0),
]

EXCEPTIONAL_S_SUPPORT = [
    (4, 0, 0),
    (5, 0, 1),
    (6, 0, 2),
    (6, 1, 0),
    (7, 0, 3),
    (7, 1, 1),
    (8, 1, 2),
    (8, 2, 0),
    (9, 1, 3),
    (9, 2, 1),
    (10, 2, 2),
    (10, 3, 0),
    (11, 2, 3),
    (11, 3, 1),
    (12, 3, 2),
    (13, 3, 3),
    (13, 4, 1),
    (14, 4, 2),
    (14, 5, 0),
    (15, 5, 1),
    (16, 6, 0),
]

EXCEPTIONAL_T_SUPPORT = [
    (5, 0, 0),
    (6, 0, 1),
    (7, 0, 2),
    (7, 1, 0),
    (8, 1, 1),
    (9, 1, 2),
    (9, 2, 0),
    (10, 2, 1),
    (11, 2, 2),
    (11, 3, 0),
    (12, 3, 1),
    (13, 3, 2),
    (14, 4, 1),
]


def clear_coefficient(value, denominator, fraction_field, polynomial_ring):
    """Clear ``denominator`` and return a polynomial-ring coefficient."""

    cleared = fraction_field(denominator) * value
    if cleared.denom.degree() != 0:
        raise AssertionError(f"uncancelled denominator: {cleared.denom}")
    unit = cleared.denom[(0, 0)]
    return cleared.numer * (polynomial_ring.domain.one / unit)


def verify_generic() -> dict[str, object]:
    field = QQ.frac_field("a", "tau")
    a, tau = field.gens
    ring = a.numer.ring

    started = perf_counter()
    S, T = degree_five_family(
        field,
        a,
        tau,
        verify_canonical=False,
    )

    # These are common denominators, not necessarily least common
    # denominators.  Bilinearity gives the cleared equation
    #
    # Ds*{S_2,B} + Dt*{A,T_2} = -Pi^3(A,B)/24,
    #
    # where A=Ds*S and B=Dt*T.
    a_poly = a.numer
    denominator_s = a_poly**4 * (a_poly + ring.one) ** 10
    denominator_t = a_poly**4 * (a_poly + ring.one) ** 8
    A = {
        monomial: clear_coefficient(
            coefficient,
            denominator_s,
            field,
            ring,
        )
        for monomial, coefficient in S.items()
    }
    B = {
        monomial: clear_coefficient(
            coefficient,
            denominator_t,
            field,
            ring,
        )
        for monomial, coefficient in T.items()
    }

    columns = [
        scale(poisson({monomial: ring.one}, B), denominator_s)
        for monomial in S_SUPPORT
    ]
    columns += [
        scale(poisson(A, {monomial: ring.one}), denominator_t)
        for monomial in T_SUPPORT
    ]
    rhs = scale(
        pi_power(A, B, 3),
        -ring.domain.one / ring.domain(24),
    )

    field_columns = [
        {monomial: field(coefficient) for monomial, coefficient in column.items()}
        for column in columns
    ]
    field_rhs = {
        monomial: field(coefficient) for monomial, coefficient in rhs.items()
    }
    particular, kernel, rank = solve_affine(
        field_columns,
        field_rhs,
        field,
    )

    residual = dict(field_rhs)
    for column_index, coefficient in particular.items():
        residual = add(
            residual,
            field_columns[column_index],
            -coefficient,
        )
    if residual:
        raise AssertionError(f"nonzero residual with {len(residual)} terms")
    if rank != 42 or kernel or len(particular) != 42:
        raise AssertionError(
            "unexpected sparse system dimensions: "
            f"rank={rank}, kernel={len(kernel)}, support={len(particular)}"
        )
    for coefficient in particular.values():
        cleared = field(denominator_s) * coefficient
        if cleared.denom.degree() != 0:
            raise AssertionError(
                "the hbar^3 lift has an unexpected parameter pole: "
                f"{coefficient.denom}"
            )

    return {
        "chart": "generic",
        "field": "QQ(a,tau)",
        "columns": len(columns),
        "output_equations": len(field_rhs),
        "rank": rank,
        "kernel_dimension": len(kernel),
        "solution_support": len(particular),
        "residual_terms": len(residual),
        "numerator_term_range": (
            min(len(coefficient.numer) for coefficient in particular.values()),
            max(len(coefficient.numer) for coefficient in particular.values()),
        ),
        "denominator_degrees": sorted(
            {coefficient.denom.degree() for coefficient in particular.values()}
        ),
        "additional_parameter_poles": 0,
        "seconds": round(perf_counter() - started, 2),
    }


def verify_exceptional() -> dict[str, object]:
    field = QQ.frac_field("tau")
    tau = field.gens[0]
    ring = tau.numer.ring
    started = perf_counter()

    S, T = degree_five_exceptional_family(
        field,
        tau,
        verify_canonical=False,
    )

    def polynomial_coefficient(coefficient):
        if coefficient.denom.degree() != 0:
            raise AssertionError(
                f"unexpected exceptional denominator: {coefficient.denom}"
            )
        unit = coefficient.denom[(0,)]
        return coefficient.numer * (ring.domain.one / unit)

    A = {
        monomial: polynomial_coefficient(coefficient)
        for monomial, coefficient in S.items()
    }
    B = {
        monomial: polynomial_coefficient(coefficient)
        for monomial, coefficient in T.items()
    }
    columns = [
        exceptional_poisson({monomial: ring.one}, B)
        for monomial in EXCEPTIONAL_S_SUPPORT
    ]
    columns += [
        exceptional_poisson(A, {monomial: ring.one})
        for monomial in EXCEPTIONAL_T_SUPPORT
    ]
    rhs = scale(
        exceptional_pi_power(A, B, 3),
        -ring.domain.one / ring.domain(24),
    )

    field_columns = [
        {monomial: field(coefficient) for monomial, coefficient in column.items()}
        for column in columns
    ]
    field_rhs = {
        monomial: field(coefficient) for monomial, coefficient in rhs.items()
    }
    particular, kernel, rank = solve_affine(
        field_columns,
        field_rhs,
        field,
    )
    residual = dict(field_rhs)
    for column_index, coefficient in particular.items():
        residual = add(
            residual,
            field_columns[column_index],
            -coefficient,
        )
    if residual:
        raise AssertionError(f"nonzero residual with {len(residual)} terms")
    if rank != 34 or kernel or len(particular) != 34:
        raise AssertionError(
            "unexpected exceptional sparse system dimensions: "
            f"rank={rank}, kernel={len(kernel)}, support={len(particular)}"
        )

    return {
        "chart": "kappa=-1 replacement",
        "field": "QQ(tau)",
        "columns": len(columns),
        "output_equations": len(field_rhs),
        "rank": rank,
        "kernel_dimension": len(kernel),
        "solution_support": len(particular),
        "residual_terms": len(residual),
        "numerator_term_range": (
            min(len(coefficient.numer) for coefficient in particular.values()),
            max(len(coefficient.numer) for coefficient in particular.values()),
        ),
        "denominator_degrees": sorted(
            {coefficient.denom.degree() for coefficient in particular.values()}
        ),
        "seconds": round(perf_counter() - started, 2),
    }


def main() -> None:
    print("degree-five hbar^3 function-field certificates")
    for result in (verify_generic(), verify_exceptional()):
        print(f"[{result['chart']}]")
        for key, value in result.items():
            if key != "chart":
                print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
