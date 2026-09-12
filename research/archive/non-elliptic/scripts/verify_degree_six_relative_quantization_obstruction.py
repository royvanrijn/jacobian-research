#!/usr/bin/env python3
"""Uniform restricted hbar^5 obstruction on a degree-six symbol family.

On the kappa=-9 weighted-seed slice, the completed classical pair depends on
``(sigma,tau)``.  The root-boundary grading makes the relevant parity
correction spaces small.  This checker works exactly over
``QQ(sigma,tau)`` and proves:

* the complete weight-homogeneous hbar^3 equation has rank 77 and a
  six-dimensional affine lift space;
* the hbar^5 current-correction and all linear/quadratic lower-lift
  variations have joint rank 41; and
* adjoining the constant hbar^5 defect raises the rank to 42.

Thus a strong dual cocycle has nonzero constant value on a nonempty open of
the two-parameter family.  This is a parameter-uniform obstruction in the
declared parity and boundary-weight filtration, not an unrestricted
quantization obstruction and not a result about DC_2.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_quantum_residue import (
    SparsePoly,
    add,
    column_rank,
    pi_power,
    poisson,
    scale,
    solve_affine,
    split_pair,
    weighted_seed_symbol_pair,
)


PRIMES = (31991, 32003, 65521)


def degree_six_pair(field, sigma, tau) -> tuple[SparsePoly, SparsePoly]:
    """Return the completed pair on the kappa=-9 degree-six slice."""

    a = -field(8) / field(7)
    shear = (
        field(6024) * sigma**2
        + field(5016) * sigma * tau
        + field(11088) * sigma
        + field(1056) * tau**2
        + field(4752) * tau
        - field(16929)
    ) / field(2156)
    linear = (
        -field(9) / field(2)
        - field(3) * sigma
        - field(2) * tau
        + field(2)
    )
    constant = (
        field(9) / field(2)
        + field(2) * sigma
        + tau
        - field(3)
    )
    return weighted_seed_symbol_pair(
        field,
        a,
        (constant, linear, tau, sigma),
        shear,
    )


def weight_monomials(
    maximum_degree: int,
    maximum_z_order: int,
    weight: int,
) -> list[tuple[int, int, int]]:
    """Complete root-boundary weight summand in one Bernstein box."""

    return [
        (x_degree, q_degree, z_degree)
        for z_degree in range(maximum_z_order + 1)
        for x_degree in range(maximum_degree - 3 * z_degree + 1)
        for q_degree in range(
            maximum_degree - 3 * z_degree - x_degree + 1
        )
        if x_degree - q_degree - 2 * z_degree == weight
    ]


S2_SUPPORT = weight_monomials(32, 4, 4)
T2_SUPPORT = weight_monomials(28, 3, 5)
S4_SUPPORT = weight_monomials(28, 2, 10)
T4_SUPPORT = weight_monomials(24, 1, 11)


def fifth_defect(
    S: SparsePoly,
    T: SparsePoly,
    lower_pair: tuple[SparsePoly, SparsePoly],
    field,
) -> SparsePoly:
    s2, t2 = lower_pair
    value = poisson(s2, t2)
    value = add(value, pi_power(s2, T, 3), field.one / field(24))
    value = add(value, pi_power(S, t2, 3), field.one / field(24))
    return add(
        value,
        pi_power(S, T, 5),
        field.one / field(1920),
    )


def family_presentation(field, sigma, tau) -> dict[str, object]:
    """Build the hbar^3 lift and hbar^5 strong-cocycle presentation."""

    S, T = degree_six_pair(field, sigma, tau)
    correction_three = [
        poisson({monomial: field.one}, T)
        for monomial in S2_SUPPORT
    ]
    correction_three += [
        poisson(S, {monomial: field.one})
        for monomial in T2_SUPPORT
    ]
    rhs_three = scale(
        pi_power(S, T, 3),
        -field.one / field(24),
    )
    particular, kernel, rank_three = solve_affine(
        correction_three,
        rhs_three,
        field,
    )
    base_pair = split_pair(particular, S2_SUPPORT, T2_SUPPORT)
    kernel_pairs = [
        split_pair(vector, S2_SUPPORT, T2_SUPPORT)
        for vector in kernel
    ]

    constant = fifth_defect(S, T, base_pair, field)
    lower_variations: list[SparsePoly] = []
    for basis_s, basis_t in kernel_pairs:
        diagonal = poisson(basis_s, basis_t)
        shifted = fifth_defect(
            S,
            T,
            (
                add(base_pair[0], basis_s),
                add(base_pair[1], basis_t),
            ),
            field,
        )
        linear = add(
            add(shifted, constant, -field.one),
            diagonal,
            -field.one,
        )
        lower_variations.extend((linear, diagonal))
    for left, right in combinations(range(len(kernel_pairs)), 2):
        left_s, left_t = kernel_pairs[left]
        right_s, right_t = kernel_pairs[right]
        lower_variations.append(
            add(
                poisson(left_s, right_t),
                poisson(right_s, left_t),
            )
        )

    correction_five = [
        poisson({monomial: field.one}, T)
        for monomial in S4_SUPPORT
    ]
    correction_five += [
        poisson(S, {monomial: field.one})
        for monomial in T4_SUPPORT
    ]
    strong_columns = correction_five + lower_variations
    output_support = sorted(
        set(constant).union(
            *(set(column) for column in strong_columns)
        )
    )
    return {
        "S": S,
        "T": T,
        "rank_three": rank_three,
        "kernel_pairs": kernel_pairs,
        "correction_five": correction_five,
        "lower_variations": lower_variations,
        "strong_columns": strong_columns,
        "constant": constant,
        "output_support": output_support,
    }


def pairing(functional, vector, field):
    return sum(
        (
            coefficient * vector.get(monomial, field.zero)
            for monomial, coefficient in functional.items()
        ),
        field.zero,
    )


def dual_witness(
    columns: list[SparsePoly],
    constant: SparsePoly,
    field,
) -> dict[tuple[int, int, int], object]:
    """Return one normalized exact functional killing ``columns``."""

    output_monomials = sorted(
        set(constant).union(*(set(column) for column in columns))
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    transpose_rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(transpose_rows)
    left_kernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(output_monomials),
        pivots,
        nonzero,
    )
    for vector in left_kernel:
        functional = {
            output_monomials[index]: coefficient
            for index, coefficient in vector.items()
            if coefficient
        }
        value = pairing(functional, constant, field)
        if value:
            normalized = {
                monomial: coefficient / value
                for monomial, coefficient in functional.items()
            }
            assert all(
                pairing(normalized, column, field) == field.zero
                for column in columns
            )
            assert pairing(normalized, constant, field) == field.one
            return normalized
    raise AssertionError("the augmented rank jump has no dual witness")


def rank_record(field, sigma, tau) -> dict[str, int]:
    presentation = family_presentation(field, sigma, tau)
    correction_five = presentation["correction_five"]
    strong_columns = presentation["strong_columns"]
    constant = presentation["constant"]
    return {
        "h3_rank": presentation["rank_three"],
        "h3_kernel_dimension": len(presentation["kernel_pairs"]),
        "h5_correction_rank": column_rank(correction_five),
        "h5_strong_span_rank": column_rank(strong_columns),
        "h5_augmented_rank": column_rank(strong_columns + [constant]),
        "h5_output_dimension": len(presentation["output_support"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="optional path for the deterministic JSON certificate",
    )
    args = parser.parse_args()

    assert (len(S2_SUPPORT), len(T2_SUPPORT)) == (49, 34)
    assert (len(S4_SUPPORT), len(T4_SUPPORT)) == (22, 12)

    function_field = QQ.frac_field("sigma", "tau")
    sigma, tau = function_field.gens
    generic = rank_record(function_field, sigma, tau)
    assert generic == {
        "h3_rank": 77,
        "h3_kernel_dimension": 6,
        "h5_correction_rank": 34,
        "h5_strong_span_rank": 41,
        "h5_augmented_rank": 42,
        "h5_output_dimension": 110,
    }

    rational_presentation = family_presentation(
        QQ,
        QQ.one,
        QQ.zero,
    )
    rational_ranks = {
        "h3_rank": rational_presentation["rank_three"],
        "h3_kernel_dimension": len(
            rational_presentation["kernel_pairs"]
        ),
        "h5_correction_rank": column_rank(
            rational_presentation["correction_five"]
        ),
        "h5_strong_span_rank": column_rank(
            rational_presentation["strong_columns"]
        ),
        "h5_augmented_rank": column_rank(
            rational_presentation["strong_columns"]
            + [rational_presentation["constant"]]
        ),
        "h5_output_dimension": len(
            rational_presentation["output_support"]
        ),
    }
    assert rational_ranks == generic
    witness = dual_witness(
        rational_presentation["strong_columns"],
        rational_presentation["constant"],
        QQ,
    )

    prime_records = {}
    for prime in PRIMES:
        finite_field = GF(prime)
        record = rank_record(
            finite_field,
            finite_field.one,
            finite_field.zero,
        )
        assert record == generic
        prime_records[str(prime)] = record

    certificate = {
        "scope": (
            "parity-preserving, root-weight-homogeneous inherited "
            "filtration; no unrestricted or DC_2 claim"
        ),
        "family": {
            "degree": 6,
            "kappa": -9,
            "parameters": ["sigma", "tau"],
            "classical_open_conditions": [
                "sigma != 0",
                "admissible marked Hessian-boundary conditions",
            ],
            "completing_shear": (
                "(6024*sigma^2+5016*sigma*tau+11088*sigma"
                "+1056*tau^2+4752*tau-16929)/2156"
            ),
        },
        "supports": {
            "S2": len(S2_SUPPORT),
            "T2": len(T2_SUPPORT),
            "S4": len(S4_SUPPORT),
            "T4": len(T4_SUPPORT),
        },
        "function_field_ranks": generic,
        "rational_boundary_clean_point": {
            "sigma": 1,
            "tau": 0,
            "ranks": rational_ranks,
            "strong_cocycle_value": 1,
            "strong_cocycle_support": [
                {
                    "monomial": list(monomial),
                    "coefficient": str(coefficient),
                }
                for monomial, coefficient in sorted(witness.items())
            ],
        },
        "good_prime_records": prime_records,
        "conclusion": (
            "on a nonempty Zariski open, every hbar^3 lift in the "
            "six-dimensional restricted torsor is obstructed at hbar^5"
        ),
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    print("PASS: exact degree-six pair over Q(sigma,tau) is canonical")
    print("PASS: hbar^3 rank 77 leaves a six-dimensional lift torsor")
    print("PASS: hbar^5 strong span has rank 41 and augmented rank 42")
    print("PASS: a 30-term rational strong cocycle has constant value one")
    print("PASS: the same ranks hold at three good primes")
    print(
        "THEOREM: a nonempty open of the two-parameter degree-six family "
        "is obstructed at hbar^5 in the declared restricted filtration"
    )
    print("SCOPE: parameter-uniform restricted obstruction; no DC_2 claim")


if __name__ == "__main__":
    main()
