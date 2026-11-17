#!/usr/bin/env python3
"""Restricted order-five presentation for a degree-seven marked-root slice.

The classical seed is

    H = w^2 (w-1) P,
    P = w^4 + sigma*w^3 + tau*w^2
        + (-13/2-3*sigma-2*tau)*w + 9/2+2*sigma+tau.

It has ``H'(1)=-1`` and ``H''(1)=-9``.  The coefficient of ``w^4`` in
``P`` is fixed to one, leaving an explicit two-parameter, exact-degree-seven
slice of the normalized weighted family.  The completing Q^2 shear was
reconstructed by exact quadratic interpolation of the uniform residue
functional; ``derive_degree_seven_marked_root_shear.py`` independently
replays those homotopy calculations.

This checker constructs the full inherited root-weight summands through
order five.  Its conclusions are restricted to that parity-preserving PBW
filtration and make no unrestricted quantization or DC_2 claim.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

from sympy.polys.domains import GF, QQ

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


GOOD_PRIMES = (17, 19, 23, 29, 31)


def degree_seven_shear(field, sigma, tau):
    """Unique classical completing shear on the selected monic slice."""

    numerator = field(3) * (
        field(26104) * sigma**2
        + field(21736) * sigma * tau
        + field(134160) * sigma
        + field(4576) * tau**2
        + field(56160) * tau
        + field(75285)
    )
    return numerator / field(28028)


def degree_seven_pair(field, sigma, tau) -> tuple[SparsePoly, SparsePoly]:
    """Return the completed canonical pair on the monic degree-seven slice."""

    a = -field(8) / field(7)
    factor_coefficients = (
        field(9) / field(2) + field(2) * sigma + tau,
        -field(13) / field(2) - field(3) * sigma - field(2) * tau,
        tau,
        sigma,
        field.one,
    )
    return weighted_seed_symbol_pair(
        field,
        a,
        factor_coefficients,
        degree_seven_shear(field, sigma, tau),
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


# The classical bounds are (deg_Z,deg_B)=(7,43) for S and (6,39) for T.
# Each even PBW step lowers differential order by two and Bernstein degree by
# four, while the root-boundary weight rises by six.
S2_SUPPORT = weight_monomials(39, 5, 4)
T2_SUPPORT = weight_monomials(35, 4, 5)
S4_SUPPORT = weight_monomials(35, 3, 10)
T4_SUPPORT = weight_monomials(31, 2, 11)


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
    return add(value, pi_power(S, T, 5), field.one / field(1920))


def relative_presentation(
    field,
    S,
    T,
    s2_support,
    t2_support,
    s4_support,
    t4_support,
) -> dict[str, object]:
    """Build the inherited order-three/order-five data for one symbol pair."""

    correction_three = [
        poisson({monomial: field.one}, T) for monomial in s2_support
    ]
    correction_three += [
        poisson(S, {monomial: field.one}) for monomial in t2_support
    ]
    rhs_three = scale(pi_power(S, T, 3), -field.one / field(24))
    particular, kernel, rank_three = solve_affine(
        correction_three,
        rhs_three,
        field,
    )
    base_pair = split_pair(particular, s2_support, t2_support)
    kernel_pairs = [
        split_pair(vector, s2_support, t2_support) for vector in kernel
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
        poisson({monomial: field.one}, T) for monomial in s4_support
    ]
    correction_five += [
        poisson(S, {monomial: field.one}) for monomial in t4_support
    ]
    strong_columns = correction_five + lower_variations
    output_support = sorted(
        set(constant).union(*(set(column) for column in strong_columns))
    )
    return {
        "S": S,
        "T": T,
        "rank_three": rank_three,
        "base_pair": base_pair,
        "kernel_pairs": kernel_pairs,
        "correction_five": correction_five,
        "lower_variations": lower_variations,
        "strong_columns": strong_columns,
        "constant": constant,
        "output_support": output_support,
    }


def family_presentation(field, sigma, tau) -> dict[str, object]:
    """Build the degree-seven inherited order-three/order-five data."""

    S, T = degree_seven_pair(field, sigma, tau)
    return relative_presentation(
        field,
        S,
        T,
        S2_SUPPORT,
        T2_SUPPORT,
        S4_SUPPORT,
        T4_SUPPORT,
    )


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


def pairing(functional, vector, field):
    """Pair one sparse output functional with one sparse defect vector."""

    return sum(
        (
            coefficient * vector.get(monomial, field.zero)
            for monomial, coefficient in functional.items()
        ),
        field.zero,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    assert tuple(map(len, (S2_SUPPORT, T2_SUPPORT))) == (72, 54)
    assert tuple(map(len, (S4_SUPPORT, T4_SUPPORT))) == (38, 25)

    # This rational point is enough to certify that the generic augmented
    # rank is at least 69.  The modular probes guard against bad reduction and
    # provide the baseline for the exhaustive Fitting scans.
    rational = rank_record(QQ, QQ.one, QQ.zero)
    expected = {
        "h3_rank": 118,
        "h3_kernel_dimension": 8,
        "h5_correction_rank": 61,
        "h5_strong_span_rank": 68,
        "h5_augmented_rank": 69,
        "h5_output_dimension": 178,
    }
    assert rational == expected

    prime_records = {}
    for prime in GOOD_PRIMES:
        field = GF(prime)
        record = rank_record(field, field.one, field.zero)
        assert record == expected
        prime_records[str(prime)] = record

    certificate = {
        "scope": (
            "parity-preserving inherited root-weight filtration; "
            "no unrestricted quantization or DC_2 claim"
        ),
        "family": {
            "degree": 7,
            "kappa": -9,
            "parameters": ["sigma", "tau"],
            "factor": (
                "w^4+sigma*w^3+tau*w^2"
                "+(-13/2-3*sigma-2*tau)*w+9/2+2*sigma+tau"
            ),
            "completing_shear": (
                "3*(26104*sigma^2+21736*sigma*tau+134160*sigma"
                "+4576*tau^2+56160*tau+75285)/28028"
            ),
        },
        "classical_symbol_bounds": {
            "S": {"z_order": 7, "bernstein_degree": 43},
            "T": {"z_order": 6, "bernstein_degree": 39},
        },
        "supports": {
            "S2": len(S2_SUPPORT),
            "T2": len(T2_SUPPORT),
            "S4": len(S4_SUPPORT),
            "T4": len(T4_SUPPORT),
        },
        "rational_point": {"sigma": 1, "tau": 0, "ranks": rational},
        "good_prime_records": prime_records,
        "generic_rank_lower_bound": expected,
        "conclusion": (
            "the selected two-parameter family is obstructed at order five "
            "on a nonempty open in the declared filtration"
        ),
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    print("PASS: exact degree-seven pair is canonical at a rational point")
    print("PASS: order-three rank 118 leaves an eight-dimensional torsor")
    print("PASS: order-five ranks are 61, 68, and 69")
    print("PASS: the same signature holds at five good primes")
    print("SCOPE: restricted filtration only; no unrestricted or DC_2 claim")


if __name__ == "__main__":
    main()
