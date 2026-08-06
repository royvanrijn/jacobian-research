#!/usr/bin/env python3
"""Restricted order-five presentation for a degree-nine marked-root slice."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sympy.polys.domains import GF, QQ

from explore_degree_five_quantum_residue import (
    SparsePoly,
    column_rank,
    weighted_seed_symbol_pair,
)
from verify_degree_seven_relative_quantization_obstruction import (
    relative_presentation,
    weight_monomials,
)


GOOD_PRIMES = (23, 29, 31, 37)


def degree_nine_shear(field, sigma, tau):
    numerator = field(3) * (
        field(2437664) * sigma**2
        + field(3435360) * sigma * tau
        + field(8147828) * sigma
        + field(1214208) * tau**2
        + field(5770548) * tau
        + field(5171465)
    )
    return numerator / field(476476)


def degree_nine_pair(field, sigma, tau) -> tuple[SparsePoly, SparsePoly]:
    factor_coefficients = (
        field(13) / field(2) + field(4) * sigma + field(3) * tau,
        -field(17) / field(2) - field(5) * sigma - field(4) * tau,
        field.zero,
        field.zero,
        tau,
        sigma,
        field.one,
    )
    return weighted_seed_symbol_pair(
        field,
        -field(8) / field(7),
        factor_coefficients,
        degree_nine_shear(field, sigma, tau),
    )


S2_SUPPORT = weight_monomials(53, 7, 4)
T2_SUPPORT = weight_monomials(49, 6, 5)
S4_SUPPORT = weight_monomials(49, 5, 10)
T4_SUPPORT = weight_monomials(45, 4, 11)
S6_SUPPORT = weight_monomials(45, 3, 16)
T6_SUPPORT = weight_monomials(41, 2, 17)


def family_presentation(field, sigma, tau):
    S, T = degree_nine_pair(field, sigma, tau)
    return relative_presentation(
        field,
        S,
        T,
        S2_SUPPORT,
        T2_SUPPORT,
        S4_SUPPORT,
        T4_SUPPORT,
    )


def rank_record(field, sigma, tau):
    presentation = family_presentation(field, sigma, tau)
    correction = presentation["correction_five"]
    strong = presentation["strong_columns"]
    constant = presentation["constant"]
    correction_rank = column_rank(correction)
    return {
        "h3_rank": presentation["rank_three"],
        "h3_kernel_dimension": len(presentation["kernel_pairs"]),
        "h5_correction_rank": correction_rank,
        "h5_correction_kernel_dimension": len(correction) - correction_rank,
        "h5_strong_span_rank": column_rank(strong),
        "h5_augmented_rank": column_rank(strong + [constant]),
        "h5_output_dimension": len(presentation["output_support"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    assert tuple(map(len, (S2_SUPPORT, T2_SUPPORT))) == (132, 107)
    assert tuple(map(len, (S4_SUPPORT, T4_SUPPORT))) == (84, 64)
    assert tuple(map(len, (S6_SUPPORT, T6_SUPPORT))) == (46, 31)
    expected = {
        "h3_rank": 227,
        "h3_kernel_dimension": 12,
        "h5_correction_rank": 142,
        "h5_correction_kernel_dimension": 6,
        "h5_strong_span_rank": 149,
        "h5_augmented_rank": 150,
        "h5_output_dimension": 371,
    }
    rational = rank_record(QQ, QQ.one, QQ.zero)
    assert rational == expected
    prime_records = {}
    for prime in GOOD_PRIMES:
        field = GF(prime)
        record = rank_record(field, field.one, field.zero)
        assert record == expected
        prime_records[str(prime)] = record
    bad_field = GF(19)
    bad_reduction = rank_record(bad_field, bad_field.one, bad_field.zero)
    assert bad_reduction["h3_rank"] == 226

    certificate = {
        "scope": (
            "degree-nine rank discovery in the parity-preserving inherited "
            "root-weight filtration; no vanishing-locus or DC_2 claim"
        ),
        "family": {
            "degree": 9,
            "kappa": -9,
            "parameters": ["sigma", "tau"],
            "factor": (
                "w^6+sigma*w^5+tau*w^4"
                "+(-17/2-5*sigma-4*tau)*w+13/2+4*sigma+3*tau"
            ),
            "completing_shear": (
                "3*(2437664*sigma^2+3435360*sigma*tau+8147828*sigma"
                "+1214208*tau^2+5770548*tau+5171465)/476476"
            ),
        },
        "supports": {
            "S2": len(S2_SUPPORT),
            "T2": len(T2_SUPPORT),
            "S4": len(S4_SUPPORT),
            "T4": len(T4_SUPPORT),
            "S6": len(S6_SUPPORT),
            "T6": len(T6_SUPPORT),
        },
        "rational_point": {"sigma": 1, "tau": 0, "ranks": rational},
        "good_prime_records": prime_records,
        "excluded_bad_reduction": {"19": bad_reduction},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print("PASS: exact degree-nine ranks extend the degree-six-through-eight pattern")
    print("RATIONAL:", rational)
    for prime, record in prime_records.items():
        print(f"GF({prime}):", record)
    print("SCOPE: rank discovery only; scan the obstruction section next")


if __name__ == "__main__":
    main()
