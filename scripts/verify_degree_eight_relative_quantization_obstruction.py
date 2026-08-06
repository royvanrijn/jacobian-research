#!/usr/bin/env python3
"""Restricted order-five presentation for a degree-eight marked-root slice."""

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


GOOD_PRIMES = (19, 23, 29, 31)


def degree_eight_shear(field, sigma, tau):
    numerator = field(3) * (
        field(71424) * sigma**2
        + field(86112) * sigma * tau
        + field(279300) * sigma
        + field(26104) * tau**2
        + field(169352) * tau
        + field(176269)
    )
    return numerator / field(28028)


def degree_eight_pair(field, sigma, tau) -> tuple[SparsePoly, SparsePoly]:
    a = -field(8) / field(7)
    factor_coefficients = (
        field(11) / field(2) + field(3) * sigma + field(2) * tau,
        -field(15) / field(2) - field(4) * sigma - field(3) * tau,
        field.zero,
        tau,
        sigma,
        field.one,
    )
    return weighted_seed_symbol_pair(
        field,
        a,
        factor_coefficients,
        degree_eight_shear(field, sigma, tau),
    )


# For degree n, the classical bounds are
# S:(deg_Z,deg_B)=(n,7n-6), T:(n-1,7n-10).  Each even PBW step lowers
# deg_Z by two, deg_B by four, and raises root weight by six.
S2_SUPPORT = weight_monomials(46, 6, 4)
T2_SUPPORT = weight_monomials(42, 5, 5)
S4_SUPPORT = weight_monomials(42, 4, 10)
T4_SUPPORT = weight_monomials(38, 3, 11)


def family_presentation(field, sigma, tau):
    S, T = degree_eight_pair(field, sigma, tau)
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
    return {
        "h3_rank": presentation["rank_three"],
        "h3_kernel_dimension": len(presentation["kernel_pairs"]),
        "h5_correction_rank": column_rank(correction),
        "h5_correction_kernel_dimension": len(correction)
        - column_rank(correction),
        "h5_strong_span_rank": column_rank(strong),
        "h5_augmented_rank": column_rank(strong + [constant]),
        "h5_output_dimension": len(presentation["output_support"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    assert tuple(map(len, (S2_SUPPORT, T2_SUPPORT))) == (100, 78)
    assert tuple(map(len, (S4_SUPPORT, T4_SUPPORT))) == (59, 42)
    rational = rank_record(QQ, QQ.one, QQ.zero)
    expected = {
        "h3_rank": 168,
        "h3_kernel_dimension": 10,
        "h5_correction_rank": 97,
        "h5_correction_kernel_dimension": 4,
        "h5_strong_span_rank": 104,
        "h5_augmented_rank": 105,
        "h5_output_dimension": 269,
    }
    assert rational == expected
    prime_records = {}
    for prime in GOOD_PRIMES:
        field = GF(prime)
        record = rank_record(
            field, field.one, field.zero
        )
        assert record == expected
        prime_records[str(prime)] = record
    bad_field = GF(17)
    bad_reduction = rank_record(bad_field, bad_field.one, bad_field.zero)
    assert bad_reduction["h3_rank"] == 167
    certificate = {
        "scope": (
            "degree-eight rank discovery in the parity-preserving inherited "
            "root-weight filtration; no vanishing-locus or DC_2 claim"
        ),
        "family": {
            "degree": 8,
            "kappa": -9,
            "parameters": ["sigma", "tau"],
            "factor": (
                "w^5+sigma*w^4+tau*w^3"
                "+(-15/2-4*sigma-3*tau)*w+11/2+3*sigma+2*tau"
            ),
            "completing_shear": (
                "3*(71424*sigma^2+86112*sigma*tau+279300*sigma"
                "+26104*tau^2+169352*tau+176269)/28028"
            ),
        },
        "supports": {
            "S2": len(S2_SUPPORT),
            "T2": len(T2_SUPPORT),
            "S4": len(S4_SUPPORT),
            "T4": len(T4_SUPPORT),
        },
        "rational_point": {"sigma": 1, "tau": 0, "ranks": rational},
        "good_prime_records": prime_records,
        "excluded_bad_reduction": {"17": bad_reduction},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print("PASS: exact degree-eight ranks match the degree-six/seven pattern")
    print("RATIONAL:", rational)
    for prime, record in prime_records.items():
        print(f"GF({prime}):", record)
    print("SCOPE: rank discovery only; scan the obstruction section next")


if __name__ == "__main__":
    main()
