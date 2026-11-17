#!/usr/bin/env python3
"""Eliminate the Q-fibre variable s6 before applying the dense s3 pivot.

On the residual Q component, mu4 is quadratic and mu5 is cubic in the
rescaled fibre variable s6.  Writing

    mu4 = a*s6^2 + f1*s6 + f0,
    mu5 = d*s6^3 + g2*s6^2 + g1*s6 + g0,

two division-free pseudo-remainder steps give

    r2 = a*mu5 - d*s6*mu4,
    L  = a*r2 - coeff(r2,s6^2)*mu4 = P*s6 + Q.

The s6-resultant is then

    H = a*Q^2 - f1*P*Q + f0*P^2.

This ordering deliberately keeps s3 polynomial.  Applying the dense
linear pivot s3=-B/(6A) before these steps causes a much larger
intermediate expression.  The output is a modular research profile, not
a characteristic-zero certificate or a treatment of the exceptional
leading-coefficient loci.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

from flint import nmod_mpoly_ctx

from research_two_pair_sic_bidegree33_t0_Q_flint import (
    ROOT,
    configure_coefficient_ring,
    flint_input_data,
    parse_moments,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_raw_prs_mod1000003.json"
)
RAW_VARIABLES = ("T", "s1", "lam", "v", "s3", "s5")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1000003)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--through",
        choices=(
            "linear",
            "resultant",
            "border",
            "subresultant",
            "escape",
        ),
        default="linear",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def profile(value) -> dict[str, object]:
    return {
        "terms": len(value.to_dict()),
        "degrees": {
            variable: int(degree)
            for variable, degree in zip(
                RAW_VARIABLES,
                value.degrees(),
                strict=True,
            )
        },
        "total_degree": int(value.total_degree()),
        "sha256": hashlib.sha256(str(value).encode()).hexdigest(),
    }


def coefficient_in_s6(raw, degree: int, context):
    terms = {}
    for (s6_power, s5_power, s3_power), coefficient in raw.items():
        if s6_power != degree:
            continue
        for exponent, scalar in coefficient.to_dict().items():
            lifted_exponent = exponent + (s3_power, s5_power)
            terms[lifted_exponent] = (
                terms.get(lifted_exponent, 0) + int(scalar)
            )
    return context.from_dict(terms)


def variable_degree(value, index: int) -> int:
    return max(
        (exponent[index] for exponent in value.to_dict()),
        default=-1,
    )


def variable_coefficient(value, index: int, degree: int, context):
    terms = {}
    for exponent, scalar in value.to_dict().items():
        if exponent[index] != degree:
            continue
        coefficient_exponent = (
            exponent[:index] + (0,) + exponent[index + 1 :]
        )
        terms[coefficient_exponent] = (
            terms.get(coefficient_exponent, 0) + int(scalar)
        )
    return context.from_dict(terms)


def pseudo_remainder(value, divisor, variable, variable_index, context):
    """Division-free remainder in one selected variable."""

    divisor_degree = variable_degree(divisor, variable_index)
    divisor_leading = variable_coefficient(
        divisor,
        variable_index,
        divisor_degree,
        context,
    )
    remainder = value
    steps = 0
    while variable_degree(remainder, variable_index) >= divisor_degree:
        remainder_degree = variable_degree(remainder, variable_index)
        remainder_leading = variable_coefficient(
            remainder,
            variable_index,
            remainder_degree,
            context,
        )
        shift = remainder_degree - divisor_degree
        remainder = (
            divisor_leading * remainder
            - remainder_leading * variable**shift * divisor
        )
        steps += 1
    return remainder, steps


def main() -> None:
    arguments = parse_arguments()
    if arguments.prime in (0, 2, 3, 5, 7, 13):
        raise ValueError("choose a prime avoiding the displayed denominators")
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")

    started = time.monotonic()
    configure_coefficient_ring(arguments.prime)
    _, _, _, input_metadata = flint_input_data()
    moments, moment_profiles = parse_moments(singular, arguments.timeout)
    context = nmod_mpoly_ctx.get(
        RAW_VARIABLES,
        modulus=arguments.prime,
    )
    f2, f1, f0 = (
        coefficient_in_s6(moments[4], degree, context)
        for degree in (2, 1, 0)
    )
    g3, g2, g1, g0 = (
        coefficient_in_s6(moments[5], degree, context)
        for degree in (3, 2, 1, 0)
    )

    linear_started = time.monotonic()
    r2_2 = f2 * g2 - g3 * f1
    r2_1 = f2 * g1 - g3 * f0
    r2_0 = f2 * g0
    linear_p = f2 * r2_1 - r2_2 * f1
    linear_q = f2 * r2_0 - r2_2 * f0
    linear_seconds = time.monotonic() - linear_started

    payload: dict[str, object] = {
        "format": "two-pair-sic-bidegree33-t0-Q-raw-prs-v1",
        "status": (
            f"exact finite-field calculation modulo {arguments.prime}; "
            "not a characteristic-zero certificate"
        ),
        "prime": arguments.prime,
        "variables": list(RAW_VARIABLES),
        "through": arguments.through,
        "input": input_metadata,
        "moments": moment_profiles,
        "quadratic_coefficients": {
            "a": profile(f2),
            "f1": profile(f1),
            "f0": profile(f0),
        },
        "cubic_coefficients": {
            "d": profile(g3),
            "g2": profile(g2),
            "g1": profile(g1),
            "g0": profile(g0),
        },
        "first_pseudo_remainder": {
            "s6^2": profile(r2_2),
            "s6": profile(r2_1),
            "constant": profile(r2_0),
        },
        "linear_relation": {
            "P": profile(linear_p),
            "Q": profile(linear_q),
            "seconds": round(linear_seconds, 6),
            "meaning": "P*s6+Q belongs to (mu4,mu5)",
        },
    }
    if arguments.through in (
        "resultant",
        "border",
        "subresultant",
        "escape",
    ):
        resultant_started = time.monotonic()
        resultant = (
            f2 * linear_q**2
            - f1 * linear_p * linear_q
            + f0 * linear_p**2
        )
        payload["resultant"] = {
            **profile(resultant),
            "seconds": round(
                time.monotonic() - resultant_started,
                6,
            ),
            "meaning": "Res_s6(mu4,mu5) before content removal",
        }
        if arguments.through in ("border", "subresultant", "escape"):
            border_started = time.monotonic()
            border = linear_p.resultant(resultant, "s5")
            payload["border"] = {
                **profile(border),
                "seconds": round(
                    time.monotonic() - border_started,
                    6,
                ),
                "meaning": (
                    "Res_s5(P,Res_s6(mu4,mu5)); away from this "
                    "polynomial, P is invertible in the length-five "
                    "fibre algebra. This is an additional "
                    "P-invertibility locus, not the canonical "
                    "leading-coefficient Q-border Delta."
                ),
            }
            if arguments.through in ("subresultant", "escape"):
                subresultant_started = time.monotonic()
                bad_root_relation, steps = pseudo_remainder(
                    resultant,
                    linear_p,
                    context.gens()[5],
                    5,
                    context,
                )
                payload["linear_subresultant"] = {
                    **profile(bad_root_relation),
                    "steps": steps,
                    "seconds": round(
                        time.monotonic() - subresultant_started,
                        6,
                    ),
                    "meaning": (
                        "division-free remainder of H by P; on the "
                        "additional P-invertibility locus its linear "
                        "factor identifies the root that must be "
                        "checked for escape to infinity"
                    ),
                }
                if arguments.through == "escape":
                    escape_started = time.monotonic()
                    escape_value, escape_steps = pseudo_remainder(
                        linear_q,
                        bad_root_relation,
                        context.gens()[5],
                        5,
                        context,
                    )
                    if variable_degree(escape_value, 5) > 0:
                        raise AssertionError(
                            "linear escape pseudo-remainder retained s5"
                        )
                    payload["escape_test"] = {
                        **profile(escape_value),
                        "steps": escape_steps,
                        "seconds": round(
                            time.monotonic() - escape_started,
                            6,
                        ),
                        "meaning": (
                            "fraction-free resultant of the bad-root "
                            "linear relation with Q; nonvanishing "
                            "excludes that root from the affine "
                            "equation P*s6+Q=0"
                        ),
                    }
    payload["seconds"] = round(time.monotonic() - started, 6)
    payload["reproduction_command"] = " ".join(sys.argv)

    output = arguments.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
