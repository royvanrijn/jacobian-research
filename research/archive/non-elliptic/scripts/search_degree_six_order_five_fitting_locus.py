#!/usr/bin/env python3
"""Modular discovery scan for the degree-six order-five Fitting locus.

This script specializes the exact relative presentation from
``verify_degree_six_relative_quantization_obstruction.py`` at every point of
one finite parameter plane.  It records:

* failure or rank change of the order-three lift equation;
* rank drop of the order-five strong matrix; and
* points where adjoining the constant order-five defect does not raise rank.

The last condition is the fiberwise consistency condition for the strong
order-five system.  It is only a modular component-discovery screen.  In
particular, points found over one prime are not characteristic-zero
components, and the degree-drop divisor ``sigma=0`` is kept separate.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

from sympy.polys.domains import GF

from verify_degree_six_relative_quantization_obstruction import rank_record


GENERIC_SIGNATURE = (77, 6, 34, 41, 42)


def reconstructed_rank_drop_divisor(
    sigma_value: int,
    tau_value: int,
    prime: int,
) -> int:
    """Seven-prime reconstruction candidate for the rank-40 divisor.

    The primitive polynomial is normalized to have ``sigma^4`` coefficient
    2563.  Using it to restrict a scan is an optimization for modular
    discovery, not a characteristic-zero determinantal certificate.
    """

    sigma = sigma_value
    tau = tau_value
    value = (
        2563 * sigma**4
        + 3954 * sigma**3 * tau
        + 2319 * sigma**2 * tau**2
        + 608 * sigma * tau**3
        + 60 * tau**4
        + 7240 * sigma**3
        + 7200 * sigma**2 * tau
        + 2280 * sigma * tau**2
        + 200 * tau**3
        + 6970 * sigma**2
        + 3400 * sigma * tau
        + 250 * tau**2
        + 2250 * sigma
    )
    return value % prime


def point_record(task: tuple[int, int, int]) -> dict[str, object]:
    """Compute one exact finite-field rank record."""

    prime, sigma_value, tau_value = task
    field = GF(prime)
    try:
        ranks = rank_record(
            field,
            field(sigma_value),
            field(tau_value),
        )
    except (ArithmeticError, AssertionError, ValueError, ZeroDivisionError) as error:
        return {
            "sigma": sigma_value,
            "tau": tau_value,
            "error": f"{type(error).__name__}: {error}",
        }

    signature = (
        ranks["h3_rank"],
        ranks["h3_kernel_dimension"],
        ranks["h5_correction_rank"],
        ranks["h5_strong_span_rank"],
        ranks["h5_augmented_rank"],
    )
    return {
        "sigma": sigma_value,
        "tau": tau_value,
        "degree_drop": sigma_value == 0,
        "signature": list(signature),
        "strong_rank_drop": ranks["h5_strong_span_rank"] < 41,
        "order_five_consistent": (
            ranks["h5_augmented_rank"]
            == ranks["h5_strong_span_rank"]
        ),
    }


def scan_prime(
    prime: int,
    jobs: int,
    restrict_to_candidate_divisor: bool,
) -> dict[str, object]:
    tasks = [
        (prime, sigma_value, tau_value)
        for sigma_value in range(prime)
        for tau_value in range(prime)
        if not restrict_to_candidate_divisor
        or reconstructed_rank_drop_divisor(
            sigma_value,
            tau_value,
            prime,
        )
        == 0
    ]
    if jobs == 1:
        records = [point_record(task) for task in tasks]
    else:
        context = mp.get_context("spawn")
        with context.Pool(processes=jobs) as pool:
            records = pool.map(point_record, tasks, chunksize=4)

    errors = [record for record in records if "error" in record]
    exceptional = [
        record
        for record in records
        if "signature" in record
        and tuple(record["signature"]) != GENERIC_SIGNATURE
    ]
    consistent = [
        record
        for record in records
        if record.get("order_five_consistent")
    ]
    interior_consistent = [
        record for record in consistent if not record["degree_drop"]
    ]
    return {
        "prime": prime,
        "parameter_points": len(tasks),
        "restricted_to_reconstructed_rank_drop_divisor": (
            restrict_to_candidate_divisor
        ),
        "generic_signature": list(GENERIC_SIGNATURE),
        "errors": errors,
        "exceptional_points": exceptional,
        "order_five_consistent_points": consistent,
        "interior_order_five_consistent_points": interior_consistent,
        "counts": {
            "errors": len(errors),
            "exceptional": len(exceptional),
            "order_five_consistent": len(consistent),
            "interior_order_five_consistent": len(interior_consistent),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--candidate-divisor",
        action="store_true",
        help=(
            "scan only points on the seven-prime reconstructed rank-drop "
            "quartic"
        ),
    )
    args = parser.parse_args()

    if args.prime <= 7:
        raise SystemExit("choose a prime greater than 7")
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive")

    result = scan_prime(args.prime, args.jobs, args.candidate_divisor)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    counts = result["counts"]
    print(
        "PASS: scanned GF({})^2; exceptional={}, consistent={}, "
        "interior-consistent={}".format(
            args.prime,
            counts["exceptional"],
            counts["order_five_consistent"],
            counts["interior_order_five_consistent"],
        )
    )
    for record in result["interior_order_five_consistent_points"]:
        print(
            "SURVIVOR: sigma={} tau={} signature={}".format(
                record["sigma"],
                record["tau"],
                record["signature"],
            )
        )
    print(
        "SCOPE: modular discovery only; no characteristic-zero component "
        "is certified"
    )


if __name__ == "__main__":
    main()
