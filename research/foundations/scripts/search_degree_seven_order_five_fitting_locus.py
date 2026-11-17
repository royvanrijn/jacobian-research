#!/usr/bin/env python3
"""Modular Fitting and consistency scan for degree-seven/nine marked slices.

Every point of one finite parameter plane is specialized into the complete
inherited order-five presentation.  The output separates strong-matrix rank
drop from fibrewise consistency after adjoining the constant defect.  This
is modular component discovery, not a characteristic-zero component proof.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

from sympy.polys.domains import GF

from verify_degree_seven_relative_quantization_obstruction import (
    rank_record as degree_seven_rank_record,
)
from verify_degree_eight_relative_quantization_obstruction import (
    rank_record as degree_eight_rank_record,
)
from verify_degree_nine_relative_quantization_obstruction import (
    rank_record as degree_nine_rank_record,
)


GENERIC_SIGNATURES = {
    7: (118, 8, 61, 68, 69),
    8: (168, 10, 97, 104, 105),
    9: (227, 12, 142, 149, 150),
}


RANK_RECORDS = {
    7: degree_seven_rank_record,
    8: degree_eight_rank_record,
    9: degree_nine_rank_record,
}


def point_record(task: tuple[int, int, int, int]) -> dict[str, object]:
    degree, prime, sigma_value, tau_value = task
    field = GF(prime)
    rank_record = RANK_RECORDS[degree]
    try:
        ranks = rank_record(field, field(sigma_value), field(tau_value))
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
        "signature": list(signature),
        "strong_rank_drop": (
            ranks["h5_strong_span_rank"] < GENERIC_SIGNATURES[degree][3]
        ),
        "order_five_consistent": (
            ranks["h5_augmented_rank"]
            == ranks["h5_strong_span_rank"]
        ),
    }


def scan_prime(degree: int, prime: int, jobs: int) -> dict[str, object]:
    tasks = [
        (degree, prime, sigma_value, tau_value)
        for sigma_value in range(prime)
        for tau_value in range(prime)
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
        and tuple(record["signature"]) != GENERIC_SIGNATURES[degree]
    ]
    consistent = [
        record for record in records if record.get("order_five_consistent")
    ]
    return {
        "prime": prime,
        "degree": degree,
        "parameter_points": len(tasks),
        "generic_signature": list(GENERIC_SIGNATURES[degree]),
        "errors": errors,
        "exceptional_points": exceptional,
        "order_five_consistent_points": consistent,
        "counts": {
            "errors": len(errors),
            "exceptional": len(exceptional),
            "order_five_consistent": len(consistent),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, choices=(7, 8, 9), default=7)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.prime <= 13:
        raise SystemExit("choose a prime greater than 13")
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive")

    result = scan_prime(args.degree, args.prime, args.jobs)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    counts = result["counts"]
    print(
        "PASS: scanned GF({})^2; exceptional={}, consistent={}".format(
            args.prime,
            counts["exceptional"],
            counts["order_five_consistent"],
        )
    )
    for record in result["order_five_consistent_points"]:
        print(
            "SURVIVOR: sigma={} tau={} signature={}".format(
                record["sigma"], record["tau"], record["signature"]
            )
        )
    print("SCOPE: modular discovery only; no characteristic-zero component")


if __name__ == "__main__":
    main()
