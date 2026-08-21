#!/usr/bin/env python3
"""Pin two exact rank-19 split-infinity Mestre frontiers.

The two specializations were found by bounded point searches.  This script
replays those cached searches (or reruns them), verifies the exact mod-3
independence certificates and exact PARI conductor data, and evaluates the
audited Delta=11/5 explicit-formula diagnostic.  The rank lower bounds are
unconditional.  The analytic upper bounds are conditional on GRH.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import mpmath
from pathlib import Path
import shutil
from typing import Any

from explicit_formula_mestre_rank15_490_9_delta22 import (
    DELTA,
    NUMERICAL_ALLOWANCE,
    run_prime_sum,
)
from explicit_formula_rank20_t5081_delta22 import (
    archimedean_closed_form,
    explicit_formula_upper,
)
import search_mestre_dsquare_four as search


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elliptic_mestre_dsquare_rank19_frontiers.json"
)
FRONTIERS = (
    {
        "label": "family2_u483",
        "family_index": 2,
        "numerator": 483,
        "denominator": 1,
        "height": 100_000,
        "denominator_bound": 2_000,
        "raw_root": ROOT / "artifacts/local/elliptic-curves/mestre-family2-integer-pilot-h100k/ratpoints-raw",
        "expected_basis_sha256": "cfcdf33e95f20b448c822668140331b5b9e16597a62e8f82c0fd3a63b4d7c176",
        "expected_upper_interval": ("20.30", "20.31"),
    },
    {
        "label": "family3_u660",
        "family_index": 3,
        "numerator": 660,
        "denominator": 1,
        "height": 2_000_000,
        "denominator_bound": 13_000,
        "raw_root": ROOT / "artifacts/local/elliptic-curves/mestre-family3-full-h2m/ratpoints-raw",
        "expected_basis_sha256": "6ea2924633ed408bf9409024bfc71fe66e6b9e8ac13fe90233cf8d068f734c4d",
        "expected_upper_interval": ("20.41", "20.42"),
    },
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def build_record(spec: dict[str, Any], *, timeout: float) -> dict[str, Any]:
    candidate = {
        "family_index": spec["family_index"],
        "numerator": spec["numerator"],
        "denominator": spec["denominator"],
    }
    search.POINT_HEIGHT = spec["height"]
    search.POINT_DENOMINATOR_BOUND = spec["denominator_bound"]
    result = search.point_worker(candidate, str(spec["raw_root"]), timeout)
    certificate = result["finite_reduction_certificate"]
    if certificate["combined_exact_rank_over_F3"] != 19:
        raise AssertionError(f"{spec['label']} lost its exact rank-19 certificate")
    if certificate["independent_subset_sha256"] != spec["expected_basis_sha256"]:
        raise AssertionError(f"{spec['label']} changed its independent basis")

    global_result = search.conductor_worker(candidate, timeout)
    global_curve = global_result["global_curve"]
    if not global_result["below_strict_log_conductor_182_72"]:
        raise AssertionError(f"{spec['label']} crossed the conductor target")
    if global_curve["root_number"] != -1:
        raise AssertionError(f"{spec['label']} changed root-number parity")

    prime_sum = run_prime_sum(
        global_curve["minimal_model"], timeout=timeout, stack_bytes=search.STACK_BYTES
    )
    if prime_sum["conductor"] != global_curve["conductor"]:
        raise AssertionError(f"{spec['label']} conductor replay changed")
    with mpmath.workdps(60):
        _, archimedean = archimedean_closed_form(DELTA)
        upper = explicit_formula_upper(
            log_conductor=str(prime_sum["log_conductor"]),
            prime_sum=str(prime_sum["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=archimedean,
            numerical_allowance=NUMERICAL_ALLOWANCE,
        )
        lower_expected, upper_expected = map(mpmath.mpf, spec["expected_upper_interval"])
        if not lower_expected < upper < upper_expected or not upper < 21:
            raise AssertionError(f"{spec['label']} explicit-formula gate changed")
        upper_text = mpmath.nstr(upper, 60)

    family = search.FAMILIES[spec["family_index"]]
    parameter_u = Q(spec["numerator"], spec["denominator"])
    parameter_t = search.base_parameter(family, parameter_u)
    coefficients = family.construction.primitive_jacobian_coefficients(parameter_t)
    point_search = dict(result["point_search"])
    point_search["charts"] = [
        {
            key: value
            for key, value in chart.items()
            if key not in ("cached", "wall_seconds")
        }
        for chart in point_search["charts"]
    ]
    stable_prime_sum = {
        key: value
        for key, value in prime_sum.items()
        if key != "pari_milliseconds_after_curve_setup"
    }
    return {
        "label": spec["label"],
        "family_index": spec["family_index"],
        "roots": list(family.roots),
        "u": rational_text(parameter_u),
        "T": rational_text(parameter_t),
        "short_weierstrass_coefficients": [rational_text(Q(value)) for value in coefficients],
        "minimal_model": global_curve["minimal_model"],
        "conductor": global_curve["conductor"],
        "log_conductor": global_curve["log_conductor"],
        "root_number": global_curve["root_number"],
        "point_search": point_search,
        "pool_point_count_modulo_inverse": result["pool_point_count_modulo_inverse"],
        "pool_point_sha256": result["pool_point_sha256"],
        "exact_rank_certificate": certificate,
        "explicit_formula_delta_11_over_5": {
            **stable_prime_sum,
            "numerical_allowance": str(NUMERICAL_ALLOWANCE),
            "conservative_upper": upper_text,
            "under_grh_analytic_rank_at_most_19": True,
            "parity_reason": "root number -1 forces odd analytic order and the upper is below 21",
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--check", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    ratpoints = shutil.which("ratpoints")
    if ratpoints is None:
        raise SystemExit("ratpoints is required")
    search.RATPOINTS = Path(ratpoints)
    records = [build_record(spec, timeout=args.timeout) for spec in FRONTIERS]
    payload = {
        "schema_version": 1,
        "status": "two exact low-conductor rank-19 frontiers certified",
        "claim": {
            "unconditional": "each displayed curve has algebraic rank at least 19 and log conductor below 182.72",
            "conditional_on_grh": "each displayed curve has analytic rank at most 19",
            "conditional_on_grh_and_bsd": "each displayed curve has algebraic rank exactly 19",
            "rank21_target_achieved": False,
        },
        "frontiers": records,
    }
    payload["result_sha256"] = search.canonical_digest(payload)
    if args.check:
        existing = json.loads(args.output.read_text())
        if existing.get("result_sha256") != payload["result_sha256"]:
            raise AssertionError("the pinned rank-19 frontier artifact changed")
        print(f"PASS {args.output} sha256={payload['result_sha256']}")
    else:
        search.atomic_json(args.output, payload)
        print(f"WROTE {args.output} sha256={payload['result_sha256']}")


if __name__ == "__main__":
    main()
