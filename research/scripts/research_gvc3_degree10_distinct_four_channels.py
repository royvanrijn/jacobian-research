#!/usr/bin/env python3
"""Modular frontier for pairwise-distinct degree-ten four-channel profiles.

For balanced degree ten, the positive even harmonic degrees are
``2,4,6,8,10``.  This script checks every four-of-five profile on the
pairwise-distinct direction and coefficient torus.  It compiles invariant
Reynolds moments and uses quotient saturation modulo several primes.

The output is discovery only.  A unit ideal in several finite
characteristics is not promoted to a characteristic-zero theorem.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import sympy as sp

from research_gvc3_many_coherent_channels import (
    configuration_discriminant,
    moment,
    modular_saturation_cutoff,
    primitive_polynomial,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_degree10_distinct_four_channels_modular.json"
)
PROFILES = tuple(itertools.combinations((2, 4, 6, 8, 10), 4))
GROUPS = ((0,), (1,), (2,), (3,))
PRIMES = (101, 103, 107)


def compile_profile(degrees: tuple[int, ...], max_order: int) -> dict[str, object]:
    coefficients = sp.symbols("a0:4")
    a0, a1, a2, a3 = coefficients
    lam = sp.Symbol("lam")
    variables = (a1, a2, a3, lam)
    equations: dict[int, sp.Expr] = {}
    term_counts: dict[str, int] = {}
    moment_sha256: dict[str, str] = {}
    for order in range(2, max_order + 1):
        expression = moment(degrees, GROUPS, order, coefficients, (lam,))
        polynomial = primitive_polynomial(expression.subs(a0, 1), variables)
        term_counts[str(order)] = (
            len(sp.Poly(polynomial, *variables).terms())
            if polynomial != 0
            else 0
        )
        if polynomial != 0:
            equations[order] = polynomial
            moment_sha256[str(order)] = hashlib.sha256(
                str(polynomial).encode()
            ).hexdigest()
    return {
        "degrees": degrees,
        "variables": variables,
        "equations": equations,
        "saturation": a1
        * a2
        * a3
        * configuration_discriminant(4, (lam,)),
        "moment_term_counts": term_counts,
        "moment_sha256": moment_sha256,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=7)
    parser.add_argument("--primes", nargs="+", type=int, default=list(PRIMES))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    compiled = [
        compile_profile(profile, arguments.max_order) for profile in PROFILES
    ]
    tasks = {}
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        for profile_index, data in enumerate(compiled):
            for prime in arguments.primes:
                future = executor.submit(
                    modular_saturation_cutoff,
                    data["equations"],
                    data["saturation"],
                    data["variables"],
                    prime,
                    arguments.singular,
                    arguments.timeout,
                )
                tasks[future] = (profile_index, prime)
        results: list[list[dict[str, object]]] = [
            [] for _ in compiled
        ]
        for future in as_completed(tasks):
            profile_index, prime = tasks[future]
            result = future.result()
            results[profile_index].append(result)
            print(
                "DONE",
                ",".join(map(str, compiled[profile_index]["degrees"])),
                f"p={prime}",
                f"status={result.get('status')}",
                f"unit={result.get('unit')}",
                f"seconds={result.get('seconds', 0):.3f}",
                flush=True,
            )

    profiles = []
    for data, profile_results in zip(compiled, results, strict=True):
        profile_results.sort(key=lambda result: int(result["prime"]))
        all_units = all(
            result.get("status") == "completed" and result.get("unit") == 1
            for result in profile_results
        )
        profiles.append(
            {
                "harmonic_degrees": list(data["degrees"]),
                "moment_term_counts": data["moment_term_counts"],
                "moment_sha256": data["moment_sha256"],
                "modular_saturations": profile_results,
                "all_declared_primes_unit": all_units,
            }
        )

    all_profiles_unit = all(
        profile["all_declared_primes_unit"] for profile in profiles
    )
    artifact = {
        "format": "gvc3-degree10-distinct-four-channels-modular-v1",
        "status": "bounded modular discovery; not a characteristic-zero theorem",
        "balanced_degree": 10,
        "laplacian_power": 5,
        "scope": (
            "all five four-of-five positive-even harmonic profiles, with "
            "pairwise-distinct directions and nonzero coefficients"
        ),
        "normalization": "directions infinity, zero, one, lam; a0=1",
        "configuration_saturation": "lam*(lam-1)",
        "coefficient_saturation": "a1*a2*a3",
        "max_moment_order": arguments.max_order,
        "primes": list(arguments.primes),
        "workers": arguments.workers,
        "per_saturation_timeout_seconds": arguments.timeout,
        "saturation_method": "successive ideal quotients (Singular sat_with_exp)",
        "moment_sequence_promotion_gate": {
            "type": "finite proper-hypergeometric Wick/occupation sum",
            "reference": "gvc3_four_coherent_channels.json",
        },
        "profiles": profiles,
        "all_profiles_unit_at_all_declared_primes": all_profiles_unit,
        "conclusion": (
            "every pairwise-distinct degree-ten four-channel profile is a modular unit candidate through moment seven"
            if all_profiles_unit
            else "at least one declared profile/prime did not return a unit within the declared run"
        ),
        "not_checked_here": [
            "characteristic-zero promotion",
            "direction-collision strata",
            "the five-channel degree-ten profile",
            "multiplier survival, because no pure-moment survivor was produced",
        ],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    )
    if all_profiles_unit:
        print("PASS all five profiles are units at all declared primes")
    else:
        print("OPEN at least one profile/prime did not return a unit")


if __name__ == "__main__":
    main()
