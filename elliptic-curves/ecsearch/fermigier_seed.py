"""First CRT discriminant-shaping seed in the Fermigier adapter family."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Any

from .crt_lattice import (
    crt,
    evaluate_polynomial,
    first_rationals_by_height,
    gauss_reduce_congruence_lattice,
    hensel_lift_roots,
    p_adic_valuation,
)
from .fermigier import (
    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
    fermigier_canonical_coefficients,
    fermigier_discriminant_factor,
)
from .local_data import weierstrass_local_data


CONSTRAINTS = (
    {"prime": 89, "exponent": 2},
    {"prime": 131, "exponent": 2},
    {"prime": 137, "exponent": 2},
)


def homogenized_discriminant_factor(numerator: int, denominator: int) -> int:
    """Return ``b^20*Phi(a/b)`` for the canonical adapter discriminant."""

    value = fermigier_discriminant_factor(Fraction(numerator, denominator))
    homogeneous = value * denominator**20
    if homogeneous.denominator != 1:
        raise ArithmeticError("homogenized factor is unexpectedly nonintegral")
    return homogeneous.numerator


def _format_fraction(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _pari_local_replay(
    numerator: int, denominator: int, primes: tuple[int, ...]
) -> dict[str, Any]:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP executable 'gp' is required")
    coefficients = fermigier_canonical_coefficients(
        Fraction(numerator, denominator)
    )
    model = ",".join(f"({_format_fraction(value)})" for value in coefficients)
    lines = ["setrand(1);", f"E=ellinit([{model}]);"]
    for prime in primes:
        lines.extend(
            (
                f"L=elllocalred(E,{prime});",
                f'print("LOCAL_{prime}=",L[1],",",L[2],",",L[4],",",valuation(E.disc,{prime}),",",ellap(E,{prime}));',
            )
        )
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "***" in combined:
        raise RuntimeError(combined)
    local: dict[str, dict[str, int]] = {}
    for line in completed.stdout.splitlines():
        match = re.fullmatch(
            r"LOCAL_(\d+)=(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+)",
            line,
        )
        if match is None:
            continue
        (
            prime,
            conductor_exponent,
            kodaira_code,
            tamagawa,
            discriminant_valuation,
            local_euler_coefficient,
        ) = map(int, match.groups())
        local[str(prime)] = {
            "conductor_exponent": conductor_exponent,
            "kodaira_code": kodaira_code,
            "kodaira_symbol": f"I_{kodaira_code - 4}",
            "tamagawa_number": tamagawa,
            "minimal_discriminant_valuation": discriminant_valuation,
            "local_euler_coefficient": local_euler_coefficient,
        }
    if set(local) != {str(prime) for prime in primes}:
        raise RuntimeError(f"PARI local replay was incomplete: {combined}")
    version = subprocess.run(
        [gp, "--version-short"],
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    return {"version": version, "setrand_seed": 1, "local_data": local}


def build_fermigier_seed(*, maximum_height: int = 2**23) -> dict[str, Any]:
    """Build the deterministic three-prime high-family seed manifest."""

    derivative = tuple(
        degree * FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS[degree]
        for degree in range(1, len(FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS))
    )
    local_constraints: list[dict[str, Any]] = []
    for item in CONSTRAINTS:
        prime = item["prime"]
        exponent = item["exponent"]
        roots_mod_prime = hensel_lift_roots(
            FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS, prime, 1
        )
        roots = hensel_lift_roots(
            FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS, prime, exponent
        )
        reductions = []
        for root in roots_mod_prime:
            coefficients = tuple(
                int(value) for value in fermigier_canonical_coefficients(root)
            )
            local = weierstrass_local_data(coefficients, prime)
            reductions.append(
                {
                    "root": root,
                    "derivative_mod_prime": evaluate_polynomial(
                        derivative, root, prime
                    ),
                    "reduction": local.reduction,
                    "local_euler_coefficient": local.local_euler_coefficient,
                }
            )
        local_constraints.append(
            {
                "prime": prime,
                "exponent": exponent,
                "modulus": prime**exponent,
                "roots_mod_prime": roots_mod_prime,
                "roots_mod_prime_power": roots,
                "root_reductions": reductions,
            }
        )

    moduli = [item["modulus"] for item in local_constraints]
    combinations: list[dict[str, Any]] = []
    for roots in product(
        *(item["roots_mod_prime_power"] for item in local_constraints)
    ):
        residue, modulus = crt(zip(roots, moduli, strict=True))
        basis = gauss_reduce_congruence_lattice(
            residue, modulus, weights=(1, 1)
        )
        representatives, checked_height = first_rationals_by_height(
            residue,
            modulus,
            count=1,
            maximum_height=maximum_height,
            weights=(1, 1),
        )
        representative = representatives[0]
        homogeneous = homogenized_discriminant_factor(
            representative.numerator, representative.denominator
        )
        combinations.append(
            {
                "roots": list(roots),
                "crt_residue": residue,
                "crt_modulus": modulus,
                "reduced_basis": [list(basis[0]), list(basis[1])],
                "numerator": representative.numerator,
                "denominator": representative.denominator,
                "height": representative.height,
                "weighted_norm": representative.weighted_norm,
                "exhaustive_height_checked": checked_height,
                "forced_valuations": {
                    str(item["prime"]): p_adic_valuation(
                        homogeneous, item["prime"]
                    )
                    for item in local_constraints
                },
            }
        )
    combinations.sort(
        key=lambda item: (
            item["height"],
            item["weighted_norm"],
            abs(item["numerator"]),
            item["denominator"],
            item["numerator"] < 0,
            item["roots"],
        )
    )
    best = combinations[0]
    numerator = best["numerator"]
    denominator = best["denominator"]
    homogeneous = homogenized_discriminant_factor(numerator, denominator)
    forced_product = 1
    for item in local_constraints:
        forced_product *= item["prime"] ** item["exponent"]
    if homogeneous % forced_product:
        raise ArithmeticError("forced prime powers do not divide the factor")
    coefficients = fermigier_canonical_coefficients(
        Fraction(numerator, denominator)
    )
    primes = tuple(item["prime"] for item in local_constraints)
    pari = _pari_local_replay(numerator, denominator, primes)
    return {
        "schema": "elliptic-curves.fermigier-crt-seed.v1",
        "claim_level": "exact_local_constraint_seed",
        "generator": "elliptic-curves/scripts/run_fermigier_crt_seed.py",
        "canonical_pinned_command": (
            "python3 elliptic-curves/scripts/run_fermigier_crt_seed.py "
            "--output artifacts/generated-results/elliptic-curves/"
            "fermigier_crt_seed_v1.json"
        ),
        "randomness": "none; all root combinations and height boxes are exhaustive",
        "family": {
            "adapter_parameter": "u=s/2",
            "weierstrass_model": "[1,a2(u),1,a4(u),a6(u)]",
            "discriminant_factor_coordinate": "canonical adapter u",
            "published_generic_rank_lower_bound": 12,
            "generic_independence_reproduced_here": False,
            "normalization_warning": (
                "the paper's printed E22 shift and the exact reconstruction "
                "differ by a still-unexplained factor two"
            ),
        },
        "search": {
            "maximum_height": maximum_height,
            "lattice_weights": [1, 1],
            "representative_count_per_residue": 1,
            "constraints": local_constraints,
            "root_combinations_tested": len(combinations),
            "combinations": combinations,
        },
        "best_seed": {
            "numerator": numerator,
            "denominator": denominator,
            "height": best["height"],
            "crt_residue": best["crt_residue"],
            "crt_modulus": best["crt_modulus"],
            "homogenized_discriminant_factor": homogeneous,
            "forced_prime_power_product": forced_product,
            "uncontrolled_cofactor": homogeneous // forced_product,
            "weierstrass_coefficients": [
                _format_fraction(value) for value in coefficients
            ],
            "pari_gp": pari,
        },
        "limitations": {
            "global_conductor": "not computed; the uncontrolled cofactor is not factored",
            "rank": "no specialization points or independence certificate computed",
            "record_status": "not a target candidate",
        },
    }


def dump_fermigier_seed(result: dict[str, Any], output: Path) -> None:
    """Write a fresh manifest and refuse to replace an existing artifact."""

    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing artifact: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
