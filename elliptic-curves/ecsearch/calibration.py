"""End-to-end calibration of prime-power discriminant engineering.

The calibration family is

    E_t: y^2 = x^3 - t^2*x + t^2,

with independent sections ``Q=(t,t)`` and ``R=(-t,t)`` (certified by the
rank-two specialization at ``t=5``).  The additional visible point ``(0,t)``
is ``-(Q+R)``.  For ``t=a/b`` the integral model used here is

    Y^2 = X^3 - a^2*b^2*X + a^2*b^4,

and its discriminant is

    -16*a^4*b^6*(27*b^2 - 4*a^2).
"""

from __future__ import annotations

import ast
import json
import re
import shutil
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

from .crt_lattice import (
    crt,
    first_rationals_by_height,
    gauss_reduce_congruence_lattice,
    hensel_lift_roots,
    p_adic_valuation,
)


DISCRIMINANT_FACTOR_COEFFICIENTS = (27, 0, -4)
CONSTRAINTS = (
    {"prime": 23, "exponent": 3},
    {"prime": 47, "exponent": 2},
    {"prime": 73, "exponent": 2},
)


def binary_discriminant_factor(numerator: int, denominator: int) -> int:
    return 27 * denominator * denominator - 4 * numerator * numerator


def integral_model(numerator: int, denominator: int) -> tuple[int, int]:
    a4 = -(numerator * numerator) * (denominator * denominator)
    a6 = (numerator * numerator) * (denominator**4)
    return a4, a6


def model_discriminant(numerator: int, denominator: int) -> int:
    return (
        -16
        * numerator**4
        * denominator**6
        * binary_discriminant_factor(numerator, denominator)
    )


def _legendre_symbol(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    symbol = pow(value, (prime - 1) // 2, prime)
    return -1 if symbol == prime - 1 else symbol


def _parse_factor_matrix(value: str) -> list[list[int]]:
    pairs = re.findall(r"(-?\d+)\s*,\s*(-?\d+)", value)
    return [[int(prime), int(exponent)] for prime, exponent in pairs]


def _run_pari(numerator: int, denominator: int) -> dict[str, Any]:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP executable 'gp' is required")
    a4, a6 = integral_model(numerator, denominator)
    point_y = numerator * denominator * denominator
    program = f"""
setrand(1);
E=ellinit([0,0,0,{a4},{a6}]);
P=[0,{point_y}];
if(!ellisoncurve(E,P),error("calibration point is not on the curve"));
G=ellglobalred(E);
T=elltors(E);
R=ellrank(E,1);
print("CONDUCTOR=",G[1]);
print("GLOBAL_CHANGE=",G[2]);
print("TAMAGAWA_PRODUCT=",G[3]);
print("CONDUCTOR_FACTORS=",G[4]);
print("TORSION_ORDER=",T[1]);
print("POINT_ORDER=",ellorder(E,P));
print("RANK_LOWER=",R[1]);
print("RANK_UPPER=",R[2]);
print("DISCRIMINANT=",E.disc);
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "***" in combined:
        raise RuntimeError(combined)
    values: dict[str, str] = {}
    for line in completed.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    required = {
        "CONDUCTOR",
        "GLOBAL_CHANGE",
        "TAMAGAWA_PRODUCT",
        "CONDUCTOR_FACTORS",
        "TORSION_ORDER",
        "POINT_ORDER",
        "RANK_LOWER",
        "RANK_UPPER",
        "DISCRIMINANT",
    }
    missing = required - values.keys()
    if missing:
        raise RuntimeError(f"PARI output omitted {sorted(missing)}: {combined}")
    version = subprocess.run(
        [gp, "--version-short"],
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    return {
        "version": version,
        "conductor": int(values["CONDUCTOR"]),
        "global_minimal_change": ast.literal_eval(values["GLOBAL_CHANGE"]),
        "tamagawa_product": int(values["TAMAGAWA_PRODUCT"]),
        "conductor_factorization": _parse_factor_matrix(
            values["CONDUCTOR_FACTORS"]
        ),
        "torsion_order": int(values["TORSION_ORDER"]),
        "known_point_order": int(values["POINT_ORDER"]),
        "rank_bounds": [int(values["RANK_LOWER"]), int(values["RANK_UPPER"])],
        "discriminant": int(values["DISCRIMINANT"]),
    }


def build_calibration(*, maximum_height: int = 262144) -> dict[str, Any]:
    """Run the complete fixed calibration and return a deterministic manifest."""

    local_constraints: list[dict[str, Any]] = []
    for constraint in CONSTRAINTS:
        prime = constraint["prime"]
        exponent = constraint["exponent"]
        modulus = prime**exponent
        roots = hensel_lift_roots(
            DISCRIMINANT_FACTOR_COEFFICIENTS, prime, exponent
        )
        local_constraints.append(
            {
                "prime": prime,
                "exponent": exponent,
                "modulus": modulus,
                "roots": roots,
                "reduction_goal": "split_multiplicative",
                "split_test_legendre_6": _legendre_symbol(6, prime),
            }
        )
    combinations: list[dict[str, Any]] = []
    moduli = [constraint["modulus"] for constraint in local_constraints]
    root_lists = [constraint["roots"] for constraint in local_constraints]
    for roots in product(*root_lists):
        residue, modulus = crt(zip(roots, moduli, strict=True))
        reduced_basis = gauss_reduce_congruence_lattice(
            residue, modulus, weights=(1, 1)
        )
        representatives, checked_height = first_rationals_by_height(
            residue,
            modulus,
            count=1,
            maximum_height=maximum_height,
            weights=(1, 1),
        )
        if not representatives:
            raise RuntimeError("bounded lattice neighbourhood has no admissible pair")
        representative = representatives[0]
        factor = binary_discriminant_factor(
            representative.numerator, representative.denominator
        )
        combinations.append(
            {
                "roots": list(roots),
                "crt_residue": residue,
                "crt_modulus": modulus,
                "reduced_basis": [list(reduced_basis[0]), list(reduced_basis[1])],
                "numerator": representative.numerator,
                "denominator": representative.denominator,
                "height": representative.height,
                "weighted_norm": representative.weighted_norm,
                "exhaustive_height_checked": checked_height,
                "binary_discriminant_factor": factor,
                "forced_valuations": {
                    str(constraint["prime"]): p_adic_valuation(
                        factor, constraint["prime"]
                    )
                    for constraint in local_constraints
                },
            }
        )
    combinations.sort(
        key=lambda item: (
            item["height"],
            item["weighted_norm"],
            abs(item["numerator"]),
            item["denominator"],
            item["numerator"],
            item["roots"],
        )
    )
    best = combinations[0]
    numerator = best["numerator"]
    denominator = best["denominator"]
    a4, a6 = integral_model(numerator, denominator)
    pari = _run_pari(numerator, denominator)
    discriminant = model_discriminant(numerator, denominator)
    if pari["discriminant"] != discriminant:
        raise RuntimeError("PARI and direct discriminant calculations disagree")
    return {
        "schema": "elliptic-curves.crt-lattice-calibration.v1",
        "claim_level": "exact_computation",
        "generator": "elliptic-curves/scripts/run_crt_lattice_calibration.py",
        "canonical_pinned_command": (
            "python3 elliptic-curves/scripts/run_crt_lattice_calibration.py "
            "--output artifacts/generated-results/elliptic-curves/"
            "crt_lattice_calibration_v1.json"
        ),
        "randomness": {
            "search": "deterministic exhaustive enumeration",
            "pari_gp_setrand_seed": 1,
        },
        "family": {
            "generic_model": "y^2 = x^3 - t^2*x + t^2",
            "generic_rank_lower_bound": 2,
            "known_independent_sections": [["t", "t"], ["-t", "t"]],
            "dependent_section": {
                "point": [0, "t"],
                "relation": "P=-(Q+R)",
            },
            "integral_specialization": (
                "Y^2 = X^3 - a^2*b^2*X + a^2*b^4 for t=a/b"
            ),
            "discriminant": "-16*a^4*b^6*(27*b^2 - 4*a^2)",
        },
        "search": {
            "polynomial_coefficients_low_to_high": list(
                DISCRIMINANT_FACTOR_COEFFICIENTS
            ),
            "maximum_height": maximum_height,
            "lattice_weights": [1, 1],
            "representative_count_per_residue": 1,
            "constraints": local_constraints,
            "root_combinations_tested": len(combinations),
            "combinations": combinations,
        },
        "best_candidate": {
            "numerator": numerator,
            "denominator": denominator,
            "crt_residue": best["crt_residue"],
            "crt_modulus": best["crt_modulus"],
            "height": best["height"],
            "binary_discriminant_factor": best["binary_discriminant_factor"],
            "factor_identity": "27*b^2 - 4*a^2 = 23^3*47^2*73^2",
            "weierstrass_coefficients": [0, 0, 0, a4, a6],
            "known_integral_point": [0, numerator * denominator * denominator],
            "discriminant": discriminant,
            "pari_gp": pari,
        },
        "interpretation": {
            "result": (
                "The selected multiplicative primes occur to exponents 3,2,2 "
                "in the minimal discriminant factor and only to exponent one "
                "in the conductor."
            ),
            "limitation": (
                "This calibration family has generic rank at least two but is "
                "not a candidate for either record target."
            ),
        },
    }


def dump_calibration(result: dict[str, Any], output: Path) -> None:
    """Write a fresh compact manifest; never overwrite a pinned artifact."""

    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing artifact: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
