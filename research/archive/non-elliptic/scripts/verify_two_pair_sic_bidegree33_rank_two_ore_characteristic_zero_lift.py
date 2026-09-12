#!/usr/bin/env python3
"""Verify the reconstructed characteristic-zero order-14 Ore operator.

This checker validates the simultaneous projective reconstruction against
every cached prime, checks the primitive integer normalization and the known
forward factor, and tests exact rational moments.  These are exact finite
checks of the lifted operator; they do not replace a characteristic-zero
all-order telescoping certificate.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_sic_bidegree33_rank_two_interior_cyclic_split import (  # noqa: E402
    exact_normalized_moments,
)


LIFT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json"
)
CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift_verification.json"
)
ORDER = 14
M_DEGREE = 58
FORWARD_SHIFTS = (32, 34, 35, 37, 38, 40, 41, 43)
EXACT_MAXIMUM_MOMENT = 40


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(coefficients: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def main() -> None:
    lift = json.loads(LIFT.read_text())
    cache = json.loads(CACHE.read_text())
    if lift["status"] != "stable simultaneous rational reconstruction":
        raise RuntimeError("operator lift is not stable")
    if lift["kind"] != "common":
        raise RuntimeError("operator lift is not the common order-14 factor")
    if lift["operator"] != {
        "order": ORDER,
        "m_degree": M_DEGREE,
        "coefficient_count": (ORDER + 1) * (M_DEGREE + 1),
    }:
        raise RuntimeError("unexpected lifted-operator shape")

    denominator = int(lift["common_denominator"])
    coefficients = [
        [int(value) for value in polynomial]
        for polynomial in lift["primitive_integer_coefficients"]
    ]
    if len(coefficients) != ORDER + 1:
        raise RuntimeError("wrong lifted shift order")
    if any(len(polynomial) != M_DEGREE + 1 for polynomial in coefficients):
        raise RuntimeError("wrong lifted m-degree width")
    if any(polynomial[-1] == 0 for polynomial in coefficients):
        raise RuntimeError("a lifted coefficient lost degree 58")
    content = 0
    for polynomial in coefficients:
        for coefficient in polynomial:
            content = gcd(content, abs(coefficient))
    if content != 1:
        raise RuntimeError("lifted integer operator is not primitive")
    if coefficients[-1][-1] != denominator:
        raise RuntimeError("lifted normalization is not monic at S^14*m^58")

    images = {
        int(prime): [int(value) for value in image]
        for prime, image in cache["images"].items()
    }
    requested_primes = [
        int(prime)
        for prime in lift["build_primes"] + lift["holdout_primes"]
    ]
    flat = [
        coefficient
        for polynomial in coefficients
        for coefficient in polynomial
    ]
    for prime in requested_primes:
        if prime not in images:
            raise RuntimeError(f"missing cached prime {prime}")
        if denominator % prime == 0:
            raise RuntimeError(f"normalization denominator vanishes at {prime}")
        inverse = pow(denominator, -1, prime)
        reduction = [
            coefficient * inverse % prime for coefficient in flat
        ]
        if reduction != images[prime]:
            raise RuntimeError(f"lifted operator mismatch at prime {prime}")

    m = sp.symbols("m")
    forward = sp.Poly(
        sum(
            sp.Integer(coefficient) * m**degree
            for degree, coefficient in enumerate(coefficients[-1])
        ),
        m,
        domain=sp.ZZ,
    )
    fixed = sp.Poly(
        sp.prod(3 * m + shift for shift in FORWARD_SHIFTS),
        m,
        domain=sp.ZZ,
    )
    quotient, remainder = sp.div(forward, fixed, domain=sp.ZZ)
    if not remainder.is_zero or quotient.degree() != 50:
        raise RuntimeError("exact forward-factor audit failed")

    moments = exact_normalized_moments(EXACT_MAXIMUM_MOMENT)
    exact_rows = EXACT_MAXIMUM_MOMENT - ORDER + 1
    for moment_index in range(exact_rows):
        residual = sum(
            Fraction(evaluate(polynomial, moment_index))
            * moments[moment_index + shift]
            for shift, polynomial in enumerate(coefficients)
        )
        if residual:
            raise RuntimeError(
                f"exact rational recurrence mismatch at m={moment_index}"
            )

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "ore-characteristic-zero-lift-verification-v1"
        ),
        "status": (
            "exact finite verification of a stable characteristic-zero "
            "order-14 operator lift; all-order telescoping remains open"
        ),
        "point": 0,
        "operator": {"order": ORDER, "m_degree": M_DEGREE},
        "primitive_integer_operator": True,
        "maximum_primitive_coefficient_bits": lift[
            "successful_lattice"
        ]["maximum_coefficient_bits"],
        "reconstruction_prime_count": len(lift["build_primes"]),
        "independent_holdout_prime_count": len(lift["holdout_primes"]),
        "all_cached_images_replayed": len(requested_primes),
        "forward_fixed_factor": (
            "product_(k in {32,34,35,37,38,40,41,43})(3m+k)"
        ),
        "forward_residual_degree": quotient.degree(),
        "exact_rational_moment_rows": exact_rows,
        "all_order_status": (
            "not yet proved over characteristic zero; the first exact "
            "interior divergence level closes, while a five-level direct "
            "expanded run exceeded 1800 seconds"
        ),
        "files_sha256": {
            str(LIFT.relative_to(ROOT)): sha256(LIFT),
            str(CACHE.relative_to(ROOT)): sha256(CACHE),
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS primitive order-{ORDER}, m-degree-{M_DEGREE} integer lift")
    print(f"PASS replayed {len(requested_primes)} exact modular images")
    print("PASS five fresh holdout primes")
    print("PASS exact forward factor and degree-50 residual")
    print(f"PASS {exact_rows} exact rational moment identities")
    print("PASS all-order characteristic-zero certificate remains open")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
