#!/usr/bin/env python3
"""Reconstruct the sampled fixed-fiber order-14 Ore operator over QQ.

This is a bounded research calculation, not a verifier.  It recomputes the
primitive order-14 common right factor from the order-18 and order-27 modular
operators at the first integral rank-two point, combines coefficient images
by CRT, and tests balanced rational reconstruction at an independent prime.

The purpose is to produce a prescribed characteristic-zero target for a
direct divergence-certificate calculation without computing the full
D-module pushforward ideal.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd, isqrt
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from sympy import isprime, nextprime

from verify_two_pair_sic_bidegree33_rank_two_ore_gcd import (
    EXPECTED_COMMON_DEGREE,
    EXPECTED_COMMON_ORDER,
    EXPECTED_REMAINDER_ORDERS,
    MAXIMUM_MOMENT,
    ORDER_DEGREES,
    ROOT,
    SOURCE,
    ShiftOreField,
    run_moments,
    run_operator,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_research.json"
)
IMAGE_CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json"
)
DEFAULT_PRIMES = (
    1000003,
    1000033,
    1000037,
    1000039,
)
COEFFICIENT_COUNT = (EXPECTED_COMMON_ORDER + 1) * (
    EXPECTED_COMMON_DEGREE + 1
)


def primitive_operator(
    executable: Path,
    prime: int,
    point: int,
) -> list[object]:
    ore = ShiftOreField(prime)
    order_18 = ore.parse_operator(
        run_operator(executable, prime, point, *ORDER_DEGREES[0]),
        *ORDER_DEGREES[0],
    )
    order_27 = ore.parse_operator(
        run_operator(executable, prime, point, *ORDER_DEGREES[1]),
        *ORDER_DEGREES[1],
    )
    common, remainder_orders = ore.greatest_common_right_divisor(
        order_27,
        order_18,
    )
    assert tuple(remainder_orders) == EXPECTED_REMAINDER_ORDERS
    primitive = ore.primitive_polynomial_operator(common)
    assert len(primitive) == EXPECTED_COMMON_ORDER + 1
    assert all(
        coefficient.degree() == EXPECTED_COMMON_DEGREE
        for coefficient in primitive
    )
    sequence = run_moments(executable, prime, point)
    assert len(sequence) == MAXIMUM_MOMENT + 1
    assert ore.verify_recurrence(primitive, sequence) == 487
    return primitive


def flatten_operator(operator: list[object], prime: int) -> list[int]:
    flattened: list[int] = []
    for polynomial in operator:
        coefficients = [0] * (EXPECTED_COMMON_DEGREE + 1)
        for monomial, coefficient in polynomial.to_dict().items():
            coefficients[monomial[0]] = int(coefficient) % prime
        flattened.extend(coefficients)
    assert len(flattened) == COEFFICIENT_COUNT
    return flattened


def crt_merge(
    residues: list[int],
    modulus: int,
    image: list[int],
    prime: int,
) -> tuple[list[int], int]:
    if not residues:
        return image[:], prime
    assert len(residues) == len(image)
    inverse = pow(modulus, -1, prime)
    merged = [
        residue
        + modulus * (((value - residue) % prime) * inverse % prime)
        for residue, value in zip(residues, image, strict=True)
    ]
    return merged, modulus * prime


def rational_reconstruct(
    residue: int,
    modulus: int,
) -> Fraction | None:
    residue %= modulus
    bound = isqrt(modulus // 2)
    old_remainder, remainder = modulus, residue
    old_coefficient, coefficient = 0, 1
    while abs(remainder) >= bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_coefficient, coefficient = (
            coefficient,
            old_coefficient - quotient * coefficient,
        )
    if (
        coefficient == 0
        or abs(coefficient) >= bound
        or gcd(remainder, coefficient) != 1
    ):
        return None
    if coefficient < 0:
        remainder = -remainder
        coefficient = -coefficient
    if (residue * coefficient - remainder) % modulus:
        return None
    return Fraction(remainder, coefficient)


def matches_mod_prime(
    candidate: Fraction | None,
    residue: int,
    prime: int,
) -> bool:
    if candidate is None or candidate.denominator % prime == 0:
        return False
    value = (
        candidate.numerator
        * pow(candidate.denominator, -1, prime)
    ) % prime
    return value == residue


def nested_operator(
    candidates: list[Fraction | None],
) -> list[list[list[int] | None]]:
    result: list[list[list[int] | None]] = []
    width = EXPECTED_COMMON_DEGREE + 1
    for shift in range(EXPECTED_COMMON_ORDER + 1):
        block = candidates[shift * width : (shift + 1) * width]
        result.append(
            [
                (
                    [value.numerator, value.denominator]
                    if value is not None
                    else None
                )
                for value in block
            ]
        )
    return result


def nested_residues(values: list[int]) -> list[list[int]]:
    result: list[list[int]] = []
    width = EXPECTED_COMMON_DEGREE + 1
    for shift in range(EXPECTED_COMMON_ORDER + 1):
        result.append(values[shift * width : (shift + 1) * width])
    return result


def load_image_cache(path: Path, point: int) -> dict[int, list[int]]:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    if payload.get("format") != (
        "two-pair-sic-bidegree33-rank-two-"
        "ore-reconstruct-images-v1"
    ):
        raise ValueError("unsupported Ore-reconstruction image cache")
    if int(payload.get("point", -1)) != point:
        raise ValueError("Ore-reconstruction image-cache point mismatch")
    images = {
        int(prime): [int(value) for value in values]
        for prime, values in payload["images"].items()
    }
    if not all(len(values) == COEFFICIENT_COUNT for values in images.values()):
        raise ValueError("invalid cached Ore-operator image width")
    return images


def write_image_cache(
    path: Path,
    point: int,
    images: dict[int, list[int]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "ore-reconstruct-images-v1"
        ),
        "status": (
            "exact modular primitive order-14 operator images for "
            "resumable CRT reconstruction"
        ),
        "point": point,
        "coefficient_count": COEFFICIENT_COUNT,
        "images": {
            str(prime): values
            for prime, values in sorted(images.items())
        },
    }
    path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--point", type=int, default=0)
    parser.add_argument(
        "--primes",
        type=int,
        nargs="+",
        default=DEFAULT_PRIMES,
    )
    parser.add_argument(
        "--prime-count",
        type=int,
        help="use this many consecutive primes strictly above --prime-start",
    )
    parser.add_argument("--prime-start", type=int, default=1_000_000)
    parser.add_argument(
        "--holdout-count",
        type=int,
        default=1,
        help="reserve this many final primes as independent holdouts",
    )
    parser.add_argument(
        "--image-cache",
        type=Path,
        default=IMAGE_CACHE,
        help="resumable exact modular-image cache",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.prime_count is not None:
        if arguments.prime_count <= arguments.holdout_count:
            raise ValueError("prime-count must exceed holdout-count")
        arguments.primes = []
        prime = arguments.prime_start
        for _ in range(arguments.prime_count):
            prime = int(nextprime(prime))
            arguments.primes.append(prime)

    if not 1 <= arguments.holdout_count < len(arguments.primes):
        raise ValueError(
            "holdout-count must be positive and smaller than prime count"
        )
    if len(arguments.primes) < 2:
        raise ValueError("at least one reconstruction and one holdout prime")
    if len(set(arguments.primes)) != len(arguments.primes):
        raise ValueError("primes must be distinct")
    if not all(isprime(prime) for prime in arguments.primes):
        raise ValueError("every modulus must be prime")
    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required")

    cached_images = load_image_cache(
        arguments.image_cache,
        arguments.point,
    )
    images: list[tuple[int, list[int]]] = []
    with tempfile.TemporaryDirectory(
        prefix="sic33-ore-reconstruct-"
    ) as path:
        executable = Path(path) / "recurrence-probe"
        subprocess.run(
            [
                compiler,
                "-O3",
                "-std=c++17",
                str(SOURCE),
                "-o",
                str(executable),
            ],
            check=True,
            timeout=30,
        )
        for prime in arguments.primes:
            image = cached_images.get(prime)
            if image is None:
                operator = primitive_operator(
                    executable,
                    prime,
                    arguments.point,
                )
                image = flatten_operator(operator, prime)
                cached_images[prime] = image
                write_image_cache(
                    arguments.image_cache,
                    arguments.point,
                    cached_images,
                )
                print(f"PASS computed prime {prime}", flush=True)
            else:
                print(f"PASS cached prime {prime}", flush=True)
            images.append((prime, image))

    build_images = images[: -arguments.holdout_count]
    holdout_images = images[-arguments.holdout_count :]
    residues: list[int] = []
    modulus = 1
    for prime, image in build_images:
        residues, modulus = crt_merge(
            residues,
            modulus,
            image,
            prime,
        )
    candidates = [
        rational_reconstruct(residue, modulus)
        for residue in residues
    ]
    reconstructed = sum(value is not None for value in candidates)
    per_holdout_matches = {
        str(prime): sum(
            matches_mod_prime(candidate, residue, prime)
            for candidate, residue in zip(
                candidates,
                holdout,
                strict=True,
            )
        )
        for prime, holdout in holdout_images
    }
    stable = (
        reconstructed == COEFFICIENT_COUNT
        and all(
            matches == COEFFICIENT_COUNT
            for matches in per_holdout_matches.values()
        )
    )

    payload = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "ore-reconstruct-research-v1"
        ),
        "status": (
            "bounded modular reconstruction target; not a telescoping "
            "certificate"
        ),
        "point": arguments.point,
        "operator": {
            "order": EXPECTED_COMMON_ORDER,
            "m_degree": EXPECTED_COMMON_DEGREE,
            "coefficient_order": (
                "shift-major; each polynomial is low-to-high in m"
            ),
            "coefficient_count": COEFFICIENT_COUNT,
        },
        "reconstruction_primes": [
            prime for prime, _ in build_images
        ],
        "crt_modulus": str(modulus),
        "rational_reconstruction_bound": str(isqrt(modulus // 2)),
        "holdout_primes": [
            prime for prime, _ in holdout_images
        ],
        "modular_operator": {
            "prime": images[0][0],
            "coefficients": nested_residues(images[0][1]),
        },
        "balanced_reconstructions": reconstructed,
        "holdout_matches": per_holdout_matches,
        "fully_stable": stable,
        "candidate_coefficients": nested_operator(candidates),
        "remaining_gate": (
            "add primes until every coefficient is stable at independent "
            "holdouts, then solve and independently verify a polynomial "
            "relative-divergence certificate for the prescribed operator"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"PASS reconstructed {reconstructed}/{COEFFICIENT_COUNT}; "
        f"holdout matches {per_holdout_matches}"
    )
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
