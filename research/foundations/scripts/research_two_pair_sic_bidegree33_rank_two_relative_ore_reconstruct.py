#!/usr/bin/env python3
"""Reconstruct the fixed-fiber order-18 relative Ore operator over QQ.

The order-18 recurrence is the natural zero-boundary target for the
18-dimensional relative logarithmic critical algebra.  This resumable
multi-prime calculation reconstructs its normalized coefficients and tests
them at independent holdout primes.  It remains a bounded reconstruction
until an all-order relative-divergence certificate is supplied.
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path
import shutil
import subprocess
import tempfile

from sympy import isprime, nextprime

from research_two_pair_sic_bidegree33_rank_two_ore_reconstruct import (
    crt_merge,
    matches_mod_prime,
    rational_reconstruct,
)
from verify_two_pair_sic_bidegree33_rank_two_ore_gcd import (
    MAXIMUM_MOMENT,
    ROOT,
    SOURCE,
    ShiftOreField,
    run_moments,
    run_operator,
)


ORDER = 18
M_DEGREE = 18
COEFFICIENT_COUNT = (ORDER + 1) * (M_DEGREE + 1)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_relative_ore_reconstruct.json"
)
IMAGE_CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_relative_ore_images.json"
)
DEFAULT_PRIMES = (
    1000003,
    1000033,
    1000037,
    1000039,
    1000081,
    1000099,
    1000117,
    1000121,
)


def flatten_operator(operator: list[object], prime: int) -> list[int]:
    flattened: list[int] = []
    for polynomial in operator:
        coefficients = [0] * (M_DEGREE + 1)
        for monomial, coefficient in polynomial.to_dict().items():
            coefficients[monomial[0]] = int(coefficient) % prime
        flattened.extend(coefficients)
    if len(flattened) != COEFFICIENT_COUNT:
        raise RuntimeError("relative operator coefficient-count mismatch")
    return flattened


def relative_operator(
    executable: Path,
    prime: int,
    point: int,
) -> list[int]:
    ore = ShiftOreField(prime)
    parsed = ore.parse_operator(
        run_operator(
            executable,
            prime,
            point,
            ORDER,
            M_DEGREE,
        ),
        ORDER,
        M_DEGREE,
    )
    operator = []
    for coefficient in parsed:
        if coefficient.denom.degree() != 0:
            raise RuntimeError("fitted relative operator is not polynomial")
        denominator = int(coefficient.denom.LC) % prime
        inverse = ore.polynomial_ring.domain.convert(
            pow(denominator, -1, prime)
        )
        operator.append(coefficient.numer * inverse)
    sequence = run_moments(executable, prime, point)
    if ore.verify_recurrence(operator, sequence) != (
        MAXIMUM_MOMENT - ORDER + 1
    ):
        raise RuntimeError("relative recurrence holdout failure")
    return flatten_operator(operator, prime)


def load_cache(path: Path, point: int) -> dict[int, list[int]]:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    if payload.get("format") != (
        "two-pair-sic-bidegree33-rank-two-relative-ore-images-v1"
    ):
        raise ValueError("unsupported relative-Ore image cache")
    if int(payload.get("point", -1)) != point:
        raise ValueError("relative-Ore image-cache point mismatch")
    images = {
        int(prime): [int(value) for value in values]
        for prime, values in payload["images"].items()
    }
    if not all(len(image) == COEFFICIENT_COUNT for image in images.values()):
        raise ValueError("invalid relative-Ore cached image width")
    return images


def write_cache(
    path: Path,
    point: int,
    images: dict[int, list[int]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-relative-ore-images-v1"
        ),
        "status": (
            "exact modular normalized order-18 relative-operator images"
        ),
        "point": point,
        "coefficient_count": COEFFICIENT_COUNT,
        "images": {
            str(prime): image
            for prime, image in sorted(images.items())
        },
    }
    path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


def nested_residues(values: list[int]) -> list[list[int]]:
    width = M_DEGREE + 1
    return [
        values[shift * width : (shift + 1) * width]
        for shift in range(ORDER + 1)
    ]


def nested_candidates(
    candidates: list[object | None],
) -> list[list[list[int] | None]]:
    width = M_DEGREE + 1
    return [
        [
            (
                [value.numerator, value.denominator]
                if value is not None
                else None
            )
            for value in candidates[shift * width : (shift + 1) * width]
        ]
        for shift in range(ORDER + 1)
    ]


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
    parser.add_argument("--holdout-count", type=int, default=2)
    parser.add_argument("--image-cache", type=Path, default=IMAGE_CACHE)
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
    if len(set(arguments.primes)) != len(arguments.primes):
        raise ValueError("primes must be distinct")
    if not all(isprime(prime) for prime in arguments.primes):
        raise ValueError("every modulus must be prime")
    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required")

    cache = load_cache(arguments.image_cache, arguments.point)
    images: list[tuple[int, list[int]]] = []
    with tempfile.TemporaryDirectory(
        prefix="sic33-relative-ore-reconstruct-"
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
            image = cache.get(prime)
            if image is None:
                image = relative_operator(
                    executable,
                    prime,
                    arguments.point,
                )
                cache[prime] = image
                write_cache(arguments.image_cache, arguments.point, cache)
                print(f"PASS computed prime {prime}", flush=True)
            else:
                print(f"PASS cached prime {prime}", flush=True)
            images.append((prime, image))

    build_images = images[: -arguments.holdout_count]
    holdout_images = images[-arguments.holdout_count :]
    residues: list[int] = []
    modulus = 1
    for prime, image in build_images:
        residues, modulus = crt_merge(residues, modulus, image, prime)
    candidates = [
        rational_reconstruct(residue, modulus)
        for residue in residues
    ]
    reconstructed = sum(candidate is not None for candidate in candidates)
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
            "relative-ore-reconstruct-v1"
        ),
        "status": (
            "stable rational order-18 relative-operator reconstruction"
            if stable
            else "bounded modular relative-operator reconstruction"
        ),
        "point": arguments.point,
        "operator": {
            "order": ORDER,
            "m_degree": M_DEGREE,
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
        "candidate_coefficients": nested_candidates(candidates),
        "remaining_gate": (
            "construct and replay an exact characteristic-zero "
            "relative-divergence certificate"
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
