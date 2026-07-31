#!/usr/bin/env python3
"""Research the compact relative Picard--Fuchs relation at point zero.

The beta substitution

    x = u,  y = u*t/(1-t)

turns the generating two-form into

    (x+y) dx dy / ((x+y)^3-z*Phi(x,y)).

The closed-cycle Griffiths--Dwork calculation has differential order eight.
For the interval period, this script independently fits the corresponding
inhomogeneous differential relation modulo primes, converts its homogeneous
tail to a shift operator, and compares that operator with the stored order-14
Ore factor.  Everything in this file is a finite modular calculation; the
fitted relation is not an all-order characteristic-zero certificate.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

from flint import nmod_mat
from sympy import Poly, cancel, expand, nextprime, symbols

from verify_two_pair_sic_bidegree33_rank_two_ore_gcd import (
    ROOT,
    SOURCE,
    ShiftOreField,
)


ORDER = 8
Z_DEGREE = 72
COEFFICIENT_COUNT = (ORDER + 1) * (Z_DEGREE + 1)
MAXIMUM_MOMENT = 850
FIT_EXTRA_ROWS = 40
REFERENCE_PRIME = 1_000_003
TAIL_SHIFT = 64
EXPECTED_RESIDUAL_DEGREE = 55
EXPECTED_SHIFT_ORDER = 64
EXPECTED_SHIFT_M_DEGREE = 8
EXPECTED_RIGHT_FACTOR_ORDER = 14
EXPECTED_LEFT_QUOTIENT_ORDER = 50

IMAGE_CACHE = (
    ROOT
    / "artifacts"
    / "local"
    / "two_pair_sic_bidegree33_rank_two_compact_relative_pf_images.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json"
)
CLOSED_OPERATOR = (
    ROOT
    / "artifacts"
    / "local"
    / "two_pair_sic_bidegree33_rank_two_compact_picard_fuchs.ore"
)
COMMON_IMAGE_CACHE = (
    ROOT
    / "artifacts"
    / "local"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json"
)


def falling_mod(value: int, order: int, prime: int) -> int:
    answer = 1
    for offset in range(order):
        answer = answer * (value - offset) % prime
    return answer


def parse_moments(output: str, maximum: int) -> list[int]:
    lines = output.strip().splitlines()
    if lines[0] != f"MOMENTS maximum={maximum} point=0":
        raise ValueError("unexpected moment-probe header")
    sequence = []
    for expected, line in enumerate(lines[1:]):
        index, value = (int(entry) for entry in line.split())
        if index != expected:
            raise ValueError("nonconsecutive moment output")
        sequence.append(value)
    if len(sequence) != maximum + 1:
        raise ValueError("incomplete moment output")
    return sequence


def compute_moments(
    executable: Path,
    prime: int,
    maximum: int,
) -> list[int]:
    completed = subprocess.run(
        [
            str(executable),
            "--moments",
            str(prime),
            str(maximum),
            "0",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=180,
    )
    return parse_moments(completed.stdout, maximum)


def differential_matrix(
    sequence: list[int],
    prime: int,
    order: int,
    degree: int,
    extra_rows: int = FIT_EXTRA_ROWS,
) -> nmod_mat:
    columns = (order + 1) * (degree + 1)
    rows = []
    for coefficient_index in range(
        degree,
        degree + columns + extra_rows,
    ):
        row = []
        for derivative_order in range(order + 1):
            for z_exponent in range(degree + 1):
                moment_index = (
                    coefficient_index - z_exponent + derivative_order
                )
                row.append(
                    falling_mod(
                        moment_index,
                        derivative_order,
                        prime,
                    )
                    * sequence[moment_index]
                    % prime
                )
        rows.append(row)
    return nmod_mat(rows, prime)


def fit_differential_operator(
    sequence: list[int],
    prime: int,
) -> tuple[list[int], int, int]:
    matrix = differential_matrix(
        sequence,
        prime,
        ORDER,
        Z_DEGREE,
    )
    basis, nullity = matrix.nullspace()
    if nullity != 1:
        raise ValueError(f"expected nullity one, got {nullity}")
    vector = [int(basis[row, 0]) % prime for row in range(COEFFICIENT_COUNT)]
    normalization = vector[-1]
    if normalization == 0:
        raise ValueError("top differential coefficient cannot normalize")
    inverse = pow(normalization, -1, prime)
    vector = [value * inverse % prime for value in vector]
    if vector[-1] != 1:
        raise AssertionError("failed to normalize differential operator")
    return vector, matrix.rank(), nullity


def nested_differential(flattened: list[int]) -> list[list[int]]:
    width = Z_DEGREE + 1
    return [
        flattened[index * width : (index + 1) * width]
        for index in range(ORDER + 1)
    ]


def differential_support(flattened: list[int]) -> list[dict[str, int]]:
    result = []
    for derivative_order, coefficients in enumerate(
        nested_differential(flattened)
    ):
        support = [
            exponent
            for exponent, coefficient in enumerate(coefficients)
            if coefficient
        ]
        result.append(
            {
                "derivative_order": derivative_order,
                "minimum_z_degree": min(support),
                "maximum_z_degree": max(support),
                "term_count": len(support),
            }
        )
    return result


def differential_residuals(
    flattened: list[int],
    sequence: list[int],
    prime: int,
) -> list[int]:
    coefficients = nested_differential(flattened)
    residuals = []
    for coefficient_index in range(len(sequence) - ORDER):
        total = 0
        for derivative_order, row in enumerate(coefficients):
            for z_exponent, value in enumerate(row):
                if value == 0 or z_exponent > coefficient_index:
                    continue
                moment_index = (
                    coefficient_index - z_exponent + derivative_order
                )
                total += (
                    value
                    * falling_mod(moment_index, derivative_order, prime)
                    * sequence[moment_index]
                )
        residuals.append(total % prime)
    return residuals


def multiply_linear_mod(
    polynomial: list[int],
    constant: int,
    prime: int,
) -> list[int]:
    result = [0] * (len(polynomial) + 1)
    for exponent, coefficient in enumerate(polynomial):
        result[exponent] = (
            result[exponent] + constant * coefficient
        ) % prime
        result[exponent + 1] = (
            result[exponent + 1] + coefficient
        ) % prime
    return result


def differential_to_shift(
    flattened: list[int],
    prime: int,
) -> tuple[list[list[int]], int]:
    coefficients = nested_differential(flattened)
    delta = max(
        z_exponent - derivative_order
        for derivative_order, row in enumerate(coefficients)
        for z_exponent, value in enumerate(row)
        if value
    )
    shift_coefficients = [
        [0] * (ORDER + 1) for _ in range(delta + 1)
    ]
    for derivative_order, row in enumerate(coefficients):
        for z_exponent, value in enumerate(row):
            if value == 0:
                continue
            shift = delta - z_exponent + derivative_order
            falling = [1]
            for offset in range(derivative_order):
                falling = multiply_linear_mod(
                    falling,
                    (shift - offset) % prime,
                    prime,
                )
            for exponent, factor in enumerate(falling):
                shift_coefficients[shift][exponent] = (
                    shift_coefficients[shift][exponent]
                    + value * factor
                ) % prime
    return shift_coefficients, delta


def polynomial_operator(
    ore: ShiftOreField,
    coefficient_rows: list[list[int]],
) -> list[Any]:
    domain = ore.polynomial_ring.domain
    return [
        ore.polynomial_ring.from_dict(
            {
                (exponent,): domain.convert(value)
                for exponent, value in enumerate(row)
                if value % ore.prime
            }
        )
        for row in coefficient_rows
    ]


def load_common_operator_image(prime: int) -> list[int]:
    payload = json.loads(COMMON_IMAGE_CACHE.read_text())
    image = payload["images"].get(str(prime))
    if image is None:
        raise ValueError(f"stored order-14 image is missing prime {prime}")
    return [int(value) for value in image]


def reference_analysis(
    flattened: list[int],
    sequence: list[int],
    prime: int,
) -> dict[str, object]:
    residuals = differential_residuals(flattened, sequence, prime)
    nonzero_residuals = [
        index for index, value in enumerate(residuals) if value
    ]
    residual_degree = max(nonzero_residuals)
    if residual_degree != EXPECTED_RESIDUAL_DEGREE:
        raise AssertionError("unexpected inhomogeneous residual degree")
    if any(residuals[residual_degree + 1 :]):
        raise AssertionError("differential tail did not vanish")

    shift_rows, delta = differential_to_shift(flattened, prime)
    if delta != TAIL_SHIFT or len(shift_rows) - 1 != EXPECTED_SHIFT_ORDER:
        raise AssertionError("unexpected shift conversion")
    shift_degrees = [
        max(exponent for exponent, value in enumerate(row) if value)
        for row in shift_rows
    ]
    if shift_degrees != [EXPECTED_SHIFT_M_DEGREE] * len(shift_rows):
        raise AssertionError("unexpected shift coefficient degrees")

    ore = ShiftOreField(prime)
    shift_polynomial = polynomial_operator(ore, shift_rows)
    recurrence_checks = 0
    for moment_index in range(len(sequence) - len(shift_rows) + 1):
        total = sum(
            ore.evaluate_polynomial(coefficient, moment_index)
            * sequence[moment_index + shift]
            for shift, coefficient in enumerate(shift_polynomial)
        )
        if total % prime:
            raise AssertionError(
                f"shift recurrence failed at m={moment_index}"
            )
        recurrence_checks += 1

    common_flat = load_common_operator_image(prime)
    common_rows = [
        common_flat[index * 59 : (index + 1) * 59]
        for index in range(15)
    ]
    common_polynomial = polynomial_operator(ore, common_rows)
    shift_field = [ore.field(value) for value in shift_polynomial]
    common_field = [ore.field(value) for value in common_polynomial]
    quotient, remainder = ore.left_division(shift_field, common_field)
    if remainder:
        raise AssertionError("compact shift relation has nonzero remainder")
    if len(quotient) - 1 != EXPECTED_LEFT_QUOTIENT_ORDER:
        raise AssertionError("unexpected compact left quotient order")
    common_gcrd, euclidean_orders = ore.greatest_common_right_divisor(
        shift_field,
        common_field,
    )
    if len(common_gcrd) - 1 != EXPECTED_RIGHT_FACTOR_ORDER:
        raise AssertionError("unexpected greatest common right divisor")

    quotient_degree_pairs = [
        [coefficient.numer.degree(), coefficient.denom.degree()]
        for coefficient in quotient
    ]
    if quotient_degree_pairs[0] != [0, 50]:
        raise AssertionError("unexpected trailing quotient coefficient")
    if quotient_degree_pairs[-1] != [0, 50]:
        raise AssertionError("unexpected forward quotient coefficient")

    immediate_boxes = []
    for order, degree in ((8, 71), (7, 72), (7, 80)):
        matrix = differential_matrix(
            sequence,
            prime,
            order,
            degree,
        )
        rank = matrix.rank()
        immediate_boxes.append(
            {
                "differential_order": order,
                "z_degree": degree,
                "rank": rank,
                "column_count": matrix.ncols(),
                "nullity": matrix.ncols() - rank,
            }
        )
        if rank != matrix.ncols():
            raise AssertionError("unexpected relation in a lower box")

    return {
        "prime": prime,
        "differential_operator": {
            "order": ORDER,
            "z_degree": Z_DEGREE,
            "coefficient_count": COEFFICIENT_COUNT,
            "normalization": "coefficient of z^72*d_z^8 is one",
            "support": differential_support(flattened),
            "residual_polynomial_degree": residual_degree,
            "zero_tail_starts_at_z_degree": residual_degree + 1,
            "coefficient_identities_checked": len(residuals),
            "immediate_lower_boxes": immediate_boxes,
        },
        "shift_conversion": {
            "coefficient_index_offset": delta,
            "order": len(shift_rows) - 1,
            "coefficient_m_degrees": shift_degrees,
            "moment_equations_checked": recurrence_checks,
        },
        "ore_comparison": {
            "identity": "R_64,8 = Q_50 * G_14,58",
            "left_quotient_order": len(quotient) - 1,
            "zero_remainder": True,
            "greatest_common_right_divisor_order": len(common_gcrd) - 1,
            "euclidean_remainder_orders": euclidean_orders,
            "left_quotient_rational_degree_pairs": quotient_degree_pairs,
            "forward_quotient_numerator_degree": (
                quotient[-1].numer.degree()
            ),
            "forward_quotient_denominator_degree": (
                quotient[-1].denom.degree()
            ),
        },
    }


def verify_birational_compression() -> dict[str, object]:
    u, t, x, y = symbols("u t x y")
    q = (
        216 - 648*t + 648*t**2 - 216*t**3
        + 91*u - 5*u*t - 263*u*t**2 + 177*u*t**3
        + 23*u**2 + 44*u**2*t + 145*u**2*t**2 - 212*u**2*t**3
        + 11*u**3 - 4*u**3*t + 102*u**3*t**2 + 245*u**3*t**3
        + 13*u**4*t + 5*u**4*t**2 + 131*u**4*t**3
        + 17*u**5*t**2 + 20*u**5*t**3
        + 19*u**6*t**3
    )
    phi = (
        19*x**3*y**3 + 17*x**3*y**2 + 13*x**3*y + 11*x**3
        + 37*x**2*y**3 + 31*x**2*y**2 + 29*x**2*y + 23*x**2
        + 149*x*y**3 + 127*x*y**2 + 113*x*y + 91*x
        + 354*y**3 + 302*y**2 + 268*y + 216
    )
    substituted = cancel(
        (x + y) ** 3
        * q.subs({u: x, t: y / (x + y)})
        / x**3
    )
    if expand(substituted - phi) != 0:
        raise AssertionError("birational beta compression failed")
    return {
        "substitution": "x=u, y=u*t/(1-t)",
        "inverse": "u=x, t=y/(x+y)",
        "potential_identity": "Q/u^3 = Phi/(x+y)^3",
        "measure_identity": "(du/u)dt = dx*dy/(x+y)^2",
        "compact_form": "(x+y)dxdy/((x+y)^3-z*Phi)",
        "phi_term_count": len(Poly(phi, x, y).terms()),
        "phi_total_degree": Poly(phi, x, y).total_degree(),
        "verified_by_exact_symbolic_expansion": True,
    }


def closed_operator_metadata() -> dict[str, object]:
    text = CLOSED_OPERATOR.read_text()
    body = text.removeprefix("PICARD_FUCHS_BEGIN\n").removesuffix(
        "PICARD_FUCHS_END\n"
    )
    powers = [int(value) for value in re.findall(r"\)dz\^(\d+)", body)]
    if powers != list(range(8, 1, -1)):
        raise ValueError("unexpected closed Picard--Fuchs operator syntax")
    if body.count(")dz +") != 1:
        raise ValueError("closed Picard--Fuchs first derivative missing")
    coefficient_chunks = re.split(r"\)dz(?:\^\d+)? \+ \(", body[1:-1])
    z_degrees = []
    for chunk in coefficient_chunks:
        exponents = [int(value) for value in re.findall(r"z\^(\d+)", chunk)]
        if re.search(r"(?:^|[ +\-])\d+\*z(?:[ +\-]|$)", chunk):
            exponents.append(1)
        z_degrees.append(max(exponents, default=0))
    if len(z_degrees) != 9:
        raise ValueError("expected nine closed operator coefficients")
    return {
        "producer": "MultivariateCreativeTelescoping.jl 0.1.3 CRT",
        "cycle_type": "closed-cycle projective Picard--Fuchs operator",
        "differential_order": 8,
        "coefficient_z_degrees_in_descending_derivative_order": z_degrees,
        "whole_file_bytes": CLOSED_OPERATOR.stat().st_size,
        "interpretation_limit": (
            "the interval period has boundary forcing and is not "
            "annihilated by this closed-cycle operator"
        ),
    }


def load_image_cache(path: Path) -> dict[int, list[int]]:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    if payload.get("format") != (
        "two-pair-sic-bidegree33-rank-two-compact-relative-pf-images-v1"
    ):
        raise ValueError("unsupported compact-relative image cache")
    images = {
        int(prime): [int(value) for value in values]
        for prime, values in payload["images"].items()
    }
    if any(len(image) != COEFFICIENT_COUNT for image in images.values()):
        raise ValueError("invalid compact-relative image width")
    return images


def write_image_cache(path: Path, images: dict[int, list[int]]) -> None:
    payload = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-compact-relative-"
            "pf-images-v1"
        ),
        "status": (
            "exact finite modular images of a fitted inhomogeneous "
            "interval differential relation; not an all-order "
            "characteristic-zero certificate"
        ),
        "point": 0,
        "differential_order": ORDER,
        "z_degree": Z_DEGREE,
        "normalization": "coefficient of z^72*d_z^8 is one",
        "coefficient_count": COEFFICIENT_COUNT,
        "maximum_moment": MAXIMUM_MOMENT,
        "fit_extra_rows": FIT_EXTRA_ROWS,
        "images": {
            str(prime): image for prime, image in sorted(images.items())
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


def prime_list(count: int, start: int) -> list[int]:
    primes = []
    prime = start
    for _ in range(count):
        prime = int(nextprime(prime))
        primes.append(prime)
    return primes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-count", type=int, default=1)
    parser.add_argument("--prime-start", type=int, default=1_000_000)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--image-cache", type=Path, default=IMAGE_CACHE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.prime_count < 1:
        raise ValueError("prime-count must be positive")
    if arguments.jobs < 1:
        raise ValueError("jobs must be positive")

    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required")
    requested_primes = prime_list(
        arguments.prime_count,
        arguments.prime_start,
    )
    images = load_image_cache(arguments.image_cache)

    with tempfile.TemporaryDirectory(prefix="sic33-compact-relative-") as path:
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
        missing = [prime for prime in requested_primes if prime not in images]
        if missing:
            with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
                sequences = executor.map(
                    lambda prime: compute_moments(
                        executable,
                        prime,
                        MAXIMUM_MOMENT,
                    ),
                    missing,
                )
                for prime, sequence in zip(missing, sequences, strict=True):
                    image, rank, nullity = fit_differential_operator(
                        sequence,
                        prime,
                    )
                    images[prime] = image
                    write_image_cache(arguments.image_cache, images)
                    print(
                        f"PASS prime={prime} rank={rank} "
                        f"nullity={nullity}",
                        flush=True,
                    )
        reference_prime = requested_primes[0]
        reference_sequence = compute_moments(
            executable,
            reference_prime,
            MAXIMUM_MOMENT,
        )

    reference_image = images[reference_prime]
    fitted_again, fit_rank, fit_nullity = fit_differential_operator(
        reference_sequence,
        reference_prime,
    )
    if fitted_again != reference_image:
        raise AssertionError("cached reference image does not replay")
    analysis = reference_analysis(
        reference_image,
        reference_sequence,
        reference_prime,
    )
    analysis["differential_operator"].update(
        {"fit_rank": fit_rank, "fit_nullity": fit_nullity}
    )

    artifact = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-compact-relative-"
            "pf-research-v1"
        ),
        "status": (
            "exact symbolic compression and finite modular relative "
            "Picard--Fuchs/Ore-factor calculation; not an all-order "
            "characteristic-zero certificate"
        ),
        "point": 0,
        "birational_compression": verify_birational_compression(),
        "closed_cycle_operator": closed_operator_metadata(),
        "reference_modular_analysis": analysis,
        "image_cache": str(arguments.image_cache.relative_to(ROOT)),
        "cached_primes": sorted(images),
        "interpretation": (
            "The interval relation is inhomogeneous: L_8 F is a "
            "polynomial of degree 55. Its homogeneous coefficient tail "
            "is the order-64, m-degree-8 shift relation R. At the "
            "reference prime, R is exactly a left multiple Q_50 G_14 "
            "with zero remainder. Characteristic-zero reconstruction, "
            "a reduction certificate for L_8 F, and the finite initial-"
            "value/nonvanishing argument remain necessary."
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS exact birational beta compression")
    print("PASS closed-cycle operator has differential order 8")
    print("PASS interval residual has degree 55")
    print("PASS R_64,8 = Q_50 * G_14,58 with zero remainder")
    print("PASS greatest common right divisor has order 14")
    print(f"PASS wrote {arguments.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
