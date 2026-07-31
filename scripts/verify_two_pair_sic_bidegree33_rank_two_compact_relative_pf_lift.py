#!/usr/bin/env python3
"""Verify the characteristic-zero compact relative PF/Ore bridge.

This checker replays the simultaneous reconstruction, converts the lifted
order-eight differential relation to its order-64 shift tail, verifies exact
rational moments, and proves the exact right factorization

    R_64,8 = Q_50 G_14,58

in QQ(m)[S; S*m=(m+1)*S].  The factorization plus 50 exact initial zeros and
the positive forward denominator of Q prove that an all-order certificate for
R would imply the stored order-14 recurrence.  This checker does not itself
prove the relative divergence/boundary identity for R.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, factorial, gcd
from pathlib import Path
import sys
from typing import Any

from sympy import Poly, ZZ, symbols
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_sic_bidegree33_rank_two_ore_gcd import (  # noqa: E402
    ShiftOreField,
)


DIFFERENTIAL_LIFT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_pf_"
        "characteristic_zero_lift.json"
    )
)
IMAGE_CACHE = (
    ROOT
    / "artifacts"
    / "local"
    / "two_pair_sic_bidegree33_rank_two_compact_relative_pf_images.json"
)
COMMON_LIFT = (
    ROOT
    / "artifacts"
    / "local"
    / "two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json"
)
SHIFT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_compact_relative_shift_operator.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_pf_"
        "characteristic_zero_verification.json"
    )
)

DIFFERENTIAL_ORDER = 8
Z_DEGREE = 72
SHIFT_ORDER = 64
M_DEGREE = 8
COMMON_ORDER = 14
COMMON_M_DEGREE = 58
QUOTIENT_ORDER = 50
EXACT_MAXIMUM_MOMENT = 100
EXPECTED_RESIDUAL_DEGREE = 55


class RationalShiftOreField(ShiftOreField):
    """Exact QQ(m) specialization of the small shift-Ore engine."""

    def __init__(self) -> None:
        self.prime = 0
        self.field = QQ.frac_field("m")
        self.m = self.field.gens[0]
        self.polynomial_ring = self.m.numer.ring
        self.polynomial_m = self.polynomial_ring.gens[0]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply_cubic(
    polynomial: list[int],
    cubic: list[int],
) -> list[int]:
    result = [0] * (len(polynomial) + 3)
    for left_degree, left_value in enumerate(polynomial):
        for right_degree, right_value in enumerate(cubic):
            result[left_degree + right_degree] += left_value * right_value
    return result


def exact_normalized_moments(maximum: int) -> list[Fraction]:
    """Rank-two formula for nu_m, avoiding a full bivariate power."""

    a1 = [1, 0, 2, 5]
    a2 = [0, 1, 3, 7]
    p1 = [11, 13, 17, 19]
    p2 = [23, 29, 31, 37]
    dual = [[1]]
    coordinate = [[1]]
    sequence = []
    for order in range(maximum + 1):
        factorials = [factorial(index) for index in range(3 * order + 2)]
        moment = 0
        for inner in range(order + 1):
            pairing = sum(
                dual[inner][degree]
                * coordinate[inner][degree]
                * factorials[degree]
                * factorials[3 * order - degree]
                for degree in range(3 * order + 1)
            )
            moment += comb(order, inner) * pairing
        sequence.append(Fraction(moment, factorials[3 * order + 1]))
        if order == maximum:
            break
        next_dual = [None] * (order + 2)
        next_coordinate = [None] * (order + 2)
        for inner in range(order + 1):
            next_dual[inner] = multiply_cubic(dual[inner], a2)
            next_coordinate[inner] = multiply_cubic(
                coordinate[inner],
                p2,
            )
        next_dual[order + 1] = multiply_cubic(dual[order], a1)
        next_coordinate[order + 1] = multiply_cubic(
            coordinate[order],
            p1,
        )
        dual = next_dual
        coordinate = next_coordinate
    return sequence


def multiply_linear(polynomial: list[int], constant: int) -> list[int]:
    result = [0] * (len(polynomial) + 1)
    for exponent, coefficient in enumerate(polynomial):
        result[exponent] += constant * coefficient
        result[exponent + 1] += coefficient
    return result


def differential_to_shift(
    coefficients: list[list[int]],
) -> list[list[int]]:
    delta = max(
        z_exponent - derivative_order
        for derivative_order, row in enumerate(coefficients)
        for z_exponent, value in enumerate(row)
        if value
    )
    if delta != SHIFT_ORDER:
        raise RuntimeError(f"unexpected coefficient offset {delta}")
    shift_coefficients = [
        [0] * (M_DEGREE + 1) for _ in range(SHIFT_ORDER + 1)
    ]
    for derivative_order, row in enumerate(coefficients):
        for z_exponent, value in enumerate(row):
            if value == 0:
                continue
            shift = delta - z_exponent + derivative_order
            falling = [1]
            for offset in range(derivative_order):
                falling = multiply_linear(falling, shift - offset)
            for exponent, factor in enumerate(falling):
                shift_coefficients[shift][exponent] += value * factor
    if any(row[-1] == 0 for row in shift_coefficients):
        raise RuntimeError("a shift coefficient lost m-degree eight")
    return shift_coefficients


def evaluate(polynomial: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def falling(value: int, order: int) -> int:
    result = 1
    for offset in range(order):
        result *= value - offset
    return result


def polynomial_operator(
    ore: RationalShiftOreField,
    coefficient_rows: list[list[int]],
) -> list[Any]:
    domain = ore.polynomial_ring.domain
    return [
        ore.field(
            ore.polynomial_ring.from_dict(
                {
                    (exponent,): domain.convert(value)
                    for exponent, value in enumerate(row)
                    if value
                }
            )
        )
        for row in coefficient_rows
    ]


def write_shift_operator(
    coefficients: list[list[int]],
    differential_hash: str,
) -> None:
    payload = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-compact-relative-"
            "shift-operator-v1"
        ),
        "status": (
            "exact characteristic-zero coefficient conversion from the "
            "reconstructed compact differential operator; its all-order "
            "relative divergence and endpoint identities require a "
            "separate certificate"
        ),
        "point": 0,
        "operator": {
            "order": SHIFT_ORDER,
            "m_degree": M_DEGREE,
            "coefficient_count": (SHIFT_ORDER + 1) * (M_DEGREE + 1),
        },
        "ore_convention": "S*f(m)=f(m+1)*S",
        "coefficient_index_conversion": (
            "z^e*d_z^k maps at n=m+64 to "
            "S^(64-e+k)*(m+64-e+k)_falling_k"
        ),
        "primitive_integer_coefficients": coefficients,
        "source_differential_lift_sha256": differential_hash,
    }
    SHIFT_OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-exact-ore-division",
        action="store_true",
        help="run the fast exact checks but omit the roughly eight-minute division",
    )
    arguments = parser.parse_args()

    lift = json.loads(DIFFERENTIAL_LIFT.read_text())
    cache = json.loads(IMAGE_CACHE.read_text())
    common_lift = json.loads(COMMON_LIFT.read_text())
    if lift.get("status") != "stable simultaneous rational reconstruction":
        raise RuntimeError("compact differential lift is not stable")
    if lift.get("kind") != "compact":
        raise RuntimeError("wrong reconstructed operator kind")
    if lift.get("operator") != {
        "order": DIFFERENTIAL_ORDER,
        "m_degree": Z_DEGREE,
        "coefficient_count": (DIFFERENTIAL_ORDER + 1) * (Z_DEGREE + 1),
    }:
        raise RuntimeError("unexpected differential operator shape")

    denominator = int(lift["common_denominator"])
    differential = [
        [int(value) for value in row]
        for row in lift["primitive_integer_coefficients"]
    ]
    flat = [value for row in differential for value in row]
    if differential[-1][-1] != denominator:
        raise RuntimeError("differential normalization mismatch")
    content = 0
    for value in flat:
        content = gcd(content, abs(value))
    if content != 1:
        raise RuntimeError("differential operator is not primitive")

    images = {
        int(prime): [int(value) for value in image]
        for prime, image in cache["images"].items()
    }
    requested_primes = [
        int(prime)
        for prime in lift["build_primes"] + lift["holdout_primes"]
    ]
    for prime in requested_primes:
        inverse = pow(denominator, -1, prime)
        reduction = [value * inverse % prime for value in flat]
        if reduction != images.get(prime):
            raise RuntimeError(f"differential image mismatch at {prime}")

    shift = differential_to_shift(differential)
    shift_content = 0
    for row in shift:
        for value in row:
            shift_content = gcd(shift_content, abs(value))
    if shift_content != 1:
        raise RuntimeError("shift operator is not primitive")
    differential_hash = sha256(DIFFERENTIAL_LIFT)
    write_shift_operator(shift, differential_hash)

    moments = exact_normalized_moments(EXACT_MAXIMUM_MOMENT)
    residuals = []
    for coefficient_index in range(
        EXACT_MAXIMUM_MOMENT - DIFFERENTIAL_ORDER + 1
    ):
        total = Fraction(0)
        for derivative_order, row in enumerate(differential):
            for z_exponent, coefficient in enumerate(row):
                if coefficient == 0 or z_exponent > coefficient_index:
                    continue
                moment_index = (
                    coefficient_index - z_exponent + derivative_order
                )
                total += (
                    coefficient
                    * falling(moment_index, derivative_order)
                    * moments[moment_index]
                )
        residuals.append(total)
    residual_support = [
        index for index, value in enumerate(residuals) if value
    ]
    if residual_support != list(range(EXPECTED_RESIDUAL_DEGREE + 1)):
        raise RuntimeError("unexpected exact differential residual support")

    shift_rows = EXACT_MAXIMUM_MOMENT - SHIFT_ORDER + 1
    for moment_index in range(shift_rows):
        total = sum(
            Fraction(evaluate(polynomial, moment_index))
            * moments[moment_index + shift_index]
            for shift_index, polynomial in enumerate(shift)
        )
        if total:
            raise RuntimeError(
                f"exact compact shift mismatch at m={moment_index}"
            )

    common = [
        [int(value) for value in row]
        for row in common_lift["primitive_integer_coefficients"]
    ]
    if len(common) != COMMON_ORDER + 1 or any(
        len(row) != COMMON_M_DEGREE + 1 for row in common
    ):
        raise RuntimeError("unexpected common operator shape")
    for moment_index in range(QUOTIENT_ORDER):
        total = sum(
            Fraction(evaluate(polynomial, moment_index))
            * moments[moment_index + shift_index]
            for shift_index, polynomial in enumerate(common)
        )
        if total:
            raise RuntimeError(
                f"exact common initial row mismatch at m={moment_index}"
            )

    m = symbols("m")
    shift_forward = Poly(
        sum(value * m**degree for degree, value in enumerate(shift[-1])),
        m,
        domain=ZZ,
    )
    common_forward_shifted = Poly(
        sum(
            value * (m + QUOTIENT_ORDER) ** degree
            for degree, value in enumerate(common[-1])
        ),
        m,
        domain=ZZ,
    )
    forward_gcd = shift_forward.gcd(common_forward_shifted)
    forward_numerator = shift_forward.exquo(forward_gcd)
    forward_denominator = common_forward_shifted.exquo(forward_gcd)
    if forward_numerator.degree() != 0 or forward_denominator.degree() != 50:
        raise RuntimeError("unexpected forward quotient degrees")
    denominator_coefficients = [
        int(value) for value in forward_denominator.all_coeffs()
    ]
    if not all(value > 0 for value in denominator_coefficients):
        raise RuntimeError("forward quotient denominator is not positive")

    quotient_degree_pairs = None
    exact_factorization = False
    if not arguments.skip_exact_ore_division:
        ore = RationalShiftOreField()
        quotient, remainder = ore.left_division(
            polynomial_operator(ore, shift),
            polynomial_operator(ore, common),
        )
        if remainder or len(quotient) - 1 != QUOTIENT_ORDER:
            raise RuntimeError("exact shift-Ore factorization failed")
        quotient_degree_pairs = [
            [coefficient.numer.degree(), coefficient.denom.degree()]
            for coefficient in quotient
        ]
        expected_pairs = [[0, 50]] + [[50, 100]] * 49 + [[0, 50]]
        if quotient_degree_pairs != expected_pairs:
            raise RuntimeError("unexpected exact quotient degree pattern")
        exact_factorization = True

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-compact-relative-"
            "pf-characteristic-zero-verification-v1"
        ),
        "status": (
            "exact characteristic-zero compact relative PF/Ore bridge"
            if exact_factorization
            else "fast exact checks; full rational Ore division skipped"
        ),
        "point": 0,
        "differential_operator": {
            "order": DIFFERENTIAL_ORDER,
            "z_degree": Z_DEGREE,
            "maximum_primitive_coefficient_bits": lift[
                "successful_lattice"
            ]["maximum_coefficient_bits"],
            "build_prime_count": len(lift["build_primes"]),
            "holdout_prime_count": len(lift["holdout_primes"]),
            "exact_residual_support": [0, EXPECTED_RESIDUAL_DEGREE],
            "exact_coefficient_rows_checked": len(residuals),
        },
        "shift_operator": {
            "order": SHIFT_ORDER,
            "m_degree": M_DEGREE,
            "maximum_primitive_coefficient_bits": max(
                abs(value).bit_length() for row in shift for value in row
            ),
            "exact_moment_rows_checked": shift_rows,
            "artifact": str(SHIFT_OUTPUT.relative_to(ROOT)),
        },
        "ore_factorization": {
            "identity": "R_64,8 = Q_50 * G_14,58",
            "exact_rational_division_performed": exact_factorization,
            "zero_remainder": exact_factorization,
            "left_quotient_order": QUOTIENT_ORDER,
            "left_quotient_degree_pairs": quotient_degree_pairs,
            "forward_coefficient": {
                "numerator_degree": forward_numerator.degree(),
                "denominator_degree": forward_denominator.degree(),
                "denominator_all_51_coefficients_positive": True,
                "nonzero_for_every_integer_m_at_least_zero": True,
            },
        },
        "initial_value_closure": {
            "exact_G14_rows": [0, QUOTIENT_ORDER - 1],
            "exact_initial_zero_count": QUOTIENT_ORDER,
            "consequence": (
                "if R annihilates all moments, the order-50 forward "
                "recurrence for G14(moment) and its 50 initial zeros "
                "force G14(moment)=0 for every m>=0"
            ),
        },
        "remaining_gate": (
            "prove the exact relative divergence identity for R and its "
            "two endpoint traces; operator reconstruction, exact "
            "factorization, forward nonvanishing, and initial values close"
        ),
        "files_sha256": {
            str(DIFFERENTIAL_LIFT.relative_to(ROOT)): differential_hash,
            str(SHIFT_OUTPUT.relative_to(ROOT)): sha256(SHIFT_OUTPUT),
            str(COMMON_LIFT.relative_to(ROOT)): sha256(COMMON_LIFT),
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS replayed 95 reconstruction primes and five holdouts")
    print("PASS exact residual is supported in degrees 0 through 55")
    print(f"PASS {shift_rows} exact R_64,8 moment rows")
    print("PASS 50 exact initial G_14 rows")
    print("PASS Q_50 forward denominator has 51 positive coefficients")
    if exact_factorization:
        print("PASS exact R_64,8 = Q_50 * G_14,58 with zero remainder")
    else:
        print("PASS exact Ore division deliberately skipped")
    print(f"PASS wrote {SHIFT_OUTPUT.relative_to(ROOT)}")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
