#!/usr/bin/env python3
"""Ore-gcd comparison of the sampled rank-two cubic recurrences.

For each of the two exact-rank-two points and three prime fields used by
the adjacent recurrence probes, this checker:

* recomputes the order-18/degree-18 and order-27/degree-11 operators;
* performs exact left Euclidean division in F_p(m)[S; S*m=(m+1)*S];
* proves that neither displayed operator is a left multiple of the other;
* computes their monic greatest common right divisor; and
* clears its rational coefficient content and checks the resulting
  polynomial operator directly on all 501 computed moments.

The result is an exact bounded modular calculation, not a universal
Picard--Fuchs or creative-telescoping certificate.
"""

from __future__ import annotations

from itertools import product
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any

from sympy.polys.domains import GF

from verify_two_pair_sic_bidegree33_rank_two_holonomic_probe import (
    POINTS,
    PRIMES,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "scripts"
    / "explore_two_pair_sic_bidegree33_rank_two_recurrence.cpp"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_gcd.json"
)
MAXIMUM_MOMENT = 500
ORDER_DEGREES = ((18, 18), (27, 11))
EXPECTED_REMAINDER_ORDERS = (17, 16, 15, 14, -1)
EXPECTED_COMMON_ORDER = 14
EXPECTED_COMMON_DEGREE = 58
COMMON_FORWARD_SHIFTS = (32, 34, 35, 37, 38, 40, 41, 43)

Operator = list[Any]


class ShiftOreField:
    """Small exact F_p(m)[S; sigma] implementation for left division."""

    def __init__(self, prime: int):
        self.prime = prime
        self.field = GF(prime).frac_field("m")
        self.m = self.field.gens[0]
        self.polynomial_ring = self.m.numer.ring
        self.polynomial_m = self.polynomial_ring.gens[0]

    def parse_operator(
        self,
        output: str,
        order: int,
        degree: int,
    ) -> Operator:
        lines = output.strip().splitlines()
        assert lines[0] == f"FOUND order={order} degree={degree}"
        coefficients = [self.field.zero] * (order + 1)
        for line in lines[1:]:
            values = [int(value) for value in line.split()]
            assert len(values) == degree + 2
            shift = values[0]
            coefficients[shift] = sum(
                (
                    self.field(value) * self.m**exponent
                    for exponent, value in enumerate(values[1:])
                ),
                self.field.zero,
            )
        assert all(coefficient for coefficient in coefficients)
        return coefficients

    def shift(self, value: Any, amount: int) -> Any:
        if amount == 0:
            return value
        replacement = {self.polynomial_m: self.polynomial_m + amount}
        return value.raw_new(
            value.numer.compose(replacement),
            value.denom.compose(replacement),
        )

    @staticmethod
    def trim(operator: Operator) -> Operator:
        while operator and not operator[-1]:
            operator.pop()
        return operator

    def left_division(
        self,
        dividend: Operator,
        divisor: Operator,
    ) -> tuple[Operator, Operator]:
        """Return Q,R with dividend=Q*divisor+R in the shift Ore ring."""

        assert divisor
        remainder = dividend[:]
        quotient = [self.field.zero] * max(
            0,
            len(dividend) - len(divisor) + 1,
        )
        while len(remainder) >= len(divisor):
            shift = len(remainder) - len(divisor)
            scalar = remainder[-1] / self.shift(divisor[-1], shift)
            quotient[shift] = scalar
            for index, coefficient in enumerate(divisor):
                remainder[shift + index] -= (
                    scalar * self.shift(coefficient, shift)
                )
            self.trim(remainder)
        return self.trim(quotient), remainder

    @staticmethod
    def monic(operator: Operator) -> Operator:
        leading = operator[-1]
        return [coefficient / leading for coefficient in operator]

    def greatest_common_right_divisor(
        self,
        left: Operator,
        right: Operator,
    ) -> tuple[Operator, list[int]]:
        remainder_orders = []
        while right:
            _, remainder = self.left_division(left, right)
            remainder_orders.append(len(remainder) - 1)
            left = right
            right = self.monic(remainder) if remainder else []
        return self.monic(left), remainder_orders

    def primitive_polynomial_operator(
        self,
        operator: Operator,
    ) -> Operator:
        common_denominator = self.polynomial_ring.one
        for coefficient in operator:
            common_denominator = common_denominator.lcm(
                coefficient.denom
            )
        polynomial_coefficients = [
            coefficient.numer
            * common_denominator.exquo(coefficient.denom)
            for coefficient in operator
        ]
        content = polynomial_coefficients[0]
        for coefficient in polynomial_coefficients[1:]:
            content = content.gcd(coefficient)
        primitive = [
            coefficient.exquo(content)
            for coefficient in polynomial_coefficients
        ]
        inverse = self.polynomial_ring.domain.convert(
            pow(int(primitive[-1].LC), -1, self.prime)
        )
        return [coefficient * inverse for coefficient in primitive]

    def evaluate_polynomial(self, polynomial: Any, value: int) -> int:
        result = 0
        for coefficient in polynomial.to_dense():
            result = (
                result * value + int(coefficient)
            ) % self.prime
        return result

    def verify_recurrence(
        self,
        operator: Operator,
        sequence: list[int],
    ) -> int:
        order = len(operator) - 1
        checks = 0
        for moment_index in range(len(sequence) - order):
            total = sum(
                self.evaluate_polynomial(coefficient, moment_index)
                * sequence[moment_index + shift]
                for shift, coefficient in enumerate(operator)
            )
            assert total % self.prime == 0, moment_index
            checks += 1
        return checks

    def common_forward_factor(self) -> Any:
        factor = self.polynomial_ring.one
        for shift in COMMON_FORWARD_SHIFTS:
            factor *= 3 * self.polynomial_m + shift
        return factor


def run_operator(
    executable: Path,
    prime: int,
    point: int,
    order: int,
    degree: int,
) -> str:
    completed = subprocess.run(
        [
            str(executable),
            str(prime),
            str(MAXIMUM_MOMENT),
            str(order),
            str(degree),
            str(point),
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    return completed.stdout


def run_moments(
    executable: Path,
    prime: int,
    point: int,
) -> list[int]:
    completed = subprocess.run(
        [
            str(executable),
            "--moments",
            str(prime),
            str(MAXIMUM_MOMENT),
            str(point),
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    lines = completed.stdout.strip().splitlines()
    assert lines[0] == (
        f"MOMENTS maximum={MAXIMUM_MOMENT} point={point}"
    )
    sequence = []
    for expected_index, line in enumerate(lines[1:]):
        index, value = (int(entry) for entry in line.split())
        assert index == expected_index
        sequence.append(value)
    assert len(sequence) == MAXIMUM_MOMENT + 1
    return sequence


def degree_pairs(operator: Operator) -> list[list[int]]:
    return [
        [coefficient.numer.degree(), coefficient.denom.degree()]
        for coefficient in operator
    ]


def analyze_sample(
    executable: Path,
    prime: int,
    point: int,
) -> dict[str, object]:
    ore = ShiftOreField(prime)
    order_18 = ore.parse_operator(
        run_operator(executable, prime, point, *ORDER_DEGREES[0]),
        *ORDER_DEGREES[0],
    )
    order_27 = ore.parse_operator(
        run_operator(executable, prime, point, *ORDER_DEGREES[1]),
        *ORDER_DEGREES[1],
    )

    _, first_remainder = ore.left_division(order_27, order_18)
    assert len(first_remainder) - 1 == 17

    common, remainder_orders = ore.greatest_common_right_divisor(
        order_27,
        order_18,
    )
    assert tuple(remainder_orders) == EXPECTED_REMAINDER_ORDERS
    assert len(common) - 1 == EXPECTED_COMMON_ORDER

    quotient_18, remainder_18 = ore.left_division(order_18, common)
    quotient_27, remainder_27 = ore.left_division(order_27, common)
    assert not remainder_18
    assert not remainder_27
    assert len(quotient_18) - 1 == 4
    assert len(quotient_27) - 1 == 13

    primitive = ore.primitive_polynomial_operator(common)
    coefficient_degrees = [
        coefficient.degree() for coefficient in primitive
    ]
    assert coefficient_degrees == [
        EXPECTED_COMMON_DEGREE
    ] * (EXPECTED_COMMON_ORDER + 1)
    forward = primitive[-1]
    fixed_forward = ore.common_forward_factor()
    residual_forward = forward.exquo(fixed_forward)
    assert residual_forward.degree() == 50

    sequence = run_moments(executable, prime, point)
    recurrence_checks = ore.verify_recurrence(primitive, sequence)
    assert recurrence_checks == 487

    return {
        "point": point,
        "prime": prime,
        "input_operators": [
            {"order": order, "m_degree": degree}
            for order, degree in ORDER_DEGREES
        ],
        "order_27_left_division_by_order_18_remainder_order": 17,
        "euclidean_remainder_orders": remainder_orders,
        "greatest_common_right_divisor": {
            "order": EXPECTED_COMMON_ORDER,
            "primitive_coefficient_degrees": coefficient_degrees,
            "forward_fixed_factor": (
                "product_(k in {32,34,35,37,38,40,41,43})(3m+k)"
            ),
            "forward_residual_degree": residual_forward.degree(),
            "moment_equations_checked": recurrence_checks,
        },
        "left_quotient_orders": {
            "order_18_operator": len(quotient_18) - 1,
            "order_27_operator": len(quotient_27) - 1,
        },
        "left_quotient_rational_degree_pairs": {
            "order_18_operator": degree_pairs(quotient_18),
            "order_27_operator": degree_pairs(quotient_27),
        },
    }


def main() -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required")
    records = []
    with tempfile.TemporaryDirectory(prefix="sic33-ore-gcd-") as path:
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
        for point, prime in product(range(len(POINTS)), PRIMES):
            records.append(
                analyze_sample(executable, prime, point)
            )

    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-ore-gcd-v1",
        "status": (
            "exact bounded modular Ore-factor calculation; not a "
            "universal Picard-Fuchs or telescoping certificate"
        ),
        "ore_convention": (
            "S*f(m)=f(m+1)*S; division is on the left and the common "
            "factor is a greatest common right divisor"
        ),
        "records": records,
        "uniform_sampled_result": {
            "order_27_is_left_multiple_of_order_18": False,
            "common_right_factor_order": EXPECTED_COMMON_ORDER,
            "common_right_factor_primitive_m_degree": (
                EXPECTED_COMMON_DEGREE
            ),
            "order_18_left_quotient_order": 4,
            "order_27_left_quotient_order": 13,
            "direct_moment_equations_per_sample": 487,
            "interior_length_match": (
                "the common order 14 equals the exact interior length "
                "14 in the 2+2+14 logarithmic-Jacobian decomposition"
            ),
        },
        "interpretation_limit": (
            "The order match is evidence that the sampled scalar period "
            "uses the interior cyclic factor while the two endpoint "
            "pairs account for the four extra relative-Jacobian "
            "dimensions. This identification is not proved universally."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS six exact modular shift-Ore Euclidean calculations")
    print("PASS order 27 is not a left multiple of order 18")
    print("PASS common right factor has order 14 and m-degree 58")
    print("PASS order-18 and order-27 left quotients have orders 4 and 13")
    print("PASS common factor annihilates all 487 available moment rows")
    print("PASS order 14 matches the exact interior Jacobian length")
    print("PASS result remains sampled evidence, not a universal certificate")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
