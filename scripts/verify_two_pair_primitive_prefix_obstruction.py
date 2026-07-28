#!/usr/bin/env python3
"""Exact replay of the all-degree primitive finite-prefix obstruction."""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial
from pathlib import Path

from verify_two_pair_image_mathieu_counterexample import (
    ZERO,
    add,
    contraction,
    double_factorial_odd,
    monomial,
    multiply,
    power,
    scale,
    witness,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_primitive_prefix_obstruction.json"
)
DEGREES = range(4, 9)
PHASE_TWO_NUMERATOR_COEFFICIENTS = (
    225434243432448,
    2353938430951424,
    -34631484170043392,
    -967600078304837632,
    -10586540331031855104,
    -73015690816373391360,
    -360309321821207396352,
    -1348856283083292540928,
    -3958268002248099233792,
    -9285501721191809486848,
    -17615580262001968017408,
    -27188834549765293099008,
    -34199757839826016106496,
    -34988488014787406053120,
    -28953216357761535637760,
    -19199991306428683987072,
    -10063262895913229883456,
    -4086905609184770420496,
    -1249692217173815594688,
    -275573767406411075760,
    -40876542617863385700,
    -3589300834376709375,
    -137104131892612500,
)


def evaluate_integer_polynomial(
    coefficients: tuple[int, ...],
    value: int,
) -> int:
    result = 0
    for coefficient in coefficients:
        result = result * value + coefficient
    return result


def phase_two_angular_residual(height_half: int) -> Fraction:
    h = height_half
    numerator_polynomial = evaluate_integer_polynomial(
        PHASE_TWO_NUMERATOR_COEFFICIENTS,
        h,
    )
    numerator = (
        -16
        * (4 * h + 1)
        * (4 * h + 3)
        * (4 * h + 5)
        * numerator_polynomial
    )
    denominator = (
        81
        * (2 * h + 1) ** 3
        * (2 * h + 3) ** 4
        * (2 * h + 5) ** 4
        * (2 * h + 7) ** 2
        * (2 * h + 11)
        * (4 * h + 9)
        * (4 * h + 11)
        * (4 * h + 13)
        * (6 * h + 5)
        * (6 * h + 7)
        * (6 * h + 11)
        * (8 * h + 1)
        * (8 * h + 3)
        * (8 * h + 5)
        * (8 * h + 7)
        * (8 * h + 9)
        * (8 * h * h + 18 * h + 29) ** 2
    )
    return Fraction(numerator, denominator)


def main() -> None:
    f, _, generators = witness()
    r = generators["R"]
    z = generators["Z"]
    replay: dict[str, dict[str, str]] = {}
    triangular_replay: dict[str, dict[str, str]] = {}

    for degree in DEGREES:
        radial_seed = multiply(power(r, degree - 4), f)
        primitive = add(radial_seed, power(z, degree))
        primitive_power = monomial(ZERO)
        replay[str(degree)] = {}
        for order in range(1, degree + 2):
            primitive_power = multiply(primitive_power, primitive)
            actual = contraction(primitive_power)
            if order <= degree:
                assert actual == {}
                replay[str(degree)][str(order)] = "0"
            else:
                expected = Fraction(
                    (degree + 1)
                    * factorial(degree * (degree + 1) + 1)
                    * factorial(degree),
                    double_factorial_odd(2 * degree + 1),
                )
                assert actual == {ZERO: expected}
                replay[str(degree)][str(order)] = str(expected)

    for degree in range(4, 8):
        radial_seed = multiply(power(r, degree - 4), f)
        triangular_replay[str(degree)] = {}
        for phase in range(1, degree + 1):
            correction = multiply(
                power(r, degree - phase),
                power(z, phase),
            )
            candidate = add(radial_seed, correction)
            candidate_power = monomial(ZERO)
            for order in range(1, phase + 2):
                candidate_power = multiply(candidate_power, candidate)
                actual = contraction(candidate_power)
                if order <= phase:
                    assert actual == {}
                else:
                    expected = Fraction(
                        (phase + 1)
                        * factorial(degree * (phase + 1) + 1)
                        * factorial(phase),
                        double_factorial_odd(2 * phase + 1),
                    )
                    assert actual == {ZERO: expected}
                    triangular_replay[str(degree)][str(phase)] = str(
                        expected
                    )

    degree_five_seed = multiply(r, f)
    z5 = power(z, 5)
    w5 = power(generators["W"], 5)
    extreme_constants = {
        "mu2_ab": int(2 * contraction(multiply(z5, w5))[ZERO]),
        "mu4_b": int(
            4 * contraction(multiply(power(degree_five_seed, 3), w5))[ZERO]
        ),
        "mu6_a": int(
            6 * contraction(multiply(power(degree_five_seed, 5), z5))[ZERO]
        ),
    }
    assert extreme_constants == {
        "mu2_ab": 921600,
        "mu4_b": -12009388769280000,
        "mu6_a": 569547266090245736292679680000000,
    }

    odd_height_replay: dict[str, str] = {}
    t = generators["T"]
    for degree in range(4, 8):
        radial_seed = multiply(power(r, degree - 4), f)
        correction = multiply(power(z, degree - 1), t)
        candidate = add(radial_seed, correction)
        candidate_power = monomial(ZERO)
        for order in range(1, 2 * degree + 1):
            candidate_power = multiply(candidate_power, candidate)
            actual = contraction(candidate_power)
            if order < 2 * degree:
                assert actual == {}
            else:
                expected = Fraction(
                    (2 * degree * (2 * degree - 1) // 2)
                    * factorial(2 * degree * degree + 1)
                    * factorial(2 * degree - 2),
                    double_factorial_odd(4 * degree - 1),
                )
                assert actual == {ZERO: expected}
                odd_height_replay[str(degree)] = str(expected)

    opposite_odd_height_replay: dict[str, dict[str, list[str]]] = {}
    w = generators["W"]
    for degree in range(4, 8):
        radial_seed = multiply(power(r, degree - 4), f)
        opposite_odd_height_replay[str(degree)] = {}
        for height in range(1, degree, 2):
            phase = degree - height
            positive = multiply(power(z, phase), power(t, height))
            negative = multiply(power(w, phase), power(t, height))
            cross = 2 * contraction(multiply(positive, negative))[ZERO]
            positive_order = 2 * phase + 2
            positive_value = (
                (positive_order * (positive_order - 1) // 2)
                * contraction(
                    multiply(
                        power(radial_seed, 2 * phase),
                        power(positive, 2),
                    )
                )[ZERO]
            )
            negative_order = phase + 2
            negative_value = (
                (negative_order * (negative_order - 1) // 2)
                * contraction(
                    multiply(
                        power(radial_seed, phase),
                        power(negative, 2),
                    )
                )[ZERO]
            )
            expected_cross = Fraction(
                2
                * factorial(2 * degree + 1)
                * factorial(phase)
                * double_factorial_odd(2 * height - 1),
                double_factorial_odd(2 * degree + 1),
            )
            expected_positive = Fraction(
                (positive_order * (positive_order - 1) // 2)
                * factorial(degree * positive_order + 1)
                * factorial(2 * phase)
                * double_factorial_odd(2 * height - 1),
                double_factorial_odd(4 * phase + 2 * height + 1),
            )
            expected_negative = Fraction(
                (negative_order * (negative_order - 1) // 2)
                * factorial(degree * negative_order + 1)
                * (-1) ** phase
                * factorial(2 * phase)
                * double_factorial_odd(2 * degree - 1),
                2**phase
                * double_factorial_odd(4 * phase + 2 * degree + 1),
            )
            assert cross == expected_cross
            assert positive_value == expected_positive
            assert negative_value == expected_negative
            opposite_odd_height_replay[str(degree)][str(height)] = [
                str(cross),
                str(positive_value),
                str(negative_value),
            ]

    opposite_even_height_replay: dict[str, dict[str, list[str]]] = {}
    for degree in range(4, 8):
        radial_seed = multiply(power(r, degree - 4), f)
        opposite_even_height_replay[str(degree)] = {}
        for height in range(0, degree + 1, 2):
            phase = degree - height
            if phase < 3:
                continue
            positive = multiply(power(z, phase), power(t, height))
            negative = multiply(power(w, phase), power(t, height))
            cross = 2 * contraction(multiply(positive, negative))[ZERO]
            positive_order = phase + 1
            positive_value = positive_order * contraction(
                multiply(power(radial_seed, phase), positive)
            )[ZERO]
            half_phase = (phase + 1) // 2
            negative_order = half_phase + 1
            negative_value = negative_order * contraction(
                multiply(power(radial_seed, half_phase), negative)
            )[ZERO]
            gamma = 1 if phase % 2 == 0 else 3 * half_phase
            expected_cross = Fraction(
                2
                * factorial(2 * degree + 1)
                * factorial(phase)
                * double_factorial_odd(2 * height - 1),
                double_factorial_odd(2 * degree + 1),
            )
            expected_positive = Fraction(
                positive_order
                * factorial(degree * positive_order + 1)
                * factorial(phase)
                * double_factorial_odd(height - 1),
                double_factorial_odd(2 * phase + height + 1),
            )
            expected_negative = Fraction(
                negative_order
                * factorial(degree * negative_order + 1)
                * gamma
                * (-1) ** half_phase
                * factorial(phase)
                * double_factorial_odd(height + 2 * half_phase - 1),
                2**half_phase
                * double_factorial_odd(
                    height + 2 * half_phase + 2 * phase + 1
                ),
            )
            assert cross == expected_cross
            assert positive_value == expected_positive
            assert negative_value == expected_negative
            opposite_even_height_replay[str(degree)][str(height)] = [
                str(cross),
                str(positive_value),
                str(negative_value),
            ]

    phase_one_positive = multiply(z, power(t, 4))
    phase_one_negative = multiply(w, power(t, 4))
    phase_one_terms = (degree_five_seed, phase_one_positive, phase_one_negative)
    phase_one_moments: dict[str, dict[str, int]] = {}
    expected_phase_one = {
        2: (
            11520,
            {(1, 1): 70, (1, 0): 198, (0, 1): -165},
        ),
        3: (
            -1393459200,
            {
                (1, 1): 490,
                (1, 0): 858,
                (0, 2): 21,
                (0, 1): -845,
            },
        ),
        4: (
            10534551552000,
            {
                (2, 2): 8580,
                (2, 1): 35112,
                (2, 0): 45220,
                (1, 2): -36036,
                (1, 1): 45220,
                (1, 0): 100776,
                (0, 2): 52269,
                (0, 1): -142120,
            },
        ),
    }
    for order, (scalar, expected_coefficients) in expected_phase_one.items():
        actual_coefficients: dict[tuple[int, int], int] = {}
        for positive_count in range(order + 1):
            for negative_count in range(order - positive_count + 1):
                seed_count = order - positive_count - negative_count
                value = contraction(
                    multiply(
                        multiply(
                            power(phase_one_terms[0], seed_count),
                            power(phase_one_terms[1], positive_count),
                        ),
                        power(phase_one_terms[2], negative_count),
                    )
                ).get(ZERO, Fraction(0))
                if not value:
                    continue
                multinomial = Fraction(
                    factorial(order),
                    factorial(seed_count)
                    * factorial(positive_count)
                    * factorial(negative_count),
                )
                coefficient = multinomial * value
                assert coefficient.denominator == 1
                actual_coefficients[(positive_count, negative_count)] = int(
                    coefficient
                )
        assert actual_coefficients == {
            exponent: scalar * coefficient
            for exponent, coefficient in expected_coefficients.items()
        }
        phase_one_moments[str(order)] = {
            f"a^{positive_count}b^{negative_count}": coefficient
            for (positive_count, negative_count), coefficient in sorted(
                actual_coefficients.items()
            )
        }
    phase_one_resultant = (
        735 * 54155455**2
        - 12929 * 155631189 * 54155455
        - 12870 * 155631189**2
    )
    assert phase_one_resultant == -418538718730248905250

    phase_two_modulus = 47
    phase_two_residues = [
        evaluate_integer_polynomial(
            PHASE_TWO_NUMERATOR_COEFFICIENTS,
            residue,
        )
        % phase_two_modulus
        for residue in range(phase_two_modulus)
    ]
    assert all(phase_two_residues)

    phase_two_replay: dict[str, dict[str, str]] = {}
    for degree in range(6, 10, 2):
        height_half = (degree - 2) // 2
        radial_seed = multiply(power(r, degree - 4), f)
        positive = multiply(power(z, 2), power(t, 2 * height_half))
        negative = multiply(power(w, 2), power(t, 2 * height_half))
        positive_coefficient = Fraction(
            (4 * height_half + 1)
            * (4 * height_half + 3)
            * (4 * height_half + 5),
            2
            * (2 * height_half + 3)
            * (2 * height_half + 5)
            * (2 * height_half + 7),
        )
        negative_coefficient = Fraction(
            -2
            * (2 * height_half + 9)
            * (4 * height_half + 1)
            * (4 * height_half + 3)
            * (4 * height_half + 5)
            * (4 * height_half + 7),
            9
            * (2 * height_half + 1)
            * (2 * height_half + 3)
            * (2 * height_half + 5)
            * (
                8 * height_half * height_half
                + 18 * height_half
                + 29
            ),
        )
        exceptional_candidate = add(
            radial_seed,
            scale(positive_coefficient, positive),
            scale(negative_coefficient, negative),
        )
        assert contraction(power(exceptional_candidate, 2)) == {}
        assert contraction(power(exceptional_candidate, 3)) == {}
        fourth_moment = contraction(power(exceptional_candidate, 4))[ZERO]
        expected_fourth_moment = (
            factorial(4 * degree + 1)
            * phase_two_angular_residual(height_half)
        )
        assert fourth_moment == expected_fourth_moment
        assert fourth_moment
        phase_two_replay[str(degree)] = {
            "a": str(positive_coefficient),
            "b": str(negative_coefficient),
            "moment_4": str(fourth_moment),
        }

    artifact = {
        "format": "two-pair-primitive-prefix-obstruction-v1",
        "field": "characteristic zero",
        "family": "G_(d,lambda)=R^(d-4)*F+lambda*Z^d",
        "all_order_statement": {
            "orders_1_through_d": "E_2(G_(d,lambda)^m)=0",
            "order_d_plus_1": (
                "(d+1)*lambda*(d*(d+1)+1)!*d!/(2*d+1)!!"
            ),
        },
        "direct_sparse_replay": replay,
        "positive_phase_triangular_replay": triangular_replay,
        "degree_five_extreme_phase_obstruction": extreme_constants,
        "odd_height_prefix_replay": odd_height_replay,
        "opposite_odd_height_replay": opposite_odd_height_replay,
        "opposite_even_height_phase_at_least_three_replay": (
            opposite_even_height_replay
        ),
        "degree_five_phase_one_obstruction": {
            "moments": phase_one_moments,
            "resultant": phase_one_resultant,
        },
        "uniform_even_height_phase_two_obstruction": {
            "numerator_coefficients_descending": (
                PHASE_TWO_NUMERATOR_COEFFICIENTS
            ),
            "modulus": phase_two_modulus,
            "nonzero_residues": phase_two_residues,
            "direct_sparse_replay": phase_two_replay,
        },
        "written_source": (
            "extended-geometry/TWO_PAIR_PRIMITIVE_PREFIX_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print(
        "PASS primitive prefixes: moments 1..d vanish and moment d+1 "
        "has the claimed nonzero value for 4<=d<=8"
    )
    print(
        "PASS triangular corrections: least phase s is detected exactly "
        "at moment s+1 for 4<=d<=7"
    )
    print(
        "PASS degree-five extreme phases: moments 2,4,6 force a=b=0 "
        "in RF+aZ^5+bW^5"
    )
    print(
        "PASS odd-height prefixes: moments below 2d vanish and moment "
        "2d has the claimed value for 4<=d<=7"
    )
    print(
        "PASS opposite odd heights: mixed moment 2 and both branch "
        "certificates agree for 4<=d<=7"
    )
    print(
        "PASS opposite even heights of phase at least 3: mixed moment 2 "
        "and both linear branch certificates agree for 4<=d<=7"
    )
    print(
        "PASS degree-five phase one: moments 2,3,4 and the nonzero "
        "elimination resultant force a=b=0"
    )
    print(
        "PASS uniform even-height phase two: the exceptional branch "
        "fails at moment 4 and its numerator has no root modulo 47"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
