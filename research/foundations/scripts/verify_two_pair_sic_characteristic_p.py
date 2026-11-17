#!/usr/bin/env python3
"""Exact characteristic-p audit of the denominator-cleared SIC2C4 seed.

Put Fbar=2F, where F is the bidegree-(4,4) two-pair witness from
verify_two_pair_image_mathieu_counterexample.py.  This checker verifies:

* the integral polynomial, quadric, coefficient-rank, and exceptional-prime
  data;
* ordinary contractions of Fbar^m and Z*Fbar^m;
* the Legendre/Lucas/Kummer nonvanishing criterion for the mixed moment;
* the naive Hasse (divided-derivative) contraction formulas; and
* exact Hilbert--Mumford chart calculations at p=2 and p=3.

The bounded loops replay the formulas.  The all-order arguments and the
positive-characteristic Image-kernel statement are written in
extended-geometry/TWO_PAIR_SIC_CHARACTERISTIC_P.md.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import sympy as sp

from verify_two_pair_image_mathieu_counterexample import (
    ZERO,
    Polynomial,
    coefficient_matrix,
    contraction,
    determinant,
    double_factorial_odd,
    hopf_profile_integral,
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
    / "two_pair_sic_characteristic_p.json"
)
MOMENT_CUTOFF = 8
PRIME_CUTOFF = 101
PROPAGATED_DEGREE_CUTOFF = 20


def primes_through(bound: int) -> list[int]:
    result = []
    for candidate in range(2, bound + 1):
        if all(candidate % divisor for divisor in range(2, int(candidate**0.5) + 1)):
            result.append(candidate)
    return result


def matrix_rank_mod_p(matrix: list[list[Fraction]], prime: int) -> int:
    work = [
        [
            (entry.numerator * pow(entry.denominator, -1, prime)) % prime
            for entry in row
        ]
        for row in matrix
    ]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                row
                for row in range(rank, len(work))
                if work[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(entry * inverse) % prime for entry in work[rank]]
        for row in range(len(work)):
            if row == rank:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[rank])
                ]
        rank += 1
    return rank


def balanced_coefficient_matrix(
    polynomial: Polynomial,
    degree: int,
) -> list[list[Fraction]]:
    """Coefficient matrix in xi1^i xi2^(d-i), z1^j z2^(d-j)."""
    matrix = [
        [Fraction(0) for _column in range(degree + 1)]
        for _row in range(degree + 1)
    ]
    for (xi1, xi2, z1, z2), coefficient in polynomial.items():
        assert xi1 + xi2 == z1 + z2 == degree
        matrix[xi1][z1] = coefficient
    return matrix


def p_adic_factorial(order: int, prime: int) -> int:
    result = 0
    quotient = order
    while quotient:
        quotient //= prime
        result += quotient
    return result


def p_adic_integer(value: int, prime: int) -> int:
    assert value > 0
    result = 0
    current = value
    while current % prime == 0:
        current //= prime
        result += 1
    return result


def mixed_valuation(order: int, prime: int) -> int:
    """v_p((4m+2)! m!/(2m+1)!!) for an odd prime."""
    assert prime % 2
    return (
        p_adic_factorial(4 * order + 2, prime)
        + 2 * p_adic_factorial(order, prime)
        - p_adic_factorial(2 * order + 1, prime)
    )


def mixed_floor_valuation(order: int, prime: int) -> int:
    """Equivalent floor-sum form of mixed_valuation."""
    assert prime % 2
    result = 0
    prime_power = prime
    while prime_power <= 4 * order + 2:
        result += (
            (4 * order + 2) // prime_power
            - (2 * order + 1) // prime_power
            + 2 * (order // prime_power)
        )
        prime_power *= prime
    return result


def base_p_digit_sum(value: int, prime: int) -> int:
    result = 0
    current = value
    while current:
        result += current % prime
        current //= prime
    return result


def mixed_digit_valuation(order: int, prime: int) -> int:
    """Legendre digit-sum form of mixed_valuation."""
    numerator = (
        4 * order
        + 1
        - base_p_digit_sum(4 * order + 2, prime)
        - 2 * base_p_digit_sum(order, prime)
        + base_p_digit_sum(2 * order + 1, prime)
    )
    assert numerator % (prime - 1) == 0
    return numerator // (prime - 1)


def propagated_mixed_integer(
    degree: int,
    seed_power: int,
    order: int,
) -> int:
    """Mixed moment for R^k*Fbar^r, where degree=4r+k."""
    assert seed_power >= 1
    assert degree >= 4 * seed_power
    angular_order = seed_power * order
    return (
        4**angular_order
        * factorial(degree * order + 2)
        * factorial(angular_order) ** 2
        // factorial(2 * angular_order + 1)
    )


def propagated_valuation(
    degree: int,
    seed_power: int,
    order: int,
    prime: int,
) -> int:
    """Odd-prime valuation of the propagated cleared mixed moment."""
    angular_order = seed_power * order
    return (
        p_adic_factorial(degree * order + 2, prime)
        + 2 * p_adic_factorial(angular_order, prime)
        - p_adic_factorial(2 * angular_order + 1, prime)
    )


def propagated_digit_valuation(
    degree: int,
    seed_power: int,
    order: int,
    prime: int,
) -> int:
    angular_order = seed_power * order
    numerator = (
        degree * order
        + 1
        - base_p_digit_sum(degree * order + 2, prime)
        - 2 * base_p_digit_sum(angular_order, prime)
        + base_p_digit_sum(2 * angular_order + 1, prime)
    )
    assert numerator % (prime - 1) == 0
    return numerator // (prime - 1)


def binomial_or_zero(top: int, bottom: int) -> int:
    if bottom < 0 or bottom > top:
        return 0
    return comb(top, bottom)


def first_lift_convolution_matrix(radial_power: int) -> list[list[Fraction]]:
    """Coefficient matrix of R^k*Fbar from its four diagonal symbols."""
    degree = 4 + radial_power
    matrix = [
        [Fraction(0) for _column in range(degree + 1)]
        for _row in range(degree + 1)
    ]
    for row in range(degree + 1):
        if row + 1 <= degree:
            matrix[row][row + 1] = 4 * binomial_or_zero(
                radial_power + 3,
                row,
            )
        matrix[row][row] = -2 * (
            binomial_or_zero(radial_power + 2, row - 2)
            - 4 * binomial_or_zero(radial_power + 2, row - 1)
            + binomial_or_zero(radial_power + 2, row)
        )
        if row >= 1:
            matrix[row][row - 1] = -3 * (
                binomial_or_zero(radial_power + 1, row - 3)
                - 2 * binomial_or_zero(radial_power + 1, row - 2)
                + binomial_or_zero(radial_power + 1, row - 1)
            )
        if row >= 2:
            matrix[row][row - 2] = -(
                binomial_or_zero(radial_power, row - 4)
                - 2 * binomial_or_zero(radial_power, row - 3)
                + binomial_or_zero(radial_power, row - 2)
            )
    return matrix


def lower_hessenberg_determinant(
    matrix: list[list[Fraction]],
) -> Fraction:
    """Continuant for a matrix with one upper and two lower diagonals."""
    leading = [Fraction(1)]
    for size in range(1, len(matrix) + 1):
        row = size - 1
        value = matrix[row][row] * leading[size - 1]
        if size >= 2:
            value -= (
                matrix[row][row - 1]
                * matrix[row - 1][row]
                * leading[size - 2]
            )
        if size >= 3:
            value += (
                matrix[row][row - 2]
                * matrix[row - 2][row - 1]
                * matrix[row - 1][row]
                * leading[size - 3]
            )
        leading.append(value)
    return leading[-1]


def lucas_binomial_mod(top: int, bottom: int, prime: int) -> int:
    result = 1
    upper = top
    lower = bottom
    while upper or lower:
        upper_digit = upper % prime
        lower_digit = lower % prime
        if lower_digit > upper_digit:
            return 0
        result = result * comb(upper_digit, lower_digit) % prime
        upper //= prime
        lower //= prime
    return result


def audit_hasse_image_intertwining() -> None:
    """Check the universal divided-raising/Hasse-derivative identity."""
    for coordinate_order in range(13):
        for dual_order in range(13):
            for hasse_order in range(1, 13):
                left = (
                    binomial_or_zero(coordinate_order, hasse_order)
                    * binomial_or_zero(
                        coordinate_order - hasse_order,
                        dual_order,
                    )
                    if coordinate_order >= hasse_order
                    else 0
                )
                right = binomial_or_zero(
                    dual_order + hasse_order,
                    hasse_order,
                ) * binomial_or_zero(
                    coordinate_order,
                    dual_order + hasse_order,
                )
                assert left == right


def audit_p_typical_generation() -> None:
    """Check the no-carry units used to reduce to p-power operators."""
    for prime in primes_through(43):
        for order in range(1, 256):
            partial = 0
            quotient = order
            digit_place = 1
            while quotient:
                digit = quotient % prime
                for _copy in range(digit):
                    coefficient = comb(
                        partial + digit_place,
                        digit_place,
                    )
                    assert coefficient % prime != 0
                    partial += digit_place
                quotient //= prime
                digit_place *= prime
            assert partial == order


def audit_frobenius_sic_collapse() -> None:
    """Audit E((xi^a z^b)^p)=0 for a!=0 and z^(pb) for a=0."""

    def falling_mod(top: int, order: int, prime: int) -> int:
        result = 1
        for offset in range(order):
            result = result * (top - offset) % prime
        return result

    for prime in primes_through(19):
        assert factorial(prime - 1) % prime != 0
        assert factorial(prime) % prime == 0
        for xi1 in range(5):
            for xi2 in range(5):
                for z1 in range(5):
                    for z2 in range(5):
                        coefficient = (
                            falling_mod(prime * z1, prime * xi1, prime)
                            * falling_mod(prime * z2, prime * xi2, prime)
                        ) % prime
                        if xi1 <= z1 and xi2 <= z2:
                            expected = 1 if xi1 == xi2 == 0 else 0
                            assert coefficient == expected
                        else:
                            assert coefficient == 0


def profile_two_integral(order: int) -> Fraction:
    denominator = 1
    for index in range(order + 1):
        denominator *= 4 * index + 1
    return Fraction(4**order * factorial(order), denominator)


def profile_mixed_integer(height: int, order: int) -> int:
    value = (
        Fraction(factorial(4 * height * order + 2))
        * hopf_profile_integral(height, order)
    )
    assert value.denominator == 1
    return value.numerator


def profile_two_mixed_integer(order: int) -> int:
    value = Fraction(factorial(8 * order + 2)) * profile_two_integral(order)
    assert value.denominator == 1
    return value.numerator


def audit_non_power_profiles() -> None:
    """Audit the Phi_2 closed form and general necessary prime cutoff."""
    for order in range(1, 81):
        assert hopf_profile_integral(2, order) == profile_two_integral(order)

    for height in range(1, 9):
        for order in range(1, 13):
            mixed = profile_mixed_integer(height, order)
            bound = 4 * height * order + 2
            for prime in primes_through(bound):
                assert mixed % prime == 0

    for prime in primes_through(PRIME_CUTOFF):
        for order in range(1, 2 * prime + 1):
            mixed = profile_two_mixed_integer(order)
            assert (mixed % prime != 0) == (8 * order + 2 < prime)

    # General profiles can have additional numerator-prime holes above the
    # factorial cutoff; these are the first small exact examples.
    assert profile_mixed_integer(6, 1) % 47 == 0
    assert 47 > 4 * 6 + 2
    assert profile_mixed_integer(4, 5) % 89 == 0
    assert 89 > 4 * 4 * 5 + 2


def audit_prime_power_survival() -> None:
    """Audit radial monotonicity and non-radial re-entry modulo p^a."""
    radial_values = {
        order: propagated_mixed_integer(4, 1, order)
        for order in range(1, 201)
    }
    for order in range(1, 200):
        quotient = (
            16
            * (4 * order + 3)
            * (4 * order + 5)
            * (order + 1) ** 2
        )
        assert radial_values[order + 1] == radial_values[order] * quotient

    # Every radial-power sequence is a subsequence A_(r*m).
    for seed_power in range(1, 6):
        for prime in primes_through(43):
            valuations = [
                propagated_valuation(
                    4 * seed_power,
                    seed_power,
                    order,
                    prime,
                )
                for order in range(1, 101)
            ]
            assert valuations == sorted(valuations)

    # The general consecutive-order recurrence is local in linear factors.
    for degree in range(4, 17):
        for seed_power in range(1, degree // 4 + 1):
            for prime in primes_through(31):
                for order in range(1, 41):
                    delta = 2 * seed_power if prime == 2 else 0
                    delta += sum(
                        p_adic_integer(degree * order + offset, prime)
                        for offset in range(3, degree + 3)
                    )
                    delta += 2 * sum(
                        p_adic_integer(seed_power * order + offset, prime)
                        for offset in range(1, seed_power + 1)
                    )
                    delta -= sum(
                        p_adic_integer(2 * seed_power * order + offset, prime)
                        for offset in range(2, 2 * seed_power + 2)
                    )
                    next_valuation = propagated_valuation(
                        degree,
                        seed_power,
                        order + 1,
                        prime,
                    )
                    current_valuation = propagated_valuation(
                        degree,
                        seed_power,
                        order,
                        prime,
                    )
                    if prime == 2:
                        next_valuation += 2 * seed_power * (order + 1)
                        current_valuation += 2 * seed_power * order
                    assert delta == (
                        next_valuation - current_valuation
                    )

    # The first non-radial lift already has prime-power re-entry.
    fourth = propagated_mixed_integer(5, 1, 4)
    fifth = propagated_mixed_integer(5, 1, 5)
    assert propagated_valuation(5, 1, 4, 11) == 2
    assert propagated_valuation(5, 1, 5, 11) == 1
    assert fourth % 121 == 0
    assert fifth % 121 == 22


def audit_first_lift_tensor_ranks(
    cleared: Polynomial,
    radial: Polynomial,
) -> dict[str, dict[str, object]]:
    """Compute the complete coefficient-rank phase diagram for R^k*Fbar."""
    expected = {
        0: (1536, {2: 4, 3: 4}),
        1: (-17408, {2: 4, 17: 5}),
        2: (376832, {2: 4, 23: 6}),
        3: (-15298560, {2: 4, 3: 7, 5: 7, 83: 7}),
        4: (362086400, {2: 8, 5: 7, 13: 8, 17: 8}),
    }
    result = {}
    recurrence_determinants = {}
    for radial_power in range(21):
        degree = 4 + radial_power
        polynomial = multiply(power(radial, radial_power), cleared)
        expanded_matrix = balanced_coefficient_matrix(polynomial, degree)
        convolution_matrix = first_lift_convolution_matrix(radial_power)
        assert expanded_matrix == convolution_matrix
        recurrence_determinant = lower_hessenberg_determinant(
            convolution_matrix
        )
        assert recurrence_determinant == determinant(expanded_matrix)
        recurrence_determinants[str(radial_power)] = int(
            recurrence_determinant
        )

    # In characteristic two, put n=k+2.  Lucas parity reduces the matrix
    # to a bidiagonal one of rank 2^(1+s_2(floor(n/2))).
    for radial_power in range(129):
        expected_rank_two = 2 ** (
            1 + ((radial_power + 2) // 2).bit_count()
        )
        assert matrix_rank_mod_p(
            first_lift_convolution_matrix(radial_power),
            2,
        ) == expected_rank_two

    for radial_power, (expected_determinant, exceptional_ranks) in (
        expected.items()
    ):
        degree = 4 + radial_power
        polynomial = multiply(power(radial, radial_power), cleared)
        matrix = balanced_coefficient_matrix(polynomial, degree)
        actual_determinant = determinant(matrix)
        assert actual_determinant == expected_determinant
        for prime in primes_through(PRIME_CUTOFF):
            rank = matrix_rank_mod_p(matrix, prime)
            expected_rank = exceptional_ranks.get(prime, degree + 1)
            assert rank == expected_rank
        result[str(radial_power)] = {
            "degree": degree,
            "determinant": expected_determinant,
            "exceptional_ranks": {
                str(prime): rank
                for prime, rank in exceptional_ranks.items()
            },
        }
    result["recurrence_audit"] = {
        "radial_powers": [0, 20],
        "determinants": recurrence_determinants,
    }
    return result


def audit_hasse_sic_counterexample() -> None:
    """Audit f=xi*z^p, g=z against the Hasse Image in one pair."""
    for prime in primes_through(43):
        for order in range(1, 501):
            assert lucas_binomial_mod(
                prime * order,
                order,
                prime,
            ) == 0
        geometric_sum = 0
        prime_power = 1
        for _exponent in range(1, 11):
            geometric_sum += prime_power
            prime_power *= prime
            assert lucas_binomial_mod(
                prime * geometric_sum + 1,
                geometric_sum,
                prime,
            ) == 1


def hasse_contraction_scalar(polynomial: Polynomial) -> Fraction:
    """Balanced Hasse contraction: sum of diagonal coefficients."""
    result = Fraction(0)
    for (xi1, xi2, z1, z2), coefficient in polynomial.items():
        if xi1 == z1 and xi2 == z2:
            result += coefficient
    return result


def sympy_seed() -> tuple[
    sp.Expr,
    tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
]:
    xi1, xi2, z1, z2 = sp.symbols("xi1 xi2 z1 z2")
    radial = xi1 * z1 + xi2 * z2
    phase = xi1 * z2
    opposite = 2 * xi2 * z1
    torus = xi1 * z1 - xi2 * z2
    cleared = sp.expand(
        (radial + phase)
        * (
            2 * radial**2 * opposite
            - (2 * radial + phase) * torus**2
        )
    )
    return cleared, (xi1, xi2, z1, z2)


def transformed_nonpositive_coefficients() -> tuple[
    list[sp.Expr],
    tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
]:
    """Equations for a conjugate supported only in weights j-i>0."""
    cleared, (xi1, xi2, z1, z2) = sympy_seed()
    a, b, c, d = sp.symbols("a b c d")
    transformed = sp.expand(
        cleared.subs(
            {
                xi1: d * xi1 - c * xi2,
                xi2: -b * xi1 + a * xi2,
                z1: a * z1 + b * z2,
                z2: c * z1 + d * z2,
            },
            simultaneous=True,
        )
    )
    polynomial = sp.Poly(transformed, xi1, xi2, z1, z2)
    equations = []
    for (xi1_order, _xi2_order, z1_order, _z2_order), coefficient in (
        polynomial.terms()
    ):
        if z1_order <= xi1_order:
            equations.append(coefficient)
    return equations, (a, b, c, d)


def one_sided_after_matrix(
    prime: int,
    matrix: tuple[int, int, int, int],
) -> bool:
    equations, (a, b, c, d) = transformed_nonpositive_coefficients()
    substitution = dict(zip((a, b, c, d), matrix))
    determinant_equation = a * d - b * c - 1
    return all(
        int(equation.subs(substitution)) % prime == 0
        for equation in equations + [determinant_equation]
    )


def main() -> None:
    seed, multiplier, generators = witness()
    cleared = scale(2, seed)
    radial = generators["R"]
    phase = generators["Z"]
    opposite = generators["W"]
    torus = generators["T"]

    # Integral geometry and the two exceptional primes.
    assert all(coefficient.denominator == 1 for coefficient in cleared.values())
    quadric_defect = {}
    for polynomial, coefficient in (
        (power(torus, 2), 1),
        (power(radial, 2), -1),
        (multiply(phase, opposite), 2),
    ):
        for exponent, value in polynomial.items():
            quadric_defect[exponent] = (
                quadric_defect.get(exponent, Fraction(0))
                + coefficient * value
            )
    assert not {exponent: value for exponent, value in quadric_defect.items() if value}

    cleared_matrix = coefficient_matrix(cleared)
    assert determinant(cleared_matrix) == 1536
    assert matrix_rank_mod_p(cleared_matrix, 2) == 4
    assert matrix_rank_mod_p(cleared_matrix, 3) == 4
    for prime in primes_through(PRIME_CUTOFF):
        if prime not in (2, 3):
            assert matrix_rank_mod_p(cleared_matrix, prime) == 5

    factorial_diagonal = sp.diag(24, 6, 4, 6, 24)
    sympy_matrix = sp.Matrix(
        [[sp.Integer(entry) for entry in row] for row in cleared_matrix]
    )
    invariant_i2 = sp.trace((factorial_diagonal * sympy_matrix) ** 2)
    assert invariant_i2 == 4608

    # Ordinary contractions and the integer mixed formula.
    cleared_power = monomial(ZERO)
    replay = {}
    for order in range(1, MOMENT_CUTOFF + 1):
        cleared_power = multiply(cleared_power, cleared)
        assert contraction(cleared_power) == {}
        base_moment = (
            factorial(4 * order + 2)
            * factorial(order)
            // double_factorial_odd(2 * order + 1)
        )
        mixed_moment = 2**order * base_moment
        assert mixed_moment == (
            4**order
            * comb(4 * order + 2, 2 * order + 1)
            * factorial(2 * order + 1)
            * factorial(order) ** 2
        )
        assert contraction(multiply(multiplier, cleared_power)) == {
            ZERO: Fraction(mixed_moment)
        }

        hasse_pure = hasse_contraction_scalar(cleared_power)
        hasse_mixed = hasse_contraction_scalar(
            multiply(multiplier, cleared_power)
        )
        assert hasse_pure == 16**order
        assert hasse_mixed == 2 * order * 16**order
        replay[str(order)] = {
            "ordinary_mixed": mixed_moment,
            "hasse_pure": int(hasse_pure),
            "hasse_mixed": int(hasse_mixed),
        }

    # Legendre, Lucas, and Kummer forms agree.  For odd p, the moment is
    # nonzero exactly when the doubling (2m+1)+(2m+1) has no base-p carry
    # and the factorial factors have not yet met p, equivalently 4m+2<p.
    phase_diagram = {}
    for prime in primes_through(PRIME_CUTOFF):
        nonzero_orders = []
        for order in range(1, 2 * prime + 1):
            base_moment = (
                factorial(4 * order + 2)
                * factorial(order)
                // double_factorial_odd(2 * order + 1)
            )
            cleared_moment = 2**order * base_moment
            if prime == 2:
                assert cleared_moment % prime == 0
                continue
            valuation = mixed_valuation(order, prime)
            assert valuation == mixed_floor_valuation(order, prime)
            assert valuation == mixed_digit_valuation(order, prime)
            assert (cleared_moment % prime != 0) == (valuation == 0)
            assert (valuation == 0) == (4 * order + 2 < prime)
            if valuation == 0:
                nonzero_orders.append(order)
        expected = list(range(1, (prime - 3) // 4 + 1))
        assert nonzero_orders == expected
        phase_diagram[str(prime)] = nonzero_orders

    # The same arithmetic becomes simpler for every propagated witness
    # R^k*Fbar^r of balanced degree d=4r+k.  After cancellation, the mixed
    # integer is 4^(rm)*(rm)!^2 times the consecutive product from
    # 2rm+2 through dm+2.  Since d>=4r, it avoids p exactly when dm+2<p.
    for degree in range(4, PROPAGATED_DEGREE_CUTOFF + 1):
        for seed_power in range(1, degree // 4 + 1):
            for prime in primes_through(PRIME_CUTOFF):
                nonzero_orders = []
                for order in range(1, 2 * prime + 1):
                    mixed_integer = propagated_mixed_integer(
                        degree,
                        seed_power,
                        order,
                    )
                    if prime == 2:
                        assert mixed_integer % prime == 0
                        continue
                    valuation = propagated_valuation(
                        degree,
                        seed_power,
                        order,
                        prime,
                    )
                    assert valuation == propagated_digit_valuation(
                        degree,
                        seed_power,
                        order,
                        prime,
                    )
                    assert (mixed_integer % prime != 0) == (valuation == 0)
                    assert (valuation == 0) == (
                        degree * order + 2 < prime
                    )
                    if valuation == 0:
                        nonzero_orders.append(order)
                expected = list(range(1, (prime - 3) // degree + 1))
                assert nonzero_orders == expected

    # The full Hasse-compatible Image system uses every positive
    # divided-power order.  The binomial identity checked here is the
    # termwise intertwining that makes its common image equal ker(H).
    audit_hasse_image_intertwining()
    audit_p_typical_generation()
    audit_frobenius_sic_collapse()
    audit_non_power_profiles()
    audit_prime_power_survival()
    first_lift_tensor_ranks = audit_first_lift_tensor_ranks(cleared, radial)
    audit_hasse_sic_counterexample()

    # The characteristic-two reduction is explicitly one-sided.
    assert one_sided_after_matrix(2, (0, 1, 1, 0))

    # At p=3 the exact Hilbert--Mumford chart ideal is the unit ideal, so
    # no algebraic-closure conjugate has support only in weights j-i>0.
    equations, (a, b, c, d) = transformed_nonpositive_coefficients()
    groebner_three = sp.groebner(
        equations + [a * d - b * c - 1],
        a,
        b,
        c,
        d,
        modulus=3,
        order="grevlex",
    )
    assert groebner_three.contains(sp.Integer(1))
    assert len(groebner_three.polys) == 1

    sympy_cleared, (xi1, xi2, z1, z2) = sympy_seed()
    radial_expr = xi1 * z1 + xi2 * z2
    phase_expr = xi1 * z2
    opposite_expr = 2 * xi2 * z1
    assert sp.Poly(
        sympy_cleared - phase_expr * radial_expr**2 * (radial_expr + phase_expr),
        xi1,
        xi2,
        z1,
        z2,
        modulus=2,
    ).is_zero
    assert sp.Poly(
        sympy_cleared
        - (radial_expr + phase_expr)
        * (
            radial_expr**3
            - radial_expr**2 * opposite_expr
            - radial_expr**2 * phase_expr
            + radial_expr * phase_expr * opposite_expr
            - phase_expr**2 * opposite_expr
        ),
        xi1,
        xi2,
        z1,
        z2,
        modulus=3,
    ).is_zero

    artifact = {
        "format": "two-pair-sic-characteristic-p-v1",
        "integral_seed": "Fbar=2F",
        "coefficient_matrix_determinant": 1536,
        "full_matrix_rank_primes": "all primes except 2 and 3",
        "quadric_nondegenerate_primes": "all odd primes",
        "ordinary_contraction": {
            "E_2(Fbar^m)": "0",
            "E_2(Z*Fbar^m)": "2^m*(4m+2)!*m!/(2m+1)!!",
            "mod_p_nonzero": "p odd and 4m+2<p",
            "fixed_p_pattern": "m=1,...,floor((p-3)/4)",
        },
        "valuation_formulas_for_odd_p": {
            "Legendre": (
                "v_p((4m+2)!)+2*v_p(m!)-v_p((2m+1)!)"
            ),
            "floor_sum": (
                "sum_j(floor((4m+2)/p^j)-floor((2m+1)/p^j)"
                "+2*floor(m/p^j))"
            ),
            "digit_sum": (
                "(4m+1-s_p(4m+2)-2s_p(m)+s_p(2m+1))/(p-1)"
            ),
            "Kummer_factorization": (
                "2^m*binom(4m+2,2m+1)*(2m+1)!*(m!)^2"
            ),
        },
        "hasse_contraction": {
            "H_2(Fbar^m)": "16^m",
            "H_2(Z*Fbar^m)": "2m*16^m",
            "conclusion": "naive Hasse replacement destroys pure vanishing for odd p",
        },
        "propagated_family": {
            "form": "Gbar_(r,k)=R^k*Fbar^r, d=4r+k",
            "mixed_moment": (
                "4^(r*m)*(d*m+2)!*((r*m)!)^2/(2*r*m+1)!"
            ),
            "mod_p_nonzero": "p odd and d*m+2<p",
            "fixed_p_pattern": "m=1,...,floor((p-3)/d)",
            "p_adic_valuation": (
                "(d*m+1-s_p(d*m+2)-2s_p(r*m)"
                "+s_p(2*r*m+1))/(p-1)"
            ),
            "replay": {
                "degrees": [4, PROPAGATED_DEGREE_CUTOFF],
                "all_seed_powers_r_with_4r<=d": True,
                "primes_through": PRIME_CUTOFF,
                "orders_through": "2*p",
            },
        },
        "complete_hasse_image_system": {
            "raising": (
                "X^[a](xi^c)=binom(c+a,a)*xi^(c+a)"
            ),
            "operators": "D_[a]=partial_z^[a]-X_xi^[a] for every a>0",
            "p_typical_generators": (
                "D_(i,p^j) for 1<=i<=r and j>=0 suffice"
            ),
            "kernel": "sum_a image(D_[a])=kernel(H)",
            "intertwining_identity": (
                "binom(b,a)*binom(b-a,c)"
                "=binom(c+a,a)*binom(b,c+a)"
            ),
            "audit_order_cutoff": 12,
        },
        "full_ordinary_positive_characteristic_sic": {
            "frobenius_identity": "E_(r,p)(f^p)=f(0,z)^p",
            "single_moment_premise": "E_(r,p)(f^p)=0",
            "conclusion_cutoff": "E_(r,p)(g*f^m)=0 for every m>=p",
            "consequence": "ordinary SIC(r) holds for every r and p>0",
            "sharpness": (
                "f=xi1, g=z1^(p-1): all pure moments vanish but "
                "E(g*f^(p-1))=(p-1)! is nonzero"
            ),
        },
        "non_power_hopf_profiles": {
            "necessary_nonvanishing_condition": "p>4*h*m+2",
            "height_two_integral": (
                "C_(2,m)=4^m*m!/product_(j=0)^m(4j+1)"
            ),
            "height_two_nonzero": "p odd and 8*m+2<p",
            "additional_hole_examples": [
                {"height": 6, "order": 1, "prime": 47},
                {"height": 4, "order": 5, "prime": 89},
            ],
            "general_conclusion": (
                "above the factorial cutoff, numerator primes of C_(h,m) "
                "can create holes"
            ),
        },
        "prime_power_survival": {
            "radial_family": "d=4*r, so M_(r,m)=A_(r*m)",
            "radial_ratio": (
                "A_(s+1)/A_s="
                "16*(4*s+3)*(4*s+5)*(s+1)^2"
            ),
            "radial_consequence": (
                "v_p(M_(r,m)) is nondecreasing in m for every prime"
            ),
            "general_increment": (
                "v_p(M_(m+1))-v_p(M_m) is the signed valuation sum "
                "of the three consecutive-factor blocks in M_(m+1)/M_m"
            ),
            "non_radial_reentry": {
                "family": "R*Fbar (d=5,r=1)",
                "prime_power": "11^2",
                "order_4_valuation": 2,
                "order_5_valuation": 1,
                "order_4_residue": 0,
                "order_5_residue": 22,
            },
            "conclusion": (
                "outside d=4*r, prime-power survival need not be "
                "an initial interval"
            ),
        },
        "first_lift_tensor_ranks": {
            "family": "R^k*Fbar for k>=0",
            "data": first_lift_tensor_ranks,
            "diagonal_symbols": {
                "upper": "4*(1+x)^(k+3)",
                "diagonal": "-2*(1+x)^(k+2)*(x^2-4*x+1)",
                "lower_1": "-3*x*(x-1)^2*(1+x)^(k+1)",
                "lower_2": "-x^2*(x-1)^2*(1+x)^k",
            },
            "determinant_method": (
                "three-term lower-Hessenberg continuant, audited for "
                "0<=k<=20"
            ),
            "characteristic_two_rank": (
                "2^(1+s_2(floor((k+2)/2))), audited for 0<=k<=128"
            ),
            "conclusion": (
                "higher degree introduces exceptional odd rank primes, "
                "but semistability still propagates"
            ),
        },
        "hasse_sic_counterexample": {
            "pair_dimension": 1,
            "f": "xi*z^p",
            "g": "z",
            "pure": "H(f^m)=binom(p*m,m)*z^((p-1)*m)=0 for all m>0",
            "mixed_orders": "m_e=(p^e-1)/(p-1)",
            "mixed": (
                "H(g*f^m_e)=z^((p-1)*m_e+1), since "
                "binom(p*m_e+1,m_e)=1 mod p"
            ),
            "consequence": "the p-typical Hasse Image is not Mathieu",
        },
        "exceptional_geometry": {
            "p=2": {
                "matrix_rank": 4,
                "reduction": "Z*R^2*(R+Z)",
                "displayed_seed_in_nullcone": True,
            },
            "p=3": {
                "matrix_rank": 4,
                "reduction": (
                    "(R+Z)*(R^3-R^2*W-R^2*Z+R*Z*W-Z^2*W)"
                ),
                "hilbert_mumford_chart_groebner_basis": ["1"],
                "displayed_seed_in_nullcone": False,
            },
            "p>=5": {
                "matrix_rank": 5,
                "I_2(Fbar)": 4608,
                "displayed_seed_in_nullcone": False,
            },
        },
        "replay_cutoff": MOMENT_CUTOFF,
        "prime_cutoff": PRIME_CUTOFF,
        "replay": replay,
        "phase_diagram": phase_diagram,
        "written_proof": "extended-geometry/TWO_PAIR_SIC_CHARACTERISTIC_P.md",
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS SIC2C4 mod p: integral seed and exceptional primes 2,3")
    print(
        "PASS SIC2C4 mod p: mixed moment is nonzero exactly when 4m+2<p"
    )
    print("PASS SIC2C4 mod p: Legendre/Lucas/Kummer criteria agree")
    print("PASS SIC2C4 mod p: all propagated degrees have cutoff d*m+2<p")
    print("PASS SIC2C4 mod p: naive Hasse moments are 16^m and 2m*16^m")
    print("PASS SIC2C4 mod p: complete Hasse Image intertwining")
    print("PASS characteristic-p SIC: E(f^p)=f(0,z)^p and full collapse")
    print("PASS Hopf profiles mod p: Phi_2 exact; higher profiles have holes")
    print("PASS prime powers: radial monotonicity; non-radial re-entry mod 11^2")
    print("PASS tensor ranks: exact first-lift recurrence through k=20")
    print("PASS Hasse-SIC: one-pair counterexample in every characteristic")
    print("PASS SIC2C4 mod p: p=2 nullcone, p=3 semistable")
    print(f"PASS SIC2C4 mod p: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
