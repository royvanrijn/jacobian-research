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


def p_adic_factorial(order: int, prime: int) -> int:
    result = 0
    quotient = order
    while quotient:
        quotient //= prime
        result += quotient
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
    print("PASS SIC2C4 mod p: naive Hasse moments are 16^m and 2m*16^m")
    print("PASS SIC2C4 mod p: p=2 nullcone, p=3 semistable")
    print(f"PASS SIC2C4 mod p: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
