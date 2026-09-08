#!/usr/bin/env python3
"""Certify a geometric rank upper bound for Nagao's section-7 K3 family.

The elliptic K3 surface attached to the roots ``(346,260,255,146,55,0)``
has twenty finite ``I1`` fibers and one split ``I4`` fiber at infinity.  The
trivial lattice and the twelve independently certified generic sections give
seventeen divisor classes defined over ``Q``.

This verifier counts the good reduction at ``p=29`` over both ``F_p`` and
``F_{p^2}``.  Removing the seventeen known Frobenius eigenvalues leaves a
degree-five reciprocal factor.  Its only eigenvalue of the form ``p*zeta``
with ``zeta`` a root of unity is ``-p``.  Specialization of Neron--Severi
groups and the cycle-class eigenvalue condition therefore give

    rho(S over Qbar) <= 18,

and Shioda--Tate gives geometric generic Mordell--Weil rank at most thirteen.

More strongly, the residual factor has no eigenvalue ``+p``.  Every divisor
defined over ``Q`` reduces to a class in the ``+p`` eigenspace, so the
seventeen displayed rational divisor classes exhaust ``NS(S/Q)``.  Hence the
arithmetic generic Mordell--Weil rank over ``Q(T)`` is exactly twelve.  The
geometric generic rank over ``Qbar(T)`` remains in the interval 12 through 13.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

import sympy as sp

from verify_nagao_section7_linear_sections import _quartic


GOOD_PRIME = 29
KNOWN_RATIONAL_NS_RANK = 17
TRIVIAL_LATTICE_RANK = 5
EXPECTED_SURFACE_POINT_COUNTS = {1: 1212, 2: 723600}
EXPECTED_H2_TRACES = {1: 370, 2: 16318}
EXPECTED_RESIDUAL_TRACES = {1: -123, 2: 2021}
ROOT_OF_UNITY_ORDERS_OF_DEGREE_AT_MOST_FOUR = (1, 2, 3, 4, 5, 6, 8, 10, 12)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_nagao_section7_picard_bound.py"
)


def section7_short_polynomials() -> tuple[sp.Poly, sp.Poly, sp.Poly]:
    """Return the integral short Weierstrass ``A``, ``B`` and discriminant."""

    parameter, x = sp.symbols("T X")
    quartic = sp.Poly(_quartic(parameter, x), x)
    e, d, c, b, a = (
        quartic.coeff_monomial(x**degree) for degree in range(5)
    )
    invariant_i = sp.expand(12 * a * e - 3 * b * d + c**2)
    invariant_j = sp.expand(
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    coefficient_a = sp.Poly(sp.expand(-27 * invariant_i), parameter, domain=sp.ZZ)
    coefficient_b = sp.Poly(sp.expand(-27 * invariant_j), parameter, domain=sp.ZZ)
    discriminant = sp.Poly(
        sp.expand(-16 * (4 * coefficient_a.as_expr() ** 3 + 27 * coefficient_b.as_expr() ** 2)),
        parameter,
        domain=sp.ZZ,
    )
    if (coefficient_a.degree(), coefficient_b.degree(), discriminant.degree()) != (
        8,
        12,
        20,
    ):
        raise AssertionError("the section-7 short model degrees changed")
    return coefficient_a, coefficient_b, discriminant


def verify_good_reduction(prime: int = GOOD_PRIME) -> dict[str, Any]:
    """Check that the fiber configuration has good reduction at ``prime``."""

    if prime != GOOD_PRIME:
        raise ValueError("this pinned certificate is for p=29")
    coefficient_a, coefficient_b, discriminant = section7_short_polynomials()
    reduced_a = sp.Poly(coefficient_a, modulus=prime)
    reduced_b = sp.Poly(coefficient_b, modulus=prime)
    reduced_delta = sp.Poly(discriminant, modulus=prime)
    if (reduced_a.degree(), reduced_b.degree(), reduced_delta.degree()) != (8, 12, 20):
        raise AssertionError("a leading coefficient vanished modulo p")
    if sp.gcd(reduced_delta, reduced_delta.diff()).degree() != 0:
        raise AssertionError("the finite discriminant is not squarefree modulo p")
    if sp.gcd(reduced_delta, reduced_a).degree() != 0:
        raise AssertionError("a finite discriminant root is not type I1 modulo p")

    # The infinity cubic is (X-108)^2(X+216).  Its two roots remain distinct,
    # and its tangent square 18^2 remains nonzero, so the I4 fiber stays split.
    if (324 % prime) == 0 or ((108 - (-216)) % prime) == 0:
        raise AssertionError("the split I4 fiber degenerated modulo p")
    return {
        "prime": prime,
        "short_coefficient_degrees_mod_p": [8, 12],
        "finite_discriminant_degree_mod_p": 20,
        "finite_discriminant_squarefree_mod_p": True,
        "finite_fibers": "20 geometrically simple I1 fibers",
        "infinity_fiber": "split I4",
        "good_reduction_of_resolved_k3": True,
    }


def _ascending_coefficients_mod(polynomial: sp.Poly, prime: int) -> tuple[int, ...]:
    parameter = polynomial.gens[0]
    return tuple(
        int(polynomial.coeff_monomial(parameter**degree)) % prime
        for degree in range(polynomial.degree() + 1)
    )


def _evaluate_prime(coefficients: tuple[int, ...], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def _least_nonsquare(prime: int) -> int:
    squares = {value * value % prime for value in range(prime)}
    return next(value for value in range(2, prime) if value not in squares)


def finite_affine_character_sum(prime: int, extension_degree: int) -> int:
    """Return ``sum_{t,x} chi(x^3+A(t)x+B(t))`` exactly.

    For degree two, write ``F_{p^2}=F_p[w]/(w^2-d)``.  The quadratic
    character of a nonzero element equals the Legendre symbol of its norm,
    which keeps the count dependency-free and fully integral.
    """

    coefficient_a, coefficient_b, _ = section7_short_polynomials()
    coefficients_a = _ascending_coefficients_mod(coefficient_a, prime)
    coefficients_b = _ascending_coefficients_mod(coefficient_b, prime)
    nonzero_squares = {value * value % prime for value in range(1, prime)}

    if extension_degree == 1:
        total = 0
        for parameter in range(prime):
            value_a = _evaluate_prime(coefficients_a, parameter, prime)
            value_b = _evaluate_prime(coefficients_b, parameter, prime)
            for x in range(prime):
                right = (x * x * x + value_a * x + value_b) % prime
                if right:
                    total += 1 if right in nonzero_squares else -1
        return total
    if extension_degree != 2:
        raise ValueError("only F_p and F_{p^2} are used by this certificate")

    nonsquare = _least_nonsquare(prime)
    elements = tuple((real, imag) for real in range(prime) for imag in range(prime))

    def multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return (
            (left[0] * right[0] + nonsquare * left[1] * right[1]) % prime,
            (left[0] * right[1] + left[1] * right[0]) % prime,
        )

    def evaluate(
        coefficients: tuple[int, ...], value: tuple[int, int]
    ) -> tuple[int, int]:
        result = (0, 0)
        for coefficient in reversed(coefficients):
            product = multiply(result, value)
            result = ((product[0] + coefficient) % prime, product[1])
        return result

    x_data = []
    for x in elements:
        square = multiply(x, x)
        cube = multiply(square, x)
        x_data.append((x[0], x[1], cube[0], cube[1]))

    total = 0
    for parameter in elements:
        value_a = evaluate(coefficients_a, parameter)
        value_b = evaluate(coefficients_b, parameter)
        for x_real, x_imag, cube_real, cube_imag in x_data:
            product_real = (
                value_a[0] * x_real + nonsquare * value_a[1] * x_imag
            ) % prime
            product_imag = (
                value_a[0] * x_imag + value_a[1] * x_real
            ) % prime
            right_real = (cube_real + product_real + value_b[0]) % prime
            right_imag = (cube_imag + product_imag + value_b[1]) % prime
            if right_real == 0 and right_imag == 0:
                continue
            norm = (right_real * right_real - nonsquare * right_imag * right_imag) % prime
            total += 1 if norm in nonzero_squares else -1
    return total


def surface_point_count(prime: int, extension_degree: int) -> dict[str, int]:
    """Count the resolved K3 surface over ``F_{p^extension_degree}``."""

    q = prime**extension_degree
    character_sum = finite_affine_character_sum(prime, extension_degree)
    # The q finite fibers contribute q(q+1)+character_sum.  The split I4
    # fiber at infinity contributes 4q points.
    point_count = q * q + 5 * q + character_sum
    h2_trace = point_count - 1 - q * q
    return {
        "field_size": q,
        "finite_affine_character_sum": character_sum,
        "surface_point_count": point_count,
        "h2_frobenius_trace": h2_trace,
    }


def reconstruct_residual_frobenius(
    prime: int, trace_degree_one: int, trace_degree_two: int
) -> dict[str, Any]:
    """Reconstruct and analyze the reciprocal residual degree-five factor."""

    residual_trace_one = trace_degree_one - KNOWN_RATIONAL_NS_RANK * prime
    residual_trace_two = trace_degree_two - KNOWN_RATIONAL_NS_RANK * prime**2
    valid_signs = []
    sign_records = []
    pair_square_sum = residual_trace_two + 3 * prime**2
    for sign in (-1, 1):
        pair_sum = residual_trace_one - sign * prime
        pair_product_numerator = pair_sum**2 - pair_square_sum
        if pair_product_numerator % 2:
            valid = False
            discriminant = None
            pair_product = None
        else:
            pair_product = pair_product_numerator // 2
            discriminant = pair_sum**2 - 4 * pair_product
            valid = discriminant >= 0
            if valid:
                lower_rhs = pair_sum + 4 * prime
                upper_rhs = 4 * prime - pair_sum
                valid = (
                    lower_rhs >= 0
                    and upper_rhs >= 0
                    and discriminant <= lower_rhs**2
                    and discriminant <= upper_rhs**2
                )
        sign_records.append(
            {
                "real_eigenvalue_sign": sign,
                "pair_sum_U_plus_V": pair_sum,
                "pair_product_U_times_V": pair_product,
                "pair_sum_discriminant": discriminant,
                "satisfies_exact_Weil_interval": valid,
            }
        )
        if valid:
            valid_signs.append((sign, pair_sum, pair_product, discriminant))
    if valid_signs != [(-1, -94, 2146, 252)]:
        raise AssertionError(f"the residual Weil reconstruction changed: {valid_signs}")

    sign, pair_sum, pair_product, discriminant = valid_signs[0]
    variable = sp.symbols("X")
    quartic = sp.Poly(
        variable**4
        - pair_sum * variable**3
        + (pair_product + 2 * prime**2) * variable**2
        - prime**2 * pair_sum * variable
        + prime**4,
        variable,
        domain=sp.ZZ,
    )
    residual = sp.Poly(
        sp.expand((variable - sign * prime) * quartic.as_expr()),
        variable,
        domain=sp.ZZ,
    )
    residual_at_positive_p = int(residual.eval(prime))
    if residual_at_positive_p != 534704436:
        raise AssertionError("the residual +p eigenvalue exclusion changed")
    normalized_variable = sp.symbols("Z")
    normalized_quartic = sp.Poly(
        sp.expand(quartic.as_expr().subs(variable, prime * normalized_variable) / prime**4),
        normalized_variable,
        domain=sp.QQ,
    )
    cyclotomic_gcds = []
    for order in ROOT_OF_UNITY_ORDERS_OF_DEGREE_AT_MOST_FOUR:
        cyclotomic = sp.Poly(
            sp.cyclotomic_poly(order, normalized_variable),
            normalized_variable,
            domain=sp.QQ,
        )
        common = sp.gcd(normalized_quartic, cyclotomic)
        cyclotomic_gcds.append(
            {"order": order, "gcd_degree": common.degree()}
        )
        if common.degree() != 0:
            raise AssertionError("the residual quartic gained a root of unity")

    if residual_trace_one != EXPECTED_RESIDUAL_TRACES[1] or residual_trace_two != EXPECTED_RESIDUAL_TRACES[2]:
        raise AssertionError("the pinned residual traces changed")
    return {
        "known_rational_Neron_Severi_rank": KNOWN_RATIONAL_NS_RANK,
        "residual_traces": {"degree_1": residual_trace_one, "degree_2": residual_trace_two},
        "reciprocal_pair_reconstruction": sign_records,
        "unique_real_eigenvalue": str(sign * prime),
        "pair_sums_U_V": ["-47+3*sqrt(7)", "-47-3*sqrt(7)"],
        "residual_quartic": str(quartic.as_expr()),
        "residual_characteristic_polynomial": str(residual.as_expr()),
        "residual_characteristic_polynomial_at_positive_p": residual_at_positive_p,
        "normalized_residual_quartic": str(normalized_quartic.as_expr()),
        "cyclotomic_gcd_checks": cyclotomic_gcds,
        "residual_p_times_root_of_unity_multiplicity_upper_bound": 1,
        "residual_positive_p_eigenvalue_multiplicity": 0,
        "rational_Neron_Severi_rank": KNOWN_RATIONAL_NS_RANK,
        "arithmetic_generic_Mordell_Weil_rank_over_Q_of_T": (
            KNOWN_RATIONAL_NS_RANK - TRIVIAL_LATTICE_RANK
        ),
        "geometric_Picard_rank_upper_bound": KNOWN_RATIONAL_NS_RANK + 1,
        "geometric_generic_Mordell_Weil_rank_upper_bound": (
            KNOWN_RATIONAL_NS_RANK + 1 - TRIVIAL_LATTICE_RANK
        ),
    }


def verify_picard_bound() -> dict[str, Any]:
    """Run the complete exact point-count and Frobenius certificate."""

    good_reduction = verify_good_reduction()
    counts = {
        degree: surface_point_count(GOOD_PRIME, degree) for degree in (1, 2)
    }
    for degree in (1, 2):
        if counts[degree]["surface_point_count"] != EXPECTED_SURFACE_POINT_COUNTS[degree]:
            raise AssertionError("the K3 surface point count changed")
        if counts[degree]["h2_frobenius_trace"] != EXPECTED_H2_TRACES[degree]:
            raise AssertionError("the H2 Frobenius trace changed")
    residual = reconstruct_residual_frobenius(
        GOOD_PRIME,
        counts[1]["h2_frobenius_trace"],
        counts[2]["h2_frobenius_trace"],
    )
    return {
        "good_reduction": good_reduction,
        "point_counts": {str(degree): counts[degree] for degree in (1, 2)},
        "residual_frobenius": residual,
        "inference": (
            "smooth proper specialization injects the characteristic-zero Neron--Severi group; "
            "a geometric divisor class has Frobenius eigenvalue p times a root of unity. "
            "The exact residual factor has only the eigenvalue -p of that form, so rho<=18. "
            "It has no further +p eigenvalue, so the 17 explicit Q-defined divisor classes "
            "exhaust NS over Q. Shioda--Tate for 20 I1 plus split I4 therefore gives "
            "arithmetic generic rank exactly 12 and geometric generic rank at most 13."
        ),
        "proved_arithmetic_generic_rank_over_Q_of_T": 12,
        "proved_generic_rank_interval": [12, 13],
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic-curves/elliptic_nagao_section7_picard_bound.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    result = verify_picard_bound()
    script_path = Path(__file__).resolve()
    linear_verifier = script_path.with_name("verify_nagao_section7_linear_sections.py")
    artifact = {
        "schema_version": 1,
        "status": "exact finite-field computation and arithmetic generic-rank theorem",
        **result,
        "scope_limits": [
            "does not determine whether the geometric generic rank over Qbar(T) is 12 or 13",
            "does not imply rank 21 at any specialization",
        ],
        "target_hit": False,
        "software": {"python": platform.python_version(), "sympy": sp.__version__},
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "rank12_certificate_verifier_sha256": hashlib.sha256(
            linear_verifier.read_bytes()
        ).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print("#S(F_29)=1212; #S(F_29^2)=723600")
    print("proved arithmetic generic Mordell--Weil rank over Q(T) equals 12")
    print("proved 12 <= geometric generic Mordell--Weil rank <= 13")
    print(args.output)


if __name__ == "__main__":
    main()
