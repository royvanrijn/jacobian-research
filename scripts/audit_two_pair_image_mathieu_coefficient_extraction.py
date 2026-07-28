#!/usr/bin/env python3
"""Independent exact audit of the direct SIC2C4 coefficient proof.

This script does not use Gaussian moments, radial variables, Hopf
coordinates, or the Taylor-jet calculation from the first proof.  It checks
the algebraic chart identity and the finite-difference / recurrence
certificates in Section 4 of

    extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md.

The loops are regressions, not a replacement for the all-order proof.  The
universal certificate consists of:

* the monomial constant-term formula (audited over a range of bidegrees);
* degree < m for the pure finite-difference polynomial;
* divisibility of A_m(X)-1 by 2X+1 with quotient degree <= m-2; and
* the general denominator-remainder invariance of the quotient transform; and
* the repeated-pole jet formula and triangular generalized-beta recurrence;
* the termwise recurrence
      (2m+1) B_m - 2m B_(m-1) = sum_k (-1)^k binom(m,k).
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial


Poly = list[Fraction]  # coefficients in ascending powers of X
Bivariate = dict[tuple[int, int], Fraction]  # Laurent exponents of (x, v)
CERTIFICATE_CUTOFF = 160
MONOMIAL_CUTOFF = 12
REMAINDER_CUTOFF = 40
CHART_EXPANSION_CUTOFF = 10
REPEATED_POLE_CUTOFF = 40
MAX_POLE_ORDER = 6


def trim(poly: Poly) -> Poly:
    result = poly[:]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_scale(poly: Poly, scalar: Fraction | int) -> Poly:
    value = Fraction(scalar)
    return trim([value * coefficient for coefficient in poly])


def poly_add(left: Poly, right: Poly) -> Poly:
    result = [Fraction(0)] * max(len(left), len(right))
    for degree, coefficient in enumerate(left):
        result[degree] += coefficient
    for degree, coefficient in enumerate(right):
        result[degree] += coefficient
    return trim(result)


def poly_multiply(left: Poly, right: Poly) -> Poly:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    return trim(result)


def poly_power(poly: Poly, exponent: int) -> Poly:
    result = [Fraction(1)]
    base = poly
    value = exponent
    while value:
        if value % 2:
            result = poly_multiply(result, base)
        base = poly_multiply(base, base)
        value //= 2
    return result


def poly_evaluate(poly: Poly, value: Fraction | int) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def poly_divmod(dividend: Poly, divisor: Poly) -> tuple[Poly, Poly]:
    """Exact polynomial division over the rationals."""
    numerator = trim(dividend)
    denominator = trim(divisor)
    assert denominator != [Fraction(0)]
    if len(numerator) < len(denominator):
        return [Fraction(0)], numerator
    quotient = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    remainder = numerator[:]
    while (
        remainder != [Fraction(0)]
        and len(remainder) >= len(denominator)
    ):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[shift] += coefficient
        for degree, denominator_coefficient in enumerate(denominator):
            remainder[degree + shift] -= coefficient * denominator_coefficient
        remainder = trim(remainder)
    return trim(quotient), remainder


def divide_by_two_x_plus_one(poly: Poly) -> tuple[Poly, Fraction]:
    """Return quotient and remainder for division by 2*X+1."""
    work = trim(poly)
    if len(work) == 1:
        return [Fraction(0)], work[0]
    quotient = [Fraction(0)] * (len(work) - 1)
    remainder_work = work[:]
    for degree in range(len(work) - 1, 0, -1):
        coefficient = remainder_work[degree] / 2
        quotient[degree - 1] = coefficient
        remainder_work[degree] -= 2 * coefficient
        remainder_work[degree - 1] -= coefficient
    assert all(value == 0 for value in remainder_work[1:])
    return trim(quotient), remainder_work[0]


def finite_difference_sum(order: int, poly: Poly) -> Fraction:
    return sum(
        Fraction((-1) ** index * comb(order, index))
        * poly_evaluate(poly, index)
        for index in range(order + 1)
    )


def b_add(*polynomials: Bivariate) -> Bivariate:
    result: Bivariate = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def b_scale(polynomial: Bivariate, scalar: Fraction | int) -> Bivariate:
    value = Fraction(scalar)
    return {
        exponent: value * coefficient
        for exponent, coefficient in polynomial.items()
        if value * coefficient
    }


def b_multiply(left: Bivariate, right: Bivariate) -> Bivariate:
    result: Bivariate = {}
    for (left_x, left_v), left_coefficient in left.items():
        for (right_x, right_v), right_coefficient in right.items():
            exponent = (left_x + right_x, left_v + right_v)
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def b_power(polynomial: Bivariate, exponent: int) -> Bivariate:
    result: Bivariate = {(0, 0): Fraction(1)}
    base = polynomial
    value = exponent
    while value:
        if value % 2:
            result = b_multiply(result, base)
        base = b_multiply(base, base)
        value //= 2
    return result


def odd_double_factorial(order: int) -> int:
    result = 1
    for value in range(1, order + 1, 2):
        result *= value
    return result


def formal_beta(first: int, second: int) -> Fraction:
    """Apply I(t^j)=1/(j+1) after expanding t^first*(1-t)^second."""
    return sum(
        Fraction((-1) ** index * comb(second, index), first + index + 1)
        for index in range(second + 1)
    )


def audit_chart_identity() -> None:
    one: Bivariate = {(0, 0): Fraction(1)}
    x: Bivariate = {(1, 0): Fraction(1)}
    x_inverse: Bivariate = {(-1, 0): Fraction(1)}
    v: Bivariate = {(0, 1): Fraction(1)}
    v_squared = b_power(v, 2)
    w = b_scale(
        b_multiply(x_inverse, b_add(one, b_scale(v_squared, -1))),
        Fraction(1, 2),
    )

    # Build F=(R+Z)(R^2 W-(2R+Z)T^2/2) with R=1, Z=x, T=v.
    built = b_multiply(
        b_add(one, x),
        b_add(
            w,
            b_scale(
                b_multiply(b_add(b_scale(one, 2), x), v_squared),
                Fraction(-1, 2),
            ),
        ),
    )
    target = b_scale(
        b_multiply(
            b_multiply(b_add(one, x), x_inverse),
            b_add(
                one,
                b_scale(
                    b_multiply(v_squared, b_power(b_add(one, x), 2)),
                    -1,
                ),
            ),
        ),
        Fraction(1, 2),
    )
    assert built == target


def audit_chart_constant_terms() -> None:
    """Expand the chart itself and check both constant-term formulas."""
    one: Bivariate = {(0, 0): Fraction(1)}
    x: Bivariate = {(1, 0): Fraction(1)}
    x_inverse: Bivariate = {(-1, 0): Fraction(1)}
    v_squared: Bivariate = {(0, 2): Fraction(1)}
    chart_f = b_scale(
        b_multiply(
            b_multiply(b_add(one, x), x_inverse),
            b_add(
                one,
                b_scale(
                    b_multiply(v_squared, b_power(b_add(one, x), 2)),
                    -1,
                ),
            ),
        ),
        Fraction(1, 2),
    )
    for order in range(1, CHART_EXPANSION_CUTOFF + 1):
        powered = b_power(chart_f, order)
        pure_constant_term = {
            exponent: coefficient
            for exponent, coefficient in powered.items()
            if exponent[0] == 0
        }
        mixed_constant_term = {
            exponent: coefficient
            for exponent, coefficient in b_multiply(x, powered).items()
            if exponent[0] == 0
        }
        expected_pure = {
            (0, 2 * index): Fraction(
                (-1) ** index
                * comb(order, index)
                * comb(order + 2 * index, order),
                2**order,
            )
            for index in range(order + 1)
        }
        expected_mixed = {
            (0, 2 * index): Fraction(
                (-1) ** index
                * comb(order, index)
                * comb(order + 2 * index, order - 1),
                2**order,
            )
            for index in range(order + 1)
        }
        assert pure_constant_term == expected_pure
        assert mixed_constant_term == expected_mixed


def audit_monomial_coefficient_formula() -> None:
    # For xi1^a xi2^(n-a) z1^b z2^(n-b), CT_u forces a=b.
    # Independently expand the beta polynomial under I(t^j)=1/(j+1).
    for degree in range(MONOMIAL_CUTOFF + 1):
        for dual_index in range(degree + 1):
            for coordinate_index in range(degree + 1):
                direct = (
                    factorial(dual_index) * factorial(degree - dual_index)
                    if dual_index == coordinate_index
                    else 0
                )
                extracted = (
                    factorial(degree + 1)
                    * formal_beta(coordinate_index, degree - coordinate_index)
                    if dual_index == coordinate_index
                    else 0
                )
                assert direct == extracted
                if dual_index == coordinate_index:
                    beta = formal_beta(
                        coordinate_index,
                        degree - coordinate_index,
                    )
                    assert beta == Fraction(
                        factorial(coordinate_index)
                        * factorial(degree - coordinate_index),
                        factorial(degree + 1),
                    )
                    if coordinate_index:
                        assert (
                            (degree - coordinate_index + 1) * beta
                            == coordinate_index
                            * formal_beta(
                                coordinate_index - 1,
                                degree - coordinate_index + 1,
                            )
                        )


def audit_all_order_certificates() -> None:
    previous_b = Fraction(1)
    for order in range(1, CERTIFICATE_CUTOFF + 1):
        # Equation (4.7): degree m-1, hence its m-th difference is zero.
        pure_poly: Poly = [Fraction(1)]
        for shift in range(2, order + 1):
            pure_poly = poly_multiply(pure_poly, [Fraction(shift), Fraction(2)])
        pure_poly = poly_scale(pure_poly, Fraction(1, factorial(order)))
        assert len(pure_poly) - 1 <= order - 1
        assert finite_difference_sum(order, pure_poly) == 0

        # Equations (4.8)-(4.9): A_m(-1/2)=1, so A_m-1 is divisible
        # by 2X+1 and the quotient has degree at most m-2.
        a_poly: Poly = [Fraction(1)]
        for shift in range(2, order + 1):
            a_poly = poly_multiply(a_poly, [Fraction(shift), Fraction(2)])
        a_poly = poly_scale(a_poly, Fraction(1, factorial(order - 1)))
        a_minus_one = a_poly[:]
        a_minus_one[0] -= 1
        quotient, remainder = divide_by_two_x_plus_one(a_minus_one)
        assert remainder == 0
        assert quotient == [Fraction(0)] or len(quotient) - 1 <= order - 2
        assert finite_difference_sum(order, quotient) == 0

        # Equation (4.10), followed by the termwise recurrence (4.11).
        current_b = sum(
            Fraction((-1) ** index * comb(order, index), 2 * index + 1)
            for index in range(order + 1)
        )
        alternating_binomial = sum(
            (-1) ** index * comb(order, index)
            for index in range(order + 1)
        )
        assert (
            (2 * order + 1) * current_b - 2 * order * previous_b
            == alternating_binomial
            == 0
        )
        assert current_b == Fraction(
            2**order * factorial(order),
            odd_double_factorial(2 * order + 1),
        )

        # Equations (4.4)-(4.5), evaluated from their displayed sums.
        pure_sum = Fraction(0)
        mixed_sum = Fraction(0)
        for index in range(order + 1):
            signed_choose = (-1) ** index * comb(order, index)
            pure_sum += Fraction(
                signed_choose * comb(order + 2 * index, order),
                2**order * (2 * index + 1),
            )
            mixed_sum += Fraction(
                signed_choose * comb(order + 2 * index, order - 1),
                2**order * (2 * index + 1),
            )
        assert pure_sum == 0
        assert mixed_sum == Fraction(
            factorial(order),
            odd_double_factorial(2 * order + 1),
        )
        previous_b = current_b


def audit_rank_one_quotient_transform() -> None:
    """Check the strengthened endpoint-residue identity (4.8)."""
    for order in range(1, 61):
        # A deterministic degree-(m-1) test polynomial with mixed signs.
        polynomial = [
            Fraction((-1) ** degree * (degree + 2), degree + 1)
            for degree in range(order)
        ]
        transformed = sum(
            Fraction((-1) ** index * comb(order, index), 2 * index + 1)
            * poly_evaluate(polynomial, index)
            for index in range(order + 1)
        )
        beta_sum = sum(
            Fraction((-1) ** index * comb(order, index), 2 * index + 1)
            for index in range(order + 1)
        )
        assert transformed == poly_evaluate(polynomial, Fraction(-1, 2)) * beta_sum


def quotient_transform(order: int, denominator: Poly, numerator: Poly) -> Fraction:
    return sum(
        Fraction((-1) ** index * comb(order, index))
        * poly_evaluate(numerator, index)
        / poly_evaluate(denominator, index)
        for index in range(order + 1)
    )


def repeated_beta_sum(order: int, pole_order: int) -> Fraction:
    return sum(
        Fraction(
            (-1) ** index * comb(order, index),
            (2 * index + 1) ** pole_order,
        )
        for index in range(order + 1)
    )


def generalized_beta_product_jet(order: int, jet_length: int) -> list[Fraction]:
    """Taylor coefficients at a=1 of 2^m*m!/prod_{ell=0}^m(a+2*ell)."""
    denominator = [Fraction(1)]
    for index in range(order + 1):
        denominator = poly_multiply(
            denominator,
            [Fraction(2 * index + 1), Fraction(1)],
        )[:jet_length]
    reciprocal = [Fraction(0)] * jet_length
    reciprocal[0] = 1 / denominator[0]
    for degree in range(1, jet_length):
        reciprocal[degree] = -sum(
            (
                denominator[index] if index < len(denominator) else 0
            )
            * reciprocal[degree - index]
            for index in range(1, degree + 1)
        ) / denominator[0]
    return [
        Fraction(2**order * factorial(order)) * coefficient
        for coefficient in reciprocal
    ]


def audit_general_remainder_transform() -> None:
    """Check the degree-r denominator-remainder invariance (4.8a)-(4.8b)."""
    for order in range(1, REMAINDER_CUTOFF + 1):
        for rank in range(1, 5):
            # L=(X+1/2)^r is monic and nonzero at every node k >= 0.
            denominator = poly_power([Fraction(1, 2), Fraction(1)], rank)
            remainder = [
                Fraction((-1) ** degree * (order + degree + 1), degree + 2)
                for degree in range(rank)
            ]
            quotient = [
                Fraction((-1) ** degree * (rank + degree + 1), degree + 1)
                for degree in range(order)
            ]
            numerator = poly_add(
                remainder,
                poly_multiply(denominator, quotient),
            )
            recovered_quotient, recovered_remainder = poly_divmod(
                numerator,
                denominator,
            )
            assert recovered_quotient == quotient
            assert recovered_remainder == trim(remainder)
            assert len(quotient) - 1 < order
            assert len(numerator) - 1 < order + rank
            assert quotient_transform(order, denominator, numerator) == (
                quotient_transform(order, denominator, remainder)
            )


def audit_repeated_pole_jet_transform() -> None:
    """Check the triangular recurrence and explicit jet formula (4.8c)-(4.8d)."""
    local_coordinate = [Fraction(1), Fraction(2)]  # lambda=2X+1
    for order in range(1, REPEATED_POLE_CUTOFF + 1):
        assert repeated_beta_sum(order, 0) == 0
        product_jet = generalized_beta_product_jet(order, MAX_POLE_ORDER)
        for pole_order in range(1, MAX_POLE_ORDER + 1):
            current = repeated_beta_sum(order, pole_order)
            assert (2 * order + 1) * current == (
                repeated_beta_sum(order, pole_order - 1)
                + 2 * order * repeated_beta_sum(order - 1, pole_order)
            )
            assert current == (
                (-1) ** (pole_order - 1) * product_jet[pole_order - 1]
            )

        for pole_order in range(1, 5):
            coefficients = [
                Fraction(
                    (-1) ** degree * (order + pole_order + degree),
                    degree + 1,
                )
                for degree in range(pole_order)
            ]
            remainder = [Fraction(0)]
            for degree, coefficient in enumerate(coefficients):
                remainder = poly_add(
                    remainder,
                    poly_scale(
                        poly_power(local_coordinate, degree),
                        coefficient,
                    ),
                )
            quotient = [
                Fraction((-1) ** degree * (degree + 2), degree + 1)
                for degree in range(order)
            ]
            denominator = poly_power(local_coordinate, pole_order)
            numerator = poly_add(
                remainder,
                poly_multiply(denominator, quotient),
            )
            expected = sum(
                coefficient
                * repeated_beta_sum(order, pole_order - degree)
                for degree, coefficient in enumerate(coefficients)
            )
            assert quotient_transform(order, denominator, numerator) == expected


def main() -> None:
    audit_chart_identity()
    audit_chart_constant_terms()
    audit_monomial_coefficient_formula()
    audit_all_order_certificates()
    audit_rank_one_quotient_transform()
    audit_general_remainder_transform()
    audit_repeated_pole_jet_transform()
    print("PASS SIC2C4 coefficient proof: formal contraction chart identity")
    print(
        "PASS SIC2C4 coefficient proof: chart constant-term expansions "
        f"through m={CHART_EXPANSION_CUTOFF}"
    )
    print(
        "PASS SIC2C4 coefficient proof: monomial CT-beta formula "
        f"through bidegree {MONOMIAL_CUTOFF}"
    )
    print(
        "PASS SIC2C4 coefficient proof: finite-difference divisibility and "
        f"recurrence certificates through m={CERTIFICATE_CUTOFF}"
    )
    print(
        "PASS SIC2C4 coefficient proof: rank-one quotient transform "
        "through m=60"
    )
    print(
        "PASS SIC2C4 coefficient proof: degree-r remainder transform "
        f"for r<=4 through m={REMAINDER_CUTOFF}"
    )
    print(
        "PASS SIC2C4 coefficient proof: repeated-pole jet recurrence/product "
        f"for s<={MAX_POLE_ORDER} through m={REPEATED_POLE_CUTOFF}"
    )
    print("PASS SIC2C4 coefficient proof: no Gaussian or Hopf input")


if __name__ == "__main__":
    main()
