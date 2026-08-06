#!/usr/bin/env python3
"""Certify quotient-module descent through mu_7 on the cubic Hurwitz chart.

After eliminating ``a3`` with ``mu_2``, write ``mu_3`` as

    A*a2^2 + B*a2 + C.

On ``A=0`` and the ``M01`` channel-minor open one has ``b0 != 0``.  Put
``q=b1/b0`` and ``z=b0*a2``; then ``A`` eliminates ``b2`` and the transformed
``mu_3`` is ``L=B*z+C``.

This checker proves two exact irreducibility statements by primitive reduction
modulo 29 followed by Singular factorization:

* on ``B != 0``, eliminating ``z`` between ``L`` and ``mu_4`` gives an
  irreducible polynomial in ``(a1,b0,q,lambda)``;
* on the degree-drop boundary ``B=C=0``, reducing the cubic ``C`` modulo the
  constant-leading quadratic ``B`` in ``a1`` gives ``U+V*a1``.  Its norm is an
  irreducible polynomial in ``(b0,q,lambda)``.  The descended ``mu_4`` is
  cubic in ``z``; cancelling the cubic lead against descended ``mu_5`` gives
  the next irreducible quadratic-in-``z`` subresultant and exposes its exact
  leading-coefficient boundaries;
* the cubic--quadratic norm gives an irreducible degree-118 base equation,
  and the descended mu_6 and mu_7 module determinants give irreducible
  degree-130 and degree-162 base equations.

The claims are over characteristic zero.  The modular calculations are not
used as evidence about fibres: equality with the primitive exact polynomials
is checked coefficient by coefficient before irreducibility is inferred.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess

from flint import fmpz_mpoly_ctx, nmod_mpoly_ctx

from research_two_pair_sic_bidegree33_rank_two_hurwitz import (
    EXACT_LIFT_PRIME,
    RationalParameterPolynomial,
    active_parameter_polynomial_string,
    add_parameter_polynomial,
    add_rational_polynomial,
    base_polynomials,
    channel_minor_polynomial,
    coefficient_groups,
    divide_rational_polynomials,
    exact_moment_polynomials,
    exact_parameter_polynomial,
    linear_substitution_numerator,
    moment,
    multiply_parameter_polynomials,
    multiply_rational_polynomials,
    parameter_polynomial_power,
    powers,
    primitive_integer_coefficients,
    rational_linear_substitution_numerator,
    rational_polynomial,
    rational_polynomial_power,
    split_linear_parameter_polynomial,
    split_linear_rational_polynomial,
    strip_parameter_factor,
    strip_rational_factor,
    substitute_parameter_variable,
    substitute_rational_variable,
    zero_exponent,
)
from research_two_pair_sic_bidegree33_rank_two_hurwitz_linear_incidence import (
    finite_ratio_replacements,
    rational_coefficient_groups,
    rational_ratio_replacements,
    rescale_a2_by_b0,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_hurwitz_module_descent.json"
)
PRIME = 29
ORDERS = (2, 3, 4, 5, 6, 7)


def polynomial_degree(polynomial: dict[tuple[int, ...], object]) -> int:
    return max((sum(exponent) for exponent in polynomial), default=-1)


def variable_degrees(
    polynomial: dict[tuple[int, ...], object], variables: tuple[int, ...]
) -> dict[str, int]:
    names = ("lambda", "a1", "z", "a3", "b0", "q", "b2")
    return {
        names[variable]: max(
            (exponent[variable] for exponent in polynomial), default=-1
        )
        for variable in variables
    }


def plain_coefficient_groups(
    polynomial: dict[tuple[int, ...], int], variable: int
) -> dict[int, dict[tuple[int, ...], int]]:
    groups: dict[int, dict[tuple[int, ...], int]] = {}
    for exponent, coefficient in polynomial.items():
        reduced = list(exponent)
        power = reduced[variable]
        reduced[variable] = 0
        groups.setdefault(power, {})[tuple(reduced)] = coefficient
    return groups


def flint_coefficient_polynomial(polynomial, context):
    coefficients = {}
    for exponent, coefficient in polynomial.items():
        if any(exponent[index] for index in (1, 2, 3, 6)):
            raise AssertionError("unexpected variable in a base coefficient")
        coefficients[(exponent[4], exponent[5], exponent[0])] = int(
            coefficient
        )
    return context.from_dict(coefficients)


def flint_linear_norm_descent(
    cubic: dict[tuple[int, ...], int],
    quadratic: dict[tuple[int, ...], int],
    context,
):
    cubic_groups = plain_coefficient_groups(cubic, 2)
    quadratic_groups = plain_coefficient_groups(quadratic, 2)
    assert set(cubic_groups) == {0, 1, 2, 3}
    assert set(quadratic_groups) == {0, 1, 2}
    a0, a1, a2, a3 = (
        flint_coefficient_polynomial(cubic_groups[power], context)
        for power in range(4)
    )
    b0, b1, b2 = (
        flint_coefficient_polynomial(quadratic_groups[power], context)
        for power in range(3)
    )
    remainder_constant = (
        b2**2 * a0 - b2 * a2 * b0 + a3 * b1 * b0
    )
    remainder_linear = (
        b2**2 * a1
        - b2 * a2 * b1
        + a3 * (b1**2 - b2 * b0)
    )
    norm = (
        b2 * remainder_constant**2
        - b1 * remainder_constant * remainder_linear
        + b0 * remainder_linear**2
    )
    known_factor = (
        context.gen(0) ** 2 * a3**2 * b2**2
    )
    residual, residual_remainder = divmod(norm, known_factor)
    assert not residual_remainder
    return {
        "remainder_constant": remainder_constant,
        "remainder_linear": remainder_linear,
        "norm": norm,
        "residual": residual,
    }


def flint_quadratic_module_remainder(
    polynomial: dict[tuple[int, ...], int],
    quadratic: dict[tuple[int, ...], int],
    context,
):
    groups = plain_coefficient_groups(polynomial, 2)
    quadratic_groups = plain_coefficient_groups(quadratic, 2)
    assert set(quadratic_groups) == {0, 1, 2}
    degree = max(groups)
    q0, q1, q2 = (
        flint_coefficient_polynomial(quadratic_groups[power], context)
        for power in range(3)
    )
    zero = context.constant(0)
    one = context.constant(1)
    remainders = [(one, zero), (zero, one), (-q0, -q1)]
    q2_q0 = q2 * q0
    for power in range(3, degree + 1):
        previous_constant, previous_linear = remainders[power - 1]
        earlier_constant, earlier_linear = remainders[power - 2]
        remainders.append(
            (
                -q1 * previous_constant - q2_q0 * earlier_constant,
                -q1 * previous_linear - q2_q0 * earlier_linear,
            )
        )
    q2_powers = [q2**power for power in range(degree)]
    answer_constant = zero
    answer_linear = zero
    for power, group in groups.items():
        coefficient = flint_coefficient_polynomial(group, context)
        multiplier = q2_powers[
            degree - 1 if power <= 1 else degree - power
        ]
        reduced_constant, reduced_linear = remainders[power]
        answer_constant += coefficient * multiplier * reduced_constant
        answer_linear += coefficient * multiplier * reduced_linear
    return answer_constant, answer_linear


def flint_later_base_equation(
    polynomial: dict[tuple[int, ...], int],
    cubic: dict[tuple[int, ...], int],
    quadratic: dict[tuple[int, ...], int],
    linear_norm: dict[str, object],
    context,
    b0_power: int | None,
):
    module_constant, module_linear = flint_quadratic_module_remainder(
        polynomial, quadratic, context
    )
    determinant = (
        linear_norm["remainder_linear"] * module_constant
        - linear_norm["remainder_constant"] * module_linear
    )
    cubic_lead = flint_coefficient_polynomial(
        plain_coefficient_groups(cubic, 2)[3], context
    )
    quadratic_lead = flint_coefficient_polynomial(
        plain_coefficient_groups(quadratic, 2)[2], context
    )
    if b0_power is None:
        residual = determinant
        removed_factors = {}
        for name, factor in (
            ("b0", context.gen(0)),
            ("lc_z_N4", cubic_lead),
            ("lc_z_S5", quadratic_lead),
        ):
            exponent = 0
            while True:
                quotient, remainder = divmod(residual, factor)
                if remainder:
                    break
                residual = quotient
                exponent += 1
            removed_factors[name] = exponent
    else:
        known_factor = (
            context.gen(0) ** b0_power
            * cubic_lead
            * quadratic_lead**2
        )
        residual, remainder = divmod(determinant, known_factor)
        assert not remainder
        removed_factors = {
            "b0": b0_power,
            "lc_z_N4": 1,
            "lc_z_S5": 2,
        }
    return {
        "module_constant": module_constant,
        "module_linear": module_linear,
        "determinant": determinant,
        "residual": residual,
        "removed_factors": removed_factors,
    }


def compare_flint_reduction(exact, modular, prime: int) -> int:
    reduction = {
        exponent: int(coefficient) % prime
        for exponent, coefficient in exact.to_dict().items()
        if int(coefficient) % prime
    }
    modular_dict = {
        exponent: int(coefficient) % prime
        for exponent, coefficient in modular.to_dict().items()
        if int(coefficient) % prime
    }
    assert set(reduction) == set(modular_dict)
    lead = max(reduction)
    scalar = (
        modular_dict[lead] * pow(reduction[lead], -1, prime) % prime
    )
    assert {
        exponent: coefficient * scalar % prime
        for exponent, coefficient in reduction.items()
    } == modular_dict
    assert exact.total_degree() == modular.total_degree()
    return scalar


def full_parameter_polynomial_from_flint(polynomial) -> dict[tuple[int, ...], int]:
    answer = {}
    for exponent, coefficient in polynomial.to_dict().items():
        b0_power, q_power, lambda_power = exponent
        answer[(lambda_power, 0, 0, 0, b0_power, q_power, 0)] = int(
            coefficient
        )
    return answer


def flint_profile(polynomial) -> dict[str, object]:
    degrees = tuple(int(value) for value in polynomial.degrees())
    return {
        "terms": len(polynomial.to_dict()),
        "total_degree": int(polynomial.total_degree()),
        "variable_degrees": {
            "b0": degrees[0],
            "q": degrees[1],
            "lambda": degrees[2],
        },
    }


def rational_quadratic_remainder_numerator(
    polynomial: RationalParameterPolynomial,
    variable: int,
    quadratic: tuple[
        RationalParameterPolynomial,
        RationalParameterPolynomial,
        RationalParameterPolynomial,
    ],
) -> tuple[RationalParameterPolynomial, RationalParameterPolynomial]:
    """Reduce modulo A*x^2+B*x+C while clearing the minimal power of A."""

    a_polynomial, b_polynomial, c_polynomial = quadratic
    groups = rational_coefficient_groups(polynomial, variable)
    maximum = max(groups, default=0)
    zero: RationalParameterPolynomial = {}
    one: RationalParameterPolynomial = {zero_exponent(): Fraction(1)}
    if maximum == 0:
        return groups.get(0, {}), zero
    if maximum == 1:
        return groups.get(0, {}), groups.get(1, {})

    # R_k=A^(k-1)*x^k modulo A*x^2+B*x+C.  The k=2 base is
    # (-C,-B); only the k>=3 recurrence contains A*C.
    remainders: list[
        tuple[RationalParameterPolynomial, RationalParameterPolynomial]
    ] = [(one, zero), (zero, one)]
    remainders.append(
        (
            {exponent: -coefficient for exponent, coefficient in c_polynomial.items()},
            {exponent: -coefficient for exponent, coefficient in b_polynomial.items()},
        )
    )
    a_times_c = multiply_rational_polynomials(a_polynomial, c_polynomial)
    for _ in range(3, maximum + 1):
        previous_constant, previous_linear = remainders[-1]
        earlier_constant, earlier_linear = remainders[-2]
        constant = multiply_rational_polynomials(
            b_polynomial, previous_constant
        )
        add_rational_polynomial(
            constant,
            multiply_rational_polynomials(a_times_c, earlier_constant),
        )
        constant = {
            exponent: -coefficient
            for exponent, coefficient in constant.items()
        }
        linear = multiply_rational_polynomials(b_polynomial, previous_linear)
        add_rational_polynomial(
            linear,
            multiply_rational_polynomials(a_times_c, earlier_linear),
        )
        linear = {
            exponent: -coefficient for exponent, coefficient in linear.items()
        }
        remainders.append((constant, linear))

    a_powers = [
        rational_polynomial_power(a_polynomial, power)
        for power in range(maximum)
    ]
    answer_constant: RationalParameterPolynomial = {}
    answer_linear: RationalParameterPolynomial = {}
    for power, group in groups.items():
        multiplier = a_powers[maximum - 1 if power <= 1 else maximum - power]
        group_multiplier = multiply_rational_polynomials(group, multiplier)
        reduced_constant, reduced_linear = remainders[power]
        add_rational_polynomial(
            answer_constant,
            multiply_rational_polynomials(
                group_multiplier, reduced_constant
            ),
        )
        add_rational_polynomial(
            answer_linear,
            multiply_rational_polynomials(group_multiplier, reduced_linear),
        )
    return answer_constant, answer_linear


def shifted_rational_polynomial(
    polynomial: RationalParameterPolynomial, variable: int, power: int
) -> RationalParameterPolynomial:
    answer: RationalParameterPolynomial = {}
    for exponent, coefficient in polynomial.items():
        shifted = list(exponent)
        shifted[variable] += power
        answer[tuple(shifted)] = coefficient
    return answer


def verify_rational_quadratic_remainder(
    polynomial: RationalParameterPolynomial,
    variable: int,
    quadratic: tuple[
        RationalParameterPolynomial,
        RationalParameterPolynomial,
        RationalParameterPolynomial,
    ],
    constant: RationalParameterPolynomial,
    linear: RationalParameterPolynomial,
) -> None:
    """Check A^(d-1)*P == constant+linear*x modulo the quadratic."""

    a_polynomial, b_polynomial, c_polynomial = quadratic
    degree = max(
        (exponent[variable] for exponent in polynomial), default=0
    )
    clearing = rational_polynomial_power(
        a_polynomial, max(0, degree - 1)
    )
    difference = multiply_rational_polynomials(polynomial, clearing)
    add_rational_polynomial(difference, constant, Fraction(-1))
    add_rational_polynomial(
        difference,
        shifted_rational_polynomial(linear, variable, 1),
        Fraction(-1),
    )
    divisor = shifted_rational_polynomial(a_polynomial, variable, 2)
    add_rational_polynomial(
        divisor, shifted_rational_polynomial(b_polynomial, variable, 1)
    )
    add_rational_polynomial(divisor, c_polynomial)
    divide_rational_polynomials(difference, divisor)


def modular_quadratic_remainder_numerator(
    polynomial: dict[tuple[int, ...], int],
    variable: int,
    quadratic: tuple[
        dict[tuple[int, ...], int],
        dict[tuple[int, ...], int],
        dict[tuple[int, ...], int],
    ],
    prime: int,
) -> tuple[dict[tuple[int, ...], int], dict[tuple[int, ...], int]]:
    """Local copy of the rank-two recurrence, kept explicit for the audit."""

    a_polynomial, b_polynomial, c_polynomial = quadratic
    groups = coefficient_groups(polynomial, variable, prime)
    maximum = max(groups, default=0)
    zero: dict[tuple[int, ...], int] = {}
    one = {zero_exponent(): 1}
    if maximum == 0:
        return groups.get(0, {}), zero
    if maximum == 1:
        return groups.get(0, {}), groups.get(1, {})
    remainders = [(one, zero), (zero, one)]
    remainders.append(
        (
            {
                exponent: (-coefficient) % prime
                for exponent, coefficient in c_polynomial.items()
            },
            {
                exponent: (-coefficient) % prime
                for exponent, coefficient in b_polynomial.items()
            },
        )
    )
    a_times_c = multiply_parameter_polynomials(
        a_polynomial, c_polynomial, prime
    )
    for _ in range(3, maximum + 1):
        previous_constant, previous_linear = remainders[-1]
        earlier_constant, earlier_linear = remainders[-2]
        constant = multiply_parameter_polynomials(
            b_polynomial, previous_constant, prime
        )
        add_parameter_polynomial(
            constant,
            multiply_parameter_polynomials(
                a_times_c, earlier_constant, prime
            ),
            1,
            prime,
        )
        constant = {
            exponent: (-coefficient) % prime
            for exponent, coefficient in constant.items()
        }
        linear = multiply_parameter_polynomials(
            b_polynomial, previous_linear, prime
        )
        add_parameter_polynomial(
            linear,
            multiply_parameter_polynomials(
                a_times_c, earlier_linear, prime
            ),
            1,
            prime,
        )
        linear = {
            exponent: (-coefficient) % prime
            for exponent, coefficient in linear.items()
        }
        remainders.append((constant, linear))
    a_powers = [
        parameter_polynomial_power(a_polynomial, power, prime)
        for power in range(maximum)
    ]
    answer_constant: dict[tuple[int, ...], int] = {}
    answer_linear: dict[tuple[int, ...], int] = {}
    for power, group in groups.items():
        multiplier = a_powers[maximum - 1 if power <= 1 else maximum - power]
        group_multiplier = multiply_parameter_polynomials(
            group, multiplier, prime
        )
        reduced_constant, reduced_linear = remainders[power]
        add_parameter_polynomial(
            answer_constant,
            multiply_parameter_polynomials(
                group_multiplier, reduced_constant, prime
            ),
            1,
            prime,
        )
        add_parameter_polynomial(
            answer_linear,
            multiply_parameter_polynomials(
                group_multiplier, reduced_linear, prime
            ),
            1,
            prime,
        )
    return answer_constant, answer_linear


def exact_descent() -> dict[str, object]:
    moments = {
        order: rational_polynomial(polynomial)
        for order, polynomial in exact_moment_polynomials(ORDERS).items()
    }
    pivot, rest = split_linear_rational_polynomial(moments[2], 3)
    eliminated = {
        order: rational_linear_substitution_numerator(
            moments[order], 3, pivot, rest
        )
        for order in ORDERS
        if order > 2
    }

    pre_ratio_groups = rational_coefficient_groups(eliminated[3], 2)
    assert set(pre_ratio_groups) == {0, 1, 2}
    quadratic_lead = primitive_integer_coefficients(pre_ratio_groups[2])
    expected_lead = {
        (2, 0, 0, 0, 2, 0, 0): 3,
        (1, 0, 0, 0, 2, 0, 0): 1,
        (0, 0, 0, 0, 2, 0, 0): 6,
        (1, 0, 0, 0, 1, 1, 0): 1,
        (0, 0, 0, 0, 1, 1, 0): 5,
        (0, 0, 0, 0, 1, 0, 1): -8,
        (0, 0, 0, 0, 1, 0, 0): -8,
        (0, 0, 0, 0, 0, 2, 0): 4,
    }
    assert quadratic_lead == expected_lead
    minor = exact_parameter_polynomial(
        channel_minor_polynomial("01", EXACT_LIFT_PRIME)
    )
    assert minor == {
        (0, 0, 0, 0, 0, 1, 0): 1,
        (0, 1, 0, 0, 1, 0, 0): -1,
    }

    b1_replacement, b2_replacement = rational_ratio_replacements()

    def transform(polynomial: RationalParameterPolynomial):
        answer = substitute_rational_variable(
            polynomial, 5, b1_replacement
        )
        answer = substitute_rational_variable(
            answer, 6, b2_replacement
        )
        return rescale_a2_by_b0(answer)[0]

    transformed = {order: transform(polynomial) for order, polynomial in eliminated.items()}
    incidence_groups = rational_coefficient_groups(transformed[3], 2)
    assert set(incidence_groups) == {0, 1}
    constant, linear = incidence_groups[0], incidence_groups[1]
    projection = rational_linear_substitution_numerator(
        transformed[4], 2, linear, constant
    )

    quadratic = rational_coefficient_groups(linear, 1)
    assert set(quadratic) == {0, 1, 2}
    assert len(quadratic[2]) == 1
    leading_value = next(iter(quadratic[2].values()))
    assert leading_value != 0 and next(iter(quadratic[2])) == zero_exponent()
    u_polynomial, v_polynomial = rational_quadratic_remainder_numerator(
        constant,
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
    )
    verify_rational_quadratic_remainder(
        constant,
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
        u_polynomial,
        v_polynomial,
    )
    norm = multiply_rational_polynomials(
        quadratic[2], multiply_rational_polynomials(u_polynomial, u_polynomial)
    )
    add_rational_polynomial(
        norm,
        multiply_rational_polynomials(
            quadratic[1],
            multiply_rational_polynomials(u_polynomial, v_polynomial),
        ),
        Fraction(-1),
    )
    add_rational_polynomial(
        norm,
        multiply_rational_polynomials(
            quadratic[0],
            multiply_rational_polynomials(v_polynomial, v_polynomial),
        ),
    )
    # The a3 elimination is localized at P1.  Its ratio-chart form is
    # 3*lambda+3+4*q, and denominator clearing contributes its square to
    # this resultant.  Remove exactly that invertible factor.
    first_pivot_ratio: RationalParameterPolynomial = {
        (1, 0, 0, 0, 0, 0, 0): Fraction(3),
        (0, 0, 0, 0, 0, 1, 0): Fraction(4),
        zero_exponent(): Fraction(3),
    }
    norm, pivot_exponent = strip_rational_factor(norm, first_pivot_ratio)
    assert pivot_exponent == 2
    mu4_constant, mu4_linear = rational_quadratic_remainder_numerator(
        transformed[4],
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
    )
    verify_rational_quadratic_remainder(
        transformed[4],
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
        mu4_constant,
        mu4_linear,
    )
    descended_mu4 = multiply_rational_polynomials(
        v_polynomial, mu4_constant
    )
    add_rational_polynomial(
        descended_mu4,
        multiply_rational_polynomials(u_polynomial, mu4_linear),
        Fraction(-1),
    )
    descended_mu4, mu4_pivot_exponent = strip_rational_factor(
        descended_mu4, first_pivot_ratio
    )
    assert mu4_pivot_exponent == 3
    mu5_constant, mu5_linear = rational_quadratic_remainder_numerator(
        transformed[5],
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
    )
    verify_rational_quadratic_remainder(
        transformed[5],
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
        mu5_constant,
        mu5_linear,
    )
    descended_mu5 = multiply_rational_polynomials(
        v_polynomial, mu5_constant
    )
    add_rational_polynomial(
        descended_mu5,
        multiply_rational_polynomials(u_polynomial, mu5_linear),
        Fraction(-1),
    )
    descended_mu5, mu5_pivot_exponent = strip_rational_factor(
        descended_mu5, first_pivot_ratio
    )
    assert mu5_pivot_exponent == 4
    mu4_by_z = rational_coefficient_groups(descended_mu4, 2)
    mu5_by_z = rational_coefficient_groups(descended_mu5, 2)
    assert max(mu4_by_z) == max(mu5_by_z) == 3
    subresultant = multiply_rational_polynomials(
        mu4_by_z[3], descended_mu5
    )
    add_rational_polynomial(
        subresultant,
        multiply_rational_polynomials(mu5_by_z[3], descended_mu4),
        Fraction(-1),
    )
    subresultant_by_z = rational_coefficient_groups(subresultant, 2)
    assert set(subresultant_by_z) == {0, 1, 2}
    later_descended: dict[int, RationalParameterPolynomial] = {}
    later_removed_powers: dict[str, int] = {}
    for order, expected_power in ((6, 4), (7, 5)):
        later_constant, later_linear = (
            rational_quadratic_remainder_numerator(
                transformed[order],
                1,
                (quadratic[2], quadratic[1], quadratic[0]),
            )
        )
        verify_rational_quadratic_remainder(
            transformed[order],
            1,
            (quadratic[2], quadratic[1], quadratic[0]),
            later_constant,
            later_linear,
        )
        descended = multiply_rational_polynomials(
            v_polynomial, later_constant
        )
        add_rational_polynomial(
            descended,
            multiply_rational_polynomials(u_polynomial, later_linear),
            Fraction(-1),
        )
        descended, removed_power = strip_rational_factor(
            descended, first_pivot_ratio
        )
        assert removed_power == expected_power, (
            order,
            removed_power,
            expected_power,
        )
        later_descended[order] = descended
        later_removed_powers[str(order)] = removed_power
    return {
        "quadratic_lead": quadratic_lead,
        "projection": primitive_integer_coefficients(projection),
        "norm": primitive_integer_coefficients(norm),
        "descended_mu4": primitive_integer_coefficients(descended_mu4),
        "mu4_z_lead": primitive_integer_coefficients(mu4_by_z[3]),
        "mu5_subresultant": primitive_integer_coefficients(subresultant),
        "subresultant_z_lead": primitive_integer_coefficients(
            subresultant_by_z[2]
        ),
        "descended_mu6": primitive_integer_coefficients(
            later_descended[6]
        ),
        "descended_mu7": primitive_integer_coefficients(
            later_descended[7]
        ),
        "incidence_terms": len(transformed[3]),
        "quadratic_terms": {
            str(power): len(polynomial)
            for power, polynomial in quadratic.items()
        },
        "remainder_terms": {"U": len(u_polynomial), "V": len(v_polynomial)},
        "norm_removed_P1_power": pivot_exponent,
        "mu4_removed_P1_power": mu4_pivot_exponent,
        "mu5_removed_P1_power": mu5_pivot_exponent,
        "later_removed_P1_powers": later_removed_powers,
        "exact_quadratic_remainder_identities": True,
    }


def modular_descent(
    prime: int, later_orders: tuple[int, ...] = (6, 7)
) -> dict[str, object]:
    assert later_orders and min(later_orders) >= 6
    orders = tuple(range(2, max(later_orders) + 1))
    b_base, d_base = base_polynomials(prime)
    b_powers = powers(b_base, max(orders), prime)
    d_powers = powers(d_base, max(orders), prime)
    moments = {
        order: moment(order, b_powers, d_powers, prime)
        for order in orders
    }
    pivot, rest = split_linear_parameter_polynomial(moments[2], 3, prime)
    eliminated = {
        order: linear_substitution_numerator(
            moments[order], 3, pivot, rest, prime
        )
        for order in orders
        if order > 2
    }
    b1_replacement, b2_replacement = finite_ratio_replacements(prime)

    def transform(polynomial: dict[tuple[int, ...], int]):
        answer = substitute_parameter_variable(
            polynomial, 5, b1_replacement, prime
        )
        answer = substitute_parameter_variable(
            answer, 6, b2_replacement, prime
        )
        return rescale_a2_by_b0(answer, prime)[0]

    transformed = {order: transform(polynomial) for order, polynomial in eliminated.items()}
    incidence_groups = coefficient_groups(transformed[3], 2, prime)
    assert set(incidence_groups) == {0, 1}
    constant, linear = incidence_groups[0], incidence_groups[1]
    projection = linear_substitution_numerator(
        transformed[4], 2, linear, constant, prime
    )
    quadratic = coefficient_groups(linear, 1, prime)
    assert set(quadratic) == {0, 1, 2}
    u_polynomial, v_polynomial = modular_quadratic_remainder_numerator(
        constant,
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
        prime,
    )
    norm = multiply_parameter_polynomials(
        quadratic[2],
        multiply_parameter_polynomials(u_polynomial, u_polynomial, prime),
        prime,
    )
    add_parameter_polynomial(
        norm,
        multiply_parameter_polynomials(
            quadratic[1],
            multiply_parameter_polynomials(u_polynomial, v_polynomial, prime),
            prime,
        ),
        -1,
        prime,
    )
    add_parameter_polynomial(
        norm,
        multiply_parameter_polynomials(
            quadratic[0],
            multiply_parameter_polynomials(v_polynomial, v_polynomial, prime),
            prime,
        ),
        1,
        prime,
    )
    first_pivot_ratio = {
        (1, 0, 0, 0, 0, 0, 0): 3 % prime,
        (0, 0, 0, 0, 0, 1, 0): 4 % prime,
        zero_exponent(): 3 % prime,
    }
    norm, pivot_exponent = strip_parameter_factor(
        norm, first_pivot_ratio, prime
    )
    assert pivot_exponent == 2
    mu4_constant, mu4_linear = modular_quadratic_remainder_numerator(
        transformed[4],
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
        prime,
    )
    descended_mu4 = multiply_parameter_polynomials(
        v_polynomial, mu4_constant, prime
    )
    add_parameter_polynomial(
        descended_mu4,
        multiply_parameter_polynomials(
            u_polynomial, mu4_linear, prime
        ),
        -1,
        prime,
    )
    descended_mu4, mu4_pivot_exponent = strip_parameter_factor(
        descended_mu4, first_pivot_ratio, prime
    )
    assert mu4_pivot_exponent == 3
    mu5_constant, mu5_linear = modular_quadratic_remainder_numerator(
        transformed[5],
        1,
        (quadratic[2], quadratic[1], quadratic[0]),
        prime,
    )
    descended_mu5 = multiply_parameter_polynomials(
        v_polynomial, mu5_constant, prime
    )
    add_parameter_polynomial(
        descended_mu5,
        multiply_parameter_polynomials(
            u_polynomial, mu5_linear, prime
        ),
        -1,
        prime,
    )
    descended_mu5, mu5_pivot_exponent = strip_parameter_factor(
        descended_mu5, first_pivot_ratio, prime
    )
    assert mu5_pivot_exponent == 4
    mu4_by_z = coefficient_groups(descended_mu4, 2, prime)
    mu5_by_z = coefficient_groups(descended_mu5, 2, prime)
    assert max(mu4_by_z) == max(mu5_by_z) == 3
    subresultant = multiply_parameter_polynomials(
        mu4_by_z[3], descended_mu5, prime
    )
    add_parameter_polynomial(
        subresultant,
        multiply_parameter_polynomials(
            mu5_by_z[3], descended_mu4, prime
        ),
        -1,
        prime,
    )
    subresultant_by_z = coefficient_groups(subresultant, 2, prime)
    assert set(subresultant_by_z) == {0, 1, 2}
    later_descended: dict[int, dict[tuple[int, ...], int]] = {}
    later_removed_powers: dict[str, int] = {}
    for order in later_orders:
        expected_power = {6: 4, 7: 5}.get(order)
        later_constant, later_linear = (
            modular_quadratic_remainder_numerator(
                transformed[order],
                1,
                (quadratic[2], quadratic[1], quadratic[0]),
                prime,
            )
        )
        descended = multiply_parameter_polynomials(
            v_polynomial, later_constant, prime
        )
        add_parameter_polynomial(
            descended,
            multiply_parameter_polynomials(
                u_polynomial, later_linear, prime
            ),
            -1,
            prime,
        )
        descended, removed_power = strip_parameter_factor(
            descended, first_pivot_ratio, prime
        )
        if expected_power is not None:
            assert removed_power == expected_power, (
                order,
                removed_power,
                expected_power,
            )
        later_descended[order] = descended
        later_removed_powers[str(order)] = removed_power
    answer = {
        "projection": projection,
        "norm": norm,
        "degree_drop_U": u_polynomial,
        "degree_drop_V": v_polynomial,
        "degree_drop_quadratic": quadratic,
        "descended_mu4": descended_mu4,
        "mu4_z_lead": mu4_by_z[3],
        "mu5_subresultant": subresultant,
        "subresultant_z_lead": subresultant_by_z[2],
    }
    answer.update(
        {
            f"descended_mu{order}": later_descended[order]
            for order in later_orders
        }
    )
    answer["later_removed_P1_powers"] = later_removed_powers
    return answer


def compare_primitive_reduction(
    exact: dict[tuple[int, ...], int],
    modular: dict[tuple[int, ...], int],
    prime: int,
) -> int:
    reduction = {
        exponent: coefficient % prime
        for exponent, coefficient in exact.items()
        if coefficient % prime
    }
    assert set(reduction) == set(modular)
    lead = max(reduction)
    scalar = modular[lead] * pow(reduction[lead], -1, prime) % prime
    assert {
        exponent: coefficient * scalar % prime
        for exponent, coefficient in reduction.items()
    } == modular
    assert polynomial_degree(exact) == polynomial_degree(modular)
    return scalar


def singular_factorization(
    polynomial: dict[tuple[int, ...], int],
    variables: tuple[str, ...],
    active: tuple[int, ...],
    label: str,
) -> dict[str, int]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required")
    serializer_variables = ("r", *variables)
    expression = active_parameter_polynomial_string(
        polynomial, PRIME, serializer_variables, active
    )
    source = "\n".join(
        [
            f"ring R={PRIME},({','.join(variables)}),dp;",
            f"poly F={expression};",
            "ideal factors=factorize(F,1);",
            f'print("{label}_BEGIN");',
            'print("factor_count="+string(size(factors)));',
            'print("polynomial_degree="+string(deg(F)));',
            'print("factor_degree="+string(deg(factors[1])));',
            f'print("{label}_END");',
            "quit;",
        ]
    )
    completed = subprocess.run(
        [executable, "-q"],
        input=source,
        text=True,
        capture_output=True,
        timeout=300,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            "Singular factorization failed:\n"
            + completed.stdout[-2000:]
            + completed.stderr[-2000:]
        )
    match = re.search(
        label
        + r"_BEGIN\s+factor_count=(\d+)\s+polynomial_degree=(\d+)\s+"
        + r"factor_degree=(\d+)\s+"
        + label
        + r"_END",
        completed.stdout,
    )
    if match is None:
        raise AssertionError("could not parse Singular output: " + completed.stdout[-2000:])
    factor_count, degree, factor_degree = map(int, match.groups())
    assert factor_count == 1 and degree == factor_degree
    return {
        "factor_count": factor_count,
        "degree": degree,
        "factor_degree": factor_degree,
        "source_sha256": sha256(source.encode()).hexdigest(),
    }


def profile(polynomial: dict[tuple[int, ...], int], variables: tuple[int, ...]) -> dict[str, object]:
    return {
        "terms": len(polynomial),
        "total_degree": polynomial_degree(polynomial),
        "variable_degrees": variable_degrees(polynomial, variables),
    }


def main() -> None:
    exact = exact_descent()
    modular = modular_descent(PRIME)
    projection = exact["projection"]
    norm = exact["norm"]
    assert isinstance(projection, dict) and isinstance(norm, dict)
    projection_scalar = compare_primitive_reduction(
        projection, modular["projection"], PRIME
    )
    norm_scalar = compare_primitive_reduction(norm, modular["norm"], PRIME)
    descended_mu4 = exact["descended_mu4"]
    assert isinstance(descended_mu4, dict)
    mu4_scalar = compare_primitive_reduction(
        descended_mu4, modular["descended_mu4"], PRIME
    )
    mu4_z_lead = exact["mu4_z_lead"]
    mu5_subresultant = exact["mu5_subresultant"]
    subresultant_z_lead = exact["subresultant_z_lead"]
    assert isinstance(mu4_z_lead, dict)
    assert isinstance(mu5_subresultant, dict)
    assert isinstance(subresultant_z_lead, dict)
    mu4_z_lead_scalar = compare_primitive_reduction(
        mu4_z_lead, modular["mu4_z_lead"], PRIME
    )
    subresultant_scalar = compare_primitive_reduction(
        mu5_subresultant, modular["mu5_subresultant"], PRIME
    )
    subresultant_z_lead_scalar = compare_primitive_reduction(
        subresultant_z_lead,
        modular["subresultant_z_lead"],
        PRIME,
    )
    exact_context = fmpz_mpoly_ctx.get(["b0", "q", "lambda"])
    modular_context = nmod_mpoly_ctx.get(
        ["b0", "q", "lambda"], PRIME
    )
    exact_linear_norm = flint_linear_norm_descent(
        descended_mu4, mu5_subresultant, exact_context
    )
    modular_linear_norm = flint_linear_norm_descent(
        modular["descended_mu4"],
        modular["mu5_subresultant"],
        modular_context,
    )
    exact_norm_residual = exact_linear_norm["residual"]
    assert int(exact_norm_residual.content()) == 1
    remainder_constant_scalar = compare_flint_reduction(
        exact_linear_norm["remainder_constant"],
        modular_linear_norm["remainder_constant"],
        PRIME,
    )
    remainder_linear_scalar = compare_flint_reduction(
        exact_linear_norm["remainder_linear"],
        modular_linear_norm["remainder_linear"],
        PRIME,
    )
    assert remainder_constant_scalar == remainder_linear_scalar
    linear_norm_residual_scalar = compare_flint_reduction(
        exact_norm_residual,
        modular_linear_norm["residual"],
        PRIME,
    )
    modular_norm_residual = full_parameter_polynomial_from_flint(
        modular_linear_norm["residual"]
    )
    later_base_equations = {}
    modular_later_residuals = {}
    for order, b0_power in ((6, 1), (7, 2)):
        exact_later = flint_later_base_equation(
            exact[f"descended_mu{order}"],
            descended_mu4,
            mu5_subresultant,
            exact_linear_norm,
            exact_context,
            b0_power,
        )
        modular_later = flint_later_base_equation(
            modular[f"descended_mu{order}"],
            modular["descended_mu4"],
            modular["mu5_subresultant"],
            modular_linear_norm,
            modular_context,
            b0_power,
        )
        module_constant_scalar = compare_flint_reduction(
            exact_later["module_constant"],
            modular_later["module_constant"],
            PRIME,
        )
        module_linear_scalar = compare_flint_reduction(
            exact_later["module_linear"],
            modular_later["module_linear"],
            PRIME,
        )
        assert module_constant_scalar == module_linear_scalar
        determinant_scalar = compare_flint_reduction(
            exact_later["determinant"],
            modular_later["determinant"],
            PRIME,
        )
        residual_content = int(exact_later["residual"].content())
        assert residual_content > 0
        primitive_residual = exact_later["residual"] // residual_content
        residual_scalar = compare_flint_reduction(
            primitive_residual,
            modular_later["residual"],
            PRIME,
        )
        modular_later_residuals[order] = full_parameter_polynomial_from_flint(
            modular_later["residual"]
        )
        later_base_equations[order] = {
            "exact": exact_later,
            "modular": modular_later,
            "primitive_residual": primitive_residual,
            "module_reduction_scalar": module_constant_scalar,
            "determinant_reduction_scalar": determinant_scalar,
            "residual_content": residual_content,
            "residual_reduction_scalar": residual_scalar,
        }
    projection_factorization = singular_factorization(
        modular["projection"],
        ("a1", "b0", "q", "lambda"),
        (1, 4, 5, 0),
        "PROJECTION",
    )
    norm_factorization = singular_factorization(
        modular["norm"],
        ("b0", "q", "lambda"),
        (4, 5, 0),
        "NORM",
    )
    mu4_factorization = singular_factorization(
        modular["descended_mu4"],
        ("z", "b0", "q", "lambda"),
        (2, 4, 5, 0),
        "DESCENDED_MU4",
    )
    mu4_z_lead_factorization = singular_factorization(
        modular["mu4_z_lead"],
        ("b0", "q", "lambda"),
        (4, 5, 0),
        "MU4_Z_LEAD",
    )
    subresultant_factorization = singular_factorization(
        modular["mu5_subresultant"],
        ("z", "b0", "q", "lambda"),
        (2, 4, 5, 0),
        "MU5_SUBRESULTANT",
    )
    subresultant_z_lead_factorization = singular_factorization(
        modular["subresultant_z_lead"],
        ("b0", "q", "lambda"),
        (4, 5, 0),
        "SUBRESULTANT_Z_LEAD",
    )
    linear_norm_residual_factorization = singular_factorization(
        modular_norm_residual,
        ("b0", "q", "lambda"),
        (4, 5, 0),
        "LINEAR_NORM_RESIDUAL",
    )
    later_base_factorizations = {
        order: singular_factorization(
            modular_later_residuals[order],
            ("b0", "q", "lambda"),
            (4, 5, 0),
            f"MU{order}_BASE_RESIDUAL",
        )
        for order in (6, 7)
    }

    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-hurwitz-module-descent-v1",
        "field": "characteristic zero",
        "chart": "A=0 on the M01-open cubic Hurwitz chart",
        "coordinates": "q=b1/b0, z=b0*a2, with b2 eliminated by A",
        "exact_b0_open_proof": {
            "A": (
                "3*b0^2*lambda^2+b0^2*lambda+6*b0^2+"
                "b0*b1*lambda+5*b0*b1-8*b0*b2-8*b0+4*b1^2"
            ),
            "M01": "b1-a1*b0",
            "argument": "at b0=0, A=4*b1^2 while M01=b1",
        },
        "prime": PRIME,
        "prime_guard": "29>3*7, so no factorial in moments through mu_7 vanishes",
        "incidence_terms": exact["incidence_terms"],
        "open_projection": {
            "construction": "numerator of mu4 after z=-C/B",
            "profile": profile(projection, (0, 1, 4, 5)),
            "reduction_scalar": projection_scalar,
            "factorization": projection_factorization,
            "conclusion": (
                "the B-open incidence-plus-mu4 projection is irreducible over QQ"
            ),
        },
        "degree_drop_norm": {
            "quadratic_terms_by_a1_power": exact["quadratic_terms"],
            "remainder_terms": exact["remainder_terms"],
            "exact_remainder_identities": exact[
                "exact_quadratic_remainder_identities"
            ],
            "removed_localizer_factor": {
                "factor": "3*lambda+3+4*q=P1/(3*b0)",
                "power": exact["norm_removed_P1_power"],
            },
            "construction": "A2*U^2-A1*U*V+A0*V^2",
            "profile": profile(norm, (0, 4, 5)),
            "reduction_scalar": norm_scalar,
            "factorization": norm_factorization,
            "conclusion": (
                "B=C=0 has one irreducible characteristic-zero component; "
                "its V-open is birational to the norm hypersurface"
            ),
            "first_later_moment": {
                "equation": "V*(mu4 mod B)_0-U*(mu4 mod B)_1",
                "removed_P1_power": exact["mu4_removed_P1_power"],
                "profile": profile(descended_mu4, (0, 2, 4, 5)),
                "reduction_scalar": mu4_scalar,
                "factorization": mu4_factorization,
                "conclusion": (
                    "the descended mu4 equation is irreducible over QQ"
                ),
                "z_leading_coefficient": {
                    "profile": profile(mu4_z_lead, (0, 4, 5)),
                    "reduction_scalar": mu4_z_lead_scalar,
                    "factorization": mu4_z_lead_factorization,
                },
            },
            "second_later_moment": {
                "removed_mu5_P1_power": exact["mu5_removed_P1_power"],
                "construction": (
                    "lc_z(N4)*N5-lc_z(N5)*N4; quadratic in z"
                ),
                "profile": profile(mu5_subresultant, (0, 2, 4, 5)),
                "reduction_scalar": subresultant_scalar,
                "factorization": subresultant_factorization,
                "z_leading_coefficient": {
                    "profile": profile(subresultant_z_lead, (0, 4, 5)),
                    "reduction_scalar": subresultant_z_lead_scalar,
                    "factorization": subresultant_z_lead_factorization,
                },
                "conclusion": (
                    "on lc_z(N4)!=0, mu5 descends to an irreducible "
                    "quadratic-in-z subresultant; both displayed leading "
                    "coefficient divisors are irreducible over QQ"
                ),
                "linear_remainder": {
                    "construction": (
                        "pseudo-remainder of cubic N4 modulo quadratic S5"
                    ),
                    "constant_profile": flint_profile(
                        exact_linear_norm["remainder_constant"]
                    ),
                    "linear_profile": flint_profile(
                        exact_linear_norm["remainder_linear"]
                    ),
                    "common_reduction_scalar": remainder_constant_scalar,
                },
                "linear_norm": {
                    "construction": "S2*R0^2-S1*R0*R1+S0*R1^2",
                    "raw_profile": flint_profile(
                        exact_linear_norm["norm"]
                    ),
                    "removed_factors": {
                        "b0": 2,
                        "lc_z_N4": 2,
                        "lc_z_S5": 2,
                    },
                    "residual_profile": flint_profile(
                        exact_norm_residual
                    ),
                    "modular_residual_terms": len(
                        modular_linear_norm["residual"].to_dict()
                    ),
                    "reduction_scalar": linear_norm_residual_scalar,
                    "factorization": linear_norm_residual_factorization,
                    "conclusion": (
                        "on b0*lc_z(N4)*lc_z(S5)!=0, the common-root "
                        "projection of N4 and N5 is the irreducible "
                        "degree-118 hypersurface K=0 over QQ"
                    ),
                },
                "later_base_equations": {
                    str(order): {
                        "construction": (
                            "R1*(mu_order mod S5)_0-"
                            "R0*(mu_order mod S5)_1"
                        ),
                        "removed_P1_power_before_z_reduction": exact[
                            "later_removed_P1_powers"
                        ][str(order)],
                        "module_constant_profile": flint_profile(
                            later_base_equations[order]["exact"][
                                "module_constant"
                            ]
                        ),
                        "module_linear_profile": flint_profile(
                            later_base_equations[order]["exact"][
                                "module_linear"
                            ]
                        ),
                        "common_module_reduction_scalar": (
                            later_base_equations[order][
                                "module_reduction_scalar"
                            ]
                        ),
                        "raw_determinant_profile": flint_profile(
                            later_base_equations[order]["exact"][
                                "determinant"
                            ]
                        ),
                        "residual_content_before_primitive_normalization": (
                            later_base_equations[order][
                                "residual_content"
                            ]
                        ),
                        "raw_determinant_reduction_scalar": (
                            later_base_equations[order][
                                "determinant_reduction_scalar"
                            ]
                        ),
                        "removed_factors": {
                            "b0": 1 if order == 6 else 2,
                            "lc_z_N4": 1,
                            "lc_z_S5": 2,
                        },
                        "residual_profile": flint_profile(
                            later_base_equations[order][
                                "primitive_residual"
                            ]
                        ),
                        "modular_residual_terms": len(
                            later_base_equations[order]["modular"][
                                "residual"
                            ].to_dict()
                        ),
                        "residual_reduction_scalar": (
                            later_base_equations[order][
                                "residual_reduction_scalar"
                            ]
                        ),
                        "factorization": later_base_factorizations[order],
                        "conclusion": (
                            "the primitive base equation is irreducible over "
                            "QQ; on V*R1!=0 it is equivalent to imposing "
                            "this later moment after N4=S5=0, while V=0 "
                            "and R0=R1=0 remain separate boundaries"
                        ),
                    }
                    for order in (6, 7)
                },
            },
        },
        "logic": (
            "each primitive exact polynomial has degree-preserving reduction "
            "equal up to a nonzero scalar to its independently constructed "
            "GF(29) polynomial; irreducibility modulo 29 implies "
            "irreducibility over QQ"
        ),
        "scope": (
            "component and subresultant structure through mu7, including exact "
            "irreducible base equations K, J6, and J7; their common zero set, "
            "the V=0 and R0=R1=0 boundaries, and the localization boundaries "
            "P1=0, Delta=0, and M01=0 remain separate"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print("PASS exact ratio chart and b0-open implication")
    print("PASS irreducible B-open mu4 projection over QQ")
    print("PASS irreducible B=C degree-drop component over QQ")
    print("PASS irreducible descended mu4 equation over QQ")
    print("PASS irreducible cubic-cubic mu5 subresultant over QQ")
    print("PASS irreducible degree-118 linear-remainder norm over QQ")
    print("PASS irreducible degree-130 mu6 base equation over QQ")
    print("PASS irreducible degree-162 mu7 base equation over QQ")


if __name__ == "__main__":
    main()
