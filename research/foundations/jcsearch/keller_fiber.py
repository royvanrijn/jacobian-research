"""Public quadratic-gauge compiler for prescribed finite étale fibers."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Sequence

import sympy as sp


@dataclass(frozen=True)
class StableMultiplicityCertificate:
    """Exact invariant data supplied by the proved multiplicity families.

    This is a reproducible parameter-and-value record; the stable
    functoriality that makes the value separating remains a written theorem.
    """

    family_parameter: int
    gauge_exponent: int
    separation_invariant: str
    separation_value: int
    fitting_support: tuple[tuple[int, int], ...] = ()
    boundary_prime_count: int | None = None
    boundary_ramification_index: int | None = None


@dataclass(frozen=True)
class KellerFiberCompilation:
    """Exact data for one polynomial compiled into a complete Keller fiber."""

    polynomial: sp.Poly
    translation: sp.Rational
    inverse_variable: sp.Symbol
    source_variables: tuple[sp.Symbol, sp.Symbol, sp.Symbol]
    seed: sp.Expr
    determinant_minus_two_map: tuple[sp.Expr, sp.Expr, sp.Expr]
    determinant_one_map: tuple[sp.Expr, sp.Expr, sp.Expr]
    target: tuple[sp.Expr, sp.Expr, sp.Expr]
    inverse_polynomial: sp.Expr
    geometric_degree: int
    coordinate_degrees: tuple[int, int, int]
    lifted_seed: sp.Expr | None = None
    stable_multiplicity: StableMultiplicityCertificate | None = None


def _seed_coefficients(
    seed: sp.Expr,
    inverse_variable: sp.Symbol,
) -> tuple[int, dict[int, sp.Expr]]:
    polynomial = sp.Poly(sp.expand(seed), inverse_variable, domain=sp.QQ)
    degree = polynomial.degree()
    coefficients = {
        index: polynomial.coeff_monomial(inverse_variable**index)
        for index in range(1, degree + 1)
    }
    g1 = coefficients.get(1, 0)
    g3 = coefficients.get(3, 0)
    if degree < 3 or not g1 or not g3 or not coefficients.get(degree, 0):
        raise ValueError(
            "seed must have degree at least three and g_1 g_3 g_N nonzero"
        )
    return degree, coefficients


def _quadratic_gauge_data(
    seed: sp.Expr,
    inverse_variable: sp.Symbol,
    source_variables: Sequence[sp.Symbol],
    stable_parameter: int | None,
) -> tuple[tuple[sp.Expr, sp.Expr, sp.Expr], sp.Expr]:
    if len(source_variables) != 3:
        raise ValueError("the quadratic-gauge compiler needs three source variables")
    x, y, z = source_variables
    if not all(isinstance(symbol, sp.Symbol) for symbol in (x, y, z)):
        raise ValueError("source_variables must contain three symbols")
    if len({x, y, z}) != 3 or inverse_variable in {x, y, z}:
        raise ValueError("source and inverse variables must be distinct")
    degree, coefficients = _seed_coefficients(seed, inverse_variable)
    g1 = coefficients[1]
    g2 = coefficients.get(2, 0)
    g3 = coefficients[3]

    t = 1 + x * y
    q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
    pi = t * q
    second = (
        y
        + 3 * (g3 / g1) * x * q
        + 2 * (g2 / g1) * t * q
    )
    third = x * (5 - 3 * t) - (g3 / g1) * x**3 * z
    if stable_parameter is not None and degree == 3:
        exponent = stable_parameter + 4
        second += 3 * (g3 / g1) * (
            t ** (exponent - 1) * x * q**exponent - t**2 * x * q**3
        )
        third -= (g3 / g1) * (
            t ** (exponent - 3) * x**3 * q**exponent - x**3 * q**3
        )
        lifted_seed = (
            g1 * inverse_variable
            + g2 * pi * inverse_variable**2
            + g3
            * pi
            * (1 + pi ** (exponent - 1) - pi**2)
            * inverse_variable**3
        )
    else:
        shift = stable_parameter or 0
        for index in range(4, degree + 1):
            second += (
                index
                * (coefficients[index] / g1)
                * t ** (shift + 2)
                * x ** (index - 2)
                * q ** (index + shift)
            )
            third -= (
                (index - 2)
                * (coefficients[index] / g1)
                * t**shift
                * x**index
                * q ** (index + shift)
            )
        lifted_seed = (
            g1 * inverse_variable
            + pi * (g2 * inverse_variable**2 + g3 * inverse_variable**3)
            + sum(
                coefficients[index]
                * pi ** (index + shift)
                * inverse_variable**index
                for index in range(4, degree + 1)
            )
        )
    mapping = tuple(sp.cancel(component) for component in (pi, second, third))
    return mapping, lifted_seed


def quadratic_gauge_map(
    seed: sp.Expr,
    inverse_variable: sp.Symbol,
    source_variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the determinant-minus-two minimal quadratic gauge for ``seed``."""
    mapping, _ = _quadratic_gauge_data(
        seed, inverse_variable, source_variables, stable_parameter=None
    )
    return mapping


def stable_multiplicity_gauge_map(
    seed: sp.Expr,
    inverse_variable: sp.Symbol,
    source_variables: Sequence[sp.Symbol],
    family_parameter: int,
) -> tuple[tuple[sp.Expr, sp.Expr, sp.Expr], sp.Expr]:
    """Return a fiber-invisible multiplicity lift and its deformed seed.

    In degree three, ``family_parameter=k`` selects the cubic lift with
    exponent ``n=k+4``.  In every degree at least four it is the common
    power shift ``m=k``.
    """
    parameter = _as_nonnegative_integer(family_parameter)
    degree, coefficients = _seed_coefficients(seed, inverse_variable)
    if degree >= 4 and any(
        not coefficients.get(index, 0) for index in range(4, degree + 1)
    ):
        raise ValueError(
            "stable power shifts require every coefficient g_4,...,g_N nonzero"
        )
    return _quadratic_gauge_data(
        seed, inverse_variable, source_variables, stable_parameter=parameter
    )


def jacobian_one_normalization(
    mapping: Sequence[sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Apply the target normalization ``diag(1,-1/2,1)``."""
    if len(mapping) != 3:
        raise ValueError("a Keller map must have three coordinates")
    return (mapping[0], sp.cancel(-mapping[1] / 2), mapping[2])


def total_coordinate_degrees(
    mapping: Sequence[sp.Expr],
    source_variables: Sequence[sp.Symbol],
) -> tuple[int, int, int]:
    """Return exact total degrees of three polynomial map coordinates."""
    if len(mapping) != 3 or len(source_variables) != 3:
        raise ValueError("expected three map coordinates and three source variables")
    return tuple(
        int(
            sp.Poly(
                sp.expand(component),
                *source_variables,
                domain=sp.QQ,
            ).total_degree()
        )
        for component in mapping
    )


def _as_rational(value: int | Fraction | sp.Rational) -> sp.Rational:
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    rational = sp.sympify(value)
    if not rational.is_Rational:
        raise ValueError("translation must be rational")
    return sp.Rational(rational)


def _as_nonnegative_integer(value: int | sp.Integer) -> int:
    if isinstance(value, bool):
        raise ValueError("stable_parameter must be a nonnegative integer")
    integer = sp.sympify(value)
    if not integer.is_Integer or integer < 0:
        raise ValueError("stable_parameter must be a nonnegative integer")
    return int(integer)


def choose_admissible_translation(
    polynomial: sp.Poly,
    derivative_orders: Sequence[int] = (1, 3),
) -> sp.Rational:
    """Choose the first integer where all requested derivatives are nonzero."""
    variable = polynomial.gens[0]
    normalized_orders = []
    for order in derivative_orders:
        if isinstance(order, bool):
            raise ValueError("derivative_orders must contain nonnegative integers")
        integer = sp.sympify(order)
        if not integer.is_Integer or integer < 0:
            raise ValueError("derivative_orders must contain nonnegative integers")
        normalized_orders.append(int(integer))
    orders = tuple(dict.fromkeys(normalized_orders))
    if not orders:
        raise ValueError("derivative_orders must contain nonnegative integers")
    derivatives = [
        sp.diff(polynomial.as_expr(), variable, order) for order in orders
    ]
    if any(derivative == 0 for derivative in derivatives):
        raise ValueError("a requested derivative vanishes identically")
    radius = 0
    while True:
        candidates = (0,) if radius == 0 else (radius, -radius)
        for candidate in candidates:
            if all(
                derivative.subs(variable, candidate) for derivative in derivatives
            ):
                return sp.Rational(candidate)
        radius += 1


def _stable_certificate(
    degree: int,
    family_parameter: int,
) -> StableMultiplicityCertificate:
    if degree == 3:
        exponent = family_parameter + 4
        return StableMultiplicityCertificate(
            family_parameter=family_parameter,
            gauge_exponent=exponent,
            separation_invariant="geometric_boundary_target_components",
            separation_value=exponent,
        )
    support = ((0, 0), (1, 2)) + tuple(
        (index + family_parameter, index - 1)
        for index in range(4, degree + 1)
    )
    common_divisor = gcd(degree - 3, family_parameter + 2)
    return StableMultiplicityCertificate(
        family_parameter=family_parameter,
        gauge_exponent=family_parameter,
        separation_invariant="normalized_fitting_newton_area",
        separation_value=2 * degree - 3 + (degree - 2) * family_parameter,
        fitting_support=support,
        boundary_prime_count=common_divisor,
        boundary_ramification_index=(degree - 3) // common_divisor,
    )


def compile_polynomial_to_keller_fiber(
    polynomial: sp.Poly | sp.Expr,
    variable: sp.Symbol,
    *,
    translation: int | Fraction | sp.Rational | None = None,
    inverse_variable: sp.Symbol | None = None,
    source_variables: Sequence[sp.Symbol] | None = None,
    stable_parameter: int | sp.Integer | None = None,
) -> KellerFiberCompilation:
    """Compile a squarefree rational polynomial into a complete Keller fiber.

    The returned map is the determinant-one normalization of the
    root-engineered quadratic gauge.  The theorem identifies its fiber over
    ``target`` with ``Spec(Q[T]/(P))`` after the displayed translation.

    Supplying ``stable_parameter=k>=0`` chooses a certified infinite-family
    lift without changing that fiber: cubic exponent ``n=k+4`` in degree
    three, and common power shift ``m=k`` in every degree at least four.
    """
    rational_polynomial = sp.Poly(polynomial, variable, domain=sp.QQ)
    degree = rational_polynomial.degree()
    if degree < 3:
        raise ValueError("complete quadratic-gauge fibers require degree at least three")
    if rational_polynomial.gcd(rational_polynomial.diff()).degree() != 0:
        raise ValueError("the prescribed polynomial must be squarefree")

    parameter = (
        None if stable_parameter is None else _as_nonnegative_integer(stable_parameter)
    )
    derivative_orders = (
        (1, 3)
        if parameter is None or degree == 3
        else (1, 3, *range(4, degree + 1))
    )
    a = (
        choose_admissible_translation(rational_polynomial, derivative_orders)
        if translation is None
        else _as_rational(translation)
    )
    derivative_values = {
        order: sp.diff(rational_polynomial.as_expr(), variable, order).subs(
            variable, a
        )
        for order in derivative_orders
    }
    if any(not value for value in derivative_values.values()):
        orders = ",".join(str(order) for order in derivative_orders)
        raise ValueError(
            f"translation must have nonzero polynomial derivatives of orders {orders}"
        )
    first_at_a = derivative_values[1]

    inverse_symbol = inverse_variable or sp.Symbol("S")
    sources = tuple(source_variables or sp.symbols("x y z"))
    if len(sources) != 3 or not all(
        isinstance(symbol, sp.Symbol) for symbol in sources
    ):
        raise ValueError("source_variables must contain three symbols")
    if len(set(sources)) != 3 or inverse_symbol in set(sources):
        raise ValueError("source and inverse variables must be distinct")

    translated = sp.expand(
        rational_polynomial.as_expr().subs(variable, a + inverse_symbol)
    )
    value_at_a = rational_polynomial.eval(a)
    seed = sp.expand(translated - value_at_a)
    if parameter is None:
        minus_two, lifted_seed = _quadratic_gauge_data(
            seed, inverse_symbol, sources, stable_parameter=None
        )
        stable_certificate = None
    else:
        minus_two, lifted_seed = stable_multiplicity_gauge_map(
            seed, inverse_symbol, sources, parameter
        )
        stable_certificate = _stable_certificate(degree, parameter)
    determinant_one = jacobian_one_normalization(minus_two)
    target = (
        sp.Integer(1),
        sp.Integer(0),
        sp.cancel(-2 * value_at_a / first_at_a),
    )
    seed_polynomial = sp.Poly(seed, inverse_symbol, domain=sp.QQ)
    g1 = seed_polynomial.coeff_monomial(inverse_symbol)
    inverse_polynomial = sp.expand(seed - g1 * target[2] / 2)
    assert sp.expand(inverse_polynomial - translated) == 0

    return KellerFiberCompilation(
        polynomial=rational_polynomial,
        translation=a,
        inverse_variable=inverse_symbol,
        source_variables=sources,
        seed=seed,
        lifted_seed=lifted_seed,
        determinant_minus_two_map=minus_two,
        determinant_one_map=determinant_one,
        target=target,
        inverse_polynomial=inverse_polynomial,
        geometric_degree=degree,
        coordinate_degrees=total_coordinate_degrees(determinant_one, sources),
        stable_multiplicity=stable_certificate,
    )
