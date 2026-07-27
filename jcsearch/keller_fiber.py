"""Public quadratic-gauge compiler for prescribed finite étale fibers."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

import sympy as sp


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


def quadratic_gauge_map(
    seed: sp.Expr,
    inverse_variable: sp.Symbol,
    source_variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the determinant-minus-two quadratic gauge attached to ``seed``."""
    if len(source_variables) != 3:
        raise ValueError("the quadratic-gauge compiler needs three source variables")
    x, y, z = source_variables
    polynomial = sp.Poly(sp.expand(seed), inverse_variable, domain=sp.QQ)
    degree = polynomial.degree()
    coefficients = {
        index: polynomial.coeff_monomial(inverse_variable**index)
        for index in range(1, degree + 1)
    }
    g1 = coefficients.get(1, 0)
    g3 = coefficients.get(3, 0)
    if degree < 3 or not g1 or not g3 or not coefficients.get(degree, 0):
        raise ValueError("seed must have degree at least three and g_1 g_3 g_N nonzero")

    t = 1 + x * y
    q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
    pi = t * q
    second = (
        y
        + 3 * (g3 / g1) * x * q
        + 2 * (coefficients.get(2, 0) / g1) * t * q
    )
    third = x * (5 - 3 * t) - (g3 / g1) * x**3 * z
    for index in range(4, degree + 1):
        second += (
            index
            * (coefficients[index] / g1)
            * t**2
            * x ** (index - 2)
            * q**index
        )
        third -= (
            (index - 2)
            * (coefficients[index] / g1)
            * (x * q) ** index
        )
    return tuple(sp.cancel(component) for component in (pi, second, third))


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


def choose_admissible_translation(polynomial: sp.Poly) -> sp.Rational:
    """Choose the first integer ``a`` with ``P'(a)P'''(a) != 0``."""
    variable = polynomial.gens[0]
    first = sp.diff(polynomial.as_expr(), variable)
    third = sp.diff(polynomial.as_expr(), variable, 3)
    radius = 0
    while True:
        candidates = (0,) if radius == 0 else (radius, -radius)
        for candidate in candidates:
            if first.subs(variable, candidate) and third.subs(variable, candidate):
                return sp.Rational(candidate)
        radius += 1


def compile_polynomial_to_keller_fiber(
    polynomial: sp.Poly | sp.Expr,
    variable: sp.Symbol,
    *,
    translation: int | Fraction | sp.Rational | None = None,
    inverse_variable: sp.Symbol | None = None,
    source_variables: Sequence[sp.Symbol] | None = None,
) -> KellerFiberCompilation:
    """Compile a squarefree rational polynomial into a complete Keller fiber.

    The returned map is the determinant-one normalization of the
    root-engineered quadratic gauge.  The theorem identifies its fiber over
    ``target`` with ``Spec(Q[T]/(P))`` after the displayed translation.
    """
    rational_polynomial = sp.Poly(polynomial, variable, domain=sp.QQ)
    degree = rational_polynomial.degree()
    if degree < 3:
        raise ValueError("complete quadratic-gauge fibers require degree at least three")
    if rational_polynomial.gcd(rational_polynomial.diff()).degree() != 0:
        raise ValueError("the prescribed polynomial must be squarefree")

    a = (
        choose_admissible_translation(rational_polynomial)
        if translation is None
        else _as_rational(translation)
    )
    first_at_a = sp.diff(rational_polynomial.as_expr(), variable).subs(variable, a)
    third_at_a = sp.diff(
        rational_polynomial.as_expr(), variable, 3
    ).subs(variable, a)
    if not first_at_a or not third_at_a:
        raise ValueError("translation must satisfy P'(a)P'''(a) != 0")

    inverse_symbol = inverse_variable or sp.Symbol("S")
    sources = tuple(source_variables or sp.symbols("x y z"))
    if len(sources) != 3 or not all(isinstance(symbol, sp.Symbol) for symbol in sources):
        raise ValueError("source_variables must contain three symbols")

    translated = sp.expand(
        rational_polynomial.as_expr().subs(variable, a + inverse_symbol)
    )
    value_at_a = rational_polynomial.eval(a)
    seed = sp.expand(translated - value_at_a)
    minus_two = quadratic_gauge_map(seed, inverse_symbol, sources)
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
        determinant_minus_two_map=minus_two,
        determinant_one_map=determinant_one,
        target=target,
        inverse_polynomial=inverse_polynomial,
        geometric_degree=degree,
        coordinate_degrees=total_coordinate_degrees(determinant_one, sources),
    )
