#!/usr/bin/env python3
"""Exact finite-residual replays for a slice LND image.

Let D=d/ds on Q[x,s].  The complete intersection used below is the radical
ideal of six points arranged in three vertical fibers.  For h in Q[x,s],

    h in D(I)

is equivalent to the zero-constant s-primitive of h taking one common value
on every vertical fiber.  The script checks this equivalence on a generic
degree window, verifies the reduced complete-intersection certificate, and
replays pure and mixed powers for a polynomial which vanishes on precisely
the two nontrivial fibers.  It also checks two nonreduced length-two residual
schemes carried by q=s^2 and both branches of the carrier-free primary
slice-or-vanishing dichotomy.

The generic degree window and the displayed powers are regressions.  The
all-degree arguments are written in
extended-geometry/LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md.
"""

from __future__ import annotations

from collections.abc import Iterable

import sympy as sp


x, s = sp.symbols("x s")

# Two collision fibers and one singleton fiber.
fibers: dict[sp.Rational, tuple[sp.Rational, ...]] = {
    sp.Rational(0): (sp.Rational(0), sp.Rational(1), sp.Rational(2)),
    sp.Rational(2): (sp.Rational(-1), sp.Rational(3)),
    sp.Rational(5): (sp.Rational(4),),
}


def lagrange_idempotent(a: sp.Rational) -> sp.Expr:
    """Return the Lagrange idempotent for x=a on the three-point base."""

    numerator = sp.prod(x - other for other in fibers if other != a)
    denominator = sp.prod(a - other for other in fibers if other != a)
    return sp.cancel(numerator / denominator)


base_polynomial = sp.expand(sp.prod(x - a for a in fibers))
fiber_polynomials = {
    a: sp.expand(sp.prod(s - b for b in roots))
    for a, roots in fibers.items()
}
vertical_polynomial = sp.expand(
    sum(
        lagrange_idempotent(a) * fiber_polynomials[a]
        for a in fibers
    )
)
assert sp.Poly(vertical_polynomial, x, s, domain=sp.QQ).as_expr() == (
    vertical_polynomial
)

# The ideal I=(base_polynomial, vertical_polynomial).
ideal_groebner = sp.groebner(
    [base_polynomial, vertical_polynomial],
    s,
    x,
    order="lex",
    domain=sp.QQ,
)


def reduce_mod_ideal(polynomial: sp.Expr) -> sp.Expr:
    """Reduce a polynomial modulo I."""

    return sp.expand(ideal_groebner.reduce(sp.expand(polynomial))[1])


def primitive(polynomial: sp.Expr) -> sp.Expr:
    """Take the unique s-primitive with zero constant term."""

    integrated = sp.integrate(sp.expand(polynomial), s)
    return sp.expand(integrated - integrated.subs(s, 0))


def fiber_obstructions(polynomial: sp.Expr) -> tuple[sp.Expr, ...]:
    """Return primitive differences on all nontrivial vertical fibers."""

    antiderivative = primitive(polynomial)
    equations: list[sp.Expr] = []
    for a, roots in fibers.items():
        reference = roots[0]
        reference_value = antiderivative.subs({x: a, s: reference})
        for root in roots[1:]:
            equations.append(
                sp.expand(
                    antiderivative.subs({x: a, s: root})
                    - reference_value
                )
            )
    return tuple(equations)


def interpolation_constant(polynomial: sp.Expr) -> sp.Expr:
    """Choose c(x) cancelling the primitive at one point of each fiber."""

    antiderivative = primitive(polynomial)
    return sp.expand(
        -sum(
            lagrange_idempotent(a)
            * antiderivative.subs({x: a, s: roots[0]})
            for a, roots in fibers.items()
        )
    )


def coefficient_matrix(
    expressions: Iterable[sp.Expr],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    """Return the matrix of linear expressions in the given parameters."""

    rows = [
        [sp.expand(expression).coeff(parameter) for parameter in parameters]
        for expression in expressions
        if sp.expand(expression) != 0
    ]
    if not rows:
        return sp.zeros(0, len(parameters))
    return sp.Matrix(rows)


# The complete intersection has exactly the displayed six reduced points.
point_count = sum(len(roots) for roots in fibers.values())
assert point_count == 6
for a, roots in fibers.items():
    assert base_polynomial.subs(x, a) == 0
    assert sp.expand(vertical_polynomial.subs(x, a) - fiber_polynomials[a]) == 0
    base_derivative = sp.diff(base_polynomial, x).subs(x, a)
    for root in roots:
        assert vertical_polynomial.subs({x: a, s: root}) == 0
        fiber_derivative = sp.diff(vertical_polynomial, s).subs(
            {x: a, s: root}
        )
        # This is the Jacobian determinant up to sign.  Its nonvanishing
        # makes every point reduced.
        assert base_derivative * fiber_derivative != 0

# Check the exact slice-membership criterion on a generic 3-by-5 window.
parameters = sp.symbols("c0:15")
generic_h = sp.expand(
    sum(
        parameters[5 * i + j] * x**i * s**j
        for i in range(3)
        for j in range(5)
    )
)
generic_primitive = primitive(generic_h)
generic_lift = sp.expand(
    generic_primitive + interpolation_constant(generic_h)
)
generic_ideal_groebner = sp.groebner(
    [base_polynomial, vertical_polynomial],
    s,
    x,
    order="lex",
    domain=sp.QQ.frac_field(*parameters),
)
generic_remainder = sp.expand(
    generic_ideal_groebner.reduce(generic_lift)[1]
)

remainder_polynomial = sp.Poly(generic_remainder, s, x)
remainder_coefficients = [
    coefficient for _, coefficient in remainder_polynomial.terms()
]
remainder_matrix = coefficient_matrix(remainder_coefficients, parameters)
obstruction_matrix = coefficient_matrix(
    fiber_obstructions(generic_h), parameters
)
stacked_matrix = remainder_matrix.col_join(obstruction_matrix)
assert remainder_matrix.rank() == 3
assert obstruction_matrix.rank() == 3
assert stacked_matrix.rank() == 3

# A seed may be arbitrary on singleton fibers, but it must vanish identically
# on every collision fiber.  This f does exactly that.
f = x * (x - 2)
g = 1 + x * s + s**3
assert f.subs(x, 5) != 0
for a, roots in fibers.items():
    if len(roots) > 1:
        assert sp.expand(f.subs(x, a)) == 0

for exponent in range(1, 13):
    for candidate in (sp.expand(f**exponent), sp.expand(g * f**exponent)):
        assert fiber_obstructions(candidate) == (0, 0, 0)
        lift = sp.expand(
            primitive(candidate) + interpolation_constant(candidate)
        )
        assert reduce_mod_ideal(lift) == 0
        assert sp.diff(lift, s) == candidate

# The constraints are substantive: s already fails a pure-power condition.
bad_seed = s
bad_obstructions = [
    fiber_obstructions(sp.expand(bad_seed**exponent))
    for exponent in range(1, 4)
]
assert any(any(value != 0 for value in row) for row in bad_obstructions)

# The first carrier/residual regressions use I=s^2*J with nonreduced
# length-two residual schemes.  Membership is equivalent to
#
#     primitive(h)/s^2 in J,
#
# because every element of I vanishes on s=0.
carrier = s**2


def divide_primitive_by_carrier(polynomial: sp.Expr) -> sp.Expr:
    """Return primitive(polynomial)/s^2 and require exact divisibility."""

    quotient, remainder = sp.div(primitive(polynomial), carrier, s)
    assert sp.expand(remainder) == 0
    return sp.expand(quotient)


off_carrier_ideal = sp.groebner(
    [x**2, s - 1],
    s,
    x,
    order="lex",
    domain=sp.QQ,
)
on_carrier_ideal = sp.groebner(
    [x**2, s],
    s,
    x,
    order="lex",
    domain=sp.QQ,
)


def residual_remainder(
    polynomial: sp.Expr, residual_ideal: sp.GroebnerBasis
) -> sp.Expr:
    """Reduce the carrier quotient modulo a residual ideal."""

    quotient = divide_primitive_by_carrier(polynomial)
    return sp.expand(residual_ideal.reduce(quotient)[1])


# Off the carrier, the first pure condition is a genuine interval
# cancellation.  Higher powers enter (x^2,s-1) through their x-adic order.
off_carrier_seed = x * s * (2 - 3 * s)
off_carrier_multiplier = 1 + x + s**2
for exponent in range(1, 13):
    assert residual_remainder(
        sp.expand(off_carrier_seed**exponent), off_carrier_ideal
    ) == 0
for exponent in range(2, 13):
    assert residual_remainder(
        sp.expand(
            off_carrier_multiplier * off_carrier_seed**exponent
        ),
        off_carrier_ideal,
    ) == 0

# On the carrier, powers gain s-adic order after integration and division.
on_carrier_seed = s**2
on_carrier_multiplier = 1 + x + s**3
for exponent in range(1, 13):
    for candidate in (
        sp.expand(on_carrier_seed**exponent),
        sp.expand(on_carrier_multiplier * on_carrier_seed**exponent),
    ):
        assert residual_remainder(candidate, on_carrier_ideal) == 0

# Carrier-free primary ideals split into two local cases.  The ideal
# (x^2,s-x) contains the slice s-x, so its image is the whole ring.
primary_with_slice = sp.groebner(
    [x**2, s - x],
    s,
    x,
    order="lex",
    domain=sp.QQ,
)
assert primary_with_slice.reduce(s - x)[1] == 0
assert sp.diff(s - x, s) == 1

# For the double point m^2, pure membership forces the seed to vanish at
# the support.  Integration from the support raises local order, so the
# displayed pure and mixed powers lie in D(m^2).
primary_without_slice = sp.groebner(
    [x**2, x * s, s**2],
    s,
    x,
    order="lex",
    domain=sp.QQ,
)
primary_seed = x + s
primary_multiplier = 1 + x + s**2
for exponent in range(1, 13):
    for candidate in (
        sp.expand(primary_seed**exponent),
        sp.expand(primary_multiplier * primary_seed**exponent),
    ):
        lift = primitive(candidate)
        assert primary_without_slice.reduce(lift)[1] == 0
        assert sp.diff(lift, s) == candidate

print("PASS: I is the reduced complete intersection of six displayed points")
print("PASS: generic slice membership has exactly three fiber obstructions")
print("PASS: primitive-plus-interpolation and fiber tests have the same kernel")
print("PASS: pure and mixed powers replay exactly through exponent twelve")
print("PASS: two nonreduced carrier/residual charts replay through order twelve")
print("PASS: the carrier-free primary dichotomy replays through order twelve")
print("NOTE: the all-degree Mathieu argument is in the accompanying note")
