"""Exact pieces of the Fermigier--Mestre family.

Fermigier fixes six roots and forms ``p6(x-s)*p6(x+s)``.  Its polynomial
square part leaves a quartic remainder.  Here ``fermigier_quartic(s)`` uses the
literal shift ``s`` from that displayed construction.

There is an unresolved normalization discrepancy in the published E22
specialization: the paper prints ``s=19754/39``, whereas exact reconstruction
of its displayed E22 model requires ``s=39508/39``.  The canonical integral
model below therefore uses an explicitly named adapter parameter ``u=s/2``;
at ``u=19754/39`` it reproduces E22.  No claim is made that this factor of two
is explained by the paper.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


Polynomial = tuple[Fraction, ...]
FERMIGIER_ROOTS = (0, 55, 314, 378, 1007, 1036)
FERMIGIER_REPORTED_PARAMETER = Fraction(19754, 39)
FERMIGIER_E22_RECONSTRUCTION_SHIFT = 2 * FERMIGIER_REPORTED_PARAMETER
FERMIGIER_SOURCE_PARAMETERS = (3, 5)
FERMIGIER_THIRTEENTH_X_INTERCEPT = Fraction(1256, 5)
FERMIGIER_THIRTEENTH_X_SLOPE = Fraction(-17, 35)
FERMIGIER_THIRTEENTH_Y_SCALE = Fraction(50616, 1225)

# Primitive discriminant in the literal symmetric shift ``s`` used by the
# paper's displayed product.  It satisfies Psi(2*u) = 16*Phi(u), where Phi is
# the canonical adapter discriminant below.
FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS = (
    178751927593521952355531210961215773383808056363395226527454720000,
    0,
    -18656650379213797389666658813347865821574479769276654052771500,
    0,
    529759946663021294253196844469955118570259147664813128103,
    0,
    -2794585724783814536180462100993070634116448921747400,
    0,
    -74195486391707382079823737680731330188240761295,
    0,
    1197167655825977922608410702757393878691200,
    0,
    -7102859996723716046015855835033267744,
    0,
    19079756640074226760165837126400,
    0,
    -17584994369128502534584064,
    0,
    -11924229808341504000,
    0,
    20185251840000,
)

# Primitive discriminant of the canonical [1,a2,1,a4,a6] model, low to high
# in the adapter parameter u=s/2.  Odd coefficients vanish.
FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS = (
    11171995474595122022220700685075985836488003522712201657965920000,
    0,
    -4664162594803449347416664703336966455393619942319163513192875,
    0,
    529759946663021294253196844469955118570259147664813128103,
    0,
    -11178342899135258144721848403972282536465795686989600,
    0,
    -1187127782267318113277179802891701283011852180720,
    0,
    76618729972862587046938284976473208236236800,
    0,
    -1818332159161271307780059093768516542464,
    0,
    19537670799436008202409817217433600,
    0,
    -72028136935950346381656326144,
    0,
    -195366581179867201536000,
    0,
    1322860664586240000,
)


@dataclass(frozen=True)
class FermigierQuartic:
    shift: Fraction
    product: Polynomial
    square_part: Polynomial
    quartic: Polynomial


def multiply_polynomials(left: Sequence[Fraction], right: Sequence[Fraction]) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return tuple(result)


def evaluate_polynomial(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def _shifted_root_polynomial(shift: Fraction) -> Polynomial:
    result: Polynomial = (Fraction(1),)
    for root in FERMIGIER_ROOTS:
        result = multiply_polynomials(
            result,
            (shift - root, Fraction(1)),
        )
    return result


def _polynomial_square_part(product_polynomial: Polynomial) -> Polynomial:
    if len(product_polynomial) != 13 or product_polynomial[12] != 1:
        raise ValueError("the Fermigier product must be monic of degree twelve")
    square_part = [Fraction(0)] * 7
    square_part[6] = Fraction(1)
    for degree in range(11, 5, -1):
        coefficient_index = degree - 6
        current_square = multiply_polynomials(square_part, square_part)
        square_part[coefficient_index] = (
            product_polynomial[degree] - current_square[degree]
        ) / 2
    return tuple(square_part)


def fermigier_quartic(shift: Fraction | int) -> FermigierQuartic:
    """Return the quartic genus-one model at the literal product shift."""

    shift = Fraction(shift)
    left = _shifted_root_polynomial(-shift)
    right = _shifted_root_polynomial(shift)
    product_polynomial = multiply_polynomials(left, right)
    square_part = _polynomial_square_part(product_polynomial)
    square = multiply_polynomials(square_part, square_part)
    remainder = tuple(
        square[degree] - product_polynomial[degree] for degree in range(13)
    )
    if any(remainder[degree] for degree in range(5, 13)):
        raise ArithmeticError("the Fermigier remainder did not reduce to a quartic")
    quartic = remainder[:5]
    if quartic[4] == 0:
        raise ArithmeticError("degenerate Fermigier specialization")
    return FermigierQuartic(
        shift=shift,
        product=product_polynomial,
        square_part=square_part,
        quartic=quartic,
    )


def twelve_visible_points(model: FermigierQuartic) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return the twelve points with abscissas ``root +/- shift``.

    Fermigier also displays a thirteenth point ``x=A+B*s``; it is returned by
    :func:`thirteenth_visible_point` rather than this backwards-compatible
    twelve-point helper.
    """

    points: list[tuple[Fraction, Fraction]] = []
    for root in FERMIGIER_ROOTS:
        for sign in (-1, 1):
            x_coordinate = Fraction(root) + sign * model.shift
            y_coordinate = evaluate_polynomial(model.square_part, x_coordinate)
            assert y_coordinate * y_coordinate == evaluate_polynomial(
                model.quartic, x_coordinate
            )
            points.append((x_coordinate, y_coordinate))
    return tuple(points)


def thirteenth_visible_point(
    model: FermigierQuartic,
    *,
    y_sign: int = 1,
) -> tuple[Fraction, Fraction]:
    """Return Fermigier's material thirteenth point on the fixed-root quartic.

    In the source two-parameter construction, the fixed roots arise at
    ``(u,v)=(3,5)`` after the affine root change ``X=x/2+465``.  Substitution
    in Fermigier's displayed ``A+B*t`` formula gives, in the fixed-root
    coordinate, ``x=1256/5-(17/35)*s``.  Exact expansion supplies the
    rational ordinate

    ``y=(50616/1225)*s*(936*s^3-254422*s^2-283436139*s+34925066050)``.

    Its negative represents the other point over the same abscissa.
    """

    if y_sign not in (-1, 1):
        raise ValueError("y_sign must be either -1 or 1")
    x_coordinate = (
        FERMIGIER_THIRTEENTH_X_INTERCEPT
        + FERMIGIER_THIRTEENTH_X_SLOPE * model.shift
    )
    shift = model.shift
    y_coordinate = y_sign * FERMIGIER_THIRTEENTH_Y_SCALE * shift * (
        936 * shift**3
        - 254422 * shift**2
        - 283436139 * shift
        + 34925066050
    )
    assert y_coordinate * y_coordinate == evaluate_polynomial(
        model.quartic, x_coordinate
    )
    return x_coordinate, y_coordinate


def thirteen_visible_points(
    model: FermigierQuartic,
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return the twelve root points and Fermigier's thirteenth point."""

    return twelve_visible_points(model) + (thirteenth_visible_point(model),)


def quartic_point_to_canonical_point(
    model: FermigierQuartic,
    point: tuple[Fraction | int, Fraction | int],
) -> tuple[Fraction, Fraction]:
    """Map a non-ramification quartic point to the canonical family model.

    This is the classical binary-quartic covariant map to
    ``Y^2=X^3-27*I*X-27*J``, followed by the exact invariant bridge to the
    repository's canonical ``[1,a2,1,a4,a6]`` model.  The covariant map is
    undefined when the supplied quartic ordinate is zero.
    """

    x_coordinate, y_coordinate = map(Fraction, point)
    if model.shift == 0:
        raise ArithmeticError("the canonical invariant bridge degenerates at shift zero")
    if y_coordinate == 0:
        raise ArithmeticError("the quartic covariant map is undefined at y=0")
    if y_coordinate * y_coordinate != evaluate_polynomial(
        model.quartic, x_coordinate
    ):
        raise ValueError("the supplied point is not on the Fermigier quartic")

    e, d, c, b, a = model.quartic
    g4 = (
        (3 * b * b - 8 * a * c) * x_coordinate**4
        + 4 * (b * c - 6 * a * d) * x_coordinate**3
        + 2 * (2 * c * c - 24 * a * e - 3 * b * d) * x_coordinate**2
        + 4 * (c * d - 6 * b * e) * x_coordinate
        + 3 * d * d
        - 8 * c * e
    )
    g6 = (
        (b**3 + 8 * a * a * d - 4 * a * b * c) * x_coordinate**6
        + 2
        * (16 * a * a * e + 2 * a * b * d - 4 * a * c * c + b * b * c)
        * x_coordinate**5
        + 5 * (8 * a * b * e + b * b * d - 4 * a * c * d) * x_coordinate**4
        + 20 * (b * b * e - a * d * d) * x_coordinate**3
        - 5 * (8 * a * d * e + b * d * d - 4 * b * c * e) * x_coordinate**2
        - 2
        * (16 * a * e * e + 2 * b * d * e - 4 * c * c * e + c * d * d)
        * x_coordinate
        - d**3
        - 8 * b * e * e
        + 4 * c * d * e
    )
    raw_x = 3 * g4 / (2 * y_coordinate) ** 2
    raw_y = 27 * g6 / (2 * y_coordinate) ** 3

    adapter_parameter = model.shift / 2
    scale = 101232 * adapter_parameter
    canonical = fermigier_canonical_coefficients(adapter_parameter)
    a1, a2, a3, _, _ = canonical
    b2 = a1 * a1 + 4 * a2
    canonical_x = (raw_x / scale**2 - 3 * b2) / 36
    canonical_y = (raw_y / (108 * scale**3) - a1 * canonical_x - a3) / 2
    return canonical_x, canonical_y


def fermigier_canonical_coefficients(
    adapter_parameter: Fraction | int,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
    """Return the canonical family model in ``u = shift/2``.

    The coefficients lie in ``Z[u]`` and their discriminant is exactly
    ``FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS`` evaluated at ``u``.  This
    removes the ``u^12`` and twelfth-power content present in the raw quartic
    conversion, both of which are nonminimal coordinate artifacts.
    """

    u = Fraction(adapter_parameter)
    a2 = -8 * u**4 + 1718550 * u**2 + 298803565660
    a4 = (
        -64 * u**8
        - 18151200 * u**6
        + 1028009011008 * u**4
        + 317946481466562025 * u**2
        + 27856983036916830925012
    )
    a6 = (
        512 * u**12
        + 35222400 * u**10
        - 31029122010320 * u**8
        + 573566223236700400 * u**6
        + 109667821527431621677482 * u**4
        + 14321756340366264921294086000 * u**2
        + 829998138277457737118423455411406
    )
    return Fraction(1), a2, Fraction(1), a4, a6


def weierstrass_discriminant(
    coefficients: Sequence[Fraction | int],
) -> Fraction:
    """Return the exact discriminant of ``[a1,a2,a3,a4,a6]``."""

    if len(coefficients) != 5:
        raise ValueError("five Weierstrass coefficients are required")
    a1, a2, a3, a4, a6 = map(Fraction, coefficients)
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    )
    return -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def weierstrass_c_invariants(
    coefficients: Sequence[Fraction | int],
) -> tuple[Fraction, Fraction]:
    """Return exact ``(c4,c6)`` for ``[a1,a2,a3,a4,a6]``."""

    if len(coefficients) != 5:
        raise ValueError("five Weierstrass coefficients are required")
    a1, a2, a3, a4, a6 = map(Fraction, coefficients)
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    return b2 * b2 - 24 * b4, -b2**3 + 36 * b2 * b4 - 216 * b6


def fermigier_discriminant_factor(adapter_parameter: Fraction | int) -> Fraction:
    """Evaluate the primitive degree-twenty canonical discriminant factor."""

    return evaluate_polynomial(
        tuple(map(Fraction, FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS)),
        Fraction(adapter_parameter),
    )
