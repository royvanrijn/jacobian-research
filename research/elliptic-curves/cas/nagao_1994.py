#!/usr/bin/env python3
"""Exact data for Nagao's 1994 Mestre families.

The primary source is Koh-ichi Nagao, *Construction of high-rank elliptic
curves*, Kobe J. Math. 11 (1994), 211--219.  This module separates three
different kinds of information:

* exact polynomial identities for the six-root Mestre construction;
* exact point-on-curve checks for Nagao's printed specializations; and
* the paper's cited independence results, which are not promoted to a new
  repository-local proof merely because modern PARI reproduces the numerical
  height determinants.

``SixRootMestreConstruction`` removes the fixed square content from the
quartic.  All quartic points returned here use that primitive normalization.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd, lcm
from typing import Iterable, Sequence

from mestre_root_tuples import SixRootMestreConstruction


Q = Fraction

PRIMARY_SOURCE = "https://da.lib.kobe-u.ac.jp/da/kernel/E0003610/E0003610.pdf"

RANK13_ROOTS = (148, 116, 104, 57, 25, 0)
RANK21_ROOTS = (399, 380, 352, 47, 4, 0)

RANK13_CONSTRUCTION = SixRootMestreConstruction(tuple(Q(root) for root in RANK13_ROOTS))
RANK21_CONSTRUCTION = SixRootMestreConstruction(tuple(Q(root) for root in RANK21_ROOTS))

# Nagao writes E_{14721/376}.  Direct replay of his printed minimal model in
# the q(X-T)q(X+T) constructor uses T=14721/188.  Keeping both names explicit
# prevents this factor of two from becoming a silent convention.
RANK21_PUBLISHED_PARAMETER = Q(14721, 376)
RANK21_CONSTRUCTOR_PARAMETER = 2 * RANK21_PUBLISHED_PARAMETER

RANK13_BASE_CHANGE_CONSTANT = 23550

RANK13_PUBLISHED_MODEL = (
    1,
    0,
    0,
    -1970473859866423938027563293202211,
    33666977357380599346718366106269137257495662819841,
)
RANK13_PUBLISHED_CONDUCTOR = int(
    "683806043817733548135842157722731703695131359344925852223201804394372410"
)
RANK13_CONDUCTOR_FACTORIZATION = (
    (2, 1),
    (3, 1),
    (5, 1),
    (7, 1),
    (13, 1),
    (17, 2),
    (19, 1),
    (23, 1),
    (31, 1),
    (67, 1),
    (4597, 1),
    (241261, 1),
    (274751, 1),
    (3133672472267157655375122808011344461248585991, 1),
)

RANK21_PUBLISHED_MODEL = (
    1,
    1,
    1,
    -215843772422443922015169952702159835,
    -19474361277787151947255961435459054151501792241320535,
)
RANK21_PUBLISHED_CONDUCTOR = int(
    "26112076915897777815571388664430310998157918697219343275140810790098571234096793308930"
)
RANK21_CONDUCTOR_FACTORIZATION = (
    (2, 1),
    (3, 1),
    (5, 1),
    (7, 1),
    (13, 1),
    (17, 1),
    (23, 1),
    (47, 1),
    (4507, 1),
    (
        115482611374267602141168398241396608699381902319617225736297616061235976719,
        1,
    ),
)


def _points(raw: Iterable[tuple[str, str]]) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((Q(x), Q(y)) for x, y in raw)


# Pages 216 and 218 of the paper.  They are stored as rational strings so the
# transcription remains readable and denominator cubes are not obscured.
RANK13_PUBLISHED_POINTS = _points(
    (
        ("25386174421432494", "67109523606414317896233"),
        ("26509492766907566", "-245653984040295182797783"),
        ("26399379919965810", "-214794859310156325577539"),
        ("24966084178226490", "182907525210079265587329"),
        ("24682088239127826", "260827678931979907473297"),
        ("25869170713328826", "-66811158987551134918203"),
        ("25873709433784550", "68073651220906477508921"),
        ("24768067209581670", "-237268561874887509232071"),
        ("24875645484381876", "-207753103174362148710903"),
        ("26306960164950894", "188927676301479656226153"),
        ("26597669117443566", "270396628393543357053417"),
        ("25390343800982316", "-65957033567927528804583"),
        ("230674528980533950/9", "-1830275649909410733133/27"),
    )
)

RANK21_PUBLISHED_POINTS = _points(
    (
        ("800843008889340065933/16", "22662214190910903990783584765347/64"),
        (
            "10610541066763914590637/2209",
            "1087744114825178454840094794778034/103823",
        ),
        ("907186946780634143", "728916386168451830641677698"),
        (
            "196833201085564442194083107/227919409",
            "2277807398930440819587410184793923763894/3440899317673",
        ),
        (
            "185463474139064652528000075/366301321",
            "225699857838583242849473830466481978146/7010640982619",
        ),
        (
            "-12485261071234691432503/123904",
            "1543303353428939982282171752702539/43614208",
        ),
        ("-59703014087684747037/361", "741881245094154068525036126962/6859"),
        ("-73270463404799613067/361", "866878137858638792891117943482/6859"),
        ("-360733396398627565", "106985840484096728947883974"),
        ("-389445180957906897", "74288355118790673852542098"),
        (
            "-1474458350349858512665407/14205361",
            "2278493401578368084310409028259332632/53540005609",
        ),
        (
            "-114305856035468892691779277/278589481",
            "16972779768877136292841029639987095378/4649937027371",
        ),
        ("-21972533600828202797/81", "100790786584963504563876005302/729"),
        (
            "-25047938415396324842058977/71216721",
            "68347192566984943007522052612937752062/600997908519",
        ),
        (
            "3434828081885118352213715284707/5137262501809",
            "4279912483838925044234939165329697576812433846/11643877735262694377",
        ),
        ("-227656313261676647", "133660024327268949095297798"),
        (
            "-4098089434105992137835293/12552849",
            "5660088413991351759301403659890889706/44474744007",
        ),
        (
            "2657828735869178020212617/1495729",
            "4174499731549997186596131721273201376/1829276567",
        ),
        (
            "883965004314243424124994323/850947241",
            "23250077986002214917145041708721276812178/24822981967211",
        ),
        (
            "37543938954172817209003/73441",
            "1224097915991280099903835490020298/19902511",
        ),
        (
            "19165312347502458410162233/17214201",
            "75593839815741485450348997055551694952/71421719949",
        ),
    )
)


def primitive_quartic_coefficients(
    construction: SixRootMestreConstruction, parameter: Fraction
) -> tuple[Fraction, ...]:
    """Return ascending ``(e,d,c,b,a)`` coefficients of the primitive quartic."""

    return construction.primitive_quartic_coefficients(Q(parameter))


def rank13_published_quartic_coefficients(
    parameter: Fraction,
) -> tuple[Fraction, ...]:
    """Return Nagao's printed rank-13 quartic, in ascending order."""

    parameter = Q(parameter)
    return (
        9 * parameter**6
        - 159200 * parameter**4
        + 891699592 * parameter**2
        + 4156297690000,
        2700 * parameter**4
        - 29575350 * parameter**2
        - 284435346600,
        -18 * parameter**4
        + 396150 * parameter**2
        + 6706476489,
        -2700 * parameter**2 - 63901710,
        9 * parameter**2 + 211950,
    )


def quartic_value(coefficients: Sequence[Fraction], x: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(x) + Q(coefficient)
    return answer


def primitive_visible_points(
    construction: SixRootMestreConstruction, parameter: Fraction
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return the twelve Mestre points after fixed-square content removal."""

    scale = construction.quartic_square_scale
    points = tuple(
        (x, y / scale) for x, y in construction.visible_points(Q(parameter))
    )
    coefficients = primitive_quartic_coefficients(construction, parameter)
    if any(y**2 != quartic_value(coefficients, x) for x, y in points):
        raise AssertionError("a primitive visible point failed exactly")
    return points


def binary_invariants(
    coefficients: Sequence[Fraction],
) -> tuple[Fraction, Fraction]:
    """Return classical ``I,J`` for ascending binary-quartic coefficients."""

    if len(coefficients) != 5:
        raise ValueError("a binary quartic has five coefficients")
    e, d, c, b, a = (Q(value) for value in coefficients)
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    return invariant_i, invariant_j


def short_jacobian_coefficients(
    construction: SixRootMestreConstruction, parameter: Fraction
) -> tuple[Fraction, ...]:
    invariant_i, invariant_j = binary_invariants(
        primitive_quartic_coefficients(construction, parameter)
    )
    return (Q(0), Q(0), Q(0), -27 * invariant_i, -27 * invariant_j)


def quartic_covariants_at(
    coefficients: Sequence[Fraction], x: Fraction
) -> tuple[Fraction, Fraction]:
    """Evaluate the standard binary-quartic covariants ``g,h`` at ``(x,1)``."""

    e, d, c, b, a = (Q(value) for value in coefficients)
    x = Q(x)
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    g2 = c**2 / 12 - b * d / 8 - a * e
    g3 = c * d / 12 - b * e / 2
    g4 = d**2 / 16 - c * e / 6
    u_x = 4 * a * x**3 + 3 * b * x**2 + 2 * c * x + d
    u_y = b * x**3 + 2 * c * x**2 + 3 * d * x + 4 * e
    g_value = g0 * x**4 + g1 * x**3 + g2 * x**2 + g3 * x + g4
    g_x = 4 * g0 * x**3 + 3 * g1 * x**2 + 2 * g2 * x + g3
    g_y = g1 * x**3 + 2 * g2 * x**2 + 3 * g3 * x + 4 * g4
    h_value = (u_x * g_y - u_y * g_x) / 8
    return g_value, h_value


def quartic_point_to_short_jacobian(
    construction: SixRootMestreConstruction,
    parameter: Fraction,
    point: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    """Map a checked affine quartic point to ``y^2=x^3-27Ix-27J``."""

    coefficients = primitive_quartic_coefficients(construction, parameter)
    x, y = (Q(value) for value in point)
    if y == 0:
        raise ValueError("the affine covariant map requires a nonzero ordinate")
    if y**2 != quartic_value(coefficients, x):
        raise ValueError("the supplied point is not on the primitive quartic")
    g_value, h_value = quartic_covariants_at(coefficients, x)
    jacobian_point = (36 * g_value / y**2, 108 * h_value / y**3)
    _, _, _, coefficient_a, coefficient_b = short_jacobian_coefficients(
        construction, parameter
    )
    if jacobian_point[1] ** 2 != (
        jacobian_point[0] ** 3
        + coefficient_a * jacobian_point[0]
        + coefficient_b
    ):
        raise AssertionError("the binary-quartic covariant identity failed")
    return jacobian_point


def rank13_extra_point(parameter: Fraction) -> tuple[Fraction, Fraction]:
    """Nagao's thirteenth section on the primitive rank-13 quartic."""

    parameter = Q(parameter)
    x = (parameter + 703) / 15
    y = (
        -224 * parameter**3
        - 844 * parameter**2
        + 900484 * parameter
        + 2161725
    ) / 75
    coefficients = primitive_quartic_coefficients(RANK13_CONSTRUCTION, parameter)
    if y**2 != quartic_value(coefficients, x):
        raise AssertionError("Nagao's thirteenth section failed exactly")
    return x, y


def rank13_known_quartic_points(
    parameter: Fraction,
) -> tuple[tuple[Fraction, Fraction], ...]:
    return primitive_visible_points(RANK13_CONSTRUCTION, parameter) + (
        rank13_extra_point(parameter),
    )


def rank13_base_parameter(parameter_u: Fraction) -> Fraction:
    """Return ``T=(-u^2+23550)/(2u)`` from Nagao's quadratic base change."""

    parameter_u = Q(parameter_u)
    if parameter_u == 0:
        raise ValueError("the base-change parameter u must be nonzero")
    return (RANK13_BASE_CHANGE_CONSTANT - parameter_u**2) / (2 * parameter_u)


def rank13_leading_square(parameter_u: Fraction) -> Fraction:
    """Return the square root of the base-changed quartic leading coefficient."""

    parameter_u = Q(parameter_u)
    if parameter_u == 0:
        raise ValueError("the base-change parameter u must be nonzero")
    return 3 * (parameter_u**2 + RANK13_BASE_CHANGE_CONSTANT) / (2 * parameter_u)


def rank13_base_changed_short_jacobian_coefficients(
    parameter_u: Fraction,
) -> tuple[Fraction, ...]:
    return short_jacobian_coefficients(
        RANK13_CONSTRUCTION, rank13_base_parameter(parameter_u)
    )


def rank21_short_jacobian_coefficients(
    published_parameter: Fraction = RANK21_PUBLISHED_PARAMETER,
) -> tuple[Fraction, ...]:
    """Return the short Jacobian using the record's explicit factor-two replay."""

    return short_jacobian_coefficients(
        RANK21_CONSTRUCTION, 2 * Q(published_parameter)
    )


def point_on_extended_weierstrass(
    model: Sequence[int | Fraction], point: tuple[Fraction, Fraction]
) -> bool:
    if len(model) != 5:
        raise ValueError("an extended Weierstrass model has five coefficients")
    a1, a2, a3, a4, a6 = (Q(value) for value in model)
    x, y = (Q(value) for value in point)
    return y**2 + a1 * x * y + a3 * y == x**3 + a2 * x**2 + a4 * x + a6


def polynomial_content(coefficients: Sequence[Fraction]) -> Fraction:
    """Return the positive rational content of a coefficient tuple."""

    values = tuple(Q(value) for value in coefficients if value)
    if not values:
        return Q(0)
    denominator = 1
    for value in values:
        denominator = lcm(denominator, value.denominator)
    numerators = [abs((value * denominator).numerator) for value in values]
    return Q(gcd(*numerators), denominator)


def even_discriminant_polynomial(
    construction: SixRootMestreConstruction,
) -> tuple[Fraction, ...]:
    """Return the degree-ten polynomial ``F(U)`` with ``disc(T)=F(T^2)``."""

    coefficients = construction.primitive_discriminant_polynomial
    if any(coefficients[index] for index in range(1, len(coefficients), 2)):
        raise AssertionError("the Mestre discriminant unexpectedly ceased to be even")
    return tuple(coefficients[2 * index] for index in range((len(coefficients) + 1) // 2))


def rank13_base_changed_discriminant_numerator() -> tuple[Fraction, ...]:
    """Return ``(2u)^20 disc((23550-u^2)/(2u))`` in ascending powers of u."""

    in_u_squared = even_discriminant_polynomial(RANK13_CONSTRUCTION)
    answer = [Q(0)] * 41
    constant = RANK13_BASE_CHANGE_CONSTANT
    for power, coefficient in enumerate(in_u_squared):
        exponent = 2 * power
        scalar = coefficient * 2 ** (20 - exponent)
        base_degree = 20 - exponent
        for chosen in range(exponent + 1):
            degree = base_degree + 2 * chosen
            answer[degree] += (
                scalar
                * comb(exponent, chosen)
                * constant ** (exponent - chosen)
                * (-1) ** chosen
            )
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def factorization_product(factorization: Sequence[tuple[int, int]]) -> int:
    answer = 1
    for prime, exponent in factorization:
        answer *= prime**exponent
    return answer
