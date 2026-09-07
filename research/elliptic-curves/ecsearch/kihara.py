"""Exact replay of Kihara's arithmetic rank-at-least-14 family.

Kihara starts with twelve symmetric roots ``b_i = a_i +/- u`` and writes

``prod_i (x-b_i) = G(x)^2 - r(x)``,

where ``r`` is a quartic.  The paper gives a one-parameter specialization and
three additional abscissas.  At ``t=2`` the first fourteen points are proved
independent in the paper by a canonical-height determinant, with the fifteenth
point used as the origin of the quartic group law.

This module reconstructs that specialization using rational arithmetic and
maps it to a Weierstrass model.  The repository's pinned finite-reduction
certificate is an independent exact replay of the fourteen-point lower bound;
it is not evidence for rank 30.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from math import comb, isqrt
from typing import Mapping, Sequence

from .fermigier import evaluate_polynomial, multiply_polynomials
from .rank_certification import (
    IndependenceCertificate,
    build_independence_certificate,
    verify_independence_certificate,
)


Polynomial = tuple[Fraction, ...]
AffinePoint = tuple[Fraction, Fraction]
WeierstrassModel = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]

KIHARA_SOURCE_DOI = "10.3792/pjaa.77.50"
KIHARA_SOURCE_URL = (
    "https://projecteuclid.org/journals/"
    "proceedings-of-the-japan-academy-series-a-mathematical-sciences/"
    "volume-77/issue-4/On-an-elliptic-curve-over-mathbf-Qt-of-rank-geq/"
    "10.3792/pjaa.77.50.pdf"
)
KIHARA_CERTIFICATE_PARAMETER = Fraction(2)


@dataclass(frozen=True)
class KiharaQuartic:
    """Kihara's quartic construction at one rational parameter."""

    parameter: Fraction
    p: Fraction
    q: Fraction
    u: Fraction
    roots: tuple[Fraction, ...]
    product: Polynomial
    square_part: Polynomial
    quartic: Polynomial


@dataclass(frozen=True)
class KiharaRank14Replay:
    """The fourteen points mapped to a Weierstrass curve."""

    quartic: KiharaQuartic
    quartic_points: tuple[AffinePoint, ...]
    quartic_origin: AffinePoint
    weierstrass_coefficients: WeierstrassModel
    weierstrass_points: tuple[AffinePoint, ...]


def _square_part(product: Polynomial) -> Polynomial:
    if len(product) != 13 or product[12] != 1:
        raise ValueError("the Kihara product must be monic of degree twelve")
    result = [Fraction(0)] * 7
    result[6] = Fraction(1)
    for degree in range(11, 5, -1):
        coefficient_index = degree - 6
        current_square = multiply_polynomials(result, result)
        result[coefficient_index] = (
            product[degree] - current_square[degree]
        ) / 2
    return tuple(result)


def _rational_square_root(value: Fraction) -> Fraction:
    if value < 0:
        raise ArithmeticError("a negative rational is not a rational square")
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if (
        numerator * numerator != value.numerator
        or denominator * denominator != value.denominator
    ):
        raise ArithmeticError("the rational value is not a square")
    return Fraction(numerator, denominator)


def kihara_quartic(parameter: Fraction | int) -> KiharaQuartic:
    """Return the quartic in Kihara's rank-14 one-parameter family."""

    t = Fraction(parameter)
    if t == 0:
        raise ArithmeticError("Kihara's displayed parameterization has a pole at t=0")
    p = t**2 * (8 + 3 * t**2)
    q = -6 * (2 + t**2) * (4 + t**2)
    u = (
        4
        * (2 + t**2)
        * (2304 + 2400 * t**2 + 928 * t**4 + 150 * t**6 + 9 * t**8)
        * (1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8)
        / t
    )
    a_values = (
        Fraction(0),
        (2 * p**2 + p * q + 2 * q**2) ** 2,
        2 * (p + q) ** 2 * (2 * p**2 + p * q + q**2),
        q**2 * (4 * p**2 - p * q + 4 * q**2),
        p * (2 * p - q) * (2 * p**2 + 4 * p * q + 5 * q**2),
        4 * p**4 + 8 * p**3 * q + 9 * p**2 * q**2 - 2 * p * q**3 + 2 * q**4,
    )
    roots = tuple(u + value for value in a_values) + tuple(
        -u + value for value in a_values
    )
    product: Polynomial = (Fraction(1),)
    for root in roots:
        product = multiply_polynomials(product, (-root, Fraction(1)))
    square_part = _square_part(product)
    square = multiply_polynomials(square_part, square_part)
    remainder = tuple(square[index] - product[index] for index in range(13))
    if any(remainder[index] for index in range(5, 13)):
        raise ArithmeticError("the Kihara remainder did not reduce to a quartic")
    quartic = remainder[:5]
    if quartic[4] == 0:
        raise ArithmeticError("degenerate Kihara specialization")
    return KiharaQuartic(
        parameter=t,
        p=p,
        q=q,
        u=u,
        roots=roots,
        product=product,
        square_part=square_part,
        quartic=quartic,
    )


def _point_from_abscissa(model: KiharaQuartic, x_coordinate: Fraction) -> AffinePoint:
    rhs = evaluate_polynomial(model.quartic, x_coordinate)
    y_coordinate = _rational_square_root(rhs)
    return x_coordinate, y_coordinate


def fifteen_visible_points(model: KiharaQuartic) -> tuple[AffinePoint, ...]:
    """Return Kihara's ``P1,...,P15`` at a rational specialization.

    The paper specifies ``y=G(b_i)`` for the first twelve points and only the
    abscissas of the final three.  For those three this replay consistently
    takes the nonnegative rational square root.
    """

    points: list[AffinePoint] = [
        (root, evaluate_polynomial(model.square_part, root))
        for root in model.roots
    ]
    for point in points:
        assert point[1] ** 2 == evaluate_polynomial(model.quartic, point[0])

    p, q, u, t = model.p, model.q, model.u, model.parameter
    common_denominator = 2 * p**2 + 2 * p * q + 3 * q**2
    p13_polynomial = (
        8 * p**6
        + 28 * p**5 * q
        + 58 * p**4 * q**2
        + 69 * p**3 * q**3
        + 76 * p**2 * q**4
        + 40 * p * q**5
        + 22 * q**6
    )
    x13 = (
        (2 * p**2 + 4 * p * q + 5 * q**2) * u + p13_polynomial
    ) / common_denominator

    shared_factor = 1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8
    p14_polynomial = (
        10616832
        - 18579456 * t
        + 33619968 * t**2
        - 51535872 * t**3
        + 45895680 * t**4
        - 61848576 * t**5
        + 35397888 * t**6
        - 41945856 * t**7
        + 16968640 * t**8
        - 17591104 * t**9
        + 5232272 * t**10
        - 4675248 * t**11
        + 1035180 * t**12
        - 769824 * t**13
        + 126252 * t**14
        - 71874 * t**15
        + 8559 * t**16
        - 2916 * t**17
        + 243 * t**18
    )
    x14 = (
        -4 * shared_factor * p14_polynomial
        / (t * (2304 + 3168 * t**2 + 1580 * t**4 + 339 * t**6 + 27 * t**8))
    )
    x15 = (
        4
        * (-48 + 24 * t - 34 * t**2 + 16 * t**3 - 6 * t**4 + 3 * t**5)
        * (96 + 80 * t**2 + 4 * t**3 + 18 * t**4 + 3 * t**5)
        * shared_factor
        / t
    )
    points.extend(_point_from_abscissa(model, value) for value in (x13, x14, x15))
    return tuple(points)


def _translate_quartic(
    coefficients: Sequence[Fraction], x_coordinate: Fraction
) -> Polynomial:
    translated = [Fraction(0)] * len(coefficients)
    for degree, coefficient in enumerate(coefficients):
        for new_degree in range(degree + 1):
            translated[new_degree] += (
                coefficient
                * comb(degree, new_degree)
                * x_coordinate ** (degree - new_degree)
            )
    return tuple(translated)


def quartic_points_to_weierstrass(
    coefficients: Sequence[Fraction],
    origin: AffinePoint,
    points: Sequence[AffinePoint],
) -> tuple[WeierstrassModel, tuple[AffinePoint, ...]]:
    """Map a quartic with a rational origin to a Weierstrass model.

    After translating the origin to ``(0,q)``, write
    ``v^2=a*u^4+b*u^3+c*u^2+d*u+q^2``.  The returned model is

    ``Y^2=X^3+c*X^2+(d*b-4*q^2*a)*X+(q^2*b^2+a*d^2-4*q^2*a*c)``.

    The map is evaluated only away from the selected origin, which maps to the
    point at infinity.
    """

    quartic = tuple(map(Fraction, coefficients))
    if len(quartic) != 5:
        raise ValueError("a quartic must have five low-to-high coefficients")
    origin_x, origin_y = map(Fraction, origin)
    if origin_y == 0:
        raise ArithmeticError("the selected quartic origin must have nonzero ordinate")
    if origin_y**2 != evaluate_polynomial(quartic, origin_x):
        raise ValueError("the selected origin is not on the quartic")
    constant, d, c, b, a = _translate_quartic(quartic, origin_x)
    assert constant == origin_y**2
    a4 = d * b - 4 * origin_y**2 * a
    a6 = origin_y**2 * b**2 + a * d**2 - 4 * origin_y**2 * a * c
    weierstrass: WeierstrassModel = (Fraction(0), c, Fraction(0), a4, a6)
    images: list[AffinePoint] = []
    for point in points:
        x_coordinate, y_coordinate = map(Fraction, point)
        if y_coordinate**2 != evaluate_polynomial(quartic, x_coordinate):
            raise ValueError("a supplied point is not on the quartic")
        local_x = x_coordinate - origin_x
        if local_x == 0:
            raise ArithmeticError(
                "the quartic-to-Weierstrass map is undefined at the origin"
            )
        X = (
            2 * origin_y * (y_coordinate + origin_y) + d * local_x
        ) / local_x**2
        Y = (
            2 * (X**2 - 4 * origin_y**2 * a) * local_x
            - 2 * d * X
            - 4 * origin_y**2 * b
        ) / (4 * origin_y)
        assert Y**2 == X**3 + c * X**2 + a4 * X + a6
        images.append((X, Y))
    return weierstrass, tuple(images)


def kihara_rank14_replay(
    parameter: Fraction | int = KIHARA_CERTIFICATE_PARAMETER,
) -> KiharaRank14Replay:
    """Reconstruct ``P1,...,P14`` with ``P15`` as the group origin."""

    quartic = kihara_quartic(parameter)
    visible = fifteen_visible_points(quartic)
    coefficients, points = quartic_points_to_weierstrass(
        quartic.quartic,
        visible[14],
        visible[:14],
    )
    return KiharaRank14Replay(
        quartic=quartic,
        quartic_points=visible[:14],
        quartic_origin=visible[14],
        weierstrass_coefficients=coefficients,
        weierstrass_points=points,
    )


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _arithmetic_digests(replay: KiharaRank14Replay) -> dict[str, str]:
    payloads = {
        "quartic_coefficients": [
            _fraction_text(value) for value in replay.quartic.quartic
        ],
        "weierstrass_coefficients": [
            _fraction_text(value) for value in replay.weierstrass_coefficients
        ],
        "weierstrass_points": [
            [_fraction_text(x_coordinate), _fraction_text(y_coordinate)]
            for x_coordinate, y_coordinate in replay.weierstrass_points
        ],
    }
    return {
        name: hashlib.sha256(
            json.dumps(value, separators=(",", ":")).encode("ascii")
        ).hexdigest()
        for name, value in payloads.items()
    }


def build_kihara_rank14_manifest(
    *, maximum_reduction_prime: int = 700
) -> dict[str, object]:
    """Build the compact exact certificate manifest for ``t=2``."""

    replay = kihara_rank14_replay()
    certificate = build_independence_certificate(
        replay.weierstrass_coefficients,
        replay.weierstrass_points,
        relation_prime=5,
        maximum_reduction_prime=maximum_reduction_prime,
    )
    return {
        "schema": "elliptic-curves.kihara-rank14-certificate.v1",
        "claim": (
            "the Kihara t=2 specialization has fourteen independent rational "
            "points; this replays the rank-at-least-14 family baseline"
        ),
        "target_status": (
            "baseline only: this is neither a rank-30 candidate nor progress "
            "toward thirty independent points"
        ),
        "source": {
            "doi": KIHARA_SOURCE_DOI,
            "url": KIHARA_SOURCE_URL,
            "paper_height_determinant": "221792776617402574.10",
        },
        "specialization": {
            "parameter_t": "2",
            "p": _fraction_text(replay.quartic.p),
            "q": _fraction_text(replay.quartic.q),
            "u": _fraction_text(replay.quartic.u),
            "point_order": "P1,...,P14 with P15 as the quartic origin",
        },
        "arithmetic_sha256": _arithmetic_digests(replay),
        "independence_certificate": certificate.to_json_object(),
        "generation": {
            "command": (
                "python3 elliptic-curves/scripts/run_kihara_rank14.py --output "
                "artifacts/generated-results/elliptic-curves/"
                "kihara_rank14_t2_v1.json"
            ),
            "maximum_reduction_prime": maximum_reduction_prime,
            "arithmetic": "Python Fraction and exact finite-field group operations",
            "external_software": "none",
        },
        "scope": {
            "specialization": (
                "unconditional rank lower bound 14 from exact rational points "
                "and a finite-reduction infinite-descent certificate"
            ),
            "generic_family": (
                "the paper supplies these as rational-function points; "
                "independence at one defined specialization rules out a "
                "generic relation"
            ),
            "not_computed": [
                "the exact rank of the t=2 specialization",
                "a conductor-optimized specialization",
                "any fifteenth point independent of P1,...,P14",
                "any rank-30 curve",
            ],
        },
    }


def verify_kihara_rank14_manifest(manifest: Mapping[str, object]) -> None:
    """Replay a pinned manifest without searching for replacement rows."""

    assert manifest["schema"] == "elliptic-curves.kihara-rank14-certificate.v1"
    specialization = manifest["specialization"]
    if not isinstance(specialization, Mapping):
        raise AssertionError("malformed Kihara specialization metadata")
    assert specialization["parameter_t"] == "2"
    replay = kihara_rank14_replay(Fraction(str(specialization["parameter_t"])))
    assert specialization["p"] == _fraction_text(replay.quartic.p)
    assert specialization["q"] == _fraction_text(replay.quartic.q)
    assert specialization["u"] == _fraction_text(replay.quartic.u)
    assert manifest["arithmetic_sha256"] == _arithmetic_digests(replay)
    raw_certificate = manifest["independence_certificate"]
    if not isinstance(raw_certificate, Mapping):
        raise AssertionError("malformed Kihara independence certificate")
    certificate = IndependenceCertificate.from_json_object(raw_certificate)
    verify_independence_certificate(
        replay.weierstrass_coefficients,
        replay.weierstrass_points,
        certificate,
    )
