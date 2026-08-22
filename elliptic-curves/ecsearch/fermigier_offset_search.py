"""Exact denominator-aware offset search on Fermigier's normalized quartic."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from typing import Iterable

from .fermigier import evaluate_polynomial, fermigier_quartic


Q = Fraction
SQUARE_MODULI = (64, 63, 65, 11, 13)
SQUARE_RESIDUES = {
    modulus: frozenset((value * value) % modulus for value in range(modulus))
    for modulus in SQUARE_MODULI
}


@dataclass(frozen=True)
class DenominatorOffsetPoint:
    sign: int
    offset_numerator: int
    offset_denominator: int
    x: Fraction
    normalized_y: Fraction
    raw_y: Fraction


def normalized_quartic_integer_value(
    adapter_u: Fraction | int,
    sign: int,
    offset_numerator: int,
    offset_denominator: int,
) -> tuple[int, int, int]:
    """Return ``(V,A,B^3*d^2)`` with ``R_s(x)=V/(B^6*d^4)``.

    Here ``s=A/B=2u`` is reduced and
    ``x=sign*s+offset_numerator/offset_denominator``.  The clearing factor
    ``B^6*d^4`` is a square, so the normalized quartic value is rationally a
    square exactly when the returned integer ``V`` is a square.
    """

    if sign not in (-1, 1):
        raise ValueError("sign must be -1 or 1")
    if offset_denominator < 1:
        raise ValueError("offset denominator must be positive")
    if gcd(abs(offset_numerator), offset_denominator) != 1:
        raise ValueError("offset must be reduced")
    shift = 2 * Q(adapter_u)
    a_value, b_value = shift.numerator, shift.denominator
    b2 = b_value * b_value
    a2 = a_value * a_value
    a4 = a2 * a2
    b4 = b2 * b2
    x_numerator = sign * a_value * offset_denominator + offset_numerator * b_value
    d = offset_denominator

    quartic4 = a2 + 1_149_050 * b2
    quartic3 = -30 * (62 * a2 + 68_377_393 * b2)
    quartic2 = -(2 * a4 - 1_718_550 * a2 * b2 - 1_195_214_262_641 * b4)
    quartic1 = 30 * (
        62 * a4 - 21_690_305 * a2 * b2 - 8_594_794_400_346 * b4
    )
    quartic0 = (
        a4 * a2
        - 879_500 * a4 * b2
        + 102_302_344_648 * a2 * b4
        + 18_103_855_887_324_900 * b4 * b2
    )
    value = (
        quartic4 * x_numerator**4
        + quartic3 * x_numerator**3 * b_value * d
        + quartic2 * x_numerator**2 * d**2
        + quartic1 * x_numerator * b_value * d**3
        + quartic0 * d**4
    )
    return value, a_value, b_value**3 * d**2


def _integer_square_root(value: int) -> int | None:
    if value < 0:
        return None
    for modulus in SQUARE_MODULI:
        if value % modulus not in SQUARE_RESIDUES[modulus]:
            return None
    root = isqrt(value)
    return root if root * root == value else None


def denominator_offset_points(
    adapter_u: Fraction | int,
    *,
    maximum_denominator: int,
    maximum_abs_numerator: int,
    minimum_denominator: int = 2,
) -> tuple[DenominatorOffsetPoint, ...]:
    """Exhaust the declared reduced offsets ``n/d`` on both translated charts."""

    if not 1 <= minimum_denominator <= maximum_denominator:
        raise ValueError("invalid denominator interval")
    if maximum_abs_numerator < 0:
        raise ValueError("numerator bound must be nonnegative")
    shift = 2 * Q(adapter_u)
    model = fermigier_quartic(shift)
    points: dict[Fraction, DenominatorOffsetPoint] = {}
    for denominator in range(minimum_denominator, maximum_denominator + 1):
        for numerator in range(-maximum_abs_numerator, maximum_abs_numerator + 1):
            if gcd(abs(numerator), denominator) != 1:
                continue
            for sign in (-1, 1):
                value, _, square_denominator = normalized_quartic_integer_value(
                    adapter_u, sign, numerator, denominator
                )
                root = _integer_square_root(value)
                if root is None:
                    continue
                normalized_y = Q(root, square_denominator)
                x_value = sign * shift + Q(numerator, denominator)
                raw_y = 50_616 * shift * normalized_y
                if raw_y * raw_y != evaluate_polynomial(model.quartic, x_value):
                    raise AssertionError("the cleared denominator point failed exact replay")
                points.setdefault(
                    x_value,
                    DenominatorOffsetPoint(
                        sign,
                        numerator,
                        denominator,
                        x_value,
                        normalized_y,
                        raw_y,
                    ),
                )
    return tuple(sorted(points.values(), key=lambda point: (point.x, point.sign)))


def point_stream_digest(points: Iterable[DenominatorOffsetPoint]) -> str:
    import hashlib

    digest = hashlib.sha256()
    for point in points:
        digest.update(
            (
                f"{point.sign}|{point.offset_numerator}/{point.offset_denominator}|"
                f"{point.x}|{point.raw_y}\n"
            ).encode()
        )
    return digest.hexdigest()
