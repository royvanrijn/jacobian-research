#!/usr/bin/env python3
"""Exact six-root Mestre constructions, without rank claims.

For a monic polynomial ``q(X)=prod_i (X-a_i)`` of degree six, put

``P_T(X) = q(X-T) q(X+T)``.

There is a unique monic degree-six polynomial ``g_T`` for which
``g_T^2-P_T`` has degree at most five.  Mestre's elliptic-curve construction
requires the degree-five coefficient to vanish.  When it does, the remainder
is divisible by ``T^2`` and

``Y^2 = (g_T(X)^2-P_T(X))/T^2``

is a binary quartic with twelve displayed rational points.  This module
checks those identities exactly.  It deliberately makes no assertion that
the displayed points are independent.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property
from itertools import combinations
from math import gcd, isqrt, lcm
from typing import Iterable, Sequence


Q = Fraction


def _multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return tuple(answer)


def _evaluate(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def _add(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    length = max(len(left), len(right))
    return tuple(
        (left[index] if index < len(left) else Q(0))
        + (right[index] if index < len(right) else Q(0))
        for index in range(length)
    )


def monic_polynomial_from_roots(
    roots: Iterable[Fraction],
) -> tuple[Fraction, ...]:
    """Return ascending coefficients of the monic polynomial with ``roots``."""

    answer = (Q(1),)
    for root in roots:
        root = Q(root)
        answer = _multiply(answer, (-root, Q(1)))
    return answer


def mestre_quartic_condition(
    polynomial: Sequence[Fraction],
) -> Fraction:
    """Return Mestre's degree-five obstruction for a monic sextic.

    Writing ``q=X^6+a1*X^5+...+a6``, the coefficient of ``X^5`` in
    ``(g_T^2-q(X-T)q(X+T))/T^2`` is the negative of the returned value.
    Thus the remainder is quartic exactly when this value is zero.
    """

    polynomial = tuple(Q(value) for value in polynomial)
    if len(polynomial) != 7 or polynomial[6] != 1:
        raise ValueError("the Mestre condition requires a monic sextic")
    a1, a2, a3, a4, a5 = (
        polynomial[5],
        polynomial[4],
        polynomial[3],
        polynomial[2],
        polynomial[1],
    )
    return (
        a1**5
        - 6 * a1**3 * a2
        + 7 * a1**2 * a3
        + 8 * a1 * a2**2
        - 8 * a1 * a4
        - 12 * a2 * a3
        + 24 * a5
    )


def normalize_integer_root_tuple(roots: Iterable[int]) -> tuple[int, ...]:
    """Normalize an integral six-root configuration up to affine symmetry.

    Translation sends the least root to zero, the common gcd of all
    differences is removed, and reflection in the midpoint is resolved by
    retaining the lexicographically smaller orientation.
    """

    roots = tuple(roots)
    if len(roots) != 6:
        raise ValueError("exactly six roots are required")
    if not all(isinstance(root, int) for root in roots):
        raise TypeError("integer-root normalization requires integers")
    ordered = tuple(sorted(roots))
    if len(set(ordered)) != 6:
        raise ValueError("the six roots must be distinct")
    translated = tuple(root - ordered[0] for root in ordered)
    scale = gcd(*translated[1:])
    primitive = tuple(root // scale for root in translated)
    diameter = primitive[-1]
    reflected = tuple(diameter - root for root in reversed(primitive))
    return min(primitive, reflected)


def affine_normalized_integer_root_tuples(
    max_root: int,
) -> tuple[tuple[int, ...], ...]:
    """Enumerate primitive normalized tuples with diameter at most ``max_root``."""

    if max_root < 5:
        raise ValueError("max_root must be at least 5")
    answer: list[tuple[int, ...]] = []
    for diameter in range(5, max_root + 1):
        for middle in combinations(range(1, diameter), 4):
            roots = (0, *middle, diameter)
            if gcd(*roots[1:]) != 1:
                continue
            reflected = tuple(diameter - root for root in reversed(roots))
            if roots <= reflected:
                answer.append(roots)
    return tuple(answer)


@dataclass(frozen=True)
class VisiblePointDegeneracy:
    """Specialized coincidences among the twelve displayed affine points."""

    distinct_abscissae: int
    collision_loss: int
    zero_ordinates: int


@dataclass(frozen=True)
class SixRootMestreConstruction:
    """An exact six-root polynomial-square construction."""

    roots: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        roots = tuple(sorted(Q(root) for root in self.roots))
        if len(roots) != 6 or len(set(roots)) != 6:
            raise ValueError("exactly six distinct rational roots are required")
        object.__setattr__(self, "roots", roots)

    @property
    def polynomial(self) -> tuple[Fraction, ...]:
        return monic_polynomial_from_roots(self.roots)

    @property
    def quartic_condition(self) -> Fraction:
        return mestre_quartic_condition(self.polynomial)

    @property
    def is_quartic_family(self) -> bool:
        return self.quartic_condition == 0

    @property
    def is_reflection_symmetric(self) -> bool:
        total = self.roots[0] + self.roots[-1]
        return all(
            self.roots[index] + self.roots[-1 - index] == total
            for index in range(3)
        )

    def product_coefficients(self, parameter: Fraction) -> tuple[Fraction, ...]:
        """Return ascending coefficients of ``q(X-T)q(X+T)``."""

        parameter = Q(parameter)
        shifted_roots = tuple(root + parameter for root in self.roots) + tuple(
            root - parameter for root in self.roots
        )
        return monic_polynomial_from_roots(shifted_roots)

    def square_approximant_coefficients(
        self, parameter: Fraction
    ) -> tuple[Fraction, ...]:
        """Return the unique monic ``g_T`` matching degrees 12 through 6."""

        product = self.product_coefficients(parameter)
        approximant = [Q(0)] * 7
        approximant[6] = Q(1)
        for index in range(5, -1, -1):
            square = _multiply(approximant, approximant)
            degree = 6 + index
            approximant[index] = (product[degree] - square[degree]) / 2
        return tuple(approximant)

    def remainder_coefficients(
        self, parameter: Fraction
    ) -> tuple[Fraction, ...]:
        """Return ascending coefficients of ``g_T^2-q(X-T)q(X+T)``."""

        product = self.product_coefficients(parameter)
        approximant = self.square_approximant_coefficients(parameter)
        square = _multiply(approximant, approximant)
        return tuple(left - right for left, right in zip(square, product))

    def quartic_coefficients(
        self, parameter: Fraction
    ) -> tuple[Fraction, ...]:
        """Return ascending coefficients of the normalized quartic ``R_T``."""

        parameter = Q(parameter)
        if parameter == 0:
            raise ValueError("the direct R_T/T^2 formula requires T != 0")
        if not self.is_quartic_family:
            raise ValueError("this root tuple leaves a degree-five remainder")
        remainder = self.remainder_coefficients(parameter)
        if any(remainder[index] for index in range(5, len(remainder))):
            raise AssertionError("the exact Mestre quartic condition failed")
        answer = tuple(value / parameter**2 for value in remainder[:5])
        if all(value == 0 for value in answer):
            raise AssertionError("the normalized quartic vanished identically")
        return answer

    @cached_property
    def quartic_content(self) -> Fraction:
        """Return the fixed rational-square content of the quartic family.

        Every quartic coefficient has degree at most six in ``T``.  The gcd of
        its values on seven consecutive integers therefore gives its fixed
        integer-valued content.  Combining all five coefficients makes the
        normalization canonical for this API.  In the configurations used by
        this construction the content is a rational square; this is checked.
        """

        values = [
            coefficient
            for parameter in range(1, 8)
            for coefficient in self.quartic_coefficients(Q(parameter))
        ]
        common_denominator = 1
        for value in values:
            common_denominator = lcm(common_denominator, value.denominator)
        integer_values = [
            abs((value * common_denominator).numerator) for value in values
        ]
        content = Q(gcd(*integer_values), common_denominator)
        if (
            isqrt(content.numerator) ** 2 != content.numerator
            or isqrt(content.denominator) ** 2 != content.denominator
        ):
            raise AssertionError("the fixed quartic content is not a rational square")
        return content

    @cached_property
    def quartic_square_scale(self) -> Fraction:
        """Return the positive square root of :attr:`quartic_content`."""

        return Q(
            isqrt(self.quartic_content.numerator),
            isqrt(self.quartic_content.denominator),
        )

    def primitive_quartic_coefficients(
        self, parameter: Fraction
    ) -> tuple[Fraction, ...]:
        """Remove the tuple's fixed rational-square quartic content."""

        return tuple(
            value / self.quartic_content
            for value in self.quartic_coefficients(parameter)
        )

    def primitive_binary_invariants(
        self, parameter: Fraction
    ) -> tuple[Fraction, Fraction]:
        """Return ``I,J`` after fixed-square quartic normalization."""

        invariant_i, invariant_j = self.binary_invariants(parameter)
        return (
            invariant_i / self.quartic_content**2,
            invariant_j / self.quartic_content**3,
        )

    def primitive_jacobian_coefficients(
        self, parameter: Fraction
    ) -> tuple[Fraction, ...]:
        """Return ``[0,0,0,-27I,-27J]`` for the primitive quartic."""

        invariant_i, invariant_j = self.primitive_binary_invariants(parameter)
        return (Q(0), Q(0), Q(0), -27 * invariant_i, -27 * invariant_j)

    def quartic_value(self, parameter: Fraction, x: Fraction) -> Fraction:
        return _evaluate(self.quartic_coefficients(parameter), Q(x))

    def binary_invariants(
        self, parameter: Fraction
    ) -> tuple[Fraction, Fraction]:
        """Return classical ``I,J`` for the normalized binary quartic."""

        e, d, c, b, a = self.quartic_coefficients(parameter)
        invariant_i = 12 * a * e - 3 * b * d + c**2
        invariant_j = (
            72 * a * c * e
            + 9 * b * c * d
            - 27 * a * d**2
            - 27 * b**2 * e
            - 2 * c**3
        )
        return invariant_i, invariant_j

    def quartic_discriminant(self, parameter: Fraction) -> Fraction:
        """Return the exact binary-quartic discriminant ``(4I^3-J^2)/27``."""

        invariant_i, invariant_j = self.binary_invariants(parameter)
        return (4 * invariant_i**3 - invariant_j**2) / 27

    def primitive_quartic_discriminant(self, parameter: Fraction) -> Fraction:
        """Return the discriminant after fixed-square content removal."""

        return self.quartic_discriminant(parameter) / self.quartic_content**6

    @cached_property
    def primitive_discriminant_polynomial(self) -> tuple[Fraction, ...]:
        """Return the primitive quartic discriminant as a polynomial in ``T``.

        Its degree is at most 20.  We recover it from 21 consecutive exact
        values in the Newton binomial basis, verify the next finite difference
        vanishes, and convert to ascending power-basis coefficients.
        """

        values = [
            self.primitive_quartic_discriminant(Q(parameter))
            for parameter in range(1, 23)
        ]
        differences = values
        newton_coefficients: list[Fraction] = []
        while differences:
            newton_coefficients.append(differences[0])
            differences = [
                right - left
                for left, right in zip(differences, differences[1:])
            ]
        if newton_coefficients[21] != 0:
            raise AssertionError("the discriminant exceeded its degree-20 bound")

        answer = (Q(0),)
        # basis_j = binomial(T-1,j)
        basis = (Q(1),)
        for index, coefficient in enumerate(newton_coefficients[:21]):
            answer = _add(
                answer,
                tuple(coefficient * value for value in basis),
            )
            next_factor = (Q(-(index + 1)), Q(1))
            basis = tuple(
                value / (index + 1) for value in _multiply(basis, next_factor)
            )
        while len(answer) > 1 and answer[-1] == 0:
            answer = answer[:-1]
        if _evaluate(answer, Q(23)) != self.primitive_quartic_discriminant(Q(23)):
            raise AssertionError("discriminant interpolation failed its check value")
        return answer

    def primitive_discriminant_value(self, parameter: Fraction) -> Fraction:
        """Evaluate the cached degree-at-most-20 discriminant polynomial."""

        return _evaluate(self.primitive_discriminant_polynomial, Q(parameter))

    def visible_points(
        self, parameter: Fraction
    ) -> tuple[tuple[Fraction, Fraction], ...]:
        """Return and exactly check the twelve displayed quartic points."""

        parameter = Q(parameter)
        if parameter == 0:
            raise ValueError("the normalized visible-point formula requires T != 0")
        approximant = self.square_approximant_coefficients(parameter)
        points: list[tuple[Fraction, Fraction]] = []
        for root in self.roots:
            for sign in (-1, 1):
                x = root + sign * parameter
                y = _evaluate(approximant, x) / parameter
                if y**2 != self.quartic_value(parameter, x):
                    raise AssertionError("a displayed Mestre point failed exactly")
                points.append((x, y))
        return tuple(points)

    def visible_point_degeneracy(
        self, parameter: Fraction
    ) -> VisiblePointDegeneracy:
        points = self.visible_points(parameter)
        distinct_abscissae = len({point[0] for point in points})
        return VisiblePointDegeneracy(
            distinct_abscissae=distinct_abscissae,
            collision_loss=12 - distinct_abscissae,
            zero_ordinates=sum(point[1] == 0 for point in points),
        )

    def collision_parameters(self) -> tuple[Fraction, ...]:
        """Return positive ``T`` where two displayed abscissae coincide."""

        return tuple(
            sorted(
                {
                    (right - left) / 2
                    for left, right in combinations(self.roots, 2)
                }
            )
        )
