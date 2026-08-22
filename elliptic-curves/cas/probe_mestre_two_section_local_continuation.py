#!/usr/bin/env python3
"""Local two-section continuation at the rank-13 Mestre D-square seed.

This is deliberately a *jet* computation.  It reconstructs ``R`` by the
monic-square recursion and carries only values and first derivatives of its
coefficients.  In particular, it never builds any of the degree 40, 50, or 60
expanded residual polynomials from ``mestre_affine_section_elimination.sing``.

At the affine-normalized root point for ``(0,25,95,143,168,205)`` we impose
the two nonvisible sections with normalized abscissae

    -95/23 + (37/23) T,   583/115 + (13/23) T.

The seven equations are ``M`` and the three recursive residuals for each
section.  The labelled incidence has rank-six Jacobian and a two-dimensional
tangent excess at the displayed seeds.  The script gives finite-field and
p-adic Hensel witnesses for that unresolved local germ, not a smooth
one-dimensional component.  Those witnesses do not provide a rational
parameterization, a global component equation, section intersections, a
Shioda Gram matrix, saturation, or a rank-14 claim.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import dataclass
from fractions import Fraction
from math import comb, gcd, isqrt
from pathlib import Path
from typing import Iterable, Sequence


Q = Fraction
VARIABLES = ("c1", "c2", "c3", "c4", "x01", "x11", "x02", "x12")
EQUATION_NAMES = ("M", "E2_1", "E1_1", "E0_1", "E2_2", "E1_2", "E0_2")
# Coefficients of
# (X-19/5)(X-143/25)(X-168/25)(X-41/5), in descending order after 1.
# They are reconstructed from the roots used by the independent Singular
# verifier, rather than copied from a secondary display in the prose note.
SEED_MODULI = (
    Q(-611, 25), Q(136799, 625), Q(-530557, 625), Q(18714696, 15625)
)
SECTIONS = (
    (Q(-95, 23), Q(37, 23)),
    (Q(583, 115), Q(13, 23)),
)
COMPARISON_MODULI = (
    Q(-23), Q(4489, 23), Q(-8810207, 12167), Q(277065600, 279841)
)
COMPARISON_SECTIONS = (
    (Q(619, 69), Q(31, 21)),
    (Q(619, 69), Q(-31, 21)),
    (Q(304, 69), Q(19, 21)),
    (Q(304, 69), Q(-19, 21)),
    (Q(115, 69), Q(1, 21)),
    (Q(115, 69), Q(-1, 21)),
)
DEFAULT_PRIMES = (17, 29, 37)


@dataclass(frozen=True)
class Field:
    """Either Q or one prime field, with explicit rational reduction."""

    modulus: int | None = None
    tangent: bool = True
    series_order: int = 0

    def coerce(self, value: int | Fraction) -> int | Fraction:
        value = Q(value)
        if self.modulus is None:
            return value
        denominator = value.denominator % self.modulus
        if denominator == 0:
            raise ZeroDivisionError("a seed denominator is not invertible")
        return (value.numerator % self.modulus) * pow(
            denominator, -1, self.modulus
        ) % self.modulus

    def add(self, left: int | Fraction, right: int | Fraction) -> int | Fraction:
        answer = left + right
        return answer if self.modulus is None else answer % self.modulus

    def neg(self, value: int | Fraction) -> int | Fraction:
        return -value if self.modulus is None else (-value) % self.modulus

    def mul(self, left: int | Fraction, right: int | Fraction) -> int | Fraction:
        answer = left * right
        return answer if self.modulus is None else answer % self.modulus

    def inverse(self, value: int | Fraction) -> int | Fraction:
        if self.modulus is None:
            return Q(1, 1) / value
        return pow(int(value) % self.modulus, -1, self.modulus)


@dataclass(frozen=True)
class Jet:
    """A scalar together with its derivatives in the eight input variables."""

    field: Field
    value: int | Fraction
    gradient: tuple[int | Fraction, ...]

    @classmethod
    def constant(cls, field: Field, value: int | Fraction) -> "Jet":
        dimension = len(VARIABLES) if field.tangent else 0
        return cls(field, field.coerce(value), (field.coerce(0),) * dimension)

    @classmethod
    def variable(
        cls, field: Field, value: int | Fraction, index: int
    ) -> "Jet":
        gradient = [field.coerce(0)] * (len(VARIABLES) if field.tangent else 0)
        if field.tangent:
            gradient[index] = field.coerce(1)
        return cls(field, field.coerce(value), tuple(gradient))

    def __add__(self, other: "Jet") -> "Jet":
        return Jet(
            self.field,
            self.field.add(self.value, other.value),
            tuple(
                self.field.add(left, right)
                for left, right in zip(self.gradient, other.gradient)
            ),
        )

    def __neg__(self) -> "Jet":
        return Jet(
            self.field,
            self.field.neg(self.value),
            tuple(self.field.neg(value) for value in self.gradient),
        )

    def __sub__(self, other: "Jet") -> "Jet":
        return self + (-other)

    def __mul__(self, other: "Jet") -> "Jet":
        field = self.field
        return Jet(
            field,
            field.mul(self.value, other.value),
            tuple(
                field.add(
                    field.mul(left, other.value), field.mul(self.value, right)
                )
                for left, right in zip(self.gradient, other.gradient)
            ),
        )

    def power(self, exponent: int) -> "Jet":
        answer = type(self).constant(self.field, 1)
        base = self
        while exponent:
            if exponent & 1:
                answer = answer * base
            base = base * base
            exponent //= 2
        return answer

    def is_zero(self) -> bool:
        return self.value == 0 and all(value == 0 for value in self.gradient)


@dataclass(frozen=True)
class OrderTwoJet:
    """A degree-two series in one formal continuation parameter."""

    field: Field
    value: Fraction
    linear: Fraction
    quadratic: Fraction

    @classmethod
    def constant(cls, field: Field, value: int | Fraction) -> "OrderTwoJet":
        return cls(field, Q(value), Q(0), Q(0))

    @classmethod
    def seed(cls, field: Field, value: Fraction, linear: Fraction) -> "OrderTwoJet":
        return cls(field, Q(value), Q(linear), Q(0))

    def __add__(self, other: "OrderTwoJet") -> "OrderTwoJet":
        return OrderTwoJet(
            self.field,
            self.value + other.value,
            self.linear + other.linear,
            self.quadratic + other.quadratic,
        )

    def __neg__(self) -> "OrderTwoJet":
        return OrderTwoJet(self.field, -self.value, -self.linear, -self.quadratic)

    def __sub__(self, other: "OrderTwoJet") -> "OrderTwoJet":
        return self + (-other)

    def __mul__(self, other: "OrderTwoJet") -> "OrderTwoJet":
        return OrderTwoJet(
            self.field,
            self.value * other.value,
            self.value * other.linear + self.linear * other.value,
            self.value * other.quadratic + self.linear * other.linear + self.quadratic * other.value,
        )

    def power(self, exponent: int) -> "OrderTwoJet":
        answer = OrderTwoJet.constant(self.field, 1)
        base = self
        while exponent:
            if exponent & 1:
                answer = answer * base
            base = base * base
            exponent //= 2
        return answer

    def is_zero(self) -> bool:
        return self.value == self.linear == self.quadratic == 0


@dataclass(frozen=True)
class OrderThreeJet:
    """A degree-three series in one formal continuation parameter."""

    field: Field
    value: Fraction
    linear: Fraction
    quadratic: Fraction
    cubic: Fraction

    @classmethod
    def constant(cls, field: Field, value: int | Fraction) -> "OrderThreeJet":
        return cls(field, Q(value), Q(0), Q(0), Q(0))

    @classmethod
    def seed(
        cls, field: Field, value: Fraction, linear: Fraction, quadratic: Fraction
    ) -> "OrderThreeJet":
        return cls(field, Q(value), Q(linear), Q(quadratic), Q(0))

    def __add__(self, other: "OrderThreeJet") -> "OrderThreeJet":
        return OrderThreeJet(
            self.field,
            self.value + other.value,
            self.linear + other.linear,
            self.quadratic + other.quadratic,
            self.cubic + other.cubic,
        )

    def __neg__(self) -> "OrderThreeJet":
        return OrderThreeJet(
            self.field, -self.value, -self.linear, -self.quadratic, -self.cubic
        )

    def __sub__(self, other: "OrderThreeJet") -> "OrderThreeJet":
        return self + (-other)

    def __mul__(self, other: "OrderThreeJet") -> "OrderThreeJet":
        return OrderThreeJet(
            self.field,
            self.value * other.value,
            self.value * other.linear + self.linear * other.value,
            self.value * other.quadratic + self.linear * other.linear + self.quadratic * other.value,
            self.value * other.cubic + self.linear * other.quadratic
            + self.quadratic * other.linear + self.cubic * other.value,
        )

    def power(self, exponent: int) -> "OrderThreeJet":
        answer = OrderThreeJet.constant(self.field, 1)
        base = self
        while exponent:
            if exponent & 1:
                answer = answer * base
            base = base * base
            exponent //= 2
        return answer

    def is_zero(self) -> bool:
        return self.value == self.linear == self.quadratic == self.cubic == 0


@dataclass(frozen=True)
class FormalSeries:
    """A truncated Q[[t]] coefficient vector for recursive implicit lifting."""

    field: Field
    coefficients: tuple[Fraction, ...]

    @property
    def value(self) -> Fraction:
        return self.coefficients[0]

    @classmethod
    def constant(cls, field: Field, value: int | Fraction) -> "FormalSeries":
        return cls(field, (Q(value),) + (Q(0),) * field.series_order)

    @classmethod
    def seed(cls, field: Field, coefficients: Sequence[Fraction]) -> "FormalSeries":
        if len(coefficients) != field.series_order + 1:
            raise ValueError("wrong formal-series truncation")
        return cls(field, tuple(map(Q, coefficients)))

    def _coerce(self, other: "FormalSeries | int | Fraction") -> "FormalSeries":
        if isinstance(other, FormalSeries):
            if other.field != self.field:
                raise ValueError("formal series belong to different fields")
            return other
        return FormalSeries.constant(self.field, other)

    def __add__(self, other: "FormalSeries | int | Fraction") -> "FormalSeries":
        other = self._coerce(other)
        return FormalSeries(
            self.field,
            tuple(left + right for left, right in zip(self.coefficients, other.coefficients)),
        )

    def __neg__(self) -> "FormalSeries":
        return FormalSeries(self.field, tuple(-value for value in self.coefficients))

    def __sub__(self, other: "FormalSeries | int | Fraction") -> "FormalSeries":
        return self + (-other)

    def __radd__(self, other: "int | Fraction") -> "FormalSeries":
        return self + other

    def __rsub__(self, other: "int | Fraction") -> "FormalSeries":
        return self._coerce(other) - self

    def __mul__(self, other: "FormalSeries | int | Fraction") -> "FormalSeries":
        other = self._coerce(other)
        order = self.field.series_order
        return FormalSeries(
            self.field,
            tuple(
                sum(
                    (self.coefficients[left] * other.coefficients[degree - left]
                     for left in range(degree + 1)),
                    Q(0),
                )
                for degree in range(order + 1)
            ),
        )

    def __rmul__(self, other: "int | Fraction") -> "FormalSeries":
        return self * other

    def inverse(self) -> "FormalSeries":
        """Return the inverse of a series with nonzero constant coefficient."""

        if not self.value:
            raise ZeroDivisionError("a formal series with zero constant term is not a unit")
        coefficients = [Q(1, 1) / self.value]
        for degree in range(1, self.field.series_order + 1):
            coefficients.append(
                -sum(
                    (
                        self.coefficients[index] * coefficients[degree - index]
                        for index in range(1, degree + 1)
                    ),
                    Q(0),
                )
                / self.value
            )
        return FormalSeries.seed(self.field, coefficients)

    def __truediv__(self, other: "FormalSeries | int | Fraction") -> "FormalSeries":
        return self * self._coerce(other).inverse()

    def __rtruediv__(self, other: "int | Fraction") -> "FormalSeries":
        return self._coerce(other) / self

    def power(self, exponent: int) -> "FormalSeries":
        answer = FormalSeries.constant(self.field, 1)
        base = self
        while exponent:
            if exponent & 1:
                answer = answer * base
            base = base * base
            exponent //= 2
        return answer

    def __pow__(self, exponent: int) -> "FormalSeries":
        return self.power(exponent)

    def is_zero(self) -> bool:
        return all(value == 0 for value in self.coefficients)


@dataclass(frozen=True)
class FormalBivariate:
    """A truncated Q[[u,v]] series, cut by total degree."""

    field: Field
    coefficients: tuple[tuple[tuple[int, int], Fraction], ...]

    @property
    def value(self) -> Fraction:
        return dict(self.coefficients).get((0, 0), Q(0))

    @classmethod
    def constant(cls, field: Field, value: int | Fraction) -> "FormalBivariate":
        return cls.seed(field, {(0, 0): Q(value)})

    @classmethod
    def seed(
        cls, field: Field, coefficients: dict[tuple[int, int], Fraction]
    ) -> "FormalBivariate":
        order = field.series_order
        return cls(
            field,
            tuple(
                sorted(
                    (monomial, Q(value))
                    for monomial, value in coefficients.items()
                    if sum(monomial) <= order and value
                )
            ),
        )

    def dictionary(self) -> dict[tuple[int, int], Fraction]:
        return dict(self.coefficients)

    def coefficient(self, monomial: tuple[int, int]) -> Fraction:
        return self.dictionary().get(monomial, Q(0))

    def __add__(self, other: "FormalBivariate") -> "FormalBivariate":
        answer = self.dictionary()
        for monomial, value in other.coefficients:
            answer[monomial] = answer.get(monomial, Q(0)) + value
        return FormalBivariate.seed(self.field, answer)

    def __neg__(self) -> "FormalBivariate":
        return FormalBivariate.seed(
            self.field, {monomial: -value for monomial, value in self.coefficients}
        )

    def __sub__(self, other: "FormalBivariate") -> "FormalBivariate":
        return self + (-other)

    def __mul__(self, other: "FormalBivariate") -> "FormalBivariate":
        order = self.field.series_order
        answer: dict[tuple[int, int], Fraction] = {}
        for (left_u, left_v), left_value in self.coefficients:
            for (right_u, right_v), right_value in other.coefficients:
                monomial = (left_u + right_u, left_v + right_v)
                if sum(monomial) <= order:
                    answer[monomial] = answer.get(monomial, Q(0)) + left_value * right_value
        return FormalBivariate.seed(self.field, answer)

    def power(self, exponent: int) -> "FormalBivariate":
        answer = FormalBivariate.constant(self.field, 1)
        base = self
        while exponent:
            if exponent & 1:
                answer = answer * base
            base = base * base
            exponent //= 2
        return answer

    def is_zero(self) -> bool:
        return not self.coefficients


Poly = dict[tuple[int, int], Jet]


def jet_is_zero(jet: Jet) -> bool:
    """Do not discard a vanishing value whose first derivative is nonzero."""

    return jet.is_zero()


def poly_add(left: Poly, right: Poly) -> Poly:
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer[key] + value if key in answer else value
    return {key: value for key, value in answer.items() if not jet_is_zero(value)}


def poly_neg(poly: Poly) -> Poly:
    return {key: -value for key, value in poly.items()}


def poly_mul(left: Poly, right: Poly) -> Poly:
    answer: Poly = {}
    for (x_left, t_left), left_value in left.items():
        for (x_right, t_right), right_value in right.items():
            key = (x_left + x_right, t_left + t_right)
            term = left_value * right_value
            answer[key] = answer[key] + term if key in answer else term
    return {key: value for key, value in answer.items() if not jet_is_zero(value)}


def shifted_univariate(coefficients: Sequence[Jet], sign: int) -> Poly:
    """Return q(X + sign*T), with q stored in ascending X degree."""

    answer: Poly = {}
    for degree, coefficient in enumerate(coefficients):
        for x_degree in range(degree + 1):
            t_degree = degree - x_degree
            scalar = comb(degree, x_degree) * sign**t_degree
            answer[(x_degree, t_degree)] = coefficient * type(coefficient).constant(
                coefficient.field, scalar
            )
    return answer


def substitute_line(poly: Poly, intercept: Jet, slope: Jet) -> list[Jet]:
    """Return coefficients of poly(intercept+slope*T,T) in ascending T."""

    coefficients = [type(intercept).constant(intercept.field, 0) for _ in range(7)]
    for (x_degree, t_degree), value in poly.items():
        for slope_degree in range(x_degree + 1):
            degree = t_degree + slope_degree
            if degree >= len(coefficients):
                coefficients.extend(
                    type(intercept).constant(intercept.field, 0)
                    for _ in range(degree + 1 - len(coefficients))
                )
            term = (
                value
                * intercept.power(x_degree - slope_degree)
                * slope.power(slope_degree)
                * type(intercept).constant(intercept.field, comb(x_degree, slope_degree))
            )
            coefficients[degree] = coefficients[degree] + term
    return coefficients


def residuals(coords: Sequence[int | Fraction], field: Field) -> list[Jet]:
    """Evaluate M and two copies of (E2,E1,E0) in recursive form."""

    if len(coords) != len(VARIABLES):
        raise ValueError("wrong coordinate count")
    jets = [Jet.variable(field, value, index) for index, value in enumerate(coords)]
    return residuals_from_jets(jets)


def residuals_from_jets(jets: Sequence[Jet | OrderTwoJet]) -> list[Jet | OrderTwoJet]:
    """Evaluate the recursive residuals for any compatible truncated jet."""

    if len(jets) != len(VARIABLES):
        raise ValueError("wrong coordinate count")
    c1, c2, c3, c4, x01, x11, x02, x12 = jets
    field = c1.field
    kind = type(c1)
    zero = kind.constant(field, 0)
    one = kind.constant(field, 1)
    # q=X(X-1)(X^4+c1 X^3+c2 X^2+c3 X+c4), ascending in X.
    q = [zero] * 7
    for index, coefficient in enumerate((c4, c3, c2, c1, one)):
        q[index + 2] = q[index + 2] + coefficient
        q[index + 1] = q[index + 1] - coefficient
    product = poly_mul(shifted_univariate(q, -1), shifted_univariate(q, 1))
    g: Poly = {(6, 0): one}
    half = kind.constant(field, field.inverse(field.coerce(2)))
    for lower_degree in range(5, -1, -1):
        target_degree = 6 + lower_degree
        square = poly_mul(g, g)
        product_part = {(0, t): value for (x, t), value in product.items() if x == target_degree}
        square_part = {(0, t): value for (x, t), value in square.items() if x == target_degree}
        correction = poly_mul(poly_add(product_part, poly_neg(square_part)), {(0, 0): half})
        g = poly_add(g, {(lower_degree, t): value for (_, t), value in correction.items()})
    numerator = poly_add(poly_mul(g, g), poly_neg(product))
    remainder = {(x, t - 2): value for (x, t), value in numerator.items() if t >= 2}

    a1, a2, a3, a4, a5 = q[5], q[4], q[3], q[2], q[1]
    mestre = (
        a1.power(5) - a1.power(3) * a2 * kind.constant(field, 6)
        + a1.power(2) * a3 * kind.constant(field, 7)
        + a1 * a2.power(2) * kind.constant(field, 8)
        - a1 * a4 * kind.constant(field, 8)
        - a2 * a3 * kind.constant(field, 12)
        + a5 * kind.constant(field, 24)
    )
    leading = (
        a1.power(4) * kind.constant(field, 5)
        - a1.power(2) * a2 * kind.constant(field, 24)
        + a1 * a3 * kind.constant(field, 32)
        + a2.power(2) * kind.constant(field, 16)
        - a4 * kind.constant(field, 64)
    )

    answer = [mestre]
    for intercept, slope in ((x01, x11), (x02, x12)):
        f = substitute_line(remainder, intercept, slope)
        s2 = (one - slope.power(2)).power(2) * leading * kind.constant(field, Q(1, 4))
        n1 = s2 * f[4] * kind.constant(field, 4) - f[5].power(2)
        n0 = s2.power(2) * f[3] * kind.constant(field, 8) - n1 * f[5]
        e2 = s2.power(3) * f[2] * kind.constant(field, 64) - n1.power(2) - n0 * f[5] * kind.constant(field, 4)
        e1 = s2.power(4) * f[1] * kind.constant(field, 64) - n0 * n1
        e0 = s2.power(5) * f[0] * kind.constant(field, 256) - n0.power(2)
        answer.extend((e2, e1, e0))
    return answer


def row_reduce(matrix: Sequence[Sequence[int | Fraction]], field: Field) -> tuple[int, list[int]]:
    work = [list(row) for row in matrix]
    row = 0
    pivots: list[int] = []
    if not work:
        return 0, pivots
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work)) if work[index][column] != 0), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = field.inverse(work[row][column])
        work[row] = [field.mul(value, inverse) for value in work[row]]
        for index in range(len(work)):
            if index == row or work[index][column] == 0:
                continue
            scalar = work[index][column]
            work[index] = [
                field.add(value, field.neg(field.mul(scalar, pivot_value)))
                for value, pivot_value in zip(work[index], work[row])
            ]
        pivots.append(column)
        row += 1
        if row == len(work):
            break
    return row, pivots


def solve_square(matrix: Sequence[Sequence[int]], right: Sequence[int], prime: int) -> list[int]:
    augmented = [list(row) + [value] for row, value in zip(matrix, right)]
    size = len(matrix)
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column] % prime), None)
        if pivot is None:
            raise AssertionError("chosen transverse Jacobian minor is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column] % prime, -1, prime)
        augmented[column] = [(value * inverse) % prime for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scalar = augmented[row][column] % prime
            augmented[row] = [
                (value - scalar * pivot_value) % prime
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[-1] % prime for row in augmented]


def solve_square_over_q(
    matrix: Sequence[Sequence[Fraction]], right: Sequence[Fraction]
) -> list[Fraction]:
    """Exact Gauss solve used only for the two tangent-kernel directions."""

    size = len(matrix)
    augmented = [[Q(value) for value in row] + [Q(target)] for row, target in zip(matrix, right)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column]), None)
        if pivot is None:
            raise AssertionError("the exact six-by-six tangent minor is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = Q(1, 1) / augmented[column][column]
        augmented[column] = [value * inverse for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scalar = augmented[row][column]
            augmented[row] = [
                value - scalar * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[-1] for row in augmented]


def tangent_quadratic_data() -> dict[str, object]:
    """Compute the exact quadratic obstruction on the two-dimensional tangent plane."""

    seed = SEED_MODULI + SECTIONS[0] + SECTIONS[1]
    jacobian = [list(jet.gradient) for jet in residuals(seed, Field())]
    rows = list(range(6))
    pivot_columns = list(range(6))
    free_columns = [6, 7]
    pivot_matrix = [[jacobian[row][column] for column in pivot_columns] for row in rows]
    directions = []
    for free in free_columns:
        right = [-jacobian[row][free] for row in rows]
        direction = solve_square_over_q(pivot_matrix, right)
        full = direction + [Q(0), Q(0)]
        full[free] = Q(1)
        directions.append(full)
    for direction in directions:
        if any(sum(entry * derivative for entry, derivative in zip(direction, row)) != 0 for row in jacobian):
            raise AssertionError("computed vector is not tangent to all seven equations")

    def quadratic(direction: Sequence[Fraction]) -> list[Fraction]:
        jets = [
            OrderTwoJet.seed(Field(), value, derivative)
            for value, derivative in zip(seed, direction)
        ]
        return [jet.quadratic for jet in residuals_from_jets(jets)]

    first = quadratic(directions[0])
    second = quadratic(directions[1])
    combined = quadratic([left + right for left, right in zip(directions[0], directions[1])])
    cross = [value - left - right for value, left, right in zip(combined, first, second)]
    seventh = jacobian[6]

    def corrected_obstruction(quadratic_terms: Sequence[Fraction]) -> Fraction:
        correction = solve_square_over_q(
            pivot_matrix, [-quadratic_terms[row] for row in rows]
        )
        return quadratic_terms[6] + sum(
            coefficient * value
            for coefficient, value in zip(seventh[:6], correction)
        )

    obstruction = {
        "u2": corrected_obstruction(first),
        "uv": corrected_obstruction(cross),
        "v2": corrected_obstruction(second),
    }

    def cubic_obstruction(direction: Sequence[Fraction]) -> Fraction:
        second_order = quadratic(direction)
        correction_two = solve_square_over_q(
            pivot_matrix, [-second_order[row] for row in rows]
        ) + [Q(0), Q(0)]
        jets = [
            OrderThreeJet.seed(Field(), value, linear, quadratic_value)
            for value, linear, quadratic_value in zip(seed, direction, correction_two)
        ]
        cubic_terms = [jet.cubic for jet in residuals_from_jets(jets)]
        correction_three = solve_square_over_q(
            pivot_matrix, [-cubic_terms[row] for row in rows]
        )
        return cubic_terms[6] + sum(
            coefficient * value
            for coefficient, value in zip(seventh[:6], correction_three)
        )

    cubic_directions = {
        "u": directions[0],
        "v": directions[1],
        "u_plus_v": [left + right for left, right in zip(directions[0], directions[1])],
        "u_plus_2v": [left + 2 * right for left, right in zip(directions[0], directions[1])],
    }
    cubic = {name: cubic_obstruction(direction) for name, direction in cubic_directions.items()}
    forms = {
        name: {"u2": str(left), "uv": str(middle), "v2": str(right)}
        for name, left, middle, right in zip(EQUATION_NAMES, first, cross, second)
    }
    return {
        "free_tangent_coordinates": [VARIABLES[index] for index in free_columns],
        "quadratic_forms": forms,
        "all_quadratic_terms_vanish": all(
            value == 0 for values in zip(first, cross, second) for value in values
        ),
        "residual_E0_2_after_solving_the_other_six_to_second_order": {
            key: str(value) for key, value in obstruction.items()
        },
        "quadratic_obstruction_vanishes": all(value == 0 for value in obstruction.values()),
        "cubic_obstruction_samples": {key: str(value) for key, value in cubic.items()},
        "all_cubic_obstruction_samples_vanish": all(value == 0 for value in cubic.values()),
    }


def formal_transverse_line(
    order: int = 8, direction: tuple[int, int] = (1, 0)
) -> dict[str, object]:
    """Solve six equations in Q[[t]] along a chosen line in (x02,x12)."""

    if order < 1:
        raise ValueError("formal order must be positive")
    seed = SEED_MODULI + SECTIONS[0] + SECTIONS[1]
    jacobian = [list(jet.gradient) for jet in residuals(seed, Field())]
    rows = list(range(6))
    pivot_columns = list(range(6))
    pivot_matrix = [[jacobian[row][column] for column in pivot_columns] for row in rows]
    coefficients = [[Q(0)] * (order + 1) for _ in VARIABLES]
    for index, value in enumerate(seed):
        coefficients[index][0] = value
    coefficients[6][1], coefficients[7][1] = map(Q, direction)
    field = Field(series_order=order)
    seventh_coefficients = []
    for degree in range(1, order + 1):
        jets = [FormalSeries.seed(field, row) for row in coefficients]
        values = residuals_from_jets(jets)
        correction = solve_square_over_q(
            pivot_matrix, [-values[row].coefficients[degree] for row in rows]
        )
        for column, value in zip(pivot_columns, correction):
            coefficients[column][degree] = value
        checked = residuals_from_jets([FormalSeries.seed(field, row) for row in coefficients])
        if any(checked[row].coefficients[degree] != 0 for row in rows):
            raise AssertionError("formal implicit solve failed for a selected equation")
        seventh_coefficients.append(checked[6].coefficients[degree])
    return {
        "parameter": (
            f"x02 = 583/115 + ({direction[0]})*t; "
            f"x12 = 13/23 + ({direction[1]})*t"
        ),
        "order": order,
        "all_E0_2_coefficients_through_order_vanish": all(
            value == 0 for value in seventh_coefficients
        ),
        "E0_2_coefficients": [str(value) for value in seventh_coefficients],
        "coordinate_series": {
            name: [str(value) for value in row]
            for name, row in zip(VARIABLES, coefficients)
        },
    }


def small_pade_models(
    coordinate_series: dict[str, Sequence[str]], max_degree: int = 5
) -> list[dict[str, object]]:
    """Recognize only rigorously checked low-degree rational functions in t."""

    models = []
    for name, rendered in coordinate_series.items():
        if name in {"x02", "x12"}:
            continue
        series = list(map(Q, rendered))
        for numerator_degree in range(max_degree + 1):
            for denominator_degree in range(1, max_degree + 1):
                equations = []
                right = []
                for degree in range(numerator_degree + 1, len(series)):
                    equations.append(
                        [
                            series[degree - offset] if degree >= offset else Q(0)
                            for offset in range(1, denominator_degree + 1)
                        ]
                    )
                    right.append(-series[degree])
                # Require at least one coefficient beyond those used to solve
                # the denominator.  Equality only constructs a Padé
                # interpolant and is not recognition evidence.
                if len(equations) <= denominator_degree:
                    continue
                try:
                    denominator_tail = solve_square_over_q(
                        equations[:denominator_degree], right[:denominator_degree]
                    )
                except AssertionError:
                    continue
                if all(
                    sum(left * right for left, right in zip(row, denominator_tail)) == target
                    for row, target in zip(equations, right)
                ):
                    models.append(
                        {
                            "coordinate": name,
                            "numerator_degree": numerator_degree,
                            "denominator_degree": denominator_degree,
                        }
                    )
    return models


def formal_bivariate_germ(
    order: int = 5, include_coordinates: bool = False
) -> dict[str, object]:
    """Solve the six-row implicit system in Q[[u,v]] and audit E0_2."""

    if order < 1:
        raise ValueError("formal order must be positive")
    seed = SEED_MODULI + SECTIONS[0] + SECTIONS[1]
    jacobian = [list(jet.gradient) for jet in residuals(seed, Field())]
    rows = list(range(6))
    pivot_columns = list(range(6))
    pivot_matrix = [[jacobian[row][column] for column in pivot_columns] for row in rows]
    coefficients = [{(0, 0): value} for value in seed]
    coefficients[6][(1, 0)] = Q(1)
    coefficients[7][(0, 1)] = Q(1)
    field = Field(series_order=order)
    seventh = {}
    for degree in range(1, order + 1):
        for u_degree in range(degree + 1):
            monomial = (u_degree, degree - u_degree)
            values = residuals_from_jets(
                [FormalBivariate.seed(field, row) for row in coefficients]
            )
            correction = solve_square_over_q(
                pivot_matrix,
                [-values[row].coefficient(monomial) for row in rows],
            )
            for column, value in zip(pivot_columns, correction):
                coefficients[column][monomial] = value
            checked = residuals_from_jets(
                [FormalBivariate.seed(field, row) for row in coefficients]
            )
            if any(checked[row].coefficient(monomial) != 0 for row in rows):
                raise AssertionError("bivariate implicit solve failed for a selected equation")
            seventh[monomial] = checked[6].coefficient(monomial)
    result = {
        "parameters": "x02=583/115+u; x12=13/23+v",
        "order": order,
        "all_E0_2_coefficients_through_total_order_vanish": all(
            value == 0 for value in seventh.values()
        ),
        "E0_2_coefficients": {
            f"u^{u_degree}v^{v_degree}": str(value)
            for (u_degree, v_degree), value in sorted(seventh.items())
        },
    }
    if include_coordinates:
        result["coordinate_coefficients"] = {
            name: {
                f"u^{u_degree}v^{v_degree}": str(value)
                for (u_degree, v_degree), value in sorted(row.items())
            }
            for name, row in zip(VARIABLES, coefficients)
        }
    return result


def pair_affine_intersection_audit() -> dict[str, object]:
    """Audit the finite affine intersection of the two pinned seed sections."""

    collision_t = Q(-529, 60)
    ordinate_one = (
        -Q(44352, 2645) * collision_t**3
        + Q(94026768, 330625) * collision_t**2
        - Q(18751906008, 8265625) * collision_t
        + Q(87599424, 15625)
    )
    ordinate_two = (
        Q(19008, 2645) * collision_t**3
        - Q(3274128, 330625) * collision_t**2
        - Q(585453528, 8265625) * collision_t
        + Q(19541808, 78125)
    )
    if ordinate_one == ordinate_two or ordinate_one == -ordinate_two:
        raise AssertionError("the selected sections unexpectedly meet in the affine chart")
    return {
        "unique_abscissa_collision_parameter_T": str(collision_t),
        "y1_minus_y2_at_collision": str(ordinate_one - ordinate_two),
        "y1_plus_y2_at_collision": str(ordinate_one + ordinate_two),
        "finite_affine_intersection_count": 0,
        "scope_limit": (
            "does not include the infinity fiber, component corrections, or a Shioda pairing"
        ),
    }


def all_companion_pair_tangent_audit() -> dict[str, object]:
    """Determine whether a different pair avoids the seed's tangent excess."""

    def ranks_at(
        moduli: tuple[Fraction, ...], companions: tuple[tuple[Fraction, Fraction], ...]
    ) -> dict[str, int]:
        answer = {}
        for left in range(len(companions)):
            for right in range(left + 1, len(companions)):
                values = residuals(moduli + companions[left] + companions[right], Field())
                if any(value.value != 0 for value in values):
                    raise AssertionError("a displayed companion pair is not an exact residual zero")
                rank, _ = row_reduce([list(value.gradient) for value in values], Field())
                answer[f"{left + 1},{right + 1}"] = rank
        return answer

    ranks = ranks_at(SEED_MODULI, (
        (Q(-95, 23), Q(37, 23)), (Q(-95, 23), Q(-37, 23)),
        (Q(583, 115), Q(13, 23)), (Q(583, 115), Q(-13, 23)),
        (Q(3444, 575), Q(7, 23)), (Q(3444, 575), Q(-7, 23)),
    ))
    comparison_ranks = ranks_at(COMPARISON_MODULI, COMPARISON_SECTIONS)
    if set(ranks.values()) != {6} or set(comparison_ranks.values()) != {6}:
        raise AssertionError("a companion pair escaped the common rank-six tangent pattern")
    return {
        "pair_count": len(ranks),
        "jacobian_rank_histogram": {"6": len(ranks)},
        "comparison_seed_roots": [0, 23, 93, 128, 133, 175],
        "comparison_seed_pair_count": len(comparison_ranks),
        "comparison_seed_jacobian_rank_histogram": {"6": len(comparison_ranks)},
        "interpretation": (
            "every pair at both displayed six-companion seeds has the same two-dimensional "
            "tangent space; neither seed can supply the naively transverse rank-seven pair"
        ),
    }


def rational_mod(value: Fraction, modulus: int) -> int:
    return (value.numerator % modulus) * pow(value.denominator % modulus, -1, modulus) % modulus


def rational_reconstruction(value: int, modulus: int) -> Fraction | None:
    """Return a conventional low-height reconstruction, if one is certified."""

    bound = isqrt(modulus // 2)
    r0, r1 = modulus, value % modulus
    t0, t1 = 0, 1
    while abs(r1) > bound:
        quotient = r0 // r1
        r0, r1 = r1, r0 - quotient * r1
        t0, t1 = t1, t0 - quotient * t1
    if (
        t1 == 0
        or abs(t1) > bound
        or gcd(r1, t1) != 1
        or (r1 - value * t1) % modulus
    ):
        return None
    return Q(r1, t1)


def selected_minor(matrix: Sequence[Sequence[int]], prime: int) -> tuple[list[int], list[int]]:
    """Use c4 and x12 as free coordinates when their six-by-six minor works."""

    preferred_free = ((3, 7),) + tuple(
        free for free in itertools.combinations(range(len(VARIABLES)), 2) if free != (3, 7)
    )
    for free in preferred_free:
        columns = [index for index in range(len(VARIABLES)) if index not in free]
        for rows in itertools.combinations(range(7), 6):
            minor = [[matrix[row][column] for column in columns] for row in rows]
            if row_reduce(minor, Field(prime))[0] == 6:
                return list(rows), columns
    raise AssertionError("no two-dimensional Hensel slice is nonsingular")


def lift_selected_rows(
    coords: Sequence[int], *, prime: int, exponent: int, rows: Sequence[int], columns: Sequence[int],
    free_digits: Sequence[int], jacobian: Sequence[Sequence[int]],
) -> list[int] | None:
    """Lift six independent rows, then retain only a full seven-row solution."""

    next_modulus = prime ** (exponent + 1)
    trial = list(coords)
    free = [index for index in range(len(VARIABLES)) if index not in columns]
    for index, digit in zip(free, free_digits):
        trial[index] = (trial[index] + prime**exponent * digit) % next_modulus
    values = residuals(trial, Field(next_modulus, tangent=False))
    if any(int(values[row].value) % (prime**exponent) for row in rows):
        raise AssertionError("attempted to lift a nonsolution")
    minor = [[jacobian[row][column] for column in columns] for row in rows]
    right = [(-int(values[row].value) // (prime**exponent)) % prime for row in rows]
    correction = solve_square(minor, right, prime)
    for column, digit in zip(columns, correction):
        trial[column] = (trial[column] + prime**exponent * digit) % next_modulus
    return trial if all(
        jet.value == 0 for jet in residuals(trial, Field(next_modulus, tangent=False))
    ) else None


def hensel_witness(prime: int, precision: int) -> dict[str, object]:
    """Continue a nonzero tangent direction while auditing the seventh row."""

    seed = SEED_MODULI + SECTIONS[0] + SECTIONS[1]
    coords = [rational_mod(value, prime) for value in seed]
    jacobian = [list(map(int, jet.gradient)) for jet in residuals(coords, Field(prime))]
    rows, columns = selected_minor(jacobian, prime)
    first_lifts = []
    for c4_digit in range(prime):
        for x12_digit in range(prime):
            lifted = lift_selected_rows(
                coords, prime=prime, exponent=1, rows=rows, columns=columns,
                free_digits=(c4_digit, x12_digit), jacobian=jacobian,
            )
            if lifted is not None:
                first_lifts.append(((c4_digit, x12_digit), lifted))
    # The next lift is where a redundant differential can become a nonlinear
    # obstruction.  Free digits do not affect that obstruction at this order,
    # so zero digits suffice to test each first-order tangent direction.
    survivors = []
    for digits, first in first_lifts:
        lifted = lift_selected_rows(
            first, prime=prime, exponent=2, rows=rows, columns=columns,
            free_digits=(0, 0), jacobian=jacobian,
        )
        if lifted is not None:
            survivors.append((digits, lifted))
    nonzero = next((item for item in survivors if item[0] != (0, 0)), None)
    if nonzero is None:
        raise AssertionError("no nontrivial tangent direction survived to p^3")
    digits, current = nonzero
    exponent = 3
    while exponent < precision:
        next_lift = None
        for c4_digit in range(prime):
            for x12_digit in range(prime):
                next_lift = lift_selected_rows(
                    current, prime=prime, exponent=exponent, rows=rows, columns=columns,
                    free_digits=(c4_digit, x12_digit), jacobian=jacobian,
                )
                if next_lift is not None:
                    break
            if next_lift is not None:
                break
        if next_lift is None:
            raise AssertionError(f"chosen nonzero branch stopped at p^{exponent + 1}")
        current = next_lift
        exponent += 1
    modulus = prime**precision
    free = [index for index in range(len(VARIABLES)) if index not in columns]
    reconstructions = [rational_reconstruction(value, modulus) for value in current]
    exact_reconstruction = False
    if all(value is not None for value in reconstructions):
        exact_reconstruction = all(
            jet.value == 0 for jet in residuals(reconstructions, Field())
        )
    return {
        "prime": prime,
        "precision": precision,
        "modulus": modulus,
        "free_coordinates": [VARIABLES[index] for index in free],
        "independent_equation_rows": [EQUATION_NAMES[index] for index in rows],
        "first_order_lift_count": len(first_lifts),
        "second_order_surviving_tangent_direction_count": len(survivors),
        "second_order_all_tangent_directions_survive": len(survivors) == prime * prime,
        "chosen_nonzero_tangent_direction": list(digits),
        "coordinates_mod_prime_power": dict(zip(VARIABLES, map(int, current))),
        "all_seven_residuals_zero_mod_prime_power": True,
        "low_height_rational_reconstruction": {
            "all_coordinates_reconstructed": all(value is not None for value in reconstructions),
            "exactly_satisfies_the_seven_Q_equations": exact_reconstruction,
        },
    }


def run(
    primes: Iterable[int], precision: int, bivariate_order: int | None = None,
    line_order: int = 12, pade_degree: int = 5, line_direction: tuple[int, int] = (1, 0),
) -> dict[str, object]:
    seed = SEED_MODULI + SECTIONS[0] + SECTIONS[1]
    exact = residuals(seed, Field())
    if any(jet.value != 0 for jet in exact):
        raise AssertionError("the pinned two-section seed is not on the incidence")
    jacobian = [list(jet.gradient) for jet in exact]
    exact_rank, exact_pivots = row_reduce(jacobian, Field())
    if exact_rank != 6:
        raise AssertionError(f"expected exact rank 6, found {exact_rank}")
    quadratic_data = tangent_quadratic_data()
    if not quadratic_data["quadratic_obstruction_vanishes"]:
        raise AssertionError("the seventh residual has a nonzero quadratic obstruction")
    formal_data = formal_transverse_line(line_order, direction=line_direction)
    if not formal_data["all_E0_2_coefficients_through_order_vanish"]:
        raise AssertionError("the seventh residual obstructs the formal transverse line")
    pade_models = small_pade_models(formal_data["coordinate_series"], max_degree=pade_degree)
    intersection_audit = pair_affine_intersection_audit()
    all_pair_audit = all_companion_pair_tangent_audit()
    modular = []
    witnesses = []
    for prime in primes:
        reduced = residuals([rational_mod(value, prime) for value in seed], Field(prime))
        rank, pivots = row_reduce([list(jet.gradient) for jet in reduced], Field(prime))
        if rank != 6:
            raise AssertionError(f"rank dropped modulo {prime}")
        modular.append({"prime": prime, "jacobian_rank": rank, "pivot_columns": [VARIABLES[index] for index in pivots]})
        witnesses.append(hensel_witness(prime, precision))
    result = {
        "status": "two-dimensional tangent space and bounded nontrivial p-adic continuations recorded",
        "seed": {
            "normalized_moduli": [str(value) for value in SEED_MODULI],
            "sections": [[str(value) for value in section] for section in SECTIONS],
        },
        "equations": "M plus (E2,E1,E0) for each section, evaluated recursively as jets",
        "expanded_residuals_materialized": False,
        "exact_jacobian_rank": exact_rank,
        "exact_pivot_columns": [VARIABLES[index] for index in exact_pivots],
        "tangent_space_dimension": len(VARIABLES) - exact_rank,
        "exact_second_order_obstruction": {
            "free_tangent_coordinates": quadratic_data["free_tangent_coordinates"],
            "residual_E0_2_after_solving_the_other_six_to_second_order": quadratic_data[
                "residual_E0_2_after_solving_the_other_six_to_second_order"
            ],
            "quadratic_obstruction_vanishes": quadratic_data[
                "quadratic_obstruction_vanishes"
            ],
        },
        "exact_third_order_obstruction": {
            "samples": quadratic_data["cubic_obstruction_samples"],
            "all_spanning_samples_vanish": quadratic_data[
                "all_cubic_obstruction_samples_vanish"
            ],
        },
        "characteristic_zero_formal_transverse_line": {
            "parameter": formal_data["parameter"],
            "order": formal_data["order"],
            "all_E0_2_coefficients_through_order_vanish": formal_data[
                "all_E0_2_coefficients_through_order_vanish"
            ],
            "low_degree_pade_search": {
                "maximum_numerator_and_denominator_degree": pade_degree,
                "recognized_models": pade_models,
            },
        },
        "pair_affine_intersection_audit": intersection_audit,
        "all_companion_pair_tangent_audit": all_pair_audit,
        "finite_field_tangent_checks": modular,
        "hensel_continuations": witnesses,
        "not_established": [
            "the actual local dimension (the tangent excess must still be resolved)",
            "rational parameterization or a global plane model",
            "pair intersections or a Shioda Gram matrix",
            "saturation or independence from the generic rank-13 subgroup",
            "generic rank at least 14",
        ],
    }
    if bivariate_order is not None:
        bivariate = formal_bivariate_germ(bivariate_order)
        if not bivariate["all_E0_2_coefficients_through_total_order_vanish"]:
            raise AssertionError("the seventh residual obstructs the bivariate formal germ")
        result["characteristic_zero_formal_bivariate_germ"] = {
            "parameters": bivariate["parameters"],
            "order": bivariate["order"],
            "all_E0_2_coefficients_through_total_order_vanish": bivariate[
                "all_E0_2_coefficients_through_total_order_vanish"
            ],
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", nargs="+", type=int, default=DEFAULT_PRIMES)
    parser.add_argument("--precision", type=int, default=4)
    parser.add_argument(
        "--bivariate-order", type=int,
        help="also solve the full Q[[u,v]] germ through this total degree",
    )
    parser.add_argument("--line-order", type=int, default=12)
    parser.add_argument("--pade-degree", type=int, default=5)
    parser.add_argument("--line-direction", nargs=2, type=int, default=(1, 0))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(
        args.primes, args.precision, args.bivariate_order,
        line_order=args.line_order, pade_degree=args.pade_degree,
        line_direction=tuple(args.line_direction),
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
