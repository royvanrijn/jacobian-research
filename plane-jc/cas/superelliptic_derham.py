#!/usr/bin/env python3
"""Exact Hermite/de Rham reduction for superelliptic weighted Wronskians.

The basic identity is

    d(D/y**b) = (a*A*D' - b*A'*D)/(a*y**(a+b)) dt,
    y**a = A(t).

The reducer works character by character for the action y -> zeta*y.  It is
intended as a small exact-algebra kernel for Newton-chain leading blocks, not
as a general-purpose algebraic-curves package.
"""

from dataclasses import dataclass
from math import gcd
from typing import Iterable

import sympy as sp


@dataclass(frozen=True)
class ExactTerm:
    """A primitive term whose differential was removed.

    ``kind == "pole"`` represents numerator/y**exponent.  ``kind ==
    "infinity"`` represents numerator*y**exponent.
    """

    kind: str
    numerator: sp.Expr
    exponent: int


@dataclass(frozen=True)
class ReductionResult:
    character: int
    exact_terms: tuple[ExactTerm, ...]
    remainder: sp.Expr
    affine_coordinates: tuple[sp.Expr, ...]
    residue: sp.Expr
    compact_coordinates: tuple[sp.Expr, ...] | None

    @property
    def is_exact(self) -> bool:
        return all(value == 0 for value in self.affine_coordinates)


@dataclass(frozen=True)
class QuotientDifferentialResult:
    """The invariant-character differential descended to the rational base."""

    character: int
    differential: sp.Expr

    @property
    def is_exact(self) -> bool:
        # Exactness on the punctured rational quotient requires an ordinary
        # rational Hermite/residue calculation, kept separate from the
        # nontrivial superelliptic eigenspace reducer.
        return self.differential == 0


class SuperellipticDeRham:
    """Hermite reduction on the smooth model of ``y**a = A(t)``.

    ``A`` must be squarefree and have positive degree.  Nontrivial characters
    are reduced in superelliptic eigenspaces.  Character zero is returned as
    an ordinary rational differential on the quotient coordinate ``t``.
    """

    def __init__(
        self, t: sp.Symbol, A: sp.Expr, a: int, *, check_squarefree: bool = True
    ):
        if a < 2:
            raise ValueError("the covering exponent a must be at least two")
        self.t = t
        self.a = int(a)
        self.A = sp.Poly(A, t)
        if self.A.degree() < 1:
            raise ValueError("A must have positive degree")
        if check_squarefree and sp.gcd(self.A, self.A.diff()).degree() != 0:
            raise ValueError("A must be squarefree")
        self.n = self.A.degree()
        self.delta = gcd(self.a, self.n)

    @property
    def genus(self) -> int:
        return ((self.a - 1) * (self.n - 1) - (self.delta - 1)) // 2

    @property
    def affine_h1_dimension(self) -> int:
        return (self.a - 1) * (self.n - 1)

    @property
    def compact_h1_dimension(self) -> int:
        return 2 * self.genus

    def character_dimension(self, character: int, compact: bool = True) -> int:
        r = character % self.a
        if r == 0:
            return 0
        if not compact:
            return self.n - 1
        return self.n - 2 if self._has_infinity_residue(r) else self.n - 1

    def _has_infinity_residue(self, r: int) -> bool:
        # The d points over infinity carry the augmentation representation of
        # mu_d.  A residue occurs precisely when a/d divides the character.
        return r % (self.a // self.delta) == 0

    def residue_coefficients(self, character: int) -> tuple[sp.Expr, ...]:
        """Return the normalized residue functional on t^i dt/y^r.

        The tuple has length n-1.  It is zero if this character has no
        logarithmic direction at infinity.  Otherwise its pivot coefficient
        is one.  Vanishing of its pairing with the affine remainder is the
        second-kind condition needed to descend to compact de Rham cohomology.
        """

        r = character % self.a
        if r == 0:
            raise ValueError("character zero is not handled")
        zeros = [sp.S.Zero] * (self.n - 1)
        if not self._has_infinity_residue(r):
            return tuple(zeros)

        e = self.a // self.delta
        N = self.n // self.delta
        pivot = (r // e) * N - 1
        z = sp.Dummy("z")
        lc = self.A.LC()
        unit = sp.expand(self.A.as_expr().subs(self.t, 1 / z) * z**self.n / lc)
        order = self.n - 1 - pivot
        inverse_power = sp.series(
            unit ** sp.Rational(-r, self.a), z, 0, order
        ).removeO()
        for i in range(pivot, self.n - 1):
            zeros[i] = sp.expand(inverse_power).coeff(z, i - pivot)
        assert zeros[pivot] == 1
        return tuple(zeros)

    def compact_basis(self, character: int) -> tuple[sp.Expr, ...]:
        """Canonical polynomial numerators for the compact character basis."""

        r = character % self.a
        residues = self.residue_coefficients(r)
        if not any(residues):
            return tuple(self.t**i for i in range(self.n - 1))
        pivot = next(i for i, value in enumerate(residues) if value != 0)
        return tuple(
            self.t**i - residues[i] * self.t**pivot
            for i in range(self.n - 1)
            if i != pivot
        )

    def reduce(
        self, numerator: sp.Expr, denominator_exponent: int
    ) -> ReductionResult | QuotientDifferentialResult:
        """Reduce ``numerator*dt/y**denominator_exponent`` exactly.

        The output remainder is in the affine basis t^i dt/y^r with
        0 <= i <= n-2.  If its infinity residue vanishes, compact coordinates
        are also returned (the pivot/logarithmic coordinate is omitted).
        """

        m = int(denominator_exponent)
        if m <= 0:
            raise ValueError("the denominator exponent must be positive")
        r = m % self.a
        if r == 0:
            quotient_power = m // self.a
            return QuotientDifferentialResult(
                character=0,
                differential=sp.cancel(numerator / self.A.as_expr() ** quotient_power),
            )

        P = sp.Poly(sp.cancel(numerator), self.t)
        Aprime = self.A.diff()
        inverse_aprime = sp.invert(Aprime, self.A)
        exact_terms: list[ExactTerm] = []

        # Finite-pole Hermite reduction.  At a root of A, choose S so that
        # P + (m-a)/a A'S vanishes; the quotient has denominator y^(m-a).
        while m > self.a:
            k = m - self.a
            S = (-sp.Rational(self.a, k) * P * inverse_aprime).rem(self.A)
            lifted = P - self.A * S.diff() + sp.Rational(k, self.a) * Aprime * S
            P, remainder = sp.div(lifted, self.A)
            if not remainder.is_zero:
                raise ArithmeticError("finite-pole Hermite division was not exact")
            exact_terms.append(ExactTerm("pole", S.as_expr(), k))
            m = k

        # Polynomial reduction at infinity using
        # d(Q*y^(a-m)) = (A Q' + (a-m)/a A'Q) dt/y^m.
        infinity_primitive = sp.Poly(0, self.t)
        while not P.is_zero and P.degree() >= self.n - 1:
            j = P.degree() - self.n + 1
            monomial = sp.Poly(self.t**j, self.t)
            operator = self.A * monomial.diff() + sp.Rational(
                self.a - m, self.a
            ) * Aprime * monomial
            factor = P.LC() / operator.LC()
            P = P - factor * operator
            infinity_primitive += factor * monomial
        if not infinity_primitive.is_zero:
            exact_terms.append(
                ExactTerm("infinity", infinity_primitive.as_expr(), self.a - m)
            )

        coords = tuple(sp.cancel(P.nth(i)) for i in range(self.n - 1))
        residue_coeffs = self.residue_coefficients(r)
        residue = sp.cancel(sum(c * x for c, x in zip(residue_coeffs, coords)))
        if residue == 0:
            if any(residue_coeffs):
                pivot = next(i for i, value in enumerate(residue_coeffs) if value != 0)
                compact = tuple(value for i, value in enumerate(coords) if i != pivot)
            else:
                compact = coords
        else:
            compact = None
        return ReductionResult(
            character=r,
            exact_terms=tuple(exact_terms),
            remainder=P.as_expr(),
            affine_coordinates=coords,
            residue=residue,
            compact_coordinates=compact,
        )

    def reduce_weighted_wronskian(
        self, R: sp.Expr, b: int
    ) -> ReductionResult | QuotientDifferentialResult:
        """Reduce the differential attached to a*A*D' - b*A'D = R."""

        return self.reduce(sp.Rational(1, self.a) * R, self.a + int(b))


def weighted_wronskian_compatibility(
    t: sp.Symbol,
    A: sp.Expr,
    R: sp.Expr,
    a: int,
    b: int,
    D_coefficients: Iterable[sp.Symbol],
    D_exponents: Iterable[int],
) -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Expr, ...]]:
    """Triangularly solve a support-constrained weighted Wronskian block.

    Returns the solved D coefficients and the residual compatibility
    equations.  No division by parameter expressions is performed: a pivot is
    accepted only when its coefficient is a nonzero constant.
    """

    variables = tuple(D_coefficients)
    exponents = tuple(D_exponents)
    if len(variables) != len(exponents):
        raise ValueError("D coefficients and exponents must have equal length")
    D = sum(value * t**power for value, power in zip(variables, exponents))
    equation = sp.Poly(
        sp.expand(a * A * sp.diff(D, t) - b * sp.diff(A, t) * D - R), t
    )
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    residual: list[sp.Expr] = []
    for degree in range(max(0, equation.degree()) + 1):
        coefficient = sp.cancel(equation.nth(degree).subs(substitutions))
        if coefficient == 0:
            continue
        live = [variable for variable in variables if coefficient.has(variable)]
        solved = False
        if len(live) == 1 and sp.Poly(coefficient, live[0]).degree() == 1:
            variable = live[0]
            pivot = sp.diff(coefficient, variable)
            if not pivot.has(*variables) and not pivot.has(*A.free_symbols - {t}):
                substitutions[variable] = sp.cancel(
                    -(coefficient - pivot * variable) / pivot
                )
                solved = True
        if not solved:
            residual.append(sp.factor(coefficient))
    return substitutions, tuple(residual)


if __name__ == "__main__":
    t = sp.symbols("t")
    A = t**8 + 2 * t**7 - t**3 + t + 1
    curve = SuperellipticDeRham(t, A, 2)
    D = t**12 - 3 * t**7 + 2 * t**2
    R = sp.expand(2 * A * sp.diff(D, t) - 3 * sp.diff(A, t) * D)
    result = curve.reduce_weighted_wronskian(R, 3)
    assert result.is_exact
    print("genus", curve.genus)
    print("character dimensions (affine, compact)", (7, curve.character_dimension(1)))
    print("weighted-Wronskian exactness", result.is_exact)
