"""Exact data for the global degree-seven Davenport--Sunada pair.

The two covers are the conjugate Davenport polynomials over
``Q(a)/(a^2+a+2)``.  This module deliberately keeps the symbolic cover data
and the finite Fano-plane permutation data together: the former identifies
the common branch base, while the latter certifies equality of every
unramified fiber zeta function.
"""
from __future__ import annotations

from itertools import product

import sympy as sp


Y, Z, T, U, A = sp.symbols("Y Z T U a")
MINIMAL_A = A**2 + A + 2


def reduce_a(expression: sp.Expr, *parameters: sp.Symbol) -> sp.Expr:
    """Reduce a rational expression modulo ``a^2+a+2``.

    The denominator may involve ``a``; it is inverted in the quadratic
    coefficient field.  ``parameters`` generate the surrounding rational
    function field.
    """
    domain = sp.QQ.frac_field(*parameters) if parameters else sp.QQ
    numerator, denominator = sp.fraction(sp.cancel(expression))
    modulus = sp.Poly(MINIMAL_A, A, domain=domain)
    numerator_poly = sp.Poly(sp.expand(numerator), A, domain=domain)
    denominator_poly = sp.Poly(sp.expand(denominator), A, domain=domain)
    inverse_denominator = sp.invert(denominator_poly, modulus)
    return sp.rem(numerator_poly * inverse_denominator, modulus).as_expr()


def conjugate_a(expression: sp.Expr, *parameters: sp.Symbol) -> sp.Expr:
    """Apply the nontrivial automorphism ``a -> -1-a``."""
    return reduce_a(expression.subs(A, -1 - A), *parameters)


def davenport_polynomial(variable: sp.Symbol = Y) -> sp.Expr:
    """Return Cassou-Nogues--Couveignes' degree-seven polynomial ``g_T``."""
    return (
        sp.Rational(1, 7) * variable**7
        + (1 + A) * T * variable**5
        + (1 + A) * T * variable**4
        - (3 - 2 * A) * T**2 * variable**3
        - 2 * (1 - 2 * A) * T**2 * variable**2
        - sp.Rational(1, 28)
        * (5 + 3 * A)
        * (28 * T - 2 - 11 * A)
        * T**2
        * variable
        - (1 + A) * T**3
    )


def davenport_pair() -> tuple[sp.Expr, sp.Expr]:
    """Return ``g_T(Y)`` and its conjugate ``h_T(Z)``."""
    g = davenport_polynomial(Y)
    h = conjugate_a(davenport_polynomial(Z), T, Z)
    return g, h


def correspondence() -> sp.Expr:
    """Return the cubic component of ``g_T(Y)-h_T(Z)=0``."""
    return (
        Y**3
        + A * Y**2 * Z
        + Y**2 * Z
        + A * Y * Z**2
        - Z**3
        + (5 + 3 * A) * T * Y
        + (-2 + 3 * A) * T * Z
        + (2 * A + 1) * T
    )


def branch_cubic() -> sp.Expr:
    """Return the common reduced finite-branch divisor in the ``(T,U)`` base."""
    return (
        17920 * T**10
        + 10472 * T**9
        - 2464 * T**8
        - 1792 * T**7 * U
        - 1728 * T**7
        + 11956 * T**6 * U
        + 7056 * T**5 * U
        - 4802 * T**3 * U**2
        + 343 * U**3
    )


Matrix = tuple[int, ...]
Permutation = tuple[int, ...]


def matrix_vector(matrix: Matrix, vector: int) -> int:
    """Apply a 3-by-3 binary matrix to a nonzero vector encoded by 1..7."""
    bits = tuple((vector >> column) & 1 for column in range(3))
    output = tuple(
        sum(matrix[3 * row + column] * bits[column] for column in range(3)) % 2
        for row in range(3)
    )
    return output[0] | (output[1] << 1) | (output[2] << 2)


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        sum(left[3 * row + k] * right[3 * k + column] for k in range(3)) % 2
        for row in range(3)
        for column in range(3)
    )


def transpose(matrix: Matrix) -> Matrix:
    return tuple(matrix[3 * column + row] for row in range(3) for column in range(3))


def gl32() -> tuple[Matrix, ...]:
    """Enumerate ``GL(3,2)`` exactly."""
    matrices = []
    for entries in product((0, 1), repeat=9):
        images = {matrix_vector(entries, vector) for vector in range(1, 8)}
        if len(images) == 7 and 0 not in images:
            matrices.append(entries)
    return tuple(matrices)


IDENTITY: Matrix = (1, 0, 0, 0, 1, 0, 0, 0, 1)


def inverses(matrices: tuple[Matrix, ...]) -> dict[Matrix, Matrix]:
    return {
        matrix: next(
            candidate
            for candidate in matrices
            if matrix_product(matrix, candidate) == IDENTITY
        )
        for matrix in matrices
    }


def point_permutation(matrix: Matrix) -> Permutation:
    """Natural action on the seven nonzero vectors."""
    return tuple(matrix_vector(matrix, vector) - 1 for vector in range(1, 8))


def line_permutation(matrix: Matrix, inverse: Matrix) -> Permutation:
    """Contragredient action on the seven nonzero covectors."""
    dual = transpose(inverse)
    return tuple(matrix_vector(dual, covector) - 1 for covector in range(1, 8))


def cycle_partition(permutation: Permutation) -> tuple[int, ...]:
    """Return the sorted Frobenius orbit lengths of a permutation."""
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def zeta_denominator(permutation: Permutation) -> tuple[tuple[int, int], ...]:
    """Encode ``prod_d (1-u^d)^multiplicity`` for a finite etale fiber."""
    counts: dict[int, int] = {}
    for length in cycle_partition(permutation):
        counts[length] = counts.get(length, 0) + 1
    return tuple(sorted(counts.items()))
