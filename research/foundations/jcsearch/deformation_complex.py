"""Exact linear algebra for three-term deformation complexes.

The common bounded calculation in several parts of the repository is

    gauge directions -> corrections -> relation defects.

This module provides the field-level core of that calculation.  Polynomial
parameter rings, Fitting strata, saturation, completion, and conductor
descent remain separate module/scheme calculations; specializing such a
family to an exact field point produces an instance of ``ThreeTermComplex``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterable

import sympy as sp


def total_degree_monomials(
    variables: tuple[sp.Symbol, ...], maximum_degree: int
) -> tuple[sp.Expr, ...]:
    """Return monomials ordered first by total degree."""

    if not variables:
        return (sp.Integer(1),)
    return tuple(
        sp.prod(
            variable**exponent
            for variable, exponent in zip(variables, exponents)
        )
        for total_degree in range(maximum_degree + 1)
        for exponents in product(
            range(total_degree + 1), repeat=len(variables)
        )
        if sum(exponents) == total_degree
    )


def polynomial_coordinate_column(
    polynomial: sp.Expr,
    basis: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    """Write a polynomial as a column in a declared monomial basis."""

    expanded = sp.Poly(sp.expand(polynomial), *variables)
    support = set(expanded.monoms())
    basis_support = {
        sp.Poly(monomial, *variables).monoms()[0] for monomial in basis
    }
    outside = support - basis_support
    if outside:
        raise ValueError(f"polynomial has monomials outside the basis: {outside}")
    return sp.Matrix(
        [expanded.coeff_monomial(monomial) for monomial in basis]
    )


def polynomial_action_matrix(
    domain_basis: tuple[sp.Expr, ...],
    codomain_basis: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    action,
) -> sp.Matrix:
    """Build the matrix of a polynomial linear differential operator."""

    return sp.Matrix.hstack(
        *(
            polynomial_coordinate_column(
                action(element), codomain_basis, variables
            )
            for element in domain_basis
        )
    )


def _column(vector: sp.MatrixBase, dimension: int, name: str) -> sp.Matrix:
    result = sp.Matrix(vector)
    if result.shape == (1, dimension):
        result = result.T
    if result.shape != (dimension, 1):
        raise ValueError(
            f"{name} must have shape ({dimension}, 1), got {result.shape}"
        )
    return result


def rank_mod_prime(matrix: sp.MatrixBase, prime: int) -> int:
    """Return the rank of a rational matrix over ``GF(prime)``.

    Denominators divisible by ``prime`` are rejected rather than silently
    treated as zero.  This makes bad specializations explicit.
    """

    if prime <= 1 or not sp.ntheory.primetest.isprime(prime):
        raise ValueError(f"{prime} is not prime")

    rows = []
    for row in sp.Matrix(matrix).tolist():
        modular_row = []
        for value in row:
            rational = sp.Rational(value)
            numerator = int(rational.p) % prime
            denominator = int(rational.q) % prime
            if denominator == 0:
                raise ValueError(
                    f"denominator of {value} vanishes modulo {prime}"
                )
            modular_row.append(
                numerator * pow(denominator, -1, prime) % prime
            )
        rows.append(modular_row)

    row_count = len(rows)
    column_count = len(rows[0]) if rows else sp.Matrix(matrix).cols
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if rows[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [
            entry * inverse % prime for entry in rows[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            multiplier = rows[row][column] % prime
            if multiplier:
                rows[row] = [
                    (entry - multiplier * pivot_entry) % prime
                    for entry, pivot_entry in zip(
                        rows[row], rows[pivot_row]
                    )
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


@dataclass(frozen=True)
class ThreeTermComplex:
    """A finite exact complex ``C0 -> C1 -> C2`` over a field.

    Matrices use column-vector conventions.  Thus ``d0`` has shape
    ``(dim C1, dim C0)`` and ``d1`` has shape ``(dim C2, dim C1)``.
    """

    d0: sp.MatrixBase
    d1: sp.MatrixBase
    name: str = "deformation complex"

    def __post_init__(self) -> None:
        d0 = sp.Matrix(self.d0)
        d1 = sp.Matrix(self.d1)
        if d0.rows != d1.cols:
            raise ValueError(
                "incompatible middle dimensions: "
                f"d0 has {d0.rows} rows but d1 has {d1.cols} columns"
            )
        object.__setattr__(self, "d0", d0)
        object.__setattr__(self, "d1", d1)
        composite = d1 * d0
        if composite != sp.zeros(d1.rows, d0.cols):
            raise ValueError(f"{self.name}: d1*d0 is not zero")

    @property
    def dimensions(self) -> tuple[int, int, int]:
        return self.d0.cols, self.d0.rows, self.d1.rows

    @property
    def ranks(self) -> tuple[int, int]:
        return self.d0.rank(), self.d1.rank()

    @property
    def cohomology_dimensions(self) -> tuple[int, int, int]:
        """Return dimensions of ``H0``, ``H1``, and ``H2``."""

        c0, c1, c2 = self.dimensions
        rank_d0, rank_d1 = self.ranks
        return (
            c0 - rank_d0,
            c1 - rank_d0 - rank_d1,
            c2 - rank_d1,
        )

    def dual_obstruction_cocycles(self) -> tuple[sp.Matrix, ...]:
        """Return a basis of functionals annihilating ``im(d1)``.

        Each returned column has length ``dim C2`` and acts on a defect by
        transpose multiplication.  These vectors identify ``H2.dual``.
        """

        return tuple(sp.Matrix(vector) for vector in self.d1.T.nullspace())

    def obstruction_coordinates(
        self, defect: sp.MatrixBase
    ) -> tuple[sp.Expr, ...]:
        """Evaluate a defect on a basis of dual obstruction cocycles."""

        vector = _column(defect, self.d1.rows, "defect")
        return tuple(
            sp.expand((functional.T * vector)[0])
            for functional in self.dual_obstruction_cocycles()
        )

    def defect_is_correctable(self, defect: sp.MatrixBase) -> bool:
        """Test whether ``d1(correction) = -defect`` is solvable."""

        return all(
            coordinate == 0
            for coordinate in self.obstruction_coordinates(defect)
        )

    def correction_solution(self, defect: sp.MatrixBase) -> sp.Matrix:
        """Return one correction solving ``d1(correction) = -defect``.

        Free parameters in the affine solution torsor are set to zero.  Its
        translation space is available as ``d1.nullspace()``.
        """

        vector = _column(defect, self.d1.rows, "defect")
        if not self.defect_is_correctable(vector):
            raise ValueError(f"{self.name}: defect has a nonzero H2 class")
        solution, parameters = self.d1.gauss_jordan_solve(-vector)
        zero_parameters = {
            parameter: sp.Integer(0) for parameter in parameters
        }
        return sp.Matrix(solution.subs(zero_parameters))

    def prime_rank_profile(
        self, primes: Iterable[int]
    ) -> dict[int, tuple[int, int]]:
        """Return ``(rank d0, rank d1)`` at each good prime."""

        return {
            prime: (
                rank_mod_prime(self.d0, prime),
                rank_mod_prime(self.d1, prime),
            )
            for prime in primes
        }


@dataclass(frozen=True)
class RelativeFiberConnectionComplex:
    """The bounded canonical fiber-plus-connection complex and its raw data.

    In canonical coordinates ``(S,T,R)``, corrections are triples ``(s,t,a)``
    and raw defects are

        (s_S+t_T, s_R-a_T, t_R+a_S).

    The returned three-term complex uses coordinates in the kernel of the
    Bianchi map ``(F,G,H) -> F_R-G_S-H_T`` as its final module.  Keeping the
    raw matrices here makes boundary-lattice restrictions reproducible
    without rebuilding the differential in each census script.
    """

    complex: ThreeTermComplex
    raw_d1: sp.Matrix
    bianchi: sp.Matrix
    closed_basis: sp.Matrix


def relative_fiber_connection_complex(
    gauge_degree: int = 4,
) -> RelativeFiberConnectionComplex:
    """Construct the canonical bounded relative complex.

    The correction, defect, and Bianchi degrees are respectively one, two,
    and three less than ``gauge_degree``.  The default is the exact
    ``35 -> 60 -> 26`` complex used by the rank-two incidence programme.
    """

    if gauge_degree < 3:
        raise ValueError("gauge_degree must be at least three")

    S, T, R = sp.symbols("S T R")
    variables = (S, T, R)
    gauge_basis = total_degree_monomials(variables, gauge_degree)
    scalar_corrections = total_degree_monomials(
        variables, gauge_degree - 1
    )
    scalar_defects = total_degree_monomials(variables, gauge_degree - 2)
    bianchi_targets = total_degree_monomials(
        variables, gauge_degree - 3
    )

    def block_column(*columns: sp.Matrix) -> sp.Matrix:
        return sp.Matrix.vstack(*columns)

    gauge_columns = []
    for hamiltonian in gauge_basis:
        gauge_columns.append(
            block_column(
                polynomial_coordinate_column(
                    -sp.diff(hamiltonian, T),
                    scalar_corrections,
                    variables,
                ),
                polynomial_coordinate_column(
                    sp.diff(hamiltonian, S),
                    scalar_corrections,
                    variables,
                ),
                polynomial_coordinate_column(
                    -sp.diff(hamiltonian, R),
                    scalar_corrections,
                    variables,
                ),
            )
        )
    d0 = sp.Matrix.hstack(*gauge_columns)

    raw_defect_columns = []
    for component in range(3):
        for correction in scalar_corrections:
            if component == 0:
                defects = (
                    sp.diff(correction, S),
                    sp.diff(correction, R),
                    sp.Integer(0),
                )
            elif component == 1:
                defects = (
                    sp.diff(correction, T),
                    sp.Integer(0),
                    sp.diff(correction, R),
                )
            else:
                defects = (
                    sp.Integer(0),
                    -sp.diff(correction, T),
                    sp.diff(correction, S),
                )
            raw_defect_columns.append(
                block_column(
                    *(
                        polynomial_coordinate_column(
                            defect,
                            scalar_defects,
                            variables,
                        )
                        for defect in defects
                    )
                )
            )
    raw_d1 = sp.Matrix.hstack(*raw_defect_columns)

    bianchi_blocks = (
        polynomial_action_matrix(
            scalar_defects,
            bianchi_targets,
            variables,
            lambda polynomial: sp.diff(polynomial, R),
        ),
        polynomial_action_matrix(
            scalar_defects,
            bianchi_targets,
            variables,
            lambda polynomial: -sp.diff(polynomial, S),
        ),
        polynomial_action_matrix(
            scalar_defects,
            bianchi_targets,
            variables,
            lambda polynomial: -sp.diff(polynomial, T),
        ),
    )
    bianchi = sp.Matrix.hstack(*bianchi_blocks)
    if bianchi * raw_d1 != sp.zeros(bianchi.rows, raw_d1.cols):
        raise AssertionError("raw defects violate the Bianchi identity")

    closed_basis = sp.Matrix.hstack(*bianchi.nullspace())
    d1_closed, parameters = closed_basis.gauss_jordan_solve(raw_d1)
    if list(parameters):
        raise AssertionError("closed-defect coordinates are not unique")
    if closed_basis * d1_closed != raw_d1:
        raise AssertionError("closed-defect coordinate reconstruction failed")

    complex_ = ThreeTermComplex(
        d0,
        d1_closed,
        "relative fiber-plus-connection complex",
    )
    return RelativeFiberConnectionComplex(
        complex=complex_,
        raw_d1=raw_d1,
        bianchi=bianchi,
        closed_basis=closed_basis,
    )
