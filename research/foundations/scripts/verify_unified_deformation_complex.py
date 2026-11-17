#!/usr/bin/env python3
"""Regression examples for the unified three-term deformation complex.

The examples instantiate

    gauge directions -> corrections -> relation defects

for a commutative left--right slice, an exact symplectic pair, an exact Weyl
pair, and a constant-Hessian constraint.  They verify only the common
field-level linear algebra.  Parameter components, higher Kuranishi maps,
boundary saturation, and conductor descent are deliberately outside this
small regression.
"""

from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.deformation_complex import ThreeTermComplex  # noqa: E402


def polynomial_basis(
    variables: tuple[sp.Symbol, ...], maximum_degree: int
) -> tuple[sp.Expr, ...]:
    """Return monomials in total-degree order."""

    if len(variables) == 1:
        return tuple(variables[0] ** degree for degree in range(maximum_degree + 1))
    if len(variables) != 2:
        raise ValueError("the regression helper supports one or two variables")
    left, right = variables
    return tuple(
        left**left_degree * right ** (total_degree - left_degree)
        for total_degree in range(maximum_degree + 1)
        for left_degree in range(total_degree + 1)
    )


def coordinate_column(
    polynomial: sp.Expr,
    basis: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    expanded = sp.Poly(sp.expand(polynomial), *variables)
    return sp.Matrix(
        [expanded.coeff_monomial(monomial) for monomial in basis]
    )


def matrix_from_action(
    domain: tuple[sp.Expr, ...],
    codomain: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    action,
) -> sp.Matrix:
    return sp.Matrix.hstack(
        *(
            coordinate_column(action(element), codomain, variables)
            for element in domain
        )
    )


def commutative_left_right_complex() -> ThreeTermComplex:
    """A determinant-one left--right slice at the identity of A1."""

    (x,) = sp.symbols("x,")
    corrections = polynomial_basis((x,), 2)
    defects = polynomial_basis((x,), 1)

    # Constant source and target translations are the admissible
    # divergence-free infinitesimal left--right gauges in this slice.
    d0 = sp.Matrix([[1, 1], [0, 0], [0, 0]])
    d1 = matrix_from_action(
        corrections,
        defects,
        (x,),
        lambda correction: sp.diff(correction, x),
    )
    return ThreeTermComplex(d0, d1, "commutative left-right")


def exact_symplectic_complex() -> ThreeTermComplex:
    """Hamiltonians -> pair corrections -> first symplectic defect."""

    q, p = sp.symbols("q p")
    hamiltonians = polynomial_basis((q, p), 3)
    scalar_corrections = polynomial_basis((q, p), 2)
    defects = polynomial_basis((q, p), 1)

    correction_columns = []
    for hamiltonian in hamiltonians:
        left = -sp.diff(hamiltonian, p)
        right = sp.diff(hamiltonian, q)
        correction_columns.append(
            sp.Matrix.vstack(
                coordinate_column(left, scalar_corrections, (q, p)),
                coordinate_column(right, scalar_corrections, (q, p)),
            )
        )
    d0 = sp.Matrix.hstack(*correction_columns)

    defect_columns = []
    for side in range(2):
        for correction in scalar_corrections:
            defect = (
                sp.diff(correction, q)
                if side == 0
                else sp.diff(correction, p)
            )
            defect_columns.append(
                coordinate_column(defect, defects, (q, p))
            )
    d1 = sp.Matrix.hstack(*defect_columns)
    return ThreeTermComplex(d0, d1, "exact symplectic")


def weyl_multiply(
    left: tuple[int, int], right: tuple[int, int]
) -> dict[tuple[int, int], sp.Integer]:
    """Multiply q^a p^b q^c p^d with [p,q]=1 in PBW order."""

    a, b = left
    c, d = right
    result: dict[tuple[int, int], sp.Integer] = {}
    for contractions in range(min(b, c) + 1):
        coefficient = (
            sp.factorial(contractions)
            * comb(b, contractions)
            * comb(c, contractions)
        )
        exponent = (
            a + c - contractions,
            b + d - contractions,
        )
        result[exponent] = result.get(exponent, 0) + coefficient
    return result


def weyl_commutator(
    left: dict[tuple[int, int], sp.Expr],
    right: dict[tuple[int, int], sp.Expr],
) -> dict[tuple[int, int], sp.Expr]:
    def multiply(
        first: dict[tuple[int, int], sp.Expr],
        second: dict[tuple[int, int], sp.Expr],
    ) -> dict[tuple[int, int], sp.Expr]:
        result: dict[tuple[int, int], sp.Expr] = {}
        for left_exponent, left_coefficient in first.items():
            for right_exponent, right_coefficient in second.items():
                for exponent, coefficient in weyl_multiply(
                    left_exponent, right_exponent
                ).items():
                    result[exponent] = (
                        result.get(exponent, 0)
                        + left_coefficient * right_coefficient * coefficient
                    )
        return {
            exponent: sp.expand(coefficient)
            for exponent, coefficient in result.items()
            if coefficient != 0
        }

    left_right = multiply(left, right)
    right_left = multiply(right, left)
    support = set(left_right) | set(right_left)
    return {
        exponent: sp.expand(
            left_right.get(exponent, 0) - right_left.get(exponent, 0)
        )
        for exponent in support
        if left_right.get(exponent, 0) != right_left.get(exponent, 0)
    }


def weyl_basis(maximum_degree: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (q_degree, total_degree - q_degree)
        for total_degree in range(maximum_degree + 1)
        for q_degree in range(total_degree + 1)
    )


def weyl_coordinate_column(
    element: dict[tuple[int, int], sp.Expr],
    basis: tuple[tuple[int, int], ...],
) -> sp.Matrix:
    return sp.Matrix([element.get(exponent, 0) for exponent in basis])


def exact_weyl_complex() -> ThreeTermComplex:
    """Inner gauges -> PBW corrections -> exact commutator defects."""

    hamiltonians = weyl_basis(3)
    scalar_corrections = weyl_basis(2)
    defects = weyl_basis(1)
    q = {(1, 0): sp.Integer(1)}
    p = {(0, 1): sp.Integer(1)}

    gauge_columns = []
    for exponent in hamiltonians:
        hamiltonian = {exponent: sp.Integer(1)}
        q_variation = weyl_commutator(hamiltonian, q)
        p_variation = weyl_commutator(hamiltonian, p)
        gauge_columns.append(
            sp.Matrix.vstack(
                weyl_coordinate_column(q_variation, scalar_corrections),
                weyl_coordinate_column(p_variation, scalar_corrections),
            )
        )
    d0 = sp.Matrix.hstack(*gauge_columns)

    defect_columns = []
    for side in range(2):
        for exponent in scalar_corrections:
            correction = {exponent: sp.Integer(1)}
            # The relation is [P,Q]=1, so the first variation is
            # [delta P,Q]+[P,delta Q].
            defect = (
                weyl_commutator(p, correction)
                if side == 0
                else weyl_commutator(correction, q)
            )
            defect_columns.append(
                weyl_coordinate_column(defect, defects)
            )
    d1 = sp.Matrix.hstack(*defect_columns)
    return ThreeTermComplex(d0, d1, "exact Weyl")


def hessian_determinant_complex() -> ThreeTermComplex:
    """A restricted cubic-correction constant-Hessian complex."""

    x, y = sp.symbols("x y")
    gauges = polynomial_basis((x, y), 1)
    corrections = polynomial_basis((x, y), 3)
    defects = polynomial_basis((x, y), 2)

    d0 = sp.Matrix.hstack(
        *(
            coordinate_column(gauge, corrections, (x, y))
            for gauge in gauges
        )
    )
    d1 = matrix_from_action(
        corrections,
        defects,
        (x, y),
        lambda correction: (
            sp.diff(correction, x, 2)
            + sp.diff(correction, y, 2)
        ),
    )
    return ThreeTermComplex(d0, d1, "Hessian determinant")


def main() -> None:
    complexes = (
        commutative_left_right_complex(),
        exact_symplectic_complex(),
        exact_weyl_complex(),
        hessian_determinant_complex(),
    )

    for complex_ in complexes:
        assert complex_.prime_rank_profile((101, 103, 107)) == {
            prime: complex_.ranks for prime in (101, 103, 107)
        }
        for gauge_column in range(complex_.d0.cols):
            defect = complex_.d1 * complex_.d0[:, gauge_column]
            assert complex_.defect_is_correctable(defect)
        print(
            "PASS:",
            complex_.name,
            f"dims={complex_.dimensions}",
            f"ranks={complex_.ranks}",
            f"cohomology={complex_.cohomology_dimensions}",
        )

    symplectic = complexes[1]
    weyl = complexes[2]
    assert symplectic.ranks == weyl.ranks
    assert symplectic.cohomology_dimensions == weyl.cohomology_dimensions

    # A quadratic defect cannot be removed by the restricted cubic
    # potential-correction space: its Laplacian has degree at most one.
    hessian = complexes[3]
    quadratic_defect = sp.Matrix([0, 0, 0, 0, 0, 1])
    assert not hessian.defect_is_correctable(quadratic_defect)
    assert any(
        coordinate != 0
        for coordinate in hessian.obstruction_coordinates(quadratic_defect)
    )
    correctable_defect = hessian.d1[:, 5]
    correction = hessian.correction_solution(correctable_defect)
    assert hessian.d1 * correction == -correctable_defect
    print("PASS: the canonical-pair symplectic and exact PBW complexes agree")
    print("PASS: a restricted Hessian quadratic defect survives in H2")
    print("PASS: a correctable defect has an exact particular correction")
    print("SCOPE: field fibers only; no component or global descent theorem")


if __name__ == "__main__":
    main()
