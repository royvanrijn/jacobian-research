#!/usr/bin/env python3
"""Exact rank-one constant-Hessian-pencil recognition for HC4 potentials.

For a four-variable potential ``psi`` with Hessian ``H``, a constant
covector ``ell`` gives the quadratic direction

    A = (ell.x)^2 / 2,       Hess(A) = ell*ell^T.

The rank-one determinant identity shows that this is a constant-Hessian
pencil direction exactly when

    ell^T adj(H) ell = 0

as a polynomial in the source variables.  Coefficient extraction gives a
finite homogeneous quadratic scheme in P^3.  Its four standard affine charts
decide nonemptiness over the algebraic closure by exact Groebner bases.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import sympy as sp


@dataclass(frozen=True)
class ConstantNullCovectorSystem:
    potential: sp.Expr
    variables: tuple[sp.Symbol, ...]
    covector: tuple[sp.Symbol, ...]
    hessian: sp.Matrix
    hessian_determinant: sp.Expr
    metric_numerator: sp.Expr
    coefficient_equations: tuple[sp.Expr, ...]

    @property
    def has_constant_nonzero_hessian_determinant(self) -> bool:
        return (
            self.hessian_determinant != 0
            and not (set(self.variables) & self.hessian_determinant.free_symbols)
        )


@dataclass(frozen=True)
class ProjectiveNullChart:
    chart_index: int
    normalized_covector: sp.Symbol
    is_empty: bool
    groebner_basis: tuple[sp.Expr, ...]


def constant_null_covector_system(
    potential: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    covector: Sequence[sp.Symbol] | None = None,
) -> ConstantNullCovectorSystem:
    """Compile the coefficient quadrics cutting out constant null covectors."""

    variables_tuple = tuple(variables)
    if len(variables_tuple) != 4 or len(set(variables_tuple)) != 4:
        raise ValueError("HC4 recognition requires four distinct variables")
    covector_tuple = tuple(covector or sp.symbols("ell0:4"))
    if len(covector_tuple) != 4 or len(set(covector_tuple)) != 4:
        raise ValueError("the covector requires four distinct symbols")

    hessian = sp.hessian(potential, variables_tuple)
    determinant = sp.factor(hessian.det(method="domain-ge"))
    adjugate = hessian.adjugate(method="domain-ge")
    ell = sp.Matrix(covector_tuple)
    metric_numerator = sp.expand((ell.T * adjugate * ell)[0])
    polynomial = sp.Poly(metric_numerator, *variables_tuple)
    equations = tuple(
        dict.fromkeys(
            sp.factor(coefficient)
            for coefficient in polynomial.coeffs()
            if coefficient != 0
        )
    )
    return ConstantNullCovectorSystem(
        potential=sp.expand(potential),
        variables=variables_tuple,
        covector=covector_tuple,
        hessian=hessian,
        hessian_determinant=determinant,
        metric_numerator=metric_numerator,
        coefficient_equations=equations,
    )


def projective_null_covector_charts(
    system: ConstantNullCovectorSystem,
) -> tuple[ProjectiveNullChart, ...]:
    """Decide the four standard projective charts over the algebraic closure."""

    if not system.has_constant_nonzero_hessian_determinant:
        raise ValueError("the potential does not have constant nonzero Hessian")

    charts = []
    for chart_index, normalized in enumerate(system.covector):
        remaining = tuple(
            variable
            for index, variable in enumerate(system.covector)
            if index != chart_index
        )
        equations = tuple(
            sp.expand(equation.subs(normalized, 1))
            for equation in system.coefficient_equations
        )
        if not equations:
            charts.append(
                ProjectiveNullChart(
                    chart_index=chart_index,
                    normalized_covector=normalized,
                    is_empty=False,
                    groebner_basis=(),
                )
            )
            continue
        basis = sp.groebner(equations, *remaining, order="grevlex")
        charts.append(
            ProjectiveNullChart(
                chart_index=chart_index,
                normalized_covector=normalized,
                is_empty=basis.contains(sp.Integer(1)),
                groebner_basis=tuple(
                    sp.factor(polynomial.as_expr()) for polynomial in basis.polys
                ),
            )
        )
    return tuple(charts)


def has_rank_one_pencil_over_algebraic_closure(
    system: ConstantNullCovectorSystem,
) -> bool:
    """Return whether the projective constant-null scheme is nonempty."""

    return any(not chart.is_empty for chart in projective_null_covector_charts(system))


def certify_constant_null_covector(
    system: ConstantNullCovectorSystem,
    values: Sequence[sp.Expr],
) -> bool:
    """Verify one explicit nonzero constant covector exactly."""

    values_tuple = tuple(values)
    if len(values_tuple) != 4 or all(value == 0 for value in values_tuple):
        return False
    substitution = dict(zip(system.covector, values_tuple, strict=True))
    return all(
        sp.expand(equation.subs(substitution)) == 0
        for equation in system.coefficient_equations
    )

