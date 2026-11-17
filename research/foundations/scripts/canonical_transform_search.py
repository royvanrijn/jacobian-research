#!/usr/bin/env python3
"""Shared exact primitives for polynomial canonical-transformation searches.

The HC4 and DC2 searches have different admission gates, but they use the
same low-level operations: constant Poisson brackets, exact Hamiltonian
shears, word composition, pullback, and symplectic verification.  This module
keeps those operations convention-explicit.  It does not decide whether a
map preserves a Hessian determinant, retains a collision, has a polynomial
inverse beyond the declared shear word, or admits a Weyl lift.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Sequence

import sympy as sp


@dataclass(frozen=True)
class MixedLineHamiltonian:
    """One reusable ``coefficient*(q_i+epsilon*p_j)^degree`` letter."""

    source_index: int
    dual_index: int
    epsilon: int
    degree: int
    coefficient: int
    hamiltonian: sp.Expr


def mixed_line_hamiltonians(
    q: Sequence[sp.Symbol],
    p: Sequence[sp.Symbol],
    *,
    degrees: Sequence[int],
    coefficients: Sequence[int],
    epsilons: Sequence[int] = (-1, 1),
) -> tuple[MixedLineHamiltonian, ...]:
    """Generate the common mixed source--dual Hamiltonian alphabet."""

    if len(q) != len(p) or not q:
        raise ValueError("q and p must contain the same positive pair count")
    if any(degree < 2 for degree in degrees):
        raise ValueError("mixed-line shear degrees must be at least two")
    if any(epsilon not in (-1, 1) for epsilon in epsilons):
        raise ValueError("epsilons must be -1 or 1")
    return tuple(
        MixedLineHamiltonian(
            source_index=source_index,
            dual_index=dual_index,
            epsilon=epsilon,
            degree=degree,
            coefficient=coefficient,
            hamiltonian=(
                coefficient
                * (q[source_index] + epsilon * p[dual_index]) ** degree
            ),
        )
        for source_index, dual_index, epsilon, degree, coefficient in product(
            range(len(q)),
            range(len(p)),
            epsilons,
            degrees,
            coefficients,
        )
    )


def canonical_poisson_matrix(
    pair_count: int,
    *,
    q_p_bracket: int = 1,
) -> sp.Matrix:
    """Return Pi with ``{q_i,p_j}=q_p_bracket*delta_ij``.

    Variables are ordered ``(q_1,...,q_n,p_1,...,p_n)``.  The explicit sign
    parameter accommodates the opposite bracket conventions in the HC4 and
    DC2 notes without a hidden coordinate permutation.
    """

    if pair_count <= 0:
        raise ValueError("pair_count must be positive")
    if q_p_bracket not in (-1, 1):
        raise ValueError("q_p_bracket must be -1 or 1")
    identity = q_p_bracket * sp.eye(pair_count)
    return sp.zeros(pair_count).row_join(identity).col_join(
        (-identity).row_join(sp.zeros(pair_count))
    )


def poisson_bracket(
    left: sp.Expr,
    right: sp.Expr,
    variables: Sequence[sp.Symbol],
    poisson_matrix: sp.Matrix,
) -> sp.Expr:
    """Evaluate the constant Poisson bracket exactly."""

    dimension = len(variables)
    if poisson_matrix.shape != (dimension, dimension):
        raise ValueError("Poisson matrix and variable dimensions disagree")
    left_gradient = sp.Matrix([sp.diff(left, variable) for variable in variables])
    right_gradient = sp.Matrix(
        [sp.diff(right, variable) for variable in variables]
    )
    return sp.expand((left_gradient.T * poisson_matrix * right_gradient)[0])


def hamiltonian_vector_field(
    hamiltonian: sp.Expr,
    variables: Sequence[sp.Symbol],
    poisson_matrix: sp.Matrix,
) -> tuple[sp.Expr, ...]:
    """Return ``Pi*grad(H)`` in the declared convention."""

    gradient = sp.Matrix(
        [sp.diff(hamiltonian, variable) for variable in variables]
    )
    return tuple(sp.expand(value) for value in poisson_matrix * gradient)


def exact_invariant_shear(
    hamiltonian: sp.Expr,
    variables: Sequence[sp.Symbol],
    poisson_matrix: sp.Matrix,
) -> tuple[sp.Expr, ...]:
    """Return the exact time-one map ``z -> z+Pi*grad(H)``.

    This finite formula is valid when every velocity component is invariant
    under the Hamiltonian flow.  The condition is verified rather than
    assumed, preventing accidental use of a truncated Lie series.
    """

    velocity = hamiltonian_vector_field(
        hamiltonian,
        variables,
        poisson_matrix,
    )
    for component in velocity:
        assert poisson_bracket(
            hamiltonian,
            component,
            variables,
            poisson_matrix,
        ) == 0
    return tuple(
        sp.expand(variable + component)
        for variable, component in zip(variables, velocity, strict=True)
    )


def compose_maps(
    outer: Sequence[sp.Expr],
    inner: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return ``outer o inner`` by simultaneous substitution."""

    if len(outer) != len(variables) or len(inner) != len(variables):
        raise ValueError("map and variable dimensions disagree")
    substitution = dict(zip(variables, inner, strict=True))
    return tuple(
        sp.expand(expression.subs(substitution, simultaneous=True))
        for expression in outer
    )


def pullback_polynomial(
    polynomial: sp.Expr,
    coordinate_map: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> sp.Expr:
    """Pull back one polynomial by a coordinate map."""

    if len(coordinate_map) != len(variables):
        raise ValueError("map and variable dimensions disagree")
    substitution = dict(zip(variables, coordinate_map, strict=True))
    return sp.expand(polynomial.subs(substitution, simultaneous=True))


def verify_symplectic_map(
    outputs: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    poisson_matrix: sp.Matrix,
    *,
    require_unit_determinant: bool = True,
) -> sp.Matrix:
    """Verify exact preservation of ``Pi`` and return the Jacobian."""

    if len(outputs) != len(variables):
        raise ValueError("map and variable dimensions disagree")
    jacobian = sp.Matrix(outputs).jacobian(variables)
    assert (jacobian * poisson_matrix * jacobian.T).applyfunc(
        sp.expand
    ) == poisson_matrix
    if require_unit_determinant:
        assert sp.expand(jacobian.det()) == 1
    return jacobian
