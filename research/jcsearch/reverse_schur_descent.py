"""Exact equation builders for reverse Schur descent.

The scalar family is

    Phi(t, x) = lambda*t**2/2 + t*A(x) + B(x).

The module deliberately separates three conditions which are easy to
conflate in coefficient searches:

* singularity (or a prescribed corank) of ``Hess(B+s*A)``;
* constancy of the full bordered Hessian determinant; and
* equality of two gradients after Schur reduction.

It also implements the simultaneous matrix-pivot identities.  The routines
only build and verify exact polynomial equations; they do not promote a
bounded coefficient search to a classification theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence

import sympy as sp


def _symbols(values: Sequence[sp.Symbol], *, name: str) -> tuple[sp.Symbol, ...]:
    result = tuple(values)
    if not result or any(not isinstance(value, sp.Symbol) for value in result):
        raise ValueError(f"{name} must be a nonempty sequence of symbols")
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must contain distinct symbols")
    return result


def coefficient_equations(
    expressions: Sequence[sp.Expr],
    polynomial_variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return the distinct nonzero coefficient equations of expressions.

    Symbols not listed in ``polynomial_variables`` remain coefficient
    parameters.  The order is deterministic and follows SymPy's monomial
    ordering expression by expression.
    """

    variables = _symbols(polynomial_variables, name="polynomial_variables")
    equations: list[sp.Expr] = []
    seen: set[sp.Expr] = set()
    for expression in expressions:
        polynomial = sp.Poly(sp.expand(expression), *variables)
        for _, coefficient in polynomial.terms():
            coefficient = sp.factor(coefficient)
            if coefficient != 0 and coefficient not in seen:
                seen.add(coefficient)
                equations.append(coefficient)
    return tuple(equations)


def hessian_integrability_residuals(
    matrix: sp.Matrix,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return symmetry and third-derivative residuals for a Hessian matrix.

    A polynomial matrix ``H`` is a Hessian, up to affine ambiguity in the
    potential, precisely when it is symmetric and each row one-form is
    closed.  The latter equations are

    ``d_k H_ij - d_j H_ik = 0``.
    """

    x = _symbols(variables, name="variables")
    if matrix.shape != (len(x), len(x)):
        raise ValueError("matrix size must match the variable count")
    residuals = [
        sp.expand(matrix[left, right] - matrix[right, left])
        for left in range(len(x))
        for right in range(left + 1, len(x))
    ]
    residuals.extend(
        sp.expand(
            sp.diff(matrix[row, left], x[right])
            - sp.diff(matrix[row, right], x[left])
        )
        for row in range(len(x))
        for left in range(len(x))
        for right in range(left + 1, len(x))
    )
    return tuple(residuals)


def adjugate_divergence_residuals(
    matrix: sp.Matrix,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return the row divergences of ``adj(matrix)``.

    For an exact Hessian these vanish by the polynomial Piola identity.
    Keeping the identity executable is useful when a reverse search begins
    with a corank-one matrix rather than an already integrated potential.
    """

    x = _symbols(variables, name="variables")
    if matrix.shape != (len(x), len(x)):
        raise ValueError("matrix size must match the variable count")
    adjugate = matrix.adjugate(method="domain-ge")
    return tuple(
        sp.expand(
            sum(
                sp.diff(adjugate[row, column], x[column])
                for column in range(len(x))
            )
        )
        for row in range(len(x))
    )


def kernel_line_piola_residuals(
    kernel_vector: sp.Matrix,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return ``D_v(v) + div(v)*v`` for a polynomial kernel vector.

    If a corank-one Hessian has ``adj(H)=epsilon*v*v.T`` with ``epsilon`` an
    ``x``-constant unit, the adjugate Piola identity forces these residuals
    to vanish.
    """

    x = _symbols(variables, name="variables")
    if kernel_vector.shape != (len(x), 1):
        raise ValueError("kernel-vector size must match the variable count")
    jacobian = kernel_vector.jacobian(x)
    divergence = sp.trace(jacobian)
    return tuple(
        sp.expand(entry)
        for entry in jacobian * kernel_vector + divergence * kernel_vector
    )


def rank_at_most_equations(
    matrix: sp.Matrix,
    rank_bound: int,
) -> tuple[sp.Expr, ...]:
    """Return all minors imposing ``rank(matrix) <= rank_bound``."""

    rows, columns = matrix.shape
    if rank_bound < 0:
        raise ValueError("rank_bound must be nonnegative")
    minor_size = rank_bound + 1
    if minor_size > min(rows, columns):
        return ()
    return tuple(
        sp.factor(matrix.extract(row_set, column_set).det(method="domain-ge"))
        for row_set in combinations(range(rows), minor_size)
        for column_set in combinations(range(columns), minor_size)
    )


def corank_one_adjugate_scalar(
    matrix: sp.Matrix,
    kernel_vector: sp.Matrix,
) -> sp.Expr:
    """Certify ``adj(matrix) = q*v*v.T`` and return the exact scalar ``q``.

    This is the rank-one adjugate factor used by the scalar reverse-Schur
    classifier.  A ``ValueError`` is raised if the supplied vector is not a
    kernel vector or if the adjugate has a different rank-one factorization.
    """

    rows, columns = matrix.shape
    if rows != columns or kernel_vector.shape != (rows, 1):
        raise ValueError("matrix and kernel-vector sizes are incompatible")
    if any(sp.factor(entry) != 0 for entry in matrix * kernel_vector):
        raise ValueError("the supplied vector is not in the kernel")

    adjugate = matrix.adjugate(method="domain-ge")
    scalar: sp.Expr | None = None
    for left in range(rows):
        for right in range(rows):
            product = kernel_vector[left] * kernel_vector[right]
            if product != 0:
                scalar = sp.cancel(adjugate[left, right] / product)
                break
        if scalar is not None:
            break
    if scalar is None:
        raise ValueError("the kernel vector is zero")
    expected = kernel_vector * kernel_vector.T * scalar
    if any(
        sp.factor(adjugate[left, right] - expected[left, right]) != 0
        for left in range(rows)
        for right in range(rows)
    ):
        raise ValueError("the adjugate does not factor through this kernel")
    if scalar == 0:
        raise ValueError("the matrix does not have generic corank one")
    return sp.factor(scalar)


@dataclass(frozen=True)
class ScalarPivotSchurFamily:
    """A scalar reverse-Schur family ``t*A+B`` in exact coordinates."""

    variables: tuple[sp.Symbol, ...]
    pivot_variable: sp.Symbol
    a: sp.Expr
    b: sp.Expr

    def __post_init__(self) -> None:
        variables = _symbols(self.variables, name="variables")
        if self.pivot_variable in variables:
            raise ValueError("the pivot variable must be separate from x")

    @property
    def gradient_a(self) -> sp.Matrix:
        return sp.Matrix([sp.diff(self.a, variable) for variable in self.variables])

    def pencil_potential(self, parameter: sp.Expr) -> sp.Expr:
        return sp.expand(self.b + parameter * self.a)

    def pencil_hessian(self, parameter: sp.Expr) -> sp.Matrix:
        return sp.hessian(self.pencil_potential(parameter), self.variables)

    def bordered_hessian(self, parameter: sp.Expr) -> sp.Matrix:
        gradient = self.gradient_a
        pencil = self.pencil_hessian(parameter)
        return sp.Matrix.vstack(
            sp.Matrix.hstack(sp.zeros(1, 1), gradient.T),
            sp.Matrix.hstack(gradient, pencil),
        )

    def bordered_determinant(self, parameter: sp.Expr) -> sp.Expr:
        return sp.factor(self.bordered_hessian(parameter).det(method="domain-ge"))

    def parent_potential(self, pivot_coefficient: sp.Expr) -> sp.Expr:
        t = self.pivot_variable
        return sp.expand(pivot_coefficient * t**2 / 2 + t * self.a + self.b)

    def parent_hessian_determinant(self, pivot_coefficient: sp.Expr) -> sp.Expr:
        variables = (self.pivot_variable, *self.variables)
        return sp.factor(
            sp.hessian(self.parent_potential(pivot_coefficient), variables).det(
                method="domain-ge"
            )
        )

    def schur_descendant(self, kappa: sp.Expr, mu: sp.Expr) -> sp.Expr:
        """Return ``B+kappa*A^2/2+mu*A``.

        For a genuine quadratic parent with pivot coefficient ``lambda`` and
        fixed pivot-gradient value ``y``, take ``kappa=-1/lambda`` and
        ``mu=y/lambda``.
        """

        return sp.expand(self.b + kappa * self.a**2 / 2 + mu * self.a)

    def quadratic_pivot_reduction(
        self,
        pivot_coefficient: sp.Expr,
        pivot_gradient_value: sp.Expr,
    ) -> sp.Expr:
        return sp.expand(
            self.b
            - (self.a - pivot_gradient_value) ** 2 / (2 * pivot_coefficient)
        )

    def quadratic_schur_identity_residual(
        self,
        pivot_coefficient: sp.Expr,
        pivot_gradient_value: sp.Expr,
    ) -> sp.Expr:
        """Return the exact block-Schur determinant residual."""

        critical_pivot = (
            pivot_gradient_value - self.a
        ) / pivot_coefficient
        parent_determinant = self.parent_hessian_determinant(pivot_coefficient)
        reduced = self.quadratic_pivot_reduction(
            pivot_coefficient,
            pivot_gradient_value,
        )
        reduced_determinant = sp.hessian(reduced, self.variables).det(
            method="domain-ge"
        )
        return sp.factor(
            parent_determinant.subs(self.pivot_variable, critical_pivot)
            - pivot_coefficient * reduced_determinant
        )

    def singular_pencil_equations(self, parameter: sp.Symbol) -> tuple[sp.Expr, ...]:
        return coefficient_equations(
            (self.pencil_hessian(parameter).det(method="domain-ge"),),
            (parameter, *self.variables),
        )

    def constant_parent_equations(
        self,
        pivot_coefficient: sp.Expr,
        target_constant: sp.Expr,
    ) -> tuple[sp.Expr, ...]:
        return coefficient_equations(
            (
                self.parent_hessian_determinant(pivot_coefficient)
                - target_constant,
            ),
            (self.pivot_variable, *self.variables),
        )

    def collision_equations(
        self,
        kappa: sp.Expr,
        mu: sp.Expr,
        left_point: Sequence[sp.Expr],
        right_point: Sequence[sp.Expr],
    ) -> tuple[sp.Expr, ...]:
        """Return exact equal-gradient equations for two reduced points."""

        if len(left_point) != len(self.variables) or len(right_point) != len(
            self.variables
        ):
            raise ValueError("collision points must match the variable count")
        descendant = self.schur_descendant(kappa, mu)
        gradient = [sp.diff(descendant, variable) for variable in self.variables]
        left = dict(zip(self.variables, left_point, strict=True))
        right = dict(zip(self.variables, right_point, strict=True))
        return tuple(
            sp.factor(coordinate.subs(left) - coordinate.subs(right))
            for coordinate in gradient
        )


@dataclass(frozen=True)
class MatrixPivotSchurFamily:
    """Simultaneous matrix-pivot reverse-Schur equation builder."""

    variables: tuple[sp.Symbol, ...]
    pivot_variables: tuple[sp.Symbol, ...]
    a: tuple[sp.Expr, ...]
    b: sp.Expr
    pivot_matrix: sp.Matrix

    def __post_init__(self) -> None:
        variables = _symbols(self.variables, name="variables")
        pivots = _symbols(self.pivot_variables, name="pivot_variables")
        if set(variables) & set(pivots):
            raise ValueError("pivot and reduced variables must be disjoint")
        if len(self.a) != len(pivots):
            raise ValueError("one A_i is required per pivot")
        if self.pivot_matrix.shape != (len(pivots), len(pivots)):
            raise ValueError("pivot matrix has the wrong size")
        if self.pivot_matrix != self.pivot_matrix.T:
            raise ValueError("pivot matrix must be symmetric")
        if sp.factor(self.pivot_matrix.det(method="domain-ge")) == 0:
            raise ValueError("pivot matrix must be invertible")

    @property
    def jacobian_a(self) -> sp.Matrix:
        return sp.Matrix(
            [
                [sp.diff(expression, variable) for variable in self.variables]
                for expression in self.a
            ]
        )

    def pencil_potential(self, parameters: Sequence[sp.Expr]) -> sp.Expr:
        if len(parameters) != len(self.a):
            raise ValueError("one pencil parameter is required per pivot")
        return sp.expand(
            self.b
            + sum(
                parameter * expression
                for parameter, expression in zip(parameters, self.a, strict=True)
            )
        )

    def pencil_hessian(self, parameters: Sequence[sp.Expr]) -> sp.Matrix:
        return sp.hessian(self.pencil_potential(parameters), self.variables)

    @property
    def parent_potential(self) -> sp.Expr:
        t = sp.Matrix(self.pivot_variables)
        a = sp.Matrix(self.a)
        return sp.expand((t.T * self.pivot_matrix * t)[0] / 2 + (t.T * a)[0] + self.b)

    def reduced_potential(self, pivot_gradient_values: Sequence[sp.Expr]) -> sp.Expr:
        if len(pivot_gradient_values) != len(self.a):
            raise ValueError("one fixed pivot-gradient value is required per pivot")
        a_minus_y = sp.Matrix(self.a) - sp.Matrix(pivot_gradient_values)
        inverse = self.pivot_matrix.inv()
        return sp.expand(self.b - (a_minus_y.T * inverse * a_minus_y)[0] / 2)

    def schur_identity_residual(
        self,
        pivot_gradient_values: Sequence[sp.Expr],
    ) -> sp.Expr:
        """Return ``det Hess(Phi)|critical - det(Lambda)*det Hess(psi)``."""

        if len(pivot_gradient_values) != len(self.a):
            raise ValueError("one fixed pivot-gradient value is required per pivot")
        critical = self.pivot_matrix.inv() * (
            sp.Matrix(pivot_gradient_values) - sp.Matrix(self.a)
        )
        parent_variables = (*self.pivot_variables, *self.variables)
        parent_determinant = sp.hessian(
            self.parent_potential,
            parent_variables,
        ).det(method="domain-ge")
        substitution = dict(zip(self.pivot_variables, critical, strict=True))
        reduced_determinant = sp.hessian(
            self.reduced_potential(pivot_gradient_values),
            self.variables,
        ).det(method="domain-ge")
        return sp.factor(
            parent_determinant.subs(substitution, simultaneous=True)
            - self.pivot_matrix.det(method="domain-ge") * reduced_determinant
        )

    def pencil_corank_equations(
        self,
        parameters: Sequence[sp.Symbol],
        corank_at_least: int,
    ) -> tuple[sp.Expr, ...]:
        if len(parameters) != len(self.a):
            raise ValueError("one pencil parameter is required per pivot")
        if not 1 <= corank_at_least <= len(self.variables):
            raise ValueError("invalid requested corank")
        rank_bound = len(self.variables) - corank_at_least
        minors = rank_at_most_equations(
            self.pencil_hessian(parameters),
            rank_bound,
        )
        return coefficient_equations(
            minors,
            (*parameters, *self.variables),
        )
