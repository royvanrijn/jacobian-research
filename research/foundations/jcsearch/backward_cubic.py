"""Exact primitives for backward cubic reduction.

The forward rank-compressed homogenization of

    F(x) = x + Q(x) + B c(x)

introduces companion variables ``y`` and a homogenizing variable ``t``.
On the invariant slice ``t=1`` the companion block cancels by polynomial
left--right equivalence.  This module makes that reverse operation, linear
fixed-covector slices, and pairwise collision bookkeeping first-class.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Sequence

import sympy as sp


@dataclass(frozen=True)
class CompanionCancellation:
    """Exact data for cancelling one rank-factor companion block."""

    variables: tuple[sp.Symbol, ...]
    companion_variables: tuple[sp.Symbol, ...]
    base_map: tuple[sp.Expr, ...]
    specialized_map: tuple[sp.Expr, ...]
    source_shear: tuple[sp.Expr, ...]
    factored_map: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class ParametricCompanionCancellation:
    """Relative cancellation over the homogenizing line."""

    variables: tuple[sp.Symbol, ...]
    companion_variables: tuple[sp.Symbol, ...]
    parameter: sp.Symbol
    base_map: tuple[sp.Expr, ...]
    scaled_base_family: tuple[sp.Expr, ...]
    parent_map: tuple[sp.Expr, ...]
    source_shear: tuple[sp.Expr, ...]
    factored_map: tuple[sp.Expr, ...]
    special_fiber: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class FixedCovectorRestriction:
    """Restriction of a map to one fixed affine-linear level."""

    covector: sp.Matrix
    level: sp.Expr
    pivot: int
    variables: tuple[sp.Symbol, ...]
    embedding: tuple[sp.Expr, ...]
    restricted_map: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class BackwardTerminalProfile:
    """The two dimension objectives attached to a quadratic--cubic endpoint."""

    base_dimension: int
    cubic_output_rank: int

    @property
    def homogeneous_dimension(self) -> int:
        return self.base_dimension + self.cubic_output_rank + 1

    @property
    def direct_cubic_key(self) -> tuple[int, int, int]:
        """Optimize arbitrary degree-three dimension before homogenization."""

        return (
            self.base_dimension,
            self.homogeneous_dimension,
            self.cubic_output_rank,
        )

    @property
    def homogeneous_key(self) -> tuple[int, int, int]:
        """Optimize the rank-compressed cubic-homogeneous endpoint."""

        return (
            self.homogeneous_dimension,
            self.base_dimension,
            self.cubic_output_rank,
        )


def coefficient_matrix(
    expressions: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> sp.Matrix:
    """Return component rows against the union of polynomial monomials."""

    polynomials = [
        sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
        for expression in expressions
    ]
    monomials = sorted(
        {
            exponents
            for polynomial in polynomials
            for exponents, coefficient in polynomial.terms()
            if coefficient
        }
    )
    if not monomials:
        return sp.zeros(len(expressions), 0)
    return sp.Matrix(
        [
            [polynomial.coeff_monomial(exponents) for exponents in monomials]
            for polynomial in polynomials
        ]
    )


def fixed_linear_covectors(
    corrections: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Matrix, ...]:
    """Compute all ``ell`` with ``ell(corrections)=0``."""

    return tuple(coefficient_matrix(corrections, variables).T.nullspace())


def collision_compatible_covectors(
    covectors: Sequence[sp.Matrix], points: Sequence[sp.Matrix]
) -> tuple[sp.Matrix, ...]:
    """Restrict a fixed-covector space to common levels on the given points."""

    if not covectors:
        return ()
    if len(points) < 2:
        return tuple(covectors)
    basis = sp.Matrix.hstack(*covectors)
    constraints = sp.Matrix.vstack(
        *((point - points[0]).T * basis for point in points[1:])
    )
    return tuple(basis * vector for vector in constraints.nullspace())


def surviving_collision_pairs(
    points: Sequence[sp.Matrix], projection: sp.Matrix | None = None
) -> tuple[tuple[int, int], ...]:
    """Return every pair that remains distinct after an optional projection."""

    projected = (
        list(points)
        if projection is None
        else [projection * point for point in points]
    )
    return tuple(
        (left, right)
        for left, right in combinations(range(len(projected)), 2)
        if projected[left] != projected[right]
    )


def companion_cancellation(
    variables: Sequence[sp.Symbol],
    companion_variables: Sequence[sp.Symbol],
    quadratic: Sequence[sp.Expr],
    matrix_b: sp.Matrix,
    cubic_basis: Sequence[sp.Expr],
) -> CompanionCancellation:
    """Construct and verify the specialization/cancellation identity.

    With ``F=x+Q+B*c``, the specialized companion map

        M=(x+Q+B*y, y-c)

    satisfies ``M=A_B o (F x I) o S_c`` where
    ``S_c(x,y)=(x,y-c)`` and ``A_B(u,v)=(u+B*v,v)``.
    """

    x = tuple(variables)
    y = tuple(companion_variables)
    q = tuple(sp.expand(value) for value in quadratic)
    c = tuple(sp.expand(value) for value in cubic_basis)
    n = len(x)
    rank = len(y)
    if len(q) != n or len(c) != rank:
        raise ValueError("incompatible quadratic or companion dimensions")
    if matrix_b.shape != (n, rank):
        raise ValueError("B has incompatible shape")

    b_c = matrix_b * sp.Matrix(c)
    b_y = matrix_b * sp.Matrix(y)
    base = tuple(sp.expand(x[i] + q[i] + b_c[i]) for i in range(n))
    source_second = tuple(sp.expand(y[j] - c[j]) for j in range(rank))
    specialized = tuple(
        [sp.expand(x[i] + q[i] + b_y[i]) for i in range(n)]
        + list(source_second)
    )
    source_shear = tuple(list(x) + list(source_second))

    # Apply F x I after S_c, then the target shear A_B.
    factored_first = [
        sp.expand(base[i] + sum(matrix_b[i, j] * source_second[j] for j in range(rank)))
        for i in range(n)
    ]
    factored = tuple(factored_first + list(source_second))
    if any(sp.expand(left - right) != 0 for left, right in zip(specialized, factored)):
        raise AssertionError("companion cancellation identity failed")
    return CompanionCancellation(
        variables=x,
        companion_variables=y,
        base_map=base,
        specialized_map=specialized,
        source_shear=source_shear,
        factored_map=factored,
    )


def parametric_companion_cancellation(
    variables: Sequence[sp.Symbol],
    companion_variables: Sequence[sp.Symbol],
    parameter: sp.Symbol,
    quadratic: Sequence[sp.Expr],
    matrix_b: sp.Matrix,
    cubic_basis: Sequence[sp.Expr],
) -> ParametricCompanionCancellation:
    """Cancel companions relatively and classify every parameter fiber.

    For homogeneous ``Q`` of degree two and ``c`` of degree three, set

        F=x+Q+B*c,
        E_t=x+t*Q+t^2*B*c,
        V=(x+t*Q+t^2*B*y, y-c, t).

    Then ``V=A_{t^2 B} o (E_t x I) o S_c`` over the full affine parameter
    line.  Moreover ``E_t(x)=t^-1 F(tx)`` over ``t != 0``, while the special
    fiber is the triangular automorphism ``(x,y-c(x))``.
    """

    x = tuple(variables)
    y = tuple(companion_variables)
    q = tuple(sp.expand(value) for value in quadratic)
    c = tuple(sp.expand(value) for value in cubic_basis)
    n = len(x)
    rank = len(y)
    if parameter in x or parameter in y:
        raise ValueError("parameter must be distinct from all map variables")
    if len(q) != n or len(c) != rank:
        raise ValueError("incompatible quadratic or companion dimensions")
    if matrix_b.shape != (n, rank):
        raise ValueError("B has incompatible shape")

    def is_homogeneous_of_degree(expression: sp.Expr, degree: int) -> bool:
        polynomial = sp.Poly(expression, *x, domain=sp.QQ)
        return all(
            sum(exponents) == degree
            for exponents, coefficient in polynomial.terms()
            if coefficient
        )

    if not all(is_homogeneous_of_degree(value, 2) for value in q):
        raise ValueError("Q must be homogeneous of degree two")
    if not all(is_homogeneous_of_degree(value, 3) for value in c):
        raise ValueError("c must be homogeneous of degree three")

    b_c = matrix_b * sp.Matrix(c)
    b_y = matrix_b * sp.Matrix(y)
    base = tuple(sp.expand(x[i] + q[i] + b_c[i]) for i in range(n))
    family = tuple(
        sp.expand(x[i] + parameter * q[i] + parameter**2 * b_c[i])
        for i in range(n)
    )
    parent_first = tuple(
        sp.expand(x[i] + parameter * q[i] + parameter**2 * b_y[i])
        for i in range(n)
    )
    source_second = tuple(sp.expand(y[j] - c[j]) for j in range(rank))
    parent = tuple(
        list(parent_first) + list(source_second) + [parameter]
    )
    source_shear = tuple(list(x) + list(source_second) + [parameter])
    factored_first = tuple(
        sp.expand(
            family[i]
            + parameter**2
            * sum(
                matrix_b[i, j] * source_second[j]
                for j in range(rank)
            )
        )
        for i in range(n)
    )
    factored = tuple(
        list(factored_first) + list(source_second) + [parameter]
    )
    if any(
        sp.expand(left - right) != 0
        for left, right in zip(parent, factored)
    ):
        raise AssertionError("parametric companion cancellation failed")

    scaling = {variable: parameter * variable for variable in x}
    scaled_base = tuple(
        sp.cancel(value.subs(scaling, simultaneous=True) / parameter)
        for value in base
    )
    if any(
        sp.expand(left - right) != 0
        for left, right in zip(family, scaled_base)
    ):
        raise AssertionError("nonzero-fiber scaling identity failed")

    special = tuple(
        sp.expand(value.subs(parameter, 0)) for value in parent[:-1]
    )
    expected_special = tuple(list(x) + list(source_second))
    if special != expected_special:
        raise AssertionError("special fiber is not the triangular shear")

    return ParametricCompanionCancellation(
        variables=x,
        companion_variables=y,
        parameter=parameter,
        base_map=base,
        scaled_base_family=family,
        parent_map=parent,
        source_shear=source_shear,
        factored_map=factored,
        special_fiber=special,
    )


def lift_point_to_nonzero_companion_slice(
    point: sp.Matrix,
    variables: Sequence[sp.Symbol],
    cubic_basis: Sequence[sp.Expr],
    parameter_value: sp.Expr,
) -> sp.Matrix:
    """Lift a base point to the corresponding nonzero parent slice."""

    value = sp.cancel(parameter_value)
    if value == 0:
        raise ValueError("the special fiber has no scaled base point")
    x = tuple(variables)
    if point.shape != (len(x), 1):
        raise ValueError("point dimension does not match the base variables")
    scaled = sp.Matrix([sp.cancel(entry / value) for entry in point])
    substitution = dict(zip(x, scaled))
    companions = [
        sp.expand(expression.subs(substitution, simultaneous=True))
        for expression in cubic_basis
    ]
    return scaled.col_join(sp.Matrix(companions)).col_join(
        sp.Matrix([value])
    )


def restrict_fixed_covector(
    expressions: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    covector: sp.Matrix,
    level: sp.Expr,
    *,
    prefix: str = "slice_x",
) -> FixedCovectorRestriction:
    """Restrict ``F`` to ``covector*x=level`` when that covector is fixed."""

    x = tuple(variables)
    f = tuple(sp.expand(value) for value in expressions)
    n = len(x)
    if len(f) != n or covector.shape != (n, 1):
        raise ValueError("map and covector dimensions do not match")
    correction = sp.Matrix(f) - sp.Matrix(x)
    if sp.expand((covector.T * correction)[0]) != 0:
        raise ValueError("covector is not fixed by the map")

    pivot = next(
        (index for index in range(n - 1, -1, -1) if covector[index] != 0),
        None,
    )
    if pivot is None:
        raise ValueError("zero covector does not define a hyperplane")
    retained = [index for index in range(n) if index != pivot]
    slice_variables = tuple(sp.symbols(f"{prefix}0:{n - 1}"))
    embedding: list[sp.Expr] = [sp.Integer(0)] * n
    for variable, index in zip(slice_variables, retained):
        embedding[index] = variable
    embedding[pivot] = sp.cancel(
        (
            level
            - sum(covector[index] * embedding[index] for index in retained)
        )
        / covector[pivot]
    )
    substitution = dict(zip(x, embedding))
    restricted = tuple(
        sp.expand(f[index].subs(substitution, simultaneous=True))
        for index in retained
    )
    if sp.expand(
        (covector.T * sp.Matrix(f))[0].subs(substitution, simultaneous=True)
        - level
    ) != 0:
        raise AssertionError("restricted map left its fixed level")
    return FixedCovectorRestriction(
        covector=covector,
        level=sp.expand(level),
        pivot=pivot,
        variables=slice_variables,
        embedding=tuple(embedding),
        restricted_map=restricted,
    )


def project_point_to_restriction(
    point: sp.Matrix, restriction: FixedCovectorRestriction
) -> sp.Matrix:
    """Drop the solved pivot coordinate from a point on the fixed level."""

    if sp.expand((restriction.covector.T * point)[0] - restriction.level) != 0:
        raise ValueError("point is not on the fixed level")
    return sp.Matrix(
        [point[index] for index in range(point.rows) if index != restriction.pivot]
    )


def profile_from_cubic_components(
    dimension: int,
    cubic_components: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> BackwardTerminalProfile:
    """Build both search objectives from an exact cubic output vector."""

    rank = coefficient_matrix(cubic_components, variables).rank()
    return BackwardTerminalProfile(dimension, rank)


def common_image(
    expressions: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    points: Iterable[sp.Matrix],
) -> sp.Matrix | None:
    """Return the common exact image, or ``None`` if the images differ."""

    images = [
        sp.Matrix(expressions).subs(dict(zip(variables, point)))
        for point in points
    ]
    if not images or any(image != images[0] for image in images[1:]):
        return None
    return images[0]
