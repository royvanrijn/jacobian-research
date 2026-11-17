"""Projective-degree/Segre utilities for affine-map compactifications.

For an affine polynomial map ``F: A^n -> A^n`` of polynomial degree ``m``,
the compactification used here is

    [X0^m : F_1^h : ... : F_n^h] : P^n --> P^n.

This is deliberately not the full polar map of a homogenized potential.
If ``sigma_k`` are the signed degrees of the pushed-forward Segre class of
the base scheme, the projective degrees ``g_i`` are related by the
triangular binomial transform implemented below.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from typing import Iterable, Sequence

import sympy as sp


def _integer_tuple(values: Iterable[int], *, name: str) -> tuple[int, ...]:
    result = tuple(values)
    if any(not isinstance(value, int) for value in result):
        raise TypeError(f"{name} must contain Python integers")
    return result


def projective_degrees_from_segre(
    map_degree: int,
    segre_degrees: Sequence[int],
) -> tuple[int, ...]:
    """Return ``(g_0,...,g_n)`` from ``(sigma_1,...,sigma_n)``.

    The convention is

    ``g_i = m^i - sum_{k=1}^i binom(i,k)m^(i-k)sigma_k``.
    """

    if not isinstance(map_degree, int) or map_degree < 1:
        raise ValueError("map_degree must be a positive integer")
    sigmas = _integer_tuple(segre_degrees, name="segre_degrees")
    return (1,) + tuple(
        map_degree**index
        - sum(
            comb(index, codimension)
            * map_degree ** (index - codimension)
            * sigmas[codimension - 1]
            for codimension in range(1, index + 1)
        )
        for index in range(1, len(sigmas) + 1)
    )


def segre_degrees_from_projective(
    map_degree: int,
    projective_degrees: Sequence[int],
) -> tuple[int, ...]:
    """Return ``(sigma_1,...,sigma_n)`` from ``(g_0,...,g_n)``."""

    if not isinstance(map_degree, int) or map_degree < 1:
        raise ValueError("map_degree must be a positive integer")
    degrees = _integer_tuple(
        projective_degrees,
        name="projective_degrees",
    )
    if not degrees or degrees[0] != 1:
        raise ValueError("projective_degrees must start with g_0=1")

    sigmas: list[int] = []
    for index in range(1, len(degrees)):
        known = sum(
            comb(index, codimension)
            * map_degree ** (index - codimension)
            * sigmas[codimension - 1]
            for codimension in range(1, index)
        )
        sigmas.append(map_degree**index - known - degrees[index])
    return tuple(sigmas)


def total_segre_correction(
    map_degree: int,
    *,
    ambient_dimension: int,
    affine_degree: int,
) -> int:
    """Return the weighted top Segre correction ``m^n-g_n``."""

    if ambient_dimension < 1:
        raise ValueError("ambient_dimension must be positive")
    if affine_degree < 0:
        raise ValueError("affine_degree must be nonnegative")
    return map_degree**ambient_dimension - affine_degree


def equal_degree_complete_intersection_segre(
    *,
    ambient_dimension: int,
    codimension: int,
    generator_degree: int,
) -> tuple[int, ...]:
    """Return the pushed-forward Segre degrees of an equal-degree CI.

    For a complete intersection of codimension ``c`` cut out by degree
    ``d`` forms in projective ``n``-space, this expands

    ``(d H)^c / (1+d H)^c``

    through codimension ``n``.
    """

    if ambient_dimension < 1:
        raise ValueError("ambient_dimension must be positive")
    if not 1 <= codimension <= ambient_dimension:
        raise ValueError("codimension must lie between 1 and n")
    if generator_degree < 1:
        raise ValueError("generator_degree must be positive")
    return tuple(
        0
        if index < codimension
        else (
            (-1) ** (index - codimension)
            * comb(index - 1, codimension - 1)
            * generator_degree**index
        )
        for index in range(1, ambient_dimension + 1)
    )


def is_log_concave(projective_degrees: Sequence[int]) -> bool:
    """Test the elementary log-concavity inequalities for a degree list."""

    degrees = _integer_tuple(
        projective_degrees,
        name="projective_degrees",
    )
    return all(
        degrees[index] ** 2
        >= degrees[index - 1] * degrees[index + 1]
        for index in range(1, len(degrees) - 1)
    )


def homogenize_to_degree(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    homogenizing_variable: sp.Symbol,
    degree: int,
) -> sp.Expr:
    """Homogenize ``polynomial`` to the explicitly supplied total degree."""

    if degree < 0:
        raise ValueError("degree must be nonnegative")
    poly = sp.Poly(sp.expand(polynomial), *variables)
    if poly.total_degree() > degree:
        raise ValueError("polynomial has degree larger than requested degree")
    result = sp.S.Zero
    for exponents, coefficient in poly.terms():
        term_degree = sum(exponents)
        monomial = sp.prod(
            variable**exponent
            for variable, exponent in zip(variables, exponents)
        )
        result += (
            coefficient
            * homogenizing_variable ** (degree - term_degree)
            * monomial
        )
    return sp.expand(result)


def affine_map_compactification(
    coordinates: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    homogenizing_variable: sp.Symbol,
    *,
    map_degree: int | None = None,
) -> tuple[sp.Expr, ...]:
    """Construct ``[X0^m:F_1^h:...:F_n^h]`` for an affine map."""

    if len(coordinates) != len(variables):
        raise ValueError("an affine self-map needs one coordinate per variable")
    actual_degree = max(
        sp.Poly(sp.expand(coordinate), *variables).total_degree()
        for coordinate in coordinates
    )
    degree = actual_degree if map_degree is None else map_degree
    if degree < 1:
        raise ValueError("the affine map must have positive polynomial degree")
    if degree < actual_degree:
        raise ValueError("map_degree is below the actual polynomial degree")
    return (homogenizing_variable**degree,) + tuple(
        homogenize_to_degree(
            coordinate,
            variables,
            homogenizing_variable,
            degree,
        )
        for coordinate in coordinates
    )


def affine_gradient_compactification(
    potential: sp.Expr,
    variables: Sequence[sp.Symbol],
    homogenizing_variable: sp.Symbol,
) -> tuple[sp.Expr, ...]:
    """Construct the actual affine-gradient compactification ``Gamma_Psi``."""

    gradient = tuple(sp.diff(potential, variable) for variable in variables)
    return affine_map_compactification(
        gradient,
        variables,
        homogenizing_variable,
    )


def full_polar_map(
    potential: sp.Expr,
    variables: Sequence[sp.Symbol],
    homogenizing_variable: sp.Symbol,
) -> tuple[sp.Expr, ...]:
    """Construct the full polar map of the homogenized potential.

    This function is intentionally separate from
    :func:`affine_gradient_compactification`.
    """

    potential_degree = sp.Poly(sp.expand(potential), *variables).total_degree()
    homogenized = homogenize_to_degree(
        potential,
        variables,
        homogenizing_variable,
        potential_degree,
    )
    return tuple(
        sp.expand(sp.diff(homogenized, variable))
        for variable in (homogenizing_variable, *variables)
    )


def homogeneous_leading_forms(
    coordinates: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    map_degree: int | None = None,
) -> tuple[sp.Expr, ...]:
    """Return the common top-degree forms cutting out the infinity support."""

    actual_degree = max(
        sp.Poly(sp.expand(coordinate), *variables).total_degree()
        for coordinate in coordinates
    )
    degree = actual_degree if map_degree is None else map_degree
    if degree < actual_degree:
        raise ValueError("map_degree is below the actual polynomial degree")
    result = []
    for coordinate in coordinates:
        poly = sp.Poly(sp.expand(coordinate), *variables)
        result.append(
            sp.expand(
                sum(
                    coefficient
                    * sp.prod(
                        variable**exponent
                        for variable, exponent in zip(variables, exponents)
                    )
                    for exponents, coefficient in poly.terms()
                    if sum(exponents) == degree
                )
            )
        )
    return tuple(result)


def integrability_residuals(
    coordinates: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return all curls ``dF_i/dx_j-dF_j/dx_i``."""

    if len(coordinates) != len(variables):
        raise ValueError("integrability needs one coordinate per variable")
    return tuple(
        sp.expand(
            sp.diff(coordinates[left], variables[right])
            - sp.diff(coordinates[right], variables[left])
        )
        for left in range(len(variables))
        for right in range(left + 1, len(variables))
    )


def integrate_homogeneous_gradient(
    leading_forms: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> sp.Expr:
    """Recover a homogeneous potential from an integrable homogeneous tuple.

    Raises ``ValueError`` if the forms do not have one common positive degree
    or fail the curl equations.  Euler's identity gives the unique potential
    with zero constant term.
    """

    if len(leading_forms) != len(variables):
        raise ValueError("one leading form is required per variable")
    degrees = {
        sp.Poly(sp.expand(form), *variables).total_degree()
        for form in leading_forms
        if form != 0
    }
    if len(degrees) != 1:
        raise ValueError("nonzero leading forms must have one common degree")
    if not degrees:
        raise ValueError("the zero tuple has no positive leading degree")
    degree = degrees.pop()
    if degree < 1:
        raise ValueError("leading forms must have positive degree")
    if any(integrability_residuals(leading_forms, variables)):
        raise ValueError("leading forms are not an integrable gradient tuple")
    potential = sp.expand(
        sum(
            variable * form
            for variable, form in zip(variables, leading_forms)
        )
        / (degree + 1)
    )
    if any(
        sp.expand(sp.diff(potential, variable) - form) != 0
        for variable, form in zip(variables, leading_forms)
    ):
        raise ValueError("Euler reconstruction failed")
    return potential


def equal_degree_ci_hilbert_function(
    *,
    generator_degree: int,
    codimension: int,
) -> tuple[int, ...]:
    """Return coefficients of ``(1+z+...+z^(d-1))^c``.

    This is the Hilbert function of an Artinian complete intersection of
    ``codimension`` forms, all of degree ``generator_degree``.
    """

    if generator_degree < 1:
        raise ValueError("generator_degree must be positive")
    if codimension < 1:
        raise ValueError("codimension must be positive")
    base = (1,) * generator_degree
    result = (1,)
    for _ in range(codimension):
        product = [0] * (len(result) + len(base) - 1)
        for left_index, left_value in enumerate(result):
            for right_index, right_value in enumerate(base):
                product[left_index + right_index] += (
                    left_value * right_value
                )
        result = tuple(product)
    return result


def filtered_missing_generator_drop(
    *,
    map_degree: int,
    epsilon_order: int,
    cyclic_ideal_dimension: int,
) -> int:
    """Lower-bound the length removed by one missing gradient component.

    If its initial form in the smooth-essential normal slice is
    ``epsilon^q*s``, this returns

    ``(map_degree-q) * dim(B*s)``.
    """

    if map_degree < 2:
        raise ValueError("map_degree must be at least two")
    if not 1 <= epsilon_order < map_degree:
        raise ValueError("epsilon_order must lie between 1 and m-1")
    if cyclic_ideal_dimension < 1:
        raise ValueError("cyclic_ideal_dimension must be positive")
    return (
        map_degree - epsilon_order
    ) * cyclic_ideal_dimension


def truncated_dvr_module_length(
    *,
    truncation_order: int,
    generic_rank: int,
    torsion_orders: Sequence[int],
) -> int:
    """Return the length after truncating a finite DVR module.

    For

    ``M = R^rho + direct_sum_j R/(epsilon^a_j)``,

    with ``R=K[[epsilon]]``, the quotient by ``epsilon^m`` has length

    ``m*rho + sum_j min(m,a_j)``.
    """

    if truncation_order < 1:
        raise ValueError("truncation_order must be positive")
    if generic_rank < 0:
        raise ValueError("generic_rank must be nonnegative")
    orders = _integer_tuple(torsion_orders, name="torsion_orders")
    if any(order < 1 for order in orders):
        raise ValueError("torsion orders must be positive")
    return (
        truncation_order * generic_rank
        + sum(min(truncation_order, order) for order in orders)
    )


@dataclass(frozen=True)
class SmoothEssentialGradientNormalSlice:
    """Dimension-free normal slice for a smooth essential gradient top.

    The top potential depends on ``essential_rank`` variables and defines a
    smooth hypersurface in their projective space.  At the generic point of
    the kernel vertex, the active Jacobian algebra is an equal-degree
    complete intersection and the compactifying parameter is truncated at
    ``map_degree``.
    """

    ambient_dimension: int
    map_degree: int
    essential_rank: int

    def __post_init__(self) -> None:
        if self.ambient_dimension < 2:
            raise ValueError("ambient_dimension must be at least two")
        if self.map_degree < 2:
            raise ValueError("map_degree must be at least two")
        if not 1 <= self.essential_rank < self.ambient_dimension:
            raise ValueError(
                "essential_rank must lie between 1 and ambient_dimension-1"
            )

    @property
    def kernel_dimension(self) -> int:
        return self.ambient_dimension - self.essential_rank

    @property
    def support_dimension(self) -> int:
        return self.kernel_dimension - 1

    @property
    def base_codimension(self) -> int:
        return self.essential_rank + 1

    @property
    def jacobian_hilbert_function(self) -> tuple[int, ...]:
        return equal_degree_ci_hilbert_function(
            generator_degree=self.map_degree,
            codimension=self.essential_rank,
        )

    @property
    def jacobian_length(self) -> int:
        return self.map_degree**self.essential_rank

    @property
    def jacobian_socle_degree(self) -> int:
        return self.essential_rank * (self.map_degree - 1)

    @property
    def truncated_active_length(self) -> int:
        return self.map_degree ** (self.essential_rank + 1)

    @property
    def unit_penultimate_segre_degree(self) -> int:
        """Leading Segre degree when one epsilon-order-one term is a unit."""

        return self.jacobian_length

    def missing_generator_drop(
        self,
        *,
        epsilon_order: int,
        cyclic_ideal_dimension: int,
    ) -> int:
        if cyclic_ideal_dimension > self.jacobian_length:
            raise ValueError("cyclic ideal cannot be larger than B")
        return filtered_missing_generator_drop(
            map_degree=self.map_degree,
            epsilon_order=epsilon_order,
            cyclic_ideal_dimension=cyclic_ideal_dimension,
        )

    def leading_segre_upper_bound(
        self,
        *,
        epsilon_order: int,
        cyclic_ideal_dimension: int,
    ) -> int:
        """Upper-bound ``sigma_(r+1)`` from one missing component."""

        return self.truncated_active_length - self.missing_generator_drop(
            epsilon_order=epsilon_order,
            cyclic_ideal_dimension=cyclic_ideal_dimension,
        )

    def isolated_vertex_affine_degree_lower_bound(
        self,
        *,
        epsilon_order: int,
        cyclic_ideal_dimension: int,
    ) -> int:
        """Lower-bound the affine degree when the kernel vertex is a point."""

        if self.kernel_dimension != 1:
            raise ValueError("the kernel vertex is not zero-dimensional")
        return self.missing_generator_drop(
            epsilon_order=epsilon_order,
            cyclic_ideal_dimension=cyclic_ideal_dimension,
        )


@dataclass(frozen=True)
class SingularEssentialGradientNormalSlice:
    """Generic normal slice along one singular component of the top.

    The top potential uses ``essential_rank`` variables.  Its selected
    projective singular component has dimension
    ``singular_locus_dimension``, degree ``component_degree``, and generic
    transverse Jacobian length ``transverse_jacobian_length``.

    Lower potential layers are encoded by the finite
    ``K[[epsilon]]``-module profile of the active gradient quotient:
    ``generic_rank`` and ``torsion_orders``.
    """

    ambient_dimension: int
    map_degree: int
    essential_rank: int
    singular_locus_dimension: int
    transverse_jacobian_length: int
    component_degree: int = 1

    def __post_init__(self) -> None:
        if self.ambient_dimension < 3:
            raise ValueError("ambient_dimension must be at least three")
        if self.map_degree < 2:
            raise ValueError("map_degree must be at least two")
        if not 2 <= self.essential_rank < self.ambient_dimension:
            raise ValueError(
                "essential_rank must lie between 2 and ambient_dimension-1"
            )
        if not (
            0
            <= self.singular_locus_dimension
            <= self.essential_rank - 2
        ):
            raise ValueError(
                "singular-locus dimension must lie between 0 and r-2"
            )
        if self.transverse_jacobian_length < 1:
            raise ValueError("transverse Jacobian length must be positive")
        if self.component_degree < 1:
            raise ValueError("component degree must be positive")

    @property
    def kernel_dimension(self) -> int:
        return self.ambient_dimension - self.essential_rank

    @property
    def joined_support_dimension(self) -> int:
        return (
            self.kernel_dimension
            + self.singular_locus_dimension
        )

    @property
    def base_codimension(self) -> int:
        return (
            self.essential_rank
            - self.singular_locus_dimension
        )

    @property
    def flat_active_truncated_length(self) -> int:
        return (
            self.map_degree
            * self.transverse_jacobian_length
        )

    def active_truncated_length(
        self,
        *,
        generic_rank: int,
        torsion_orders: Sequence[int],
    ) -> int:
        """Return the active length for a validated DVR profile."""

        orders = _integer_tuple(torsion_orders, name="torsion_orders")
        if (
            generic_rank + len(orders)
            != self.transverse_jacobian_length
        ):
            raise ValueError(
                "generic rank plus torsion summands must equal "
                "the transverse Jacobian length"
            )
        return truncated_dvr_module_length(
            truncation_order=self.map_degree,
            generic_rank=generic_rank,
            torsion_orders=orders,
        )

    def leading_segre_contribution_bounds(
        self,
        *,
        generic_rank: int,
        torsion_orders: Sequence[int],
    ) -> tuple[int, int]:
        """Bound this component's leading Segre contribution.

        Missing kernel-gradient components are epsilon-divisible.  Thus the
        final transverse quotient still has special fiber of length ``mu``,
        while it is a quotient of the active truncated algebra.
        """

        active_length = self.active_truncated_length(
            generic_rank=generic_rank,
            torsion_orders=torsion_orders,
        )
        return (
            self.component_degree * self.transverse_jacobian_length,
            self.component_degree * active_length,
        )

    @property
    def order_one_active_segre_contribution(self) -> int:
        """Exact contribution when the active DVR profile has order one."""

        return (
            self.component_degree
            * self.transverse_jacobian_length
        )

    @property
    def unit_kernel_gradient_segre_contribution(self) -> int:
        """Exact contribution when a missing component is epsilon times a unit."""

        return self.order_one_active_segre_contribution

    def flat_missing_generator_upper_bound(
        self,
        *,
        epsilon_order: int,
        cyclic_ideal_dimension: int,
    ) -> int:
        """Upper-bound the contribution when the active DVR module is flat."""

        if cyclic_ideal_dimension > self.transverse_jacobian_length:
            raise ValueError("cyclic ideal cannot be larger than B")
        drop = filtered_missing_generator_drop(
            map_degree=self.map_degree,
            epsilon_order=epsilon_order,
            cyclic_ideal_dimension=cyclic_ideal_dimension,
        )
        return self.component_degree * (
            self.flat_active_truncated_length - drop
        )


@dataclass(frozen=True)
class ProjectiveGradientSegreRecord:
    """Validated invariant record for one actual affine compactification."""

    ambient_dimension: int
    map_degree: int
    projective_degrees: tuple[int, ...]
    segre_degrees: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.ambient_dimension < 1:
            raise ValueError("ambient_dimension must be positive")
        if len(self.projective_degrees) != self.ambient_dimension + 1:
            raise ValueError("projective degree list has the wrong length")
        if len(self.segre_degrees) != self.ambient_dimension:
            raise ValueError("Segre degree list has the wrong length")
        if projective_degrees_from_segre(
            self.map_degree,
            self.segre_degrees,
        ) != self.projective_degrees:
            raise ValueError("projective and Segre degrees are inconsistent")

    @classmethod
    def from_projective_degrees(
        cls,
        *,
        map_degree: int,
        projective_degrees: Sequence[int],
    ) -> "ProjectiveGradientSegreRecord":
        degrees = _integer_tuple(
            projective_degrees,
            name="projective_degrees",
        )
        if len(degrees) < 2:
            raise ValueError("a positive ambient dimension is required")
        sigmas = segre_degrees_from_projective(map_degree, degrees)
        return cls(len(degrees) - 1, map_degree, degrees, sigmas)

    @classmethod
    def from_segre_degrees(
        cls,
        *,
        map_degree: int,
        segre_degrees: Sequence[int],
    ) -> "ProjectiveGradientSegreRecord":
        sigmas = _integer_tuple(segre_degrees, name="segre_degrees")
        if not sigmas:
            raise ValueError("a positive ambient dimension is required")
        degrees = projective_degrees_from_segre(map_degree, sigmas)
        return cls(len(sigmas), map_degree, degrees, sigmas)

    @property
    def affine_degree(self) -> int:
        """Top projective degree of the actual affine-map compactification."""

        return self.projective_degrees[-1]

    @property
    def weighted_correction(self) -> int:
        """Weighted top Segre correction."""

        return self.map_degree**self.ambient_dimension - self.affine_degree

    @property
    def has_fixed_divisor(self) -> bool:
        """Whether the degree list detects a nonzero codimension-one class."""

        return bool(self.segre_degrees[0])
