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
