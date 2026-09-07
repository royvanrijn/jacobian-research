#!/usr/bin/env python3
"""Exact numerical gates for intrinsic SNC boundaries of the affine plane.

The boundary-lattice prefilter checks the localization sequence from boundary
classes.  This module uses the full intersection matrix.  For a smooth SNC
completion ``X \\ A2 = D_1 + ... + D_r``, adjunction reconstructs the
canonical class in the boundary basis, and rationality forces

    K_X^2 + rho(X) = 10.

For a resolved polynomial Keller map ``f: X -> P2``, the pole vector of the
target line at infinity then determines the ordinary and logarithmic
ramification vectors.  All computations are exact over the integers and
rationals.

These are necessary conditions.  Passing them does not construct a surface
or a Keller map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import sympy as sp

from finite_normalization_signatures import (
    ResidualDifferentAudit,
    ResidualDifferentComponent,
    audit_residual_different,
)


def _integral_matrix(matrix: Sequence[Sequence[int]] | sp.Matrix) -> sp.Matrix:
    result = sp.Matrix(matrix)
    if result.rows == 0 or result.cols == 0:
        raise ValueError("an intrinsic boundary matrix must be nonempty")
    if result.rows != result.cols:
        raise ValueError("an intrinsic boundary intersection matrix must be square")
    if result != result.T:
        raise ValueError("an intrinsic boundary intersection matrix must be symmetric")
    if any(value.is_Integer is not True for value in result):
        raise ValueError("boundary intersection numbers must be integers")
    return result


def symmetric_inertia(
    matrix: Sequence[Sequence[int]] | sp.Matrix,
) -> tuple[int, int, int]:
    """Return the exact positive, negative, and zero inertia.

    Symmetric Gaussian elimination is performed by congruences.  A nonzero
    diagonal entry contributes a one-dimensional pivot.  If every diagonal
    entry is zero, a nonzero off-diagonal entry contributes the invertible
    block ``[[0,b],[b,0]]``, whose inertia is ``(1,1,0)``.
    """

    work = sp.Matrix(matrix).applyfunc(sp.Rational)
    if work.rows != work.cols or work != work.T:
        raise ValueError("inertia requires a symmetric square matrix")

    positive = negative = zero = 0
    while work.rows:
        size = work.rows
        diagonal = next(
            (index for index in range(size) if work[index, index] != 0),
            None,
        )
        if diagonal is not None:
            order = [diagonal] + [i for i in range(size) if i != diagonal]
            work = work.extract(order, order)
            pivot = work[0, 0]
            if pivot > 0:
                positive += 1
            else:
                negative += 1
            if size == 1:
                break
            column = work[1:, 0]
            work = (
                work[1:, 1:] - column * column.T / pivot
            ).applyfunc(sp.cancel)
            continue

        off_diagonal = next(
            (
                (i, j)
                for i in range(size)
                for j in range(i + 1, size)
                if work[i, j] != 0
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break

        first, second = off_diagonal
        order = [first, second] + [
            i for i in range(size) if i not in (first, second)
        ]
        work = work.extract(order, order)
        positive += 1
        negative += 1
        if size == 2:
            break
        pivot = work[:2, :2]
        cross = work[:2, 2:]
        work = (
            work[2:, 2:] - cross.T * pivot.inv() * cross
        ).applyfunc(sp.cancel)

    return positive, negative, zero


def _tree_edges(matrix: sp.Matrix) -> tuple[tuple[int, int], ...]:
    return tuple(
        (i, j)
        for i in range(matrix.rows)
        for j in range(i + 1, matrix.cols)
        if matrix[i, j] != 0
    )


def _is_connected(vertex_count: int, edges: Sequence[tuple[int, int]]) -> bool:
    adjacency = [set() for _ in range(vertex_count)]
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    seen = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency[current] - seen:
            seen.add(neighbor)
            frontier.append(neighbor)
    return len(seen) == vertex_count


@dataclass(frozen=True)
class IntrinsicA2Boundary:
    """The numerical skeleton of a proposed smooth SNC boundary."""

    names: tuple[str, ...]
    intersection_matrix: sp.Matrix
    genera: tuple[int, ...]

    def __post_init__(self) -> None:
        matrix = _integral_matrix(self.intersection_matrix)
        object.__setattr__(self, "intersection_matrix", matrix)
        if len(self.names) != matrix.rows:
            raise ValueError("boundary names must index the intersection matrix")
        if len(set(self.names)) != len(self.names):
            raise ValueError("boundary names must be distinct")
        if len(self.genera) != matrix.rows:
            raise ValueError("one geometric genus is required per boundary component")
        if any(not isinstance(genus, int) or genus < 0 for genus in self.genera):
            raise ValueError("boundary genera must be nonnegative integers")

    @property
    def rank(self) -> int:
        return self.intersection_matrix.rows

    @property
    def edges(self) -> tuple[tuple[int, int], ...]:
        return _tree_edges(self.intersection_matrix)

    @property
    def valencies(self) -> tuple[int, ...]:
        result = [0] * self.rank
        for first, second in self.edges:
            result[first] += 1
            result[second] += 1
        return tuple(result)

    @property
    def adjunction_vector(self) -> sp.Matrix:
        """The intersections ``K_X . D_i`` forced by adjunction."""

        matrix = self.intersection_matrix
        return sp.Matrix(
            [
                2 * self.genera[index] - 2 - matrix[index, index]
                for index in range(self.rank)
            ]
        )

    @property
    def canonical_coefficients(self) -> sp.Matrix:
        """Coefficients of ``K_X`` in the boundary-prime basis."""

        return self.intersection_matrix.inv() * self.adjunction_vector

    @property
    def log_canonical_coefficients(self) -> sp.Matrix:
        """Coefficients of ``K_X + D`` in the boundary-prime basis."""

        return self.canonical_coefficients + sp.ones(self.rank, 1)

    @property
    def canonical_square(self) -> sp.Rational:
        coefficients = self.canonical_coefficients
        return sp.cancel(
            (coefficients.T * self.intersection_matrix * coefficients)[0]
        )


@dataclass(frozen=True)
class A2BoundaryAudit:
    boundary: IntrinsicA2Boundary
    determinant: int
    inertia: tuple[int, int, int]
    canonical_coefficients: tuple[sp.Rational, ...]
    canonical_square: sp.Rational
    expected_canonical_square: int
    failures: tuple[str, ...]

    @property
    def passes(self) -> bool:
        return not self.failures

    @property
    def nonnegative_canonical_components(self) -> tuple[str, ...]:
        return tuple(
            name
            for name, coefficient in zip(
                self.boundary.names, self.canonical_coefficients
            )
            if coefficient >= 0
        )

    @property
    def can_support_nonproper_keller_boundary(self) -> bool:
        """A nonproper Keller resolution needs a dicritical with ``k_i >= 0``."""

        return self.passes and bool(self.nonnegative_canonical_components)

    def as_dict(self) -> dict[str, object]:
        return {
            "names": list(self.boundary.names),
            "intersection_matrix": [
                [int(value) for value in row]
                for row in self.boundary.intersection_matrix.tolist()
            ],
            "genera": list(self.boundary.genera),
            "determinant": self.determinant,
            "inertia": list(self.inertia),
            "canonical_coefficients": [
                str(value) for value in self.canonical_coefficients
            ],
            "canonical_square": str(self.canonical_square),
            "expected_canonical_square": self.expected_canonical_square,
            "nonnegative_canonical_components": list(
                self.nonnegative_canonical_components
            ),
            "passes": self.passes,
            "failures": list(self.failures),
        }


def audit_a2_boundary(boundary: IntrinsicA2Boundary) -> A2BoundaryAudit:
    """Apply exact necessary conditions for ``X \\ D`` to be ``A2``."""

    matrix = boundary.intersection_matrix
    rank = boundary.rank
    failures: list[str] = []

    if any(genus != 0 for genus in boundary.genera):
        failures.append("the boundary is not a union of rational curves")

    invalid_crossings = tuple(
        (boundary.names[i], boundary.names[j], int(matrix[i, j]))
        for i in range(rank)
        for j in range(i + 1, rank)
        if matrix[i, j] not in (0, 1)
    )
    if invalid_crossings:
        failures.append(
            "the boundary is not a simple rational tree: "
            f"invalid crossings {invalid_crossings}"
        )

    edges = boundary.edges
    if len(edges) != rank - 1 or not _is_connected(rank, edges):
        failures.append("the boundary dual graph is not a tree")

    determinant = int(matrix.det())
    expected_determinant = (-1) ** (rank - 1)
    if determinant != expected_determinant:
        failures.append(
            f"det(Q)={determinant}, expected {expected_determinant}"
        )

    inertia = symmetric_inertia(matrix)
    if inertia != (1, rank - 1, 0):
        failures.append(
            f"inertia(Q)={inertia}, expected {(1, rank - 1, 0)}"
        )

    canonical: tuple[sp.Rational, ...]
    canonical_square: sp.Rational
    if determinant == 0:
        canonical = ()
        canonical_square = sp.nan
        failures.append("adjunction cannot reconstruct K_X from a singular Q")
    else:
        canonical_matrix = boundary.canonical_coefficients
        canonical = tuple(sp.cancel(value) for value in canonical_matrix)
        if any(value.q != 1 for value in canonical):
            failures.append(
                f"the adjunction solution for K_X is nonintegral: {canonical}"
            )
        canonical_square = boundary.canonical_square
        expected_square = 10 - rank
        if canonical_square != expected_square:
            failures.append(
                f"K_X^2={canonical_square}, expected 10-r={expected_square}"
            )

    return A2BoundaryAudit(
        boundary=boundary,
        determinant=determinant,
        inertia=inertia,
        canonical_coefficients=canonical,
        canonical_square=canonical_square,
        expected_canonical_square=10 - rank,
        failures=tuple(failures),
    )


@dataclass(frozen=True)
class KellerPoleAudit:
    """Numerical canonical/ramification audit for one target pole vector."""

    boundary_audit: A2BoundaryAudit
    pole_vector: tuple[int, ...]
    hyperplane_intersections: tuple[int, ...]
    geometric_degree: int
    ramification_coefficients: tuple[int, ...]
    log_ramification_coefficients: tuple[int, ...]
    dicritical_candidates: tuple[str, ...]
    failures: tuple[str, ...]

    @property
    def passes(self) -> bool:
        return not self.failures

    def as_dict(self) -> dict[str, object]:
        return {
            "names": list(self.boundary_audit.boundary.names),
            "pole_vector": list(self.pole_vector),
            "hyperplane_intersections": list(self.hyperplane_intersections),
            "geometric_degree": self.geometric_degree,
            "ramification_coefficients": list(self.ramification_coefficients),
            "log_ramification_coefficients": list(
                self.log_ramification_coefficients
            ),
            "dicritical_candidates": list(self.dicritical_candidates),
            "passes": self.passes,
            "failures": list(self.failures),
        }


@dataclass(frozen=True)
class KellerDicriticalDatum:
    """Clean finite-chart input for one intrinsic dicritical component.

    The pole audit determines the transverse ramification index and the
    product ``residue_degree * target_curve_degree``.  Thus only the degree
    of the image curve and its companion-sheet intersection must be supplied.
    """

    name: str
    target_curve_degree: int
    companion_intersection: int

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a Keller dicritical datum needs a component name")
        if self.target_curve_degree <= 0:
            raise ValueError("the target curve degree must be positive")
        if self.companion_intersection < 0:
            raise ValueError("the companion intersection must be nonnegative")


@dataclass(frozen=True)
class KellerDicriticalBudget:
    """Intrinsic reconstruction of one clean dicritical's finite-cover data."""

    name: str
    boundary_index: int
    ramification_index: int
    residue_degree: int
    available_residual_intersection: int
    forced_companion_intersection: int

    @property
    def feasible(self) -> bool:
        return self.forced_companion_intersection >= 0


@dataclass(frozen=True)
class BoundaryContractionAudit:
    """Exact intersection form after contracting a negative-definite block."""

    original_names: tuple[str, ...]
    contracted_names: tuple[str, ...]
    surviving_names: tuple[str, ...]
    contracted_inertia: tuple[int, int, int]
    pullback_correction: sp.Matrix
    finite_intersection_matrix: sp.Matrix
    failures: tuple[str, ...]

    @property
    def passes(self) -> bool:
        return not self.failures


@dataclass(frozen=True)
class FiniteModelDicriticalProjectionBudget:
    """Projection/different comparison on the finite Stein model."""

    name: str
    target_curve_degree: int
    ramification_index: int
    residue_degree: int
    finite_self_intersection: sp.Rational
    surface_different_degree: sp.Rational
    residual_forced_companion_intersection: sp.Rational
    projection_forced_companion_intersection: sp.Rational
    target_normalization_correction: int
    corrected_budget_gap: sp.Rational
    contraction: BoundaryContractionAudit

    @property
    def budgets_match(self) -> bool:
        return self.corrected_budget_gap == 0


def audit_keller_pole_vector(
    boundary_audit: A2BoundaryAudit,
    pole_vector: Sequence[int],
    *,
    require_nonproper: bool = True,
) -> KellerPoleAudit:
    """Audit ``f^*L`` for a resolved Keller map ``f: X -> P2``.

    If ``p`` is the coefficient vector of ``f^*L`` and ``k`` is that of
    ``K_X``, then

    ``Qp`` records intersections with a general target line,
    ``p^T Q p`` is the geometric degree,
    ``k+3p`` is the ordinary ramification vector, and
    ``k+1+2p`` is the logarithmic ramification vector.

    A component is dicritical over the affine target exactly when
    ``p_i=0`` and ``(Qp)_i>0``.
    """

    boundary = boundary_audit.boundary
    failures = list(boundary_audit.failures)
    if len(pole_vector) != boundary.rank:
        raise ValueError("the pole vector must have one entry per boundary prime")
    if any(not isinstance(value, int) for value in pole_vector):
        raise ValueError("pole orders must be integers")

    p = sp.Matrix(pole_vector)
    if any(value < 0 for value in p):
        failures.append("the target-line pullback is not effective")
    intersections = boundary.intersection_matrix * p
    if any(value < 0 for value in intersections):
        failures.append(
            f"f^*L is not nef on the boundary: Qp={tuple(intersections)}"
        )

    degree_value = sp.cancel(
        (p.T * boundary.intersection_matrix * p)[0]
    )
    if degree_value.q != 1:
        failures.append(f"the geometric degree is nonintegral: {degree_value}")
    degree = int(degree_value)
    if degree <= 0:
        failures.append(f"the geometric degree is not positive: {degree}")

    canonical = sp.Matrix(boundary_audit.canonical_coefficients)
    if canonical.rows != boundary.rank:
        ramification = sp.zeros(boundary.rank, 1)
        log_ramification = sp.zeros(boundary.rank, 1)
    else:
        ramification = canonical + 3 * p
        log_ramification = canonical + sp.ones(boundary.rank, 1) + 2 * p
        if any(value < 0 for value in ramification):
            failures.append(
                "the ordinary ramification vector k+3p is not effective: "
                f"{tuple(ramification)}"
            )
        if any(value < 0 for value in log_ramification):
            failures.append(
                "the logarithmic ramification vector k+1+2p is not effective: "
                f"{tuple(log_ramification)}"
            )

    dicritical_indices = tuple(
        index
        for index in range(boundary.rank)
        if p[index] == 0 and intersections[index] > 0
    )
    dicritical = tuple(boundary.names[index] for index in dicritical_indices)
    if require_nonproper and not dicritical:
        failures.append(
            "no boundary prime is dicritical over the affine target "
            "(p_i=0 and (Qp)_i>0)"
        )
    for index in dicritical_indices:
        if canonical.rows == boundary.rank and canonical[index] < 0:
            failures.append(
                f"dicritical {boundary.names[index]} has negative "
                f"canonical coefficient {canonical[index]}"
            )

    return KellerPoleAudit(
        boundary_audit=boundary_audit,
        pole_vector=tuple(int(value) for value in p),
        hyperplane_intersections=tuple(int(value) for value in intersections),
        geometric_degree=degree,
        ramification_coefficients=tuple(int(value) for value in ramification),
        log_ramification_coefficients=tuple(
            int(value) for value in log_ramification
        ),
        dicritical_candidates=dicritical,
        failures=tuple(failures),
    )


def contract_boundary_curves(
    boundary: IntrinsicA2Boundary,
    contracted_names: Sequence[str],
) -> BoundaryContractionAudit:
    """Contract a boundary block and compute the Mumford intersection form.

    If ``Q_VV`` is the intersection block of the contracted curves, the
    pullback of the surviving primes is corrected by
    ``-Q_VV^{-1}Q_VH`` and their intersection matrix is

    ``Q_HH-Q_HV Q_VV^{-1}Q_VH``.
    """

    requested = tuple(contracted_names)
    if len(set(requested)) != len(requested):
        raise ValueError("contracted boundary names must be distinct")
    unknown = tuple(name for name in requested if name not in boundary.names)
    if unknown:
        raise ValueError(f"unknown contracted boundary names: {unknown}")
    contracted_set = set(requested)
    vertical = tuple(
        index
        for index, name in enumerate(boundary.names)
        if name in contracted_set
    )
    horizontal = tuple(
        index
        for index, name in enumerate(boundary.names)
        if name not in contracted_set
    )
    if not horizontal:
        raise ValueError("at least one boundary curve must survive")
    surviving_names = tuple(boundary.names[index] for index in horizontal)
    contracted_ordered = tuple(boundary.names[index] for index in vertical)
    q_hh = boundary.intersection_matrix.extract(horizontal, horizontal)
    failures: list[str] = []
    if not vertical:
        return BoundaryContractionAudit(
            original_names=boundary.names,
            contracted_names=(),
            surviving_names=surviving_names,
            contracted_inertia=(0, 0, 0),
            pullback_correction=sp.zeros(0, len(horizontal)),
            finite_intersection_matrix=q_hh,
            failures=(),
        )

    q_vv = boundary.intersection_matrix.extract(vertical, vertical)
    inertia = symmetric_inertia(q_vv)
    if inertia != (0, len(vertical), 0):
        failures.append(
            "the contracted intersection block is not negative definite: "
            f"inertia={inertia}"
        )
    if q_vv.det() == 0:
        failures.append("the contracted intersection block is singular")
        correction = sp.zeros(len(vertical), len(horizontal))
        finite = q_hh
    else:
        q_vh = boundary.intersection_matrix.extract(vertical, horizontal)
        q_hv = q_vh.T
        correction = (-q_vv.inv() * q_vh).applyfunc(sp.cancel)
        finite = (q_hh - q_hv * q_vv.inv() * q_vh).applyfunc(sp.cancel)
    return BoundaryContractionAudit(
        original_names=boundary.names,
        contracted_names=contracted_ordered,
        surviving_names=surviving_names,
        contracted_inertia=inertia,
        pullback_correction=correction,
        finite_intersection_matrix=finite,
        failures=tuple(failures),
    )


def contract_keller_vertical_boundary(
    pole_audit: KellerPoleAudit,
) -> BoundaryContractionAudit:
    """Contract exactly the boundary curves killed by ``f^*L``."""

    boundary = pole_audit.boundary_audit.boundary
    vertical = tuple(
        name
        for name, intersection in zip(
            boundary.names,
            pole_audit.hyperplane_intersections,
        )
        if intersection == 0
    )
    return contract_boundary_curves(boundary, vertical)


def infer_keller_dicritical_budget(
    pole_audit: KellerPoleAudit,
    name: str,
    target_curve_degree: int,
) -> KellerDicriticalBudget:
    """Infer ``(e,f,M.E)`` from ``Q``, ``p``, ``k``, and ``deg(C)``."""

    if target_curve_degree <= 0:
        raise ValueError("the target curve degree must be positive")
    boundary = pole_audit.boundary_audit.boundary
    name_to_index = {
        component_name: index
        for index, component_name in enumerate(boundary.names)
    }
    if name not in name_to_index:
        raise ValueError(f"unknown boundary component {name!r}")
    if name not in set(pole_audit.dicritical_candidates):
        raise ValueError(f"boundary component {name!r} is not dicritical")
    index = name_to_index[name]
    line_intersection = pole_audit.hyperplane_intersections[index]
    residue_degree, remainder = divmod(
        line_intersection, target_curve_degree
    )
    if remainder:
        raise ValueError(
            f"(f^*L).{name}={line_intersection} is not "
            f"divisible by deg(C)={target_curve_degree}"
        )
    if residue_degree <= 0:
        raise ValueError("a dicritical residue degree must be positive")
    matrix = boundary.intersection_matrix
    coefficients = pole_audit.ramification_coefficients
    available = sum(
        coefficients[j] * int(matrix[j, index])
        for j in range(boundary.rank)
        if j != index
    )
    return KellerDicriticalBudget(
        name=name,
        boundary_index=index,
        ramification_index=coefficients[index] + 1,
        residue_degree=residue_degree,
        available_residual_intersection=available,
        forced_companion_intersection=available - 2 * residue_degree + 2,
    )


def infer_finite_model_dicritical_projection_budget(
    pole_audit: KellerPoleAudit,
    name: str,
    target_curve_degree: int,
) -> FiniteModelDicriticalProjectionBudget:
    """Compare companion incidence on the finite Stein model.

    The uncorrected residual-different calculation assumes a smooth target
    curve.  A rational plane curve of degree ``c`` contributes normalization
    correction ``f(c-1)(c-2)``.  The returned gap checks that this correction
    exactly reconciles the residual and projection-formula budgets.
    """

    local = infer_keller_dicritical_budget(
        pole_audit, name, target_curve_degree
    )
    contraction = contract_keller_vertical_boundary(pole_audit)
    if not contraction.passes:
        raise ValueError(
            "the H-null boundary block is not contractible: "
            + "; ".join(contraction.failures)
        )
    finite_index = contraction.surviving_names.index(name)
    q_finite = contraction.finite_intersection_matrix
    self_intersection = sp.cancel(q_finite[finite_index, finite_index])
    original_names = pole_audit.boundary_audit.boundary.names
    surviving_indices = tuple(
        original_names.index(component_name)
        for component_name in contraction.surviving_names
    )
    canonical = tuple(
        pole_audit.boundary_audit.canonical_coefficients[index]
        for index in surviving_indices
    )
    ramification = tuple(
        pole_audit.ramification_coefficients[index]
        for index in surviving_indices
    )
    canonical_intersection = sp.cancel(
        sum(
            canonical[j] * q_finite[j, finite_index]
            for j in range(len(surviving_indices))
        )
    )
    surface_different = sp.cancel(
        canonical_intersection + self_intersection + 2
    )
    available = sp.cancel(
        sum(
            ramification[j] * q_finite[j, finite_index]
            for j in range(len(surviving_indices))
            if j != finite_index
        )
    )
    residual_companion = sp.cancel(
        available
        - 2 * local.residue_degree
        + 2
        - surface_different
    )
    projection_companion = sp.cancel(
        local.residue_degree * target_curve_degree**2
        - local.ramification_index * self_intersection
    )
    target_correction = (
        local.residue_degree
        * (target_curve_degree - 1)
        * (target_curve_degree - 2)
    )
    corrected_gap = sp.cancel(
        residual_companion
        + target_correction
        - projection_companion
    )
    return FiniteModelDicriticalProjectionBudget(
        name=name,
        target_curve_degree=target_curve_degree,
        ramification_index=local.ramification_index,
        residue_degree=local.residue_degree,
        finite_self_intersection=self_intersection,
        surface_different_degree=surface_different,
        residual_forced_companion_intersection=residual_companion,
        projection_forced_companion_intersection=projection_companion,
        target_normalization_correction=target_correction,
        corrected_budget_gap=corrected_gap,
        contraction=contraction,
    )


def audit_keller_residual_different(
    pole_audit: KellerPoleAudit,
    components: Sequence[KellerDicriticalDatum],
) -> ResidualDifferentAudit:
    """Apply the clean residual-different identity to intrinsic dicriticals.

    For a component ``E_i`` mapping to a target curve ``C`` of degree ``c``,

    ``(f^*L).E_i = f_i c`` and ``r_i = e_i-1``.

    The intrinsic boundary package therefore supplies ``e_i`` and ``f_i``.
    This adapter rejects non-dicritical names and nonintegral residue degrees
    before forwarding the complete boundary graph to the exact weighted-edge
    audit.
    """

    boundary = pole_audit.boundary_audit.boundary
    residual_components: list[ResidualDifferentComponent] = []
    seen_names: set[str] = set()
    for component in components:
        if component.name in seen_names:
            raise ValueError("each dicritical component may be audited only once")
        seen_names.add(component.name)
        budget = infer_keller_dicritical_budget(
            pole_audit,
            component.name,
            component.target_curve_degree,
        )
        residual_components.append(
            ResidualDifferentComponent(
                name=component.name,
                boundary_index=budget.boundary_index,
                ramification_index=budget.ramification_index,
                residue_degree=budget.residue_degree,
                companion_intersection=component.companion_intersection,
            )
        )

    return audit_residual_different(
        boundary.intersection_matrix.tolist(),
        pole_audit.ramification_coefficients,
        residual_components,
    )
