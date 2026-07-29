"""Cellular linear algebra for coefficient-decorated Ritt complexes.

The coefficient-free Ritt complex only records vertices, moves, and Coxeter
cells.  This module supplies the small piece of the Hessian--Ritt cotangent
prototype that is already justified by the degree-30 and degree-42
calculations:

* a constant coefficient block on the reduced Dickson component; and
* coefficient blocks supported on a half-braid relative to its endpoints.

For a three-edge path ``v0-v1-v2-v3``, relative cellular cochains modulo the
endpoints are

    D^{v1,v2} -> D^{e01,e12,e23}.

Their only cohomology is one copy of ``D`` in degree one.  This is the
cellular location of a path-to-boundary conormal module.  The construction is
valid for any coefficient module ``D``; finite-dimensional matrices below
are used for fibers and exact Artin jets.

This is an associated-graded prototype.  Direct-summing two relative blocks
does not assert that the corresponding completed cotangent transitivity
triangle splits.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from jcsearch.deformation_complex import ThreeTermComplex
from jcsearch.ritt_complex import RittMoveComplex, Word


@dataclass(frozen=True)
class CellularCochainModel:
    """A three-term cellular cochain complex with labelled summands."""

    complex: ThreeTermComplex
    degree_zero_labels: tuple[str, ...]
    degree_one_labels: tuple[str, ...]
    degree_two_labels: tuple[str, ...]

    def __post_init__(self) -> None:
        expected = self.complex.dimensions
        actual = (
            len(self.degree_zero_labels),
            len(self.degree_one_labels),
            len(self.degree_two_labels),
        )
        if actual != expected:
            raise ValueError(
                f"label dimensions {actual} do not match complex {expected}"
            )


def ritt_cellular_coboundaries(
    complex_: RittMoveComplex,
) -> CellularCochainModel:
    """Return the ordinary cellular cochains of a Ritt two-complex."""

    vertices = tuple(sorted(complex_.words))
    vertex_index = {word: index for index, word in enumerate(vertices)}
    edges = tuple(complex_.edges)
    edge_index = {
        frozenset(edge.endpoints): index for index, edge in enumerate(edges)
    }

    d0 = sp.zeros(len(edges), len(vertices))
    for row, edge in enumerate(edges):
        left, right = edge.endpoints
        d0[row, vertex_index[left]] = -1
        d0[row, vertex_index[right]] = 1

    d1 = sp.zeros(len(complex_.two_cells), len(edges))
    for row, cell in enumerate(complex_.two_cells):
        first, second = cell.paths
        oriented_boundary = list(zip(first, first[1:])) + list(
            zip(reversed(second), reversed(second[:-1]))
        )
        for left, right in oriented_boundary:
            column = edge_index[frozenset((left, right))]
            stored_left, stored_right = edges[column].endpoints
            d1[row, column] += (
                1 if (left, right) == (stored_left, stored_right) else -1
            )

    return CellularCochainModel(
        ThreeTermComplex(d0, d1, "Ritt cellular cochains"),
        tuple(f"vertex:{word_label(word)}" for word in vertices),
        tuple(
            "move:"
            + word_label(edge.endpoints[0])
            + "->"
            + word_label(edge.endpoints[1])
            for edge in edges
        ),
        tuple(
            f"{cell.relation}-cell:{index}"
            for index, cell in enumerate(complex_.two_cells)
        ),
    )


def word_label(word: Word) -> str:
    """Use the repository's compact outer-to-inner word notation."""

    return "".join(map(str, word))


def tensor_with_vector_space(
    model: CellularCochainModel,
    dimension: int,
    coefficient_name: str,
) -> CellularCochainModel:
    """Tensor a cellular model with an exact vector space over ``Q``."""

    if dimension < 0:
        raise ValueError("coefficient dimension must be nonnegative")
    identity = sp.eye(dimension)
    d0 = sp.kronecker_product(model.complex.d0, identity)
    d1 = sp.kronecker_product(model.complex.d1, identity)

    def expanded(labels: tuple[str, ...]) -> tuple[str, ...]:
        return tuple(
            f"{label}|{coefficient_name}[{index}]"
            for label in labels
            for index in range(dimension)
        )

    return CellularCochainModel(
        ThreeTermComplex(
            d0,
            d1,
            f"{model.complex.name} tensor {coefficient_name}",
        ),
        expanded(model.degree_zero_labels),
        expanded(model.degree_one_labels),
        expanded(model.degree_two_labels),
    )


def relative_half_path_model(
    coefficient_dimension: int,
    coefficient_name: str,
    path_name: str,
) -> CellularCochainModel:
    """Return cochains of a three-edge half-braid relative to its endpoints."""

    if coefficient_dimension < 0:
        raise ValueError("coefficient dimension must be nonnegative")
    # The scalar relative coboundary sends (v1,v2) to
    # (v1, -v1+v2, -v2).  It is injective with one-dimensional cokernel.
    scalar_d0 = sp.Matrix(
        [
            [1, 0],
            [-1, 1],
            [0, -1],
        ]
    )
    identity = sp.eye(coefficient_dimension)
    d0 = sp.kronecker_product(scalar_d0, identity)
    d1 = sp.zeros(0, 3 * coefficient_dimension)

    vertex_labels = tuple(
        f"{path_name}:interior-vertex-{vertex}|{coefficient_name}[{index}]"
        for vertex in (1, 2)
        for index in range(coefficient_dimension)
    )
    edge_labels = tuple(
        f"{path_name}:move-{edge}|{coefficient_name}[{index}]"
        for edge in range(3)
        for index in range(coefficient_dimension)
    )
    return CellularCochainModel(
        ThreeTermComplex(d0, d1, f"relative half-path {path_name}"),
        vertex_labels,
        edge_labels,
        (),
    )


def direct_sum(
    models: tuple[CellularCochainModel, ...],
    name: str,
) -> CellularCochainModel:
    """Form a direct sum while preserving the cohomological degrees."""

    if not models:
        empty = sp.zeros(0, 0)
        return CellularCochainModel(
            ThreeTermComplex(empty, empty, name), (), (), ()
        )

    c0 = sum(model.complex.dimensions[0] for model in models)
    c1 = sum(model.complex.dimensions[1] for model in models)
    c2 = sum(model.complex.dimensions[2] for model in models)
    d0 = sp.zeros(c1, c0)
    d1 = sp.zeros(c2, c1)

    c0_offset = c1_offset = c2_offset = 0
    for model in models:
        m0, m1, m2 = model.complex.dimensions
        d0[
            c1_offset : c1_offset + m1,
            c0_offset : c0_offset + m0,
        ] = model.complex.d0
        d1[
            c2_offset : c2_offset + m2,
            c1_offset : c1_offset + m1,
        ] = model.complex.d1
        c0_offset += m0
        c1_offset += m1
        c2_offset += m2

    return CellularCochainModel(
        ThreeTermComplex(d0, d1, name),
        tuple(
            label for model in models for label in model.degree_zero_labels
        ),
        tuple(
            label for model in models for label in model.degree_one_labels
        ),
        tuple(
            label for model in models for label in model.degree_two_labels
        ),
    )


def braid_totalization(
    braid: RittMoveComplex,
    base_dimension: int,
    defect_dimensions: tuple[int, ...],
    defect_names: tuple[str, ...],
    name: str,
) -> CellularCochainModel:
    """Build the reduced-base plus relative-defect cellular totalization."""

    if len(defect_dimensions) != len(defect_names):
        raise ValueError("each defect dimension needs a name")
    base = tensor_with_vector_space(
        ritt_cellular_coboundaries(braid),
        base_dimension,
        "Dickson-base",
    )
    defects = tuple(
        relative_half_path_model(dimension, defect_name, f"half-{index}")
        for index, (dimension, defect_name) in enumerate(
            zip(defect_dimensions, defect_names)
        )
    )
    return direct_sum((base,) + defects, name)
