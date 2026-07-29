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
triangle splits.  The finite-module classes below retain the extension tower
itself: commuting coordinate actions, equivariant quotient maps, their exact
kernel layers, and compatible-section tests.
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


@dataclass(frozen=True)
class FiniteModuleRepresentation:
    """A finite rational module represented by commuting coordinate actions."""

    variable_names: tuple[str, ...]
    actions: tuple[sp.Matrix, ...]
    name: str

    def __post_init__(self) -> None:
        if len(self.variable_names) != len(self.actions):
            raise ValueError("each coordinate needs one action matrix")
        normalized = tuple(sp.Matrix(action) for action in self.actions)
        dimensions = {
            action.rows
            for action in normalized
            if action.rows == action.cols
        }
        if any(action.rows != action.cols for action in normalized):
            raise ValueError("module actions must be square")
        if len(dimensions) > 1:
            raise ValueError("all module actions must have the same size")
        object.__setattr__(self, "actions", normalized)
        if not self.actions:
            return
        for index, left in enumerate(self.actions):
            for right in self.actions[index + 1 :]:
                if left * right != right * left:
                    raise ValueError(
                        f"coordinate actions do not commute on {self.name}"
                    )

    @property
    def dimension(self) -> int:
        """Return the underlying rational vector-space dimension."""

        return self.actions[0].rows if self.actions else 0

    @classmethod
    def zero(
        cls,
        variable_names: tuple[str, ...],
        name: str,
    ) -> "FiniteModuleRepresentation":
        """Return the zero module over the declared coordinate algebra."""

        return cls(
            variable_names,
            tuple(sp.zeros(0, 0) for _ in variable_names),
            name,
        )


@dataclass(frozen=True)
class EquivariantModuleSurjection:
    """An equivariant surjection between finite module representations."""

    source: FiniteModuleRepresentation
    target: FiniteModuleRepresentation
    projection: sp.Matrix
    name: str

    def __post_init__(self) -> None:
        projection = sp.Matrix(self.projection)
        object.__setattr__(self, "projection", projection)
        if self.source.variable_names != self.target.variable_names:
            raise ValueError("source and target coordinate names differ")
        expected = (self.target.dimension, self.source.dimension)
        if projection.shape != expected:
            raise ValueError(
                f"projection shape {projection.shape} does not match {expected}"
            )
        if projection.rank() != self.target.dimension:
            raise ValueError(f"{self.name} is not surjective")
        for source_action, target_action in zip(
            self.source.actions, self.target.actions
        ):
            if projection * source_action != target_action * projection:
                raise ValueError(f"{self.name} is not equivariant")

    @property
    def kernel_basis(self) -> sp.Matrix:
        """Return a column basis for the invariant kernel."""

        nullspace = self.projection.nullspace()
        if not nullspace:
            return sp.zeros(self.source.dimension, 0)
        return sp.Matrix.hstack(*nullspace)

    @property
    def kernel(self) -> FiniteModuleRepresentation:
        """Return the exact kernel with its restricted coordinate actions."""

        basis = self.kernel_basis
        restricted = []
        for action in self.source.actions:
            if basis.cols == 0:
                restricted.append(sp.zeros(0, 0))
                continue
            restricted.append(
                sp.Matrix.hstack(
                    *(
                        basis.gauss_jordan_solve(
                            action * basis[:, column]
                        )[0]
                        for column in range(basis.cols)
                    )
                )
            )
        return FiniteModuleRepresentation(
            self.source.variable_names,
            tuple(restricted),
            f"kernel({self.name})",
        )

    def compatible_section(self) -> sp.Matrix | None:
        """Return one equivariant section, or ``None`` when none exists."""

        source_dimension = self.source.dimension
        target_dimension = self.target.dimension
        if target_dimension == 0:
            return sp.zeros(source_dimension, 0)
        variables = sp.symbols(
            f"section_0:{source_dimension * target_dimension}"
        )
        section = sp.Matrix(
            source_dimension,
            target_dimension,
            variables,
        )
        equations = list(
            self.projection * section - sp.eye(target_dimension)
        )
        for source_action, target_action in zip(
            self.source.actions, self.target.actions
        ):
            equations.extend(
                source_action * section - section * target_action
            )
        coefficient_matrix, right_hand_side = sp.linear_eq_to_matrix(
            equations, variables
        )
        try:
            solution, parameters = coefficient_matrix.gauss_jordan_solve(
                right_hand_side
            )
        except ValueError:
            return None
        values = solution.subs(
            {parameter: sp.Integer(0) for parameter in parameters}
        )
        return sp.Matrix(
            source_dimension,
            target_dimension,
            list(values),
        )

    @property
    def splits(self) -> bool:
        """Whether the surjection has an equivariant linear section."""

        return self.compatible_section() is not None


@dataclass(frozen=True)
class PostnikovModuleTower:
    """A composable tower of finite equivariant module surjections."""

    maps: tuple[EquivariantModuleSurjection, ...]
    name: str

    def __post_init__(self) -> None:
        for left, right in zip(self.maps, self.maps[1:]):
            if (
                left.target.variable_names
                != right.source.variable_names
                or left.target.dimension != right.source.dimension
                or left.target.actions != right.source.actions
            ):
                raise ValueError(
                    f"noncomposable adjacent maps in {self.name}"
                )

    @property
    def module_dimensions(self) -> tuple[int, ...]:
        """Dimensions of all modules, including the terminal target."""

        if not self.maps:
            return ()
        return (
            self.maps[0].source.dimension,
            *(map_.target.dimension for map_ in self.maps),
        )

    @property
    def layer_dimensions(self) -> tuple[int, ...]:
        """Dimensions of the exact kernel layers."""

        return tuple(map_.kernel.dimension for map_ in self.maps)

    @property
    def split_profile(self) -> tuple[bool, ...]:
        """Compatible-section status for every adjacent extension."""

        return tuple(map_.splits for map_ in self.maps)


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


def postnikov_braid_totalization(
    braid: RittMoveComplex,
    base_dimension: int,
    tower: PostnikovModuleTower,
    layer_names: tuple[str, ...],
    name: str,
) -> CellularCochainModel:
    """Totalize the associated graded of an arbitrary Postnikov tower."""

    if len(layer_names) != len(tower.maps):
        raise ValueError("each Postnikov layer needs a cellular label")
    return braid_totalization(
        braid,
        base_dimension,
        tower.layer_dimensions,
        layer_names,
        name,
    )
