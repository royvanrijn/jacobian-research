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


@dataclass(frozen=True)
class FourTermCellularCochainComplex:
    """A cellular cochain complex through genuine three-cells."""

    d0: sp.Matrix
    d1: sp.Matrix
    d2: sp.Matrix
    name: str

    def __post_init__(self) -> None:
        d0 = sp.Matrix(self.d0)
        d1 = sp.Matrix(self.d1)
        d2 = sp.Matrix(self.d2)
        object.__setattr__(self, "d0", d0)
        object.__setattr__(self, "d1", d1)
        object.__setattr__(self, "d2", d2)
        if d0.rows != d1.cols or d1.rows != d2.cols:
            raise ValueError("adjacent cellular differentials do not compose")
        if d1 * d0 != sp.zeros(d1.rows, d0.cols):
            raise ValueError("d1*d0 is nonzero")
        if d2 * d1 != sp.zeros(d2.rows, d1.cols):
            raise ValueError("d2*d1 is nonzero")

    @property
    def dimensions(self) -> tuple[int, int, int, int]:
        """Return dimensions in cohomological degrees zero through three."""

        return (self.d0.cols, self.d0.rows, self.d1.rows, self.d2.rows)

    @property
    def ranks(self) -> tuple[int, int, int]:
        """Return exact rational differential ranks."""

        return (self.d0.rank(), self.d1.rank(), self.d2.rank())

    @property
    def cohomology_dimensions(self) -> tuple[int, int, int, int]:
        """Return exact rational cohomology dimensions."""

        c0, c1, c2, c3 = self.dimensions
        r0, r1, r2 = self.ranks
        return (
            c0 - r0,
            c1 - r0 - r1,
            c2 - r1 - r2,
            c3 - r2,
        )


@dataclass(frozen=True)
class FiniteChainComplex:
    """A finite rational homological chain complex.

    ``boundaries[n - 1]`` is the differential from degree ``n`` to degree
    ``n - 1``.  This small exact-linear-algebra model is used to certify
    comparisons between normalized bar chains and cellular chains.
    """

    dimensions: tuple[int, ...]
    boundaries: tuple[sp.Matrix, ...]
    name: str

    def __post_init__(self) -> None:
        if not self.dimensions:
            raise ValueError("a chain complex needs at least degree zero")
        if len(self.boundaries) != len(self.dimensions) - 1:
            raise ValueError("one boundary is required between adjacent degrees")
        normalized = tuple(sp.Matrix(boundary) for boundary in self.boundaries)
        object.__setattr__(self, "boundaries", normalized)
        for degree, boundary in enumerate(normalized, start=1):
            expected = (self.dimensions[degree - 1], self.dimensions[degree])
            if boundary.shape != expected:
                raise ValueError(
                    f"degree-{degree} boundary has shape {boundary.shape}, "
                    f"expected {expected}"
                )
        for lower, upper in zip(normalized, normalized[1:]):
            if lower * upper != sp.zeros(lower.rows, upper.cols):
                raise ValueError("adjacent chain boundaries do not compose to zero")

    @property
    def homology_dimensions(self) -> tuple[int, ...]:
        """Return exact rational homology dimensions in every degree."""

        ranks = tuple(boundary.rank() for boundary in self.boundaries)
        return tuple(
            dimension
            - (ranks[degree - 1] if degree else 0)
            - (ranks[degree] if degree < len(ranks) else 0)
            for degree, dimension in enumerate(self.dimensions)
        )

    def tensor_with_vector_space(
        self,
        dimension: int,
        coefficient_name: str,
    ) -> "FiniteChainComplex":
        """Tensor the complex by an exact finite rational coefficient."""

        if dimension < 0:
            raise ValueError("coefficient dimension must be nonnegative")
        identity = sp.eye(dimension)
        return FiniteChainComplex(
            tuple(size * dimension for size in self.dimensions),
            tuple(
                sp.kronecker_product(boundary, identity)
                for boundary in self.boundaries
            ),
            f"{self.name} tensor {coefficient_name}",
        )


@dataclass(frozen=True)
class FiniteChainComparison:
    """A chain map together with an exact mapping-cone quasi-isomorphism test."""

    source: FiniteChainComplex
    target: FiniteChainComplex
    maps: tuple[sp.Matrix, ...]
    name: str

    def __post_init__(self) -> None:
        if len(self.source.dimensions) != len(self.target.dimensions):
            raise ValueError("source and target must have the same degree range")
        if len(self.maps) != len(self.source.dimensions):
            raise ValueError("one comparison map is required in every degree")
        normalized = tuple(sp.Matrix(map_) for map_ in self.maps)
        object.__setattr__(self, "maps", normalized)
        for degree, map_ in enumerate(normalized):
            expected = (
                self.target.dimensions[degree],
                self.source.dimensions[degree],
            )
            if map_.shape != expected:
                raise ValueError(
                    f"degree-{degree} comparison has shape {map_.shape}, "
                    f"expected {expected}"
                )
        for degree in range(1, len(normalized)):
            if (
                self.target.boundaries[degree - 1] * normalized[degree]
                != normalized[degree - 1]
                * self.source.boundaries[degree - 1]
            ):
                raise ValueError(f"{self.name} is not a chain map in degree {degree}")

    @property
    def mapping_cone(self) -> FiniteChainComplex:
        """Return the homological mapping cone of the comparison map."""

        source_dimensions = self.source.dimensions
        target_dimensions = self.target.dimensions
        maximum_degree = len(source_dimensions)
        dimensions = tuple(
            (target_dimensions[degree] if degree < maximum_degree else 0)
            + (source_dimensions[degree - 1] if degree > 0 else 0)
            for degree in range(maximum_degree + 1)
        )

        def target_boundary(degree: int) -> sp.Matrix:
            if 1 <= degree < maximum_degree:
                return self.target.boundaries[degree - 1]
            return sp.zeros(
                target_dimensions[degree - 1] if degree > 0 else 0,
                target_dimensions[degree] if degree < maximum_degree else 0,
            )

        def source_boundary(degree: int) -> sp.Matrix:
            if 1 <= degree < maximum_degree:
                return self.source.boundaries[degree - 1]
            return sp.zeros(
                source_dimensions[degree - 1] if degree > 0 else 0,
                source_dimensions[degree] if degree < maximum_degree else 0,
            )

        boundaries = []
        for degree in range(1, maximum_degree + 1):
            target_left = (
                target_dimensions[degree - 1]
                if degree - 1 < maximum_degree
                else 0
            )
            target_right = (
                target_dimensions[degree]
                if degree < maximum_degree
                else 0
            )
            source_left = (
                source_dimensions[degree - 2] if degree >= 2 else 0
            )
            source_right = source_dimensions[degree - 1]
            upper = target_boundary(degree).row_join(
                self.maps[degree - 1]
            )
            lower = sp.zeros(source_left, target_right).row_join(
                -source_boundary(degree - 1)
            )
            boundary = upper.col_join(lower)
            assert boundary.shape == (
                target_left + source_left,
                target_right + source_right,
            )
            boundaries.append(boundary)
        return FiniteChainComplex(
            dimensions,
            tuple(boundaries),
            f"cone({self.name})",
        )

    @property
    def is_quasi_isomorphism(self) -> bool:
        """Whether the exact rational mapping cone is acyclic."""

        return all(
            dimension == 0
            for dimension in self.mapping_cone.homology_dimensions
        )

    def tensor_with_vector_space(
        self,
        dimension: int,
        coefficient_name: str,
    ) -> "FiniteChainComparison":
        """Tensor a chain comparison by a finite rational coefficient."""

        if dimension < 0:
            raise ValueError("coefficient dimension must be nonnegative")
        identity = sp.eye(dimension)
        return FiniteChainComparison(
            self.source.tensor_with_vector_space(dimension, coefficient_name),
            self.target.tensor_with_vector_space(dimension, coefficient_name),
            tuple(
                sp.kronecker_product(map_, identity) for map_ in self.maps
            ),
            f"{self.name} tensor {coefficient_name}",
        )


def minimal_top_cell_completion(
    model: CellularCochainModel,
    name: str,
) -> FourTermCellularCochainComplex:
    """Attach algebraic top cells whose coboundary kills all two-skeleton H2.

    The rows of ``d2`` form the exact annihilator of the image of ``d1``.
    For a known regular CW boundary this recovers the oriented top-cell
    attachment.  Without such geometric input it is only a minimal algebraic
    completion, not an assertion that the required cells exist.
    """

    annihilator = model.complex.d1.T.nullspace()
    d2 = (
        sp.Matrix.vstack(*(vector.T for vector in annihilator))
        if annihilator
        else sp.zeros(0, model.complex.dimensions[2])
    )
    return FourTermCellularCochainComplex(
        model.complex.d0,
        model.complex.d1,
        d2,
        name,
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
