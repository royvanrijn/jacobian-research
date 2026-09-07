#!/usr/bin/env python3
"""Verify the skeletal linear algebra in Hessian--Ritt cotangent descent."""

from __future__ import annotations

import json
import sys
from itertools import permutations
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.hessian_ritt_cellular import (  # noqa: E402
    FiniteChainComparison,
    FiniteChainComplex,
    braid_totalization,
    minimal_top_cell_completion,
    ritt_cellular_coboundaries,
    tensor_with_vector_space,
)
from jcsearch.ritt_complex import (  # noqa: E402
    MoveType,
    compose_factors,
    dickson,
    dickson_vertex_factors,
    permutation_ritt_complex,
    symmetric_braid_complex,
)


def face_poset_bar_complex(
    cells: tuple[str, ...],
    dimensions: dict[str, int],
    proper_faces: dict[str, tuple[str, ...]],
    *,
    relative_cells: frozenset[str] = frozenset(),
    name: str,
) -> tuple[FiniteChainComplex, tuple[tuple[tuple[str, ...], ...], ...]]:
    """Return normalized chains of a finite face-poset bar construction.

    A degree-``q`` basis element is a strict chain of ``q+1`` comparable
    cells.  Quotienting by ``relative_cells`` removes the simplices entirely
    contained in that face subcomplex.
    """

    order = {cell: index for index, cell in enumerate(cells)}

    def comparable(lower: str, upper: str) -> bool:
        return lower in proper_faces.get(upper, ())

    maximum_dimension = max(dimensions.values())
    all_chains: list[list[tuple[str, ...]]] = [
        [(cell,) for cell in cells]
    ]
    for _ in range(maximum_dimension):
        next_chains = {
            chain + (upper,)
            for chain in all_chains[-1]
            for upper in cells
            if dimensions[upper] > dimensions[chain[-1]]
            and comparable(chain[-1], upper)
        }
        all_chains.append(
            sorted(
                next_chains,
                key=lambda chain: tuple(order[cell] for cell in chain),
            )
        )
    while all_chains and not all_chains[-1]:
        all_chains.pop()
    chains = [
        [
            chain
            for chain in degree_chains
            if not all(cell in relative_cells for cell in chain)
        ]
        for degree_chains in all_chains
    ]
    bases = tuple(tuple(degree_chains) for degree_chains in chains)
    indices = tuple(
        {chain: index for index, chain in enumerate(degree_chains)}
        for degree_chains in bases
    )
    boundaries = []
    for degree in range(1, len(bases)):
        boundary = sp.zeros(len(bases[degree - 1]), len(bases[degree]))
        for column, chain in enumerate(bases[degree]):
            for removed in range(len(chain)):
                face = chain[:removed] + chain[removed + 1 :]
                if all(cell in relative_cells for cell in face):
                    continue
                boundary[indices[degree - 1][face], column] += (
                    -1 if removed % 2 else 1
                )
        boundaries.append(boundary)
    return (
        FiniteChainComplex(
            tuple(len(basis) for basis in bases),
            tuple(boundaries),
            name,
        ),
        bases,
    )


def relative_path_subdivision_comparison() -> FiniteChainComparison:
    """Compare a three-edge relative cell path with its face-poset bar."""

    vertices = tuple(f"v{index}" for index in range(4))
    edges = tuple(f"e{index}" for index in range(3))
    cells = vertices + edges
    dimensions = {
        **{vertex: 0 for vertex in vertices},
        **{edge: 1 for edge in edges},
    }
    proper_faces = {
        edge: (vertices[index], vertices[index + 1])
        for index, edge in enumerate(edges)
    }
    target, bases = face_poset_bar_complex(
        cells,
        dimensions,
        proper_faces,
        relative_cells=frozenset((vertices[0], vertices[-1])),
        name="normalized face-poset bar of a relative half-braid",
    )

    d1 = sp.Matrix(
        [
            [1, -1, 0],
            [0, 1, -1],
        ]
    )
    source = FiniteChainComplex(
        (2, 3),
        (d1,),
        "cellular chains of a relative half-braid",
    )
    degree_zero_index = {
        chain: index for index, chain in enumerate(bases[0])
    }
    degree_one_index = {
        chain: index for index, chain in enumerate(bases[1])
    }
    map_zero = sp.zeros(target.dimensions[0], source.dimensions[0])
    for column, vertex in enumerate(vertices[1:-1]):
        map_zero[degree_zero_index[(vertex,)], column] = 1
    map_one = sp.zeros(target.dimensions[1], source.dimensions[1])
    for column, edge in enumerate(edges):
        left, right = proper_faces[edge]
        map_one[degree_one_index[(left, edge)], column] = 1
        map_one[degree_one_index[(right, edge)], column] = -1
    return FiniteChainComparison(
        source,
        target,
        (map_zero, map_one),
        "relative half-braid cellular subdivision",
    )


def filled_braid_subdivision_comparison() -> FiniteChainComparison:
    """Compare the filled braid CW disk with its face-poset bar chains."""

    vertices = tuple(f"v{index}" for index in range(6))
    edges = tuple(f"e{index}" for index in range(6))
    face = "f"
    cells = vertices + edges + (face,)
    dimensions = {
        **{vertex: 0 for vertex in vertices},
        **{edge: 1 for edge in edges},
        face: 2,
    }
    edge_faces = {
        edge: (vertices[index], vertices[(index + 1) % len(vertices)])
        for index, edge in enumerate(edges)
    }
    proper_faces = {
        **edge_faces,
        face: vertices + edges,
    }
    target, bases = face_poset_bar_complex(
        cells,
        dimensions,
        proper_faces,
        name="normalized face-poset bar of the filled braid",
    )

    d1 = sp.zeros(6, 6)
    for column in range(6):
        d1[column, column] = -1
        d1[(column + 1) % 6, column] = 1
    d2 = sp.ones(6, 1)
    source = FiniteChainComplex(
        (6, 6, 1),
        (d1, d2),
        "cellular chains of the filled braid",
    )
    indices = tuple(
        {chain: index for index, chain in enumerate(basis)}
        for basis in bases
    )
    map_zero = sp.zeros(target.dimensions[0], 6)
    for column, vertex in enumerate(vertices):
        map_zero[indices[0][(vertex,)], column] = 1
    map_one = sp.zeros(target.dimensions[1], 6)
    for column, edge in enumerate(edges):
        left, right = edge_faces[edge]
        map_one[indices[1][(left, edge)], column] = 1
        map_one[indices[1][(right, edge)], column] = -1
    map_two = sp.zeros(target.dimensions[2], 1)
    for edge in edges:
        left, right = edge_faces[edge]
        map_two[indices[2][(left, edge, face)], 0] = 1
        map_two[indices[2][(right, edge, face)], 0] = -1
    return FiniteChainComparison(
        source,
        target,
        (map_zero, map_one, map_two),
        "filled braid cellular subdivision",
    )


def degree_42_hessian_tangent_image(
    word: tuple[int, int, int],
) -> sp.Matrix:
    """Return the normalized chart tangent image in Hessian coefficients."""

    variable = sp.Symbol("W")
    factors = dickson_vertex_factors(
        word,
        variable,
        sp.Integer(1),
        sp.Integer(0),
    )
    columns = []
    for index, degree in enumerate(word):
        outer = (
            compose_factors(factors[:index], variable)
            if index
            else variable
        )
        inner = (
            compose_factors(factors[index + 1 :], variable)
            if index + 1 < len(factors)
            else variable
        )
        full_argument = sp.expand(factors[index].subs(variable, inner))
        multiplier = sp.expand(
            sp.diff(outer, variable).subs(variable, full_argument)
        )
        for power in range(1, degree):
            variation = sp.Poly(
                sp.expand(multiplier * inner**power),
                variable,
            )
            # Hessian projection forgets the degree-one coefficient.
            columns.append(
                sp.Matrix(
                    [variation.nth(target) for target in range(2, 42)]
                )
            )
    return sp.Matrix.hstack(*columns)


def subspace_intersection(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    """Return an exact column basis for two rational subspaces."""

    kernel = left.row_join(-right).nullspace()
    if not kernel:
        return sp.zeros(left.rows, 0)
    candidates = sp.Matrix.hstack(
        *(left * vector[: left.cols, :] for vector in kernel)
    )
    basis = candidates.columnspace()
    return (
        sp.Matrix.hstack(*basis)
        if basis
        else sp.zeros(left.rows, 0)
    )


def audit_degree_42_actual_tangent_incidence() -> dict[str, object]:
    """Compute the actual six-chart Hessian tangent incidence diagram."""

    words = tuple(sorted(permutations((2, 3, 7))))
    tangent_images = {
        word: degree_42_hessian_tangent_image(word) for word in words
    }
    vertex_ranks = {
        "".join(map(str, word)): image.rank()
        for word, image in tangent_images.items()
    }
    assert set(vertex_ranks.values()) == {9}

    edge_intersection_ranks = {}
    for word in words:
        for swapped_index in range(2):
            target_list = list(word)
            target_list[swapped_index], target_list[swapped_index + 1] = (
                target_list[swapped_index + 1],
                target_list[swapped_index],
            )
            target = tuple(target_list)
            if word >= target:
                continue
            intersection = subspace_intersection(
                tangent_images[word], tangent_images[target]
            )
            label = (
                "".join(map(str, word))
                + "-"
                + "".join(map(str, target))
            )
            edge_intersection_ranks[label] = intersection.cols
    assert edge_intersection_ranks == {
        "237-327": 8,
        "237-273": 5,
        "273-723": 6,
        "327-372": 6,
        "372-732": 5,
        "723-732": 8,
    }

    common = tangent_images[words[0]]
    for word in words[1:]:
        common = subspace_intersection(common, tangent_images[word])
    assert common.cols == 3

    variable, translation, parameter = sp.symbols("W translation parameter")
    target = sp.Poly(
        sp.expand(
            dickson(42, variable + translation, parameter)
            - dickson(42, translation, parameter)
        ),
        variable,
    )
    coefficients = sp.Matrix(
        [target.nth(power) for power in range(2, 42)]
    )
    dickson_tangent = coefficients.jacobian(
        (translation, parameter)
    ).subs({translation: 1, parameter: 0})
    assert dickson_tangent.rank() == 2
    assert common.row_join(dickson_tangent).rank() == common.rank()

    excess_polynomial = sp.Poly((variable + 1) ** 36 - 1, variable)
    excess = sp.Matrix(
        [excess_polynomial.nth(power) for power in range(2, 42)]
    )
    common_model = dickson_tangent.row_join(excess)
    assert common_model.rank() == 3
    assert common.row_join(common_model).rank() == 3

    braid = symmetric_braid_complex((2, 3, 7), MoveType.CHEBYSHEV)
    vertex_by_word = {vertex.word: vertex for vertex in braid.vertices}
    adjacency = {word: [] for word in words}
    for edge in braid.edges:
        left, right = edge.endpoints
        adjacency[left].append(right)
        adjacency[right].append(left)

    def half_braid_paths(
        word: tuple[int, int, int],
    ) -> tuple[tuple[tuple[int, int, int], ...], ...]:
        endpoint = tuple(reversed(word))
        paths = []

        def extend(path: tuple[tuple[int, int, int], ...]) -> None:
            if len(path) == 4:
                if path[-1] == endpoint:
                    paths.append(path)
                return
            for neighbor in adjacency[path[-1]]:
                if neighbor not in path:
                    extend(path + (neighbor,))

        extend((word,))
        assert len(paths) == 2
        return tuple(paths)

    def intersection_dimension(
        path: tuple[tuple[int, int, int], ...],
    ) -> int:
        intersection = tangent_images[path[0]]
        for vertex in path[1:]:
            intersection = subspace_intersection(
                intersection, tangent_images[vertex]
            )
        return intersection.cols

    all_cuts = frozenset((2, 3, 6, 7, 14, 21))
    composite_cuts = frozenset((6, 14, 21))
    sector_flags = {}
    for representative in ((2, 7, 3), (2, 3, 7), (3, 2, 7)):
        path_data = []
        for path in half_braid_paths(representative):
            cuts = frozenset().union(
                *(
                    frozenset(vertex_by_word[vertex].cuts)
                    for vertex in path
                )
            )
            omitted = all_cuts - cuts
            assert len(omitted) == 1
            path_data.append(
                (
                    next(iter(omitted)),
                    path,
                    intersection_dimension(path),
                )
            )
        thick = next(
            item for item in path_data if item[0] in composite_cuts
        )
        thin = next(
            item for item in path_data if item[0] not in composite_cuts
        )
        assert (thick[2], thin[2], common.cols, dickson_tangent.rank()) == (
            4,
            3,
            3,
            2,
        )
        label = "".join(map(str, representative))
        sector_flags[label] = {
            "opposite_chart": "".join(map(str, reversed(representative))),
            "thick_composite_omission": thick[0],
            "thin_prime_omission": thin[0],
            "tangent_dimensions": {
                "thick_path": thick[2],
                "thin_path": thin[2],
                "full_boundary": common.cols,
                "reduced_dickson_graph": dickson_tangent.rank(),
            },
            "conormal_ranks": {
                "thick_path": 9 - thick[2],
                "thin_path": 9 - thin[2],
                "full_boundary": 9 - common.cols,
                "reduced_dickson_graph": 9 - dickson_tangent.rank(),
            },
        }
    return {
        "ambient": "degree-42 Hessian coefficients W^2 through W^41",
        "vertex_tangent_ranks": vertex_ranks,
        "move_intersection_ranks": edge_intersection_ranks,
        "common_six_chart_tangent_rank": common.cols,
        "reduced_dickson_tangent_rank": dickson_tangent.rank(),
        "common_excess_generator": "(W+1)^36-1 after Hessian projection",
        "common_excess_rank": common_model.rank() - dickson_tangent.rank(),
        "three_sector_first_conormal_flags": sector_flags,
        "conclusion": (
            "the actual degree-42 tangent face data are coherent as "
            "ambient subspaces but are not a constant rank-two system: "
            "move ranks are 5, 6, or 8, all three labelled sectors have "
            "conormal ranks 5<6<7, and the six-chart intersection has one "
            "nilpotent tangent direction beyond the Dickson plane"
        ),
    }


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_ritt_cotangent_descent.json"
)


def audit_filled_braid() -> dict[str, object]:
    """The three-factor filled braid already has no missing top cell."""

    braid = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    model = ritt_cellular_coboundaries(braid)
    completion = minimal_top_cell_completion(
        model, "filled three-factor braid"
    )
    assert completion.dimensions == (6, 6, 1, 0)
    assert completion.ranks == (5, 1, 0)
    assert completion.cohomology_dimensions == (1, 0, 0, 0)
    return {
        "dimensions": completion.dimensions,
        "ranks": completion.ranks,
        "cohomology_dimensions": completion.cohomology_dimensions,
        "top_cells_added": completion.dimensions[3],
    }


def audit_permutohedron_completion() -> dict[str, object]:
    """The four-factor Coxeter boundary needs its genuine three-cell."""

    boundary = permutation_ritt_complex(
        (2, 3, 5, 7), MoveType.CHEBYSHEV
    )
    model = ritt_cellular_coboundaries(boundary)
    assert model.complex.cohomology_dimensions == (1, 0, 1)
    completion = minimal_top_cell_completion(
        model, "three-dimensional permutohedron"
    )
    assert completion.dimensions == (24, 36, 14, 1)
    assert completion.ranks == (23, 13, 1)
    assert completion.cohomology_dimensions == (1, 0, 0, 0)
    oriented_boundary = tuple(int(entry) for entry in completion.d2.tolist()[0])
    assert oriented_boundary == (
        1,
        -1,
        -1,
        1,
        -1,
        -1,
        1,
        1,
        1,
        1,
        -1,
        -1,
        1,
        1,
    )

    coefficient_model = tensor_with_vector_space(
        model, 2, "constant-perfect-coefficient"
    )
    coefficient_completion = minimal_top_cell_completion(
        coefficient_model,
        "permutohedron with a two-dimensional perfect coefficient",
    )
    assert coefficient_completion.dimensions == (48, 72, 28, 2)
    assert coefficient_completion.cohomology_dimensions == (2, 0, 0, 0)
    return {
        "two_skeleton_cohomology": (
            *model.complex.cohomology_dimensions,
            0,
        ),
        "completed_dimensions": completion.dimensions,
        "completed_ranks": completion.ranks,
        "completed_cohomology": completion.cohomology_dimensions,
        "oriented_three_cell_boundary": oriented_boundary,
        "two_dimensional_coefficient_cohomology": (
            coefficient_completion.cohomology_dimensions
        ),
    }


def audit_low_degree_skeletal_invariance() -> dict[str, object]:
    """Adding genuine three-cells changes H2 but not H0 or H1."""

    braid = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    degree_30 = braid_totalization(
        braid,
        base_dimension=2,
        defect_dimensions=(1,),
        defect_names=("sector",),
        name="degree-30 first Postnikov cellular model",
    )
    completion = minimal_top_cell_completion(
        degree_30, "degree-30 first Postnikov completion"
    )
    assert degree_30.complex.cohomology_dimensions == (2, 1, 0)
    assert completion.cohomology_dimensions == (2, 1, 0, 0)
    return {
        "two_skeleton": degree_30.complex.cohomology_dimensions,
        "with_all_minimal_top_cells": completion.cohomology_dimensions,
        "conclusion": (
            "the complete two-skeleton determines H0 and H1; H2 is "
            "trustworthy only after the next cellular differential"
        ),
    }


def audit_actual_bar_to_cellular_comparison() -> dict[str, object]:
    """Certify the normalized face-bar to cellular subdivision maps."""

    path = relative_path_subdivision_comparison()
    disk = filled_braid_subdivision_comparison()
    assert path.source.dimensions == (2, 3)
    assert path.target.dimensions == (5, 6)
    assert path.source.homology_dimensions == (0, 1)
    assert path.target.homology_dimensions == (0, 1)
    assert path.is_quasi_isomorphism
    assert disk.source.dimensions == (6, 6, 1)
    assert disk.target.dimensions == (13, 24, 12)
    assert disk.source.homology_dimensions == (1, 0, 0)
    assert disk.target.homology_dimensions == (1, 0, 0)
    assert disk.is_quasi_isomorphism

    coefficient_tests = {}
    for dimension in (2, 4, 6):
        path_with_coefficients = path.tensor_with_vector_space(
            dimension, f"D{dimension}"
        )
        disk_with_coefficients = disk.tensor_with_vector_space(
            dimension, f"D{dimension}"
        )
        assert path_with_coefficients.is_quasi_isomorphism
        assert disk_with_coefficients.is_quasi_isomorphism
        coefficient_tests[str(dimension)] = {
            "relative_path_homology": (
                path_with_coefficients.source.homology_dimensions
            ),
            "filled_braid_homology": (
                disk_with_coefficients.source.homology_dimensions
            ),
        }
    return {
        "relative_half_braid": {
            "cellular_dimensions": path.source.dimensions,
            "normalized_bar_dimensions": path.target.dimensions,
            "mapping_cone_homology": (
                path.mapping_cone.homology_dimensions
            ),
        },
        "filled_braid": {
            "cellular_dimensions": disk.source.dimensions,
            "normalized_bar_dimensions": disk.target.dimensions,
            "mapping_cone_homology": (
                disk.mapping_cone.homology_dimensions
            ),
        },
        "coefficient_dimensions_tested": coefficient_tests,
        "conclusion": (
            "the canonical cellular subdivision maps are "
            "quasi-isomorphisms; tensoring preserves the comparison for "
            "every exact perfect coefficient block"
        ),
    }


def main() -> None:
    output = {
        "schema": "hessian-ritt-cotangent-descent.v2",
        "status": "exact skeletal descent regression",
        "filled_braid": audit_filled_braid(),
        "four_factor_permutohedron": audit_permutohedron_completion(),
        "low_degree_skeletal_invariance": (
            audit_low_degree_skeletal_invariance()
        ),
        "bar_to_cellular_subdivision": (
            audit_actual_bar_to_cellular_comparison()
        ),
        "degree_42_actual_tangent_incidence": (
            audit_degree_42_actual_tangent_incidence()
        ),
        "formal_theorem": (
            "cotangent complexes commute with homotopy colimits of "
            "connective algebra diagrams; for perfect local cotangent "
            "modules, duality converts the cellular chain colimit into "
            "the cellular cochain totalization. Once the actual bar "
            "coefficient diagram factors coherently through the Ritt face "
            "category, homotopy cofinality is the remaining presentation "
            "criterion; face-poset bar to cellular subdivision is formal"
        ),
        "proved_comparison": (
            "the full simplicial bar diagram of the actual derived "
            "Hessian intersection satisfies cotangent descent"
        ),
        "remaining_effectivity_condition": (
            "construct coherent completed coefficient transport from the "
            "actual bar diagram to the Ritt face category and prove the "
            "bar-to-face functor homotopy cofinal on every completed "
            "component; degree 30 supplies all conormal-fiber braid "
            "charts, while degree 42 currently supplies one exact flag "
            "chart"
        ),
        "theorem_boundary": (
            "the checker proves the cellular and skeletal linear algebra; "
            "the canonical note proves formal cotangent-colimit descent. "
            "It does not assert all-degree coefficient effectivity or "
            "filtered algebraic H2 vanishing."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_hessian_ritt_cotangent_descent.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: the filled braid is already complete through dimension two")
    print("PASS: the oriented permutohedron three-cell kills topological H2")
    print("PASS: perfect constant coefficients duplicate the top-cell relation")
    print("PASS: genuine three-cells preserve cellular H0 and H1")
    print("PASS: the relative path cellular chains match the face-poset bar")
    print("PASS: the filled braid cellular chains match the face-poset bar")
    print("PASS: subdivision comparison is natural for perfect coefficients")
    print("PASS: all six degree-42 Hessian vertex tangent images have rank nine")
    print("PASS: degree-42 move intersections have ranks 5, 6, and 8")
    print("PASS: the common tangent is Dickson rank two plus one excess line")
    print("PASS: all three degree-42 sectors have conormal ranks 5 < 6 < 7")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
