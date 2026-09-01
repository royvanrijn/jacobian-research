#!/usr/bin/env sage -python
"""Calibrate basis-free exact matching on the active R17 consensus core.

The input core is supervised, but every rebasings test hides its integral
coordinates behind a nontrivial GL(17,Z) map.  Sage proposes an incidence-
hypergraph bijection; ``latent_lattice`` independently lifts signs and accepts
only exact unimodular ray maps.  Unequal-cloud common-subgraph matching is not
claimed here.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC))

from sage.all import Graph, Matrix, ZZ  # noqa: E402
from sage.version import version as sage_version  # noqa: E402

from latent_lattice import (  # noqa: E402
    build_relation_complex,
    canonical_unoriented,
    lift_relation_vertex_bijection,
    lift_relation_vertex_injection,
    rational_rank,
    row_embedding_smith_invariant_factors,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
CONSENSUS = ARTIFACTS / "latent_lattice_relation_consensus_v1.json"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_hypergraph_matcher_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def incidence_graph(complex_):
    graph = Graph(loops=True, multiedges=False, sparse=True)
    for index in range(len(complex_.vertices)):
        graph.add_edge(("vertex", index), ("vertex", index), "vertex")
    for edge_index, edge in enumerate(complex_.ternary_relations):
        node = ("edge", edge_index)
        graph.add_edge(node, node, "edge")
        for vertex in edge:
            graph.add_edge(("vertex", vertex), node, "incidence")
    return graph


def transform(vectors, matrix):
    return tuple(
        tuple(
            sum(int(vector[row]) * int(matrix[row, column]) for row in range(17))
            for column in range(17)
        )
        for vector in vectors
    )


def transformations():
    identity = Matrix.identity(ZZ, 17)
    shear = Matrix(identity)
    shear[0, 1] = 1
    permutation = Matrix(ZZ, 17, 17)
    for index in range(17):
        permutation[index, 16 - index] = 1
    signs = Matrix.diagonal(ZZ, [(-1 if index % 3 == 0 else 1) for index in range(17)])
    mixed = shear * permutation * signs
    second = Matrix(identity)
    second[4, 9] = -2
    second[11, 2] = 1
    return {
        "identity": identity,
        "elementary_shear": shear,
        "reverse_permutation": permutation,
        "diagonal_signs": signs,
        "mixed_product": mixed,
        "two_commuting_shears": second,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = json.loads(CONSENSUS.read_text())
    truth_document = json.loads(TRUTH.read_text())
    support_three = next(
        record
        for record in document["coefficient_support_cores"]
        if record["minimum_fibre_support"] == 3
    )["relation_complex"]
    vertices = tuple(
        tuple(map(int, vector))
        for vector, degree in zip(
            support_three["vertices"], support_three["additive_degrees"]
        )
        if int(degree) > 0
    )
    source = build_relation_complex(vertices)
    if len(source.vertices) != 134 or rational_rank(source.vertices) != 17:
        raise ArithmeticError("active R17 consensus core changed")
    source_graph = incidence_graph(source)
    vector_nodes = [("vertex", index) for index in range(len(source.vertices))]
    edge_nodes = [
        ("edge", index) for index in range(len(source.ternary_relations))
    ]
    automorphisms = tuple(
        source_graph.automorphism_group(
            partition=[vector_nodes, edge_nodes], edge_labels=True
        )
    )
    trials = []
    for label, matrix in transformations().items():
        if abs(matrix.det()) != 1:
            raise ArithmeticError("calibration transformation is not unimodular")
        target = build_relation_complex(transform(source.vertices, matrix))
        isomorphic, certificate = source_graph.is_isomorphic(
            incidence_graph(target), certificate=True, edge_labels=True
        )
        if not isomorphic:
            raise ArithmeticError("isomorphic rebasings control was rejected")
        lifts = ()
        automorphisms_tested = 0
        for automorphism in automorphisms:
            automorphisms_tested += 1
            vertex_map = tuple(
                int(certificate[automorphism(("vertex", index))][1])
                for index in range(len(source.vertices))
            )
            try:
                lifts = lift_relation_vertex_bijection(source, target, vertex_map)
            except ValueError:
                lifts = ()
            if lifts:
                break
        if not lifts:
            raise ArithmeticError(
                f"hypergraph certificate for {label} had no exact GL lift"
            )
        exact_replay = all(
            abs(Matrix(ZZ, lift).det()) == 1 for lift in lifts
        )
        if not exact_replay:
            raise ArithmeticError("non-unimodular map escaped the exact lift")
        trials.append(
            {
                "label": label,
                "input_matrix_rows": [list(map(int, row)) for row in matrix.rows()],
                "incidence_isomorphic": isomorphic,
                "source_hypergraph_automorphisms_tested": automorphisms_tested,
                "exact_gl_lift_count": len(lifts),
                "exact_gl_lift_matrix_rows": [
                    [list(row) for row in lift] for lift in lifts
                ],
                "input_map_recovered_up_to_global_sign": (
                    tuple(tuple(map(int, row)) for row in matrix.rows()) in lifts
                    or tuple(tuple(-int(value) for value in row) for row in matrix.rows())
                    in lifts
                ),
            }
        )
    if not all(record["input_map_recovered_up_to_global_sign"] for record in trials):
        raise ArithmeticError("a declared rebasings map was not recovered")

    # A changed incidence count is a cheap negative control for accidental
    # color-blind graph matching.
    corrupted = incidence_graph(source)
    edge_node = next(vertex for vertex in corrupted if vertex[0] == "edge")
    neighbour = next(vertex for vertex in corrupted.neighbors(edge_node) if vertex != edge_node)
    corrupted.delete_edge(edge_node, neighbour, "incidence")
    negative_isomorphic = source_graph.is_isomorphic(corrupted, edge_labels=True)
    if negative_isomorphic:
        raise ArithmeticError("corrupted incidence graph passed isomorphism")

    # Supervised rectangular validator controls.  The held-out truth map is
    # used only to propose the ray injection; exact relations, integral replay,
    # and Smith factors independently certify the resulting 17 x r matrix.
    rectangular_trials = []
    truth_by_label = {
        record["label"]: record for record in truth_document["positive_controls"]
    }
    control_by_label = {
        record["label"]: record for record in document["controls"]
    }
    for held in document["leave_one_fibre_out"]:
        label = held["held_out_label"]
        active_full = build_relation_complex(held["held_out_visible_vectors"])
        active_vectors = tuple(
            vector
            for vector, degree in zip(
                active_full.vertices, active_full.additive_degrees
            )
            if degree > 0
        )
        active = build_relation_complex(active_vectors)
        ambient = build_relation_complex(
            tuple(
                tuple(map(int, row.split()))
                for row in control_by_label[label][
                    "ambient_short_vector_coordinate_rows"
                ]
            )
        )
        matrix = tuple(
            tuple(map(int, row))
            for row in truth_by_label[label]["embedding_matrix_columns"]
        )
        ambient_index = {
            vector: index for index, vector in enumerate(ambient.vertices)
        }
        vertex_map = tuple(
            ambient_index[
                canonical_unoriented(
                    tuple(
                        sum(
                            int(vector[row]) * matrix[row][column]
                            for row in range(17)
                        )
                        for column in range(len(matrix[0]))
                    )
                )
            ]
            for vector in active.vertices
        )
        lifts = lift_relation_vertex_injection(active, ambient, vertex_map)
        if matrix not in lifts and tuple(
            tuple(-value for value in row) for row in matrix
        ) not in lifts:
            raise ArithmeticError("supervised rectangular truth map was not lifted")
        smith = row_embedding_smith_invariant_factors(lifts[0])
        if smith != (1,) * 17:
            raise ArithmeticError("R17 rectangular control is not primitive")

        training_vectors = tuple(
            tuple(map(int, vector))
            for vector in held["training_two_of_three_vectors"]
        )
        visible_training = []
        for vector in training_vectors:
            image = canonical_unoriented(
                tuple(
                    sum(
                        int(vector[row]) * matrix[row][column]
                        for row in range(17)
                    )
                    for column in range(len(matrix[0]))
                )
            )
            if image in ambient_index:
                visible_training.append(vector)
        rectangular_trials.append(
            {
                "label": label,
                "ambient_rank": len(matrix[0]),
                "active_supervised_ray_count": len(active.vertices),
                "active_supervised_rank": rational_rank(active.vertices),
                "active_supervised_relation_count": len(active.ternary_relations),
                "exact_rectangular_lift_count": len(lifts),
                "smith_invariant_factors": list(smith),
                "primitive": True,
                "training_core_visible_in_full_cloud_ray_count": len(
                    visible_training
                ),
                "training_core_visible_in_full_cloud_rank": rational_rank(
                    tuple(visible_training)
                ),
            }
        )

    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-hypergraph-matcher.v1",
        "status": "PASS_BASIS_FREE_REBASING_CONTROLS",
        "scope": "Phase-0 R17 consensus core only; no wgxli target is loaded",
        "active_core": {
            "ray_count": len(source.vertices),
            "rank": rational_rank(source.vertices),
            "ternary_relation_count": len(source.ternary_relations),
            "incidence_automorphism_group_order": len(automorphisms),
            "canonical_digest": source.canonical_digest,
        },
        "trials": trials,
        "supervised_rectangular_validator_trials": rectangular_trials,
        "corrupted_incidence_negative_control_isomorphic": negative_isomorphic,
        "gate_decision": (
            "CLOSED. Exact rebasings, rectangular lifts, and primitive Smith replay "
            "pass when a ray injection is supplied; blind maximum-common-subgraph "
            "recovery between unequal fibre clouds still fails its bounded control."
        ),
        "proof_boundary": (
            "Incidence isomorphism is exact. Every proposed ray bijection is lifted "
            "through exact signed ternary constraints, required to have determinant "
            "plus or minus one, and replayed on every integer ray. The active core "
            "itself comes from supervised control alignments."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CONSENSUS,
                TRUTH,
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {
            "python": platform.python_version(),
            "sage": sage_version,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice hypergraph-matcher artifact is stale")
        print(
            f"LATENTMATCH|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTMATCH|status={payload['status']}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
