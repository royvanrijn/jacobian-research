#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of N(2A9+D6).

Recover the two A9 and one D6 root components and their combined index-twenty
root basis intrinsically. Exhaust all sixteen A9 component/diagram and D6
diagram chamber maps. A candidate preserves the Niemeier glue exactly when
its lift in the supplied ambient basis is integral. Certify the retained
matrix group, fixed lattices, and conjugacy classes.
"""

from __future__ import annotations

import argparse
import itertools
import json
import runpy
from collections import Counter
from pathlib import Path

from sage.all import Graph, QQ, ZZ, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-2a9-d6-residual-group-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")

rows = COMMON["rows"]
digest = COMMON["digest"]
matrix_key = COMMON["matrix_key"]
row_module_basis = COMMON["row_module_basis"]
primitive_closure_index = COMMON["primitive_closure_index"]
signed_roots = COMMON["signed_roots"]


def component_roots(gram):
    roots = signed_roots(gram)
    unseen = set(range(len(roots)))
    result = []
    while unseen:
        component = {min(unseen)}
        pending = list(component)
        unseen.difference_update(component)
        while pending:
            current = pending.pop()
            adjacent = {
                index
                for index in unseen
                if roots[current] * gram * roots[index] != 0
            }
            component.update(adjacent)
            unseen.difference_update(adjacent)
            pending.extend(sorted(adjacent))
        vectors = [roots[index] for index in sorted(component)]
        result.append((matrix(ZZ, vectors).rank(), len(vectors), vectors))
    assert sorted((rank, count) for rank, count, _roots in result) == [
        (6, 60),
        (9, 90),
        (9, 90),
    ]
    return result


def dynkin_graph(basis, gram):
    rank = basis.nrows()
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(rank))
    graph.add_edges(
        (row, column)
        for row in range(rank)
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def simple_roots(gram, roots, rank):
    for trial in range(1, 100):
        chamber = vector(
            ZZ,
            [
                (index + 1) ** 2 + trial * (index + 1) + trial**2
                for index in range(gram.nrows())
            ],
        )
        values = [root * gram * chamber for root in roots]
        if all(value != 0 for value in values):
            break
    else:
        raise AssertionError("failed to choose a regular chamber")
    positive = [root for root, value in zip(roots, values) if value > 0]
    positive_set = {tuple(map(int, root)) for root in positive}
    basis = matrix(
        ZZ,
        [
            root
            for root in positive
            if not any(
                tuple(map(int, root - other)) in positive_set
                for other in positive
            )
        ],
    )
    cartan = basis * gram * basis.transpose()
    assert basis.nrows() == basis.rank() == rank
    assert all(cartan[index, index] == 2 for index in range(rank))
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(rank)
        for column in range(rank)
        if row != column
    )
    graph = dynkin_graph(basis, gram)
    if rank == 9:
        assert sorted(graph.degree(vertex) for vertex in graph) == [1, 1] + [2] * 7
    else:
        assert rank == 6
        assert sorted(graph.degree(vertex) for vertex in graph) == [
            1, 1, 1, 2, 2, 3
        ]
    assert graph.automorphism_group().order() == 2
    return basis


def graph_variants(left, right, gram):
    rank = left.nrows()
    left_graph = dynkin_graph(left, gram)
    right_graph = dynkin_graph(right, gram)
    isomorphic, initial = left_graph.is_isomorphic(right_graph, certificate=True)
    assert isomorphic
    variants = [
        [automorphism(initial[index]) for index in range(rank)]
        for automorphism in right_graph.automorphism_group()
    ]
    variants = sorted(variants, key=lambda value: value != list(range(rank)))
    assert len(variants) == 2
    if left is right:
        assert variants[0] == list(range(rank))
    return variants


def residue_vector(values):
    return tuple(QQ(value - value.floor()) for value in values)


def add_residues(left, right):
    return residue_vector([a + b for a, b in zip(left, right)])


def quotient_cosets(root_basis_inverse):
    generators = []
    for ambient_index in range(24):
        ambient = vector(ZZ, [int(index == ambient_index) for index in range(24)])
        generators.append(residue_vector(ambient * root_basis_inverse))
    cosets = {(QQ(0),) * 24}
    for generator in generators:
        while True:
            enlarged = cosets | {add_residues(value, generator) for value in cosets}
            if enlarged == cosets:
                break
            cosets = enlarged
    assert len(cosets) == 20
    return generators, cosets


def matrix_order(action):
    identity = identity_matrix(ZZ, 24)
    power = identity
    for order in range(1, 17):
        power *= action
        if power == identity:
            return order
    raise AssertionError("unexpected residual action order")


def candidate_actions(a9_bases, d6_basis, root_basis_inverse, gram):
    a9_variants = {
        (source, target): graph_variants(a9_bases[source], a9_bases[target], gram)
        for source in range(2)
        for target in range(2)
    }
    d6_variants = graph_variants(d6_basis, d6_basis, gram)
    candidates = []
    elements = []
    seen = set()
    for permutation in itertools.permutations(range(2)):
        for a9_selections in itertools.product(range(2), repeat=2):
            for d6_selection in range(2):
                target_basis = matrix(
                    QQ,
                    [
                        a9_bases[permutation[source]][
                            a9_variants[source, permutation[source]][
                                a9_selections[source]
                            ][row]
                        ]
                        for source in range(2)
                        for row in range(9)
                    ]
                    + [d6_basis[d6_variants[d6_selection][row]] for row in range(6)],
                )
                rational_action = root_basis_inverse * target_basis
                integral = all(
                    entry.denominator() == 1 for entry in rational_action.list()
                )
                candidates.append(
                    {
                        "a9_component_permutation_zero_based": list(permutation),
                        "a9_component_diagram_variant_indices": list(a9_selections),
                        "d6_diagram_variant_index": d6_selection,
                        "integral_ambient_lift": integral,
                    }
                )
                if not integral:
                    continue
                action = matrix(ZZ, rational_action)
                assert abs(action.det()) == 1
                assert action * gram * action.transpose() == gram
                key = matrix_key(action)
                assert key not in seen
                seen.add(key)
                elements.append(
                    {
                        "a9_component_permutation_zero_based": tuple(permutation),
                        "a9_component_diagram_variant_indices": tuple(a9_selections),
                        "d6_diagram_variant_index": d6_selection,
                        "matrix": action,
                        "order": matrix_order(action),
                        "fixed_rank": int(
                            24 - (action - identity_matrix(ZZ, 24)).rank()
                        ),
                    }
                )
    assert len(candidates) == 16
    assert len(elements) == 4
    return candidates, elements


def certify_closure_and_classes(elements, gram):
    by_key = {matrix_key(row["matrix"]): index for index, row in enumerate(elements)}
    identity = identity_matrix(ZZ, 24)
    identity_index = by_key[matrix_key(identity)]
    multiplication = {}
    inverses = {}
    for left_index, left in enumerate(elements):
        inverse_key = matrix_key(left["matrix"].inverse().change_ring(ZZ))
        assert inverse_key in by_key
        inverses[left_index] = by_key[inverse_key]
        for right_index, right in enumerate(elements):
            product_key = matrix_key(left["matrix"] * right["matrix"])
            assert product_key in by_key
            multiplication[left_index, right_index] = by_key[product_key]
    generator_index = next(
        index for index, row in enumerate(elements) if row["order"] == 4
    )
    generated = set()
    current = identity_index
    for _power in range(4):
        generated.add(current)
        current = multiplication[current, generator_index]
    assert generated == set(range(4))

    classes = []
    for element_index, element in enumerate(elements):
        action = element["matrix"]
        fixed = row_module_basis(
            (action - identity)
            .transpose()
            .right_kernel_matrix()
            .change_ring(ZZ)
        )
        assert fixed.nrows() == element["fixed_rank"]
        if fixed.nrows():
            assert primitive_closure_index(fixed) == 1
            fixed_determinant = int((fixed * gram * fixed.transpose()).det())
        else:
            fixed_determinant = 1
        classes.append(
            {
                "class_size": 1,
                "action_order": element["order"],
                "fixed_rank": element["fixed_rank"],
                "fixed_determinant": fixed_determinant,
                "representative_a9_component_permutation_zero_based": list(
                    element["a9_component_permutation_zero_based"]
                ),
                "representative_a9_component_diagram_variant_indices": list(
                    element["a9_component_diagram_variant_indices"]
                ),
                "representative_d6_diagram_variant_index": element[
                    "d6_diagram_variant_index"
                ],
                "representative_matrix": rows(action),
                "element_indices_zero_based": [element_index],
            }
        )
    classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
            row["representative_a9_component_permutation_zero_based"],
            row["representative_a9_component_diagram_variant_indices"],
            row["representative_d6_diagram_variant_index"],
        )
    )
    for index, row in enumerate(classes, start=1):
        row["class_id"] = f"2A9_D6-C{index:02d}"
    return generator_index, classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "2A9_D6"
    )
    gram = matrix(ZZ, ambient_row["gram"])
    components = component_roots(gram)
    a9_roots = [roots for rank, _count, roots in components if rank == 9]
    d6_roots = next(roots for rank, _count, roots in components if rank == 6)
    a9_bases = [simple_roots(gram, roots, 9) for roots in a9_roots]
    d6_basis = simple_roots(gram, d6_roots, 6)
    root_basis = matrix(
        ZZ,
        [root for basis in a9_bases for root in basis.rows()]
        + list(d6_basis.rows()),
    )
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 20
    root_basis_inverse = root_basis.inverse()
    quotient_generators, quotient = quotient_cosets(root_basis_inverse)
    candidates, elements = candidate_actions(
        a9_bases, d6_basis, root_basis_inverse, gram
    )
    generator_index, classes = certify_closure_and_classes(elements, gram)
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    assert order_distribution == Counter({1: 1, 2: 1, 4: 2})
    assert fixed_rank_distribution == Counter({24: 1, 16: 1, 10: 2})
    assert len(classes) == 4
    assert [row["fixed_determinant"] for row in classes] == [1, 256, 1024, 1024]
    return {
        "schema": "elkies-k3.2a9-d6-residual-group.v1",
        "status": "PASS_EXACT_2A9_D6_GLUE_AND_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "intrinsic recovery of the 2A9+D6 root components and order-20 "
                "root quotient; all sixteen chamber candidates; the complete "
                "integral residual lift group, cyclic generation, fixed lattices, "
                "and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "2A9_D6",
        "root_lattice": {
            "a9_simple_root_bases_in_ambient": [rows(basis) for basis in a9_bases],
            "d6_simple_root_basis_in_ambient": rows(d6_basis),
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "root_lattice_quotient": {
            "order": len(quotient),
            "ambient_basis_generators_mod_root_basis": [
                [str(value) for value in generator]
                for generator in quotient_generators
            ],
            "cosets_mod_root_basis": [
                [str(value) for value in coset] for coset in sorted(quotient)
            ],
        },
        "candidate_chamber_maps": candidates,
        "residual_group": {
            "order": len(elements),
            "a9_component_permutation_image_order": len(
                {row["a9_component_permutation_zero_based"] for row in elements}
            ),
            "candidate_maps_tested": len(candidates),
            "generator_count": 1,
            "generator_element_indices_zero_based": [generator_index],
            "order_distribution": {
                str(key): value for key, value in sorted(order_distribution.items())
            },
            "fixed_rank_distribution": {
                str(key): value
                for key, value in sorted(fixed_rank_distribution.items())
            },
            "conjugacy_classes": classes,
            "elements": [
                {
                    "a9_component_permutation_zero_based": list(
                        row["a9_component_permutation_zero_based"]
                    ),
                    "a9_component_diagram_variant_indices": list(
                        row["a9_component_diagram_variant_indices"]
                    ),
                    "d6_diagram_variant_index": row["d6_diagram_variant_index"],
                    "action_order": row["order"],
                    "fixed_rank": row["fixed_rank"],
                    "matrix": rows(row["matrix"]),
                }
                for row in elements
            ],
        },
        "literature_cross_check": {
            "citation": (
                "Chenevier, The Characteristic Masses of Niemeier Lattices "
                "(2020), 2A9+D6 case"
            ),
            "arxiv": "2002.03707",
            "expected_structure": "cyclic residual group of order four",
            "role": (
                "independent naming/order cross-check only; the artifact derives "
                "the group from the supplied lattice"
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build(json.loads(arguments.catalog.read_text()))
    payload["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(arguments.catalog),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_2a9_d6_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("2A9+D6 residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "2A9_D6_GROUP|quotient={}|group={}|classes={}|status=PASS_EXACT".format(
            payload["root_lattice_quotient"]["order"],
            payload["residual_group"]["order"],
            len(payload["residual_group"]["conjugacy_classes"]),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
