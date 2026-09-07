#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of the Niemeier lattice N(2D12).

Recover the two D12 root components and the order-four glue code in
((Z/2)^2)^2 intrinsically. Exhaust all 2^2 * 2! diagram/permutation maps,
retain exactly the glue stabilizer, lift it to the supplied integral Gram
model, and certify matrix closure and conjugacy classes.
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
    ROOT / "artifacts/generated-results/elkies-k3-2d12-glue-residual-group-v1.json"
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


def residue_key(values):
    return tuple(QQ(value - value.floor()) for value in values)


def add_words(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def generated_code(generators):
    code = {(0,) * 2}
    for generator in generators:
        code |= {add_words(word, generator) for word in list(code)}
    return code


def binary_code_basis(words):
    pivots = {}
    selected = []
    for word in sorted(set(words)):
        encoded = sum(int(value) << (2 * index) for index, value in enumerate(word))
        reduced = encoded
        for pivot in sorted(pivots, reverse=True):
            if (reduced >> pivot) & 1:
                reduced ^= pivots[pivot]
        if reduced == 0:
            continue
        pivot = reduced.bit_length() - 1
        pivots[pivot] = reduced
        selected.append(word)
    return selected


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
        assert matrix(ZZ, vectors).rank() == 12
        assert len(vectors) == 264
        result.append(vectors)
    assert len(result) == 2
    return result


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(12))
    graph.add_edges(
        (row, column)
        for row in range(12)
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def simple_roots(gram, roots):
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
    assert basis.nrows() == basis.rank() == 12
    assert all(cartan[index, index] == 2 for index in range(12))
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(12)
        for column in range(12)
        if row != column
    )
    graph = dynkin_graph(basis, gram)
    assert sorted(graph.degree(vertex) for vertex in graph) == (
        [1, 1, 1] + [2] * 8 + [3]
    )
    assert graph.automorphism_group().order() == 2
    return basis


def graph_variants(left, right, gram):
    left_graph = dynkin_graph(left, gram)
    right_graph = dynkin_graph(right, gram)
    isomorphic, initial = left_graph.is_isomorphic(right_graph, certificate=True)
    assert isomorphic
    variants = [
        [automorphism(initial[index]) for index in range(12)]
        for automorphism in right_graph.automorphism_group()
    ]
    assert len(variants) == 2
    return variants


def component_label_data(bases, gram):
    label_maps = []
    representatives = []
    for basis in bases:
        cartan = basis * gram * basis.transpose()
        graph = dynkin_graph(basis, gram)
        leaves = sorted(vertex for vertex in graph if graph.degree(vertex) == 1)
        first = cartan.inverse().row(leaves[0])
        second = cartan.inverse().row(leaves[1])
        reps = [vector(QQ, [0] * 12), first, second, first + second]
        mapping = {residue_key(value): index for index, value in enumerate(reps)}
        assert len(mapping) == 4
        label_maps.append(mapping)
        representatives.append(reps)
    return label_maps, representatives


def recover_glue_code(bases, gram):
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 4
    root_basis_inverse = root_basis.inverse()
    label_maps, representatives = component_label_data(bases, gram)
    ambient_generator_words = []
    for ambient_index in range(24):
        ambient_vector = vector(
            ZZ, [int(index == ambient_index) for index in range(24)]
        )
        root_coordinates = ambient_vector * root_basis_inverse
        word = []
        for component_index in range(2):
            coordinates = root_coordinates[
                12 * component_index : 12 * component_index + 12
            ]
            word.append(label_maps[component_index][residue_key(coordinates)])
        ambient_generator_words.append(tuple(word))
    code = generated_code(ambient_generator_words)
    generators = binary_code_basis(ambient_generator_words)
    assert len(code) == 4
    assert len(generators) == 2
    assert generated_code(generators) == code
    return (
        root_basis,
        root_basis_inverse,
        label_maps,
        representatives,
        code,
        generators,
        ambient_generator_words,
    )


def variant_data(bases, gram, label_maps, representatives):
    result = {}
    for source in range(2):
        for target in range(2):
            pair = []
            for selection in graph_variants(bases[source], bases[target], gram):
                local_map = []
                for source_coordinates in representatives[source]:
                    target_coordinates = vector(QQ, [0] * 12)
                    for row, image in enumerate(selection):
                        target_coordinates[image] = source_coordinates[row]
                    local_map.append(
                        label_maps[target][residue_key(target_coordinates)]
                    )
                pair.append(
                    {
                        "selection": selection,
                        "local_map": tuple(local_map),
                    }
                )
            assert len({row["local_map"] for row in pair}) == 2
            result[source, target] = pair
    return result


def transform_word(word, permutation, selections, variants):
    output = [0] * 2
    for source, target in enumerate(permutation):
        local_map = variants[source, target][selections[source]]["local_map"]
        output[target] = local_map[word[source]]
    return tuple(output)


def matrix_order(action):
    identity = identity_matrix(ZZ, 24)
    power = identity
    for order in range(1, 9):
        power *= action
        if power == identity:
            return order
    raise AssertionError("unexpected residual action order")


def lift_stabilizer(code, bases, gram, root_basis_inverse, variants):
    elements = []
    seen = set()
    for permutation in itertools.permutations(range(2)):
        for selections in itertools.product(range(2), repeat=2):
            if {
                transform_word(word, permutation, selections, variants)
                for word in code
            } != code:
                continue
            target_basis = matrix(
                QQ,
                [
                    bases[permutation[source]][
                        variants[source, permutation[source]][selections[source]][
                            "selection"
                        ][row]
                    ]
                    for source in range(2)
                    for row in range(12)
                ],
            )
            action = root_basis_inverse * target_basis
            assert all(entry.denominator() == 1 for entry in action.list())
            action = matrix(ZZ, action)
            assert abs(action.det()) == 1
            assert action * gram * action.transpose() == gram
            key = matrix_key(action)
            assert key not in seen
            seen.add(key)
            elements.append(
                {
                    "component_permutation_zero_based": tuple(permutation),
                    "component_diagram_variant_indices": tuple(selections),
                    "component_discriminant_maps": tuple(
                        variants[source, permutation[source]][selections[source]][
                            "local_map"
                        ]
                        for source in range(2)
                    ),
                    "matrix": action,
                    "order": matrix_order(action),
                    "fixed_rank": int(
                        24 - (action - identity_matrix(ZZ, 24)).rank()
                    ),
                }
            )
    assert len(elements) == 2
    return elements


def certify_closure_and_classes(elements, gram):
    by_key = {matrix_key(row["matrix"]): index for index, row in enumerate(elements)}
    identity = identity_matrix(ZZ, 24)
    assert matrix_key(identity) in by_key
    for left in elements:
        assert matrix_key(left["matrix"].inverse().change_ring(ZZ)) in by_key
        for right in elements:
            assert matrix_key(left["matrix"] * right["matrix"]) in by_key
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
        assert primitive_closure_index(fixed) == 1
        fixed_gram = fixed * gram * fixed.transpose()
        classes.append(
            {
                "class_size": 1,
                "action_order": element["order"],
                "fixed_rank": element["fixed_rank"],
                "fixed_determinant": int(fixed_gram.det()),
                "representative_component_permutation_zero_based": list(
                    element["component_permutation_zero_based"]
                ),
                "representative_component_diagram_variant_indices": list(
                    element["component_diagram_variant_indices"]
                ),
                "representative_matrix": rows(action),
                "element_indices_zero_based": [element_index],
            }
        )
    classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
        )
    )
    for index, row in enumerate(classes, start=1):
        row["class_id"] = f"2D12-C{index:02d}"
    return classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "2D12"
    )
    gram = matrix(ZZ, ambient_row["gram"])
    bases = [simple_roots(gram, roots) for roots in component_roots(gram)]
    (
        root_basis,
        root_basis_inverse,
        label_maps,
        representatives,
        code,
        code_generators,
        ambient_generator_words,
    ) = recover_glue_code(bases, gram)
    variants = variant_data(bases, gram, label_maps, representatives)
    elements = lift_stabilizer(code, bases, gram, root_basis_inverse, variants)
    classes = certify_closure_and_classes(elements, gram)
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    component_image_order = len(
        {row["component_permutation_zero_based"] for row in elements}
    )
    kernel_order = sum(
        row["component_permutation_zero_based"] == tuple(range(2))
        for row in elements
    )
    assert order_distribution == Counter({1: 1, 2: 1})
    assert fixed_rank_distribution == Counter({24: 1, 12: 1})
    assert component_image_order == 2
    assert kernel_order == 1
    assert len(classes) == 2
    return {
        "schema": "elkies-k3.2d12-glue-residual-group.v1",
        "status": "PASS_EXACT_2D12_GLUE_AND_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "the recovered order-four N(2D12)/D12^2 glue code; its complete "
                "diagram/permutation stabilizer among all eight candidates; "
                "integral lifts, matrix closure, and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "2D12",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "glue_code": {
            "alphabet": "(Z/2)^2 with labels 0,1,2,3 and XOR addition",
            "length": 2,
            "order": len(code),
            "binary_dimension": len(code_generators),
            "binary_generators": [list(word) for word in code_generators],
            "ambient_basis_generator_words": [
                list(word) for word in ambient_generator_words
            ],
            "codewords": [list(word) for word in sorted(code)],
        },
        "residual_group": {
            "order": len(elements),
            "component_permutation_image_order": component_image_order,
            "component_kernel_order": kernel_order,
            "diagram_permutation_candidates_tested": 2**2 * 2,
            "generator_count": 1,
            "generator_element_indices_zero_based": [
                index for index, row in enumerate(elements) if row["order"] == 2
            ],
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
                    "component_permutation_zero_based": list(
                        row["component_permutation_zero_based"]
                    ),
                    "component_diagram_variant_indices": list(
                        row["component_diagram_variant_indices"]
                    ),
                    "component_discriminant_maps": [
                        list(value) for value in row["component_discriminant_maps"]
                    ],
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
                "(2020), 2D12 case"
            ),
            "arxiv": "2002.03707",
            "expected_structure": (
                "natural S2 component-permutation subgroup of the "
                "hyperoctahedral diagram group H2"
            ),
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
        "elkies-k3/scripts/probe_2d12_glue_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("2D12 residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "2D12GROUP|glue={}|group={}|image={}|kernel={}|classes={}|status=PASS_EXACT".format(
            payload["glue_code"]["order"],
            payload["residual_group"]["order"],
            payload["residual_group"]["component_permutation_image_order"],
            payload["residual_group"]["component_kernel_order"],
            len(payload["residual_group"]["conjugacy_classes"]),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
