#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of the Niemeier lattice N(3D8).

Recover the three D8 root components and the order-eight glue code in
((Z/2)^2)^3 intrinsically. Exhaust all 2^3 * 3! diagram/permutation maps,
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
    ROOT / "artifacts/generated-results/elkies-k3-3d8-glue-residual-group-v1.json"
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
    code = {(0,) * 3}
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
        assert matrix(ZZ, vectors).rank() == 8
        assert len(vectors) == 112
        result.append(vectors)
    assert len(result) == 3
    return result


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(8))
    graph.add_edges(
        (row, column)
        for row in range(8)
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
    assert basis.nrows() == basis.rank() == 8
    assert all(cartan[index, index] == 2 for index in range(8))
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(8)
        for column in range(8)
        if row != column
    )
    graph = dynkin_graph(basis, gram)
    assert sorted(graph.degree(vertex) for vertex in graph) == [
        1, 1, 1, 2, 2, 2, 2, 3
    ]
    assert graph.automorphism_group().order() == 2
    return basis


def graph_variants(left, right, gram):
    left_graph = dynkin_graph(left, gram)
    right_graph = dynkin_graph(right, gram)
    isomorphic, initial = left_graph.is_isomorphic(right_graph, certificate=True)
    assert isomorphic
    variants = [
        [automorphism(initial[index]) for index in range(8)]
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
        reps = [vector(QQ, [0] * 8), first, second, first + second]
        mapping = {residue_key(value): index for index, value in enumerate(reps)}
        assert len(mapping) == 4
        label_maps.append(mapping)
        representatives.append(reps)
    return label_maps, representatives


def recover_glue_code(bases, gram):
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 8
    root_basis_inverse = root_basis.inverse()
    label_maps, representatives = component_label_data(bases, gram)
    ambient_generator_words = []
    for ambient_index in range(24):
        ambient_vector = vector(
            ZZ, [int(index == ambient_index) for index in range(24)]
        )
        root_coordinates = ambient_vector * root_basis_inverse
        word = []
        for component_index in range(3):
            coordinates = root_coordinates[
                8 * component_index : 8 * component_index + 8
            ]
            word.append(label_maps[component_index][residue_key(coordinates)])
        ambient_generator_words.append(tuple(word))
    code = generated_code(ambient_generator_words)
    generators = binary_code_basis(ambient_generator_words)
    assert len(code) == 8
    assert len(generators) == 3
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
    for source in range(3):
        for target in range(3):
            pair = []
            for selection in graph_variants(bases[source], bases[target], gram):
                local_map = []
                for source_coordinates in representatives[source]:
                    target_coordinates = vector(QQ, [0] * 8)
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
    output = [0] * 3
    for source, target in enumerate(permutation):
        local_map = variants[source, target][selections[source]]["local_map"]
        output[target] = local_map[word[source]]
    return tuple(output)


def matrix_order(action):
    identity = identity_matrix(ZZ, 24)
    power = identity
    for order in range(1, 49):
        power *= action
        if power == identity:
            return order
    raise AssertionError("unexpected residual action order")


def lift_stabilizer(code, bases, gram, root_basis_inverse, variants):
    elements = []
    seen = set()
    for permutation in itertools.permutations(range(3)):
        for selections in itertools.product(range(2), repeat=3):
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
                    for source in range(3)
                    for row in range(8)
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
                        for source in range(3)
                    ),
                    "matrix": action,
                    "order": matrix_order(action),
                    "fixed_rank": int(
                        24 - (action - identity_matrix(ZZ, 24)).rank()
                    ),
                }
            )
    assert len(elements) == 6
    return elements


def certify_closure_generators_and_classes(elements, gram):
    by_key = {matrix_key(row["matrix"]): index for index, row in enumerate(elements)}
    identity = identity_matrix(ZZ, 24)
    identity_key = matrix_key(identity)
    assert identity_key in by_key
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

    generated = {by_key[identity_key]}
    generators = []

    def close():
        pending = list(generated)
        while pending:
            current = pending.pop()
            for generator in generators:
                for product in (
                    multiplication[current, generator],
                    multiplication[generator, current],
                ):
                    if product not in generated:
                        generated.add(product)
                        pending.append(product)

    for index in range(len(elements)):
        if index in generated:
            continue
        generators.append(index)
        close()
    assert len(generated) == len(elements)

    remaining = set(range(len(elements)))
    classes = []
    while remaining:
        representative_index = min(remaining)
        class_indices = {
            multiplication[
                multiplication[inverses[conjugator], representative_index],
                conjugator,
            ]
            for conjugator in range(len(elements))
        }
        assert class_indices <= remaining
        remaining.difference_update(class_indices)
        orders = {elements[index]["order"] for index in class_indices}
        fixed_ranks = {elements[index]["fixed_rank"] for index in class_indices}
        assert len(orders) == len(fixed_ranks) == 1
        representative = elements[representative_index]
        action = representative["matrix"]
        fixed = row_module_basis(
            (action - identity)
            .transpose()
            .right_kernel_matrix()
            .change_ring(ZZ)
        )
        fixed_rank = next(iter(fixed_ranks))
        assert fixed.nrows() == fixed_rank
        assert primitive_closure_index(fixed) == 1
        fixed_gram = fixed * gram * fixed.transpose()
        classes.append(
            {
                "class_size": len(class_indices),
                "action_order": next(iter(orders)),
                "fixed_rank": fixed_rank,
                "fixed_determinant": int(fixed_gram.det()),
                "representative_component_permutation_zero_based": list(
                    representative["component_permutation_zero_based"]
                ),
                "representative_component_diagram_variant_indices": list(
                    representative["component_diagram_variant_indices"]
                ),
                "representative_matrix": rows(action),
                "element_indices_zero_based": sorted(class_indices),
            }
        )
    classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
            row["representative_component_permutation_zero_based"],
            row["representative_component_diagram_variant_indices"],
        )
    )
    for index, row in enumerate(classes, start=1):
        row["class_id"] = f"3D8-C{index:02d}"
    return generators, classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row for row in catalog["rooted_niemeier_lattices"] if row["label"] == "3D8"
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
    generators, classes = certify_closure_generators_and_classes(elements, gram)
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    component_image_order = len(
        {row["component_permutation_zero_based"] for row in elements}
    )
    kernel_order = sum(
        row["component_permutation_zero_based"] == tuple(range(3))
        for row in elements
    )
    assert order_distribution == Counter({1: 1, 2: 3, 3: 2})
    assert fixed_rank_distribution == Counter({24: 1, 16: 3, 8: 2})
    assert component_image_order == 6
    assert kernel_order == 1
    assert len(generators) == 2
    assert len(classes) == 3
    return {
        "schema": "elkies-k3.3d8-glue-residual-group.v1",
        "status": "PASS_EXACT_3D8_GLUE_AND_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "the recovered order-eight N(3D8)/D8^3 glue code; its complete "
                "diagram/permutation stabilizer among all 48 candidates; integral "
                "lifts, matrix closure, generators, and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "3D8",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "glue_code": {
            "alphabet": "(Z/2)^2 with labels 0,1,2,3 and XOR addition",
            "length": 3,
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
            "diagram_permutation_candidates_tested": 2**3 * 6,
            "generator_count": len(generators),
            "generator_element_indices_zero_based": generators,
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
                "(2020), 3D8 case"
            ),
            "arxiv": "2002.03707",
            "expected_structure": (
                "natural S3 component-permutation subgroup of the "
                "hyperoctahedral diagram group H3"
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
        "elkies-k3/scripts/probe_3d8_glue_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("3D8 residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "3D8GROUP|glue={}|group={}|image={}|kernel={}|classes={}|status=PASS_EXACT".format(
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
