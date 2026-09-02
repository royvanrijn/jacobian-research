#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of the Niemeier lattice N(6D4).

Recover the six D4 root components intrinsically and identify N/D4^6 as a
six-dimensional binary code in ((Z/2)^2)^6.  Exhaust the complete triality
wreath product (S3)^6 semidirect S6 on a binary basis of that code, lift every
stabilizer to the supplied integral Niemeier Gram model, and certify group
closure and matrix conjugacy classes.

The recovered code is the hexacode after choosing F4 labels.  That name and
the expected group structure are literature cross-checks only: all exact
group data in the artifact are derived from the supplied lattice.
"""

from __future__ import annotations

import argparse
import itertools
import json
import runpy
from collections import Counter
from pathlib import Path

import numpy as np
from sage.all import Graph, QQ, ZZ, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-6d4-hexacode-residual-group-v1.json"
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


LOCAL_MAPS = [
    (0,) + permutation for permutation in itertools.permutations((1, 2, 3))
]
LOCAL_MAP_INDEX = {value: index for index, value in enumerate(LOCAL_MAPS)}
LOCAL_MAP_ARRAY = np.asarray(LOCAL_MAPS, dtype=np.uint8)


def residue_key(values):
    return tuple(QQ(value - value.floor()) for value in values)


def encode_word(word):
    return sum(int(value) << (2 * index) for index, value in enumerate(word))


def decode_word(encoded, length=6):
    return tuple((encoded >> (2 * index)) & 3 for index in range(length))


def add_words(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def generated_code(generators):
    code = {(0,) * 6}
    for generator in generators:
        code |= {add_words(word, generator) for word in list(code)}
    return code


def binary_code_basis(words):
    pivots = {}
    selected = []
    for word in sorted(set(words)):
        reduced = encode_word(word)
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
        assert matrix(ZZ, vectors).rank() == 4
        assert len(vectors) == 24
        result.append(vectors)
    assert len(result) == 6
    return result


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
    assert basis.nrows() == basis.rank() == 4
    assert all(cartan[index, index] == 2 for index in range(4))
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(4)
        for column in range(4)
        if row != column
    )
    graph = dynkin_graph(basis, gram)
    assert sorted(graph.degree(vertex) for vertex in graph) == [1, 1, 1, 3]
    return basis


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(4))
    graph.add_edges(
        (row, column)
        for row in range(4)
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def graph_variants(left, right, gram):
    left_graph = dynkin_graph(left, gram)
    right_graph = dynkin_graph(right, gram)
    isomorphic, initial = left_graph.is_isomorphic(
        right_graph, certificate=True
    )
    assert isomorphic
    variants = [
        [automorphism(initial[index]) for index in range(4)]
        for automorphism in right_graph.automorphism_group()
    ]
    assert len(variants) == 6
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
        reps = [vector(QQ, [0] * 4), first, second, first + second]
        mapping = {residue_key(value): index for index, value in enumerate(reps)}
        assert len(mapping) == 4
        label_maps.append(mapping)
        representatives.append(reps)
    return label_maps, representatives


def recover_glue_code(bases, gram):
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 64
    root_basis_inverse = root_basis.inverse()
    label_maps, representatives = component_label_data(bases, gram)
    ambient_generator_words = []
    for ambient_index in range(24):
        ambient_vector = vector(
            ZZ, [int(index == ambient_index) for index in range(24)]
        )
        root_coordinates = ambient_vector * root_basis_inverse
        word = []
        for component_index in range(6):
            component_coordinates = root_coordinates[
                4 * component_index : 4 * component_index + 4
            ]
            word.append(
                label_maps[component_index][residue_key(component_coordinates)]
            )
        ambient_generator_words.append(tuple(word))
    code = generated_code(ambient_generator_words)
    assert len(code) == 64
    minimal_generators = binary_code_basis(ambient_generator_words)
    assert len(minimal_generators) == 6
    assert generated_code(minimal_generators) == code
    return (
        root_basis,
        root_basis_inverse,
        label_maps,
        representatives,
        code,
        minimal_generators,
        ambient_generator_words,
    )


def transform_word(word, element):
    permutation, local_maps = element
    output = [0] * 6
    for source, target in enumerate(permutation):
        output[target] = LOCAL_MAPS[local_maps[source]][word[source]]
    return tuple(output)


def triality_permutation_stabilizer(code, generators):
    code_lookup = np.zeros(1 << 12, dtype=bool)
    for word in code:
        code_lookup[encode_word(word)] = True
    assignments = np.indices((6,) * 6, dtype=np.uint8).reshape(6, -1).T
    assert assignments.shape == (46656, 6)
    elements = []
    for permutation in itertools.permutations(range(6)):
        alive = np.ones(len(assignments), dtype=bool)
        for generator in generators:
            encoded = np.zeros(len(assignments), dtype=np.uint16)
            for source, target in enumerate(permutation):
                encoded |= (
                    LOCAL_MAP_ARRAY[assignments[:, source], generator[source]]
                    .astype(np.uint16)
                    << (2 * target)
                )
            alive &= code_lookup[encoded]
            if not np.any(alive):
                break
        for local_maps in assignments[alive]:
            element = tuple(permutation), tuple(map(int, local_maps))
            assert {transform_word(word, element) for word in code} == code
            elements.append(element)
    assert len(elements) == 2160
    return elements


def variant_local_maps(bases, gram, label_maps, representatives):
    result = {}
    for source in range(6):
        for target in range(6):
            by_map = {}
            for selection in graph_variants(bases[source], bases[target], gram):
                local_map = []
                for source_coordinates in representatives[source]:
                    target_coordinates = vector(QQ, [0] * 4)
                    for row, image in enumerate(selection):
                        target_coordinates[image] = source_coordinates[row]
                    local_map.append(
                        label_maps[target][residue_key(target_coordinates)]
                    )
                local_map = tuple(local_map)
                assert local_map in LOCAL_MAP_INDEX
                local_index = LOCAL_MAP_INDEX[local_map]
                assert local_index not in by_map
                by_map[local_index] = selection
            assert set(by_map) == set(range(6))
            result[source, target] = by_map
    return result


def element_key(element):
    return tuple(element[0]), tuple(element[1])


LOCAL_COMPOSE = {}
LOCAL_INVERSE = {}
for left_index, left in enumerate(LOCAL_MAPS):
    for right_index, right in enumerate(LOCAL_MAPS):
        composed = tuple(right[left[value]] for value in range(4))
        LOCAL_COMPOSE[left_index, right_index] = LOCAL_MAP_INDEX[composed]
    inverse = [0] * 4
    for source, target in enumerate(left):
        inverse[target] = source
    LOCAL_INVERSE[left_index] = LOCAL_MAP_INDEX[tuple(inverse)]


def compose_elements(left, right):
    left_permutation, left_maps = left
    right_permutation, right_maps = right
    return (
        tuple(
            right_permutation[left_permutation[index]] for index in range(6)
        ),
        tuple(
            LOCAL_COMPOSE[
                left_maps[index], right_maps[left_permutation[index]]
            ]
            for index in range(6)
        ),
    )


def inverse_element(element):
    permutation, local_maps = element
    inverse_permutation = [0] * 6
    inverse_maps = [0] * 6
    for source, target in enumerate(permutation):
        inverse_permutation[target] = source
        inverse_maps[target] = LOCAL_INVERSE[local_maps[source]]
    return tuple(inverse_permutation), tuple(inverse_maps)


def element_order(element):
    identity = tuple(range(6)), (0,) * 6
    power = identity
    for order in range(1, 121):
        power = compose_elements(power, element)
        if power == identity:
            return order
    raise AssertionError("unexpected triality-permutation order")


def lift_group(
    abstract_elements,
    bases,
    gram,
    root_basis_inverse,
    variants_by_map,
):
    result = []
    seen_matrices = set()
    for element in abstract_elements:
        permutation, local_maps = element
        target_basis = matrix(
            QQ,
            [
                bases[permutation[source]][
                    variants_by_map[source, permutation[source]][
                        local_maps[source]
                    ][row]
                ]
                for source in range(6)
                for row in range(4)
            ],
        )
        action = root_basis_inverse * target_basis
        assert all(entry.denominator() == 1 for entry in action.list())
        action = matrix(ZZ, action)
        assert abs(action.det()) == 1
        assert action * gram * action.transpose() == gram
        key = matrix_key(action)
        assert key not in seen_matrices
        seen_matrices.add(key)
        result.append(
            {
                "element": element,
                "matrix": action,
                "order": element_order(element),
                "fixed_rank": int(
                    24 - (action - identity_matrix(ZZ, 24)).rank()
                ),
            }
        )
    return result


def certify_closure_and_classes(elements, gram):
    by_key = {element_key(row["element"]): row for row in elements}
    identity = tuple(range(6)), (0,) * 6
    generated = {element_key(identity): identity}
    generators = []

    def close():
        pending = list(generated.values())
        while pending:
            current = pending.pop()
            for generator in generators:
                for product in (
                    compose_elements(current, generator),
                    compose_elements(generator, current),
                ):
                    key = element_key(product)
                    assert key in by_key
                    if key not in generated:
                        generated[key] = product
                        pending.append(product)

    for row in elements:
        key = element_key(row["element"])
        if key in generated:
            continue
        generators.append(row["element"])
        close()
        if len(generated) == len(elements):
            break
    assert set(generated) == set(by_key)

    abstract_group = [row["element"] for row in elements]
    remaining = set(by_key)
    conjugacy_classes = []
    while remaining:
        representative_key = min(remaining)
        representative = by_key[representative_key]
        class_keys = {
            element_key(
                compose_elements(
                    compose_elements(
                        inverse_element(conjugator),
                        representative["element"],
                    ),
                    conjugator,
                )
            )
            for conjugator in abstract_group
        }
        assert class_keys <= remaining
        remaining.difference_update(class_keys)
        orders = {by_key[key]["order"] for key in class_keys}
        fixed_ranks = {by_key[key]["fixed_rank"] for key in class_keys}
        assert len(orders) == len(fixed_ranks) == 1
        action = representative["matrix"]
        fixed = row_module_basis(
            (action - identity_matrix(ZZ, 24))
            .transpose()
            .right_kernel_matrix()
            .change_ring(ZZ)
        )
        fixed_rank = next(iter(fixed_ranks))
        assert fixed.nrows() == fixed_rank
        if fixed_rank:
            assert primitive_closure_index(fixed) == 1
            fixed_gram = fixed * gram * fixed.transpose()
            fixed_determinant = int(fixed_gram.det())
        else:
            fixed_determinant = 1
        permutation, local_maps = representative["element"]
        conjugacy_classes.append(
            {
                "class_size": len(class_keys),
                "action_order": next(iter(orders)),
                "fixed_rank": fixed_rank,
                "fixed_determinant": fixed_determinant,
                "representative_triality_permutation": {
                    "component_permutation_zero_based": list(permutation),
                    "component_triality_maps": [
                        list(LOCAL_MAPS[index]) for index in local_maps
                    ],
                },
                "representative_matrix": rows(action),
                "element_keys": [
                    [list(key[0]), list(key[1])] for key in sorted(class_keys)
                ],
            }
        )
    conjugacy_classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
            row["representative_triality_permutation"][
                "component_permutation_zero_based"
            ],
            row["representative_triality_permutation"][
                "component_triality_maps"
            ],
        )
    )
    for index, row in enumerate(conjugacy_classes, start=1):
        row["class_id"] = f"6D4-C{index:02d}"
    return generators, conjugacy_classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "6D4"
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
    abstract_elements = triality_permutation_stabilizer(code, code_generators)
    variants_by_map = variant_local_maps(
        bases, gram, label_maps, representatives
    )
    elements = lift_group(
        abstract_elements,
        bases,
        gram,
        root_basis_inverse,
        variants_by_map,
    )
    generators, conjugacy_classes = certify_closure_and_classes(elements, gram)
    component_image_order = len({row["element"][0] for row in elements})
    kernel_order = sum(row["element"][0] == tuple(range(6)) for row in elements)
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    assert len(code) == 64
    assert len(elements) == 2160
    assert component_image_order == 720
    assert kernel_order == 3
    assert len(generators) == 6
    assert len(conjugacy_classes) == 16
    assert order_distribution == Counter(
        {1: 1, 2: 135, 3: 242, 4: 360, 5: 144, 6: 810, 12: 180, 15: 288}
    )
    assert fixed_rank_distribution == Counter(
        {4: 828, 8: 804, 10: 270, 12: 167, 16: 90, 24: 1}
    )
    return {
        "schema": "elkies-k3.6d4-hexacode-residual-group.v1",
        "status": "PASS_EXACT_6D4_HEXACODE_AND_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "the complete triality-permutation stabilizer of the recovered "
                "N(6D4)/D4^6 glue code, its integral lifts to the supplied "
                "Niemeier Gram model, matrix-group closure, and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "6D4",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "glue_code": {
            "alphabet": "(Z/2)^2 with labels 0,1,2,3 and XOR addition",
            "length": 6,
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
            "triality_permutation_candidates_tested": 720 * 6**6,
            "generator_count": len(generators),
            "generators": [
                {
                    "component_permutation_zero_based": list(value[0]),
                    "component_triality_maps": [
                        list(LOCAL_MAPS[index]) for index in value[1]
                    ],
                }
                for value in generators
            ],
            "order_distribution": {
                str(key): value for key, value in sorted(order_distribution.items())
            },
            "fixed_rank_distribution": {
                str(key): value
                for key, value in sorted(fixed_rank_distribution.items())
            },
            "conjugacy_classes": conjugacy_classes,
            "elements": [
                {
                    "component_permutation_zero_based": list(row["element"][0]),
                    "component_triality_map_indices": list(row["element"][1]),
                    "component_triality_maps": [
                        list(LOCAL_MAPS[index]) for index in row["element"][1]
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
                "(2020), D4^6 case"
            ),
            "arxiv": "2002.03707",
            "expected_structure": (
                "hexacode stabilizer with triality kernel order 3 and "
                "component image S6"
            ),
            "role": (
                "independent naming/order cross-check only; the artifact "
                "derives the group from the supplied lattice"
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
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(
            arguments.catalog
        ),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_6d4_hexacode_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("6D4 hexacode residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "6D4GROUP|glue={}|group={}|image={}|kernel={}|classes={}|status=PASS_EXACT".format(
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
