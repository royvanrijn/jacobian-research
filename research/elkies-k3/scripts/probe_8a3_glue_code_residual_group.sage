#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of the Niemeier lattice N(8A3).

Recover the eight A3 root components intrinsically and identify the glue
N/R as a subgroup of (Z/4)^8.  Exhaust the full signed-coordinate group
2^8 semidirect S8 on the four glue generators, lift every code stabilizer to
the supplied integral Niemeier Gram model, and certify closure and matrix
conjugacy classes.

This is exact residual-group infrastructure.  It does not enumerate
rank-seven auxiliaries or quotient ambient embeddings by the Weyl group.
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
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-glue-code-residual-group-v1.json"
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


def encode_word(word):
    return sum(int(value) << (2 * index) for index, value in enumerate(word))


def decode_word(encoded, length=8):
    return tuple((encoded >> (2 * index)) & 3 for index in range(length))


def add_words(left, right):
    return tuple((a + b) % 4 for a, b in zip(left, right))


def scalar_word(scalar, word):
    return tuple((scalar * value) % 4 for value in word)


def generated_code(generators):
    code = {(0,) * 8}
    for generator in generators:
        code = {
            add_words(word, scalar_word(scalar, generator))
            for word in code
            for scalar in range(4)
        }
    return code


def simple_root_components(gram):
    roots = signed_roots(gram)
    unseen = set(range(len(roots)))
    components = []
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
        component_roots = [roots[index] for index in sorted(component)]
        assert matrix(ZZ, component_roots).rank() == 3
        assert len(component_roots) == 12
        components.append(component_roots)
    assert len(components) == 8

    bases = []
    for component_roots in components:
        for trial in range(1, 100):
            chamber = vector(
                ZZ,
                [
                    (index + 1) ** 2
                    + trial * (index + 1)
                    + trial**2
                    for index in range(24)
                ],
            )
            values = [root * gram * chamber for root in component_roots]
            if all(value != 0 for value in values):
                break
        positive = [
            root
            for root, value in zip(component_roots, values)
            if value > 0
        ]
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
        assert basis.nrows() == basis.rank() == 3
        assert sorted(cartan.list()) == [-1, -1, -1, -1, 0, 0, 2, 2, 2]
        bases.append(basis)
    return bases


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(3))
    graph.add_edges(
        (row, column)
        for row in range(3)
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
        [automorphism(initial[index]) for index in range(3)]
        for automorphism in right_graph.automorphism_group()
    ]
    assert len(variants) == 2
    return variants


def component_label_maps(bases, gram):
    label_maps = []
    generators = []
    for basis in bases:
        cartan = basis * gram * basis.transpose()
        graph = dynkin_graph(basis, gram)
        leaf = min(vertex for vertex in graph if graph.degree(vertex) == 1)
        generator = cartan.inverse().row(leaf)
        mapping = {
            residue_key(scalar * generator): scalar
            for scalar in range(4)
        }
        assert len(mapping) == 4
        label_maps.append(mapping)
        generators.append(generator)
    return label_maps, generators


def recover_glue_code(bases, gram):
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 256
    root_basis_inverse = root_basis.inverse()
    label_maps, component_generators = component_label_maps(bases, gram)
    ambient_generator_words = []
    for ambient_index in range(24):
        ambient_vector = vector(
            ZZ, [int(index == ambient_index) for index in range(24)]
        )
        root_coordinates = ambient_vector * root_basis_inverse
        word = []
        for component_index in range(8):
            component_coordinates = root_coordinates[
                3 * component_index : 3 * component_index + 3
            ]
            word.append(
                label_maps[component_index][
                    residue_key(component_coordinates)
                ]
            )
        ambient_generator_words.append(tuple(word))
    code = generated_code(ambient_generator_words)
    assert len(code) == 256
    unique_ambient_words = sorted(set(ambient_generator_words))
    minimal_generators = None
    for count in range(1, 9):
        for candidates in itertools.combinations(unique_ambient_words, count):
            if generated_code(candidates) == code:
                minimal_generators = list(candidates)
                break
        if minimal_generators is not None:
            break
    assert minimal_generators is not None
    return (
        root_basis,
        root_basis_inverse,
        component_generators,
        code,
        minimal_generators,
        ambient_generator_words,
    )


def permute_word(word, permutation):
    output = [0] * 8
    for source, target in enumerate(permutation):
        output[target] = word[source]
    return tuple(output)


EXPANDED_SIGN_MASKS = [
    sum(((mask >> index) & 1) << (2 * index + 1) for index in range(8))
    for mask in range(256)
]


def odd_coordinate_mask(encoded):
    return sum(
        (((encoded >> (2 * index)) & 1) << index)
        for index in range(8)
    )


def signed_encoded_word(encoded, output_sign_mask):
    return encoded ^ EXPANDED_SIGN_MASKS[
        output_sign_mask & odd_coordinate_mask(encoded)
    ]


def transform_word(word, signed_permutation):
    permutation, signs = signed_permutation
    output = [0] * 8
    for source, target in enumerate(permutation):
        output[target] = (signs[source] * word[source]) % 4
    return tuple(output)


def signed_permutation_stabilizer(code, generators):
    code_encoded = {encode_word(word) for word in code}
    elements = []
    permutations_tested = 0
    signed_candidates_tested = 0
    for permutation in itertools.permutations(range(8)):
        permutations_tested += 1
        valid_output_masks = list(range(256))
        for generator in generators:
            permuted = encode_word(permute_word(generator, permutation))
            valid_output_masks = [
                mask
                for mask in valid_output_masks
                if signed_encoded_word(permuted, mask) in code_encoded
            ]
            if not valid_output_masks:
                break
        signed_candidates_tested += 256
        for output_mask in valid_output_masks:
            signs = tuple(
                -1 if (output_mask >> permutation[source]) & 1 else 1
                for source in range(8)
            )
            signed_permutation = (tuple(permutation), signs)
            assert {
                transform_word(word, signed_permutation) for word in code
            } == code
            elements.append(signed_permutation)
    assert permutations_tested == 40320
    assert signed_candidates_tested == 10321920
    return elements


def variant_signs(bases, gram, component_generators):
    label_maps, unused_generators = component_label_maps(bases, gram)
    assert unused_generators == component_generators
    result = {}
    for source in range(8):
        for target in range(8):
            variants = graph_variants(bases[source], bases[target], gram)
            by_sign = {}
            for selection in variants:
                target_coordinates = vector(QQ, [0, 0, 0])
                for row, image in enumerate(selection):
                    target_coordinates[image] = component_generators[source][row]
                label = label_maps[target][residue_key(target_coordinates)]
                assert label in (1, 3)
                sign = 1 if label == 1 else -1
                assert sign not in by_sign
                by_sign[sign] = selection
            assert set(by_sign) == {-1, 1}
            result[source, target] = by_sign
    return result


def signed_key(signed_permutation):
    permutation, signs = signed_permutation
    return tuple(permutation), tuple(signs)


def compose_signed(left, right):
    left_permutation, left_signs = left
    right_permutation, right_signs = right
    return (
        tuple(
            right_permutation[left_permutation[index]]
            for index in range(8)
        ),
        tuple(
            left_signs[index]
            * right_signs[left_permutation[index]]
            for index in range(8)
        ),
    )


def inverse_signed(value):
    permutation, signs = value
    inverse_permutation = [0] * 8
    inverse_signs = [0] * 8
    for source, target in enumerate(permutation):
        inverse_permutation[target] = source
        inverse_signs[target] = signs[source]
    return tuple(inverse_permutation), tuple(inverse_signs)


def signed_order(value):
    identity = tuple(range(8)), (1,) * 8
    power = identity
    for order in range(1, 49):
        power = compose_signed(power, value)
        if power == identity:
            return order
    raise AssertionError("unexpected signed-permutation order")


def lift_group(
    signed_elements,
    bases,
    gram,
    root_basis_inverse,
    variants_by_sign,
):
    result = []
    seen_matrices = set()
    for signed_permutation in signed_elements:
        permutation, signs = signed_permutation
        target_basis = matrix(
            QQ,
            [
                bases[permutation[source]][
                    variants_by_sign[source, permutation[source]][
                        signs[source]
                    ][row]
                ]
                for source in range(8)
                for row in range(3)
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
        fixed_rank = int(24 - (action - identity_matrix(ZZ, 24)).rank())
        result.append(
            {
                "signed_permutation": signed_permutation,
                "matrix": action,
                "order": signed_order(signed_permutation),
                "fixed_rank": fixed_rank,
            }
        )
    return result


def certify_closure_and_classes(elements, gram):
    by_signed_key = {
        signed_key(row["signed_permutation"]): row for row in elements
    }
    identity = tuple(range(8)), (1,) * 8
    generated = {signed_key(identity): identity}
    generators = []

    def close():
        pending = list(generated.values())
        while pending:
            current = pending.pop()
            for generator in generators:
                for product in (
                    compose_signed(current, generator),
                    compose_signed(generator, current),
                ):
                    key = signed_key(product)
                    assert key in by_signed_key
                    if key not in generated:
                        generated[key] = product
                        pending.append(product)

    for row in elements:
        key = signed_key(row["signed_permutation"])
        if key in generated:
            continue
        generators.append(row["signed_permutation"])
        close()
        if len(generated) == len(elements):
            break
    assert set(generated) == set(by_signed_key)

    remaining = set(by_signed_key)
    conjugacy_classes = []
    signed_group = [row["signed_permutation"] for row in elements]
    while remaining:
        representative_key = min(remaining)
        representative = by_signed_key[representative_key]
        class_keys = {
            signed_key(
                compose_signed(
                    compose_signed(inverse_signed(conjugator), representative["signed_permutation"]),
                    conjugator,
                )
            )
            for conjugator in signed_group
        }
        assert class_keys <= remaining
        remaining.difference_update(class_keys)
        orders = {by_signed_key[key]["order"] for key in class_keys}
        fixed_ranks = {
            by_signed_key[key]["fixed_rank"] for key in class_keys
        }
        assert len(orders) == len(fixed_ranks) == 1
        action = representative["matrix"]
        fixed = row_module_basis(
            (action - identity_matrix(ZZ, 24))
            .transpose()
            .right_kernel_matrix()
            .change_ring(ZZ)
        )
        assert primitive_closure_index(fixed) == 1
        fixed_gram = fixed * gram * fixed.transpose()
        conjugacy_classes.append(
            {
                "class_size": len(class_keys),
                "action_order": next(iter(orders)),
                "fixed_rank": next(iter(fixed_ranks)),
                "fixed_determinant": int(fixed_gram.det()),
                "representative_signed_permutation": {
                    "component_permutation_zero_based": list(
                        representative["signed_permutation"][0]
                    ),
                    "component_diagram_signs": list(
                        representative["signed_permutation"][1]
                    ),
                },
                "representative_matrix": rows(action),
                "element_signed_keys": [
                    [list(key[0]), list(key[1])]
                    for key in sorted(class_keys)
                ],
            }
        )
    conjugacy_classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
            row["representative_signed_permutation"][
                "component_permutation_zero_based"
            ],
            row["representative_signed_permutation"][
                "component_diagram_signs"
            ],
        )
    )
    for index, row in enumerate(conjugacy_classes, start=1):
        row["class_id"] = f"8A3-C{index:02d}"
    return generators, conjugacy_classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "8A3"
    )
    gram = matrix(ZZ, ambient_row["gram"])
    bases = simple_root_components(gram)
    (
        root_basis,
        root_basis_inverse,
        component_generators,
        code,
        code_generators,
        ambient_generator_words,
    ) = recover_glue_code(bases, gram)
    signed_elements = signed_permutation_stabilizer(code, code_generators)
    variants_by_sign = variant_signs(bases, gram, component_generators)
    elements = lift_group(
        signed_elements,
        bases,
        gram,
        root_basis_inverse,
        variants_by_sign,
    )
    generators, conjugacy_classes = certify_closure_and_classes(elements, gram)
    component_image_order = len(
        {row["signed_permutation"][0] for row in elements}
    )
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    assert len(code) == 256
    assert len(code_generators) == 4
    assert len(elements) == 2688
    assert component_image_order == 1344
    assert len(generators) == 7
    assert order_distribution == Counter(
        {1: 1, 2: 99, 3: 224, 4: 588, 6: 672, 7: 384, 8: 336, 14: 384}
    )
    assert fixed_rank_distribution == Counter(
        {4: 720, 6: 1000, 8: 308, 10: 336, 12: 238, 16: 85, 24: 1}
    )
    assert len(conjugacy_classes) == 16
    return {
        "schema": "elkies-k3.8a3-glue-code-residual-group.v1",
        "status": "PASS_EXACT_8A3_GLUE_CODE_AND_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "the complete signed-coordinate stabilizer of the recovered "
                "N(8A3)/A3^8 glue code, its integral lifts to the supplied "
                "Niemeier Gram model, matrix-group closure, and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding "
                "orbits, K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "8A3",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "glue_code": {
            "modulus": 4,
            "length": 8,
            "order": len(code),
            "minimal_generator_count": len(code_generators),
            "minimal_generators": [list(word) for word in code_generators],
            "ambient_basis_generator_words": [
                list(word) for word in ambient_generator_words
            ],
            "codewords": [list(word) for word in sorted(code)],
        },
        "residual_group": {
            "order": len(elements),
            "component_permutation_image_order": component_image_order,
            "signed_coordinate_candidates_tested": 10321920,
            "generator_count": len(generators),
            "generators": [
                {
                    "component_permutation_zero_based": list(value[0]),
                    "component_diagram_signs": list(value[1]),
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
                    "component_permutation_zero_based": list(
                        row["signed_permutation"][0]
                    ),
                    "component_diagram_signs": list(
                        row["signed_permutation"][1]
                    ),
                    "action_order": row["order"],
                    "fixed_rank": row["fixed_rank"],
                    "matrix": rows(row["matrix"]),
                }
                for row in elements
            ],
        },
        "literature_cross_check": {
            "citation": (
                "Cheng-Duncan-Harvey, Umbral Moonshine and the Niemeier "
                "Lattices, Research in the Mathematical Sciences 1 (2014)"
            ),
            "arxiv": "1307.5793",
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
        "elkies-k3/scripts/probe_8a3_glue_code_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("8A3 glue-code residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "8A3GROUP|glue={}|group={}|image={}|classes={}|status=PASS_EXACT".format(
            payload["glue_code"]["order"],
            payload["residual_group"]["order"],
            payload["residual_group"]["component_permutation_image_order"],
            len(payload["residual_group"]["conjugacy_classes"]),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
