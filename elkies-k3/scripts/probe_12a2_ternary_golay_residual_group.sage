#!/usr/bin/env sage-python
"""Recover the ternary Golay glue and full residual group of N(12A2).

Recover the twelve A2 root components and their index-729 root basis from the
supplied Niemeier Gram matrix. Express N/A2^12 as an intrinsic six-dimensional
ternary code, certify its Golay parameters, and compute its complete monomial
automorphism group. Lift a deterministic generator set and every conjugacy-
class representative to exact integral ambient isometries, then certify fixed
lattices and class masses.

The 2^12 * 12! chamber envelope is never enumerated. Code automorphisms are
computed before any rank-seven shell search.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import runpy
from collections import Counter
from pathlib import Path

from sage.all import (
    GF,
    Graph,
    LinearCode,
    Permutation,
    PermutationGroup,
    QQ,
    ZZ,
    identity_matrix,
    matrix,
    vector,
)
from sage.libs.gap.libgap import libgap


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-12a2-ternary-golay-residual-group-v1.json"
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

F3 = GF(3)


def residue_key(values):
    return tuple(QQ(value - value.floor()) for value in values)


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
        assert matrix(ZZ, vectors).rank() == 2
        assert len(vectors) == 6
        result.append(vectors)
    assert len(result) == 12
    return result


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(2))
    if cartan[0, 1] == -1:
        graph.add_edge(0, 1)
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
    assert basis.nrows() == basis.rank() == 2
    assert basis * gram * basis.transpose() == matrix(ZZ, [[2, -1], [-1, 2]])
    assert dynkin_graph(basis, gram).automorphism_group().order() == 2
    return basis


def graph_variants(left, right, gram):
    left_graph = dynkin_graph(left, gram)
    right_graph = dynkin_graph(right, gram)
    isomorphic, initial = left_graph.is_isomorphic(right_graph, certificate=True)
    assert isomorphic
    variants = [
        [automorphism(initial[index]) for index in range(2)]
        for automorphism in right_graph.automorphism_group()
    ]
    assert len(variants) == 2
    return variants


def component_label_data(bases, gram):
    label_maps = []
    representatives = []
    for basis in bases:
        cartan = basis * gram * basis.transpose()
        generator = cartan.inverse().row(0)
        reps = [vector(QQ, [0, 0]), generator, 2 * generator]
        mapping = {residue_key(value): index for index, value in enumerate(reps)}
        assert len(mapping) == 3
        label_maps.append(mapping)
        representatives.append(reps)
    return label_maps, representatives


def recover_glue_code(bases, gram):
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 729
    root_basis_inverse = root_basis.inverse()
    label_maps, representatives = component_label_data(bases, gram)
    ambient_words = []
    for ambient_index in range(24):
        ambient = vector(ZZ, [int(index == ambient_index) for index in range(24)])
        coordinates = ambient * root_basis_inverse
        word = []
        for component_index in range(12):
            local = coordinates[2 * component_index : 2 * component_index + 2]
            word.append(label_maps[component_index][residue_key(local)])
        ambient_words.append(word)
    generator_matrix = matrix(F3, ambient_words).row_space().basis_matrix()
    code = LinearCode(generator_matrix)
    weight_distribution = [0] * 13
    for coefficients in itertools.product(range(3), repeat=6):
        word = vector(F3, coefficients) * generator_matrix
        weight_distribution[sum(bool(value) for value in word)] += 1
    assert code.length() == 12
    assert code.dimension() == 6
    assert code.cardinality() == 729
    assert generator_matrix * generator_matrix.transpose() == matrix(F3, 6)
    assert weight_distribution == [
        1,
        0,
        0,
        0,
        0,
        0,
        264,
        0,
        0,
        440,
        0,
        0,
        24,
    ]
    return (
        root_basis,
        root_basis_inverse,
        label_maps,
        representatives,
        ambient_words,
        code,
        weight_distribution,
    )


def variant_multiplier_data(bases, gram, label_maps, representatives):
    result = {}
    for source in range(12):
        for target in range(12):
            by_multiplier = {}
            variants = graph_variants(bases[source], bases[target], gram)
            for selection_index, selection in enumerate(variants):
                target_coordinates = vector(QQ, [0, 0])
                source_coordinates = representatives[source][1]
                for row, image in enumerate(selection):
                    target_coordinates[image] = source_coordinates[row]
                multiplier = label_maps[target][residue_key(target_coordinates)]
                assert multiplier in (1, 2)
                assert multiplier not in by_multiplier
                by_multiplier[multiplier] = {
                    "selection_index": selection_index,
                    "selection": selection,
                }
            assert set(by_multiplier) == {1, 2}
            result[source, target] = by_multiplier
    return result


def semimonomial_matrix(transformation):
    result = matrix(F3, transformation * identity_matrix(F3, 12))
    assert all(sum(bool(value) for value in row) == 1 for row in result.rows())
    assert all(sum(bool(result[row, column]) for row in range(12)) == 1 for column in range(12))
    return result


def monomial_data(action):
    permutation = []
    multipliers = []
    for source in range(12):
        targets = [target for target in range(12) if action[source, target]]
        assert len(targets) == 1
        target = targets[0]
        permutation.append(target)
        multipliers.append(int(action[source, target]))
    assert sorted(permutation) == list(range(12))
    assert all(value in (1, 2) for value in multipliers)
    return tuple(permutation), tuple(multipliers)


def signed_permutation(action):
    permutation, multipliers = monomial_data(action)
    images = [
        target + 1 if multiplier == 1 else target + 13
        for target, multiplier in zip(permutation, multipliers)
    ]
    images += [image + 12 if image <= 12 else image - 12 for image in images]
    result = Permutation(images)
    assert result.size() <= 24
    return result


def signed_permutation_matrix(element):
    action = matrix(F3, 12)
    for source in range(12):
        image = int(element(source + 1))
        if image <= 12:
            target = image - 1
            multiplier = 1
        else:
            target = image - 13
            multiplier = 2
        action[source, target] = multiplier
        assert int(element(source + 13)) == (
            target + 13 if multiplier == 1 else target + 1
        )
    return action


def ambient_action(
    monomial_action,
    bases,
    root_basis_inverse,
    variants_by_multiplier,
    gram,
):
    permutation, multipliers = monomial_data(monomial_action)
    selections = []
    target_rows = []
    for source, (target, multiplier) in enumerate(zip(permutation, multipliers)):
        variant = variants_by_multiplier[source, target][multiplier]
        selections.append(variant["selection_index"])
        target_rows.extend(
            bases[target][variant["selection"][row]] for row in range(2)
        )
    rational_action = root_basis_inverse * matrix(QQ, target_rows)
    assert all(value.denominator() == 1 for value in rational_action.list())
    action = matrix(ZZ, rational_action)
    assert abs(action.det()) == 1
    assert action * gram * action.transpose() == gram
    return action, permutation, multipliers, tuple(selections)


def fixed_lattice_data(action, gram):
    identity = identity_matrix(ZZ, 24)
    fixed = row_module_basis(
        (action - identity).transpose().right_kernel_matrix().change_ring(ZZ)
    )
    fixed_rank = fixed.nrows()
    if fixed_rank:
        assert primitive_closure_index(fixed) == 1
        fixed_determinant = int((fixed * gram * fixed.transpose()).det())
    else:
        fixed_determinant = 1
    return fixed_rank, fixed_determinant


def recover_residual_group(
    code,
    bases,
    root_basis_inverse,
    variants_by_multiplier,
    gram,
):
    libgap.set_seed(0)
    transformations, order = code.automorphism_group_gens(equivalence="linear")
    assert int(order) == 190080
    assert len(transformations) == 3
    monomial_generators = [semimonomial_matrix(row) for row in transformations]
    for generator in monomial_generators:
        assert LinearCode(code.generator_matrix() * generator) == code
    component_permutation_group = PermutationGroup(
        [
            Permutation([target + 1 for target in monomial_data(row)[0]])
            for row in monomial_generators
        ]
    )
    component_image_order = int(component_permutation_group.order())
    assert component_image_order == 95040
    assert int(order) // component_image_order == 2
    permutation_generators = [signed_permutation(row) for row in monomial_generators]
    permutation_group = PermutationGroup(permutation_generators)
    assert permutation_group.order() == int(order)

    generator_rows = []
    ambient_generator_keys = set()
    for monomial, permutation_element in zip(
        monomial_generators, permutation_generators
    ):
        action, component_permutation, multipliers, selections = ambient_action(
            monomial,
            bases,
            root_basis_inverse,
            variants_by_multiplier,
            gram,
        )
        key = matrix_key(action)
        assert key not in ambient_generator_keys
        ambient_generator_keys.add(key)
        assert action.multiplicative_order() == permutation_element.order()
        generator_rows.append(
            {
                "component_permutation_zero_based": list(component_permutation),
                "component_discriminant_multipliers_mod_3": list(multipliers),
                "component_diagram_variant_indices": list(selections),
                "signed_coordinate_permutation_on_24_points_one_based": list(
                    map(int, (permutation_element(index) for index in range(1, 25)))
                ),
                "action_order": int(permutation_element.order()),
                "matrix": rows(action),
            }
        )

    classes = []
    class_mass = 0
    for representative in permutation_group.conjugacy_classes_representatives():
        monomial = signed_permutation_matrix(representative)
        assert LinearCode(code.generator_matrix() * monomial) == code
        action, component_permutation, multipliers, selections = ambient_action(
            monomial,
            bases,
            root_basis_inverse,
            variants_by_multiplier,
            gram,
        )
        fixed_rank, fixed_determinant = fixed_lattice_data(action, gram)
        class_size = int(permutation_group.conjugacy_class(representative).cardinality())
        class_mass += class_size
        classes.append(
            {
                "class_size": class_size,
                "action_order": int(representative.order()),
                "fixed_rank": fixed_rank,
                "fixed_determinant": fixed_determinant,
                "representative_component_permutation_zero_based": list(
                    component_permutation
                ),
                "representative_component_discriminant_multipliers_mod_3": list(
                    multipliers
                ),
                "representative_component_diagram_variant_indices": list(
                    selections
                ),
                "representative_signed_coordinate_permutation_on_24_points_one_based": list(
                    map(int, (representative(index) for index in range(1, 25)))
                ),
                "representative_matrix": rows(action),
            }
        )
    assert class_mass == int(order)
    assert len(classes) == 26
    classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
            row["representative_component_permutation_zero_based"],
            row["representative_component_discriminant_multipliers_mod_3"],
        )
    )
    for index, row in enumerate(classes, start=1):
        row["class_id"] = f"12A2-C{index:02d}"
    return int(order), component_image_order, generator_rows, classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "12A2"
    )
    gram = matrix(ZZ, ambient_row["gram"])
    bases = [simple_roots(gram, roots) for roots in component_roots(gram)]
    (
        root_basis,
        root_basis_inverse,
        label_maps,
        representatives,
        ambient_words,
        code,
        weight_distribution,
    ) = recover_glue_code(bases, gram)
    variants_by_multiplier = variant_multiplier_data(
        bases, gram, label_maps, representatives
    )
    order, component_image_order, generators, classes = recover_residual_group(
        code,
        bases,
        root_basis_inverse,
        variants_by_multiplier,
        gram,
    )
    order_distribution = Counter()
    fixed_rank_class_distribution = Counter()
    for row in classes:
        order_distribution[row["action_order"]] += row["class_size"]
        fixed_rank_class_distribution[row["fixed_rank"]] += 1
    assert order_distribution == Counter(
        {
            1: 1,
            2: 991,
            3: 4400,
            4: 12672,
            5: 9504,
            6: 36080,
            8: 47520,
            10: 9504,
            11: 17280,
            12: 15840,
            20: 19008,
            22: 17280,
        }
    )
    assert fixed_rank_class_distribution == Counter(
        {24: 1, 16: 1, 12: 3, 10: 1, 8: 4, 6: 5, 4: 6, 2: 5}
    )
    return {
        "schema": "elkies-k3.12a2-ternary-golay-residual-group.v1",
        "status": "PASS_EXACT_12A2_TERNARY_GOLAY_AND_FULL_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "intrinsic recovery of the twelve A2 root components, index-729 "
                "root quotient, and [12,6,6] self-dual ternary Golay glue code; "
                "the complete order-190080 monomial code automorphism group; "
                "integral ambient lifts of deterministic generators and every "
                "conjugacy-class representative; class masses and fixed lattices"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "12A2",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "glue_code": {
            "alphabet": "GF(3), with A2 diagram reversal acting by -1",
            "length": int(code.length()),
            "dimension": int(code.dimension()),
            "order": int(code.cardinality()),
            "minimum_weight": int(next(
                index
                for index, count in enumerate(weight_distribution)
                if index and count
            )),
            "self_dual": True,
            "weight_distribution": weight_distribution,
            "generator_matrix": rows(code.generator_matrix()),
            "ambient_basis_generator_words": ambient_words,
        },
        "residual_group": {
            "order": order,
            "structure_cross_check": "2.M12",
            "component_permutation_image_order": component_image_order,
            "central_diagram_kernel_order": order // component_image_order,
            "ambient_chamber_envelope_order_not_enumerated": (
                2**12 * math.factorial(12)
            ),
            "ambient_chamber_envelope_factorization": "2^12 * 12!",
            "generator_count": len(generators),
            "generators": generators,
            "conjugacy_class_count": len(classes),
            "order_distribution": {
                str(key): value for key, value in sorted(order_distribution.items())
            },
            "fixed_rank_class_distribution": {
                str(key): value
                for key, value in sorted(fixed_rank_class_distribution.items())
            },
            "conjugacy_classes": classes,
        },
        "literature_cross_check": {
            "citation": (
                "Chenevier, The Characteristic Masses of Niemeier Lattices "
                "(2020), 12A2 case"
            ),
            "arxiv": "2002.03707",
            "expected_structure": (
                "the ternary Golay monomial automorphism group 2.M12 of order "
                "190080"
            ),
            "role": (
                "independent naming/order cross-check only; the artifact recovers "
                "the code and group from the supplied lattice"
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
        "elkies-k3/scripts/probe_12a2_ternary_golay_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("12A2 residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "12A2GROUP|glue={}|dimension={}|group={}|classes={}|status=PASS_EXACT".format(
            payload["glue_code"]["order"],
            payload["glue_code"]["dimension"],
            payload["residual_group"]["order"],
            payload["residual_group"]["conjugacy_class_count"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
