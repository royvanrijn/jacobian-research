#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of the Niemeier lattice N(2A12).

Recover the two A12 root components and their combined index-thirteen root
basis intrinsically. Exhaust all 2^2 * 2! component/diagram chamber maps; a
candidate preserves the Niemeier glue exactly when its lift in the supplied
ambient basis is integral. Certify the retained matrix group, fixed lattices,
and conjugacy classes.
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
    ROOT / "artifacts/generated-results/elkies-k3-2a12-residual-group-v1.json"
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
        assert matrix(ZZ, vectors).rank() == 12
        assert len(vectors) == 156
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
    assert sorted(graph.degree(vertex) for vertex in graph) == [1, 1] + [2] * 10
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
    variants = sorted(variants, key=lambda value: value != list(range(12)))
    assert len(variants) == 2
    if left is right:
        assert variants[0] == list(range(12))
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
        generated = set(cosets)
        current = generator
        for _multiple in range(1, 14):
            generated |= {add_residues(value, current) for value in cosets}
            current = add_residues(current, generator)
        cosets = generated
    assert len(cosets) == 13
    return generators, cosets


def matrix_order(action):
    identity = identity_matrix(ZZ, 24)
    power = identity
    for order in range(1, 9):
        power *= action
        if power == identity:
            return order
    raise AssertionError("unexpected residual action order")


def candidate_actions(bases, root_basis_inverse, gram):
    variants = {
        (source, target): graph_variants(bases[source], bases[target], gram)
        for source in range(2)
        for target in range(2)
    }
    candidates = []
    elements = []
    seen = set()
    for permutation in itertools.permutations(range(2)):
        for selections in itertools.product(range(2), repeat=2):
            target_basis = matrix(
                QQ,
                [
                    bases[permutation[source]][
                        variants[source, permutation[source]][selections[source]][row]
                    ]
                    for source in range(2)
                    for row in range(12)
                ],
            )
            rational_action = root_basis_inverse * target_basis
            integral = all(
                entry.denominator() == 1 for entry in rational_action.list()
            )
            candidates.append(
                {
                    "component_permutation_zero_based": list(permutation),
                    "component_diagram_variant_indices": list(selections),
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
                    "component_permutation_zero_based": tuple(permutation),
                    "component_diagram_variant_indices": tuple(selections),
                    "matrix": action,
                    "order": matrix_order(action),
                    "fixed_rank": int(
                        24 - (action - identity_matrix(ZZ, 24)).rank()
                    ),
                }
            )
    assert len(candidates) == 8
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
    generated = {identity_index}
    generator_index = next(
        index for index, row in enumerate(elements) if row["order"] == 4
    )
    current = identity_index
    for _power in range(4):
        generated.add(current)
        current = multiplication[current, generator_index]
    assert generated == set(range(4))

    remaining = set(range(4))
    classes = []
    while remaining:
        representative_index = min(remaining)
        class_indices = {
            multiplication[
                multiplication[inverses[conjugator], representative_index],
                conjugator,
            ]
            for conjugator in range(4)
        }
        assert class_indices <= remaining
        remaining.difference_update(class_indices)
        representative = elements[representative_index]
        orders = {elements[index]["order"] for index in class_indices}
        fixed_ranks = {elements[index]["fixed_rank"] for index in class_indices}
        assert len(orders) == len(fixed_ranks) == 1
        action = representative["matrix"]
        fixed = row_module_basis(
            (action - identity)
            .transpose()
            .right_kernel_matrix()
            .change_ring(ZZ)
        )
        fixed_rank = next(iter(fixed_ranks))
        assert fixed.nrows() == fixed_rank
        if fixed_rank:
            assert primitive_closure_index(fixed) == 1
            fixed_determinant = int((fixed * gram * fixed.transpose()).det())
        else:
            fixed_determinant = 1
        classes.append(
            {
                "class_size": len(class_indices),
                "action_order": next(iter(orders)),
                "fixed_rank": fixed_rank,
                "fixed_determinant": fixed_determinant,
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
        row["class_id"] = f"2A12-C{index:02d}"
    return generator_index, classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row for row in catalog["rooted_niemeier_lattices"] if row["label"] == "2A12"
    )
    gram = matrix(ZZ, ambient_row["gram"])
    bases = [simple_roots(gram, roots) for roots in component_roots(gram)]
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 13
    root_basis_inverse = root_basis.inverse()
    quotient_generators, quotient = quotient_cosets(root_basis_inverse)
    candidates, elements = candidate_actions(bases, root_basis_inverse, gram)
    generator_index, classes = certify_closure_and_classes(elements, gram)
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    assert order_distribution == Counter({1: 1, 2: 1, 4: 2})
    assert fixed_rank_distribution == Counter({24: 1, 12: 1, 6: 2})
    assert len(classes) == 4
    return {
        "schema": "elkies-k3.2a12-residual-group.v1",
        "status": "PASS_EXACT_2A12_GLUE_AND_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "intrinsic recovery of the two A12 root components and order-13 "
                "root quotient; all eight component/diagram chamber candidates; "
                "the complete integral residual lift group, cyclic generation, "
                "fixed lattices, and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "2A12",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
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
            "component_permutation_image_order": len(
                {row["component_permutation_zero_based"] for row in elements}
            ),
            "component_diagram_kernel_order": sum(
                row["component_permutation_zero_based"] == tuple(range(2))
                for row in elements
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
                    "component_permutation_zero_based": list(
                        row["component_permutation_zero_based"]
                    ),
                    "component_diagram_variant_indices": list(
                        row["component_diagram_variant_indices"]
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
                "Chenevier, The Characteristic Masses of Niemeier Lattices "
                "(2020), 2A12 case"
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
        "elkies-k3/scripts/probe_2a12_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("2A12 residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "2A12GROUP|quotient={}|group={}|classes={}|status=PASS_EXACT".format(
            payload["root_lattice_quotient"]["order"],
            payload["residual_group"]["order"],
            len(payload["residual_group"]["conjugacy_classes"]),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
