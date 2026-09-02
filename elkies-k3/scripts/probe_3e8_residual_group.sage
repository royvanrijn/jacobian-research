#!/usr/bin/env sage-python
"""Recover the exact chamber residual group of the Niemeier lattice N(3E8).

Recover the three E8 root components intrinsically.  Because E8 is
unimodular and has no Dynkin-diagram automorphisms, every chamber-preserving
isometry is a component permutation.  Lift all six permutations to the
supplied integral Niemeier Gram model and certify closure and matrix
conjugacy classes without assuming a catalog coordinate decomposition.
"""

from __future__ import annotations

import argparse
import itertools
import json
import runpy
from collections import Counter
from pathlib import Path

from sage.all import Graph, ZZ, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-3e8-residual-group-v1.json"
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
        assert matrix(ZZ, vectors).rank() == 8
        assert len(vectors) == 240
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
    assert graph.automorphism_group().order() == 1
    return basis


def graph_variant(left, right, gram):
    isomorphic, certificate = dynkin_graph(left, gram).is_isomorphic(
        dynkin_graph(right, gram), certificate=True
    )
    assert isomorphic
    assert dynkin_graph(right, gram).automorphism_group().order() == 1
    return [certificate[index] for index in range(8)]


def compose(left, right):
    return tuple(right[left[index]] for index in range(3))


def inverse(permutation):
    result = [0] * 3
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def permutation_order(permutation):
    identity = tuple(range(3))
    power = identity
    for order in range(1, 7):
        power = compose(power, permutation)
        if power == identity:
            return order
    raise AssertionError("unexpected S3 element order")


def lift_group(bases, gram):
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    assert abs(root_basis.det()) == 1
    root_basis_inverse = root_basis.inverse()
    variants = {
        (source, target): graph_variant(bases[source], bases[target], gram)
        for source in range(3)
        for target in range(3)
    }
    elements = []
    seen = set()
    for permutation in itertools.permutations(range(3)):
        target_basis = matrix(
            ZZ,
            [
                bases[permutation[source]][
                    variants[source, permutation[source]][row]
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
                "permutation": tuple(permutation),
                "matrix": action,
                "order": permutation_order(permutation),
                "fixed_rank": int(
                    24 - (action - identity_matrix(ZZ, 24)).rank()
                ),
            }
        )
    assert len(elements) == 6
    return root_basis, elements


def certify_closure_and_classes(elements, gram):
    by_permutation = {row["permutation"]: row for row in elements}
    identity = tuple(range(3))
    assert identity in by_permutation
    for left in by_permutation:
        assert inverse(left) in by_permutation
        for right in by_permutation:
            product = compose(left, right)
            assert product in by_permutation
            assert (
                by_permutation[left]["matrix"]
                * by_permutation[right]["matrix"]
                == by_permutation[product]["matrix"]
            )

    remaining = set(by_permutation)
    classes = []
    while remaining:
        representative_key = min(remaining)
        representative = by_permutation[representative_key]
        class_keys = {
            compose(compose(inverse(conjugator), representative_key), conjugator)
            for conjugator in by_permutation
        }
        assert class_keys <= remaining
        remaining.difference_update(class_keys)
        orders = {by_permutation[key]["order"] for key in class_keys}
        fixed_ranks = {by_permutation[key]["fixed_rank"] for key in class_keys}
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
        assert primitive_closure_index(fixed) == 1
        fixed_gram = fixed * gram * fixed.transpose()
        classes.append(
            {
                "class_size": len(class_keys),
                "action_order": next(iter(orders)),
                "fixed_rank": fixed_rank,
                "fixed_determinant": int(fixed_gram.det()),
                "representative_component_permutation_zero_based": list(
                    representative_key
                ),
                "representative_matrix": rows(action),
                "element_keys": [list(key) for key in sorted(class_keys)],
            }
        )
    classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["fixed_determinant"],
            row["representative_component_permutation_zero_based"],
        )
    )
    for index, row in enumerate(classes, start=1):
        row["class_id"] = f"3E8-C{index:02d}"
    return classes


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row for row in catalog["rooted_niemeier_lattices"] if row["label"] == "3E8"
    )
    gram = matrix(ZZ, ambient_row["gram"])
    bases = [simple_roots(gram, roots) for roots in component_roots(gram)]
    root_basis, elements = lift_group(bases, gram)
    classes = certify_closure_and_classes(elements, gram)
    order_distribution = Counter(row["order"] for row in elements)
    fixed_rank_distribution = Counter(row["fixed_rank"] for row in elements)
    assert order_distribution == Counter({1: 1, 2: 3, 3: 2})
    assert fixed_rank_distribution == Counter({24: 1, 16: 3, 8: 2})
    assert [row["fixed_determinant"] for row in classes] == [1, 256, 6561]
    return {
        "schema": "elkies-k3.3e8-residual-group.v1",
        "status": "PASS_EXACT_3E8_RESIDUAL_GROUP",
        "proof_scope": {
            "proved": (
                "intrinsic recovery of the three E8 root components; the complete "
                "component-permutation chamber residual group, its integral lifts "
                "to the supplied Niemeier Gram model, closure, and conjugacy classes"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding orbits, "
                "K3 admissibility, or determinant-band completeness"
            ),
        },
        "ambient_label": "3E8",
        "root_lattice": {
            "simple_root_bases_in_ambient": [rows(basis) for basis in bases],
            "combined_basis_in_ambient": rows(root_basis),
            "index_in_niemeier": abs(int(root_basis.det())),
        },
        "residual_group": {
            "order": len(elements),
            "component_permutation_image_order": len(elements),
            "component_kernel_order": 1,
            "generator_count": 2,
            "generators_component_permutations_zero_based": [
                [1, 0, 2],
                [1, 2, 0],
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
                    "component_permutation_zero_based": list(row["permutation"]),
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
                "(2020), 3E8 case"
            ),
            "arxiv": "2002.03707",
            "expected_structure": "S3 component-permutation group",
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
        "elkies-k3/scripts/probe_3e8_residual_group.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("3E8 residual-group artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "3E8GROUP|group={}|image={}|kernel={}|classes={}|status=PASS_EXACT".format(
            payload["residual_group"]["order"],
            payload["residual_group"]["component_permutation_image_order"],
            payload["residual_group"]["component_kernel_order"],
            len(payload["residual_group"]["conjugacy_classes"]),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
