#!/usr/bin/env sage-python
"""Recover the exact chamber residual sections of N(4A6) and N(4E6).

For each backend, recover the four rank-six root components intrinsically,
enumerate all 4!*2^4 component/diagram maps of the root lattice, and retain
exactly those extending integrally to the Niemeier overlattice.  The output
records the complete matrix group, its component-permutation image, a closure
certificate, and fixed-lattice data by action profile.  This is group
infrastructure only; it does not enumerate rank-seven auxiliaries.
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
    / "artifacts/generated-results/elkies-k3-4a6-4e6-residual-sections-v1.json"
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


def finite_order(action):
    identity = identity_matrix(ZZ, action.nrows())
    power = identity
    for order in range(1, 25):
        power *= action
        if power == identity:
            return order
    raise AssertionError("unexpected residual action order")


def cycle_type(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        current = min(unseen)
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


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
        assert matrix(ZZ, component_roots).rank() == 6
        components.append(component_roots)
    assert len(components) == 4

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
        assert basis.nrows() == basis.rank() == 6
        assert all(cartan[index, index] == 2 for index in range(6))
        assert sum(
            cartan[row, column] == -1
            for row in range(6)
            for column in range(row)
        ) == 5
        bases.append(basis)
    return bases, [len(component) for component in components]


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(6))
    graph.add_edges(
        (row, column)
        for row in range(6)
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
        [automorphism(initial[index]) for index in range(6)]
        for automorphism in right_graph.automorphism_group()
    ]
    assert len(variants) == 2
    return variants


def enumerate_section(label, gram):
    bases, root_counts = simple_root_components(gram)
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    root_basis_inverse = root_basis.inverse()
    variants = {
        (source, target): graph_variants(
            bases[source], bases[target], gram
        )
        for source in range(4)
        for target in range(4)
    }

    elements = []
    seen = set()
    candidates = 0
    for permutation in itertools.permutations(range(4)):
        for selections in itertools.product(
            *[
                variants[source, permutation[source]]
                for source in range(4)
            ]
        ):
            candidates += 1
            target_basis = matrix(
                QQ,
                [
                    bases[permutation[source]][selections[source][row]]
                    for source in range(4)
                    for row in range(6)
                ],
            )
            action = root_basis_inverse * target_basis
            if not all(entry.denominator() == 1 for entry in action.list()):
                continue
            action = matrix(ZZ, action)
            if abs(action.det()) != 1:
                continue
            assert action * gram * action.transpose() == gram
            key = matrix_key(action)
            assert key not in seen
            seen.add(key)
            fixed = row_module_basis(
                (action - identity_matrix(ZZ, 24))
                .transpose()
                .right_kernel_matrix()
                .change_ring(ZZ)
            )
            assert primitive_closure_index(fixed) == 1
            fixed_gram = fixed * gram * fixed.transpose()
            elements.append(
                {
                    "matrix": action,
                    "component_permutation": tuple(permutation),
                    "cycle_type": cycle_type(permutation),
                    "order": finite_order(action),
                    "fixed_rank": fixed.nrows(),
                    "fixed_determinant": int(fixed_gram.det()),
                }
            )
    assert candidates == 384
    by_key = {matrix_key(row["matrix"]): row for row in elements}

    identity = identity_matrix(ZZ, 24)
    generated = {matrix_key(identity): identity}
    generators = []

    def close():
        pending = list(generated.values())
        while pending:
            current = pending.pop()
            for generator in generators:
                for product in (current * generator, generator * current):
                    key = matrix_key(product)
                    assert key in by_key
                    if key not in generated:
                        generated[key] = product
                        pending.append(product)

    for row in elements:
        if matrix_key(row["matrix"]) in generated:
            continue
        generators.append(row["matrix"])
        close()
        if len(generated) == len(elements):
            break
    assert set(generated) == set(by_key)

    remaining = set(by_key)
    conjugacy_classes = []
    group_matrices = [row["matrix"] for row in elements]
    while remaining:
        representative_key = min(remaining)
        representative = by_key[representative_key]
        class_keys = {
            matrix_key(
                conjugator.inverse()
                * representative["matrix"]
                * conjugator
            )
            for conjugator in group_matrices
        }
        assert class_keys <= set(by_key)
        assert class_keys <= remaining
        remaining.difference_update(class_keys)
        profiles_in_class = {
            (
                by_key[key]["cycle_type"],
                by_key[key]["order"],
                by_key[key]["fixed_rank"],
                by_key[key]["fixed_determinant"],
            )
            for key in class_keys
        }
        assert len(profiles_in_class) == 1
        profile_in_class = next(iter(profiles_in_class))
        conjugacy_classes.append(
            {
                "component_cycle_type": list(profile_in_class[0]),
                "action_order": profile_in_class[1],
                "fixed_rank": profile_in_class[2],
                "fixed_determinant": profile_in_class[3],
                "class_size": len(class_keys),
                "representative_matrix": rows(representative["matrix"]),
                "element_matrix_keys": sorted(class_keys),
            }
        )
    conjugacy_classes.sort(
        key=lambda row: (
            row["action_order"],
            -row["fixed_rank"],
            row["component_cycle_type"],
            row["fixed_determinant"],
            row["element_matrix_keys"][0],
        )
    )
    for index, row in enumerate(conjugacy_classes, start=1):
        row["class_id"] = f"{label}-C{index:02d}"

    profile = Counter(
        (
            row["cycle_type"],
            row["order"],
            row["fixed_rank"],
            row["fixed_determinant"],
        )
        for row in elements
    )
    profiles = [
        {
            "component_cycle_type": list(key[0]),
            "action_order": key[1],
            "fixed_rank": key[2],
            "fixed_determinant": key[3],
            "number_of_elements": count,
        }
        for key, count in sorted(profile.items())
    ]
    return {
        "ambient_label": label,
        "candidate_component_diagram_maps": candidates,
        "residual_group_order": len(elements),
        "component_permutation_image_order": len(
            {row["component_permutation"] for row in elements}
        ),
        "root_counts_by_component": root_counts,
        "generator_count": len(generators),
        "closure_certified": True,
        "conjugacy_classes": conjugacy_classes,
        "action_profiles": profiles,
        "elements": [
            {
                "component_permutation_zero_based": list(
                    row["component_permutation"]
                ),
                "component_cycle_type": list(row["cycle_type"]),
                "action_order": row["order"],
                "fixed_rank": row["fixed_rank"],
                "fixed_determinant": row["fixed_determinant"],
                "matrix": rows(row["matrix"]),
            }
            for row in elements
        ],
    }


def build(catalog):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_by_label = {
        row["label"]: row
        for row in catalog["rooted_niemeier_lattices"]
    }
    backends = []
    for label in ("4A6", "4E6"):
        gram = matrix(ZZ, ambient_by_label[label]["gram"])
        backends.append(enumerate_section(label, gram))
    assert all(
        row["candidate_component_diagram_maps"] == 384
        for row in backends
    )
    assert [
        (
            row["ambient_label"],
            row["residual_group_order"],
            row["component_permutation_image_order"],
            len(row["conjugacy_classes"]),
        )
        for row in backends
    ] == [("4A6", 24, 12, 7), ("4E6", 48, 24, 8)]
    return {
        "schema": "elkies-k3.4a6-4e6-residual-sections.v1",
        "status": "PASS_EXACT_RESIDUAL_SECTIONS",
        "proof_scope": {
            "proved": (
                "complete chamber-preserving component/diagram lift groups "
                "for the supplied N(4A6) and N(4E6) Gram models"
            ),
            "not_proved": (
                "rank-seven auxiliary enumeration, full Weyl embedding "
                "orbits, determinant-band completeness, or K3 admissibility"
            ),
        },
        "backends": backends,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/probe_4a6_4e6_residual_sections.sage"
        ),
    }


parser = argparse.ArgumentParser()
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
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == encoded
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)
print(
    "A6E6SECTION|{}|status=PASS_EXACT".format(
        "|".join(
            "{}:group={},image={}".format(
                row["ambient_label"],
                row["residual_group_order"],
                row["component_permutation_image_order"],
            )
            for row in payload["backends"]
        )
    )
)
