#!/usr/bin/env sage
"""Enumerate D5-root anchors in the surviving Niemeier lattices.

status: ACTIVE_PROOF
claim: The pinned auxiliary D5 root sublattice has sixteen anchor orbits under
  the full automorphism groups of the thirteen D5-admissible rooted Niemeier
  lattices.  An explicit D5 basis is retained for every orbit.
inputs: artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json
outputs: artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json
supersedes/superseded-by: none
"""

import argparse
import itertools
import json
from collections import deque
from pathlib import Path

from sage.all import CartanMatrix, Graph, QQ, ZZ, identity_matrix, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json"
D5_GRAM = CartanMatrix(["D", 5])

TYPE_FROM_INVARIANT = {
    (rank, rank * (rank + 1)): ("A", rank) for rank in range(1, 25)
}
TYPE_FROM_INVARIANT.update(
    {(rank, 2 * rank * (rank - 1)): ("D", rank) for rank in range(4, 25)}
)
TYPE_FROM_INVARIANT.update({(6, 72): ("E", 6), (7, 126): ("E", 7), (8, 240): ("E", 8)})


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def signed_roots_and_components(gram):
    minimum_data = pari(gram).qfminim(2)
    representatives = matrix(ZZ, minimum_data[2].sage()).transpose()
    roots = [vector(ZZ, row) for row in representatives.rows()]
    roots += [-root for root in roots]
    unseen = set(range(len(roots)))
    components = []
    while unseen:
        component = {min(unseen)}
        frontier = list(component)
        unseen.difference_update(component)
        while frontier:
            current = frontier.pop()
            neighbours = {
                index
                for index in unseen
                if roots[current] * gram * roots[index] != 0
            }
            component.update(neighbours)
            unseen.difference_update(neighbours)
            frontier.extend(sorted(neighbours))
        component_roots = [roots[index] for index in sorted(component)]
        rank = matrix(ZZ, component_roots).rank()
        invariant = (int(rank), len(component_roots))
        components.append(
            {
                "roots": component_roots,
                "invariant": invariant,
                "type": TYPE_FROM_INVARIANT[invariant],
            }
        )
    return roots, components


def simple_roots(gram, roots):
    for trial in range(1, 100):
        chamber_vector = vector(
            ZZ,
            [
                (index + 1) ** 2 + trial * (index + 1) + trial**2
                for index in range(gram.nrows())
            ],
        )
        values = [root * gram * chamber_vector for root in roots]
        if all(value != 0 for value in values):
            break
    else:
        raise AssertionError("failed to choose a regular chamber vector")
    positive = [root for root, value in zip(roots, values) if value > 0]
    positive_set = {tuple(map(int, root)) for root in positive}
    result = matrix(
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
    cartan = result * gram * result.transpose()
    assert result.rank() == result.nrows()
    assert all(cartan[index, index] == 2 for index in range(cartan.nrows()))
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(cartan.nrows())
        for column in range(cartan.ncols())
        if row != column
    )
    return result


def dynkin_graph(cartan):
    graph = Graph()
    graph.add_vertices(range(cartan.nrows()))
    graph.add_edges(
        (row, column)
        for row in range(cartan.nrows())
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def graph_isometry_variants(left_basis, right_basis, gram):
    left = dynkin_graph(left_basis * gram * left_basis.transpose())
    right = dynkin_graph(right_basis * gram * right_basis.transpose())
    isomorphic, initial = left.is_isomorphic(right, certificate=True)
    assert isomorphic
    return [
        [automorphism(initial[index]) for index in range(left_basis.nrows())]
        for automorphism in right.automorphism_group()
    ]


def lifted_component_permutations(gram, bases, component_types):
    root_basis = matrix(ZZ, [row for basis in bases for row in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == 24
    groups = {}
    for index, component_type in enumerate(component_types):
        groups.setdefault(component_type, []).append(index)

    lifted = []
    type_permutations = itertools.product(
        *[list(itertools.permutations(indices)) for indices in groups.values()]
    )
    for choices in type_permutations:
        permutation = list(range(len(bases)))
        for indices, images in zip(groups.values(), choices):
            for source, target in zip(indices, images):
                permutation[source] = target

        variants = [
            graph_isometry_variants(bases[source], bases[target], gram)
            for source, target in enumerate(permutation)
        ]
        for selections in itertools.product(*variants):
            target_basis = matrix(
                QQ,
                [
                    bases[permutation[source]][selections[source][row]]
                    for source in range(len(bases))
                    for row in range(bases[source].nrows())
                ],
            )
            automorphism = root_basis.inverse() * target_basis
            if not all(entry.denominator() == 1 for entry in automorphism.list()):
                continue
            automorphism = matrix(ZZ, automorphism)
            if abs(automorphism.det()) != 1:
                continue
            assert automorphism * gram * automorphism.transpose() == gram
            lifted.append(permutation)
            break
    return sorted(set(map(tuple, lifted)))


def component_orbits(permutations, eligible):
    unseen = set(eligible)
    result = []
    while unseen:
        orbit = {min(unseen)}
        frontier = list(orbit)
        while frontier:
            current = frontier.pop()
            for permutation in permutations:
                image = permutation[current]
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    return result


def first_d5_basis(gram, roots):
    degrees = [
        sum(D5_GRAM[index, other] == -1 for other in range(5))
        for index in range(5)
    ]
    order = sorted(range(5), key=lambda index: (-degrees[index], index))
    selected = {order[0]: roots[0]}

    def extend(position):
        if position == 5:
            result = matrix(ZZ, [selected[index] for index in range(5)])
            assert result * gram * result.transpose() == D5_GRAM
            return result
        target = order[position]
        for root in roots:
            if all(
                root * gram * selected[index] == D5_GRAM[target, index]
                for index in selected
            ):
                selected[target] = root
                result = extend(position + 1)
                if result is not None:
                    return result
                del selected[target]
        return None

    result = extend(1)
    assert result is not None
    return result


def hnf_key(value):
    return tuple(tuple(map(int, row)) for row in value.hermite_form().rows())


def enumerate_exceptional_d5_orbit(family_rank):
    gram = CartanMatrix(["E", family_rank])
    minimum_data = pari(gram).qfminim(2)
    representatives = matrix(ZZ, minimum_data[2].sage()).transpose()
    roots = [vector(ZZ, row) for row in representatives.rows()]
    roots += [-root for root in roots]
    fixed_root = roots[0]
    degrees = [
        sum(D5_GRAM[index, other] == -1 for other in range(5))
        for index in range(5)
    ]
    order = sorted(range(5), key=lambda index: (-degrees[index], index))
    selected = {order[0]: fixed_root}
    subsystems = set()

    def extend(position):
        if position == 5:
            basis = matrix(ZZ, [selected[index] for index in range(5)])
            subsystems.add(hnf_key(basis))
            return
        target = order[position]
        for root in roots:
            if all(
                root * gram * selected[index] == D5_GRAM[target, index]
                for index in selected
            ):
                selected[target] = root
                extend(position + 1)
                del selected[target]

    extend(1)
    reflections = []
    identity = identity_matrix(ZZ, family_rank)
    for root in roots:
        if root * gram * fixed_root == 0:
            reflections.append(identity - gram * root.column() * root.row())

    unseen = set(subsystems)
    orbit_count = 0
    while unseen:
        orbit_count += 1
        orbit = {min(unseen)}
        frontier = deque(orbit)
        while frontier:
            basis = matrix(ZZ, frontier.popleft())
            for reflection in reflections:
                image = hnf_key(basis * reflection)
                assert image in subsystems
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        unseen.difference_update(orbit)
    return {
        "fixed_root_D5_subsystems": len(subsystems),
        "fixed_root_stabilizer_orbits": orbit_count,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument(
    "--deep",
    action="store_true",
    help="Re-enumerate the exceptional E6/E7/E8 D5 subsystems.",
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

catalog_payload = json.loads(CATALOG.read_text())
admissible_labels = set(
    catalog_payload["accounting"]["classes_requiring_embedding_enumeration"]
)

if arguments.deep:
    exceptional = {
        str(rank): enumerate_exceptional_d5_orbit(rank) for rank in (6, 7, 8)
    }
else:
    exceptional = json.loads(arguments.output.read_text())[
        "internal_D5_weyl_orbits"
    ]["exceptional_enumeration"]

assert exceptional == {
    "6": {"fixed_root_D5_subsystems": 15, "fixed_root_stabilizer_orbits": 1},
    "7": {"fixed_root_D5_subsystems": 120, "fixed_root_stabilizer_orbits": 1},
    "8": {"fixed_root_D5_subsystems": 1260, "fixed_root_stabilizer_orbits": 1},
}

ambient_rows = []
anchor_rows = []
for entry in catalog_payload["rooted_niemeier_lattices"]:
    if entry["label"] not in admissible_labels:
        continue
    gram = matrix(ZZ, entry["gram"])
    unused_roots, components = signed_roots_and_components(gram)
    bases = [simple_roots(gram, component["roots"]) for component in components]
    component_types = [component["type"] for component in components]
    permutations = lifted_component_permutations(gram, bases, component_types)
    eligible = [
        index
        for index, component in enumerate(components)
        if (component["type"][0] == "D" and component["type"][1] >= 5)
        or component["type"][0] == "E"
    ]
    orbits = component_orbits(permutations, eligible)
    ambient_anchors = []
    for local_index, orbit in enumerate(orbits, start=1):
        representative = orbit[0]
        d5_basis = first_d5_basis(gram, components[representative]["roots"])
        anchor = {
            "anchor_index": local_index,
            "component_orbit": orbit,
            "component_type": list(components[representative]["type"]),
            "D5_basis_in_ambient": rows(d5_basis),
        }
        ambient_anchors.append(anchor)
        anchor_rows.append({"niemeier": entry["label"], **anchor})
    ambient_rows.append(
        {
            "niemeier": entry["label"],
            "root_component_types": [list(value) for value in component_types],
            "lifted_type_preserving_component_permutations": [
                list(permutation) for permutation in permutations
            ],
            "D5_anchor_orbit_count": len(orbits),
            "D5_anchor_orbits": ambient_anchors,
        }
    )

assert len(ambient_rows) == 13
assert len(anchor_rows) == 16

result = {
    "schema": "elkies-k3.niemeier-d5-anchor-orbits.v1",
    "status": "PASS_EXACT_D5_ANCHOR_ORBITS_NOT_FULL_AUXILIARY_EMBEDDINGS",
    "classification_scope": {
        "proved": (
            "The thirteen D5-admissible rooted Niemeier lattices contain "
            "sixteen D5 root-sublattice anchor orbits under their full "
            "automorphism groups. An explicit D5 basis is retained for each."
        ),
        "not_proved": (
            "The norm-12 and norm-24 auxiliary generators extending each D5 "
            "anchor have not yet been enumerated modulo the anchor stabilizer."
        ),
    },
    "accounting": {
        "D5_admissible_niemeier_classes": len(ambient_rows),
        "D5_anchor_orbits": len(anchor_rows),
        "excluded_niemeier_classes_from_catalog_gate": 11,
    },
    "internal_D5_weyl_orbits": {
        "D_family": (
            "One orbit in every D_n for n>=5 by the signed-coordinate "
            "description: a D5 subsystem selects five coordinate axes and "
            "the Weyl group is transitive on such selections."
        ),
        "exceptional_enumeration": exceptional,
    },
    "ambients": ambient_rows,
    "anchors": anchor_rows,
}

payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == payload
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(payload)

print(
    "NIEMEIERD5|ambients={}|anchors={}|exceptional=1,1,1|"
    "extensions_complete=0|status=PASS".format(len(ambient_rows), len(anchor_rows))
)
