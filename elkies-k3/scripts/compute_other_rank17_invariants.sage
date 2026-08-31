#!/usr/bin/env sage
"""Compute exact intrinsic invariants of the alternate rootless rank-17 frame.

The theta expansion is deliberately bounded and labelled as such.  Every
reported coefficient is exact, but no claim is made that four coefficients
determine the full theta series.
"""

import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, identity_matrix, lcm, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-other-rank17-invariants.json"


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def canonical_line(value):
    result = vector(ZZ, value)
    first = next(entry for entry in result if entry)
    if first < 0:
        result = -result
    return tuple(map(int, result))


def matrix_group(generators, form):
    identity = identity_matrix(ZZ, form.nrows())
    found = {tuple(identity.list()): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator
            key = tuple(candidate.list())
            if key in found:
                continue
            assert candidate.transpose() * form * candidate == form
            found[key] = candidate
            frontier.append(candidate)
    return tuple(found.values())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

source = json.loads(SOURCE.read_text())
gram = matrix(ZZ, source["rootless_frame"])
assert gram.nrows() == gram.ncols() == 17
assert gram.is_positive_definite() and gram.det() == 948

reduction = gram.LLL_gram()
reduced = reduction.transpose() * gram * reduction
assert abs(reduction.det()) == 1
assert reduced.det() == 948
assert set(reduced.diagonal()) == {4}

# PARI returns one representative from each sign pair in the third component.
minimum_data = pari(gram).qfminim(4)
assert int(minimum_data[0]) == 2626
minimum_columns = matrix(ZZ, minimum_data[2])
minimal_lines = sorted(
    canonical_line(minimum_columns.column(index))
    for index in range(minimum_columns.ncols())
)
assert len(minimal_lines) == len(set(minimal_lines)) == 1313
minimal_matrix = matrix(ZZ, minimal_lines)
pairings = minimal_matrix * gram * minimal_matrix.transpose()

# The minimal-vector graph uses the standard additive relation: two norm-four
# lines are adjacent exactly when their absolute pairing is two, equivalently
# when one choice of signs has a norm-four sum or difference.
adjacency = [[] for _ in minimal_lines]
pairing_distribution = Counter()
adjacency_hasher = hashlib.sha256()
for left in range(len(minimal_lines)):
    for right in range(left + 1, len(minimal_lines)):
        absolute_pairing = abs(int(pairings[left, right]))
        pairing_distribution[absolute_pairing] += 1
        adjacency_hasher.update(bytes((absolute_pairing,)))
        if absolute_pairing == 2:
            adjacency[left].append(right)
            adjacency[right].append(left)

degree_distribution = Counter(map(len, adjacency))
unseen = set(range(len(minimal_lines)))
component_sizes = []
while unseen:
    start = min(unseen)
    unseen.remove(start)
    queue = deque([start])
    size = 0
    while queue:
        current = queue.popleft()
        size += 1
        for following in adjacency[current]:
            if following in unseen:
                unseen.remove(following)
                queue.append(following)
    component_sizes.append(size)

automorphism_data = pari(gram).qfauto()
automorphism_order = int(automorphism_data[0])
automorphism_generators = tuple(
    matrix(ZZ, generator) for generator in automorphism_data[1]
)
automorphisms = matrix_group(automorphism_generators, gram)
assert len(automorphisms) == automorphism_order == 4

line_index = {line: index for index, line in enumerate(minimal_lines)}
line_permutations = []
for automorphism in automorphisms:
    permutation = tuple(
        line_index[canonical_line(automorphism * vector(ZZ, line))]
        for line in minimal_lines
    )
    line_permutations.append(permutation)
induced_line_group_order = len(set(line_permutations))
line_orbits = []
unseen = set(range(len(minimal_lines)))
while unseen:
    start = min(unseen)
    orbit = {permutation[start] for permutation in line_permutations}
    unseen.difference_update(orbit)
    line_orbits.append(len(orbit))

smith_diagonal = list(map(int, gram.smith_form()[0].diagonal()))
assert smith_diagonal == [1] * 16 + [948]
inverse_gram = gram.inverse()
discriminant_generator_column = next(
    index
    for index in range(gram.nrows())
    if lcm(entry.denominator() for entry in inverse_gram.column(index)) == 948
)
discriminant_generator = inverse_gram.column(discriminant_generator_column)
discriminant_quadratic_value = (
    discriminant_generator * gram * discriminant_generator
)
discriminant_quadratic_value_mod_two = (
    discriminant_quadratic_value.numerator()
    % (2 * discriminant_quadratic_value.denominator())
) / discriminant_quadratic_value.denominator()
discriminant_actions = []
for generator in automorphism_generators:
    image = generator * discriminant_generator
    multipliers = [
        multiplier
        for multiplier in range(948)
        if all(
            entry.denominator() == 1
            for entry in image - multiplier * discriminant_generator
        )
    ]
    assert len(multipliers) == 1
    assert (
        (multipliers[0] ** 2 - 1)
        * discriminant_quadratic_value_mod_two
        / 2
    ).denominator() == 1
    discriminant_actions.append(multipliers[0])
assert set(discriminant_actions) == {473, 947}

# Exact cumulative qfminim counts, converted to individual theta
# coefficients.  The convention is theta_L(q)=sum_v q^(v.v/2).
cumulative = {0: 1}
for norm in (2, 4, 6, 8):
    cumulative[norm] = 1 + int(pari(gram).qfminim(norm)[0])
theta_by_norm = {0: 1}
previous = 1
for norm in (2, 4, 6, 8):
    theta_by_norm[norm] = cumulative[norm] - previous
    previous = cumulative[norm]
assert theta_by_norm == {0: 1, 2: 0, 4: 2626, 6: 53290, 8: 460360}

payload = {
    "schema": "elkies-k3.other-rank17-invariants.v1",
    "status": "PASS_EXACT_ALTERNATE_RANK17_INTRINSIC_INVARIANTS",
    "frame": {
        "rank": 17,
        "determinant": 948,
        "minimum_squared_norm": 4,
        "norm_two_root_count": 0,
        "smith_invariants": smith_diagonal,
        "reduced_gram": rows(reduced),
        "reduced_basis_columns_in_source_basis": rows(reduction),
    },
    "theta_series": {
        "convention": "theta_L(q)=sum_{v in L} q^(v.v/2)",
        "exact_through_squared_norm": 8,
        "coefficients_by_squared_norm": {
            str(norm): count for norm, count in theta_by_norm.items()
        },
        "expansion": "1 + 2626*q^2 + 53290*q^3 + 460360*q^4 + O(q^5)",
        "proof_boundary": (
            "The displayed truncation is exact. It is not asserted to determine "
            "the infinite theta series."
        ),
    },
    "minimal_vector_graph": {
        "vertices_unoriented_norm_four_lines": len(minimal_lines),
        "edge_rule": "absolute pairing equals 2",
        "edges": sum(map(len, adjacency)) // 2,
        "connected_component_sizes": sorted(component_sizes, reverse=True),
        "degree_distribution": {
            str(degree): count for degree, count in sorted(degree_distribution.items())
        },
        "absolute_pairing_distribution_on_distinct_lines": {
            str(pairing): count
            for pairing, count in sorted(pairing_distribution.items())
        },
        "lexicographic_pair_label_sha256": adjacency_hasher.hexdigest(),
        "lattice_automorphism_orbit_sizes_on_lines": sorted(line_orbits),
    },
    "automorphism_group": {
        "order": automorphism_order,
        "generator_count": len(automorphism_generators),
        "generators_column_convention": [rows(item) for item in automorphism_generators],
        "induced_group_order_on_unoriented_minimal_lines": induced_line_group_order,
    },
    "discriminant_form_action": {
        "group": "Z/948Z",
        "dual_generator_source_basis_column": discriminant_generator_column,
        "dual_generator": [str(entry) for entry in discriminant_generator],
        "quadratic_value_of_generator_mod_2Z": str(
            discriminant_quadratic_value_mod_two
        ),
        "automorphism_generator_multipliers_mod_948": discriminant_actions,
        "full_image_multipliers_mod_948": sorted(
            {1, *discriminant_actions, (discriminant_actions[0] * discriminant_actions[1]) % 948}
        ),
        "image_order": 4,
    },
    "inputs": {
        str(SOURCE.relative_to(ROOT)): hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "rootless_frame_sha256": source["rootless_frame_sha256"],
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compute_other_rank17_invariants.sage"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != serialized:
        raise SystemExit("alternate rank-17 invariant artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(serialized)

print(
    "OTHERR17INV|theta_norm8=exact|minimal_lines=1313|graph_edges={}|"
    "aut_order=4|disc_image=1,473,475,947|status=PASS".format(
        payload["minimal_vector_graph"]["edges"]
    ),
    flush=True,
)
