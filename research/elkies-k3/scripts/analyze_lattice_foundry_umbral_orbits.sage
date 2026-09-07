#!/usr/bin/env sage-python
"""Pilot umbral-orbit analysis for selected lattice-foundry complements.

Status: exact finite computation, except for the explicitly sampled degree-three data.
Claim: recovers stabilizer images and orbit-resolved short-vector/coset counts.
Inputs: rooted Niemeier catalogue, foundry catalogue, and rootless-J2 catalogue.
Output: artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json.
Supersedes: no earlier computation.

The ambient is the Niemeier lattice with root system ``2A7+2D5``.  Its
umbral group is realized by the canonical chamber-preserving section
``Aut(N,Phi)`` used in Cheng--Duncan--Harvey.  For a stored auxiliary
embedding ``K -> N``, this script computes the *literal section stabilizer*

    {g in Aut(N,Phi) : g(K) = K},

and also the full ambient stabilizer by testing compatible
``Aut(K) x Aut(M)`` pairs against the primitive Niemeier gluing.  The latter's
umbral image and induced action on ``M=K^perp`` are used for orbit/fixed-point
data on norm-four vectors, rational-bisection cosets in ``M/2M``, and a
deterministic group-invariant sample in ``M/3M``.

The literal stabilizer is deliberately retained as a diagnostic because it
can be narrower than the image in ``G^X`` of the full ambient stabilizer: an
outer class may stabilize ``K`` only after multiplication by a Weyl element.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import itertools
import json
import random
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from sage.all import (
    Graph,
    QQ,
    ZZ,
    block_diagonal_matrix,
    ceil,
    identity_matrix,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
ROOTLESS = ROOT / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json"

SELECTED = (
    ("published_R17", "NS0001-F001"),
    ("alternate_Q80", "NS0001-F002"),
    ("NS0024", "NS0024-F005"),
    ("NS0032-F011", "NS0032-F011"),
    ("NS0028-F005", "NS0028-F005"),
    ("NS0033-F026", "NS0033-F026"),
)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def matrix_key(value):
    return tuple(map(int, value.list()))


def signed_roots_and_components(gram):
    result = pari(gram).qfminim(2)
    representatives = matrix(ZZ, result[2].sage()).transpose()
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
        signed_count = len(component_roots)
        if signed_count == rank * (rank + 1):
            family = "A"
        elif signed_count == 2 * rank * (rank - 1):
            family = "D"
        else:
            raise AssertionError((rank, signed_count))
        components.append({"roots": component_roots, "type": (family, int(rank))})
    return components


def simple_roots(gram, roots):
    for trial in range(1, 100):
        chamber = vector(
            ZZ,
            [(index + 1) ** 2 + trial * (index + 1) + trial**2 for index in range(24)],
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
            if not any(tuple(map(int, root - other)) in positive_set for other in positive)
        ],
    )
    assert basis.rank() == basis.nrows()
    return basis


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(basis.nrows()))
    graph.add_edges(
        (row, column)
        for row in range(basis.nrows())
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def graph_isometry_variants(left_basis, right_basis, gram):
    left = dynkin_graph(left_basis, gram)
    right = dynkin_graph(right_basis, gram)
    isomorphic, initial = left.is_isomorphic(right, certificate=True)
    assert isomorphic
    return [
        [automorphism(initial[index]) for index in range(left_basis.nrows())]
        for automorphism in right.automorphism_group()
    ]


def multiplicative_order(value):
    identity = identity_matrix(ZZ, value.nrows())
    power = identity
    for order in range(1, 17):
        power *= value
        if power == identity:
            return order
    raise AssertionError("unexpected automorphism order")


def umbral_class(component_types, permutation, automorphism):
    fixed = {
        family: sum(
            permutation[index] == index
            for index, (component_family, unused_rank) in enumerate(component_types)
            if component_family == family
        )
        for family in ("A", "D")
    }
    order = multiplicative_order(automorphism)
    if order == 1:
        label = "1A"
    elif order == 4:
        label = "4A"
    elif fixed == {"A": 2, "D": 2}:
        label = "2A"
    elif fixed == {"A": 0, "D": 2}:
        label = "2B"
    elif fixed == {"A": 2, "D": 0}:
        label = "2C"
    else:
        raise AssertionError((order, fixed, permutation))
    return label, order, fixed


def umbral_group_section(gram):
    components = signed_roots_and_components(gram)
    bases = [simple_roots(gram, component["roots"]) for component in components]
    component_types = [component["type"] for component in components]
    root_basis = matrix(ZZ, [row for basis in bases for row in basis.rows()])
    assert abs(root_basis.det()) > 0
    groups = {}
    for index, component_type in enumerate(component_types):
        groups.setdefault(component_type, []).append(index)

    found = {}
    type_permutations = itertools.product(
        *[list(itertools.permutations(indices)) for indices in groups.values()]
    )
    for choices in type_permutations:
        permutation = list(range(len(bases)))
        for indices, images in zip(groups.values(), choices):
            for source, target in zip(indices, images):
                permutation[source] = target
        variants = [
            graph_isometry_variants(bases[source], bases[permutation[source]], gram)
            for source in range(len(bases))
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
            key = matrix_key(automorphism)
            label, order, fixed = umbral_class(
                component_types, tuple(permutation), automorphism
            )
            found[key] = {
                "matrix": automorphism,
                "class": label,
                "order": order,
                "fixed_components": fixed,
                "component_permutation": list(map(int, permutation)),
            }
    result = sorted(found.values(), key=lambda row: (row["class"], matrix_key(row["matrix"])))
    assert len(result) == 8
    assert Counter(row["class"] for row in result) == Counter(
        {"1A": 1, "2A": 1, "2B": 2, "2C": 2, "4A": 2}
    )
    return result, bases


def integral_automorphism_group(gram):
    """Enumerate a positive-definite lattice automorphism group in row convention."""

    data = pari(gram).qfauto()
    claimed_order = int(data[0])
    generators = [matrix(ZZ, item).transpose() for item in data[1]]
    identity = identity_matrix(ZZ, gram.nrows())
    found = {matrix_key(identity): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator
            key = matrix_key(candidate)
            if key not in found:
                assert candidate * gram * candidate.transpose() == gram
                found[key] = candidate
                frontier.append(candidate)
    assert len(found) == claimed_order
    return list(found.values())


def is_weyl_element(automorphism, component_bases, ambient_gram):
    """Test the Weyl kernel via its trivial action on component discriminants."""

    for basis in component_bases:
        try:
            action = basis.transpose().solve_right(
                (basis * automorphism).transpose()
            ).transpose()
        except ValueError:
            return False
        if not all(entry.denominator() == 1 for entry in action.list()):
            return False
        action = matrix(ZZ, action)
        cartan = basis * ambient_gram * basis.transpose()
        if action * cartan * action.transpose() != cartan:
            return False
        discriminant_difference = cartan.inverse() * action - cartan.inverse()
        if not all(entry.denominator() == 1 for entry in discriminant_difference.list()):
            return False
    return True


def outer_group_element(ambient_automorphism, group_rows, component_bases, ambient_gram):
    matches = []
    for index, row in enumerate(group_rows):
        quotient = ambient_automorphism * row["matrix"].inverse()
        if is_weyl_element(quotient, component_bases, ambient_gram):
            matches.append(index)
    assert len(matches) == 1
    return matches[0]


def full_ambient_stabilizer(auxiliary, complement, ambient_gram, group_rows, component_bases):
    """Enumerate compatible Aut(K) x Aut(M) pairs extending to the Niemeier lattice."""

    auxiliary_gram = auxiliary * ambient_gram * auxiliary.transpose()
    complement_gram = complement * ambient_gram * complement.transpose()
    auxiliary_automorphisms = integral_automorphism_group(auxiliary_gram)
    complement_automorphisms = integral_automorphism_group(complement_gram)
    sublattice_basis = auxiliary.stack(complement)
    sublattice_inverse = sublattice_basis.inverse()
    extended_count = 0
    outer_indices = set()
    action_outer_indices = {}
    for auxiliary_action in auxiliary_automorphisms:
        for complement_action in complement_automorphisms:
            block_action = block_diagonal_matrix(auxiliary_action, complement_action)
            ambient_action = sublattice_inverse * block_action * sublattice_basis
            if not all(entry.denominator() == 1 for entry in ambient_action.list()):
                continue
            ambient_action = matrix(ZZ, ambient_action)
            assert abs(ambient_action.det()) == 1
            assert ambient_action * ambient_gram * ambient_action.transpose() == ambient_gram
            extended_count += 1
            outer_index = outer_group_element(
                ambient_action, group_rows, component_bases, ambient_gram
            )
            outer_indices.add(outer_index)
            action_outer_indices.setdefault(matrix_key(complement_action), set()).add(outer_index)
    actions = {
        matrix_key(action): action
        for action in complement_automorphisms
        if matrix_key(action) in action_outer_indices
    }
    assert extended_count > 0
    assert matrix_key(identity_matrix(ZZ, 17)) in actions
    return {
        "ambient_order": extended_count,
        "auxiliary_automorphism_order": len(auxiliary_automorphisms),
        "complement_automorphism_order": len(complement_automorphisms),
        "outer_indices": sorted(outer_indices),
        "actions": [actions[key] for key in sorted(actions)],
        "action_outer_indices": {
            key: sorted(action_outer_indices[key]) for key in sorted(actions)
        },
    }


def parity_mask(value) -> int:
    return sum((int(entry) % 2) << index for index, entry in enumerate(value))


def mask_vector(mask, dimension):
    return vector(ZZ, [(mask >> index) & 1 for index in range(dimension)])


def exact_rational_bisection_cosets(gram):
    change = gram.LLL_gram().transpose()
    reduced = change * gram * change.transpose()
    assert abs(change.det()) == 1
    result = pari(reduced).qfminim(10)
    masks_by_norm = {norm: set() for norm in range(2, 11, 2)}
    signed_counts = Counter()
    for column in matrix(ZZ, result[2].sage()).columns():
        value = vector(ZZ, column)
        norm = int(value * reduced * value)
        signed_counts[norm] += 2
        masks_by_norm[norm].add(parity_mask(value * change))
    rational = masks_by_norm[10] - masks_by_norm[6] - masks_by_norm[2]
    return rational, {str(norm): signed_counts[norm] for norm in range(2, 11, 2)}


def ldl_data(gram):
    dimension = gram.nrows()
    lower = matrix(QQ, dimension, dimension)
    diagonal = []
    for index in range(dimension):
        lower[index, index] = 1
        value = QQ(gram[index, index]) - sum(
            lower[index, prior] ** 2 * diagonal[prior] for prior in range(index)
        )
        assert value > 0
        diagonal.append(value)
        for row in range(index + 1, dimension):
            lower[row, index] = (
                QQ(gram[row, index])
                - sum(
                    lower[row, prior] * lower[index, prior] * diagonal[prior]
                    for prior in range(index)
                )
            ) / value
    return lower, diagonal


def minimum_coset_norm_through_bound(gram, lower, diagonal, degree, residue, bound):
    dimension = gram.nrows()
    shift = [QQ(entry) / degree for entry in residue]
    coordinates = [ZZ.zero()] * dimension
    best = None

    def visit(index, used):
        nonlocal best
        if best is not None and used > QQ(best) / (degree * degree):
            return
        if index < 0:
            scaled = vector(
                ZZ, [int(residue[i] + degree * coordinates[i]) for i in range(dimension)]
            )
            exact = int(scaled * gram * scaled)
            assert exact <= bound
            if best is None or exact < best:
                best = exact
            return
        remaining = QQ(bound) / (degree * degree) - used
        if remaining < 0:
            return
        center = shift[index] + sum(
            lower[row, index] * (coordinates[row] + shift[row])
            for row in range(index + 1, dimension)
        )
        radius_squared = remaining / diagonal[index]
        floor_radius = ZZ(radius_squared.floor()).isqrt()
        radius_ceiling = floor_radius + (QQ(floor_radius**2) < radius_squared)
        coordinate_bound = int(radius_ceiling + ceil(abs(center)))
        for entry in range(-coordinate_bound, coordinate_bound + 1):
            contribution = diagonal[index] * (QQ(entry) + center) ** 2
            if contribution <= remaining:
                coordinates[index] = ZZ(entry)
                visit(index - 1, used + contribution)

    visit(dimension - 1, QQ.zero())
    return best


def tuple_mod(value, modulus):
    return tuple(int(entry) % modulus for entry in value)


def orbit_histogram(points, transforms, image):
    unseen = set(points)
    histogram = Counter()
    while unseen:
        seed = min(unseen)
        orbit = {image(seed, transform) for transform in transforms}
        assert orbit <= set(points)
        unseen.difference_update(orbit)
        histogram[len(orbit)] += 1
    return {str(size): count for size, count in sorted(histogram.items())}


def fixed_counts(points, transforms, image):
    return [sum(image(point, transform) == point for point in points) for transform in transforms]


def induced_complement_action(complement, ambient_automorphism):
    coefficients = complement.transpose().solve_right(
        (complement * ambient_automorphism).transpose()
    ).transpose()
    assert all(entry.denominator() == 1 for entry in coefficients.list())
    coefficients = matrix(ZZ, coefficients)
    assert coefficients * complement == complement * ambient_automorphism
    return coefficients


def norm_four_vectors(gram):
    result = pari(gram).qfminim(4)
    vectors = set()
    for column in matrix(ZZ, result[2].sage()).columns():
        value = vector(ZZ, column)
        if int(value * gram * value) == 4:
            vectors.add(tuple(map(int, value)))
            vectors.add(tuple(map(int, -value)))
    return vectors


def selected_embeddings(foundry, rootless):
    frames = {
        frame["frame_id"]: frame
        for ns_row in foundry["ns_classes"]
        for frame in ns_row["frames"]
    }
    published = next(row for row in rootless["rootless_classes"] if row["matches_published_R17"])
    alternate = next(row for row in rootless["rootless_classes"] if row["matches_alternate_Q80"])
    controls = {
        "NS0001-F001": published["representative_embedding"],
        "NS0001-F002": alternate["representative_embedding"],
    }
    result = []
    for display_name, frame_id in SELECTED:
        frame = frames[frame_id]
        if frame_id in controls:
            embedding = controls[frame_id]
            provenance = "complete rootless-J2 representative embedding"
        else:
            assert frame["embedding_count_in_declared_shell"] == 1
            embedding = frame["embeddings"][0]
            provenance = "unique embedding in one-root-control-shell-v1"
        result.append((display_name, frame_id, frame, embedding, provenance))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--frame-id", action="append", default=[])
    parser.add_argument("--d3-orbit-seeds", type=int, default=512)
    parser.add_argument("--pari-stack-gb", type=int, default=4)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--isolated-child", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()

    requested_frames = arguments.frame_id or [frame_id for unused_name, frame_id in SELECTED]
    if len(requested_frames) > 1 and not arguments.isolated_child:
        with tempfile.TemporaryDirectory(prefix="umbral-orbits-") as temporary_directory:
            payloads = []
            for index, frame_id in enumerate(requested_frames):
                child_output = Path(temporary_directory) / f"target-{index:02d}.json"
                command = [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--isolated-child",
                    "--frame-id",
                    frame_id,
                    "--d3-orbit-seeds",
                    str(arguments.d3_orbit_seeds),
                    "--pari-stack-gb",
                    str(arguments.pari_stack_gb),
                    "--output",
                    str(child_output),
                ]
                subprocess.run(command, cwd=ROOT, check=True)
                payloads.append(json.loads(child_output.read_text()))
        merged = payloads[0]
        merged["targets"] = [target for payload in payloads for target in payload["targets"]]
        assert [target["frame_id"] for target in merged["targets"]] == requested_frames
        for payload in payloads[1:]:
            assert payload["inputs"] == merged["inputs"]
            assert payload["umbral_case"] == merged["umbral_case"]
            assert payload["group_section"] == merged["group_section"]
        merged["execution"] = (
            "Targets are run in isolated child processes so PARI norm-shell memory is "
            "released between frames."
        )
        merged["reproduce"] = (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/analyze_lattice_foundry_umbral_orbits.sage "
            f"--d3-orbit-seeds {arguments.d3_orbit_seeds} --pari-stack-gb "
            f"{arguments.pari_stack_gb}"
            + "".join(f" --frame-id {frame_id}" for frame_id in arguments.frame_id)
        )
        serialized = json.dumps(merged, indent=2, sort_keys=True) + "\n"
        output_path = arguments.output.resolve()
        if arguments.check:
            if output_path.read_text() != serialized:
                raise SystemExit("umbral-orbit artifact is stale")
        else:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(serialized)
        print(f"UMBRALORBIT|targets={len(merged['targets'])}|status=PASS_MERGED", flush=True)
        return

    pari.allocatemem(arguments.pari_stack_gb * 1024**3)

    catalog = json.loads(CATALOG.read_text())
    foundry = json.loads(FOUNDRY.read_text())
    rootless = json.loads(ROOTLESS.read_text())
    ambient_row = next(
        row for row in catalog["rooted_niemeier_lattices"] if row["label"] == "2A7_2D5"
    )
    ambient_gram = matrix(ZZ, ambient_row["gram"])
    group, component_bases = umbral_group_section(ambient_gram)

    selected_rows = selected_embeddings(foundry, rootless)
    if arguments.frame_id:
        wanted = set(arguments.frame_id)
        selected_rows = [row for row in selected_rows if row[1] in wanted]
        missing = wanted - {row[1] for row in selected_rows}
        if missing:
            parser.error(f"unknown pilot frame ids: {sorted(missing)}")
    targets = []
    for display_name, frame_id, frame, embedding, provenance in selected_rows:
        auxiliary = matrix(ZZ, embedding["auxiliary_basis_in_ambient"])
        complement = matrix(ZZ, embedding["complement_basis_in_ambient"])
        assert auxiliary * ambient_gram * complement.transpose() == 0
        complement_gram = complement * ambient_gram * complement.transpose()
        assert pari(complement_gram).qfisom(pari(matrix(ZZ, frame["gram"]))) != 0

        literal_stabilizer_rows = []
        literal_actions = []
        for group_row in group:
            automorphism = group_row["matrix"]
            if auxiliary.row_module() != (auxiliary * automorphism).row_module():
                continue
            action = induced_complement_action(complement, automorphism)
            assert action * complement_gram * action.transpose() == complement_gram
            literal_stabilizer_rows.append(group_row)
            literal_actions.append(action)
        assert literal_stabilizer_rows

        full_stabilizer = full_ambient_stabilizer(
            auxiliary, complement, ambient_gram, group, component_bases
        )
        actions = full_stabilizer["actions"]
        outer_rows = [group[index] for index in full_stabilizer["outer_indices"]]

        shell = norm_four_vectors(complement_gram)
        shell_image = lambda point, action: tuple(map(int, vector(ZZ, point) * action))
        shell_fixed = fixed_counts(shell, actions, shell_image)
        shell_orbits = orbit_histogram(shell, actions, shell_image)
        unoriented_shell = {min(point, tuple(-entry for entry in point)) for point in shell}
        pair_image = lambda point, action: min(
            tuple(map(int, vector(ZZ, point) * action)),
            tuple(map(int, -vector(ZZ, point) * action)),
        )
        pair_orbits = orbit_histogram(unoriented_shell, actions, pair_image)

        bisections, signed_shell_counts = exact_rational_bisection_cosets(complement_gram)
        d2_image = lambda point, action: parity_mask(mask_vector(point, 17) * action)
        d2_fixed = fixed_counts(bisections, actions, d2_image)
        d2_orbits = orbit_histogram(bisections, actions, d2_image)

        change = complement_gram.LLL_gram().transpose()
        reduced = change * complement_gram * change.transpose()
        change_inverse = matrix(ZZ, change.inverse())
        reduced_actions = [change * action * change_inverse for action in actions]
        assert all(action.change_ring(ZZ) == action for action in reduced_actions)
        lower, diagonal = ldl_data(reduced)
        seed_value = int(hashlib.sha256(frame_id.encode()).hexdigest()[:16], 16)
        rng = random.Random(20260901 ^ seed_value)
        d3_sample = {(0,) * 17}
        while len(d3_sample) < arguments.d3_orbit_seeds:
            residue = tuple(rng.randrange(3) for unused in range(17))
            d3_sample.update(
                tuple_mod(vector(ZZ, residue) * action, 3) for action in reduced_actions
            )
        d3_minima = {}
        d3_rational = set()
        for residue in sorted(d3_sample):
            residue_vector = vector(ZZ, residue)
            minimum = minimum_coset_norm_through_bound(
                reduced, lower, diagonal, 3, residue_vector, 24
            )
            d3_minima[residue] = minimum
            residue_norm = int(residue_vector * reduced * residue_vector)
            if minimum is not None and minimum >= 20 and (residue_norm - 2) % 6 == 0:
                d3_rational.add(residue)
        d3_image = lambda point, action: tuple_mod(vector(ZZ, point) * action, 3)
        d3_fixed = fixed_counts(d3_rational, reduced_actions, d3_image)
        d3_orbits = orbit_histogram(d3_rational, reduced_actions, d3_image)
        minima_histogram = Counter(
            "above_24" if value is None else str(value) for value in d3_minima.values()
        )

        targets.append(
            {
                "display_name": display_name,
                "frame_id": frame_id,
                "determinant": int(complement_gram.det()),
                "niemeier_root_system": "A7^2 D5^2",
                "embedding_provenance": provenance,
                "stored_ambient_label": embedding.get("niemeier", embedding.get("ambient")),
                "literal_section_stabilizer": {
                    "order": len(literal_actions),
                    "class_multiset": dict(
                        sorted(Counter(row["class"] for row in literal_stabilizer_rows).items())
                    ),
                    "classes": [row["class"] for row in literal_stabilizer_rows],
                    "induced_action_matrices": [rows(action) for action in literal_actions],
                },
                "full_ambient_stabilizer": {
                    "order": full_stabilizer["ambient_order"],
                    "auxiliary_automorphism_group_order": full_stabilizer[
                        "auxiliary_automorphism_order"
                    ],
                    "complement_automorphism_group_order": full_stabilizer[
                        "complement_automorphism_order"
                    ],
                    "induced_complement_action_order": len(actions),
                    "umbral_image_order": len(full_stabilizer["outer_indices"]),
                    "umbral_image_class_multiset": dict(
                        sorted(Counter(row["class"] for row in outer_rows).items())
                    ),
                    "umbral_image_classes": [row["class"] for row in outer_rows],
                    "induced_actions": [
                        {
                            "matrix": rows(action),
                            "compatible_umbral_classes": sorted(
                                {
                                    group[index]["class"]
                                    for index in full_stabilizer["action_outer_indices"][
                                        matrix_key(action)
                                    ]
                                }
                            ),
                            "fixed_norm_four_signed_vectors": shell_fixed[action_index],
                            "fixed_degree_two_rational_cosets": d2_fixed[action_index],
                            "fixed_sampled_degree_three_rational_cosets": d3_fixed[
                                action_index
                            ],
                        }
                        for action_index, action in enumerate(actions)
                    ],
                },
                "norm_four": {
                    "signed_vectors": len(shell),
                    "unoriented_pairs": len(unoriented_shell),
                    "signed_orbit_size_histogram": shell_orbits,
                    "unoriented_orbit_size_histogram": pair_orbits,
                },
                "degree_two_rational_cosets": {
                    "count": len(bisections),
                    "orbit_size_histogram": d2_orbits,
                    "signed_vector_shell_counts_through_norm_ten": signed_shell_counts,
                },
                "degree_three_rational_sample": {
                    "group_invariant_sample_size": len(d3_sample),
                    "orbit_seed_floor": arguments.d3_orbit_seeds,
                    "qualifying_count": len(d3_rational),
                    "minimum_norm_histogram": dict(sorted(minima_histogram.items())),
                    "orbit_size_histogram": d3_orbits,
                    "status": "EXACT_CVP_INSIDE_DETERMINISTIC_GROUP_INVARIANT_SAMPLE",
                },
            }
        )
        print(
            f"UMBRALORBIT|frame={frame_id}|literal_stab={len(literal_actions)}|"
            f"ambient_stab={full_stabilizer['ambient_order']}|"
            f"umbral_image={len(full_stabilizer['outer_indices'])}|M_action={len(actions)}|"
            f"norm4={len(shell)}|d2={len(bisections)}|d3_sample={len(d3_sample)}|"
            f"d3_rational={len(d3_rational)}|status=PASS",
            flush=True,
        )
        del shell, bisections, d3_minima, d3_sample, d3_rational
        gc.collect()

    output = {
        "schema": "elkies-k3.lattice-foundry-umbral-orbits.v1",
        "status": "PASS_EXACT_AMBIENT_STABILIZERS_D2_ORBITS_AND_SAMPLED_D3_ORBITS",
        "scope": {
            "proved": (
                "All six stored embeddings lie in N(A7^2 D5^2). The eight-element "
                "chamber-preserving section is enumerated exactly and has Dih_4 class "
                "distribution 1A,2A,2B,2C,4A. Full ambient stabilizers are recovered "
                "from compatible Aut(K) x Aut(M) pairs across the primitive gluing. "
                "Their umbral images, induced complement actions, norm-four orbits, "
                "and all rational bisection-coset orbits are exact."
            ),
            "not_proved": (
                "Degree-three data are sampled. No equality or graded correspondence "
                "with an umbral module is claimed, and no geometric or arithmetic "
                "consequence is inferred from the sampled degree-three cosets."
            ),
        },
        "umbral_case": {
            "X": "A7^2 D5^2",
            "coxeter_number": 8,
            "lambency": 8,
            "G_X": "Dih_4 (dihedral group of order 8)",
            "G_X_order": 8,
            "class_distribution": dict(sorted(Counter(row["class"] for row in group).items())),
            "reference": "Cheng--Duncan--Harvey, Umbral Moonshine and the Niemeier Lattices, Tables 2, 18, 38--44",
            "low_grade_trace_controls": {
                "H_g_1_discriminant_31": {"1A": 2, "2A": 2, "2BC": -2, "4A": 2},
                "H_g_2_discriminant_28": {"1A": 4, "2A": -4, "2BC": 0, "4A": 0},
            },
        },
        "twining_convention": (
            "A group element permutes a finite vector/coset shell. Its permutation "
            "character is the number of fixed shell elements. The expression tr(g|v) "
            "for an individual lattice vector is not canonical; fixed-point theta "
            "coefficients are the concrete convention used here."
        ),
        "group_section": [
            {
                "class": row["class"],
                "order": row["order"],
                "fixed_components": row["fixed_components"],
                "component_permutation": row["component_permutation"],
                "ambient_matrix": rows(row["matrix"]),
            }
            for row in group
        ],
        "targets": targets,
        "inputs": {
            relative(CATALOG): digest(CATALOG),
            relative(FOUNDRY): digest(FOUNDRY),
            relative(ROOTLESS): digest(ROOTLESS),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/analyze_lattice_foundry_umbral_orbits.sage "
            f"--d3-orbit-seeds {arguments.d3_orbit_seeds}"
            + "".join(f" --frame-id {frame_id}" for frame_id in arguments.frame_id)
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("umbral-orbit artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(f"UMBRALORBIT|targets={len(targets)}|status=PASS", flush=True)


if __name__ == "__main__":
    main()
