#!/usr/bin/env sage-python
"""Certify the prime-local bridge-mutation normal form and bad-prime census.

For a rank-four relative position of two primitive hyperbolic planes, the
script does three small exact tasks.

* It dissects the seven stored edges with gcd(det(C), |det(NS)|) > 1.
* It enumerates every marked maximal graph at the shared prime, transports
  that graph across the rank-four bridge, and records the exact ADE change.
* It repeats the enumeration for the non-cyclic R17 Z/4 + Z/8 control.

The prime-to-shared-prime part of each historical glue is held at its stored
value.  This is the finite experiment required by Theorem H-1d; it is not an
enumeration of unrelated global changes to the common core or good-prime
decorations.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import floor
from pathlib import Path

from sage.env import SAGE_VERSION
from sage.all import (
    FreeQuadraticModule,
    Genus,
    QQ,
    QuadraticForm,
    ZZ,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    matrix,
    pari,
    vector,
)
from sage.modules.torsion_quadratic_module import TorsionQuadraticModule


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-prime-local-bridge-mutation-v1.json"
)
BRIDGE_CORPUS = Path(
    "artifacts/generated-results/"
    "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
)
RELATIVE_CORPUS = Path(
    "artifacts/generated-results/"
    "elkies-k3-relative-u-bridge-lifting-regression-v1.json"
)
R17_CERTIFICATE = Path(
    "artifacts/generated-results/"
    "elkies-k3-r17-local-bridge-mutation-v1.json"
)
R17_SHORT_GRAM = Path("elkies-k3/data/lattice/short_vector_basis_gram.txt")
HYPERBOLIC_PLANE = matrix(ZZ, [[0, 1], [1, 0]])


def load_json(relative_path):
    return json.loads((ROOT / relative_path).read_text())


def load_matrix(relative_path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in (ROOT / relative_path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def sha256(relative_path):
    return hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def discriminant_module(gram):
    ambient = FreeQuadraticModule(ZZ, gram.nrows(), gram)
    dual = ambient.span_of_basis(gram.inverse().rows())
    return TorsionQuadraticModule(dual, ambient)


def finite_form_key(form):
    normal = form.normal_form()
    return (
        tuple(map(int, normal.invariants())),
        tuple(
            tuple(str(entry) for entry in row)
            for row in normal.gram_matrix_quadratic().rows()
        ),
    )


def finite_form_record(form):
    normal = form.normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rational_rows(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def local_form_record(gram, prime, sign=1):
    form = Genus(gram).discriminant_form().primary_part(prime)
    if sign == -1:
        form = form.twist(-1)
    return finite_form_record(form)


def fractional_part(value):
    value = QQ(value)
    return value - floor(value)


def element_key(value):
    return tuple(fractional_part(entry) for entry in value)


def cyclic_subgroup_key(generator, order):
    return tuple(sorted(element_key(index * generator) for index in range(order)))


def subgroup_key(generators, orders):
    assert generators and len(generators) == len(orders)
    values = {tuple(QQ(0) for _ in range(len(generators[0])))}
    for generator, order in zip(generators, orders):
        values = {
            element_key(vector(QQ, value) + index * generator)
            for value in values
            for index in range(order)
        }
    return tuple(sorted(values))


def even_overlattice_gram(split_gram, glue_lifts):
    ambient = FreeQuadraticModule(ZZ, split_gram.nrows(), split_gram)
    lattice = ambient.span(
        list(identity_matrix(ZZ, split_gram.nrows()).rows()) + glue_lifts
    )
    basis = lattice.basis_matrix()
    gram = basis * split_gram * basis.transpose()
    assert gram in matrix(ZZ, gram.nrows(), gram.ncols()).parent()
    gram = gram.change_ring(ZZ)
    assert all(entry % 2 == 0 for entry in gram.diagonal())
    return gram


def signed_roots(gram):
    result = pari(gram).qfminim(2)
    positive = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2].sage()).columns()
    ]
    answer = positive + [-root for root in positive]
    assert len(answer) == int(result[0])
    return answer


def root_system_type(gram):
    roots = signed_roots(gram)
    seen = set()
    components = []
    for start in range(len(roots)):
        if start in seen:
            continue
        seen.add(start)
        frontier = [start]
        indices = []
        while frontier:
            left = frontier.pop()
            indices.append(left)
            for right in range(len(roots)):
                if (
                    right not in seen
                    and roots[left] * gram * roots[right] != 0
                ):
                    seen.add(right)
                    frontier.append(right)
        component_rank = matrix(
            ZZ, [roots[index] for index in indices]
        ).rank()
        root_count = len(indices)
        if root_count == component_rank * (component_rank + 1):
            name = f"A{component_rank}"
        elif (
            component_rank >= 4
            and root_count == 2 * component_rank * (component_rank - 1)
        ):
            name = f"D{component_rank}"
        elif (component_rank, root_count) in ((6, 72), (7, 126), (8, 240)):
            name = f"E{component_rank}"
        else:
            raise ArithmeticError(
                f"unrecognized simply-laced component ({component_rank}, {root_count})"
            )
        components.append((component_rank, name))

    multiplicities = Counter(name for _, name in components)
    ordered_names = []
    for _, name in sorted(set(components)):
        multiplicity = multiplicities[name]
        ordered_names.append(f"{multiplicity if multiplicity > 1 else ''}{name}")
    return "+".join(ordered_names) if ordered_names else "0"


def root_rank(root_type):
    if root_type == "0":
        return 0
    answer = 0
    for component in root_type.split("+"):
        split = 0
        while split < len(component) and component[split].isdigit():
            split += 1
        multiplicity = int(component[:split]) if split else 1
        answer += multiplicity * int(component[split + 1 :])
    return answer


def row_isometry(source, target):
    witness = QuadraticForm(ZZ, source).is_globally_equivalent_to(
        QuadraticForm(ZZ, target), return_matrix=True
    )
    assert witness is not False
    witness = matrix(ZZ, witness)
    inverse = witness.inverse()
    candidates = [witness, witness.transpose()]
    if inverse in matrix(ZZ, inverse.nrows(), inverse.ncols()).parent():
        inverse = inverse.change_ring(ZZ)
        candidates.extend([inverse, inverse.transpose()])
    for candidate in candidates:
        if candidate * source * candidate.transpose() == target:
            assert abs(candidate.det()) == 1
            return candidate
    raise ArithmeticError("no verified row-convention integral isometry")


def row_automorphisms(gram):
    generators = [
        matrix(ZZ, item).transpose()
        for item in pari(gram).qfauto()[1].sage()
    ]
    identity = identity_matrix(ZZ, gram.nrows())
    known = {tuple(identity.list()): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator
            assert candidate * gram * candidate.transpose() == gram
            key = tuple(candidate.list())
            if key not in known:
                known[key] = candidate
                frontier.append(candidate)
    assert len(known) == int(pari(gram).qfauto()[0])
    return list(known.values())


def raw_bridge_transport(cross):
    old_raw = cross.transpose() * HYPERBOLIC_PLANE * cross - HYPERBOLIC_PLANE
    new_raw = cross * HYPERBOLIC_PLANE * cross.transpose() - HYPERBOLIC_PLANE
    new_u_in_old = (cross.transpose() * HYPERBOLIC_PLANE).augment(
        identity_matrix(ZZ, 2)
    )
    reverse_bridge_in_old = (
        identity_matrix(ZZ, 2)
        - cross
        * HYPERBOLIC_PLANE
        * cross.transpose()
        * HYPERBOLIC_PLANE
    ).augment(-cross * HYPERBOLIC_PLANE)
    change = new_u_in_old.stack(reverse_bridge_in_old)
    assert abs(change.det()) == 1
    assert (
        change
        * block_diagonal_matrix(HYPERBOLIC_PLANE, -old_raw)
        * change.transpose()
        == block_diagonal_matrix(HYPERBOLIC_PLANE, -new_raw)
    )
    return old_raw, new_raw, change


def bridge_transport_witnesses(
    cross,
    old_bridge,
    new_bridge,
    actual_old_generators,
    actual_new_key,
    actual_orders,
):
    old_raw, new_raw, change = raw_bridge_transport(cross)
    old_seed = row_isometry(old_raw, old_bridge)
    new_seed = row_isometry(new_raw, new_bridge)
    core_rank = len(actual_old_generators[0]) - 2
    witnesses = []
    for old_automorphism in row_automorphisms(old_bridge):
        old_identification = old_automorphism * old_seed
        for new_automorphism in row_automorphisms(new_bridge):
            new_identification = new_automorphism * new_seed

            def transport(
                value,
                old_identification=old_identification,
                new_identification=new_identification,
            ):
                old_coordinates = vector(
                    QQ,
                    [0, 0]
                    + list(value[-2:] * old_identification),
                )
                new_coordinates = old_coordinates * change.inverse()
                assert all(entry in ZZ for entry in new_coordinates[:2])
                new_bridge_coordinates = (
                    new_coordinates[-2:] * new_identification.inverse()
                )
                return vector(
                    QQ,
                    list(value[:core_rank]) + list(new_bridge_coordinates),
                )

            transported_actual_key = subgroup_key(
                [transport(item) for item in actual_old_generators],
                actual_orders,
            )
            if transported_actual_key == actual_new_key:
                witnesses.append(transport)
    assert witnesses
    return old_raw, new_raw, witnesses


def stored_glue_generator(frame_record):
    generators = frame_record["glue_generators"]
    assert len(generators) == 1
    generator = generators[0]
    return int(generator["order"]), vector(
        QQ, [QQ(value) for value in generator["K_plus_C_dual_coordinates"]]
    )


def graph_generators_from_split(split_coordinates):
    diagonal_matrix, _, smith_right = split_coordinates.smith_form()
    orders = [abs(int(value)) for value in diagonal_matrix.diagonal()]
    inverse_right = smith_right.inverse().change_ring(ZZ)
    inverse_split = split_coordinates.inverse()
    answer = []
    for coordinate_index, order in enumerate(orders):
        if order <= 1:
            continue
        smith_generator = vector(
            ZZ,
            identity_matrix(ZZ, split_coordinates.nrows()).row(coordinate_index),
        )
        frame_coordinates = smith_generator * inverse_right
        lift = frame_coordinates * inverse_split
        assert order * lift in ZZ ** split_coordinates.nrows()
        answer.append((order, lift))
    return answer


def transition_histogram_record(histogram):
    return [
        {
            "old_root_system": old_type,
            "new_root_system": new_type,
            "marked_graph_count": count,
        }
        for (old_type, new_type), count in sorted(histogram.items())
    ]


def historical_edge_census():
    bridge_corpus = load_json(BRIDGE_CORPUS)
    relative_corpus = load_json(RELATIVE_CORPUS)
    assert bridge_corpus["status"] == "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES"
    assert relative_corpus["status"] == "PASS_EXACT_RELATIVE_U_BRIDGE_LIFTING_REGRESSION"
    relative_by_edge = {
        (edge["corridor"], int(edge["edge_index"])): edge
        for edge in relative_corpus["edges"]
    }

    answer = []
    for edge in bridge_corpus["edges"]:
        old_frame = edge["old_frame"]
        new_frame = edge["new_frame"]
        common_bridge_determinant = int(
            old_frame["bridge_determinant_absolute"]
        )
        common_glue_order = int(old_frame["K_plus_C_index_in_W"])
        core_determinant = int(edge["core"]["determinant_absolute"])
        ambient_determinant = (
            core_determinant
            * common_bridge_determinant
            // common_glue_order**2
        )
        shared = gcd(common_bridge_determinant, ambient_determinant)
        if shared == 1:
            continue
        assert len(ZZ(shared).prime_divisors()) == 1
        prime = int(ZZ(shared).prime_divisors()[0])
        bridge_exponent = int(ZZ(common_bridge_determinant).valuation(prime))
        ambient_exponent = int(ZZ(ambient_determinant).valuation(prime))
        core_exponent = int(ZZ(core_determinant).valuation(prime))
        prime_power = prime**bridge_exponent
        support_only_exponents = list(
            range(max(0, bridge_exponent - ambient_exponent), bridge_exponent + 1)
        )
        forced_glue_exponent_numerator = (
            core_exponent + bridge_exponent - ambient_exponent
        )
        assert forced_glue_exponent_numerator % 2 == 0
        forced_glue_exponent = forced_glue_exponent_numerator // 2
        assert forced_glue_exponent == bridge_exponent

        key = (edge["corridor"], int(edge["edge_index"]))
        relative_edge = relative_by_edge[key]
        cross = matrix(
            ZZ,
            relative_edge["relative_u"]["old_to_new"]["cross_pairing_A"],
        )
        old_raw = matrix(
            ZZ,
            relative_edge["relative_u"]["old_to_new"][
                "positive_projection_gram_G_A"
            ],
        )
        new_raw = matrix(
            ZZ,
            relative_edge["relative_u"]["new_to_old"][
                "positive_projection_gram_G_A"
            ],
        )
        assert old_raw == cross.transpose() * HYPERBOLIC_PLANE * cross - HYPERBOLIC_PLANE
        assert new_raw == cross * HYPERBOLIC_PLANE * cross.transpose() - HYPERBOLIC_PLANE
        saturation_index = int(
            relative_edge["relative_u"]["old_to_new"]["saturation_index"]
        )
        assert saturation_index == int(
            relative_edge["relative_u"]["new_to_old"]["saturation_index"]
        ) == 1

        core = matrix(ZZ, edge["core"]["gram"])
        old_bridge = matrix(ZZ, old_frame["bridge_gram"])
        new_bridge = matrix(ZZ, new_frame["bridge_gram"])
        old_order, actual_old = stored_glue_generator(old_frame)
        new_order, actual_new = stored_glue_generator(new_frame)
        assert old_order == new_order == common_bridge_determinant
        actual_old_prime = (common_bridge_determinant // prime_power) * actual_old
        actual_new_prime = (common_bridge_determinant // prime_power) * actual_new
        actual_old_key = cyclic_subgroup_key(actual_old_prime, prime_power)
        actual_new_key = cyclic_subgroup_key(actual_new_prime, prime_power)

        _, _, transport_witnesses = bridge_transport_witnesses(
            cross,
            old_bridge,
            new_bridge,
            [actual_old_prime],
            actual_new_key,
            [prime_power],
        )

        core_local = discriminant_module(core).primary_part(prime)
        old_bridge_local = discriminant_module(old_bridge).primary_part(prime)
        new_bridge_local = discriminant_module(new_bridge).primary_part(prime)
        assert finite_form_key(old_bridge_local) == finite_form_key(new_bridge_local)
        assert old_bridge_local.invariants() == (prime_power,)

        graph_candidates = {}
        for core_element in core_local:
            if core_element.order() != prime_power:
                continue
            for bridge_element in old_bridge_local:
                if bridge_element.order() != prime_power:
                    continue
                if not (core_element.q() + bridge_element.q()).is_zero():
                    continue
                lift = vector(
                    QQ,
                    list(core_element.lift()) + list(bridge_element.lift()),
                )
                graph_candidates.setdefault(
                    cyclic_subgroup_key(lift, prime_power), lift
                )
        assert actual_old_key in graph_candidates

        old_split = block_diagonal_matrix(core, old_bridge)
        new_split = block_diagonal_matrix(core, new_bridge)
        stored_old_gram = even_overlattice_gram(old_split, [actual_old])
        stored_new_gram = even_overlattice_gram(new_split, [actual_new])
        assert abs(stored_old_gram.det()) == abs(stored_new_gram.det()) == ambient_determinant
        target_frame_form_key = finite_form_key(
            Genus(stored_old_gram).discriminant_form().primary_part(prime)
        )
        assert target_frame_form_key == finite_form_key(
            Genus(stored_new_gram).discriminant_form().primary_part(prime)
        )

        histogram_results = []
        actual_results = []
        q_ns_compatible_counts = []
        for transport in transport_witnesses:
            histogram = Counter()
            actual_transition = None
            q_ns_compatible = 0
            for graph_key, old_prime_lift in graph_candidates.items():
                new_prime_lift = transport(old_prime_lift)
                old_candidate = even_overlattice_gram(
                    old_split, [prime_power * actual_old, old_prime_lift]
                )
                new_candidate = even_overlattice_gram(
                    new_split, [prime_power * actual_new, new_prime_lift]
                )
                assert abs(old_candidate.det()) == abs(new_candidate.det()) == ambient_determinant
                old_form_key = finite_form_key(
                    Genus(old_candidate).discriminant_form().primary_part(prime)
                )
                new_form_key = finite_form_key(
                    Genus(new_candidate).discriminant_form().primary_part(prime)
                )
                if old_form_key == new_form_key == target_frame_form_key:
                    q_ns_compatible += 1
                else:
                    continue
                transition = (
                    root_system_type(old_candidate),
                    root_system_type(new_candidate),
                )
                histogram[transition] += 1
                if graph_key == actual_old_key:
                    actual_transition = transition
            assert actual_transition is not None
            histogram_results.append(tuple(sorted(histogram.items())))
            actual_results.append(actual_transition)
            q_ns_compatible_counts.append(q_ns_compatible)

        assert len(set(histogram_results)) == 1
        assert len(set(actual_results)) == 1
        assert len(set(q_ns_compatible_counts)) == 1
        histogram = Counter(dict(histogram_results[0]))
        actual_transition = actual_results[0]
        q_ns_compatible_count = q_ns_compatible_counts[0]
        assert q_ns_compatible_count == len(graph_candidates)
        assert histogram[actual_transition] == 1
        assert root_rank(actual_transition[0]) == int(edge["source_root_rank"])
        assert root_rank(actual_transition[1]) == int(edge["target_root_rank"])

        answer.append(
            {
                "corridor": edge["corridor"],
                "edge_index": int(edge["edge_index"]),
                "cross_pairing_A": rows(cross),
                "det_A": int(cross.det()),
                "raw_bridge_grams": {
                    "old": rows(old_raw),
                    "new": rows(new_raw),
                },
                "saturation_index": saturation_index,
                "ambient_determinant": ambient_determinant,
                "core_determinant": core_determinant,
                "common_bridge_determinant": common_bridge_determinant,
                "shared_bad_prime": prime,
                "valuations": {
                    "v_l_D": ambient_exponent,
                    "v_l_det_K": core_exponent,
                    "v_l_det_C": bridge_exponent,
                },
                "support_theorem_only_local_glue_orders": [
                    prime**exponent for exponent in support_only_exponents
                ],
                "core_and_q_NS_forced_local_glue_order": prime**forced_glue_exponent,
                "local_modules": {
                    "A_K_l_invariants": list(map(int, core_local.invariants())),
                    "A_C_old_l": finite_form_record(old_bridge_local),
                    "A_C_new_l": finite_form_record(new_bridge_local),
                    "saturated_bridge_choice_count_per_orientation": 1,
                    "reason": "the raw bridges are saturated (m=1)",
                },
                "target_discriminant_forms": {
                    "q_frame_l": local_form_record(stored_old_gram, prime),
                    "q_NS_l": local_form_record(stored_old_gram, prime, sign=-1),
                },
                "marked_maximal_graph_count": len(graph_candidates),
                "q_NS_compatible_marked_graph_count": q_ns_compatible_count,
                "q_NS_selects_graph_label": q_ns_compatible_count < len(graph_candidates),
                "bridge_transport_coordinate_witness_count": len(transport_witnesses),
                "root_transition_histogram": transition_histogram_record(histogram),
                "distinct_root_transition_count": len(histogram),
                "actual_historical_transition": {
                    "old_root_system": actual_transition[0],
                    "new_root_system": actual_transition[1],
                    "multiplicity_among_marked_graphs": histogram[actual_transition],
                    "unique_after_declared_ADE_transition": histogram[actual_transition] == 1,
                },
                "prime_to_bad_glue_policy": "held at the stored value",
            }
        )

    assert len(answer) == 7
    assert all(edge["det_A"] == 0 for edge in answer)
    assert all(
        edge["actual_historical_transition"]["unique_after_declared_ADE_transition"]
        for edge in answer
    )
    return answer


def r17_noncyclic_control():
    pinned = load_json(R17_CERTIFICATE)
    assert pinned["status"] == "PASS_EXACT_R17_LOCAL_BRIDGE_MUTATION"
    short = load_matrix(R17_SHORT_GRAM)
    zero17 = vector(ZZ, [0] * 17)
    v = vector(ZZ, [1] + [0] * 16)
    w = vector(ZZ, [0] * 6 + [1] + [0] * 4 + [-1] + [0] * 5)
    old_bridge_basis = matrix(ZZ, [v, w])
    old_bridge = old_bridge_basis * short * old_bridge_basis.transpose()
    assert old_bridge == matrix(ZZ, [[4, 0], [0, 8]])
    r = -v + w
    r2 = -2 * v + w
    ambient = block_diagonal_matrix(HYPERBOLIC_PLANE, -short)
    new_u = matrix(ZZ, [[3, 2] + list(r), [4, 3] + list(r2)])
    assert new_u * ambient * new_u.transpose() == HYPERBOLIC_PLANE
    new_frame_basis = (new_u * ambient).right_kernel_matrix()
    new_frame = -(new_frame_basis * ambient * new_frame_basis.transpose())

    core_coordinates_old = matrix(
        ZZ, [r * short, r2 * short]
    ).right_kernel_matrix()
    core = core_coordinates_old * short * core_coordinates_old.transpose()
    core_basis = matrix(
        ZZ, [[0, 0] + list(item) for item in core_coordinates_old.rows()]
    )
    core_coordinates_new = new_frame_basis.solve_left(core_basis).change_ring(ZZ)
    new_bridge_coordinates = (
        core_coordinates_new * new_frame
    ).right_kernel_matrix()
    new_bridge = (
        new_bridge_coordinates
        * new_frame
        * new_bridge_coordinates.transpose()
    )
    assert new_bridge.det() == 32

    old_split_coordinates = core_coordinates_old.stack(old_bridge_basis)
    new_split_coordinates = core_coordinates_new.stack(new_bridge_coordinates)
    old_actual = graph_generators_from_split(old_split_coordinates)
    new_actual = graph_generators_from_split(new_split_coordinates)
    assert [order for order, _ in old_actual] == [4, 8]
    assert [order for order, _ in new_actual] == [4, 8]
    actual_old_generators = [item for _, item in old_actual]
    actual_new_generators = [item for _, item in new_actual]
    orders = [4, 8]
    actual_old_key = subgroup_key(actual_old_generators, orders)
    actual_new_key = subgroup_key(actual_new_generators, orders)

    cross = matrix(ZZ, [[2, 3], [3, 4]])
    old_raw, new_raw, transport_witnesses = bridge_transport_witnesses(
        cross,
        old_bridge,
        new_bridge,
        actual_old_generators,
        actual_new_key,
        orders,
    )
    assert old_raw == matrix(ZZ, [[12, 16], [16, 24]])
    assert new_raw == old_raw

    core_local = discriminant_module(core).primary_part(2)
    bridge_local = discriminant_module(old_bridge).primary_part(2)
    bridge_generators = bridge_local.gens()
    assert tuple(map(int, bridge_local.invariants())) == (4, 8)
    assert [int(item.order()) for item in bridge_generators] == orders
    first_bridge, second_bridge = bridge_generators

    graph_candidates = {}
    for first_core in core_local:
        if first_core.order() != 4:
            continue
        if not (first_core.q() + first_bridge.q()).is_zero():
            continue
        for second_core in core_local:
            if second_core.order() != 8:
                continue
            if not (second_core.q() + second_bridge.q()).is_zero():
                continue
            if not (
                first_core.b(second_core)
                + first_bridge.b(second_bridge)
            ).is_zero():
                continue
            first_lift = vector(
                QQ,
                list(first_core.lift()) + list(first_bridge.lift()),
            )
            second_lift = vector(
                QQ,
                list(second_core.lift()) + list(second_bridge.lift()),
            )
            key = subgroup_key([first_lift, second_lift], orders)
            if len(key) == 32:
                graph_candidates.setdefault(key, (first_lift, second_lift))
    assert actual_old_key in graph_candidates
    assert len(graph_candidates) == 32

    old_split = block_diagonal_matrix(core, old_bridge)
    new_split = block_diagonal_matrix(core, new_bridge)
    stored_old_gram = even_overlattice_gram(old_split, actual_old_generators)
    stored_new_gram = even_overlattice_gram(new_split, actual_new_generators)
    assert abs(stored_old_gram.det()) == abs(stored_new_gram.det()) == 948
    target_form_key = finite_form_key(
        Genus(stored_old_gram).discriminant_form().primary_part(2)
    )
    assert target_form_key == finite_form_key(
        Genus(stored_new_gram).discriminant_form().primary_part(2)
    )

    histogram_results = []
    actual_results = []
    q_ns_compatible_counts = []
    for transport in transport_witnesses:
        histogram = Counter()
        actual_transition = None
        q_ns_compatible = 0
        for graph_key, old_lifts in graph_candidates.items():
            new_lifts = [transport(item) for item in old_lifts]
            old_candidate = even_overlattice_gram(old_split, list(old_lifts))
            new_candidate = even_overlattice_gram(new_split, new_lifts)
            assert abs(old_candidate.det()) == abs(new_candidate.det()) == 948
            if (
                finite_form_key(
                    Genus(old_candidate).discriminant_form().primary_part(2)
                )
                == finite_form_key(
                    Genus(new_candidate).discriminant_form().primary_part(2)
                )
                == target_form_key
            ):
                q_ns_compatible += 1
            else:
                continue
            transition = (
                root_system_type(old_candidate),
                root_system_type(new_candidate),
            )
            histogram[transition] += 1
            if graph_key == actual_old_key:
                actual_transition = transition
        assert actual_transition is not None
        histogram_results.append(tuple(sorted(histogram.items())))
        actual_results.append(actual_transition)
        q_ns_compatible_counts.append(q_ns_compatible)

    assert len(set(histogram_results)) == 1
    assert len(set(actual_results)) == 1
    assert len(set(q_ns_compatible_counts)) == 1
    histogram = Counter(dict(histogram_results[0]))
    actual_transition = actual_results[0]
    q_ns_compatible_count = q_ns_compatible_counts[0]
    assert actual_transition == ("0", "4A1")
    assert histogram[actual_transition] == 4
    assert q_ns_compatible_count == len(graph_candidates) == 32

    return {
        "name": "published R17 non-cyclic degree-two control",
        "cross_pairing_A": rows(cross),
        "det_A": int(cross.det()),
        "raw_bridge_grams": {"old": rows(old_raw), "new": rows(new_raw)},
        "saturation_index": 1,
        "ambient_determinant": 948,
        "core_determinant": int(core.det()),
        "common_bridge_determinant": 32,
        "shared_bad_prime": 2,
        "local_modules": {
            "A_K_2_invariants": list(map(int, core_local.invariants())),
            "A_C_old_2": finite_form_record(bridge_local),
            "A_C_new_2": finite_form_record(
                discriminant_module(new_bridge).primary_part(2)
            ),
            "noncyclic": True,
        },
        "target_discriminant_forms": {
            "q_frame_2": local_form_record(stored_old_gram, 2),
            "q_NS_2": local_form_record(stored_old_gram, 2, sign=-1),
        },
        "forced_local_glue_group_invariants": [4, 8],
        "forced_local_glue_order": 32,
        "marked_maximal_graph_count": len(graph_candidates),
        "q_NS_compatible_marked_graph_count": q_ns_compatible_count,
        "q_NS_selects_graph_label": False,
        "bridge_transport_coordinate_witness_count": len(transport_witnesses),
        "root_transition_histogram": transition_histogram_record(histogram),
        "distinct_root_transition_count": len(histogram),
        "actual_transition": {
            "old_root_system": actual_transition[0],
            "new_root_system": actual_transition[1],
            "multiplicity_among_marked_graphs": histogram[actual_transition],
            "unique_after_declared_ADE_transition": False,
        },
        "conclusion": (
            "The maximal graph is genuinely non-cyclic. The finite quadratic "
            "form and the rootless-to-4A1 transition leave four marked graph "
            "choices, so the cyclic one-generator model and ADE uniqueness both fail."
        ),
    }


def build_result():
    historical = historical_edge_census()
    r17 = r17_noncyclic_control()
    return {
        "schema": "elkies-k3.prime-local-bridge-mutation.v1",
        "status": "PASS_EXACT_PRIME_LOCAL_BRIDGE_MUTATION_CLASSIFICATION",
        "theorem": {
            "name": "Prime-local bridge-mutation normal form H-1d",
            "raw_bridges": "G_0=A^t J A-J and G_1=A J A^t-J",
            "saturation": (
                "At each prime l, order-l^v_l(m) isotropic subgroups of "
                "A_(B_i,l) classify the local saturated bridges C_i."
            ),
            "graph_glue": (
                "For fixed K and C_i, local graph choices are isotropic "
                "H_(i,l) in A_(K,l)+A_(C_i,l), with injective projections "
                "and H_(i,l)^perp/H_(i,l) isomorphic to q_(frame,l)=-q_(NS,l)."
            ),
            "coupling": (
                "The rank-four bridge isometry tau_l couples the two sides: "
                "H_(1,l)=(id on A_K + tau_l)(H_(0,l))."
            ),
            "forced_order": (
                "2*v_l(|H_l|)=v_l(det K)+v_l(det C)-v_l(|det NS|)."
            ),
            "good_prime_rigidity": (
                "If l does not divide |det NS|, the local quotient is zero; "
                "both projections are onto the full discriminant modules and "
                "the graph is a single unmarked local orbit."
            ),
            "local_global_normal_form": (
                "Even overlattices are classified by the product of their "
                "primary saturation and graph data. A safe finite support is "
                "l dividing 2*|det NS|*det(G_0); after the saturated bridges "
                "are fixed, nontrivial graph defects occur only at l dividing "
                "gcd(det C, |det NS|)."
            ),
            "root_gate": (
                "Finite quadratic forms do not determine root births. For fixed "
                "metric K and C, each selected graph must be evaluated by the "
                "theta/coset minimum test after the primary choices are assembled."
            ),
            "literal_det_A_cutoff_is_invalid_here": (
                "All seven historical shared-prime cross matrices have det(A)=0. "
                "Thus l | 2*D*det(A) is not a finite prime set; det(G_A) is the "
                "required nonzero replacement in the rank-four case."
            ),
        },
        "historical_shared_prime_edges": historical,
        "historical_aggregate": {
            "edge_count": len(historical),
            "bad_primes": sorted(
                set(edge["shared_bad_prime"] for edge in historical)
            ),
            "all_raw_bridges_saturated": all(
                edge["saturation_index"] == 1 for edge in historical
            ),
            "all_local_bridge_modules_cyclic": all(
                len(edge["local_modules"]["A_C_old_l"]["invariants"]) == 1
                for edge in historical
            ),
            "all_orders_forced_maximal_by_core_determinant": all(
                edge["core_and_q_NS_forced_local_glue_order"]
                == edge["shared_bad_prime"]
                ** edge["valuations"]["v_l_det_C"]
                for edge in historical
            ),
            "q_NS_selects_no_marked_graph_label": all(
                not edge["q_NS_selects_graph_label"] for edge in historical
            ),
            "actual_ADE_transition_unique_on_every_edge": all(
                edge["actual_historical_transition"][
                    "unique_after_declared_ADE_transition"
                ]
                for edge in historical
            ),
        },
        "r17_noncyclic_negative_control": r17,
        "inputs": {
            "paths": [
                str(BRIDGE_CORPUS),
                str(RELATIVE_CORPUS),
                str(R17_CERTIFICATE),
                str(R17_SHORT_GRAM),
            ],
            "sha256": {
                str(path): sha256(path)
                for path in (
                    BRIDGE_CORPUS,
                    RELATIVE_CORPUS,
                    R17_CERTIFICATE,
                    R17_SHORT_GRAM,
                )
            },
        },
        "software": {
            "sage": SAGE_VERSION,
            "arithmetic": (
                "exact finite quadratic modules, Smith forms, rational "
                "overlattice bases, PARI qfminim/qfauto, and Sage integral isometry"
            ),
        },
        "proof_boundary": {
            "proved": (
                "The local normal form is finite discriminant-form/overlattice "
                "algebra. The census exhausts every marked bad-prime graph while "
                "holding the stored good-prime glue fixed, transports it to the "
                "opposite bridge, and computes every norm-two root system exactly."
            ),
            "not_proved": (
                "The census does not vary the common core, the prime-to-bad glue, "
                "or equation/Galois markings. Good-prime rigidity is an unmarked "
                "local statement; global root equivalence still depends on whether "
                "discriminant automorphisms lift to the metric lattices."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_prime_local_bridge_mutation.sage --check"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    result = build_result()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        assert output.read_text() == encoded, f"stale artifact: {output}"
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "PRIMELOCALBRIDGE|historical=7|bad_primes=3,5|"
        "historical_ADE_unique=7|r17_graphs=32|r17_ADE_multiplicity=4|"
        f"status={result['status']}|output={output}"
    )


if __name__ == "__main__":
    main()
