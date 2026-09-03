#!/usr/bin/env sage-python
"""Certify inverse ADE constraints on the first NS0024 core neighbour.

The checker compiles the desired child root system from the parent core K,
the fixed binary bridge C, its graph glue H, and the isotropic line y.  It
enumerates the affine dual layers before constructing the child core.  Only
after the predicted physical root set and its metric have been fixed does it
materialize the neighbour and compare an independent completed-frame root
enumeration.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import GF, QQ, ZZ, block_diagonal_matrix, lcm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
ROUTE_SCRIPT = ROOT / "elkies-k3/scripts/certify_ns0024_new_rootless_source_route.sage"
SIGNATURE_SCRIPT = (
    ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_root_system_signature.sage"
)
BIRTH_DEATH_SCRIPT = (
    ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_q80_defect_birth_death.sage"
)
OUTPUT = GENERATED / "elkies-k3-ns0024-inverse-ade-mutation-v1.json"


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational_vector(value):
    return [str(entry) for entry in value]


def norm(gram, value):
    return value * gram * value


def dual_vectors_through_two(gram, discriminant_class):
    """Enumerate every dual vector of norm at most two, with its class."""

    inverse = gram.inverse()
    denominator = lcm(entry.denominator() for entry in inverse.list())
    scaled_inverse = (denominator * inverse).change_ring(ZZ)
    enumeration = pari(scaled_inverse).qfminim(2 * denominator)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    zero = vector(QQ, [0] * gram.nrows())
    answer = [(zero, QQ(0), discriminant_class(zero))]
    for representative in representatives:
        for pairing in (vector(ZZ, representative), -vector(ZZ, representative)):
            value = pairing * inverse
            value_norm = norm(gram, value)
            if value_norm <= 2:
                answer.append((value, value_norm, discriminant_class(value)))
    assert len(answer) == int(enumeration[0]) + 1
    assert len({tuple(row[0]) for row in answer}) == len(answer)
    return answer


def child_dual_affine_layers(gram, prime, lift, birth_death, discriminant_class):
    """Enumerate K_ell^dual through norm two in its p affine layers."""

    pairing_kernel = birth_death["pairing_kernel_basis"](prime, lift)
    dual_basis = pairing_kernel * gram.inverse()
    dual_gram = dual_basis * gram * dual_basis.transpose()
    denominator = lcm(entry.denominator() for entry in dual_gram.list())
    scaled_gram = (denominator * dual_gram).change_ring(ZZ)
    inverse_basis = dual_basis.inverse()
    answer = []
    for layer in range(prime):
        shift = QQ(layer) * lift / prime
        target = -shift * inverse_basis
        bounded = pari(scaled_gram).qfcvp(
            pari(target), B=2 * denominator + QQ(1) / 10**8, flag=0
        )
        for candidate in matrix(ZZ, bounded[2].sage()).columns():
            value = (vector(QQ, candidate) - target) * dual_basis
            value_norm = norm(gram, value)
            if not 0 <= value_norm <= 2:
                continue
            base = value - shift
            answer.append(
                (value, value_norm, discriminant_class(base), layer)
            )
    assert len({tuple(row[0]) for row in answer}) == len(answer)
    return answer


def canonical_line(core_value, bridge_value):
    positive = tuple(core_value) + tuple(bridge_value)
    negative = tuple(-entry for entry in positive)
    return min(positive, negative)


def graph_root_witnesses(
    core_vectors,
    bridge_vectors,
    core_generator,
    bridge_generator,
    multiplier,
    order,
    core_rank,
    layer_modulus=None,
):
    """Join core and bridge dual shells along H and total norm two."""

    graph_by_bridge_class = {}
    graph_labels = {}
    for label in range(order):
        core_class = tuple(
            entry for entry in discriminant_class(label * multiplier * core_generator)
        )
        bridge_class = tuple(
            entry for entry in discriminant_class(label * bridge_generator)
        )
        assert bridge_class not in graph_by_bridge_class
        graph_by_bridge_class[bridge_class] = core_class
        graph_labels[bridge_class] = label

    core_by_cell = defaultdict(list)
    for row in core_vectors:
        value, value_norm, residue = row[:3]
        layer = row[3] if len(row) == 4 else None
        core_by_cell[(tuple(residue), value_norm)].append((value, layer))

    signed = {}
    for bridge_value, bridge_norm, bridge_class in bridge_vectors:
        bridge_class = tuple(bridge_class)
        if bridge_class not in graph_by_bridge_class:
            continue
        core_class = graph_by_bridge_class[bridge_class]
        target_norm = QQ(2) - bridge_norm
        if not 0 <= target_norm <= 2:
            continue
        for core_value, layer in core_by_cell.get((core_class, target_norm), []):
            line = canonical_line(core_value, bridge_value)
            positive = tuple(core_value) + tuple(bridge_value)
            sign = 1 if positive == line else -1
            canonical_bridge_class = discriminant_class(line[core_rank:])
            canonical_layer = layer
            if layer is not None and sign == -1:
                assert layer_modulus is not None
                canonical_layer = (-layer) % layer_modulus
            record = {
                "core": vector(QQ, line[:core_rank]),
                "bridge": vector(QQ, line[core_rank:]),
                "graph_label": graph_labels[canonical_bridge_class],
                "affine_layer": canonical_layer,
            }
            if line in signed:
                assert signed[line] == record
            signed[line] = record
    return [signed[key] for key in sorted(signed)]


def discriminant_class(value):
    return tuple(entry - entry.floor() for entry in map(QQ, value))


def root_metric_signature(witnesses, core_gram, bridge_gram):
    split = block_diagonal_matrix(core_gram, bridge_gram)
    roots = [vector(QQ, list(row["core"]) + list(row["bridge"])) for row in witnesses]
    pairings = matrix(QQ, roots) * split * matrix(QQ, roots).transpose()
    adjacency = [set() for _ in roots]
    edges = []
    for left in range(len(roots)):
        assert pairings[left, left] == 2
        for right in range(left + 1, len(roots)):
            value = ZZ(pairings[left, right])
            assert value in (-1, 0, 1)
            if value:
                adjacency[left].add(right)
                adjacency[right].add(left)
                edges.append([left, right, int(value)])

    unseen = set(range(len(roots)))
    components = []
    names = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        unseen.remove(seed)
        indices = []
        while stack:
            item = stack.pop()
            indices.append(item)
            for neighbor in adjacency[item]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        indices.sort()
        rank = int(matrix(QQ, [roots[index] for index in indices]).rank())
        signed_count = 2 * len(indices)
        if signed_count == rank * (rank + 1):
            name = f"A{rank}"
        elif rank >= 4 and signed_count == 2 * rank * (rank - 1):
            name = f"D{rank}"
        else:
            name = {(6, 72): "E6", (7, 126): "E7", (8, 240): "E8"}[
                (rank, signed_count)
            ]
        names.append(name)
        components.append(
            {
                "type": name,
                "root_line_indices": indices,
                "rank": rank,
                "signed_root_count": signed_count,
            }
        )
    counts = Counter(names)
    family_order = {"A": 0, "D": 1, "E": 2}
    labels = sorted(counts, key=lambda item: (family_order[item[0]], int(item[1:])))
    ade = "+".join(
        (str(counts[label]) if counts[label] > 1 else "") + label
        for label in labels
    ) or "rootless"
    triangular = [
        [int(pairings[left, right]) for right in range(left, len(roots))]
        for left in range(len(roots))
    ]
    return {
        "ade_type": ade,
        "root_rank": int(matrix(QQ, roots).rank()) if roots else 0,
        "signed_root_count": 2 * len(roots),
        "root_line_count": len(roots),
        "components": components,
        "nonorthogonal_root_line_edges": edges,
        "pairwise_inner_products_upper_triangular": triangular,
        "pairwise_inner_products_sha256": hashlib.sha256(
            json.dumps(triangular, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def physical_key(core_value, bridge_value):
    return canonical_line(core_value, bridge_value)


def actual_child_witnesses(signature, child_to_parent):
    answer = []
    for row in signature["physical_bridge_witnesses"]["root_lines"]:
        core_child = vector(QQ, map(QQ, row["core_dual_coordinates"]))
        bridge = vector(QQ, map(QQ, row["bridge_dual_coordinates"]))
        core_parent = core_child * child_to_parent
        line = canonical_line(core_parent, bridge)
        answer.append(line)
    assert len(set(answer)) == len(answer)
    return set(answer)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    route = runpy.run_path(str(ROUTE_SCRIPT))
    signature_tools = runpy.run_path(str(SIGNATURE_SCRIPT))
    birth_death = runpy.run_path(str(BIRTH_DEATH_SCRIPT))
    base = runpy.run_path(str(route["BASE_SCRIPT"]))
    search = runpy.run_path(str(route["SEARCH_SCRIPT"]))
    core_tools = runpy.run_path(str(route["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(route["REVERSE_SCRIPT"]))
    bridges = json.loads(route["BRIDGES"].read_text())
    theta = json.loads(route["THETA"].read_text())
    masked = json.loads(route["MASKED"].read_text())

    prepared = search["prepare_corridor"](
        "NS0024", bridges, theta, base, core_tools, reverse
    )
    base["BRIDGE_ORDER"] = prepared["order"]
    search["configure_order"](base, prepared["order"])
    masked_row = next(row for row in masked["corridors"] if row["corridor"] == "NS0024")
    bridge_index = masked_row["completion"]["bridge_class_index"]
    bridge = next(
        row for row in prepared["viable_bridges"]
        if row["bridge_class_index"] == bridge_index
    )

    order = int(prepared["order"])
    parent_core = base["lll_reduce"](prepared["seed"])
    prime, raw_line = route["PATH"][0]
    line = vector(ZZ, raw_line)
    parent_quadratic = base["quadratic_form"](parent_core)
    assert parent_quadratic(line) % prime == 0
    adjusted_lift = birth_death["adjusted_isotropic_lift"](
        parent_core, prime, vector(ZZ, list(line))
    )

    parent_choices = route["completed_frames"](
        parent_core, bridge, order, base, core_tools
    )
    parent_multiplier, parent_frame = min(
        (multiplier, frame)
        for multiplier, frame in parent_choices
        if route["root_data"](frame) == (13, 280, 4)
    )
    assert parent_multiplier == 59
    core_generator = base["primary_generator"](parent_core, order)
    bridge_generator = bridge["generator"]
    bridge_vectors = dual_vectors_through_two(
        bridge["gram"], reverse["discriminant_class"]
    )
    parent_core_vectors = dual_vectors_through_two(
        parent_core, reverse["discriminant_class"]
    )
    parent_witnesses = graph_root_witnesses(
        parent_core_vectors,
        bridge_vectors,
        core_generator,
        bridge_generator,
        parent_multiplier,
        order,
        parent_core.nrows(),
    )
    parent_signature = root_metric_signature(
        parent_witnesses, parent_core, bridge["gram"]
    )
    assert parent_signature["ade_type"] == "D5+E8"

    incidence_rows = []
    survivor_keys = set()
    for index, witness in enumerate(parent_witnesses):
        core_value = witness["core"]
        dual_pairing = core_value * parent_core
        assert dual_pairing in ZZ ** parent_core.nrows()
        residue = int(dual_pairing.dot_product(line) % prime)
        survives = residue == 0
        if survives:
            survivor_keys.add(physical_key(core_value, witness["bridge"]))
        incidence_rows.append(
            {
                "parent_root_line_index": index,
                "core_dual_coordinates": rational_vector(core_value),
                "bridge_dual_coordinates": rational_vector(witness["bridge"]),
                "graph_glue_label": witness["graph_label"],
                "modular_linear_form": [int(value % prime) for value in dual_pairing],
                "value_on_selected_line": residue,
                "constraint": "equals_zero" if survives else "nonzero",
                "fate": "survives" if survives else "dies",
            }
        )

    # Prediction boundary: enumerate all p affine layers and join them to C
    # through the transported graph H before constructing the child core.
    child_core_vectors = child_dual_affine_layers(
        parent_core,
        prime,
        adjusted_lift,
        birth_death,
        reverse["discriminant_class"],
    )
    predicted_witnesses = graph_root_witnesses(
        child_core_vectors,
        bridge_vectors,
        core_generator,
        bridge_generator,
        parent_multiplier,
        order,
        parent_core.nrows(),
        layer_modulus=prime,
    )
    predicted_signature = root_metric_signature(
        predicted_witnesses, parent_core, bridge["gram"]
    )
    assert predicted_signature["ade_type"] == "3A1+A2"
    predicted_keys = {
        physical_key(row["core"], row["bridge"]) for row in predicted_witnesses
    }
    predicted_layer_counts = Counter(row["affine_layer"] for row in predicted_witnesses)
    affine_vector_layer_counts = Counter(row[3] for row in child_core_vectors)
    affine_theta_cells = {
        (tuple(row[2]), row[1]) for row in child_core_vectors
    }
    born_keys = predicted_keys - survivor_keys
    assert survivor_keys == predicted_keys
    assert not born_keys
    assert predicted_layer_counts == {0: 6}

    prediction_record = {
        "computed_before_child_construction": True,
        "desired_ade_type": "3A1+A2",
        "metric_signature": predicted_signature,
        "affine_layer_root_line_counts": {
            str(layer): count for layer, count in sorted(predicted_layer_counts.items())
        },
        "affine_dual_enumeration": {
            "core_dual_vectors_through_norm_two": len(child_core_vectors),
            "occupied_core_theta_cells": len(affine_theta_cells),
            "vectors_by_layer": {
                str(layer): count
                for layer, count in sorted(affine_vector_layer_counts.items())
            },
            "bridge_dual_vectors_through_norm_two": len(bridge_vectors),
        },
        "surviving_parent_root_lines": len(survivor_keys),
        "born_root_lines": len(born_keys),
        "no_additional_norm_two_witnesses": True,
        "root_lines_in_parent_core_plus_bridge_coordinates": [
            {
                "core_dual_coordinates": rational_vector(row["core"]),
                "bridge_dual_coordinates": rational_vector(row["bridge"]),
                "graph_glue_label": row["graph_label"],
                "affine_layer": row["affine_layer"],
            }
            for row in predicted_witnesses
        ],
    }

    # Independent truth check after the inverse constraint record is fixed.
    transform = parent_quadratic.find_p_neighbor_from_vec(
        prime, line, return_matrix=True
    )
    child_to_parent = transform.transpose()
    child_quadratic = parent_quadratic.find_p_neighbor_from_vec(prime, line)
    child_core = child_quadratic.Hessian_matrix()
    child_choices = route["completed_frames"](
        base["lll_reduce"](child_core), bridge, order, base, core_tools
    )
    child_multiplier, expected_child_frame = min(
        (multiplier, frame)
        for multiplier, frame in child_choices
        if route["root_data"](frame) == (5, 12, 24)
    )
    assert child_multiplier == 50

    # completed_frames uses an LLL-reduced child basis.  Recover that basis
    # explicitly so its physical roots can be transported back to the parent.
    lll_transform = matrix(ZZ, pari(child_core).qflllgram()).transpose()
    reduced_child = lll_transform * child_core * lll_transform.transpose()
    assert reduced_child == base["lll_reduce"](child_core)
    reduced_child_to_parent = lll_transform * child_to_parent
    child_core_generator = base["primary_generator"](reduced_child, order)
    child_glue = vector(
        QQ,
        list(child_multiplier * child_core_generator) + list(bridge_generator),
    )
    materialized_frame, split, split_basis = signature_tools["overlattice_with_basis"](
        reduced_child, bridge["gram"], child_glue
    )
    assert materialized_frame == expected_child_frame
    materialized_signature = signature_tools["metric_root_signature"](
        materialized_frame, split, split_basis, child_glue
    )
    actual_keys = actual_child_witnesses(
        materialized_signature, reduced_child_to_parent
    )
    assert actual_keys == predicted_keys
    assert materialized_signature["ade_type"] == predicted_signature["ade_type"]

    input_paths = (
        ROUTE_SCRIPT,
        SIGNATURE_SCRIPT,
        BIRTH_DEATH_SCRIPT,
        route["MASKED"],
        route["BRIDGES"],
        route["THETA"],
        route["BASE_SCRIPT"],
        route["SEARCH_SCRIPT"],
        route["CORE_SCRIPT"],
        route["REVERSE_SCRIPT"],
    )
    payload = {
        "schema": "elkies-k3.ns0024-inverse-ade-mutation.v1",
        "status": "PASS_EXACT_INVERSE_ADE_MUTATION_CONSTRAINTS",
        "transition": {
            "parent_ade_type": parent_signature["ade_type"],
            "desired_child_ade_type": "3A1+A2",
            "good_neighbor_prime": prime,
            "selected_isotropic_line": list(map(int, line)),
            "adjusted_isotropic_lift": list(map(int, adjusted_lift)),
            "lift_norm_divided_by_2p2": int(
                norm(parent_core, adjusted_lift) / (2 * prime**2)
            ),
            "cyclic_bridge_order": order,
            "bridge_class_index": bridge_index,
            "parent_graph_multiplier": int(parent_multiplier),
            "child_graph_multiplier_in_reduced_child_basis": int(child_multiplier),
        },
        "compiled_line_constraints": {
            "isotropy": "y^2=0 mod 2*p^2 after lift adjustment",
            "old_witness_incidence": incidence_rows,
            "old_root_line_counts": {
                "total": len(parent_witnesses),
                "must_survive": len(survivor_keys),
                "must_die": len(parent_witnesses) - len(survivor_keys),
                "survival_equality_rank_over_Fp": int(
                    matrix(
                        GF(prime),
                        [
                            row["modular_linear_form"]
                            for row in incidence_rows
                            if row["fate"] == "survives"
                        ],
                    ).rank()
                ),
            },
            "birth_constraints": (
                "For every graph class (a,b), bridge vector c in C^dual+b, "
                "and 0<=j<p, enumerate M+r_a+k0+j*y/p at norm 2-c^2."
            ),
            "target_metric_constraint": (
                "The complete joined witness set must have the prescribed "
                "pairwise metric of 3A1+A2."
            ),
            "exclusion_constraint": (
                "Every other joined affine shell of total norm two is empty."
            ),
        },
        "predicted_child": prediction_record,
        "post_prediction_materialized_check": {
            "performed_after_constraint_evaluation": True,
            "physical_root_set_equals_prediction": True,
            "pairwise_metric_equals_prediction": True,
            "ade_type": materialized_signature["ade_type"],
            "root_rank": materialized_signature["root_rank"],
            "signed_root_count": materialized_signature["signed_root_count"],
        },
        "algorithm_boundary": {
            "proved": (
                "For fixed positive even K,C,H and a good odd-prime line, "
                "the displayed modular incidences plus finite affine-CVP "
                "shells determine the complete child root metric before "
                "constructing the child lattice."
            ),
            "control": (
                "On the first NS0024 edge, exactly six of 140 parent root "
                "lines survive, 134 die, no nonzero affine layer contributes "
                "a root, and the survivors have type 3A1+A2."
            ),
            "not_proved": (
                "No completeness or complexity bound for finding a line "
                "satisfying a prescribed ADE target, elliptic-neighbour "
                "interpretation, equation lift, or field of definition is claimed."
            ),
        },
        "literature_context": [
            {
                "reference": "https://arxiv.org/abs/2104.06846",
                "role": "isotropic-line construction of good-prime Kneser neighbours",
            },
            {
                "reference": "https://arxiv.org/abs/2410.18788",
                "role": (
                    "visible-root filtering by modular incidence; the affine "
                    "graph-glue layers here account for nonvisible roots"
                ),
            },
        ],
        "inputs": {relative(path): digest(path) for path in input_paths},
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_ns0024_inverse_ade_mutation.sage --check"
        ),
    }

    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"missing or stale artifact: {output}")
        print("PASS exact inverse ADE mutation constraints")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
