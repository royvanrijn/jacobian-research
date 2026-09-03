#!/usr/bin/env sage-python
"""Certify the complete dual-layer defect transition on the Q80 path.

For a good-prime neighbour N=M+Z*y/p, the checker predicts every forbidden
dual vector from the parent, the line, and affine CVP queries in M.  Only
after the prediction is complete does it construct N and compare the full
physical witness set with an independent child-dual enumeration.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, inverse_mod, lcm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
COMPLETION_SCRIPT = (
    ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_q80_defect_completion.sage"
)
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-q80-defect-birth-death-v1.json"


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def norm(gram, value):
    return value.dot_product(gram * value)


def adjusted_isotropic_lift(gram, prime, line):
    """Lift ``line`` modulo p so its norm is zero modulo 2*p^2."""

    quadratic_value = ZZ(norm(gram, line) / 2)
    assert quadratic_value % prime == 0
    covector = gram * line
    pivot = next(index for index, value in enumerate(covector) if value % prime)
    correction = (
        -ZZ(quadratic_value // prime)
        * inverse_mod(ZZ(covector[pivot]), prime)
    ) % prime
    lift = vector(ZZ, line)
    lift[pivot] += prime * correction
    assert norm(gram, lift) % (2 * prime**2) == 0
    return lift


def congruence_kernel_basis(gram, prime, lift):
    """Basis of M={z in K:<z,lift>=0 mod p} in parent coordinates."""

    covector = gram * lift
    pivot = next(index for index, value in enumerate(covector) if value % prime)
    pivot_inverse = inverse_mod(ZZ(covector[pivot]), prime)
    rows = []
    for index in range(gram.nrows()):
        if index == pivot:
            continue
        row = vector(ZZ, [0] * gram.nrows())
        row[index] = 1
        row[pivot] = (-ZZ(covector[index]) * pivot_inverse) % prime
        rows.append(row)
    row = vector(ZZ, [0] * gram.nrows())
    row[pivot] = prime
    rows.append(row)
    basis = matrix(ZZ, rows)
    assert abs(basis.det()) == prime
    assert all(
        norm(gram, row + lift) - norm(gram, row) - norm(gram, lift)
        == 2 * row.dot_product(gram * lift)
        for row in basis.rows()
    )
    assert all(row.dot_product(gram * lift) % prime == 0 for row in basis.rows())
    return basis, pivot, pivot_inverse


def affine_vectors_of_norm(gram, lattice_basis, shift, target_norm):
    """Enumerate exactly the vectors in ``lattice_basis*Z + shift`` of a norm."""

    reduced_gram = lattice_basis * gram * lattice_basis.transpose()
    target = -shift * lattice_basis.inverse()
    # PARI uses floating point internally for its bound.  The rational guard
    # prevents a boundary vector from being dropped; exact QQ filtering below
    # remains the authority.
    bounded = pari(reduced_gram).qfcvp(
        pari(target), B=QQ(target_norm) + QQ(1) / 10**8, flag=0
    )
    candidates = matrix(ZZ, bounded[2].sage()).columns()
    answer = []
    for candidate in candidates:
        value = (vector(QQ, candidate) - target) * lattice_basis
        if norm(gram, value) == target_norm:
            answer.append(value)
    return answer


def pairing_kernel_basis(prime, lift):
    """Basis of integer dual-pairing rows a satisfying a.y=0 mod p."""

    pivot = next(index for index, value in enumerate(lift) if value % prime)
    pivot_inverse = inverse_mod(ZZ(lift[pivot]), prime)
    rows = []
    for index in range(len(lift)):
        if index == pivot:
            continue
        row = vector(ZZ, [0] * len(lift))
        row[index] = 1
        row[pivot] = (-ZZ(lift[index]) * pivot_inverse) % prime
        rows.append(row)
    row = vector(ZZ, [0] * len(lift))
    row[pivot] = prime
    rows.append(row)
    basis = matrix(ZZ, rows)
    assert abs(basis.det()) == prime
    assert all(row.dot_product(lift) % prime == 0 for row in basis.rows())
    return basis


def full_dual_layer_profile(gram, prime, lift, reverse):
    """Compute the complete norm-at-most-two child theta profile by layers."""

    pairing_kernel = pairing_kernel_basis(prime, lift)
    dual_basis = pairing_kernel * gram.inverse()
    dual_gram = dual_basis * gram * dual_basis.transpose()
    denominator = lcm(value.denominator() for value in dual_gram.list())
    scaled_gram = (denominator * dual_gram).change_ring(ZZ)
    inverse_basis = dual_basis.inverse()
    profile = Counter()
    physical = set()
    layer_counts = []
    for layer in range(prime):
        shift = QQ(layer) * lift / prime
        target = -shift * inverse_basis
        bounded = pari(scaled_gram).qfcvp(
            pari(target), B=2 * denominator + QQ(1) / 10**8, flag=0
        )
        candidates = matrix(ZZ, bounded[2].sage()).columns()
        layer_count = 0
        for candidate in candidates:
            value = (vector(QQ, candidate) - target) * dual_basis
            value_norm = norm(gram, value)
            if not 0 <= value_norm <= 2:
                continue
            base_value = value - shift
            residue = reverse["discriminant_class"](base_value)
            profile[(residue, value_norm)] += 1
            physical.add(tuple(value))
            layer_count += 1
        layer_counts.append(layer_count)
    assert sum(profile.values()) == len(physical)
    return profile, physical, layer_counts


def materialized_dual_vectors_through_two(gram):
    """Independent enumeration of every child-dual vector through norm two."""

    inverse = gram.inverse()
    denominator = lcm(value.denominator() for value in inverse.list())
    scaled_inverse = (denominator * inverse).change_ring(ZZ)
    enumeration = pari(scaled_inverse).qfminim(2 * denominator)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    answer = {tuple(vector(QQ, [0] * gram.nrows()))}
    for representative in representatives:
        for pairing in (
            vector(ZZ, representative),
            -vector(ZZ, representative),
        ):
            dual_vector = pairing * inverse
            if norm(gram, dual_vector) <= 2:
                answer.add(tuple(dual_vector))
    assert len(answer) == int(enumeration[0]) + 1
    return answer


def profile_digest(profile):
    rows = [
        {
            "class": [str(value) for value in residue],
            "norm": str(value_norm),
            "multiplicity": multiplicity,
        }
        for (residue, value_norm), multiplicity in sorted(profile.items())
    ]
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode()).hexdigest(), len(rows)


def detailed_mask(gram, bridge, order, base, reverse):
    masks, _, generator = base["mask_profile"](gram, [bridge], reverse)
    assert len(masks) == 1
    multiplier = masks[0]["isotropic_multipliers"][0]
    return reverse["reverse_mask"](
        reverse["CoreCellOracle"](gram),
        bridge["theta_profile"],
        generator,
        bridge["generator"],
        multiplier,
        order,
    )


def predict_mask_layers(gram, prime, line, requirements):
    """Predict the requested child cells without constructing the child."""

    lift = adjusted_isotropic_lift(gram, prime, line)
    kernel, pivot, pivot_inverse = congruence_kernel_basis(gram, prime, lift)
    rows = []
    physical = set()
    for requirement_index, requirement in enumerate(requirements):
        residue = vector(
            QQ, [QQ(value) for value in requirement["core_discriminant_class"]]
        )
        target_norm = QQ(requirement["required_core_norm"])
        dual_pairing = residue * gram
        assert dual_pairing in ZZ**gram.nrows()

        # Choose k0 in K with <residue+k0,lift>=0 mod p.  Then every vector
        # in this child discriminant cell and layer j is in the affine coset
        # M + residue + k0 + j*lift/p.
        correction = vector(ZZ, [0] * gram.nrows())
        correction[pivot] = (
            -ZZ(dual_pairing.dot_product(lift)) * pivot_inverse
        ) % prime
        assert (residue + correction).dot_product(gram * lift) % prime == 0

        layers = []
        for layer in range(prime):
            shift = residue + correction + QQ(layer) * lift / prime
            values = affine_vectors_of_norm(gram, kernel, shift, target_norm)
            for value in values:
                physical.add(tuple(value))
            if values:
                layers.append(
                    {
                        "layer": layer,
                        "physical_witness_count": len(values),
                    }
                )
        if layers:
            rows.append(
                {
                    "requirement_index": requirement_index,
                    "graph_label": requirement["graph_label"],
                    "core_discriminant_class": requirement[
                        "core_discriminant_class"
                    ],
                    "required_core_norm": requirement["required_core_norm"],
                    "occupied_layers": layers,
                }
            )
    return {
        "adjusted_lift": lift,
        "kernel_basis": kernel,
        "occupied_requirements": rows,
        "physical_witnesses": physical,
    }


def model_neighbor_basis(kernel, lift, prime):
    """Basis of M+Z*lift/p, formed only for the post-prediction comparison."""

    generators = matrix(
        ZZ,
        [list(prime * row) for row in kernel.rows()] + [list(lift)],
    )
    numerator_basis = generators.row_module(ZZ).basis_matrix()
    assert numerator_basis.nrows() == kernel.nrows()
    return numerator_basis.change_ring(QQ) / prime


def rational_vector_record(value):
    return [str(entry) for entry in value]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    completion = runpy.run_path(str(COMPLETION_SCRIPT))
    directed = runpy.run_path(str(completion["DIRECTED_SCRIPT"]))
    base = runpy.run_path(str(directed["BASE_SCRIPT"]))
    search = runpy.run_path(str(directed["SEARCH_SCRIPT"]))
    control = runpy.run_path(str(directed["CONTROL_SCRIPT"]))
    core = runpy.run_path(str(directed["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(directed["REVERSE_SCRIPT"]))
    prepared, bridge, initial_gram, _, _, _ = directed["initial_q80"](
        base, search, control, core, reverse
    )
    current = base["quadratic_form"](initial_gram)
    transitions = []
    witness_counts = []

    initial_mask, _, initial_witnesses = directed["masked_witness_data"](
        initial_gram, bridge, prepared["order"], base, reverse
    )
    assert initial_mask["occupied_forbidden_cells"] == 2
    witness_counts.append(len(initial_witnesses))

    for step, (prime, raw_line) in enumerate(completion["DIRECTED_PATH"], start=1):
        parent_gram = current.Hessian_matrix()
        parent_mask = detailed_mask(
            parent_gram, bridge, prepared["order"], base, reverse
        )
        _, parent_cells, parent_witnesses = directed["masked_witness_data"](
            parent_gram, bridge, prepared["order"], base, reverse
        )
        line = vector(ZZ, raw_line)
        pairings = [
            int(witness.dot_product(line)) % prime
            for witness in parent_witnesses
        ]
        assert all(pairings)

        # This is the prediction boundary: no child lattice, child Gram matrix,
        # or child dual enumeration exists before this call returns.
        prediction = predict_mask_layers(
            parent_gram, prime, line, parent_mask["requirements"]
        )
        predicted_vectors = prediction["physical_witnesses"]
        predicted_cells = len(prediction["occupied_requirements"])
        full_profile, full_predicted_vectors, full_layer_counts = (
            full_dual_layer_profile(
                parent_gram,
                prime,
                prediction["adjusted_lift"],
                reverse,
            )
        )
        full_profile_hash, full_profile_cells = profile_digest(full_profile)
        assert all(
            layer["layer"] != 0
            for row in prediction["occupied_requirements"]
            for layer in row["occupied_layers"]
        )

        # Materialize the neighbour only after the affine prediction, and use
        # the existing independent child-dual enumerator as the truth set.
        transform = current.find_p_neighbor_from_vec(
            prime, line, return_matrix=True
        )
        neighbor_basis = transform.transpose()
        child = current.find_p_neighbor_from_vec(prime, line)
        child_gram = child.Hessian_matrix()
        model_basis = model_neighbor_basis(
            prediction["kernel_basis"], prediction["adjusted_lift"], prime
        )
        change = neighbor_basis * model_basis.inverse()
        assert all(value in ZZ for value in change.list())
        assert abs(change.det()) == 1

        child_mask, child_cells, child_witnesses = directed["masked_witness_data"](
            child_gram, bridge, prepared["order"], base, reverse
        )
        actual_vectors = set()
        for witness in child_witnesses:
            child_dual_vector = witness * child_gram.inverse()
            actual_vectors.add(tuple(child_dual_vector * neighbor_basis))
        assert predicted_vectors == actual_vectors
        assert predicted_cells == len(child_cells)
        assert predicted_cells == child_mask["occupied_forbidden_cells"]

        actual_child_coordinates = materialized_dual_vectors_through_two(child_gram)
        full_actual_vectors = {
            tuple(vector(QQ, value) * neighbor_basis)
            for value in actual_child_coordinates
        }
        assert full_predicted_vectors == full_actual_vectors
        discriminant_exponent = lcm(
            value.denominator() for value in parent_gram.inverse().list()
        )
        inverse_prime = inverse_mod(prime, discriminant_exponent)
        actual_profile = Counter()
        for value_tuple in full_actual_vectors:
            value = vector(QQ, value_tuple)
            residue = reverse["discriminant_class"](
                inverse_prime * prime * value
            )
            actual_profile[(residue, norm(parent_gram, value))] += 1
        assert full_profile == actual_profile

        witness_counts.append(len(actual_vectors))
        transitions.append(
            {
                "step_after_near_miss": step,
                "prime": prime,
                "line": list(map(int, raw_line)),
                "adjusted_isotropic_lift": list(
                    map(int, prediction["adjusted_lift"])
                ),
                "lift_norm_divided_by_2p2": int(
                    norm(parent_gram, prediction["adjusted_lift"])
                    / (2 * prime**2)
                ),
                "parent_occupied_cells": len(parent_cells),
                "parent_physical_witnesses": len(parent_witnesses),
                "nonzero_parent_pairings": pairings,
                "predicted_before_child_construction": True,
                "predicted_child_occupied_cells": predicted_cells,
                "predicted_child_physical_witnesses": len(predicted_vectors),
                "occupied_affine_layers": prediction["occupied_requirements"],
                "j_zero_layer_empty_on_reverse_mask": True,
                "prediction_equals_materialized_child": True,
                "full_sigma2_transition": {
                    "theta_cells": full_profile_cells,
                    "dual_vectors_through_norm_two": len(full_predicted_vectors),
                    "vectors_by_affine_layer": full_layer_counts,
                    "canonical_profile_sha256": full_profile_hash,
                    "equals_materialized_child": True,
                },
                "predicted_witnesses_in_parent_ambient_coordinates": [
                    rational_vector_record(vector(QQ, value))
                    for value in sorted(predicted_vectors)
                ],
            }
        )
        current = child

    assert witness_counts == [4, 6, 4, 4, 0]
    assert [row["predicted_child_occupied_cells"] for row in transitions] == [
        2,
        2,
        2,
        0,
    ]

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-q80-defect-birth-death.v1",
        "status": "PASS_EXACT_Q80_DEFECT_BIRTH_DEATH_LAYERS",
        "inputs": {
            relative(COMPLETION_SCRIPT): digest(COMPLETION_SCRIPT),
            relative(completion["DIRECTED_SCRIPT"]): digest(
                completion["DIRECTED_SCRIPT"]
            ),
            relative(directed["BRIDGES"]): digest(directed["BRIDGES"]),
            relative(directed["THETA"]): digest(directed["THETA"]),
            relative(directed["BASE_SCRIPT"]): digest(directed["BASE_SCRIPT"]),
            relative(directed["CONTROL_SCRIPT"]): digest(
                directed["CONTROL_SCRIPT"]
            ),
            relative(directed["CORE_SCRIPT"]): digest(directed["CORE_SCRIPT"]),
            relative(directed["REVERSE_SCRIPT"]): digest(
                directed["REVERSE_SCRIPT"]
            ),
        },
        "dual_layer_theorem": {
            "statement": (
                "For K_y^dual={x in K^dual:<x,y>=0 mod p}, "
                "K_l^dual is the disjoint union of K_y^dual+j*y/p for "
                "0<=j<p."
            ),
            "canonical_discriminant_identification": (
                "The child class of x+j*y/p maps to the parent class of x; "
                "equivalently [v] maps to p^(-1)[p*v] because p is prime "
                "to det(K)."
            ),
            "affine_cell_query": (
                "For a parent class representative r and k0 solving "
                "<r+k0,y>=0 mod p, query M+r+k0+j*y/p at the requested "
                "norm for j=0,...,p-1."
            ),
        },
        "strong_acceptance_criterion": {
            "old_death": (
                "The line is nonorthogonal to every current physical "
                "forbidden witness; equivalently every requested j=0 layer "
                "is empty."
            ),
            "no_birth": (
                "Every requested j=1,...,p-1 affine layer is empty."
            ),
            "equivalence_checked_on_all_transitions": True,
        },
        "physical_witness_count_regression": witness_counts,
        "transitions": transitions,
        "proof_boundary": (
            "The layer decomposition and finite affine-CVP transition are "
            "general for positive even good-prime neighbours. The abstract "
            "counted signature Sigma_2(K) alone is not sufficient: the line "
            "pairings and physical representatives, or an equivalent affine "
            "CVP oracle for K, are required. This checker proves exact outcome "
            "prediction on the four stored Q80 transitions, not a uniform "
            "running-time speedup."
        ),
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_q80_defect_birth_death.sage --check"
        ),
    }

    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS exact Q80 defect birth-death layers")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
