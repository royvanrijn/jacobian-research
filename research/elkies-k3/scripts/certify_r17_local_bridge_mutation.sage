#!/usr/bin/env sage-python
"""Certify the local bridge-mutation law and the R17 non-cyclic example.

The general identities are proved in Theorem H-1c.  This checker supplies
exact symbolic controls, audits the 42-edge corpus, and constructs the new
degree-two primitive U embedding in U + R17(-1).  It also performs the chamber
and physical-zero gates and compares the resulting 4A1 frame with both stored
H3 4A1 frames.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path

from sage.env import SAGE_VERSION
from sage.all import (
    GF,
    PolynomialRing,
    QuadraticForm,
    ZZ,
    block_diagonal_matrix,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-local-bridge-mutation-v1.json"
)
SHORT_GRAM = Path("elkies-k3/data/lattice/short_vector_basis_gram.txt")
SHORT_COORDS = Path("elkies-k3/data/lattice/short_vector_basis_coords.txt")
PINNED_GRAM = Path("elkies-k3/data/lattice/rank17_gram.txt")
RELATIVE_CORPUS = Path(
    "artifacts/generated-results/"
    "elkies-k3-relative-u-bridge-lifting-regression-v1.json"
)
GLUE_CORPUS = Path(
    "artifacts/generated-results/"
    "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
)
HISTORICAL_MARKING = Path(
    "artifacts/generated-results/"
    "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
)
PHYSICAL_Q8_FRAME = Path(
    "artifacts/generated-results/"
    "elkies-k3-h3-q4o164-c8-q8o376-4a1-old_zero-frame.txt"
)
HYPERBOLIC_PLANE = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(relative_path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in (ROOT / relative_path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def load_json(relative_path):
    return json.loads((ROOT / relative_path).read_text())


def sha256(relative_path):
    return hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def smith_invariants(gram):
    diagonal = gram.smith_form()[0].diagonal()
    return [abs(int(value)) for value in diagonal if abs(int(value)) > 1]


def pari_short_vectors(gram, bound):
    result = pari(gram).qfminim(bound)
    positive = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2].sage()).columns()
    ]
    signed = positive + [-item for item in positive]
    assert len(signed) == int(result[0])
    return signed


def row_isometry(source, target):
    witness = QuadraticForm(ZZ, source).is_globally_equivalent_to(
        QuadraticForm(ZZ, target), return_matrix=True
    )
    assert witness is not False
    witness = matrix(ZZ, witness)
    candidates = [witness, witness.transpose()]
    inverse = witness.inverse()
    if inverse.change_ring(ZZ) == inverse:
        inverse = inverse.change_ring(ZZ)
        candidates.extend([inverse, inverse.transpose()])
    for candidate in candidates:
        if candidate * source * candidate.transpose() == target:
            assert abs(candidate.det()) == 1
            return candidate
    raise ArithmeticError("no verified row-convention integral isometry")


def symbolic_controls():
    ring = PolynomialRing(ZZ, names=("a", "b", "c", "d"))
    a, b, c, d = ring.gens()
    relative = matrix(ring, [[a, b], [c, d]])
    hyperbolic = HYPERBOLIC_PLANE.change_ring(ring)
    old_raw = relative.transpose() * hyperbolic * relative - hyperbolic
    new_raw = relative * hyperbolic * relative.transpose() - hyperbolic
    expected_old = matrix(ring, [[2 * a * c, a * d + b * c - 1], [a * d + b * c - 1, 2 * b * d]])
    assert old_raw == expected_old
    expected_det = 2 * (a * d + b * c) - 1 - (a * d - b * c) ** 2
    assert old_raw.det() == new_raw.det() == expected_det

    parity_ring = PolynomialRing(GF(2), names=("a", "b", "c", "d"))
    aa, bb, cc, dd = parity_ring.gens()
    parity_relative = matrix(parity_ring, [[aa, bb], [cc, dd]])
    parity_hyperbolic = HYPERBOLIC_PLANE.change_ring(parity_ring)
    parity_raw = (
        parity_relative.transpose()
        * parity_hyperbolic
        * parity_relative
        - parity_hyperbolic
    )
    parity_det = aa * dd - bb * cc
    assert parity_raw == matrix(
        parity_ring, [[0, parity_det - 1], [parity_det - 1, 0]]
    )

    counterexample_A = matrix(ZZ, [[2, 3], [6, 8]])
    counterexample_gram = (
        counterexample_A.transpose()
        * HYPERBOLIC_PLANE
        * counterexample_A
        - HYPERBOLIC_PLANE
    )
    assert counterexample_gram == matrix(ZZ, [[24, 33], [33, 48]])
    assert counterexample_gram.is_positive_definite()
    assert counterexample_gram.det() == 63
    assert smith_invariants(counterexample_gram) == [3, 21]
    counterexample_ns = block_diagonal_matrix(HYPERBOLIC_PLANE, -counterexample_gram)
    counterexample_target = matrix(ZZ, [[6, 2, 1, 0], [8, 3, 0, 1]])
    assert counterexample_target * counterexample_ns * counterexample_target.transpose() == HYPERBOLIC_PLANE

    return {
        "raw_gram_old_to_new": "A^t J A - J",
        "raw_gram_new_to_old": "A J A^t - J",
        "common_raw_determinant": "2(ad+bc)-1-(det A)^2",
        "mod_two_identity": "G_i = [[0,det(A)-1],[det(A)-1,0]] mod 2",
        "corrected_parity_law": {
            "det_A_even": "the saturated bridge has trivial 2-primary discriminant",
            "det_A_odd_and_m_odd": "the 2-primary discriminant requires two generators",
            "degree_two": "st even kills the 2-primary part; st odd with m odd gives two generators",
        },
        "saturated_odd_prime_counterexample": {
            "A": rows(counterexample_A),
            "d_s_t_z": [2, 1, 4, 1],
            "st_even": True,
            "raw_and_saturated_gram": rows(counterexample_gram),
            "saturation_index": 1,
            "discriminant_group_invariants": smith_invariants(counterexample_gram),
            "conclusion": "Even det(A) does not force the full discriminant group to be cyclic.",
        },
    }


def corpus_controls():
    relative = load_json(RELATIVE_CORPUS)
    glue = load_json(GLUE_CORPUS)
    assert relative["status"] == "PASS_EXACT_RELATIVE_U_BRIDGE_LIFTING_REGRESSION"
    assert glue["status"] == "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES"
    relative_by_edge = {
        (edge["corridor"], int(edge["edge_index"])): edge
        for edge in relative["edges"]
    }
    glue_by_edge = {
        (edge["corridor"], int(edge["edge_index"])): edge
        for edge in glue["edges"]
    }
    assert relative_by_edge.keys() == glue_by_edge.keys()
    assert len(relative_by_edge) == 42

    forced_maximal = []
    shared_bad_primes = []
    orientation_count = 0
    for key in sorted(relative_by_edge):
        relative_edge = relative_by_edge[key]
        glue_edge = glue_by_edge[key]
        old = glue_edge["old_frame"]
        new = glue_edge["new_frame"]
        c_old = int(old["bridge_determinant_absolute"])
        c_new = int(new["bridge_determinant_absolute"])
        h_old = int(old["K_plus_C_index_in_W"])
        h_new = int(new["K_plus_C_index_in_W"])
        k = int(glue_edge["core"]["determinant_absolute"])
        assert c_old == c_new
        assert h_old == h_new
        c_value = c_old
        h_value = h_old
        assert k * c_value % (h_value**2) == 0
        ambient_determinant = k * c_value // (h_value**2)
        assert c_value // h_value == 1

        for direction in ("old_to_new", "new_to_old"):
            orientation_count += 1
            record = relative_edge["relative_u"][direction]
            cross = matrix(ZZ, record["cross_pairing_A"])
            gram = matrix(ZZ, record["positive_projection_gram_G_A"])
            assert cross.det() % 2 == 0
            assert int(record["saturation_index"]) == 1
            assert int(record["saturated_bridge_determinant"]) == c_value
            assert int(gram.det()) == c_value
            assert c_value % 2 == 1
            assert gram.det() % 2 == 1
            assert gcd(gcd(abs(int(gram[0, 0])), abs(int(gram[0, 1]))), abs(int(gram[1, 1]))) == 1

        edge_record = {
            "corridor": key[0],
            "edge_index": key[1],
            "ambient_determinant": ambient_determinant,
            "bridge_determinant": c_value,
            "common_glue_order": h_value,
            "gcd_bridge_ambient": gcd(c_value, ambient_determinant),
        }
        if edge_record["gcd_bridge_ambient"] == 1:
            forced_maximal.append(edge_record)
        else:
            shared_bad_primes.append(edge_record)

    assert orientation_count == 84
    assert len(forced_maximal) == 35
    assert len(shared_bad_primes) == 7
    return {
        "edge_count": 42,
        "orientation_count": 84,
        "all_saturation_indices_one": True,
        "all_relative_determinants_even": True,
        "all_bridge_determinants_odd": True,
        "all_full_discriminant_groups_cyclic_by_exact_smith_test": True,
        "maximal_glue_edges_forced_by_gcd_c_D_one": len(forced_maximal),
        "maximal_glue_edges_requiring_shared_prime_computation": len(shared_bad_primes),
        "shared_bad_prime_edges": shared_bad_primes,
    }


def r17_controls():
    short = load_matrix(SHORT_GRAM)
    short_coordinates = load_matrix(SHORT_COORDS)
    pinned = load_matrix(PINNED_GRAM)
    assert short_coordinates * pinned * short_coordinates.transpose() == short
    assert abs(short_coordinates.det()) == 1
    assert short.det() == pinned.det() == 948
    assert smith_invariants(short) == [948]
    assert all(short[index, index] == 4 for index in range(17))
    assert not pari_short_vectors(short, 2)
    assert len(pari_short_vectors(short, 4)) == 2622

    zero17 = vector(ZZ, [0] * 17)
    v = vector(ZZ, [1] + [0] * 16)
    w = vector(ZZ, [0] * 6 + [1] + [0] * 4 + [-1] + [0] * 5)
    bridge_basis = matrix(ZZ, [v, w])
    bridge_gram = bridge_basis * short * bridge_basis.transpose()
    assert bridge_gram == matrix(ZZ, [[4, 0], [0, 8]])
    assert bridge_basis.row_module(ZZ) == bridge_basis.row_module(ZZ).saturation()

    r = -v + w
    r2 = -2 * v + w
    ambient = block_diagonal_matrix(HYPERBOLIC_PLANE, -short)
    old_u = matrix(ZZ, [[1, 0] + list(zero17), [0, 1] + list(zero17)])
    new_u = matrix(ZZ, [[3, 2] + list(r), [4, 3] + list(r2)])
    new_zero = new_u.row(1) - new_u.row(0)
    assert new_u * ambient * new_u.transpose() == HYPERBOLIC_PLANE
    assert new_zero * ambient * new_zero == -2
    assert new_zero * ambient * new_u.row(0) == 1
    cross = old_u * ambient * new_u.transpose()
    assert cross == matrix(ZZ, [[2, 3], [3, 4]])

    old_raw = matrix(ZZ, [r, r2])
    old_raw_gram = old_raw * short * old_raw.transpose()
    assert old_raw_gram == matrix(ZZ, [[12, 16], [16, 24]])
    assert old_raw.row_module(ZZ).saturation() == bridge_basis.row_module(ZZ)
    assert abs(old_raw.row_module(ZZ).index_in(old_raw.row_module(ZZ).saturation())) == 1

    new_frame_basis = (new_u * ambient).right_kernel_matrix()
    assert abs(new_u.stack(new_frame_basis).det()) == 1
    new_frame = -(new_frame_basis * ambient * new_frame_basis.transpose())
    assert new_frame.is_positive_definite()
    assert new_frame.det() == 948
    assert smith_invariants(new_frame) == [948]

    core_coordinates_old = matrix(ZZ, [r * short, r2 * short]).right_kernel_matrix()
    assert core_coordinates_old.nrows() == 15
    core_basis = matrix(ZZ, [[0, 0] + list(item) for item in core_coordinates_old.rows()])
    core_gram = core_coordinates_old * short * core_coordinates_old.transpose()
    assert core_gram.det() == 30336
    assert not pari_short_vectors(core_gram, 2)
    core_coordinates_new = new_frame_basis.solve_left(core_basis).change_ring(ZZ)
    assert core_coordinates_new * new_frame_basis == core_basis
    new_bridge_coordinates = (core_coordinates_new * new_frame).right_kernel_matrix()
    new_bridge_basis = new_bridge_coordinates * new_frame_basis
    new_bridge_gram = (
        new_bridge_coordinates
        * new_frame
        * new_bridge_coordinates.transpose()
    )
    assert new_bridge_gram.det() == 32
    assert smith_invariants(new_bridge_gram) == [4, 8]
    bridge_isometry = row_isometry(new_bridge_gram, bridge_gram)

    reverse_pairings = old_u * ambient * new_u.transpose()
    reverse_projection = reverse_pairings * HYPERBOLIC_PLANE * new_u
    reverse_residual = old_u - reverse_projection
    assert reverse_residual * ambient * new_u.transpose() == 0
    reverse_raw_coordinates = new_frame_basis.solve_left(reverse_residual).change_ring(ZZ)
    reverse_raw_module = reverse_raw_coordinates.row_module(ZZ)
    reverse_saturation = reverse_raw_module.saturation()
    assert reverse_saturation == new_bridge_coordinates.row_module(ZZ)
    assert abs(reverse_raw_module.index_in(reverse_saturation)) == 1
    reverse_raw_gram = reverse_raw_coordinates * new_frame * reverse_raw_coordinates.transpose()
    assert reverse_raw_gram == cross * HYPERBOLIC_PLANE * cross.transpose() - HYPERBOLIC_PLANE

    old_glue_index = abs(int(core_coordinates_old.stack(bridge_basis).det()))
    new_glue_index = abs(int(core_coordinates_new.stack(new_bridge_coordinates).det()))
    assert old_glue_index == new_glue_index == 32
    assert short.det() == core_gram.det() * bridge_gram.det() // old_glue_index**2
    assert new_frame.det() == core_gram.det() * new_bridge_gram.det() // new_glue_index**2

    residue_generators = {
        "b5": [
            int(short[4, 0] % 4),
            int((short[4, 6] - short[4, 11]) % 8),
        ],
        "b2_plus_2b5": [
            int((short[1, 0] + 2 * short[4, 0]) % 4),
            int(
                (
                    short[1, 6]
                    - short[1, 11]
                    + 2 * (short[4, 6] - short[4, 11])
                )
                % 8
            ),
        ],
    }
    assert residue_generators == {"b5": [1, 0], "b2_plus_2b5": [0, 1]}

    coset_vectors = []
    short_through_eight = pari_short_vectors(short, 8)
    for item in short_through_eight:
        if all((item[index] - r[index]) % 2 == 0 for index in range(17)):
            coset_vectors.append(item)
    assert len(coset_vectors) == 8
    assert {int(item * short * item) for item in coset_vectors} == {8}

    equality_components = []
    for witness in coset_vectors:
        source_part = (r - witness) / 2
        assert source_part in ZZ**17
        alpha = (source_part * short * source_part - 2) / 2
        assert alpha in ZZ
        component = vector(ZZ, [alpha, 1] + list(source_part))
        assert component * ambient * component == -2
        assert component * ambient * new_u.row(0) == 0
        assert component * ambient * old_u.row(0) == 1
        equality_components.append(component)
    component_pairs = []
    for left_index, left in enumerate(equality_components):
        for right_index in range(left_index + 1, len(equality_components)):
            if left + equality_components[right_index] == new_u.row(0):
                component_pairs.append([left_index, right_index])
    assert len(component_pairs) == 4
    assert {index for pair in component_pairs for index in pair} == set(range(8))
    zero_intersections = [
        int(new_zero * ambient * component) for component in equality_components
    ]
    assert sorted(zero_intersections) == [0, 0, 0, 0, 1, 1, 1, 1]
    simple_components = [
        component
        for component, intersection in zip(equality_components, zero_intersections)
        if intersection == 0
    ]
    assert len(simple_components) == 4
    assert (
        matrix(ZZ, simple_components)
        * (-ambient)
        * matrix(ZZ, simple_components).transpose()
        == 2 * matrix.identity(ZZ, 4)
    )

    target_roots = pari_short_vectors(new_frame, 2)
    assert len(target_roots) == 8
    target_root_ambient = {
        tuple(item * new_frame_basis) for item in target_roots
    }
    assert target_root_ambient == {
        tuple(signed)
        for component in simple_components
        for signed in (component, -component)
    }
    simple_coordinates = new_frame_basis.solve_left(
        matrix(ZZ, simple_components)
    ).change_ring(ZZ)
    simple_module = simple_coordinates.row_module(ZZ)
    assert abs(simple_module.index_in(simple_module.saturation())) == 1

    historical_data = load_json(HISTORICAL_MARKING)
    historical_basis = matrix(
        ZZ,
        historical_data["current_suffix_stages"]["current_4A1"][
            "basis_in_pinned_R17"
        ],
    )
    pinned_ambient = block_diagonal_matrix(HYPERBOLIC_PLANE, -pinned)
    historical_split = historical_basis * pinned_ambient * historical_basis.transpose()
    assert historical_split[:2, :2] == HYPERBOLIC_PLANE
    assert not historical_split[:2, 2:]
    historical_frame = -historical_split[2:, 2:]
    physical_frame = load_matrix(PHYSICAL_Q8_FRAME)

    frame_comparisons = {}
    for name, candidate in (
        ("new", new_frame),
        ("historical_current_4A1", historical_frame),
        ("physical_q8_orbit376_4A1", physical_frame),
    ):
        norm_two_signed = len(pari_short_vectors(candidate, 2))
        norm_at_most_four_signed = len(pari_short_vectors(candidate, 4))
        frame_comparisons[name] = {
            "determinant": int(candidate.det()),
            "norm_two_signed": norm_two_signed,
            "norm_four_pairs": (norm_at_most_four_signed - norm_two_signed) // 2,
            "automorphism_group_order": int(pari(candidate).qfauto()[0]),
        }
    assert frame_comparisons == {
        "new": {
            "determinant": 948,
            "norm_two_signed": 8,
            "norm_four_pairs": 1301,
            "automorphism_group_order": 32,
        },
        "historical_current_4A1": {
            "determinant": 948,
            "norm_two_signed": 8,
            "norm_four_pairs": 1263,
            "automorphism_group_order": 32,
        },
        "physical_q8_orbit376_4A1": {
            "determinant": 948,
            "norm_two_signed": 8,
            "norm_four_pairs": 1337,
            "automorphism_group_order": 64,
        },
    }
    historical_isometric = QuadraticForm(ZZ, new_frame).is_globally_equivalent_to(
        QuadraticForm(ZZ, historical_frame)
    )
    physical_isometric = QuadraticForm(ZZ, new_frame).is_globally_equivalent_to(
        QuadraticForm(ZZ, physical_frame)
    )
    assert historical_isometric is False
    assert physical_isometric is False

    return {
        "source_R17": {
            "determinant": 948,
            "discriminant_group_invariants": [948],
            "norm_four_pairs": 1311,
            "short_vector_coordinate_basis_determinant": int(short_coordinates.det()),
        },
        "relative_U": {
            "cross_pairing_A": rows(cross),
            "d_s_t_z": [2, 1, 1, 0],
            "det_A": int(cross.det()),
            "target_U_basis_in_U_plus_short_R17": rows(new_u),
            "old_raw_gram": rows(old_raw_gram),
            "new_raw_gram": rows(reverse_raw_gram),
            "saturation_indices": [1, 1],
        },
        "bridge": {
            "C0_gram": rows(bridge_gram),
            "C1_gram": rows(new_bridge_gram),
            "C1_basis_in_ambient_NS": rows(new_bridge_basis),
            "integral_isometry_C1_to_C0": rows(bridge_isometry),
            "common_determinant": 32,
            "discriminant_group_invariants": [4, 8],
            "core_rank": 15,
            "core_determinant": 30336,
            "core_gram": rows(core_gram),
            "core_root_count_signed": 0,
            "glue_orders": [old_glue_index, new_glue_index],
            "glue_is_maximal": True,
            "ambient_bridge_gcd": gcd(948, 32),
            "projection_to_A_C_generators": residue_generators,
        },
        "geometric_gate": {
            "coset": "r + 2 R17",
            "coset_minimum": 8,
            "norm_eight_witnesses": len(coset_vectors),
            "norm_eight_witnesses_in_short_R17": [
                list(map(int, item)) for item in coset_vectors
            ],
            "no_norm_four_or_six_witnesses": True,
            "new_fibre_nef_in_old_chamber": True,
            "vertical_root_count_signed": len(equality_components),
            "frame_root_system": "4A1",
            "root_span_saturation_index": 1,
            "new_zero_square": -2,
            "new_zero_fibre_intersection": 1,
            "vertical_component_intersections_with_new_zero": sorted(
                zero_intersections
            ),
            "simple_components_in_ambient_NS": [
                list(map(int, item)) for item in simple_components
            ],
            "new_zero_is_physical": True,
            "mw_rank": 13,
            "mw_torsion_order": 1,
        },
        "target_frame": {
            "basis_in_ambient_NS": rows(new_frame_basis),
            "gram": rows(new_frame),
            "determinant": 948,
            "discriminant_group_invariants": [948],
        },
        "frame_comparisons": frame_comparisons,
        "new_is_historical_current_4A1": False,
        "new_is_physical_q8_orbit376_4A1": False,
        "conclusion": "This is a distinct 4A1/MW13 J2 frame, not a shortcut to either stored H3 4A1 node.",
    }


def build_result():
    inputs = [
        SHORT_GRAM,
        SHORT_COORDS,
        PINNED_GRAM,
        RELATIVE_CORPUS,
        GLUE_CORPUS,
        HISTORICAL_MARKING,
        PHYSICAL_Q8_FRAME,
    ]
    return {
        "schema": "elkies-k3.r17-local-bridge-mutation.v1",
        "status": "PASS_EXACT_R17_LOCAL_BRIDGE_MUTATION",
        "theorem": {
            "name": "Local Bridge-Mutation Theorem H-1c",
            "scope": "rank-two relative U position, two-sided saturation, glue defect support, and 2-primary parity",
        },
        "symbolic_controls": symbolic_controls(),
        "corpus_controls": corpus_controls(),
        "r17_example": r17_controls(),
        "inputs": {
            "paths": [str(item) for item in inputs],
            "sha256": {str(item): sha256(item) for item in inputs},
        },
        "software": {
            "sage": SAGE_VERSION,
            "arithmetic": "exact ZZ lattice arithmetic with PARI qfminim/qfauto and Sage integral isometry",
        },
        "proof_boundary": (
            "The structural theorem is integral lattice algebra.  The checker gives "
            "exact symbolic identities, a complete replay of the stored 42-edge "
            "corpus, and an exact lattice/chamber/divisor certificate for the new "
            "R17 degree-two fibration.  It does not construct a Weierstrass equation, "
            "identify a J1 surface-automorphism orbit, or determine Galois action."
        ),
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_local_bridge_mutation.sage --check"
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_result()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        assert output.read_text() == encoded, f"stale artifact: {output}"
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "R17BRIDGE|edges=42|forced_maximal=35|"
        "bridge=Z/4+Z/8|roots=4A1|mw=13|historical=false|"
        f"status={result['status']}|output={output}"
    )


if __name__ == "__main__":
    main()
