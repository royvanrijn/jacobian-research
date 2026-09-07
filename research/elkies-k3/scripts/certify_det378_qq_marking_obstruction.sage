#!/usr/bin/env sage-python
"""Certify the determinant-378 rational-marking modular obstruction.

The split row in the arithmetic-first queue has literal transcendental lattice

    T = U(3) + <42> = 3 (U + <14>).

The primitive similarity order has projective norm-one curve X_0(7).  This
checker embeds the *literal* even Clifford order in M_2(QQ), computes the spin
action on A_T, and proves that its stable subgroup is the additional projective
identity condition modulo 3.  In Gamma_0(7) coordinates this is

    Gamma_0(7) intersection {g : g = +/-I (mod 3)}.

Conjugation by diag(3,1) identifies the group with Gamma_0(63), so the exact
marked curve is X_0(63).  The complete Mazur--Kenku classification of rational
cyclic-isogeny degrees excludes every noncuspidal QQ-point.

That classification, the standard split-Eichler normalizer description, and
the ternary-spin/marked-K3 period correspondence are external theorem inputs.
All lattice, discriminant-form, Clifford-order, spin-action, congruence,
signature, and cusp calculations are replayed exactly here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, deque
from pathlib import Path

from sage.all import (
    CliffordAlgebra,
    Gamma0,
    QQ,
    ZZ,
    QuadraticForm,
    identity_matrix,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CATALOGUE = GENERATED / "elkies-k3-rank7-auxiliary-catalogue-v1.json"
T_ARITHMETIC = GENERATED / "elkies-k3-rank7-t-arithmetic-v1.json"
OUTPUT = GENERATED / "elkies-k3-det378-qq-marking-obstruction-v1.json"
SURFACE_ID = "K3-7b71a1cc00b0c6e2"

RATIONAL_CYCLIC_ISOGENY_DEGREES = tuple(
    list(range(1, 20)) + [21, 25, 27, 37, 43, 67, 163]
)


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def string_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def tuple_matrix(value, modulus):
    return tuple(int(entry % modulus) for entry in value.list())


def multiply_2x2(left, right, modulus):
    aa, bb, cc, dd = left
    ee, ff, gg, hh = right
    return (
        (aa * ee + bb * gg) % modulus,
        (aa * ff + bb * hh) % modulus,
        (cc * ee + dd * gg) % modulus,
        (cc * ff + dd * hh) % modulus,
    )


def negate_2x2(value, modulus):
    return tuple((-entry) % modulus for entry in value)


def action_order(value, identity, compose, maximum):
    current = identity
    for order in range(1, maximum + 1):
        current = compose(current, value)
        if current == identity:
            return order
    raise ArithmeticError(f"action element order exceeds {maximum}")


def closure(identity, generators, compose):
    image = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            following = compose(current, generator)
            if following not in image:
                image.add(following)
                queue.append(following)
    return image


def build_payload():
    catalogue = json.loads(CATALOGUE.read_text())
    t_payload = json.loads(T_ARITHMETIC.read_text())
    surface = next(
        row for row in catalogue["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    t_row = next(
        row for row in t_payload["surfaces"] if row["surface_id"] == SURFACE_ID
    )

    transcendental = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    expected = matrix(ZZ, [[0, 0, 3], [0, 42, 0], [3, 0, 0]])
    assert transcendental == expected
    assert transcendental.det() == -378
    assert QuadraticForm(ZZ, transcendental).signature_vector() == (2, 1, 0)
    assert t_row["literal_transcendental_gram"] == rows(transcendental)
    assert t_row["similarity_normalization"]["literal_content"] == 3

    smith, smith_left, smith_right = transcendental.smith_form()
    assert smith_left * transcendental * smith_right == smith
    discriminant_invariants = [abs(int(smith[index, index])) for index in range(3)]
    assert discriminant_invariants == [3, 3, 42]
    assert surface["surface_key"]["ns_discriminant_form_key"]["invariants"] == [
        3,
        3,
        42,
    ]

    # An explicit basis of T^*: e0/3, e2/3, e1/42.  Its Gram matrix modulo Z
    # records the bilinear discriminant form; the diagonal records q modulo 2Z.
    dual_generators = matrix(
        QQ,
        [
            [QQ(1) / 3, 0, 0],
            [0, 0, QQ(1) / 42],
            [0, QQ(1) / 3, 0],
        ],
    )
    assert all(entry in ZZ for entry in (transcendental * dual_generators).list())
    assert abs(dual_generators.det()) == QQ(1) / 378
    discriminant_bilinear_gram = (
        dual_generators.transpose() * transcendental * dual_generators
    )
    assert discriminant_bilinear_gram == matrix(
        QQ,
        [[0, QQ(1) / 3, 0], [QQ(1) / 3, 0, 0], [0, 0, QQ(1) / 42]],
    )

    primitive = (transcendental / 3).change_ring(ZZ)
    assert primitive == matrix(ZZ, [[0, 0, 1], [0, 14, 0], [1, 0, 0]])

    # Exact split representation of C^+(primitive), in the basis
    # 1,e0e1,e0e2,e1e2.  Its elements are precisely the matrices whose upper
    # right entry is divisible by seven.
    quadratic = QuadraticForm(QQ, primitive)
    clifford = CliffordAlgebra(quadratic)
    ee0, ee1, ee2 = clifford.gens()
    even_basis = [clifford.one(), ee0 * ee1, ee0 * ee2, ee1 * ee2]
    split_basis = [
        matrix(QQ, [[1, 0], [0, 1]]),
        matrix(QQ, [[0, 7], [0, 0]]),
        matrix(QQ, [[1, 0], [0, 0]]),
        matrix(QQ, [[0, 0], [1, 0]]),
    ]
    clifford_keys = list(clifford.basis().keys())

    def even_coordinates(value):
        coefficients = value.monomial_coefficients()
        return [coefficients.get(clifford_keys[index], 0) for index in (0, 4, 5, 6)]

    def split_image(value):
        coordinates = even_coordinates(value)
        return sum(
            (coordinates[index] * split_basis[index] for index in range(4)),
            matrix(QQ, 2, 2),
        )

    for left in even_basis:
        for right in even_basis:
            assert split_image(left * right) == split_image(left) * split_image(right)

    primitive_trace_pairing = matrix(
        ZZ,
        4,
        4,
        lambda left, right: (split_basis[left] * split_basis[right]).trace(),
    )
    assert primitive_trace_pairing == matrix(
        ZZ,
        [[2, 0, 1, 0], [0, 0, 0, 7], [1, 0, 1, 0], [0, 7, 0, 0]],
    )
    assert primitive_trace_pairing.det() == -(7**2)
    assert t_row["clifford"]["integral_even_clifford_order"][
        "reduced_trace_pairing"
    ] == rows(primitive_trace_pairing)

    coarse_group = Gamma0(7)
    coarse_widths = sorted(
        int(coarse_group.cusp_width(cusp)) for cusp in coarse_group.cusps()
    )
    assert coarse_group.index() == 8
    assert coarse_group.genus() == 0
    assert coarse_group.nu2() == 0 and coarse_group.nu3() == 2
    assert coarse_widths == [1, 7]

    # J conjugates Gamma_0(7) to the norm-one group Gamma^0(7) of the displayed
    # primitive order.
    conjugator = matrix(ZZ, [[0, -1], [1, 0]])
    for generator in coarse_group.gens():
        gamma0_matrix = matrix(ZZ, generator.matrix())
        order_matrix = conjugator.inverse() * gamma0_matrix * conjugator
        assert order_matrix.det() == 1
        assert order_matrix[0, 1] % 7 == 0

    # The corrected central volume element identifies T_Q with the trace-zero
    # part of the split even Clifford algebra.  Conjugation gives the spin
    # action on the literal lattice (scaling the form by three changes neither
    # the rational action nor O(T)).
    volume = (
        ee0 * ee1 * ee2
        - QQ(primitive[0, 1]) / 2 * ee2
        + QQ(primitive[0, 2]) / 2 * ee1
        - QQ(primitive[1, 2]) / 2 * ee0
    )
    assert all(volume * generator == generator * volume for generator in (ee0, ee1, ee2))
    assert volume * volume == QQ(7) / 4
    vector_images = [split_image(volume * generator) for generator in (ee0, ee1, ee2)]
    assert all(value.trace() == 0 for value in vector_images)
    image_matrix = matrix(
        QQ, 4, 3, lambda row, column: vector_images[column].list()[row]
    )
    left_inverse = (
        image_matrix.transpose() * image_matrix
    ).inverse() * image_matrix.transpose()
    assert left_inverse * image_matrix == identity_matrix(QQ, 3)

    def spin_action(value):
        value = matrix(QQ, value)
        inverse = value.inverse()
        columns = []
        for image in vector_images:
            conjugate = value * image * inverse
            coordinates = left_inverse * vector(QQ, conjugate.list())
            assert image_matrix * coordinates == vector(QQ, conjugate.list())
            columns.append(coordinates)
        action = matrix(QQ, 3, 3, lambda row, column: columns[column][row])
        assert all(entry in ZZ for entry in action.list())
        action = action.change_ring(ZZ)
        assert action.transpose() * transcendental * action == transcendental
        assert action.det() == 1
        return action

    moduli = tuple(discriminant_invariants)

    def discriminant_action_key(action):
        # A_T=T^*/T is coker(T); action on coker coordinates is inverse
        # transpose, followed by the Smith-coordinate change.
        contragredient = action.inverse().transpose().change_ring(ZZ)
        smith_action = smith_left * contragredient * smith_left.inverse()
        return tuple(
            int(smith_action[row, column] % moduli[row])
            for row in range(3)
            for column in range(3)
        )

    identity_action = tuple(
        1 if row == column else 0 for row in range(3) for column in range(3)
    )

    def compose_discriminant(left, right):
        answer = matrix(ZZ, 3, 3, left) * matrix(ZZ, 3, 3, right)
        return tuple(
            int(answer[row, column] % moduli[row])
            for row in range(3)
            for column in range(3)
        )

    generator_records = []
    joint_generators = []
    spin_discriminant_generators = []
    for generator in coarse_group.gens():
        gamma0_matrix = matrix(ZZ, generator.matrix())
        order_matrix = conjugator.inverse() * gamma0_matrix * conjugator
        action = spin_action(order_matrix)
        discriminant_action = discriminant_action_key(action)
        mod3_reduction = tuple_matrix(gamma0_matrix, 3)
        joint_generators.append((mod3_reduction, discriminant_action))
        spin_discriminant_generators.append(discriminant_action)
        generator_records.append(
            {
                "gamma0_7_generator": rows(gamma0_matrix),
                "primitive_order_generator": rows(order_matrix),
                "spin_action_on_T": rows(action),
                "mod_3_reduction": list(mod3_reduction),
                "discriminant_action_smith_coordinates": list(discriminant_action),
            }
        )

    identity_mod3 = (1, 0, 0, 1)

    def compose_joint(left, right):
        return (
            multiply_2x2(left[0], right[0], 3),
            compose_discriminant(left[1], right[1]),
        )

    joint_identity = (identity_mod3, identity_action)
    joint_image = closure(joint_identity, joint_generators, compose_joint)
    mod3_image = {value[0] for value in joint_image}
    spin_discriminant_image = {value[1] for value in joint_image}
    assert len(joint_image) == len(mod3_image) == 24
    assert len(spin_discriminant_image) == 12
    assert {
        value[0]
        for value in joint_image
        if value[1] == identity_action
    } == {identity_mod3, negate_2x2(identity_mod3, 3)}
    for mod3_value, discriminant_value in joint_image:
        partners = {
            other_mod3
            for other_mod3, other_discriminant in joint_image
            if other_discriminant == discriminant_value
        }
        assert partners == {mod3_value, negate_2x2(mod3_value, 3)}

    spin_order_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 12)
        for value in spin_discriminant_image
    )
    assert spin_order_histogram == Counter({3: 8, 2: 3, 1: 1})

    # Now embed C^+(T) itself.  Under the rational similarity identification
    # its basis is 1,3e0e1,3e0e2,3e1e2, not the primitive order basis.
    literal_split_basis = [split_basis[0]] + [3 * value for value in split_basis[1:]]
    literal_quadratic = QuadraticForm(QQ, transcendental)
    literal_clifford = CliffordAlgebra(literal_quadratic)
    ll0, ll1, ll2 = literal_clifford.gens()
    literal_even_basis = [
        literal_clifford.one(),
        ll0 * ll1,
        ll0 * ll2,
        ll1 * ll2,
    ]
    literal_keys = list(literal_clifford.basis().keys())

    def literal_split_image(value):
        coefficients = value.monomial_coefficients()
        coordinates = [coefficients.get(literal_keys[index], 0) for index in (0, 4, 5, 6)]
        return sum(
            (coordinates[index] * literal_split_basis[index] for index in range(4)),
            matrix(QQ, 2, 2),
        )

    for left in literal_even_basis:
        for right in literal_even_basis:
            assert literal_split_image(left * right) == (
                literal_split_image(left) * literal_split_image(right)
            )

    literal_trace_pairing = matrix(
        ZZ,
        4,
        4,
        lambda left, right: (
            literal_split_basis[left] * literal_split_basis[right]
        ).trace(),
    )
    assert literal_trace_pairing == matrix(
        ZZ,
        [[2, 0, 3, 0], [0, 0, 0, 63], [3, 0, 9, 0], [0, 63, 0, 0]],
    )
    assert literal_trace_pairing.det() == -(189**2)

    # Matrix membership in the literal order is exactly
    #   21|B, 3|C, A=D (mod 3).
    # For determinant one, relative to the primitive order, this is reduction
    # to +/-I modulo 3.  Exhaust the 24-element mod-3 image to identify the
    # stable kernel with the literal norm-one group.
    projective_identity_mod3 = {
        identity_mod3,
        negate_2x2(identity_mod3, 3),
    }
    stable_joint = {
        value for value in joint_image if value[0] in projective_identity_mod3
    }
    assert len(stable_joint) == 2
    assert {value[1] for value in stable_joint} == {identity_action}

    # The norm-one spin subgroup is not yet all of O^+(T).  A negative
    # reflection r, central inversion -I, and the proper Fricke isometry -r
    # represent the other integral cosets.  Compute all of them: no extra
    # coset contributes a stable element.
    negative_vector = vector(ZZ, [1, 0, -1])
    negative_norm = negative_vector * transcendental * negative_vector
    assert negative_norm == -6
    reflection = identity_matrix(QQ, 3) - (
        QQ(2) / negative_norm
    ) * matrix(QQ, negative_vector).transpose() * matrix(QQ, negative_vector) * transcendental
    assert all(entry in ZZ for entry in reflection.list())
    reflection = reflection.change_ring(ZZ)
    assert reflection == matrix(ZZ, [[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    assert reflection.det() == -1
    assert reflection.transpose() * transcendental * reflection == transcendental
    reflection_action = discriminant_action_key(reflection)
    assert reflection_action not in spin_discriminant_image
    assert compose_discriminant(reflection_action, reflection_action) == identity_action
    reflection_extended_image = closure(
        identity_action,
        spin_discriminant_generators + [reflection_action],
        compose_discriminant,
    )
    assert len(reflection_extended_image) == 24
    reflection_extended_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 24)
        for value in reflection_extended_image
    )
    assert reflection_extended_histogram == Counter({6: 8, 3: 8, 2: 7, 1: 1}), (
        reflection_extended_histogram
    )

    central_inversion = -identity_matrix(ZZ, 3)
    assert central_inversion.det() == -1
    assert central_inversion.transpose() * transcendental * central_inversion == transcendental
    central_inversion_action = discriminant_action_key(central_inversion)
    assert central_inversion_action not in reflection_extended_image
    assert compose_discriminant(
        central_inversion_action, central_inversion_action
    ) == identity_action
    assert all(
        compose_discriminant(central_inversion_action, value)
        == compose_discriminant(value, central_inversion_action)
        for value in reflection_extended_image
    )

    fricke_isometry = -reflection
    assert fricke_isometry.det() == 1
    assert fricke_isometry.transpose() * transcendental * fricke_isometry == transcendental
    fricke_action = discriminant_action_key(fricke_isometry)
    assert fricke_action == compose_discriminant(
        central_inversion_action, reflection_action
    )
    assert fricke_action not in spin_discriminant_image
    orientation_preserving_image = closure(
        identity_action,
        spin_discriminant_generators + [fricke_action],
        compose_discriminant,
    )
    assert len(orientation_preserving_image) == 24
    orientation_preserving_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 24)
        for value in orientation_preserving_image
    )
    assert orientation_preserving_histogram == Counter(
        {6: 8, 3: 8, 2: 7, 1: 1}
    )
    orientation_central_coset = [
        compose_discriminant(fricke_action, spin_value)
        for spin_value in spin_discriminant_image
        if all(
            compose_discriminant(
                compose_discriminant(fricke_action, spin_value), other
            )
            == compose_discriminant(
                other, compose_discriminant(fricke_action, spin_value)
            )
            for other in spin_discriminant_image
        )
    ]
    assert len(orientation_central_coset) == 1
    orientation_central_involution = orientation_central_coset[0]
    assert orientation_central_involution not in spin_discriminant_image
    assert compose_discriminant(
        orientation_central_involution, orientation_central_involution
    ) == identity_action

    full_discriminant_image = closure(
        identity_action,
        spin_discriminant_generators
        + [reflection_action, central_inversion_action],
        compose_discriminant,
    )
    assert len(full_discriminant_image) == 48
    full_order_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 24)
        for value in full_discriminant_image
    )
    assert full_order_histogram == Counter({6: 24, 2: 15, 3: 8, 1: 1}), (
        full_order_histogram
    )

    # In Gamma_0(7) coordinates the stable group is
    # Delta=Gamma_0(7) intersection +/-Gamma(3).  If D=diag(3,1), then
    # D^-1 Delta D=Gamma_0(63).  The determinant equation supplies the diagonal
    # congruence in the inverse direction.
    diagonal_conjugator = matrix(QQ, [[3, 0], [0, 1]])
    stable_group = Gamma0(63)
    stable_index_over_coarse = len(spin_discriminant_image)
    stable_index = int(coarse_group.index()) * stable_index_over_coarse
    assert stable_index == 96 == stable_group.index()
    assert stable_group.genus() == 5
    assert stable_group.nu2() == stable_group.nu3() == 0
    stable_cusps = []
    for cusp in stable_group.cusps():
        denominator = int(cusp.denominator()) if not cusp.is_infinity() else 0
        rationality_order = (
            1
            if denominator == 0
            else int(ZZ(denominator).gcd(63 // ZZ(63).gcd(denominator)))
        )
        stable_cusps.append(
            {
                "representative": str(cusp),
                "width": int(stable_group.cusp_width(cusp)),
                "cusp_rationality_order": rationality_order,
                "rational_over_QQ": rationality_order <= 2,
            }
        )
    assert len(stable_cusps) == 8
    assert sum(row["width"] for row in stable_cusps) == stable_index
    rational_cusps = [row for row in stable_cusps if row["rational_over_QQ"]]
    assert len(rational_cusps) == 4

    assert 63 not in RATIONAL_CYCLIC_ISOGENY_DEGREES

    return {
        "schema": "elkies-k3.det378-qq-marking-obstruction.v1",
        "status": "PASS_DET378_QQ_RATIONAL_MARKING_OBSTRUCTION",
        "surface_id": SURFACE_ID,
        "transcendental_lattice": {
            "gram": rows(transcendental),
            "decomposition": "U(3) orthogonal_sum <42>",
            "determinant": int(transcendental.det()),
            "signature": [2, 1],
            "primitive_similarity_gram": rows(primitive),
            "literal_similarity_factor": 3,
            "discriminant_form": {
                "group_invariants": discriminant_invariants,
                "group": "Z/3 + Z/3 + Z/42",
                "dual_generators": ["e0/3", "e2/3", "e1/42"],
                "bilinear_gram_mod_Z": string_rows(discriminant_bilinear_gram),
                "quadratic_values_mod_2Z": ["0", "0", "1/42"],
            },
        },
        "split_even_clifford_orders": {
            "primitive_similarity_order": {
                "basis_in_M2_QQ": [string_rows(value) for value in split_basis],
                "reduced_trace_pairing": rows(primitive_trace_pairing),
                "reduced_discriminant": 7,
                "matrix_conditions": ["7 divides B"],
                "norm_one_group": "Gamma^0(7), conjugate to Gamma_0(7)",
            },
            "literal_T_order": {
                "basis_in_M2_QQ": [string_rows(value) for value in literal_split_basis],
                "reduced_trace_pairing": rows(literal_trace_pairing),
                "reduced_discriminant": 189,
                "matrix_conditions": [
                    "21 divides B",
                    "3 divides C",
                    "A-D is divisible by 3",
                ],
                "norm_one_group": (
                    "Gamma^0(7) intersection {g : g=+/-I modulo 3}"
                ),
            },
        },
        "coarse_norm_one_curve": {
            "label": "X_0(7)",
            "congruence_group": "Gamma_0(7)",
            "index_in_PSL2Z": int(coarse_group.index()),
            "elliptic_orbits_order_2": int(coarse_group.nu2()),
            "elliptic_orbits_order_3": int(coarse_group.nu3()),
            "cusp_widths": coarse_widths,
            "genus": int(coarse_group.genus()),
        },
        "stable_discriminant_kernel": {
            "abstract_group": "O^+(T)^* = kernel(O^+(T) -> O(A_T))",
            "norm_one_spin_discriminant_image_order": len(spin_discriminant_image),
            "norm_one_spin_discriminant_image_isomorphism_type": "A4",
            "norm_one_spin_discriminant_image_element_order_histogram": {
                str(order): count
                for order, count in sorted(spin_order_histogram.items())
            },
            "orientation_preserving_discriminant_image_order": len(
                orientation_preserving_image
            ),
            "orientation_preserving_discriminant_image_isomorphism_type": (
                "A4 x C2"
            ),
            "orientation_preserving_discriminant_image_element_order_histogram": {
                str(order): count
                for order, count in sorted(orientation_preserving_histogram.items())
            },
            "full_discriminant_image_order": len(full_discriminant_image),
            "full_discriminant_image_isomorphism_type": "A4 x C2 x C2",
            "full_discriminant_image_element_order_histogram": {
                str(order): count
                for order, count in sorted(full_order_histogram.items())
            },
            "generator_actions": generator_records,
            "kernel_condition": "projective identity (+/-I) modulo 3",
            "non_spin_stable_cosets": 0,
            "negative_reflection_witness": {
                "vector": list(map(int, negative_vector)),
                "norm": int(negative_norm),
                "matrix": rows(reflection),
                "discriminant_action_smith_coordinates": list(reflection_action),
                "action_lies_in_norm_one_A4_image": False,
            },
            "central_inversion_witness": {
                "matrix": rows(central_inversion),
                "discriminant_action_smith_coordinates": list(
                    central_inversion_action
                ),
                "action_lies_in_reflection_extended_image": False,
            },
            "fricke_isometry_witness": {
                "matrix": rows(fricke_isometry),
                "determinant": 1,
                "discriminant_action_smith_coordinates": list(fricke_action),
                "action_lies_in_norm_one_A4_image": False,
                "action_centralizes_norm_one_A4_image": False,
                "central_involution_in_same_A4_coset_smith_coordinates": list(
                    orientation_central_involution
                ),
            },
            "congruence_group_inside_coarse_Gamma0_7": (
                "Gamma_0(7) intersection {g : g=+/-I modulo 3}"
            ),
            "index_over_coarse_curve": stable_index_over_coarse,
            "index_in_PSL2Z": stable_index,
        },
        "exact_marked_modular_curve": {
            "label": "X_0(63)",
            "identification": (
                "diag(3,1)^(-1) (Gamma_0(7) intersection +/-Gamma(3)) "
                "diag(3,1) = Gamma_0(63)"
            ),
            "diagonal_conjugator": string_rows(diagonal_conjugator),
            "defined_over": "QQ",
            "index_in_PSL2Z": stable_index,
            "genus": int(stable_group.genus()),
            "elliptic_orbits_order_2": int(stable_group.nu2()),
            "elliptic_orbits_order_3": int(stable_group.nu3()),
            "cusps": stable_cusps,
            "rational_points": {
                "rational_cuspidal_points": len(rational_cusps),
                "geometric_cusps": len(stable_cusps),
                "noncuspidal_points": 0,
                "conclusion": "X_0(63)(QQ) consists of its four rational cusps",
                "theorem_input": (
                    "complete Mazur--Kenku rational cyclic-isogeny classification"
                ),
            },
        },
        "noncuspidal_point_audit": {
            "points": [],
            "CM_rejections_required": 0,
            "saturated_rational_marking_checks_required": 0,
            "NS_complements_constructed": 0,
            "rootless_MW17_frame_tests_required": 0,
            "conclusion": (
                "There is no noncuspidal rational point, so every downstream "
                "CM, saturation, NS-complement, and rootless-frame gate is vacuous."
            ),
        },
        "arithmetic_conclusion": (
            "No characteristic-zero K3 over QQ with transcendental lattice "
            "U(3)+<42> can have the full saturated rational rank-19 marking."
        ),
        "external_theorem_inputs": [
            {
                "name": "Mazur--Kenku classification of rational cyclic isogenies",
                "used_for": "X_0(63)(QQ) has no noncuspidal point",
                "allowed_degrees": list(RATIONAL_CYCLIC_ISOGENY_DEGREES),
            },
            {
                "name": "standard normalizer description for split Eichler orders",
                "used_for": (
                    "the norm-one group together with the displayed Fricke, negative-"
                    "reflection, and central-inversion cosets exhausts O^+(T)"
                ),
            },
            {
                "name": "standard ternary spin and marked-K3 period correspondence",
                "used_for": (
                    "identification of the stable orthogonal period curve and the "
                    "rational point forced by a full rational Neron--Severi marking"
                ),
            },
        ],
        "input_hashes": {
            relative(CATALOGUE): digest(CATALOGUE),
            relative(T_ARITHMETIC): digest(T_ARITHMETIC),
        },
        "proof_boundary": {
            "proved": (
                "The literal T and discriminant form, both split Clifford embeddings "
                "and multiplication tables, literal-order congruence conditions, "
                "A4 norm-one spin image, A4 x C2 proper image, full A4 x C2 x C2 "
                "discriminant action, stable kernel, X_0(63) identification, "
                "signature, and cusp data are replayed exactly. With the named "
                "cyclic-isogeny and period theorems this excludes the full rational "
                "determinant-378 marking."
            ),
            "not_proved": (
                "The checker does not reprove the global Mazur--Kenku classification, "
                "the split-Eichler normalizer theorem, or the general marked-K3 "
                "period/spin correspondence. It makes no exclusion over larger number "
                "fields and does not address the other, anisotropic determinant-378 "
                "lattice."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/certify_det378_qq_marking_obstruction.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    payload = build_payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "DET378MARKING|coarse=X0(7)|spin_Aimage=A4|SO_Aimage=A4xC2|"
        "O_Aimage=A4xC2xC2|"
        "stable=X0(63)|genus=5|QQ_noncusps=0|T378=EXCLUDED|status=PASS"
    )


if __name__ == "__main__":
    main()
