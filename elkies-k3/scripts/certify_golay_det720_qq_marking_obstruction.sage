#!/usr/bin/env sage-python
"""Certify the determinant-720 rational-marking modular obstruction.

For the Golay rootless frame the literal transcendental lattice has Gram

    [ 8  -2   2]
    [-2  -4  10].
    [ 2  10  -4]

The script embeds both its primitive-similarity and literal even Clifford
orders in M_2(QQ).  The first norm-one group is conjugate to Gamma_0(15).
Its action on A_T = Z/2 + Z/6 + Z/60 has image S_3, and the kernel is the
additional principal mod-2 condition.  The literal order realizes that kernel
directly.  After a rational diagonal conjugation the stable marked curve is
X_0(60), which has no noncuspidal QQ-point by the complete Mazur--Kenku
classification of rational cyclic-isogeny degrees.

The classification and the standard ternary spin / marked-K3 period
correspondence are external theorem inputs.  All lattice, Clifford-order,
congruence, discriminant-action, signature, cusp, and quotient-point checks
are replayed exactly here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, deque
from itertools import product
from math import gcd
from pathlib import Path

from sage.all import (
    CliffordAlgebra,
    EllipticCurve_from_j,
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
SATURATION = GENERATED / "elkies-k3-golay-det720-3a5-saturation-rejection-v1.json"
OUTPUT = GENERATED / "elkies-k3-golay-det720-qq-marking-obstruction-v1.json"
SURFACE_ID = "K3-f43753fb154e3406"

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


def determinant_2x2(value, modulus):
    aa, bb, cc, dd = value
    return (aa * dd - bb * cc) % modulus


def sl2(modulus):
    return [
        value
        for value in product(range(modulus), repeat=4)
        if determinant_2x2(value, modulus) == 1 % modulus
    ]


def left_cosets(group, subgroup, modulus):
    unseen = set(group)
    representatives = []
    lookup = {}
    while unseen:
        representative = min(unseen)
        coset = {
            multiply_2x2(element, representative, modulus)
            for element in subgroup
        }
        index = len(representatives)
        representatives.append(representative)
        for element in coset:
            lookup[element] = index
        unseen -= coset
    return representatives, lookup


def canonical_projective_row(value, prime):
    left, right = value
    if left % prime:
        return (1, (right * pow(left, -1, prime)) % prime)
    return (0, 1)


def cycle_lengths(permutation):
    seen = set()
    result = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        result.append(length)
    return sorted(result)


def build_payload():
    catalogue = json.loads(CATALOGUE.read_text())
    t_payload = json.loads(T_ARITHMETIC.read_text())
    saturation = json.loads(SATURATION.read_text())

    surface = next(
        row for row in catalogue["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    t_row = next(
        row for row in t_payload["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    transcendental = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    expected = matrix(ZZ, [[8, -2, 2], [-2, -4, 10], [2, 10, -4]])
    assert transcendental == expected
    assert transcendental.det() == -720
    smith, smith_left, smith_right = transcendental.smith_form()
    assert smith_left * transcendental * smith_right == smith
    discriminant_invariants = [abs(int(smith[index, index])) for index in range(3)]
    assert discriminant_invariants == [2, 6, 60]
    assert t_row["literal_transcendental_gram"] == rows(transcendental)
    assert t_row["clifford"]["quaternion_discriminant"] == 1
    assert t_row["clifford"]["integral_even_clifford_order"][
        "reduced_discriminant"
    ] == 45

    primitive = transcendental / 2
    assert primitive.change_ring(ZZ) == matrix(
        ZZ, [[4, -1, 1], [-1, -2, 5], [1, 5, -2]]
    )
    primitive = primitive.change_ring(ZZ)

    # Rational U-coordinates for the primitive similarity form.  The columns
    # are two isotropic vectors pairing to one and their orthogonal complement.
    isotropic = vector(QQ, [1, 1, 0])
    companion = vector(QQ, [1, -2, 0]) / 9
    complement = vector(QQ, [1, 7, 3])
    u_change = matrix(QQ, [isotropic, companion, complement]).transpose()
    u_gram = u_change.transpose() * primitive * u_change
    assert u_gram == matrix(QQ, [[0, 1, 0], [1, 0, 0], [0, 0, 90]])

    # Exact split representation of C^+(primitive).  The basis is
    # 1,e0e1,e0e2,e1e2.  It is the rational representation induced by the
    # displayed U-coordinates, followed by conjugation by diag(9,1).
    quadratic = QuadraticForm(QQ, primitive)
    clifford = CliffordAlgebra(quadratic)
    ee0, ee1, ee2 = clifford.gens()
    even_basis = [clifford.one(), ee0 * ee1, ee0 * ee2, ee1 * ee2]
    split_basis = [
        matrix(QQ, [[1, 0], [0, 1]]),
        matrix(QQ, [[-2, 0], [0, 1]]),
        matrix(QQ, [[4, 2], [-5, -3]]),
        matrix(QQ, [[2, 1], [5, 3]]),
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
        [[2, -1, 1, 5], [-1, 5, -11, -1], [1, -11, 5, 4], [5, -1, 4, 23]],
    )
    assert primitive_trace_pairing.det() == -(45**2)

    # For [[A,B],[C,D]], membership in the primitive order is equivalent to
    #   5|C, B=C/5 (mod 3), A-D=2B (mod 3).
    # Its determinant-one subgroup is K Gamma_0(15) K^-1.
    conjugator = matrix(ZZ, [[1, 0], [5, 1]])
    group15 = sl2(15)

    def in_primitive_order_mod15(value):
        aa, bb, cc, dd = value
        if cc % 5:
            return False
        return (bb - cc // 5) % 3 == 0 and (aa - dd - 2 * bb) % 3 == 0

    primitive_norm_one_mod15 = sorted(
        value for value in group15 if in_primitive_order_mod15(value)
    )
    gamma0_15_mod15 = [value for value in group15 if value[2] % 15 == 0]
    conjugated_gamma0_15 = sorted(
        tuple_matrix(
            conjugator
            * matrix(ZZ, 2, 2, value)
            * conjugator.inverse(),
            15,
        )
        for value in gamma0_15_mod15
    )
    assert primitive_norm_one_mod15 == conjugated_gamma0_15

    coarse_group = Gamma0(15)
    coarse_widths = sorted(
        int(coarse_group.cusp_width(cusp)) for cusp in coarse_group.cusps()
    )
    assert coarse_group.index() == 24
    assert coarse_group.genus() == 1
    assert coarse_group.nu2() == coarse_group.nu3() == 0
    assert coarse_widths == [1, 3, 5, 15]

    # Recover the spin action on T.  The corrected central volume element
    # identifies T_Q with the trace-zero part of the split even Clifford algebra.
    volume = (
        ee0 * ee1 * ee2
        - QQ(primitive[0, 1]) / 2 * ee2
        + QQ(primitive[0, 2]) / 2 * ee1
        - QQ(primitive[1, 2]) / 2 * ee0
    )
    assert all(volume * generator == generator * volume for generator in (ee0, ee1, ee2))
    assert volume * volume == QQ(45) / 4
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
        action = matrix(
            QQ, 3, 3, lambda row, column: columns[column][row]
        )
        assert all(entry in ZZ for entry in action.list())
        action = action.change_ring(ZZ)
        assert action.transpose() * transcendental * action == transcendental
        assert action.det() == 1
        return action

    moduli = tuple(discriminant_invariants)

    def discriminant_action_key(action):
        # A_T = T^*/T is represented as coker(T); the induced action there is
        # action^(-transpose).  smith_left supplies the coker Smith coordinates.
        contragredient = action.inverse().transpose().change_ring(ZZ)
        smith_action = smith_left * contragredient * smith_left.inverse()
        return tuple(
            int(smith_action[row, column] % moduli[row])
            for row in range(3)
            for column in range(3)
        )

    identity_action = tuple(
        1 if row == column else 0
        for row in range(3)
        for column in range(3)
    )

    def compose_discriminant(left, right):
        answer = matrix(ZZ, 3, 3, left) * matrix(ZZ, 3, 3, right)
        return tuple(
            int(answer[row, column] % moduli[row])
            for row in range(3)
            for column in range(3)
        )

    # Sage supplies a certified generating set of Gamma_0(15).  Conjugate it
    # into the displayed order and compute both its mod-2 and A_T actions.
    generator_records = []
    joint_generators = []
    for generator in coarse_group.gens():
        gamma0_matrix = matrix(ZZ, generator.matrix())
        order_matrix = conjugator * gamma0_matrix * conjugator.inverse()
        action = spin_action(order_matrix)
        discriminant_action = discriminant_action_key(action)
        mod2_action = tuple_matrix(order_matrix, 2)
        joint_generators.append((mod2_action, discriminant_action))
        generator_records.append(
            {
                "gamma0_15_generator": rows(gamma0_matrix),
                "primitive_order_generator": rows(order_matrix),
                "spin_action_on_T": rows(action),
                "mod_2_reduction": list(mod2_action),
                "discriminant_action_smith_coordinates": list(discriminant_action),
            }
        )

    identity_mod2 = (1, 0, 0, 1)

    def compose_joint(left, right):
        return (
            multiply_2x2(left[0], right[0], 2),
            compose_discriminant(left[1], right[1]),
        )

    joint_identity = (identity_mod2, identity_action)
    joint_image = {joint_identity}
    queue = deque([joint_identity])
    while queue:
        current = queue.popleft()
        for generator in joint_generators:
            following = compose_joint(current, generator)
            if following not in joint_image:
                joint_image.add(following)
                queue.append(following)

    mod2_image = {value[0] for value in joint_image}
    discriminant_image = {value[1] for value in joint_image}
    assert len(joint_image) == len(mod2_image) == len(discriminant_image) == 6
    assert sum(value[0] == identity_mod2 for value in joint_image) == 1
    assert sum(value[1] == identity_action for value in joint_image) == 1

    def action_order(value):
        current = identity_action
        for order in range(1, 7):
            current = compose_discriminant(current, value)
            if current == identity_action:
                return order
        raise ArithmeticError("discriminant image element order exceeds six")

    order_histogram = Counter(action_order(value) for value in discriminant_image)
    assert order_histogram == Counter({2: 3, 3: 2, 1: 1})

    # The literal even Clifford order has basis 1,2e0e1,2e0e2,2e1e2 in this
    # split representation.  Verify the full multiplication table against the
    # Clifford algebra of the literal Gram, not only its trace pairing.
    literal_split_basis = [split_basis[0]] + [2 * value for value in split_basis[1:]]
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
        [[2, -2, 2, 10], [-2, 20, -44, -4], [2, -44, 20, 16], [10, -4, 16, 92]],
    )
    assert literal_trace_pairing.det() == -(360**2)

    # Literal-order membership is equivalent to the primitive-order conditions
    # plus reduction to the identity modulo two for determinant-one matrices.
    def in_literal_order_mod30(value):
        aa, bb, cc, dd = value
        if bb % 2 or cc % 10:
            return False
        return (
            (bb // 2 - cc // 10) % 3 == 0
            and (aa - dd - 2 * bb + 3 * (cc // 5)) % 6 == 0
        )

    group30 = sl2(30)
    literal_norm_one_mod30 = {
        value for value in group30 if in_literal_order_mod30(value)
    }
    stable_from_primitive = {
        value
        for value in group30
        if in_primitive_order_mod15(tuple(entry % 15 for entry in value))
        and tuple(entry % 2 for entry in value) == identity_mod2
    }
    assert literal_norm_one_mod30 == stable_from_primitive

    # A negative-vector reflection represents the determinant-minus-one
    # component of O^+(T).  Its discriminant action is outside the S_3 image
    # of SO^+(T), so that component contributes no stable element.  Thus the
    # stable O^+ and stable SO^+ have the same projective action here.
    negative_vector = vector(ZZ, [0, 1, 0])
    assert negative_vector * transcendental * negative_vector == -4
    reflection = identity_matrix(QQ, 3) - (
        QQ(2) / (negative_vector * transcendental * negative_vector)
    ) * matrix(QQ, negative_vector).transpose() * matrix(QQ, negative_vector) * transcendental
    assert all(entry in ZZ for entry in reflection.list())
    reflection = reflection.change_ring(ZZ)
    assert reflection.det() == -1
    assert reflection.transpose() * transcendental * reflection == transcendental
    reflection_action = discriminant_action_key(reflection)
    assert reflection_action not in discriminant_image

    # Delta=Gamma_0(15) intersection Gamma(2) has b even and c divisible by
    # 30.  Conjugation by diag(2,1) identifies it over QQ with Gamma_0(60).
    stable_index_over_coarse = len(discriminant_image)
    stable_index = int(coarse_group.index()) * stable_index_over_coarse
    stable_group = Gamma0(60)
    assert stable_index == 144 == stable_group.index()
    assert stable_group.genus() == 7
    assert stable_group.nu2() == stable_group.nu3() == 0
    stable_cusps = []
    for cusp in stable_group.cusps():
        denominator = int(cusp.denominator()) if not cusp.is_infinity() else 0
        rationality_order = (
            1 if denominator == 0 else gcd(denominator, 60 // gcd(60, denominator))
        )
        assert rationality_order <= 2
        stable_cusps.append(
            {
                "representative": str(cusp),
                "width": int(stable_group.cusp_width(cusp)),
                "rational_over_QQ": True,
            }
        )
    assert len(stable_cusps) == 12
    assert sum(row["width"] for row in stable_cusps) == stable_index

    # The useful sixfold quotient is X_0(15).  The complete rational-point
    # determination has four cusps and four noncuspidal points.  These four
    # j-values are checked directly against their 3- and 5-isogenies and are
    # non-CM because a rational CM j-invariant is an integer.
    quotient_j_values = [
        -QQ(5**2) / 2,
        -QQ(5**2 * 241**3) / 2**3,
        -QQ(5 * 29**3) / 2**5,
        QQ(5 * 211**3) / 2**15,
    ]
    quotient_points = []
    for j_value in quotient_j_values:
        curve = EllipticCurve_from_j(j_value).global_minimal_model()
        isogeny_counts = {
            str(prime): len(curve.isogenies_prime_degree(ZZ(prime)))
            for prime in (3, 5)
        }
        assert isogeny_counts == {"3": 1, "5": 1}
        assert j_value.denominator() > 1
        quotient_points.append(
            {
                "j": str(j_value),
                "minimal_ainvariants": list(map(int, curve.ainvs())),
                "conductor": int(curve.conductor()),
                "rational_prime_isogeny_counts": isogeny_counts,
                "non_CM": True,
            }
        )

    assert 60 not in RATIONAL_CYCLIC_ISOGENY_DEGREES
    assert saturation["status"] == (
        "PASS_EXACT_RATIONAL_POINT_REJECTED_NS_DET20_TORSION3_HALF_SECTION"
    )
    assert saturation["discriminant_form"]["explicit_overlattice_index"] == 6
    assert saturation["full_surface"]["neron_severi_determinant"] == -20

    return {
        "schema": "elkies-k3.golay-det720-qq-marking-obstruction.v1",
        "status": "PASS_GOLAY_DET720_QQ_RATIONAL_MARKING_OBSTRUCTION",
        "surface_id": SURFACE_ID,
        "transcendental_lattice": {
            "gram": rows(transcendental),
            "determinant": int(transcendental.det()),
            "signature": [2, 1],
            "discriminant_group_invariants": discriminant_invariants,
            "primitive_similarity_gram": rows(primitive),
            "rational_U_coordinates": string_rows(u_change),
            "rational_U_gram": string_rows(u_gram),
        },
        "split_even_clifford_orders": {
            "primitive_similarity_order": {
                "basis_in_M2_QQ": [string_rows(value) for value in split_basis],
                "reduced_trace_pairing": rows(primitive_trace_pairing),
                "reduced_discriminant": 45,
                "matrix_conditions": [
                    "5 divides C",
                    "B = C/5 modulo 3",
                    "A-D = 2B modulo 3",
                ],
                "norm_one_group": "K Gamma_0(15) K^-1, K=[[1,0],[5,1]]",
            },
            "literal_T_order": {
                "basis_in_M2_QQ": [string_rows(value) for value in literal_split_basis],
                "reduced_trace_pairing": rows(literal_trace_pairing),
                "reduced_discriminant": 360,
                "matrix_conditions": [
                    "2 divides B",
                    "10 divides C",
                    "B/2 = C/10 modulo 3",
                    "A-D = 2B-3C/5 modulo 6",
                ],
                "norm_one_group": (
                    "K (Gamma_0(15) intersection Gamma(2)) K^-1"
                ),
            },
        },
        "coarse_norm_one_curve": {
            "label": "X_0(15)",
            "congruence_group": "K Gamma_0(15) K^-1",
            "congruence_level": 15,
            "index_in_PSL2Z": 24,
            "elliptic_orbits_order_2": 0,
            "elliptic_orbits_order_3": 0,
            "cusp_widths": coarse_widths,
            "genus": 1,
        },
        "stable_discriminant_kernel": {
            "abstract_group": "O^+(T)^* = kernel(O^+(T) -> O(A_T))",
            "discriminant_smith_invariants": discriminant_invariants,
            "proper_discriminant_image_order": len(discriminant_image),
            "proper_discriminant_image_isomorphism_type": "S3",
            "proper_discriminant_image_element_order_histogram": {
                str(order): count for order, count in sorted(order_histogram.items())
            },
            "generator_actions": generator_records,
            "kernel_condition": "identity modulo 2 inside the primitive norm-one group",
            "determinant_minus_one_stable_elements": 0,
            "negative_reflection_witness": {
                "vector": list(map(int, negative_vector)),
                "matrix": rows(reflection),
                "discriminant_action_smith_coordinates": list(reflection_action),
                "action_lies_in_proper_S3_image": False,
            },
            "congruence_group": "K (Gamma_0(15) intersection Gamma(2)) K^-1",
            "index_over_coarse_curve": stable_index_over_coarse,
            "index_in_PSL2Z": stable_index,
        },
        "exact_marked_modular_curve": {
            "label": "X_0(60)",
            "identification": (
                "diag(2,1)^(-1) (Gamma_0(15) intersection Gamma(2)) "
                "diag(2,1) = Gamma_0(60)"
            ),
            "defined_over": "QQ",
            "index_in_PSL2Z": stable_index,
            "genus": int(stable_group.genus()),
            "elliptic_orbits_order_2": int(stable_group.nu2()),
            "elliptic_orbits_order_3": int(stable_group.nu3()),
            "cusps": stable_cusps,
            "rational_points": {
                "cuspidal_points": len(stable_cusps),
                "noncuspidal_points": 0,
                "conclusion": "X_0(60)(QQ) consists of its twelve rational cusps",
                "theorem_input": "complete Mazur--Kenku rational cyclic-isogeny classification",
            },
        },
        "useful_quotient": {
            "map": "X_0(60)-model marked curve -> X_0(15)",
            "degree": stable_index_over_coarse,
            "target_genus": 1,
            "target_rational_cusps": 4,
            "target_rational_noncuspidal_points": quotient_points,
            "target_rational_noncuspidal_count": len(quotient_points),
            "lift_conclusion": (
                "None lifts to a QQ-point of the stable marked curve; such a lift "
                "would give a rational cyclic 60-isogeny."
            ),
        },
        "known_rational_3A5_control": {
            "certificate": relative(SATURATION),
            "displayed_sublattice_determinant": 720,
            "saturation_index": 6,
            "saturated_NS_determinant": 20,
            "torsion": "Z/3",
            "rational_half_section": True,
            "interpretation": (
                "This is a nonprimitive determinant-720 sublattice inside the "
                "determinant-20 Neron--Severi lattice, not a point of the primitive "
                "determinant-720 marked moduli problem."
            ),
        },
        "arithmetic_conclusion": (
            "No characteristic-zero K3 over QQ with geometric Neron--Severi "
            "lattice equal to the Golay determinant-720 lattice can have all "
            "nineteen divisor classes rational. Consequently its rootless frame "
            "cannot have a saturated rational MW17 basis over QQ(t)."
        ),
        "external_theorem_inputs": [
            {
                "name": "Mazur--Kenku classification of rational cyclic isogenies",
                "used_for": "X_0(60)(QQ) has no noncuspidal point",
                "allowed_degrees": list(RATIONAL_CYCLIC_ISOGENY_DEGREES),
            },
            {
                "name": "standard ternary spin and marked-K3 period correspondence",
                "used_for": (
                    "identification of the stable orthogonal period curve and the "
                    "rational point forced by a full rational Neron--Severi marking"
                ),
            },
            {
                "name": "complete rational-point determination of X_0(15)",
                "used_for": "the four-point noncuspidal list on the useful quotient",
            },
        ],
        "input_hashes": {
            relative(CATALOGUE): digest(CATALOGUE),
            relative(T_ARITHMETIC): digest(T_ARITHMETIC),
            relative(SATURATION): digest(SATURATION),
        },
        "proof_boundary": {
            "proved": (
                "The split embeddings, both integral order conditions, exact "
                "congruence groups, S3 discriminant action, stable kernel, X_0(60) "
                "identification, signatures, cusps, quotient j-values, and the "
                "determinant-20 saturation control are replayed exactly. With the "
                "named cyclic-isogeny and period theorems this excludes the full "
                "rational determinant-720 marking."
            ),
            "not_proved": (
                "The checker does not reprove the global Mazur--Kenku classification, "
                "the complete X_0(15)(QQ) determination, or the general marked-K3 "
                "period/spin correspondence. It makes no exclusion over larger number "
                "fields and no claim against geometric determinant-720 K3 surfaces."
            ),
        },
        "reproduce": (
            "sage -python "
            "elkies-k3/scripts/certify_golay_det720_qq_marking_obstruction.sage"
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
        "GOLAY720MARKING|coarse=X0(15)|Aimage=S3|stable=X0(60)|"
        "genus=7|QQ_noncusps=0|NS720=EXCLUDED|status=PASS"
    )


if __name__ == "__main__":
    main()
