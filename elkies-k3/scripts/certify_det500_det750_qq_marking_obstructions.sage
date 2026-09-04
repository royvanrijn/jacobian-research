#!/usr/bin/env sage-python
"""Certify the determinant-500 and determinant-750 marking obstructions.

The two rootless-MW17 catalogue rows have literal transcendental lattices

    T_N = U(5) + <10*N> = 5 (U + <2*N>),     N in {2,3}.

The primitive-similarity even Clifford order has norm-one curve X_0(N), but
the literal order and the stable discriminant kernel retain projective
identity modulo 5.  In Gamma_0(N) coordinates the stable subgroup is

    Gamma_0(N) intersection +/-Gamma(5).

After conjugation by diag(5,1), this is Gamma_H(25*N), where H is the inverse
image of {+/-1} in (Z/5Z)^*.  The resulting marked curves have genera 4 and
9 and degree-two forgetful maps to X_0(50) and X_0(75), respectively.  The
Mazur--Kenku rational cyclic-isogeny classification excludes noncuspidal
rational points on both quotients and hence on both marked curves.

The rational-isogeny classification, the standard normalizer description of
split prime-level Eichler orders, and the ternary-spin/marked-K3 period
correspondence are external theorem inputs.  Everything else is replayed
exactly below.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, deque
from pathlib import Path

from sage.all import (
    CliffordAlgebra,
    Cusp,
    Gamma0,
    GammaH,
    Integers,
    QQ,
    ZZ,
    QuadraticForm,
    gcd,
    identity_matrix,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CATALOGUE = GENERATED / "elkies-k3-rank7-auxiliary-catalogue-v1.json"
T_ARITHMETIC = GENERATED / "elkies-k3-rank7-t-arithmetic-v1.json"
OUTPUT = GENERATED / "elkies-k3-det500-det750-qq-marking-obstructions-v1.json"

CASES = (
    {
        "surface_id": "K3-04b86146cc6b284b",
        "N": 2,
        "determinant": 500,
        "frame_id": "K3-04b86146cc6b284b-F001",
        "frame_sha256": "60c96c1674d08624622723f7d7775cc74811e7e957723291ef89c2e7efaa3ba1",
        "stable_genus": 4,
        "stable_rational_cusps": ("1/25", "3/50", "2/25", "Infinity"),
    },
    {
        "surface_id": "K3-10a14a46c14b3150",
        "N": 3,
        "determinant": 750,
        "frame_id": "K3-10a14a46c14b3150-F001",
        "frame_sha256": "faf65e8ac34588dd5e4537069bae9b1294c106fee4887166e9029ca1231ee201",
        "stable_genus": 9,
        "stable_rational_cusps": ("2/75", "1/25", "2/25", "Infinity"),
    },
)

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


def rational_cusp_records(group, level):
    units = [ZZ(value) for value in Integers(level) if value.is_unit()]
    unit_generators = [ZZ(value) for value in Integers(level).unit_gens()]
    records = []
    rational = []
    for cusp in group.cusps():
        orbit = {
            str(group.reduce_cusp(Cusp(cusp).galois_action(unit, level)))
            for unit in units
        }
        is_rational = all(
            group.are_equivalent(cusp, Cusp(cusp).galois_action(unit, level))
            for unit in unit_generators
        )
        if is_rational:
            assert len(orbit) == 1
            rational.append(str(cusp))
        records.append(
            {
                "representative": str(cusp),
                "width": int(group.cusp_width(cusp)),
                "galois_orbit_representatives": sorted(orbit),
                "galois_orbit_size": len(orbit),
                "rational_over_QQ": is_rational,
            }
        )
    return records, rational


def build_case(config, catalogue, t_payload):
    surface_id = config["surface_id"]
    level = int(config["N"])
    scale = 5
    marked_level = level * scale**2
    expected_determinant = int(config["determinant"])
    surface = next(
        row for row in catalogue["surfaces"] if row["surface_id"] == surface_id
    )
    t_row = next(
        row for row in t_payload["surfaces"] if row["surface_id"] == surface_id
    )

    rootless_frames = [
        frame
        for frame in surface["frames"]
        if frame["frame_id"] == config["frame_id"]
        and frame["gram_sha256"] == config["frame_sha256"]
        and int(frame["mw_rank_for_rho_19"]) == 17
        and frame["root_type"] == "0"
        and int(frame["root_rank"]) == 0
        and int(frame["signed_root_count"]) == 0
        and int(frame["rootless_intrinsics"]["minimum_squared_norm"]) == 4
    ]
    assert len(rootless_frames) == 1

    transcendental = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    expected = matrix(
        ZZ,
        [[0, 0, scale], [0, 2 * scale * level, 0], [scale, 0, 0]],
    )
    assert transcendental == expected
    assert transcendental.det() == -expected_determinant
    assert QuadraticForm(ZZ, transcendental).signature_vector() == (2, 1, 0)
    assert t_row["literal_transcendental_gram"] == rows(transcendental)
    assert int(t_row["similarity_normalization"]["literal_content"]) == scale
    assert int(t_row["clifford"]["integral_even_clifford_order"]["reduced_discriminant"]) == level

    smith, smith_left, smith_right = transcendental.smith_form()
    assert smith_left * transcendental * smith_right == smith
    discriminant_invariants = [
        abs(int(smith[index, index])) for index in range(3)
    ]
    assert discriminant_invariants == [scale, scale, 2 * scale * level]
    assert surface["surface_key"]["ns_discriminant_form_key"]["invariants"] == discriminant_invariants

    # In the displayed basis, these columns are e0/5, e2/5, and
    # e1/(10*N).  They give the complete discriminant bilinear and quadratic
    # forms, not only the Smith invariants.
    dual_generators = matrix(
        QQ,
        [
            [QQ(1) / scale, 0, 0],
            [0, 0, QQ(1) / (2 * scale * level)],
            [0, QQ(1) / scale, 0],
        ],
    )
    assert all(entry in ZZ for entry in (transcendental * dual_generators).list())
    assert abs(dual_generators.det()) == QQ(1) / expected_determinant
    discriminant_bilinear_gram = (
        dual_generators.transpose() * transcendental * dual_generators
    )
    assert discriminant_bilinear_gram == matrix(
        QQ,
        [
            [0, QQ(1) / scale, 0],
            [QQ(1) / scale, 0, 0],
            [0, 0, QQ(1) / (2 * scale * level)],
        ],
    )

    primitive = (transcendental / scale).change_ring(ZZ)
    assert primitive == matrix(
        ZZ, [[0, 0, 1], [0, 2 * level, 0], [1, 0, 0]]
    )

    # Exact split representation of C^+(primitive), in the basis
    # 1,e0e1,e0e2,e1e2.  It is the Eichler order Gamma^0(N).
    quadratic = QuadraticForm(QQ, primitive)
    clifford = CliffordAlgebra(quadratic)
    ee0, ee1, ee2 = clifford.gens()
    even_basis = [clifford.one(), ee0 * ee1, ee0 * ee2, ee1 * ee2]
    split_basis = [
        matrix(QQ, [[1, 0], [0, 1]]),
        matrix(QQ, [[0, level], [0, 0]]),
        matrix(QQ, [[1, 0], [0, 0]]),
        matrix(QQ, [[0, 0], [1, 0]]),
    ]
    clifford_keys = list(clifford.basis().keys())

    def even_coordinates(value):
        coefficients = value.monomial_coefficients()
        return [
            coefficients.get(clifford_keys[index], 0)
            for index in (0, 4, 5, 6)
        ]

    def split_image(value):
        coordinates = even_coordinates(value)
        return sum(
            (
                coordinates[index] * split_basis[index]
                for index in range(4)
            ),
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
        [[2, 0, 1, 0], [0, 0, 0, level], [1, 0, 1, 0], [0, level, 0, 0]],
    )
    assert primitive_trace_pairing.det() == -(level**2)
    assert t_row["clifford"]["integral_even_clifford_order"][
        "reduced_trace_pairing"
    ] == rows(primitive_trace_pairing)

    coarse_group = Gamma0(level)
    coarse_cusps = [
        {
            "representative": str(cusp),
            "width": int(coarse_group.cusp_width(cusp)),
        }
        for cusp in coarse_group.cusps()
    ]
    assert coarse_group.genus() == 0

    # J conjugates Gamma_0(N) to the determinant-one group Gamma^0(N) of the
    # displayed primitive order.
    conjugator = matrix(ZZ, [[0, -1], [1, 0]])
    for generator in coarse_group.gens():
        gamma0_matrix = matrix(ZZ, generator.matrix())
        order_matrix = conjugator.inverse() * gamma0_matrix * conjugator
        assert order_matrix.det() == 1
        assert order_matrix[0, 1] % level == 0

    # Identify T_Q with the trace-zero Clifford subspace.  Conjugation by
    # primitive-order units gives the spin action on the literal lattice.
    volume = (
        ee0 * ee1 * ee2
        - QQ(primitive[0, 1]) / 2 * ee2
        + QQ(primitive[0, 2]) / 2 * ee1
        - QQ(primitive[1, 2]) / 2 * ee0
    )
    assert all(
        volume * generator == generator * volume
        for generator in (ee0, ee1, ee2)
    )
    assert volume * volume == QQ(level) / 4
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
        for image_value in vector_images:
            conjugate = value * image_value * inverse
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
        # A_T is coker(T); actions therefore enter contragrediently, followed
        # by the Smith-coordinate change.
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
        mod5_reduction = tuple_matrix(gamma0_matrix, scale)
        joint_generators.append((mod5_reduction, discriminant_action))
        spin_discriminant_generators.append(discriminant_action)
        generator_records.append(
            {
                "gamma0_generator": rows(gamma0_matrix),
                "primitive_order_generator": rows(order_matrix),
                "spin_action_on_T": rows(action),
                "mod_5_reduction": list(mod5_reduction),
                "discriminant_action_smith_coordinates": list(discriminant_action),
            }
        )

    identity_mod5 = (1, 0, 0, 1)

    def compose_joint(left, right):
        return (
            multiply_2x2(left[0], right[0], scale),
            compose_discriminant(left[1], right[1]),
        )

    joint_identity = (identity_mod5, identity_action)
    joint_image = closure(joint_identity, joint_generators, compose_joint)
    mod5_image = {value[0] for value in joint_image}
    spin_discriminant_image = {value[1] for value in joint_image}
    assert len(joint_image) == len(mod5_image) == 120
    assert len(spin_discriminant_image) == 60
    projective_identity_mod5 = {
        identity_mod5,
        negate_2x2(identity_mod5, scale),
    }
    assert {
        value[0] for value in joint_image if value[1] == identity_action
    } == projective_identity_mod5
    for mod5_value, discriminant_value in joint_image:
        partners = {
            other_mod5
            for other_mod5, other_discriminant in joint_image
            if other_discriminant == discriminant_value
        }
        assert partners == {
            mod5_value,
            negate_2x2(mod5_value, scale),
        }

    spin_order_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 60)
        for value in spin_discriminant_image
    )
    assert spin_order_histogram == Counter({5: 24, 3: 20, 2: 15, 1: 1})

    # Embed the literal C^+(T_N).  Under rational similarity its basis is
    # 1,5e0e1,5e0e2,5e1e2.  Its norm-one matrices are precisely primitive
    # order units that are scalar, hence +/-I, modulo five.
    literal_split_basis = [split_basis[0]] + [
        scale * value for value in split_basis[1:]
    ]
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
        coordinates = [
            coefficients.get(literal_keys[index], 0)
            for index in (0, 4, 5, 6)
        ]
        return sum(
            (
                coordinates[index] * literal_split_basis[index]
                for index in range(4)
            ),
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
    expected_literal_pairing = matrix(
        ZZ,
        [
            [2, 0, scale, 0],
            [0, 0, 0, scale**2 * level],
            [scale, 0, scale**2, 0],
            [0, scale**2 * level, 0, 0],
        ],
    )
    assert literal_trace_pairing == expected_literal_pairing
    literal_reduced_discriminant = level * scale**3
    assert literal_trace_pairing.det() == -(literal_reduced_discriminant**2)

    stable_joint = {
        value for value in joint_image if value[0] in projective_identity_mod5
    }
    assert len(stable_joint) == 2
    assert {value[1] for value in stable_joint} == {identity_action}

    # The split prime-level normalizer has one additional Fricke coset.
    # Together with a negative reflection and central inversion these give all
    # components of O^+(T).  None contains a stable discriminant action.
    negative_vector = vector(ZZ, [1, 0, -1])
    negative_norm = negative_vector * transcendental * negative_vector
    assert negative_norm == -2 * scale
    reflection = identity_matrix(QQ, 3) - (
        QQ(2) / negative_norm
    ) * matrix(QQ, negative_vector).transpose() * matrix(
        QQ, negative_vector
    ) * transcendental
    assert all(entry in ZZ for entry in reflection.list())
    reflection = reflection.change_ring(ZZ)
    assert reflection == matrix(ZZ, [[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    assert reflection.det() == -1
    reflection_action = discriminant_action_key(reflection)
    assert reflection_action not in spin_discriminant_image

    central_inversion = -identity_matrix(ZZ, 3)
    central_inversion_action = discriminant_action_key(central_inversion)
    assert central_inversion_action not in spin_discriminant_image

    fricke_isometry = -reflection
    fricke_action = discriminant_action_key(fricke_isometry)
    assert fricke_action == compose_discriminant(
        central_inversion_action, reflection_action
    )
    orientation_preserving_image = closure(
        identity_action,
        spin_discriminant_generators + [fricke_action],
        compose_discriminant,
    )
    assert len(orientation_preserving_image) == 120
    orientation_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 120)
        for value in orientation_preserving_image
    )
    assert orientation_histogram == Counter(
        {4: 30, 2: 25, 5: 24, 6: 20, 3: 20, 1: 1}
    ), orientation_histogram

    full_discriminant_image = closure(
        identity_action,
        spin_discriminant_generators
        + [reflection_action, central_inversion_action],
        compose_discriminant,
    )
    assert len(full_discriminant_image) == 240
    full_histogram = Counter(
        action_order(value, identity_action, compose_discriminant, 120)
        for value in full_discriminant_image
    )
    assert full_histogram == Counter(
        {6: 60, 4: 60, 2: 51, 10: 24, 5: 24, 3: 20, 1: 1}
    ), full_histogram
    non_spin_action_cosets = [
        {
            compose_discriminant(representative, spin_value)
            for spin_value in spin_discriminant_image
        }
        for representative in (
            fricke_action,
            reflection_action,
            central_inversion_action,
        )
    ]
    assert all(len(coset) == len(spin_discriminant_image) for coset in non_spin_action_cosets)
    assert all(identity_action not in coset for coset in non_spin_action_cosets)

    # Conjugating by diag(5,1) sends the stable group to Gamma_H(25*N),
    # where H is the inverse image of +/-1 modulo five.  Unlike the mod-2 and
    # mod-3 cases, H is a proper index-two subgroup of all diamond operators;
    # the marked curve is therefore a degree-two cover of X_0(25*N), not that
    # coarse quotient itself.
    diamond_subgroup = [
        residue
        for residue in range(marked_level)
        if gcd(residue, marked_level) == 1
        and residue % scale in (1, scale - 1)
    ]
    stable_group = GammaH(marked_level, diamond_subgroup)
    stable_index_over_coarse = len(spin_discriminant_image)
    stable_index = int(coarse_group.index()) * stable_index_over_coarse
    assert stable_group.index() == stable_index
    assert stable_group.genus() == int(config["stable_genus"])
    assert stable_group.nu2() == stable_group.nu3() == 0
    stable_cusps, stable_rational_cusps = rational_cusp_records(
        stable_group, marked_level
    )
    assert len(stable_cusps) == 24
    assert sum(row["width"] for row in stable_cusps) == stable_index
    assert tuple(stable_rational_cusps) == config["stable_rational_cusps"]

    quotient_group = Gamma0(marked_level)
    quotient_cusps, quotient_rational_cusps = rational_cusp_records(
        quotient_group, marked_level
    )
    assert stable_group.index() == 2 * quotient_group.index()
    assert len(quotient_cusps) == 12
    assert len(quotient_rational_cusps) == 4
    assert marked_level not in RATIONAL_CYCLIC_ISOGENY_DEGREES

    return {
        "surface_id": surface_id,
        "determinant": expected_determinant,
        "classification": "ARITHMETICALLY_EXCLUDED",
        "rootless_candidate": {
            "frame_id": config["frame_id"],
            "frame_gram_sha256": config["frame_sha256"],
            "rank": 17,
            "minimum_squared_norm": 4,
            "signed_root_count": 0,
        },
        "transcendental_lattice": {
            "gram": rows(transcendental),
            "determinant": int(transcendental.det()),
            "signature": [2, 1],
            "smith_invariants": discriminant_invariants,
            "discriminant_bilinear_gram_in_dual_basis": string_rows(
                discriminant_bilinear_gram
            ),
            "discriminant_quadratic_diagonal_mod_2Z": [
                str(discriminant_bilinear_gram[index, index])
                for index in range(3)
            ],
            "primitive_similarity_gram": rows(primitive),
            "similarity_factor": scale,
        },
        "clifford_orders": {
            "primitive_similarity_order": {
                "basis_in_M2_QQ": [string_rows(value) for value in split_basis],
                "reduced_trace_pairing": rows(primitive_trace_pairing),
                "reduced_discriminant": level,
                "integral_matrix_conditions": [f"{level} divides B"],
                "norm_one_group": f"Gamma^0({level})",
                "conjugate_curve": f"X_0({level})",
            },
            "literal_T_order": {
                "basis_in_M2_QQ": [
                    string_rows(value) for value in literal_split_basis
                ],
                "reduced_trace_pairing": rows(literal_trace_pairing),
                "reduced_discriminant": literal_reduced_discriminant,
                "integral_matrix_conditions": [
                    f"{scale * level} divides B",
                    f"{scale} divides C",
                    f"A = D modulo {scale}",
                ],
                "norm_one_condition_inside_primitive_order": "+/-I modulo 5",
            },
        },
        "coarse_norm_one_curve": {
            "label": f"X_0({level})",
            "index_in_PSL2Z": int(coarse_group.index()),
            "genus": int(coarse_group.genus()),
            "elliptic_orbits_order_2": int(coarse_group.nu2()),
            "elliptic_orbits_order_3": int(coarse_group.nu3()),
            "cusps": coarse_cusps,
            "warning": "This primitive-similarity curve forgets the mod-5 marking.",
        },
        "stable_discriminant_kernel": {
            "spin_generators": generator_records,
            "proper_spin_image_order": len(spin_discriminant_image),
            "proper_spin_image_isomorphism_type": "PSL_2(F_5) ~= A5",
            "proper_spin_image_element_orders": {
                str(key): value for key, value in sorted(spin_order_histogram.items())
            },
            "kernel_condition": "+/-I modulo 5 inside Gamma_0(N)",
            "orientation_preserving_image_order": len(
                orientation_preserving_image
            ),
            "orientation_preserving_image_isomorphism_type": "S5",
            "full_O_plus_image_order": len(full_discriminant_image),
            "full_O_plus_image_isomorphism_type": "S5 x C2",
            "full_O_plus_image_element_orders": {
                str(key): value for key, value in sorted(full_histogram.items())
            },
            "non_spin_stable_cosets": 0,
            "index_over_coarse_curve": stable_index_over_coarse,
            "index_in_PSL2Z": stable_index,
        },
        "exact_marked_modular_curve": {
            "label": (
                f"X_H({marked_level}), H={{a in (Z/{marked_level}Z)^*: "
                "a=+/-1 mod 5}}"
            ),
            "group_before_conjugation": (
                f"Gamma_0({level}) intersection +/-Gamma(5)"
            ),
            "group_after_diag_5_1_conjugation": f"Gamma_H({marked_level})",
            "diamond_subgroup": diamond_subgroup,
            "index_in_PSL2Z": stable_index,
            "genus": int(stable_group.genus()),
            "elliptic_orbits_order_2": int(stable_group.nu2()),
            "elliptic_orbits_order_3": int(stable_group.nu3()),
            "geometric_cusp_count": len(stable_cusps),
            "cusps": stable_cusps,
            "rational_cusps": stable_rational_cusps,
            "rational_points": {
                "cuspidal_points": len(stable_rational_cusps),
                "noncuspidal_points": 0,
                "CM_points": 0,
                "non_CM_points": 0,
                "full_curve_points_realizing_only_an_overlattice": 0,
                "coarse_quotient_boundary": (
                    f"X_0({level})(QQ) is infinite because it is genus zero "
                    "with a QQ-cusp; its noncuspidal points are quotient-level "
                    "period data, not full rational markings, and none has a "
                    "noncuspidal QQ-lift to this exact marked curve."
                ),
                "conclusion": (
                    "Every QQ-point is one of the four displayed cusps; there "
                    "is no K3 period point."
                ),
            },
        },
        "quotient_maps": [
            {
                "source": f"X_H({marked_level})",
                "target": f"X_0({marked_level})",
                "degree": 2,
                "target_index_in_PSL2Z": int(quotient_group.index()),
                "target_genus": int(quotient_group.genus()),
                "target_geometric_cusp_count": len(quotient_cusps),
                "target_rational_cusps": quotient_rational_cusps,
                "target_rational_noncuspidal_points": 0,
                "exclusion": (
                    f"{marked_level} is absent from the Mazur--Kenku list of "
                    "rational cyclic-isogeny degrees."
                ),
            },
            {
                "source": f"X_H({marked_level})",
                "target": f"X_0({level})",
                "degree": stable_index_over_coarse,
                "target_genus": int(coarse_group.genus()),
                "warning": "This coarse quotient alone is not a marking decision.",
            },
        ],
        "decision": (
            f"A rational noncuspidal marked point would map to a rational "
            f"cyclic {marked_level}-isogeny, which does not exist."
        ),
    }


def build_payload():
    catalogue = json.loads(CATALOGUE.read_text())
    t_payload = json.loads(T_ARITHMETIC.read_text())
    cases = [build_case(config, catalogue, t_payload) for config in CASES]
    assert [row["classification"] for row in cases] == [
        "ARITHMETICALLY_EXCLUDED",
        "ARITHMETICALLY_EXCLUDED",
    ]
    return {
        "schema": "elkies-k3.det500-det750-qq-marking-obstructions.v1",
        "status": "PASS_DET500_DET750_QQ_RATIONAL_MARKING_OBSTRUCTIONS",
        "inputs": {
            relative(CATALOGUE): digest(CATALOGUE),
            relative(T_ARITHMETIC): digest(T_ARITHMETIC),
        },
        "cases": cases,
        "cases_by_surface_id": {row["surface_id"]: row for row in cases},
        "accounting": {
            "selected_rootless_MW17_rows": 2,
            "ARITHMETICALLY_EXCLUDED": 2,
            "ARITHMETICALLY_REALIZABLE": 0,
            "UNRESOLVED_FOR_EXPLICIT_REASON": 0,
            "equation_work_authorized": 0,
        },
        "theorem_inputs": [
            "Mazur--Kenku's complete classification of rational cyclic-isogeny degrees",
            "the standard split prime-level Eichler normalizer description",
            "the rank-three ternary-spin/fully-marked-K3 period correspondence",
        ],
        "proof_boundary": {
            "proved": (
                "The literal Clifford orders, full discriminant actions, stable "
                "congruence curves, signatures, cusp Galois orbits, quotient maps, "
                "and rational-point exclusions for the two selected rows."
            ),
            "external": (
                "The three named theorem inputs are used but not reproved by this checker."
            ),
            "not_claimed": (
                "No statement is made about other determinant-500 or determinant-750 "
                "ternary genera, or about realizability over number fields larger than QQ."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/certify_det500_det750_qq_marking_obstructions.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "DET500DET750MARKING|stable=XH50,XH75|genera=4,9|"
        "QQ_noncusps=0,0|excluded=2|status=PASS"
    )


if __name__ == "__main__":
    main()
