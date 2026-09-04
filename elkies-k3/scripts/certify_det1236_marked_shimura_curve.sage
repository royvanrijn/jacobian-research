#!/usr/bin/env sage-python
"""Identify the exact marked period curve for the determinant-1236 row.

For

    T = [[-2,0,1],[0,4,0],[1,0,154]],

the literal even Clifford order is an Eichler order of discriminant 6 and
level 103.  Its norm-one curve is X_0^6(103).  Exact normalizer elements show
that the projective stable discriminant kernel adds precisely the full
Fricke involution w_618, so the fully marked period curve is

    C_1236 = X_0^6(103)/<w_618>.

The script also computes the Atkin--Lehner quotient tower and verifies the
published genus-two quotient model.  It deliberately stops before claiming
an arithmetic realization: the rational points on C_1236 beyond its two
certified rational discriminant -3 CM points are not known here.

The Eichler normalizer theorem, Ogg's fixed-point formula, the CM residue
field formula of Gonzalez--Rotger, the published genus-two model, and the
ternary-spin/marked-K3 period correspondence are named external inputs.
All lattice, Clifford, discriminant-action, genus, model, and point
substitutions used around those inputs are replayed exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    CliffordAlgebra,
    EllipticCurve,
    HyperellipticCurve,
    ModularSymbols,
    PolynomialRing,
    QQ,
    ZZ,
    QuadraticForm,
    QuaternionAlgebra,
    factor,
    gcd,
    identity_matrix,
    kronecker,
    matrix,
    prod,
    vector,
)
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CATALOGUE = GENERATED / "elkies-k3-rank7-auxiliary-catalogue-v1.json"
T_ARITHMETIC = GENERATED / "elkies-k3-rank7-t-arithmetic-v1.json"
OUTPUT = GENERATED / "elkies-k3-det1236-marked-shimura-curve-v1.json"

SURFACE_ID = "K3-6d288cfad55e0d15"
DETERMINANT = 1236
FRAME_ID = "K3-6d288cfad55e0d15-F001"
FRAME_SHA256 = "1fb6824ca519210cbfe32af3c0c4c15fb1d3867c5847cca32878e5a4929b0682"
QUATERNION_DISCRIMINANT = 6
EICHLER_LEVEL = 103
REDUCED_ORDER_DISCRIMINANT = 618
LOW_GENUS_SOURCE_COMMIT = "6cc368fe37aa67187783118f18d149b2b1fd6230"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def squarefree_part(value):
    return ZZ(prod(prime for prime, exponent in factor(value) if exponent % 2))


def order_class_number(discriminant):
    return len(
        BinaryQF_reduced_representatives(
            ZZ(discriminant), primitive_only=True
        )
    )


def phi_from_factorization(value):
    return ZZ(value) * prod(
        1 - QQ(1) / prime for prime, _ in factor(ZZ(value))
    )


def psi_from_factorization(value):
    return prod(
        (prime + 1) * prime ** (exponent - 1)
        for prime, exponent in factor(ZZ(value))
    )


def local_embedding_count(D, N, prime, field_discriminant, conductor):
    """Ogg's local optimal-embedding factor."""
    D = ZZ(D)
    N = ZZ(N)
    prime = ZZ(prime)
    conductor = ZZ(conductor)
    symbol_order = (
        1
        if conductor % prime == 0
        else kronecker(field_discriminant, prime)
    )
    if D % prime == 0:
        return ZZ(1 - symbol_order)

    level_valuation = N.valuation(prime)
    conductor_valuation = conductor.valuation(prime)
    if level_valuation == 1:
        return ZZ(1 + symbol_order)

    def psi_prime(value, pp):
        exponent = ZZ(value).valuation(pp)
        return ZZ(1 if exponent == 0 else pp**exponent + pp ** (exponent - 1))

    field_symbol = kronecker(field_discriminant, prime)
    if level_valuation >= 2 + 2 * conductor_valuation:
        return ZZ(2 * psi_prime(conductor, prime) if field_symbol == 1 else 0)
    if level_valuation == 1 + 2 * conductor_valuation:
        if field_symbol == 1:
            return ZZ(2 * psi_prime(conductor, prime))
        if field_symbol == 0:
            return ZZ(prime**conductor_valuation)
        return ZZ(0)
    if level_valuation == 2 * conductor_valuation:
        return ZZ(
            prime ** (conductor_valuation - 1)
            * (prime + 1 + field_symbol)
        )
    if conductor**2 % (prime * N) == 0:
        if level_valuation % 2 == 0:
            return ZZ(prime**level_valuation + prime ** (level_valuation - 1))
        return ZZ(2 * prime**level_valuation)
    raise ArithmeticError("unhandled Ogg local-embedding case")


def fixed_point_record(D, N, atkin_lehner_label):
    """Ogg fixed-point count for one Atkin--Lehner involution."""
    D = ZZ(D)
    N = ZZ(N)
    label = ZZ(atkin_lehner_label)
    sf_part = squarefree_part(label)
    conductor = ZZ((label // sf_part).sqrt())
    assert conductor**2 * sf_part == label
    if label == 2:
        orders = [(-4, 1), (-8, 1)]
    elif label % 4 == 3:
        orders = [(-sf_part, conductor), (-sf_part, 2 * conductor)]
    else:
        orders = [(-4 * sf_part, conductor)]

    contributions = []
    total = ZZ(0)
    for field_discriminant, order_conductor in orders:
        local_factors = []
        for prime, _ in factor(D * N // label):
            local_factors.append(
                [
                    int(prime),
                    int(
                        local_embedding_count(
                            D,
                            N,
                            prime,
                            field_discriminant,
                            order_conductor,
                        )
                    ),
                ]
            )
        class_number = order_class_number(
            field_discriminant * order_conductor**2
        )
        contribution = ZZ(class_number) * prod(
            factor_value for _, factor_value in local_factors
        )
        total += contribution
        contributions.append(
            {
                "field_discriminant": int(field_discriminant),
                "conductor": int(order_conductor),
                "order_discriminant": int(
                    field_discriminant * order_conductor**2
                ),
                "class_number": int(class_number),
                "local_embedding_factors": local_factors,
                "contribution": int(contribution),
            }
        )
    return {
        "atkin_lehner_label": int(label),
        "fixed_points": int(total),
        "cm_contributions": contributions,
    }


def build_payload(catalogue, t_arithmetic):
    surface = next(
        row for row in catalogue["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    t_row = next(
        row for row in t_arithmetic["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    assert int(surface["determinant"]) == DETERMINANT
    frame = next(row for row in surface["frames"] if row["frame_id"] == FRAME_ID)
    assert frame["gram_sha256"] == FRAME_SHA256
    assert int(frame["mw_rank_for_rho_19"]) == 17
    assert frame["root_type"] == "0"
    assert int(frame["root_rank"]) == 0
    assert int(frame["signed_root_count"]) == 0
    assert int(frame["rootless_intrinsics"]["minimum_squared_norm"]) == 4

    transcendental = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    expected_transcendental = matrix(
        ZZ, [[-2, 0, 1], [0, 4, 0], [1, 0, 154]]
    )
    assert transcendental == expected_transcendental
    assert transcendental.det() == -DETERMINANT
    assert QuadraticForm(ZZ, transcendental).signature_vector() == (2, 1, 0)
    assert rows(transcendental) == t_row["literal_transcendental_gram"]
    assert int(t_row["similarity_normalization"]["literal_content"]) == 1
    assert matrix(
        ZZ,
        t_row["similarity_normalization"]["primitive_integral_quadratic_hessian"],
    ) == transcendental

    smith, smith_left, smith_right = transcendental.smith_form()
    assert smith_left * transcendental * smith_right == smith
    assert smith.diagonal() == [1, 1, DETERMINANT]
    smith_generator_integral = smith_left.inverse() * vector(ZZ, [0, 0, 1])
    dual_generator = transcendental.inverse() * smith_generator_integral
    assert transcendental * dual_generator == smith_generator_integral
    discriminant_bilinear_value = dual_generator * transcendental * dual_generator
    discriminant_quadratic_value = discriminant_bilinear_value / 2
    assert discriminant_bilinear_value == QQ(317) / 1236
    assert discriminant_quadratic_value == QQ(317) / 2472

    discriminant_orthogonal_units = [
        unit
        for unit in range(DETERMINANT)
        if gcd(unit, DETERMINANT) == 1
        and (unit**2 - 1) * discriminant_quadratic_value in ZZ
    ]
    assert discriminant_orthogonal_units == [
        1,
        205,
        413,
        617,
        619,
        823,
        1031,
        1235,
    ]

    quadratic = QuadraticForm(QQ, transcendental)
    clifford = CliffordAlgebra(quadratic)
    ee0, ee1, ee2 = clifford.gens()
    even_basis = [
        clifford.one(),
        ee0 * ee1,
        ee0 * ee2,
        ee1 * ee2,
    ]
    clifford_keys = list(clifford.basis().keys())
    even_indices = (0, 4, 5, 6)

    def even_coordinates(value):
        coefficients = value.monomial_coefficients()
        coordinates = vector(
            QQ,
            [coefficients.get(clifford_keys[index], 0) for index in even_indices],
        )
        reconstruction = sum(
            (
                coordinates[index] * even_basis[index]
                for index in range(4)
            ),
            clifford.zero(),
        )
        assert reconstruction == value
        return coordinates

    def scalar_part(value):
        return value.monomial_coefficients().get(clifford_keys[0], 0)

    for left in even_basis:
        for right in even_basis:
            product_coordinates = even_coordinates(left * right)
            assert all(entry in ZZ for entry in product_coordinates)

    reduced_trace_pairing = matrix(
        ZZ,
        4,
        4,
        lambda left, right: scalar_part(
            even_basis[left] * even_basis[right]
            + (even_basis[left] * even_basis[right]).clifford_conjugate()
        ),
    )
    expected_trace_pairing = matrix(
        ZZ,
        [
            [2, 0, 1, 0],
            [0, 4, 0, 2],
            [1, 0, 155, 0],
            [0, 2, 0, -308],
        ],
    )
    assert reduced_trace_pairing == expected_trace_pairing
    assert reduced_trace_pairing.det() == -(REDUCED_ORDER_DISCRIMINANT**2)
    assert t_row["clifford"]["integral_even_clifford_order"][
        "reduced_trace_pairing"
    ] == rows(reduced_trace_pairing)
    assert int(
        t_row["clifford"]["integral_even_clifford_order"][
            "reduced_discriminant"
        ]
    ) == REDUCED_ORDER_DISCRIMINANT

    quaternion = QuaternionAlgebra(QQ, 2, QQ(309) / 4)
    assert quaternion.ramified_primes() == [2, 3]
    assert quaternion.discriminant() == QUATERNION_DISCRIMINANT
    assert REDUCED_ORDER_DISCRIMINANT == (
        QUATERNION_DISCRIMINANT * EICHLER_LEVEL
    )

    volume = (
        ee0 * ee1 * ee2
        - QQ(transcendental[0, 1]) / 2 * ee2
        + QQ(transcendental[0, 2]) / 2 * ee1
        - QQ(transcendental[1, 2]) / 2 * ee0
    )
    assert volume * volume == QQ(309) / 2
    vector_images = [volume * generator for generator in (ee0, ee1, ee2)]
    image_matrix = matrix(
        QQ,
        4,
        3,
        lambda row, column: even_coordinates(vector_images[column])[row],
    )
    left_inverse = (
        image_matrix.transpose() * image_matrix
    ).inverse() * image_matrix.transpose()
    assert left_inverse * image_matrix == identity_matrix(QQ, 3)

    def even_element(coordinates):
        return sum(
            (
                QQ(coordinates[index]) * even_basis[index]
                for index in range(4)
            ),
            clifford.zero(),
        )

    def reduced_norm(value):
        product_value = value * value.clifford_conjugate()
        coordinates = even_coordinates(product_value)
        assert coordinates[1:] == vector(QQ, [0, 0, 0])
        return coordinates[0]

    def normalizer_record(label, coordinates, expected_norm, expected_unit):
        element = even_element(coordinates)
        norm = reduced_norm(element)
        assert norm == expected_norm
        inverse = element.clifford_conjugate() / norm
        order_columns = []
        for basis_value in even_basis:
            column = even_coordinates(element * basis_value * inverse)
            assert all(entry in ZZ for entry in column)
            order_columns.append(column)
        order_action = matrix(
            ZZ,
            4,
            4,
            lambda row, column: order_columns[column][row],
        )
        assert abs(order_action.det()) == 1

        t_columns = []
        for vector_image in vector_images:
            conjugate = element * vector_image * inverse
            coordinates_on_t = left_inverse * even_coordinates(conjugate)
            assert image_matrix * coordinates_on_t == even_coordinates(conjugate)
            assert all(entry in ZZ for entry in coordinates_on_t)
            t_columns.append(coordinates_on_t)
        t_action = matrix(
            ZZ,
            3,
            3,
            lambda row, column: t_columns[column][row],
        )
        assert t_action.transpose() * transcendental * t_action == transcendental
        assert t_action.det() == 1

        contragredient = t_action.inverse().transpose().change_ring(ZZ)
        smith_action = smith_left * contragredient * smith_left.inverse()
        assert all(
            smith_action[2, column] % DETERMINANT == 0
            for column in (0, 1)
        )
        discriminant_unit = int(smith_action[2, 2] % DETERMINANT)
        assert discriminant_unit == expected_unit
        return {
            "atkin_lehner_label": int(label),
            "even_order_coordinates": list(map(int, coordinates)),
            "reduced_norm": int(norm),
            "conjugation_on_even_order": rows(order_action),
            "spin_action_on_T": rows(t_action),
            "action_on_cyclic_A_T": discriminant_unit,
        }

    normalizer_records = {
        "w_2": normalizer_record(2, (-14, -3, 2, 1), 2, 619),
        "w_6": normalizer_record(6, (-10, -15, 8, -6), 6, 1031),
        "w_618": normalizer_record(618, (-4, -3, 8, 6), 618, 1235),
    }

    atkin_actions = {
        1: 1,
        2: normalizer_records["w_2"]["action_on_cyclic_A_T"],
        6: normalizer_records["w_6"]["action_on_cyclic_A_T"],
        618: normalizer_records["w_618"]["action_on_cyclic_A_T"],
    }
    atkin_actions[3] = atkin_actions[2] * atkin_actions[6] % DETERMINANT
    atkin_actions[309] = atkin_actions[618] * atkin_actions[2] % DETERMINANT
    atkin_actions[103] = atkin_actions[618] * atkin_actions[6] % DETERMINANT
    atkin_actions[206] = atkin_actions[618] * atkin_actions[3] % DETERMINANT
    atkin_actions = dict(sorted(atkin_actions.items()))
    assert list(atkin_actions) == [1, 2, 3, 6, 103, 206, 309, 618]
    assert sorted(atkin_actions.values()) == discriminant_orthogonal_units
    assert atkin_actions == {
        1: 1,
        2: 619,
        3: 413,
        6: 1031,
        103: 205,
        206: 823,
        309: 617,
        618: 1235,
    }
    projectively_stable_atkin_labels = [
        label
        for label, unit in atkin_actions.items()
        if unit in (1, DETERMINANT - 1)
    ]
    assert projectively_stable_atkin_labels == [1, 618]

    D = ZZ(QUATERNION_DISCRIMINANT)
    N = ZZ(EICHLER_LEVEL)
    mu = phi_from_factorization(D) * psi_from_factorization(N)
    elliptic_order_2 = prod(
        1 - kronecker(-4, prime) for prime, _ in factor(D)
    ) * prod(1 + kronecker(-4, prime) for prime, _ in factor(N))
    elliptic_order_3 = prod(
        1 - kronecker(-3, prime) for prime, _ in factor(D)
    ) * prod(1 + kronecker(-3, prime) for prime, _ in factor(N))
    top_genus = ZZ(
        1
        + mu / 12
        - QQ(elliptic_order_2) / 4
        - QQ(elliptic_order_3) / 3
    )
    assert (mu, elliptic_order_2, elliptic_order_3, top_genus) == (208, 0, 4, 17)

    hall_divisors = [1, 2, 3, 6, 103, 206, 309, 618]
    fixed_records = {
        label: fixed_point_record(D, N, label) for label in hall_divisors[1:]
    }
    fixed_counts = {
        label: record["fixed_points"] for label, record in fixed_records.items()
    }
    assert fixed_counts == {2: 0, 3: 4, 6: 4, 103: 0, 206: 0, 309: 12, 618: 12}

    quotient_subgroups = {
        "stable_full_fricke": [1, 618],
        "genus_two_target": [1, 2, 309, 618],
        "alternate_genus_three_w3": [1, 3, 206, 618],
        "alternate_genus_three_w6": [1, 6, 103, 618],
        "full_atkin_lehner": hall_divisors,
    }
    quotient_records = {}
    for name, subgroup in quotient_subgroups.items():
        fixed_sum = sum(fixed_counts[label] for label in subgroup if label != 1)
        quotient_genus = ZZ(
            (2 * top_genus - 2 - fixed_sum) / (2 * len(subgroup)) + 1
        )
        quotient_records[name] = {
            "atkin_lehner_subgroup": subgroup,
            "degree_from_X_0_6_103": len(subgroup),
            "fixed_point_sum": fixed_sum,
            "genus": int(quotient_genus),
        }
    assert {
        name: row["genus"] for name, row in quotient_records.items()
    } == {
        "stable_full_fricke": 6,
        "genus_two_target": 2,
        "alternate_genus_three_w3": 3,
        "alternate_genus_three_w6": 3,
        "full_atkin_lehner": 1,
    }

    # The stable quotient has twelve w_618 branch points of orbifold order 2.
    # Its four order-3 points upstairs are paired by w_618.
    stable_signature = {
        "genus": 6,
        "cusps": 0,
        "elliptic_points_order_2": 12,
        "elliptic_points_order_3": 2,
        "orbifold_mu": 104,
    }
    assert QQ(1) + QQ(104) / 12 - QQ(12) / 4 - QQ(2) / 3 == 6

    polynomial_ring = PolynomialRing(QQ, "x")
    xx = polynomial_ring.gen()
    genus_two_polynomial = 1944 * xx**6 + 441 * xx**4 - 90 * xx**2 + 9
    assert genus_two_polynomial.is_squarefree()
    genus_two_curve = HyperellipticCurve(genus_two_polynomial)
    assert genus_two_curve.genus() == 2
    rational_points_on_genus_two_target = [
        [-1, 48],
        [-1, -48],
        ["-1/3", "8/3"],
        ["-1/3", "-8/3"],
        ["-1/5", "312/125"],
        ["-1/5", "-312/125"],
        [0, 3],
        [0, -3],
        ["1/5", "312/125"],
        ["1/5", "-312/125"],
        ["1/3", "8/3"],
        ["1/3", "-8/3"],
        [1, 48],
        [1, -48],
    ]
    for x_coordinate, y_coordinate in rational_points_on_genus_two_target:
        assert QQ(y_coordinate) ** 2 == genus_two_polynomial(QQ(x_coordinate))

    cubic_model = EllipticCurve(
        [0, 441, 0, 1944 * (-90), 1944**2 * 9]
    )
    full_quotient_elliptic_curve = EllipticCurve("618f1")
    assert cubic_model.global_minimal_model() == full_quotient_elliptic_curve
    assert list(full_quotient_elliptic_curve.a_invariants()) == [1, 0, 0, -185, 1401]
    assert full_quotient_elliptic_curve.conductor() == 618
    assert full_quotient_elliptic_curve.rank() == 1
    assert full_quotient_elliptic_curve.torsion_subgroup().order() == 1
    assert full_quotient_elliptic_curve.gens() == [full_quotient_elliptic_curve(10, -29)]
    generator = full_quotient_elliptic_curve.gens()[0]
    # Quotient map from the displayed genus-two model by x -> -x.
    def elliptic_image(x_coordinate, y_coordinate):
        return full_quotient_elliptic_curve(
            54 * QQ(x_coordinate) ** 2 + 4,
            9 * QQ(y_coordinate) - 27 * QQ(x_coordinate) ** 2 - 2,
        )

    assert elliptic_image(0, 3) == 3 * generator
    assert elliptic_image(1, 48) == 4 * generator
    assert elliptic_image(QQ(1) / 3, QQ(-8) / 3) == generator
    assert elliptic_image(QQ(1) / 3, QQ(8) / 3) == -generator
    assert elliptic_image(QQ(1) / 5, QQ(-312) / 125) == 10 * generator
    assert elliptic_image(QQ(1) / 5, QQ(312) / 125) == -10 * generator

    # Jacquet--Langlands identifies the differentials on X_0^6(103) with
    # the weight-two newspace of classical level 2*3*103.  At the ramified
    # primes 2 and 3 the geometric Atkin--Lehner signs are the negatives of
    # the classical signs; at the Eichler-level prime 103 they agree.  Thus
    # the geometric w_618-invariants are exactly the classical sign-product
    # +1 orbits.  We compute all orbits so the dimension-six assertion is an
    # exact accounting rather than an inference from the known quotient.
    modular_symbols_newspace = (
        ModularSymbols(618, 2, sign=1).cuspidal_subspace().new_subspace()
    )
    assert modular_symbols_newspace.dimension() == 17
    comparison_primes = [5, 7, 11]
    rational_newform_labels = [f"618{letter}1" for letter in "abcdefg"]
    signature_to_label = {
        tuple(EllipticCurve(label).ap(prime) for prime in comparison_primes): label
        for label in rational_newform_labels
    }
    assert len(signature_to_label) == len(rational_newform_labels)

    newform_orbits = []
    for factor_space in modular_symbols_newspace.decomposition():
        classical_signs = []
        for prime in [2, 3, 103]:
            operator = factor_space.atkin_lehner_operator(prime).matrix()
            sign = ZZ(operator[0, 0])
            assert sign in [-1, 1]
            assert operator == sign * identity_matrix(
                operator.base_ring(), operator.nrows()
            )
            classical_signs.append(int(sign))

        cremona_label = None
        if factor_space.dimension() == 1:
            eigenform = factor_space.q_eigenform(12)
            ap_signature = tuple(
                ZZ(eigenform[prime]) for prime in comparison_primes
            )
            cremona_label = signature_to_label[ap_signature]

        newform_orbits.append(
            {
                "dimension": int(factor_space.dimension()),
                "classical_atkin_lehner_signs_W2_W3_W103": classical_signs,
                "cremona_label": cremona_label,
            }
        )

    assert sum(row["dimension"] for row in newform_orbits) == 17
    stable_orbits = [
        row
        for row in newform_orbits
        if prod(row["classical_atkin_lehner_signs_W2_W3_W103"]) == 1
    ]
    assert sum(row["dimension"] for row in stable_orbits) == 6
    assert all(row["dimension"] == 1 for row in stable_orbits)
    stable_factor_labels = sorted(row["cremona_label"] for row in stable_orbits)
    assert stable_factor_labels == [f"618{letter}1" for letter in "abcdef"]

    # On the stable curve the deck involution of C_1236 -> B is geometric
    # w_2=-W_2.  Its invariant part is therefore classical W_2=-1, and its
    # anti-invariant Prym is classical W_2=+1.
    genus_two_factor_labels = sorted(
        row["cremona_label"]
        for row in stable_orbits
        if row["classical_atkin_lehner_signs_W2_W3_W103"][0] == -1
    )
    prym_factor_labels = sorted(
        row["cremona_label"]
        for row in stable_orbits
        if row["classical_atkin_lehner_signs_W2_W3_W103"][0] == 1
    )
    assert genus_two_factor_labels == ["618e1", "618f1"]
    assert prym_factor_labels == ["618a1", "618b1", "618c1", "618d1"]
    stable_factor_ranks = {
        label: int(EllipticCurve(label).rank()) for label in stable_factor_labels
    }
    assert set(stable_factor_ranks.values()) == {1}
    stable_jacobian_rank = sum(stable_factor_ranks.values())
    assert stable_jacobian_rank == 6

    # Gonzalez--Rotger's residue-field formula at discriminant -3:
    # D(R)=2, N*(R)=103, m_R=3, Q=206=D(R)N*(R), and h(-3)=1.
    # Hence every such CM image on the w_618 quotient has residue field QQ.
    cm_minus_three = {
        "order_discriminant": -3,
        "class_number": order_class_number(-3),
        "D_R": 2,
        "N_star_R": 103,
        "m": 618,
        "m_R": 3,
        "Q": 206,
        "top_curve_elliptic_points": 4,
        "marked_curve_rational_points": 2,
        "residue_field": "QQ",
        "k3_interpretation": (
            "CM/rho-20 specialization; it does not realize geometric NS of "
            "rank 19 exactly."
        ),
    }
    assert cm_minus_three["class_number"] == 1
    assert cm_minus_three["Q"] == cm_minus_three["D_R"] * cm_minus_three["N_star_R"]
    cm_auxiliary_quaternion = QuaternionAlgebra(QQ, -3, 206)
    assert cm_auxiliary_quaternion.discriminant() == QUATERNION_DISCRIMINANT
    assert cm_auxiliary_quaternion.ramified_primes() == [2, 3]

    # The other w_3-fixed point downstairs comes from the four discriminant
    # -24 points fixed by w_6 upstairs.  Here D(R)=1, N*(R)=103, h(R)=2,
    # m_R=gcd(618,24)=6, and m/m_R=103.  Corollary 5.14 therefore gives a
    # quadratic, not rational, residue field on the w_618 quotient: H_R has
    # degree four over QQ and the specified involution has a degree-two fixed
    # field.  The two geometric points form one quadratic closed point.
    cm_minus_twenty_four = {
        "order_discriminant": -24,
        "class_number": order_class_number(-24),
        "D_R": 1,
        "N_star_R": 103,
        "m": 618,
        "m_R": 6,
        "Q": 103,
        "top_curve_w6_fixed_points": 4,
        "marked_curve_geometric_points": 2,
        "marked_curve_rational_points": 0,
        "ring_class_field_degree_over_QQ": 4,
        "marked_curve_residue_degree_over_QQ": 2,
        "k3_interpretation": (
            "CM/rho-20 quadratic closed point; it supplies no rational "
            "rank-19 marking."
        ),
    }
    assert cm_minus_twenty_four["class_number"] == 2
    assert cm_minus_twenty_four["Q"] == (
        cm_minus_twenty_four["D_R"] * cm_minus_twenty_four["N_star_R"]
    )
    assert cm_minus_twenty_four["m"] // cm_minus_twenty_four["m_R"] == 103
    cm_minus_twenty_four_auxiliary = QuaternionAlgebra(QQ, -6, 206)
    assert cm_minus_twenty_four_auxiliary.discriminant() == QUATERNION_DISCRIMINANT
    assert cm_minus_twenty_four_auxiliary.ramified_primes() == [2, 3]

    return {
        "schema": "elkies-k3.det1236-marked-shimura-curve.v1",
        "status": "UNRESOLVED_FOR_EXPLICIT_REASON",
        "surface_id": SURFACE_ID,
        "determinant": DETERMINANT,
        "independence_boundary": {
            "excluded_inputs": [
                "curve 356",
                "curve 385",
                "the frozen prospective experiment",
            ],
            "statement": (
                "The replay reads only the lattice catalogue and its T-arithmetic "
                "ledger; the low-genus equation is an independently published "
                "Shimura-quotient model."
            ),
        },
        "inputs": {
            relative(CATALOGUE): digest(CATALOGUE),
            relative(T_ARITHMETIC): digest(T_ARITHMETIC),
        },
        "rootless_frame_control": {
            "frame_id": FRAME_ID,
            "frame_gram_sha256": FRAME_SHA256,
            "rank": 17,
            "minimum_squared_norm": 4,
            "root_count": 0,
            "boundary": (
                "This is imported only to identify the selected arithmetic row; "
                "it is not used to infer a rational marking."
            ),
        },
        "transcendental_lattice": {
            "gram": rows(transcendental),
            "signature": [2, 1],
            "determinant": int(transcendental.det()),
            "content": 1,
            "smith_form": rows(smith),
            "discriminant_group": "Z/1236Z",
            "dual_generator": [str(entry) for entry in dual_generator],
            "discriminant_bilinear_value_mod_ZZ": str(
                discriminant_bilinear_value
            ),
            "discriminant_quadratic_value_mod_ZZ": str(
                discriminant_quadratic_value
            ),
            "orthogonal_group_units": discriminant_orthogonal_units,
        },
        "literal_even_clifford_order": {
            "basis": ["1", "e0*e1", "e0*e2", "e1*e2"],
            "reduced_trace_pairing": rows(reduced_trace_pairing),
            "reduced_discriminant": REDUCED_ORDER_DISCRIMINANT,
            "quaternion_hilbert_symbol": ["2", "309/4"],
            "quaternion_ramified_primes": [2, 3],
            "quaternion_discriminant": QUATERNION_DISCRIMINANT,
            "local_order_type": {
                "status": "PASS_EXACT_EICHLER_ORDER",
                "level": EICHLER_LEVEL,
                "explanation": (
                    "The literal order has reduced discriminant 2*3*103 in the "
                    "quaternion algebra ramified exactly at 2 and 3. It is maximal "
                    "at the ramified primes and locally the index-103 Borel/Eichler "
                    "order at the sole split discriminant prime 103."
                ),
            },
            "similarity_gap": "none: the literal content is one",
        },
        "normalizer_discriminant_action": {
            "normalizer_generators": normalizer_records,
            "atkin_lehner_actions_on_A_T": {
                f"w_{label}": unit for label, unit in atkin_actions.items()
            },
            "image_order": len(set(atkin_actions.values())),
            "image_equals_O_A_T": True,
            "central_inversion_action": 1235,
            "projectively_stable_atkin_lehner_labels": projectively_stable_atkin_labels,
        },
        "exact_marked_period_curve": {
            "status": "PASS_EXACT_PROJECTIVE_STABLE_DISCRIMINANT_KERNEL_CURVE",
            "label": "X_0^6(103)/<w_618>",
            "group": (
                "projective image of kernel(O^+(T)->O(A_T)); equivalently "
                "Gamma_0^6(103) extended by the projective class of (-1)*w_618"
            ),
            "degree_from_norm_one_curve": 2,
            "signature": stable_signature,
            "cuspidal_rational_points": 0,
            "certified_rational_cm_points": cm_minus_three,
            "certified_nonrational_cm_points": cm_minus_twenty_four,
            "rational_non_cm_points": "UNRESOLVED",
            "points_realizing_only_a_higher_picard_specialization": 2,
            "boundary": (
                "The two certified rational points are the discriminant -3 CM "
                "orbifold points. They have rho 20 and therefore do not realize "
                "the requested saturated rank-19 NS marking."
            ),
        },
        "coarse_norm_one_curve": {
            "label": "X_0^6(103)",
            "genus": int(top_genus),
            "mu": int(mu),
            "cusps": 0,
            "elliptic_points_order_2": int(elliptic_order_2),
            "elliptic_points_order_3": int(elliptic_order_3),
        },
        "ogg_fixed_point_counts": {
            f"w_{label}": record for label, record in fixed_records.items()
        },
        "quotient_tower": quotient_records,
        "genus_two_quotient": {
            "label": "X_0^6(103)/<w_2,w_309>",
            "map_from_exact_marked_curve": {
                "degree": 2,
                "deck_involution": "w_2 modulo <w_618>",
                "geometric_branch_points": 6,
                "branch_order_discriminant": -1236,
                "branch_source": (
                    "The twelve w_309-fixed points on X_0^6(103) are exactly "
                    "the order-discriminant -1236 CM points; quotienting by "
                    "w_618 pairs them into the six branch points on B."
                ),
            },
            "model": "y^2 = 1944*x^6 + 441*x^4 - 90*x^2 + 9",
            "verified_rational_points": rational_points_on_genus_two_target,
            "verified_point_images_on_618f1": {
                "(0,+/-3)": "+/-3*G",
                "(+/-1,+/-48)": "+/-4*G",
                "(+/-1/3,-/+8/3)": "+/-1*G",
                "(+/-1/5,-/+312/125)": "+/-10*G",
                "generator_G": "(10,-29)",
            },
            "w3_fixed_fiber_separation": {
                "fixed_rational_points": ["(0,3)", "(0,-3)"],
                "identification_up_to_hyperelliptic_sign": [
                    (
                        "one point is the image of the two rational "
                        "discriminant -3 CM points on C_1236"
                    ),
                    (
                        "the other is the image of the quadratic conjugate "
                        "discriminant -24 CM points on C_1236"
                    ),
                ],
                "boundary": (
                    "The published model fixes w_3 as x -> -x but does not "
                    "choose which y-sign labels the two CM image classes."
                ),
            },
            "model_source": {
                "paper": "Padurariu--Saia, Shimura curve Atkin--Lehner quotients of genus at most two",
                "arxiv": "https://arxiv.org/abs/2509.25368",
                "code": "https://github.com/fsaia/GenusAtMost2",
                "code_commit": LOW_GENUS_SOURCE_COMMIT,
            },
        },
        "jacquet_langlands_cover_precheck": {
            "status": "PASS_EXACT_NEWFORM_AND_PRYM_ACCOUNTING",
            "classical_level": 618,
            "classical_newspace_dimension": 17,
            "geometric_sign_conversion": {
                "w_2": "-W_2",
                "w_3": "-W_3",
                "w_103": "W_103",
                "w_618": "W_2*W_3*W_103",
            },
            "all_newform_orbits": newform_orbits,
            "stable_curve_jacobian_isogenous_factors": stable_factor_labels,
            "stable_curve_jacobian_rank": stable_jacobian_rank,
            "genus_two_quotient_jacobian_isogenous_factors": genus_two_factor_labels,
            "cover_prym_isogenous_factors": prym_factor_labels,
            "factor_ranks": stable_factor_ranks,
            "classical_chabauty": {
                "status": "DOES_NOT_PASS_STRICT_RANK_BOUND",
                "jacobian_rank": stable_jacobian_rank,
                "curve_genus": 6,
            },
            "quadratic_chabauty_dimension_screen": {
                "status": "PASSES_NECESSARY_DIMENSION_INEQUALITY",
                "neron_severi_rank_lower_bound": 6,
                "inequality": "rank J(QQ)=6 < genus 6 + rho(J)-1, with rho(J)>=6",
                "boundary": (
                    "This screen does not determine C_1236(QQ). It becomes an "
                    "actionable quadratic-Chabauty route only after an explicit "
                    "model of the degree-two marking cover is constructed."
                ),
            },
        },
        "full_atkin_lehner_quotient": {
            "label": "X_0^6(103)^*",
            "genus": 1,
            "cremona_label": "618f1",
            "minimal_model": "y^2 + x*y = x^3 - 185*x + 1401",
            "rank": 1,
            "torsion_order": 1,
            "generator": [10, -29],
            "map_from_genus_two_model": [
                "X = 54*x^2 + 4",
                "Y = 9*y - 27*x^2 - 2",
            ],
            "point_images": {
                "(0,3)": "3*(10,-29)",
                "(1,48)": "4*(10,-29)",
            },
            "boundary": (
                "Positive rank on this degree-four quotient of the exact marked "
                "curve does not imply a rational lift to the marked curve."
            ),
        },
        "arithmetic_certificate": {
            "status": "UNRESOLVED_FOR_EXPLICIT_REASON",
            "precise_obstruction": (
                "Determine whether C_1236(QQ) contains a non-CM point. Equivalently, "
                "construct the degree-two cover C_1236 -> "
                "X_0^6(103)/<w_2,w_309> and decide which rational points on the "
                "displayed genus-two model have rational lifts. No equation or "
                "descent class for that cover is presently certified."
            ),
            "why_lower_quotient_points_do_not_decide_it": (
                "The fourteen verified rational points on the genus-two quotient, and "
                "the positive-rank elliptic quotient 618f1, live below a nontrivial "
                "degree-two marking cover. Their existence supplies no rational point "
                "upstairs outside the two separately certified discriminant -3 CM points."
            ),
            "positive_next_gate": (
                "One exact non-CM rational lift on C_1236 would pass the period-curve "
                "gate for this content-one cyclic marking and authorize construction "
                "of NS=T^perp and the already catalogued rootless frames."
            ),
            "negative_next_gate": (
                "A complete proof that the only rational points of C_1236 are its "
                "two discriminant -3 CM points would arithmetically exclude the row."
            ),
        },
        "theorem_inputs": [
            "the local classification of squarefree-discriminant quaternion orders and the Eichler normalizer quotient",
            "Ogg's Atkin--Lehner fixed-point and Riemann--Hurwitz formulas",
            "Gonzalez--Rotger's CM residue-field formula for Atkin--Lehner quotients",
            "Padurariu--Saia's exact genus-two quotient model and elliptic quotient identification",
            "Jacquet--Langlands and the ramified-place Atkin--Lehner sign normalization",
            "the rank-three ternary-spin/fully-marked-K3 period correspondence",
        ],
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_det1236_marked_shimura_curve.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    catalogue = json.loads(CATALOGUE.read_text())
    t_arithmetic = json.loads(T_ARITHMETIC.read_text())
    payload = build_payload(catalogue, t_arithmetic)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file():
            raise FileNotFoundError(args.output)
        if args.output.read_text() != rendered:
            raise AssertionError(f"generated artifact changed: {args.output}")
        print(
            json.dumps(
                {
                    "status": "PASS_DET1236_MARKED_SHIMURA_CURVE_CHECK",
                    "output": relative(args.output),
                    "sha256": digest(args.output),
                    "arithmetic_certificate": payload["status"],
                    "exact_marked_curve": payload["exact_marked_period_curve"]["label"],
                    "marked_genus": payload["exact_marked_period_curve"]["signature"]["genus"],
                },
                sort_keys=True,
            )
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        json.dumps(
            {
                "status": "WROTE_DET1236_MARKED_SHIMURA_CURVE_CERTIFICATE",
                "output": relative(args.output),
                "sha256": digest(args.output),
                "arithmetic_certificate": payload["status"],
                "exact_marked_curve": payload["exact_marked_period_curve"]["label"],
                "marked_genus": payload["exact_marked_period_curve"]["signature"]["genus"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
