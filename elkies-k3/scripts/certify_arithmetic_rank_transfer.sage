#!/usr/bin/env sage-python
"""Certify equivariant rank transfer and the alternate-Q80 arithmetic rank.

status: ACTIVE_PROOF
claim: exact finite Galois-module rank transfer, arithmetic rank 17 for the
       norm12-orbit-11952 alternate-Q80 pencil, H3/E6 controls, and a
       fail-closed NS0024 arithmetic-promotion gate
inputs: the pinned R17 Gram, exact rational H3 sections and bisections, and
        existing H3, Q80, E6, and NS0024 certificates
outputs: elkies-k3-arithmetic-rank-transfer-controls-v1.json

The reusable verifier works with a common geometric Neron--Severi lattice,
finite integral Galois matrices, and for each fibration a marked hyperbolic
plane and geometric fibre-root basis.  It checks the action, computes fixed
ranks exactly, and verifies the representation-ring rank-transfer identity by
traces on every element of the generated finite group.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from jsonschema import Draft202012Validator
from sage.all import (
    PolynomialRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
MARKING_SCHEMA = ROOT / "elkies-k3/data/arithmetic/arithmetic-marking-v1.schema.json"
R17_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
H3_ENDPOINT = GEN / "elkies-k3-h3-q12o5867-endpoint-certificate.json"
H3_PUBLISHED_TARGET = GEN / "elkies-2026-published-r17-target.json"
H3_PUBLISHED_COORDINATE_MATCH = (
    GEN / "elkies-k3-h3-q12o5867-elkies-2026-coordinate-match.json"
)
H3_PUBLISHED_SECTIONS = (
    ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
)
E6_INCIDENCE = GEN / "elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
E6_ORBIT103 = (
    GEN / "elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json"
)
NS0024_ROUTE = GEN / "elkies-k3-ns0024-new-rootless-source-route-v1.json"
R17_BISECTION_SPLITTING = (
    GEN / "elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
R17_NORM12_CLASSIFICATION = (
    GEN / "elkies-k3-r17-norm12-isotropic-frame-classification-v1.json"
)
ALTERNATE_Q80_FRAME = GEN / "q80-alternate-fifth-q6-rootless-transport.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-arithmetic-rank-transfer-controls-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_matrix_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def matrix_key(value):
    return tuple(int(entry) for entry in value.list())


def finite_group_closure(generators, rank, maximum_order=4096):
    identity = identity_matrix(ZZ, rank)
    elements = {matrix_key(identity): identity}
    frontier = [identity]
    while frontier:
        left = frontier.pop()
        for right in generators:
            product = left * right
            key = matrix_key(product)
            if key not in elements:
                elements[key] = product
                frontier.append(product)
                if len(elements) > maximum_order:
                    raise ArithmeticError("declared Galois image exceeded closure bound")
    return [elements[key] for key in sorted(elements)]


def independent_row_indices(columns):
    if columns.ncols() == 0:
        return []
    return list(columns.transpose().pivots())


def restricted_action(action, basis, label):
    if basis.ncols() == 0:
        return matrix(QQ, 0, 0)
    rows = independent_row_indices(basis)
    square = basis.matrix_from_rows(rows)
    image = (action * basis).matrix_from_rows(rows)
    coordinates = square.inverse() * image
    if basis * coordinates != action * basis:
        raise ArithmeticError(f"{label} is not stable under the Galois action")
    return matrix(QQ, coordinates)


def common_fixed_basis(actions, rank):
    if not actions:
        return identity_matrix(ZZ, rank)
    identity = identity_matrix(ZZ, rank)
    rows = []
    for action in actions:
        rows.extend((action - identity).rows())
    equations = matrix(ZZ, rows)
    return equations.right_kernel().basis_matrix()


def common_fixed_rational_basis(actions, rank):
    if not actions:
        return identity_matrix(QQ, rank)
    identity = identity_matrix(QQ, rank)
    rows = []
    for action in actions:
        rows.extend((matrix(QQ, action) - identity).rows())
    return matrix(QQ, rows).right_kernel().basis_matrix()


def extend_to_basis(columns, ambient_rank):
    chosen = matrix(QQ, columns)
    for index in range(ambient_rank):
        if chosen.ncols() == ambient_rank:
            break
        candidate = matrix(QQ, ambient_rank, 1, lambda row, _column: int(row == index))
        trial = chosen.augment(candidate)
        if trial.rank() > chosen.rank():
            chosen = trial
    if chosen.ncols() != ambient_rank or chosen.det() == 0:
        raise ArithmeticError("could not extend the root space to a rational frame basis")
    return chosen


def validate_marking(marking):
    schema = json.loads(MARKING_SCHEMA.read_text())
    Draft202012Validator(schema).validate(marking)
    gram = matrix(ZZ, marking["gram"])
    if not gram.is_square() or gram.det() == 0:
        raise ArithmeticError(f"{marking['id']} has a degenerate NS Gram matrix")
    rank = gram.nrows()
    generators = [matrix(ZZ, item["matrix"]) for item in marking["galois_generators"]]
    for action in generators:
        if action.nrows() != rank or action.ncols() != rank:
            raise ArithmeticError(f"{marking['id']} has a wrong-size Galois matrix")
        if abs(action.det()) != 1 or action.transpose() * gram * action != gram:
            raise ArithmeticError(f"{marking['id']} has a nonintegral lattice action")
    group = finite_group_closure(generators, rank)
    fixed_ns_basis = common_fixed_basis(group, rank)
    fixed_ns_rank = fixed_ns_basis.nrows()
    fibrations = []
    character_rows = {}
    for item in marking["fibrations"]:
        u_basis = matrix(ZZ, item["u_basis_columns"]).transpose()
        root_columns = item["root_basis_columns"]
        root_basis = (
            matrix(ZZ, root_columns).transpose()
            if root_columns
            else matrix(ZZ, rank, 0)
        )
        if u_basis.nrows() != rank or u_basis.ncols() != 2:
            raise ArithmeticError(f"{item['id']} has a wrong-size U basis")
        hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
        if u_basis.transpose() * gram * u_basis != hyperbolic:
            raise ArithmeticError(f"{item['id']} does not span the marked U")
        for action in group:
            if action * u_basis != u_basis:
                raise ArithmeticError(f"{item['id']} U is not fixed pointwise")
        if root_basis.nrows() != rank:
            raise ArithmeticError(f"{item['id']} has a wrong-size root basis")
        if root_basis.ncols() and u_basis.transpose() * gram * root_basis:
            raise ArithmeticError(f"{item['id']} root basis is not orthogonal to U")
        if root_basis.rank() != root_basis.ncols():
            raise ArithmeticError(f"{item['id']} root basis is dependent")
        root_actions = [
            restricted_action(action, root_basis, f"{item['id']} root space")
            for action in group
        ]
        fixed_root_rank = common_fixed_rational_basis(
            root_actions, root_basis.ncols()
        ).nrows()
        geometric_rank = rank - 2 - root_basis.ncols()
        arithmetic_rank = fixed_ns_rank - 2 - fixed_root_rank
        w_basis = (u_basis.transpose() * gram).right_kernel().basis_matrix().transpose()
        w_actions = [
            restricted_action(action, w_basis, f"{item['id']} frame space")
            for action in group
        ]
        if root_basis.ncols():
            rows = independent_row_indices(w_basis)
            root_coordinates = (
                w_basis.matrix_from_rows(rows).inverse()
                * root_basis.matrix_from_rows(rows)
            )
            if w_basis * root_coordinates != root_basis:
                raise ArithmeticError(f"{item['id']} roots do not lie in the frame")
        else:
            root_coordinates = matrix(QQ, w_basis.ncols(), 0)
        full_frame_basis = extend_to_basis(root_coordinates, w_basis.ncols())
        root_rank = root_basis.ncols()
        quotient_actions = []
        for frame_action in w_actions:
            conjugated = full_frame_basis.inverse() * frame_action * full_frame_basis
            if root_rank and conjugated[root_rank:, :root_rank] != 0:
                raise ArithmeticError(f"{item['id']} quotient action is not well-defined")
            quotient_actions.append(conjugated[root_rank:, root_rank:])
        fixed_mw_basis = common_fixed_rational_basis(quotient_actions, geometric_rank)
        if fixed_mw_basis.nrows() != arithmetic_rank:
            raise ArithmeticError(f"{item['id']} quotient fixed-space rank disagrees with A2.2")
        trivial_basis = u_basis.augment(root_basis)
        full_ns_basis = extend_to_basis(trivial_basis, rank)
        ns_quotient_actions = []
        for action in group:
            conjugated = full_ns_basis.inverse() * action * full_ns_basis
            trivial_rank = trivial_basis.ncols()
            if trivial_rank and conjugated[trivial_rank:, :trivial_rank] != 0:
                raise ArithmeticError(f"{item['id']} NS quotient action is not well-defined")
            ns_quotient_actions.append(conjugated[trivial_rank:, trivial_rank:])
        if [action.trace() for action in ns_quotient_actions] != [
            action.trace() for action in quotient_actions
        ]:
            raise ArithmeticError(f"{item['id']} two MW quotient constructions disagree")
        section_orbits = []
        for section in item.get("sections", []):
            section_class = matrix(ZZ, rank, 1, section["ns_class"])
            coordinates = full_ns_basis.inverse() * section_class
            quotient_coordinates = vector(
                QQ, coordinates.column(0)[trivial_basis.ncols():]
            )
            orbit = {
                tuple(action * quotient_coordinates)
                for action in ns_quotient_actions
            }
            stabilizer_order = sum(
                int(action * quotient_coordinates == quotient_coordinates)
                for action in ns_quotient_actions
            )
            orbit_size = len(orbit)
            if orbit_size != section["expected_orbit_size"]:
                raise ArithmeticError(f"{item['id']} section orbit changed for {section['label']}")
            if stabilizer_order * orbit_size != len(group):
                raise ArithmeticError(f"{item['id']} section orbit-stabilizer check failed")
            section_orbits.append(
                {
                    "label": section["label"],
                    "quotient_coordinates": [str(entry) for entry in quotient_coordinates],
                    "orbit_size": orbit_size,
                    "stabilizer_order": stabilizer_order,
                    "declared_field_of_definition": section["field_of_definition"],
                }
            )
        expected = item["expected"]
        if geometric_rank != expected["geometric_mw_rank"]:
            raise ArithmeticError(f"{item['id']} geometric MW rank changed")
        if arithmetic_rank != expected["arithmetic_mw_rank"]:
            raise ArithmeticError(f"{item['id']} arithmetic MW rank changed")
        traces = []
        for action, root_action in zip(group, root_actions):
            trace_ns = int(action.trace())
            trace_root = int(root_action.trace())
            traces.append(
                {
                    "trace_ns": trace_ns,
                    "trace_root": trace_root,
                    "trace_mw": trace_ns - 2 - trace_root,
                }
            )
        character_rows[item["id"]] = [row["trace_mw"] for row in traces]
        fibrations.append(
            {
                "id": item["id"],
                "root_rank": root_basis.ncols(),
                "fixed_root_rank": fixed_root_rank,
                "geometric_mw_rank": geometric_rank,
                "arithmetic_mw_rank": arithmetic_rank,
                "rational_fixed_mw_basis_in_quotient_coordinates": rational_matrix_rows(
                    fixed_mw_basis
                ),
                "galois_action_on_root_basis": [
                    rational_matrix_rows(action) for action in root_actions
                ],
                "galois_action_on_mw_quotient": [
                    rational_matrix_rows(action) for action in ns_quotient_actions
                ],
                "component_labels": item.get("component_labels", []),
                "section_orbits": section_orbits,
                "character_traces": traces,
            }
        )
    edges = []
    by_id = {row["id"]: row for row in fibrations}
    for edge in marking.get("edges", []):
        source = by_id[edge["source"]]
        target = by_id[edge["target"]]
        geometric_delta = target["geometric_mw_rank"] - source["geometric_mw_rank"]
        arithmetic_delta = target["arithmetic_mw_rank"] - source["arithmetic_mw_rank"]
        if geometric_delta != source["root_rank"] - target["root_rank"]:
            raise ArithmeticError("geometric rank-transfer identity failed")
        if arithmetic_delta != source["fixed_root_rank"] - target["fixed_root_rank"]:
            raise ArithmeticError("arithmetic rank-transfer identity failed")
        source_character = character_rows[source["id"]]
        target_character = character_rows[target["id"]]
        source_roots = source["character_traces"]
        target_roots = target["character_traces"]
        for index in range(len(group)):
            lhs = target_character[index] - source_character[index]
            rhs = source_roots[index]["trace_root"] - target_roots[index]["trace_root"]
            if lhs != rhs:
                raise ArithmeticError("representation-ring trace identity failed")
        edges.append(
            {
                "source": source["id"],
                "target": target["id"],
                "geometric_rank_change": geometric_delta,
                "arithmetic_rank_change": arithmetic_delta,
                "representation_ring_identity": "PASS_ON_EVERY_GROUP_ELEMENT",
            }
        )
    return {
        "id": marking["id"],
        "ground_field": marking["ground_field"],
        "ns_rank": rank,
        "ns_determinant": int(gram.det()),
        "galois_image_order": len(group),
        "fixed_ns_rank": fixed_ns_rank,
        "fixed_ns_basis": matrix_rows(fixed_ns_basis),
        "fibrations": fibrations,
        "edges": edges,
    }


def h3_marking(endpoint, published_target, published_sections):
    r17 = matrix(ZZ, [
        [int(entry) for entry in line.split()]
        for line in R17_GRAM.read_text().splitlines()
        if line.strip()
    ])
    if r17.nrows() != 17 or r17.det() != 948:
        raise ArithmeticError("pinned R17 Gram changed")
    if endpoint.get("status") != "PASS_EXACT_Q12O5867_SOURCE_IDENTITY_RHO19_FULL_MW_R17":
        raise ArithmeticError("H3 endpoint certificate is not exact")
    if endpoint["picard_rank"]["geometric_picard_rank_characteristic_zero"] != 19:
        raise ArithmeticError("H3 geometric Picard rank changed")
    if endpoint["mordell_weil"]["full_geometric_mordell_weil_rank"] != 17:
        raise ArithmeticError("H3 geometric MW rank changed")
    if published_target.get("status") != "PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17":
        raise ArithmeticError("published H3 R17 certificate is not exact")
    if published_sections.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_SECTIONS_AND_CHORDS":
        raise ArithmeticError("published H3 section transcription is not exact")
    if len(published_sections["sections"]) != 17:
        raise ArithmeticError("published H3 rational section count changed")
    height = published_target["published_height_lattice"]
    equation = published_target["published_equation"]
    identification = published_target["pinned_identification"]
    if (
        height["rank"] != 17
        or height["determinant"] != 948
        or equation["published_section_identities"] != 17
        or not equation["rootless_semistable"]
        or abs(identification["basis_change_determinant"]) != 1
    ):
        raise ArithmeticError("published H3 rational R17 evidence changed")
    gram = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -r17)
    identity = identity_matrix(ZZ, 19)
    fibre = [1] + [0] * 18
    mate = [0, 1] + [0] * 17
    sections = []
    for index in range(17):
        mw_vector = [0] * 17
        mw_vector[index] = 1
        section_class = [int((r17[index, index] - 2) / 2), 1] + mw_vector
        sections.append(
            {
                "label": f"rational_pinned_R17_basis_{index + 1}",
                "ns_class": section_class,
                "expected_orbit_size": 1,
                "field_of_definition": "QQ(t)",
            }
        )
    return {
        "id": "H3_QQ_rootless_R17",
        "ground_field": "QQ",
        "gram": matrix_rows(gram),
        "galois_generators": [{"name": "identity", "matrix": matrix_rows(identity)}],
        "fibrations": [
            {
                "id": "rootless_R17",
                "u_basis_columns": [fibre, mate],
                "root_basis_columns": [],
                "component_labels": [],
                "sections": sections,
                "expected": {"geometric_mw_rank": 17, "arithmetic_mw_rank": 17},
            }
        ],
        "edges": [],
    }


def alternate_q80_arithmetic_application(
    endpoint,
    published_target,
    published_sections,
    coordinate_match,
    splitting,
    classification,
    alternate_payload,
):
    """Promote the degree-two alternate-Q80 copy in the rational R17 marking."""
    # This call checks the exact rank-19 source evidence before its rational
    # divisor basis is changed from split U+R17(-1) coordinates to
    # (F,O,Q_1,...,Q_17), where the Q_i are rational pinned-R17 sections.
    source_marking = h3_marking(endpoint, published_target, published_sections)
    split_gram = matrix(ZZ, source_marking["gram"])
    r17 = -split_gram[2:, 2:]

    if coordinate_match.get("status") != "PASS_EXACT_Q12O5867_IS_ELKIES_2026_PUBLISHED_MODEL":
        raise ArithmeticError("published R17 coordinate identification is not exact")

    if splitting.get("status") != "PASS_EXACT_SIMULTANEOUS_SPLIT_HITS_QUOTIENTED":
        raise ArithmeticError("R17 bisection source certificate is not exact")
    if classification.get("status") != "PASS_EXACT_MINIMAL_J2_ACCESSIBILITY":
        raise ArithmeticError("norm-twelve frame classification is not exact")
    if alternate_payload.get("schema") != "q80-alternate-fifth-q6-rootless-transport-v1":
        raise ArithmeticError("alternate-Q80 frame certificate is not exact")

    label = "norm12-orbit-11952"
    source_records = [
        row for row in splitting["construction"]["records"] if row["label"] == label
    ]
    classified_records = [
        row for row in classification["classification"]["records"]
        if row["label"] == label
    ]
    if len(source_records) != 1 or len(classified_records) != 1:
        raise ArithmeticError(f"{label} is not uniquely recorded")
    source_record = source_records[0]
    classified = classified_records[0]
    required_curve_gates = (
        "branch_polynomial_irreducible_over_Q",
        "branch_polynomial_squarefree",
        "branch_polynomial_coprime_to_surface_discriminant",
        "branch_polynomial_coprime_to_trace_denominator",
    )
    if not all(source_record[key] for key in required_curve_gates):
        raise ArithmeticError(f"{label} lost an exact QQ genus-one curve gate")
    if source_record["member_selection"] != "unique regular M0 member":
        raise ArithmeticError(f"{label} lost its regular member selection")
    if classified["frame_class"] != "alternate-Q80" or not classified["shared_zero"]:
        raise ArithmeticError(f"{label} is no longer the shared-zero alternate frame")

    w = vector(ZZ, source_record["pinned_rank17_w"])
    if list(w) != classified["trace_vector"] or w * r17 * w != 12:
        raise ArithmeticError(f"{label} has inconsistent pinned R17 coordinates")
    old_fibre = vector(ZZ, [1, 0] + [0] * 17)
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    fibre = vector(ZZ, [3, 2] + list(w))
    target_mate = fibre + old_zero
    target_u = matrix(ZZ, [fibre, target_mate]).transpose()
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    if target_u.transpose() * split_gram * target_u != hyperbolic:
        raise ArithmeticError(f"{label} does not span the claimed U")
    if fibre * split_gram * old_fibre != 2 or fibre * split_gram * old_zero != 1:
        raise ArithmeticError(f"{label} lost its degree-two shared-zero marking")

    complement = matrix(
        ZZ, [list(fibre * split_gram), list(target_mate * split_gram)]
    ).right_kernel_matrix()
    splitting_transport = matrix(
        ZZ,
        [list(fibre), list(target_mate)] + [list(row) for row in complement.rows()],
    )
    if abs(splitting_transport.det()) != 1:
        raise ArithmeticError(f"{label} target U is not primitive and integrally split")
    target_frame = -(complement * split_gram * complement.transpose())
    if target_frame.det() != 948 or not target_frame.is_positive_definite():
        raise ArithmeticError(f"{label} target frame has the wrong genus")
    if int(pari(target_frame).qfminim(2)[0]) != 0:
        raise ArithmeticError(f"{label} target frame is not rootless")
    alternate_frame = matrix(ZZ, alternate_payload["rootless_frame"])
    if pari(target_frame).qfisom(pari(alternate_frame)) == 0:
        raise ArithmeticError(f"{label} target frame is not alternate Q80")

    # Build a literal rational-divisor basis.  In split coordinates a section
    # Q_i with frame coordinate e_i is (a_i,1,e_i), where
    # a_i=(<Q_i,Q_i>-2)/2.  The determinant-one published-to-pinned
    # identification makes every Q_i an integral combination of the seventeen
    # published QQ(t)-sections.
    rational_columns = [old_fibre, old_zero]
    section_offsets = []
    for index in range(17):
        offset = ZZ((r17[index, index] - 2) / 2)
        section_offsets.append(int(offset))
        mw_vector = [0] * 17
        mw_vector[index] = 1
        rational_columns.append(vector(ZZ, [offset, 1] + mw_vector))
    split_from_rational = matrix(
        ZZ, 19, 19, lambda row, column: rational_columns[column][row]
    )
    if split_from_rational.det() != 1:
        raise ArithmeticError("the rational R17 divisor classes are not an integral NS basis")
    rational_gram = split_from_rational.transpose() * split_gram * split_from_rational
    rational_from_split = split_from_rational.inverse()
    rational_fibre = vector(ZZ, rational_from_split * fibre)
    rational_target_mate = vector(ZZ, rational_from_split * target_mate)
    expected_rational_fibre = vector(
        ZZ,
        [40, -1, -1, -1, 3, 0, 0, 1, -1, -1, 3, 1, 1, 0, 2, -1, 1, -2, -2],
    )
    if rational_fibre != expected_rational_fibre:
        raise ArithmeticError(f"{label} rational-divisor coordinates changed")

    source_u = [
        [1] + [0] * 18,
        [1, 1] + [0] * 17,
    ]
    sections = []
    for index in range(17):
        section_class = [0] * 19
        section_class[index + 2] = 1
        sections.append(
            {
                "label": f"rational_pinned_R17_section_Q{index + 1}",
                "ns_class": section_class,
                "expected_orbit_size": 1,
                "field_of_definition": "QQ(t)",
            }
        )
    marking = {
        "id": "published_R17_to_alternate_Q80_norm12_orbit_11952",
        "ground_field": "QQ",
        "gram": matrix_rows(rational_gram),
        "galois_generators": [
            {"name": "identity", "matrix": matrix_rows(identity_matrix(ZZ, 19))}
        ],
        "fibrations": [
            {
                "id": "published_R17",
                "u_basis_columns": source_u,
                "root_basis_columns": [],
                "component_labels": [],
                "sections": sections,
                "expected": {"geometric_mw_rank": 17, "arithmetic_mw_rank": 17},
            },
            {
                "id": "alternate_Q80_norm12_orbit_11952",
                "u_basis_columns": [
                    list(map(int, rational_fibre)),
                    list(map(int, rational_target_mate)),
                ],
                "root_basis_columns": [],
                "component_labels": [],
                "sections": [
                    {
                        "label": "shared_old_zero_O",
                        "ns_class": [0, 1] + [0] * 17,
                        "expected_orbit_size": 1,
                        "field_of_definition": "QQ(u)",
                    }
                ],
                "expected": {"geometric_mw_rank": 17, "arithmetic_mw_rank": 17},
            },
        ],
        "edges": [
            {"source": "published_R17", "target": "alternate_Q80_norm12_orbit_11952"}
        ],
    }
    gate_result = validate_marking(marking)
    by_id = {row["id"]: row for row in gate_result["fibrations"]}
    target_result = by_id["alternate_Q80_norm12_orbit_11952"]
    if (
        gate_result["fixed_ns_rank"] != 19
        or target_result["root_rank"] != 0
        or target_result["arithmetic_mw_rank"] != 17
    ):
        raise ArithmeticError("alternate-Q80 arithmetic promotion gate failed")

    branch_coefficients = source_record[
        "branch_polynomial_q_coefficients_low_to_high"
    ]
    if len(branch_coefficients) != 5:
        raise ArithmeticError(f"{label} branch polynomial is no longer quartic")
    polynomial_ring = PolynomialRing(QQ, "t")
    branch_polynomial = polynomial_ring(
        [QQ(coefficient) for coefficient in branch_coefficients]
    )
    if (
        branch_polynomial.degree() != 4
        or not branch_polynomial.is_irreducible()
        or not branch_polynomial.is_squarefree()
    ):
        raise ArithmeticError(f"{label} branch polynomial failed exact replay")
    for block, names in (
        (source_record["lifted_section"], ("x0", "x1", "y0", "y1")),
        (source_record["trace_section"], ("M0", "Nx", "Ny", "h")),
    ):
        for name in names:
            coefficients = block[f"{name}_coefficients_low_to_high"]
            if not coefficients:
                raise ArithmeticError(f"{label} lost the rational map block {name}")
            for coefficient in coefficients:
                QQ(coefficient)
    return {
        "status": "PASS_EXACT_ALTERNATE_Q80_ARITHMETIC_RANK_17_BEFORE_EQUATION_COMPILATION",
        "witness": label,
        "source_rational_ns_basis": {
            "basis": ["F", "O"] + [f"Q{index + 1}" for index in range(17)],
            "rank": 19,
            "gram_determinant": int(rational_gram.det()),
            "split_from_rational_determinant": int(split_from_rational.det()),
            "pinned_section_offsets_in_split_basis": section_offsets,
            "published_to_pinned_section_basis": published_target[
                "pinned_identification"
            ],
            "field_of_definition": "QQ",
        },
        "target_marking": {
            "split_basis": ["F", "O+F"] + [f"w{index + 1}" for index in range(17)],
            "fibre_in_split_basis": list(map(int, fibre)),
            "rational_divisor_basis": ["F", "O"]
            + [f"Q{index + 1}" for index in range(17)],
            "fibre_D_in_rational_divisor_basis": list(map(int, rational_fibre)),
            "mate_O_plus_D_in_rational_divisor_basis": list(
                map(int, rational_target_mate)
            ),
            "shared_zero_in_rational_divisor_basis": [0, 1] + [0] * 17,
            "old_fibre_degree": 2,
            "old_zero_degree": 1,
            "primitive_U_splitting_determinant": int(splitting_transport.det()),
            "frame_determinant": int(target_frame.det()),
            "frame_root_count": 0,
            "frame_class": "alternate-Q80",
        },
        "descent_evidence": {
            "source": (
                "The published fibration, zero, and seventeen rational sections form a determinant-one "
                "integral basis of the geometric rank-19 Neron-Severi lattice."
            ),
            "fibre": (
                "The exact norm12-orbit-11952 record supplies a unique regular irreducible smooth "
                "genus-one curve over QQ in class D; hence D is an effective QQ-divisor and is nef."
            ),
            "zero": (
                "The published zero O is a QQ-curve and D.O=1, so it is the QQ-rational zero section "
                "of the pencil |D|."
            ),
            "pencil": (
                "The primitive nef isotropic QQ-divisor D has h0=2 on the K3, so |D| is a QQ-defined "
                "Jacobian elliptic pencil even though its base coordinate and Weierstrass equation are not compiled."
            ),
        },
        "arithmetic_marking": marking,
        "gate_result": gate_result,
        "conclusion": {
            "geometric_picard_rank": 19,
            "fixed_ns_rank": 19,
            "geometric_root_rank": 0,
            "arithmetic_mordell_weil_rank": 17,
            "field": "QQ(u), equivalently QQ(t) after naming the new pencil coordinate t",
        },
        "proof_boundary": (
            "This proves exact arithmetic generic rank 17 for the alternate-Q80 fibration without "
            "recovering its endpoint sections or compiling its Weierstrass equation. It does not give "
            "the new base parameter, equation, individual sections, or integral Mordell-Weil Gram in "
            "an equation-side section basis."
        ),
    }


def e6_incidence_marking(incidence):
    if incidence.get("status") != "PASS_EXACT_E6_RANK4_INCIDENCE_DESCENT":
        raise ArithmeticError("E6 incidence certificate is not exact")
    ns = incidence["neron_severi"]
    mw = incidence["mordell_weil"]
    if ns["generic_picard_rank"] != 19 or mw["geometric_rank"] != 4:
        raise ArithmeticError("E6 incidence ranks changed")
    gram = matrix(ZZ, ns["integral_gram"])
    action = identity_matrix(ZZ, 19)
    action.swap_columns(15, 16)
    action.swap_columns(17, 18)
    fibre = [0, 1] + [0] * 17
    zero_plus_fibre = [1, 1] + [0] * 17
    roots = []
    for index in range(2, 15):
        root = [0] * 19
        root[index] = 1
        roots.append(root)
    section_records = []
    for label, index in zip(["P", "Q", "R1", "R2"], range(15, 19)):
        section_class = [0] * 19
        section_class[index] = 1
        section_records.append(
            {
                "label": label,
                "ns_class": section_class,
                "expected_orbit_size": 2,
                "field_of_definition": "QQ(k)(r), r^2=k^4+6*k^2+13",
            }
        )
    return {
        "id": "E6_rank4_unordered_incidence",
        "ground_field": "QQ(k)",
        "gram": matrix_rows(gram),
        "galois_generators": [{"name": "r_to_minus_r", "matrix": matrix_rows(action)}],
        "fibrations": [
            {
                "id": "two_IVstar_I2",
                "u_basis_columns": [fibre, zero_plus_fibre],
                "root_basis_columns": roots,
                "component_labels": ns["basis"][2:15],
                "sections": section_records,
                "expected": {"geometric_mw_rank": 4, "arithmetic_mw_rank": 2},
            }
        ],
        "edges": [],
    }


def orbit103_module_control(audit):
    if audit.get("status") != "PASS_EXACT_BOTH_ORBITS_ARITHMETIC_RANK2_AND_CHI_MINUS3":
        raise ArithmeticError("orbit-103 arithmetic certificate is not exact")
    row = audit["orbit103_arithmetic"]
    action = matrix(ZZ, [[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    group = finite_group_closure([action], 3)
    fixed = common_fixed_basis(group, 3)
    if row["geometric_rank"] != 3 or row["arithmetic_rank"] != fixed.nrows():
        raise ArithmeticError("orbit-103 fixed-space control changed")
    return {
        "id": "E6A1_orbit103_MW_2_plus_chi_minus3",
        "ground_field": row["base_field"],
        "geometric_mw_rank": 3,
        "arithmetic_mw_rank": fixed.nrows(),
        "galois_image_order": len(group),
        "action_on_mw_basis": matrix_rows(action),
        "fixed_basis": matrix_rows(fixed),
        "section_fields": {
            "Q_plus": row["base_field"],
            "Q_minus": row["base_field"],
            "Q_delta": row["third_geometric_direction"]["field"],
        },
        "status": "PASS_EXACT_2_TRIVIAL_PLUS_CHI_MINUS3",
    }


def abstract_transfer_marking():
    """A representation-only regression with fixed and anti-invariant root loss."""
    gram = block_diagonal_matrix(
        matrix(ZZ, [[0, 1], [1, 0]]),
        -2 * identity_matrix(ZZ, 17),
    )
    action = identity_matrix(ZZ, 19)
    action[18, 18] = -1
    fibre = [1] + [0] * 18
    mate = [0, 1] + [0] * 17

    def roots(indices):
        result = []
        for index in indices:
            root = [0] * 19
            root[index] = 1
            result.append(root)
        return result

    return {
        "id": "abstract_C2_root_character_transfer_regression",
        "ground_field": "abstract characteristic-zero field",
        "gram": matrix_rows(gram),
        "galois_generators": [{"name": "sigma", "matrix": matrix_rows(action)}],
        "fibrations": [
            {
                "id": "R_13trivial_plus_sign",
                "u_basis_columns": [fibre, mate],
                "root_basis_columns": roots(list(range(2, 15)) + [18]),
                "expected": {"geometric_mw_rank": 3, "arithmetic_mw_rank": 3},
            },
            {
                "id": "R_13trivial",
                "u_basis_columns": [fibre, mate],
                "root_basis_columns": roots(list(range(2, 15))),
                "expected": {"geometric_mw_rank": 4, "arithmetic_mw_rank": 3},
            },
            {
                "id": "R_12trivial",
                "u_basis_columns": [fibre, mate],
                "root_basis_columns": roots(list(range(2, 14))),
                "expected": {"geometric_mw_rank": 5, "arithmetic_mw_rank": 4},
            },
        ],
        "edges": [
            {"source": "R_13trivial_plus_sign", "target": "R_13trivial"},
            {"source": "R_13trivial", "target": "R_12trivial"},
        ],
    }


def ns0024_gate(route):
    if route.get("status") != "PASS_EXACT_FRAME_AND_SOURCE_ROUTE_WITH_OPEN_EQUATION_TRANSPORT":
        raise ArithmeticError("NS0024 route certificate status changed")
    stages = route["completed_core_path"]["stages"]
    geometric = [row["mw_rank_for_rho_19"] for row in stages]
    if geometric != [4, 12, 12, 17]:
        raise ArithmeticError("NS0024 geometric rank path changed")
    boundary = route["proof_boundary"]
    missing = boundary["not_proved"]
    required_phrases = [
        "equation for the D5+E8 fibration",
        "elliptic-neighbour pencils along the Kneser path",
        "field of definition",
    ]
    # The field-of-definition omission is stated by the canonical theorem note,
    # while the artifact spells out the missing equation and marked route.  Keep
    # both gates explicit rather than inferring arithmetic descent from a J2 path.
    if required_phrases[0] not in missing or required_phrases[1] not in missing:
        raise ArithmeticError("NS0024 proof boundary no longer records the equation/route gap")
    return {
        "id": "NS0024_completed_core_path",
        "geometric_mw_ranks": geometric,
        "conditional_arithmetic_ranks_if_trivial_NS_action_and_U_descent": geometric,
        "proved_arithmetic_mw_ranks": [None, None, None, None],
        "promotion_gate": {
            "rational_source_divisor_span_rank_19": "MISSING",
            "K_defined_marked_U_at_each_stage": "MISSING",
            "galois_action_on_NS": "UNKNOWN",
            "component_actions": "UNKNOWN",
            "endpoint_arithmetic_rank_17": "NOT_PROVED",
        },
        "decision": "FAIL_CLOSED_GEOMETRIC_ONLY",
        "reason": (
            "The 17,13,7 path is a completed-core Kneser path and J2 existence "
            "certificate, not a marked elliptic-neighbour corridor or a field-of-definition certificate."
        ),
    }


def build_payload():
    paths = [
        MARKING_SCHEMA,
        R17_GRAM,
        H3_ENDPOINT,
        H3_PUBLISHED_TARGET,
        H3_PUBLISHED_COORDINATE_MATCH,
        H3_PUBLISHED_SECTIONS,
        E6_INCIDENCE,
        E6_ORBIT103,
        NS0024_ROUTE,
        R17_BISECTION_SPLITTING,
        R17_NORM12_CLASSIFICATION,
        ALTERNATE_Q80_FRAME,
    ]
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
    endpoint = json.loads(H3_ENDPOINT.read_text())
    published_target = json.loads(H3_PUBLISHED_TARGET.read_text())
    coordinate_match = json.loads(H3_PUBLISHED_COORDINATE_MATCH.read_text())
    published_sections = json.loads(H3_PUBLISHED_SECTIONS.read_text())
    incidence = json.loads(E6_INCIDENCE.read_text())
    orbit103 = json.loads(E6_ORBIT103.read_text())
    ns0024 = json.loads(NS0024_ROUTE.read_text())
    splitting = json.loads(R17_BISECTION_SPLITTING.read_text())
    classification = json.loads(R17_NORM12_CLASSIFICATION.read_text())
    alternate_payload = json.loads(ALTERNATE_Q80_FRAME.read_text())
    h3_input = h3_marking(endpoint, published_target, published_sections)
    e6_input = e6_incidence_marking(incidence)
    controls = [validate_marking(h3_input), validate_marking(e6_input)]
    transfer_regression = validate_marking(abstract_transfer_marking())
    alternate_q80 = alternate_q80_arithmetic_application(
        endpoint,
        published_target,
        published_sections,
        coordinate_match,
        splitting,
        classification,
        alternate_payload,
    )
    return {
        "schema": "elkies-k3.arithmetic-rank-transfer-controls.v1",
        "status": "PASS_ALTERNATE_Q80_ARITHMETIC_RANK17_AND_EQUIVARIANT_CONTROLS",
        "inputs": {relative(path): digest(path) for path in paths},
        "arithmetic_marking_schema": {
            "required": [
                "ground_field",
                "geometric NS Gram",
                "finite integral Galois generators",
                "marked U basis for every fibration",
                "embedded geometric root basis",
            ],
            "derived": [
                "finite Galois image",
                "rational NS fixed subspace",
                "root fixed subspaces",
                "geometric and arithmetic MW ranks",
                "MW character traces",
                "representation-ring edge identities",
            ],
            "optional_integral_and_field_layer": [
                "integral action on W/R including torsion",
                "saturated fixed Mordell--Weil lattice",
                "section orbits and stabilizer fields",
            ],
            "failure_policy": (
                "Missing target-U descent or field-of-definition evidence is UNKNOWN, never rank zero or arithmetic promotion."
            ),
        },
        "full_lattice_controls": controls,
        "representation_transfer_regression": {
            "scope": (
                "Abstract integral C2-lattice regression for Theorem A2.3; not a geometric existence claim."
            ),
            "result": transfer_regression,
        },
        "module_control": orbit103_module_control(orbit103),
        "alternate_q80_application": alternate_q80,
        "ns0024_application": ns0024_gate(ns0024),
        "theorem_checks": {
            "galois_actions_integral_and_gram_preserving": True,
            "marked_U_fixed_pointwise": True,
            "root_spaces_stable": True,
            "fixed_ranks_computed_over_QQ": True,
            "representation_ring_identity_checked_by_all_group_traces": True,
        },
        "proof_boundary": {
            "proved": (
                "The exact finite-module verifier reproduces arithmetic rank 17 for the rational H3 endpoint, "
                "promotes the norm12-orbit-11952 alternate-Q80 fibration to arithmetic rank 17 over QQ, "
                "reproduces rank 2 from geometric rank 4 for the unordered E6 incidence, and checks the "
                "2+chi_-3 orbit-103 split."
            ),
            "application": (
                "The new NS0024 completed-core route remains geometric-only and is rejected by the arithmetic promotion gate."
            ),
            "not_proved": (
                "No characteristic-zero NS0024 source equation, rational rank-19 divisor span, marked elliptic-neighbour "
                "corridor, Galois action, or arithmetic MW17 endpoint is constructed."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/certify_arithmetic_rank_transfer.sage --check"
        ),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
parser.add_argument(
    "--marking",
    type=Path,
    help="validate one self-contained arithmetic-marking-v1 JSON record",
)
arguments = parser.parse_args()
if arguments.marking:
    marking_result = validate_marking(json.loads(arguments.marking.read_text()))
    print(json.dumps(marking_result, indent=2, sort_keys=True))
    raise SystemExit(0)
payload = build_payload()
rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output = arguments.output.resolve()
if arguments.check:
    if not output.exists() or output.read_text() != rendered:
        raise SystemExit("arithmetic rank-transfer certificate is stale or missing")
else:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered)
print(
    "ARITHMETICRANKTRANSFER|H3=17|ALTERNATE_Q80=17|E6_INCIDENCE=2/4|"
    "E6_ORBIT103=2/3|NS0024=GEOMETRIC_ONLY|status=PASS",
    flush=True,
)
