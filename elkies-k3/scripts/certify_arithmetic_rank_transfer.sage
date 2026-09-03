#!/usr/bin/env sage-python
"""Certify the equivariant arithmetic rank-transfer controls.

status: ACTIVE_PROOF
claim: exact finite Galois-module rank transfer, H3/E6 controls, and a
       fail-closed NS0024 arithmetic-promotion gate
inputs: the pinned R17 Gram and existing H3, E6, and NS0024 certificates
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

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
MARKING_SCHEMA = ROOT / "elkies-k3/data/arithmetic/arithmetic-marking-v1.schema.json"
R17_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
H3_ENDPOINT = GEN / "elkies-k3-h3-q12o5867-endpoint-certificate.json"
E6_INCIDENCE = GEN / "elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
E6_ORBIT103 = (
    GEN / "elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json"
)
NS0024_ROUTE = GEN / "elkies-k3-ns0024-new-rootless-source-route-v1.json"
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
        fixed_root_rank = common_fixed_basis(root_actions, root_basis.ncols()).nrows()
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


def h3_marking(endpoint):
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
    gram = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -r17)
    identity = identity_matrix(ZZ, 19)
    fibre = [1] + [0] * 18
    mate = [0, 1] + [0] * 17
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
                "expected": {"geometric_mw_rank": 17, "arithmetic_mw_rank": 17},
            }
        ],
        "edges": [],
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
        E6_INCIDENCE,
        E6_ORBIT103,
        NS0024_ROUTE,
    ]
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
    endpoint = json.loads(H3_ENDPOINT.read_text())
    incidence = json.loads(E6_INCIDENCE.read_text())
    orbit103 = json.loads(E6_ORBIT103.read_text())
    ns0024 = json.loads(NS0024_ROUTE.read_text())
    h3_input = h3_marking(endpoint)
    e6_input = e6_incidence_marking(incidence)
    controls = [validate_marking(h3_input), validate_marking(e6_input)]
    transfer_regression = validate_marking(abstract_transfer_marking())
    return {
        "schema": "elkies-k3.arithmetic-rank-transfer-controls.v1",
        "status": "PASS_EQUIVARIANT_RANK_TRANSFER_CONTROLS_AND_FAIL_CLOSED_NS0024_GATE",
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
                "rank 2 from geometric rank 4 for the unordered E6 incidence, and the 2+chi_-3 orbit-103 split."
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
    "ARITHMETICRANKTRANSFER|H3=17|E6_INCIDENCE=2/4|"
    "E6_ORBIT103=2/3|NS0024=GEOMETRIC_ONLY|status=PASS",
    flush=True,
)
