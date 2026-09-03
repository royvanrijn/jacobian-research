#!/usr/bin/env sage-python
"""Replay the marked-U realization planner controls and end-to-end R17 route.

status: REGRESSION
claim: three positive marked-U controls, one fail-closed foundry control, and
       the existing R17 equation/rank-17 route through a target-free request
inputs: the pinned relative-U, R17, equation, arithmetic, and foundry artifacts
outputs: elkies-k3-marked-u-realization-planner-controls-v1.json

The four controls are intentionally ordered:

1. published R17 -> alternate Q80 at degree two;
2. published R17 -> the new non-cyclic 4A1/MW13 marking;
3. that 4A1/MW13 marking -> a cheap rootless successor;
4. only then inspect a non-historical foundry target, failing closed because
   its route ledger does not contain a common marked ``(NS,U,W)`` source.

The final replay asks only for root rank zero (and a genuinely different frame
class), selects ``norm12-orbit-11952`` within the declared equation cost bound,
then checks the existing exact equation and arithmetic-rank certificates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
PLANNER = ROOT / "elkies-k3/scripts/plan_marked_u_realizations.sage"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
PINNED_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
CLASSIFICATION = GENERATED / "elkies-k3-r17-norm12-isotropic-frame-classification-v1.json"
ALTERNATE = GENERATED / "q80-alternate-fifth-q6-rootless-transport.json"
LOCAL_MUTATION = GENERATED / "elkies-k3-r17-local-bridge-mutation-v1.json"
DIRECT_EQUATION = GENERATED / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
ARITHMETIC = GENERATED / "elkies-k3-arithmetic-rank-transfer-controls-v1.json"
PUBLISHED_EQUATION = GENERATED / "elkies-2026-published-r17-target.json"
FOUNDRY = GENERATED / "elkies-k3-lattice-foundry-v1.json"
FOUNDRY_READINESS = GENERATED / "elkies-k3-inverse-ade-foundry-readiness-v2.json"
OUTPUT = GENERATED / "elkies-k3-marked-u-realization-planner-controls-v1.json"
J = matrix(ZZ, [[0, 1], [1, 0]])


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def load_text_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def standard_marked_source(frame, label):
    ns = block_diagonal_matrix(J, -frame)
    basis = identity_matrix(ZZ, ns.nrows())
    return {
        "label": label,
        "ns_gram": rows(ns),
        "u_basis_in_ns": rows(basis[:2, :]),
        "frame_basis_in_ns": rows(basis[2:, :]),
    }


def marked_source(ns, source_u, source_w, label):
    return {
        "label": label,
        "ns_gram": rows(ns),
        "u_basis_in_ns": rows(source_u),
        "frame_basis_in_ns": rows(source_w),
    }


def artifact(path, schema, status):
    return {
        "path": relative(path),
        "sha256": digest(path),
        "schema": schema,
        "status": status,
    }


def classification_evidence(index, gate, classification_artifact):
    bindings = [
        {
            "pointer": f"/classification/records/{index}/isotropic_fibre",
            "candidate_field": "target_fibre_in_ambient_ns",
        }
    ]
    if gate == "effective_zero":
        bindings.append(
            {
                "pointer": f"/classification/records/{index}/shared_zero",
                "equals": True,
            }
        )
    return {
        "verdict": "PASS_EXACT",
        "artifact": classification_artifact,
        "bindings": bindings,
    }


def classification_catalog(classification, classification_artifact):
    candidates = []
    for index, record in enumerate(classification["classification"]["records"]):
        trace = list(map(int, record["trace_vector"]))
        candidates.append(
            {
                "id": record["label"],
                "intersection_coordinates": [2, 1, 1, -2],
                "projected_vectors": [trace, trace],
                "metadata": {
                    "frame_class_from_independent_replay": record["frame_class"],
                    "expected_riemann_roch_ambient": 10,
                    "equation_complexity": record["equation_complexity"],
                },
                "gate_evidence": {
                    "nefness": classification_evidence(
                        index, "nefness", classification_artifact
                    ),
                    "effective_zero": classification_evidence(
                        index, "effective_zero", classification_artifact
                    ),
                },
            }
        )
    return {
        "mode": "catalog",
        "scope": (
            "the 43 exact norm-twelve irreducible smooth genus-one curves in the "
            "published-R17 bisection certificate; not every representation of G_A"
        ),
        "complete_for_declared_box": False,
        "candidates": candidates,
    }


def compact_result(result, include_selected_u=True):
    selected = result["selected_realization"]
    compact = {
        "status": result["status"],
        "target_request": result["target_request"],
        "search_order": result["search_order"],
        "search_completeness": result["search_completeness"],
        "counters": result["counters"],
        "retained_candidate_ids": [
            row["candidate_id"] for row in result["retained_realizations"]
        ],
        "selected_candidate_id": selected["candidate_id"] if selected else None,
    }
    if selected:
        compact["selected"] = {
            "intersection_coordinates": selected["intersection_coordinates"],
            "cross_pairing_A": selected["cross_pairing_A"],
            "positive_projection_gram_G_A": selected[
                "positive_projection_gram_G_A"
            ],
            "root_ADE_gate": selected["root_ADE_gate"],
            "bridge": selected["bridge"],
            "physical_gates": selected["physical_gates"],
            "equation_facing_cost": selected["equation_facing_cost"],
            "child_frame_gram_sha256": selected["child_frame_gram_sha256"],
        }
        if include_selected_u:
            compact["selected"]["primitive_u_basis_in_ambient_ns"] = selected[
                "primitive_u_basis_in_ambient_ns"
            ]
            compact["selected"]["target_zero_in_ambient_ns"] = selected[
                "target_zero_in_ambient_ns"
            ]
    return compact


def exact_engine_control(planner):
    cross = planner["physical_cross_matrix"]((2, 1, 1, 0))
    bridge_gram = cross.transpose() * J * cross - J
    source = standard_marked_source(bridge_gram, "rank-two exact-shell unit control")
    configuration = {
        "source": source,
        "target": {},
        "intersection_box": {
            "F_dot_F_prime": [2],
            "F_dot_O_prime": [1],
            "O_dot_F_prime": [1],
            "O_dot_O_prime": [0],
        },
        "prime_local_bridge_constraints": {
            "relative_rank": 2,
            "saturation_index": 1,
            "saturated_bridge_determinant": 32,
            "saturated_bridge_gram": rows(bridge_gram),
            "det_A_parity": "odd",
        },
        "representations": {"mode": "exact"},
        "physical_gates": {"required": False},
    }
    result = planner["plan"](configuration)
    assert result["status"] == "PASS_MARKED_U_REALIZATIONS_FOUND"
    assert result["search_completeness"][
        "complete_over_all_integral_representations_in_declared_box"
    ]
    assert result["counters"]["representations_enumerated"] > 0
    rejected = planner["prime_local_screen"](
        cross,
        bridge_gram,
        {
            "relative_rank": 2,
            "saturation_index": 1,
            "saturated_bridge_determinant": 31,
        },
    )
    assert not rejected["pass"]
    assert "BRIDGE_SQUARE_INDEX_OBSTRUCTION" in rejected["rejection_reasons"]
    return {
        "status": "PASS_EXACT_SHELL_AND_PREFILTER_CONTROL",
        "representations_enumerated": result["counters"][
            "representations_enumerated"
        ],
        "retained_realizations": result["counters"]["retained_realizations"],
        "complete": True,
        "prime_local_negative_control": rejected,
        "representations_enumerated_after_negative_prefilter": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    planner = runpy.run_path(str(PLANNER))
    short = load_text_matrix(SHORT_GRAM)
    pinned = load_text_matrix(PINNED_GRAM)
    published_ns = block_diagonal_matrix(J, -short)
    published_u = identity_matrix(ZZ, 19)[:2, :]
    published_zero = published_u.row(1) - published_u.row(0)

    classification = load(CLASSIFICATION)
    assert classification["status"] == "PASS_EXACT_MINIMAL_J2_ACCESSIBILITY"
    classification_artifact = artifact(
        CLASSIFICATION,
        "elkies-k3.r17-norm12-isotropic-frame-classification.v1",
        "PASS_EXACT_MINIMAL_J2_ACCESSIBILITY",
    )
    catalog = classification_catalog(classification, classification_artifact)
    alternate_payload = load(ALTERNATE)
    alternate_frame = matrix(ZZ, alternate_payload["rootless_frame"])

    common_box = {
        "F_dot_F_prime": [2],
        "F_dot_O_prime": [1],
        "O_dot_F_prime": [1],
        "O_dot_O_prime": [-2],
    }
    control_1_input = {
        "source": standard_marked_source(pinned, "published-R17 equation marking"),
        "target": {
            "root_rank": 0,
            "ade_type": "rootless",
            "frame_gram": rows(alternate_frame),
        },
        "intersection_box": common_box,
        "prime_local_bridge_constraints": {"relative_rank": 1},
        "representations": catalog,
        "physical_gates": {"required": True},
    }
    control_1 = planner["plan"](control_1_input)
    assert control_1["status"] == "PASS_MARKED_U_REALIZATIONS_FOUND"
    assert len(control_1["retained_realizations"]) == 10
    assert control_1["selected_realization"]["candidate_id"] == "norm12-orbit-11952"
    assert control_1["selected_realization"]["intersection_coordinates"] == {
        "F_dot_F_prime": 2,
        "F_dot_O_prime": 1,
        "O_dot_F_prime": 1,
        "O_dot_O_prime": -2,
    }

    local = load(LOCAL_MUTATION)
    assert local["status"] == "PASS_EXACT_R17_LOCAL_BRIDGE_MUTATION"
    local_artifact = artifact(
        LOCAL_MUTATION,
        "elkies-k3.r17-local-bridge-mutation.v1",
        "PASS_EXACT_R17_LOCAL_BRIDGE_MUTATION",
    )
    new_u = matrix(ZZ, local["r17_example"]["relative_U"]["target_U_basis_in_U_plus_short_R17"])
    new_w = matrix(ZZ, local["r17_example"]["target_frame"]["basis_in_ambient_NS"])
    local_projected = rows(new_u[:, 2:])
    local_candidate = {
        "id": "published-R17-degree2-new-noncyclic-4A1",
        "intersection_coordinates": [2, 1, 1, 0],
        "projected_vectors": local_projected,
        "metadata": {
            "expected_riemann_roch_ambient": 2,
            "equation_complexity": {},
        },
        "gate_evidence": {
            "nefness": {
                "verdict": "PASS_EXACT",
                "artifact": local_artifact,
                "bindings": [
                    {
                        "pointer": "/r17_example/relative_U/target_U_basis_in_U_plus_short_R17",
                        "candidate_field": "primitive_u_basis_in_ambient_ns",
                    },
                    {
                        "pointer": "/r17_example/geometric_gate/new_fibre_nef_in_old_chamber",
                        "equals": True,
                    },
                ],
            },
            "effective_zero": {
                "verdict": "PASS_EXACT",
                "artifact": local_artifact,
                "bindings": [
                    {
                        "pointer": "/r17_example/relative_U/target_U_basis_in_U_plus_short_R17",
                        "candidate_field": "primitive_u_basis_in_ambient_ns",
                    },
                    {
                        "pointer": "/r17_example/geometric_gate/new_zero_is_physical",
                        "equals": True,
                    },
                ],
            },
        },
    }
    noncyclic_constraints = {
        "relative_rank": 2,
        "raw_bridge_determinant": 32,
        "saturation_index": 1,
        "saturated_bridge_determinant": 32,
        "saturated_bridge_gram": [[4, 0], [0, 8]],
        "p_primary_invariants": {"2": [4, 8]},
        "det_A_parity": "odd",
    }
    control_2 = planner["plan"](
        {
            "source": standard_marked_source(short, "published-R17 short-vector marking"),
            "target": {"root_rank": 4, "ade_type": "4A1"},
            "intersection_box": {
                "F_dot_F_prime": [2],
                "F_dot_O_prime": [1],
                "O_dot_F_prime": [1],
                "O_dot_O_prime": [0],
            },
            "prime_local_bridge_constraints": noncyclic_constraints,
            "representations": {
                "mode": "catalog",
                "scope": "single exact H-1c R17 local-mutation witness",
                "complete_for_declared_box": False,
                "candidates": [local_candidate],
            },
            "physical_gates": {"required": True},
        }
    )
    assert control_2["status"] == "PASS_MARKED_U_REALIZATIONS_FOUND"
    selected_2 = control_2["selected_realization"]
    assert selected_2["primitive_u_basis_in_ambient_ns"] == rows(new_u)
    assert selected_2["root_ADE_gate"]["ade_type"] == "4A1"
    assert selected_2["root_ADE_gate"]["mordell_weil_rank_from_shioda_tate"] == 13
    assert selected_2["bridge"]["saturated_discriminant_invariants"] == [4, 8]

    new_transport = new_u.stack(new_w)
    assert abs(int(new_transport.det())) == 1
    assert new_transport * published_ns * new_transport.transpose() == block_diagonal_matrix(
        J, -matrix(ZZ, local["r17_example"]["target_frame"]["gram"])
    )
    old_in_new = published_u * new_transport.inverse()
    reverse_cross = new_u * published_ns * published_u.transpose()
    reverse_coordinates = planner["physical_coordinates"](reverse_cross)
    assert reverse_coordinates == (2, 1, 1, 0)
    assert old_in_new[:, :2] == reverse_cross.transpose() * J
    reverse_projected = old_in_new[:, 2:]

    published = load(PUBLISHED_EQUATION)
    assert published["status"] == "PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17"
    assert published["pinned_identification"]["integrally_isometric"]
    assert published["published_height_lattice"]["rank"] == 17
    assert published["published_height_lattice"]["determinant"] == 948
    published_artifact = artifact(
        PUBLISHED_EQUATION,
        "elkies-k3.elkies-2026-published-r17-target.v1",
        "PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17",
    )
    reverse_candidate = {
        "id": "new-4A1-cheap-return-to-published-rootless",
        "intersection_coordinates": list(reverse_coordinates),
        "projected_vectors": rows(reverse_projected),
        "metadata": {
            "expected_riemann_roch_ambient": 2,
            "equation_complexity": {},
        },
        "gate_evidence": {
            "nefness": {
                "verdict": "PASS_EXACT",
                "artifact": published_artifact,
                "bindings": [
                    {
                        "pointer": "/published_equation/rootless_semistable",
                        "equals": True,
                    }
                ],
                "literal_bindings": [
                    {
                        "candidate_field": "primitive_u_basis_in_ambient_ns",
                        "equals": rows(published_u),
                    }
                ],
            },
            "effective_zero": {
                "verdict": "PASS_EXACT",
                "artifact": published_artifact,
                "bindings": [
                    {
                        "pointer": "/published_equation/published_section_identities",
                        "equals": 17,
                    }
                ],
                "literal_bindings": [
                    {
                        "candidate_field": "target_zero_in_ambient_ns",
                        "equals": list(map(int, published_zero)),
                    }
                ],
            },
        },
    }
    new_frame = matrix(ZZ, local["r17_example"]["target_frame"]["gram"])
    control_3 = planner["plan"](
        {
            "source": marked_source(
                published_ns,
                new_u,
                new_w,
                "new non-cyclic 4A1/MW13 R17 marking",
            ),
            "target": {
                "root_rank": 0,
                "ade_type": "rootless",
            },
            "intersection_box": {
                "F_dot_F_prime": [2],
                "F_dot_O_prime": [1],
                "O_dot_F_prime": [1],
                "O_dot_O_prime": [0],
            },
            "prime_local_bridge_constraints": noncyclic_constraints,
            "representations": {
                "mode": "catalog",
                "scope": "reverse of the exact H-1c local-mutation witness",
                "complete_for_declared_box": False,
                "candidates": [reverse_candidate],
            },
            "physical_gates": {"required": True},
        }
    )
    assert control_3["status"] == "PASS_MARKED_U_REALIZATIONS_FOUND"
    assert control_3["target_request"]["target_frame_supplied"] is False
    assert control_3["selected_realization"]["primitive_u_basis_in_ambient_ns"] == rows(
        published_u
    )
    assert control_3["selected_realization"]["root_ADE_gate"]["ade_type"] == "rootless"
    assert control_3["selected_realization"]["root_ADE_gate"][
        "mordell_weil_rank_from_shioda_tate"
    ] == 17

    # The non-historical foundry target is deliberately inspected only after
    # the three marked positive controls.  Its bulk ledger row lacks the data
    # needed to form an explicit source, so no neighbour or representation
    # enumeration is launched.
    foundry = load(FOUNDRY)
    readiness = load(FOUNDRY_READINESS)
    target_row = next(
        row for row in foundry["rootless_targets"] if not row["is_existing_H3_control"]
    )
    assert target_row["frame_id"] == "NS0002-F007"
    readiness_row = next(
        row for row in readiness["routes"] if row["target_frame_id"] == target_row["frame_id"]
    )
    assert not readiness_row["planner_ready"]
    assert "common_marked_NS_basis_and_source_U" in readiness_row["missing_inputs"]
    foundry_control = {
        "status": "BLOCKED_MISSING_EXPLICIT_MARKED_SOURCE",
        "attempt_order": 4,
        "target_frame_id": target_row["frame_id"],
        "target_ns_id": target_row["ns_id"],
        "historical_corridor_target": target_row["is_existing_H3_control"],
        "source_frame_id": readiness_row["source_frame_id"],
        "missing_inputs": readiness_row["missing_inputs"],
        "representation_enumeration_started": False,
        "compute_campaign_started": False,
        "reason": (
            "A positive target Gram and a core route do not supply a common marked "
            "NS, source U, frame basis, or literal relative-U transport."
        ),
    }

    # End-to-end milestone: no target Gram is present.  The cost bound retains
    # the first distinct rootless frame in the certified bisection catalog.
    milestone_bound = {
        "maximum_old_fibre_degree": 2,
        "maximum_old_zero_degree": 1,
        "maximum_new_zero_old_fibre_degree": 1,
        "maximum_coefficient_l1": 13,
        "maximum_coordinate_input_bits": 3482,
    }
    milestone = planner["plan"](
        {
            "source": standard_marked_source(pinned, "published-R17 equation marking"),
            "target": {
                "root_rank": 0,
                "ade_type": "rootless",
                "exclude_source_frame_isometry": True,
            },
            "intersection_box": common_box,
            "prime_local_bridge_constraints": {"relative_rank": 1},
            "representations": catalog,
            "physical_gates": {"required": True},
            "equation_cost_bound": milestone_bound,
        }
    )
    selected_milestone = milestone["selected_realization"]
    assert milestone["target_request"]["target_frame_supplied"] is False
    assert selected_milestone["candidate_id"] == "norm12-orbit-11952"
    assert selected_milestone["root_ADE_gate"][
        "mordell_weil_rank_from_shioda_tate"
    ] == 17

    direct = load(DIRECT_EQUATION)
    arithmetic = load(ARITHMETIC)
    assert direct["status"] == "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS"
    assert direct["divisor"]["label"] == selected_milestone["candidate_id"]
    assert direct["divisor"]["class_D_in_U_plus_R17_minus"] == selected_milestone[
        "target_fibre_in_ambient_ns"
    ]
    assert direct["weierstrass_model"]["equation"] == "Y^2=X^3+A(u)*X+B(u)"
    assert direct["weierstrass_model"]["fibre_configuration"] == "24 I1"
    assert direct["sections"]["status"] == "PASS_EXACT_SATURATED_RANK17_BASIS"
    assert direct["sections"]["rank"] == 17
    arithmetic_application = arithmetic["alternate_q80_application"]
    assert arithmetic_application["witness"] == selected_milestone["candidate_id"]
    assert arithmetic_application["status"] == (
        "PASS_EXACT_ALTERNATE_Q80_ARITHMETIC_RANK_17_BEFORE_EQUATION_COMPILATION"
    )
    assert arithmetic_application["conclusion"]["arithmetic_mordell_weil_rank"] == 17

    engine_control = exact_engine_control(planner)
    payload = {
        "schema": "elkies-k3.marked-u-realization-planner-controls.v1",
        "status": "PASS_MARKED_U_PLANNER_CONTROLS_AND_R17_END_TO_END",
        "planner_separation": {
            "marked_u_planner": relative(PLANNER),
            "core_planner": "elkies-k3/scripts/plan_inverse_ade_targets.sage",
            "same_graph": False,
            "interface": (
                "A core-class proposal may supply a target request, but acceptance "
                "requires a literal primitive nef U' in one explicit marked NS."
            ),
        },
        "exact_engine_control": engine_control,
        "ordered_controls": [
            {
                "order": 1,
                "name": "published-R17 to alternate-Q80 degree-two copy",
                "result": compact_result(control_1),
            },
            {
                "order": 2,
                "name": "published-R17 to new non-cyclic 4A1/MW13",
                "result": compact_result(control_2),
            },
            {
                "order": 3,
                "name": "new non-cyclic 4A1/MW13 to cheap rootless successor",
                "result": compact_result(control_3),
            },
            {
                "order": 4,
                "name": "first non-historical foundry target readiness attempt",
                "result": foundry_control,
            },
        ],
        "major_milestone_replay": {
            "status": "PASS_EXISTING_R17_ROUTE_THROUGH_NEW_PLANNER",
            "instruction": "find a distinct rootless rank-17 fibration within cost B",
            "target_gram_supplied": False,
            "characteristic_zero_source": {
                "artifact": relative(PUBLISHED_EQUATION),
                "sha256": digest(PUBLISHED_EQUATION),
                "status": published["status"],
                "equation": published["published_equation"]["form"],
                "published_section_identities": published["published_equation"][
                    "published_section_identities"
                ],
                "integrally_identified_with_planner_source_frame": published[
                    "pinned_identification"
                ]["integrally_isometric"],
            },
            "cost_bound_B": milestone_bound,
            "planner_result": compact_result(milestone),
            "equation_compilation": {
                "artifact": relative(DIRECT_EQUATION),
                "sha256": digest(DIRECT_EQUATION),
                "status": direct["status"],
                "equation": direct["weierstrass_model"]["equation"],
                "fibre_configuration": direct["weierstrass_model"][
                    "fibre_configuration"
                ],
                "saturated_section_rank": direct["sections"]["rank"],
            },
            "arithmetic_rank_certificate": {
                "artifact": relative(ARITHMETIC),
                "sha256": digest(ARITHMETIC),
                "status": arithmetic_application["status"],
                "rank": arithmetic_application["conclusion"][
                    "arithmetic_mordell_weil_rank"
                ],
                "field": arithmetic_application["conclusion"]["field"],
            },
            "boundary": (
                "This replays an already-certified R17 equation route through the new "
                "target-free marked-U planner. It is not a claim that equation compilation "
                "has been automated for arbitrary planner hits or unready foundry targets."
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                PLANNER,
                SHORT_GRAM,
                PINNED_GRAM,
                CLASSIFICATION,
                ALTERNATE,
                LOCAL_MUTATION,
                DIRECT_EQUATION,
                ARITHMETIC,
                PUBLISHED_EQUATION,
                FOUNDRY,
                FOUNDRY_READINESS,
            )
        },
        "proof_boundary": {
            "proved": (
                "Exact marked-source validation, lexicographic cross-matrix enumeration, "
                "pre-representation prime-local screens, literal primitive-U construction, "
                "root/ADE and target-frame gates, exact evidence-bound nef/zero gates, and "
                "the four ordered controls above."
            ),
            "not_proved": (
                "Global optimality outside declared bounds/catalogs, an automatic nefness "
                "theorem without evidence, generic equation compilation, or readiness of "
                "the 936 bulk foundry routes."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_marked_u_realization_planner_controls.sage --check"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "MARKEDU|controls=3-pass+1-fail-closed|milestone=R17-rootless-rank17|"
        f"status={payload['status']}|output={relative(output)}"
    )


if __name__ == "__main__":
    main()
