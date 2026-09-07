#!/usr/bin/env sage-python
"""Certify the single curated equation-to-planner-to-MW17 foundry route.

The planner sees the marked parent core, bridge/glue state, good prime, root
survival requirement, desired root metric, and low-shell overlap fingerprint.
It does not see the stored neighbour line, child Gram matrix, target foundry
frame, target U, or physical elliptic-neighbour path.  Those objects enter only
after a hit, as independent endpoint/transport checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
ROUTE = ROOT / "elkies-k3/data/lattice-foundry/planner-ready-h3-a1-r17-v1.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-single-planner-ready-foundry-route-v1.json"
)
PLANNER_SCRIPT = ROOT / "elkies-k3/scripts/plan_inverse_ade_targets.sage"
SEARCH_SCRIPT = (
    ROOT / "elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage"
)

REQUIRED_FIELDS = [
    "common_marked_NS_basis_and_source_U",
    "rank15_core_plus_rank2_bridge_decomposition",
    "graph_glue_generator",
    "good_neighbor_prime_plan",
    "prescribed_survival_birth_templates",
    "elliptic_neighbor_transport_or_relative_U_lift",
]
FORBIDDEN_PLANNER_FIELDS = {
    "selected_isotropic_line",
    "adjusted_isotropic_lift",
    "child_core_gram",
    "target_frame_gram",
    "target_frame_id",
    "target_U_basis",
    "elliptic_neighbor_path",
}


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_matrix(value):
    return [list(map(int, row)) for row in value.rows()]


def rational_vector(value):
    return [str(entry) for entry in value]


def strip_timings(value):
    if isinstance(value, dict):
        return {
            key: strip_timings(item)
            for key, item in value.items()
            if key != "elapsed_seconds"
        }
    if isinstance(value, list):
        return [strip_timings(item) for item in value]
    return value


def all_keys(value):
    if isinstance(value, dict):
        answer = set(value)
        for item in value.values():
            answer.update(all_keys(item))
        return answer
    if isinstance(value, list):
        answer = set()
        for item in value:
            answer.update(all_keys(item))
        return answer
    return set()


def assert_hash(record):
    path = ROOT / record["path"]
    assert digest(path) == record["sha256"], path
    return json.loads(path.read_text())


def exact_shell(gram, shell_norm):
    enumeration = pari(gram).qfminim(shell_norm)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    return [value for value in representatives if value * gram * value == shell_norm]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    manifest = json.loads(ROUTE.read_text())
    assert manifest["route_count"] == 1
    assert len(manifest["routes"]) == 1
    route = manifest["routes"][0]
    assert [key for key in REQUIRED_FIELDS if key not in route] == []
    assert sum(key in route for key in REQUIRED_FIELDS) == 6

    source = route["common_marked_NS_basis_and_source_U"]
    equation = assert_hash(source["equation_source"])
    pointing = assert_hash(source["pointing_source"])
    assert equation["schema"] == source["equation_source"]["schema"]
    assert equation["status"] == source["equation_source"]["status"]
    assert pointing["status"] == source["pointing_source"]["status"]
    assert equation["child"]["ADE"] == "A1"
    assert equation["child"]["MW_rank_if_rho19"] == 16
    assert len(equation["child"]["minimal_A_coefficients_low_to_high"]) == 9
    assert len(equation["child"]["minimal_B_coefficients_low_to_high"]) == 13
    assert pointing["marking"]["unique_A1_root_bound_to_equation"]
    assert pointing["marking"]["prescribed_zero_pointed"]

    planner = runpy.run_path(str(PLANNER_SCRIPT))
    search = runpy.run_path(str(SEARCH_SCRIPT))
    base = runpy.run_path(str(search["BASE_SCRIPT"]))
    core_tools = runpy.run_path(str(search["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(search["REVERSE_SCRIPT"]))
    inverse = runpy.run_path(str(planner["INVERSE_SCRIPT"]))
    birth_death = runpy.run_path(str(inverse["BIRTH_DEATH_SCRIPT"]))
    signatures = runpy.run_path(str(planner["SIGNATURE_SCRIPT"]))
    planner["graph_state"].__globals__["search"] = search
    planner["completion_signature"].__globals__["search"] = search

    bridge_artifact_path = ROOT / source["source_frame_reconstruction"]["artifact"]
    bridge_artifact = json.loads(bridge_artifact_path.read_text())
    reconstruction = source["source_frame_reconstruction"]
    physical_edge = next(
        row
        for row in bridge_artifact["edges"]
        if row["corridor"] == reconstruction["corridor"]
        and int(row["edge_index"]) == reconstruction["edge_index"]
    )
    assert physical_edge["source_root_rank"] == 1
    assert physical_edge["target_root_rank"] == 0

    original_core = matrix(ZZ, physical_edge["core"]["gram"])
    original_bridge = matrix(
        ZZ, physical_edge[reconstruction["side"]]["bridge_gram"]
    )
    original_glue_row = physical_edge[reconstruction["side"]]["glue_generators"][0]
    original_glue = vector(
        QQ, [QQ(value) for value in original_glue_row["K_plus_C_dual_coordinates"]]
    )
    physical_source_frame = core_tools["glued_frame"](
        original_core, original_bridge, original_glue
    )
    source_transform = matrix(ZZ, pari(physical_source_frame).qflllgram()).transpose()
    source_frame = source_transform * physical_source_frame * source_transform.transpose()
    assert abs(int(source_frame.det())) == reconstruction["expected_frame_determinant"]

    decomposition = route["rank15_core_plus_rank2_bridge_decomposition"]
    bridge_basis = matrix(ZZ, decomposition["bridge_basis_rows_in_source_frame"])
    bridge_gram = bridge_basis * source_frame * bridge_basis.transpose()
    assert json_matrix(bridge_gram) == decomposition["bridge_gram"]
    assert abs(int(bridge_gram.det())) == decomposition["bridge_determinant"]
    core_basis = (bridge_basis * source_frame).right_kernel_matrix()
    core_gram_unreduced = core_basis * source_frame * core_basis.transpose()
    core_transform = matrix(ZZ, pari(core_gram_unreduced).qflllgram()).transpose()
    core_gram = core_transform * core_gram_unreduced * core_transform.transpose()
    assert json_matrix(core_gram) == decomposition["core_gram"]
    assert core_gram.nrows() == decomposition["core_rank"]
    assert abs(int(core_gram.det())) == decomposition["core_determinant"]
    assert int(pari(core_gram).qfminim(2)[0]) == 0
    assert int(pari(core_gram).qfminim(4)[0]) > 0

    theta_path = ROOT / "artifacts/generated-results/elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
    theta = json.loads(theta_path.read_text())
    h3_theta = next(row for row in theta["corridors"] if row["corridor"] == "H3")
    graph = route["graph_glue_generator"]
    bridge = next(
        row
        for row in base["bridge_data"](h3_theta, reverse)
        if row["bridge_class_index"] == graph["bridge_class_index"]
    )
    assert json_matrix(bridge["gram"]) == decomposition["bridge_gram"]
    assert rational_vector(bridge["generator"]) == graph["bridge_discriminant_generator"]
    state = planner["graph_state"](
        core_gram,
        bridge,
        graph["cyclic_order"],
        base,
        reverse,
        inverse,
    )
    assert rational_vector(state["generator"]) == graph["core_discriminant_generator"]
    assert state["multiplier"] == graph["core_multiplier"]
    source_signature = planner["completion_signature"](
        core_gram,
        bridge,
        graph["cyclic_order"],
        base,
        core_tools,
        reverse,
        signatures,
    )
    expected_source = graph["expected_source_completion"]
    for key in ("ade_type", "root_rank", "root_line_count"):
        assert source_signature[key] == expected_source[key]
    source_glue = vector(
        QQ,
        list(state["multiplier"] * state["generator"])
        + list(bridge["generator"]),
    )
    reconstructed_source_frame = core_tools["glued_frame"](
        core_gram, bridge["gram"], source_glue
    )
    assert abs(int(reconstructed_source_frame.det())) == expected_source["determinant"]
    assert pari(reconstructed_source_frame).qfisom(pari(source_frame)) != 0

    prime_plan = route["good_neighbor_prime_plan"]
    assert len(prime_plan["ordered_primes"]) == 1
    prime = ZZ(prime_plan["ordered_primes"][0])
    assert prime.is_prime() and prime % 2 and core_gram.det() % prime
    prescribed = route["prescribed_survival_birth_templates"]
    template = {
        "prime": int(prime),
        "surviving_parent_root_line_indices": prescribed[
            "surviving_parent_root_line_indices"
        ],
        "parent_root_lines": len(state["witnesses"]),
        "desired_child_signature": prescribed["desired_child_signature"],
        "desired_child_signature_key": ["rootless", 0, 0, []],
        "target_core_overlap_fingerprint": prescribed[
            "target_core_overlap_fingerprint"
        ],
    }
    public_state = {
        "core_discriminant_generator": rational_vector(state["generator"]),
        "core_multiplier": state["multiplier"],
        "parent_graph_root_witnesses": [
            {
                "core": rational_vector(row["core"]),
                "bridge": rational_vector(row["bridge"]),
                "graph_label": row["graph_label"],
            }
            for row in state["witnesses"]
        ],
    }
    planner_input = {
        "core_gram": json_matrix(core_gram),
        "bridge": {
            "class_index": bridge["bridge_class_index"],
            "gram": json_matrix(bridge["gram"]),
            "discriminant_generator": rational_vector(bridge["generator"]),
            "cyclic_order": graph["cyclic_order"],
        },
        "graph_state": public_state,
        "prime_and_template": template,
        "bounds": {
            "maximum_parameter_support": prime_plan["maximum_parameter_support"],
            "dense_probes": prime_plan["dense_probes"],
            "maximum_materialized_neighbors": prime_plan[
                "maximum_materialized_neighbors"
            ],
        },
    }
    assert not (all_keys(planner_input) & FORBIDDEN_PLANNER_FIELDS)
    planner_input_encoding = json.dumps(
        planner_input, sort_keys=True, separators=(",", ":")
    ).encode()
    planner_result = planner["plan_template"](
        core_gram,
        bridge,
        graph["cyclic_order"],
        state,
        template,
        None,
        base,
        core_tools,
        reverse,
        inverse,
        birth_death,
        signatures,
        prime_plan["maximum_parameter_support"],
        prime_plan["dense_probes"],
        False,
        prime_plan["maximum_materialized_neighbors"],
    )
    assert planner_result["status"] == "HIT"
    assert planner_result["statistics"]["materialized_neighbors"] == 1

    # The endpoint enters only here, after the planner has returned a line.
    discovered_line = vector(ZZ, planner_result["selected_isotropic_line"])
    child = base["quadratic_form"](core_gram).find_p_neighbor_from_vec(
        prime, discovered_line
    ).Hessian_matrix()
    child_core = base["lll_reduce"](child)
    child_signature = planner["completion_signature"](
        child_core,
        bridge,
        graph["cyclic_order"],
        base,
        core_tools,
        reverse,
        signatures,
    )
    assert child_signature["ade_type"] == "rootless"
    assert child_signature["root_rank"] == 0
    child_state = planner["graph_state"](
        child_core,
        bridge,
        graph["cyclic_order"],
        base,
        reverse,
        inverse,
    )
    child_glue = vector(
        QQ,
        list(child_state["multiplier"] * child_state["generator"])
        + list(bridge["generator"]),
    )
    child_frame = core_tools["glued_frame"](
        child_core, bridge["gram"], child_glue
    )
    assert abs(int(child_frame.det())) == 948

    transport = route["elliptic_neighbor_transport_or_relative_U_lift"]
    target_record = transport["target_identification"]
    foundry_path = ROOT / target_record["foundry_artifact"]
    assert digest(foundry_path) == target_record["foundry_sha256"]
    foundry = json.loads(foundry_path.read_text())
    target_frames = [
        frame
        for ns_class in foundry["ns_classes"]
        for frame in ns_class["frames"]
        if frame["frame_id"] in ("NS0001-F001", "NS0001-F002")
    ]
    target = next(
        frame for frame in target_frames if frame["frame_id"] == target_record["frame_id"]
    )
    alternate = next(
        frame for frame in target_frames if frame["frame_id"] == "NS0001-F002"
    )
    target_gram = matrix(ZZ, target["reduced_gram"])
    alternate_gram = matrix(ZZ, alternate["reduced_gram"])
    assert pari(child_frame).qfisom(pari(target_gram)) != 0
    assert pari(child_frame).qfisom(pari(alternate_gram)) == 0
    signed_norm_four = int(pari(child_frame).qfminim(4)[0])
    assert signed_norm_four == 2622

    observed_overlap = {}
    for shell_norm_text in prescribed["target_core_overlap_fingerprint"][
        "parent_shell_line_survivors"
    ]:
        shell_norm = ZZ(shell_norm_text)
        observed_overlap[shell_norm_text] = sum(
            1
            for value in exact_shell(core_gram, shell_norm)
            if (value * core_gram).dot_product(discovered_line) % prime == 0
        )
    assert observed_overlap == prescribed["target_core_overlap_fingerprint"][
        "parent_shell_line_survivors"
    ]

    relative_u_record = assert_hash(transport["relative_U_artifact"])
    relative_u_edge = next(
        row
        for row in relative_u_record["edges"]
        if row["corridor"] == transport["relative_U_artifact"]["corridor"]
        and row["edge_index"] == transport["relative_U_artifact"]["edge_index"]
    )
    physical = transport["physical_edge"]
    for key in ("q", "old_fibre_degree", "source_root_rank", "target_root_rank"):
        assert relative_u_edge[key] == physical[key]
    direction = transport["relative_U_artifact"]["direction"]
    relative_u_lift = relative_u_edge["relative_u"][direction]
    assert relative_u_lift["saturation_index"] == transport[
        "relative_U_artifact"
    ]["saturation_index"]
    assert relative_u_edge["relative_u"]["old_fibre_degree_matches_edge_record"]
    endpoint_path = ROOT / target_record["endpoint_certificate"]
    assert digest(endpoint_path) == target_record["endpoint_certificate_sha256"]

    payload = {
        "schema": "elkies-k3.single-planner-ready-foundry-route-certificate.v1",
        "status": "PASS_ONE_BLIND_PLANNER_READY_ROUTE_TO_MW17",
        "route_count": 1,
        "route_id": route["route_id"],
        "six_required_fields_present": REQUIRED_FIELDS,
        "source_equation": {
            "field": "QQ",
            "ade_type": equation["child"]["ADE"],
            "mw_rank_if_rho_19": equation["child"]["MW_rank_if_rho19"],
            "weierstrass_degrees": equation["child"]["degrees_A_B_Delta"],
            "exact_coefficients_checked": True,
            "marked_source_U_checked": True,
        },
        "planner_protocol": {
            "planner_input": planner_input,
            "planner_input_sha256": hashlib.sha256(planner_input_encoding).hexdigest(),
            "forbidden_endpoint_fields_absent": sorted(FORBIDDEN_PLANNER_FIELDS),
            "target_core_isometry_requested": False,
            "result": planner_result,
        },
        "discovered_endpoint": {
            "selected_isotropic_line": list(map(int, discovered_line)),
            "child_core_determinant": abs(int(child_core.det())),
            "completed_frame_determinant": abs(int(child_frame.det())),
            "completed_frame_ade_type": child_signature["ade_type"],
            "completed_frame_signed_norm_four_count": signed_norm_four,
            "exact_integral_isometry_to_NS0001_F001": True,
            "exact_integral_isometry_to_NS0001_F002": False,
            "observed_parent_target_shell_overlap": observed_overlap,
        },
        "relative_U_and_fibration_hop": {
            "corridor": relative_u_edge["corridor"],
            "edge_index": relative_u_edge["edge_index"],
            "q": relative_u_edge["q"],
            "orbit": physical["orbit"],
            "old_fibre_degree": relative_u_edge["old_fibre_degree"],
            "direction": direction,
            "intersection_coordinates": relative_u_lift["intersection_coordinates"],
            "saturated_bridge_gram": relative_u_lift["saturated_bridge_gram"],
            "saturation_index": relative_u_lift["saturation_index"],
            "target_root_rank": relative_u_edge["target_root_rank"],
            "target_mw_rank_if_rho_19": 17,
        },
        "chain": [
            "exact QQ A1/MW16 equation source",
            "blind inverse-ADE core planner",
            "rootless NS0001-F001 frame",
            "exact relative-U lift U'",
            "degree-2 elliptic-neighbour hop q6 orbit2247",
            "MW17",
        ],
        "inputs": {
            relative(ROUTE): digest(ROUTE),
            relative(PLANNER_SCRIPT): digest(PLANNER_SCRIPT),
            relative(bridge_artifact_path): digest(bridge_artifact_path),
            relative(theta_path): digest(theta_path),
            source["equation_source"]["path"]: source["equation_source"]["sha256"],
            source["pointing_source"]["path"]: source["pointing_source"]["sha256"],
            target_record["foundry_artifact"]: target_record["foundry_sha256"],
            transport["relative_U_artifact"]["path"]: transport[
                "relative_U_artifact"
            ]["sha256"],
            target_record["endpoint_certificate"]: target_record[
                "endpoint_certificate_sha256"
            ],
        },
        "proof_boundary": (
            "This certifies one exact lattice-planner discovery and identifies its "
            "rootless completion with the published R17 frame, followed by the "
            "already-certified exact relative-U/elliptic-neighbour transport. It "
            "does not make any of the 936 bulk foundry pairs planner-ready."
        ),
        "reproduce": (
            "sage -python "
            "elkies-k3/scripts/certify_single_planner_ready_foundry_route.sage"
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        old = json.loads(output.read_text())
        if strip_timings(old) != strip_timings(payload):
            raise SystemExit(f"stale artifact: {output}")
        print(payload["status"])
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(relative(output))


if __name__ == "__main__":
    main()
