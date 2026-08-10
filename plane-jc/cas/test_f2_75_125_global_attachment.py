#!/usr/bin/env python3
"""Regression tests for the F2 ``(75,125)`` attachment compiler."""

from copy import deepcopy

from verify_f2_75_125_global_attachment import (
    audit_candidate,
    build_payload,
    minimal_source_principal_boundary,
    topology_bounds,
)


def synthetic_squarefree_candidate() -> dict[str, object]:
    """Return an abstract gate witness, not a geometric F2 completion."""

    terminal_cycles = [
        [2, 3, 4, 5, 1, 6],
        [3, 1, 2, 6, 4, 5],
        [6, 2, 3, 1, 5, 4],
    ]
    transposition = [2, 1, 3, 4, 5, 6]
    source_graph = minimal_source_principal_boundary("squarefree_one_packet")
    terminal = source_graph["terminal_components"][0]
    components = []
    for name, self_intersection in zip(
        source_graph["component_order"],
        source_graph["self_intersections_in_component_order"],
    ):
        component = {
            "name": name,
            "self_intersection": self_intersection,
            "genus": 0,
        }
        if name == terminal:
            component.update(
                original_valuation_orders={
                    "x": -25,
                    "y_old": 5,
                    "x*y_old^5-center": 12,
                },
                carrier_center="1",
                translated_kummer_identification_certified=True,
                proximity_chain_certified=True,
            )
        components.append(component)
    return {
        "case": "squarefree_one_packet",
        "geometric_degree": 6,
        "source_boundary": {
            "complete": True,
            "open_surface": "A2",
            "smooth_completion_certified": True,
            "snc_certified": True,
            "geometric_components_certified": True,
            "components": components,
            "edges": source_graph["edges"],
            "terminal_components": source_graph["terminal_components"],
            "attachment_neighbors": source_graph["attachment_neighbors"],
        },
        "target_ledgers": [
            {
                "name": "T_terminal",
                "center": "infinity",
                "boundary_rows": [
                    {
                        "name": "E",
                        "e": 1,
                        "f": 6,
                        "principal_packet": "packet_1",
                        "source_component": terminal,
                    }
                ],
                "affine_rows": [],
                "exhaustive": True,
                "finite_flat_certified": True,
                "target_transfer_certified": True,
            },
            {
                "name": "C_affine",
                "center": "affine_nonproperness",
                "boundary_rows": [{"name": "R", "e": 2, "f": 1}],
                "affine_rows": [{"name": "A", "f": 4}],
                "exhaustive": True,
                "finite_flat_certified": True,
                "target_transfer_certified": True,
            },
        ],
        "spectator_orbits": [
            {
                "name": "S1",
                "target_ledger": "T_terminal",
                "target_branch_value": "b1",
                "transverse_index": 1,
                "residue_degree": 1,
                "inertia_permutation": transposition,
                "geometry_certified": True,
            },
            {
                "name": "S2",
                "target_ledger": "T_terminal",
                "target_branch_value": "b2",
                "transverse_index": 1,
                "residue_degree": 1,
                "inertia_permutation": transposition,
                "geometry_certified": True,
            },
        ],
        "global_meridian_systems": [
            {
                "name": "T_terminal",
                "connected": True,
                "cycles": [*terminal_cycles, transposition, transposition],
                "terminal_packets": {"packet_1": [0, 1, 2]},
            },
            {
                "name": "C_affine",
                "connected": True,
                "cycles": [
                    [2, 3, 4, 5, 6, 1],
                    [6, 1, 2, 3, 4, 5],
                ],
                "terminal_packets": {},
            }
        ],
    }


payload = build_payload()
proximity = payload["source_original_proximity_resolution"]
assert proximity["orders_on_kummer_chart"]["X"] == -5
assert proximity["orders_on_kummer_chart"]["y_translated"] == 17
assert proximity["orders_on_original_source"] == {
    "x": -25,
    "y_old": 5,
    "v_minus_1_where_v=x*y_old^5": 12,
}
assert proximity["monomial_carrier_ray"] == [-5, 1]
assert len(proximity["standard_P2_carrier_fan"]["insertion_order_coordinates"]) == 6
assert len(proximity["principal_arm_fan"]["insertion_order_coordinates"]) == 6
assert proximity["standard_P2_carrier_fan"]["boundary_intersection_determinant"] == 1
target_completion = payload["target_minimal_completion"]["fan"]
assert target_completion["boundary_self_intersections"] == [-2, -2, -1, -3, -2]
assert target_completion["boundary_intersection_determinant"] == 1
assert payload["target_valuation_uniqueness"]["uniformizer_order"] == 1
assert "same target divisor" in payload["target_valuation_uniqueness"]["double_root_consequence"]
case_records = {record["case"]: record for record in payload["cases"]}
assert case_records["double_same_target"]["geometric_degree_floor"] == 12
assert case_records["double_distinct_targets"]["status"] == "excluded_by_target_valuation_uniqueness"
assert payload["target_local_fan"]["terminal_self_intersection_in_this_local_fan"] == -1
assert len(payload["attachment_slots"]) == 5
assert sum(slot["forces_new_interior_branch"] for slot in payload["attachment_slots"]) == 3
assert topology_bounds(1)["minimum_source_boundary_components"] == 19
assert topology_bounds(1)["minimum_source_boundary_leaves"] == 6
assert topology_bounds(2)["minimum_source_boundary_components"] == 31
assert topology_bounds(2)["minimum_source_boundary_leaves"] == 10

synthetic = synthetic_squarefree_candidate()
audit = audit_candidate(synthetic)
assert audit["status"] == "passes_compiled_necessary_gates_not_an_existence_proof"
assert abs(audit["source_boundary"]["intersection_determinant"]) == 1
assert audit["source_boundary"]["canonical_coefficients_in_component_order"] is not None
assert audit["target_ledgers"]["status"] == "passes_finite_normalization_and_purity_gates"
assert audit["global_meridians"]["status"] == "passes_declared_meridian_gates"

nonunimodular = deepcopy(synthetic)
nonunimodular["source_boundary"]["components"][0]["self_intersection"] += 1
bad_audit = audit_candidate(nonunimodular)
assert bad_audit["status"] == "excluded_candidate"
assert "not unimodular" in " ".join(bad_audit["source_boundary"]["reasons"])

incomplete = deepcopy(synthetic)
del incomplete["spectator_orbits"]
incomplete_audit = audit_candidate(incomplete)
assert incomplete_audit["status"] == "incomplete"
assert incomplete_audit["spectators"]["status"] == "incomplete"

bad_meridian = deepcopy(synthetic)
bad_meridian["global_meridian_systems"][0]["cycles"].pop()
bad_meridian_audit = audit_candidate(bad_meridian)
assert bad_meridian_audit["status"] == "excluded_candidate"
assert "product one" in " ".join(bad_meridian_audit["global_meridians"]["reasons"])

print("PASS: F2 target fan orients both toric node attachments")
print("PASS: source valuation compiles as a six-blowup carrier plus principal arm")
print("PASS: one/two terminal packets force the 19/31-component proximity trees")
print("PASS: candidate mode separates passing, incomplete, and excluded ledgers")
print("PASS: class/unit/canonical and meridian gates reject malformed candidates")
