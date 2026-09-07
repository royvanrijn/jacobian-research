#!/usr/bin/env sage -python
"""
status: ACTIVE_PROOF
claim: fixed q24/orbit85 D13->D12 exact-construction contract.
inputs: q24 equation corridor, Theta discovery replay, modular component RR,
        modular orbit85 D12 signature.
outputs: artifacts/local/elkies-k3/q24-orbit85-exact-construction-manifest.json

This is not the characteristic-zero equation lift.  It is the pinned
construction manifest for that lift: exact Theta/component data plus the
already-passing modular resolved-RR regression gates that the QQ compiler must
replay.
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"

CORRIDOR = LOCAL / "q24-equation-d13-to-pinned-r17.json"
THETA = LOCAL / "q24-d12-discovery-replay-mod-100003.json"
COMPONENT_RR = LOCAL / "q24-d12-component-valuation-rr-mod-100003.json"
SIGNATURE = LOCAL / "q24-orbit85-d12-signature-mod-100003.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
OUT = LOCAL / "q24-orbit85-exact-construction-manifest.json"

for path in (CORRIDOR, THETA, COMPONENT_RR, SIGNATURE, MANIFEST):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

corridor = json.loads(CORRIDOR.read_text())
theta = json.loads(THETA.read_text())
component = json.loads(COMPONENT_RR.read_text())
signature = json.loads(SIGNATURE.read_text())
manifest = json.loads(MANIFEST.read_text())

assert corridor["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert corridor["q24"]["orbit"] == 85
assert corridor["equation_lift_manifest"]["all_manifest_transitions_match_replay"]

assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
first = manifest["forward_steps"][0]
assert first["parent"] == "D13/MW4"
assert first["child"] == "D12/MW5"
assert first["q"] == 24
assert first["orbit"] == 85
assert first["equation_lift_hint"] == "RESOLVED_RR_PREFERRED"

assert theta["status"] == "PASS_EXACT_D12_DISCOVERY_REPLAY"
td = theta["pinned_D13"]
assert td["D12_root_data"] == [12,264,4]
assert td["removed_leaf"] == 3
assert td["theta_square"] == -2
assert td["theta_old_fibre_degree"] == 2
assert td["theta_old_zero_intersection"] == 10
assert td["theta_new_fibre_intersection"] == 0
assert td["theta_old_affine_intersection"] == 0
assert td["theta_old_component_hits"] == [[12,1]]
assert theta["modular_theta_diagnostic"]["candidate_is_rational"] is False
assert theta["modular_theta_diagnostic"]["bisection_squarefree_degree"] == 18

assert component["schema"] == "elkies-k3.h3-q24-d12-component-valuation-rr-modp.v1"
assert component["status"] == "CANDIDATE_H3_Q24_EFFECTIVE_D13_D12_MODP"
assert component["geometric_trivialization"] == {
    "deterministic_fibre_twist": -7,
    "geometric_fibre_twist": -8,
    "isotropic_square": 0,
    "vertical_square": -12,
}
assert component["global_rr"] == {
    "ambient_dimension": 56,
    "smooth_collision_rank": 48,
    "post_collision_dimension": 8,
}
cluster = component["resolved_cluster"]
assert cluster["plan"] == [
    {"center":"C01","additional_order":2},
    {"center":"C02","additional_order":2},
    {"center":"C04","additional_order":2},
    {"center":"C06","additional_order":2},
    {"center":"C08","additional_order":2},
]
assert cluster["method"] == "direct_membership_in_(surface,e^threshold)_on_full_chart_cover"
assert cluster["codimension_on_post_collision"] == 6
assert cluster["kernel_dimension"] == 2
assert len(cluster["condition_ledger"]) == 12
assert cluster["condition_ledger"][-1]["cumulative_rank"] == 6
assert component["quartic_degree"] == 4
assert component["child"]["root_rank"] == 12
assert component["child"]["root_determinant"] == 4
assert component["child"]["euler"] == 24

assert signature["status"] in (
    "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
    "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
)
assert signature["source_neighbor"] == {
    "child":"D12/MW5",
    "orbit":85,
    "q":24,
    "source":"D13/MW4",
}
assert signature["rr"] == {
    "ambient":56,
    "collision_rank":48,
    "geometric_fibre_twist":-8,
    "kernel":2,
    "post_collision":8,
    "resolved_rank":6,
}
assert signature["plane_pivots"] == [0,1]
assert signature["quartic_degree"] == 4
assert signature["child_root_rank"] == 12
assert signature["child_root_det"] == 4
assert signature["child_euler"] == 24

payload = {
    "schema":"elkies-k3.h3-q24-orbit85-exact-construction-manifest.v1",
    "status":"PASS_Q24_ORBIT85_EXACT_CONSTRUCTION_MANIFEST",
    "proof_boundary":(
        "Exact characteristic-zero construction contract plus modular "
        "resolved-RR regression. This manifest fixes the q24/orbit85 divisor "
        "as the exact Theta/component construction and records the modular "
        "56->8->2 RR/quartic/D12 gates. It does not yet provide QQ kernel "
        "coefficients or a characteristic-zero child Weierstrass equation."
    ),
    "inputs":{
        "corridor":str(CORRIDOR.relative_to(ROOT)),
        "theta_discovery":str(THETA.relative_to(ROOT)),
        "component_rr_modp":str(COMPONENT_RR.relative_to(ROOT)),
        "orbit85_signature_modp":str(SIGNATURE.relative_to(ROOT)),
        "backward_manifest":str(MANIFEST.relative_to(ROOT)),
    },
    "fixed_corridor_step":{
        "parent":"D13/MW4",
        "child":"D12/MW5",
        "q":24,
        "orbit":85,
        "equation_lift_hint":first["equation_lift_hint"],
        "transition_matches_replay":True,
    },
    "exact_theta_construction":{
        "fibre_decomposition":td["fibre_decomposition"],
        "q24_fibre":td["q24_fibre"],
        "theta_class":td["theta_class"],
        "D12_highest_root_in_D13_numbering":td["D12_highest_root_in_D13_numbering"],
        "theta_square":td["theta_square"],
        "theta_old_fibre_degree":td["theta_old_fibre_degree"],
        "theta_old_zero_intersection":td["theta_old_zero_intersection"],
        "theta_new_fibre_intersection":td["theta_new_fibre_intersection"],
        "theta_old_component_hits":td["theta_old_component_hits"],
        "rejected_generic_theta_candidate":{
            "construction":theta["modular_theta_diagnostic"]["construction"],
            "squarefree_degree":theta["modular_theta_diagnostic"]["bisection_squarefree_degree"],
            "candidate_is_rational":False,
        },
    },
    "resolved_rr_contract":{
        "geometric_trivialization":component["geometric_trivialization"],
        "affine_component_coefficients":cluster["affine_component_coefficients"],
        "centre_thresholds":cluster["centre_thresholds"],
        "cluster_plan":cluster["plan"],
        "global_rr":component["global_rr"],
        "resolved_rank":cluster["codimension_on_post_collision"],
        "kernel_dimension":cluster["kernel_dimension"],
        "condition_ledger":cluster["condition_ledger"],
        "modp_kernel_basis_post_collision":cluster["kernel_basis_post_collision"],
        "modp_kernel_basis_ambient":cluster["kernel_basis_ambient"],
    },
    "modular_regression_gates":{
        "prime":signature["prime"],
        "plane_pivots":signature["plane_pivots"],
        "quartic_degree":signature["quartic_degree"],
        "quartic_coefficients":signature["quartic_coefficients"],
        "jacobian_A":signature["jacobian_A"],
        "jacobian_B":signature["jacobian_B"],
        "child_root_rank":signature["child_root_rank"],
        "child_root_det":signature["child_root_det"],
        "child_euler":signature["child_euler"],
    },
    "qq_replay_requirements":[
        "Use the exact corrected D13 equation model, not a p-adic q24 section lift.",
        "Build the 56-dimensional geometric-trivialization ambient with fibre twist -8.",
        "Impose smooth collision to rank 48 / post-collision dimension 8.",
        "Replay the exact I9* component cover with thresholds recorded in centre_thresholds.",
        "Recover an exact 2D kernel for L(D24), compile the degree-two chord quartic, and require degree 3 or 4.",
        "Compile the binary-quartic Jacobian and certify exact D12/MW5 fibre data.",
        "Regress modulo 100003 to the stored 2x56 plane, quartic, Jacobian, and D12 signature.",
    ],
}

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O85EXACTMANIFEST|"
    "theta=1|ambient=56|collision=48|post=8|resolved=6|kernel=2|"
    "quartic=4|child=D12/MW5|status=PASS_Q24_ORBIT85_EXACT_CONSTRUCTION_MANIFEST",
    flush=True,
)
