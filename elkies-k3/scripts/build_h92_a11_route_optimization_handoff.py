#!/usr/bin/env python3
"""Build the compact machine-readable A11 route-optimization handoff.

status: ACTIVE_PROOF
claim: compact lossless handoff of the promoted route and rejected searches
inputs: exact route-search and certification artifacts listed in PATHS
outputs: artifacts/generated-results/elkies-k3-h3-a11-route-optimization-handoff.json
"""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
OUTPUT = GENERATED / "elkies-k3-h3-a11-route-optimization-handoff.json"
PATHS = {
    "zero_mismatch": GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json",
    "a11_ranking": GENERATED / "elkies-k3-h3-a11-marked-target-neighbor-ranking.json",
    "crossover": GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json",
    "o1991_certificate": GENERATED / "elkies-k3-h3-a11-q8-orbit1991-lattice-certificate.json",
    "o1991_explicit_zero": GENERATED / "elkies-k3-h3-a11-q8-orbit1991-explicit-zero-frames.json",
    "o1991_continuation": GENERATED / "elkies-k3-h3-o1991-marked-target-neighbor-ranking.json",
    "o12_certificate": GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json",
    "o12_explicit_zero": GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json",
    "a5_ranking": GENERATED / "elkies-k3-h3-a5a5-marked-target-neighbor-ranking.json",
    "current_route_cost": GENERATED / "elkies-k3-h3-a5a5-current-route-equation-cost-audit.json",
    "deterministic_q6q8_gate": GENERATED / "elkies-k3-h3-a5a5-q6q8-explicit-curve-gate.json",
    "explicit_q4q6_cost": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json",
    "orbit32_certificate": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4-orbit32-lattice-certificate.json",
    "orbit32_continuation": GENERATED / "elkies-k3-h3-a5a5-q4-o32-explicit-zero-q4q6-equation-cost.json",
    "orbit3372_certificate": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit3372-lattice-certificate.json",
    "orbit3372_continuation": GENERATED / "elkies-k3-h3-a5a5-q6-o3372-explicit-zero-q4q6-equation-cost.json",
    "candidate_suffix_crossovers": GENERATED / "elkies-k3-h3-candidate-current-suffix-crossovers.json",
    "q10_gate": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q10-gate.json",
    "q10_cost": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q10-equation-cost.json",
    "q12_gate": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q12-gate.json",
    "q12_cost": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q12-equation-cost.json",
    "semistable_pinned_transport": GENERATED / "elkies-k3-h3-semistable-mw2-pinned-transport.json",
    "semistable_nef_suffix": GENERATED / "elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json",
    "semistable_a5_ranking": GENERATED / "elkies-k3-h3-a5a5-semistable-hub-survivor-ranking.json",
    "bridge_section_qq": LOCAL / "q24-a11-bridge-m-section-qq.json",
    "semistable_marking": GENERATED / "elkies-k3-h3-semistable-mw2-equation-marking.json",
    "semistable_large_frontier": GENERATED / "elkies-k3-h3-semistable-mw2-q8q10q12-equation-cost.json",
    "semistable_q12o7798_second": GENERATED / "elkies-k3-h3-semistable-mw2-q12o7798-q4q6-equation-cost.json",
    "q25_mw7_marking": GENERATED / "elkies-k3-h3-q25-mw7-equation-marking.json",
    "q25_mw7_frontier": GENERATED / "elkies-k3-h3-q25-mw7-q4q6-marked-frontier-compact.json",
    "q25_mw7_second": GENERATED / "elkies-k3-h3-q25-mw7-q6o36810-q4-marked-frontier-compact.json",
    "q25_mw4_marking": GENERATED / "elkies-k3-h3-q25-mw4-equation-marking.json",
    "q25_mw4_frontier": GENERATED / "elkies-k3-h3-q25-mw4-q4q6-marked-frontier-compact.json",
    "mw3_marking": GENERATED / "elkies-k3-h3-mw3-a5-d4-2a2-a1-equation-marking.json",
    "mw3_frontier": GENERATED / "elkies-k3-h3-mw3-a5-d4-2a2-a1-q4q6-marked-frontier-compact.json",
    "a11_degree3_gate": GENERATED / "elkies-k3-h3-a11-q6q9q12-degree3-full-nef-equation-cost.json",
    "a11_q9o1802_certificate": GENERATED / "elkies-k3-h3-a11-q9d3o1802-lattice-certificate.json",
    "a11_q9o1802_frontier": GENERATED / "elkies-k3-h3-a11-q9d3o1802-q4q6-marked-frontier-compact.json",
    "a11_q9o2800_certificate": GENERATED / "elkies-k3-h3-a11-q9d3o2800-lattice-certificate.json",
    "a11_q9o2800_frontier": GENERATED / "elkies-k3-h3-a11-q9d3o2800-q4q6-marked-frontier-compact.json",
    "a11_q9o2793_certificate": GENERATED / "elkies-k3-h3-a11-q9d3o2793-lattice-certificate.json",
    "a11_q9o2793_frontier": GENERATED / "elkies-k3-h3-a11-q9d3o2793-q4q6-marked-frontier-compact.json",
    "a11_q9o1542_certificate": GENERATED / "elkies-k3-h3-a11-q9d3o1542-lattice-certificate.json",
    "a11_q9o1542_frontier": GENERATED / "elkies-k3-h3-a11-q9d3o1542-q4q6-marked-frontier-compact.json",
    "a11_q9o956_certificate": GENERATED / "elkies-k3-h3-a11-q9d3o956-lattice-certificate.json",
    "a11_q9o956_frontier": GENERATED / "elkies-k3-h3-a11-q9d3o956-q4q6-marked-frontier-compact.json",
    "d13_marking": GENERATED / "elkies-k3-h3-equation-d13-marking.json",
    "d13_low_frontier": GENERATED / "elkies-k3-h3-equation-d13-q4q6q8-marked-frontier-adapted-compact.json",
    "d13_degree2_frontier": GENERATED / "elkies-k3-h3-equation-d13-q10to20-degree2-marked-frontier-compact.json",
    "d13_degree3_frontier": GENERATED / "elkies-k3-h3-equation-d13-q9to18-degree3-marked-frontier-compact.json",
    "d13_degree4_frontier": GENERATED / "elkies-k3-h3-equation-d13-q8to20-degree4-marked-frontier-compact.json",
    "d13_q6o42_checkpoint": GENERATED / "elkies-k3-h3-equation-d13-q6o42-marking.json",
    "d13_q6o42_second": GENERATED / "elkies-k3-h3-equation-d13-q6o42-q4q6-marked-frontier-compact.json",
    "pinned_marking": GENERATED / "elkies-k3-h3-pinned-r17-equation-marking.json",
    "pinned_q4_exact": GENERATED / "elkies-k3-h3-pinned-r17-q4-degree2-targeted-ranking.json",
    "pinned_orbit12_cvp": GENERATED / "elkies-k3-h3-pinned-r17-orbit12-targeted-shell-cvp.json",
    "pinned_q6_orbit12": GENERATED / "elkies-k3-h3-pinned-r17-q6-orbit12-cvp-lattice-certificate.json",
    "pinned_q8_orbit12": GENERATED / "elkies-k3-h3-pinned-r17-q8-orbit12-cvp-lattice-certificate.json",
    "pinned_q8_second": GENERATED / "elkies-k3-h3-pinned-r17-q8-orbit12-cvp-q4q6-cap20000-marked-frontier-compact.json",
    "orbit12_q8q18": GENERATED / "elkies-k3-h3-pinned-r17-q8q18-orbit12-cvp-lattice-certificate.json",
    "orbit12_q8q18q10": GENERATED / "elkies-k3-h3-pinned-r17-q8q18q10-orbit12-cvp-lattice-certificate.json",
    "orbit12_q8q18q10q14": GENERATED / "elkies-k3-h3-pinned-r17-q8q18q10q14-orbit12-cvp-lattice-certificate.json",
    "pair49893": GENERATED / "elkies-k3-h3-pinned-r17-q4-pair49893-lattice-certificate.json",
    "pair49893_second": GENERATED / "elkies-k3-h3-pinned-r17-q4-pair49893-q4-marked-frontier-compact.json",
    "pinned_mw3_degree3_probe": GENERATED / "elkies-k3-h3-pinned-r17-mw3-degree3-targeted-shell-cvp.json",
    "pinned_q12d3_mw3": GENERATED / "elkies-k3-h3-pinned-r17-q12d3-mw3-cvp-lattice-certificate.json",
    "q12d3_to_q25_probe": GENERATED / "elkies-k3-h3-pinned-r17-q12d3-mw3-cvp-child-q24q36-targeted-shell-cvp.json",
    "q12d3_q25_endpoint": GENERATED / "elkies-k3-h3-pinned-r17-q12d3-q28d2-q25mw7-lattice-certificate.json",
    "q25_degree3_detour": GENERATED / "elkies-k3-h3-q25mw7-pinned-r17-degree2-degree3-detour.json",
    "mw3_correct_scale_negative": GENERATED / "elkies-k3-h3-mw3-a5-d4-2a2-a1-orbit12-q24q40-targeted-shell-cvp.json",
    "equation_a11_marking": GENERATED / "elkies-k3-h3-equation-a11-marking.json",
    "equation_a11_degree2_targeted": GENERATED / "elkies-k3-h3-equation-a11-pinned-q8q16-targeted-shell-cvp.json",
    "equation_a11_degree3_targeted": GENERATED / "elkies-k3-h3-equation-a11-pinned-degree3-q15q24-targeted-shell-cvp.json",
    "equation_a11_degree4_targeted": GENERATED / "elkies-k3-h3-equation-a11-pinned-degree4-q24q40-targeted-shell-cvp.json",
    "orbit12_beam_q8_branch": GENERATED / "elkies-k3-h3-pinned-r17-q8q18-q8beam-orbit12-child-targeted-shell-cvp.json",
    "orbit12_beam_q12_branch": GENERATED / "elkies-k3-h3-pinned-r17-q8q18-q12beam-orbit12-child-targeted-shell-cvp.json",
    "pinned_q8_current_a1_identification": GENERATED / "elkies-k3-h3-pinned-r17-q8-current-route-a1-identification.json",
    "pinned_current_suffix_marking": GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json",
    "suffix_target_3A3": GENERATED / "elkies-k3-h3-pinned-r17-current_3A3-q4q10-targeted-shell-cvp.json",
    "suffix_target_A3_2A2": GENERATED / "elkies-k3-h3-pinned-r17-current_A3_2A2-q4q10-targeted-shell-cvp.json",
    "suffix_target_5A1": GENERATED / "elkies-k3-h3-pinned-r17-current_5A1-q4q10-targeted-shell-cvp.json",
    "suffix_target_4A1": GENERATED / "elkies-k3-h3-pinned-r17-current_4A1-q4q10-targeted-shell-cvp.json",
    "suffix_target_3A1": GENERATED / "elkies-k3-h3-pinned-r17-current_3A1-q4q10-targeted-shell-cvp.json",
    "suffix_target_2A1": GENERATED / "elkies-k3-h3-pinned-r17-current_2A1-q4q10-targeted-shell-cvp.json",
    "explicit_zero_detour_route": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-detour-route-certificate.json",
    "orbit3372_marking": GENERATED / "elkies-k3-h3-a5a5-q6o3372-suffix-marking.json",
    "orbit32_marking": GENERATED / "elkies-k3-h3-a5a5-q4o32-suffix-marking.json",
    "orbit3372_pinned_targeted": GENERATED / "elkies-k3-h3-a5a5-q6o3372-pinned-q4q24-targeted-shell-cvp.json",
    "orbit32_pinned_targeted": GENERATED / "elkies-k3-h3-a5a5-q4o32-pinned-q4q24-targeted-shell-cvp.json",
    "zero_loop_search": GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json",
    "orbit1307_certificate": GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit1307-lattice-certificate.json",
    "orbit1307_explicit_zero": GENERATED / "elkies-k3-h3-a5a5-q6-o1307-explicit-zero-frames.json",
    "orbit1307_return": GENERATED / "elkies-k3-h3-a5a5-q6o1307-q4-return-a5a5-certificate.json",
    "orbit1307_exit": GENERATED / "elkies-k3-h3-a5a5-q6o1307-loop-current-3a3-certificate.json",
    "promoted_route": GENERATED / "elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json",
    "q230_second_zero_search": GENERATED / "elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-presentations.json",
    "q230_promoted_route": GENERATED / "elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json",
    "a11_q8_zero_loop_search": GENERATED / "elkies-k3-h3-a11-zero-changing-q8-presentations.json",
    "d13_zero_loop_search": GENERATED / "elkies-k3-h3-d13-zero-changing-d12-presentations.json",
    "d13_promoted_route": GENERATED / "elkies-k3-h3-d13-q4o11-promoted-route-certificate.json",
    "first_q8_source_marking": GENERATED / "elkies-k3-h3-first-q8-source-marking.json",
    "first_q8_zero_loop_search": GENERATED / "elkies-k3-h3-first-q8-zero-changing-d13-presentations.json",
    "first_q8_returned_marking": GENERATED / "elkies-k3-h3-first-q8-q4o11-c1-returned-marking.json",
    "first_q8_second_loop_search": GENERATED / "elkies-k3-h3-first-q8-q4o11-c1-second-zero-changing-d13-presentations.json",
    "first_q8_q10_search": GENERATED / "elkies-k3-h3-first-q8-source-q10-zero-changing-d13-presentations.json",
    "first_q8_degree3_search": GENERATED / "elkies-k3-h3-first-q8-source-q6q9q12-degree3-zero-changing-d13-presentations.json",
    "first_q8_degree4_search": GENERATED / "elkies-k3-h3-first-q8-source-q8q12q16-degree4-zero-changing-d13-presentations.json",
    "first_q8_promoted_route": GENERATED / "elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json",
    "first_q8_landing_d13_marking": GENERATED / "elkies-k3-h3-first-q8-q4o11-landing-d13-marking.json",
    "first_q8_landing_d13_crossover": GENERATED / "elkies-k3-h3-first-q8-q4o11-landing-d13-zero-changing-d12-presentations.json",
    "physical_q4o208_promoted_route": GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json",
}


def load(name):
    return json.loads(PATHS[name].read_text())


def summary(item):
    result = {
        "candidate_id": item["candidate_id"],
        "equation_cost_score": item.get("equation_cost_score"),
        "P_dot_O": item.get("P_dot_O", item.get("horizontal", {}).get("P_dot_O")),
        "expected_RR_ambient": item.get("expected_RR_ambient"),
        "pinned_R17_degree": item.get("marked_target_degrees", {}).get("pinned_R17"),
        "explicit_degree_zero_count": len(item.get("explicit_degree_zero_curves", [])),
        "explicit_degree_one_count": len(item.get("explicit_degree_one_curves", [])),
    }
    if "child" in item:
        result["child"] = item["child"]
    return result


def compact_marked_candidate(item, target="orbit12"):
    degrees = item["marked_target_degrees"]
    return {
        "candidate_id": item["candidate_id"],
        "child": item.get("child"),
        "target_degree": degrees[target],
        "pinned_R17_degree": degrees.get("pinned_R17"),
        "minimum_section_intersection": item.get("minimum_section_intersection"),
        "P_dot_O": item.get("P_dot_O"),
        "component_degree_zero_count": item.get("component_degree_zero_count"),
        "component_degree_one_count": item.get("component_degree_one_count"),
    }


def compact_frontier(data, target="orbit12", direct_degree=None):
    ranked = data["ranked_candidates"]
    result = {
        "status": data["status"],
        "input_candidate_count": data["input_candidate_count"],
        "full_nef_candidate_count": data["full_nef_candidate_count"],
        "best": compact_marked_candidate(ranked[0], target),
    }
    if direct_degree is not None:
        genuine = next(
            (item for item in ranked if item["marked_target_degrees"][target] > direct_degree),
            None,
        )
        result["first_genuine_nonbacktrack"] = (
            compact_marked_candidate(genuine, target) if genuine else None
        )
    return result


def marked_degree(marking, target):
    return marking["target_fibres_in_root_adapted_hub"][target][1]


def compact_semistable_candidate(item):
    return {
        **summary(item),
        "orbit12_degree": item["marked_target_degrees"]["orbit12_fibre"],
        "pinned_R17_degree": item["marked_target_degrees"]["pinned_R17_fibre"],
    }


a11 = load("a11_ranking")
o1991_continuation = load("o1991_continuation")
explicit_cost = load("explicit_q4q6_cost")
o32 = load("orbit32_continuation")
o3372 = load("orbit3372_continuation")
o12_zero = load("o12_explicit_zero")
a5_ranking = load("a5_ranking")
deterministic_gate = load("deterministic_q6q8_gate")
current_route_cost = load("current_route_cost")
suffix_crossovers = load("candidate_suffix_crossovers")
bridge_section = load("bridge_section_qq")
semistable_large = load("semistable_large_frontier")
semistable_second = load("semistable_q12o7798_second")
q25_mw7_marking = load("q25_mw7_marking")
q25_mw7_frontier = load("q25_mw7_frontier")
q25_mw7_second = load("q25_mw7_second")
q25_mw4_marking = load("q25_mw4_marking")
q25_mw4_frontier = load("q25_mw4_frontier")
mw3_marking = load("mw3_marking")
mw3_frontier = load("mw3_frontier")
a11_degree3 = load("a11_degree3_gate")
a11_beam_names = ["q9o1802", "q9o2800", "q9o2793", "q9o1542", "q9o956"]
d13_marking = load("d13_marking")
d13_frontier_names = [
    "d13_low_frontier", "d13_degree2_frontier",
    "d13_degree3_frontier", "d13_degree4_frontier",
]
pinned_marking = load("pinned_marking")
pinned_orbit12_direct = pinned_marking["target_fibres_in_root_adapted_hub"]["orbit12"][1]
d13_promotion = load("d13_promoted_route")
q230_promotion = load("q230_promoted_route")
a11_q8_zero_loop = load("a11_q8_zero_loop_search")
d13_zero_loop = load("d13_zero_loop_search")
first_q8_promotion = load("first_q8_promoted_route")
first_q8_search = load("first_q8_zero_loop_search")
first_q8_second = load("first_q8_second_loop_search")
physical_q4o208 = load("physical_q4o208_promoted_route")

payload = {
    "schema": "elkies-k3.h3-a11-route-optimization-handoff.v3",
    "status": physical_q4o208["status"],
    "generated_on": "2026-08-25",
    "promotion": {
        "new_lifting_target": True,
        "active_lifting_route_unchanged": False,
        "route": physical_q4o208["promotion"],
        "full_route_q_sequence_from_H3": [6, 8, 24, 6] + physical_q4o208["full_route_q_sequence_from_A11"],
        "reason": "After the exact A11 q8/orbit12 equation reaches component-9-zero 2A5/MW7, use physical q4/orbit208. Its exact RR/Jacobian compilation has dimensions 4-to-2-to-2 and fibres 3I4+12I1; its full marked lattice route lands on canonical current 3A3 and pinned R17. Only equation-effective C5 pointing and full old-curve marking remain.",
    },
    "first_A11_q8_zero_changing_audit": {
        "status": a11_q8_zero_loop["status"],
        "declared_nef_first_edges": a11_q8_zero_loop["counts"]["first_edges"],
        "exact_first_edges": a11_q8_zero_loop["counts"]["first_exact_nef"],
        "exact_presentations": a11_q8_zero_loop["counts"]["exit_exact_nef"],
        "nonprimitive_explicit_root_spans": a11_q8_zero_loop["counts"]["nonprimitive_explicit_root_spans"],
        "direct_orbit12_inherited_explicit_score": a11_q8_zero_loop["direct_q8_orbit12"]["inherited_explicit_score"],
        "best_retained": {
            "candidate_id": a11_q8_zero_loop["ranked_presentations"][0]["first_edge_candidate_id"],
            "explicit_zero_curve": a11_q8_zero_loop["ranked_presentations"][0]["explicit_zero_curve"],
            "q_sequence": a11_q8_zero_loop["ranked_presentations"][0]["q_sequence"],
            "inherited_explicit_score": a11_q8_zero_loop["ranked_presentations"][0]["inherited_explicit_equation_cost"]["total_score"],
            "inherited_explicit_operational_score": a11_q8_zero_loop["ranked_presentations"][0]["inherited_explicit_equation_cost"]["operational_total_score"],
        },
        "strict_retained_winner_count": a11_q8_zero_loop["strict_winner_count"],
        "decision": "retain q8/orbit12; no inherited-explicit operational-cost zero loop is cheaper",
        "boundary": "All 24 nonprimitive full-root-span re-zeroings are handled by saturated unimodular frames retaining the embedded actual simple-root lattice.",
    },
    "critical_audit": {
        "quintic_bridge_status": load("zero_mismatch")["status"],
        "conclusion": "The documented A0-zero quintic AJ composition mixes the A0 and R3 D12 zero frames; do not use it as a certified equation-side bridge.",
    },
    "a11_exhaustive_marked_ranking": {
        "candidate_count": a11["candidate_count"],
        "orbit12": a11["named_orbits"]["12"],
        "orbit2162": a11["named_orbits"]["2162"],
        "conclusion": "Orbit12 is the unique marked-distance winner; branch after orbit12 rather than replacing it at A11.",
    },
    "rejected_lateral_branch": {
        "route": ["A11 q8 orbit1991", "explicit-zero A1+A2+D10/MW4"],
        "certificate_status": load("o1991_certificate")["status"],
        "continuation_candidate_count": o1991_continuation["candidate_count"],
        "best_continuation": o1991_continuation["rankings_top_100"]["pinned_R17"][0],
        "conclusion": "Every tested q4/q6/q8/q10 continuation worsens the exact pinned degree.",
    },
    "orbit12_zero_marking": {
        "certificate_status": load("o12_certificate")["status"],
        "selected_explicit_zero": o12_zero["selected_zero_curve"],
        "selected_frame_max_abs": o12_zero["selected"]["frame_max_abs"],
        "historical_route_q_in_historical_zero": a5_ranking["current_route"]["historical_q"],
        "historical_route_q_in_equation_deterministic_zero": a5_ranking["current_route"]["equation_zero_q"],
        "historical_zero_explicit_curve_matches": a5_ranking["current_route"]["historical_zero_explicit_curve_matches"],
        "deterministic_zero_q6q8_survivors": deterministic_gate["survivor_count"],
        "current_2A5_to_3A3_in_selected_explicit_zero": current_route_cost["explicit_zero_edge"],
        "current_2A5_to_3A3_cost_profile": current_route_cost.get(
            "direct_equation_cost_profile",
            current_route_cost.get("stored_nonphysical_equation_cost_profile"),
        ),
    },
    "certified_compiler_branches": [
        {
            "route": ["A11 q8 orbit12", "explicit zero old_A11_component_9", "q4 orbit32"],
            "certificate_status": load("orbit32_certificate")["status"],
            "candidate": summary(next(item for item in explicit_cost["ranked_candidates"] if item["candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 32})),
            "next_shell_pareto": [summary(item) for item in o32["pinned_degree_equation_cost_pareto_front"]],
            "endpoint_status": "not connected to pinned R17",
        },
        {
            "route": ["A11 q8 orbit12", "explicit zero old_A11_component_9", "q6 orbit3372"],
            "certificate_status": load("orbit3372_certificate")["status"],
            "candidate": summary(next(item for item in explicit_cost["ranked_candidates"] if item["candidate_id"] == {"q": 6, "old_fibre_degree": 2, "orbit_index": 3372})),
            "next_shell_pareto": [summary(item) for item in o3372["pinned_degree_equation_cost_pareto_front"]],
            "endpoint_status": "not connected to pinned R17",
        },
    ],
    "larger_q_exhaustive_gates": {
        "q10": {
            "enumeration": load("q10_gate")["summaries"][0],
            "best_equation_cost": summary(load("q10_cost")["ranked_candidates"][0]),
            "pareto_count": len(load("q10_cost")["pinned_degree_equation_cost_pareto_front"]),
            "conclusion": "Every explicit-curve/parent-affine-nef q10 candidate is compiler-cost dominated by a certified q4 or q6 branch.",
        },
        "q12": {
            "enumeration": load("q12_gate")["summaries"][0],
            "best_equation_cost": summary(load("q12_cost")["ranked_candidates"][0]),
            "pareto_count": len(load("q12_cost")["pinned_degree_equation_cost_pareto_front"]),
            "conclusion": "Every explicit-curve/parent-affine-nef q12 candidate is compiler-cost dominated by a certified q4 or q6 branch.",
        },
        "proof_boundary": "The q10/q12 enumerations and cost profiles are exact, but dominated candidates were not given child-root/all-section/transport certificates.",
    },
    "candidate_to_certified_suffix_crossovers": {
        "status": suffix_crossovers["status"],
        "best_by_candidate": suffix_crossovers["best_by_candidate"],
        "later_suffix_records": [
            {
                "candidate": item["candidate"],
                "target_stage": item["target_stage"],
                "old_fibre_degree": item["old_fibre_degree"],
                "q_in_candidate_zero": item["q_in_candidate_zero"],
                "minimum_section_intersection": item["minimum_section_intersection"],
                "negative_explicit_curves": item["negative_explicit_curves"],
            }
            for item in suffix_crossovers["records"]
            if item["target_stage_index"] > 2
        ],
        "conclusion": "The only nef crossover is the common 2A5 parent. Every genuine later-suffix target is huge and fails a known section wall, so no pinned splice is certified.",
    },
    "reverse_hub_observation": {
        "semistable_A5_A4_2A3_MW2": {
            "status": load("semistable_nef_suffix")["status"],
            "transport_status": load("semistable_pinned_transport")["status"],
            "route": load("semistable_nef_suffix")["route"],
            "reverse_degrees": [item["nef_audit"]["old_fibre_degree"] for item in load("semistable_nef_suffix")["steps"]],
            "minimum_section_intersections": [item["nef_audit"]["minimum_section_intersection"] for item in load("semistable_nef_suffix")["steps"]],
            "a11_q8_best_orbit": a11["rankings_top_100"]["mw2_a5_a4_2a3_semistable"][0],
            "orbit12_explicit_zero_direct_degree": next(
                item["direct_old_fibre_degree"] for item in load("crossover")["records"]
                if item["state"] == "q8_orbit12_explicit_zero" and item["target"] == "mw2_a5_a4_2a3_semistable"
            ),
            "best_gated_orbit12_exit": load("semistable_a5_ranking")["top_100"][0],
            "conclusion": "This is now a fully lattice-certified pinned reverse hub. Orbit12 remains the best A11 q8 approach, but no cheap certified A11-to-hub continuation has been found.",
        },
        "exact_crossover_audit": str(PATHS["crossover"].relative_to(ROOT)),
    },
    "lifting_agent_exchange": {
        "status": bridge_section["status"],
        "equation_P_dot_O": bridge_section["section"]["equation_P_dot_O"],
        "pinned_lattice_P_dot_O": bridge_section["section"]["pinned_lattice_P_dot_O"],
        "degrees_X_Y_Z": bridge_section["section"]["degrees_X_Y_Z"],
        "max_coefficient_bits": bridge_section["section"]["max_coefficient_bits"],
        "conclusion": "The exact QQ bridge section exists, but its 909707-bit coefficients make the current suffix an exceptionally expensive compiler target and strengthen the case for a route change only after a complete pinned lattice splice is found.",
    },
    "bidirectional_equation_cost_search": {
        "semistable_MW2_reverse_beam": {
            "marking_status": load("semistable_marking")["status"],
            "input_candidate_count": semistable_large["input_candidate_count"],
            "full_nef_candidate_count": semistable_large["full_nef_candidate_count"],
            "best_equation_cost": compact_semistable_candidate(semistable_large["ranked_candidates"][0]),
            "closest_to_orbit12": compact_semistable_candidate(min(
                semistable_large["ranked_candidates"],
                key=lambda item: item["marked_target_degrees"]["orbit12_fibre"],
            )),
            "closest_checkpoint_second_layer": compact_semistable_candidate(min(
                semistable_second["ranked_candidates"],
                key=lambda item: item["marked_target_degrees"]["orbit12_fibre"],
            )),
            "conclusion": "The exact q8/q10/q12 reverse beam and the closest checkpoint's q4/q6 continuation both remain billions to trillions of marked degrees from orbit12; no meeting corridor exists in these shells.",
        },
        "q25_4_4_4_reverse_hubs": {
            "q25_mw7": {
                "direct_orbit12_degree": marked_degree(q25_mw7_marking, "orbit12"),
                "first_layer": compact_frontier(q25_mw7_frontier),
                "second_layer": compact_frontier(
                    q25_mw7_second,
                    direct_degree=marked_degree(q25_mw7_marking, "orbit12"),
                ),
            },
            "q25_mw4": {
                "direct_orbit12_degree": marked_degree(q25_mw4_marking, "orbit12"),
                "first_layer": compact_frontier(q25_mw4_frontier),
            },
            "mw3_A5_D4_2A2_A1": {
                "direct_orbit12_degree": marked_degree(mw3_marking, "orbit12"),
                "first_layer": compact_frontier(mw3_frontier),
            },
            "conclusion": "All exact low-q reverse frontiers worsen the direct orbit12 degree. The only q25-MW7 second-layer return is the original hub; its first genuine nonbacktrack is worse, so the proposed q25,4,4,4 crossover is not present in the searched shells.",
        },
    },
    "a11_non_degree_two_exits": {
        "status": a11_degree3["status"],
        "input_candidate_count": a11_degree3["input_candidate_count"],
        "full_nef_candidate_count": a11_degree3["full_nef_candidate_count"],
        "rejected_apparent_cost_leader": a11_degree3["rejected_cost_leader_q6d3o385"],
        "best_certified_equation_cost": summary(a11_degree3["best_candidate"]),
        "orbit12_degree_equation_cost_pareto": [
            summary(item) for item in a11_degree3["orbit12_degree_equation_cost_pareto_front"]
        ],
        "certified_second_layer_beam": [
            {
                "first_edge_certificate_status": load(f"a11_{name}_certificate")["status"],
                "first_edge": summary(load(f"a11_{name}_certificate")),
                "second_layer": compact_frontier(load(f"a11_{name}_frontier")),
            }
            for name in a11_beam_names
        ],
        "conclusion": "Five exact degree-three Pareto exits have full marked-U, root, nef, and bidirectional unimodular transport certificates. Every complete q4/q6 second layer moves away from orbit12, and the orbit12-degree-2 exit saves only 229 cost units before paying for an additional edge; none gives a cheaper pinned route.",
    },
    "earlier_equation_D13_branch": {
        "marking_status": d13_marking["status"],
        "current_q24_D12_to_pinned_R17_degree": d13_marking["current_D12_to_pinned_R17_degree"],
        "frontiers": [compact_frontier(load(name), target="pinned_R17") for name in d13_frontier_names],
        "best_compiler_friendly_checkpoint": {
            "status": load("d13_q6o42_checkpoint")["status"],
            "candidate": compact_marked_candidate(
                load("d13_low_frontier")["ranked_candidates"][0], "pinned_R17"
            ),
            "child": load("d13_q6o42_checkpoint")["child"],
            "second_layer": compact_frontier(load("d13_q6o42_second"), target="pinned_R17"),
        },
        "zero_changing_q24_replacement": {
            "search_status": d13_zero_loop["status"],
            "search_counts": d13_zero_loop["counts"],
            "strict_winner_count": d13_zero_loop["strict_winner_count"],
            "promoted_certificate_status": d13_promotion["status"],
            "splice": d13_promotion["new_D13_splice"],
            "combined_bottleneck_comparison": d13_promotion["combined_bottleneck_comparison"],
        },
        "conclusion": "Direct-to-pinned D13 beams remain poor, but changing the D13 zero through q4/orbit11 makes the existing q24 D12 exit cheaper. The q4,q4,q24 splice is fully certified to the exact current D12 basis and pinned R17 and is now promoted.",
    },
    "target_directed_pinned_reverse_search": {
        "method_status": load("pinned_orbit12_cvp")["status"],
        "q4_exact_calibration": load("pinned_orbit12_cvp")["q4_exact_calibration"],
        "direct_pinned_to_orbit12_degree": pinned_orbit12_direct,
        "current_route_a1_identification": {
            "status": load("pinned_q8_current_a1_identification")["status"],
            "full_unimodular_marking_identification": load("pinned_q8_current_a1_identification")["full_unimodular_marking_identification"],
            "component_chamber_aligned": load("pinned_q8_current_a1_identification")["simple_and_affine_A1_component_chamber_aligned"],
            "edge_profiles": load("pinned_q8_current_a1_identification")["edge_profiles"],
            "conclusion": load("pinned_q8_current_a1_identification")["conclusion"],
        },
        "certified_greedy_chain": [
            {
                "certificate_status": load(name)["status"],
                "candidate_id": load(name)["candidate_id"],
                "child": load(name)["child"],
                "orbit12_degree": load(name)["marked_target_degrees"]["orbit12"],
            }
            for name in (
                "pinned_q8_orbit12", "orbit12_q8q18", "orbit12_q8q18q10",
                "orbit12_q8q18q10q14",
            )
        ],
        "best_certified_marked_degree": load("orbit12_q8q18q10q14")["marked_target_degrees"]["orbit12"],
        "bounded_q8_second_layer": compact_frontier(load("pinned_q8_second")),
        "non_greedy_branches": {
            "q8_branch_best_next_degree": min(
                item["marked_target_degrees"]["orbit12"]
                for search in load("orbit12_beam_q8_branch")["searches"].values()
                for item in search["rankings"]
            ),
            "q12_branch_best_next_degree": min(
                item["marked_target_degrees"]["orbit12"]
                for search in load("orbit12_beam_q12_branch")["searches"].values()
                for item in search["rankings"]
            ),
            "conclusion": "Both alternative 2A1/MW15 branches return to the same 1.49M/1.54M local cluster and do not cross the 1,461,358 greedy endpoint.",
        },
        "conclusion": (
            "Target-directed CVP is calibrated by reproducing the exhaustive pinned q4 optimum. "
            "The q8 state is exactly the current-route A1/MW16 fibre under a full unimodular "
            "marking identification. From that existing suffix state, the lateral q18,q10,q14 "
            "chain lowers the orbit12 marked degree to 1,461,358 but ends in a 7-root/MW10 "
            "local diamond; tested continuations return to earlier states. It is not a new "
            "A11 meeting route."
        ),
    },
    "new_exact_q25_corridors": {
        "degree_two_pair49893": {
            "certificate_status": load("pair49893")["status"],
            "intermediate_root_data": load("pair49893")["child"]["root_data"],
            "q25_reverse_edge": load("pair49893")["target_transports"]["q25_mw7"]["reverse_edge_profile"],
            "pinned_reverse_edge": load("pair49893")["target_profiles"]["pinned_R17"],
            "complete_q4_second_layer": compact_frontier(load("pair49893_second")),
            "conclusion": "Exact q25 -> 6A1/MW11 -> pinned corridor with degree-two edges q8/P.O2 and q16/P.O6; its complete q4 continuation does not approach orbit12.",
        },
        "degree_three_intermediate": {
            "composed_status": load("q25_degree3_detour")["status"],
            "route": load("q25_degree3_detour")["route"],
            "endpoint_identification": load("q25_degree3_detour")["endpoint_identification"],
            "comparison": load("q25_degree3_detour")["comparison"],
            "mw3_contact_degree": load("pinned_q12d3_mw3")["marked_target_degrees"]["mw3_a5_d4_2a2_a1"],
            "q25_endpoint_certificate_status": load("q12d3_q25_endpoint")["status"],
            "conclusion": "Exact physical q25 -> 2A2+2A1/MW11 -> pinned corridor with q6,d2,P.O1 then q39,d3,P.O10. It lowers the first P.O and gives MW3 contact 250, but is not clearly cheaper than the all-degree-two pair49893 corridor.",
        },
        "mw3_correct_scale_negative": {
            "status": load("mw3_correct_scale_negative")["status"],
            "q_values": load("mw3_correct_scale_negative")["search_parameters"]["q_values"],
            "primitive_nef_counts": {
                q: data["primitive_nef_fixed_norm_candidates"]
                for q, data in load("mw3_correct_scale_negative")["searches"].items()
            },
            "conclusion": "The bounded target-directed degree-two scan at the physical MW3 frame's correct q24..q40 scale found no nef candidate in its fixed component chamber.",
        },
        "promotion_decision": "NO_CHANGE_WITHOUT_CERTIFIED_A11_OR_ORBIT12_PREFIX",
    },
    "equation_a11_correct_scale_search": {
        "marking_status": load("equation_a11_marking")["status"],
        "degree_two": {
            q: {
                "nef_candidates": data["primitive_nef_fixed_norm_candidates"],
                "best_pinned_degree": (
                    data["rankings"][0]["marked_target_degrees"]["pinned_R17"]
                    if data["rankings"] else None
                ),
            }
            for q, data in load("equation_a11_degree2_targeted")["searches"].items()
        },
        "degree_three": {
            q: {
                "nef_candidates": data["primitive_nef_fixed_norm_candidates"],
                "best_pinned_degree": (
                    data["rankings"][0]["marked_target_degrees"]["pinned_R17"]
                    if data["rankings"] else None
                ),
            }
            for q, data in load("equation_a11_degree3_targeted")["searches"].items()
        },
        "degree_four": {
            q: data["primitive_nef_fixed_norm_candidates"]
            for q, data in load("equation_a11_degree4_targeted")["searches"].items()
        },
        "conclusion": (
            "At equation A11 the target-directed q8 degree-two probe reproduces orbit12 "
            "at pinned degree 29,900,919. q10/q12/q16 are at least 1.406 billion and q14 "
            "has no hit. The previously untested correct-scale degree-three q15..q24 "
            "and degree-four q24..q40 samples give no competitor; q15 degree three is "
            "already 1.446 billion and all higher tested shells are empty."
        ),
    },
    "current_suffix_bidirectional_target_search": {
        "marking_status": load("pinned_current_suffix_marking")["status"],
        "exact_overlap_checks": load("pinned_current_suffix_marking")["exact_overlap_checks"],
        "targets": {
            target: {
                q: {
                    "nef_candidates": data["primitive_nef_fixed_norm_candidates"],
                    "best_target_degree": (
                        data["rankings"][0]["marked_target_degrees"][target]
                        if data["rankings"] else None
                    ),
                }
                for q, data in load(path_name)["searches"].items()
            }
            for target, path_name in (
                ("current_3A3", "suffix_target_3A3"),
                ("current_A3_2A2", "suffix_target_A3_2A2"),
                ("current_5A1", "suffix_target_5A1"),
                ("current_4A1", "suffix_target_4A1"),
                ("current_3A1", "suffix_target_3A1"),
                ("current_2A1", "suffix_target_2A1"),
            )
        },
        "conclusion": (
            "All q4/q6/q8/q10 target-directed probes into the middle and late current "
            "suffix select the existing current-route A1 fibre as the q8 winner. Its "
            "degrees against 3A3, A3+2A2, 5A1, 4A1, 3A1, and 2A1 are respectively "
            "365373, 55374, 2761, 149, 8, and 2. No new two-edge degree-two crossover "
            "skips a material portion of the suffix in these bounded shells."
        ),
    },
    "explicit_zero_zero_changing_loop": {
        "route_certificate_status": load("explicit_zero_detour_route")["status"],
        "promotion": load("explicit_zero_detour_route")["promotion"],
        "q_sequence_from_A11": load("explicit_zero_detour_route")["full_route_q_sequence_from_A11"],
        "zero_changing_loop": load("explicit_zero_detour_route")["zero_changing_loop"],
        "endpoint": load("explicit_zero_detour_route")["endpoint"],
        "bounded_correct_scale_continuations": {
            "q4_orbit32_pinned": {
                q: data["primitive_nef_fixed_norm_candidates"]
                for q, data in load("orbit32_pinned_targeted")["searches"].items()
            },
            "q6_orbit3372_pinned": {
                q: data["primitive_nef_fixed_norm_candidates"]
                for q, data in load("orbit3372_pinned_targeted")["searches"].items()
            },
        },
        "conclusion": (
            "The q6 orbit3372 branch has an exact degree-two q6 return to the original "
            "2A5 fibration with a new zero, followed by an exact q12 edge to current 3A3. "
            "Continuing through marked suffix fibres gives a full determinant-one route to "
            "pinned R17. The three-edge loop passes exact nef gates while the direct q104 "
            "presentation rejects two named explicit curves, but its raw score 14641 is "
            "still above the direct comparator's 13518, so it is not promoted."
        ),
    },
    "promoted_q6_orbit1307_route": {
        "route_certificate_status": load("promoted_route")["status"],
        "promotion": load("promoted_route")["promotion"],
        "new_splice": load("promoted_route")["new_splice"],
        "landing_current_3A3_identification": load("promoted_route")["landing_current_3A3_identification"],
        "endpoint": load("promoted_route")["endpoint"],
        "full_route_q_sequence_from_A11": load("promoted_route")["full_route_q_sequence_from_A11"],
        "exhaustive_search": {
            "status": load("zero_loop_search")["status"],
            "counts": load("zero_loop_search")["counts"],
            "strict_cost_winner_count": load("zero_loop_search")["strict_cost_winner_count"],
            "best": load("zero_loop_search")["ranked_loops"][0],
        },
        "lifting_instruction": (
            "After the current A11 q8/orbit12 lift reaches the equation-explicit 2A5 zero, "
            "switch to q6 orbit1307, use old_A11_component_10 as the child zero, take the "
            "exact q4 return to 2A5 and q6 exit to current 3A3, then resume the existing suffix."
        ),
    },
    "promoted_q4_orbit230_q6_orbit1315_route": {
        "route_certificate_status": q230_promotion["status"],
        "promotion": q230_promotion["promotion"],
        "a11_splice": q230_promotion["a11_splice"],
        "combined_bottleneck_comparison": q230_promotion["combined_bottleneck_comparison"],
        "landing_current_3A3_identification": q230_promotion["landing_current_3A3_identification"],
        "endpoint": q230_promotion["endpoint"],
        "full_route_q_sequence_from_A11": q230_promotion["full_route_q_sequence_from_A11"],
        "full_route_q_sequence_from_H3": q230_promotion["full_route_q_sequence_from_H3"],
        "second_zero_search": {
            "status": load("q230_second_zero_search")["status"],
            "parameters": load("q230_second_zero_search")["search_parameters"],
            "counts": load("q230_second_zero_search")["counts"],
            "best": load("q230_second_zero_search")["ranked_presentations"][0],
        },
        "lifting_instruction": (
            "After A11 q8/orbit12 reaches the equation-explicit 2A5 zero, switch to "
            "q4/orbit230 with old_A11_component_10, q4 return, q6/orbit1315 with "
            "old_A5A5_component_1, q4 return, q4 exit to current 3A3, then resume the suffix."
        ),
    },
    "promoted_first_q8_q4_orbit11_route": {
        "route_certificate_status": first_q8_promotion["status"],
        "promotion": first_q8_promotion["promotion"],
        "first_q8_replacement": first_q8_promotion["first_q8_replacement"],
        "equation_D13_identification": first_q8_promotion["equation_D13_identification"],
        "combined_equation_cost_comparison": first_q8_promotion["combined_equation_cost_comparison"],
        "endpoint": first_q8_promotion["endpoint"],
        "full_route_q_sequence_from_H3": first_q8_promotion["full_route_q_sequence_from_H3"],
        "first_loop_search": {
            "status": first_q8_search["status"],
            "counts": first_q8_search["counts"],
            "strict_winner_count": first_q8_search["strict_winner_count"],
            "best": first_q8_search["ranked_presentations"][0],
        },
        "second_loop_boundary": {
            "status": first_q8_second["status"],
            "counts": first_q8_second["counts"],
            "strict_winner_count": first_q8_second["strict_winner_count"],
            "best": first_q8_second["ranked_presentations"][0],
            "decision": "No q4/q6/q8 second changed-zero loop beats the direct 1,952-point D13 exit.",
        },
        "widened_first_edge_boundary": {
            "q10_degree2": {
                "status": load("first_q8_q10_search")["status"],
                "counts": load("first_q8_q10_search")["counts"],
                "best_operational_score": load("first_q8_q10_search")["ranked_presentations"][0]["inherited_explicit_equation_cost"]["operational_total_score"],
            },
            "q6_q9_q12_degree3": {
                "status": load("first_q8_degree3_search")["status"],
                "counts": load("first_q8_degree3_search")["counts"],
                "best_operational_score": load("first_q8_degree3_search")["ranked_presentations"][0]["inherited_explicit_equation_cost"]["operational_total_score"],
            },
            "q8_q12_q16_degree4": {
                "status": load("first_q8_degree4_search")["status"],
                "counts": load("first_q8_degree4_search")["counts"],
                "best_operational_score": load("first_q8_degree4_search")["ranked_presentations"][0]["inherited_explicit_equation_cost"]["operational_total_score"],
            },
            "decision": "No widened degree-two q10, degree-three q6/q9/q12, or degree-four q8/q12/q16 presentation beats 3,961.",
        },
        "cross_stage_D13_boundary": {
            "marking_status": load("first_q8_landing_d13_marking")["status"],
            "D12_old_fibre_degree": load("first_q8_landing_d13_marking")["target_fibres_in_root_adapted_hub"]["current_0_D12"][1],
            "search_status": load("first_q8_landing_d13_crossover")["status"],
            "counts": load("first_q8_landing_d13_crossover")["counts"],
            "direct_changed_zero_D13_to_D12_score": load("first_q8_landing_d13_crossover")["direct_q24"]["score"],
            "best_zero_loop_raw_credit_score": load("first_q8_landing_d13_crossover")["ranked_presentations"][0]["total_equation_cost_score"],
            "best_zero_loop_operational_score": load("first_q8_landing_d13_crossover")["ranked_presentations"][0]["inherited_explicit_equation_cost"]["operational_total_score"],
            "decision": "Do not carry the first-q8 landing zero into the D13 splice: q60 coefficient and RR growth make it more expensive than the canonical 25,323-point D13 continuation.",
        },
        "lifting_instruction": (
            "After the initial H3 q6 reaches E8+E6/MW3, switch to q4/orbit11 with "
            "old_E8E6_component_1 as zero, q4 return, and q4 exit to equation D13; "
            "then resume the already certified D13 and A11 optimized continuations."
        ),
    },
    "requested_lifting_agent_exchange": [
        "The actual zero section chosen by the equation compiler after orbit12 and its class in equation-A11 coordinates.",
        "Any additional measured resolved-RR dimensions or coefficient-growth data beyond the exact q24 bridge section now ingested.",
    ],
    "inputs": {
        "paths": {name: str(path.relative_to(ROOT)) for name, path in PATHS.items()},
        "sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in PATHS.items()},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"A11ROUTEHANDOFF|status={payload['status']}|promotion=1|output={OUTPUT.resolve()}")
