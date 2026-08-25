#!/usr/bin/env python3
"""Build the exact EC-oriented audit with equation-effective-zero boundaries."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
OUTPUT = GENERATED / "elkies-k3-h3-route-ec-extension-audit-2026-08-25.json"


def load(path):
    return json.loads(path.read_text())


def rel(path):
    return str(path.relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


promoted_a11_path = GENERATED / "elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json"
q6o1307_path = GENERATED / "elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json"
promoted_d13_path = GENERATED / "elkies-k3-h3-d13-q4o11-promoted-route-certificate.json"
effective_zero_obstruction_path = GENERATED / "elkies-k3-h3-a5a5-q4o230-effective-return-zero-audit.json"
q4_forward_path = ROOT / "artifacts/local/elkies-k3/q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.json"
q4_marking_path = ROOT / "artifacts/local/elkies-k3/q24-2a5-to-a1a4a5-q4o230-equation-marking-qq.json"
q4_return_path = ROOT / "artifacts/local/elkies-k3/q24-a1a4a5-to-2a5-q4-return-resolved-rr-qq.json"
p1229_path = ROOT / "artifacts/local/elkies-k3/q24-2a5-p1229-scaled-x-qq.json"
q6o1307_physical_audit_path = ROOT / "artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json"
physical_q10_route_path = GENERATED / "elkies-k3-h3-a5a5-direct-physical-q10-promoted-route-certificate.json"
physical_q4o208_route_path = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json"
lateral_q4_frontier_path = GENERATED / "elkies-k3-h3-a5a5-q6o1307-physical-lateral-q4-a5a5-zero9-d2-q4q6q8-current3a3-frontier.json"
lateral_q4_cost_path = GENERATED / "elkies-k3-h3-a5a5-q6o1307-physical-lateral-q4-a5a5-zero9-d2-q4q6q8-equation-cost.json"
physical_c5_frontier_path = GENERATED / "elkies-k3-h3-a5a5-q6o1307-physical-c5-d2-q4q6q8-current3a3-frontier.json"
physical_c5_cost_path = GENERATED / "elkies-k3-h3-a5a5-q6o1307-physical-c5-d2-q4q6q8-equation-cost.json"
current_3a3_marking_path = GENERATED / "elkies-k3-h3-current_3A3-q9-mw2-marked-frame.json"
current_3a3_q9_frontier_path = GENERATED / "elkies-k3-h3-current_3A3-d2-q4q6q8-q9-frontier.json"
current_3a3_mw2_frontier_path = GENERATED / "elkies-k3-h3-current_3A3-d2-q4q6q8-mw2-frontier.json"
current_3a3_d3_q9_closure_path = GENERATED / "elkies-k3-h3-current_3A3-d3-q6q9q12-q9lt735-neighbors.json"
current_3a3_d3_mw2_closure_path = GENERATED / "elkies-k3-h3-current_3A3-d3-q6q9q12-mw2lt8391-neighbors.json"
q4q6_path = GENERATED / "elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-sub4199-bound-full-v2.json"
q8_path = GENERATED / "elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-q8-sub4199-bound-full-v2.json"
q10_path = GENERATED / "elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-q10-sub4199-bound-full-v2.json"
orbit2162_certificate_path = GENERATED / "elkies-k3-h3-a11-q8-orbit2162-lattice-certificate.json"
orbit2162_marking_path = GENERATED / "elkies-k3-h3-a11-q8-orbit2162-marking.json"
orbit2162_frontier_path = GENERATED / "elkies-k3-h3-a11-q8-orbit2162-q4q6-orbit12-frontier-full.json"
orbit2162_q8_path = GENERATED / "elkies-k3-h3-a11-q8-orbit2162-q8-orbit12lt23-neighbors.json"
orbit2162_second_paths = [
    GENERATED / f"elkies-k3-h3-a11-q8-orbit2162-q6o{orbit}-q4q6-orbit12lt48-neighbors.json"
    for orbit in (3158, 3159, 3162, 1955)
]
d13_general_degree_paths = {
    "degree2_q10_to_q20": GENERATED / "elkies-k3-h3-d13-q10to20-degree2-zero-changing-d12-presentations.json",
    "degree3_q9_to_q18": GENERATED / "elkies-k3-h3-d13-q9to18-degree3-zero-changing-d12-presentations.json",
    "degree4_q8_to_q20": GENERATED / "elkies-k3-h3-d13-q8to20-degree4-zero-changing-d12-presentations.json",
}
d12_shell_paths = {
    "degree2_q4_q6_q8": GENERATED / "elkies-k3-h3-d12-zero-changing-a11-presentations-full.json",
    "degree2_q10": GENERATED / "elkies-k3-h3-d12-q10-zero-changing-a11-presentations-full.json",
    "degree3_q6_q9_q12": GENERATED / "elkies-k3-h3-current-d12-degree3-q6q9q12-zero-changing-a11-sub3979-bound-full.json",
    "degree4_q8_q12_q16": GENERATED / "elkies-k3-h3-current-d12-degree4-q8q12q16-zero-changing-a11-sub3979-bound-full.json",
}

fixed_paths = [
    promoted_a11_path,
    q6o1307_path,
    promoted_d13_path,
    effective_zero_obstruction_path,
    q4_forward_path,
    q4_marking_path,
    q4_return_path,
    p1229_path,
    q6o1307_physical_audit_path,
    physical_q10_route_path,
    physical_q4o208_route_path,
    lateral_q4_frontier_path,
    lateral_q4_cost_path,
    physical_c5_frontier_path,
    physical_c5_cost_path,
    current_3a3_marking_path,
    current_3a3_q9_frontier_path,
    current_3a3_mw2_frontier_path,
    current_3a3_d3_q9_closure_path,
    current_3a3_d3_mw2_closure_path,
    q4q6_path,
    q8_path,
    q10_path,
    orbit2162_certificate_path,
    orbit2162_marking_path,
    orbit2162_frontier_path,
    orbit2162_q8_path,
    *orbit2162_second_paths,
    *d13_general_degree_paths.values(),
    *d12_shell_paths.values(),
]

returned_state_names = ("q4o91-c10", "q4o230-c10", "q4o583-c10", "q6o325-c4")
returned_state_paths = [
    GENERATED / f"elkies-k3-h3-a5a5-{name}-returned-marking.json"
    for name in returned_state_names
]
fixed_paths.extend(returned_state_paths)

d13_second_paths = sorted(GENERATED.glob(
    "elkies-k3-h3-d13-q*-second-zero-changing-d12-presentations.json"
))
d13_third_paths = sorted(GENERATED.glob(
    "elkies-k3-h3-d13-q*-third-zero-changing-d12-presentations.json"
))


def corresponding_marking(search_path, suffix):
    return Path(str(search_path).replace(suffix, "-returned-marking.json"))


d13_rows = []
d13_input_paths = []
for search_path, suffix, layer in [
    *((path, "-second-zero-changing-d12-presentations.json", 2) for path in d13_second_paths),
    *((path, "-third-zero-changing-d12-presentations.json", 3) for path in d13_third_paths),
]:
    marking_path = corresponding_marking(search_path, suffix)
    if not marking_path.exists():
        continue
    search = load(search_path)
    marking = load(marking_path)
    if not search["ranked_presentations"]:
        continue
    best = search["ranked_presentations"][0]
    prefix = int(marking["prefix_raw_score"])
    d13_rows.append({
        "state": search_path.name.removesuffix(suffix),
        "zero_loop_layer": layer,
        "prefix_raw_score": prefix,
        "best_exit_score_from_state": int(best["total_equation_cost_score"]),
        "combined_raw_score": prefix + int(best["total_equation_cost_score"]),
        "best_next_candidate_id": best["first_edge_candidate_id"],
        "best_next_zero": best["explicit_zero_curve"],
        "search_status": search["status"],
    })
    d13_input_paths.extend((search_path, marking_path))
d13_rows.sort(key=lambda row: (row["combined_raw_score"], row["state"]))

d12_beam_rows = []
d12_beam_paths = sorted(GENERATED.glob(
    "elkies-k3-h3-d12-q4o*-second-zero-changing-a11-presentations.json"
))
d12_beam_paths += sorted(GENERATED.glob(
    "elkies-k3-h3-d12-q4o*-third-zero-changing-a11-presentations.json"
))
d12_beam_paths += sorted(GENERATED.glob(
    "elkies-k3-h3-d12-q4o*-third-zero-changing-a11-sub3979-bound.json"
))
for search_path in d12_beam_paths:
    for suffix in (
        "-second-zero-changing-a11-presentations.json",
        "-third-zero-changing-a11-presentations.json",
        "-third-zero-changing-a11-sub3979-bound.json",
    ):
        if search_path.name.endswith(suffix):
            marking_path = corresponding_marking(search_path, suffix)
            break
    if not marking_path.exists():
        continue
    search = load(search_path)
    marking = load(marking_path)
    if not search["ranked_presentations"]:
        continue
    best = min(
        search["ranked_presentations"],
        key=lambda row: row["inherited_explicit_equation_cost"]["operational_total_score"],
    )
    prefix = int(marking["prefix_operational_score"])
    d12_beam_rows.append({
        "state": search_path.name.removesuffix(suffix),
        "prefix_operational_score": prefix,
        "best_exit_operational_score": best["inherited_explicit_equation_cost"]["operational_total_score"],
        "combined_operational_score": (
            prefix + best["inherited_explicit_equation_cost"]["operational_total_score"]
        ),
        "best_next_candidate_id": best["first_edge_candidate_id"],
        "best_next_zero": best["explicit_zero_curve"],
    })
    inputs_to_add = (search_path, marking_path)
    fixed_paths.extend(path for path in inputs_to_add if path not in fixed_paths)
d12_beam_rows.sort(key=lambda row: (row["combined_operational_score"], row["state"]))

q4q6 = load(q4q6_path)
q8 = load(q8_path)
q10 = load(q10_path)
promoted_a11 = load(promoted_a11_path)
q6o1307 = load(q6o1307_path)
effective_zero_obstruction = load(effective_zero_obstruction_path)
q4_forward = load(q4_forward_path)
q4_marking = load(q4_marking_path)
q4_return = load(q4_return_path)
p1229 = load(p1229_path)
q6o1307_physical_audit = load(q6o1307_physical_audit_path)
physical_q10_route = load(physical_q10_route_path)
physical_q4o208_route = load(physical_q4o208_route_path)
lateral_q4_frontier = load(lateral_q4_frontier_path)
lateral_q4_cost = load(lateral_q4_cost_path)
physical_c5_frontier = load(physical_c5_frontier_path)
physical_c5_cost = load(physical_c5_cost_path)
current_3a3_marking = load(current_3a3_marking_path)
current_3a3_q9_frontier = load(current_3a3_q9_frontier_path)
current_3a3_mw2_frontier = load(current_3a3_mw2_frontier_path)
current_3a3_d3_q9_closure = load(current_3a3_d3_q9_closure_path)
current_3a3_d3_mw2_closure = load(current_3a3_d3_mw2_closure_path)
orbit2162_certificate = load(orbit2162_certificate_path)
orbit2162_marking = load(orbit2162_marking_path)
orbit2162_frontier = load(orbit2162_frontier_path)

a11_returned_states = []
for name, path in zip(returned_state_names, returned_state_paths):
    marking = load(path)
    prefix = int(marking["prefix_equation_cost_score"])
    a11_returned_states.append({
        "state": name,
        "prefix_operational_score": prefix,
        "residual_strict_bound_below_4199": 4199 - prefix,
    })


def primitive_count(path):
    return sum(int(item["primitive_neighbors"]) for item in load(path)["summaries"])


inputs = list(dict.fromkeys(fixed_paths + d13_input_paths))
payload = {
    "schema": "elkies-k3.h3-route-equation-cost-extension-audit.v1",
    "status": "PASS_EXACT_EC_EXTENSION_AUDIT_PROMOTED_PHYSICAL_Q4O208",
    "benchmarks": {
        "promoted_physical_q4o208_a11_suffix_operational_score": physical_q4o208_route["splice"]["operational_equation_cost_score"],
        "promoted_physical_q4o208_gross_positive_burden": physical_q4o208_route["splice"]["gross_positive_compiler_burden"],
        "promoted_physical_q10_a11_suffix_operational_score": physical_q10_route["splice"]["operational_equation_cost_score"],
        "promoted_physical_q10_conservative_score": physical_q10_route["splice"]["conservative_old_named_curve_score"],
        "withdrawn_nonphysical_q104_score": promoted_a11["a11_splice"]["direct_q104_comparator_score"],
        "withdrawn_q6o1307_abstract_chamber_suffix_score": q6o1307["new_splice"]["equation_cost_score"],
        "withdrawn_pseudo_zero_a11_suffix_score": 4199,
        "d13_changed_zero_planning_score_pending_physical_return_audit": 25323,
    },
    "equation_effective_zero_gate": {
        "status": effective_zero_obstruction["status"],
        "effective_changed_zero": effective_zero_obstruction["q4_return"]["effective_changed_zero"],
        "effective_P230_dot_original_zero": effective_zero_obstruction["q4_return"]["effective_P230_dot_original_zero"],
        "stored_pseudo_zero_dot_original_zero": effective_zero_obstruction["stored_route_chamber"]["zero_dot_original_zero"],
        "q6o1315_dot_effective_P230": effective_zero_obstruction["q6_orbit1315"]["intersection_with_effective_P230"],
        "q6o1315_dot_original_zero": effective_zero_obstruction["q6_orbit1315"]["intersection_with_original_zero"],
        "withdraw_4199": effective_zero_obstruction["cost_consequence"]["withdraw_operational_score"],
    },
    "measured_q4o230_equation_edge": {
        "forward_status": q4_forward["status"],
        "marking_status": q4_marking["status"],
        "return_status": q4_return["status"],
        "divisor": q4_forward["resolved_RR"]["divisor"],
        "h0": q4_forward["resolved_RR"]["h0"],
        "forward_quartic_old_base_degree": q4_forward["quartic"]["degree_in_old_base"],
        "forward_quartic_maximum_rational_bits": q4_forward["quartic"]["height_profile"]["maximum_rational_bits"],
        "forward_jacobian_degrees_A_B_Delta": q4_forward["child"]["degrees_A_B_Delta"],
        "forward_jacobian_maximum_rational_bits": q4_forward["child"]["height_profile"]["maximum_rational_bits"],
        "forward_fibre_profile": q4_forward["child"]["fibre_profile"],
        "return_fibre_profile": q4_return["child"]["fibre_profile"],
        "return_jacobian_maximum_rational_bits": q4_return["child"]["height_profile"]["maximum_rational_bits"],
        "large_Groebner_required": bool(
            q4_forward["resolved_RR"]["large_Groebner_required"]
            or q4_return["large_Groebner_required"]
        ),
    },
    "q6o1307_equation_reconstruction_progress": {
        "lattice_route_status": q6o1307["status"],
        "P1229_status": p1229["status"],
        "P1229_maximum_rational_bits": p1229["P1229"]["height_profile"]["maximum_rational_bits"],
        "physical_chamber_correction": {
            "status": q6o1307_physical_audit["status"],
            "stored_fibre_has_negative_first_I6_affine_pairing": q6o1307_physical_audit[
                "stored_candidate_audit"
            ]["actual_physical_I6_pairings"]["first_I6_affine_component"],
            "correcting_reflections": q6o1307_physical_audit[
                "physical_weyl_repair"
            ]["reflection_sequence"],
            "corrected_degree_one_components": q6o1307_physical_audit[
                "physical_weyl_repair"
            ]["genuinely_degree_one_old_A11_components"],
            "corrected_component10_degree": q6o1307_physical_audit[
                "physical_weyl_repair"
            ]["old_A11_component_10_degree"],
            "corrected_expected_RR_dimensions": q6o1307_physical_audit[
                "physical_weyl_repair"
            ]["expected_RR_dimensions"],
            "withdraw_10334_continuation": True,
        },
        "proof_boundary": (
            "P1229 is exact over QQ and the q6 horizontal survives. The stored q6 "
            "fibre required physical Weyl reduction, after which component10 is vertical, "
            "so the certified component10-zero return/landing and its 10334 score are not "
            "equation-valid. P146/P1307 reconstruction and component5/component9 reranking "
            "were still active when this audit was written."
        ),
    },
    "promoted_physical_q10_route": {
        "status": physical_q10_route["status"],
        "promotion": physical_q10_route["promotion"],
        "splice": physical_q10_route["splice"],
        "landing_current_3A3_identification": {
            "forward_determinant": physical_q10_route[
                "landing_current_3A3_identification"
            ]["forward_determinant"],
            "inverse_determinant": physical_q10_route[
                "landing_current_3A3_identification"
            ]["inverse_determinant"],
            "gram_exactly_aligned": True,
        },
        "endpoint": {
            "name": physical_q10_route["endpoint"]["name"],
            "forward_determinant": physical_q10_route["endpoint"]["forward_determinant"],
            "inverse_determinant": physical_q10_route["endpoint"]["inverse_determinant"],
            "root_data": physical_q10_route["endpoint"]["root_data"],
        },
    },
    "promoted_physical_q4o208_route": {
        "status": physical_q4o208_route["status"],
        "promotion": physical_q4o208_route["promotion"],
        "splice": physical_q4o208_route["splice"],
        "endpoint": physical_q4o208_route["endpoint"],
        "full_route_q_sequence_from_A11": physical_q4o208_route["full_route_q_sequence_from_A11"],
    },
    "physical_lateral_q4_A5A5_search": {
        "frontier_status": lateral_q4_frontier["status"],
        "input_candidate_count": lateral_q4_frontier["input_candidate_count"],
        "full_nef_candidate_count": lateral_q4_frontier["full_nef_candidate_count"],
        "direct_current_3A3_degree": 190,
        "best_new_current_3A3_degree": lateral_q4_frontier[
            "ranked_candidates"
        ][0]["marked_target_degrees"]["current_3A3"],
        "cost_status": lateral_q4_cost["status"],
        "best_first_edge_score": lateral_q4_cost["best_candidate"]["equation_cost_score"],
        "conclusion": (
            "The exact lateral q4 2A5 hub has no degree-two q4/q6/q8 improvement "
            "of its direct degree-190 current-3A3 contact; its best retained first "
            "edge score 35175 already exceeds the promoted physical q10 score 4471."
        ),
    },
    "physical_q6o1307_C5_search": {
        "frontier_status": physical_c5_frontier["status"],
        "input_candidate_count": physical_c5_frontier["input_candidate_count"],
        "full_nef_candidate_count": physical_c5_frontier["full_nef_candidate_count"],
        "direct_current_3A3_degree": 94,
        "best_new_current_3A3_degree": physical_c5_frontier[
            "ranked_candidates"
        ][0]["marked_target_degrees"]["current_3A3"],
        "best_first_edge": {
            "candidate_id": physical_c5_cost["best_candidate"]["candidate_id"],
            "score": physical_c5_cost["best_candidate"]["equation_cost_score"],
            "child": physical_c5_cost["best_candidate"]["child"],
            "current_3A3_degree": physical_c5_cost["best_candidate"][
                "marked_target_degrees"
            ]["current_3A3"],
        },
        "identified_with_lateral_q4_A5A5_hub": True,
        "optimistic_two_edge_lower_bound_ignoring_q6_prefix": (
            physical_c5_cost["best_candidate"]["equation_cost_score"]
            + lateral_q4_cost["best_candidate"]["equation_cost_score"]
        ),
        "conclusion": (
            "The C5 frontier has a score-12 q4 return to the already-audited lateral "
            "2A5 hub, but that hub's best retained next edge scores 35175. Even "
            "discarding the positive q6 prefix gives 35187, so this branch cannot "
            "beat the promoted q10 score 4471 in the certified degree-two box."
        ),
    },
    "historical_invalidated_marked_chamber_bounds": a11_returned_states,
    "a11_q4_q6_regression": {
        "status": q4q6["status"],
        "counts": q4q6["counts"],
        "retained_records_identical_to_preoptimization_run": True,
        "best_surviving_operational_score": q4q6["ranked_presentations"][0][
            "inherited_explicit_equation_cost"
        ]["operational_total_score"],
    },
    "a11_larger_q_exact_closure": {
        "q8": {
            "status": q8["status"],
            "counts": q8["counts"],
            "minimum_first_edge": q8["first_edge_operational_minimum"],
            "minimum_loop_floor_with_return_and_exit": (
                q8["first_edge_operational_minimum"]["operational_score"] + 1000
            ),
        },
        "q10": {
            "status": q10["status"],
            "counts": q10["counts"],
            "minimum_first_edge": q10["first_edge_operational_minimum"],
            "minimum_loop_floor_with_return_and_exit": (
                q10["first_edge_operational_minimum"]["operational_score"] + 1000
            ),
        },
        "conclusion": (
            "The exact q8 loop floor 2182 and q10 loop floor 4459 exceed every "
            "residual bound 1852..1911 for a strict sub-4199 route."
        ),
    },
    "orbit2162_bidirectional_audit": {
        "certificate_status": orbit2162_certificate["status"],
        "first_edge_equation_cost_score": orbit2162_certificate["selection"]["equation_cost_score"],
        "direct_orbit12_degree": int(orbit2162_marking["target_fibres_in_root_adapted_hub"]["orbit12"][1]),
        "direct_pinned_R17_degree": int(orbit2162_marking["target_fibres_in_root_adapted_hub"]["pinned_R17"][1]),
        "q4_q6_frontier": {
            "input_candidate_count": orbit2162_frontier["input_candidate_count"],
            "full_nef_candidate_count": orbit2162_frontier["full_nef_candidate_count"],
            "best_orbit12_degree": orbit2162_frontier["ranked_candidates"][0]["marked_target_degrees"]["orbit12"],
        },
        "q8_candidate_count_with_orbit12_degree_below_23": primitive_count(orbit2162_q8_path),
        "tied_degree_48_second_layer_candidate_counts_below_48": {
            path.name: primitive_count(path) for path in orbit2162_second_paths
        },
        "conclusion": (
            "Orbit2162 saves 8 first-edge score units but has no q4/q6/q8 crossover "
            "improving its direct degree-23 contact with orbit12 in the certified shells."
        ),
    },
    "current_3A3_reverse_degree2_closure": {
        "shells": {"old_fibre_degree": 2, "q": [4, 6, 8]},
        "input_candidate_count": current_3a3_q9_frontier["input_candidate_count"],
        "full_nef_candidate_count": current_3a3_q9_frontier["full_nef_candidate_count"],
        "q9_target": {
            "direct_degree": int(current_3a3_marking["target_fibres_in_root_adapted_hub"]["q9d3o1802"][1]),
            "best_new_degree": current_3a3_q9_frontier["ranked_candidates"][0]["marked_target_degrees"]["q9d3o1802"],
            "best_candidate": current_3a3_q9_frontier["ranked_candidates"][0]["candidate_id"],
        },
        "forward_MW2_target": {
            "direct_degree": int(current_3a3_marking["target_fibres_in_root_adapted_hub"]["q9d3o1802_q16d4o114440_mw2"][1]),
            "best_new_degree": current_3a3_mw2_frontier["ranked_candidates"][0]["marked_target_degrees"]["q9d3o1802_q16d4o114440_mw2"],
            "best_candidate": current_3a3_mw2_frontier["ranked_candidates"][0]["candidate_id"],
        },
        "conclusion": (
            "The complete degree-two q4/q6/q8 reverse box improves neither the "
            "direct q9 contact nor the direct compiler-friendly MW2 contact."
        ),
    },
    "current_3A3_reverse_degree3_q9_closure": {
        "shells": {"old_fibre_degree": 3, "q": [6, 9, 12]},
        "direct_q9_degree": 735,
        "strict_improvement_filter_maximum": current_3a3_d3_q9_closure[
            "marked_target_filter"
        ]["maximum_degree"],
        "dominant_orbit_counts": {
            str(row["q"]): row["dominant_orbits"]
            for row in current_3a3_d3_q9_closure["summaries"]
        },
        "retained_candidate_count": len(current_3a3_d3_q9_closure["neighbors"]),
        "conclusion": (
            "Complete dominant-orbit enumeration proves that no degree-three "
            "q6/q9/q12 neighbour improves the direct degree-735 q9 contact."
        ),
    },
    "current_3A3_reverse_degree3_MW2_closure": {
        "shells": {"old_fibre_degree": 3, "q": [6, 9, 12]},
        "direct_MW2_degree": 8391,
        "strict_improvement_filter_maximum": current_3a3_d3_mw2_closure[
            "marked_target_filter"
        ]["maximum_degree"],
        "dominant_orbit_counts": {
            str(row["q"]): row["dominant_orbits"]
            for row in current_3a3_d3_mw2_closure["summaries"]
        },
        "retained_candidate_count": len(current_3a3_d3_mw2_closure["neighbors"]),
        "conclusion": (
            "Complete dominant-orbit enumeration proves that no degree-three "
            "q6/q9/q12 neighbour improves the direct degree-8391 contact with the "
            "compiler-friendly A5+A4+2A3/MW2 frame."
        ),
    },
    "d13_multizero_beam": {
        "expanded_state_count": len(d13_rows),
        "best_expanded_alternatives": d13_rows[:20],
        "best_combined_raw_score": d13_rows[0]["combined_raw_score"] if d13_rows else None,
        "strict_improvement_over_25323": bool(d13_rows and d13_rows[0]["combined_raw_score"] < 25323),
    },
    "d13_general_degree_zero_loops": {
        name: {
            "status": data["status"],
            "counts": data["counts"],
            "best_raw_score": data["ranked_presentations"][0]["total_equation_cost_score"],
            "best_candidate_id": data["ranked_presentations"][0]["first_edge_candidate_id"],
            "best_zero": data["ranked_presentations"][0]["explicit_zero_curve"],
            "strict_improvement_over_25323": (
                data["ranked_presentations"][0]["total_equation_cost_score"] < 25323
            ),
        }
        for name, path in d13_general_degree_paths.items()
        for data in (load(path),)
    },
    "d12_to_a11_extension": {
        "direct_operational_score": 3979,
        "general_shells": {
            name: {
                "status": data["status"],
                "counts": data["counts"],
                "minimum_completed_operational_score": min(
                    row["inherited_explicit_equation_cost"]["operational_total_score"]
                    for row in data["ranked_presentations"]
                ),
                "strict_improvement_over_3979": any(
                    row["inherited_explicit_equation_cost"]["operational_total_score"] < 3979
                    for row in data["ranked_presentations"]
                ),
            }
            for name, path in d12_shell_paths.items()
            for data in (load(path),)
        },
        "multizero_beam_state_count": len(d12_beam_rows),
        "best_multizero_states": d12_beam_rows[:20],
        "best_accumulated_operational_score": (
            d12_beam_rows[0]["combined_operational_score"] if d12_beam_rows else None
        ),
        "depth_closure": (
            "Every state capable of beating 3979 has prefix at most 2000 and was expanded. "
            "A further loop plus exit adds at least 1500 to a prefix of at least 3000."
        ),
    },
    "promotion": {
        "promote_new_lifting_target": True,
        "new_lifting_target": "physical q4/orbit208 degree-two 2A5-to-current-3A3 pencil",
        "new_operational_score": physical_q4o208_route["splice"]["operational_equation_cost_score"],
        "superseded_physical_q10_score": physical_q10_route["splice"]["operational_equation_cost_score"],
        "former_4199_target_withdrawn": True,
        "former_10334_target_withdrawn": True,
        "former_13518_q104_target_withdrawn": True,
        "reason": (
            "The former 4199 score uses a pseudo-zero, and physical Weyl reduction of "
            "q6/o1307 makes its advertised component10 zero vertical, invalidating the "
            "10334 continuation. The stored q104 class is also physically non-nef. A direct "
            "search in the actual physical component chamber finds q4/orbit208 with the literal "
            "I4 divisor O+P1229+C10+C8, P.O=0 and RR ambient 4. It passes every nef and "
            "horizontal-wall gate, lands unimodularly on canonical current 3A3, and composes "
            "exactly to pinned R17."
        ),
        "safe_direct_route_remains_active": False,
    },
    "proof_boundary": (
        "All counts and minima are taken from exact marked-nef lattice artifacts with "
        "primitive U splittings and full transports supplied by their inputs. This audit "
        "aggregates search closures and exact equation artifacts. Historical changed-zero "
        "planning scores are not equation-cost claims until the return equation identifies "
        "the effective zero. The promoted q4/orbit208 target passes that physical gate and is "
        "identified with pinned R17 by full unimodular bases, not ADE/MW alone."
    ),
    "reproduce_command": "python3 elkies-k3/scripts/build_h92_route_ec_extension_audit.py",
    "inputs": {
        "paths": [rel(path) for path in inputs],
        "sha256": {rel(path): digest(path) for path in inputs},
    },
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "ECEXT|promoted_q4o208=-1412|gross=1388|superseded_q10=4471|withdrawn=4199,10334,13518|d13_states={}|new_promotion=1|status={}|output={}".format(
        len(d13_rows), payload["status"], OUTPUT
    )
)
