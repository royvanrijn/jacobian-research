#!/usr/bin/env python3
"""Verify the current scoped evidence and keep the full goal explicitly open."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    verification = json.loads((HERE / "verification.json").read_text())
    assert verification["status"] == "PASS"
    for path, digest in verification["input_sha256"].items():
        assert sha(HERE / path) == digest, f"verification input changed: {path}"
    source_manifest = json.loads((HERE / "sources/manifest.json").read_text())
    assert len(source_manifest) == 18
    for source in source_manifest:
        assert "unavailable" not in source
        assert sha(HERE / "sources" / source["file"]) == source["sha256"]
    catalogue = json.loads((HERE / "sources.json").read_text())
    assert len(catalogue["sources"]) == 46
    assert len({r["id"] for r in catalogue["sources"]}) == 46
    entries = {e["id"]: e for e in json.loads((ROOT / "MATH_STATUS.json").read_text())["entries"]}
    result = entries["MRC-AFFINE-SIGNED-LIFTS-20261116"]
    objective = entries["OP-MRC"]
    assert result["state"] == "proved" and objective["state"] == "open"
    assert result["artifact_hash"] == "sha256:" + sha(HERE / "verify_signed_lifts.py")
    assert not result["closes_problems"] and not result["formal_verification"] and not result["external_review"]
    extension = entries["MRC-LIFT-DIMENSION-RADICAL-20261116"]
    assert extension["state"] == "proved" and extension["id"] in objective["dependencies"]
    assert extension["artifact_hash"] == "sha256:" + sha(HERE / "verify_lift_dimension.py")
    dimension = json.loads((HERE / "lift_dimension_verification.json").read_text())
    assert dimension["status"] == "PASS" and dimension["new_sharp_panel_certificates"] == 0
    for path, digest in dimension["input_sha256"].items():
        assert sha(HERE / path) == digest, f"lift-dimension input changed: {path}"
    assert len(dimension["frozen_panel_independent_replay"]) == 16
    for source in json.loads((HERE / "sources/wirt-hm/manifest.json").read_text())["files"]:
        assert sha(HERE / "sources/wirt-hm" / source["path"]) == source["sha256"]
    outputs = {p.name: sha(p) for p in sorted((HERE / "output/pdf").glob("*.pdf"))}
    assert set(outputs) == {"PROOF.pdf", "LITERATURE.pdf", "LIFT_DIMENSION_AND_RADICAL_TOWERS.pdf", "ALL_CHARACTERISTIC_LIFTS.pdf", "BRANCHED_COVER_LIFT_BOUNDS.pdf", "ICOSAHEDRAL_LIFTS_AND_QUADRATIC_COEFFICIENTS.pdf", "JOINT_SIMPLE_QUOTIENT_KERNELS.pdf", "AUGMENTATION_COVER_HOMOLOGY.pdf", "SIX_POINT_COVERS_AND_FORMAL_ORTHOGONALITY.pdf", "TORUS_WHITEHEAD_PERIPHERAL_CURVES.pdf", "MERIDIAN_DEFORMATIONS_AND_STABILIZATION.pdf", "SIMPLE_MODULE_BOUNDS_AND_A5_EXHAUSTION.pdf", "GL3_QUOTIENTS_AND_SINGULAR_AFFINE_CERTIFICATES.pdf", "RANK_TWO_AFFINE_EQUIVALENCE.pdf", "TREFOIL_AMALGAMS_SHARP_FOUR.pdf", "PERIODIC_TWIST_FAMILIES.pdf", "SIMPLE_COVER_RANK_ADDITIVITY.pdf", "PERMUTATION_INDEX_GATE.pdf", "TERNARY_AFFINE_AND_COVER_CERTIFICATES.pdf", "PROJECTIVE_THIRTEEN_COVER_HOMOLOGY.pdf", "MUTATION_COMMON_DOUBLE_COVERS.pdf", "COMPACT_DOUBLE_COVER_PRESENTATIONS.pdf", "FOUR_PIECE_DOUBLE_COVERS_SHARP_FOUR.pdf", "THREE_PIECE_DOUBLE_COVERS_SHARP_FOUR.pdf", "FINAL_FIVE_DOUBLE_COVERS_SHARP_FOUR.pdf", "ORIGINAL_SIXTEEN_KNOWN_GRAPH_COVER_FAMILY.pdf"}
    integer = json.loads((HERE / "integer_lift_verification.json").read_text())
    current13 = json.loads((HERE / "current13_certificate_verification.json").read_text())
    all_field = entries["MRC-ALL-CHARACTERISTIC-LIFTS-20261116"]
    assert all_field["state"] == "proved" and all_field["id"] in objective["dependencies"]
    assert all_field["artifact_hash"] == "sha256:" + sha(HERE / "verify_integer_lifts.py")
    assert integer["status"] == current13["status"] == "PASS"
    assert integer["sharp_five_certificates"] == 16
    assert integer["independently_replayed_lift_matrices"] == 512
    assert current13["certificates"] == 243
    assert sha(HERE / "sharp_five_certificates.json") == integer["sharp_certificates_sha256"]
    for receipt in [integer, current13]:
        for path, digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"all-field input changed: {path}"
    current_manifest = json.loads((HERE / "sources/current13/manifest.json").read_text())
    for item in current_manifest["files"]:
        assert sha(HERE / "sources/current13" / item["path"]) == item["sha256"]
    diagram = json.loads((HERE / "diagram_audit_verification.json").read_text())
    names = json.loads((HERE / "current13_name_audit.json").read_text())
    snapshots = json.loads((HERE / "knotinfo_snapshot_audit.json").read_text())
    parabolic = json.loads((HERE / "parabolic_panel_verification.json").read_text())
    assert diagram["status"] == parabolic["status"] == "PASS"
    assert diagram["named_certificate_transfers"] == diagram["independent_Montesinos_constructions"] == 16
    assert names["status"] == "COMPLETE" and len(names["records"]) == 243 and not names["failures"]
    assert parabolic["targets"] == 16 and parabolic["sharp_four_certificates"] == 0
    for receipt in [diagram, names, snapshots, parabolic]:
        for path, digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"identity/parabolic input changed: {path}"
    for result in snapshots["results"]:
        assert result["author536_exactly_match_uncertain_set"] and result["author243_PD_exact_matches"] == 243
        assert sha(HERE / result["subset"]) == result["subset_sha256"]
    from verify_diagram_audit import edge_pairs
    author = {c["name"]: c for c in json.loads((HERE / "sources/current13/data/certificates.json").read_text())}
    for result in names["records"]:
        a, b = edge_pairs(author[result["name"]]["pd"]), edge_pairs(result["pd"])
        witness = result["source_to_census"]
        images, direction, flip = witness["dart_images"], witness["orientation"], witness["over_under_swap"]
        assert sorted(images) == list(range(len(a))) and direction in [-1, 1] and flip in [0, 1]
        for d in a:
            image = images[d]
            assert images[a[d]] == b[image]
            assert images[4*(d//4)+(d+1)%4] == 4*(image//4)+(image+direction)%4
            assert (d % 2) ^ (image % 2) == flip
    branched = json.loads((HERE / "branched_cover_verification.json").read_text())
    core = json.loads((HERE / "nonabelian_core_verification.json").read_text())
    subgroup = json.loads((HERE / "core_subgroup_homology_verification.json").read_text())
    cover_claim = entries["MRC-BRANCHED-COVER-LIFTS-20261116"]
    assert cover_claim["state"] == "proved" and cover_claim["id"] in objective["dependencies"]
    assert cover_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_branched_cover_bounds.py")
    assert not cover_claim["closes_problems"] and not cover_claim["formal_verification"] and not cover_claim["external_review"]
    assert branched["status"] == core["status"] == subgroup["status"] == "PASS"
    assert len(branched["all_degree_field_panel"]) == core["targets"] == 16
    assert len(branched["homology_controls"]) == 9
    for receipt in [branched, core, subgroup]:
        for path, digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"branched/core input changed: {path}"
    for case in branched["all_degree_field_panel"]:
        assert [d["degree"] for d in case["degrees"]] == [2,3,4,5]
        assert case["degrees"][0]["maximum_lift_dimension"] == 2
        assert all(d["surjective_bases"] == 0 for d in case["degrees"][1:])
    assert sorted(r["normalized_A5_core_colourings"] for r in core["results"]) == [1]*4+[121]*8+[241]*4
    assert not any(r["sharp_four"] for r in core["results"])
    assert subgroup["covers"] == len(subgroup["results"]) == 12 and len(subgroup["skipped_trivial"]) == 4
    determinants = {c["name"]:c["degrees"][0]["H1_torsion_factors"][0] for c in branched["all_degree_field_panel"]}
    for case in subgroup["results"]:
        assert case["A5_image_order"] == 60 and case["degree"] == 5
        assert case["H1_free_rank"] == 0 and case["H1_torsion_factors"] == [3*determinants[case["name"]]]
        assert case["maximum_H1_dimension_all_fields"] == 1 and not case["sharp_four"]
    assert len(subgroup["controls"]) == 2 and subgroup["invalid_relator_rejected"]
    quadratic = json.loads((HERE / "icosahedral_lift_verification.json").read_text())
    quadratic_claim = entries["MRC-QUADRATIC-COEFFICIENT-LIFTS-20261116"]
    assert quadratic_claim["state"] == "proved" and quadratic_claim["id"] in objective["dependencies"]
    assert quadratic_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_icosahedral_lifts.py")
    assert not quadratic_claim["closes_problems"] and not quadratic_claim["formal_verification"] and not quadratic_claim["external_review"]
    assert quadratic["status"] == "PASS" and quadratic["fixed_H3_bases"] == 12
    for path,digest in quadratic["input_sha256"].items():
        assert sha(HERE / path) == digest, f"quadratic input changed: {path}"
    assert quadratic["integral_split_factorizations"] == 12
    assert quadratic["every_field_fixed_base_affine_image_order_at_most"] == 120
    assert quadratic["characteristic2_no_fixed_point_control_orbits"] == [12,12,20,20]
    assert len(quadratic["quadratic_field_controls"]) == 4
    for case in quadratic["results"]:
        assert case["expanded_Smith_factors"] == [1]*20
        assert case["generic_dimension"] == case["characteristic5_dimension"] == case["maximum_all_fields"] == 3
        assert not case["sharp_four"]
    joint = json.loads((HERE / "joint_a5_verification.json").read_text())
    joint_claim = entries["MRC-JOINT-SIMPLE-QUOTIENT-KERNELS-20261116"]
    assert joint_claim["state"] == "proved" and joint_claim["id"] in objective["dependencies"]
    assert joint_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_joint_a5.py")
    assert not joint_claim["closes_problems"] and not joint_claim["formal_verification"] and not joint_claim["external_review"]
    assert joint["status"] == "PASS" and joint["cases"] == 293
    for path,digest in joint["input_sha256"].items():
        assert sha(HERE / path) == digest, f"joint-kernel input changed: {path}"
    assert joint["generating_pairs"] == 2280 and joint["automorphisms"] == 120
    assert joint["pair_kernel_capacity"] == 19 and joint["triple_kernel_capacity"] == 1668
    assert joint["same_knot_kernel_histogram"] == {"0":151,"1":118,"2":22,"3":1,"5":1}
    assert len(joint["results"]) == 293 and not any(r["sharp_four"] for r in joint["results"])
    assert joint["three_trefoil_sum_control_kernels"] == 49 and joint["two_trefoil_sum_control_kernels"] == 1
    assert len(joint["orbifold_patterns"]) == 40 and joint["orbifold_patterns_reaching_twenty"] == 39
    pretzel = joint["known_pretzel_control"]
    assert pretzel["A5_quotient_kernels"] == 23 and pretzel["Hom_to_A5"] == 2901
    assert pretzel["raw_lower"] == 3 and pretzel["joint_lower"] == pretzel["bridge_upper"] == 4
    geometry_probe = json.loads((HERE / "current13_montesinos_results.json").read_text())
    assert geometry_probe["status"] == "COMPLETE" and geometry_probe["completed_candidates"] == 30
    assert geometry_probe["diagrams_checked"] == 240 and not geometry_probe["matches"]
    for path,digest in geometry_probe["input_sha256"].items():
        assert sha(HERE / path) == digest, f"geometry-probe input changed: {path}"
    augmentation = json.loads((HERE / "augmentation_cover_verification.json").read_text())
    whitehead = json.loads((HERE / "whitehead_dependency_audit.json").read_text())
    augmentation_claim = entries["MRC-AUGMENTATION-COVER-HOMOLOGY-20261116"]
    assert augmentation_claim["state"] == "proved" and augmentation_claim["id"] in objective["dependencies"]
    assert augmentation_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_augmentation_covers.py")
    assert not augmentation_claim["closes_problems"] and not augmentation_claim["formal_verification"] and not augmentation_claim["external_review"]
    assert augmentation["status"] == "PASS" and augmentation["knots"] == augmentation["A5_kernel_classes"] == 16
    assert augmentation["knots_with_quotients"] == 12 and augmentation["all_index_five_subgroups_covered"]
    assert augmentation["maximum_gain_all_primes_except5"] == 1 and augmentation["new_sharp_cases"] == 0
    assert len(augmentation["skipped_no_index_five_subgroup"]) == 4
    assert [c["maximum_gain_all_primes_except5"] for c in augmentation["controls"]] == [4,2,5]
    assert augmentation["characteristic5_countercontrol"] == {"V4_invariant_dimension":1,"H1_G_V4":5,"b1_H_minus_b1_G":4}
    assert whitehead["status"] == "PASS scoped source and dependency audit"
    assert whitehead["all16_targets_still_in_open_list"] and not whitehead["own_proved_entries_using_DUT22"]
    assert whitehead["Dutra_original_v1_gap"]["common_nontrivial_reduced_word"] == [2,-1,-2,1]
    for receipt in [augmentation,whitehead]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"augmentation/source input changed: {path}"
    six = json.loads((HERE / "six_point_cover_verification.json").read_text())
    generator_probe = json.loads((HERE / "core_generator_probe_verification.json").read_text())
    formal_claim = entries["MRC-FORMAL-ORTHOGONALITY-20261116"]
    assert formal_claim["state"] == "proved" and formal_claim["id"] in objective["dependencies"]
    assert formal_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_six_point_covers.py")
    assert not formal_claim["closes_problems"] and not formal_claim["formal_verification"] and not formal_claim["external_review"]
    assert six["status"] == "PASS" and six["covers"] == len(six["results"]) == 16
    assert six["maximum_gain"] == 1 and six["new_sharp_cases"] == 0
    assert not six["orbifold_nonprincipal_classes"] and six["V4_V5_cover_cocycle_rank_bound_at_most2_all_fields"]
    for multiplier in [4,28]:
        assert sum(r["H1_H_torsion_factors"] == [multiplier*determinants[r["name"]]] for r in six["results"]) == 8
    assert [r["gain"] for r in six["controls"]] == [5,3,7]
    for row in six["results"]:
        assert all(c["orbifold_H1_V5"] == 0 for c in row["field_checks"])
        assert [(c["invariants"],c["H1_dimension"],c["Z1_dimension"]) for c in row["V5_nonsplit_characteristics"]] == [(1,2,6),(1,1,5)]
        assert row["V4_characteristic5"]["H1_dimension"] == row["V4_characteristic5"]["invariants"] == 1
    assert generator_probe["status"].startswith("PASS") and generator_probe["cases"] == 16
    assert generator_probe["bounded_candidates_per_case"] == 444 and generator_probe["actual_two_generator_certificates"] == 0
    assert not generator_probe["rank_conclusions"]
    for receipt in [six,generator_probe]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"six-point/formal input changed: {path}"
    torus = json.loads((HERE / "torus_whitehead_verification.json").read_text())
    torus_claim = entries["MRC-TORUS-WHITEHEAD-CURVES-20261116"]
    assert torus_claim["state"] == "proved" and torus_claim["id"] in objective["dependencies"]
    assert torus_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_torus_whitehead.py")
    assert not torus_claim["formal_verification"] and not torus_claim["external_review"]
    assert not torus_claim["independent_replay"] and not torus_claim["closes_problems"]
    assert torus["status"] == "PASS" and len(torus["torus_finite_controls"]) == 13
    assert sum(len(r["samples"]) for r in torus["torus_finite_controls"]) == 39
    assert torus["symbolic_T35"]["specialization_group_order"] == 120
    assert torus["finite_whitehead"]["triples_examined"] == 2976
    assert [r["valid_blocks"] for r in torus["finite_whitehead"]["valid"]] == [0,4,64,88]
    assert all(torus["symbolic_whitehead"][k] for k in ["boundary_commutator", "longitude_word_eigenvalue", "rank_one_meridian"])
    for path,digest in torus["input_sha256"].items():
        assert sha(HERE / path) == digest, f"torus/Whitehead input changed: {path}"
    tangent = json.loads((HERE / "meridian_tangent_verification.json").read_text())
    jets = json.loads((HERE / "meridian_jet_verification.json").read_text())
    deformation_claim = entries["MRC-MERIDIAN-DEFORMATIONS-STABILIZATION-20261116"]
    assert deformation_claim["state"] == "proved" and deformation_claim["id"] in objective["dependencies"]
    assert deformation_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_meridian_tangents.py")
    assert not deformation_claim["formal_verification"] and not deformation_claim["external_review"]
    assert not deformation_claim["closes_problems"]
    assert tangent["status"] == jets["status"] == "PASS"
    assert tangent["tangent_bases"] == len(tangent["tangents"]) == 16 and tangent["restricted_modules"] == 32
    assert len({r["name"] for r in tangent["tangents"]}) == 12
    for row in tangent["tangents"]:
        assert (row["normalized_tangent_dimension"],row["residual_conjugation_dimension"],row["modulo_conjugation_dimension"]) == (5,4,1)
        assert row["nonzero_meridian_eigenvalue_derivative"]
    for row in tangent["restricted_cocycles"]:
        assert len(row["modules"]) == 2
        assert all((m["restricted_cocycle_dimension"],m["principal_dimension"],m["restricted_H1"]) == (3,3,0) for m in row["modules"])
    assert jets["jets"] == 32 and jets["bases"] == 16 and jets["prime"] == 11 and jets["maximum_degree"] == 3
    for receipt in [tangent,jets]+[json.loads((HERE/p).read_text()) for p in ["meridian_tangent_results.json","stabilization_results.json","meridian_jet_results.json"]]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"deformation/stabilization input changed: {path}"
    all_a5 = json.loads((HERE / "a5_all_modules_verification.json").read_text())
    all_a5_claim = entries["MRC-SIMPLE-MODULE-A5-EXHAUSTION-20261116"]
    assert all_a5_claim["state"] == "proved" and all_a5_claim["id"] in objective["dependencies"]
    assert all_a5_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_a5_all_modules.py")
    assert not all_a5_claim["formal_verification"] and not all_a5_claim["external_review"] and not all_a5_claim["closes_problems"]
    assert all_a5["status"] == "PASS" and all_a5["kernel_classes"] == all_a5["new_W3_matrices"] == 16
    assert all_a5["new_binary_simple_actions"] == 32 and all_a5["new_sharp_cases"] == 0
    assert all_a5["all_A5_module_ratio_bound"] == "3/2" and all_a5["all_A5_ordinary_cocycle_rank_bounds_at_most"] == 2
    assert sorted(r["W3_integer_nonunit_factors"] for r in all_a5["results"]) == [[5]]*8+[[7,35]]*8
    assert len(all_a5["tensor_pairs"]) == 16 and set(all_a5["joint_image_orders"].values()) == {3600}
    assert all(r["Z1_dimension"] == 6 and r["H1_dimension"] == 2 for r in all_a5["tensor_pairs"])
    assert [r["rank_lower_bound"] for r in all_a5["controls"]] == [2,3]
    for receipt in [all_a5]+[json.loads((HERE/p).read_text()) for p in ["a5_all_modules_results.json","a5_tensor_pair_results.json"]]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"complete A5-module input changed: {path}"
    gl3 = json.loads((HERE / "gl3f2_verification.json").read_text())
    gl3_claim = entries["MRC-GL3-SINGULAR-AFFINE-CERTIFICATES-20261116"]
    assert gl3_claim["state"] == "proved" and gl3_claim["id"] in objective["dependencies"]
    assert gl3_claim["artifact_hash"] == "sha256:" + sha(HERE / "verify_gl3f2.py")
    assert not gl3_claim["formal_verification"] and not gl3_claim["external_review"] and not gl3_claim["closes_problems"]
    assert gl3["status"] == "PASS" and gl3["targets"] == len(gl3["results"]) == 293
    assert gl3["pair_capacity"] == 57 and gl3["triple_capacity"] == 13368
    assert gl3["kernel_histogram"] == {"0":139,"1":123,"2":29,"3":1,"7":1}
    assert len(gl3["module_checks"]) == len(gl3["reflection_checks"]) == 191
    assert gl3["binary_root_checks"] == 58 and gl3["deck_type_histogram"] == {"outer":162,"inner":29}
    assert gl3["module_Z1_histograms"] == {"natural":{"4":179,"5":12},"dual":{"4":179,"5":12},"Steinberg":{"8":165,"10":26}}
    assert gl3["absolute_simple_matrix_algebra_dimensions"] == {"natural":9,"dual":9,"Steinberg":64}
    assert len(gl3["characteristic2_regular_classes"]) == 4 and gl3["new_sharp_cases"] == 0
    assert gl3["T47_control"]["affine_group_order"] == 1344
    assert gl3["T47_control"]["global_fixed_codimension"] == 4
    assert gl3["T47_control"]["seed_pairing_rank"] == 3
    for receipt in [gl3]+[json.loads((HERE/p).read_text()) for p in ["gl3f2_results.json","outer_reflection7_results.json","gl3f2_full_panel_results.json","reflection7_full_panel_results.json","gl3_source_scope_audit.json"]]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE / path) == digest, f"GL3 source or result changed: {path}"
    gl3_catalogue = json.loads((HERE/"gl3f2_catalogue.json").read_text())
    assert gl3_catalogue["model_sha256"] == sha(HERE/"gl3f2.py")
    for row in json.loads((HERE/"gl3_source_scope_audit.json").read_text())["files"]:
        assert sha(HERE/row["path"]) == row["sha256"]
    rank_two = json.loads((HERE / "rank_two_equivalence_verification.json").read_text())
    rank_two_claim = entries["MRC-RANK-TWO-AFFINE-EQUIVALENCE-20261117"]
    assert rank_two_claim["state"] == "proved" and rank_two_claim["id"] in objective["dependencies"]
    assert rank_two_claim["artifact_hash"] == "sha256:" + sha(HERE/"verify_rank_two_equivalence.py")
    assert not rank_two_claim["formal_verification"] and not rank_two_claim["external_review"] and not rank_two_claim["closes_problems"]
    assert rank_two["status"] == "PASS" and rank_two["targets"] == len(rank_two["results"]) == 293
    assert rank_two["seed_orbits"] == {"2":1,"3":13,"4":39,"5":25}
    assert rank_two["deck_histogram"] == {"inner":170}
    assert rank_two["symmetric_base_counts"] == {"2":293,"3":26,"4":21,"5":0}
    assert rank_two["S3_binary_lift_dimensions"] == {"2":5,"3":21}
    assert rank_two["original_common_form_no_representation_targets"] == 272
    assert rank_two["original_common_form_maximum_bound"] == 3 and rank_two["new_open_cases"] == 0
    assert [r["image_order"] for r in rank_two["positive_known_pretzel_controls"]] == [24,192,1536]
    assert rank_two["binary3_linear_image_classification"]["eligible_image_subgroups_by_order"] == {"24":14,"6":28,"2":21,"168":1}
    assert rank_two["binary3_linear_image_classification"]["all_binary3_root_bounds_at_most"] == 3
    nonlinear = json.loads((HERE/"three_sheet_nonabelian_verification.json").read_text())
    assert nonlinear["status"] == "PASS" and nonlinear["covers"] == 26
    assert nonlinear["cyclic_orders"] == {"2":17,"4":1,"1":1} and nonlinear["remaining_two_generator_groups"] == 7
    assert nonlinear["coset_followup"] == {"CAPPED_INCOMPLETE":4,"CYCLIC_IMAGE_ONLY":3} and nonlinear["new_sharp_cases"] == 0
    for receipt in [rank_two,nonlinear]+[json.loads((HERE/p).read_text()) for p in ["small_symmetric_full_panel_results.json","rank_two_control_diagrams.json","three_sheet_nonabelian_results.json","three_sheet_coset_results.json"]]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path) == digest, f"rank-two/symmetric/nonlinear input changed: {path}"
    splices = json.loads((HERE/"three_sheet_splice_verification.json").read_text())
    splice_claim = entries["MRC-TREFOIL-AMALGAM-SHARP-FOUR-20261117"]
    assert splice_claim["state"] == "proved" and splice_claim["id"] in objective["dependencies"]
    assert splice_claim["artifact_hash"] == "sha256:" + sha(HERE/"verify_three_sheet_splices.py")
    assert not splice_claim["formal_verification"] and not splice_claim["external_review"] and not splice_claim["closes_problems"]
    assert splices["status"] == "PASS exact sharp-four evidence" and splices["sharp_named_knots"] == 4
    assert [r["name"] for r in splices["results"]] == ["13a_3145","13a_3147","13a_3148","13a_3149"]
    assert [(r["model"],r["epsilon"]) for r in splices["results"]] == [("A",1),("A",-1),("B",1),("B",-1)]
    assert all(r["cover_group_rank"]==2 and r["meridional_lower"]==r["bridge_upper"]==r["mr_equals_bridge"]==4 for r in splices["results"])
    assert splices["proved_cyclic_covers"]==22 and splices["proved_rank2_covers"]==4 and splices["unresolved_cover_groups"]==0
    assert all(r["abelian_cover"]["universal_abelian_cover_H1"]=="0" and r["abelian_cover"]["all_solvable_quotients"]=="cyclic" for r in splices["results"])
    structure = json.loads((HERE/"three_sheet_structure_verification.json").read_text())
    assert structure["status"]=="PASS" and structure["new_proved_cyclic"]==3
    provenance = json.loads((HERE/"splice_knot_provenance.json").read_text())
    assert provenance["status"]=="PASS exact named PD provenance" and len(provenance["results"])==4
    assert all(r["literal_PD_matches_after_adding_one"] and r["in_pinned_and_fetched_project_open_list"] and r["absent_from243_deposited_certificates"] for r in provenance["results"])
    four = json.loads((HERE/"four_sheet_group_verification.json").read_text())
    assert four["status"]=="PASS exact four-sheet group triviality" and four["covers"]==21 and four["all_groups_trivial"]
    for receipt in [splices,structure,provenance,four]+[json.loads((HERE/p).read_text()) for p in ["three_sheet_structural_results.json","three_sheet_metabelian_results.json","three_sheet_splice_results.json","four_sheet_group_results.json"]]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"trefoil-amalgam proof input changed: {path}"
    finite = json.loads((HERE/"finite_cover_quotient_verification.json").read_text())
    finite_claim = entries["MRC-FINITE-COVER-COMMUTATOR-CERTIFICATES-20261117"]
    assert finite_claim["state"]=="proved" and finite_claim["id"] in objective["dependencies"]
    assert finite_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_finite_cover_quotients.py")
    assert not finite_claim["formal_verification"] and not finite_claim["external_review"] and not finite_claim["closes_problems"]
    assert finite["status"]=="PASS independent finite commutator certificates"
    assert finite["finite_certificate_count"]==4 and finite["degree"]==25 and finite["rejected_corruptions"]==4
    assert finite["new_named_cases_beyond_retained_four"]==0
    assert finite["vertex_image_order"]==300 and finite["peripheral_translation_group_order"]==25
    assert finite["exact_cyclic_controls"]==[["13a_3291",2],["13a_3346",2],["13n_4637",4]]
    assert [r["name"] for r in finite["cases"]]==[r["name"] for r in splices["results"]]
    assert all(r["original_two_generator_relators_checked"]==6 and r["original_full_cover_relators_checked"]==39 and r["full_cover_generators"]==11 and r["all_relators_fix_every_point"] for r in finite["cases"])
    for receipt in [finite,json.loads((HERE/"finite_cover_full_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"finite commutator input changed: {path}"
    periodic = json.loads((HERE/"periodic_twist_family_verification.json").read_text())
    periodic_source = json.loads((HERE/"periodic_twist_source_audit.json").read_text())
    periodic_claim = entries["MRC-PERIODIC-TWIST-FAMILIES-20261117"]
    assert periodic_claim["state"]=="proved" and periodic_claim["id"] in objective["dependencies"]
    assert periodic_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_periodic_twist_families.py")
    assert not periodic_claim["formal_verification"] and not periodic_claim["external_review"] and not periodic_claim["closes_problems"]
    assert periodic["status"]=="PASS independent periodic twist family inputs"
    assert periodic["family_count"]==4 and periodic["explicit_expanded_diagrams"]==12 and periodic["rejected_corruptions"]==4
    assert (periodic["knot_action_degree"],periodic["blocks"],periodic["block_size"])==(75,3,25)
    assert [r["name"] for r in periodic["families"]]==[r["name"] for r in finite["cases"]]
    assert [(r["determinant_constant"],r["determinant_slope_in_N"]) for r in periodic["families"]]==[(543,1776),(537,1680),(543,1776),(537,1680)]
    assert all(r["period_added_crossings"]==6 and r["meridional_lower"]==r["bridge_upper"]==4 for r in periodic["families"])
    assert periodic_source["status"]=="PASS scoped primary source records" and periodic_source["catalogue_distinct_ids"]==len(catalogue["sources"])
    assert [r["id"] for r in periodic_source["records"]]==["STOIM07","JP25"]
    for receipt in [periodic,periodic_source,json.loads((HERE/"periodic_twist_family_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"periodic-twist input changed: {path}"
    sums=json.loads((HERE/"cover_sum_verification.json").read_text())
    sum_sources=json.loads((HERE/"cover_sum_source_audit.json").read_text())
    sum_claim=entries["MRC-SIMPLE-COVER-RANK-ADDITIVITY-20261117"]
    assert sum_claim["state"]=="proved" and sum_claim["id"] in objective["dependencies"]
    assert sum_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_cover_sums.py")
    assert not sum_claim["formal_verification"] and not sum_claim["external_review"] and not sum_claim["closes_problems"]
    assert sums["status"]=="PASS independent simple-cover sum inputs" and sums["control_count"]==8 and sums["rejected_corruptions"]==4
    assert [(r["degree"],r["sphere_graph_cycles"],r["conditional_cover_rank"],r["meridional_lower"]) for r in sums["controls"]]==[(2,0,2,3),(3,1,5,7),(4,0,4,7),(3,2,8,10),(5,0,6,10),(3,0,3,5),(4,0,5,8),(3,0,2,4)]
    assert sums["T47_boundary"]["double_cover_rank"]==2 and sums["T47_boundary"]["maximum_simple_cover_rank_bound"]==3 and sums["T47_boundary"]["meridional_rank"]==4
    assert sum_sources["status"]=="PASS scoped cover-sum source audit" and sum_sources["catalogue_distinct_ids"]==len(catalogue["sources"])
    assert [r["id"] for r in sum_sources["new_records"]]==["STALL71","SCHULT03"]
    for receipt in [sums,sum_sources,json.loads((HERE/"cover_sum_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"cover-sum input changed: {path}"
    permutation=json.loads((HERE/"double_transposition_verification.json").read_text())
    permutation_sources=json.loads((HERE/"double_transposition_source_audit.json").read_text())
    permutation_claim=entries["MRC-DOUBLE-TRANSPOSITION-INDEX-GATE-20261117"]
    assert permutation_claim["state"]=="proved" and permutation_claim["id"] in objective["dependencies"]
    assert permutation_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_double_transpositions.py")
    assert permutation["status"]=="PASS independent degree8/9 double-transposition gate"
    assert [(r["degree"],r["catalogue"]["transitive_orbits"],r["catalogue"]["raw_normalized_tuples"]) for r in permutation["degrees"]]==[(8,18450,9261000),(9,6561,54010152)]
    assert all(len(d["targets"])==16 and all(r["hits"]==0 for r in d["targets"]) for d in permutation["degrees"])
    assert [(c["pairs"],c["image_orders"]) for c in permutation["controls"]]==[(14,{"1344":14}),(15,{"20160":15})]
    assert len(permutation["rejected_corruptions"])==6
    assert permutation_sources["catalogue_distinct_ids"]==len(catalogue["sources"]) and permutation_sources["new_record"]["id"]=="ATLASA8"
    for receipt in [permutation,permutation_sources]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"permutation-index input changed: {path}"
    ternary=json.loads((HERE/"ternary_certificate_verification.json").read_text())
    ternary_bases=json.loads((HERE/"ternary_base_verification.json").read_text())
    ternary_covers=json.loads((HERE/"ternary_cover_module_verification.json").read_text())
    ternary_claim=entries["MRC-TERNARY-AFFINE-COVER-20261117"]
    assert ternary_claim["state"]=="proved" and ternary_claim["id"] in objective["dependencies"]
    assert ternary_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_ternary_certificates.py")
    assert ternary["status"]=="PASS scoped complete ternary affine and cover package"
    assert ternary["targets"]==16 and ternary["tuple_orbits"]==[21632,34084] and ternary["raw_tuples"]==[1124864,1601613]
    assert ternary["all_rank_one_admissible_affine_lifts_principal"] and ternary["rejected_corruptions"]==10
    assert ternary["full_GL3_F3_representation_orbits"]==16 and ternary["knots_with_full_GL3_F3_image"]==8
    assert ternary["ordinary_cover_module_tests"]==32 and ternary["ordinary_Z1_dimension"]==4 and ternary["required_Z1_dimension"]==7
    assert ternary["projective_index13_control"]["homology_gain"]==15 and ternary["projective_index13_control"]["sharp_gain_threshold"]==13
    assert ternary["new_target_closures"]==0 and not ternary["full_MRC_proved"]
    for receipt in [ternary,ternary_bases,ternary_covers]+[json.loads((HERE/p).read_text()) for p in [
        "ternary_transvection_control.json","ternary_reflection_control.json","ternary_transvection_results.json",
        "ternary_reflection_results.json","ternary_cover_module_results.json",
        "ternary_transvection_control_verification.json","ternary_reflection_control_verification.json"]]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"ternary package input changed: {path}"
    projective=json.loads((HERE/"projective13_cover_verification.json").read_text())
    projective_sources=json.loads((HERE/"projective13_source_audit.json").read_text())
    degree_control=json.loads((HERE/"characteristic_degree_control_verification.json").read_text())
    projective_claim=entries["MRC-SL3-PROJECTIVE13-HOMOLOGY-20261117"]
    assert projective_claim["state"]=="proved" and projective_claim["id"] in objective["dependencies"]
    assert projective_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_projective13_covers.py")
    assert projective["status"]=="PASS independent projective13 all-field cover test"
    assert (projective["source_actions"],projective["distinct_source_covers"],projective["knots"])==(32,16,8)
    assert projective["maximum_Betti_gain_except13"]==4 and projective["maximum_augmentation_Z1_all_fields"]==16 and projective["required_augmentation_Z1"]==25
    assert all(r["H1_free_rank"]==2 and r["H1_torsion_factors"] in [[2,882878],[2,3128602]] and r["augmentation_Z1_F13"]==14 for r in projective["results"])
    assert projective["retained_quotient_kernels_per_source_among_these_maps"]==1 and projective["new_target_closures"]==0
    assert len(projective["rejected_corruptions"])==4
    assert degree_control["status"]=="PASS independent positive degree-characteristic control"
    assert (degree_control["degree"],degree_control["characteristic"],degree_control["gain"],degree_control["augmentation_Z1"],degree_control["split_prediction"])==(13,13,23,36,35)
    assert len(degree_control["rejected_corruptions"])==2
    assert projective_sources["status"]=="PASS pinned projective13 structural source cues" and not projective_sources["common_cover_or_mutation_proved"]
    assert [r["double_cover_H1_order"] for r in projective_sources["groups"]]==[391,401]
    for receipt in [projective,projective_sources,degree_control,json.loads((HERE/"characteristic_degree_control.json").read_text()),json.loads((HERE/"projective13_cover_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"projective13 input changed: {path}"
    mutation=json.loads((HERE/"mutation_cover_verification.json").read_text())
    mutation_sources=json.loads((HERE/"mutation_source_audit.json").read_text())
    mutation_claim=entries["MRC-MUTATION-COMMON-DOUBLE-COVERS-20261117"]
    assert mutation_claim["state"]=="proved" and mutation_claim["id"] in objective["dependencies"]
    assert mutation_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_mutation_covers.py")
    assert mutation["status"]=="PASS independent mutation/common-cover inputs"
    assert (len(mutation["diagrams"]),len(mutation["pairs"]),len(mutation["controls"]),len(mutation["corruptions_rejected"]))==(8,6,2,4)
    assert [r["switches"] for r in mutation["pairs"]]==[1,1,0,1,0,1]
    assert all(r["cycle_vectors"]==128 and r["edge_subsets_checked"]==8192 for r in mutation["pairs"])
    assert not mutation["actual_rank3_proved"]
    assert mutation_sources["status"]=="PASS scoped mutation source audit" and mutation_sources["catalogue_distinct_ids"]==len(catalogue["sources"])
    assert [r["id"] for r in mutation_sources["new_records"]]==["GRE13","MOR14"]
    for receipt in [mutation,mutation_sources,json.loads((HERE/"mutation_cover_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"mutation package input changed: {path}"
    compact=json.loads((HERE/"compact_cover_verification.json").read_text())
    ports=json.loads((HERE/"port_assembly_verification.json").read_text())
    compact_claim=entries["MRC-COMPACT-FOUR-PORT-COVER-GROUPS-20261117"]
    assert compact_claim["state"]=="proved" and compact_claim["id"] in objective["dependencies"]
    assert compact_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_port_assemblies.py")
    assert compact["status"]=="PASS independent compact common-cover presentations"
    assert (len(compact["targets"]),len(compact["controls"]),len(compact["corruptions_rejected"]))==(2,3,5)
    assert [[len(w) for w in r["presentation"]["relators"]] for r in compact["targets"]]==[[9,13,15],[12,15,11]]
    assert [r["abelianization_smith"] for r in compact["targets"]]==[[1,1,391],[1,1,401]]
    assert all(r["nonabelian_quotient"]["order"]==5616 and r["actual_rank_interval"]==[2,3] for r in compact["targets"])
    assert [r["graph_network"]["variable_branch_T_F"] for r in compact["targets"]]==[[19,14],[14,19]]
    assert ports["status"]=="PASS independent four-port core assemblies" and len(ports["corruptions_rejected"])==4
    assert ports["common_formal_inner_and_outer_groups"] and ports["outer_ribbon_map_orientation"]==-1
    assert all(r["inner_relator_length"]==11 and r["outer_relator_length"]==7 and r["unresolved_arcs_after_outer_generation"]==[0,1,2,3,4] for r in ports["targets"])
    assert not ports["rank2_proved"] and not ports["rank3_proved"] and ports["new_target_closures"]==0
    geometry=json.loads((HERE/"compact_cover_geometry_cues.json").read_text())
    recognition=json.loads((HERE/"compact_cover_recognition_cues.json").read_text())
    assert geometry["status"]==recognition["status"]=="COMPLETE_CUES"
    assert recognition["results"][0]["manifold"]=="L(5,2)"
    assert all(r["recognition"] is None and r["compact_algebra_simplification"]["generators"]==3 for r in recognition["results"][1:])
    for receipt in [compact,ports,geometry,recognition,json.loads((HERE/"compact_cover_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"compact/port input changed: {path}"
    four_piece=json.loads((HERE/"four_piece_cover_verification.json").read_text())
    four_sources=json.loads((HERE/"four_piece_source_audit.json").read_text())
    four_claim=entries["MRC-FOUR-PIECE-DOUBLE-COVER-SHARP-FOUR-20261117"]
    assert four_claim["state"]=="proved" and four_claim["id"] in objective["dependencies"]
    assert four_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_four_piece_covers.py")
    assert four_piece["status"]=="PASS independent four-piece peripheral and upper4 inputs"
    assert (len(four_piece["roots"]),len(four_piece["upper4"]),len(four_piece["corruptions_rejected"]))==(2,8,5)
    assert all(r["fiber_intersections"]==[1,1,1] and r["proposed_JSJ_pieces"]==4 for r in four_piece["roots"])
    assert [r["source_H1_order"] for r in four_piece["roots"]]==[391,401]
    new_equalities=[r["name"] for r in four_piece["upper4"]]
    assert sorted(new_equalities)==['13a_129','13a_135','13a_136','13a_137','13a_194','13a_196','13a_197','13a_198']
    assert four_sources["status"]=="PASS scoped four-piece source audit" and four_sources["catalogue_distinct_ids"]==46
    assert [r["id"] for r in four_sources["new_records"]]==["BW05","HAT00","MW16"]
    assert four_sources["classical_rank2_maximum_JSJ_pieces"]==3 and max(four_sources["rank2_necessary_type_piece_counts"])==3
    for receipt in [four_piece,four_sources,json.loads((HERE/"four_piece_cover_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():
            assert sha(HERE/path)==digest, f"four-piece input changed: {path}"
    three_piece=json.loads((HERE/"three_piece_cover_verification.json").read_text())
    three_sources=json.loads((HERE/"three_piece_source_audit.json").read_text())
    three_claim=entries["MRC-THREE-PIECE-DOUBLE-COVER-SHARP-FOUR-20261117"]
    assert three_claim["state"]=="proved" and three_claim["id"] in objective["dependencies"]
    assert three_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_three_piece_covers.py")
    assert three_piece["status"]=="PASS independent three-root peripheral and upper4 inputs; five transfers unproved"
    assert [r["name"] for r in three_piece["roots"]]==['13a_233','13a_239','13a_253']
    assert [r["H1_order"] for r in three_piece["roots"]]==[403,361,337]
    assert [r["central_filling_slope_a_f"] for r in three_piece["roots"]]==[[3,-2],[-2,-1],[-3,2]]
    assert (len(three_piece["mutation_pairs"]),len(three_piece["unproved_transfers"]),len(three_piece["upper4"]),len(three_piece["corruptions_rejected"]))==(0,5,8,5)
    assert sum(r["fragment_reductions"] for r in three_piece["small_source_audit"])==14
    assert all(r['source_determinant']!=r['target_determinant'] for r in three_piece['unproved_transfers'])
    assert three_sources["status"]=="PASS scoped three-piece case6 source audit" and three_sources["rank_source"]["only_three_piece_case"]==6
    for receipt in [three_piece,three_sources,json.loads((HERE/"three_piece_cover_results.json").read_text()),json.loads((HERE/"remaining_cover_structure_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():assert sha(HERE/path)==digest,f"three-piece input changed: {path}"
    final_five=json.loads((HERE/"final_five_cover_verification.json").read_text())
    final_sources=json.loads((HERE/"final_five_source_audit.json").read_text())
    final_claim=entries["MRC-ORIGINAL-SIXTEEN-DOUBLE-COVER-SHARP-FOUR-20261117"]
    assert final_claim["state"]=="proved" and final_claim["id"] in objective["dependencies"]
    assert final_claim["artifact_hash"]=="sha256:"+sha(HERE/"verify_final_five_covers.py")
    assert final_five["status"]=="PASS independent final-five source, peripheral, filling and upper4 inputs"
    assert [r["name"] for r in final_five["targets"]]==['13a_234','13a_242','13a_243','13a_244','13a_265']
    assert [r["source_H1_order"] for r in final_five["targets"]]==[389,359,337,311,311]
    assert [r["central_filling_slope_c_f"] for r in final_five["targets"]]==[[3,-2],[-2,1],[-3,2],[-3,2],[-3,2]]
    assert [r["inner_fiber_sign"] for r in final_five["targets"]]==[1,1,1,1,-1]
    assert all(r['nonabelian_filling_certificate']['image_order']==6 for r in final_five['targets'])
    assert (len(final_five["source_audit"]),sum(r['fragment_reductions'] for r in final_five['source_audit']),len(final_five["positive_controls"]),len(final_five["upper4"]),len(final_five["corruptions_rejected"]))==(8,32,3,5,6)
    assert final_five["actual_rank_transfers"]==0
    assert final_sources["status"]=="PASS scoped final-five classical source reuse" and final_sources["catalogue_distinct_ids"]==46
    for receipt in [final_five,final_sources,json.loads((HERE/"final_five_structure_results.json").read_text()),json.loads((HERE/"final_five_peripheral_results.json").read_text())]:
        for path,digest in receipt["input_sha256"].items():assert sha(HERE/path)==digest,f"final-five input changed: {path}"
    original_equalities=sorted(new_equalities+[r['name'] for r in three_piece['roots']]+[r['name'] for r in final_five['targets']])
    assert original_equalities==sorted(r['name'] for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases'])
    assert len(set(original_equalities))==16
    family=json.loads((HERE/"graph_cover_family_scope_verification.json").read_text())
    family_sources=json.loads((HERE/"graph_cover_family_source_audit.json").read_text())
    family_claim=entries["MRC-ORIGINAL-SIXTEEN-KNOWN-GRAPH-COVER-FAMILY-20261117"]
    assert family_claim['state']=='proved' and family_claim['id'] in objective['dependencies']
    assert family_claim['artifact_hash']=='sha256:'+sha(HERE/'audit_graph_cover_family_scope.py')
    assert family['status']=='PASS sixteen exact case/source bindings for known graph-cover family'
    assert (family['existing_family_count'],family['new_MRC_family_cases'],len(family['cases']))==(16,0,16)
    assert [r['name'] for r in family['cases']]==original_equalities
    assert all(r['group_rank']==3 and r['upper4_steps']==9 and r['existing_MRC_family']=='BJW18 published Corollary1.3' for r in family['cases'])
    assert family['topology_proved_by_this_script'] is False and family['full_MRC_proved'] is False
    assert family_sources['catalogue_distinct_ids']==46 and family_sources['new_ids']==['AFW15']
    assert {'AFW15','BJW18'}=={r['id'] for r in family_sources['records']}
    assert all(record['catalogue_at_original_checkpoint']==45 and record['catalogue_distinct_ids']==46 for record in [four_sources,three_sources,final_sources])
    for receipt in [family,family_sources]:
        for path,digest in receipt['input_sha256'].items():assert sha(HERE/path)==digest,f'known-family input changed: {path}'
    integration = json.loads((HERE / "integration.json").read_text())
    report = {
        "date": "2026-11-17",
        "goal_status": "active",
        "goal_complete": False,
        "previous_goal_turn_classification": "progress",
        "previous_progress_evidence": [
            "FINAL_FIVE_DOUBLE_COVERS_SHARP_FOUR.tex", "final_five_cover_verification.json", "final_five_source_audit.json"
        ],
        "this_turn_classification": "progress",
        "completion_audit": {
            "full_MRC_proof_or_counterexample": "NOT ACHIEVED; OP-MRC is open",
            "exhaustive_literature_and_source_audit": "INCOMPLETE; broad35-source map and18 primary snapshots, now including primary ATLAS A5 and L32 representation lists. Classification theorems and standard cross-characteristic reduction are explicit imports; supplied modular models are checked by full matrix-algebra spans. Original Dutra gap remains checked with no affected known-family dependency. Pinned external GL4 independent-direction/det guard and bounded random three-seed search inspected; this concerns coverage, not validity of the243 certificates. Cornwell torus/maximal-braid-rank proofs and needed GGLW one-step/peripheral/iterable arguments are read; full recent papers, original Hall/LM/orbifold and broad priority audits remain open",
            "construct_representations_and_matching_upper_bounds": "ACHIEVED for known-family controls,4earlier named simple3cover cases and all16original-target double-cover rank3 cases with individual upper4 traces. Total20named panel equalities proved here. These use explicit classical topology/rank imports and are unformalized/externally unreviewed. FullMRC remains open",
            "new_previously_unresolved_knot_case": "NOT ESTABLISHED. All16original targets are now proved members of the existing BJW18 graph-cover MRC family using AFW15 closed rigidity; they are not cases outside known MRC families. Twenty explicit named equalities remain certified here. Exact named bridge/rank computational priority and the4earlier sharp seeds known-family scope remain UNESTABLISHED",
            "bridge_unpublished_quadratic_argument": "ACHIEVED through signed lifts, all-field nullity, radical towers and cover homology/rank. The rank2-form part is now characterized exactly: for a fixed S3 colouring and simple3cover Y, dimension k exists iff b1(Y;F2)>=k-1. Two S4 kernels over the same S3 kernel force lower4. Every transitive rank2 image has order24*8^(k-2). This gives a construction and exact stopping condition, but no new open knot or universal MRC proof. Earlier formal icosahedral deformation boundaries retained",
            "test_if_affine_criterion_can_solve_full_goal": "UNIVERSALITY EXCLUDED on the stated infinite torus-knot family; other methods remain necessary"
        },
        "verified_controls": {k: verification[k] for k in [
            "exhaustive_three_by_three_binary_matrices", "relator_word_sign_cases_through_length_six",
            "exhaustive_connected_simple_signed_graphs_n3_n4", "parallel_edge_signings",
            "independent_group_enumerations", "invalid_symmetric_relator_rejected"]},
        "primary_snapshots_hash_checked": len(source_manifest),
        "official_repository_inputs_hash_checked": 12,
        "lift_dimension_progress": {
            "uniform_counting_lower_bound": "meridional rank at least full signed-lift dimension ell",
            "S4_affine_realization": "2^(ell-1) points, form rank2, image order24*2^(3*(ell-3))",
            "independent_small_group_orders": dimension["independent_group_orders"],
            "frozen_panel_cases_independently_replayed": 16,
            "panel_bounds": "Binary coefficients only: lower4, upper5 under every S3/S4 colouring. The all-characteristic extension now reaches sharp5 over F5; see all_characteristic_progress",
            "new_previously_unresolved_knot_case": False,
            "scope": "The panel's lack of sharp finite Coxeter hits does not establish otherwise-unresolved MRC status. Uniform theorems are written proofs, not formally verified."
        },
        "PDF_sha256": outputs,
        "cross_programme_handoff": {
            "review": "CROSS_PROGRAMME_SCOPE.json",
            "centralizer_lemma": "Written group argument audited: natural H_q centralizer is exactly constant radical translations (C2)^q",
            "surface_exclusion": "Retained conditionally with extra Euler, exact-volume and whole-plane hypotheses; full geometric proof not independently adopted here",
            "knot_targets_removed": False
        },
        "all_characteristic_progress": {
            "uniform_bound": "mr at least dim Z over every field; exact all-characteristic maximum is number of columns minus unit Smith factors",
            "sharp_source_table_diagrams": 16,
            "sharp_bounds": "mr=bridge=5 via F5 lifts and five-seed upper certificates",
            "sharp_case_novelty": "EXCLUDED for this16-case panel: all independently identified five-noninteger-tangle Montesinos, covered by Boileau-Zieschang1985",
            "independent_lift_matrices": integer["independently_replayed_lift_matrices"],
            "smith_gcd_of_minors_controls": integer["smith_two_by_two_gcd_of_minors_controls"],
            "current13_panel": "Same first16 now have all symmetric degrees/all fields excluded as a sharp root-lift route: maximum bound2. The original binary common-form affine criterion is impossible in every dimension; see branched_cover_progress",
            "October_release_algebra_and_upper_certificates_independently_replayed": 243,
            "October_release_named_census_matches": 243,
            "KnotInfo_snapshots_checked": ["2026.10.5", "2026.11.2"],
            "Montesinos_named_certificate_transfers": 16,
            "independent_Montesinos_tangle_constructions": 16,
            "parabolic_panel": "Same16 unresolved13crossing targets, complete normalized noncommuting SL2 transvection classes in fields3,5,7: only root-lift dimensions2 or3, no sharp4. Commuting families have dimension1. Separate matrix/dense-affine replay passes; PD parser shared. Other fields untested",
            "binary_realization_of_F5_lifts": "NOT ASSERTED"
        },
        "branched_cover_progress": {
            "uniform_rank_bound": "mr(K)>=r-1+rank pi1(Y_phi) for the simple r-sheet transposition branched cover",
            "uniform_homology_identity": "dim Z_phi(F)=r-1+dim H1(Y_phi;F); nonunit integer Smith factors equal H1 torsion factors, including characteristics dividing r",
            "nonabelian_count": "admissible A wreath lifts=|A|^(r-1)*|Hom(pi1 Y,A)|; r2 Kelly-Wada core-group precedent credited",
            "sharpness_gate": "rank pi1(Y)>=b-r+1 suffices with b-bridge upper certificate; rank=genus and cover maximal genus not assumed",
            "full_directed_edge_homology_controls": 9,
            "known_sharp_control": "15a1381 S4 cover H1=Z+Z/5 explains characteristic5 improvement",
            "current13_all_degrees_fields": "Bridge<=4 forces degree<=5. S3/S4/S5 have no base; S2 Smith gives maximum nullity2 over every field on all16 cases",
            "current13_original_binary_all_dimensions": "EXCLUDED: k+1<=mr<=b<=4 implies k<=3, nonzero alternating form rank2, forcing absent S3 quotient",
            "correction": "Earlier suggestion of higher-rank binary forms on these targets withdrawn",
            "A5_core_counts": {"1":4,"121":8,"241":4},
            "A5_sharp_four_threshold": "Raw count: >3600 normalized colourings; no hit. Stronger same-cover kernel criterion: >=20 distinct A5 quotient kernels; see joint_simple_kernel_progress",
            "finite_index_extension": "For an index-d subgroup J of pi1(Y), mr>=r+ceil((dim H1(J;F)-1)/d), by Schreier; index5 double-cover threshold for sharp4 is7",
            "first_A5_witness_index5_covers": "Twelve fixed first nonconstant witnesses, all image60. All covers have H1=Z/(3D), maximum field dimension1. Independent Fox-prefix/Smith replay and free/Kurosh controls pass. Other witnesses and covers untested",
            "verification_boundary": "Separate full-chain construction with elementary Smith replay and separate PSL2(F5) coefficient-group replay. PD parser and Smith producer/replay are shared. Uniform working proofs are unformalized and externally unreviewed",
            "new_previously_unresolved_knot_cases": 0
        },
        "quadratic_coefficient_progress": {
            "uniform_field_reduction": "O=Z[tau], tau^2=tau+1. Regular integer expansion doubles generic rank. Nonzero Smith factors isolate rank-drop primes away from5; test split/inert residues and5 separately. Every coefficient field covered after extension",
            "uniform_cocycle_comparison": "For meridian defect rank d>0, mr>=ceil(ell_root/d)>=1+ceil((dim Z1-dim V)/d). Full cocycles cannot improve the restricted-root estimate on the same representation",
            "original_nondegenerate_bound_recovered": "Transitivity forces direction span V; nondegenerate B gives V^G=0 and k principal cocycles. The transitive affine cocycle is nonprincipal, so ell>=k+1. Degenerate case still needs the retained stronger quadratic proof",
            "structural_base_gate": "All12 existing first nonconstant A5 core maps have unique double-transposition deck implementers, giving onto H3 reflection bases. No new map or knot substituted",
            "fixed_bases": 12,
            "integer_expanded_Smith": "twenty unit factors per matrix; generic nullity3, separate characteristic5 nullity3, hence exactly3 over every field",
            "finite_image_factorization": "Universal H3 root lifts have3 parameters in every field. Integral split maps LP=I and MP=0 plus meridian words for simple generators show every lift of these fixed bases factors through H3, image<=120",
            "characteristic2_boundary": "F4 action on64 points has no common fixed point and image120, but four orbits12,12,20,20. Underlying binary meridian rank2, so original rank-one binary criterion not satisfied",
            "replay": "Independent ordinary integer6x6 group and reverse affine-point evaluation, elementary Smith replay, residue bases, split factorization and bit-vector orbit replay; shared PD parser and Smith checker disclosed",
            "source_audit": "GGLW2026 Lemma6.2 full proof/Remark6.3 checked; entire paper and general priority remain unaudited",
            "new_previously_unresolved_knot_cases": 0,
            "boundary": "Only fixed first-map H3 bases excluded. Other maps, representations, nonlinear deformations and geometric upper improvements remain open"
        },
        "joint_simple_kernel_progress": {
            "user_proposed_application": "Twenty distinct A5 quotient kernels of one double-cover group force rank>=3 and mr>=4; with bridge<=4 this is sharp",
            "classical_source": "Hall1936 via Burness-Crestani2013 p2 explicit delta(A5)=19;2280 generating pairs divided by120 automorphisms",
            "uniform_bound": "kappa_S(H)>|Epi(F_d,S)|/|Aut S| forces rank H>=d+1 and mr>=r+d for a simple degree-r branched cover",
            "joint_surjection": "For nonabelian simple S, distinct kernels on one H jointly surject H onto S^kappa",
            "single_domain_guard": "Retained source identifiers, Aut-orbit canonicalization and rejection controls; maps on different knots are never pooled",
            "complete_frozen_uncertainty_panel": 293,
            "same_knot_kernel_histogram": joint["same_knot_kernel_histogram"],
            "maximum_one_knot": 5,
            "new_sharp_uncertainty_cases": 0,
            "known_strict_improvement_control": {"knot":"P(2,3,3,5)","kernels":23,"Hom":2901,"raw_lower":3,"joint_lower":4,"bridge_upper":4},
            "known_family": "Four-tangle Montesinos with denominators divisible by2,3,3,5, via one-cover pullback from Gamma(2,3,3,5); already Boileau-Zieschang1985",
            "finite_denominator_patterns":40,
            "patterns_reaching_twenty":39,
            "exception":"Gamma(2,3,3,3) has12; no universal A5 detection claim",
            "trefoil_sum_controls":{"two":1,"three":49},
            "replay":"Separate PSL2/PGL2(F5), exhaustive orbit marking and crossing closure reproduce all counts; direct permutation relations and two independent tangle constructions checked. PD parser shared",
            "proof_boundary":"Uniform written proofs unformalized and externally unreviewed; no novel knot or general MRC result"
        },
        "geometric_probe": {"candidates":30,"marked_diagrams":240,"matches":0,"scope":"Exploratory bounded PD probe only; no Montesinos nonmembership, flype or isotopy exclusion; no new mathematical claim registered"},
        "augmentation_cover_progress": {
            "user_proposed_bound":"For index-n H<=G and characteristic not dividing n, Delta=b1(H)-b1(G)=dim H1(G,V_(n-1)); rank G>=1+ceil(Delta/(n-1)); mr>=r+ceil(Delta/(n-1))",
            "sharp_four_gate":"In the double cover with index5, gain>=5 suffices and requires no knot-group root extension",
            "frozen_targets":16,"all_A5_kernel_classes_tested":16,"knots_with_maps":12,
            "all_connected_index5_covers_covered":True,
            "completeness_reason":"All base H1 are cyclic of order coprime to30. Any transitive S5 image has trivial sign, then every proper nontrivial A5 image would force an absent abelian quotient at2,3,5. Thus every index5 cover has A5 monodromy and is in the Aut-kernel panel",
            "integral_result":"Every base H1=Z/D and every tested subgroup H1=Z/(3D), with gcd(D,30)=1",
            "all_prime_gains":"One at3; zero at every other prime. No gain5 and no new sharp case",
            "old_exclusion_comparison":"Different group and module from fixed knot-root representations; the prior restricted/full comparison on one representation did not exclude this test",
            "controls":{"F2_gain":4,"two_trefoil_gain":2,"three_trefoil_gain":5},
            "characteristic5_countercontrol":augmentation["characteristic5_countercontrol"],
            "replay":"Independent Fox-prefix integer construction and third direct augmentation-cocycle matrices; Smith-step replay and all-prime classification. PD parser and Smith replay shared",
            "peripheral_gate":"For a diagonalizable rank-one meridian preserving a nondegenerate symmetric form, the commuting longitude has lambda^2=1, forbidden by the GGLW step. A useful deformation must escape orthogonality and come from actual augmentation points",
            "remaining":"Other cover degrees, modules and actual peripheral curves remain open; uniform proofs unformalized and externally unreviewed"
        },
        "whitehead_dependency_progress": {
            "all16_in_current_external_open_list":True,"current_blob":whitehead["current_open_list_blob"],
            "independently_checked_Dutra_gap":"Original p7 nontrivial free word [2,-1,-2,1] lies in both U_alpha,U_omega, contradicting Lemma3.1 at g1; original p12 invokes it in Lemma3.8(2)",
            "own_established_family_dependencies_on_affected_2022_theorem":0,
            "disposition":"Keep broader claim flagged, main conclusion not disproved. Preserve distinct2017 graph-manifold/braid-pattern result. GGLW full proof and a repair of the broader theorem remain unaudited"
        },
        "six_point_and_formal_progress": {
            "mathematical_gate":"Six-point A5 augmentation module has irreducible dimension5 character, orthogonal to V4; index6 permits characteristic5",
            "frozen_A5_kernel_classes":16,"cover_matrix_shape":[84,73],
            "integral_homology":"Eight source classes have H1=Z/(4D), eight have Z/(28D); all D coprime to210",
            "sharp_four_threshold":6,"maximum_gain":1,"new_sharp_cases":0,
            "all_degree6_covers_excluded":False,
            "scope":"All A5-induced six-point covers on the retained classes only; A6/other monodromy untested",
            "coefficient_group_replay":"Independent PSL2(F5) action on P1(F5), explicit group isomorphism and sheet bijection, transported Fox-prefix/Smith replay and direct augmentation cocycles",
            "all_field_module_boundary":"V4 at5 has (invariants,H1,Z1)=(1,1,4); V5 at2 has(1,2,6), at3(1,1,5). With split cases, neither augmentation cocycle estimate forces ordinary cover rank3 in any field",
            "orbifold_descent":"All V5 Q cohomology is zero at p>=5; the eight nonzero mod7 cover classes do not descend",
            "uniform_form_lifting":"In characteristic different from2, if rho0 preserves B0 and H1(Q,Sym2 W*)=0, every formal representation rho_t reducing to rho0 preserves a varying nondegenerate symmetric form. At each order the invariance defect is a cocycle killed by a form correction",
            "complex_application":"Sym2 of either standard icosahedral character equals the six-point permutation character. Shapiro uses finite H1(H6) and index2 restriction to give the form-lifting vanishing at every tested H3 point",
            "formal_peripheral_boundary":"Keeping mu^2=1 keeps a rank-one meridian eigensummand; its commuting longitude preserves the lifted form, so lambda(t)^2=1. Varying meridian eigenvalues, other components/dimensions remain open",
            "primitive_word_probe":"All16 old three-generator core presentations,444 candidate primitive consequences per knot, no certificate. Substitutions/inverse automorphisms and controls checked; failed heuristic not independently rerun or promoted to any actual group-rank conclusion",
            "uniform_proof_assurance":"Written working proofs, unformalized and externally unreviewed; finite checks establish specific input hypotheses, not full MRC"
        },
        "torus_whitehead_progress": {
            "result":"For coprime2<=p<q, ar=mr=b=p*2^N for every finite untwisted Whitehead iterate of T(p,q), for the GGLW pattern/conventions",
            "status":"Source-derived corollary with explicit November2026 GGLW preprint dependency; no novelty claim or independent full-source certification",
            "actual_curve":"Companion-matrix rank-p irreducible KCH family realizes Cornwell's lambda=(-1)^(p+1)*mu^(-(p-1)q); eta=mu^(-2(p-1)q) has iterable orders(-1,2(p-1)q)",
            "broader_gate":"An n-strand closure with rank-n augmentation and nonzero writhe w has an actual curve by Cornwell2017 Theorem1.10 proof; orders(-1,2w) give equality n*2^N for all iterates",
            "T35_control":"mu=s^3,lambda=s^-30; s=-1 gives H3 order120; actual variation escapes orthogonality",
            "finite_torus_controls":39,"whitehead_triples":2976,"valid_whitehead_blocks":156,
            "panel_closures":0,"panel_varying_meridian_tangents_executed":True,
            "not_certified":"Full MRC, broad priority, specialist/formal verification, entire GGLW microlocal/two-bridge theory and actual peripheral curves for the retained16 targets"
        },
        "meridian_deformation_progress": {
            "frozen_bases":16,"knots_with_points":12,"four_no_A5_targets_not_covered":True,
            "normalized_tangent_matrix_shape":[144,28],"exact_rank":23,"nullity":5,
            "residual_conjugation_dimension":4,"modulo_conjugation_dimension":1,
            "eigenvalue_differential":"Nonzero in every case, over Q(tau) and both complex embeddings",
            "fixed_eigenvalue_rigidity":"An order-by-order conjugation proof makes every formal rank3 fixed-involution family conjugate to the original point, strengthening orthogonality",
            "integrability_gate":"After transverse gauge, the local germ is either an actual curve with parameter s=(mu+1)/2, or an isolated thick point cut out by a finite-order univariate residual. The alternative is not decided",
            "stabilization_theorem":"If both restricted W/Wdual H1 vanish, every formal rank-one-meridian family through rho0 plus a finite trivial summand splits, even allowing varying eigenvalue",
            "restricted_modules_checked":32,"restricted_cocycle_dimension":3,"principal_dimension":3,
            "formal_higher_rank_boundary":"No irreducible higher-rank family through these stabilized special fibres. Other special fibres, representations and components remain open",
            "finite_jets":{"embeddings":32,"prime":11,"degree":3,"all_pass":True,"scope":"Residue-field finite jets only; no characteristic-zero curve or all-order lifting claim"},
            "replay":"Separate Fraction arithmetic, oriented dual-number relators, all minors, gauge and affine-cocycle columns. Jets use separate entrywise polynomial adjugates/scalar inverses. PD parser and retained H3 coefficient inputs shared",
            "new_sharp_cases":0,"formal_or_specialist_review":False
        },
        "complete_A5_module_progress": {
            "uniform_reduction":"The maximum dim Z1/dim V through a fixed finite quotient occurs on a simple module; composition-series inequality includes nonsplit extensions and commutes with scalar extension",
            "classification":"Absolute simple lists imported from primary ATLAS; modular models independently identified with full matrix-algebra spans and conjugate traces",
            "new_ordinary_W3_matrices":16,"new_characteristic2_simple_actions":32,
            "integer_shape":[84,78],"nonunit_Smith_patterns":"Eight have[5], eight[7,35]; invariant factors1^4,2,2. Separate ramified5 and all etale embeddings checked",
            "all_module_bound":"dim Z1(G,V)<=3*dim V/2 over every field and every finite-dimensional module through any single A5 homomorphism on the first16 targets, including indecomposable modules",
            "ordinary_rank_estimate_maximum":2,"actual_rank2_proved":False,
            "all_A5_induced_covers":"For H=phi^-1(L), n=[G:H], b1(H;k)<=floor(n/2)+1; no Schreier or split augmentation-gain threshold can force rank3. Other monodromy not covered",
            "rational_homology":"All A5-induced covers have b1=0 over Q, including all60-sheet kernel covers",
            "joint_product_test":{"knots":4,"joint_image_order":3600,"natural_tensor_actions":16,"field":"F4","module_dimension":4,"Z1_dimension":6,"H1_dimension":2,"threshold_H1":5,"scope":"Natural2x2 tensors only, other product modules and fields untested"},
            "independent_replay":"Integer6x6 coefficient model and reverse affine-point matrices, Smith replay, residue ranks, quotient identifications, conjugate traces, matrix-algebra spans and nonsplit/freeF2/F3 controls; source PD/maps/Smith replay shared",
            "new_sharp_cases":0,"formal_or_external_review":False
        },
        "GL3_and_singular_affine_progress": {
            "uniform_affine_lemma":"An invariant-free d-dimensional linear action and nonprincipal admissible cocycle homogenize to a rank-one meridian certificate with zero fixed space in dimension d+1",
            "known_sharp_family":"T(4,7*(2r+1)), r>=0, mr=bridge=4 by explicit binary matrices and torus meridian congruences; already known torus MRC, no new knot or priority claim",
            "positive_control":"T47 split AGL3(F2), order1344, transitive8; no invariant alternating form on natural3; four conjugate meridians have independent functionals and pairing rank3",
            "external_search_boundary":"Pinned glsearch requires independent4seed directions and det(pairing)!=0, excluding all homogeneous affine3 extensions. Valid criterion and243 certificates unchanged; w3 random moves are not an exhaustive3bridge exclusion",
            "frozen_targets":293,"quotient_maps":191,"knots_with_maps":154,
            "kernel_histogram":gl3["kernel_histogram"],"same_source_sharp_threshold":58,
            "deck_types":gl3["deck_type_histogram"],"ordinary_module_Z1":gl3["module_Z1_histograms"],
            "characteristic2_all_module_bound":"dim Z1<=5*dim V/3 through every retained epimorphism, all coefficient dimensions/fields of characteristic2 including nonsplit extensions; not all proper images or other characteristics",
            "binary_inner_deck_root_actions":58,"F7_both_deck_root_actions":191,"all_root_dimensions":3,
            "independent_replay":"Mobius model and full product transport, tuple orbit marking, iterative diagram closure, reverse affine matrices, Steinberg/adjoint8 intertwiner, dense adjoint7 and positive affine control; PD parser/source conventions shared",
            "new_open_case_closures":0,"actual_rank2_proved":False,"formal_or_external_review":False
        },
        "rank_two_equivalence_progress": {
            "uniform_radical_lemma":"Trivial-radical extensions with surviving meridian rank-one image lines preserve admissible cocycle spaces over every field, using knot H2=0 and meridian H1 evaluation",
            "exact_original_criterion":"Over a fixed S3 transposition colouring, h=b1(simple3branched cover;F2) gives dimension-k transitive rank2-form certificate iff h>=k-1; maximal k=h+1 and bound h+2 when h>0, no action when h=0",
            "relative_kernel_gate":"Exactly2^h-1 S4 transposition quotient kernels lie above the fixed S3 kernel; two distinct kernels on the same knot group force lower4 and construct a transitive eight-point original certificate",
            "proof_mechanism":"Dependence gives an explicit invariant quadratic function. Independence restricts to full W^h on the S3 kernel, and commutators supply every radical translation. Image order24*8^(k-2)",
            "mod2_cover_identity":"b1(Y3;F2)=1+b1(Y4;F2) for an S4 transposition colouring and its S3 quotient; no integral/all-field identity",
            "frozen_targets":293,"symmetric_bases":rank_two["symmetric_base_counts"],
            "S3_binary_dimensions":rank_two["S3_binary_lift_dimensions"],"A5_decks":rank_two["deck_histogram"],
            "all_symmetric_degree_field_maximum_bounds":{"1":6,"2":251,"3":36},
            "original_criterion_targets_with_no_representation_any_dimension":272,"targets_with_only_dimension2_bound3":21,
            "all_binary3_linear_bases":"Complete subgroup/normal-generation and module reduction closes proper-image gap; prior58 full GL3 bases reused, every root dimension<=3",
            "known_positive_controls":"Actual P233/P2333/P23333 diagrams give binary dimensions3/4/5 and transitive image orders24/192/1536, all basis-pair relators and dependent invariant checked",
            "nonlinear_three_sheet_checkpoint":{"covers":26,"initial_proved_cyclic_orders":nonlinear["cyclic_orders"],"initial_unresolved_two_generator_groups":7,"finite_targets":["S3","A4","D10","D14","A5"],"nonabelian_images_found":0,"coset_outcomes":nonlinear["coset_followup"],"superseded_by":"The trefoil-amalgam continuation below now determines all26 actual ranks; the old finite/coset records are retained checkpoints, and quotient images alone still do not prove cyclicity"},
            "replay_assurance":"Independent tuple orbit/trace closure, affine-point matrices and Smith replay; separate PSL2F5 deck model, full subgroup lattice and module normal forms. Independent signed edge words/Tietze steps and generated finite targets; shared PD parser and older routines disclosed",
            "new_open_cases":0,"uniform_proof_formal_or_specialist_review":False
        },
        "trefoil_amalgam_progress": {
            "named_equalities":["13a_3145:mr=bridge=4","13a_3147:mr=bridge=4","13a_3148:mr=bridge=4","13a_3149:mr=bridge=4"],
            "cover_models":["A+","A-","B+","B-"],"cover_group_rank":2,"groups_infinite_and_nonabelian":True,
            "uniform_proof":"The trefoil meridian/fibre subgroup injects as Z2, by infinite fibre and infinite-order meridian in C2*C3. Determinant+-1 edge maps form actual amalgams; classical normal form embeds the nonabelian vertex group. Exact two-generator presentations give rank2",
            "lower_bound":"mr>=r-1+rank pi1Y, with r3 gives4 for every meridional generating set",
            "upper_bound":"Four-seed diagram propagation plus published Wirtinger-number=bridge theorem gives b<=4",
            "independent_provenance":"Four PDs literally match both archived KnotInfo releases and independent HT census(+1 labels). Planar under-to-under arcs independently regenerate branch words/upper traces; old parser is used only for index transport",
            "all_three_sheet_ranks":{"C2":19,"C4":2,"trivial":1,"rank2_trefoil_amalgams":4},
            "all_four_sheet_ranks":{"trivial":21},"degree5_bases":0,"higher_simple_degrees_excluded_by_upper4":True,
            "simple_cover_rank_route_degrees_at_least3":"Frozen degree>=3 route gives exactly4sharp named cases; the new double-cover theorem adds8other names below",
            "abelianization_orders":[37,35,37,35],"universal_abelian_cover_H1":"0 over Z, exact Schreier/Smith replay",
            "solvable_quotient_boundary":"Each commutator subgroup is nontrivial/perfect, hence every solvable quotient is cyclic; enlarging solvable target searches cannot detect these nonabelian groups",
            "new_named_certificates":4,"worldwide_previously_unresolved_or_priority_claim":False,
            "source_status":"All4 in fetched external open list17Nov2026 and absent from its243 certificates; archived bridge intervals[3,4] checked. Prior-family and global novelty not established",
            "simple3cover_original16_targets_closed":0,"checkpoint_other289_panel_names_unsettled_here":True,
            "uniform_proof_formal_or_specialist_review":False,"approximate_manifold_recognition_used":False
        },
        "finite_commutator_progress": {
            "targets":[r["name"] for r in finite["cases"]],"finite_certificates":4,"degree":25,
            "checkpoint_relators_per_case":6,"full_cover_relators_per_case":39,"full_cover_generators_per_case":11,
            "commutator_witnesses":[r["commutator_witness"] for r in finite["cases"]],
            "construction":"Trefoil A5xC5 action on F5^2 with independent peripheral translations, conjugated by linear gluing maps and transported to original generators",
            "vertex_image_order":300,"full_two_generator_image_order":"not computed; not needed",
            "algebra_replay":"Standalone stdlib checker plus independent reverse-letter pointwise evaluation on all25 points. Full branched words rebuilt from independent planar arc traces; parser used only for old indices",
            "rejected_corruptions":4,"exact_cyclic_controls":finite["exact_cyclic_controls"],
            "new_named_cases_beyond_retained_four":0,
            "proof_boundary":"Additional finite proofs of the same4 mr=bridge4 equalities. Finite nonabelianity needs no residual-finiteness or amalgam normal-form input. General cover-rank lemma, named PD and bridge upper certificates remain dependencies. Global priority and full MRC unclaimed; unformalized and externally unreviewed"
        },
        "periodic_twist_family_progress": {
            "families":periodic["families"],"uniform_result":"For all N>=0, replacing the first source crossing by1+6N alternating crossings gives mr=bridge4",
            "finite_action":{"degree":75,"blocks":3,"block_size":25,"full_image_order":"not needed or computed"},
            "finite_lower_lemma":"A transitive r-block group generated by d exchanges fixing the other blocks pointwise has block-stabilizer label rank<=d-r+1, by spanning-tree gauges; nonabelian stabilizer forces d>=r+1",
            "periodicity_proof":"T(x,y)=(y,yxy), T^2 conjugates by(xy)^-1 and(xy)^3=1, so T^6=1. All original generator images remain outside, preserving the full finite image",
            "uniform_upper":"A fixed4outside-seed trace uses a two-input integer-tangle move, proved for arbitrary length by crossing-by-crossing induction",
            "distinctness":"Parallel-edge tree counts give strictly positive determinant slopes1776/1680 in N, so each family contains infinitely many distinct knot types. Family mutual disjointness unclaimed",
            "independent_replay":"Pointwise source/expanded meridian words, fixed-block/transitive-block checks and stabilizer commutators; BFS outside arcs; corner regions/direct tree enumeration; separate rational Fox determinants on12PDs;4corruptions reject",
            "all_N_assurance":"Written uniform graph/Hurwitz/propagation proofs required; no extrapolation from finite samples. Unformalized and externally unreviewed",
            "source_additions":["STOIM07: published determinant/tree lemma and explicit6star arborescence warning","JP25: classical finite-cycle and periodic twist construction in preprint Section3; published metadata checked, full published proof comparison open"],
            "catalogue_distinct_sources":37,"new_supplementary_snapshots":2,"frozen_main_manifest_unchanged":True,
            "worldwide_novelty_or_prior_family_exclusion":False,"full_MRC_proved":False
        },
        "simple_cover_additivity_progress": {
            "uniform_invariants":"C(K)=max(r-1+rank pi1Y) and A_k(K)=max(r-1+b1(Y;k)) over all connected simple transposition covers; well-defined finite maxima",
            "all_knots_result":"Both C and fixed-field A_k are(-1)additive under connected sum. Every restricted coloring has one active support by knot normal generation; converse compatible supports realize individual maxima",
            "topology":"Y=#Y_i #^beta(S1xS2), beta=sum(r_i-2)-(R-2), from every lifted sphere and capped component. Grushko adds ordinary ranks; H1 adds directly",
            "sharpness":"C=mr implies conditional meridional additivity. With matching bridge certificates, sums stay sharp. A sum of t periodic-family members and s trefoils has mr=bridge=3t+s+1",
            "defects":"Schubert bridge additivity makes nonnegative b-C and b-A_k additive; C=b is closed under connected sum and taking summands",
            "strict_boundary":sums["T47_boundary"],
            "unbounded_gap":"For J_n=#^nT47, C=2n+1 and A_k=1+n[char7], while known graph-manifold MRC and Schubert give mr=bridge=3n+1. No S3 coloring exists, excluding only rank-two-form original criterion in any dimension",
            "controls":sums["controls"],"rejected_corruptions":4,
            "verification_boundary":"Eight finite PD/word/Smith controls and separate torus matrices check inputs. Uniform all-colorings topology and Grushko step remain written classical arguments, unformalized and externally unreviewed",
            "source_additions":["STALL71 primary Grushko-Wagner statement and official1965 proof metadata; original full proof not audited","SCHULT03 primary bridge-additivity theorem and connected-sum argument, Schubert credit retained"],
            "catalogue_distinct_sources":39,"new_supplementary_snapshots":4,
            "general_meridional_additivity_or_full_MRC_proved":False,"worldwide_novelty_claim":False
        },
        "permutation_index_progress": {
            "uniform_bound":"A transitive degree-r quotient with meridian index delta gives mr>=ceil((r-1)/delta), by connected transposition-edge counting",
            "known_full_A8_control":"X=(0751)(2643),Y=(01234)(567),X^-1Y^4=(06)(25),image order20160; same quotient gives sharp known T(4,60j+15) for every j>=0",
            "class_boundary":permutation["class_boundary"],
            "negative_scope":"Exactly original16 PDs, transitive double-transposition actions in degrees8/9; only such degrees where elementary graph lower reaches4 under certified upper4",
            "degrees":permutation["degrees"],"controls":permutation["controls"],
            "raw_tuple_coverage":"Full centralizer orbit marking is disjoint and covers every normalized tail; no subgroup classification assumption",
            "independent_replay":"Permutation filtering, geometric arcs, one-pass upper trace and direct pointwise crossing tests; source parser reused only for old index names",
            "rejected_corruptions":permutation["rejected_corruptions"],
            "catalogue_at_checkpoint":40,"new_primary_snapshot":"ATLASA8; identification/order/list read; class counts proved independently",
            "this_permutation_gate_original16_targets_closed":0,"worldwide_novelty_claim":False,"full_MRC_proved":False
        },
        "ternary_affine_cover_progress": {
            "base_scope":"Both nonidentity rank-one GL3(F3) meridian classes on the unchanged original16;104 order3 transvections and117 order2 reflections",
            "tuple_orbits":[21632,34084],"raw_normalized_tuples":[1124864,1601613],
            "affine_boundary":"All admissible lifts of every rank-one base are principal; trivial linear bases give cyclic translations. No ternary affine3 transitive27point or homogeneous4 zero-common-fixed-space certificate on the16",
            "order3_result":"Only one constant cyclic base on every target, root/principal dimension1",
            "reflection_result":"Two fullGL3(F3) representation orbits on each first8 targets with root/principal dimension3; one constant root-dimension1 base on every target. Representation orbits not claimed distinct quotient kernels",
            "cover_quotients":"Each fullGL3 image induces a fullSL3(F3) double-cover quotient by determinant parity, order5616",
            "cover_modules":"All32 natural/dual3 ordinary F3 cover modules have Z1=4,H1=1; no root/meridian-extension condition imposed; threshold7 not reached",
            "known_controls":"Explicit order3/order2 homogeneous rank-one matrices for the known triple-trefoil sum, full functional rank4, nine-crossing PD/four-seed upper, root dimension4, affine image orders151632/648",
            "projective_index13_control":ternary["projective_index13_control"],
            "verification":"All3^9 matrix filtering, complete actual raw-orbit marking, geometric signed traces, affine-point coefficient reconstruction, exact kernels/principal columns and permutation image orders;4 independent replays reject10 corruptions",
            "uniform_assurance":"Written cocycle/fixed-space/core-group and Kurosh arguments separate from finite inputs; unformalized and externally unreviewed",
            "new_target_closures":0,"all_fields_or_all_SL3_modules_excluded":False,"target13point_homology_at_original_checkpoint":False,
            "current_target13point_homology_record":"MRC-SL3-PROJECTIVE13-HOMOLOGY-20261117",
            "full_MRC_proved":False,"worldwide_novelty_claim":False
        },
        "projective13_homology_progress": {
            "uniform_all_characteristic_bound":"For every finite-index r>=2 subgroup H of a finitely generated G and every field, rankG>=1+ceil((b1H-b1G)/(r-1)); cocycle sum map has rank<=b1G, so no splitting is required",
            "prime_restriction_boundary":"Only the earlier split cohomology identity needs char not dividingr; gain>=r rank3 implication is valid at every characteristic",
            "positive_degree_characteristic_control":degree_control,
            "uniform_known_control_family":"G0=Cp*Cp*Cp, odd primep, firstfactorpcycle and otherfactorsidentity; direct Schreier rewriting H0=*2p Cp gives gain2p-3>=p atp and augmentationZ1=3(p-1), one more than split prediction",
            "source_actions":32,"distinct_labelled_covers":16,"knots":8,
            "integral_homology":"All H1=Z^2+C2+C882878 in block391 or Z^2+C2+C3128602 in block401",
            "rational_first_Betti":2,"maximum_Betti_gain_except13":4,
            "augmentation_maximum_Z1_all_fields":16,"augmentation_Z1_characteristic13":14,"required_Z1":25,
            "quotient_kernel_correction":"Among the retained2 representations per source, one faithfulSL3 quotient kernel; point/dualGsets remain inequivalent, and no cross-knot pooling or all-quotient completeness claim",
            "verification":"Independent plane actions, paired group5616 and stabilizer432, geometric core words, Fox/Smith replay and13coordinate augmentation kernels;4 main+2 nonsplit-control corruptions reject",
            "structural_source_cues":projective_sources["groups"],"common_cover_or_mutation_proved":True,
            "structural_checkpoint_supersession":"Original source-cue audit proves no mutation; subsequent exact mutation theorem below proves common covers for exactly these2blocks",
            "this_homology_test_target_closures":0,"actual_cover_ranks_determined_by_this_homology_test":False,"all_SL3_modules_or_subgroups_excluded":False,
            "full_MRC_proved":False,"worldwide_novelty_claim":False,
            "uniform_assurance":"All-field inequality and every-prime control require written cochain/Schreier proofs; finite controls check inputs, not all groups. Unformalized and externally unreviewed"
        },
        "mutation_common_cover_progress": {
            "blocks":mutation["blocks"],"source_diagrams":8,"root_to_member_pairs":6,
            "maximum_Whitney_switches_per_pair":1,"cycle_vectors_per_named_graph":128,
            "edge_subsets_per_graph":8192,"positive_controls":2,"corruptions_rejected":4,
            "primary_imports":"Greene Cor4.5 with classical graph/plane-drawing inputs; Morton Thm3 torus-lift proof credited to Viro1976",
            "conclusion":"Actual Conway mutation and common double covers within each4knot block;8 ordinary cover-rank questions reduce to2",
            "abelianizations":[391,401],"rank_interval_at_mutation_checkpoint":[2,3],"rank3_proved_at_mutation_checkpoint":False,"current_shared_rank":3,
            "conditional_knot_consequence":"Rank3 for one common group settles all4members using each retained upper4; no mutation invariance of bridge or meridional rank assumed",
            "mutation_alone_target_closures":0,"other8_original_targets_untouched_at_checkpoint":True,
            "assurance":"Independent finite graph replay; classical topology imported, literal planar mutation movie/Tietze maps not supplied. Uniform proof unformalized and externally unreviewed; no worldwide novelty claim",
            "catalogue_at_mutation_checkpoint":42
        },
        "compact_four_port_progress": {
            "source_roots":["13a_129","13a_135"],"closed_presentations":[[9,13,15],[12,15,11]],
            "exact_moves":"Single-occurrence Tietze, inverse free automorphisms and normal-closure relator replacements; all accepted moves and generator maps replay",
            "formal_pieces":"Common2generator1relator groups P(length11),Q(length7)",
            "assembly":"G_X=(P*Q*free<z>)/normal(p_e*z*q_e^-1 for4ports); right base change and17to13arc gluing give an exact algebraic identity",
            "boundary_generation":"Outside normalized words stall exactly on5full arcs0..4 in both roots; no rank2 or rank3 implication",
            "verification":"Two roots,3known algebraic controls,4fragment and2assembly reductions;5compact+4port corruptions reject; H1=391/401 and transported nonabelian group order5616",
            "graph_cues":"K4 with adjacent7edge and2edge networks; branch tree/forest19/14 or14/19, total tree count11T+13F. Conductance not certified rational-tangle slope",
            "geometry_cues":"Candidate fillings match H1 but have degenerate numerical shapes and no single-pass Regina standard recognition;4_1controlL(5,2). No peripheral/core transport or geometric classification theorem",
            "rank_interval_at_compact_checkpoint":[2,3],"current_shared_rank":3,"compact_presentation_alone_target_closures":0,"other8_original_targets_untouched_at_checkpoint":True,
            "new_classical_source_imports":0,"catalogue_at_compact_checkpoint":42,
            "uniform_assurance":"Elementary algebraic assembly proof written; formal/external review absent. No injectivity/JSJ/Seifert/gluing rank theorem imported"
        },
        "four_piece_sharp_progress": {
            "new_original_target_equalities":sorted(new_equalities),"double_cover_group_rank":3,"meridional_rank_and_bridge":4,
            "exact_models":"P andQ each split into D(2,3) and A(2) with explicit primitive peripheral sections; ports reduce to an injective torus amalgam",
            "graph_manifold_JSJ_chain":["D(2,3)","A(2)","A(2)","D(2,3)"],
            "edge_fiber_intersections":[1,1,1],"unimodular_edge_maps":3,
            "topological_gate":"Constructed model closed orientable irreducible; negative-base fiber uniqueness makes all3tori canonical. Exact group isomorphism suffices; no source-cover homeomorphism from numerical triangulations assumed",
            "classical_rank_input":"Boileau-Weidmann2005 Theorem2 necessary list allows at most3JSJpieces for rank2, hence model rank>=3; exact3generator upper gives equality",
            "finite_verification":"Independent cyclic-amalgam normal forms, free-basis/port identities, slope determinants and8individual upper4 traces;5further corruptions reject",
            "uniform_assurance":"Seifert/JSJ and classification are written/imported classical arguments, unformalized and externally unreviewed; finite receipts certify inputs only",
            "remaining_original_targets_at_checkpoint":["13a_233","13a_234","13a_239","13a_242","13a_243","13a_244","13a_253","13a_265"],
            "panel_names_certified_at_checkpoint":12,"panel_names_without_certificate_at_checkpoint":281,
            "catalogue_at_checkpoint":45,"current_refinement":"Millichap-Worden2016 Cor5.4 removes original irregular-cover type; context read, full hidden-symmetry proof not audited or needed for eight-case bound","worldwide_previously_unresolved_status_at_checkpoint":"UNESTABLISHED","full_MRC_proved":False
        },
        "three_piece_sharp_progress": {
            "new_original_target_equalities":[r["name"] for r in three_piece["roots"]],"double_cover_group_rank":3,"meridional_rank_and_bridge":4,
            "graph_manifold_JSJ_chain":"D(2,3)-A(m)-D(2,3),m=2,3,2","central_filling_slopes":[r["central_filling_slope_a_f"] for r in three_piece["roots"]],
            "rank_gate":"BW05 Theorem2 case6 is the only3piece type, requires central2bridge link meridians. Filling the actual attached fiber has nonabelian C_m*C_abs(r) quotient ontoS3, impossible for the solid torus from a2bridge-component meridian filling",
            "finite_verification":"Five full core and14fragment reductions; physical arc classes/all balanced4port cuts/17to13joins, complete cyclic-amalgam normal forms, primitive slopes and3S3quotients. Eight uppers replay;5corruptions reject",
            "failed_transfers":three_piece["unproved_transfers"],"actual_rank_transfers":0,
            "remaining_original_targets_at_checkpoint":["13a_234","13a_242","13a_243","13a_244","13a_265"],
            "original_targets_certified_at_checkpoint":11,"panel_names_certified_at_checkpoint":15,"panel_names_without_certificate_at_checkpoint":278,
            "catalogue_at_checkpoint":45,"uniform_assurance":"Written Seifert/JSJ/classification/filling argument, explicit classical imports; unformalized and externally unreviewed. No source-cover homeomorphism from numerical geometry assumed",
            "worldwide_previously_unresolved_status_at_checkpoint":"UNESTABLISHED","full_MRC_proved":False
        },
        "original_sixteen_completion_progress": {
            "new_original_target_equalities":[r["name"] for r in final_five["targets"]],
            "all_original_target_equalities":original_equalities,"original_targets_certified_here":16,"remaining_original_targets":[],
            "double_cover_group_rank":3,"meridional_rank_and_bridge":4,
            "new_exact_models":"Five direct D(2,3)-A(m)-D(2,3) source group models; Q_R(2/3),Q_L,Q_S root presentations, common conjugations, ks/sk outer base paths and inner fiber inversion at265",
            "central_filling_slopes":[r['central_filling_slope_c_f'] for r in final_five['targets']],
            "finite_verification":"Eight source/full and32fragment reductions, all balanced4port cuts and17to13joins, complete normal forms, primitive gluing matrices,3old controls,5actual S3 filling images,5individual uppers and6corruptions",
            "rank_gate":"Actual nonabelian central filling excludes the only three-piece BW05 rank2 JSJ type; three pieces alone do not suffice. Written classical-import realization/JSJ and exact source upper3 required",
            "cross_knot_rank_transfers":0,"failed_transfer_checkpoint_retained":True,
            "panel_names_certified_here":20,"panel_names_without_certificate_here":273,"catalogue_at_checkpoint":45,
            "uniform_assurance":"Explicit working proofs and classical imports; unformalized and externally unreviewed. No numerical source-cover identification or rank=genus assumption",
            "worldwide_previously_unresolved_status_at_checkpoint":"UNESTABLISHED","current_known_family_status":"ESTABLISHED BJW18 graph-cover family by closed rigidity","exact_named_computational_priority":"INCOMPLETE","full_MRC_proved":False
        },
        "known_graph_cover_family_progress": {
            "original_targets_known_family":16,"new_cases_outside_known_MRC_families":0,
            "actual_double_covers":"Closed prime graph manifolds homeomorphic up to orientation to the exact constructed models",
            "rigidity_input":"AFW15 final-author Theorems2.1.2-2.1.3, arbitrary closed classification; infinite prime model removes lens/prime-factor exceptions and proves source primeness",
            "existing_MRC_family":"BJW18 published Theorem1.2 and Cor1.3, actual prime graph-manifold double branched cover and exact bridge4",
            "upper_four_scope":"Known meridional ranks1/2 and rank3 graph-cover result already imply MRC for this class with b<=4; exact cover rank3 separately determines b=mr4",
            "finite_binding":"All16literal named PDs and individual upper4 traces, retained exact model/lower seals and first8mutation aliases; topology not proved by finite checker",
            "unchanged_exact_equalities":20,"other_sharp_seed_family_audit_names":family['other_sharp_seed_family_audit_names'],
            "exact_named_bridge_rank_priority":"UNESTABLISHED","catalogue_current_sources":46,
            "source_audit":"Final AFW author pp30-31/PDF34-35 visually checked and v3 scopes compared; BJW publishedpp61-64 reread. Full constituent proofs and published-version comparison remain imported",
            "provenance_refresh":"29cheap already frozen source/provenance/replay recipes; no broad293/base/module campaign",
            "uniform_assurance":"Written classical-rigidity consequence, unformalized and externally unreviewed. No branch-involution transport, boundary-manifold shortcut or numerical homeomorphism assumed",
            "full_MRC_proved":False
        },
        "PDF_visual_review": "Earlier proofs retained. New actual-cover/known-family scope proof and updated literature map compiled twice, rendered and final previews inspected; relevant AFW/BJW primary pages visually checked; overfull boxes and unresolved references checked.",
        "repository_navigation": {
            "render_navigation": "FAILED on pre-existing registry validation error (see exact_failure)",
            "check_navigation": "FAILED on the same pre-existing registry validation error",
            "exact_failure": integration["inherited_registry_failure"],
            "earlier_prior_work_search_failure": "METHOD-EC-CENTRE-POLICY missing artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json",
            "response": "Retained failures; did not reconstruct unrelated missing artifacts or hand-edit generated status. Scoped new records and source hashes checked."
        },
        "next_gate": "All16original targets are established BJW18 known-family cases by actual-cover closed rigidity; their exact rank3/b4 certificates stand, but no outside-known-family MRC case follows. Next audit the4earlier sharp simple3cover seeds3145,3147,3148,3149 for actual double-cover/other known-family scope under a small frozen plan, retaining the16as known controls. For prospective new-MRC targets beyond BJW18, justify actual double-cover topology outside the graph-manifold class and check other known families; hyperbolic knot exterior or non-arborescence alone is insufficient. A rank>=3 quotient may give cover lower4 without actual graph-cover identification, but such an outside-family construction remains prospective. Keep20named equalities and273panel names without certificates here separate from worldwide novelty. Preserve failures and do not rerun original16/293/base/module campaigns as discovery. Continue original/recent full-source, named-priority and fullMRC work.",
        "live_background_process": None
    }
    report["completion_audit"]["exhaustive_literature_and_source_audit"] = report["completion_audit"]["exhaustive_literature_and_source_audit"].replace("broad35-source map and18 primary snapshots", "broad46-reference map,18 frozen main snapshots and separately sealed supplementary source reads; original16known-family membership now established through AFW15/BJW18")
    report["completion_audit"]["construct_representations_and_matching_upper_bounds"] += " The four seeds now also give four explicit infinite1+6N twist families, with finite75point block actions, uniform4seed upper induction and strictly growing determinants. The uniform proof is distinct from the12 finite input samples."
    (HERE / "progress.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS scoped evidence seals", "goal_complete": False,
                      "control_counts": report["verified_controls"],
                      "navigation_failure_retained": integration["inherited_registry_failure"]}, indent=2))


if __name__ == "__main__":
    main()
