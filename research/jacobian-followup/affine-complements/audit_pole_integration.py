"""Scoped adoption audit; does not replace the repository-wide navigation checks."""
from pathlib import Path
import gzip,hashlib,json,re,sys,py_compile,tempfile

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint,validate_lessons

checks=[]
def check(name,value):
    assert value,name
    checks.append(name)

out=HERE/"pole_validation.json"
# Make the result target available while checking its documentation link.
if not out.exists():
    out.write_text('{"focused_outcome":"RUNNING"}\n')

index=json.loads((ROOT/"MATH_STATUS.json").read_text())
ids=[x["id"] for x in index["entries"]]
check("unique mathematical-status IDs",len(ids)==len(set(ids)))
records={x["id"]:x for x in index["entries"] if x["id"] in {"PF2RULED1","PF2RULED2","PF2RULED3","PF2RULED4","PF2RULED5","PF2RULED6","PF2RULED7","PF2RULED8","PF2RULED9","PF2RULED10","PF2RULED11","PF2RULED12","PF2RULED13","PF2RULED14","PF2RULED15","PF2RULED16","PF2RULED17","PF2RULED18","PF2RULED19","PF2RULED20","PF2RULED21","PF2RULED22","PF2RULED23","PF2RULED24","PF2RULED25","PF2RULED26","PF2RULED27","PF2RULED28","PF2RULED29","PF2RULED30","PF2RULED31","PF2RULED32","PF2RULED33","PF2RULED34","PF2RULED35","PF2RULED36","PF2RULED37","PF2RULED38","PF2RULED39","PF2RULED40","PF2RULED41","PF2RULED42","PF2RULED43","PF2RULED44","PF2RULED45","PF2RULED46","PF2RULED47","PF2RULED48","PF2RULED49","PF2RULED50","PF2RULED51","PF2RULED52","PF2RULED53","PF2RULED54","PF2RULED55","PF2RULED56","PF2RULED57","PF2RULED58","PF2RULED59","PF2RULED60","PF2RULED61","PF2RULED62","PF2RULED63","PF2RULED64","PF2RULED65","PF2RULED66","PF2RULED67","PF2RULED68","PF2RULED69","PF2RULED70","PF2RULED71","PF2RULED72","PF2RULED73","PF2RULED74","PF2RULED75","PF2RULED76","PF2RULED77","PF2RULED78","PF2RULED79","PF2RULED80","PF2RULED81","PF2RULED82","PF2RULED83","PF2RULED84","PF2RULED85","PF2RULED86","PF2RULED87","PF2RULED88","PF2RULED89","PF2RULED90","PF2RULED91","PF2RULED92","PF2RULED93","PF2RULED94","PF2RULED95","PF2RULED96","PF2RULED97","PF2RULED98","PF2RULED99","PF2RULED100","PF2RULED101","PF2RULED102","PF2RULED103","PF2RULED104","PF2RULED105","PF2RULED106","PF2RULED107","PF2RULED108","PF2RULED109","PF2RULED110","PF2RULED111","PF2RULED112","PF2RULED113","PF2RULED114","PF2RULED115","PF2RULED116","PF2RULED117","PF2RULED118","PF2RULED119","PF2RULED120","PF2RULED121","PF2RULED122","PF2RULED123","PF2RULED124","PF2RULED125","PF2RULED126","PF2RULED127","PF2RULED128","PF2RULED129","PF2RULED130","PF2RULED131","PF2RULED132","PF2RULED133","PF2RULED134","PF2RULED135","PF2RULED136","PF2RULED137","PF2RULED138","PF2RULED139","PF2RULED140","PF2RULED141","PF2RULED142","PF2COORD1","PF2FACTOR1","PF2PRIM1","PF2PERIOD1","PF2WRIGHT1","PF2GRADE1","PF2TRACE1"}}
check("all one hundred forty-nine scoped claims registered",len(records)==149)
for id_,entry in records.items():
    check("schema "+id_,REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS)
    check("programme isolation "+id_,entry["programme_status"]=="archived")
    check("dependencies "+id_,all(d in ids or d.startswith("external: ") for d in entry["dependencies"]))
    check("checker source hash "+id_,
          entry["artifact_hash"]=="sha256:"+hashlib.sha256((ROOT/entry["checker"]).read_bytes()).hexdigest())
    for path in [entry["canonical_source"],entry["checker"]]+entry["software_lock"]:
        check("present "+id_+" "+path,(ROOT/path).is_file())

lessons=json.loads((ROOT/"knowledge/lessons.json").read_text())["lessons"]
selected=[x for x in lessons if x["id"] in {
 "METHOD-AFFINE-COMPLEMENT-RULED","METHOD-RULED-TARGET-POLE-REDUCTION",
 "METHOD-PRIME-RAMIFICATION-POLYNOMIAL-RECEIVER",
 "METHOD-GENERIC-FIBRE-POLE-CENSUS","METHOD-PRIME-RAMIFICATION-FIELD-TOWER",
 "METHOD-EXPLICIT-COMPLEMENT-DEGREE-TRANSFER","METHOD-GLOBAL-VOLUME-EXACTNESS","METHOD-GRADING-DESCENT-TRACE-BOUNDARY","METHOD-WEIGHTED-LOCAL-NONCUBE","METHOD-FORMAL-CENTRALIZER-JETS","METHOD-BOUNDARY-LEADING-PERIOD","METHOD-FIRST-NEWTON-EDGE-POLE-BUDGET","METHOD-CUBE-CLUSTER-DISCRIMINANT-BUDGET","METHOD-CUSP-JET-LEADING-POWER","METHOD-DOUBLE-FILTRATION-CUSP-JETS","METHOD-PRIMITIVE-CHART-LEADING-FORM","METHOD-RESONANT-CHART-FRONTIER","METHOD-KUMMER-BOUNDARY-PERIOD","METHOD-GENERIC-POLAR-CUSP-BUDGET","METHOD-GENERIC-QUARTIC-POLAR-GAP","METHOD-GENERIC-CRITICAL-VALUE-UNIBRANCH","METHOD-COMPLETE-CUSP-JET-LIFTING","METHOD-HIGHEST-JACOBIAN-KERNEL-BLOCK","METHOD-AFFINE-JACOBIAN-FIBER-SURJECTIVITY","METHOD-UNIVERSAL-KERNEL-DETERMINANT-OBSTRUCTION","METHOD-SQUARE-BOUNDARY-CRITICAL-REMAINDER","METHOD-GLOBAL-POLYNOMIAL-GRAPH-ELIMINATION","METHOD-UNIT-ROW-GLOBAL-CONSTANT-RECOVERY","METHOD-QUADRATIC-COEFFICIENT-UNFOLDING","METHOD-SIGN-FIXED-FIELD-CERTIFICATE-DESCENT","METHOD-FINITE-REDUCTION-CIRCUIT-INVERSION","METHOD-GLOBAL-EXCEPTIONAL-CURVE-GATE","METHOD-RESIDUAL-DISCRIMINANT-COVER","METHOD-MARKED-END-SEIFERT-OBSTRUCTION","METHOD-QUARTIC-MARKED-CABLING-DESCENT","METHOD-INTRINSIC-ONE-BOUNDARY-END-GATE","METHOD-CYCLIC-MULTIPLE-FIBER-VOLUME","METHOD-COORDINATE-POLE-REGULAR-SECTOR","METHOD-ODD-LENS-CORE-QUARTIC-DESCENT","METHOD-REGULAR-ROOT-CYCLIC-CHART-RESIDUE","METHOD-ATOROIDAL-LENS-QUARTIC-BOUNDARY","METHOD-COPRIME-REGULAR-SEIFERT-END","METHOD-UNIVERSAL-REGULAR-AFFINE-END","METHOD-REGULAR-QUARTIC-JSJ-FIBER-DESCENT","METHOD-CANONICAL-KEYFORM-EXACT-COORDINATE","METHOD-UNRAMIFIED-PACKET-OREVKOV-QUOTIENT","METHOD-LEFSCHETZ-MERIDIAN-CENTRALIZER","METHOD-UNRAMIFIED-SHEET-PAIRWISE-COMPLETE-CURVE","METHOD-TWO-SIDED-NORM-BAND-CUT","METHOD-GLOBAL-COORDINATE-PENCIL-MATRIX","METHOD-BOUNDARY-FIBRE-SCHUR-GRADED-MODELS","METHOD-SOURCE-FIELD-SYMMETRY-ESCAPE","METHOD-TRACE-DUAL-PENCIL-EXCESS","METHOD-QUARTIC-FOX-LIFT-SIGN-TRIPLE","METHOD-SIGN-PICARD-SQUARE-LINKING","METHOD-POWER-DIFFERENCE-NONMONOGENIC-RESOLVENT","METHOD-CANONICAL-CONTENT-BINARY-CONDUCTOR","METHOD-CANONICAL-QUADRICS-FULL-SOURCE-RESIDUAL","METHOD-NODAL-DEFORMATION-NONPOSITIVE-NORMAL","METHOD-PERIPHERAL-LONGITUDE-ZERO-NORMAL-CONTROL","METHOD-FULL-SOURCE-END-LATTICE-MARKING","METHOD-MONODROMY-ORBIT-TREE-EXACT-PERIOD","METHOD-ALL-BINOMIAL-PERIODIC-GATE-EXISTENCE","METHOD-ODD-DIHEDRAL-FORMAL-BRANCH-INDEX","METHOD-CYCLIC-COMPANION-MARKED-END-EXCLUSION","METHOD-ORIENTED-LINKING-COMMON-SOURCE-GATE","METHOD-COMMON-PENCIL-GRAM-NONZERO-MARK","METHOD-CYCLIC-COMPANION-PARITY-FREE-MARK","METHOD-SEIFERT-FIBRE-ORBITS-AND-BASE-GENUS","METHOD-PERIPHERAL-CONE-TORSION-INDIVIDUAL-MARK-BUDGET","METHOD-FOLD-STRICT-HALF-LINKING-AND-INTEGRAL-S4-PIECES","METHOD-EXACT-SEIFERT-FIBRE-ORDER-UNIT-TRANSFER","METHOD-LIFTED-SURFACE-BALANCES-AND-FULL-END-GROUP","METHOD-RATIONAL-COMMON-END-NO-ORDERABLE-IMAGE","METHOD-ACTUAL-SEIFERT-SECTIONS-WINDING-COLLAPSE","METHOD-MATCHED-MERIDIAN-ENDPOINT-COMMON-CLOSURE","METHOD-COMPLEMENTARY-CORNER-ROOTED-JSJ-ENDPOINT","METHOD-COMPLEMENTARY-FIBRE-FUSION-AND-TREFOIL-CAPS","METHOD-TERMINAL-CONE-MARKED-CAPS-BINOMIAL-CLASS-CLOSURE","METHOD-GENERAL-SPLICE-TRANSFER-ODD-FIBRE-PERIOD","METHOD-MERIDIAN-BALANCE-LSPACE-GENERAL-TWO-STAGE-CLOSURE","METHOD-ITERATED-PERIPHERAL-STATES-AND-ABSORBING-TAILS","METHOD-PERSISTENT-ROOTED-JSJ-AND-TERMINAL-MODULUS","METHOD-CYCLIC-LOOP-GEOMETRY-AND-FIRST-ARRIVAL","METHOD-NORMALIZED-TERMINAL-CONE-TWISTS","METHOD-MERIDIAN-CYCLE-FRAMING-AND-S4-INCOMING","METHOD-CYCLIC-PATTERN-MARK-TRANSFER-AND-A4-REDUCTION","METHOD-UNIMODULAR-SOLID-PREFIX-AND-INITIAL-S4-CLOSURE","METHOD-NONUNIMODULAR-A4-CAPS-AND-THREE-TERMINAL-PAIRS","METHOD-UNFILLED-CYCLIC-PERIPHERY-AND-S3-SOLID-HISTORY","METHOD-GENUS-FRAMING-AND-FINAL-S3-CLOSURE","METHOD-RECURSIVE-PERIPHERAL-LINKING-AND-ENDPOINT-TAILS","METHOD-INTEGER-SLOPE-GRID-AND-COMMON-ITERATED-CLOSURE","METHOD-SOLID-CAP-GENUS-AND-CYCLIC-PRIME-MARK","METHOD-TERMINAL-BRIDGE-AND-S3-SOLID-A4-CLOSURE","METHOD-PREFIX-EXCESS-AND-GENUINE-ITERATED-QUARTIC-CLOSURE","METHOD-UNIT-PERIPHERAL-REDUCTION-AND-TARGET-COVERAGE","METHOD-PROPER-FINITE-MARKED-END-AND-PERIPHERAL-LATTICES","METHOD-MONIC-KEY-JACOBIAN-AND-ALL-INDEX-CHART-NECESSITY","METHOD-POLYNOMIAL-KEY-CHAIN-AND-LAURENT-TERMINAL-NECESSITY","METHOD-FULL-PERIPHERAL-RATIONAL-LINKING-AND-LONGITUDE-JOINING-RECURSION","METHOD-PRIMITIVE-END-CANONICAL-CONTACT-AND-TRIPLE-TORUS-EXCLUSION","METHOD-MULTIPLE-BOUNDARY-CANONICAL-CONTACT-AND-A-TYPE-END-EXCLUSION","METHOD-FIBRE-TRANSITIVE-SHEAR-RECURSION-AND-SPIN-BUDGET-EXCLUSION","METHOD-TRIPLE-BOUNDARY-RESIDUE-ONE-AND-NORMAL-ITERATED-EXCLUSION","METHOD-CYCLIC-QUOTIENT-ORDER-CONDUCTOR-AND-MORSE-RESIDUE-BUDGET","METHOD-COPRIME-LEADING-FORMS-DEGENERATE-RESIDUE-BUDGET","METHOD-ALL-ADE-END-CANONICAL-ENERGY-EXCLUSION","METHOD-NORMAL-RESOLUTION-CHERN-GATE-AND-ADE-END-EXCLUSION","METHOD-PLANAR-CHERN-OBSTRUCTION-AND-CIRCLE-BUNDLE-SOURCE-REDUCTION","METHOD-WHOLE-PLANE-THOM-LATTICE-AND-COMPLETE-PERIPHERAL-OBSTRUCTION","METHOD-HIGH-INDEX-COLLISION-EXCESS-AND-NEAR-MAXIMAL-EXCLUSION","METHOD-LOCAL-MONODROMY-CAPACITY-AND-WEIGHTED-FULL-BRANCH-EXCLUSION","METHOD-TORUS-CUSP-SOLID-TORUS-EXCESS-AND-PROPORTIONAL-INDEX-BOUND","METHOD-RESIDUAL-MILNOR-IDENTITY-AND-MINIMAL-EXCESS","METHOD-LOGARITHMIC-PULLBACK-AND-SHARP-DESSIN-CUSP-COVERS","METHOD-NILRADICAL-TRACE-RESTRICTION-AND-CYCLIC-DESCENT","METHOD-HALF-SUPPORT-INERTIA-AND-PLANE-ACTION-EXCLUSION","METHOD-SIMPLE-RAMIFICATION-MULTIPLICITY-AND-HALF-SUPPORT-CUSPS","METHOD-PRIMITIVE-PRODUCT-ACTION-EXCLUSION-AND-PRIME-REDUCTION","METHOD-DISJOINT-BRANCH-SUPPORTS-AND-PRIMITIVE-REDUCTION","METHOD-HEREDITARY-BLOCK-QUOTIENT-MONODROMY-REDUCTION","METHOD-INTERMEDIATE-CRITICAL-AXIS-DESCENT-AND-AFFINE-CUSP-PROFILES","METHOD-INTERMEDIATE-CHERN-BETTI-AND-END-OBSTRUCTIONS","METHOD-BINARY-AFFINE-LONGITUDE-AND-FILLING-PROFILES","METHOD-SEIFERT-FRAMING-CANCELLATION-AND-RESIDUAL-TWO-COVERS","METHOD-ALL-INDEX-QUOTIENT-CORRECTED-PERIPHERAL-SHEAR","METHOD-MINIMAL-DEGREE-PSEUDOPLANE-COMPRESSION-AND-DECK-DICHOTOMY","METHOD-HEREDITARY-EULER-EXCESS-AND-QUOTIENT-DEGREE-BOUNDS","METHOD-RELATIVE-EULER-DEFECT-AND-RULED-CHART-COMPRESSION","METHOD-PROPER-BINARY-LINEAR-GROUPS-AND-TRIVIAL-LONGITUDE","METHOD-LOW-EULER-CYCLICITY-AND-DECORATED-RADICAL-TOWER-EXCLUSION","METHOD-DECK-NORMALIZED-EULER-AND-UNRESTRICTED-RADICAL-TOWER-EXCLUSION","METHOD-FINITE-ORBIFOLD-GROUPS-AND-UNIVERSAL-COMPLETE-RULINGS","METHOD-COTANGENT-DECK-NORMAL-FORM-AND-RELATIVE-TRACE","METHOD-COMMON-OMITTED-ORBITS-AND-ALMOST-SURJECTIVE-DECK-QUOTIENTS","METHOD-SEPARATED-VARIABLE-VOLUME-AND-SINGLE-ROOT-OBSTRUCTION","METHOD-TRANSLATION-ATLAS-EXACTNESS-AND-UNBOUNDED-MOVING-LABEL-FIBRES","METHOD-REPLICA-EQUIVARIANT-DEGREE-DESCENT-AND-SQUARE-FREE-OBSTRUCTION","METHOD-COMPLETELY-MISSED-INTERMEDIATE-LINES-AND-SOURCE-IMAGE-EULER-BUDGET","METHOD-GENERIC-SOURCE-VELOCITY-MULTIPLIERS-AND-TRACE-POLE-CHARGES","METHOD-PARTIAL-KUMMER-PLANE-CHARTS-AND-CYCLIC-EQUIVARIANT-FACTORS"}]
validate_lessons(selected,index["entries"],check_reviews=True)
check("all one hundred forty-two related lessons pass shared schema/path/scope validator",len(selected)==142)
work=json.loads((ROOT/"knowledge/work_items.json").read_text())["actions"]
action=next(x for x in work if x["id"]=="WORK-QUARTIC-RULED-COMPATIBILITY")
check("work gate preserves all one hundred forty-nine claims",set(records)<=set(action["claims"]))
check("work gate scope review",{id_:action["reviewed_claims"].get(id_) for id_ in records}=={
 id_:claim_fingerprint(e) for id_,e in records.items()})
check("work gate sources exist",all((ROOT/x).is_file() for x in action["sources"]))
all_records={entry["id"]:entry for entry in index["entries"]}
check("shared work gate retains current reviews for additional registered claims",
      set(action["claims"])==set(action["reviewed_claims"]) and
      all(id_ in all_records and action["reviewed_claims"][id_]==claim_fingerprint(all_records[id_])
          for id_ in action["claims"]))


for document in sorted(HERE.glob("*.md"))+[ROOT/"jacobian-followup/README.md"]:
    for target in re.findall(r"\[[^]]+\]\(([^)]+)\)",document.read_text()):
        if "://" in target or target.startswith(("#","mailto:")):
            continue
        check("link "+document.name+" -> "+target,
              (document.parent/target.split("#",1)[0]).resolve().exists())
with tempfile.TemporaryDirectory() as temp:
    for source in sorted(HERE.glob("*.py")):
        py_compile.compile(str(source),cfile=str(Path(temp)/(source.stem+".pyc")),doraise=True)
        check("compile "+source.name,True)

receipts={}
for filename in ("pole_reduction_receipt.json","pole_four_receipt.json","pole4_exclusion_receipt.json",
                 "low_coordinate_receipt.json","primitive_cover_receipt.json","pole35_receipt.json","source_coordinates_receipt.json","volume_period_receipt.json","wright_grading_receipt.json","pole_three_all_height_receipt.json","pole_four_companions_receipt.json","pole_four_nonsquare_receipt.json","boundary_leading_period_receipt.json","pole58_receipt.json","pole68_cube_receipt.json","cusp_leading_receipt.json","contact_filtration_receipt.json","height_eight_receipt.json","height_nine_frontier_receipt.json","cubic_cover_period_receipt.json","polar_cusp_sector_receipt.json","quartic_polar_gap_receipt.json","unibranch_polar33_receipt.json","full_cusp_lift_receipt.json","cusp_high_jacobian_receipt.json","full_next_normal_receipt.json","second_normal_obstruction_receipt.json","square_boundary_gap_receipt.json","degree8_affine_reduction_receipt.json","degree8_unit_row_receipt.json","degree8_quadratic_receipt.json","degree8_quartic_descent_receipt.json","degree8_linear_remainder_receipt.json","global_boundary_receipt.json","residual_cover_receipt.json","torus_infinity_receipt.json","quartic_cabling_receipt.json","intrinsic_normalization_receipt.json","multiple_fiber_cyclic_receipt.json","coordinate_regular_sector_receipt.json","regular_sector_end_receipt.json","higher_regular_coordinate_receipt.json","borderline_regular_end_receipt.json","coprime_regular_end_receipt.json","universal_regular_end_receipt.json","regular_quartic_exclusion_receipt.json","exact_coordinate_bridge_receipt.json","unramified_orevkov_receipt.json","lefschetz_centralizer_receipt.json","unramified_sheet_pair_receipt.json","two_sided_norm_receipt.json","coordinate_pencil_matrix_receipt.json","boundary_fibre_models_receipt.json","source_symmetry_gate_receipt.json","trace_dual_pencil_receipt.json","quartic_colouring_receipt.json","sign_square_receipt.json","power_difference_receipt.json","canonical_formats_receipt.json","canonical_quadrics_receipt.json","nodal_deformation_receipt.json","peripheral_zero_receipt.json","marked_end_lattice_receipt.json","orbit_tree_receipt.json","binomial_periodic_receipt.json","formal_branch_index_receipt.json","cyclic_companion_receipt.json","oriented_pattern_linking_receipt.json","common_pencil_mark_receipt.json","cyclic_parity_receipt.json","seifert_base_genus_receipt.json","peripheral_cone_budget_receipt.json","fold_schur_s4_receipt.json","fibre_order_transfer_receipt.json","lifted_surface_s3_receipt.json","rational_end_orderability_receipt.json","seifert_sections_winding_three_receipt.json","matched_meridian_common_exclusion_receipt.json","complementary_corner_jsj_receipt.json","trefoil_complementary_fusion_receipt.json","terminal_cone_binomial_closure_receipt.json","general_two_stage_receipt.json","meridian_balance_two_stage_receipt.json","iterated_peripheral_tail_receipt.json","persistent_jsj_terminal_modulus_receipt.json","cyclic_loop_first_arrival_receipt.json","terminal_cone_twists_receipt.json","meridian_cycle_s4_receipt.json","cyclic_pattern_a4_receipt.json","unimodular_initial_s4_receipt.json","nonunimodular_a4_receipt.json","global_ruling_unfilled_cyclic_receipt.json","genus_final_s3_receipt.json","recursive_linking_endpoints_receipt.json","integer_grid_common_receipt.json","cyclic_solid_cap_genus_receipt.json","terminal_bridge_closure_receipt.json","prefix_excess_closure_receipt.json","unit_peripheral_coverage_receipt.json","proper_finite_marked_end_receipt.json","monic_key_chart_audit_receipt.json","polynomial_key_chain_receipt.json","full_peripheral_linking_receipt.json","primitive_end_contact_receipt.json","multiple_boundary_contact_receipt.json","fibre_transitive_spin_receipt.json","triple_boundary_residue_receipt.json","cyclic_quotient_morse_receipt.json","coprime_leading_forms_receipt.json","ade_end_exclusion_receipt.json","normal_resolution_chern_receipt.json","planar_circle_bundle_receipt.json","whole_plane_peripheral_receipt.json","high_index_budget_receipt.json","local_monodromy_branch_receipt.json","torus_cusp_excess_receipt.json","residual_milnor_receipt.json","logarithmic_pullback_receipt.json","nilradical_trace_receipt.json","half_support_inertia_receipt.json","simple_ramification_ledger_receipt.json","primitive_product_action_receipt.json","disjoint_branch_supports_receipt.json","block_quotient_monodromy_receipt.json","intermediate_critical_axis_receipt.json","intermediate_chern_betti_receipt.json","binary_affine_longitude_receipt.json","seifert_framing_cancellation_receipt.json","quotient_corrected_shear_receipt.json","pseudoplane_compression_receipt.json","hereditary_euler_excess_receipt.json","ruled_chart_compression_receipt.json","proper_binary_longitude_receipt.json","radical_tower_exclusion_receipt.json","deck_normalized_euler_receipt.json","finite_orbifold_groups_receipt.json","cotangent_deck_trace_receipt.json","common_omitted_orbits_receipt.json","separated_variable_volume_receipt.json","translation_atlas_exactness_receipt.json","replica_degree_descent_receipt.json","missed_intermediate_lines_receipt.json","source_velocity_multipliers_receipt.json","partial_kummer_charts_receipt.json"):
    payload=json.loads((HERE/filename).read_text())
    check("receipt outcome "+filename,payload["outcome"]=="PASS")
    for path,digest in payload["source_sha256"].items():
        check("receipt seal "+filename+" "+path,
              hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    if "first_open_reduced_pairs" in payload and filename not in ("height_eight_receipt.json","height_nine_frontier_receipt.json","cubic_cover_period_receipt.json","polar_cusp_sector_receipt.json","quartic_polar_gap_receipt.json","unibranch_polar33_receipt.json","full_cusp_lift_receipt.json","cusp_high_jacobian_receipt.json","full_next_normal_receipt.json","second_normal_obstruction_receipt.json","square_boundary_gap_receipt.json","degree8_affine_reduction_receipt.json","degree8_unit_row_receipt.json","degree8_quadratic_receipt.json","degree8_quartic_descent_receipt.json","degree8_linear_remainder_receipt.json","global_boundary_receipt.json","residual_cover_receipt.json","torus_infinity_receipt.json","quartic_cabling_receipt.json"):
        check("old frontier explicitly historical "+filename,
              payload["frontier_superseded_by"]=="PF2RULED16" and
              payload["frontier_handoff_scope"].startswith("historical pre-PF2RULED16"))
    receipts[filename]=payload["check_count"]
check("earlier control receipt retains its narrow assurance",
      json.loads((HERE/"pole_four_receipt.json").read_text())["assurance"]["full_pole_four_layer_excluded"] is False)
check("research goal remains active","Status: active." in (HERE/"RESEARCH_GOAL.md").read_text())
full=json.loads((HERE/"pole4_exclusion_receipt.json").read_text())
check("all five marked charts are covered in the full receipt",
      full["leading_atlas"]["mixed_charts"]==3 and full["leading_atlas"]["triple_charts"]==2)
check("actual quartic packet exclusion remains unclaimed",
      full["assurance"]["quartic_packets_excluded"] is False)
check("historical pole-four receipt retains its intrinsic higher-pole and comparison limits",
      full["assurance"]["higher_poles_excluded"] is False and
      full["assurance"]["actual_ruled_compatibility_proved"] is False)
low=json.loads((HERE/"low_coordinate_receipt.json").read_text())
prim=json.loads((HERE/"primitive_cover_receipt.json").read_text())
check("all-height coordinate bound retains open frontier",
      low["first_open_height_five_pairs"]==[])
check("coordinate argument retains geometry and full-packet assurance limits",
      low["assurance"]["complete_generic_fibre_geometry_formally_certified"] is False and
      low["assurance"]["quartic_packets_excluded"] is False)
check("finite quartic subgroup census is complete",
      prim["s4_subgroups"]==dict(total=30,transitive=9,primitive_orders=[12,24],
                                primitive_with_transposition_orders=[24]))
check("primitive field does not exclude the sign field in its closure",
      prim["assurance"]["sign_field_in_galois_closure_excluded"] is False)
check("purity pin exists and is externally imported",
      (HERE/"sources/stacks-purity-0BMB.html").is_file() and
      "external" in prim["assurance"]["purity"])
check("current full-atlas handoff omits the now-closed pole-two branch",
      full["open_pole_pairs_at_height_five"]==[])
p35=json.loads((HERE/"pole35_receipt.json").read_text())
check("complete pole-(3,5) receipt preserves the open frontier",
      p35["first_open_height_five_pairs"]==[] and p35["mixed_chart_count"]==3)
check("pole-(3,5) does not overclaim full height five or packet exclusion",
      p35["assurance"]["complete_pole_height_five_excluded"] is False and
      p35["assurance"]["quartic_packets_excluded"] is False and
      p35["assurance"]["actual_ruled_compatibility_proved"] is False)
coords=json.loads((HERE/"source_coordinates_receipt.json").read_text())
check("explicit source coordinates cover both marked contact positions",
      coords["source_charts"]==2 and coords["section_degrees"]==[16,12,8,4])
check("source transfer closes all heights through seven with the exact next frontier",
      coords["excluded_pole_height_upper_bound"]==7 and
      coords["first_open_reduced_pairs"]==[[6,8],[7,8]])
check("source degree theorem remains external and the2026 conditional upgrade unused",
      "external" in coords["assurance"]["plane_degree_theorem"] and
      coords["assurance"]["conditional_2026_degree125_upgrade_used"] is False)
check("source transfer retains actual compatibility and full-packet gaps",
      coords["assurance"]["actual_ruled_compatibility_proved"] is False and
      coords["assurance"]["quartic_packets_excluded"] is False)
period=json.loads((HERE/"volume_period_receipt.json").read_text())
check("period obstruction retains one complex boundary-pair orbit",
      period["surviving_complex_boundary_orbits"]==1 and
      period["surviving_boundary"]==dict(B=["t","s**2"],E=["t**2-s**2","s**2*t"]))
check("Wright comparison retains actual source inclusion and omitted values",
      "inside C[X,Y]" in period["marked_ring"] and
      period["omitted_source_values"]==dict(X=-1,Y=2))
check("zero-period model remains unexcluded and no finite map is fabricated",
      period["assurance"]["remaining_zero_period_model_excluded"] is False and
      period["assurance"]["finite_map_constructed"] is False)
check("period argument does not prove actual compatible embedding or full quartic exclusion",
      period["assurance"]["actual_ruled_compatibility_proved"] is False and
      period["assurance"]["quartic_packets_excluded"] is False)
grading=json.loads((HERE/"wright_grading_receipt.json").read_text())
check("grading gate retains the mixed-weight and actual compatibility gaps",
      grading["assurance"]["unrestricted_mixed_weight_pairs_excluded"] is False and
      grading["assurance"]["actual_ruled_compatibility_proved"] is False)
check("trace tangent-field inference is not silently adopted",
      grading["assurance"]["actual_wright_trace_tangency_proved"] is False and
      grading["trace_countercontrol"]["branch_tangency"] is False)
check("countercontrol is a finite quartic control, not a fabricated Keller map",
      grading["trace_countercontrol"]["finite_degree"]==4 and
      grading["trace_countercontrol"]["is_wright_candidate"] is False and
      grading["trace_countercontrol"]["is_keller_counterexample"] is False)
check("invariant-target exclusion records its exact external input",
      "external" in grading["assurance"]["invariant_target_algebra_obstruction"])
p3all=json.loads((HERE/"pole_three_all_height_receipt.json").read_text())
check("intrinsic all-height pole-three exclusion gives minimum coordinate pole four",
      p3all["minimum_target_coordinate_pole"]==4 and p3all["first_open_reduced_pairs"]==[[6,8],[7,8]])
check("pole-three proof uses neither external degree bounds nor grading descent",
      p3all["assurance"]["external_2022_degree_theorem_used"] is False and
      p3all["assurance"]["plane_torus_linearization_used"] is False)
check("pole-three claim preserves remaining global gates",
      p3all["assurance"]["unrestricted_higher_minimum_poles_excluded"] is False and
      p3all["assurance"]["actual_ruled_compatibility_proved"] is False and
      p3all["assurance"]["quartic_packets_excluded"] is False)
p4companions=json.loads((HERE/"pole_four_companions_receipt.json").read_text())
check("pole-four odd exclusion preserves the even branch and exact four sextic partitions",
      p4companions["remaining_pole_four_companion_residue"]==dict(modulus=4,residue=2,minimum=10) and
      p4companions["remaining_sextic_root_partitions"]==[[6],[4,2],[3,3],[2,2,2]])
check("pole-four result does not inflate the coordinate bound or close height eight",
      p4companions["minimum_target_coordinate_pole"]==4 and
      p4companions["first_open_reduced_pairs"]==[[6,8],[7,8]] and
      p4companions["assurance"]["even_companion_layer_excluded"] is False and
      p4companions["assurance"]["all_pole_four_coordinates_excluded"] is False)
check("pole-four proof distinguishes the external low-height base from the intrinsic high odd step",
      p4companions["assurance"]["external_2022_degree_theorem_used_for_height_seven_base"] is True and
      p4companions["assurance"]["odd_companions_n_ge_9_use_external_degree_theorem"] is False)
check("pole-four claim preserves actual compatibility and full-packet gaps",
      p4companions["assurance"]["actual_ruled_compatibility_proved"] is False and
      p4companions["assurance"]["quartic_packets_excluded"] is False)
p4nonsquare=json.loads((HERE/"pole_four_nonsquare_receipt.json").read_text())
check("PF10 handoff retains exactly three cubic-square partitions before PF11",
      p4nonsquare["excluded_sextic_root_partitions"]==[[3,3]] and
      p4nonsquare["remaining_sextic_root_partitions"]==[[6],[4,2],[2,2,2]] and
      p4nonsquare["assurance"]["remaining_partitions_are_PF2RULED10_handoff_before_PF2RULED11"] is True)
check("small-root lemma strengthens but does not exclude the square-free cubic case",
      p4nonsquare["square_free_cubic_case_additional_form"]=="P=S^2+H, S in V2, H in V1" and
      p4nonsquare["assurance"]["remaining_cubic_square_cases_excluded"] is False)
check("formal-jet theorem does not raise minimum coordinate pole or close height eight",
      p4nonsquare["minimum_target_coordinate_pole"]==4 and
      p4nonsquare["first_open_reduced_pairs"]==[[6,8],[7,8]] and
      p4nonsquare["assurance"]["all_pole_four_coordinates_excluded"] is False)
check("finite coefficient controls are not misrepresented as the all-height proof",
      p4nonsquare["assurance"]["finite_z1_controls_are_not_the_all_height_proof"] is True and
      p4nonsquare["assurance"]["formal_verification"] is False and
      p4nonsquare["assurance"]["formal_centralizer_and_uniform_jet_coverage"]=="written proof")
check("nonsquare high-height step is intrinsic and global compatibility remains open",
      p4nonsquare["assurance"]["all_n_ge_10_nonsquare_step_uses_external_degree_theorem"] is False and
      p4nonsquare["assurance"]["actual_ruled_compatibility_proved"] is False and
      p4nonsquare["assurance"]["quartic_packets_excluded"] is False)
leadingperiod=json.loads((HERE/"boundary_leading_period_receipt.json").read_text())
check("leading-period classification covers all cubic and sextic contact multiplicities",
      leadingperiod["rational_primitive_classification"][0]["allowed_away_root_multiplicities"]==[0,1,2,3] and
      leadingperiod["rational_primitive_classification"][1]["allowed_away_root_multiplicities"]==[0,2,3,4,5,6])
check("leading period excludes squarefree cubics but preserves two pole-four types",
      leadingperiod["excluded_pole_four_sextic_root_partitions"]==[[2,2,2]] and
      leadingperiod["remaining_pole_four_sextic_root_partitions"]==[[6],[4,2]] and
      leadingperiod["minimum_target_coordinate_pole"]==4)
check("leading rational primitive is intrinsic without assuming a finite map or external degree bound",
      leadingperiod["assurance"]["root_period_proof_uses_finiteness"] is False and
      leadingperiod["assurance"]["root_period_proof_uses_external_degree_theorem"] is False)
check("source exactness is extended only through contact order3",
      leadingperiod["exact_source_degree_contact_orders"]==[0,1,2,3])
check("height-eight corollary retains its precise published exception and open frontier",
      leadingperiod["assurance"]["height_eight_numeric_corollary_uses_external_2022_degree_theorem"] is True and
      leadingperiod["assurance"]["published_72_108_exception_retained"] is True and
      leadingperiod["assurance"]["conditional_2026_degree125_upgrade_used"] is False and
      leadingperiod["first_open_reduced_pairs"]==[[6,8],[7,8]])
check("meromorphic positive controls are not mistaken for globally regular finite realizations",
      all(x["finite_map"] is False and x["globally_regular_pair"] is False
          for x in leadingperiod["meromorphic_leading_models"]))
check("Newton-corner comparison and full quartic completion are not silently inferred",
      leadingperiod["assurance"]["standard_Newton_corner_identification_proved"] is False and
      leadingperiod["assurance"]["height_eight_completely_excluded"] is False and
      leadingperiod["assurance"]["all_pole_four_coordinates_excluded"] is False and
      leadingperiod["assurance"]["actual_ruled_compatibility_proved"] is False and
      leadingperiod["assurance"]["quartic_packets_excluded"] is False)
pole58=json.loads((HERE/"pole58_receipt.json").read_text())
check("complete pole58 exclusion gives exact two-pair current frontier",
      pole58["excluded_pole_pairs"]==[[5,8]] and
      pole58["first_open_reduced_pairs"]==[[6,8],[7,8]] and
      pole58["minimum_target_coordinate_pole"]==4)
check("noncontact pole78 excluded with the exact contact exception retained",
      pole58["excluded_noncontact_pole_pairs"]==[[5,8],[7,8]] and
      pole58["remaining_pole78_requires"]==dict(contact_symbol=True,ordinary_degrees=[72,108]) and
      pole58["assurance"]["coalescent_Newton_faces_retained"] is True)
check("all six complete first-edge pole bounds exceed quartic degree",
      [x["Q_pole_lower_bound_per_extra_ramification"] for x in pole58["complete_coprime_height_eight_Newton_tables"]]==[12,17,10,15,20,25])
check("older checker current frontiers explicitly adopt the later theorem",
      all(json.loads((HERE/name).read_text()).get("frontier_adoption","").startswith("PF2RULED12")
          for name in ("pole4_exclusion_receipt.json","low_coordinate_receipt.json","pole35_receipt.json",
                       "source_coordinates_receipt.json","volume_period_receipt.json","wright_grading_receipt.json",
                       "pole_three_all_height_receipt.json","pole_four_companions_receipt.json",
                       "pole_four_nonsquare_receipt.json","boundary_leading_period_receipt.json")))
check("period table retains its necessary pre-finiteness coverage",
      "before the PF2RULED12" in leadingperiod["height_eight_frontier_table_scope"])
check("pole58 imports neither superseded gcd25 nor a global gcd16 theorem",
      pole58["assurance"]["global_gcd16_exclusion_imported"] is False and
      pole58["assurance"]["gcd25_from_superseded_2017_v1_imported"] is False and
      pole58["assurance"]["conditional_2026_degree125_upgrade_used"] is False)
check("new pair exclusions preserve all remaining mathematical gates",
      pole58["assurance"]["complete_height_eight_excluded"] is False and
      pole58["assurance"]["all_pole_four_coordinates_excluded"] is False and
      pole58["assurance"]["actual_ruled_compatibility_proved"] is False and
      pole58["assurance"]["quartic_packets_excluded"] is False)
cube68=json.loads((HERE/"pole68_cube_receipt.json").read_text())
check("noncontact pole68 cube reduction excludes slope4 and retains exact three profiles",
      cube68["excluded_perfect_cube_first_slopes"]==[4] and cube68["remaining_first_slopes"]==[5,6] and
      [p["high_discriminant_order"] for p in cube68["necessary_profiles"]]==[61,42,40] and
      [p["Q_poles"] for p in cube68["necessary_profiles"]]==[[3,1],[1,1,1,1],[1,1,1,1]])
check("cube profiles and local curves do not fabricate Q or finite realization",
      cube68["assurance"]["profiles_are_necessary_not_realizability_certificates"] is True and
      cube68["assurance"]["local_curve_controls_construct_Q"] is False and
      cube68["assurance"]["required_fibre_residues_proved_zero"] is False and
      cube68["assurance"]["any_profile_realizable_as_global_finite_map_proved"] is False and
      all(x["Q_constructed"] is False and x["is_global_finite_Wright_map"] is False
          for x in cube68["local_curve_controls"]))
check("historical cube reduction is intrinsic and retains its own pre-PF14 coverage",
      cube68["assurance"]["new_noncontact_cube_analysis_uses_external_degree_theorem"] is False and
      cube68["assurance"]["noncontact_pole68_excluded"] is False and
      cube68["first_open_reduced_pairs"]==[[6,8],[7,8]] and cube68["minimum_target_coordinate_pole"]==4)
check("cube-profile theorem preserves all global completion gaps",
      cube68["assurance"]["complete_height_eight_excluded"] is False and
      cube68["assurance"]["all_pole_four_coordinates_excluded"] is False and
      cube68["assurance"]["actual_ruled_compatibility_proved"] is False and
      cube68["assurance"]["quartic_packets_excluded"] is False)
cusp=json.loads((HERE/"cusp_leading_receipt.json").read_text())
check("historical profile gate is explicitly superseded",
      cube68["profile_gate_superseded_by"]=="PF2RULED14" and
      records["PF2RULED13"]["replaced_by"]==["PF2RULED14"] and
      records["PF2RULED14"]["supersedes"]==["PF2RULED13"])
check("global cusp obstruction covers all heights and noncoordinate target elements",
      cusp["all_height_one_root_noncontact_excluded"] is True and
      cusp["excludes_noncoordinate_target_elements_with_indicated_symbol"] is True and
      cusp["noncontact_pole68_excluded"] is True and
      cusp["supersedes_profile_realization_gate"]=="PF2RULED13")
check("historical PF14 height-eight handoff was contact72,108 only",
      cusp["first_open_reduced_pairs"]==[[6,8],[7,8]] and
      cusp["remaining_height_eight"]==dict(all_contact=True,ordinary_degrees=[72,108]))
check("contact pole-four degree bound is not an automatic Keller exclusion",
      cusp["pole_four_source_degree_upper_bound"]==48 and
      cusp["minimum_target_coordinate_pole"]==4 and
      cusp["assurance"]["arbitrary_source_component_degree48_excluded"] is False and
      cusp["assurance"]["all_pole_four_coordinates_excluded"] is False)
check("cusp algebra uses the actual marked target without an intermediate-plane inference",
      cusp["assurance"]["target_algebra_identity_uses_actual_marked_inclusion"] is True and
      cusp["assurance"]["auxiliary_lift_put_in_target_algebra"] is False and
      cusp["assurance"]["quartic_cusp_and_Euler_input"]=="PF2OMIT9 written theorem" and
      "PF2OMIT9" in records["PF2RULED14"]["dependencies"])
cert=json.loads((HERE/"cusp_leading_certificates.json").read_text())
check("exact incompatible-coefficient unit certificate is retained",
      cert["unit_certificate"]=="2*(2*a**2-1)-(4*a**2-3)=1" and
      cert["construct_finite_map"] is False)
check("global cusp proof preserves contact, higher gcds and full compatibility gaps",
      all(cusp["assurance"][key] is False for key in
          ["contact_height_eight_excluded","higher_gcds_classified",
           "conditional_2026_degree125_upgrade_used","actual_ruled_compatibility_proved",
           "quartic_packets_excluded","formal_verification","independent_full_reproof","external_review"]))
contact=json.loads((HERE/"contact_filtration_receipt.json").read_text())
contactcert=json.loads((HERE/"contact_filtration_certificates.json").read_text())
check("historical PF15 contact filtration retains exactly three necessary strata",
      contact["remaining_contact_multiplicities"]=={"6,8":[2],"7,8":[2,3]} and
      contact["first_open_reduced_pairs"]==[[6,8],[7,8]] and
      contact["ordinary_degrees"]==[72,108])
check("actual boundary restriction degrees and cusp multiplier slope are fixed",
      contact["boundary_restriction_degrees"]==[18,27] and
      contact["common_cusp_multiplier_degree_upper_bound"]==1 and
      contact["common_cusp_multiplier_slope"]=="20/27")
check("all finite source section spaces and rational basis transformations are certified",
      contact["assurance"]["complete_finite_section_spaces"] is True and
      contact["assurance"]["exact_rational_basis_change"] is True and
      [x["dimension"] for x in contactcert["filtered_summaries"]]==[63,77,112] and
      [x["top_form"] for x in contactcert["filtered_summaries"]]==["X^54*Y^18","X^54*Y^18","X^81*Y^27"])
check("linear lower-pole first-jet freedom is retained as a negative control",
      contact["pole7_first_jet_obstruction_from_unrestricted_lower_pole_lifts"] is False and
      [(x["rank"],x["jet_monomials"]) for x in contactcert["jet_ranks"]]==[(35,38),(38,38),(52,56)])
check("contact cusp unit and pole-six coefficient constraints are retained",
      contactcert["cusp_multiplier"]["delta"]=="-7*q0_26/(729*q0_27)" and
      contactcert["cusp_multiplier"]["delta_nonzero"] is True and
      contactcert["cusp_multiplier"]["pole6_symbol_order6"]=="-p0_18/3")
check("contact filtration does not silently import a Newton-corner identification or completion",
      all(contact["assurance"][key] is False for key in
          ["source_top_terms_assigned_to_published_Newton_corner","conditional_2026_degree125_upgrade_used",
           "complete_height_eight_excluded","all_pole_four_coordinates_excluded",
           "actual_ruled_compatibility_proved","quartic_packets_excluded","formal_verification",
           "independent_full_reproof","external_review"]))
height8=json.loads((HERE/"height_eight_receipt.json").read_text())
height8cert=json.loads((HERE/"height_eight_certificates.json").read_text())
check("complete height-eight exclusion closes every historical contact stratum",
      height8["complete_pole_height_eight_excluded"] is True and
      height8["all_three_PF2RULED15_contact_strata_excluded"] is True and
      height8["excluded_complete_pole_pairs"]==[[6,8],[7,8]] and
      records["PF2RULED15"]["replaced_by"]==["PF2RULED16"] and
      records["PF2RULED16"]["supersedes"]==["PF2RULED15"])
check("current frontier is height9 with four precise reduced pairs",
      height8["minimum_compatible_quartic_pole_height"]==9 and
      height8["minimum_target_coordinate_pole"]==4 and
      height8["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]])
check("complete chart-leading projection and lower-coordinate bound are retained",
      height8["leading_chart_companion"]=="c*v^7*(u+v)^20" and
      height8["lower_coordinate_chart_degree_upper_bound"]==20 and
      height8cert["dimensions"]=={"W_7_72":77,"W_8_108":112} and
      height8cert["universal_companion_high_total_identity"]=="Q_[>=27]=c*v^7*(u+v)^20-s4*u*v^7*(u+v)^20-s5*u*v^7*(u+v)^19")
check("all diagonal eigenvalues in the complete degree range are nonzero",
      set(height8cert["eigenvalues"])==set(map(str,range(1,21))) and
      all(len(values)==int(d)+1 and all(value!=0 for value in values)
          for d,values in height8cert["eigenvalues"].items()) and
      height8cert["resonant_degree27_control"] is True)
check("intrinsic leading obstruction separates core algebra from the imported quartic base",
      all(height8["assurance"][key] is False for key in
          ["core_chart_obstruction_uses_finiteness","core_chart_obstruction_uses_cusp_Euler_input",
           "core_chart_obstruction_uses_external_degree_theorem"]) and
      height8["assurance"]["quartic_height_eight_corollary_uses_retained_2022_degree_alternative"] is True and
      height8["assurance"]["uses_shared_polynomial_utilities_from_PF2RULED15"] is True)
check("height-eight theorem preserves the original higher-pole and compatibility obligations",
      all(height8["assurance"][key] is False for key in
          ["published_72_108_exception_excluded_in_all_pole_heights",
           "source_top_terms_assigned_to_published_Newton_corner","conditional_2026_degree125_upgrade_used",
           "height_nine_excluded","all_pole_four_coordinates_excluded","actual_ruled_compatibility_proved",
           "quartic_packets_excluded","formal_verification","independent_full_reproof","external_review"]))
height9=json.loads((HERE/"height_nine_frontier_receipt.json").read_text())
height9cert=json.loads((HERE/"height_nine_frontier_certificates.json").read_text())
check("height-nine refinement preserves all four open pairs and minimum height",
      height9["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      height9["minimum_compatible_quartic_pole_height"]==9)
rows9=height9["gcd_one_height_nine_necessary_rows"]
check("ten exact necessary gcd-one rows retained without a finite-map inference",
      len(rows9)==10 and len([x for x in rows9 if x["ordinary_degrees"][1]==128])==8 and
      len([x for x in rows9 if x["ordinary_degrees"]==[72,108]])==2 and
      height9["assurance"]["all_retained_linear_rows_realizable_as_finite_maps"] is False)
check("pole-five row has exact ordinary degree and contact restrictions",
      [(x["ordinary_degrees"],x["primitive_cubic_contact_multiplicities"]) for x in rows9 if x["poles"]==[5,9]]==[([48,128],[2,3])])
check("normal exceptional resonance survives while its reverse is excluded",
      height9["exceptional_normal_orientation_retained"] is True and
      height9["exceptional_reverse_orientation_excluded"] is True and
      height9cert["exceptional_resonant_chart_forms"]==["c1*v^6*(u+v)^12","c2*v^9*(u+v)^18"] and
      height9cert["multiplicity_gcds"]==[1,1,1,1,1,1,1,9])
check("distinct cusp multiplier slopes and source shear limitations are explicit",
      height9["normal_exception_cusp_multiplier_slope"]=="2/3" and
      height9["degree128_cusp_multiplier_slope"]=="3/4" and
      height9["assurance"]["ordinary_shear_preserves_pole_height"] is False and
      height9["assurance"]["ordinary_shear_preserves_Keller_and_field_degree"] is True and
      height9cert["excluded_divisible_source_partner_exponents"]==[1,2,4])
check("complete contact kernels retain all lower-pole terms",
      height9["assurance"]["complete_finite_section_spaces"] is True and
      height9["assurance"]["all_lower_pole_terms_retained"] is True and
      {k:v["dimension"] for k,v in height9cert["spaces"].items()}=={"5":51,"7":92,"8":117,"9":145} and
      all(x["constraint_rank"]+x["kernel_dimension"]==x["input_dimension"] for x in height9cert["kernel_certificates"]))
check("gcd3 and original compatibility goal remain unresolved after the height-nine refinement",
      all(height9["assurance"][key] is False for key in
          ["gcd3_pair_6_9_classified","all_gcd_one_height_nine_pairs_excluded",
           "unrestricted_72_108_exception_excluded","conditional_2026_degree125_upgrade_used",
           "actual_ruled_compatibility_proved","quartic_packets_excluded","formal_verification",
           "independent_full_reproof","external_review"]))
cubic=json.loads((HERE/"cubic_cover_period_receipt.json").read_text())
cubiccert=json.loads((HERE/"cubic_cover_period_certificates.json").read_text())
check("all-height gcd3 period has the correct cover exactness and per-input polynomial test",
      cubic["all_height_gcd3_leading_period_criterion_proved"] is True and
      cubic["period_test_necessary_and_sufficient_for_cover_exactness"] is True and
      cubiccert["period_polynomial_identity"]=="3*s*A_prime+(3*s_prime-2*s*h_prime/h)*A=1")
check("low-degree root positions are classified without claiming every degree9 locus",
      cubic["degree6_multi_root_positions_classified"] is True and
      cubic["all_period_root_positions_classified"] is False and
      [x["criterion"] for x in cubiccert["degree6_three_root_classification"]]==["2*t**2 + t + 2","2*t + 1","t**2 - t + 1"])
check("cover primitive degree and extremal root partition controls are retained",
      cubiccert["primitive_pole_degree"]=="2*D-3*r+number_of_simple_roots" and
      cubiccert["primitive_infinity_contact_order"]=="2*D/3-1" and
      [x["partition"] for x in cubiccert["root_partition_census"] if x["degree"]==9 and
       x["root_count"]==6 and not x["cube"] and x["passes_primitive_pole_budget"]]==[[4,1,1,1,1,1]])
check("positive-genus and noncontact period controls are not fabricated quartic maps",
      cubic["noncontact_period_controls_retained"] is True and
      cubic["assurance"]["contact_forced_for_every_gcd3_pair"] is False and
      cubic["assurance"]["period_solution_constructs_actual_PQ"] is False and
      cubiccert["finite_map_constructed"] is False and
      any(x["h"]=="(w**3-1)**2" for x in cubiccert["positive_period_controls"]))
frontier69=cubic["pole69_necessary_frontier"]
check("historical PF18 numerical groups preceded the generic-polar filter",
      cubic["complete_source_degree_cases_for_pole69"] is True and
      cubic["pole69_contact_multiplicity2_excluded"] is True and
      [x["ordinary_degrees"] for x in frontier69]==[[96,144],[88,132],[80,128],[48,128],[72,108]] and
      [x["contact_multiplicities"] for x in frontier69]==[[0],[1],list(range(3,10)),[7,8,9],[3]] and
      [x["cusp_multiplier_slope"] for x in frontier69]==["3/4","8/11","3/4","3/4","2/3"])
check("pole69 source kernels retain exact minima and all lower-pole terms",
      {k:v["dimension"] for k,v in cubiccert["source_spaces"].items()}=={"6":70,"9":145} and
      [x["minimum_symbol_order"] for x in cubiccert["constrained_partner_spaces"]]==[13,5,6] and
      all(x["rank"]+x["kernel_dimension"]==x["input_dimension"] for x in cubiccert["kernel_certificates"]))
check("auxiliary cover and intrinsic period are separate from the actual finite target field",
      cubic["assurance"]["auxiliary_cover_is_target_intermediate_field"] is False and
      cubic["assurance"]["period_criterion_uses_finiteness"] is False and
      cubic["assurance"]["period_criterion_uses_external_degree_theorem"] is False and
      cubic["assurance"]["source_degree_frontier_uses_retained_2022_alternative"] is True)
check("cubic-cover refinement preserves every open pair and the original goal limits",
      cubic["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      cubic["minimum_compatible_quartic_pole_height"]==9 and
      all(cubic["assurance"][key] is False for key in
          ["all_pole69_pairs_excluded","all_degree9_period_positions_classified",
           "actual_ruled_compatibility_proved","quartic_packets_excluded","conditional_2026_degree125_upgrade_used",
           "formal_verification","independent_full_reproof","external_review"]))
polar=json.loads((HERE/"polar_cusp_sector_receipt.json").read_text())
polarcert=json.loads((HERE/"polar_cusp_sector_certificates.json").read_text())
check("all-height generic-polar local budget includes every boundary contact correction",
      polar["all_height_generic_polar_budget_proved"] is True and
      polar["all_height_common_root_multiplicity_rule_proved"] is True and
      polar["root_support_bound_for_finite_quartic"]==4 and
      polarcert["generic_polar_formula"]=="sum nu=I(F,A_w)-(m-2)*l-I(F,E)<=2*l-m-I(F,E)" and
      polarcert["common_symbol_root_rule"]=="2*k>d at every root of common section in L^d")
check("old five-group source stage is explicitly superseded while its period theorem remains",
      cubic["source_frontier_superseded_by"]=="PF2RULED19" and
      cubic["source_frontier_scope"].startswith("PF2RULED18 before generic-polar"))
check("current pole69 source frontier has three exact groups and excludes72,108",
      polar["current_pole69_frontier"]==[
        {"contacts":[0],"ordinary_degrees":[96,144]},
        {"contacts":[3,4,5,6,7,9],"ordinary_degrees":[80,128]},
        {"contacts":[7,9],"ordinary_degrees":[48,128]}] and
      polar["pole69_72_108_excluded"] is True and polarcert["excluded_contacts"]==[1,2,8])
check("historical PF19 noncontact handoff had five necessary partitions",
      polar["current_noncontact_nonic_partitions"]==[[7,2],[5,4],[5,2,2],[4,3,2],[3,2,2,2]])
check("complete high cusp-jet block checks both full pole spaces",
      polar["assurance"]["full_section_high_jet_coverage"] is True and
      polarcert["complete_section_jet_counts"]=={"6":70,"9":145} and
      polarcert["noncontact_marked_cusp"]["required_nonzero"]=="a")
sectors=polarcert["two_root_marked_sectors"]
check("two-root leading sectors retain sixteen exact marked points without actual map claims",
      polar["marked_two_root_sector_exact_point_count"]==16 and
      [(x["multiplicities"],x["distinct_marked_coefficient_points"]) for x in sectors]==[([2,7],8),([4,5],8)] and
      all(x["complete_leading_period_and_high_jet_sector"] is True and x["global_finite_map"] is False and
          set(x["units"])=={"a","root_separation","x","y"} for x in sectors) and
      polar["assurance"]["marked_two_root_classification_is_only_leading_tests"] is True and
      polar["assurance"]["sixteen_points_construct_actual_maps"] is False)
profile=polarcert["elliptic_contact3_necessary_profile"]
check("historical PF19 elliptic profile is retained with its own pre-gap coverage",
      profile["Q_poles"]==[1,1,1,1] and profile["away_polar_intersections"]==[17,17,17] and
      profile["contact_polar_intersection"]==29 and profile["global_PQ_constructed"] is False and
      profile["profile_excluded"] is False and polar["assurance"]["elliptic_contact3_profile_excluded"] is False)
check("polar and cusp sector filters preserve all whole-pair and actual-compatibility gaps",
      polar["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      polar["minimum_compatible_quartic_pole_height"]==9 and
      all(polar["assurance"][key] is False for key in
          ["period_solution_implies_finite_map","all_noncontact_pole69_excluded","all_pole69_excluded",
           "full_degree9_period_position_classification","all_height_actual_72_108_exception_excluded",
           "actual_ruled_compatibility_proved","quartic_packets_excluded","conditional_2026_degree125_upgrade_used",
           "formal_verification","independent_full_reproof","external_review"]))
quarticgap=json.loads((HERE/"quartic_polar_gap_receipt.json").read_text())
quarticgapcert=json.loads((HERE/"quartic_polar_gap_certificates.json").read_text())
check("generic quartic discriminant gap gives the full cone and exact local pole costs",
      quarticgap["generic_quartic_polar17_excluded"] is True and
      quarticgap["forced_generic_polar_order"]==18 and
      quarticgap["forced_quartic_orders"]=={"p":3,"q":5,"r":6} and
      quarticgap["finite_local_branch_count"]==2 and quarticgap["finite_local_Q_poles"]==[1,1] and
      len(quarticgapcert["cone_monomial_bounds"])==13 and
      all(x["minimum_z_order"]>=18 for x in quarticgapcert["cone_monomial_bounds"]))
check("genericity and all-order normal-form scope are preserved",
      quarticgap["assurance"]["higher_terms_truncated_in_proof"] is False and
      quarticgap["assurance"]["normal_form_coordinate_change_depends_on_c"] is False and
      quarticgap["assurance"]["claims_no_polar17_for_all_special_fibres"] is False and
      any(x.get("special_polar_order")==17 and x.get("generic_polar_order")==14 and
          x["global_finite_map"] is False for x in quarticgapcert["controls"]))
check("elliptic realization gate and old partition handoff are explicitly superseded",
      polar["elliptic_profile_superseded_by"]=="PF2RULED20" and
      polar["noncontact_partition_frontier_superseded_by"]=="PF2RULED20" and
      profile["profile_stage"]=="PF2RULED19 necessary data before PF2RULED20" and
      quarticgap["elliptic_contact3_pole69_excluded"] is True and
      quarticgapcert["minimum_elliptic_contact3_Q_pole_degree"]==7)
check("historical PF20 stage closed two sectors and handed off three necessary types",
      quarticgap["excluded_noncontact_nonic_partitions"]==[[5,2,2],[3,2,2,2]] and
      quarticgap["remaining_noncontact_nonic_partitions"]==[[7,2],[5,4],[4,3,2]] and
      quarticgap["current_contact3_away_types"]==[[6],[4,2]])
nextnoncontact=quarticgapcert["noncontact432_necessary_profile"]
nextcontact=quarticgapcert["contact42_necessary_profile"]
check("historical PF20 higher-polar profiles retain their own stage coverage",
      nextnoncontact["polar_intersections_by_root"]==[33,25,18] and
      nextnoncontact["missing_points_by_root"]==[1,1,2] and
      nextcontact["polar_intersections_away"]==[33,18] and
      nextcontact["polar_intersection_at_contact"]==29 and
      all(x["global_PQ_constructed"] is False and x["profile_excluded"] is False
          for x in (nextnoncontact,nextcontact)))
check("quartic-gap theorem does not close sixteen leading candidates or a whole height-nine pair",
      quarticgap["marked_two_root_leading_candidates_retained"]==16 and
      quarticgap["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      quarticgap["minimum_compatible_quartic_pole_height"]==9 and
      all(quarticgap["assurance"][key] is False for key in
          ["local_generic_cone_constructs_global_finite_map","noncontact432_profile_excluded",
           "contact42_profile_excluded","all_pole69_excluded","actual_ruled_compatibility_proved",
           "quartic_packets_excluded","conditional_2026_degree125_upgrade_used","formal_verification",
           "independent_full_reproof","external_review"]))
unibranch=json.loads((HERE/"unibranch_polar33_receipt.json").read_text())
unibranchcert=json.loads((HERE/"unibranch_polar33_certificates.json").read_text())
check("PF20 higher-polar handoff is explicitly superseded by PF21",
      quarticgap["higher_polar_profiles_superseded_by"]=="PF2RULED21" and
      quarticgap["noncontact_partition_frontier_superseded_by"]=="PF2RULED21" and
      quarticgapcert["higher_polar_profiles_superseded_by"]=="PF2RULED21" and
      quarticgap["profile_scope"]=="PF2RULED20 necessary data before PF2RULED21")
check("generic single-branch degree8 polar33 is excluded by exhaustive chains and flow",
      unibranch["generic_unibranch_degree8_polar33_excluded"] is True and
      unibranchcert["exhaustive_degree8_polar33_chains"]==[[[2,2],[21,1]]] and
      unibranchcert["forced_original_w_order"]==2 and
      unibranchcert["polar33_required_N"]==27 and unibranchcert["polar33_flow_order"]==24 and
      unibranchcert["polar33_bound"]==30)
check("all-height obstruction requires the whole generic germ and original coordinates",
      unibranch["all_height_whole_unibranch_w_order2_excluded_for_poles_at_least6"] is True and
      unibranchcert["all_height_exclusion_m_minimum"]==6 and
      unibranchcert["critical_value_identity"]=="partial_c K=-phi_c^m" and
      unibranch["assurance"]["whole_generic_germ_unibranch_required"] is True and
      unibranch["assurance"]["uses_parameter_dependent_ambient_change_before_flow"] is False)
check("both old polar33 profiles close without assigning an exact multiplicity4 branch count",
      unibranch["noncontact432_excluded"] is True and unibranch["contact3_away42_excluded"] is True and
      unibranchcert["common_multiplicity4_Q_pole_lower_bound"]==2 and
      unibranchcert["lower_bound_sharpness_proved"] is False and
      unibranch["assurance"]["claims_exact_two_branches_at_common_multiplicity4"] is False)
check("remaining noncontact types are72/54 and contact3 has only cube6",
      unibranch["remaining_noncontact_nonic_partitions"]==[[7,2],[5,4]] and
      unibranch["current_contact3_away_types"]==[[6]] and
      unibranch["current_pole69_frontier"]==polar["current_pole69_frontier"])
check("every retained pole69 common nonic has at most two support points",
      unibranch["all_current_pole69_common_nonic_support_at_most2"] is True and
      [x["contact"] for x in unibranchcert["current_contact_root_forms"]]==[3,4,5,6,7,9] and
      all(len(x["away_multiplicities"])<=1 for x in unibranchcert["current_contact_root_forms"]))
check("current auxiliary components are rational and period gate already passes",
      unibranch["all_current_auxiliary_components_rational"] is True and
      unibranch["all_current_root_forms_pass_leading_period"] is True and
      [x["k"] for x in unibranchcert["auxiliary_rational_parametrizations"]]==[2,4,5,7] and
      all(3*x["u"]+x["k"]*x["v"]==1 for x in unibranchcert["auxiliary_rational_parametrizations"]) and
      unibranch["assurance"]["auxiliary_cover_is_target_intermediate_field"] is False and
      unibranch["assurance"]["auxiliary_Mobius_change_is_marked_source_automorphism"] is False)
control33=unibranchcert["special33_control"]
check("special polar33 control preserves the generic and whole-germ distinctions",
      control33["special_polar_order"]==33 and control33["generic_total_polar"]==30 and
      control33["generic_branch_polar_orders"]==[15,15] and control33["finite_map"] is False and
      unibranch["assurance"]["claims_no_selected_fibre_polar33"] is False and
      unibranch["assurance"]["claims_all_multibranch_generic_polar33_impossible"] is False)
check("sixteen points stay leading candidates and the full original goal stays open",
      unibranch["marked_two_root_leading_candidates_retained"]==16 and
      unibranch["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      unibranch["minimum_compatible_quartic_pole_height"]==9 and
      unibranchcert["actual_PQ_constructed"] is False and
      all(unibranch["assurance"][key] is False for key in
          ["sixteen_points_construct_actual_maps","all_degree9_period_loci_classified",
           "all_pole69_excluded","actual_ruled_compatibility_proved","quartic_packets_excluded",
           "conditional_2026_degree125_upgrade_used","formal_verification","independent_full_reproof",
           "external_review","research_goal_completed"]))
fullcusp=json.loads((HERE/"full_cusp_lift_receipt.json").read_text())
fullcuspcert=json.loads((HERE/"full_cusp_lift_certificates.json").read_text())
check("all sixteen marked sectors lift through the complete first cusp jet",
      fullcusp["all_sixteen_full_first_jet_lifts_exist"] is True and
      fullcusp["complete_noncontact_first_cusp_jet_locus_has_sixteen_marked_points"] is True and
      fullcusp["marked_two_root_leading_candidates_retained"]==16 and
      fullcuspcert["marked_candidate_count"]==16 and
      len(fullcuspcert["sectors"])==2 and
      all(x["distinct_marked_points"]==8 and x["full_first_jet_lifts_exist_at_every_point"] is True
          for x in fullcuspcert["sectors"]))
check("triangular right inverse covers every required pole level and complex embedding",
      fullcusp["determinant_unit_certificates"]==16 and
      [x["n"] for x in fullcuspcert["triangular_matrices"]]==list(range(1,9)) and
      all([x["n"] for x in sec["determinant_units"]]==list(range(1,9))
          for sec in fullcuspcert["sectors"]) and
      fullcuspcert["full_section_coverage"]=={"3":22,"5":51,"6":70,"8":117,"9":145})
check("full cusp-lift dimensions include translations and are not finite-map dimension claims",
      fullcusp["full_cusp_lift_affine_dimensions"]=={"P":31,"Q":85,"pair":116} and
      fullcusp["assurance"]["dimensions_quotient_target_translations"] is False and
      fullcusp["assurance"]["dimension_count_implies_realizability"] is False)
check("clean singleton ordinary cusp and boundary birationality are available on an open locus",
      fullcusp["clean_singleton_cusp_and_boundary_birational_open_nonempty"] is True and
      [x["pole_upper_bound"] for x in fullcuspcert["cusp_modifications"]]==[4,5] and
      all(x["clean_singleton_ordinary_cusp_open_nonempty"] is True and
          x["boundary_birational_open_nonempty"] is True for x in fullcuspcert["sectors"]))
check("next actual normal Jacobian gate is explicit with bounded coefficient count",
      fullcusp["next_Jacobian_parameter_count"]==116 and
      fullcusp["next_Jacobian_quotient_degree_upper_bound"]==52 and
      fullcusp["next_Jacobian_quadratic_equation_count_upper_bound"]==53 and
      fullcuspcert["next_Jacobian_u1"]=="2*h*(P2*q-p*Q2)+lambda^2*(p*q_prime-p_prime*q)" and
      fullcuspcert["next_Jacobian_equations_solved"] is False)
check("actual-ring full cusp control rejects a source Keller inference",
      fullcuspcert["countercontrol"]["full_first_jet_cusp"] is True and
      fullcuspcert["countercontrol"]["Keller"] is False and
      fullcuspcert["countercontrol"]["finite_surface_map"] is False and
      fullcusp["assurance"]["first_jet_lifts_prove_source_Keller_pairs"] is False)
check("complete first-jet result preserves all higher Jacobian, finiteness and compatibility gaps",
      fullcusp["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      fullcusp["minimum_compatible_quartic_pole_height"]==9 and
      all(fullcusp["assurance"][key] is False for key in
          ["next_Jacobian_equations_solved","all_Jacobian_coefficients_solved",
           "field_degree_four_realization_constructed","global_surface_finiteness_proved",
           "all_noncontact_pole69_excluded","all_pole69_excluded","contact_cases_classified_by_this_lift",
           "actual_ruled_compatibility_proved","quartic_packets_excluded","conditional_2026_degree125_upgrade_used",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
normalfreedom=fullcuspcert["normal_freedom"]
check("normal freedoms use complete exact first-two-jet kernels and full second-coefficient images",
      [(x["pole"],x["first_jet_kernel_dimension"],x["second_jet_image_degree"],
        x["second_jet_image_dimension"],x["invisible_through_second_coefficient"]) for x in normalfreedom]==
      [(5,14,10,11,3),(8,56,22,23,33)] and
      all(len(x["nullspace_columns"])==x["first_jet_kernel_dimension"] for x in normalfreedom))
check("next normal gate accounts for all116 parameters before elimination",
      fullcusp["boundary_parameter_count"]==46 and fullcusp["second_coefficient_control_count"]==34 and
      fullcusp["invisible_through_second_coefficient"]==36 and 46+34+36==116)
check("eight immutable highest coefficients and eleven membership conditions are the reduced next gate",
      fullcusp["next_quotient_immutable_top_degrees"]==list(range(45,53)) and
      fullcuspcert["truncated_multiplication_matrix_shape"]==[45,34] and
      fullcusp["clean_open_truncated_multiplication_rank"]==34 and
      fullcusp["low_coefficient_consistency_count"]==11 and
      fullcusp["local_rational_consistency_count_upper_bound"]==19 and
      fullcusp["assurance"]["nineteen_reduced_conditions_claimed_quadratic"] is False and
      fullcusp["assurance"]["all_eight_high_coefficient_equations_solved"] is False)
high=json.loads((HERE/"cusp_high_jacobian_receipt.json").read_text())
highcert=json.loads((HERE/"cusp_high_jacobian_certificates.json").read_text())
check("highest next-Jacobian identities cover the full116-dimensional lift family and all16 points",
      high["marked_points"]==16 and high["full_affine_lift_dimension"]==116 and
      high["assurance"]["full_affine_space_covered"] and high["assurance"]["all_complex_embeddings_covered"] and
      high["automatic_quotient_degrees"]==list(range(48,53)))
check("three remaining high equations have exact all-point rank witnesses and solutions",
      high["active_high_quotient_degrees"]==[45,46,47] and high["fixed_Q_high_rank"]==3 and
      high["fixed_Q_affine_solution_dimension"]==28 and high["all_eight_high_equations_solved_at_all_marked_points"] and
      all(len(x["rank_minor_kernel_columns"])==3 and x["solved_P_lower_pole_correction"] for x in highcert["sectors"]))
check("high-block chart leaves43 boundary parameters and eleven low membership conditions",
      high["rank_three_chart_dimension"]==113 and high["rank_three_chart_boundary_parameters"]==43 and
      high["remaining_low_membership_conditions"]==11 and
      high["clean_singleton_cusp_open_nonempty_after_high_block"])
check("selected high-block witnesses fail low membership without excluding other boundary data",
      all(x["witness_low_matrix_rank"]==34 and x["witness_low_augmented_rank"]==35 for x in highcert["sectors"]) and
      high["assurance"]["no_candidate_removed"])
check("high-block result preserves all remaining global obligations",
      high["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(high["assurance"][key] is False for key in
          ["lower_membership_equations_solved","full_next_Jacobian_solved","all_higher_Jacobian_coefficients_solved",
           "field_degree_four_realization_constructed","global_finiteness_proved","contact_cases_excluded",
           "actual_ruled_compatibility_proved","quartic_packets_excluded","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]))
nextnormal=json.loads((HERE/"full_next_normal_receipt.json").read_text())
nextcert=json.loads((HERE/"full_next_normal_certificates.json").read_text())
check("all sixteen candidates lift through the full next-normal equation",
      nextnormal["marked_points"]==16 and nextnormal["full_next_normal_solutions_exist_at_all_marked_points"] and
      nextnormal["generic_fixed_P_next_normal_rank"]==49 and
      nextnormal["assurance"]["all_eleven_previous_low_conditions_have_simultaneous_lifts"])
check("residue units and irreducibility support the stated written characteristic-zero implication",
      nextnormal["residue_characteristic"]==1009 and
      nextnormal["assurance"]["characteristic_zero_existence_uses_written_unit_lifting"] and
      nextnormal["assurance"]["modular_rank_used_as_unit_certificate_not_search_evidence"] and
      nextnormal["assurance"]["all_complex_embeddings_covered"] and
      [x["next_normal_determinant_modp"] for x in nextcert["sectors"]]==[796,285] and
      [x["singleton_cusp_resultant_modp"] for x in nextcert["sectors"]]==[38,240])
check("full next-normal solution charts retain the clean singleton cusp and hidden dimensions",
      nextnormal["fixed_P_normalized_Q_affine_dimension"]==36 and
      nextnormal["generic_unnormalized_solution_chart_dimension"]==68 and
      nextnormal["generic_boundary_parameter_count"]==32 and
      nextnormal["clean_singleton_cusp_open_nonempty"])
check("following normal controls use both complete rational third-jet images",
      [(x["pole"],x["third_image_rank"],x["invisible_through_third_coefficient"]) for x in nextcert["following_normal_controls"]]==[(5,3,0),(8,18,15)] and
      nextnormal["following_normal_control_rank"]==21 and nextnormal["following_normal_correction_degree"]==40 and
      nextnormal["following_normal_immutable_high_degrees"]==list(range(41,44)) and
      nextnormal["following_normal_cusp_value_condition"] and nextnormal["following_normal_low_membership_conditions"]==19)
check("intrinsic Jacobian pole bound sharpens the next test and supplies a finite full-identity gate",
      nextnormal["source_Jacobian_pole_bound"]==12 and nextnormal["following_normal_degree_bound"]==43 and
      nextnormal["following_normal_consistency_count_upper_bound"]==23 and
      nextnormal["source_Jacobian_full_identity_follows_from_chart_orders"]==list(range(2,10)) and
      nextnormal["first_nine_V12_jets_rank"]==247)
check("following-normal failures of selected lifts do not close any full marked candidate or global gate",
      all(x["next_Jacobian_u2_modp"]!="0" for x in nextcert["sectors"]) and
      nextnormal["assurance"]["no_marked_candidate_removed"] and
      nextnormal["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(nextnormal["assurance"][key] is False for key in
          ["rank_drop_strata_excluded","following_normal_equations_solved","all_higher_Jacobian_coefficients_solved",
           "field_degree_four_realization_constructed","global_surface_finiteness_proved",
           "actual_ruled_compatibility_proved","contact_cases_classified","quartic_packets_excluded",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
second=json.loads((HERE/"second_normal_obstruction_receipt.json").read_text())
secondcert=json.loads((HERE/"second_normal_obstruction_certificates.json").read_text())
check("generic first-coordinate exclusion permits every Jacobian constant and covers all16 candidates",
      second["marked_points"]==16 and second["generic_P_excluded_even_with_variable_kappa"] and
      second["assurance"]["all_values_of_kappa_allowed_in_obstruction"] and second["assurance"]["all_complex_embeddings_covered"])
check("global determinant uses universal rank bounds and exact all-point unit minors",
      second["universal_Q_rank_upper_bound"]==68 and second["universal_variable_kappa_rank_upper_bound"]==69 and
      second["necessary_global_P_polynomial_degree_upper_bound"]==69 and
      [x["generic_chart_minor_determinant_modp"] for x in secondcert["sectors"]]==[280,326] and
      [x["obstruction_determinant_modp"] for x in secondcert["sectors"]]==[145,885])
check("global necessary polynomial retains rank-drop solutions while local elimination is complete on its stated chart",
      second["all_rank_drop_solutions_still_satisfy_global_obstruction"] and
      second["assurance"]["global_obstruction_includes_rank_drop_strata"] and
      second["generic_fixed_P_compatibility_count"]==24 and second["generic_Q_kappa_fiber_if_nonempty_dimension"]==17 and
      all(len(x["fixed_P_compatibility_rows"])==24 and len(x["obstruction_rows"])==70 and len(x["obstruction_columns"])==70 for x in secondcert["sectors"]))
check("reversible target normalization leaves30 P parameters and15 hidden companion directions",
      second["normalized_P_parameter_count"]==30 and second["normalized_Q_fiber_if_nonempty_dimension"]==15)
check("zero-kappa power controls do not establish any actual nonzero-kappa realization",
      secondcert["relaxed_kappa0_determinantal_locus_nonempty"] and
      second["assurance"]["commuting_power_control_only_has_kappa0"] and
      all(x["relaxed_kappa0_solution_constructed"] and not x["nonzero_kappa_solution_constructed"] for x in secondcert["sectors"]))
check("second-normal determinant leaves the nonzero locus and full research obligations open",
      second["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(second["assurance"][key] is False for key in
          ["local24_conditions_asserted_independent","nonzero_kappa_solution_locus_nonempty_proved",
           "second_normal_equations_solved_on_zero_locus","marked_candidate_excluded",
           "full_source_Jacobian_identity_solved","field_degree_four_realization_constructed",
           "global_surface_finiteness_proved","actual_ruled_compatibility_proved","contact_cases_classified",
           "quartic_packets_excluded","formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
gap=json.loads((HERE/"square_boundary_gap_receipt.json").read_text())
gapcert=json.loads((HERE/"square_boundary_gap_certificates.json").read_text())
gapcheckpoint=json.loads(gzip.decompress((HERE/"square_boundary_degree8_checkpoint.json.gz").read_bytes()))
check("rooted boundary degree-seven exclusion covers all sixteen candidates and ten exact field cases",
      gap["marked_points"]==16 and gap["boundary_difference_degree_excluded_through"]==7 and
      gap["nonzero_first_normal_Jacobian_requires_boundary_difference_degree_at_least"]==8 and
      gap["exact_number_field_unit_certificates"]==10 and
      all(x["full_S_affine_dimension"]==4 and [c["degree"] for c in x["exact_unit_cases"]]==[3,4,5,6,7] for x in gapcert["sectors"]))
check("actual-root condition and all normal freedoms remain explicit",
      gap["full_actual_S_pole3_cusp_family_covered"] and
      gap["assurance"]["condition_on_boundary_root_lift_not_automatic_for_all_P"] and
      gap["assurance"]["complete_S_family_and_all_normal_freedoms_retained"] and
      gap["assurance"]["exact_Bezout_identities_replayed"])
check("degree-eight checkpoint is retained without a completion or live-process claim",
      gap["full_degree8_ideal_solved"] is False and gap["degree8_high_remainder_dimension"]==[5,5] and
      all(x["full_saturated_ideal_solved"] is False for x in gapcheckpoint["sectors"]) and
      gapcheckpoint["earlier_exploratory_attempt"]["process_terminal_verified"])
check("boundary-gap theorem preserves the general determinant, higher and finite-map obligations",
      gap["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(gap["assurance"][key] is False for key in
          ["arbitrary_determinantal_zero_locus_excluded","degree8_nonzero_kappa_locus_solved",
           "general_nonzero_kappa_second_normal_locus_solved","whole_marked_candidate_excluded",
           "full_source_Jacobian_pair_constructed","field_degree_four_realization_constructed",
           "global_surface_finiteness_proved","actual_ruled_compatibility_proved","contact_cases_classified",
           "quartic_packets_excluded","formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
affine8=json.loads((HERE/"degree8_affine_reduction_receipt.json").read_text())
affine8cert=json.loads(gzip.decompress((HERE/"degree8_affine_reduction_certificates.json.gz").read_bytes()))
check("degree8 high ideal is globally polynomial affine5 at all sixteen marked points",
      affine8["marked_points"]==16 and affine8["high_locus_global_polynomial_isomorphism"]=="A5" and
      affine8["assurance"]["global_parameterization_has_no_variable_denominator"] and
      affine8["assurance"]["high_ideal_isomorphism_replayed"] and affine8["assurance"]["all_complex_embeddings_covered"])
check("exact old divisibility ideal becomes six linear root-value equations over four parameters",
      affine8["remaining_core_parameter_count"]==4 and affine8["full_linearized_divisibility_equation_count"]==6 and
      affine8["assurance"]["old_divisibility_ideal_preserved"] and affine8["assurance"]["nonzero_kappa_condition_retained"] and
      all(x["E_arithmetic_circuit"]=="remainder_v(M*(W/(v-s)), L)" for x in affine8cert["sectors"]))
check("all coefficient-pivot cases including the full zero stratum are retained",
      affine8["nonzero_D_chart_count"]==5 and affine8["zero_D_stratum_retained"] and
      affine8["assurance"]["all_zero_pivot_strata_retained"] and
      all(len(x["nonzero_D_charts"])==5 and x["zero_D_stratum"]["excluded"] is False for x in affine8cert["sectors"]))
check("affine reduction does not close any actual nonzero-kappa or finite-map obligation",
      affine8["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(affine8["assurance"][key] is False for key in
          ["zero_D_stratum_excluded","full_degree8_nonzero_kappa_locus_solved","general_determinant_locus_solved",
           "actual_global_root_hypothesis_removed","whole_marked_candidate_excluded","source_Keller_pair_constructed",
           "field_degree_four_realization_constructed","global_surface_finiteness_proved","actual_ruled_compatibility_proved",
           "contact_cases_classified","quartic_packets_excluded","formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
unitrow=json.loads((HERE/"degree8_unit_row_receipt.json").read_text())
recovery=json.loads((HERE/"degree8_global_recovery.json").read_text())
check("exact degree3 unit rows close the complete all-zero-coefficient stratum at all16 points",
      unitrow["marked_points"]==16 and unitrow["exact_Bezout_multiplier_degree"]==3 and
      unitrow["complete_zero_D_stratum_excluded"] and
      unitrow["assurance"]["full_polynomial_identity_checked_in_characteristic_zero"] and
      unitrow["assurance"]["all_complex_embeddings_covered"])
check("global polynomial root recovery preserves all coefficient-zero cases and nonzero saturation",
      unitrow["global_root_constant_polynomial_recovery"] and
      unitrow["assurance"]["no_parameter_denominator_introduced"] and
      unitrow["assurance"]["all_individual_coefficient_zero_cases_retained"] and
      unitrow["assurance"]["elimination_and_saturation_ideals_preserved"])
check("four-parameter core retains exact recovery circuits and finite degree bounds",
      unitrow["core_parameter_count"]==4 and unitrow["saturated_variable_count"]==5 and
      unitrow["root_constant_degree_upper_bound"]==18 and unitrow["compatibility_degree_upper_bound"]==22 and
      all(x["recovered_root_constant_circuit"]=="Zstar=sum(C[i]*E[i],i=0..4)" and x["nonzero_condition"]=="M(s)*Zstar != 0" for x in recovery["sectors"]))
check("unit-row result preserves the full necessary-core and actual finite-map obligations",
      unitrow["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(unitrow["assurance"][key] is False for key in
          ["four_parameter_nonzero_kappa_locus_solved","degree8_actual_section_lifts_constructed",
           "general_determinant_locus_solved","actual_global_root_hypothesis_removed","full_source_Keller_pair_constructed",
           "field_degree_four_realization_constructed","global_surface_finiteness_proved","actual_ruled_compatibility_proved",
           "contact_cases_classified","whole_marked_candidate_excluded","quartic_packets_excluded",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
quadratic=json.loads((HERE/"degree8_quadratic_receipt.json").read_text())
quadratic_obs=json.loads((HERE/"degree8_native_saturation_observations.json").read_text())
quadratic_models=json.loads(gzip.decompress((HERE/"degree8_quadratic_presentations.json.gz").read_bytes()))
check("exact quadratic presentation retains all16 embeddings with38 variables39 equations and27 recoverable auxiliaries",
      quadratic["marked_points"]==16 and quadratic["exact_quadratic_variable_count"]==38 and
      quadratic["exact_quadratic_equation_count"]==39 and quadratic["maximum_equation_degree"]==2 and
      quadratic["uniquely_recoverable_auxiliary_count"]==27 and quadratic["base_high_variable_count"]==10 and
      all(len(x["variables"])==38 and len(x["generators"])==39 for x in quadratic_models["models"]))
check("quadratic equivalence has no variable-denominator or parameter rank-drop loss",
      quadratic["nonzero_condition_preserved"] and all(quadratic["assurance"][key] for key in
          ["exact_presentation_covers_all_complex_embeddings","monic_elimination_has_no_variable_denominator",
           "complete_original_necessary_system_preserved","all_parameter_rank_drop_cases_retained"]))
check("native residue units match both full primary inputs and remain modular observations",
      quadratic["modular_unit_observations"]==2 and quadratic["modular_prime"]==1073741827 and
      quadratic["assurance"]["modular_results_not_promoted_to_characteristic_zero"] and
      len(quadratic_obs["observations"])==2 and all(x["actual_input_matches_primary_core"] and
          x["unit_basis_recorded"] and x["characteristic_zero_unit_identity_obtained"] is False
          for x in quadratic_obs["observations"]) and quadratic_obs["characteristic_zero_exclusion_proved"] is False)
check("quadratic presentation preserves open characteristic-zero and full global obligations",
      quadratic["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(quadratic["assurance"][key] is False for key in
          ["characteristic_zero_degree8_unit_certificate_obtained","full_degree8_nonzero_kappa_core_solved",
           "actual_global_section_lifts_constructed","general_determinant_locus_solved","actual_root_hypothesis_removed",
           "full_source_Keller_pair_constructed","field_degree_four_realization_constructed","global_surface_finiteness_proved",
           "actual_ruled_compatibility_proved","contact_cases_classified","whole_marked_candidate_excluded",
           "quartic_packets_excluded","formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
descent=json.loads((HERE/"degree8_quartic_descent_receipt.json").read_text())
descentcert=json.loads(gzip.decompress((HERE/"degree8_quartic_descent_certificates.json.gz").read_bytes()))
check("quartic fixed-field descent covers all16 octic embeddings with8 quartic embedding problems",
      descent["marked_octic_embeddings"]==16 and descent["quartic_embedding_problems"]==8 and
      descent["original_coefficient_field_degree"]==8 and descent["descended_coefficient_field_degree"]==4 and
      descent["rational_unknowns_per_field_unknown_before"]==8 and descent["rational_unknowns_per_field_unknown_after"]==4)
check("complete quadratic and core transport retain all parameters and nonzero conditions",
      descent["exact_quadratic_variable_count"]==38 and descent["exact_quadratic_equation_count"]==39 and
      descent["maximum_equation_degree"]==2 and descent["core_parameter_count"]==4 and
      all(descent["assurance"][key] for key in
          ["all_sixteen_octic_embeddings_covered","all_polynomial_coefficients_transported_exactly",
           "invertible_constant_changes_only","all_parameter_strata_retained","nonzero_kappa_preserved",
           "original_cusp_and_root_hypotheses_preserved","source_sign_automorphism_preserves_finite_map_obligations"]) and
      len(descentcert["sectors"])==2 and all(len(x["quadratic_generators"])==39 for x in descentcert["sectors"]))
check("descended unit rows are exact characteristic-zero identities without modular inference",
      descent["descended_unit_multiplier_degree"]==3 and
      descent["assurance"]["descended_Bezout_identity_checked_in_characteristic_zero"] and
      descent["assurance"]["descent_uses_no_modular_inference"] and
      all(x["unit_identity_verified"] and x["core_solved"] is False for x in descentcert["sectors"]))
check("coefficient-field reduction closes no full-core or finite-map obligation",
      descent["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(descent["assurance"][key] is False for key in
          ["characteristic_zero_unit_certificate_for_full_core_obtained","full_degree8_nonzero_kappa_core_solved",
           "actual_global_section_lifts_constructed","general_determinant_locus_solved","actual_root_hypothesis_removed",
           "full_source_Keller_pair_constructed","field_degree_four_realization_constructed","global_surface_finiteness_proved",
           "actual_ruled_compatibility_proved","contact_cases_classified","whole_marked_candidate_excluded",
           "quartic_packets_excluded","formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
linearM=json.loads((HERE/"degree8_linear_remainder_receipt.json").read_text())
check("constant and linear M including zero are excluded at all16 marked points by exact inverses",
      linearM["marked_points"]==16 and linearM["coefficient_fields"]==2 and
      linearM["excluded_critical_remainder_degrees"]==[0,1] and linearM["zero_polynomial_M_excluded"] and
      linearM["necessary_minimum_critical_remainder_degree"]==2 and linearM["necessary_minimum_D_polynomial_degree"]==1)
check("finite-rule memberships and full primary E5 circuits and inverse products are exact characteristic zero",
      linearM["exact_reduction_rule_count_per_field"]==35 and linearM["inverse_polynomial_degree"]==3 and
      all(linearM["assurance"][key] for key in ["full_characteristic_zero_polynomial_identities_checked",
          "complete_primary_necessary_E5_circuit_checked","all_sixteen_complex_markings_covered",
          "all_parameters_and_normal_freedoms_retained","no_parameter_denominator_introduced",
          "exclusion_precedes_nonzero_kappa_condition","modular_unit_flags_not_used_as_proof"]) and
      all(x["all_rule_memberships_verified"] and x["full_E5_circuit_checked"] and
          x["full_inverse_identity_checked_in_characteristic_zero"] for x in linearM["summaries"]))
check("linear-M exclusion preserves remaining degree8 and global finite-map obligations",
      linearM["first_open_reduced_pairs"]==[[5,9],[6,9],[7,9],[8,9]] and
      all(linearM["assurance"][key] is False for key in
          ["actual_global_root_hypothesis_removed","full_degree8_nonzero_kappa_core_solved",
           "general_determinant_locus_solved","actual_global_section_lifts_constructed","full_source_Keller_pair_constructed",
           "field_degree_four_realization_constructed","global_surface_finiteness_proved","actual_ruled_compatibility_proved",
           "contact_cases_classified","whole_marked_candidate_excluded","quartic_packets_excluded",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
global_boundary=json.loads((HERE/"global_boundary_receipt.json").read_text())
check("new gate has all-degree scope and explicitly imported primary theorems",
      global_boundary["scope"]=="all pole heights and every finite field degree on the actual marked Wright model" and
      global_boundary["assurance"]["explicit_imported_Chau_theorems"] and
      global_boundary["assurance"]["imported_results_independently_reproved"] is False and
      global_boundary["assurance"]["all_degree_scope_is_written_not_inferred_from_finite_controls"])
check("global projection, boundary collision, symmetric monodromy and quartic image conclusions are retained",
      all(global_boundary["global_conclusions"].values()))
check("uniform all-order control contracts a retained line and does not weaken fixed-pole cutoffs",
      all(global_boundary["uniform_jet_control"][key] for key in
          ["all_nonnegative_orders_covered_by_written_induction","source_Jacobian_vanishes_on_retained_X_zero",
           "every_approximant_contracts_the_retained_X_zero_line","polynomial_mate_to_D_excluded_in_every_degree",
           "uniform_finite_jet_sufficiency_refuted_without_pole_bound","fixed_pole_space_completion_cutoff_unaffected"]))
check("updated user priority is global and quadratic coefficient calculation is only historical",
      "general all-degree geometry" in action["next_step"] and
      global_boundary["assurance"]["global_result_uses_no_new_Groebner_search"] and
      global_boundary["assurance"]["quadratic_remainder_computation_checkpointed_as_historical"] and
      json.loads((HERE/"sources/quadratic-remainder-checkpoint/STATUS.json").read_text())["quadratic_M_stratum_excluded"] is False)
check("global gate leaves the colliding-image and compatible-completion problem open",
      all(global_boundary["assurance"][key] is False for key in
          ["actual_finite_map_constructed","all_mixed_weight_Keller_pairs_excluded","all_boundary_collisions_excluded",
           "quartic_surjectivity_proved","arbitrary_normalization_compatible_completion_proved","quartic_packets_excluded",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
residual=json.loads((HERE/"residual_cover_receipt.json").read_text())
check("all-degree residual conormal orientation and canonical sign-cover arguments are retained",
      residual["all_degree_identity"]=="[u^2]H(P,Q)=-kappa*rho/(2*h)" and
      residual["all_degree_residual_discriminant_squareclass"]=="beta=[u^2]H(P,Q)" and
      residual["necessary_finite_field_degree_at_least"]==4 and
      all(residual["assurance"][key] for key in
          ["all_degree_conormal_and_orientation_arguments_are_written",
           "all_degree_sign_pullback_canonical_sheaf_trivial",
           "all_degree_sign_map_unramified_in_codimension_one",
           "primary_orientation_and_Hurwitz_results_imported"]))
check("quartic residual passport covers all collision orders without node or pole bounds",
      residual["quartic_scope"]=="every pole height, collision count and positive tangency multiplicity" and
      all(residual["assurance"][key] for key in ["quartic_parity_cases_cover_all_positive_collision_orders",
          "no_transverse_node_assumption","no_pole_bound_or_coefficient_stratum_used",
          "residual_cover_is_not_an_intermediate_field_of_the_quartic_surface_extension"]))
passport=residual["quartic_passport"]
check("residual normalization genus holes source poles and divisor class remain exact",
      passport["connected_iff"]=="o>=1" and passport["compact_genus_if_connected"]=="o-1" and
      passport["punctures_if_connected"]=="4*n-2*o+3" and passport["total_Euler_characteristic"]=="1-4*n" and
      passport["finite_X_pole_orders"]=={"cusp_missing_point":2,"odd_collision_parameter":"m_i","even_collision_each_point":"m_i/2"} and
      passport["residual_divisor_class"]=="T=-2*K_N")
sign=residual["quartic_sign_cover"]
check("sign-cover finite degree codimension-one singularity and Euler constraints are explicit",
      sign["canonical_sheaf_trivial"] and sign["degree_over_target_sign_cover"]==4 and
      sign["unramified_in_codimension_one"] and
      sign["finite_singularities"]=="one A2; two A_(m_i-1) per collision, A0 smooth" and
      sign["Euler_characteristic"]=="2+2*n" and sign["minimal_resolution_Euler_characteristic"]=="4+2*sum(m_i)")
check("residual global reduction leaves actual embedding and whole goal open",
      residual["assurance"]["abstract_passport_is_not_a_global_realization"] and
      all(residual["assurance"][key] is False for key in
          ["actual_plane_embedding_constructed","actual_finite_map_constructed",
           "all_colliding_Wright_covers_excluded","quartic_packets_excluded",
           "arbitrary_normalization_compatible_completion_proved","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]) and
      "PF2RULED33" in action["claims"])
torus=json.loads((HERE/"torus_infinity_receipt.json").read_text())
check("torus-infinity exclusion preserves all-degree scope and actual marked end",
      torus["scope"]=="every pole height and every possible finite field degree on the actual marked Wright surface" and
      torus["marked_end"]=={"H1":"Z/3","ramification_class":"4*mu=mu!=0","B_square":3,"E_dot_B":4})
check("general Seifert lemma compares complete homology and orbifold quotient",
      torus["topological_lemma"]["degree_lower_bound"]==3 and
      torus["topological_lemma"]["monodromy"]=="S_d" and
      torus["topological_lemma"]["finite_H1_required"] and
      torus["topological_lemma"]["source_orders"]==["A","s*k","r*k"] and
      torus["topological_lemma"]["first_homology_order"]==torus["topological_lemma"]["orbifold_abelianization_order"]=="k" and
      torus["topological_lemma"]["ramification_knot_class"]==0)
check("torus and coprime-degree exclusions imply composite coordinate boundary degrees uniformly",
      all(torus["consequences"][key] for key in ["unknot_at_infinity_excluded",
          "all_torus_knots_at_infinity_excluded",
          "coprime_boundary_normalization_degrees_excluded_in_every_target_presentation",
          "every_prime_boundary_coordinate_degree_excluded"]) and
      torus["consequences"]["minimum_boundary_coordinate_degree"]==4)
check("unmarked lens control and balanced weighted Euler route prevent unsupported exclusions",
      torus["unmarked_T4_9_control"]=={"monodromy":"S4","H1":"Z/3","ramification_class":0,"is_actual_Wright_cover":False} and
      torus["rejected_euler_route"]["weighted_multiplicativity_balances_for_all_orders"] and
      torus["rejected_euler_route"]["exclusion_obtained"] is False)
check("all-degree end proof separates written topology and imports from finite controls",
      all(torus["assurance"][key] for key in
          ["actual_marked_ring_finite_map_and_constant_source_Jacobian_are_hypotheses",
           "all_degree_proof_is_written_not_inferred_from_finite_controls",
           "primary_topology_results_explicitly_imported",
           "branched_base_Euler_naturality_proved_by_section_obstruction",
           "finite_controls_do_not_replace_topological_proof","uses_no_new_coefficient_or_Groebner_search",
           "boundary_coordinate_degrees_distinguished_from_source_degrees_and_B_poles"]))
check("historical torus-only receipt retains its original cabling and full-objective limits",
      all(torus["assurance"][key] is False for key in
          ["all_iterated_cable_knots_excluded","all_colliding_Wright_covers_excluded",
           "actual_finite_map_constructed","arbitrary_normalization_compatible_completion_proved",
           "whole_height_nine_pair_excluded","quartic_packets_excluded",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]) and
      "PF2RULED34" in action["claims"])
cabling=json.loads((HERE/"quartic_cabling_receipt.json").read_text())
check("new global quartic exclusion has all-pole-height scope on the actual marked ring",
      cabling["scope"]=="all pole heights on the actual marked Wright surface; finite field degree exactly four" and
      cabling["actual_marked_quartic_Keller_maps_excluded"] and
      cabling["compatible_ruled_quartic_realizations_excluded"] and
      cabling["actual_marked_quartic_realization_frontier"]==[])
check("actual end core and residual-infinity preferred longitude are retained",
      cabling["marked_end"]=={"H1":"Z/3","ramification_exterior":"solid torus","ramification_class":"generator"} and
      cabling["peripheral"]["S4_preferred_longitude"]=="identity" and
      cabling["peripheral"]["A4_preferred_longitude"]=="identity")
layer=cabling["finite_peripheral_layer"]
check("finite permutation layer is exhaustive modulo12 with all companion image states",
      layer["exhaustive_S4_subgroups"]==30 and layer["exponent_modulus"]==12 and
      layer["S_fiber_trivial_state_count"]==9 and layer["A_fiber_trivial_state_count"]==4 and
      layer["S_residue_controls"]==208 and layer["A_residue_controls"]==96)
check("written cabling descent tracks unbounded depth and actual filling lattices",
      all(cabling["assurance"][key] for key in
          ["all_positive_winding_and_signed_coprime_cable_parameters_covered_by_written_proof",
           "finite_cabling_depth_unbounded_and_handled_by_descent",
           "residual_infinity_used_to_prove_preferred_longitude_identity",
           "compression_ball_alternative_excluded_using_nonzero_marked_class",
           "actual_peripheral_filling_lattices_tracked","transparent_step_preserves_lens_core_and_H1_three",
           "primary_topology_results_explicitly_imported","end_lens_recognition_explicitly_imported","finite_permutation_layer_exhaustive_but_topological_proof_written"]))
check("global realization obstruction does not fabricate coefficient unit ideals or full packet closure",
      all(cabling["assurance"][key] is False for key in
          ["coefficient_search_used","all_finite_field_degrees_on_Wright_excluded",
           "arbitrary_normalization_compatible_completion_proved",
           "arbitrary_quartic_normalizations_excluded","quartic_packets_excluded",
           "truncated_degree8_necessary_core_unit_ideal_proved",
           "formal_verification","independent_full_reproof","external_review","research_goal_completed"]))
check("historical marked quartic handoff remains retained without reopening its coefficient frontier",
      "PF2RULED35" in action["claims"] and cabling["assurance"]["quartic_packets_excluded"] is False)
intrinsic=json.loads((HERE/"intrinsic_normalization_receipt.json").read_text())
check("all-degree one-prime valuation theorem separates local nilpotence and uniform annihilation",
      intrinsic["valuation_module"]==dict(annihilator="symbolic power p^(-v(b))",
          associated_primes=["p"],support="V(p)",uniform_annihilator=0,
          p_action_locally_nilpotent=True,prime_hypothesis_automatic=True))
check("intrinsic quartic residual infinity gives longitude without a Wright marking",
      intrinsic["intrinsic_quartic"]["residual_infinity_points"]==2 and
      intrinsic["intrinsic_quartic"]["residual_infinity_unramified"] and
      intrinsic["intrinsic_quartic"]["preferred_longitude"]=="identity" and
      intrinsic["intrinsic_quartic"]["Wright_marking_used"] is False)
check("general odd three-divisible lens obstruction keeps the primitive core hypothesis",
      intrinsic["excluded_end"]["homology_order"]=="n odd and divisible by 3" and
      intrinsic["excluded_end"]["marked_class"]=="primitive generator" and
      intrinsic["excluded_end"]["all_finite_cabling_depths"])
check("reduced ruling canonical obstruction distinguishes multiple fibers",
      all(intrinsic["ruling"][k] for k in
         ["multiplicity_one_quartic_excluded","negative_coordinate_pole_one_excluded",
          "all_degree_reduced_simple_boundary_torus_end_excluded"]) and
      intrinsic["ruling"]["canonical_linear_bundle_degree"]==-3)
check("intrinsic reduction retains full comparison, multiplicity and assurance gaps",
      all(intrinsic["assurance"][k] is False for k in
          ["full_three_manifold_proof_implemented","coordinate_control_is_finite_Keller_extension",
           "higher_index_single_prime_quartic_excluded","arbitrary_quartic_normalizations_excluded",
           "multiple_fiber_rulings_excluded","all_finite_Wright_field_degrees_excluded",
           "general_lens_core_end_proved","coefficient_search_used","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]) and
      intrinsic["assurance"]["reduced_boundary_not_confused_with_reduced_fiber"])
check("historical intrinsic multiplicity and end handoff remains retained",
      "PF2RULED36" in action["claims"])
cyclic=json.loads((HERE/"multiple_fiber_cyclic_receipt.json").read_text())
check("cyclic reduction preserves ramification and fixes all-degree torsor index",
      cyclic["general_smooth_geometry"]["linear_torsor_index"]=="q+e" and
      cyclic["general_smooth_geometry"]["ramification_of_cyclic_cover"]=="(q-1)L_prime" and
      cyclic["general_smooth_geometry"]["canonical_class"]=="(q+e-2)*sum R_i" and
      cyclic["general_smooth_geometry"]["split_boundary_lines"]=="q")
check("one shared residue detects global volume independently of multiplicity",
      cyclic["volume"]["residue"]=="-q*[z^e]g_i" and
      cyclic["volume"]["deck_shared_scalar"] and
      cyclic["volume"]["exact_if_and_only_if_shared_coefficient_zero"] and
      cyclic["volume"]["actual_index_two_necessary_coefficient"]=="[z^2]g_i=0" and
      cyclic["volume"]["tested_exponent_independent_of_q"])
check("extra immersed branch curve and wreath monodromy prevent a simple-cover transfer",
      cyclic["actual_finite_cover"]["monodromy"]=="C_q wreath S_d" and
      cyclic["actual_finite_cover"]["group_order"]=="q^d*d!" and
      cyclic["actual_finite_cover"]["additional_curve_distinct_from_C"] and
      cyclic["actual_finite_cover"]["additional_curve_immersed_normalization"] and
      cyclic["actual_finite_cover"]["source_Jacobian"]=="c*q*z^(q-1)")
check("explicit affine controls retain exact volume at every multiple fiber",
      all(cyclic["affine_control_family"][k] for k in
          ["all_positive_integer_multiplicities","smooth","one_omitted_A1",
           "exact_volume_for_every_q_at_least_two"]) and
      cyclic["affine_control_family"]["finite_Keller_extension_constructed"] is False and
      cyclic["affine_control_family"]["general_normalization_classification"] is False)
check("cyclic proof and controls keep full finite-algebra and end gates open",
      cyclic["assurance"]["all_q_geometry_proved_in_writing"] and
      cyclic["assurance"]["ordinary_ramification_coefficients_preserved"] and
      all(cyclic["assurance"][k] is False for k in
          ["whole_multiple_fiber_locus_excluded","cyclic_cover_mistaken_for_Keller",
           "restriction_to_one_piece_asserted_finite","actual_end_lens_or_core_proved",
           "arbitrary_quartic_packets_excluded","all_finite_Wright_degrees_excluded",
           "coefficient_or_Groebner_search_used","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]))
check("historical all-degree cyclic gluing handoff remains retained",
      "PF2RULED37" in action["claims"])
coordinate=json.loads((HERE/"coordinate_regular_sector_receipt.json").read_text())
check("coordinate minimum is intrinsic over all source polynomial automorphisms",
      coordinate["coordinate_invariant"]["regular_coordinates"]=="a*y+b with a nonzero" and
      coordinate["coordinate_invariant"]["negative_coordinate_poles"]=="q*positive integers" and
      coordinate["coordinate_invariant"]["minimum_pole"]=="q" and
      coordinate["coordinate_invariant"]["all_polynomial_source_automorphisms_covered_by_written_normal_form"] and
      coordinate["coordinate_invariant"]["bounded_word_controls_are_not_the_proof"])
check("whole simple regular-coordinate sector has an affine marked normal form",
      coordinate["regular_sector"]["complete_marked_normal_form"] and
      coordinate["regular_sector"]["source_chart"]=="x=u*t+h(u), y=u" and
      coordinate["regular_sector"]["marked_equivalence"]=="same q and c_j_prime=b*alpha^j*c_j")
check("the full principal-part volume criterion has a constructive uniform primitive",
      coordinate["regular_sector"]["exact_volume_iff"]=="c_1=0" and
      coordinate["regular_sector"]["primitive_constructed_for_all_c1_zero"] and
      coordinate["regular_sector"]["negative_part_inverse_determinant"]=="c_q^(q-1)")
check("cyclic geometric controls have no hidden reduced coordinate ruling",
      all(coordinate["sharpened_controls"][k] for k in
          ["every_q_at_least_two","no_alternative_source_coordinate_pole_one","no_marked_Wright_comparison"]) and
      coordinate["sharpened_controls"]["finite_Keller_extension_constructed"] is False)
check("actual finite algebra requires the monic norm and horizontal conductor",
      coordinate["finite_algebra"]["actual_finiteness_is_a_hypothesis"] and
      coordinate["finite_algebra"]["primitive_integral_coordinate"]=="K(y)=C(N)" and
      coordinate["finite_algebra"]["monic_norm_degree"]=="d" and
      coordinate["finite_algebra"]["every_D_beta_self_colliding"] and
      coordinate["finite_algebra"]["target_plus_y_algebra_not_normal"] and
      coordinate["finite_algebra"]["conductor_has_y_horizontal_component"] and
      coordinate["finite_algebra"]["norm_family_is_a_linear_pencil"] is False)
check("marked regular-sector classification preserves all unclosed global cases",
      all(coordinate["assurance"][k] for k in
          ["all_q_and_all_degree_proofs_written","primary_normal_form_and_fibration_imports_explicit",
           "one_line_injectivity_import_explicit","exact_volume_does_not_construct_finite_Keller_map"]) and
      all(coordinate["assurance"][k] is False for k in
          ["actual_regular_order_one_sector_excluded","regular_coordinate_higher_vanishing_classified",
           "no_regular_coordinate_case_classified","arbitrary_quartic_packets_excluded",
           "all_finite_Wright_degrees_excluded","generic_end_lens_or_core_proved",
           "coefficient_or_Groebner_search_used","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]))
check("historical coordinate-invariant finite-algebra handoff remains retained",
      "PF2RULED38" in action["claims"])
end=json.loads((HERE/"regular_sector_end_receipt.json").read_text())
check("whole marked simple-regular end is a lens with a primitive core",
      end["end"]["H1"]=="Z/(2q+1)" and
      end["end"]["ramification_exterior"]=="solid torus" and
      end["end"]["ramification_class"]=="primitive generator" and
      end["end"]["pole_recovered_from_abstract_end"]=="q=(|H1|-1)/2" and
      end["end"]["q_at_least_two_abstractly_not_Wright"])
check("new cabling proof covers every odd order with alternating descent",
      end["strengthened_topology"]["homology_order"]=="every odd n>1" and
      all(end["strengthened_topology"][k] for k in
          ["simple_S4_identity_longitude_lens_core_type_excluded",
           "alternating_A4_lens_core_type_excluded",
           "fourth_letter_pattern_excluded_by_marked_exterior",
           "new_A_to_A_removal_preserves_actual_boundary_homology",
           "depth_unbounded_written_descent"]) and
      end["strengthened_topology"]["n_one_excluded_by_this_theorem"] is False)
check("new global obstruction closes every simple-regular quartic principal part",
      all(end["algebraic_consequences"][k] for k in
          ["whole_simple_regular_coordinate_quartic_sector_excluded","all_q_all_principal_parts",
           "all_degree_torus_infinity_excluded","all_degree_prime_boundary_coordinate_degrees_excluded"]) and
      end["algebraic_consequences"]["minimum_boundary_coordinate_degree"]==4 and
      end["algebraic_consequences"]["minimum_possible_finite_field_degree"]==5 and
      end["algebraic_consequences"]["degree_five_realization_claimed"] is False)
check("locally valid odd-order removals are retained instead of falsely excluded",
      end["locally_consistent_removals"]==[
          dict(kind="A",p=2,a=1,n=5,fiber=[1,6],filling=[1,5],marked_class=2,is_actual_global_cover=False),
          dict(kind="S",p=3,a=1,n=5,fiber=[1,6],filling=[1,5],marked_class=3,is_actual_global_cover=False)])
check("punctured topology and actual filling lattices are explicit proof inputs",
      all(end["assurance"][k] for k in
          ["whole_marked_compactification_checked_in_writing",
           "actual_omitted_knot_exterior_tracked_by_punctured_graph",
           "primary_plumbing_diagrams_visually_inspected",
           "not_just_intersection_determinant_or_homology_recognition",
           "primary_plumbing_import_explicit",
           "actual_radial_end_recognition_uses_retained_primary_imports",
           "all_odd_cabling_proof_written",
           "finite_permutation_layer_unchanged_and_exhaustive",
           "ball_alternative_excluded_by_essential_marked_class",
           "all_replacements_track_actual_boundary_lattice",
           "equal_H1_not_used_to_identify_closed_manifolds"]))
check("whole simple-regular closure preserves every remaining full-goal limit",
      all(end["assurance"][k] is False for k in
          ["all_finite_field_degrees_excluded","higher_regular_coordinate_order_excluded",
           "no_regular_coordinate_case_excluded","arbitrary_quartic_packets_excluded",
           "coefficient_or_Groebner_search_used","finite_Keller_map_constructed",
           "formal_verification","independent_full_reproof","external_review",
           "research_goal_completed"]))
check("historical simple-regular quartic handoff remains retained",
      "PF2RULED39" in action["claims"])
higher=json.loads((HERE/"higher_regular_coordinate_receipt.json").read_text())
check("higher regular-coordinate form uses entire planes on an existing affine cyclic root",
      higher["cyclic_form"]["root"]=="y=z^r" and
      higher["cyclic_form"]["smooth_affine"] and
      higher["cyclic_form"]["entire_plane_chart_count"]=="r+1" and
      higher["cyclic_form"]["linear_exponent"]=="2-r" and
      higher["cyclic_form"]["minimum_pole_inequality"]=="q>=r-1" and
      higher["cyclic_form"]["support_separation"]=="gcd(r,support(h))=1")
check("one common residue excludes the whole diagonal in every finite degree",
      higher["volume"]["common_residue"]=="r*c_r" and
      higher["volume"]["exact_iff"]=="c_r=0" and
      higher["volume"]["all_finite_degree_diagonal_excluded"]=="q=r" and
      higher["volume"]["written_affine_Cech_sufficiency"] and
      higher["volume"]["invariant_primitive_descends_across_branch_line"])
check("extended source translation has its exact fixed divisor without target invariance",
      higher["additive_action"]["locally_nilpotent_on_normalization"] and
      higher["additive_action"]["free_when"]=="r=2" and
      higher["additive_action"]["zero_divisor_when_r_at_least_three"]=="(r-2)R" and
      higher["additive_action"]["target_algebra_invariance_proved"] is False)
check("actual finite algebra retains the repeated special norm and cyclic wreath lift",
      higher["finite_algebra"]["actual_one_index_two_prime_is_an_extra_hypothesis"] and
      higher["finite_algebra"]["special_norm"]=="H_C^r*H_0" and
      higher["finite_algebra"]["nonzero_norm"]=="H_beta reduced irreducible" and
      higher["finite_algebra"]["primitive_integral_coordinate"]=="K(y)=C(N)" and
      higher["finite_algebra"]["conductor_has_y_horizontal_component"] and
      higher["finite_algebra"]["cyclic_lift_monodromy"]=="C_r wreath S_d" and
      higher["finite_algebra"]["norm_family_is_a_linear_pencil"] is False and
      higher["finite_algebra"]["lifted_cover_is_Keller"] is False and
      higher["finite_algebra"]["deleting_boundary_components_preserves_finiteness"] is False)
check("new reduction retains the formal borderline and every global assurance limit",
      higher["controls"]["exact_support_cases"]==567 and
      any(c["q"]==c["r"]-1 and c["passes_necessary_chart_and_volume_tests"] for c in
          higher["controls"]["surviving_formal_data"]) and
      all(higher["assurance"][k] for k in
          ["all_r_geometric_proof_written","primary_fibration_import_explicit",
           "existing_affine_normalization_is_a_hypothesis","exact_controls_not_all_r_geometric_proof"]) and
      all(higher["assurance"][k] is False for k in
          ["affine_converse_for_arbitrary_chart_data","finite_Keller_map_constructed",
           "whole_higher_regular_sector_excluded","no_regular_coordinate_case_excluded",
           "arbitrary_quartic_packets_excluded","all_finite_Wright_degrees_excluded",
           "coefficient_or_Groebner_search_used","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]))
check("every borderline r has an explicit affine smooth quotient with regular primitive",
      higher["borderline_controls"]["all_r_affine_construction_written"] and
      higher["borderline_controls"]["smooth_quotient"] and
      higher["borderline_controls"]["minimum_coordinate_pole"]=="q=r-1" and
      higher["borderline_controls"]["volume_divisor"]=="R" and
      higher["borderline_controls"]["regular_primitive"]=="x dy" and
      higher["borderline_controls"]["hidden_pole_one_when_r_at_least_three"] is False and
      higher["borderline_controls"]["finite_Keller_target_constructed"] is False)
check("historical all-r regular-coordinate handoff remains retained",
      "PF2RULED40" in action["claims"])
borderline=json.loads((HERE/"borderline_regular_end_receipt.json").read_text())
check("whole borderline compactification tracks the internal puncture and actual surgery slope",
      borderline["end"]["all_r_entire_compactification_written"] and
      borderline["end"]["H1"]=="Z/(r(r-1)+1)" and
      borderline["end"]["puncture"]=="first negative-two vertex" and
      borderline["end"]["exterior"]=="T(r-1,r) knot exterior" and
      borderline["end"]["filling_slope"]=="(r(r-1)+1)mu+lambda" and
      borderline["end"]["marked_class_primitive"] and
      borderline["end"]["Heegaard_core_iff"]=="r=2")
check("radial exterior comparison uses the primary rigidity theorem with full hypotheses",
      borderline["end"]["actual_radial_exterior_recognized_using_Waldhausen"] and
      all(borderline["assurance"][k] for k in
          ["primary_rigidity_import_explicit","primary_rigidity_statement_visually_inspected",
           "actual_marked_exterior_not_inferred_from_H1"]))
check("quartic odd-lens theorem extends from cores to primitive atoroidal markings",
      borderline["strengthened_topology"]["homology_order"]=="every odd n>1" and
      all(borderline["strengthened_topology"][k] for k in
          ["full_S4","simple_branch","identity_longitude","source_lens",
           "primitive_ramification","atoroidal_ramification_exterior",
           "such_cover_excluded","exact_old_permutation_layer_reused",
           "primitive_regular_fiber_quotient_argument",
           "genuine_cable_torus_boundary_parallel_exception_checked",
           "final_S3_torus_companion_forced_trefoil",
           "actual_trefoil_and_pattern_boundary_lattices_tracked"]) and
      borderline["strengthened_topology"]["Heegaard_core_required"] is False)
check("both mirrors in the final companion case have a noncyclic third cone",
      len(borderline["controls"]["signed_trefoil_last_cases"])==4 and
      all(c["third_cone"] in [3,5] and c["lens_possible"] is False
          for c in borderline["controls"]["signed_trefoil_last_cases"]))
check("real order-seven cable cover is retained with its wrong marked exterior",
      borderline["countercontrol"]["actual_topological_simple_S4_cover_constructed"] and
      borderline["countercontrol"]["lens_H1_order"]==7 and
      borderline["countercontrol"]["primitive_marked_class"]==2 and
      borderline["countercontrol"]["marked_exterior_extra_cone"]==5 and
      borderline["countercontrol"]["atoroidal"] is False and
      borderline["countercontrol"]["marked_N3_end_realized"] is False and
      borderline["countercontrol"]["finite_Keller_map_constructed"] is False)
check("whole unbounded borderline quartic family is closed and the remaining strip is explicit",
      all(borderline["consequences"][k] for k in
          ["whole_borderline_quartic_sector_excluded","no_coefficient_or_pole_bound",
           "all_degree_torus_infinity_excluded","all_degree_prime_boundary_coordinate_degree_excluded"]) and
      borderline["consequences"]["remaining_actual_quartic_regular_sector"]=="r>=2 and q>=r+1" and
      borderline["consequences"]["higher_finite_degrees_excluded"] is False)
check("atoroidal marked-end proof preserves all wider full-goal limits",
      all(borderline["assurance"][k] for k in
          ["all_r_and_all_depth_proofs_written","no_target_equivariance_assumed",
           "finite_controls_not_full_geometric_proof"]) and
      all(borderline["assurance"][k] is False for k in
          ["affine_converse_for_arbitrary_Laurent_data","whole_higher_regular_sector_excluded",
           "no_regular_coordinate_case_excluded","extra_omitted_primes_excluded",
           "arbitrary_quartic_packets_excluded","all_finite_Wright_degrees_excluded",
           "finite_Keller_map_constructed","formal_verification","independent_full_reproof",
           "external_review","research_goal_completed"]))
check("historical closed borderline family and its original countercontrols remain retained",
      "PF2RULED41" in action["claims"])
coprime=json.loads((HERE/"coprime_regular_end_receipt.json").read_text())
check("whole mixed Laurent chart and marked contact are retained in the coprime sector",
      coprime["compactification"]["all_lower_Laurent_coefficients_retained"] and
      coprime["compactification"]["free_blowups"]=="m=q-r+2" and
      coprime["compactification"]["full_omitted_line_checked"] and
      coprime["compactification"]["arithmetic_control_is_not_unbounded_proof"])
check("coprime marked end is the actual surgery pair with both atoroidality properties",
      coprime["end"]["surgery"]=="n=qr+q-r+2 surgery on T(q,r)" and
      coprime["end"]["H1"]=="Z/n" and coprime["end"]["marked_class_primitive"] and
      coprime["end"]["closed_Seifert_orders"]==["q","r","q-r+2"] and
      coprime["end"]["closed_source_atoroidal"] and coprime["end"]["marked_exterior_atoroidal"])
check("actual radial source is recognized from the marked exterior without closed-group classification",
      coprime["end"]["radial_irreducibility_from_nontrivial_center"] and
      coprime["end"]["actual_radial_pair_recognized_via_peripheral_Waldhausen"] and
      coprime["end"]["group_only_closed_Seifert_classification_assumed"] is False)
check("new quartic theorem uses closed atoroidality and the strict homology threshold",
      coprime["strengthened_topology"]["closed_source_hypothesis"]=="irreducible and atoroidal" and
      coprime["strengthened_topology"]["source_lens_required"] is False and
      coprime["strengthened_topology"]["source_H1"]=="Z/n, n>3 odd" and
      all(coprime["strengthened_topology"][k] for k in
          ["full_S4","simple_branch","identity_longitude","primitive_ramification",
           "atoroidal_ramification_exterior","such_cover_excluded",
           "unchanged_exhaustive_permutation_layer","repeated_cone_quotient_all_fillings",
           "replacement_homology_preserved_without_distance_one_assumption"]) and
      coprime["strengthened_topology"]["last_trefoil_cases_H1_orders"]==[1,3])
check("an actual order-three topological cover prevents a false stronger theorem",
      coprime["sharpness_control"]["actual_topological_cover"] and
      coprime["sharpness_control"]["source_H1_order"]==3 and
      coprime["sharpness_control"]["source_Seifert_orders"]==[2,3,3] and
      coprime["sharpness_control"]["closed_source_atoroidal"] and
      coprime["sharpness_control"]["marked_exterior_atoroidal"] and
      coprime["sharpness_control"]["source_lens"] is False and
      coprime["sharpness_control"]["polynomial_infinity_realization_claimed"] is False and
      coprime["sharpness_control"]["finite_Keller_map_constructed"] is False)
check("all coprime monomial controls are actual affine geometry with a regular Euler primitive",
      coprime["affine_controls"]["entire_monomial_coprime_family_constructed"] and
      coprime["affine_controls"]["smooth_affine_quotient"] and
      coprime["affine_controls"]["global_regular_primitive"]=="(-q*x dy-r*y dx)/(r-q)" and
      coprime["affine_controls"]["no_target_equivariance_assumed"] and
      coprime["affine_controls"]["finite_Keller_map_constructed"] is False)
check("all coprime quartic pole/order pairs close with the exact arithmetic next region",
      coprime["consequences"]["whole_coprime_quartic_sector_excluded"] and
      coprime["consequences"]["all_q_r_all_Laurent_parts"] and
      coprime["consequences"]["remaining_quartic_regular_region"]=="r composite, q>r, gcd(q,r)>1" and
      coprime["consequences"]["sharpened_minimum_pole"]=="q>=r+smallest_prime_divisor(r)" and
      coprime["consequences"]["minimum_gap"]==2 and
      coprime["consequences"]["prime_r_requires"]=="q>=2r" and
      coprime["consequences"]["all_degree_torus_infinity_excluded"] and
      coprime["consequences"]["all_degree_prime_boundary_degree_excluded"] and
      coprime["consequences"]["all_finite_degrees_excluded"] is False)
check("coprime proof retains every global unresolved assurance boundary",
      all(coprime["assurance"][k] for k in
          ["all_parameter_and_all_depth_proofs_written",
           "primary_closed_horizontal_surface_criterion_visually_checked",
           "primary_imports_explicit","actual_peripheral_data_tracked",
           "finite_controls_not_full_geometric_proof"]) and
      all(coprime["assurance"][k] is False for k in
          ["affine_converse_for_arbitrary_Laurent_data","noncoprime_regular_sector_excluded",
           "no_regular_coordinate_case_excluded","extra_omitted_primes_excluded",
           "arbitrary_quartic_packets_excluded","all_finite_Wright_degrees_excluded",
           "finite_Keller_map_constructed","formal_verification",
           "independent_full_reproof","external_review","research_goal_completed"]))
check("integral leading powers retain the original boundary and close every prime order",
      all(coprime["single_characteristic_extension"][k] for k in
          ["integral_power_removal_written","all_lower_terms_retained",
           "original_deleted_boundary_attachment_tracked",
           "all_prime_boundary_orders_quartic_excluded",
           "all_single_characteristic_quartic_models_excluded"]) and
      coprime["single_characteristic_extension"]["remaining_boundary_order"]=="r composite, r>=4" and
      coprime["single_characteristic_extension"]["remaining_reference_germ"]=="at least two characteristic pairs" and
      coprime["single_characteristic_extension"]["remaining_translation_zero_order_at_least"]==2)
check("historical coprime marked-end handoff remains retained without reopening the closed noncoprime case",
      "PF2RULED42" in action["claims"])
universal=json.loads((HERE/"universal_regular_end_receipt.json").read_text())
check("universal affine converse closes every support using both transition directions",
      all(universal["affine_classification"][k] for k in
          ["necessary_and_sufficient","all_Laurent_supports","both_directions_all_pairs_polar",
           "smooth_affine_quotient","all_polynomial_source_shears_extend"]) and
      universal["affine_classification"]["Dubouloz_section"]=="2.8" and
      universal["affine_classification"]["marked_equivalence"]=="c'_j=b*a^j*c_j" and
      universal["affine_classification"]["volume_exact_iff"]=="c_r=0" and
      universal["affine_classification"]["target_algebra_invariance_assumed"] is False)
check("full characteristic-chain end gives primitive odd marked homology without false recognition",
      universal["universal_end"]["surgery"]=="n=2delta+2q+1" and
      all(universal["universal_end"][k] for k in
          ["n_odd","marked_class_primitive","plumbing_pair_identified",
           "stable_peripheral_data_transfer_to_radial_section"]) and
      universal["universal_end"]["arbitrary_multicharacteristic_radial_homeomorphism_claimed"] is False)
check("terminal distance-two gap and the distance-one companion bound are explicit",
      universal["terminal_structure"]["distance_two_impossible"] and
      universal["terminal_structure"]["multi_distance_one_companion_filling"]=="n/p_last^2" and
      universal["terminal_structure"]["next_fiber_distance_at_least"]==7 and
      universal["terminal_structure"]["no_consecutive_distance_one_collapse"] and
      universal["terminal_structure"]["actual_radial_JSJ_recognition_claimed"] is False)
check("all-degree target boundary consequences cover every regular-coordinate support",
      all(universal["all_degree_consequences"][k] for k in
          ["entire_regular_coordinate_sector","target_torus_or_unknot_infinity_excluded",
           "target_coprime_boundary_degrees_excluded","target_prime_boundary_degrees_excluded"]) and
      universal["all_degree_consequences"]["all_Keller_maps_excluded"] is False)
check("new scope closes geometry while preserving finite-map and assurance gaps",
      all(universal["assurance"][k] for k in
          ["written_unbounded_proofs","primary_imports_explicit",
           "finite_checks_not_complete_geometric_proof","general_affine_converse_closed"]) and
      all(universal["assurance"][k] is False for k in
          ["multi_characteristic_quartic_sector_excluded","no_regular_coordinate_case_excluded",
           "extra_omitted_primes_excluded","finite_Keller_map_constructed",
           "formal_verification","independent_full_reproof","external_review",
           "novelty_priority","research_goal_completed"]))
check("historical universal affine and plumbing handoff remains retained",
      "PF2RULED43" in action["claims"])

regular44=json.loads((HERE/"regular_quartic_exclusion_receipt.json").read_text())
check("the full one-index-two-prime regular-coordinate quartic sector is now excluded",
      all(regular44["theorem"][k] for k in
          ["all_regular_source_coordinate_quartic_models_excluded","one_omitted_prime",
           "finite_quartic_Keller_map","all_Laurent_coefficients","all_cabling_depths"]) and
      regular44["theorem"]["generic_index"]==2 and
      regular44["theorem"]["remaining_one_index_two_prime_requires"]=="every source coordinate has a pole at R")
check("actual radial recognition is proved before transferring source fibers or JSJ",
      all(regular44["all_degree_geometry"][k] for k in
          ["radial_closed_pair_recognized","radial_marked_exterior_recognized",
           "irreducibility_proved_before_Waldhausen",
           "full_peripheral_system_and_omitted_meridian_preserved",
           "source_fiber_slopes_and_JSJ_transfer_valid"]) and
      regular44["all_degree_geometry"]["H1_only_recognition_assumed"] is False and
      regular44["all_degree_geometry"]["arbitrary_abstract_group_recognition_assumed"] is False)
check("local simple-ramification multiplicity bound does not assume trivial longitude",
      regular44["all_degree_fiber_bound"]["simple_ramification_multiplicity"]==[1,2] and
      regular44["all_degree_fiber_bound"]["preferred_longitude_identity_required"] is False and
      regular44["all_degree_fiber_bound"]["nontrivial_fiber_branched_orbit_need_not_have_size_two"] and
      regular44["all_degree_fiber_bound"]["essential_nonproduct_piece_forces_kappa_one"])
check("all-depth auxiliary cover lemmas retain growth and companion alternatives",
      regular44["exterior_lemmas"]["no_unbranched_degree_four_between_nontrivial_iterated_torus_knot_exteriors"] and
      regular44["exterior_lemmas"]["source_chain_leaf_reflection_explicitly_excluded"] and
      regular44["exterior_lemmas"]["S3_C2_companions_retained"] and
      regular44["exterior_lemmas"]["finite_H1_not_assumed_for_solid_partial_three_sheeted_lemma"] and
      regular44["exterior_lemmas"]["solid_partial_three_sheeted_cover_under_absolute_algebraic_growth"]=="trefoil or mirror")
check("every quartic peripheral state is covered and the A4 replacement has its actual meridian",
      regular44["quartic_reduction"]["terminal_distance_at_least_three_excluded"] and
      regular44["quartic_reduction"]["distance_one_all_cases_excluded"] and
      regular44["quartic_reduction"]["product_exception_excluded_by_solid_three_sheeted_lemma"] and
      len(regular44["quartic_reduction"]["S4_decisions"])==9 and
      regular44["quartic_reduction"]["A4_intrinsic_fiber_distance_bound"]==[1,3] and
      regular44["quartic_reduction"]["A4_replacement_minimum_distance"]==7 and
      regular44["quartic_reduction"]["identity_companion_forces"]=="p_target=g_source")
check("source peeling constructs geometry and rational slopes without descending a Keller algebra",
      regular44["source_peeling"]["all_degree"] and
      regular44["source_peeling"]["shortened_model_actual_affine"] and
      regular44["source_peeling"]["shortened_terminal_distance_at_least"]==3 and
      regular44["source_peeling"]["new_companion_meridian"]=="n*mu_prime+g^2*lambda_prime" and
      regular44["source_peeling"]["finite_Keller_map_descended"] is False)
check("new actual controls preserve the necessity of growth and the marked exterior",
      regular44["sharpness_controls"]["non_growth_S3"]["actual_topological_cover"] and
      regular44["sharpness_controls"]["non_growth_S3"]["partial_exterior_solid"] and
      regular44["sharpness_controls"]["non_growth_S3"]["algebraic_growth"] is False and
      regular44["sharpness_controls"]["retained_S4"]["actual_topological_cover"] and
      regular44["sharpness_controls"]["retained_S4"]["iterated_S3_knot_exterior_condition"] is False)
check("whole regular quartic exclusion preserves all other global obligations",
      all(regular44["assurance"][k] for k in
          ["written_unbounded_geometric_and_all_case_proofs","primary_imports_explicit",
           "finite_controls_not_full_topological_proof",
           "source_reference_and_target_infinity_knots_distinguished",
           "compression_ball_alternatives_handled","all_solid_and_product_degenerations_handled",
           "no_automatic_target_algebra_invariance_assumed"]) and
      all(regular44["assurance"][k] is False for k in
          ["coefficient_search_used","no_regular_coordinate_case_excluded",
           "additional_omitted_primes_excluded","higher_index_packets_excluded",
           "all_finite_Keller_degrees_excluded","finite_Keller_map_constructed",
           "formal_verification","independent_full_reproof","external_review",
           "novelty_priority","research_goal_completed"]))
check("historical coordinate-negative handoff is retained while the new bridge closes it",
      "PF2RULED44" in action["claims"] and "PF2RULED45" in action["claims"] and
      "higher finite field degrees" in action["next_step"])

bridge45=json.loads((HERE/"exact_coordinate_bridge_receipt.json").read_text())
check("general canonical defect applies to every index without primitive compactification",
      bridge45["general_index_bound"]["all_indices"] and
      bridge45["general_index_bound"]["normal_form_coprime_defect"]=="g*(a*b-a-b)<=e-1" and
      bridge45["general_index_bound"]["maximum_reduced_pole"]=="3*(e-1)" and
      bridge45["general_index_bound"]["minimum_reduced_pole"]=="2*(e-1)" and
      bridge45["general_index_bound"]["e1_regular_coordinate_automatic"] and
      bridge45["general_index_bound"]["higher_index_cases_excluded"] is False and
      bridge45["assurance"]["primitive_compactification_assumed"] is False)
check("all-degree exact coordinate bridge compares residue with the entire source contraction ring",
      bridge45["exact_coordinate_bridge"]["regular_source_coordinate_exists"] and
      bridge45["exact_coordinate_bridge"]["all_finite_field_degrees"] and
      bridge45["exact_coordinate_bridge"]["global_exactness_required"] and
      bridge45["exact_coordinate_bridge"]["only_coordinate_negative_initial_pair"]==[3,2] and
      bridge45["exact_coordinate_bridge"]["essential_key_sequence"]==[3,2,-1] and
      bridge45["exact_coordinate_bridge"]["residue"]=="3*b" and
      bridge45["exact_coordinate_bridge"]["exactness_forces"]=="b=0" and
      bridge45["exact_coordinate_bridge"]["affineness_forces"]=="b!=0" and
      bridge45["exact_coordinate_bridge"]["affinization_ring"]=="C[G,x*G^3,y*G^2]")
check("the nonexact coordinate-negative family is actual smooth affine geometry",
      all(bridge45["nonexact_control_family"][k] for k in
          ["actual_smooth_affine_models","all_a_c_complex","b_nonzero",
           "one_boundary_A1","every_source_coordinate_polar",
           "proper_genus_two_family_with_one_section_removed",
           "relative_ample_section_proof_written",
           "complete_parametrized_normal_form_with_source_choices"]) and
      bridge45["nonexact_control_family"]["global_volume_exact"] is False and
      bridge45["nonexact_control_family"]["finite_Keller_map_constructed"] is False)
check("the entire single-index-two-prime quartic packet is excluded without removing another prime",
      all(bridge45["Keller_consequences"][k] for k in
          ["one_index_two_prime_no_regular_coordinate_gap_closed",
           "complete_regular_coordinate_geometric_classification_all_degrees",
           "source_additive_action","all_degree_primitive_odd_marked_end",
           "all_degree_target_torus_and_prime_boundary_degree_exclusions",
           "all_quartic_single_index_two_prime_normalizations_excluded",
           "individual_e1_e2_primes_have_regular_coordinate_on_affine_open"]) and
      bridge45["Keller_consequences"]["source_additive_kernel"]=="C[y]" and
      bridge45["Keller_consequences"]["target_algebra_invariance_assumed"] is False and
      bridge45["Keller_consequences"]["common_coordinate_for_multiple_primes_claimed"] is False and
      bridge45["Keller_consequences"]["removing_other_primes_preserves_finiteness"] is False)
check("new bridge retains higher-index, multi-prime, higher-degree and assurance limits",
      all(bridge45["assurance"][k] for k in
          ["written_unbounded_canonical_and_global_proofs","general_signed_key_formula_imported",
           "primary_statements_and_complete_page_audits",
           "finite_controls_not_complete_geometric_proof","no_deep_coefficient_search"]) and
      all(bridge45["assurance"][k] is False for k in
          ["higher_index_normalizations_excluded","multiple_prime_quartic_packets_excluded",
           "higher_finite_Keller_degrees_excluded","finite_Keller_map_constructed",
           "formal_verification","independent_full_reproof","external_review",
           "novelty_priority","research_goal_completed"]))
check("next gate retains the finite algebra and distinguishes quartic from higher-degree indices",
      "PF2RULED45" in action["next_step"] and
      "two-prime quartic packet" in action["next_step"] and
      "QPK4 already excludes higher-index quartic cases" in action["next_step"] and
      "exact volume forces an actual regular source coordinate" in action["next_step"] and
      "higher finite field degrees" in action["next_step"])

packet46=json.loads((HERE/"unramified_orevkov_receipt.json").read_text())
check("new quartic reduction leaves only the two-prime zero-excess row",
      packet46["quartic_reduction"]["only_two_primes"] and
      packet46["quartic_reduction"]["indices"]==[2,1] and
      packet46["quartic_reduction"]["residue_degrees"]==[1,1] and
      packet46["quartic_reduction"]["excesses"]==[0,0] and
      packet46["quartic_reduction"]["two_prime_case_excluded"] is False)
check("quartic distinct images and double ramified value force non-surjectivity",
      all(packet46["quartic_reduction"][k] for k in
          ["smooth_N","distinct_immersive_images","ramified_image_has_double_value",
           "necessarily_non_surjective"]) and
      packet46["quartic_reduction"]["missing_pairs"]==[[0,4],[1,2],[2,0]] and
      packet46["quartic_reduction"]["singleton_pairs"]==[[0,3],[1,1]])
check("all-degree Orevkov quotient retains the full finite cover and unramified primes",
      all(packet46["all_degree_orevkov"][k] for k in
          ["full_finite_normalization_required","sole_ramified_prime",
           "simple_index_two","immersive_image_required","unramified_primes_retained",
           "positive_braid_infinity_excluded"]) and
      packet46["all_degree_orevkov"]["monodromy"]=="S_d" and
      packet46["all_degree_orevkov"]["quotient"]=="O(C) ->> S_d" and
      packet46["all_degree_orevkov"]["unramified_prime_count"]=="d-3")
check("nodal and tangential Nori scopes preserve the actual contact cost",
      packet46["all_degree_orevkov"]["nodal_bound"]=="a_n<=1" and
      packet46["all_degree_orevkov"]["quartic_tangency_bound"]=="a_n<=1+2*T" and
      packet46["all_degree_orevkov"]["T_definition"]=="sum intersection multiplicities m>=2, not number of tangencies" and
      packet46["assurance"]["explicit_tangency_blowup_cost"])
check("unramified affine classification covers every datum and boundary order",
      all(packet46["unramified_geometry"][k] for k in
          ["all_boundary_orders","regular_coordinate_automatic","gcd_condition",
           "all_admissible_data_actual_smooth_affine","all_polynomial_source_shears",
           "actual_exact_controls_all_r"]) and
      packet46["unramified_geometry"]["full_support"]=="r<=j<=q" and
      packet46["unramified_geometry"]["exactness_iff"]=="c_r=0" and
      packet46["unramified_geometry"]["boundary_LND_order"]=="r-1")
check("new reduction does not assume common coordinates, finite deletion or target invariance",
      all(packet46["unramified_geometry"][k] is False for k in
          ["target_algebra_invariance_assumed","removing_other_primes_preserves_finiteness",
           "common_coordinate_assumed"]))
check("Orevkov factorization retains cusp, multiple-ramification and positivity boundaries",
      packet46["all_degree_orevkov"]["cusp_factorization_claimed"] is False and
      packet46["all_degree_orevkov"]["multiple_ramified_prime_symmetric_factorization_claimed"] is False and
      packet46["assurance"]["all_polynomial_links_assumed_positive"] is False and
      packet46["assurance"]["withdrawn_Jacobian_reduction_used"] is False)
check("new packet reduction keeps full exclusion and assurance limits",
      all(packet46["assurance"][k] is False for k in
          ["full_quartic_exclusion","higher_degree_exclusion","formal_verification",
           "independent_full_reproof","external_review","novelty_priority","research_goal_completed"]))
check("all-index multi-prime zero-excess packets use the full branch union and actual nonabelian group",
      all(packet46["general_zero_excess_orevkov"][k] for k in
          ["all_field_degrees","all_generic_indices","any_number_ramified_primes",
           "all_ramified_residue_degrees_one_required","all_ramified_image_normalizations_immersive_required",
           "full_branch_union_required","G_nonabelian","abelian_transitive_group_regular",
           "Galois_obstruction_written","positive_braid_full_branch_union_excluded",
           "NN_criterion_on_all_simultaneous_infinity_paths"]) and
      packet46["general_zero_excess_orevkov"]["quotient"]=="O(B) ->> G" and
      all(packet46["general_zero_excess_orevkov"][k] is False for k in
          ["G_symmetric_assumed","separate_component_diagrams_suffice","positive_excess_packets_excluded"]))
check("current work and navigation adopt the actual sole remaining quartic packet",
      "PF2RULED46" in action["next_step"] and
      "S4 quotient of O(C)" in action["next_step"] and
      all("PF2RULED46" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

lefschetz47=json.loads((HERE/"lefschetz_centralizer_receipt.json").read_text())
check("componentwise quantitative WLT keeps all degrees and actual group indices",
      all(lefschetz47["quantitative_theorem"][k] for k in
          ["all_field_degrees","all_generic_indices","componentwise","actual_G_required",
           "selected_normalization_immersive_required","selected_support_zero_excess_required",
           "full_finite_normalization_retained","simultaneous_target_infinity_model_required",
           "no_moduli_connectivity_required","bare_image_tube_replaced_by_normalization_tube"]) and
      lefschetz47["quantitative_theorem"]["G_symmetric_assumed"] is False and
      lefschetz47["quantitative_theorem"]["conclusion"]=="(b-1)*eta<=2*delta" and
      lefschetz47["quantitative_theorem"]["positivity_required"]=="eta>0")
check("other ramified images may have cusp monodromy without invalidating the selected cyclic factor",
      lefschetz47["quantitative_theorem"]["all_other_supports_zero_excess_required"] is False and
      lefschetz47["quantitative_theorem"]["other_local_groups_may_be_nonabelian"] and
      lefschetz47["assurance"]["automatic_global_Orevkov_factorization"] is False)
check("the full normalization neighborhood image is checked against the meridian centralizer",
      lefschetz47["quantitative_theorem"]["tube_image"]=="rho(K) subset C_G(sigma)" and
      lefschetz47["quantitative_theorem"]["WLT_index"]=="[Pi:K]<=Cbar^2/eta" and
      lefschetz47["assurance"]["written_neighborhood_centralization_proof"] and
      lefschetz47["assurance"]["precise_Nori_WLT_A_C_import"])
check("symmetric all-degree and quartic singularity bounds retain positivity and the open packet",
      lefschetz47["symmetric_consequences"]["transposition_class"]=="binomial(d,2)" and
      lefschetz47["symmetric_consequences"]["field_degree_bound"]=="d*(d-1)<=2+4*delta/eta" and
      lefschetz47["symmetric_consequences"]["quartic_bound"]=="a_n<=1+min(2*T,floor(2*(k+T)/5))" and
      lefschetz47["symmetric_consequences"]["positive_normal_degree_required"] and
      lefschetz47["symmetric_consequences"]["quartic_case_excluded"] is False)
check("a single fold bounds field degree even before global monodromy is identified",
      all(lefschetz47["unidentified_group_fold_bound"][k] for k in
          ["any_transitive_G","single_transposition_required","prime_degree_forces_S_d",
           "eta_positive_required","other_inertia_unrestricted"]) and
      lefschetz47["unidentified_group_fold_bound"]["normal_transposition_block_size"]=="m>=2; m divides d" and
      lefschetz47["unidentified_group_fold_bound"]["exact_class"]=="d*(m-1)/2" and
      lefschetz47["unidentified_group_fold_bound"]["explicit_field_degree_bound"]=="d<=2+4*delta/eta" and
      lefschetz47["unidentified_group_fold_bound"]["intermediate_polynomial_plane_inferred"] is False)
check("new degree-six gate removes the old finite-incidence restrictions only for immersive quintic rows",
      lefschetz47["degree_six_application"]["oriented_rows"]==["A","C","N"] and
      lefschetz47["degree_six_application"]["second_image_degree_hypothesis_preserved"] and
      lefschetz47["degree_six_application"]["selected_second_image_immersive_required"] and
      all(lefschetz47["degree_six_application"][k] is False for k in
          ["selected_second_image_nodal_required","finite_singular_locus_avoidance_required",
           "finite_tangency_bound_required","nonimmersive_second_image_row_B_covered"]) and
      lefschetz47["degree_six_application"]["nonimmersed_first_image_allowed_in_N"])
check("new contact lower bound uses exact simultaneous infinity normal degree and conjugacy sizes",
      lefschetz47["degree_six_application"]["normal_degree"]=="7-(I-20)/2" and
      lefschetz47["degree_six_application"]["contact_lower_bound"]==34 and
      lefschetz47["degree_six_application"]["conjugacy_sizes"]==[15,40] and
      lefschetz47["degree_six_application"]["remaining_contacts"]==[34,38,40] and
      lefschetz47["degree_six_application"]["remaining_contacts_realized_as_Keller_images"] is False)
check("new quantitative gate preserves nonpositive-normal-degree and proof-assurance boundaries",
      all(lefschetz47["assurance"][k] is False for k in
          ["proof_from_bounded_search","eta_nonpositive_cases_excluded",
           "all_inertia_assumed_noncentral","full_quartic_exclusion",
           "unrestricted_second_image_degree_derived","formal_verification",
           "independent_full_reproof","external_review","novelty_priority","research_goal_completed"]))
r7_action=next(a for a in work if a["id"]=="WORK-PLANE-R7-20261113")
check("both related work gates adopt the new componentwise result without obsolete finite-tangency escape",
      "PF2RULED47" in action["claims"] and "PF2RULED47" in r7_action["claims"] and
      r7_action["reviewed_claims"]["PF2RULED47"]==claim_fingerprint(all_records["PF2RULED47"]) and
      "Finite tangencies cannot rescue contacts20 through32" in r7_action["next_step"] and
      "(b-1)*eta<=2*delta" in action["next_step"] and
      all("PF2RULED47" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

pair48=json.loads((HERE/"unramified_sheet_pair_receipt.json").read_text())
check("unramified normalization tube fixes a marked sheet rather than assuming all local inertia trivial",
      all(pair48["marked_sheet"][k] for k in
          ["all_field_degrees","full_finite_normalization_required","unramified_omitted_prime",
           "residue_one_immersive_image","shared_ramified_image_allowed",
           "other_local_groups_may_be_nonabelian","positive_normal_degree_required",
           "retained_unramified_coordinate_line_corollary"]) and
      pair48["marked_sheet"]["tube_image"]=="L subset G_a" and
      pair48["marked_sheet"]["bound"]=="(d-1)*eta_H<=2*delta_H" and
      pair48["marked_sheet"]["eta_nonpositive_excluded"] is False)
check("whole complementary simple-zero pair is excluded with the simultaneous coordinate constant retained",
      all(pair48["simple_pair"][k] for k in
          ["all_pole_orders","all_lower_coefficients_retained","simultaneous_constant_A_retained",
           "source_coordinates_complementary_required","R_regular_coordinate_simple",
           "D_regular_coordinate_simple","exact_residues_required","whole_sector_excluded",
           "other_primes_allowed","finite_map_not_required_for_curve_obstruction",
           "complete_curve_constructed"]) and
      pair48["simple_pair"]["indices"]==[2,1])
check("higher-order norm construction covers all generic indices under explicit smoothness and coordinate hypotheses",
      all(pair48["higher_pair"][k] for k in
          ["all_generic_indices","R_component_open_smooth_required","regular_coordinate_exists_required",
           "source_coordinates_complementary_required","D_regular_coordinate_simple_required",
           "full_norm_resultant","entire_top_column_repaired","complete_curve_constructed"]) and
      pair48["higher_pair"]["R_order"]=="r>=max(2,e)" and
      pair48["higher_pair"]["D_pole"]=="p>=r" and
      pair48["higher_pair"]["sufficient_defect"]=="r*q-S>=e")
check("index-two pairwise corollaries include full coprime and borderline sectors without excluding failed hypotheses",
      pair48["higher_pair"]["index_two_all_coprime_q_r_with_q_at_least_two_and_p_at_least_r_excluded"] and
      pair48["higher_pair"]["index_two_borderline_q_rminus1_with_r_at_least_three_and_p_at_least_r_excluded"] and
      all(pair48["higher_pair"][k] is False for k in
          ["p_less_than_r_excluded","index_two_defect_at_most_one_excluded","r2_q1_uniformly_excluded"]) and
      pair48["higher_pair"]["failed_defect_control"]["affine_realization_claimed"] is False)
check("pairwise geometric proof retains whole source charts and all infinity branches without target positivity",
      all(pair48["assurance"][k] for k in
          ["written_unbounded_marked_sheet_proof","written_unbounded_complete_curve_proof",
           "all_infinity_branches_accounted_for","actual_source_chart_embeddings_retained",
           "no_deep_coefficient_search","finite_controls_not_universal_proof"]))
check("new pairwise result keeps the original coordinate, deletion, full-exclusion and assurance boundaries",
      all(pair48["assurance"][k] is False for k in
          ["common_complementary_coordinates_automatic","deleting_other_primes_preserves_finiteness",
           "all_D_vanishing_orders_excluded","full_quartic_exclusion","higher_degree_exclusion",
           "formal_verification","independent_full_reproof","external_review",
           "novelty_priority","research_goal_completed"]))
check("current work and navigation adopt pairwise source geometry and the marked-sheet gate",
      "PF2RULED48" in action["next_step"] and
      "simple zero on the index-two R is impossible" in action["next_step"] and
      "rq-S<=1" in action["next_step"] and
      "(d-1)eta_H<=2delta_H" in action["next_step"] and
      all("PF2RULED48" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

norm49=json.loads((HERE/"two_sided_norm_receipt.json").read_text())
check("two full norms allow arbitrary unramified boundary order with explicit coordinate and smoothness hypotheses",
      all(norm49["two_sided"][k] for k in
          ["all_field_degrees","all_indices_on_R","other_primes_retained",
           "complementary_coordinates_required","smooth_affine_component_opens_required",
           "both_full_principal_parts_retained","exact_trace_residue_used",
           "prescribed_norm_band_cut","all_infinity_branches_accounted_for","complete_curve_constructed"]) and
      norm49["two_sided"]["D_simple_required"] is False and
      norm49["two_sided"]["target_positivity_required"] is False and
      norm49["two_sided"]["R_order"]=="r>=max(2,e)" and
      norm49["two_sided"]["R_pole"]=="q>=max(r+1,s)" and
      norm49["two_sided"]["D_pole"]=="p>=r" and
      norm49["two_sided"]["R_defect"]=="Delta_R>=e" and
      norm49["two_sided"]["D_defect"]=="Delta_D>=s+1")
check("simple index-two R now allows arbitrary unramified order without losing the simultaneous constant",
      norm49["simple_R"]["index"]==2 and norm49["simple_R"]["R_order"]==1 and
      norm49["simple_R"]["D_any_order"] and norm49["simple_R"]["constant_A_retained"] and
      norm49["simple_R"]["R_pole"]=="q>=max(2,s)" and
      norm49["simple_R"]["D_defect"]=="Delta_D>=s+1")
check("common-ruling controls are actual whole affine exact-volume geometry and do not supply finite Keller maps",
      all(norm49["common_ruling"][k] for k in
          ["unbounded_smooth_affine_exact_volume_controls","both_component_opens_affine",
           "whole_translation_criterion_required","global_explicit_primitive"]) and
      norm49["common_ruling"]["orders"]==[1,0] and
      norm49["common_ruling"]["finite_Keller_realization_claimed"] is False and
      all(row["torsor_multiindex"]==[1,0,1] for row in norm49["common_ruling"]["controls"]))
check("new norm implementation bypasses raw resultant sign while preserving older seals",
      norm49["norm_software"]["method"]=="quotient multiplication determinant" and
      norm49["norm_software"]["independent_algorithm_check"]=="reduced monic resultant" and
      norm49["norm_software"]["correct_root_product"]==875 and
      norm49["norm_software"]["older_sealed_PF48_inputs_modified"] is False and
      norm49["norm_software"]["older_raw_resultant_helper_safe_for_arbitrary_new_inputs"] is False)
check("two-sided written proof does not promote finite controls or close unsupported/full-goal cases",
      norm49["assurance"]["written_unbounded_proof"] and
      norm49["assurance"]["finite_controls_not_universal_proof"] and
      all(norm49["assurance"][k] is False for k in
          ["common_or_complementary_coordinates_automatic","insufficient_defect_excluded",
           "full_quartic_exclusion","full_higher_degree_exclusion","formal_verification",
           "independent_full_reproof","external_review","novelty_priority","research_goal_completed"]))
check("current work adopts higher D order and common-ruling finite-algebra gates",
      "PF2RULED49" in action["next_step"] and "Delta_D>=s+1" in action["next_step"] and
      "Common-ruling smooth affine exact-volume controls" in action["next_step"] and
      all("PF2RULED49" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

matrix50=json.loads((HERE/"coordinate_pencil_matrix_receipt.json").read_text())
check("Xie matrix uses primitive duals and gives only a necessary affine test",
      matrix50["matrix"]["primary"]=="Xie Theorem1.2" and
      all(matrix50["matrix"][k] for k in
          ["unnormalized_primitive_duals","all_other_primes_retained","necessary_for_affine_N",
           "all_principal_parts_retained","all_source_automorphisms_covered"]) and
      matrix50["matrix"]["passing_test_proves_affineness"] is False and
      matrix50["matrix"]["diagonal"]=="-(S+2q-r+e)" and
      matrix50["matrix"]["distinct_coordinate_pencil_cross"]=="m*q*p" and
      matrix50["matrix"]["strict_necessary_condition"]=="n_R*n_D>m^2*q^2*p^2")
check("all-degree global unramified pencil and common-complementary index-two dichotomy are recorded",
      all(matrix50["global_conclusions"][k] for k in
          ["all_Keller_field_degrees","all_unramified_primes_share_one_intrinsic_coordinate_pencil",
           "boundary_constants_need_not_agree","index_two_prime_pencil_common_or_complementary_to_unramified",
           "noncomplementary_pencils_excluded_without_pole_or_defect_restrictions"]) and
      matrix50["global_conclusions"]["unramified_bound"]=="0<n_D<=p^2" and
      matrix50["global_conclusions"]["index_two_bound"]=="0<n_R<=3q^2" and
      matrix50["global_conclusions"]["common_coordinate_case_excluded"] is False)
check("sharp complementary reduction removes every higher and multicharacteristic ramified principal part",
      all(matrix50["complementary"][k] for k in
          ["all_D_orders","exact_product_test_retained","simple_R_forced_q2",
           "higher_R_forced_q_rminus1","higher_R_full_principal_part_single_monomial",
           "all_higher_and_multicharacteristic_R_poles_excluded",
           "simple_R_requires_D_order_at_least_three","simple_D_R_order_at_least_three_forced_p2"]) and
      matrix50["complementary"]["simple_D_R_order_two_poles"]==[2,3,4,5])
check("matrix controls preserve common pencils and distinguish exactness and affineness",
      matrix50["controls_boundary"]["common_ruling_matrix"]==[[-5,-2],[-2,-4]] and
      all(matrix50["controls_boundary"][k] for k in
          ["nonexact_residue_control","passing_matrix_nonaffineness_control",
           "extractor_uses_exact_numeric_curvettes","shared_initial_centres_merged"]))
check("matrix theorem keeps imported proof and full-goal assurance boundaries",
      matrix50["assurance"]["written_unbounded_matrix_proof"] and
      matrix50["assurance"]["finite_controls_not_universal_proof"] and
      all(matrix50["assurance"][k] is False for k in
          ["imported_primary_theorem_independently_reproved","full_quartic_exclusion",
           "higher_degree_problem_excluded","target_Ga_action_or_ruling_inferred",
           "formal_verification","independent_full_reproof","external_review",
           "novelty_priority","research_goal_completed"]))
check("current work adopts the global pencil and minimal complementary gate",
      "PF2RULED50" in action["next_step"] and
      "The noncomplementary-coordinate gate is closed" in action["next_step"] and
      "q=r-1 and a single Laurent monomial" in action["next_step"] and
      all("PF2RULED50" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

models51=json.loads((HERE/"boundary_fibre_models_receipt.json").read_text())
budget51=models51["fibre_budget"]
check("full grouped matrix budget retains all field degrees, orders and primes",
      all(budget51[k] for k in ["all_field_degrees","full_matrix_required","all_D_orders_and_poles",
          "other_primes_retained","complementary_index_two_required","distinct_values_cross_zero"]) and
      budget51["strict_budget"]=="n_R/q^2>sum_lambda E_lambda" and
      budget51["group_energy_lower_bound"]=="E_lambda>=max_i(p_i^2/n_i)>=1")
check("one-value rule, exceptional two-value budget and conditional component counts are distinct",
      budget51["one_boundary_value_except_r2q1"] and budget51["r2q1_at_most_two_boundary_values"] and
      budget51["unconditional_prime_count_bound"] is False and
      budget51["conditional_all_simple_D_counts"]==dict(R_order3=6,R_order4=2,R_order_at_least5=1) and
      budget51["common_ruling_case_excluded"] is False)
family51=models51["affine_family"]
check("unbounded affine models have a global smooth presentation, marked source and exact volume",
      all(family51[k] for k in ["all_n_at_least_one","entire_affine_presentation",
          "two_smooth_principal_opens_cover","full_marked_A2_open","complementary_R_r2q1",
          "all_D_simple_pole2","full_weight_modules","source_Cstar_action"]) and
      family51["global_exact_primitive"]=="-2y*dx-x*dy" and
      family51["matrix_determinant"]=="2^n*(n+3)" and
      family51["conditional_zero_excess_degree"]=="n+3")
check("actual affine controls supply no finite Keller map or target action",
      family51["finite_Keller_map_supplied"] is False and family51["target_Cstar_action_inferred"] is False)
finite51=models51["finite_algebra"]
check("primitive monomial argument uses the full finite map and all-boundary verticality",
      all(finite51[k] for k in ["primitive_source_monomial_general_a_b",
          "global_regularity_and_all_boundary_vertical_required","full_finite_Keller_map_required",
          "generic_Cstar_image_residue_one","fixed_origin_deck_obstruction","primitive_element",
          "monic_norm_degree_is_full_Keller_degree","coincident_global_image_supports_allowed"]) and
      finite51["normal_monogenic_algebra_asserted"] is False and
      finite51["specializations"]==["chi(0)=H_x^2*H_y","chi(d_i)=J_i*H_i^2","chi(C)=G_R*J_C^f, f=1 or2"])
check("failed compact-curve shortcut retains every additional infinity end",
      all(models51["failures"][k] for k in ["resonance_C_equals_di_excluded",
          "extra_infinity_ends_retained","selected_chart_tests_not_all_ends"]) and
      all(row["additional_ends"]==row["r"]-2 and row["complete_curve_exclusion_valid"] is False
          for row in models51["failures"]["failed_shortcut_controls"]))
check("new grouped and affine results retain proof and full-goal assurance boundaries",
      models51["assurance"]["written_unbounded_proofs"] and
      models51["assurance"]["finite_controls_not_universal_proof"] and
      all(models51["assurance"][k] is False for k in ["formal_verification","independent_full_reproof",
          "external_review","novelty_priority","full_quartic_exclusion","higher_degree_exclusion",
          "research_goal_completed"]))
check("current work adopts the full matrix, actual unbounded models and genuine finite-map gate",
      "PF2RULED51" in action["next_step"] and "primitive monomial" in action["next_step"] and
      "no unconditional component-count bound" in action["next_step"] and
      all("PF2RULED51" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

symmetry52=json.loads((HERE/"source_symmetry_gate_receipt.json").read_text())
check("target-field intersection uses quasi-finite divisor domination",
      symmetry52["intersection"]["target_field_intersect_source_polynomials_equals_target_algebra"] and
      symmetry52["intersection"]["quasi_finite_divisor_domination_required"] and
      symmetry52["intersection"]["generic_finiteness_alone_suffices"] is False)
check("all-degree source symmetry descent is closed even at field level",
      all(symmetry52["descent"][k] for k in ["all_Keller_field_degrees",
          "all_polynomial_source_coordinate_systems","translation_derivation_field_stability_forces_automorphism",
          "nonzero_integral_diagonal_Euler_field_stability_forces_automorphism","all_weight_sign_patterns",
          "invariant_algebra_not_homogeneous_generators"]) and
      all(symmetry52["descent"][k] is False for k in ["target_action_linearization_required",
          "automatic_descent_from_finiteness","source_action_existence_excluded"]))
marked52=symmetry52["marked_models"]
check("explicit regular weight components escape the target field and retain primitive hypotheses",
      all(marked52[k] for k in ["all_n","all_other_primes_retained","both_components_outside_target_FIELD",
          "at_least_one_Euler_coordinate_derivative_outside_target_field",
          "primitive_components_require_no_intermediate_fields","retained_zero_excess_Sd_hypothesis_suffices"]) and
      marked52["every_target_polynomial_coordinate_max_weight_at_least"]==2 and
      marked52["every_target_polynomial_coordinate_min_weight_at_most"]==-4 and
      marked52["normal_monogenic_algebra_claimed"] is False)
check("finite trace and rational norm gate retains precise integrality and image-support limits",
      symmetry52["finite_gate"]["rational_norm_x"]=="H_x*product(H_i)/G_R" and
      symmetry52["finite_gate"]["rational_norm_y"]=="H_y*G_R^2/product(H_i)^2" and
      symmetry52["finite_gate"]["trace_dual"]=="Tr(x*B_n) subset A" and
      symmetry52["finite_gate"]["residue_degree_one_required_for_displayed_norms"] and
      symmetry52["finite_gate"]["global_image_support_cancellations_allowed"] and
      symmetry52["finite_gate"]["trace_x_squared_always_regular"] is False and
      symmetry52["finite_gate"]["trace_dual_implies_x_integral"] is False)
check("source-field symmetry proof retains all six failed-inference controls",
      len(symmetry52["negative_controls"])==6 and all(symmetry52["negative_controls"].values()))
check("source-field symmetry result retains imported proof and full-goal boundaries",
      symmetry52["assurance"]["written_unbounded_proofs"] and symmetry52["assurance"]["finite_controls_not_universal_proof"] and
      all(symmetry52["assurance"][k] is False for k in ["imported_one_line_theorem_independently_reproved",
          "finite_Keller_pair_supplied","formal_verification","independent_full_reproof","external_review",
          "novelty_priority","full_quartic_exclusion","higher_degree_exclusion","research_goal_completed"]))
check("current work rejects source-symmetry descent and adopts escaping primitive finite-algebra gates",
      "PF2RULED52" in action["next_step"] and "target FIELD" in action["next_step"] and
      "Tr(x*B_n) subset A" in action["next_step"] and
      all("PF2RULED52" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

pencil53=json.loads((HERE/"trace_dual_pencil_receipt.json").read_text())
check("trace-dual theorem covers the whole sole-ramified r2q1 row without source symmetry",
      all(pencil53["scope"][k] for k in ["all_field_degrees","whole_complementary_r2q1_row",
          "sole_index_two_ramified_prime_required","all_other_primes_unramified_and_retained",
          "arbitrary_D_orders_poles_values_and_coefficients"]) and
      pencil53["scope"]["source_symmetry_required"] is False and pencil53["scope"]["zero_excess_required"] is False)
check("entire trace dual and block surjection are exact, not generic-rank approximations",
      pencil53["trace_dual"]["equality"]=="B^vee=B+x*B" and
      pencil53["trace_dual"]["whole_R_nonvanishing_simple_pole_required"] and
      pencil53["trace_dual"]["quotient"]=="O_R(R)" and
      pencil53["trace_dual"]["target_basis_only_local_required"] and
      pencil53["trace_dual"]["generic_block_rank_suffices"] is False)
check("radical-rank formula recovers all local and global ramification excess",
      pencil53["excess"]["pointwise"]=="rank(G1|ker G0)=sum_(q in R over b)(ell_q-2)" and
      pencil53["excess"]["degree"]=="d=n+3+sum_b rank(G1|ker G0)" and
      all(pencil53["excess"][k] for k in ["all_nonimmersive_lengths","multiple_folds_not_nonimmersive",
          "zero_excess_iff_radical_totally_isotropic_everywhere","quartic_two_prime_requires_zero_rank_everywhere"]))
check("primitive inverse retains all finite-chart unramified primes and support coincidences",
      pencil53["inverse"]["every_source_coordinate_primitive_without_Sd"] and
      pencil53["inverse"]["leading"]=="G_R" and
      pencil53["inverse"]["degree_mod_generic_GR"]=="d-2" and
      pencil53["inverse"]["all_coefficient_ideal_is_unit"] and
      pencil53["inverse"]["coincident_global_supports_allowed"] and
      pencil53["inverse"]["finite_x_chart_normalization"]=="N minus R, including every D_i" and
      pencil53["inverse"]["normal_monogenic_algebra_asserted"] is False)
check("conductor compatibility keeps sealed old formulas and exact derivative pole",
      pencil53["conductor"]["retained_PF2OMIT4_formulas_not_claimed_new"] and
      pencil53["conductor"]["exact_R_derivative_pole"]=="d-3=n+E_R" and
      pencil53["conductor"]["generic_R_conductor"]==0 and
      pencil53["conductor"]["conductor_may_contain_D_i"] and
      pencil53["conductor"]["unregistered_draft_used"] is False)
check("trace pencil preserves non-Keller controls and every full-goal assurance limit",
      pencil53["assurance"]["written_unbounded_proofs"] and pencil53["assurance"]["finite_controls_not_universal_proof"] and
      pencil53["assurance"]["finite_non_Keller_controls_pass"] and
      all(pencil53["assurance"][k] is False for k in ["source_A2_from_controls","pencil_conditions_sufficient",
          "imported_duality_independently_reproved","formal_verification","independent_full_reproof",
          "external_review","novelty_priority","full_quartic_exclusion","higher_degree_exclusion","research_goal_completed"]))
check("current work adopts full dual generation and radical-rank finite-algebra gate",
      "PF2RULED53" in action["next_step"] and "B^vee=B+xB" in action["next_step"] and
      "radical isotropy everywhere" in action["next_step"] and
      all("PF2RULED53" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

colour54=json.loads((HERE/"quartic_colouring_receipt.json").read_text())
check("colouring exclusion covers every source order and target degree in the whole quartic packet",
      all(colour54["scope"][k] for k in ["actual_surviving_quartic_packet","all_source_pole_orders",
          "all_target_image_degrees","all_common_or_complementary_pencils","other_unramified_prime_retained",
          "actual_global_S4_required"]) and colour54["scope"]["determinant_gate_asserted_for_degrees_at_least_five"] is False)
check("actual infinity S4 quotient forces determinant divisible by3 through every cable stage",
      all(colour54["infinity"][k] for k in ["primary_epimorphism_without_general_position_at_infinity",
          "actual_S4_survives","nonconstant_Fox3_required","prime_to_three_knots_excluded","all_finite_cabling_lengths"]) and
      colour54["infinity"]["quotient"]=="S4/V4=S3" and
      colour54["infinity"]["necessary"]=="3 divides abs(Delta_K(-1))")
check("linear parabolic lift is stronger than three-colourability and retains global relations",
      all(colour54["lifting"][k] for k in ["quotient_by_coboundaries",
          "surjective_S4_lift_iff_nonzero_parabolic_H1","H1_S3_V_zero","linear_once_phi_fixed"]) and
      colour54["lifting"]["knot_infinity_solution_implies_global_solution"] is False and
      colour54["lifting"]["trefoil"]["parabolic_H1_dimension"]==0 and
      colour54["lifting"]["torus34"]["parabolic_H1_dimension"]==1)
check("sign triple is etale also at singularities and retains the Kummer units term",
      all(colour54["sign_surface"][k] for k in ["connected_etale_cyclic_three_cover_required",
          "etale_at_finite_singularities","finite_double_branch_local_image_C2_after_S3_quotient",
          "nonzero_antiinvariant_H1_et_mu3_required","unit_term_in_Kummer_sequence_retained",
          "Pic3_nonzero_requires_units_all_cubes"]) and
      colour54["sign_surface"]["cubic_subfield_of_quartic_asserted"] is False and
      colour54["sign_surface"]["simply_connected_from_smoothness"] is False)
check("all colouring and smooth-sign negative controls are retained",
      len(colour54["negative_controls"])==5 and all(colour54["negative_controls"].values()))
check("new quartic exclusion retains unbounded proof and all remaining assurance limits",
      colour54["assurance"]["written_unbounded_proofs"] and colour54["assurance"]["finite_controls_not_universal_proof"] and
      all(colour54["assurance"][k] is False for k in ["imported_primary_theorems_independently_reproved",
          "formal_verification","independent_full_reproof","external_review","novelty_priority",
          "finite_Keller_pair_supplied","full_quartic_exclusion","higher_degree_exclusion","research_goal_completed"]))
check("current work adopts determinant exclusion, exact lifts and etale sign triple",
      "PF2RULED54" in action["next_step"] and "determinant prime to3" in action["next_step"] and
      "anti-invariant H1_et(mu3)" in action["next_step"] and
      all("PF2RULED54" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

square55=json.loads((HERE/"sign_square_receipt.json").read_text())
check("square determinant gate covers the entire quartic with unrestricted source, target and ruling data",
      all(square55["scope"][k] for k in ["entire_actual_quartic_packet","all_source_orders",
          "all_target_degrees","all_rulings","both_omitted_primes_retained",
          "underlying_sign_surface_theorem_all_branch_degrees"]) and
      square55["scope"]["required_triple_supplied_for_higher_Keller_degrees"] is False)
check("irreducible prime-cover units close the Picard3 alternative even on the singular sign surface",
      all(square55["units"].values()) and square55["controls"]["Oka_Picard_invariants"]==[3] and
      square55["controls"]["cone_Cl_invariants"]==[3] and square55["controls"]["cone_Picard_torsion_empty"])
check("square formula retains singular affine freeness and the complete boundary sequence",
      all(square55["topology"][k] for k in ["singular_affine_Hamm_theorem_used",
          "exceptional_rational_trees","H1_resolution_equals_H1_sign","H2_resolution_torsion_free",
          "H3_resolution_zero","knot_double_boundary_rational_homology_sphere",
          "H1_boundary_surjects_H1_filling","H1_sign_finite"]) and
      square55["topology"]["negative_definiteness_required"] is False and
      square55["topology"]["rational_compactification_required"] is False)
check("whole-quartic gate requires9 and an isotropic3 subgroup through every cable stage",
      square55["determinant"]["identity"]=="D=abs(det Q)*abs(H1(S,Z))^2" and
      square55["determinant"]["necessary_quartic"]=="9 divides abs(Delta_K(-1))" and
      all(square55["determinant"][k] for k in ["Pic_torsion_order_squared_divides_D",
          "excludes_v3_zero_and_one","nonzero_isotropic_order3_subgroup_required","all_cable_stages"]) and
      square55["controls"]["linking_order9_anisotropic"]==[[0,0]])
check("T34 control passes PF54 lift but fails the new global sign extension gate",
      colour54["lifting"]["torus34"]["parabolic_H1_dimension"]==1 and
      square55["determinant"]["cable_controls"][0]["passes_old_three_gate"] and
      square55["determinant"]["cable_controls"][0]["passes_new_nine_gate"] is False)
check("square determinant proof retains all universal and independent assurance limits",
      square55["assurance"]["written_unbounded_proofs"] and square55["assurance"]["finite_controls_not_universal_proof"] and
      all(v is False for k,v in square55["assurance"].items() if k not in ["written_unbounded_proofs","finite_controls_not_universal_proof"]))
check("work and navigation adopt square determinant and actual Cartier Picard requirements",
      all(t in action["next_step"] for t in ["PF2RULED55","9 dividing","Pic(S)[3]","isotropic3"]) and
      all("PF2RULED55" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

power56=json.loads((HERE/"power_difference_receipt.json").read_text())
check("power-difference theorem retains all degrees, exponents and both coefficient image dimensions",
      all(power56["power_difference"][k] for k in ["all_polynomial_degrees","all_coprime_exponents_at_least_two",
          "irreducible_reduced_G_required","A1_normalization_required","dependent_image_case_proved",
          "dominant_image_case_proved","all_positive_valuation_rays_proved","retained_interior_exceptional_checked"]) and
      power56["power_difference"]["log_subadditivity_conjecture_used"] is False)
check("actual cubic quotient retains the whole packet and everywhere-etale singular local partition",
      all(power56["actual_quartic"][k] for k in ["entire_packet","all_source_orders","all_target_degrees",
          "all_rulings","both_omitted_primes_retained","cubic_in_S4_closure",
          "sign_triple_etale_at_singularities_required","trace_discriminant_reduced_G",
          "locally_monogenic_everywhere","globally_monogenic_excluded"]) and
      power56["actual_quartic"]["all_branch_fibres_partition"]==[2,1] and
      power56["actual_quartic"]["cubic_subfield_of_quartic_asserted"] is False)
check("essential cubic index is retained without forcing global monogenicity or distinct branch support",
      power56["index"]["primitive_integral_element_has_nonconstant_index"] and
      power56["index"]["discriminant_identity"]=="disc(F_t)=constant*G*I(t)^2" and
      power56["index"]["support_allowed_to_coincide_with_G"] and
      power56["index"]["free_module_basis_implies_algebra_generator"] is False and
      power56["index"]["monogenicity_forced_by_marked_source"] is False)
check("immersive A1 A5 and Cstar controls prevent dropping cover or end hypotheses",
      all(power56["controls"].values()))
check("power-difference obstruction retains all universal assurance limits",
      power56["assurance"]["written_unbounded_proofs"] and power56["assurance"]["finite_controls_not_universal_proof"] and
      all(v is False for k,v in power56["assurance"].items() if k not in ["written_unbounded_proofs","finite_controls_not_universal_proof"]))
check("work and navigation adopt the all-degree theorem and essential cubic index without full exclusion",
      all(t in action["next_step"] for t in ["PF2RULED56","coprime power-difference","nonconstant index I(t)","module freeness"]) and
      all("PF2RULED56" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

formats57=json.loads((HERE/"canonical_formats_receipt.json").read_text())
check("canonical and binary gates retain all degrees, ramification and residue profiles",
      all(formats57["scope"].values()) and formats57["canonical"]["nonzero_nontorsion"] and
      formats57["canonical"]["global_affine_complete_intersection_excluded"] and
      formats57["canonical"]["local_complete_intersections_excluded"] is False)
check("canonical-content format bound remains conditional and normal binary exclusion is all-degree",
      formats57["projective"]["s_divides_g"] and
      formats57["projective"]["conditional_degree_bound"]=="d<=2^(g+1)" and
      formats57["projective"]["primitive_formats"]==[[3],[2,2]] and
      formats57["projective"]["primitive_degree_at_least_five_CI_excluded"] and
      formats57["projective"]["all_higher_degree_algebras_excluded"] is False and
      formats57["projective"]["normal_binary_embedding_excluded_all_nontrivial_Keller_degrees"] and
      formats57["projective"]["maximal_index_remark_required"])
check("positive divisorial surface conductor is retained without misidentifying class coordinates or curve conductors",
      formats57["conductor"]["class_formula"]=="C_cond=(d-2)L-K_N" and
      all(formats57["conductor"][k] for k in ["all_admissible_sourcefinite_coordinates",
          "class_coefficients_strictly_positive_for_d_at_least_four","divisorial_conductor_nonzero","target_index_nonconstant"]) and
      all(formats57["conductor"][k] is False for k in ["curve_point_conductor_identified_with_surface_conductor",
          "class_vector_identified_with_boundary_multiplicities","codimension_two_questions_closed"]))
check("quartic binary and cubic share the same index with arbitrary r2q1 D pole and a distinct index curve",
      all(formats57["quartic"][k] for k in ["binary_and_cubic_normalization_indices_equal",
          "normalization_of_resolvent_order_is_actual_C3","whole_r2q1_D_pole_unrestricted",
          "generic_r2q1_index_coprime_both_omitted_images","distinct_index_curve_in_etale_source_locus_required",
          "actual_full_two_quadric_pencil_has_no_global_everywhere_smooth_member"]) and
      formats57["quartic"]["r2q1_conductor_class"]=="R+2*p*D" and
      formats57["quartic"]["distinct_support_means_disjoint_curves"] is False and
      formats57["quartic"]["global_matrix_model_constructed"] is False)
check("canonical format proof retains all full-exclusion and independent assurance limits",
      formats57["assurance"]["written_unbounded_proofs"] and formats57["assurance"]["finite_controls_not_universal_proof"] and
      all(v is False for k,v in formats57["assurance"].items() if k not in ["written_unbounded_proofs","finite_controls_not_universal_proof"]))
check("work and navigation adopt all-degree canonical gates and the full shared conductor",
      all(t in action["next_step"] for t in ["PF2RULED57","ALL nontrivial Keller field degrees",
          "SAME target index J","arbitrary D pole p","complete two-quadric pencil"]) and
      all("PF2RULED57" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

quadrics58=json.loads((HERE/"canonical_quadrics_receipt.json").read_text())
check("all-degree canonical construction keeps its Gorenstein hypothesis and proves the actual quadratic ideal",
      all(quadrics58["general"][k] for k in ["all_degrees_at_least_three","Gorenstein_hypothesis_explicit",
          "canonical_uses_dual_annihilator_of_one","canonical_embedding_proved","arithmetic_Gorensteinness_proved",
          "quadratic_generation_for_all_d_at_least_four","modules_free_by_Quillen_Suslin"]) and
      quadrics58["general"]["quadratic_rank"]=="d(d-3)/2" and
      quadrics58["general"]["all_normalizations_assumed_Gorenstein"] is False and
      quadrics58["general"]["non_Gorenstein_cases_excluded"] is False)
check("full global quartic pencil now exists in all rulings and orders without claiming exclusion",
      all(v for k,v in quadrics58["quartic"].items() if k not in ["full_fibre_types","finite_Keller_pair_constructed"]) and
      quadrics58["quartic"]["full_fibre_types"]==["1+1+1+1","2+1+1","2+2"] and
      quadrics58["quartic"]["finite_Keller_pair_constructed"] is False)
check("source residual retains full omitted sheets, normal cubic pullback and source-side descent limits",
      all(quadrics58["source"][k] for k in ["full_normalization_fibre_product_retained","diagonal_component_split",
          "residual_is_actual_normal_cubic_pullback","residual_to_N_etale_quasifinite",
          "source_side_projective_and_pencil_changes","universal_residual_identity","pullback_branch_smooth"]) and
      quadrics58["source"]["residual_smooth_finite_flat_degree"]==3 and
      all(quadrics58["source"][k] is False for k in ["ordinary_source_offdiagonal_automatically_finite",
          "changes_asserted_to_descend_to_target","pullback_branch_normalizations_assumed_A1",
          "residual_assumed_rational_or_plane_or_monogenic"]))
check("2026 primary correction is handled by a directly checked canonical fibre property",
      quadrics58["primary_correction"]["year"]==2026 and
      quadrics58["primary_correction"]["arbitrary_nondegenerate_embedding_not_sufficient"] and
      quadrics58["primary_correction"]["canonical_fibre_property_checked_directly"] and
      quadrics58["primary_correction"]["original_unrestricted_uniqueness_imported"] is False)
check("canonical quadrics retain every full-exclusion and independent assurance limit",
      quadrics58["assurance"]["written_unbounded_proofs"] and quadrics58["assurance"]["finite_controls_not_universal_proof"] and
      all(v is False for k,v in quadrics58["assurance"].items() if k not in ["written_unbounded_proofs","finite_controls_not_universal_proof"]))
check("work and navigation adopt unconditional full quartic pencil and complete residual",
      all(t in action["next_step"] for t in ["PF2RULED58","explicit Gorenstein hypothesis",
          "ENTIRE surviving quartic","full source residual","2026 corrected"]) and
      all("PF2RULED58" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

nodal59=json.loads((HERE/"nodal_deformation_receipt.json").read_text())
check("nodal deformation verifies the fixed completion, full end and LCP without global moduli connectedness",
      all(nodal59["deformation"][k] for k in ["immersed_normalizations_required","unobstructed_for_eta_at_least_one",
          "fixed_compactification_and_boundary","all_infinity_intersection_numbers_retained",
          "proper_end_isotopy_verified","normalizations_remain_A1","LCP_verified",
          "unrestricted_multiple_values_require_eta_at_least_two","eta1_two_point_use_requires_initial_double_values"]) and
      nodal59["deformation"]["global_moduli_connectedness_assumed"] is False and
      nodal59["deformation"]["eta_nonpositive_excluded"] is False)
check("only Orevkov invariant is transported and nodal Nori is applied after deformation",
      all(nodal59["topology"][k] for k in ["Orevkov_invariant_preserved","Nori_applied_only_after_nodal_deformation",
          "simultaneous_full_branch_union","whole_union_eta_at_least_two_implies_abelian",
          "whole_union_eta_at_least_one_and_individual_double_values_implies_abelian"]) and
      nodal59["topology"]["original_complement_group_preserved"] is False and
      nodal59["topology"]["ordinary_original_tangential_complement_asserted_abelian"] is False)
check("whole quartic positive eta range is excluded for every contact while all-degree hypotheses remain explicit",
      all(v for k,v in nodal59["Keller"].items() if k not in ["all_packets_assumed_zero_excess","unramified_image_H_treated_as_branch"]) and
      nodal59["Keller"]["all_packets_assumed_zero_excess"] is False and
      nodal59["Keller"]["unramified_image_H_treated_as_branch"] is False)
check("new normal-degree theorem retains predecessor scopes and its nonabelian original-complement control",
      all(nodal59["relation_to_retained"].values()))
check("nodal normal-degree proof retains every full-exclusion and independent assurance limit",
      nodal59["assurance"]["written_unbounded_proofs"] and nodal59["assurance"]["finite_controls_not_universal_proof"] and
      all(v is False for k,v in nodal59["assurance"].items() if k not in ["written_unbounded_proofs","finite_controls_not_universal_proof"]))
check("work and navigation focus the remaining quartic on eta nonpositive",
      all(t in action["next_step"] for t in ["PF2RULED59","ENTIRE positive-eta quartic range",
          "a_n<=1","FIXED target compactification","eta1 requires a triple value","Focus quartic work on eta<=0"]) and
      all("PF2RULED59" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

peripheral60=json.loads((HERE/"peripheral_zero_receipt.json").read_text())
check("whole quartic longitude is identity via the full paired residual cover with framing retained",
      all(peripheral60["peripheral"][k] for k in ["whole_surviving_quartic","all_source_rulings_poles_and_contacts",
          "full_residual_rank_two_cover_retained","finite_branch_parameters_paired","normalization_infinity_split",
          "preferred_longitude_identity","all_three_degrees_one","unramified_D_retained"]) and
      peripheral60["peripheral"]["preimage_branch_knot_components"]==3 and
      all(peripheral60["peripheral"][k] is False for k in ["meridian_framing_assumed_zero",
          "cuspidal_one_prime_model_substituted","longitude_asserted_in_arbitrary_higher_degrees"]))
check("longitude-zero linear gate keeps the distinction between group controls and affine or knot monodromy",
      peripheral60["linear"]["longitude_equation"]=="z(lambda)=0" and
      all(peripheral60["linear"][k] for k in ["coboundary_invariant","nonzero_constrained_H1_iff_S4_peripheral_lift",
          "all_finite_presentations_and_cable_lengths"]) and peripheral60["linear"]["Coxeter_dimensions"]==[1,0] and
      peripheral60["linear"]["T34_dimensions"][1]>0 and
      peripheral60["linear"]["control_Coxeter_is_knot_group"] is False and
      peripheral60["linear"]["passing_knot_lift_implies_affine_lift"] is False)
check("unbounded Wronskian normal-degree equality retains the actual primitive binomial hypothesis",
      all(v for k,v in peripheral60["binomial"].items() if k not in ["general_key_assumed_binomial",
          "all_two_stage_keys_excluded","target_field_symmetry_inferred"]) and
      all(peripheral60["binomial"][k] is False for k in ["general_key_assumed_binomial",
          "all_two_stage_keys_excluded","target_field_symmetry_inferred"]))
check("realized eta0 six-node curve passes the peripheral infinity gates without supplying a Keller cover",
      peripheral60["control"]["normal_degree"]==0 and peripheral60["control"]["canonical_a2"]==1 and
      peripheral60["control"]["nodes"]==6 and peripheral60["control"]["collision_parameters"]==12 and
      all(peripheral60["control"][k]==0 for k in ["triple_values","cusps","tangencies"]) and
      peripheral60["control"]["infinity_determinant"]==9 and
      all(peripheral60["control"][k] for k in ["finite_birational_A1_normalization","immersive",
          "infinity_S4_peripheral_representation","infinity_preferred_longitude_identity"]) and
      all(peripheral60["control"][k] is False for k in ["affine_complement_S4_representation_supplied",
          "sign_Picard_or_linking_gate_certified","finite_quartic_cover_supplied","marked_plane_source_supplied",
          "control_declared_Keller"]))
check("peripheral and eta0 proof retains every full-exclusion and independent assurance limit",
      peripheral60["assurance"]["written_unbounded_proofs"] and peripheral60["assurance"]["finite_controls_not_universal_proof"] and
      all(v is False for k,v in peripheral60["assurance"].items() if k not in ["written_unbounded_proofs","finite_controls_not_universal_proof"]))
check("work and navigation retain preferred longitude and the genuine eta0 control boundary",
      all(t in action["next_step"] for t in ["PF2RULED60","preferred longitude","z(lambda)=0",
          "primitive binomial key","eta=0","six-node eta0","affine complement"]) and
      all("PF2RULED60" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

end61=json.loads((HERE/"marked_end_lattice_receipt.json").read_text())
check("whole smooth affine source has its full integral end lattice and canonical mark in all degrees",
      all(end61["source_lattice"].values()) and end61["quartic"]["canonical_covector"]==[1,0] and
      all(v for k,v in end61["quartic"].items() if k!="canonical_covector"))
check("C34 comparison retains every peripheral representation and every meridian filling cycle",
      end61["control"]["3"]["total_representations"]==72 and
      len(end61["control"]["3"]["orbits"])==3 and
      all(x["orbit_size"]==24 and x["H1"]==x["quotient_by_ramification"] and
          x["H1"]["free"]==0 and x["H1"]["torsion"] in [[],[3]]
          for x in end61["control"]["3"]["orbits"]))
check("marked end gate excludes unbounded source sectors without claiming every source or target",
      all(v for k,v in end61["exclusions"].items() if k not in ["remaining_complementary_r_mod6",
          "all_C34_source_sectors_excluded","all_w8_representations_excluded"]) and
      end61["exclusions"]["remaining_complementary_r_mod6"]==[0,4] and
      end61["exclusions"]["all_C34_source_sectors_excluded"] is False and
      end61["exclusions"]["all_w8_representations_excluded"] is False)
check("source end gate is strictly stronger than the peripheral lift and keeps its assurance limits",
      end61["control"]["8"]["total_representations"]==168 and
      len(end61["control"]["8"]["orbits"])==7 and
      sum(x["H1"]==dict(free=1,torsion=[26]) for x in end61["control"]["8"]["orbits"])==3 and
      all(v is False for k,v in end61["assurance"].items() if k not in ["written_unbounded_topology",
          "finite_group_enumeration_exhaustive_for_two_fixed_knots","matrices_are_not_independent_topological_review"]))
check("work and navigation adopt full source lattice with target-specific source-sector boundaries",
      all(t in action["next_step"] for t in ["PF2RULED61","FULL source end lattice","canonical knot marking",
          "ALL peripheral S4 classes","minimal r2q1","r=0 or4 mod6"]) and
      all("PF2RULED61" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

orbit62=json.loads((HERE/"orbit_tree_receipt.json").read_text())
check("whole orbit-graph quotient covers every smooth degree and every finite cabling tree",
      all(v for k,v in orbit62["graph_theorem"].items() if k!="tree_sufficient_for_affine_cover") and
      orbit62["graph_theorem"]["tree_sufficient_for_affine_cover"] is False and
      all(v for k,v in orbit62["consequences"].items() if k!="only_dcycles_asserted_in_higher_degrees") and
      orbit62["consequences"]["only_dcycles_asserted_in_higher_degrees"] is False)
check("unbounded winding reduction comes from exact S4 exponent and keeps actual key hypotheses",
      all(orbit62["periodic"][k] for k in ["primitive_actual_binomial_trefoil_key_required",
          "quartic_eta_zero_forces_parameter_w_plus_one","period_from_exact_group_exponent","all_windings_at_least_two"]) and
      orbit62["periodic"]["period"]==12 and orbit62["periodic"]["combined_period"]==36 and
      orbit62["periodic"]["allowed_mod36"]==[3,8,9,15,21,26,27,33] and
      orbit62["periodic"]["all_passing_infinity_types_realized"] is False)
check("complete residue census retains residue2 and all tree survivors while rejecting cycles",
      all((orbit62["periodic"]["residue_certificates"][r]["conjugacy_classes"],
           orbit62["periodic"]["residue_certificates"][r]["tree_classes"])==expected
          for r,expected in {"2":(3,3),"3":(3,3),"8":(7,4),"9":(6,6)}.items()) and
      orbit62["periodic"]["w8_mod12_cycle_classes"]==3 and
      orbit62["periodic"]["even_w26_mod36_tree_classes"]==3)
check("orbit-tree result retains full-exclusion and independent assurance limits",
      orbit62["assurance"]["written_unbounded_proofs"] and
      orbit62["assurance"]["exhaustive_finite_group_certificate_after_period_reduction"] and
      all(v is False for k,v in orbit62["assurance"].items() if k not in ["written_unbounded_proofs",
          "exhaustive_finite_group_certificate_after_period_reduction"]))
check("work and navigation adopt all-degree whole graph and exact unbounded winding gates",
      all(t in action["next_step"] for t in ["PF2RULED62","WHOLE orbit graph tree","every cable length",
          "w=3 mod6 OR w=8 mod18","separate one-edge passes are insufficient","full source lattice/linking/canonical mark"]) and
      all("PF2RULED62" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

binomial63=json.loads((HERE/"binomial_periodic_receipt.json").read_text())
check("all primitive binomial parameters use exact exponent and consistent primary peripheral equations",
      all(v for k,v in binomial63["reduction"].items() if k not in ["general_keys_assumed_binomial","peripheral_period"]) and
      binomial63["reduction"]["general_keys_assumed_binomial"] is False and binomial63["reduction"]["peripheral_period"]==12)
check("complete three-parameter census retains all images and separate determinant period",
      all(binomial63["census"][k]==v for k,v in {"admissible_mod12_cells":1152,"peripheral_cells":244,
          "tree_cells":244,"representations":14304,"classes":596,"tree_classes":484,"cycle_classes":112,
          "combined_period":36,"admissible_mod36_cells":31104,"combined_surviving_cells":2772}.items()) and
      len(binomial63["census"]["complete_cells"])==1152 and binomial63["census"]["residue_counts_not_curve_counts"])
check("primary Zannier import supplies all eta0 immersions without inventing double-value or affine covers",
      all(binomial63["existence"][k] for k in ["Zannier_Theorem1_imported","all_coprime_p_q_and_r",
          "complex_coefficients","finite_birational_A1_normalization","immersive"]) and
      binomial63["existence"]["normal_degree"]==0 and
      all(binomial63["existence"][k] is False for k in ["finite_multiple_value_types_classified",
          "all_controls_nodal","quartic_double_value_bound_supplied"]) and
      all(binomial63["odd_family"][k] is False for k in ["sign_Picard_or_linking_supplied",
          "affine_complement_S4_extension_supplied","marked_plane_finite_cover_supplied"]))
check("infinite all-odd passing family retains primary and full-exclusion assurance limits",
      all(binomial63["odd_family"][k] for k in ["all_coprime_odd_p_q_at_least_three",
          "torus_S4_images_cyclic","all_r9mod36_explicit_peripheral_S4_and_tree",
          "polynomial_immersions_exist_for_all_these_parameters"]) and
      all(v is False for k,v in binomial63["assurance"].items() if k not in ["written_unbounded_reduction",
          "finite_census_exhaustive_after_period_proof","Clay_Watson_primary_peripherals_used"]))
check("work and navigation adopt all-parameter existence with the actual remaining affine obligation",
      all(t in action["next_step"] for t in ["PF2RULED63","ALL primitive actual-binomial p,q,r",
          "Zannier Theorem1","EVERY coprime p,q,r","r=9 mod36","quartic double-value condition",
          "genuine finite target algebra"]) and
      all("PF2RULED63" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

formal64=json.loads((HERE/"formal_branch_index_receipt.json").read_text())
check("odd dihedral torsor splitting retains singular gluing and all finite thickenings",
      all(formal64["theorem"][k] for k in ["all_odd_cyclic_degrees","dihedral_inversion_required",
          "etale_at_singularities_required","entire_ramification_curve_split","singular_gluing_retained",
          "all_finite_thickenings_split","compatible_formal_splitting","reflection_quotient_product",
          "branch_neighbourhood_monogenic"]) and
      all(formal64["theorem"][k] is False for k in ["normalization_A1_required_for_formal_theorem",
          "arbitrary_Sd_covers_covered","even_degree_covered"]))
check("actual cubic index avoids branch at every finite order but does not lift all orders globally",
      all(formal64["quartic"][k] for k in ["entire_packet","all_source_orders_and_rulings",
          "all_branch_degrees","every_finite_order_has_global_primitive_lift","index_disjoint_from_branch",
          "index_nonconstant_by_PF56","different_lift_per_order"]) and
      formal64["quartic"]["normalized_relation"]=="a^3-b^2=3*G*I^2" and
      formal64["quartic"]["one_global_lift_of_all_orders_exists"] is False)
check("formal branch result retains global torsor and pure power-difference boundaries",
      all(formal64["controls"][k] for k in ["Oka_nontrivial_Pic3_with_split_branch",
          "quasietale_A5_control_excluded","index_square_factor_retained"]) and
      all(formal64["controls"][k] is False for k in ["finite_jet_increase_supplies_exclusion",
          "Picard_restriction_injective_asserted","formal_split_implies_global_split"]))
check("formal splitting proof retains primary-import and full-exclusion assurance limits",
      all(formal64["assurance"][k] for k in ["written_unbounded_proofs",
          "primary_nilpotent_invariance_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in formal64["assurance"].items() if k not in ["written_unbounded_proofs",
          "primary_nilpotent_invariance_imported","finite_controls_not_universal_proof"]))
check("work and navigation require a global next gate after universal branch formal splitting",
      all(t in action["next_step"] for t in ["PF2RULED64","ALL odd dihedral",
          "ENTIRE ramification curve","DISJOINT from the branch","Different orders",
          "genuinely GLOBAL obstruction","shared normal cubic/quartic conductor"]) and
      all("PF2RULED64" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

cyclic65=json.loads((HERE/"cyclic_companion_receipt.json").read_text())
check("all-degree cyclic companion replacement retains full groups and primitive canonical mark",
      all(cyclic65["general"][k] for k in ["all_finite_cover_degrees","cyclic_companion_restrictions_required",
          "all_orbit_lengths","integral_direct_sum_replacement","preferred_longitude_zero_in_each_cyclic_piece",
          "outer_mark_preserved","all_smooth_source_degrees_with_disjoint_A1",
          "primitive_canonical_vector_required_for_integral_generator_bound"]) and
      all(cyclic65["general"][k] is False for k in ["arbitrary_noncyclic_pieces_replaced",
          "hypotheses_automatic_in_higher_degrees"]))
check("entire all-odd C3 case with3 dividing first pair is excluded by zero mark and noncyclic homology",
      all(cyclic65["all_odd"][k] for k in ["all_coprime_odd_p_q_at_least_three","all_windings",
          "actual_primitive_binomial_key_required","all_peripheral_S4_classes_covered","ramification_class_zero",
          "entire_C3_case_if_3_divides_pq_excluded","p3_q5_all_r32mod36_excluded"]) and
      cyclic65["all_odd"]["refined_necessary_gate"]=="r=9 mod36 OR (4|r AND9|k AND3 does not divide pq)" and
      cyclic65["all_odd"]["full_all_odd_class_excluded"] is False)
check("remaining all-odd ends retain exact order mark and unrestricted source-order consequences",
      cyclic65["all_odd"]["C4_group"]=="Z/(r/3)" and
      cyclic65["all_odd"]["C3_group"]=="Z/(k/3) direct-sum H1(Sigma3(T(p,q)))" and
      all(cyclic65["all_odd"][k] for k in ["remaining_D_odd_multiple_three","remaining_source_mark_zero",
          "minimal_source_r2q1_excluded_in_entire_target_class"]) and
      cyclic65["all_odd"]["complementary_source_orders_mod6"]==[0,4])
check("new end exclusion has exact symbolic proofs and retains affine and topological assurance limits",
      all(cyclic65["assurance"][k] for k in ["written_unbounded_proofs","exact_symbolic_exponent_reduction",
          "finite_controls_not_universal_proof","retained_primary_peripherals_imported",
          "Hillman_Sakuma_primary_integral_split_crosschecked"]) and
      all(v is False for k,v in cyclic65["assurance"].items() if k not in ["written_unbounded_proofs",
          "exact_symbolic_exponent_reduction","finite_controls_not_universal_proof",
          "retained_primary_peripherals_imported","Hillman_Sakuma_primary_integral_split_crosschecked"]))
check("work and navigation adopt the whole-family exclusion without claiming remaining affine compatibility",
      all(t in action["next_step"] for t in ["PF2RULED65","ALL-DEGREE cyclic companion",
          "ENTIRE C3 case with3|pq","ALL r32mod36","zero ramification class forces CYCLIC",
          "full linking form","not automatically source-coordinate resolvents"]) and
      all("PF2RULED65" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

oriented66=json.loads((HERE/"oriented_pattern_linking_receipt.json").read_text())
check("all-degree homology torus substitution retains orientation and marks without contact or homeomorphism claims",
      all(oriented66["general"][k] for k in ["all_finite_cover_degrees","integral_homology_solid_torus_required",
          "cyclic_companion_branched_H1_zero_required","oriented_linking_preserved","outer_mark_preserved"]) and
      all(oriented66["general"][k] is False for k in ["replacement_implies_homeomorphism",
          "replacement_implies_contact_equivalence","arbitrary_nontrivial_companion_homology_deleted"]))
check("explicit pattern orientation and unrestricted linking formula are fixed by complex and A-chain controls",
      all(oriented66["pattern"][k] for k in ["complex_outward_orientation_fixed",
          "A_chain_independent_orientation_control","unrestricted_Schur_and_rescaling_proved"]) and
      oriented66["pattern"]["model"]=="w^4+u^t*w+v^s=0" and
      oriented66["pattern"]["euler"]=="-1/(3st)" and oriented66["pattern"]["linking"]=="<-t/s>")
check("entire all-odd complementary source is excluded without pole or source-order limits",
      all(oriented66["all_odd"][k] for k in ["all_coprime_odd_first_pairs",
          "all_primitive_actual_binomial_windings","full_remaining_peripheral_classes_covered",
          "ramification_mark_zero","source_nR_congruent2_mod3","ENTIRE_complementary_source_excluded",
          "all_source_orders_and_unramified_poles","common_source_pencil_required",
          "common_source_boundary_value_coincides"]) and
      oriented66["all_odd"]["distinct_target_images_identified"] is False and
      oriented66["all_odd"]["common_source_configuration_excluded"] is False)
check("oriented reduction imports primary conventions but no filling classification or full exclusion",
      all(oriented66["assurance"][k] for k in ["written_unbounded_proofs",
          "primary_Euler_and_plumbing_conventions_imported","finite_controls_not_universal_proof",
          "orientation_sign_checked"]) and
      all(v is False for k,v in oriented66["assurance"].items() if k not in ["written_unbounded_proofs",
          "primary_Euler_and_plumbing_conventions_imported","finite_controls_not_universal_proof",
          "orientation_sign_checked"]))
check("work and navigation direct the all-odd programme to the coincident common sector",
      all(t in action["next_step"] for t in ["CURRENT PF2RULED66","ENTIRE complementary source",
          "COMMON intrinsic pencil","SAME boundary value","n_R=2 mod3",
          "genuine finite target algebra/shared conductor"]) and
      all("PF2RULED66" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

common67=json.loads((HERE/"common_pencil_mark_receipt.json").read_text())
check("all-degree common-pencil theorem keeps whole smooth disjoint boundary and coordinate map hypotheses",
      all(common67["source_theorem"][k] for k in ["all_field_degrees","any_unramified_count",
          "any_positive_ramified_volume_order","smooth_affine_extension_required",
          "entire_disjoint_A1_boundary_required","common_regular_source_coordinate_required",
          "complementary_x_global_P1_morphism_required","sole_positive_volume_prime_required",
          "other_volume_orders_zero_required","values_need_not_coincide"]))
check("unbounded leaf ancestry Gram and strict Schur force nonzero R self-linking",
      all(common67["source_theorem"][k] for k in ["whole_lattice_negative_Gram",
          "retained_components_free_leaves","unramified_leaves_incomparable",
          "no_unramified_leaf_descends_R","Schur_strictly_greater_than_one",
          "R_self_linking_strictly_between_zero_and_one","R_mark_nonzero"]) and
      common67["source_theorem"]["canonical_mark_nonzero_asserted_for_arbitrary_index"] is False)
check("entire all-odd actual-binomial class closes with no representation or source sector omitted",
      all(common67["all_odd"][k] for k in ["entire_actual_primitive_binomial_two_stage_class_excluded",
          "all_coprime_odd_p_q_at_least_three","all_windings_and_source_orders_poles",
          "all_peripheral_representation_classes_covered","complementary_closed_by_PF66",
          "common_closed_by_nonzero_R_mark","PF63_polynomial_immersions_still_exist"]) and
      all(common67["all_odd"][k] is False for k in ["polynomial_curves_declared_nonexistent",
          "general_keys_assumed_binomial","even_first_pairs_excluded","longer_sequences_excluded"]))
check("sharp Schur equality control explicitly violates the whole disjoint leaf hypothesis",
      common67["countercontrol"]["Gram"]==[[5,4],[4,4]] and
      common67["countercontrol"]["Schur"]==1 and common67["countercontrol"]["R_mark_zero"] and
      common67["countercontrol"]["leaf_degrees"]==[1,2] and
      common67["countercontrol"]["valid_separated_A1_extension"] is False)
check("new class closure retains primary factorization and full-goal assurance limits",
      all(common67["assurance"][k] for k in ["written_unbounded_proofs",
          "primary_surface_factorization_imported","finite_controls_not_universal_tree_census"]) and
      all(v is False for k,v in common67["assurance"].items() if k not in ["written_unbounded_proofs",
          "primary_surface_factorization_imported","finite_controls_not_universal_tree_census"]))
check("current work and navigation leave the closed all-odd class and preserve the full frontier",
      "CURRENT PF2RULED67" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE coprime-odd","TWO-STAGE",
          "NONZERO R self-linking","EVEN first pairs","longer sequences and higher degree",
          "PF63 polynomial immersions still exist"]) and
      all("PF2RULED67" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

parity68=json.loads((HERE/"cyclic_parity_receipt.json").read_text())
check("parity-free cyclic companion theorem preserves its actual-key and cyclic restriction",
      all(parity68["theorem"][k] for k in ["all_coprime_first_pairs",
          "all_windings_and_source_orders_poles","actual_primitive_binomial_two_stage_required",
          "cyclic_companion_restriction_required","cyclic_sector_only_C3_C4",
          "all_cyclic_sector_R_marks_zero","cyclic_cover3_formula_parity_free"]))
check("all C3 representations and specified C4 source sectors close across parity",
      all(parity68["theorem"][k] for k in ["C3_all_first_pair_parities_excluded",
          "C4_common_source_excluded","C4_even_exponent_divisible4_excluded"]) and
      parity68["theorem"]["C4_remaining_even_exponent_mod4"]==2 and
      parity68["theorem"]["C4_remaining_source"]=="complementary" and
      parity68["theorem"]["C4_remaining_coprimality"]=="gcd(odd first exponent,r/3)=1")
check("new cyclic closure preserves open noncyclic and C4 sectors and actual curves",
      all(parity68["theorem"][k] is False for k in ["C4_complementary_sector_closed",
          "noncyclic_companion_sectors_closed","all_even_first_pairs_excluded"]) and
      parity68["theorem"]["PF63_polynomial_immersions_still_exist"])
check("noncyclic nonzero mark countercontrol prevents an invalid universal zero-mark extension",
      parity68["noncyclic_nonzero_mark_control"]=={
          "H1":{"free":0,"torsion":[24]},"quotient":{"free":0,"torsion":[6]}})
check("new parity module and pattern proofs retain assurance limits",
      all(parity68["assurance"][k] for k in ["written_unbounded_proofs",
          "exact_symbolic_new_pattern","retained_primary_imports","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in parity68["assurance"].items() if k not in ["written_unbounded_proofs",
          "exact_symbolic_new_pattern","retained_primary_imports","finite_controls_not_universal_proof"]))
check("current parity-free work gate keeps representation and source distinctions",
      "CURRENT PF2RULED68" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL cyclic-order3","BOTH first-pair parities",
          "NONCYCLIC","NONZERO R class","COMPLEMENTARY","gcd(odd,r/3)=1",
          "other keys/longer sequences and higher degree"]) and
      all("PF2RULED68" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

genus69=json.loads((HERE/"seifert_base_genus_receipt.json").read_text())
check("all-degree base genus quotient retains explicit orientability smoothness and boundary assumptions",
      all(genus69["general"][k] for k in ["all_finite_cover_degrees","all_cabling_lengths",
          "orientable_selected_Seifert_bases_required","arbitrary_boundary_gluings_and_fillings",
          "cutting_tori_need_not_be_incompressible","minimal_JSJ_not_required",
          "smooth_actual_extension_and_disjoint_A1_required","all_lifted_Seifert_bases_must_be_planar"]))
check("base sheets retain nonidentity fibre monodromy and all-degree companion permutation tree",
      all(genus69["general"][k] for k in ["regular_fibre_orbits_are_base_sheets",
          "nontrivial_fibre_monodromy_retained",
          "arbitrary_degree_centreless_transitive_companion_and_transitive_cable_force_dessin_tree",
          "global_Sd_not_assumed_to_imply_piece_Sd"]) and
      genus69["general"]["companion_tree_identity"]=="cycles(a)+cycles(b)=d+1")
check("complete retained census has exactly32 new genus1 tree exclusions and452 genus zero trees",
      all(genus69["quartic"][k]==v for k,v in {
          "complete_period12_cells":1152,"periodic_classes":596,"old_tree_classes":484,
          "new_positive_genus_tree_exclusions":32,"remaining_genus_zero_tree_classes":452,
          "whole_residue_cells_closed":0}.items()) and
      genus69["quartic"]["S4_S4_generator_order3_order4_sector_excluded"] and
      genus69["quartic"]["all_S4_S4_representations_excluded"] is False and
      genus69["quartic"]["other_prior_gates_already_applied_to_452"] is False)
check("raw-sheet false-positive and genus-zero nonzero-mark controls are retained",
      genus69["nonidentity_fibre_finite_control"]=={
          "H1":{"free":0,"torsion":[3]},"quotient":{"free":0,"torsion":[3]}} and
      genus69["remaining_nonzero_mark_control"]=={
          "H1":{"free":0,"torsion":[24]},"quotient":{"free":0,"torsion":[6]}})
check("new general topology and source review retain all assurance limits",
      all(genus69["assurance"][k] for k in ["written_unbounded_Mayer_Vietoris_and_fibre_quotient_proofs",
          "complete_existing_census_audited","primary_Seifert_model_visually_read",
          "finite_controls_not_universal_topology_proof"]) and
      all(v is False for k,v in genus69["assurance"].items() if k not in [
          "written_unbounded_Mayer_Vietoris_and_fibre_quotient_proofs",
          "complete_existing_census_audited","primary_Seifert_model_visually_read",
          "finite_controls_not_universal_topology_proof"]))
check("current work uses general base planarity and retains genus-zero frontier",
      "CURRENT PF2RULED69" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL-DEGREE","FIBRE ORBITS","PLANAR",
          "32 further","452 genus-zero","NO whole residue cell","GENUS-ZERO noncyclic",
          "longer keys/higher degree","Full research goal is open"]) and
      all("PF2RULED69" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

cone70=json.loads((HERE/"peripheral_cone_budget_receipt.json").read_text())
check("all-degree peripheral cone quotient preserves hypotheses and is not a homology splitting",
      all(cone70["general"][k] for k in ["all_finite_field_degrees","all_cabling_lengths",
          "orientable_selected_Seifert_bases_required","arbitrary_boundary_gluings_and_fillings",
          "protected_group_is_quotient_not_direct_summand","all_filling_core_marks_killed",
          "smooth_disjoint_entire_A1_source_hypotheses_required","no_cyclic_companion_assumption",
          "actual_cone_orders_not_period12"]))
check("individual omitted-prime marks give b-s bound without canonical primitivity",
      all(cone70["general"][k] for k in ["individual_omitted_prime_marks_give_b_minus_s_bound",
          "individual_mark_bound_needs_no_canonical_primitivity",
          "all_actual_ramified_omitted_primes_are_filling_cores",
          "unramified_omitted_primes_not_automatically_filling_cores",
          "primitive_canonical_vector_required_for_canonical_only_all_prime_improvement",
          "canonical_only_nonprimitive_bad_primes_keep_b_bound"]))
check("quartic noncyclic reductions retain exact exponent and gcd boundaries",
      cone70["quartic"]["marked_end_quotient_cyclic"] and
      cone70["quartic"]["all_planar_noncyclic_cycle_classes_audited"]==272 and
      cone70["quartic"]["S4_S4_all_even_first_exponents_greater2_excluded"] and
      cone70["quartic"]["S4_S4_remaining_even_exponent"]==2 and
      cone70["quartic"]["S4_S4_remaining_gcd"]=="gcd(o/3,r/3)=1" and
      cone70["quartic"]["A4_S3_necessary_gcd"]=="gcd(e*o/6,r/2)=1" and
      cone70["quartic"]["S4_C2_necessary_gcd"]=="gcd(o/3,r)=1" and
      cone70["quartic"]["S3_S3_new_restriction"] is False and
      cone70["quartic"]["all_S4_S4_representations_excluded"] is False)
check("new marked cone argument retains source review and full-goal assurance limits",
      all(cone70["assurance"][k] for k in ["written_unbounded_quotient_prime_budget_and_cone_list_proofs",
          "primary_presentation_visually_read","finite_controls_not_universal_proof",
          "independent_spanning_tree_full_cover_predecessor_reused"]) and
      all(v is False for k,v in cone70["assurance"].items() if k not in [
          "written_unbounded_quotient_prime_budget_and_cone_list_proofs",
          "primary_presentation_visually_read","finite_controls_not_universal_proof",
          "independent_spanning_tree_full_cover_predecessor_reused"]))
control70=next(x for x in cone70["full_controls"] if (x["p"],x["q"],x["r"])==(2,3,9))
check("passing nonzero R mark remains distinct from protected quotient",
      control70["protected"]=={"free":0,"torsion":[6]} and
      control70["full"]=={"H1":{"free":0,"torsion":[24]},"quotient":{"free":0,"torsion":[6]}} and
      control70["excluded"] is False)
check("current gate uses actual individual marks and leaves passing sectors open",
      "CURRENT PF2RULED70" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL-DEGREE","OMITTED-PRIME KNOTS",
          "no canonical primitivity","EXACTLY e2","NOT period12","QUOTIENT",
          "NONZERO R","full goal stay OPEN"]) and
      all("PF2RULED70" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

fold71=json.loads((HERE/"fold_schur_s4_receipt.json").read_text())
check("strict fold bound is all-degree/count with sole order1 and whole common source hypotheses",
      all(fold71["source"][k] for k in ["all_finite_field_degrees","any_unramified_count",
          "sole_volume_order_exactly1_required","smooth_affine_and_entire_disjoint_A1_required",
          "common_regular_f_and_global_P1_x_required","retained_free_leaf_hypothesis_required",
          "no_positive_negative_boundary_edge","no_zero_or_new_free_b1_descends_positive",
          "b1_parent_squared_distance_strictly_greater1","fold_Schur_strictly_greater2",
          "R_self_linking_strictly_between0_and_half"]) and
      fold71["source"]["arbitrary_higher_volume_order_bound_asserted"] is False and
      fold71["source"]["complementary_source_bound_asserted"] is False)
check("integral S4 S4 peripheral proof determines full end group and exact R marking",
      all(fold71["target"][k] for k in ["actual_primitive_binomial_two_stage_quartic_S4_S4_required",
          "even_first_exponent2_required_by_PF70","s_t_positive_odd","filled_cable_m_zero",
          "R_has_no_odd_primary_component","two_primary_group_cyclic_order8"]) and
      fold71["target"]["end_H1"]=="Z/(8t) plus Z/s" and
      fold71["target"]["marked_quotient"]=="Z/(2t) plus Z/s" and
      fold71["target"]["R_exact_order"]==4 and fold71["target"]["R_self_linking"]=="1/2")
check("entire common S4 S4 sector closes while complementary remains open",
      fold71["target"]["entire_S4_S4_common_source_sector_excluded"] and
      fold71["target"]["all_first_exponents_windings_source_orders_poles"] and
      fold71["target"]["complementary_source_sector_excluded"] is False and
      fold71["target"]["remaining_complementary_gate"]=="even first exponent2 and gcd(s,t)=1")
check("new fold and peripheral proofs retain all realization and assurance limits",
      all(fold71["assurance"][k] for k in ["written_unbounded_log_order_ancestry_and_projection_proof",
          "written_integral_peripheral_homology_proof","retained_primary_audits_reused",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in fold71["assurance"].items() if k not in [
          "written_unbounded_log_order_ancestry_and_projection_proof",
          "written_integral_peripheral_homology_proof","retained_primary_audits_reused",
          "finite_controls_not_universal_proof"]))
check("current work leaves closed common S4 S4 and keeps actual-mark obligations",
      "CURRENT PF2RULED71" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL common S4/S4","ALL-DEGREE","Schur S>2",
          "EXACT order4","NO odd-primary part","Complementary S4/S4 remains OPEN",
          "ACTUAL R marks","other noncyclic common sectors","full goal remain OPEN"]) and
      all("PF2RULED71" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

fibre72=json.loads((HERE/"fibre_order_transfer_receipt.json").read_text())
check("general fibre order and unit transfer preserve all degree-independent hypotheses",
      all(fibre72["general"][k] for k in ["all_finite_cover_degrees",
          "all_cabling_lengths_when_hypotheses_hold","closed_orientable_Seifert_base_required",
          "nonzero_Euler_required","unit_transfer_two_cones_base_D_required",
          "cone_u_divides_v_required","actual_boundary_basis_and_pushforward_required",
          "unit_A_mod_u_required","integral_meridian_zero_proved",
          "peripheral_and_marked_torsion_coordinates_zero"]) and
      fibre72["general"]["regular_fibre_exact_order"]=="abs(Euler*lcm(cone orders))" and
      fibre72["general"]["arbitrary_noncylclic_pieces_replaced"] is False and
      fibre72["general"]["linking_form_splitting_asserted"] is False and
      fibre72["general"]["two_remaining_boundary_tori_covered_by_unit_lemma"] is False)
check("two complete noncyclic common sectors close with full groups and zero R marking",
      all(fibre72["quartic"][k] for k in ["actual_primitive_binomial_two_stage_required",
          "all_first_exponents_windings_source_orders_poles","all144_periodic_A4_S3_and_S4_C2_classes_covered",
          "A4_S3_R_zero","A4_S3_common_entirely_excluded","S4_C2_R_zero","S4_C2_common_entirely_excluded"]) and
      fibre72["quartic"]["A4_S3_end_H1"]=="Z/(eo/3) plus Z/(r/2)" and
      fibre72["quartic"]["S4_C2_end_H1"]=="Z/(o/3) plus Z/r" and
      fibre72["quartic"]["A4_S3_necessary_gcd"]=="gcd(eo/3,r/2)=1" and
      fibre72["quartic"]["A4_S3_necessary_r_mod4"]==2)
check("entire trefoil r8mod36 closes while common S3 and complementary passing rows remain open",
      fibre72["quartic"]["trefoil_all_r8mod36_excluded_in_both_source_pencils"] and
      fibre72["quartic"]["trefoil_remaining_even_winding_mod36"]==26 and
      fibre72["quartic"]["remaining_common_sector"]=="S3/S3 only" and
      fibre72["quartic"]["S3_S3_common_entirely_excluded"] is False and
      fibre72["quartic"]["passing_complementary_sectors_excluded"] is False)
check("S3 countercontrol retains nonzero R rather than an invalid transfer extension",
      any(x["r"]==15 and x["H1"]=={"free":0,"torsion":[3]} and
          x["quotient"]=={"free":0,"torsion":[]} for x in fibre72["S3_nonzero_mark_countercontrols"]))
check("fibre and transfer result preserves primary audits and complete research-goal limits",
      all(fibre72["assurance"][k] for k in ["written_unbounded_fibre_order_and_integral_transfer_proofs",
          "primary_Euler_naturality_complete_pages_rechecked","retained_Seifert_presentation_audit_reused",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in fibre72["assurance"].items() if k not in [
          "written_unbounded_fibre_order_and_integral_transfer_proofs",
          "primary_Euler_naturality_complete_pages_rechecked","retained_Seifert_presentation_audit_reused",
          "finite_controls_not_universal_proof"]))
check("current work narrows common frontier without assuming a two-boundary zero mark",
      "CURRENT PF2RULED72" in action["next_step"] and
      all(t in action["next_step"] for t in ["ONLY S3/S3","ACTUAL TWO-BOUNDARY",
          "NONZERO R","ENTIRE trefoil r8mod36","COMPLEMENTARY",
          "full goal remain OPEN","PF63 polynomial immersions still exist"]) and
      all("PF2RULED72" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

surface73=json.loads((HERE/"lifted_surface_s3_receipt.json").read_text())
check("all-degree lifted surfaces and linking replacement retain exact peripheral hypotheses",
      all(surface73["general"][k] for k in ["all_finite_cover_degrees",
          "oriented_meridian_disk_balance","oriented_Seifert_surface_balance",
          "torsion_linking_orthogonal_replacement","actual_primitive_boundary_maps_required",
          "boundary_torsion_coordinates_zero_required","interior_mark_torsion_retained",
          "peripheral_compatible_piece_epimorphism_induces_glued_group_quotient",
          "knot_exterior_longitude_filling_group_quotient"]) and
      surface73["general"]["individual_longitude_vanishing_asserted"] is False and
      surface73["general"]["group_equality_for_replacement_asserted"] is False)
check("full S3 marked form and common trefoil gate adopted without entire-sector closure",
      surface73["quartic"]["all64_S3_S3_classes"] and
      surface73["quartic"]["full_H1"]=="Z/D plus Z/b" and
      surface73["quartic"]["oriented_linking"]=="negative inverse(s)/D orthogonal negative t/b" and
      surface73["quartic"]["ramification_mark"]=="(sk*g,0)" and
      surface73["quartic"]["new_necessary_gate"]=="gcd(D,b)=1" and
      surface73["quartic"]["trefoil_common_necessary_r_mod48"]==15 and
      surface73["quartic"]["trefoil_common_necessary_r_min"]==63 and
      surface73["quartic"]["entire_S3_S3_common_sector_excluded"] is False)
check("unbounded linking-compatible paired family closes by full groups, not false homeomorphism",
      all(surface73["quartic"][k] for k in ["unbounded_full_marked_linking_matching_cluster_family",
          "same_explicit_family_excluded_by_full_end_group",
          "actual_target_group_surjects_virtual_filling_group"]) and
      surface73["quartic"]["virtual_filling_claimed_homeomorphic_to_actual_target"] is False and
      surface73["quartic"]["all_sources_for_those_target_keys_excluded"] is False and
      len(surface73["cluster_and_group_controls"])==3)
check("new end-group result preserves exact scope and open research goal",
      all(surface73["assurance"][k] for k in ["written_unrestricted_surface_linking_and_group_proofs",
          "prior_complete_primary_page_audits_reused","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in surface73["assurance"].items() if k not in [
          "written_unrestricted_surface_linking_and_group_proofs",
          "prior_complete_primary_page_audits_reused","finite_controls_not_universal_proof"]))
check("current frontier uses full groups after completed two-boundary linking comparison",
      "CURRENT PF2RULED73" in action["next_step"] and
      all(t in action["next_step"] for t in ["FULL marked linking","FINITE CYCLIC",
          "INFINITE target triangle quotient","PAIRED","OTHER common configurations",
          "GROUP SURJECTION","COMPLEMENTARY","full goal remain OPEN"]) and
      all("PF2RULED73" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

rational74=json.loads((HERE/"rational_end_orderability_receipt.json").read_text())
check("all-degree source rational decomposition keeps whole common-coordinate hypotheses",
      all(rational74["source_theorem"][k] for k in ["all_field_degrees",
          "all_ramification_volume_orders","any_omitted_line_count",
          "smooth_affine_extension_required","entire_disjoint_A1s_required",
          "common_regular_coordinate_finite_values_required",
          "x_global_P1_morphism_with_poles_required",
          "end_connected_sum_of_rational_negative_plumbing_links","no_extra_S1S2_summand",
          "unbounded_integer_chi_blowup_proof",
          "every_homomorphism_to_left_orderable_group_trivial"]) and
      rational74["source_theorem"]["source_end_assumed_irreducible"] is False and
      rational74["source_theorem"]["all_source_end_groups_finite_or_cyclic_asserted"] is False)
check("all S3 targets receive the actual Seifert quotient without false homeomorphism",
      all(rational74["quotient"][k] for k in ["all_actual_binomial_two_stage_quartic_S3_S3",
          "actual_to_Seifert_group_surjection","simultaneous_companion_abelianizations",
          "nontrivial_LO_quotient_excludes_common_source"]) and
      rational74["quotient"]["group_equality_or_homeomorphism_asserted"] is False and
      rational74["quotient"]["rationality_classified_for_all_first_pairs"] is False)
check("entire trefoil common closure retains complementary and target-key limits",
      all(rational74["trefoil"][k] for k in ["all_remaining_S3_S3_common_sources_excluded",
          "all_source_configurations_and_orders","all_windings_common_pencil_excluded",
          "all_common_representation_sectors_excluded","actual_target_has_nontrivial_LO_quotient"]) and
      rational74["trefoil"]["entire_target_keys_in_both_pencils_excluded"] is False and
      rational74["trefoil"]["complementary_passing_sources_excluded"] is False and
      rational74["trefoil"]["unbounded_nonrational_floor_witness"]=="l=s+1")
check("rational-end result explicitly imports deep topology and leaves full goal open",
      all(rational74["assurance"][k] for k in ["written_unbounded_source_and_target_proofs",
          "complete_primary_pages_visually_audited","deep_rational_link_orderability_theorems_imported",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in rational74["assurance"].items() if k not in [
          "written_unbounded_source_and_target_proofs","complete_primary_pages_visually_audited",
          "deep_rational_link_orderability_theorems_imported","finite_controls_not_universal_proof"]))
check("current work closes all trefoil common configurations and redirects structural comparison",
      "CURRENT PF2RULED74" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE TREFOIL COMMON","EVERY homomorphism",
          "LEFT-ORDERABLE","ACTUAL group SURJECTION","OTHER first pairs",
          "TRUE normalized Seifert data","COMPLEMENTARY","full goal remain OPEN"]) and
      all("PF2RULED74" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

sections75=json.loads((HERE/"seifert_sections_winding_three_receipt.json").read_text())
check("actual relative Seifert sections normalize even with shared cone primes",
      all(sections75["general"][k] for k in ["PF74_all_degree_source_obstruction_retained",
          "actual_peripheral_section_maps_required","closed_group_quotient_suffices"]) and
      sections75["general"]["all_higher_degree_packets_assumed_common"] is False and
      all(sections75["quartic"][k] for k in ["all64_S3_S3_incidence_types",
          "actual_annulus_relative_Euler_zero","root_and_boundary_sections_retained"]) and
      sections75["quartic"]["pairwise_coprime_cone_orders_required"] is False and
      sections75["quartic"]["twists_inferred_from_H1"] is False)
check("one-step nonrational quotient forces common winding exactly3 without closing residual",
      sections75["quartic"]["non_trefoil_s_at_least3_common_entirely_excluded"] and
      sections75["quartic"]["universal_Laufer_first_intersection"]==2 and
      sections75["quartic"]["universal_Pinkham_witness"]==1 and
      sections75["quartic"]["entire_actual_two_stage_common_requires_r"]==3 and
      sections75["quartic"]["trefoil_common_entirely_excluded"] and
      sections75["quartic"]["r3_non_trefoil_common_entirely_excluded"] is False and
      sections75["quartic"]["passing_complementary_sectors_excluded"] is False)
check("residual r3 retains ACTUAL pieces and primitive marked gluing",
      sections75["residual"]["virtual_quotient_group"]=="Z/D" and
      all(sections75["residual"][k] for k in ["actual_main_cable_homeomorphic_T2_times_I",
          "actual_companion_pieces_retained"]) and
      sections75["residual"]["actual_gluing_m"]==[-1,0] and
      sections75["residual"]["actual_gluing_L"]=="(k/2,1)" and
      sections75["residual"]["actual_R_mark"]=="(k,3)" and
      sections75["residual"]["exact_JSJ_classification_asserted"] is False and
      sections75["residual"]["affine_source_realization_asserted"] is False)
check("section winding theorem keeps imported topology and complete scope limits",
      all(sections75["assurance"][k] for k in ["written_unbounded_section_and_winding_proofs",
          "prior_complete_primary_page_audits_reused","imported_rational_link_orderability_results",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in sections75["assurance"].items() if k not in [
          "written_unbounded_section_and_winding_proofs","prior_complete_primary_page_audits_reused",
          "imported_rational_link_orderability_results","finite_controls_not_universal_proof"]))
check("current structural frontier is exact winding3 actual companion gluing",
      "CURRENT PF2RULED75" in action["next_step"] and
      all(t in action["next_step"] for t in ["rEXACTLY3","SHARED cone primes",
          "ACTUAL pieces","T2xI","primitive R slope(k,3)","TWO-SEIFERT",
          "COMPLEMENTARY","full goal remain OPEN"]) and
      all("PF2RULED75" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

matched76=json.loads((HERE/"matched_meridian_common_exclusion_receipt.json").read_text())
check("degree-free matched endpoint lemmas retain all actual graph/source hypotheses",
      all(matched76["general"][k] for k in ["no_field_degree_bound",
          "matched_actual_slopes_outside_both_L_space_interiors","Seifert_endpoint_lemma",
          "normalized_negative_Euler_required","at_least_two_genuine_cones_required",
          "central_plus_cone_count_at_least1_required","non_solid_torus_pieces_required",
          "actual_graph_gluing_required","irreducibility_verified_before_LO_implication",
          "PF74_explicit_common_source_hypotheses_required","nu2_handled_without_extending_Pinkham"]))
check("entire actual-binomial two-stage common sector closes at all parameters",
      all(matched76["quartic"][k] for k in ["actual_primitive_binomial_two_stage_required",
          "all_first_exponents_windings_source_orders_poles",
          "all_residual_r3_non_trefoil_common_sources_excluded","entire_common_source_sector_excluded",
          "companion_meridian_not_in_L_space_interior","ordinary_meridian_not_in_L_space_interior",
          "actual_meridians_matched","actual_r3_end_not_L_space",
          "actual_r3_end_group_nontrivial_left_orderable"]) and
      matched76["quartic"]["every_remaining_source_pencil"]=="COMPLEMENTARY only" and
      matched76["quartic"]["passing_complementary_sources_excluded"] is False and
      matched76["quartic"]["entire_target_keys_in_both_pencils_excluded"] is False)
check("matched-slope controls include rational endpoints and genuinely nonrational fillings",
      any(x["meridian_status"]=="endpoint" for x in matched76["exact_controls"]) and
      any(x["meridian_status"]=="outside_L_space_interval" for x in matched76["exact_controls"]) and
      any(x["critical_sum"]==1 and x["Laufer"]["rational"] is False for x in matched76["exact_controls"]))
check("matched closure keeps imported topology, actual companions and full goal limits",
      all(matched76["assurance"][k] for k in ["written_unbounded_endpoint_normalization_and_gluing_proofs",
          "complete_Rasmussen_and_Nemethi_primary_pages_visually_audited",
          "deep_interval_gluing_and_orderability_theorems_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in matched76["assurance"].items() if k not in [
          "written_unbounded_endpoint_normalization_and_gluing_proofs",
          "complete_Rasmussen_and_Nemethi_primary_pages_visually_audited",
          "deep_interval_gluing_and_orderability_theorems_imported","finite_controls_not_universal_proof"]))
check("current work redirects the remaining class to actual complementary source geometry",
      "CURRENT PF2RULED76" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE COMMON-source","ACTUAL-BINOMIAL TWO-STAGE QUARTIC",
          "BOTH OUTSIDE L-space INTERIORS","COMPLEMENTARY ONLY","ACTUAL COMPLEMENTARY",
          "full goal remain OPEN","PF63 polynomial immersions still exist"]) and
      all("PF2RULED76" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

corner77=json.loads((HERE/"complementary_corner_jsj_receipt.json").read_text())
check("degree-free complementary source theorem keeps exact two-line chart hypotheses",
      all(corner77["source"][k] for k in ["no_field_degree_bound",
          "two_entire_omitted_lines_required","one_fold_one_unramified_required",
          "complementary_actual_coordinate_pair_required","exact_source_volume_required",
          "PF50_ramified_pole_collapse_used","higher_order_R_side_actual_solid_torus",
          "arbitrary_D_characteristic_support_retained","terminal_Klein_exception_included"]) and
      corner77["source"]["arbitrary_omitted_count_endpoint_theorem_asserted"] is False and
      corner77["source"]["canonical_R_boundary_JSJ_adjacent_tori_at_most"]==1 and
      corner77["source"]["R_actual_unit_cable_slope"]=="mu+q*ell")
check("non-trefoil S3 target boundary middle closes both pencils with actual peripheral recognition",
      all(corner77["target"][k] for k in ["actual_binomial_two_stage_quartic_S3_S3_required",
          "all_first_exponents_windings_source_orders_poles","all64_incidence_types_checked",
          "actual_fixed_sheets_retained","actual_pieces_negative_orbifold_Euler",
          "nonparallel_edge_fibres_verified","full_peripheral_recognition_verified",
          "non_trefoil_S3_S3_complementary_entirely_excluded",
          "non_trefoil_S3_S3_sector_excluded_in_both_pencils"]) and
      corner77["target"]["non_trefoil_R_boundary_JSJ_adjacent_tori"]==2 and
      corner77["target"]["trefoil_S3_S3_complementary_entirely_excluded"] is False and
      corner77["target"]["other_representations_of_same_keys_excluded"] is False)
check("source corner controls retain arbitrary-support depth and exact actual cable framings",
      len(corner77["R_side_controls"])==5 and
      any(x["node_count"]>=4 for x in corner77["source_cluster_controls"]) and
      all(x["root_node_degree"]<=1 for x in corner77["source_cluster_controls"]) and
      all(x["R_slope_in_meridian_longitude"]==[1,x["q"]] for x in corner77["R_side_controls"]))
check("corner JSJ closure preserves deep imports and complete goal limits",
      all(corner77["assurance"][k] for k in ["written_unbounded_corner_spine_and_recognition_proofs",
          "complete_Hatcher_Waldhausen_primary_pages_rechecked","deep_rigidity_JSJ_theorems_imported",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in corner77["assurance"].items() if k not in [
          "written_unbounded_corner_spine_and_recognition_proofs",
          "complete_Hatcher_Waldhausen_primary_pages_rechecked","deep_rigidity_JSJ_theorems_imported",
          "finite_controls_not_universal_proof"]))
check("current work removes non-trefoil S3 while retaining solid-companion trefoil and other sectors",
      "CURRENT PF2RULED77" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE NON-TREFOIL S3/S3","BOTH source pencils",
          "TWO-line","KLEIN INCLUDED","AT MOST ONE","EXACTLY TWO","TREFOIL S3/S3",
          "MERIDIAN-FIBRE","full goal remain OPEN"]) and
      all("PF2RULED77" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

fusion78=json.loads((HERE/"trefoil_complementary_fusion_receipt.json").read_text())
check("unbounded complementary fibre reduction keeps two-line whole-source hypotheses",
      all(fusion78["source"][k] for k in ["no_degree_order_support_depth_bound",
          "two_entire_omitted_lines_required","full_smooth_exact_complementary_source_and_whole_charts_required",
          "internal_D_fibres_distinct_by_negative_chain_determinant","rho2_product_collar_removed",
          "Klein_alternative_allowed_and_excluded_by_terminal_piece_type"]) and
      fusion78["source"]["Klein_fibre_uniqueness_assumed"] is False and
      fusion78["source"]["fused_annulus_cones"]==["rho-1","rho"])
check("actual trefoil caps close every S3 representation in both pencils",
      all(fusion78["target"][k] for k in ["actual_primitive_binomial_two_stage_quartic_S3_S3_required",
          "s_ge3_requires_consecutive_cones_and_hence_s3","all_trefoil_complementary_windings_excluded",
          "entire_S3_S3_sector_excluded_in_both_pencils"]) and
      fusion78["target"]["other_representations_of_same_target_keys_excluded"] is False and
      fusion78["target"]["s1_R_meridian_fibre_distance"]==5 and
      fusion78["target"]["s1_source_marked_distance"]==3 and
      fusion78["target"]["s3_source_forced"]["H1_order"]==212 and
      fusion78["target"]["s3_target_end_H1_order"]==2)
check("fibre-fusion closure keeps actual recognition imports and all scope limits",
      all(fusion78["assurance"][k] for k in ["written_unbounded_Euclidean_bridge_Klein_and_cone_proofs",
          "full_peripheral_Waldhausen_recognition","s1_sufficient_largeness_checked_with_essential_annulus",
          "complete_primary_pages_rechecked","deep_topology_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in fusion78["assurance"].items() if k not in [
          "written_unbounded_Euclidean_bridge_Klein_and_cone_proofs","full_peripheral_Waldhausen_recognition",
          "s1_sufficient_largeness_checked_with_essential_annulus","complete_primary_pages_rechecked",
          "deep_topology_imported","finite_controls_not_universal_proof"]))
check("current frontier advances beyond all S3 S3 to other complementary sectors",
      "CURRENT PF2RULED78" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE S3/S3","BOTH source pencils","TWO-line",
          "CONSECUTIVE","KLEIN ALTERNATIVE","A4/S3,S4/S4,S4/C2","even cyclic4",
          "full goal remain OPEN"]) and
      all("PF2RULED78" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

terminal79=json.loads((HERE/"terminal_cone_binomial_closure_receipt.json").read_text())
check("terminal and marked-root criteria preserve all explicit source and target hypotheses",
      all(terminal79["general"][k] for k in ["no_degree_order_support_depth_bound",
          "two_entire_omitted_lines_required","full_smooth_exact_complementary_source_and_whole_charts_required",
          "terminal_opposite_R_piece_exactly_two_cones_if_JSJ_edge","negative_Euler_orientable_target_pieces_required",
          "root_two_cones_must_be_consecutive","single_piece_source_disk_two_cones_or_consecutive_three_cones",
          "Klein_survival_excluded_by_virtual_abelianness","Klein_alternate_merger_excluded_by_centre"]) and
      terminal79["general"]["Klein_other_alignment_assumed_absent"] is False and
      terminal79["general"]["smooth_D_diagonal"]=="2*p_D" and
      terminal79["general"]["rho2_one_cone_marked_distance"].startswith("abs(Q-3b)"))
check("entire actual binomial two-stage quartic target class closes all representations both pencils",
      all(terminal79["quartic"][k] for k in ["actual_primitive_binomial_two_stage_required",
          "PF46_50_actual_source_packet_hypotheses_required","all_coprime_first_pairs_ge2_all_windings_ge2_all_source_data",
          "S4_S4_complementary_entirely_excluded","A4_S3_complementary_entirely_excluded",
          "S4_C2_complementary_entirely_excluded","remaining_even_C4_complementary_entirely_excluded",
          "all_representations_and_both_pencils_in_stated_target_class_excluded",
          "entire_actual_binomial_two_stage_quartic_class_excluded","actual_fixed_sheet_caps_retained",
          "actual_solid_topology_proved_before_homological_kernel_slopes","both_solid_cap_signs_retained",
          "A4_fibre_parallel_cap_contradicts_boundary_injection"]) and
      sum(terminal79["quartic"]["complete596_class_exhaustion"].values())==596)
check("terminal class closure retains full recognition imports and general-key higher-degree limits",
      all(terminal79["assurance"][k] for k in ["written_unbounded_terminal_cap_kernel_and_contact_proofs",
          "complete_existing_census_exhausted","full_peripheral_recognition_hypotheses_checked",
          "compressible_target_excluded_without_Waldhausen","prior_complete_Hatcher_Waldhausen_page_audits_reused",
          "deep_JSJ_fibre_rigidity_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in terminal79["assurance"].items() if k not in [
          "written_unbounded_terminal_cap_kernel_and_contact_proofs","complete_existing_census_exhausted",
          "full_peripheral_recognition_hypotheses_checked","compressible_target_excluded_without_Waldhausen",
          "prior_complete_Hatcher_Waldhausen_page_audits_reused","deep_JSJ_fibre_rigidity_imported",
          "finite_controls_not_universal_proof"]))
check("current frontier advances beyond ALL actual binomial sectors to genuine lower-weight keys",
      "CURRENT PF2RULED79" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE ACTUAL PRIMITIVE-BINOMIAL TWO-STAGE QUARTIC TARGET CLASS",
          "BOTH source pencils","EVERY representation","TWO-line","ALL596",
          "GENERAL LOWER-WEIGHT/NON-BINOMIAL TWO-STAGE","Do NOT assume that identification",
          "full quartic and full goal remain OPEN"]) and
      all("PF2RULED79" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

general80=json.loads((HERE/"general_two_stage_receipt.json").read_text())
check("general eta0 transfer uses actual splice with lower keys and full preferred peripherals",
      all(general80["general"][k] for k in ["actual_reduced_positive_two_stage_splice_required",
          "genuine_first_and_last_pairs_required","lower_weight_key_terms_allowed",
          "eta0_all_general_two_stage_quartic_packets_excluded",
          "transfer_preserves_preferred_meridian_longitude_and_all_fillings","cable_Bezout_independence"]) and
      general80["general"]["literal_binomial_required"] is False and
      general80["general"]["global_polynomial_key_normalization_asserted"] is False and
      general80["general"]["general_peripheral_period"]==12)
check("odd product retains nonidentity fibre and closes both source pencils with full mark",
      all(general80["odd_sector"][k] for k in ["all_general_two_stage_odd_product_quartic_packets_excluded_in_both_pencils",
          "M_equals_mu_equals_z","L_identity","R_mark_zero","common_closed_by_PF67","complementary_closed_with_both_cap_signs"]) and
      general80["odd_sector"]["fibre_image"]=="transposition" and
      general80["odd_sector"]["generic_cable_fibre_degree"]==2 and
      general80["odd_sector"]["R_meridian_fibre_distance"]==2 and
      general80["odd_sector"]["filled_main_cable_cones"]==["r",2])
check("general complete census and negative-even frontier remain exact and fail closed",
      general80["census"]["admissible_cells"]==9216 and general80["census"]["conjugacy_classes"]==2752 and
      general80["census"]["odd_classes"]==256 and general80["census"]["remaining_even_planar_tree_types"]==1920 and
      general80["census"]["old1152_eta0_slice_exactly_recovered"] and
      general80["frontier"]["all_surviving_genuine_two_stage_quartic_packets_require_eta_negative_even"] and
      general80["frontier"]["admissible_negative_even_cells"]==6144 and
      general80["frontier"]["all_negative_even_two_stage_cases_excluded"] is False)
check("general splice theorem retains complete primary audit and full goal limitations",
      all(general80["assurance"][k] for k in ["written_unbounded_splice_transfer_period_and_odd_fibre_proofs",
          "complete_Neumann_Norbury_Clay_Watson_primary_pages_visually_audited",
          "deep_topology_imported_with_retained_PF79_hypotheses","full_peripheral_recognition_before_slope_comparison",
          "complete_finite_census_after_unbounded_period_proof","bounded_source_controls_not_universal_proof"]) and
      all(v is False for k,v in general80["assurance"].items() if k not in [
          "written_unbounded_splice_transfer_period_and_odd_fibre_proofs","complete_Neumann_Norbury_Clay_Watson_primary_pages_visually_audited",
          "deep_topology_imported_with_retained_PF79_hypotheses","full_peripheral_recognition_before_slope_comparison",
          "complete_finite_census_after_unbounded_period_proof","bounded_source_controls_not_universal_proof"]))
check("current work advances to general negative-even keys with independent cabling parameter",
      "CURRENT PF2RULED80" in action["next_step"] and
      all(t in action["next_step"] for t in ["GENERAL genuine TWO-STAGE","STRICTLY NEGATIVE EVEN",
          "ACTUAL MARKED SPLICE TRANSFER","TRANSPOSITION fibre","9216","1920",
          "NEGATIVE-EVEN","independent of rc+1","full quartic and full goal remain OPEN"]) and
      all("PF2RULED80" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

balance81=json.loads((HERE/"meridian_balance_two_stage_receipt.json").read_text())
check("degree-free nonunit meridian balance propagates actual missing interior slopes",
      all(balance81["general"][k] for k in ["no_field_degree_or_source_depth_bound",
          "nonunit_two_boundary_balance_allowed","balanced_slope_becomes_rational_longitude_after_other_filling",
          "nonsolid_boundary_incompressible_actual_graph_pieces_required",
          "actual_matched_slopes_outside_L_space_interiors_required",
          "two_step_Rasmussen_gluing_propagates_missing_interior",
          "actual_irreducibility_before_graph_LO_implication","PF74_common_source_no_LO_image_hypotheses_required"]))
check("all-k S3 actual sections and two explicit trefoil denominators close both pencils",
      all(balance81["S3"][k] for k in ["actual_all_k_section_presentation",
          "relative_section_j_zero_proved_from_integral_disk","non_trefoil_s1_matched_meridians",
          "non_trefoil_s_ge3_balanced_propagation","trefoil_n_outside_interval_interior",
          "full_common_sector_excluded","full_complementary_sector_excluded"]) and
      balance81["S3"]["trefoil_two_exact_denominators"]==["r","r-k"])
check("entire GENERAL genuine two-stage quartic class closes all arbitrary-k representations",
      all(balance81["quartic"][k] for k in ["general_genuine_two_stage_actual_PF46_50_packet_required",
          "arbitrary_lower_weight_key_terms_allowed","k_independent_of_rc_plus1",
          "all_coprime_first_pairs_windings_coefficients_and_source_data",
          "all_six_even_sectors_excluded_in_both_pencils","all1920_even_planar_types_audited",
          "general_eta0_and_odd_product_PF80_closures_retained",
          "entire_general_genuine_two_stage_target_class_excluded",
          "source_full_fusion_third_cone_exactly2_used","all_odd_C4_s1_closed_by_unbounded_divisor16_identity",
          "all_actual_fixed_sheets_caps_and_R_marks_retained","both_solid_cap_signs_retained"]))
check("meridian balance proof keeps imported deep topology and longer/full-goal limits",
      all(balance81["assurance"][k] for k in ["written_unbounded_balance_section_interval_and_sector_proofs",
          "complete_Rasmussen_Nemethi_primary_pages_visually_rechecked",
          "deep_L_space_gluing_graph_LO_JSJ_and_rigidity_imported","finite_controls_not_universal_proof",
          "interval_interiors_kept_distinct_from_endpoints"]) and
      all(v is False for k,v in balance81["assurance"].items() if k not in [
          "written_unbounded_balance_section_interval_and_sector_proofs","complete_Rasmussen_Nemethi_primary_pages_visually_rechecked",
          "deep_L_space_gluing_graph_LO_JSJ_and_rigidity_imported","finite_controls_not_universal_proof",
          "interval_interiors_kept_distinct_from_endpoints"]))
check("current frontier advances beyond every genuine two-stage quartic sector",
      "CURRENT PF2RULED81" in action["next_step"] and
      all(t in action["next_step"] for t in ["ENTIRE GENERAL GENUINE TWO-STAGE QUARTIC TARGET CLASS",
          "BOTH pencils","MERIDIAN BALANCE","INTERIOR","DIVISOR16","All1920",
          "THREE OR MORE genuine target","full quartic and full goal remain OPEN"]) and
      all("PF2RULED81" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

tail82=json.loads((HERE/"iterated_peripheral_tail_receipt.json").read_text())
check("all-length ten-state finite saturation retains the whole graph and nonidentity fibres",
      tail82["saturation"]["all_states"]==10 and tail82["saturation"]["edges"]==1596 and
      len(tail82["saturation"]["transition_profiles"])==29 and
      all(tail82["general"][k] for k in ["all_finite_degrees","all_finite_cabling_lengths",
          "exact_group_exponent_period","whole_component_graph_contraction",
          "nonidentity_fibre_quotients_retained","multi_boundary_nonunit_balance_lemma",
          "every_ordered_partial_filling_nonsolid_incompressible_required"]))
check("actual arbitrary absorbing tails preserve intrinsic common-source obstructions",
      all(tail82["quartic"].values()) and
      len(tail82["ordinary_full_cover_controls"])>=18)
check("tail theorem leaves first-entry complementary and full-goal obligations open",
      all(tail82["assurance"][k] for k in ["written_unbounded_state_period_and_tail_proofs",
          "complete_finite_saturation_after_period_proof","retained_primary_page_audits_reused",
          "actual_all_meridian_cycle_fillings_in_full_cover_controls","deep_topology_imported",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in tail82["assurance"].items() if k not in [
          "written_unbounded_state_period_and_tail_proofs","complete_finite_saturation_after_period_proof",
          "retained_primary_page_audits_reused","actual_all_meridian_cycle_fillings_in_full_cover_controls",
          "deep_topology_imported","finite_controls_not_universal_proof"]))
check("current frontier is pre-final intrinsic invariants and complementary rooted JSJ",
      "CURRENT PF2RULED82" in action["next_step"] and
      all(t in action["next_step"] for t in ["TEN STATES","ABSORBING","ODD","COMMON",
          "BEFORE FIRST ARRIVAL","NONTRIVIAL LONGITUDE","COMPLEMENTARY","full goal remain OPEN"]) and
      all("PF2RULED82" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

persistent83=json.loads((HERE/"persistent_jsj_terminal_modulus_receipt.json").read_text())
check("source-wide last-gcd modulus and zero-mark factor gate retain all actual hypotheses",
      all(persistent83["source"].values()) and
      len(persistent83["terminal_modulus_controls"])==24)
check("arbitrary negative-even simple prefix tails close BOTH pencils through actual persistent JSJ",
      all(persistent83["quartic"].values()) and
      len(persistent83["ordinary_full_cover_controls"])==16 and
      sum(persistent83["even_planar_sector_counts"].values())==1920)
check("persistent closure retains imported topology and open first-arrival ranges",
      all(persistent83["assurance"][k] for k in [
          "written_unbounded_gcd_contact_modulus_and_JSJ_proofs",
          "all_fixed_sheet_and_branch_cycle_fillings_retained",
          "complete_predecessor_primary_audits_reused","deep_JSJ_and_rigidity_imported",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in persistent83["assurance"].items() if k not in [
          "written_unbounded_gcd_contact_modulus_and_JSJ_proofs",
          "all_fixed_sheet_and_branch_cycle_fillings_retained",
          "complete_predecessor_primary_audits_reused","deep_JSJ_and_rigidity_imported",
          "finite_controls_not_universal_proof"]))
check("current frontier is other first-arrival states and prefix eta ranges",
      "CURRENT PF2RULED83" in action["next_step"] and
      all(t in action["next_step"] for t in ["BOTH PENCILS","NO D WINDING","DIVIDES10",
          "ALL-DEPTH TERMINAL MODULUS","FIRST ARRIVAL","NONTRIVIAL LONGITUDE",
          "prefix eta from final eta","full goal remain OPEN"]) and
      all("PF2RULED83" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

cyclic84=json.loads((HERE/"cyclic_loop_first_arrival_receipt.json").read_text())
check("all-degree cyclic fibre annuli and exact composed cellular paths retain actual hypotheses",
      all(cyclic84["general"].values()) and
      cyclic84["quartic_self_loop_choices"]==688 and
      cyclic84["noncyclic_cable_bound"]==2 and cyclic84["nonself_state_change_bound"]==3)
check("all-depth cyclic incoming classes close with explicitly retained residual terminals",
      all(cyclic84["quartic"].values()) and
      len(cyclic84["ordinary_deep_companion_controls"])==6 and
      cyclic84["larger_terminal_forced_candidates"]==[])
check("cyclic first-arrival proof keeps open cases and exact matrix rather than expanded-word extrapolation",
      all(cyclic84["assurance"][k] for k in [
          "written_unbounded_cyclic_geometry_and_first_arrival_proofs",
          "finite_DAG_certificate_after_exact_period_proof",
          "exact_compressed_cellular_paths_cross_checked_against_expanded_CW",
          "all_actual_meridian_cycle_fillings_retained","retained_primary_audits_reused",
          "deep_JSJ_and_peripheral_rigidity_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in cyclic84["assurance"].items() if k not in [
          "written_unbounded_cyclic_geometry_and_first_arrival_proofs",
          "finite_DAG_certificate_after_exact_period_proof",
          "exact_compressed_cellular_paths_cross_checked_against_expanded_CW",
          "all_actual_meridian_cycle_fillings_retained","retained_primary_audits_reused",
          "deep_JSJ_and_peripheral_rigidity_imported","finite_controls_not_universal_proof"]))
check("current frontier is noncyclic incoming states and restricted final cyclic4r3",
      "CURRENT PF2RULED84" in action["next_step"] and
      all(t in action["next_step"] for t in ["TWO NONCYCLIC","CYCLIC3","BOTH PENCILS",
          "NO PREFIX ETA","FINAL winding3","(2,25)","SMALLER","TERMINAL MODULUS",
          "EXACT COMPOSED CW PATH MATRICES","full goal remain OPEN"]) and
      all("PF2RULED84" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

twists85=json.loads((HERE/"terminal_cone_twists_receipt.json").read_text())
check("all-degree terminal normalized residues and correct fibre-cover direction retain actual hypotheses",
      all(twists85["general"].values()) and
      len(twists85["full_terminal_cluster_controls"])==5)
check("all cyclic incoming first-arrival quartic classes close at arbitrary depth and tails in BOTH pencils",
      all(twists85["quartic"].values()))
check("terminal twist closure keeps both signs primary audits and noncyclic obligations explicit",
      all(twists85["assurance"][k] for k in [
          "written_unbounded_local_slope_free_tail_and_congruence_proofs",
          "complete_Hatcher25_26_27_and_NR3_4_11_12_visually_inspected",
          "complete_primary_bytes_pinned","deep_JSJ_fibre_uniqueness_and_rigidity_imported",
          "both_signs_used_as_weaker_necessary_match","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in twists85["assurance"].items() if k not in [
          "written_unbounded_local_slope_free_tail_and_congruence_proofs",
          "complete_Hatcher25_26_27_and_NR3_4_11_12_visually_inspected",
          "complete_primary_bytes_pinned","deep_JSJ_fibre_uniqueness_and_rigidity_imported",
          "both_signs_used_as_weaker_necessary_match","finite_controls_not_universal_proof"]))
check("current frontier is noncyclic incoming cores and initial fullS4 trajectories",
      "CURRENT PF2RULED85" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL CYCLIC INCOMING","BOTH PENCILS",
          "NORMALIZED CONE TWIST","DIVIDES16OR8","NONCYCLIC","INITIAL STATE3",
          "full goal remain OPEN"]) and
      all("PF2RULED85" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

framing86=json.loads((HERE/"meridian_cycle_s4_receipt.json").read_text())
check("all-degree transitive-meridian framing retains actual surface finite-H1 and pattern hypotheses",
      all(framing86["general"].values()))
check("all S4 incoming4/5 quartic classes close in both pencils at every companion depth and tail length",
      all(framing86["quartic"].values()) and
      len(framing86["ordinary_full_cover_controls"])==16 and
      len(framing86["actual_branching_controls"])==4)
check("framing closure retains surface rather than group-order proof and precise open obligations",
      all(framing86["assurance"][k] for k in [
          "written_unbounded_surface_cap_group_and_JSJ_proofs",
          "full_unfilled_CW_balance_and_filled_mark_controls",
          "complete_retained_primary_page_audits_reused",
          "deep_JSJ_and_full_peripheral_rigidity_imported",
          "self_linking_computed_from_surface_not_CW_group_order",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in framing86["assurance"].items() if k not in [
          "written_unbounded_surface_cap_group_and_JSJ_proofs",
          "full_unfilled_CW_balance_and_filled_mark_controls",
          "complete_retained_primary_page_audits_reused",
          "deep_JSJ_and_full_peripheral_rigidity_imported",
          "self_linking_computed_from_surface_not_CW_group_order",
          "finite_controls_not_universal_proof"]))
check("current frontier retains S3 A4 fixed-sheet caps initialstate3 and full open goal",
      "CURRENT PF2RULED86" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL S4 INCOMING","BOTH PENCILS",
          "ALL-DEGREE MERIDIAN-CYCLE FRAMING","ONE-BOUNDARY UNBRANCHED",
          "incoming2(S3),7(A4)","ACTUAL FIXED-SHEET CAPS","INITIAL STATE3",
          "full goal remain OPEN"]) and
      all("PF2RULED86" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

pattern87=json.loads((HERE/"cyclic_pattern_a4_receipt.json").read_text())
check("all-degree arbitrary cyclic pattern transfer retains coprime winding integral marks and positive-rank modules",
      all(pattern87["general"].values()) and
      pattern87["knotted_core_pattern_controls"][-1]["closed"]["free"]==2)
check("all A4 common incoming close and complementary A4 reduces to initial trefoil at arbitrary lengths",
      all(pattern87["quartic"].values()) and
      len(pattern87["A4_loop_controls"])==25 and
      len(pattern87["full_target_controls"])==75 and
      len(pattern87["actual_A4_birth_branching"])==2)
check("cyclic-pattern A4 reduction retains written proofs imported topology and explicit solid-cap residual",
      all(pattern87["assurance"][k] for k in [
          "written_unbounded_surface_MV_mark_history_and_terminal_proofs",
          "full_CW_integral_mark_and_all_fill_controls","complete_retained_primary_audits_reused",
          "deep_JSJ_and_full_peripheral_rigidity_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in pattern87["assurance"].items() if k not in [
          "written_unbounded_surface_MV_mark_history_and_terminal_proofs",
          "full_CW_integral_mark_and_all_fill_controls","complete_retained_primary_audits_reused",
          "deep_JSJ_and_full_peripheral_rigidity_imported","finite_controls_not_universal_proof"]))
check("current frontier retains complementary A4 trefoil solid-cap histories S3 initialstate3 and full open goal",
      "CURRENT PF2RULED87" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL A4 INCOMING COMMON","ALL-DEGREE CYCLIC PATTERN MARK TRANSFER",
          "COMPLEMENTARY A4 LATER BIRTHS","ALL INITIAL NON-TREFOIL","TERMINAL MODULUS9",
          "ONLY INITIAL TREFOIL","SOLID FIXED-SHEET CAP","incoming2(S3)","INITIAL STATE3",
          "full goal remain OPEN"]) and
      all("PF2RULED87" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

solid88=json.loads((HERE/"unimodular_initial_s4_receipt.json").read_text())
check("unbounded unimodular solid-prefix and whole SINGLE-node source congruence retain actual hypotheses",
      all(solid88["general"].values()))
check("ALL initial simple-S4 torus trajectories close BOTH pencils at arbitrary tail length without prefixeta",
      all(solid88["quartic"].values()) and
      len(solid88["full_initial_trajectory_controls"])==36 and
      len(solid88["composed_solid_prefix_controls"])==4)
check("initial-S4 closure retains actual cap proofs deep primary imports and nonunimodular residual",
      all(solid88["assurance"][k] for k in [
          "written_unbounded_source_cap_first_nonsolid_and_terminal_proofs",
          "actual_full_filled_and_unfilled_CW_mark_controls","complete_retained_primary_audits_reused",
          "deep_JSJ_and_full_peripheral_rigidity_imported",
          "solid_target_excluded_by_group_without_Waldhausen","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in solid88["assurance"].items() if k not in [
          "written_unbounded_source_cap_first_nonsolid_and_terminal_proofs",
          "actual_full_filled_and_unfilled_CW_mark_controls","complete_retained_primary_audits_reused",
          "deep_JSJ_and_full_peripheral_rigidity_imported",
          "solid_target_excluded_by_group_without_Waldhausen","finite_controls_not_universal_proof"]))
check("current frontier is incoming2 and complementary A4 trefoil NONUNIMODULAR cap with no initialstate3 residual",
      "CURRENT PF2RULED88" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL INITIAL STATE3","BOTH PENCILS","NO PREFIX ETA",
          "UNIMODULAR SOLID-PREFIX","ARBITRARY INTEGER A","SINGLE-NODE",
          "PRIORITY ONLY incoming2(S3)","NONUNIMODULAR CAP A*m+2L","COMPOSED CW ONLY",
          "full goal remain OPEN"]) and
      all("PF2RULED88" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

cap89=json.loads((HERE/"nonunimodular_a4_receipt.json").read_text())
check("all-depth exact local last-gcd strip and coefficient2 caps retain actual source hypotheses",
      all(cap89["general"].values()))
check("all final-stage A4 arrivals close and later tails have complete ordered three-pair reduction",
      all(cap89["quartic"].values()) and cap89["ordered_terminal_pairs"]==[[2,4],[2,8],[4,10]] and
      len(cap89["full_target_controls"])==39 and len(cap89["actual_solid_C3_controls"])==3)
check("nonunimodular reduction retains deep imports local countercontrols and open later-tail obligations",
      all(cap89["assurance"][k] for k in [
          "written_unbounded_cap_strip_root_and_three_pair_proofs","actual_full_CW_and_exact_contact_controls",
          "complete_retained_primary_audits_reused","deep_JSJ_and_full_peripheral_rigidity_imported",
          "finite_divisors_follow_proved_unbounded_congruences",
          "local_countercontrols_not_affine_or_full_topological_matches","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in cap89["assurance"].items() if k not in [
          "written_unbounded_cap_strip_root_and_three_pair_proofs","actual_full_CW_and_exact_contact_controls",
          "complete_retained_primary_audits_reused","deep_JSJ_and_full_peripheral_rigidity_imported",
          "finite_divisors_follow_proved_unbounded_congruences",
          "local_countercontrols_not_affine_or_full_topological_matches","finite_controls_not_universal_proof"]))
check("current frontier retains three ordered A4 terminal types tail restrictions and incomingS3",
      "CURRENT PF2RULED89" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL FINAL-STAGE A4 ARRIVALS CLOSED","EXACT LAST-GCD STRIP",
          "(2,4),(2,8),(4,10)","ATLEAST TWO SIMPLE TAILS","LAST TAIL WINDING5MOD6",
          "PRIORITY fullorderededge/peripheral","incoming2(S3)","full goal remain OPEN"]) and
      all("PF2RULED89" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

ruling90=json.loads((HERE/"global_ruling_unfilled_cyclic_receipt.json").read_text())
check("whole source ruling unfilled cyclic periphery and integral-zero-core endpoint retain all-degree scope",
      all(ruling90["general"].values()))
check("both nonsolid S3 pencils close at arbitrary lengths leaving exact trefoil-square solid history",
      all(ruling90["quartic"].values()) and
      len(ruling90["full_target_controls"])==32 and
      len(ruling90["actual_solid_S3_controls"])==3 and
      len(ruling90["actual_S3_birth_controls"])==2 and
      ruling90["noncable_positive_rank_controls"][-1]["added"]["free"]==2)
check("new S3 reduction separates imported Floer inputs finite controls and still open solid histories",
      all(ruling90["assurance"][k] for k in [
          "written_unbounded_ruling_MV_endpoint_and_history_proofs",
          "actual_full_filled_and_unfilled_CW_controls","complete_relevant_primary_HTML_import_audits",
          "retained_complete_topology_audits_reused","new_deep_Floer_theorems_imported",
          "deep_JSJ_and_full_peripheral_rigidity_imported","finite_controls_not_universal_proof"]) and
      all(v is False for k,v in ruling90["assurance"].items() if k not in [
          "written_unbounded_ruling_MV_endpoint_and_history_proofs",
          "actual_full_filled_and_unfilled_CW_controls","complete_relevant_primary_HTML_import_audits",
          "retained_complete_topology_audits_reused","new_deep_Floer_theorems_imported",
          "deep_JSJ_and_full_peripheral_rigidity_imported","finite_controls_not_universal_proof"]) and
      ruling90["dihedral_trefoil"]["mutual_linking"]==2 and
      ruling90["dihedral_trefoil"]["full_algorithm_implemented"] is False and
      ruling90["dihedral_trefoil"]["paper_diagonal_framing_identified"] is False)
check("current frontier is solid trefoil-square S3 and PF89 laterA4 rather than all incomingS3",
      "CURRENT PF2RULED90" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL NONSOLID CAPPED COMPANION HISTORIES CLOSED",
          "GENERAL WHOLE-SOURCE COORDINATE RULING","GENERAL ALL-DEGREE UNFILLED CYCLIC PERIPHERY",
          "INTEGRALLY ZERO","SOLID TREFOIL-SQUARE HISTORY","ACTUAL ordinaryE(J)",
          "(2,4),(2,8),(4,10)","full goal remain OPEN"]) and
      all("PF2RULED90" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

genus91=json.loads((HERE/"genus_final_s3_receipt.json").read_text())
check("all-depth genus framing strict infinity and solid last-minus classification retain exact scope",
      all(genus91["general"].values()))
check("all final S3 arrivals close both pencils at arbitrary companion depth with explicit later-tail scope",
      all(genus91["quartic"].values()) and
      genus91["complete_unmerged_root_orders"]==[2,4] and
      genus91["complete_collar_root_pairs"]==[[4,1],[5,2],[8,3],[10,3]] and
      len(genus91["actual_solid_histories"])==6 and
      len(genus91["exact_interval_controls"])==36 and
      len(genus91["full_target_controls"])==108)
check("genus endpoint proof imports deep inputs and retains actual noninvariant tail linking countercontrol",
      all(genus91["assurance"][k] for k in [
          "written_unbounded_genus_interval_root_and_lattice_proofs",
          "actual_composed_full_CW_and_integral_mark_controls",
          "fresh_complete_canonical_PDF_pages_audited",
          "new_primary_Floer_cabling_and_rank_inputs_imported",
          "retained_deep_JSJ_gluing_LO_full_peripheral_inputs",
          "finite_controls_not_universal_proof",
          "endpoint_Lspace_seed_linking_not_assumed_tail_invariant"]) and
      all(v is False for k,v in genus91["assurance"].items() if k not in [
          "written_unbounded_genus_interval_root_and_lattice_proofs",
          "actual_composed_full_CW_and_integral_mark_controls",
          "fresh_complete_canonical_PDF_pages_audited",
          "new_primary_Floer_cabling_and_rank_inputs_imported",
          "retained_deep_JSJ_gluing_LO_full_peripheral_inputs",
          "finite_controls_not_universal_proof",
          "endpoint_Lspace_seed_linking_not_assumed_tail_invariant"]) and
      genus91["endpoint_tail_countercontrol"]["linking"]=="8/17" and
      genus91["endpoint_tail_countercontrol"]["not_affine_realization"] is True)
check("current frontier requires later S3 tails and retains user dihedral linking suggestion without false application",
      "CURRENT PF2RULED91" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL FINAL-STAGE S3 ARRIVALS CLOSED BOTH PENCILS",
          "GENERAL ALL-DEPTH GENUS/FRAMING","MINUS loop is LAST","MATCHED GENUS ENDPOINT",
          "ATLEAST ONE LATER SIMPLE TAIL","H1Z17 linking8/17",
          "USER-SUGGESTED fallback","doNOTcallfullS4/A4dihedral","full goal remain OPEN"]) and
      all("PF2RULED91" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

link92=json.loads((HERE/"recursive_linking_endpoints_receipt.json").read_text())
check("all-degree integral-null-core reconstruction D3 square recursion and endpoint seed lemma retain hypotheses",
      all(link92["general"].values()))
check("common endpoint tails close and last-minus arrival frontier is reduced without unsupported endpoint transport",
      all(link92["quartic"].values()) and len(link92["D3_composed_controls"])==7 and
      len(link92["actual_birth_controls"])==8 and len(link92["endpoint_seed_tail_controls"])==4 and
      len(link92["exact_HNF_reduction_controls"])==6 and
      link92["positive_rank_control"]["closed_H1"]["free"]==2)
check("independent author-code comparison retains coloring contrast and higher-degree preferred framing boundaries",
      len(link92["independent_author_comparison"]["examples"])==9 and
      len(link92["independent_author_comparison"]["higher_degree_controls"])==2 and
      link92["independent_author_comparison"]["different_coloring_control"]["mutual"]==2 and
      link92["independent_author_comparison"]["different_coloring_control"]["cyclic_loop_coloring_same_knot_mutual"]==18 and
      link92["independent_author_comparison"]["unmodified_author_code"]["commit"]=="bcb80912ad35ed0adb437aff02139dfd6371fbc4")
check("recursive linking closure retains deep inputs exact integer evidence and open positive-arrival obligations",
      all(link92["assurance"][k] for k in [
          "written_unbounded_peripheral_linking_endpoint_frontier_proofs",
          "all_unfilled_boundaries_composed_and_integral_torsion_checked",
          "external_author_code_unmodified_commit_pinned",
          "colored_inputs_selected_by_monodromy_not_linking_value",
          "deep_graph_gluing_LO_rigidity_and_Floer_inputs_retained",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in link92["assurance"].items() if k not in [
          "written_unbounded_peripheral_linking_endpoint_frontier_proofs",
          "all_unfilled_boundaries_composed_and_integral_torsion_checked",
          "external_author_code_unmodified_commit_pinned",
          "colored_inputs_selected_by_monodromy_not_linking_value",
          "deep_graph_gluing_LO_rigidity_and_Floer_inputs_retained",
          "finite_controls_not_universal_proof"]))
check("current frontier records closed endpoint seed all-length colored recursion and explicit surviving S3 types",
      "CURRENT PF2RULED92" in action["next_step"] and
      all(t in action["next_step"] for t in ["ALL COMMON S3 ARRIVAL ETA<=0 CLOSED",
          "GENERAL ALL-DEGREE UNFILLED LINKING RECONSTRUCTION","externalcomparison COMPLETED",
          "ALL LAST-MINUS S3 survivors BOTH pencils REQUIRE ARRIVALr3",
          "0<=t-c<=(w-3)/2","t=c,D>1extraannulus remainsOPEN",
          "EXACTk2Ar+/-2","HNF exactrowlatticebasis BEFORE Smith","full goal remain OPEN"]) and
      all("PF2RULED92" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

grid93=json.loads((HERE/"integer_grid_common_receipt.json").read_text())
check("integer grid and finite graph denominator certificates retain all-degree hypotheses",
      all(grid93["general"].values()))
check("all common genuine iterated quartic class closes through actual normal section and exhaustive incoming graph",
      all(grid93["quartic"].values()) and
      len(grid93["signed_root_certificates"])==70 and
      len(grid93["cyclic_grid_controls"])==7 and
      len(grid93["product_peripheral_controls"])==4 and
      len(grid93["actual_product_seed_tail_controls"])==12 and
      set(grid93["complete_first_arrival_exhaustion"])=={"2","4","5","6","7","9"})
check("new common closure imports graph interval theorem and preserves complementary full goal assurance",
      all(grid93["assurance"][k] for k in [
          "written_unbounded_grid_section_witness_and_exhaustion_proofs",
          "new_primary_graph_interval_theorem_import_explicit",
          "integral_boundary_and_composed_actual_cover_controls",
          "retained_deep_Floer_graph_LO_and_full_peripheral_hypotheses",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in grid93["assurance"].items() if k not in [
          "written_unbounded_grid_section_witness_and_exhaustion_proofs",
          "new_primary_graph_interval_theorem_import_explicit",
          "integral_boundary_and_composed_actual_cover_controls",
          "retained_deep_Floer_graph_LO_and_full_peripheral_hypotheses",
          "finite_controls_not_universal_proof"]))
check("current frontier closes all common class with retained hypotheses and targets complementary S3 A4 only",
      "CURRENT PF2RULED93" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "ENTIRE COMMON GENUINE POSITIVE ITERATED-TORUS QUARTIC CLASS CLOSED AT ALL FINITE LENGTHS",
          "SIMPLEMERIDIAN FULLFINAL-S4","GENERAL DEGREE-FREE INTEGER GRID",
          "NEWRasmussenThm4.6","PRIORITY COMPLEMENTARY S3",
          "EXACTk2Ar+/-2","(2,4),(2,8),(4,10)",
          "KIG28fullgenericA4doesNOTcloseintermediateA4","full goal remain OPEN"]) and
      all("PF2RULED93" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

cap94=json.loads((HERE/"cyclic_solid_cap_genus_receipt.json").read_text())
check("all-degree cyclic solid cap genus and marked prime laws retain exact hypotheses",
      all(cap94["general"].values()))
check("complementary A4 allplus reset and lastminus5 exception with S3 exact3 primary reduction retain scope",
      all(cap94["quartic"].values()) and
      len(cap94["all_degree_cyclic_cap_controls"])==11 and
      len(cap94["actual_A4_solid_histories"])==9 and
      len(cap94["full_actual_target_controls"])==28 and
      len(cap94["prime_mark_lattice_controls"])==7 and
      len(cap94["S3_exact_3primary_controls"])==38 and
      all(x["labelled_pairs"]==24 for x in cap94["exhaustive_initial_A4_orbits"]))
check("cap prime reductions preserve noncyclic cancellation and complementary full goal boundaries",
      cap94["prime_mark_lattice_controls"][-1]["equal_valuation_case"] and
      cap94["prime_mark_lattice_controls"][-1]["valuation"]==3 and
      all(cap94["assurance"][k] for k in [
          "written_unbounded_cap_genus_infinity_and_frontier_proofs",
          "finite_initial_A4_coloring_classification_EXHAUSTIVE",
          "exact_integral_peripherals_and_composed_cover_controls",
          "retained_deep_topology_and_primary_Floer_inputs_imported",
          "finite_history_controls_not_universal_proof"]) and
      all(v is False for k,v in cap94["assurance"].items() if k not in [
          "written_unbounded_cap_genus_infinity_and_frontier_proofs",
          "finite_initial_A4_coloring_classification_EXHAUSTIVE",
          "exact_integral_peripherals_and_composed_cover_controls",
          "retained_deep_topology_and_primary_Floer_inputs_imported",
          "finite_history_controls_not_universal_proof"]))
check("current frontier uses square caps allplus reset exception and exact cyclic marked factor3",
      "CURRENT PF2RULED94" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "COMPLEMENTARY S3/A4 FRONTIER REDUCED BY TWO GENERAL LAWS",
          "ALL-DEGREE actualcyclicunitcap","MINUS ONLY LAST",
          "MUST ALLPLUS","ONLY LAST-MINUS WINDING5",
          "EXACT SINGLE FACTOR3 IN END ORDER",
          "EVERY SOLID-LOOP AND LATER-TAIL WINDING PRIME TO3",
          "EQUALcasecancellationnotclaimed","PF93ENTIRECOMMON",
          "full goal remain OPEN"]) and
      all("PF2RULED94" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

bridge95=json.loads((HERE/"terminal_bridge_closure_receipt.json").read_text())
check("general terminal bridge cone and signed divisor law retain local and zero-constant boundaries",
      all(bridge95["general"].values()))
check("ALL S3 and stillsolid A4 close both pencils while earlierreset A4 remains explicit",
      all(bridge95["quartic"].values()) and
      bridge95["divisor_constants"]=={"still_solid":[-8,-12],"reset_2_8":[-6,-16],"reset_4_10":[-2,-12]} and
      bridge95["postreset_first_winding_lists"]==[[2,4,8,16],[2,4]] and
      len(bridge95["complete_integer_source_graphs"])==4 and
      len(bridge95["symbolic_old_helper_full_row_controls"])==3 and
      len(bridge95["full_actual_target_controls"])==6)
check("terminal bridge proof retains primary convention deep topology and open reset assurance",
      all(bridge95["assurance"][k] for k in [
          "written_unbounded_terminal_slope_and_divisor_proofs",
          "actual_oriented_slopes_with_independent_signs_allowed",
          "exact_full_graph_and_composed_target_controls",
          "primary_characteristic_convention_and_retained_deep_topology_explicit",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in bridge95["assurance"].items() if k not in [
          "written_unbounded_terminal_slope_and_divisor_proofs",
          "actual_oriented_slopes_with_independent_signs_allowed",
          "exact_full_graph_and_composed_target_controls",
          "primary_characteristic_convention_and_retained_deep_topology_explicit",
          "finite_controls_not_universal_proof"]))
check("current frontier is ONLY earlierreset A4 with finite first postreset lists and all S3 closed",
      "CURRENT PF2RULED95" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "ALL S3 AND ALL STILL-SOLID A4 TRAJECTORIES CLOSED BOTH PENCILS",
          "TERMINAL BRIDGE LAW","TWO FIXED DIVISOR CONSTANTS",
          "ONLY REMAINING GENUINE CLASS","EARLIER-NONSOLID-RESET A4",
          "FIRST WINDING LISTS (2,8):2/4/8/16; (4,10):2/4",
          "SubsequentC3loops/later tailsNOTclosedorbounded",
          "NOpostresetloopcase remainsOPENwith>=2simpletails",
          "KIG28fullgenericA4NOTintermediateidentification","full goal remain OPEN"]) and
      all("PF2RULED95" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

prefix96=json.loads((HERE/"prefix_excess_closure_receipt.json").read_text())
check("all-degree local prefix positive determinant excess absorption and low determinant laws retain hypotheses",
      all(prefix96["general"].values()))
check("entire genuine positive iterated quartic class closes BOTH pencils all finite lengths",
      all(prefix96["quartic"].values()) and
      prefix96["winding2_arrival_alternatives"]==[[2,8,2],[4,2,-10]] and
      len(prefix96["balanced_local_prefix_controls"])==5 and
      len(prefix96["balanced_local_prefix_controls"][-1]["odd_windings"])==12 and
      len(prefix96["actual_full_reset_target_controls"])==8 and
      len(prefix96["actual_full_reset_target_controls"][-1]["tail"])==14)
check("prefix source countercontrols retain essential index parity local and final excess limits",
      prefix96["illegal_index_formal_countercontrol"]["D"]==-5 and
      prefix96["negative_determinant_valid_local_extension"]["delta"]==-94 and
      prefix96["parity_failure_old_delta2_counterchart"][-1]["delta"]==2 and
      prefix96["parity_failure_old_delta2_counterchart"][1]["delta"]==11 and
      prefix96["parity_failure_old_delta2_counterchart"][1]["kappa"]==6 and
      prefix96["missing_negative_final_excess_counterchart"]["delta"]==65)
check("genuine class closure retains open degenerate higher degree and full finite completion objective",
      all(prefix96["assurance"][k] for k in [
          "written_unbounded_prefix_inequality_and_complete_exhaustion_proofs",
          "index_order_and_even_bridge_parity_essential_countercontrols_retained",
          "actual_composed_integer_target_covers_with_full_marks",
          "primary_contact_JSJ_peripheral_inputs_reused_not_independently_reproved",
          "finite_controls_not_universal_proof"]) and
      all(v is False for k,v in prefix96["assurance"].items() if k not in [
          "written_unbounded_prefix_inequality_and_complete_exhaustion_proofs",
          "index_order_and_even_bridge_parity_essential_countercontrols_retained",
          "actual_composed_integer_target_covers_with_full_marks",
          "primary_contact_JSJ_peripheral_inputs_reused_not_independently_reproved",
          "finite_controls_not_universal_proof"]))
check("current task leaves closed A4 reset campaigns and targets actual genuine stage coverage",
      "CURRENT PF2RULED96" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "ENTIRE GENUINE POSITIVE ITERATED-TORUS QUARTIC CLASS CLOSED BOTH PENCILS AT ALL FINITE LENGTHS",
          "GENERAL DEGREE-FREE PREFIX EXCESS","LOW-DETERMINANT TERMINALITY",
          "ZEROorUNBOUNDEDpostresetC3loops","PRIORITY NOW GENUINE-STAGE COVERAGE AND DEGENERACIES",
          "fullcompatibleFINITE-MAP/COMPLETIONobjective","full goal remain OPEN"]) and
      all("PF2RULED96" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

unit97=json.loads((HERE/"unit_peripheral_coverage_receipt.json").read_text())
check("all-degree unit winding preserves both preferred generators signed slopes and colored linking",
      all(unit97["general"].values()) and len(unit97["signed_unit_winding_cover_controls"])==24 and
      {x["degree"] for x in unit97["signed_unit_winding_cover_controls"]}=={3,4,5,7})
check("unit parameter on nontrivial companion is not erased exact linking2 versus18",
      unit97["parameter1_countercontrol"]["before"]["matrix"]==[[-1,2],[2,-4]] and
      all(x["matrix"]==[[-9,18],[18,-36]] for x in unit97["parameter1_countercontrol"]["after"]))
check("primary reduced target coverage retains whole packet and explicit import boundaries",
      all(unit97["coverage"].values()) and
      unit97["assurance"]["primary_classification_independently_reproved"] is False and
      unit97["assurance"]["full_quartic_Keller_exclusion_registered"] is False and
      unit97["assurance"]["full_compatible_completion_goal_solved"] is False and
      unit97["assurance"]["earlier_whole_source_chain_and_finite_completion_audit_required"] is True)
check("current priority is whole source finite completion implication audit not fictional unit escape",
      "CURRENT PF2RULED97" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "ALL-DEGREE UNIT-WINDING BOTH-PERIPHERAL REDUCTION",
          "NO ESCAPE UNDER ALL WHOLE-PACKET HYPOTHESES",
          "PRIORITY WHOLE-SOURCE AND FINITE-COMPLETION IMPLICATION AUDIT",
          "full goal remains OPEN"]) and
      all("PF2RULED97" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

proper98=json.loads((HERE/"proper_finite_marked_end_receipt.json").read_text())
check("all-degree original proper finite marked-end and component-open geometric gates retained",
      all(proper98["geometry"].values()) and all(proper98["peripheral"].values()))
check("general full peripheral orbit lattice checked independently through degree5 with shear",
      [x["commuting_pairs"] for x in proper98["commuting_action_census"]]==[1,4,18,120,840] and
      proper98["nonzero_shear_control"][0]["lattice"]==[[5,0],[3,1]] and
      all(len(x["boundaries"])==1 and x["boundaries"][0]["longitudinal_degree"]==x["residue_degree"]
          for x in proper98["actual_finite_power_map_Hopf_controls"]) and
      proper98["actual_finite_power_map_Hopf_controls"][0]["image_order"]==4)
check("proper finite end theorem retains realization and whole-chain audit limits",
      all(proper98["assurance"][k] for k in [
          "written_all_degree_geometric_and_lattice_proofs",
          "finite_controls_not_universal_topology_or_realization_proof",
          "original_proper_finite_map_essential_countercontrol_retained",
          "complete_primary_PDF_pages_rendered_and_visually_read",
          "earlier_whole_chart_key_form_passport_and_exclusion_chain_audit_required"]) and
      all(v is False for k,v in proper98["assurance"].items() if k not in [
          "written_all_degree_geometric_and_lattice_proofs",
          "finite_controls_not_universal_topology_or_realization_proof",
          "original_proper_finite_map_essential_countercontrol_retained",
          "complete_primary_PDF_pages_rendered_and_visually_read",
          "earlier_whole_chart_key_form_passport_and_exclusion_chain_audit_required"]))
check("current priority advances from verified component-open actual-end to key-form whole-chart audit",
      "CURRENT PF2RULED98" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "COMPONENT-OPEN AFFINENESS VERIFIED AND ACTUAL MARKED END RECOGNIZED IN EVERY DEGREE",
          "GENERAL PERIPHERAL LATTICE(e,0),(a,f)",
          "PRIORITY KEY-FORM AND WHOLE-CHART NECESSITY AUDIT",
          "full quartic/Keller theorem registered", "Full goal remains OPEN"]) and
      all("PF2RULED98" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

monic99=json.loads((HERE/"monic_key_chart_audit_receipt.json").read_text())
check("all-degree monic key factor retains every exact leading and minimal degree hypothesis",
      all(monic99["monic_key"].values()) and len(monic99["monic_key_controls"])==5 and
      monic99["actual_affine_constant_residue_countercontrol"]["key_residue"]==1 and
      monic99["actual_affine_constant_residue_countercontrol"]["reduction_mod_key"]=="-X")
check("all positive index conditional whole chart support residue and minimal additive multiplier audited",
      all(monic99["whole_charts"].values()) and all(monic99["audit"].values()) and
      any(x["index"]==5 and x["boundary_order"]==3 and x["signed_support"][0]==-1 and x["minimal_multiplier"]==1
          for x in monic99["all_index_chart_controls"]))
published_quartic=next(x for x in index["entries"] if x["id"]=="HC4MYKR0")
check("published generic degree4 exclusion and marked etale finite completion application preserved",
      published_quartic["state"]=="proved" and
      claim_fingerprint(published_quartic)==monic99["published_quartic"]["reviewed_scope_fingerprint"] and
      all(monic99["published_quartic"][k] for k in [
          "unrestricted_generic_degree4_plane_Keller_exclusion_retained",
          "marked_dense_plane_etale_finite_degree4_completions_excluded",
          "quartic_realization_campaign_retired"]) and
      monic99["published_quartic"]["internal_quartic_chain_certified_by_external_theorem"] is False and
      monic99["published_quartic"]["arbitrary_finite_quartic_maps_with_interior_ramification_excluded"] is False)
check("current Wright priority preserves concurrent retirement and advances to higher degree",
      "WRIGHT CURRENT PF2RULED99" in action["next_step"] and
      all(t in action["next_step"] for t in [
          "PUBLISHED HC4MYKR0 QUARTIC RETIREMENT RETAINED",
          "PRIORITY HIGHER-DEGREE WRIGHT", "Full all-degree Wright goal remains OPEN"]) and
      all("PF2RULED99" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

chain100=json.loads((HERE/"polynomial_key_chain_receipt.json").read_text())
check("full polynomial chains and Laurent terminal necessity all degree index other primes",
      all(chain100["theorem"].values()) and
      [x["depth"] for x in chain100["deep_integer_circuits"]]==[64,1000] and
      all(x["expanded"] is False for x in chain100["deep_integer_circuits"]))
check("actual affine intermediate constant residue is not overexcluded and positive gap retained",
      chain100["actual_affine_control"]["passes_necessary_filter"] is True and
      chain100["actual_affine_control"]["Keller_realization_claimed"] is False and
      chain100["positive_prefix_gap"]["recurrences"][1]["beta"]==[-1,2])
check("unrestricted negative cone source inference rejected with actual closed wedge",
      chain100["cone_scope"]["closed_wedge_generators"]==[[1,-1],[1,2]] and
      chain100["cone_scope"]["terminal_weight"]==-4 and
      chain100["cone_scope"]["unrestricted_negative_range_cone_inference_imported"] is False and
      chain100["cone_scope"]["non_finite_generation_claimed"] is False)
check("current full theorem finish line closure first with external review deferred",
      "WRIGHT CURRENT PF2RULED100" in action["next_step"] and
      "FULL POLYNOMIAL KEY CHAINS EXCLUDED IN EVERY DEGREE" in action["next_step"] and
      "PRIORITY GENUINELY LAURENT TERMINAL CANCELLATIONS" in action["next_step"] and
      "Closure before external-review expense" in action["next_step"] and
      all("PF2RULED100" in (HERE/name).read_text() and "last missing implication" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]) and
      chain100["assurance"]["arbitrary_Laurent_boundary_excluded"] is False and
      chain100["assurance"]["research_goal_completed"] is False)

full101=json.loads((HERE/"full_peripheral_linking_receipt.json").read_text())
check("full longitude joining shear and rational linking retain exact integral torsion defects",
      all(full101["general"].values()) and
      full101["shear_and_torsion_control"]["data"]["matrix"]==[["1/2"]] and
      full101["shear_and_torsion_control"]["data"]["core_orders"]==[4] and
      full101["shear_and_torsion_control"]["data"]["cleared_relation_torsion_orders"]==[2] and
      full101["shear_and_torsion_control"]["data"]["closed_H1"]==dict(free=0,torsion=[8]))
check("all parameter all length theorem uses unimodular Zt proof rather than sample extrapolation",
      len(full101["universal_cable_theorem"]["exact_Zt_unimodular_certificates"])==2 and
      all(x["coefficient_ring"]=="Z[t]" and x["basis_determinant"] in [-1,1]
          and x["parameter_samples_used_as_proof"] is False
          for x in full101["universal_cable_theorem"]["exact_Zt_unimodular_certificates"]) and
      full101["universal_cable_theorem"]["finite_samples_used_as_universal_proof"] is False and
      any(len(x["parameters"])==12 for x in full101["all_length_family_controls"]))
check("new nondihedral seeds and sealed author overlap keep their actual scopes",
      [x["degree"] for x in full101["nontrivial_longitude_seeds"]]==[8,24] and
      [x["degree"] for x in full101["dihedral_overlap"]]==[3,5,7] and
      all(x["new_author_execution_claimed"] is False and x["author_diagonal_framing_identified"] is False
          for x in full101["dihedral_overlap"]) and
      full101["assurance"]["complete_filled_torsion_linking_form_computed"] is False and
      full101["assurance"]["new_covers_identified_as_dihedral"] is False)
check("current Wright goal remains actual higher degree exclusion or global finite algebra theorem",
      "WRIGHT CURRENT PF2RULED101" in action["next_step"] and
      "NEXT COMPLETE GOAL UNCHANGED" in action["next_step"] and
      "Closure before external-review expense" in action["next_step"] and
      all("PF2RULED101" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]) and
      full101["assurance"]["new_higher_degree_Keller_exclusion_proved"] is False and
      full101["assurance"]["research_goal_completed"] is False)

contact102=json.loads((HERE/"primitive_end_contact_receipt.json").read_text())
check("actual primitive one-line source canonical contact and torus parity hold all degree index",
      all(contact102["general"].values()) and all(contact102["primary_audit"].values()) and
      contact102["affine_index1_quadric_countercontrol"]["allowed"] is True and
      contact102["affine_index1_quadric_countercontrol"]["finite_Keller_target_claimed"] is False)
check("complete actual smooth six-sheet triple torus class closed at all parameters",
      all(v is True for k,v in contact102["six_sheet_triple_torus"].items()
          if k not in ["full_monodromy","unbounded_multiplier_order_table","finite_image_immersion_or_zero_excess_assumed"]) and
      contact102["six_sheet_triple_torus"]["full_monodromy"]=="A6" and
      contact102["six_sheet_triple_torus"]["finite_image_immersion_or_zero_excess_assumed"] is False and
      [x["fiber_order"] for x in contact102["six_sheet_triple_torus"]["unbounded_multiplier_order_table"]]==[3,2,1])
check("actual torus controls preserve full shear and exclude Keller compatibility rather than the covers",
      [x["torus"] for x in contact102["actual_A6_torus_cover_controls"]]==[[5,4],[5,3],[4,3]] and
      contact102["actual_A6_torus_cover_controls"][0]["result"]["shears"]==[2,0,0,0] and
      all(x["topological_cover_itself_excluded"] is False for x in contact102["actual_A6_torus_cover_controls"]) and
      contact102["nonprimitive_even_order_twist_control"]["core_orders"][0]==2)
check("new actual exclusion retains all singular multiple-prime iterated and global limits",
      contact102["assurance"]["actual_new_higher_degree_boundary_exclusion"] is True and
      all(contact102["assurance"][k] is False for k in [
          "singular_normalizations_excluded","multiple_prime_normalizations_excluded",
          "all_genuine_iterated_triple_knots_excluded","all_six_sheet_Keller_maps_excluded",
          "global_Laurent_finite_algebra_gap_closed","formal_verification","external_review","research_goal_completed"]))
check("current frontier now genuinely iterated triple knots singular multiple packets or actual finite algebra",
      "WRIGHT CURRENT PF2RULED102" in action["next_step"] and
      "ACTUAL HIGHER-DEGREE TORUS PACKET CLOSED" in action["next_step"] and
      "GENUINE ITERATED triple infinity knots" in action["next_step"] and
      "Closure before external-review expense" in action["next_step"] and
      all("PF2RULED102" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

multiple103=json.loads((HERE/"multiple_boundary_contact_receipt.json").read_text())
check("all-degree arbitrary-prime actual smooth A-type lens end exclusion",
      all(multiple103["assurance"][k] is True for k in [
       "written_universal_proof","primary_contact_and_filling_theorems_imported",
       "all_oriented_A_type_lens_ends_excluded_for_actual_smooth_Keller_normalization",
       "arbitrary_many_omitted_primes","all_field_degrees","external_review_deferred"]))
check("A-type exclusion preserves source scope and assurance limits",
      all(multiple103["assurance"][k] is False for k in [
       "singular_normalizations_excluded","other_lens_ends_excluded",
       "finite_checks_prove_universal_claim","formal_verification","external_review","global_JC2_closed"]))
check("zero-canonical quadric survives and formal blowup shows essential rank input",
      multiple103["controls"][0]["packet"]["canonical_nonzero"] is False and
      multiple103["controls"][2]["packet"]["theta"]=="-1" and
      multiple103["controls"][2]["packet"]["rank"]==2)
check("current all-prime A-type closure frontier and shared priorities retained",
      "WRIGHT CURRENT PF2RULED103" in action["next_step"] and
      "COMPLETE ALL-DEGREE A-TYPE END CLASS CLOSED" in action["next_step"] and
      all("PF2RULED103" in (HERE/name).read_text() for name in ["README.md","RESEARCH_GOAL.md"]))

spin104=json.loads((HERE/"fibre_transitive_spin_receipt.json").read_text())
check("all-degree full-shear independent cable recursion and integral closed transfer",
      all(spin104["universal"][k] is True for k in [
       "arbitrary_degree","arbitrary_finite_length","independent_coprime_parameters",
       "stable_nonzero_shear","filled_H1_preserved_integrally","core_class_multiplied_by_R"]) and
      spin104["universal"]["finite_samples_used_as_proof"] is False)
check("complete actual colored iterated exclusions retain extra unramified primes",
      all(spin104["actual_exclusions"].values()) and
      spin104["spin_source"]["natural_eight_budget"]["maximum_omitted_primes"]==2 and
      spin104["spin_source"]["natural_eight_budget"]["necessary_rank_interval_nonempty"] is False)
check("fibre hypothesis and genuine finite-map wrong-source countercontrols retained",
      spin104["nontransitive_torsion_countercontrol"]["H1"]==dict(free=2,torsion=[4,4]) and
      spin104["finite_degree8_polynomial_countercontrol"]["genuine"] is True and
      spin104["finite_degree8_polynomial_countercontrol"]["finite_Keller_realization"] is False and
      spin104["E8_formal_countercontrol"]["rank"]==8)
check("new actual exclusions preserve all source/color/global assurance limits",
      spin104["assurance"]["actual_higher_degree_genuine_iterated_colored_family_excluded"] is True and
      all(spin104["assurance"][k] is False for k in [
       "all_representations_of_underlying_knots_excluded","singular_natural_eight_sources_excluded",
       "all_degree_eight_Keller_maps_excluded","all_higher_degree_Keller_maps_excluded",
       "global_Laurent_finite_algebra_closed","topological_covers_declared_nonexistent",
       "formal_verification","external_review","research_goal_completed"]))
check("current PF104 colored iterated closure and complete goal frontier integrated",
      "WRIGHT CURRENT PF2RULED104" in action["next_step"] and
      "GENUINE ITERATED COLORED FAMILIES CLOSED" in action["next_step"] and
      all("PF2RULED104" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

triple105=json.loads((HERE/"triple_boundary_residue_receipt.json").read_text())
check("all-degree normal low-index residue theorem and primitive singular triple point",
      all(triple105["general"].values()))
check("PF104 colored iterated family now closes for arbitrary normal source and any triple knot joining",
      all(triple105["normal_iterated_application"].values()))
check("primary normalized power LOCATION and smooth attribution retained explicitly",
      all(triple105["attribution"].values()))
check("genuine unbounded higher-index finite counterfamily retains wrong-source distinction",
      triple105["genuine_unbounded_finite_counterfamily"]["formula_proved_for_all_positive_s"] is True and
      triple105["genuine_unbounded_finite_counterfamily"]["finite_Keller_realization"] is False and
      [(x["transverse_index"],x["residue_degree"],x["quotient_degree"])
       for x in triple105["genuine_unbounded_finite_counterfamily"]["representative_controls"]]==[(3,1,4),(5,2,16),(7,3,36)])
check("new normal-source theorem does not become small-degree or full-JC exclusion",
      all(triple105["assurance"][k] is False for k in [
       "all_eight_sheet_Keller_maps_excluded","all_indices_residue_one_claimed",
       "all_singular_triple_points_excluded","all_representation_of_underlying_knots_excluded",
       "global_Laurent_finite_algebra_closed","formal_verification","external_review","research_goal_completed"]))
check("current all-degree residue and normal iterated closure integrated with global goal",
      "WRIGHT CURRENT PF2RULED105" in action["next_step"] and
      "ALL-DEGREE LOW-INDEX RESIDUE ONE; NORMAL ITERATED FAMILY CLOSED" in action["next_step"] and
      all("PF2RULED105" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

quotient106=json.loads((HERE/"cyclic_quotient_morse_receipt.json").read_text())
check("all-index actual quotient order conductor and extremal residue constraints",
      all(quotient106["unconditional"][k] is True for k in [
       "all_field_degrees","normal_singular_source","arbitrary_other_primes",
       "quotient_order_at_most_d_minus2","all_extremal_quotients_residue1",
       "all_omitted_residue1_through_degree10","exact_conductor_balance_all_indices"]) and
      quotient106["unconditional"]["residue1_indices"]==[1,2,3,4,6])
check("scalar Morse quadratic budget explicitly CONDITIONAL",
      quotient106["scalar_Morse"]["hypothesis_required"] is True and
      quotient106["scalar_Morse"]["automatically_satisfied_claimed"] is False and
      quotient106["scalar_Morse"]["low_degree_residue_failure_requires_degenerate_quadratic_pencil"] is True)
check("sharp genuine local countercontrols fail actual whole-plane field budget",
      [(x["local_degree"],x["actual_Keller_degree_at_least"],x["realization"])
       for x in quotient106["sharp_local_countercontrols"]]==[(4,5,False),(16,17,False),(36,37,False)] and
      all(x["actual_Keller_source"] is False for x in quotient106["non_Morse_countercontrols"]))
check("through-ten residue theorem preserves full goal and non-Morse gap",
      all(quotient106["assurance"][k] is False for k in [
       "all_ten_or_lower_degree_maps_excluded","general_local_residue_bound_without_Morse_proved",
       "all_indices_residue1","degenerate_pencil_actual_Keller_realization","global_finite_algebra_closed",
       "formal_verification","external_review","research_goal_completed"]))
check("current quotient conductor and conditional quadratic frontier integrated",
      "WRIGHT CURRENT PF2RULED106" in action["next_step"] and
      "ALL-INDEX QUOTIENT ORDER/CONDUCTOR AND CONDITIONAL QUADRATIC RESIDUE BUDGET" in action["next_step"] and
      all("PF2RULED106" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

leading107=json.loads((HERE/"coprime_leading_forms_receipt.json").read_text())
check("complete coprime-leading-form degenerate sector retains explicit hypothesis",
      leading107["theorem"]["coprime_leading_forms_hypothesis_required"] is True and
      leading107["theorem"]["automatic_hypothesis"] is False and
      leading107["theorem"]["nontrivial_residue_minimum_degree"]==49)
check("projective pencil gives distinct plus-one divisor and scalar quotient",
      leading107["theorem"]["residue_divides_e_plus1"] is True and
      leading107["theorem"]["scalar_marked_quotient"] is True)
check("genuine non-Gorenstein passport controls preserve wrong source open",
      any(x["non_Gorenstein"] and x["residue_degree"]>1 for x in leading107["passport_family_controls"]) and
      all(x["complement"]=="C* times C" and not x["actual_Keller_source"] for x in leading107["passport_family_controls"]))
check("coprime initial sector preserves unresolved common tangent and global goal",
      all(leading107["assurance"][k] is False for k in [
       "all_degenerate_pencils_classified","general_local_residue_bound_proved",
       "all_indices_residue1","actual_Keller_controls","global_finite_algebra_closed",
       "formal_verification","external_review","research_goal_completed"]))
check("current coprime-leading class integrated preserving quartic retirement",
      "WRIGHT CURRENT PF2RULED107" in action["next_step"] and
      all("PF2RULED107" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

ade108=json.loads((HERE/"ade_end_exclusion_receipt.json").read_text())
check("complete all-degree actual smooth oriented ADE and all-minus-two-tree end class excluded",
      all(ade108["theorem"][k] is True for k in [
       "all_field_degrees","arbitrary_omitted_prime_count","arbitrary_indices",
       "smooth_actual_normalization_required","whole_plane_open_required",
       "oriented_negative_ADE_link_required","all_ADE_ends_excluded",
       "all_negative_definite_minus_two_tree_plumbing_ends_excluded"]))
check("ADE exclusion preserves primary import and source orientation limits",
      ade108["assurance"]["deep_contact_and_filling_theorems_imported"] is True and
      all(ade108["assurance"][k] is False for k in [
       "independent_whole_classification_proofs","singular_normalizations_excluded",
       "all_spherical_ends_excluded","opposite_orientations_excluded",
       "intersection_form_identifies_end","global_finite_algebra_closed",
       "formal_verification","external_review","research_goal_completed"]))
check("ADE canonical energy controls retain zero-Chern references and essential minimality",
      all(x["reference"]["canonical_nonzero"] is False and
          x["blowup"]["theta"]==x["reference"]["theta"] and
          x["blowup"]["rank"]==x["reference"]["rank"]+1 for x in ade108["exact_controls"]))
check("current complete smooth ADE frontier integrated preserving predecessor residue results",
      "WRIGHT CURRENT PF2RULED108" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED107" in action["next_step"] and
      all("PF2RULED108" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

normal109=json.loads((HERE/"normal_resolution_chern_receipt.json").read_text())
check("complete all-NORMAL ADE and all-minus-two-tree end class includes singular sources",
      all(normal109["theorem"][k] is True for k in [
       "all_field_degrees","all_normal_actual_sources","arbitrary_omitted_prime_count",
       "arbitrary_indices_and_residues","whole_marked_plane_required",
       "all_negative_ADE_ends_excluded","all_negative_definite_minus_two_TREE_ends_excluded",
       "oriented_end_identification_required"]))
check("minimal resolution Stein deformation carries a NONZERO rational Chern class",
      all(normal109["theorem"][k] is True for k in [
       "algebraic_Picard_Chern_injectivity_written","minimal_resolution_admits_Stein_deformation",
       "deformed_filling_simply_connected","rational_Chern_nonzero",
       "b2_at_least_all_resolution_boundary_components"]))
check("normal source bridge does not falsely preserve target map or original exactness",
      all(normal109["assurance"][k] is False for k in [
       "independent_whole_deformation_proof","original_resolution_Kahler_form_exact",
       "deformed_finite_target_map_holomorphic_claimed","deformed_whole_plane_holomorphic_claimed",
       "arbitrary_nonminimal_resolution_Stein_deformation_claimed","singularity_smoothing_claimed",
       "all_spherical_ends_excluded","reversed_orientations_excluded",
       "homology_or_lattice_identifies_end","global_finite_algebra_closed",
       "formal_verification","external_review","research_goal_completed"]))
check("genuine finite singular A1 ADE-end control survives its wrong plane-open hypothesis",
      normal109["genuine_wrong_open_control"]["finite_degree"]==12 and
      normal109["genuine_wrong_open_control"]["canonical_vector"]==[8,4] and
      normal109["genuine_wrong_open_control"]["principal_unit_relation"]==[2,1] and
      normal109["genuine_wrong_open_control"]["actual_Keller_source"] is False)
check("current ALLNORMAL ADE closure integrated without dropping prior scopes",
      "WRIGHT CURRENT PF2RULED109" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED108" in action["next_step"] and
      all("PF2RULED109" in (HERE/name).read_text() and "Domrina2000" in (HERE/name).read_text()
          for name in ["README.md","RESEARCH_GOAL.md"]))

circle110=json.loads((HERE/"planar_circle_bundle_receipt.json").read_text())
check("universal planar and complete normal circle-bundle source reduction",
      all(circle110["theorem"][k] is True for k in [
       "all_field_degrees","all_normal_actual_sources","whole_marked_plane_required",
       "planar_contact_requires_integral_boundary_Chern_nonzero","all_S3_ends_excluded",
       "all_negative_SU2_quotient_ends_excluded","all_Lp1_singular_sources_excluded",
       "all_Lp1_multiple_prime_sources_excluded","surviving_Lp1_source_smooth_sole_prime",
       "sharp_index_range_and_parity","sharp_all_p_affine_exact_source_controls",
       "oriented_actual_end_identification_required"]))
check("circle-bundle reduction retains integral torsion and finite-map realization limits",
      all(circle110["assurance"][k] is False for k in [
       "rational_boundary_Chern_zero_suffices","graph_base_genus_implies_contact_planarity",
       "nonstandard_p4_fill_uniqueness_required","original_mcduff_audited",
       "all_Lp1_ends_excluded","all_spherical_ends_excluded",
       "all_integral_homology_sphere_ends_excluded","reversed_orientations_excluded",
       "source_controls_have_finite_Keller_targets","homology_identifies_end",
       "global_finite_algebra_closed","formal_verification","external_review",
       "research_goal_completed"]))
check("genuine all-p source controls preserve whole-plane exact volume and no finite target",
      len(circle110["genuine_source_controls"])==6 and
      all(x["source_affine"] and x["whole_plane_open"] and x["globally_exact_volume"] and
          x["numerical_index"]==x["p"]-1 and not x["finite_target_supplied"]
          for x in circle110["genuine_source_controls"]))
check("current PF110 integrated while retaining normal ADE closure",
      "WRIGHT CURRENT PF2RULED110" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED109" in action["next_step"] and
      all("PF2RULED110" in (HERE/n).read_text() and "PF2RULED109" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

peripheral111=json.loads((HERE/"whole_plane_peripheral_receipt.json").read_text())
check("universal all-normal Thom lattice and complete marked peripheral reconstruction",
      all(peripheral111["theorem"][k] is True for k in [
       "all_field_degrees","all_normal_actual_sources","whole_marked_plane_required",
       "integral_boundary_Picard_Chern_isomorphism","exact_resolved_rank_r_plus_m",
       "all_normal_source_intersection_lattices_negative","all_actual_end_H1_finite",
       "strict_omitted_end_knots_generate_all_H1","all_omitted_knots_covered_exterior_torsion_free",
       "complete_F_recovers_integral_peripheral_module_if_cores_generate",
       "partial_coverage_outside_prime_budget","arbitrary_b_fibre_transitive_cable_sector_written",
       "complete_six_sheet_colored_iterated_family_excluded","arbitrary_cable_length_and_stated_parameters"]))
check("complete peripheral result preserves coverage and actual-representation scope",
      all(peripheral111["assurance"][k] is False for k in [
       "finite_sweep_proves_universal_theorem","new_deep_classification_imported",
       "arbitrary_cover_F_determines_full_torsion","zero_row_defects_certify_generation",
       "coverage_hypothesis_omitted","additional_unramified_images_assumed_absent",
       "smooth_source_assumed_for_six_sheet_exclusion","six_sheet_application_is_new_exclusion_class","all_six_sheet_Keller_maps_excluded",
       "uncolored_knots_alone_excluded","all_Lp1_source_controls_excluded",
       "global_finite_algebra_closed","formal_verification","external_review",
       "research_goal_completed"]))
check("genuine six-sheet seed has hidden homology that the Euler budget cannot repair",
      peripheral111["exact_seed"]["generation"]["closed_H1"]==dict(free=0,torsion=[4]) and
      peripheral111["exact_seed"]["generation"]["unseen_quotient"]==dict(free=0,torsion=[2]) and
      all(x["unbounded_formula"]["required_outside_prime_count"]==1 and
          x["unbounded_formula"]["remaining_Euler_budget"]==0
          for x in peripheral111["cable_controls"]))
check("current PF111 integrated with all prior Wright checkpoints retained",
      "WRIGHT CURRENT PF2RULED111" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED110" in action["next_step"] and
      all("PF2RULED111" in (HERE/n).read_text() and "PF2RULED110" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

high112=json.loads((HERE/"high_index_budget_receipt.json").read_text())
check("universal all-normal high-index costs and complete near-maximal class exclusion",
      all(high112["theorem"][k] is True for k in [
       "all_normal_actual_sources","all_field_degrees_main_high_index_theorem",
       "whole_marked_plane_and_all_omitted_primes_required",
       "high_index_residue_one_and_injective_normalization",
       "high_index_requires_positive_own_excess","high_index_requires_distinct_nonproper_image",
       "high_index_e_plus_r_at_most_d_minus1",
       "arbitrary_residue_high_image_mass_normalization_injective",
       "high_image_mass_plus_unit_residue_count_at_most_d_minus2",
       "irreducible_whole_nonproper_TOTAL_mass_half_degree_bound",
       "irreducible_whole_nonproper_normalized_excess_lower_bound",
       "all_omitted_indices_at_most_d_minus3_for_d_at_least5",
       "complete_near_maximal_e_equals_d_minus2_class_excluded",
       "irreducible_whole_nonproper_set_half_degree_bound",
       "complete_next_strip_necessary_profile","actual_quotient_order_at_most_d_minus4",
       "degree5_normalization_smooth","degree6_only_A1_quotients",
       "stronger_partial_link_excess_reserve",
       "complete_degree17_stated_colored_iterated_family_excluded",
       "arbitrary_cable_length_and_stated_parameters"]))
check("high-index proof retains the distinct Chau scopes and open next strip",
      high112["assurance"]["old_maximal_index_exclusion_credited"] is True and
      all(high112["assurance"][k] is False for k in [
       "independently_reproved_Chau_and_Orevkov","finite_sweep_proves_theorem",
       "singular_normal_sources_excluded_by_assumption","contractible_component_alone_forbidden",
       "whole_nonproper_set_confused_with_branch_curve","next_strip_all_models_excluded",
       "numerical_profile_realized","all_degree5_or_degree6_maps_excluded",
       "all_degree17_maps_excluded","uncolored_knots_without_representation_excluded",
       "arbitrary_finite_etale_open_suffices","source_Cstar_action_or_target_descent_assumed",
       "withdrawn_Jelonek_2020_claim_imported","global_finite_algebra_closed",
       "formal_verification","external_review","research_goal_completed"]))
check("degree17 seed genuinely has zero row defects but a hidden quotient",
      high112["exact17_seed"]["peripheral"]["cleared_relation_torsion_orders"]==[1,1,1] and
      high112["exact17_seed"]["generation"]["unseen_quotient"]==dict(free=0,torsion=[3]) and
      high112["exact17_seed"]["generation"]["closed_H1"]==dict(free=0,torsion=[51]))
check("current all-degree near-maximal closure integrated while full objective stays open",
      "WRIGHT CURRENT PF2RULED112" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED111" in action["next_step"] and
      all("PF2RULED112" in (HERE/n).read_text() and "PF2RULED111" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

local113=json.loads((HERE/"local_monodromy_branch_receipt.json").read_text())
check("universal all-normal localcapacity and complete weighted full-branch exclusion",
      all(local113["theorem"][k] is True for k in [
       "all_field_degrees","all_normal_actual_sources",
       "original_finite_normalization_whole_plane_all_primes_required",
       "every_local_cluster_at_most_d_minus_r","complete_finite_fibre_count_lower_bound",
       "no_singleton_TOTAL_normalization_fibre","local_monodromy_never_transitive",
       "positive_weighted_FULL_branch_union_excluded","polynomial_target_rectifications_included",
       "contractible_irreducible_FULL_branch_excluded","extra_unramified_nonproper_images_allowed",
       "high_index_forces_different_RAMIFIED_image",
       "all_indices_at_most_d_minus4_for_d_at_least7",
       "complete_prior_e_equals_d_minus3_strip_excluded",
       "shared_image_mass_plus_unit_residue_count_at_most_d_minus3",
       "irreducible_FULL_branch_mass_at_most_half_degree",
       "new_two_RAMIFIED_prime_necessary_profile",
       "quotient_order_at_most_d_minus5_for_d_at_least7"]))
check("local-global proof does not confuse source fibres, fullbranch scope or source actions",
      all(local113["assurance"][k] is False for k in [
       "independent_whole_classification_proof","finite_samples_prove_theorem",
       "source_singleton_fibres_excluded","source_surjectivity_claimed",
       "arbitrary_finite_covers_excluded","one_weighted_component_suffices",
       "weighted_analytic_germ_suffices","infinity_knot_alone_certifies_global_branch_class",
       "nonpositive_weights_suffice","simultaneous_component_rectification_assumed",
       "source_Cstar_action_lift_assumed","source_symmetry_descent_assumed",
       "new_next_strip_models_excluded","numerical_profile_realized",
       "all_degree7_or_degree8_maps_excluded","homology_identifies_ADE_end",
       "global_finite_algebra_closed","formal_verification","external_review",
       "research_goal_completed"]))
check("genuine smooth and singular finite branch controls survive the wrong-plane hypothesis",
      local113["genuine_finite_controls"]["cubic_fibre_length"]==3 and
      local113["genuine_finite_controls"]["singular_A1_fibre_length"]==12 and
      local113["genuine_finite_controls"]["Keller_plane_examples"] is False)
check("current PF113 supersedes old open profile without rewriting predecessor proof bytes",
      "WRIGHT CURRENT PF2RULED113" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED112" in action["next_step"] and
      "CURRENT FRONTIER CORRECTION BY PF113" in records["PF2RULED112"]["scope"] and
      all("PF2RULED113" in (HERE/n).read_text() and "PF2RULED112" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

cusp114=json.loads((HERE/"torus_cusp_excess_receipt.json").read_text())
check("all-degree cyclic-axis torus-cusp excess and proportional cap integrated",
      all(cusp114["theorem"][k] is True for k in [
       "all_field_degrees","all_normal_actual_sources",
       "original_finite_map_whole_plane_all_primes_required",
       "cyclic_axis_torus_cusp_mass_at_most_five_excess_minus_two_fibre_order",
       "nontrivial_local_residue_and_longitude_included","exact_orbifold_Euler_identity",
       "actual_unramified_meridian_fillings_retained",
       "universal_index_cap_floor_five_d_minus17_over6_for_d_at_least9",
       "entire_prior_e_equals_d_minus4_strip_excluded","shared_image_all_residue_cusp_cost"]))
check("new local torus-cusp theorem retains its precise hypotheses and full goal",
      all(cusp114["assurance"][k] is False for k in [
       "finite_samples_prove_theorem","arbitrary_plane_branches_covered",
       "scalar_Morse_assumed","local_residue_one_assumed","regular_fibre_identity_assumed",
       "local_quotient_identifies_global_end","arithmetic_gate_realizes_cover",
       "actual_wrong_plane_finite_models_excluded","source_symmetry_descent_assumed",
       "global_finite_algebra_closed","formal_verification","external_review",
       "research_goal_completed"]))
check("current PF114 closes prior two-ramified-prime profile with sealed predecessor retained",
      "WRIGHT CURRENT PF2RULED114" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED113" in action["next_step"] and
      "CURRENT FRONTIER CORRECTION BY PF114" in records["PF2RULED113"]["scope"] and
      all("PF2RULED114" in (HERE/n).read_text() and "PF2RULED113" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

milnor115=json.loads((HERE/"residual_milnor_receipt.json").read_text())
check("arbitrary-branch residual Milnor and minimal-excess classification integrated",
      all(milnor115["theorem"][k] is True for k in [
       "all_degrees_indices_local_residues","any_singular_chosen_plane_branch",
       "smooth_or_cyclic_axis_source_pair_required","sole_critical_branch_required",
       "exact_residual_Milnor_identity","normalized_excess_at_least_local_residue",
       "minimal_excess_complete_numerical_classification",
       "index_at_least4_forces_excess_at_least_q_plus1",
       "arbitrary_global_image_nonimmersion_costs","nonprimitive_omitted_residue_strict_excess",
       "all_omitted_residues_one_through_degree15","degree15_corollary_already_follows_from_PF114"]))
check("all-branch Milnor proof retains chosen branch, nonlinear source coordinate and global limits",
      all(milnor115["assurance"][k] is False for k in [
       "one_Puiseux_pair_assumed","scalar_Morse_assumed","source_P_coordinate_assumed_linear",
       "residual_curve_assumed_smooth","global_image_assumed_contractible",
       "ordinary_preimages_counted_as_omitted_primes",
       "PF114_linear_mass_bound_extended_to_arbitrary_branch",
       "finite_samples_prove_theorem","jet_stabilization_assumed_without_Nakayama",
       "all_degree_residue_one_claimed","degree15_inversion_claimed",
       "unique_analytic_equivalence_class_claimed","arithmetic_gates_realize_covers",
       "actual_finite_algebra_closed","formal_verification","external_review",
       "research_goal_completed"]))
check("higher-characteristic and high-index genuine finite controls remain allowed",
      milnor115["higher_characteristic_control"]["image_conductor"]==16 and
      milnor115["higher_characteristic_control"]["residual_Milnor"]==60 and
      milnor115["higher_characteristic_control"]["local_finite_germ_excluded"] is False and
      milnor115["genuine_finite_high_index_control"]["residual_Milnor"]==4 and
      milnor115["genuine_finite_high_index_control"]["excluded_locally"] is False)
check("current PF115 advances residual-divisor gate and retains proportional PF114 bound",
      "WRIGHT CURRENT PF2RULED115" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED114" in action["next_step"] and
      all("PF2RULED115" in (HERE/n).read_text() and "PF2RULED114" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

log116=json.loads((HERE/"logarithmic_pullback_receipt.json").read_text())
check("actual logarithmic map and sharp unbounded finite family integrated",
      all(log116["theorem"][k] is True for k in [
       "all_finite_smooth_plane_germ_degrees","reduced_target_contains_ALL_branch_values_required",
       "exact_logarithmic_pullback_isomorphism","exact_ramification_weighted_residue_map",
       "actual_curve_Jacobian_duality_formula","target_coefficient_differential_values_scale_by_residue",
       "source_normalization_saturation_loses_differential_gaps","unbounded_genuine_finite_sharp_cusp_family",
       "local_five_eta_coefficient_sharp","non_Morse_common_leading_zero_family","all_ordinary_residual_maps_primitive"]))
check("module bridge does not identify surface trace dual or import extra preprint/global conclusions",
      all(log116["assurance"][k] is False for k in [
       "original_surface_trace_dual_identified_with_curve_residue","ramification_diagonal_dropped",
       "multibranch_gluing_dropped","target_coefficient_module_identified_with_source_saturation",
       "newer_preprint_semifield_transfer_imported","Bernstein_root_transport_proved",
       "orbifold_BMY_equality_assumed","global_Wright_bound_optimality_claimed",
       "all_coefficient_fields_claimed_rational","sharp_family_is_Keller_completion",
       "finite_samples_prove_theorem","actual_global_finite_algebra_closed",
       "formal_verification","external_review","research_goal_completed"]))
check("actual fixed analytic modules and genuine finite invariant fibre verified",
      log116["fixed_actual_module_control"]["differential_extra_values"]==[11,15] and
      log116["fixed_actual_module_control"]["Tjurina"]==14 and
      log116["exact_finite_examples"][1]["upstairs_length"]==50 and
      log116["exact_finite_examples"][1]["invariant_length"]==10)
check("PF116 current global algebra gate preserves prior all-branch reduction",
      "WRIGHT CURRENT PF2RULED116" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED115" in action["next_step"] and
      all("PF2RULED116" in (HERE/n).read_text() and "PF2RULED115" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

trace117=json.loads((HERE/"nilradical_trace_receipt.json").read_text())
check("nilradical trace restriction and actual quotient descent integrated",
      all(trace117["theorem"][k] is True for k in [
       "all_finite_smooth_plane_degrees","ALL_branch_values_required",
       "nonreduced_base_change_retained","nilradical_annihilator_equals_weighted_trace",
       "embedded_curve_codifferent_Dr_B","all_multibranch_gluing_retained",
       "actual_small_cyclic_quotient_descent","branch_permuting_small_groups_included",
       "reduced_downstairs_branch_curve_Gorenstein","actual_surface_trace_restriction_closed",
       "actual_upstairs_residue_codifferent_bridge","invariant_quotient_residue_bridge"]))
check("local trace bridge preserves singular-source and global-algebra limits",
      all(trace117["assurance"][k] is False for k in [
       "inverse_Jacobian_reduced_as_fraction","nonreduced_base_change_replaced_by_reduced",
       "curve_ramification_diagonal_discarded","normalization_gluing_discarded",
       "ambient_singular_dual_assumed_free","intrinsic_singular_Kahler_residues_identified",
       "nonplane_curve_Jacobian_duality_imported","sharp_finite_covers_excluded",
       "new_global_index_exclusion","actual_global_finite_algebra_closed",
       "global_trace_moment_regularity_proved","Bernstein_connection_transport_proved",
       "finite_samples_prove_theorem","formal_verification","external_review","research_goal_completed"]))
check("PF117 current global gate preserves sealed local predecessor",
      "WRIGHT CURRENT PF2RULED117" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED116" in action["next_step"] and
      all("PF2RULED117" in (HERE/n).read_text() and "PF2RULED116" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]) and
      trace117["controls"]["actual_nonGorenstein_m2_weights"]==[8,1,1])

half118=json.loads((HERE/"half_support_inertia_receipt.json").read_text())
check("global half-support theorem and complete plane action classes integrated",
      all(half118["theorem"][k] is True for k in [
       "all_plane_field_degrees","whole_plane_all_omitted_primes_required",
       "exact_global_inertia_support_budget","nontrivial_inertia_moves_at_most_half",
       "large_minimal_degree_actions_excluded","both_exceptional_PSL2_11_plane_actions_excluded",
       "natural_projective_linear_point_and_hyperplane_actions_excluded_in_stated_range",
       "nonbinary_affine_linear_actions_excluded","binary_affine_half_support_profile",
       "binary_profile_normalization_smooth","binary_profile_double_values_force_empty_source_fibre"]))
check("plane group exclusion retains binary equality and global finite-algebra limits",
      all(half118["assurance"][k] is False for k in [
       "subgroup_census_rerun","finite_samples_prove_theorem","arbitrary_affine_dimension_exclusion",
       "all_semilinear_actions_excluded","binary_affine_family_excluded","all_degree6_maps_excluded",
       "actual_plane_binary_realization_constructed","source_surjectivity_assumed",
       "global_trace_moment_regularity_proved","actual_global_finite_algebra_closed",
       "historical_priority_claimed","formal_verification","external_review","research_goal_completed"]))
check("PF118 current gate keeps source gluing and nonexcluded binary32 profile",
      "WRIGHT CURRENT PF2RULED118" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED117" in action["next_step"] and
      all("PF2RULED118" in (HERE/n).read_text() and "PF2RULED117" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]) and
      half118["exceptional11"]["minimum_support"]==8 and
      half118["binary_affine32"]["ramified_prime_count"]==8 and
      half118["binary_affine32"]["unramified_primes_plus_normalized_excess"]==15 and
      not half118["binary_affine32"]["binary_affine_family_excluded"])

simple119=json.loads((HERE/"simple_ramification_ledger_receipt.json").read_text())
check("all-simple exact multiplicity ledger and two half-support cases integrated",
      all(simple119["theorem"][k] is True for k in [
       "all_simple_index2_plane_normalizations","normalization_smooth_and_all_residues1",
       "local_length_equals_image_multiplicity_plus1","singular_image_pair_consecutive",
       "local_symmetric_monodromy","exact_global_multiplicity_ledger","exact_pointwise_source_count",
       "half_support_involutory_two_cases","half_support_at_most_one_ordinary_cusp",
       "higher_and_multiple_cusps_excluded_in_half_support_case","half_support_cusp_local_group_S3"]))
check("simple ramification result does not assume generic-cover, infinity or analytic closure",
      all(simple119["assurance"][k] is False for k in [
       "finite_samples_prove_theorem","global_symmetric_group_assumed",
       "Kulikov_total_ramification_birationality_assumed","Orevkov_all_preimage_bijectivity_assumed",
       "immersed_affine_curve_complement_assumed_cyclic","analytic_L4_gap_closed",
       "binary_affine_family_excluded","all_degree6_maps_excluded",
       "fixed_cusp_sheets_all_assumed_retained","source_surjectivity_assumed",
       "actual_finite_Keller_realization_constructed","actual_global_finite_algebra_closed",
       "historical_priority_claimed","formal_verification","external_review","research_goal_completed"]))
check("PF119 current gate preserves two genuine open alternatives and prior action closures",
      "WRIGHT CURRENT PF2RULED119" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED118" in action["next_step"] and
      all("PF2RULED119" in (HERE/n).read_text() and "PF2RULED118" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]) and
      [(c["total_excess"],c["unramified_primes"]) for c in simple119["half_support32"]]==[(0,15),(8,7)])

product120=json.loads((HERE/"primitive_product_action_receipt.json").read_text())
check("complete primitive product action and prime-degree reductions integrated",
      all(product120["theorem"][k] is True for k in [
       "all_plane_degrees","whole_original_normalization_all_primes_required",
       "every_inertia_support_strictly_less_five_sixths","involutory_inertia_support_strictly_less_two_thirds",
       "all_primitive_Cartesian_product_actions_base_at_least3_excluded","product_action_exclusion_CFSG_free",
       "complete_primitive_necessary_container_reduction","all_prime_degrees_force_natural_alternating_or_symmetric",
       "nonalternating_primitive_at_most_two_branch_images","nonalternating_every_inertia_support_at_most_half"]))
check("product exclusion separates published CFSG import from remaining whole-plane obligations",
      product120["assurance"]["published_CFSG_classification_imported"] is True and
      product120["assurance"]["full_PDF_table_visually_audited"] is True and
      all(product120["assurance"][k] is False for k in [
       "finite_samples_prove_theorem","independent_CFSG_classification_proved",
       "all_wreath_actions_including_imprimitive_excluded","all_intermediate_models_assumed_planes",
       "one_meridian_per_image_assumed_to_generate","all_listed_containers_realizable",
       "all_prime_degree_maps_invertible","binary_classical_families_excluded",
       "alternating_subset_families_excluded","binary_affine_family_excluded",
       "target_blowup_affine_boundary_hypotheses_preserved","actual_global_finite_algebra_closed",
       "historical_priority_claimed","formal_verification","external_review","research_goal_completed"]))
check("PF120 current gate preserves previous exact cusp ledger and shared priorities",
      "WRIGHT CURRENT PF2RULED120" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED119" in action["next_step"] and
      all("PF2RULED120" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

disjoint121=json.loads((HERE/"disjoint_branch_supports_receipt.json").read_text())
check("universal disjoint branch supports and complete subset/classical exclusions integrated",
      all(disjoint121["theorem"][k] is True for k in [
       "all_plane_degrees","whole_original_normalization_all_primes_required",
       "universal_disjoint_nonidentity_BRANCH_meridian_pair",
       "all_faithful_interior_subset_actions_excluded","subset_exclusion_CFSG_free",
       "all_three_remaining_binary_classical_primitive_actions_excluded",
       "primitive_only_natural_alternating_symmetric_or_binary_affine",
       "binary_affine_linear_part_generated_by_transvections",
       "binary_affine_exact_two_half_support_cusp_profiles_retained"]))
check("disjoint-support theorem source and remaining-goal boundaries retained",
      disjoint121["assurance"]["published_CFSG_classification_imported"] is True and
      disjoint121["assurance"]["GM_all_element_inclusive_half_primary_web_text_imported"] is True and
      all(disjoint121["assurance"][k] is False for k in [
       "original_GM_PDF_audited","original_GM_source_bytes_locally_pinned",
       "review_record_misrepresented_as_original_source","OCR_subset_parameter_bound_used",
       "independent_CFSG_classification_proved","finite_samples_prove_theorem",
       "local_support_disjointness_requires_N_smooth",
       "one_chosen_meridian_per_image_assumed_to_generate",
       "simultaneous_rectification_of_components_assumed",
       "binary_affine_family_excluded","natural_alternating_symmetric_actions_excluded",
       "all_imprimitive_actions_excluded","intermediate_sources_assumed_planes",
       "all_listed_containers_realizable","actual_global_finite_algebra_closed",
       "historical_priority_claimed","formal_verification","external_review","research_goal_completed"]))
check("PF121 current gate preserves prior exact cusp ledger and shared priorities",
      "WRIGHT CURRENT PF2RULED121" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED120" in action["next_step"] and
      all("PF2RULED121" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

block122=json.loads((HERE/"block_quotient_monodromy_receipt.json").read_text())
check("hereditary block-quotient reduction for all monodromy types integrated",
      all(block122["theorem"][k] is True for k in [
       "all_plane_degrees_all_monodromy_types","whole_original_normalization_all_primes_required",
       "every_nontrivial_block_quotient_exact_support_defect_budget",
       "every_block_quotient_full_positive_weighted_branch_excluded",
       "every_block_quotient_irreducible_contractible_full_branch_excluded",
       "every_block_quotient_disjoint_nonidentity_BRANCH_pair","quotient_degrees2_and3_excluded",
       "all_interior_subset_block_quotients_excluded","all_primitive_product_action_block_quotients_excluded",
       "classical_transvection_support_pointwise_stabilizer_trivial",
       "all_three_binary_classical_primitive_block_quotients_excluded",
       "every_primitive_block_quotient_natural_Alt_Sym_or_binary_affine",
       "binary_affine_quotient_one_branch_original_support_exact_half"]))
check("block quotient proof keeps original geometry and intermediate-source boundaries",
      block122["assurance"]["published_CFSG_classification_imported"] is True and
      block122["assurance"]["prior_GM_primary_web_import_boundary_retained"] is True and
      all(block122["assurance"][k] is False for k in [
       "original_GM_PDF_audited","independent_CFSG_proof","finite_samples_prove_theorem",
       "intermediate_models_assumed_planes","original_N_assumed_smooth",
       "all_high_mass_residue_degrees_assumed_one","quotient_involution_assumed_to_lift_to_involution",
       "arbitrary_original_disjoint_pair_assumed_nontrivial_downstairs",
       "arbitrary_abstract_quotient_representations_covered","all_subgroup_composition_factors_classified",
       "PF119_cusp_profile_transferred_to_every_intermediate",
       "degree4_plane_inversion_applied_to_intermediate","all_natural_or_affine_block_quotients_excluded",
       "all_imprimitive_monodromy_excluded","actual_global_finite_algebra_closed",
       "historical_priority_claimed","formal_verification","external_review","research_goal_completed"]))
check("PF122 current gate retains previous original cusp profiles and shared priorities",
      "WRIGHT CURRENT PF2RULED122" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED121" in action["next_step"] and
      all("PF2RULED122" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

axis123=json.loads((HERE/"intermediate_critical_axis_receipt.json").read_text())
check("all-intermediate critical-axis descent and affine profiles integrated",
      len(axis123["theorem"])==11 and all(v is True for v in axis123["theorem"].values()))
check("intermediate descent keeps imported geometry and remaining-goal boundaries",
      axis123["assurance"]["Pedersen_v3_and_published_Giraldo_Roeder_theorems_imported"] is True and
      axis123["assurance"]["intermediate_local_geometry_proved_before_PF106_PF119_application"] is True and
      axis123["assurance"]["primary_HTML_sources_pinned"] is True and
      all(v is False for k,v in axis123["assurance"].items() if k not in [
       "Pedersen_v3_and_published_Giraldo_Roeder_theorems_imported",
       "intermediate_local_geometry_proved_before_PF106_PF119_application","primary_HTML_sources_pinned"]))
check("PF123 current gate preserves prior all-map block reduction and shared priorities",
      "WRIGHT CURRENT PF2RULED123" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED122" in action["next_step"] and
      all("PF2RULED123" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

chern124=json.loads((HERE/"intermediate_chern_betti_receipt.json").read_text())
check("all-intermediate Chern and end obstructions with exact affine Betti profiles integrated",
      len(chern124["theorem"])==11 and all(v is True for v in chern124["theorem"].values()))
check("intermediate topology credits retained work and keeps global closure boundaries",
      chern124["assurance"]["prior_canonical_and_pseudo_plane_results_credited"] is True and
      chern124["assurance"]["prior_deep_filling_and_deformation_imports_explicitly_retained"] is True and
      all(v is False for k,v in chern124["assurance"].items() if k not in [
       "prior_canonical_and_pseudo_plane_results_credited",
       "prior_deep_filling_and_deformation_imports_explicitly_retained"]))
check("PF124 current gate retains exact missing-fibre interface and shared priorities",
      "WRIGHT CURRENT PF2RULED124" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED123" in action["next_step"] and
      all("PF2RULED124" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

longitude125=json.loads((HERE/"binary_affine_longitude_receipt.json").read_text())
check("PF125 all-degree binary affine longitude and filling restrictions integrated",
      len(longitude125["theorem"])==6 and all(longitude125["theorem"].values()))
check("PF125 preserves original coverage and open global/family boundaries",
      all(longitude125["assurance"][k] is True for k in [
       "prior_PF60_rank_two_argument_credited","prior_dihedral_authors_code_comparison_retained"]) and
      all(v is False for k,v in longitude125["assurance"].items() if k not in [
       "prior_PF60_rank_two_argument_credited","prior_dihedral_authors_code_comparison_retained"]))
check("PF125 current peripheral gate retains shared priorities and PF124 end checkpoint",
      "WRIGHT CURRENT PF2RULED125" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED124" in action["next_step"] and
      all("PF2RULED125" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

seifert126=json.loads((HERE/"seifert_framing_cancellation_receipt.json").read_text())
check("PF126 componentwise actual Seifert cancellation and binary residual covers integrated",
      len(seifert126["theorem"])==9 and all(seifert126["theorem"].values()))
check("PF126 preserves original hypotheses, credits and remaining global boundaries",
      all(seifert126["assurance"][k] is True for k in [
       "prior_PF60_PF125_degree4_longitude_result_credited","prior_authors_code_comparison_retained",
       "original_all_prime_and_actual_Seifert_framing_hypotheses_preserved"]) and
      all(v is False for k,v in seifert126["assurance"].items() if k not in [
       "prior_PF60_PF125_degree4_longitude_result_credited","prior_authors_code_comparison_retained",
       "original_all_prime_and_actual_Seifert_framing_hypotheses_preserved"]))
check("PF126 current gate retains shared priorities, PF125 algebra and PF124 end checkpoint",
      "WRIGHT CURRENT PF2RULED126" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED125" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED124" in action["next_step"] and
      all("PF2RULED126" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

shear127=json.loads((HERE/"quotient_corrected_shear_receipt.json").read_text())
check("PF127 all-normal all-index quotient-corrected shear and complete peripheral exclusions integrated",
      len(shear127["theorem"])==8 and all(shear127["theorem"].values()))
check("PF127 credits retained inputs and keeps individual/global cyclicity and open-goal boundaries",
      all(shear127["assurance"][k] is True for k in [
       "prior_Borisov_geometry_and_residue_inputs_credited",
       "prior_polar_conductor_import_and_genuine_cover_families_credited",
       "actual_Seifert_framing_and_all_original_primes_retained","primary_HTML_sources_pinned"]) and
      all(v is False for k,v in shear127["assurance"].items() if k not in [
       "prior_Borisov_geometry_and_residue_inputs_credited",
       "prior_polar_conductor_import_and_genuine_cover_families_credited",
       "actual_Seifert_framing_and_all_original_primes_retained","primary_HTML_sources_pinned"]))
check("PF127 current gate preserves shared priorities and its sealed framing/end predecessors",
      "WRIGHT CURRENT PF2RULED127" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED126" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED124" in action["next_step"] and
      all("PF2RULED127" in (HERE/n).read_text() and "PF2RULED119" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

compression128=json.loads((HERE/"pseudoplane_compression_receipt.json").read_text())
check("PF128 all-degree compression and exact minimal counterexample reductions integrated",
      len(compression128["theorem"])==8 and all(compression128["theorem"].values()))
check("PF128 explicit published inputs and minimal/nonminimal boundaries preserved",
      all(compression128["assurance"][k] is True for k in [
       "published_classification_and_cover_chart_explicitly_adopted",
       "primary_PDF_bytes_pinned_and_relevant_pages_visually_read",
       "written_unbounded_proofs_and_original_whole_source_retained"]) and
      all(v is False for k,v in compression128["assurance"].items() if k not in [
       "published_classification_and_cover_chart_explicitly_adopted",
       "primary_PDF_bytes_pinned_and_relevant_pages_visually_read",
       "written_unbounded_proofs_and_original_whole_source_retained"]))
check("PF128 current navigation and sealed predecessor retained",
      "WRIGHT CURRENT PF2RULED128" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED127" in action["next_step"] and
      all("PF2RULED128" in (HERE/n).read_text() and "PF2RULED127" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

euler129=json.loads((HERE/"hereditary_euler_excess_receipt.json").read_text())
check("PF129 exact hereditary Euler, low slack, index and residue bounds integrated",
      len(euler129["theorem"])==8 and all(euler129["theorem"].values()))
check("PF129 retains all-source and intermediate/original distinctions",
      all(euler129["assurance"][k] is True for k in [
       "prior_local_axis_polar_solid_torus_and_classification_inputs_credited",
       "written_unbounded_Euler_proof_includes_all_singular_points_and_image_collisions",
       "actual_whole_original_plane_and_all_critical_primes_retained"]) and
      all(v is False for k,v in euler129["assurance"].items() if k not in [
       "prior_local_axis_polar_solid_torus_and_classification_inputs_credited",
       "written_unbounded_Euler_proof_includes_all_singular_points_and_image_collisions",
       "actual_whole_original_plane_and_all_critical_primes_retained"]))
check("PF129 current navigation and sealed PF128 retained",
      "WRIGHT CURRENT PF2RULED129" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED128" in action["next_step"] and
      all("PF2RULED129" in (HERE/n).read_text() and "PF2RULED128" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

ruled130=json.loads((HERE/"ruled_chart_compression_receipt.json").read_text())
check("PF130 relative Euler, general chart compression and exact bad-fibre control integrated",
      len(ruled130["theorem"])==8 and all(ruled130["theorem"].values()))
check("PF130 explicit input and chart/ruling scope boundaries retained",
      all(ruled130["assurance"][k] is True for k in [
       "prior_axis_Euler_universal_cover_and_deck_inputs_credited",
       "general_ruling_pi1_and_bundle_structure_explicitly_adopted",
       "entire_original_plane_and_all_normal_intermediates_retained",
       "exact_volume_countercontrol_regular_primitive_algebraically_certified"]) and
      all(v is False for k,v in ruled130["assurance"].items() if k not in [
       "prior_axis_Euler_universal_cover_and_deck_inputs_credited",
       "general_ruling_pi1_and_bundle_structure_explicitly_adopted",
       "entire_original_plane_and_all_normal_intermediates_retained",
       "exact_volume_countercontrol_regular_primitive_algebraically_certified"]))
check("PF130 current navigation and sealed PF129 retained",
      "WRIGHT CURRENT PF2RULED130" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED129" in action["next_step"] and
      all("PF2RULED130" in (HERE/n).read_text() and "PF2RULED129" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

for num,name,tc,truekeys in [
 (131,"proper_binary_longitude_receipt.json",5,["classification_explicitly_adopted_via_pinned_primary_restatement","given_module_and_separate_symplectic_proposition_checked","actual_residual_cover_needed_beyond_common_form","original_plane_and_invisible_primes_retained"]),
 (132,"radical_tower_exclusion_receipt.json",6,["coordinator_coupling_deduction_credited_and_audited","classical_ruling_and_orbifold_inputs_explicitly_adopted","HC4MYML0_full_component_selection_input_essential","entire_original_plane_and_all_omitted_primes_retained","specified_branch_permutation_hypothesis_retained","no_minimal_degree_assumption"])]:
    recent=json.loads((HERE/name).read_text())
    check("PF"+str(num)+" complete theorem integrated",len(recent["theorem"])==tc and all(recent["theorem"].values()))
    check("PF"+str(num)+" exact adoption and nonclosure boundaries retained",
          all(recent["assurance"][k] is True for k in truekeys) and
          all(v is False for k,v in recent["assurance"].items() if k not in truekeys))
    check("PF"+str(num)+" current work and prior proof retained",
          "WRIGHT CURRENT PF2RULED"+str(num) in action["next_step"] and
          all("PF2RULED"+str(num) in (HERE/n).read_text() and "PF2RULED130" in (HERE/n).read_text()
              for n in ["README.md","RESEARCH_GOAL.md"]))

deck133=json.loads((HERE/"deck_normalized_euler_receipt.json").read_text())
check("PF133 exact deck budget and all-inertia natural tower extension integrated",
      len(deck133["theorem"])==6 and all(deck133["theorem"].values()))
truth133=["PF132_surface_lemma_and_coordinator_HC4_credit_retained",
          "actual_free_deck_quotient_used_instead_of_entire_intermediate_identification",
          "all_residual_fixed_block_and_invisible_support_retained",
          "entire_original_plane_and_all_omitted_primes_retained",
          "no_minimal_degree_or_inertia_class_assumption"]
check("PF133 retained hypothesis and nonclosure boundaries",
      all(deck133["assurance"][k] is True for k in truth133) and
      all(v is False for k,v in deck133["assurance"].items() if k not in truth133))
check("PF133 removes only same-natural-action inertia restriction and preserves PF132",
      "WRIGHT CURRENT PF2RULED133" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED132" in action["next_step"] and
      all("PF2RULED133" in (HERE/n).read_text() and "PF2RULED132" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

orb134=json.loads((HERE/"finite_orbifold_groups_receipt.json").read_text())
check("PF134 finite Mobius containers and actual complete universal rulings integrated",
      len(orb134["theorem"])==6 and all(orb134["theorem"].values()))
truth134=["retained_classical_inputs_broader_scope_explicitly_adopted",
          "original_plane_entire_intermediates_and_actual_free_deck_quotient_retained",
          "coordinator_PF132_and_HC4_component_selection_credit_retained",
          "generic_deck_and_full_monodromy_roles_distinguished",
          "bad_normalized_multiplicity_fibres_remain_open"]
check("PF134 source, model and nonclosure boundaries retained",
      all(orb134["assurance"][k] is True for k in truth134) and
      all(v is False for k,v in orb134["assurance"].items() if k not in truth134))
check("PF134 current work and sealed PF133 retained",
      "WRIGHT CURRENT PF2RULED134" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED133" in action["next_step"] and
      all("PF2RULED134" in (HERE/n).read_text() and "PF2RULED133" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

cot135=json.loads((HERE/"cotangent_deck_trace_receipt.json").read_text())
check("PF135 invariant cotangent field and exact relative trace integrated",
      len(cot135["theorem"])==5 and all(cot135["theorem"].values()))
truth135=["PF134_actual_equivariant_ruling_and_coordinator_HC4_credit_retained",
          "all_original_plane_and_omitted_primes_retained_in_actual_theorem",
          "rational_coordinates_poles_and_actual_lattice_boundary_explicit",
          "negative_control_has_genuine_regular_etale_target_on_its_stated_open",
          "no_new_classification_or_Razar_type_theorem_import"]
check("PF135 exact hypothesis and nonclosure boundaries retained",
      all(cot135["assurance"][k] is True for k in truth135) and
      all(v is False for k,v in cot135["assurance"].items() if k not in truth135))
check("PF135 current work and sealed PF134 retained",
      "WRIGHT CURRENT PF2RULED135" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED134" in action["next_step"] and
      all("PF2RULED135" in (HERE/n).read_text() and "PF2RULED134" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

orbit136=json.loads((HERE/"common_omitted_orbits_receipt.json").read_text())
check("PF136 actual almost-surjective quotient and refined common-orbit budget integrated",
      len(orbit136["theorem"])==5 and all(orbit136["theorem"].values()))
truth136=["actual_original_plane_and_target_map_preserved",
          "entire_original_normalization_discriminant_and_all_defects_retained",
          "auxiliary_quotient_not_mislabeled_entire_intermediate",
          "finite_missed_points_and_mixed_fibres_remain_open",
          "shared_MO0_chart_and_MA0_residual_audit_credited"]
check("PF136 actual source, shared-credit and nonclosure boundaries retained",
      all(orbit136["assurance"][k] is True for k in truth136) and
      all(v is False for k,v in orbit136["assurance"].items() if k not in truth136))
check("PF136 current work and sealed PF135 retained",
      "WRIGHT CURRENT PF2RULED136" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED135" in action["next_step"] and
      all("PF2RULED136" in (HERE/n).read_text() and "PF2RULED135" in (HERE/n).read_text()
          for n in ["README.md","RESEARCH_GOAL.md"]))

# Current cross-programme integration is metadata validation, not mathematical replay.
hc4=ROOT/"jacobian-followup/hc4-triangular-graphs"
for cid,proof,checker,context,receipt in [
    ("HC4MYMP0","COMMON_UNRAMIFIED_PENCIL_AND_AUTOMATIC_ADDITIVE_COMPLETION.md",
     "verify_automatic_additive_completion.py","automatic-additive-completion-context.json",
     "receipt-automatic-additive-completion.json"),
    ("HC4MYMQ0","CYCLIC_DECK_SOURCE_IMAGE_AND_THE_UNIQUE_MULTIPLE_FIBRE.md",
     "verify_cyclic_deck_source_image.py","cyclic-deck-source-image-context.json",
     "receipt-cyclic-deck-source-image.json")]:
    e=next(x for x in index["entries"] if x["id"]==cid)
    r=json.loads((hc4/receipt).read_text())
    c=json.loads((hc4/"sources"/context).read_text())
    check(cid+" retained proved result and receipt",e["state"]=="proved" and r["outcome"]=="PASS")
    for key,p in [("proof_sha256",hc4/proof),("checker_sha256",hc4/checker),
                  ("context_sha256",hc4/"sources"/context)]:
        check(cid+" "+key,hashlib.sha256(p.read_bytes()).hexdigest()==r[key])
    for p,h in c["retained_proof_and_source_sha256"].items():
        check(cid+" retained pin "+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
    check(cid+" no assurance upgrade",all(r[k] is False for k in
          ["independent_replay","formal_verification","external_review"]))
    check(cid+" shared work gate current",cid in action["claims"] and
          action["reviewed_claims"][cid]==claim_fingerprint(e))
check("current Wright navigation credits MP0 and MQ0",all(
      "HC4MYMP0" in (HERE/n).read_text() and "HC4MYMQ0" in (HERE/n).read_text()
      for n in ["README.md","RESEARCH_GOAL.md","DIHEDRAL_LINKING_RECURSION_ROUTE.md",
                "LAURENT_FINITE_ALGEBRA_FINISH_LINE.md"]))
check("current Wright work gate preserves KR and full open goal",
      action["next_step"].startswith("CURRENT KR0/1/2:") and
      "WRIGHT CURRENT MP0/MQ0:" in action["next_step"])

volume137=json.loads((HERE/"separated_variable_volume_receipt.json").read_text())
check("PF137 all-degree separated-variable exact-volume result integrated",
      len(volume137["theorem"])==5 and all(volume137["theorem"].values()))
truth137=["all_degree_proof_written_not_inferred_from_controls",
          "prior_classical_atlas_and_monomial_primitive_credited",
          "same_natural_tower_and_shared_MP0_MQ0_scopes_retained",
          "x_dependent_mixed_fibres_not_excluded",
          "surviving_single_root_target_and_full_trace_gate_open"]
check("PF137 no surface classification or full-target closure overstated",
      all(volume137["assurance"][k] is True for k in truth137) and
      all(v is False for k,v in volume137["assurance"].items() if k not in truth137))
check("PF137 current work and MP0/MQ0 retained",
      "WRIGHT CURRENT PF2RULED137" in action["next_step"] and
      "WRIGHT CURRENT MP0/MQ0:" in action["next_step"] and
      all("PF2RULED137" in (HERE/n).read_text() for n in ["README.md","RESEARCH_GOAL.md"]))

atlas138=json.loads((HERE/"translation_atlas_exactness_receipt.json").read_text())
check("PF138 whole-atlas exactness and unbounded moving-label controls integrated",
      len(atlas138["theorem"])==5 and all(atlas138["theorem"].values()))
truth138=["all_degree_written_proof_not_bounded_scan",
          "primary_affine_cohomology_statements_imported_with_assurance_limits",
          "PF137_constant_label_scope_and_prior_atlas_credited",
          "shared_MP0_MQ0_actual_target_distinctions_retained",
          "actual_full_finite_order_trace_energy_gate_open"]
check("PF138 no target realization, classification or assurance upgrade",
      all(atlas138["assurance"][k] is True for k in truth138) and
      all(v is False for k,v in atlas138["assurance"].items() if k not in truth138))
check("PF138 current work preserves PF137 and shared deck reductions",
      "WRIGHT CURRENT PF2RULED138" in action["next_step"] and
      "WRIGHT CURRENT PF2RULED137" in action["next_step"] and
      "WRIGHT CURRENT MP0/MQ0:" in action["next_step"] and
      all("PF2RULED138" in (HERE/n).read_text() for n in ["README.md","RESEARCH_GOAL.md"]))

descent139=json.loads((HERE/"replica_degree_descent_receipt.json").read_text())
check("PF139 actual degree descent and all-degree square-free replica restriction integrated",
      len(descent139["theorem"])==5 and all(descent139["theorem"].values()))
truth139=["written_all_degree_proof_not_partition_census",
          "shared_MR0_and_retained_PF130_chart_construction_credited",
          "selected_source_image_not_identified_with_entire_intermediate_etale_locus",
          "minimality_needed_for_exclusion_without_smaller_map_alternative",
          "cyclic_trivial_deck_and_general_full_trace_goal_open"]
check("PF139 minimality, actual target and no assurance upgrade retained",
      all(descent139["assurance"][k] is True for k in truth139) and
      all(v is False for k,v in descent139["assurance"].items() if k not in truth139))
check("PF139 current work and shared MR0 credit retained",
      "WRIGHT CURRENT PF2RULED139" in action["next_step"] and
      "HC4MYMR0" in action["claims"] and
      all("PF2RULED139" in (HERE/n).read_text() for n in ["README.md","RESEARCH_GOAL.md"]))

missed140=json.loads((HERE/"missed_intermediate_lines_receipt.json").read_text())
check("PF140 all-intermediate missed lines and exact source-image budget integrated",
      len(missed140["theorem"])==5 and all(missed140["theorem"].values()))
truth140=["actual_whole_original_plane_and_full_finite_normalizations_retained",
          "critical_axis_avoidance_handles_whole_curves_not_punctured_images",
          "prior_boundary_and_Giraldo_Roeder_import_limits_retained",
          "shared_KV0_MR0_MS0_and_PF139_scopes_credited",
          "finite_point_gaps_and_actual_full_trace_energy_goal_open"]
check("PF140 whole-curve, selected-open and no assurance-upgrade boundaries retained",
      all(missed140["assurance"][k] is True for k in truth140) and
      all(v is False for k,v in missed140["assurance"].items() if k not in truth140))
check("PF140 current work and shared MS0 credit retained",
      "WRIGHT CURRENT PF2RULED140" in action["next_step"] and
      "HC4MYMS0" in action["claims"] and
      all("PF2RULED140" in (HERE/n).read_text() for n in ["README.md","RESEARCH_GOAL.md"]))

velocity141=json.loads((HERE/"source_velocity_multipliers_receipt.json").read_text())
check("PF141 exact full-order target multiplier and bounded inverse-leading clearing integrated",
      len(velocity141["theorem"])==5 and all(velocity141["theorem"].values()))
truth141=["actual_whole_source_and_all_ramified_unramified_primes_retained",
          "target_and_source_multiplier_ideals_distinguished",
          "canonical_Lie_derivative_divergence_shift_not_omitted",
          "shared_MS0_MT0_LO1_and_full_inverse_pole_ledger_credited",
          "unscaled_finite_landing_and_full_trace_energy_goal_open"]
check("PF141 divergence, full trace list, source-frame and nonclosure boundaries retained",
      all(velocity141["assurance"][k] is True for k in truth141) and
      all(v is False for k,v in velocity141["assurance"].items() if k not in truth141))
check("PF141 current work and shared MT0 credit retained",
      "WRIGHT CURRENT PF2RULED141" in action["next_step"] and
      "HC4MYMT0" in action["claims"] and
      all("PF2RULED141" in (HERE/n).read_text() for n in ["README.md","RESEARCH_GOAL.md"]))

kummer142=json.loads((HERE/"partial_kummer_charts_receipt.json").read_text())
check("PF142 partial Kummer plane comparison and all-LND cyclic factor result integrated",
      len(kummer142["theorem"])==5 and all(kummer142["theorem"].values()))
truth142=["all_degree_written_proof_not_parameter_census",
          "entire_cover_ramification_retained_then_deleted_on_comparison_open",
          "source_image_not_entire_intermediate_etale_locus",
          "MR0_PF130_PF139_ML0_MU0_classical_limits_credited",
          "full_original_normalization_and_trace_energy_goal_retained"]
check("PF142 minimality, whole ramified-cover and finite-lattice boundaries retained",
      all(kummer142["assurance"][k] is True for k in truth142) and
      all(v is False for k,v in kummer142["assurance"].items() if k not in truth142))
check("PF142 current work, MU0 credit and mixed minimal frontier retirement retained",
      "WRIGHT CURRENT PF2RULED142" in action["next_step"] and
      "HC4MYMU0" in action["claims"] and
      all("PF2RULED142" in (HERE/n).read_text() for n in ["README.md","RESEARCH_GOAL.md"]))

visiting=set();visited=set()
def visit(id_):
    assert id_ not in visiting, "new dependency cycle: "+id_
    if id_ in visited:return
    visiting.add(id_)
    for d in records[id_]["dependencies"]:
        if d in records:visit(d)
    visiting.remove(id_);visited.add(id_)
for id_ in records:visit(id_)
check("new dependency subgraph is acyclic",len(visited)==len(records))


report=dict(as_of="2026-11-16",focused_outcome="PASS",check_count=len(checks),checks=checks,
            mathematical_control_receipts=receipts,
            repository_navigation=dict(outcome="NOT_CERTIFIED_UNRELATED_REPOSITORY_ERRORS",
               known_first_error="KDSQ6: missing document consumer artifacts/generated-results/README.md",
               last_check_navigation_error="KDSQ6: missing document consumer artifacts/generated-results/README.md",
               last_render_navigation_error="KDSQ6: missing document consumer artifacts/generated-results/README.md",
               first_missing_path_still_absent=not(ROOT/"artifacts/generated-results/README.md").exists(),
               generated_status_and_indexes_manually_edited=False,
               missing_research_outputs_reconstructed=False),
            research_goal_completed=False,
            remaining=["CURRENT PF142: partial Kummer covers with selected whole unramified components give actual etale plane comparisons without an entire universal-cover premise. Every non-power source-LND-equivariant factor descends to generic degree at most d/2. At MINIMAL generic degree ALL such proper factors are actual cyclic deck fixed fields, for arbitrary replicas, and trivial deck gives full derivative field for EVERY source LND. Former mixed/coprime2+3/2+N minimal frontier CLOSED; MU0 arbitrary-degree conditional order/trace results retained. NEXT actual cyclic factors, necessary equivariance, or full original finite S/Lambda trace-energy implication; full goal ACTIVE OPEN",
                       "CURRENT PF141: actual generic first-source-velocity pole q+e gives exact TARGET S/Lambda stabilizer ideal product H^(1+maxceil(q/e)), with Lie divergence shift regular. Firstd power traces recover all exponents/nonproper image curves, including unram escape. FULL source-finite inverse leading square clears allfour Jacobian entries in every degree, but target clearings under failure are NOTlocallyfinite source fields. Shared MS0/MT0 rawR/LO1 and old complete inverse ledger credited. NEXT unscaled finite-order/trace-energy landing, not more jet or denominator escalation; full goal ACTIVE OPEN",
                       "CURRENT PF140: all intermediate completely missed divisors are whole disjoint smooth A1s with full smooth isomorphic residue1/Eta0 power-map packets. Actual maximal source-image T deletes r lines from entire V: exact b=chiT+r+A+Eta. General finite point gaps inside T remain. Shared MR0/PF139 mixed fibres give chiT>=3/quotientb>=r+5; full finiteN and MS0 first-derivative integrality/trace-energy gate retained. No selected-open finite-cover or pi1/deck confusion. Full mixed/cyclic/trivial-deck/non-equivariant goal ACTIVE OPEN",
                       "CURRENT PF139: shared MR0 actual source-image multiplicities plus PF130 chart criterion give smaller actual plane Keller degree d*g/m. In minimal generic degree, factors cyclic or unique mixed normalized fibre degree>=5 with repeated replica root. ALL square-free replicas only actual quadratic deck/centred odd; even degrees or odd generic deck permit no proper equivariant factor. No raw field-to-lattice leap, no selected-image/entire-locus identification. Mixed/cyclic/trivial-deck/full trace-energy goal remains OPEN",
                       "CURRENT PF138: specified AFFINE whole translation atlas has exact canonical volume IFF all transition residues zero; primary affine cohomology imported. Actual MQ0 W satisfies atlas. Explicit b^m*y=z(z-(1-m)bprime) has exact volume/whole plane/scalar units/pi1=1/free source Ga and ANY prescribed number of reduced degenerate fibres, no etale target supplied. No further canonical-volume primitive obstruction after these residues; no fibre-count bound from the listed control properties alone. NEXT actual cyclic/trivial-deck full finite order/trace-energy; full goal ACTIVE OPEN",
                       "CURRENT PF137: exactness excludes all separated a(x)y=p(z) models with multiple base roots or linear a, square-free p degree>=2, by rational primitive pole-degree argument. Single-root m>=2 passes exactness, its compatible target remains OPEN. General affine-line surfaces/x-dependent mixed fibres not classified. NEXT actual cyclic/trivial-deck full finite order/trace-energy; full goal OPEN",
                       "CURRENT MP0/MQ0: all actual generic deck groups cyclic in every degree. Original plane surjects etale onto its affine cyclic source-image quotient, with one irreducible multiple fibre and reduced ordinary fibres. This is not necessarily the entire intermediate etale locus. NEXT actual cyclic/trivial-deck full finite normalization, trace-dual multiplication or marked-energy integrality; full Wright/HC4/JC2 OPEN.",
                       "Retained PF136: common omitted curve-orbit deletion and Euler budget stay true. MP0 closes every noncyclic-deck survivor; MQ0 eliminates finite point gaps in the actual source image. No index5 noncyclic mixed-fibre case remains to revisit.",
                       "Retained PF135: conditional cotangent noncyclic-deck normal form and abstract controls preserved. MP0 excludes its actual noncyclic-completion premise. Generic field traces do not settle the remaining actual finite order.",
                       "Retained PF134: finite containers for arbitrary intermediate surface pi1 remain valid; MP0 removes dihedral/platonic generic DECK alternatives. Do not transfer deck cyclicity to arbitrary intermediate pi1 or full monodromy.",
                       "Retained PF133: all-inertia exclusion of the same natural radical tower action in degrees16,32,64,... remains proved. Other actions require their own deck audit; q1 and knot representations retained. Former index5 noncyclic profile is closed by MP0.",
                       "Retained PF132: original knot/Wright/HC4 Euler<=2 exact-volume cyclicity chain explicitly retains specified lifted-transposition inertia. The q>=2 tower exclusion remains valid. PF133 strengthens the same natural action; MP0 gives universal actual-deck cyclicity. Knot complements retain their representations and lack these surface constraints.",
                       "PRIORITY PF131: every proper binary linear primitive quotient has trivial preferred longitude by explicit McLaughlin common-form import via primary Brumer-Kramer restatement. Original PDF not obtained; full-AGL/nonidentity escape remains, no all-affine-group closure",
                       "Retained PF130: exactrelativeEuler t*c1-c2=Delta+Eta_res/Eta2=tEta1+Eta_res; equalityfiniteetaleentireopens. UniversalOPENA2chartcompressionANYEuler, chartneednotontoV; affineA1rulingcomponent=gcdcriterion constructs universalchart. Minimalproperrelativefactor throughreduced-fibreaffinerulingEXCLUDED; relativeCt/noresidual/b0=t*c-1. PrimitivefullD Ct; nonprimitivefullDcyclicforc<=3 butordercanexceedt. SharedHC4MYML0complete-baseperiodproofcredited. Genuineexactunitvolume affine surface fibre2plus3 hasglobalprimitivecertificate: exactnessdoesNOTsupplycomponentselection, noactualfiniteKellertarget.64cheapcontrols+writtenall-degreeproof/broadergeneralrulingpi1adoption, noautochart/ruling/minimality/globalclosure. Fullfinite-target/naturalAS/binary/trace-energy/EMPTYORIGINALfibreOPEN; sharedpriorities/quarticPLANEretirement/goalACTIVEpreserved",
                       "Retained PF129: exactALL-intermediateEuler c=chiV=b-A-Eta>=1, supportslackB controlsEta<=B-1. Localfive-sixthsboundtransfersinINTERMEDIATEdegreeb>=9; ALLcriticalresidues1 throughintermediateb15 evenunboundedoriginaldegree. Eulersaturationc1 triggersPF128smalleractualplaneKellerorminimalfullCt; maxsupportminimalNsm/allorigcriticalf1; B<=2 eta1onlysmoothindex2orA1index3ordinarycusp. Visibleoriginalnf<=b-2, invisibleprimesoutside.66cheapcontrols+writtenunboundedproof/creditedretaineddeepinputs, no newclassification/low-sheetimport. Fullcyclicfinite-target/naturalAS/binaryaffine/trace-energy/EMPTYORIGINALfibreOPEN, noalloriginalf1fromsmallb/planeintermediate/sourceboundarybasis orarithmeticrealization. Sharedpriorities/quarticPLANEretirement/goalACTIVEpreserved",
                       "Retained PF128: chi1entireetaleintermediatecompressesactualplaneKellerdegree b*m,m=|pi1V| divides t; minimalgenericdegreeforcesuniversalcycliccover/fullgenericD=Ct/b0t-1 WITHOUTprimitivity. Binarychi1exactlydegree4onecusps=>minimalNsm/solecriticalimage/noinvisiblebranches/solvableG/uniforminertiak|t andpi1M=C_(t/k). GeneralDelta0iffrelativeetale; primitivequotientidentifiesfullnormalizeroutsidechi1.64cheapcontrols andexplicitpublisheduniversalA2chartadoption/retainedDPclassification. Passingcyclicfinite-targetcaseOPEN; noallmapsminimal/nochiV>1closure/nointermediateplane/noautomaticpi1M0/solvableinversion/globalquadraticrefinement. Fullgoal/sharedtrace-energy/EMPTYORIGINALfibre/quarticPLANEretirementpreserved",
                       "Retained PF127: all genuine finite intermediates inherit full critical power-map LOCATION (onequotpoint, branchonlythere,n1=>f1); entireaxisnormalizedindexnroot givesa=f(2delta_C-1)+n^-1 mode forALLindices/residues/normalcyclicsources. EveryINDIVIDUALjointorbitCYCLIC/gcd(e,a,f)1, ALLzero-shearNONCOPRIMEe/fclassesEXCLUDED; evene=>a+fODD closesALLeven-index/even-residue/EVENshear. No knot/end/degree/Morsecondition. Quotientnrecovery/exactshared-imageoddTCRT/exactsingular-orderfactorizationgatebeforelinking.150controls~0.22sec73204KiB andgenuineALL-parameterfinitecoversmatchshears2/6/3, newprimaryHTMLnormalizedroots/priorpolarimportpins. Fullperipheral/globalmonodromyNOTassumedcyclic; nointermediateplane/integralcoverage/fullfinitealgebraclosure. NEXT actual marked-cover contradiction/globalfinitealgebra/sharedHC4specificEMPTYORIGINALfull-sheetevenradialprofile, naturalA/S/remainingbinaryGLOBALcompatibilityOPEN, sharedpriorities/quarticPLANEretirementpreserved",
                       "Retained PF126: all involutory branch components of EVERY genuine intermediate have pointwise fixed movedletters under preferredlongitude/zero criticalshears. Gaussnormalpush linking2delta+M_C andHessianzerosM_C cancelinactualSeifertframing. Binaryaffineonlylocalizedtranslations/order<=2/allfixedoralltworesidualprofiles; PF125order4/mixedperipheralclassesEXCLUDED allb, algebraiccontrolsretainedvalid. Fullbinaryresidualcurvecovermonodromy elementary2 W=span2v_o oddcontactcolumns, infinity=sumcolumns, exactRHgenus/infinitycount; noodd=>splitA1/lambda1, oneodd=>allCstar.144 cheapcontrols pluswrittenparameterproof; commonnondegeneratealternatinglinearform=>alloddcontactcolumnsv/pairsum0=>preferredlambdaIDENTITYinALLsuchbinarydegrees, residualb/4 identicalquadraticcovers/genusv_o-1/twoinfinitypoints, noentiremonodromyfamilyclosure/plane/integralcoverageassumption/noauthorscodererun. NEXT actual marked linking/globalfinitealgebra orsharedHC4specific EMPTY original full-sheet evenradialprofile; naturalA/S/remainingbinaryglobalcompatibilityOPEN, sharedpriorities/quarticPLANEretirementretained",
                       "Retained PF125: actual binary affine longitude fixes every critical pair; g=x+(l+1)u+m(x)v, g²=q_complement^m(u), g^4=1. Four residual orbit/filling profiles and affine critical shear counts, degree4 preferred longitude identity even with ONE ordinary cusp.177 cheap controls plus written parameter proof, no new family exclusion, no automatic integral intermediate coverage, no authors-code rerun. Correct marked linking/global compatibility or sharedHC4 EMPTY original full-sheet even radial profile remain open; full finite algebra/natural A/S/other binary affine compatibility OPEN, shared priorities/quartic PLANE retirement retained",
                       "Retained PF124: EVERYnontrivialfiniteintermediateM criticalC1classesindependentviaPF111originalThombasis/disjointpositivepullbackpackets=>c1K_MNONZEROQ. PRIORHC4MYKN1algebraicnontorsionandknownUADANVOL1/HC4MYLA1pseudo-planeexactcontrolsCREDITEDnotnew. Finitehdeg t/originalNpi1zero=>pi1Mfinite/orderdivt; resolutionchains/excision=>resolvednonzeroChern/b1Q0/b2>=critical+exceptionalcount. RetainedBODminimalresolutionSteindeformation+canonicalADEcontact/fillcomparison=>ALLINTERMEDIATEnegativeADE/S3ENDSCLOSED, singularMallowed. OSSINTEGRALboundaryc1retained; PVHM/p4onlyneedsb1Q0=>Lp1Msm/onecriticalcurve/b2one. BinaryfiniteEulerprofilechi=r(3-kappa), b2=3b/4-1immersive ORb/2-1onecusp independentdoublecontacts=>ALLb>=8/allimmersivebinaryLp1endclassesCLOSED, onlyb4ONEordinarycuspLp1intermediatecanremain. NoPicprimitive/paritytransfer/planeidentification.127controls~0.25sec67424KiB/importsretained/nofreshPDFaudit. NEXTsharedHC4specificEMPTYORIGINALfull-sheetfibrePF123 atdoubleCvalues(imprimitive/highere): excludingitwouldclosebinaryaffinequotientcontainerbutUNPROVED. NaturalAlt/Sym/otheraffineGLOBALcompatibility/invisibleoriginalimages/actualGLOBALfinitealgebraOPEN; no repeatedunmarkedcanonical/exactvolume/pseudo-controlcampaign. OriginalN/allprimes/sealedproofs/PF117closure/PF118closures/PF119profiles/PF120products/PF121classes/PF122allmapreduction/PF123descent/sharedpriorities/quarticretirement/closure-before-reviewpreserved",
                       "Retained PF123: EVERYintermediateLOCALmonodromyintransitive atALLfinitevalues, completeMfibreatleastTWO; transitivewouldforcealloriginalfactorssamevisiblebranch/uniqueCparameter=>d<=Eulerbudgetd-1. IMPORTPedersenv3finiteC2imagenormal=>quotient; liftSMALLuniformC2byHartogs/JacchainpullbackcriticalRONEsmoothaxis; IMPORTpublishedGiraldoRoederTheoremA=>Rsmooth, grouppreservesaxis+small=>CYCLIC. ALLintermediatecriticalcurvesDISJOINTSMOOTHA1/marked1/n(1,-e)/n<=e-1 viaverifiedLOCALPF106. BinaryaffineQ=>ENTIREM SMOOTH; PF119LOCALsimplecritargument=>b/4criticalA1isomorphicCnorm. OriginalN->M aboveE branchedonlysmoothaxis=>cyclic(s^k,z), originalVISIBLENsm/f1/e_i2k_i<=2t, eachpacketSUMk=t. OriginaletaC=(d/4)L/outsidecost=d/2-1-(d/4)L=>L0or1ONLY: immersedimageORONEordinarycusp. Higher/multiplecusps/visiblef>1/visiblequotientpointsCLOSEDforALLaffineblockquotients, evenimprimitiveoriginalG. Invisibleoriginalimages/singularitiesremain; outsidecostNOTonlyt0. Actualunboundedg_k hase2k/ell3k/eta k, wrongopenC*timesC; originalhigherEVENindicesNOTclosed.126controls~0.26sec68376KiB/twoexplicitgeometricimports/primaryHTMLpinned/noindependentproof/newPDFaudit. M smoothNOTA2recognition/quarticPLANEinversionnotintermediate4. NEXTnaturalAlt/Symquotients/binaryaffineGLOBALcompatibility/invisibleimages/actualGLOBALfinitealgebraOPEN; allprimes/originalN/sealedproofs/PF117closure/PF118closures/PF119profiles/PF120products/PF121classes/PF122allmapreduction/sharedpriorities/quarticretirement/closure-before-reviewpreserved",
                       "Retained PF122: hereditaryreductionforEVERYgenuinesheetBLOCKQUOTIENTofANYoriginalplaneKellerG(includingimprimitive). Exactt(b-suma)=1+t0+Eta+sum(w-ta),SUMa<=b-1/EVERYa<5b/6. QuotientweightedFULLbranchlocaltransitivity=>alloriginalfactorsoverorigincriticaloverBQ, injectiveweightedcomponents/qsumf=>fibred<=Eulerbudgetd-1 contradiction; nointermediateplane/smoothness/f1assumption. Thusirreduciblecontractiblequotientbranchclosed; TWOdisjointnonidBRANCHmeridiansforallQ=>muhalf/degrees2and3closed/subsetsclosed. Productbase>=6capclosed; base5swap4b/5+secondmin2b/5overfills, solehighmassimagecontractibleforbidden=>allproductsclosed. ClassicaltransvectionSUPPORTpointwisestabilizertrivial(Spformdiffsvperp,orthogvectorsspanV), soONEhalfmeridianGMtransvection+otherdisjointnonidcontradiction. ALLprimitiveblockquotients ONLYnaturalAlt/Sym orbinaryaffine. AffineQonebranch ORIGINALw=d/2/indicesEVENnotnecessarily2, originalotherbranchesmayQinvisible. PF119profilesONLYt1; quarticPLANEinversionnotintermediatedegree4.137controls~0.20sec65672KiB/CFSGimportsretained/GMoriginalPDFNOTaudited/noarbitraryquotientrepresentationsorcompositionfactorclassification. NEXTnatural/affinequotients/remainingimprimitivecompatibility/actualGLOBALfinitealgebraOPEN; originalN/allprimes/PF117closure/PF118closures/PF119profiles/PF120products/PF121classes/sharedpriorities/quarticretirement/closure-before-reviewpreserved",
                       "Retained PF121: universalDISJOINTnonidentityBRANCHmeridianpair forALLoriginalnormalplaneKellermaps, arbitraryindices/singularN. Distincttargetbranchesramifyindisjointlocalfactors; finiteMULTIBRANCHpointforcedbyPF113onecontractiblefullbranchclosure orChau+PalkaONEcusprectification/positivegenusshiftedfibre(nootherA1paramcurve). ALLinteriork-subsetactionsCLOSEDbyoriented-starletterseparation(noCFSG); EVERYnonaltprimitivebranchfixeshalfPF120 +explicitGMTheorem1ALL-elementinclusivehalftransvectionimport(primaryWEBTEXTread,originalPDFNOTretrieved/visuallyaudited,noOCRsubsetbound). Threebinaryclassicaltransvectionfamilieshaveoverlap2^(a-1)(2^(a-2)-1)>0ifcommuting; noncommutingcannotdisjoint=>ALLfamiliesCLOSED. PrimitiveONLYNATURALA_d/S_d orBINARYAFFINEV:H,Hirreduciblegeneratedbylineartransvections,PF119exacttwohalf-supportcuspprofilesremain.147controls~0.18sec66332KiB/nofinite-sampleproof/groupcampaign. NEXTnaturalAlt/Sym,binaryaffine,IMPRIMITIVE/actualGLOBALembeddedmultiplication/finitealgebraOPEN. Allprimes/originalN/PF117localtrace/PF118closures/PF119profiles/PF120products/sharedpriorities/quarticretirement/closure-before-reviewpreserved",
                       "Retained PF120: COMPLETEprimitiveCartesianPRODUCT-ACTIONclassCLOSEDALLdegrees: EVERYinertia<5d/6,involutions<2d/3; topcoordinateelementpure-swapinvolutionORfpr<=1/8,bothforbidden; NORMALtopquotientunramified/etaleA2trivial, contradictcoordinatepartition. CFSG-FREEexclusion, NOTimprimitivewreathaction. IMPORTEDBGTheorem4/Table2fullPDFp4-5audited=>primitiveONLYAlt/Symk-subsets,binaryaffine,threebinaryclassicalsingle-factorfamilies withtriality/U4(2).2/O6minus/A8isos. ALLprimegenericdegreesNATURALA_p/S_pnecessary, NOTinversion. Nonaltsocleprimitive3d/7<=mu<=d/2, atmostTWOdistinctFULLbranchimages, EVERYbranchinertia<=half byhighmasscuspcost36d/35contradiction.152controls/writtenunboundedproof/nogroupcampaign. NEXTsurvivingcontainers/IMPRIMITIVEcases/GLOBALembeddedmultiplication/finitealgebraOPEN; noautomaticintermediateplane/sourceactiondescent. Missing-pointtargetblowupchangesaffineboundaryhypotheses. PF119cuspcases/PF118closures/PF117localtrace/allprimes/originalN/sharedpriorities/quarticretirement/closure-before-review preserved",
                       "Retained PF119: ALL-simplee2actualnormalizationSMOOTH/f1, exactell=m_image+1/eta=m-1 fromcreditedPF105corank1, localS_(m+1)NOTglobal. Exactglobald-1=sumrj(2+Mj)+t0/pointwiseALLlengths. PF118involutoryhalf-supportd4rONEimage ONLYTWOcases:fullyimmersive(Eta0,t02r-1)orONEordinarycusp(Eta r,t0r-1); higher/multiplecuspsCLOSEDthere. Binary32:15unram/Eta0or7unram/Eta8; ENTIREfamilyNOTexcluded. Doublevalues2smoothbranches/contactarbitrary/sourceEMPTY; cuspS3profile3^r1^r, fixedsheetretentionnotassumed.107controls/writtenunboundedproof. KulikovTOTALramificationbirationalityfailsr>1, OrevkovA2EVERYpreimagebranchhypnotassumed, infinityuncontrolled/noautomaticcycliccomplement, HC4LW L4densityOPEN. NEXTentiretwohalf-supportcases/lower-supportactions/globalembeddedmultiplication/wholeplanefinitealgebraOPEN. PF118actionclosures/PF117closedlocaltrace/originalN/allprimes/sharedpriorities/quarticretirement/closure-before-review retained",
                       "Retained PF118: SOMEactualplaneinertiamoves<=half, fromsumdistinctbranchsupport+t0+Eta=d-1 andPF113FULLirreduciblebranchhalfbound. Allfaithfultransitiveactionsminimaldegree>halfCLOSEDPLANE, BOTHexceptional11PSL2(11)min8, ALLstatednaturalprojectivePGL_n(q)subgroupspoints/hyperplanes, ALLnaturalAGL_n(q)q>=3. ArbitrarydimensionalexceptionalquestionNOTclosed, nosemilineartransfer. BinaryAGL_n(2)equalityNOTexcluded: fullbranchONEC/transvectioninertia/d4ramA1e2f1/smoothN/t0+Eta=d2-1; Cnoncontractible/doublevaluesexist/atmost2smoothbranches/arbitrarycontact/UfibreEMPTYatdoublevalues/nounramprimeoverC. AGL5(2)32survives8ram/t0+Eta15.148controls/retainedKMFI95norerun/elementaryunboundedfixedspaceproof. NEXTlower-minimal-degreeactions/binaryprofile/generalarrangements/actualGLOBALembeddedmultiplication/wholeplanefinitealgebra OPEN. PF117localtracegateclosed; originalN/allprimes/sharedpriorities/quarticretirement/closure-before-review preserved",
                       "Retained PF117: LOCALsurface-trace restriction andactualCYCLICdescentCLOSED: keepNONREDUCEDfinitealgebraT, nilann=BTr_weighted andordinarycurve-codiffD_rB withALLgluing. ReynoldsdescendsinsideactualambientdualwithoutN-Gorensteinassumption, branchpermutationsincluded; reducedBdownGorensteincurve notnecessarilyplane. PF116actualresiduebridge R_B=Bvee_A*(B span g_C*R_C); downstairsinvariantquotientresiduesonly.79controls/writtenall-degreeproof/primaryStacksimports. AutomaticforsharpfinitecoverswrongC*timesC, NOnewglobalindexexclusion/strongerlocalconstant. NEXTGLOBALpolynomialmultiplication/trace-momentregularity/wholeplaneunits/Picard/integralboundarybasis/actualfinitealgebra OPEN. NoBernsteinorconnectiontransport; allprimes/originalfiniteN/sharedpriorities/quarticretirement/closure-before-review preserved",
                       "Retained PF116 (former localrestrictionOPEN is CLOSED by PF117): actual finite-plane logpullbackiso andramificationweightedcurve-residue/Jacobiandualmap proven; O_Bgluing/diagonalretained, NOTambienttrace-dualidentification. ActualhigherPuiseuxgermtau14/extravals11,15verified. GenuineUNBOUNDEDcyclicfiniteBelyifamilyattainsLOCALE5eta-2, primitiveordinaryresidualmaps/nonMorsecommonleading/NONGorenstein, butwrongopenC*timesC/unitrelation; globalWrightboundnotdeclaredoptimal. NEXTactualsurface trace-dual restrictionthroughnonreducedf(g) toBredwithconductor/Tor andWHOLEplane/globalmultiplication/finiteTARGETalgebra OPEN. NoBernsteinroottransport/newpreprintbinarytheorem/BMYequalityfromlogiso. KeeporiginalfiniteN/allprimes/sharedpriorities/quarticretirement/closure-before-review",
                       "Retained PF115: arbitrary singular chosenbranch/actual smooth-cyclicAXISsource pair gives exactresidualmu(H)=n(eta-q)(c-1),eta>=q,e>=4stricteta>=q+1. Complete minimalexcess numericalprofiles allordinarycusp, genuinefinitegerms retained. Actualnonprimitivef hasquotientq=f/singularresidualH andeta>=f+1; PF114alreadyimpliesallomittedf1throughdegree15,NOTinversion. Preserve nonlinear sourceP/arbitraryresidualbranchcount/actualordinarypreimages, no generalextensionofPF114linear5etamassbound. Next actualresidualdivisor/globalcovermultiplication/finiteTARGETalgebra remainsOPEN; no deeperindex/cablingcensus. Sharedtracepriorities,originalfiniteN/wholeplane/allprimes/quarticretirement/closure-before-review preserved",
                       "Retained PF114: local cyclic-axis toruscusp massw=e*q<=5eta-2delta with actual ordinarymeridian fillings and fullregularfibreorder; highindex e+ceil((e+2)/5)+r<=d-1 gives proportionalcapfloor((5d-17)/6) d>=9. Entire former e=d-4 twoRAMIFIED profile CLOSED. Sharedhighimage mass pays sumceil(f(e+2)/5); no scalarMorse, primitive-residue or fibre-identity assumption. Remaining noncontractibleimage arrangements/localglobalcompatibility/ACTUAL finite algebra OPEN. Do not substitute deeperindex/cabling census; preserve allprimes/wholeplane/originalfiniteN and sharedtracepriorities/quarticretirement/closure-before-review",
                       "Retained PF113 frontier (former e=d-4 OPEN is now CLOSED by PF114): original all-prime finite ledger localclusters<=d-r, COMPLETE finitefibres>=2, localmonodromy intransitive. FULL globallypositiveweightedbranch unions and irreduciblecontractibleFULLbranch excluded even with extraunram omitted images; no sourceactionlift or source-surjectivity/source-singleton claim. Universal e<=d-4 d>=7 closes olde=d-3 strip; newOPEN twoRAMIFIED profile(d-4,2),f1,oneexcess,distinctimages. Sharedhighmassw+b<=d-3; irreducibleFULLbranchmasshalfdegree, not totaloutsidebranch. Otherbrancharrangements/local-global gaps/ACTUAL finite algebra OPEN; preserve fullbranch/globalpositivity/wholeplane/allprime hypotheses and shared priorities/quartic retirement/closure-before-review",
                       "Retained PF112 frontier (former e=d-3 OPEN is CLOSED by PF113): same-image2w>d forces w+b<=d-2 (b=f1 prime count), irreducible WHOLE nonproper set TOTALmass<=floor(d/2), normalizedexcess>=d-1-floor(d/2), arbitraryresidues. Also all actual omitted e<=d-3 for d>=5, entire e=d-2 class CLOSED. Strict-high index has injective normalization, OWN excess>=1 and another DISTINCT nonproper image; two costs cannot merge. Irreducible WHOLE nonproper set =>half-degree index bound. Next e=d-3/d>=7 strip remains OPEN with exact ramified+unramified two-prime/one-excess profile; quotient n<=d-4, degree5 smooth/degree6 A1 only, not inversion. Candidate17 colored iterates CLOSED by broader index theorem, not old maximal-index gate. Other indices/arrangements/ACTUAL finite algebra OPEN; preserve whole source/all primes/classical imports/shared trace priorities/quartic retirement/closure-before-review",
                       "PRIORITY PF111: six-sheet family is already excluded by PF46/Orevkov maximal-index theorem; a consistency check, not a new exclusion. New theorem: all-normal integral Thom lattice and strict omitted end knots generate ALL H1. Complete coverage => exterior torsion-free/F complete; partial coverage hidden quotient costs outside primes. Six-sheet PSL2F5 trefoil and ALL r1/k0mod5 colored cable iterates CLOSED by hiddenZ2/Z4 versus zero Euler repair budget. Other covers/ends/common-leading germs/ACTUAL finite algebra OPEN; preserve complete-coverage and full-representation hypotheses, all source primes, PF110 genuine sources, shared trace priorities/quartic retirement/closure-before-review",
                       "PRIORITY PF110: actualPLANARcontactrequiresNONZEROINTEGRALboundaryc1;S3andallnegativeSU2quotientendsclosed. ALLL(p,1)singular/multipleprimepacketsCLOSED;survivingsmoothsoleA1,index2..p-1parity/f1. GenuineALLp>=3affinewholeplaneexactvolumesourcessaturateboundbutNOcompatiblefiniteTARGET. SurvivingLp1/otherends/commonleadinggerms/actualfinitealgebraOPEN;keeporientation/contactplanarity/integraltorsion/exceptionalp4rankproofexplicit;sharedpriorities/quarticretirementpreserved",
                       "PRIORITY PF109: completeALLNORMALnegativeADEendclass closed, includingsingularN/anydegree/primes/indices/residues. Freewhole-plane boundaryPicardbasis/rationalCherninjection=>resolvedc1nonzero; MINIMALresolutionBODSteindeformationpreservesc1/contact, ADEnegativefillforcesc1zero. Otherorientedends/commonleadingresiduegerms/ACTUALLaurentfinitealgebraOPEN;singularADEgapCLOSED. Nooriginalresolutionexactness orholomorphictargetpreservationclaim;keepunitrelationcountercontrol,actualendrecognition,sharedpriorities/quarticretirement/closure-before-review",
                       "PRIORITY PF108: completeALL-DEGREE actualSMOOTHnegativeADEendclass excluded, ANYprimecounts/indices, all-2negative-definiteTREEplumbingboundaries. KeepORIENTEDendrecognition/exactminimality anddeclareddeepcontact/fillingimports. SingularN/otherends/commonleadingresiduecases/f1packetsbeyondstatedclass/ACTUALLaurentfinitealgebraOPEN. NoQ/H1-onlyendrecognition;sharedpriorities/publishedquarticretirement/closure-before-reviewpreserved",
                       "PRIORITY PF107: coprime leading forms degenerate sector gives nf|e+1,scalarquotient1/n(1,1),e>=5nf-1,d>=6nf²+1>=12f²+1;f>1requiresd>=49. NOTautomatic. Remainingdegeneratepoints havecommonleadingzeros;f1packets/globalactualfinitealgebraOPEN. Genuine non-Gorenstein finite controls havewrongopenC*timesC;keepwhole-source/closure-before-review/publishedquarticretirement",
                       "PRIORITY PF106: unconditional all-index quotient type1/n(1,-e),coprime/orderboundn<=e-1, exactconductorintersection; scalarMorse CONDITIONAL nfdivides e-1 andd>=4f²+1. Extremalquotients forcef1; indices1,2,3,4,6 residue1atANYdegree andALL omittedf1throughgenericdegree10, notallsmall-degreeinversion. Remainingfullgoal includesf1packets, higherindexquotientresiduecovers withdegeneratequadraticpencils, otherinfinityactions andACTUAL Laurent finite algebra. NoautomaticMorse/generalnon-Morseresiduebound proved; keepwholeplane, countermodels, publishedquarticretirement andclosure-before-reviewpriority",
                       "Retained PF105: ALL actual omitted e<=3 primes have f1 at EVERY fielddegree for NORMAL N; singular triple points A1/localresidue1. All omitted f1 throughdegree8 (allodd through10), not allsmall-degree inversion. PF104 natural8 colored iterated family now closes for ALL NORMAL N, any extra primes/finite image singularities and every length/parameter, no endH1 assumption. Strong primary smooth residue consequence already usedHC4, attribution explicit. Remaining fullgoal includes f1packets, e>=4 cyclicquotient residuecovers, other colored infinity actions and ACTUAL Laurent finite algebra. Reviewdeferred, goalactive, publishedquarticretirement/sharedpriorities retained",
                       "Retained PF104: COMPLETE natural8 smooth and regular24 normal ACTUAL COLORED iterated families closed at ALL finite lengths/independent1mod6 cable parameters, with ANY extra unramified omitted primes. All-degree even-unimodular source/Euler budget and full-shear fibre-transitive integral transfer proved. Other representations, singular8, even indices/nontrivial endH1 and ACTUAL Laurent finite algebra remain OPEN. Full goal active; published quartic retirement/shared priorities and closure-before-review preserved",
                       "Retained PF103: ALL oriented A-type lens ends excluded at ALL fielddegrees and ANY omitted-prime count/index for smooth ACTUAL nonautomorphic Keller normalization. Whole-source E>0/contact theta plus exact minimality and primary filling rank bound close this complete end class. Other ends, singular corrections and actual Laurent finite algebra remain OPEN. Review deferred; published quartic retirement retained. Full goal active",
                       "Retained PF102 CLOSES actual smooth sole omitted A1/index3/fielddegree6 unknot and ALL torus infinity packets, arbitrary parameters and finite image singularities. Primitive source end generator, canonical parity and actual Stein theta remove every torus multiplier/base case; no actual topological covers are declared nonexistent. General all-degree theta=-1-(e-1)^2/h, oriented lens h>=e+1/noRP3, centre-free actual deck group and torus h|d also exclude odd-degree/odd-index packets in the stated smooth sole-line scope. NEXT COMPLETE GOAL: genuine ITERATED triple infinity knots; singular/multiple-prime lattices with proved correction terms; or global B=C[x,y] for ACTUAL normalization of C[P,Q],Jac constant. Those remain UNPROVED; all six-sheet maps are NOT excluded. Keep actual source, full peripheral shear/torsion and global finite algebra; use PF101 beyond trivial longitude without parameter campaigns. PF100 Laurent gap, countercontrols and published HC4MYKR0 quartic retirement retained. Closure before external-review expense; no outreach; full all-degree Wright goal OPEN",
                       "Higher-degree work must retain actual smooth normalization and whole source-chart hypotheses; no extension of the two-line terminal theorem to several unramified clusters is proved",
                       "For the retained rank-two residual cover use its proved z(lambda)=0 and split at infinity; do not assume trivial longitude on a general full higher-degree cover, for which PF101 uses joint orbits and shear",
                       "General lower-weight keys are not assumed binomial; higher-degree canonical quadrics require Gorensteinness, and general branch criteria require the explicit zero-excess/immersion hypotheses",
                       "Quartic Keller realizations are excluded by published HC4MYKR0; internal proof-chain review remains distinct. Higher-degree boundary exclusions remain open beyond the stated all-degree criteria; no control curve or infinity representation supplies a finite Keller map"],
            auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out.write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS {len(checks)} scoped adoption checks; research goal active; global navigation not certified")

