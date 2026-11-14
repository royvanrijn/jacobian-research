"""Scoped adoption audit; does not replace the repository-wide navigation checks."""
from pathlib import Path
import hashlib,json,re,sys,py_compile,tempfile

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
records={x["id"]:x for x in index["entries"] if x["id"] in {"PF2RULED1","PF2RULED2","PF2RULED3","PF2RULED4","PF2RULED5","PF2RULED6","PF2RULED7","PF2RULED8","PF2RULED9","PF2RULED10","PF2RULED11","PF2RULED12","PF2RULED13","PF2RULED14","PF2RULED15","PF2RULED16","PF2RULED17","PF2RULED18","PF2RULED19","PF2RULED20","PF2RULED21","PF2RULED22","PF2RULED23","PF2RULED24","PF2RULED25","PF2COORD1","PF2FACTOR1","PF2PRIM1","PF2PERIOD1","PF2WRIGHT1","PF2GRADE1","PF2TRACE1"}}
check("all thirty-two scoped claims registered",len(records)==32)
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
 "METHOD-EXPLICIT-COMPLEMENT-DEGREE-TRANSFER","METHOD-GLOBAL-VOLUME-EXACTNESS","METHOD-GRADING-DESCENT-TRACE-BOUNDARY","METHOD-WEIGHTED-LOCAL-NONCUBE","METHOD-FORMAL-CENTRALIZER-JETS","METHOD-BOUNDARY-LEADING-PERIOD","METHOD-FIRST-NEWTON-EDGE-POLE-BUDGET","METHOD-CUBE-CLUSTER-DISCRIMINANT-BUDGET","METHOD-CUSP-JET-LEADING-POWER","METHOD-DOUBLE-FILTRATION-CUSP-JETS","METHOD-PRIMITIVE-CHART-LEADING-FORM","METHOD-RESONANT-CHART-FRONTIER","METHOD-KUMMER-BOUNDARY-PERIOD","METHOD-GENERIC-POLAR-CUSP-BUDGET","METHOD-GENERIC-QUARTIC-POLAR-GAP","METHOD-GENERIC-CRITICAL-VALUE-UNIBRANCH","METHOD-COMPLETE-CUSP-JET-LIFTING","METHOD-HIGHEST-JACOBIAN-KERNEL-BLOCK","METHOD-AFFINE-JACOBIAN-FIBER-SURJECTIVITY","METHOD-UNIVERSAL-KERNEL-DETERMINANT-OBSTRUCTION"}]
validate_lessons(selected,index["entries"],check_reviews=True)
check("all twenty-five related lessons pass shared schema/path/scope validator",len(selected)==25)
work=json.loads((ROOT/"knowledge/work_items.json").read_text())["actions"]
action=next(x for x in work if x["id"]=="WORK-QUARTIC-RULED-COMPATIBILITY")
check("work gate preserves all thirty-two claims",set(records)<=set(action["claims"]))
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
                 "low_coordinate_receipt.json","primitive_cover_receipt.json","pole35_receipt.json","source_coordinates_receipt.json","volume_period_receipt.json","wright_grading_receipt.json","pole_three_all_height_receipt.json","pole_four_companions_receipt.json","pole_four_nonsquare_receipt.json","boundary_leading_period_receipt.json","pole58_receipt.json","pole68_cube_receipt.json","cusp_leading_receipt.json","contact_filtration_receipt.json","height_eight_receipt.json","height_nine_frontier_receipt.json","cubic_cover_period_receipt.json","polar_cusp_sector_receipt.json","quartic_polar_gap_receipt.json","unibranch_polar33_receipt.json","full_cusp_lift_receipt.json","cusp_high_jacobian_receipt.json","full_next_normal_receipt.json","second_normal_obstruction_receipt.json"):
    payload=json.loads((HERE/filename).read_text())
    check("receipt outcome "+filename,payload["outcome"]=="PASS")
    for path,digest in payload["source_sha256"].items():
        check("receipt seal "+filename+" "+path,
              hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    if "first_open_reduced_pairs" in payload and filename not in ("height_eight_receipt.json","height_nine_frontier_receipt.json","cubic_cover_period_receipt.json","polar_cusp_sector_receipt.json","quartic_polar_gap_receipt.json","unibranch_polar33_receipt.json","full_cusp_lift_receipt.json","cusp_high_jacobian_receipt.json","full_next_normal_receipt.json","second_normal_obstruction_receipt.json"):
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
check("higher poles and actual compatibility remain open",
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


report=dict(as_of="2026-11-14",focused_outcome="PASS",check_count=len(checks),checks=checks,
            mathematical_control_receipts=receipts,
            repository_navigation=dict(outcome="NOT_CERTIFIED_PREEXISTING_MISSING_FILES",
               known_first_error="KDSQ6: missing document consumer artifacts/generated-results/README.md",
               last_check_navigation_error="KDSQ6: missing document consumer artifacts/generated-results/README.md",
               last_render_navigation_error="KDSQ6: missing document consumer artifacts/generated-results/README.md",
               first_missing_path_still_absent=not(ROOT/"artifacts/generated-results/README.md").exists(),
               generated_status_and_indexes_manually_edited=False,
               missing_research_outputs_reconstructed=False),
            research_goal_completed=False,
            remaining=["generic P excluded by the global degree<=69 second-normal determinant in every marked noncontact family; solve Delta=0 with kappa nonzero, impose24 P-only conditions on the rank69 chart and retain all rank-drop strata, then u^3..u^9 with15 hidden directions; contact cases, gcd-one rows, field degree4, global finiteness and actual compatibility remain",
                       "all-height pole-four contact symbols and source component degree<=48; no automatic degree48 exclusion",
                       "genuinely mixed-weight finite quartic pair on the one marked zero-period Wright ring",
                       "extra Wright-specific normal-trace cancellation if branch Euler methods are used",
                       "higher pole pairs","actual finite-normalization compatibility"],
            auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out.write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS {len(checks)} scoped adoption checks; research goal active; global navigation not certified")

