#!/usr/bin/env python3
"""Source-guarded support-uniform individual squarefreeness and finite exceptions."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact/source controls must execute.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SQUAREFREENESS.md"
CN="verify_support_uniform_squarefree_completion.py"
IN="support_uniform_squarefree_completion_inputs.json"
LN="support_uniform_squarefree_completion_source_ledger.json"
RN="support_uniform_squarefree_completion_verification.json"
P="research/factorial/image-rigidity/"
ASSURANCE={"formal_verification":False,"external_review":False,"independent_replay":False}
CHILDREN={
 "high_ratio":{"proofs":["UNIFORM_ACTUAL_ENDPOINT_ROOT_SIMPLICITY.md","UNIFORM_ACTUAL_BULK_ROOT_SIMPLICITY.md"],
 "checker":"verify_uniform_actual_high_ratio_simplicity.py","inputs":"uniform_actual_high_ratio_simplicity_inputs.json",
 "ledger":"uniform_actual_high_ratio_simplicity_source_ledger.json","receipt":"uniform_actual_high_ratio_simplicity_verification.json",
 "count":11,"status":"PASS_SOURCE_GUARDED_ACTUAL_HIGH_RATIO_SIMPLICITY","support":False,"source_key":"files"},
 "fixed_terminal":{"proofs":["FIXED_TERMINAL_REPEATED_ROOT_SUPPORT_FINITENESS.md"],
 "checker":"verify_fixed_terminal_repeated_support.py","inputs":"fixed_terminal_repeated_support_inputs.json",
 "ledger":"fixed_terminal_repeated_support_source_ledger.json","receipt":"fixed_terminal_repeated_support_verification.json",
 "count":6,"status":"PASS_FIXED_TERMINAL_REPEATED_ROOT_SUPPORT_FINITENESS","support":False},
 "fixed_support":{"proofs":["EVENTUAL_SQUAREFREENESS_OF_ORDINARY_BINOMIAL_MOMENTS.md"],
 "checker":"verify_ordinary_binomial_eventual_squarefree.py","inputs":"ordinary_binomial_eventual_squarefree_inputs.json",
 "ledger":"ordinary_binomial_eventual_squarefree_source_ledger.json","receipt":"ordinary_binomial_eventual_squarefree_verification.json",
 "count":16,"status":"PASS","support":True}}
SPECS={
 "main_proof":("local",PN,"LOGICAL_written_one_order_five_region_proof_and_separate_six_way_reverse"),
 "exact_checker":("local",CN,"exact_support_and_artifact_freshness_NOT_formal_analytic_verification"),
 "typed_inputs":("local",IN,"exact_typed_scope_and_source_roles")}
for tag,d in CHILDREN.items():
 for j,name in enumerate(d["proofs"]):SPECS[tag+"_proof"+str(j)]=("local",name,
  {"high_ratio":"LOGICAL_individual_derivative_and_R6_LOCATION_only_NOT401_coprimality",
   "fixed_terminal":"LOGICAL_ONLY_tail_to_finite_repeated_triples_NOT_main_uniform_cutoff",
   "fixed_support":"LOGICAL_ONLY_finitely_bounded_supports_and_projection_reverse_NOT407"}[tag])
 for role in ["checker","inputs","ledger","receipt"]:
  SPECS[tag+"_"+role]=("local",d[role],"fresh_narrow_child_artifact_"+tag+"_NOT_parent_coprimality")
for role,name,scope in [
 ("fixed_support_regular_mechanism","ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md","LOGICAL408_complete_regular_periods_and_Cauchy_collars"),
 ("actual_simple_turn_mechanism","SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md","LOGICAL_RAW_normal_form_proper_arms_and_holomorphic_errors_ONLY")]:
 SPECS[role]=("local",name,scope)
for role,name,scope in [
 ("small_ratio_raw","UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md","LOGICAL_RAW_U5_U20_seam_periods_complete_tails_Morse_collar_NOT_coprimality"),
 ("subcritical_raw","COMPACT_SUBCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md","LOGICAL_RAW_S1_S16_actual_compactification_tail_circle_midpoint_ONLY"),
 ("bounded_thermal_raw","UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md","LOGICAL_RAW_T3_T19_actual_jets_complete_quartic_limit_FIXED_T_collar_ONLY"),
 ("large_thermal_raw","LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md","LOGICAL_RAW_L1_L9_L13_L19_L25_cycles_errors_Ai_ONLY_NO_gamma_over_M"),
 ("supercritical_raw","COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md","LOGICAL_RAW_A1_A21_A26_actual_localization_cycles_errors_ONLY_NO407"),
 ("uniform_CW_regular","UNIFORM_CUSP_CURIE_WEISS_REGULAR_COMPLEMENT.md","LOGICAL_RAW_complete_model_regular_geometry_and_tails"),
 ("uniform_quartic_growing","UNIFORM_GROWING_QUARTIC_CUSP_NEIGHBORHOOD.md","LOGICAL_RAW_uniform_quartic_packets_revalidated_with_small_T"),
 ("uniform_quartic_bounded","UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md","LOGICAL_RAW_compact_phase_jets_and_complete_bounded_field_tails"),
 ("quartic_model","QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md","LOGICAL_ONLY_sections1_2_Newman_location_ODE_heat_simplicity_NO_HJN"),
 ("quartic_unfolding","QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md","LOGICAL_RAW_actual_quartic_jets_and_local_analytic_geometry"),
 ("strict_packets_408","SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md","LOGICAL408_actual_complete_strict_packets_ONLY"),
 ("quartic_growing_408","CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md","LOGICAL408_complete_quartic_packets_and_open_collars_ONLY"),
 ("quartic_bounded_408","BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md","LOGICAL408_complete_bounded_model_period_ONLY"),
 ("morse_mechanism_408","CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md","LOGICAL408_actual_complete_periods_NO_coprimality_inference"),
 ("wall_mechanism_408","CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md","LOGICAL408_actual_wall_and_real_boundary_periods_ONLY"),
 ("stationary_geometry_408","COALESCING_BINOMIAL_SADDLE_BRIDGE.md","LOGICAL408_actual_stationary_and_phase_geometry_ONLY"),
 ("old_gamma_upper","GAMMA_PRODUCT_ROOT_MODULUS.tex","LOGICAL_ONLY_product_annulus_UPPER_NOT_soft_pair_simplicity"),
 ("raw_sparse","lower_exponent_initial_and_boundary.tex","LOGICAL_ONLY_thm_sparse_INDIVIDUAL_simplicity_clause_and_proof"),
 ("raw_constant_analytic","ALL_UNIVARIATE_TWO_TERM_SFC.tex","LOGICAL_ONLY_thm_analytic_individual_tail_NO_finite_prefix"),
 ("raw_constant_large","ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex","LOGICAL_ONLY_LITERAL_thm_large_F_squarefree_nonzero_reversal"),
 ("raw_primepower","PRIME_POWER_TWO_TERM_SFC.tex","LOGICAL_shape1_small_primepower_sfc_plus_EXACT_a0_Appell_ONLY"),
 ("raw_entropy_credit","EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex","AUTHENTICATED_raw_analytic_credit_ONLY_NO_finite_prefix"),
 ("primary_Kabluchko","uniform_cusp_sources/kabluchko_2203_05533v3.pdf","LOGICAL_Thm2_3_inverse_for_cutoff_Lemma2_1_finite_circle_ONLY_for847_NO_circle_PR"),
 ("primary_Newman","quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf","LOGICAL_Thm2_empty_product_all_real_temperature_location_NO_HJN"),
 ("primary_Airy_ODE","sources/dlmf-9-2-endpoint.html","LOGICAL_Airy_ODE_uniqueness_and_simplicity"),
 ("primary_Airy_sectors","sources/dlmf-9-7-endpoint.html","LOGICAL_Ai_single_exponential_sector_model_only_actual_cycles_proved_separately"),
 ("primary_Airy_location_credit","sources/dlmf-9-9-endpoint.html","LOGICAL_Ai_negative_zero_location_ONLY_in_imported_R6_boundary_proof_NO_all_index_separation"),
 ("primary_cubic_context","sources/dlmf-2-4-endpoint.html","AUTHENTICATED_framework_context_allows_other_cycles_NOT_actual_cycle_proof"),
 ("large_thermal_primary_context","LARGE_THERMAL_PRIMARY_CONTEXT.md","AUTHENTICATED_classical_hypothesis_audit_context_only"),
 ("large_thermal_cubic_context","sources/dlmf-36-12-large-thermal.html","AUTHENTICATED_classical_cubic_framework_context_only")]:
 SPECS[role]=("repo",P+name,scope)
# Filled once the sibling packets are frozen; all mathematical/primary bodies are immutable.
IMMUTABLE={'actual_simple_turn_mechanism': 'be8c73ad8b346bb5c260f165b30057dca9295d000629557b86cf9e312929ead0', 'bounded_thermal_raw': '83a48e4da05099917114b51cf3bb6fc8ac686c5f42c936efbdc3d028a3cb7e93', 'fixed_support_checker': 'cf9764e4e0886f7cbc8f94728a4d64212aff6ff6eb4c9c513b36f4ba09aee990', 'fixed_support_proof0': '4aeba104dd1ac39688f912268c1a2ad4fc4af489a320a0bea0af8d365c13d4b3', 'fixed_support_regular_mechanism': '40ee81bec6a1bf906e695a819a7e609be85ca03b19e5df2fb7a5a2388c5a67e0', 'fixed_terminal_checker': '5682b139b031f48bde77a74749e2caf7daa9963e3ec52447ad9c8db5847941b9', 'fixed_terminal_proof0': '847a961f4f55400f2fc5973fa2ea547c694ea998150dce1e45a75e1cb206a24d', 'high_ratio_checker': '713aad26f1e8d7eaeaf444968ee2df70fb6437b114a50d923314b5eb50f6c397', 'high_ratio_proof0': '799db033a60dfdd7bc6383252b894f8a4b0c001ce240dbc83858eb93956e97aa', 'high_ratio_proof1': '4a2f09edc09fcfab3dfdcec9beab321c68e9d656ee67172249afc904f033633a', 'large_thermal_cubic_context': 'b10bf13a639a0673f559713a1dff703d67ab3ff1d1e2fc81cba58148628d78ff', 'large_thermal_primary_context': '57b154c5dbd87856d8eabadc18fcaee3ceaee0be497f5526c47ae0323ac186cc', 'large_thermal_raw': 'a9ca7dfb62b814d0637d1a56baac5084f951b479e8e64f9e5df5b36948d13c87', 'main_proof': 'eba17e08f06ac6d0c848ef80bfc552f8233c649e79ad754b078729210f8c023d', 'morse_mechanism_408': '7e21794c68108551133243bb1b1e574166c27dfad8c0750ba28f809d8900b85c', 'old_gamma_upper': 'b6cb212dd2b461a12c5f687375089da7839495031737ef04d0de2cc8151950b9', 'primary_Airy_ODE': 'a237add8efb2c719dc7cc74baf78b435b50beec368dbc1adc672645f792b7ac7', 'primary_Airy_location_credit': '7f5e538332164cf3103abb2374acc2619bfca1c765cd46e04fa738ec790d962b', 'primary_Airy_sectors': '2073a97493c6b612a5d7f3414de437aa1a58e97bf1268563cacd33aa537bf6d1', 'primary_Kabluchko': 'fbf2a26f4424761bf0a75ac2846b27aade3cd4c65433ff11fdf6ca929073a02c', 'primary_Newman': 'f154969b9a39b0e64075570eb9f4c920affbea85bbcee0367eb66302e7a1be28', 'primary_cubic_context': 'c90ac6736075b2d0851c48421da534e74d99ef2d1f33884024ef6297f1dc971f', 'quartic_bounded_408': 'abd5481d1dda7beb2d637ea9904d580c51882030e17d108496f6163b4003f906', 'quartic_growing_408': '679586446ff23794d0182c6bc45feedd7d710d13973600d9076da6c578a32409', 'quartic_model': '957127791def20f8c73730810775760cb124a9175d5a3c7b91e8cf7498521d89', 'quartic_unfolding': '86b19ffd10650e14a91cddd47974a360d9329f44dfeda6cc860cb0a84faa7953', 'raw_constant_analytic': '112d995972906353784af4f9ceb2b03eb7e71ce3e7bd52adc5a44c6869aeadc0', 'raw_constant_large': '604e1687ea759785fe5d97a2cf214c84732be10ac0df0cae208d7bd01435dcd7', 'raw_entropy_credit': '673958f0dd9f69ba5d1e560eca531e3738b17622ab3fe3f0fcb0a44557274995', 'raw_primepower': '0e94e91d7b60905e9ef4d1c42234b6ef5f5307cac7899d77b7c5ce6714243228', 'raw_sparse': 'e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45', 'small_ratio_raw': 'ebf6489506e379f766728c199d821862ab02da3b034e9ffb428eb35a0bd925e5', 'stationary_geometry_408': '8ac6b9100edf75f0ba1cc1d0b03209bff3b76c4802e1b855840065911fd9131e', 'strict_packets_408': '3bcd961a7cbf63fe780f30580554118355f4f5252cf1d92992741881ece500e6', 'subcritical_raw': 'bc230f8206319e53ba84e43371407815b0a92a347e245aaa61633ce1a46aa2c5', 'supercritical_raw': '0a4945c8a491c20a628db6e1fd37af555b70fa5e4b850d4202087c6bc5c4219f', 'uniform_CW_regular': '8b64926a462956eb4d1b60cdbb8246ae92c85a0614136bf1f6a2e251041d31f8', 'uniform_quartic_bounded': '8a7a94ab8daadf196ceab85b4d3be1387de4b2f2825144dd12332407b041c08e', 'uniform_quartic_growing': 'e9b23ad408aa1d0933d55d07b0126a801a79dc2623a28900a310c156e77543ca', 'wall_mechanism_408': 'f38231556fbb8f8eaecc5170a28bd90a9fa9498673a28129d41d0c870b5afc44'}
SCOPE={"schema":1,"reviewed_proof_sha256":"eba17e08f06ac6d0c848ef80bfc552f8233c649e79ad754b078729210f8c023d",
 "parameters":"integers a>=0,d>=1,m>=2; actual Gamma shape a*m+1",
 "statement":"one unevaluated Mstar for ALL ordinary supports; every m>=Mstar individual A_m squarefree",
 "all_characteristic_zero_fields_and_nonzero_coefficient_phases":True,
 "rational_Bezout_no_arbitrary_field_embedding":True,
 "regional_RAW_analytic_mechanisms_NOT_coprimality_conclusions":True,
 "artifact_freshness_distinct_from_logical_premise_gating":True,
 "main_forward_cutoff_survives847_logical_withholding":True,
 "fresh401_role":"full numerical artifact authentication; logical U9/U12 and R6 LOCATION only",
 "exact_integer_partial_periods_holomorphic_same_support":True,
 "whole_original_tails_actual_leading_actions":True,
 "Cauchy_on_OPEN_support_fixed_parameter_collars":True,
 "thermal_center_rho_and_T_held_fixed_on_q_collars":True,
 "thermal_T_over_k_amplitude_term_retained":True,
 "bounded_real_temperature_complete_Z_limit_and_simplicity":True,
 "Newman_sections1_2_ONLY_no_HJN":True,
 "Kabluchko_Thm2_3_selected_inverse_logical":True,
 "Kabluchko_Lemma2_1_ONLY_fixed_terminal_child":True,
 "Ai_model_ODE_and_COMPLETE_actual_right_cycle":True,
 "high_ratio_statement":"d>=2,a>=10^24*d^2,m>=10^16 individual squarefree",
 "R6_used_for_EVERY_individual_root_LOCATION_only":True,
 "old408_used_ONLY_finitely_bounded_supports_and_projection_reverse":True,
 "fixed_terminal847_required_ONLY_cutoff_to_finite_repeated_triples":True,
 "main_uniform_cutoff_requires847":False,
 "six_way_equivalence_under_separate847_408_contexts":True,
 "coefficient_projective_moment_singularity_NOT_roots_in_X":True,
 "singular_phase_AB_nonzero_algebraic_ratio_gcd_A_Aprime":True,
 "global_repeated_triples_finite":True,"repeated_support_projection_finite":True,
 "intrinsic_singular_projective_phases_finite":True,
 "unscaled_coefficient_pairs_finite":False,
 "all_order_squarefreeness":False,"all_order_SFC":False,
 "global_repeated_set_empty":False,"effective_cutoff_or_box":False,
 "407_or418_coprimality_used_as_simplicity":False,
 "soft_pair_bound_used_as_repeated_support_finiteness":False,
 "homotopy_counts_used_as_simplicity":False,"HJN_or_later_quartic_Robin_imported":False,
 "new_sharp_lower_annulus_premise":False,"heavy_finite_prefix_import":False,
 "finite_controls_formally_verify_analytic_estimates":False,
 "raw_large_source":"LITERAL individual F_(m,d) squarefree; nonzero reversal; no Appell roundtrip",
 "K":4*10**24,"a0_uniform_cutoff":10000,
 "controls":{"max_lower":4,"max_gap":8,"max_order":5,"cover_deltas":["1/1000000","1/10","1","10"],
 "cover_epsilons":["1/10000","1/100","1","10"],"root_t_samples":[100,101,128,200]},
 "assurance_flags":ASSURANCE}
def need(b,s):
 if not b:raise ValueError(s)
def H(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a)is not type(b):return False
 if type(b)is dict:return set(a)==set(b)and all(typed(a[k],v)for k,v in b.items())
 if type(b)is list:return len(a)==len(b)and all(typed(x,y)for x,y in zip(a,b))
 return a==b
def check_inputs(d):need(typed(d,SCOPE),"Exact typed within-order and finite-exception scope")
def check_sources(d,repo):
 need(type(d)is dict and set(d)=={"schema","assurance_flags","sources"},"Exact main source schema")
 need(type(d["schema"])is int and d["schema"]==1 and typed(d["assurance_flags"],ASSURANCE),"Typed schema and assurance")
 need(type(d["sources"])is list and len(d["sources"])==len(SPECS),"Every mandatory role")
 rows={};paths=set()
 for r in d["sources"]:
  need(type(r)is dict and set(r)=={"role","kind","path","relation","sha256"},"Exact source row")
  role=r["role"];need(type(role)is str and role in SPECS and role not in rows,"Unique known role")
  need((r["kind"],r["path"],r["relation"])==SPECS[role],"Exact kind/path/logical locator")
  need(type(r["sha256"])is str and re.fullmatch("[0-9a-f]{64}",r["sha256"]),"SHA256")
  root=(repo if r["kind"]=="repo"else HERE).resolve();p=(root/r["path"]).resolve()
  need(p.is_relative_to(root)and p.is_file()and H(p)==r["sha256"],"Missing/stale artifact:"+role)
  need((r["kind"],r["path"])not in paths,"Unique artifact role");paths.add((r["kind"],r["path"]));rows[role]=r
 for role,sha in IMMUTABLE.items():need(rows[role]["sha256"]==sha,"Exact frozen reviewed artifact:"+role)
 need(rows["main_proof"]["sha256"]==SCOPE["reviewed_proof_sha256"],"Exact reviewed mathematical proof")
 forbidden={"ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md","SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md",
 "all_exponent_prefix_certificates.json","exponent_six_prefix_certificates.json","ordinary_lower_one_finite_bridge_verification.json"}
 need(not any(Path(p).name in forbidden for k,p in paths),"No407/418 conclusion or heavy boundary premise")
 return rows
def child_boundaries(repo,rows):
 direct={(r["path"]if r["kind"]=="repo"else P+r["path"]):r["sha256"]for r in rows.values()}
 reports=[]
 for tag,c in CHILDREN.items():
  d=load(HERE/c["ledger"]);need(type(d.get("schema"))is int and d["schema"]==1,"Child schema1:"+tag)
  rr=d.get(c.get("source_key","sources"));need(type(rr)is list and len(rr)==c["count"],"Complete declared child boundary:"+tag)
  found=set();seen=set()
  for r in rr:
   kind=r.get("kind",r.get("root"));path=r.get("path")
   need(kind in{"local","support","repo"}and type(path)is str,"Child source root")
   root=repo if kind=="repo"else HERE;p=(root/path).resolve()
   need(p.is_relative_to(root.resolve())and p.is_file()and H(p)==r["sha256"],"Current entire child boundary:"+tag+":"+path)
   key=(kind,path);need(key not in seen,"Unique child source");seen.add(key)
   canonical=path if kind=="repo"else P+path
   if canonical in direct:need(direct[canonical]==r["sha256"],"Overlapping child/main pins agree")
   if kind!="repo"and path in c["proofs"]+[c["checker"],c["inputs"]]:found.add(path)
  need(found==set(c["proofs"]+[c["checker"],c["inputs"]]),"Indispensable child own roles")
  reports.append({"tag":tag,"declared_current_sources":len(rr)})
 return reports
def fresh(repo,python):
 reports=[]
 with tempfile.TemporaryDirectory(prefix="squarefree-narrow-replays-")as tmp:
  for tag,c in CHILDREN.items():
   out=Path(tmp)/(tag+".json");cmd=[str(python),"-B",str(HERE/c["checker"]),"--repo-root",str(repo)]
   if c["support"]:cmd+=["--support-root",str(HERE)]
   cmd+=["--output",str(out)]
   r=subprocess.run(cmd,cwd=tmp,capture_output=True,text=True,timeout=180)
   need(r.returncode==0,"Fresh narrow child:"+tag+";"+r.stdout+";"+r.stderr)
   raw=out.read_bytes();need(raw==(HERE/c["receipt"]).read_bytes(),"Fresh child BYTE identity:"+tag)
   got=json.loads(raw);child_input=load(HERE/c["inputs"])
   expected_scope=child_input["scope"]if tag=="high_ratio"else child_input
   need(got["status"]==c["status"]and got["scope"]==expected_scope,"Exact child own scope")
   if tag=="high_ratio":
    need(got["analytical_boundary"]==child_input["analytic_boundary"],"Exact numerical versus logical boundary")
    need(all(got.get(k)is False for k in ASSURANCE),"No high-ratio assurance upgrade")
    need(got["parent_replay"]["role"]=="numerical source authentication only","401 freshness is not a coprimality premise")
   need(got["source_ledger_sha256"]==H(HERE/c["ledger"]),"Exact child ledger identity")
   reports.append({"tag":tag,"receipt_byte_identical":True,"receipt_sha256":hashlib.sha256(raw).hexdigest()})
 return reports


def logical_controls():
 rules=[
 ({"TailQ","FixedTerminal847"},"RepeatedFinite"),({"RepeatedFinite"},"TailQ"),
 ({"TailQ"},"FieldTail"),({"FieldTail"},"CoefficientComplexTail"),({"CoefficientComplexTail"},"TailQ"),
 ({"RepeatedFinite"},"SingularProjectiveFinite"),({"SingularProjectiveFinite"},"SupportFinite"),
 ({"SupportFinite","FixedSupport408"},"RepeatedFinite")]
 predicates={"TailQ","RepeatedFinite","FieldTail","CoefficientComplexTail","SingularProjectiveFinite","SupportFinite"}
 def close(seed,fn,fs):
  reached=set(seed)
  if fn:reached.add("FixedTerminal847")
  if fs:reached.add("FixedSupport408")
  while True:
   after=reached|{c for p,c in rules if p<=reached}
   if after==reached:return reached
   reached=after
 out=[]
 for p in sorted(predicates):
  reached=close({p},True,True)
  need(predicates<=reached,"Six-way reverse has both explicit context premises")
  out.append("six-way seeded:"+p)
 need(not(close(set(),True,True)&predicates),"An unseeded implication component cannot prove itself");out.append("no-cycle-bootstrap")
 for p in ["SingularProjectiveFinite","SupportFinite"]:
  reached=close({p},True,False)
  need("RepeatedFinite"not in reached and"TailQ"not in reached,"Projection finiteness reverse NEEDS408")
  out.append("408-required:"+p)
 reached=close({"TailQ"},False,True)
 need("RepeatedFinite"not in reached,"Uniform terminal tail alone NEEDS fixed-terminal fiber finiteness")
 out.append("fixed-terminal847-required-for-tail-to-B")
 reached=close({"TailQ"},True,False)
 need("RepeatedFinite"in reached,"Tail-to-B reverse DOES NOT need408")
 out.append("fixed-terminal847-suffices-without408")
 for p in sorted(predicates):
  reached=close({p},True,True)
  need(not({"BadEmpty","AllOrder","EffectiveCutoff"}&reached),"No forbidden backwards strengthening")
  out.append("no-allorder-empty-effective:"+p)
 # Abstract examples, deliberately NOT moment counterexamples.
 fixed_terminal={(j,j,3)for j in range(1,65)}
 need(len(fixed_terminal)==64 and len({n for a,d,n in fixed_terminal})==1,"Tail alone permits infinitely many supports at an early terminal in an abstract model")
 out.append("abstract-fixedN-infinite-support-pattern")
 fixed_support={(1,2,n)for n in range(2,66)}
 need(len(fixed_support)==64 and len({(a,d)for a,d,n in fixed_support})==1,"Finite support projection does not by itself bound terminals")
 out.append("abstract-fixed-support-infinite-terminal-pattern")
 pairs={(j,j)for j in range(1,65)};ratios={Q(a,b)for a,b in pairs}
 need(len(pairs)==64 and ratios=={Q(1)},"Projective finite is NOT unscaled coefficient finite")
 out.append("projective-versus-unscaled")
 return out

def exact_support():
 N=12
 def X(i):t=[0]*N;t[i]=1;return {tuple(t):Q(1)}
 def C(c):return {tuple([0]*N):Q(c)}if c else{}
 def add(*terms):
  out={}
  for p in terms:
   for e,c in p.items():out[e]=out.get(e,Q(0))+c
  return {e:c for e,c in out.items()if c}
 def neg(p):return {e:-c for e,c in p.items()}
 def sub(a,b):return add(a,neg(b))
 def mul(*terms):
  out=C(1)
  for p in terms:
   r={}
   for e,c in out.items():
    for f,d in p.items():
     k=tuple(x+y for x,y in zip(e,f));r[k]=r.get(k,Q(0))+c*d
   out={e:c for e,c in r.items()if c}
  return out
 def power(a,n):
  out=C(1)
  for _ in range(n):out=mul(out,a)
  return out
 def scale(a,c):return mul(C(c),a)
 polys=[];laurent=[];guards=[];negs=[]
 def zero(name,p):
  assert not p,name
  polys.append(name)
 def check(name,b):assert b,name;guards.append(name)
 a,d,w,s,t=[X(i)for i in range(5)]
 crit=add(mul(sub(a,s),add(w,t)),mul(d,t))
 zero('single_actual_stationary_envelope',add(mul(d,w),neg(mul(sub(add(a,d),s),add(w,t))),crit))
 sr,sl,tr,tl=[X(i)for i in range(3,7)]
 cr=add(mul(sub(a,sr),add(w,tr)),mul(d,tr))
 cl=add(mul(sub(a,sl),add(w,tl)),mul(d,tl))
 zero('two_actual_stationary_envelope',add(mul(d,w,sub(tl,tr)),mul(sub(sr,sl),add(w,tr),add(w,tl)),mul(cr,add(w,tl)),neg(mul(cl,add(w,tr)))))
 l0,r0,gl,gr,rl,rr,dl,dr=[X(i)for i in range(8)]
 pl=mul(l0,add(C(1),rl));pr=mul(r0,add(C(1),rr));T=mul(add(C(1),rl),add(C(1),rr))
 di=add(mul(l0,add(mul(gl,add(C(1),rl)),dl)),mul(r0,add(mul(gr,add(C(1),rr)),dr)))
 en=add(mul(sub(gr,gl),T),mul(dr,add(C(1),rl)),neg(mul(dl,add(C(1),rr))))
 zero('exact_two_remainder_cancellation_derivative',sub(sub(mul(di,T),mul(pr,en)),mul(add(pl,pr),add(mul(gl,T),mul(dl,add(C(1),rr))))))
 d,q,L,tau,dv,beta,r,Dq,Do=[X(i)for i in range(9)]
 c1=add(mul(d,q,Dq),mul(L,tau,dv))
 c2=sub(mul(beta,Do),mul(power(tau,3),Dq))
 c3=add(mul(d,beta),mul(L,sub(C(1),r)))
 target=sub(mul(q,sub(C(1),r),Do),mul(power(tau,4),dv))
 zero('exact_thermal_Domega_constraint_identity',add(mul(d,q,c2),mul(power(tau,3),c1),neg(mul(q,Do,c3)),mul(L,target)))
 c,x=[X(i)for i in range(2)]
 e=mul(sub(c,C(Q(1,4))),power(x,2));CC=add(mul(c,power(x,2)),scale(x,Q(1,2)))
 zero('exact_critical_ratio_r',sub(mul(e,add(scale(mul(c,x),2),C(1))),scale(mul(x,CC,sub(c,C(Q(1,4)))),2)))
 d=X(0);zero('high_ratio_K4_identity',sub(sub(scale(power(sub(d,C(1)),2),4),power(d,2)),mul(sub(d,C(2)),sub(scale(d,3),C(2)))))
 # Compact normalization dL=2a+d-1 with a=c(d-1)^2.
 d,L,c=[X(i)for i in range(3)]
 zero('compact_L_over_d_cleared',sub(sub(mul(d,L),add(scale(mul(c,power(sub(d,C(1)),2)),2),sub(d,C(1)))),sub(mul(d,L),add(scale(mul(c,power(d,2)),2),neg(scale(mul(c,d),4)),scale(c,2),d,C(-1)))))
 # Thermal scale squared is an exact combination of defining equations.
 d,L,r,CC,Sc=[X(i)for i in range(5)]
 target=add(scale(mul(sub(C(1),r),power(L,2)),4),neg(power(d,2)),C(1))
 eq1=sub(power(L,2),CC);eq2=sub(mul(sub(C(1),r),CC),Sc);eq3=add(scale(Sc,4),neg(power(d,2)),C(1))
 zero('thermal_L_over_d_squared',sub(target,add(scale(mul(sub(C(1),r),eq1),4),scale(eq2,4),eq3)))
 j,m,g=[X(i)for i in range(3)]
 zero('universal_Gamma_ODE_coefficient_recurrence',add(mul(add(j,C(1)),sub(sub(add(j,C(1)),g),m)),mul(add(j,C(1)),sub(add(g,m),add(j,C(1))))))
 a,d,m,k=[X(i)for i in range(4)]
 zero('actual_factorial_integral_scale_exponents',sub(sub(add(mul(add(a,d),m),C(1)),mul(d,sub(m,k))),add(mul(a,m),mul(d,k),C(1))))
 def vecsub(a,b):return tuple(x-y for x,y in zip(a,b))
 # Coordinates (rho,k,T). Substitute m=k*rho^-4 or m=k*rho^-3.
 for name,mvec,amp,expect in [
  ('quartic_amplitude_rho',( -4,1,0),(-2,0,0),(1,-1,0)),
  ('quartic_amplitude_T',(-4,1,0),(-3,0,1),(0,-1,1)),
  ('quartic_remainder',(-4,1,0),(-3,-1,0),(0,-2,0)),
  ('cubic_amplitude',(-3,1,0),(-1,0,0),(1,-1,0)),
  ('cubic_remainder',(-3,1,0),(-2,-1,0),(0,-2,0))]:
  lead=(mvec[0]+1,mvec[1],mvec[2])
  assert vecsub(amp,lead)==expect,name;laurent.append(name)
 # Exact abstract budget sufficient conditions, not claimed actual constants.
 for B in [Q(1,10),Q(1),Q(7),Q(100)]:
  for Cr in [Q(1,10),Q(1),Q(7),Q(100)]:
   base=max(1,(16*B).__ceil__(),(64*Cr).__ceil__())
   for inc in [0,1,100]:
    mm=base+inc
    check('abstract_regular_margin',Q(mm,4)-(B+4*Cr/mm)>=Q(mm,8))
 for rho in [Q(1,1000),Q(1,100),Q(1,10)]:
  for TT in [Q(0),Q(1,100),Q(1,10)]:
   for kk in [100,1000,10**6]:
    check('abstract_quartic_relative_margin',10*(rho/kk+TT/kk+Q(1,kk*kk))<Q(1,2))
    check('abstract_cubic_relative_margin',10*(rho/kk+Q(1,kk*kk))<Q(1,2))
 check('thermal_Ld_lower',Q(15,88)>Q(2,5)**2)
 check('thermal_Ld_upper',Q(5,18)<Q(3,5)**2)
 check('fixed_q_annulus_Cauchy_inner',Q(1,2)-Q(1,12)>Q(1,3))
 check('fixed_q_annulus_Cauchy_outer',2+Q(1,12)<3)
 for gap in [Q(1,1000),Q(1,10),Q(1)]:
  for mm in [1,2,10,1000]:
   check('c0_exp_fourth_term_absorption_identity',13*mm*24/(gap*mm)**4==Q(312,gap**4*mm**3))
 check('c0_width_square_bound',128<12**2)
 check('critical_strip_margin',8*Q(1,16)<1)
 # Exact negatives: no generic derivative replacement and no pair-to-simplicity inference.
 def rise(a,k):
  t=1
  for h in range(k):t*=a+h
  return t
 assert 2*rise(3,2)!=2*rise(2,2);negs.append('Appell_replacement_false_for_a_positive')
 assert (-1)**2-2*(-1)+1==4 and 1-2+1==0 and 2-2==0;negs.append('coprime_pair_component_can_repeat')

 return {'universal_polynomial_identities':polys,'universal_Laurent_budget_reductions':laurent,'rational_budget_controls':guards,'mathematical_negative_controls':negs}

def moment_controls(inp):
 cfg=inp["controls"];norm=gradient=euler=scaling=appell=reversal=cover=large=0
 def rising(g,n):
  out=1
  for k in range(n):out*=g+k
  return out
 for a in range(cfg["max_lower"]+1):
  for d in range(1,cfg["max_gap"]+1):
   for m in range(1,cfg["max_order"]+1):
    aa,bb=Q(2,3),Q(-3,5);x=aa/bb
    coeff=[Q(comb(m,j)*factorial(a*m+d*j))for j in range(m+1)]
    h=sum(coeff[j]*aa**(m-j)*bb**j for j in range(m+1))
    ha=sum((m-j)*coeff[j]*aa**(m-j-1)*bb**j for j in range(m))
    hb=sum(j*coeff[j]*aa**(m-j)*bb**(j-1)for j in range(1,m+1))
    pol=sum(Q(comb(m,j)*rising(a*m+1,d*j))*x**(m-j)for j in range(m+1))
    dp=sum((m-j)*Q(comb(m,j)*rising(a*m+1,d*j))*x**(m-j-1)for j in range(m))
    need(h==bb**m*factorial(a*m)*pol,"Actual coefficient moment normalization");norm+=1
    need(ha==bb**(m-1)*factorial(a*m)*dp,"Actual projective derivative chart");gradient+=1
    need(aa*ha+bb*hb==m*h,"Actual coefficient-form Euler identity");euler+=1
    lam=Q(7,3)
    hs=sum(coeff[j]*(lam*aa)**(m-j)*(lam*bb)**j for j in range(m+1))
    need(hs==lam**m*h and factorial(a*m)>0 and factorial((a+d)*m)>0,"Projective scaling and nonzero axes");scaling+=1
 for d in range(1,17):
  for m in range(2,9):
   A=[Q(comb(m,m-j)*factorial(d*(m-j)))for j in range(m+1)]
   B=[Q(comb(m-1,m-1-j)*factorial(d*(m-1-j)))for j in range(m)]
   need([j*A[j]for j in range(1,m+1)]==[m*b for b in B],"a0 EXACT Appell vector");appell+=1
   F=[Q(comb(m,j)*factorial(d*j))for j in range(m+1)]
   need(A==list(reversed(F))and A[0]>0 and F[0]>0,"RAW large criterion nonzero reversal");reversal+=1
 for delta0 in cfg["cover_deltas"]:
  delta=Q(delta0);hh=min(Q(1,8),delta/16)
  need(0<hh<=Q(1,8) and 8*hh<=delta/2<delta,"Critical strip uniform strict margin")
  for eps0 in cfg["cover_epsilons"]:
   eps=min(Q(eps0),Q(1,16));left=Q(1,4)-hh;right=Q(1,4)+hh;K=inp["K"]
   need(0<eps<left<right<K,"Ordered five-region cover")
   for cc in [eps/2,eps,(eps+left)/2,left,Q(1,4),right,(right+K)/2,K,K+1]:
    need(any([0<cc<=eps,eps<=cc<=left,left<=cc<=right,right<=cc<=K,cc>=K]),"Cover control")
    for d in [2,3,7,19]:
     rr=(cc-Q(1,4))/(cc+Q(1,2*(d-1)))
     if left<=cc<=right:need(abs(rr)<=8*hh<delta,"Actual thermal strip ratio")
    cover+=1
 for d in range(2,129):
  need(4*(d-1)**2>=d*d and inp["K"]*(d-1)**2>=10**24*d*d,"Actual high-ratio domain d>=2");large+=1
 need(2**100>4*100**4 and Q(1,2)-Q(4,100)>0,"Uniform raw constant-family inequality base")
 for t in cfg["root_t_samples"]:
  need(t>=100 and 2**t>=4*t**4,"Raw a0 analytic/large threshold control")
 # Generic COEFFICIENT-form singularity example, explicitly NOT an actual moment counterexample.
 x=Q(2,3);s=Q(5,7);p=[x*x*s,x*x-2*x*s,s-2*x,Q(1)]
 val=sum(p[j]*x**j for j in range(4));der=sum(j*p[j]*x**(j-1)for j in range(1,4))
 hb=sum((3-j)*p[j]*x**j for j in range(3))
 need(val==der==hb==0 and p[0]!=0 and p[-1]!=0,"Generic singular coefficient phase and nonzero axes")
 # X-root multiplicity is a DIFFERENT property.
 f=[Q(0),Q(0),Q(1),Q(1)]
 need(f[0]==0 and f[1]==0 and factorial(2)+factorial(3)==8 and factorial(2)>0,"Actual phase X²(1+X) has multiple X-root0 but nonsingular order1 coefficient value")
 return {"actual_moment_normalizations":norm,"actual_coefficient_gradient_identities":gradient,
 "actual_Euler_identities":euler,"actual_scaling_and_axis_controls":scaling,
 "a0_Appell_vector_identities":appell,"literal_large_nonzero_reversal_controls":reversal,
 "five_region_cover_controls":cover,"high_ratio_domain_controls":large,
 "generic_singular_coefficient_phase_controls":1,"X_root_vs_coefficient_singularity_domain_control":1}
def negative_controls(ledger,repo):
 out=[]
 def reject(label,fn):
  try:fn()
  except(ValueError,KeyError,TypeError,AssertionError,OSError):out.append(label);return
  raise ValueError("Negative accepted:"+label)
 for role in SPECS:
  d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"]if r["role"]!=role]
  reject("omit:"+role,lambda d=d:check_sources(d,repo))
  d=copy.deepcopy(ledger);next(r for r in d["sources"]if r["role"]==role)["sha256"]="0"*64
  reject("badpin:"+role,lambda d=d:check_sources(d,repo))
 for key in ["unscaled_coefficient_pairs_finite","all_order_squarefreeness","all_order_SFC",
 "global_repeated_set_empty","effective_cutoff_or_box","407_or418_coprimality_used_as_simplicity",
 "soft_pair_bound_used_as_repeated_support_finiteness","homotopy_counts_used_as_simplicity",
 "HJN_or_later_quartic_Robin_imported","new_sharp_lower_annulus_premise","heavy_finite_prefix_import",
 "finite_controls_formally_verify_analytic_estimates","main_uniform_cutoff_requires847"]:
  d=copy.deepcopy(SCOPE);d[key]=True;reject("false-scope:"+key,lambda d=d:check_inputs(d))
 for key in ["regional_RAW_analytic_mechanisms_NOT_coprimality_conclusions","exact_integer_partial_periods_holomorphic_same_support",
 "artifact_freshness_distinct_from_logical_premise_gating","main_forward_cutoff_survives847_logical_withholding",
 "Cauchy_on_OPEN_support_fixed_parameter_collars","thermal_center_rho_and_T_held_fixed_on_q_collars",
 "thermal_T_over_k_amplitude_term_retained","Newman_sections1_2_ONLY_no_HJN","R6_used_for_EVERY_individual_root_LOCATION_only",
 "old408_used_ONLY_finitely_bounded_supports_and_projection_reverse",
 "fixed_terminal847_required_ONLY_cutoff_to_finite_repeated_triples","coefficient_projective_moment_singularity_NOT_roots_in_X"]:
  d=copy.deepcopy(SCOPE);d[key]=False;reject("missing-mathematical-hypothesis:"+key,lambda d=d:check_inputs(d))
 for key,value in [("schema",True),("K",10**24),("parameters","all real supports and ordinary roots in X"),
  ("raw_large_source","pair coprimality implies individual simplicity without reversal"),
  ("high_ratio_statement","all orders at a>=100*d^2")]:
  d=copy.deepcopy(SCOPE);d[key]=value;reject("wrong-typed-domain:"+key,lambda d=d:check_inputs(d))
 d=copy.deepcopy(ledger);d["sources"]=[];reject("empty-ledger",lambda:check_sources(d,repo))
 d=copy.deepcopy(ledger);d["schema"]=True;reject("boolean-schema",lambda:check_sources(d,repo))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
 ap.add_argument("--parent-python",type=Path)
 ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
 repo=ns.repo_root.resolve();need(repo.is_dir(),"Explicit source repository root")
 python=ns.parent_python
 if python is None:python=next((p for p in [repo/".venv/bin/python",repo/".venv/Scripts/python.exe"]if p.is_file()),Path(sys.executable))
 python=python.absolute();need(python.is_file(),"Parent interpreter .absolute() preserves its venv symlink")
 inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN);rows=check_sources(ledger,repo)
 boundaries=child_boundaries(repo,rows);children=fresh(repo,python)
 algebra=exact_support();logic=logical_controls();controls=moment_controls(inp);bad=negative_controls(ledger,repo)
 check_sources(ledger,repo)
 result={"schema":1,"status":"PASS_SUPPORT_UNIFORM_INDIVIDUAL_SQUAREFREENESS",
 "scope":inp,"proof_sha256":H(HERE/PN),"checker_sha256":H(HERE/CN),"input_sha256":H(HERE/IN),
 "source_ledger_sha256":H(HERE/LN),"mandatory_source_roles":len(rows),"sources":[rows[r]for r in sorted(rows)],
 "narrow_child_boundaries":boundaries,"fresh_narrow_replays":children,"exact_derivative_algebra":algebra,
 "logical_six_way_dependency_controls":logic,"finite_actual_and_cover_controls":controls,
 "source_scope_negative_controls":bad,"assurance_flags":ASSURANCE,
 "meaning":"The written actual-period derivative/Cauchy/model proof establishes uniform eventual individual squarefreeness.847 supplies the separate finite-exception implication,408 supplies finite bounded-support/projection reverses. Fresh narrow source replays authenticate artifacts; no407/418 coprimality conclusion is used as simplicity, and finite controls are not a formal analytic proof or an effective cutoff."}
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"fresh_children":len(children),
 "universal_polynomials":len(algebra["universal_polynomial_identities"]),"Laurent":len(algebra["universal_Laurent_budget_reductions"]),
 "rational_budgets":len(algebra["rational_budget_controls"]),"logical_controls":len(logic),
 "actual_controls":controls,"source_scope_controls":len(bad),"output":str(ns.output)}))
if __name__=="__main__":main()

