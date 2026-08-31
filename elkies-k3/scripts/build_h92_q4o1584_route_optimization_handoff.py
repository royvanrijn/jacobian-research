#!/usr/bin/env python3
"""Build the durable q4/o1584 route-optimization handoff."""
import hashlib, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]; GEN=ROOT/'artifacts/generated-results'
FILES={
 'o1584_qq':ROOT/'artifacts/local/elkies-k3/q4o208-physical-q4o1584-rr-qq.json',
 'o164_equation_marking':ROOT/'artifacts/local/elkies-k3/q4o1584-second-affine-equation-marking-qq.json',
 'o164_qq':ROOT/'artifacts/local/elkies-k3/q4o1584-physical-q4o164-rr-qq.json',
 'q323_mod131_obstruction':ROOT/'artifacts/local/elkies-k3/q4o208-physical-q4o323-horizontal-mod131.json',
 'q323_lifted_shell_correction':ROOT/'artifacts/local/elkies-k3/q4o208-q4o323-horizontal-marking-qq.json',
 'q1599_qq':ROOT/'artifacts/local/elkies-k3/q4o208-q4o1599-a3-2a2-rr-qq.json',
 'q1599_cert':GEN/'elkies-k3-h3-q4o208-physical-q4o1599-a3-2a2-certificate.json',
 'q323_q207_candidate5887':ROOT/'artifacts/local/elkies-k3/q4o323-q207-deflated-horizontal-qq.json',
 'q323_q207_candidate5887_rr':ROOT/'artifacts/local/elkies-k3/q4o323-q207-smooth-rr-qq.json',
 'o1584_cert':GEN/'elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json',
 'o1584_zeros':GEN/'elkies-k3-h3-q4o208-physical-q4o1584-effective-zero-markings.json',
 'o164_cert':GEN/'elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json',
 'o164_zeros':GEN/'elkies-k3-h3-q4o208-q4o1584-q4o164-effective-zero-markings.json',
 'corrected_q8':GEN/'elkies-k3-h3-q4o208-corrected-a3-2a2-to-5a1-physical-q8-c10-certificate.json',
 'physical_q12':GEN/'elkies-k3-h3-physical-5a1-to-4a1-q12-certificate.json',
 'a2_cert':GEN/'elkies-k3-h3-physical-4a1-q4o21633-a2-certificate.json',
 'a2_marking':GEN/'elkies-k3-h3-physical-a2-effective-zero-marking.json',
 'a2_q4_shell':GEN/'elkies-k3-h3-physical-a2-q4d2-rootless-neighbors.json',
 'pinned_suffix':GEN/'elkies-k3-h3-q4o208-alternate-q4q8q12q4q8q8-pinned-r17-route-certificate.json',
 'a2_a1_cost':GEN/'elkies-k3-h3-physical-a2-q8d2-cap10000-growth-equation-cost.json',
 'a1_r17_cost':GEN/'elkies-k3-h3-physical-a2-q8o2102-a1-q8d2-cap10000-rootless-equation-cost.json',
 'q323_free_route':GEN/'elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q4o24401-q8o939-pinned-r17-route-certificate.json',
 'q12_direct_route':GEN/'elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12o4484-pinned-r17-route-certificate.json',
 'q12_rootless_cert':GEN/'elkies-k3-h3-q4o164-q8o376-q12o4484-rootless-certificate.json',
 'q12_rootless_neighbors':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap10000-all-neighbors.json',
 'q12_rootless_cost':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap10000-rootless-equation-cost.json',
 'q12_p0_word':GEN/'elkies-k3-h3-q4o164-q8o376-q12o4484-four-p0-section-word.json',
 'rootless_p0_frontier':GEN/'elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json',
 'q12_o5867_cert':GEN/'elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json',
 'q12_o5867_route':GEN/'elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12o5867-pinned-r17-route-certificate.json',
 'q12_expanded_pareto':GEN/'elkies-k3-h3-q4o164-q8o376-q12-expanded-compiler-pareto.json',
 'q12_low_pole_cert':GEN/'elkies-k3-h3-q4o164-q8o376-q12-fp-6-2-m13-10-rootless-certificate.json',
 'q12_low_pole_route':GEN/'elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12-fp-m13-10-pinned-r17-route-certificate.json',
 'q8_all_zero_known_horizontal':GEN/'elkies-k3-h3-q4o164-all-effective-zeros-q8-known-horizontal-audit.json',
 'q8_odd_i4_rational_sum':GEN/'elkies-k3-h3-q4o164-odd-i4-rational-sum-scan-p41.json',
 'p1_exact_source':ROOT/'artifacts/local/elkies-k3/h92-p1-oriented-lift-check.json',
 'q8_horizontal_qq':ROOT/'artifacts/local/elkies-k3/q4o164-q8o376-horizontal-crt-qq.json',
 'q8_smooth_rr_qq':ROOT/'artifacts/local/elkies-k3/q4o164-q8o376-smooth-rr-qq.json',
 'q12_o5867_degree1_mod131':ROOT/'artifacts/local/elkies-k3/q12o5867-degree1-compiler-branch-mod131.json',
 'q12_o5867_degree1_qq':ROOT/'artifacts/local/elkies-k3/q12o5867-degree1-compiler-branch-qq.json',
 'q12_o5867_shell_classification':ROOT/'artifacts/local/elkies-k3/q12o5867-p0-shell-lattice-classification-mod89.json',
 'q12_o5867_replacement_seeds':ROOT/'artifacts/local/elkies-k3/q12o5867-replacement-word-seeds-qq.json',
 'q12_o5867_horizontal_qq':ROOT/'artifacts/local/elkies-k3/q12o5867-target-horizontal-qq.json',
 'q12_o5867_smooth_rr_qq':ROOT/'artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json',
 'q12_o5867_endpoint_shell':ROOT/'artifacts/local/elkies-k3/q12o5867-rootless-p0-shell-mod131.json',
 'q12_o5867_endpoint_selection':ROOT/'artifacts/local/elkies-k3/q12o5867-rootless-mod131-selected-basis.json',
 'q12_o5867_endpoint_sections':ROOT/'artifacts/local/elkies-k3/q12o5867-rootless-selected-basis-qq.json',
 'q12_o5867_endpoint_height':ROOT/'artifacts/local/elkies-k3/q12o5867-rootless-height-basis-qq.json',
 'q8_direct_negative':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q8d2-cap10000-first-rootless-neighbor.json',
 'q10_direct_negative':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q10d2-cap10000-all-neighbors.json',
 'q20_direct_neighbors':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q14q16q18q20d2-cap10000-rootless-neighbors.json',
 'q20_direct_cost':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q20d2-cap10000-rootless-equation-cost.json',
 'q6_mw15_cost':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q6d2-cap10000-mw15-equation-cost.json',
 'q6_mw15_zeros':GEN/'elkies-k3-h3-q4o164-q8o376-q6o4013-2a1-effective-zero-markings.json',
 'q6_mw15_p1229_terminal':GEN/'elkies-k3-h3-q4o164-q8o376-q6o4013-2a1-p1229-q8d2-cap10000-rootless-equation-cost.json',
 'q6_mw15_second_terminal':GEN/'elkies-k3-h3-q4o164-q8o376-q6o4013-2a1-secondmissing-q8d2-cap10000-rootless-equation-cost.json',
 'q6_mw15_horizontal_terminal':GEN/'elkies-k3-h3-q4o164-q8o376-q6o4013-2a1-q6o4013_horizontal-q4q6q8d2-cap10000-rootless-neighbors.json',
 'q164_direct_high_norm_negative':GEN/'elkies-k3-h3-q4o164-c8-q10q12q14q16d2-cap10000-rootless-neighbors.json',
 'q8o10712_q12_negative':GEN/'elkies-k3-h3-q4o164-q8o10712-2a2-q12d2-cap10000-all-neighbors.json',
 'degree1_negative':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q4q6q8q10q12d1-cap10000-rootless-neighbors.json',
 'degree3_negative':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q9q12q15q18d3-cap10000-rootless-neighbors.json',
 'degree4_negative':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q8q12q16q20d4-cap10000-rootless-neighbors.json',
 'o164_o376_cost':GEN/'elkies-k3-h3-q4o164-c8-q8d2-cap10000-growth-equation-cost.json',
 'o376_o24401_cost':GEN/'elkies-k3-h3-q4o164-q8o376-4a1-p1229-q4d2-growth-equation-cost.json',
 'o24401_o939_cost':GEN/'elkies-k3-h3-q4o164-q8o376-q4o24401-a1-firstmissing-q8-cap10000-rootless-equation-cost.json',
 'old_q323_cost':GEN/'elkies-k3-h3-q4o208-corrected-q4o323-selected-equation-cost.json',
 'old_q8_cost':GEN/'elkies-k3-h3-corrected-a3-2a2-q8-5a1-selected-equation-cost.json',
 'old_q12_cost':GEN/'elkies-k3-h3-physical-5a1-q12-4a1-selected-equation-cost.json',
}
O1584_QQ_SCRIPT=ROOT/'elkies-k3/scripts/certify_h92_q4o208_physical_q4o1584_rr_qq.sage'
HASH_FILES=tuple(FILES.values())+(O1584_QQ_SCRIPT,)
data={k:json.loads(p.read_text()) for k,p in FILES.items()}
assert data['o1584_qq']['status']=='PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN'
assert data['o164_equation_marking']['status']=='PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING'
assert data['o164_qq']['status']=='PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN'
assert data['q323_mod131_obstruction']['status']=='PASS_EXACT_Q4O323_POLYNOMIAL_SECTION_SUBGROUP_OBSTRUCTION'
assert data['q323_lifted_shell_correction']['status']=='PASS_EXACT_QQ_Q4O323_LIFTED_SHELL_EXCLUDES_TARGET'
assert data['q1599_qq']['status']=='PASS_EXACT_QQ_Q4O1599_A3_2A2_RR_AND_JACOBIAN'
assert data['q1599_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['q323_q207_candidate5887']['status']=='PASS_EXACT_QQ_Q4O323_DEFLATED_CANDIDATE_SECTION'
assert data['q323_q207_candidate5887_rr']['status']=='REJECTED_EXACT_QQ_Q4O323_CANDIDATE5887_WRONG_6A1_CHILD'
assert data['o1584_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['o164_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['corrected_q8']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['physical_q12']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['a2_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['a2_marking']['status']=='PASS_EXACT_PHYSICAL_A2_EFFECTIVE_ZERO_MARKING'
assert not [n for n in data['a2_q4_shell']['neighbors'] if n['child_mw_rank']==17]
assert data['pinned_suffix']['status']=='PASS_EXACT_Q4O208_ALTERNATE_SUFFIX_TO_PINNED_R17'
assert data['q323_free_route']['status']=='PASS_EXACT_Q4O1584_Q4O164_Q8O376_Q4O24401_Q8O939_TO_PINNED_R17'
assert data['q12_direct_route']['status']=='PASS_EXACT_Q4O1584_Q4O164_Q8O376_Q12O4484_TO_PINNED_R17'
assert data['q12_rootless_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['q12_rootless_cost']['best_candidate']['candidate_id']=={'q':12,'old_fibre_degree':2,'orbit_index':4484}
assert data['q12_p0_word']['status']=='PASS_EXACT_Q12O4484_FOUR_P0_SECTION_WORD'
assert data['q12_p0_word']['selected_word']['minimum_q4o164_parent_degree_sum']==10
assert data['q12_p0_word']['selected_word']['minimum_q4o164_parent_degree_max']==3
assert data['rootless_p0_frontier']['status']=='PASS_EXACT_ROOTLESS_P0_SECTION_WORD_FRONTIER'
o5867_fibre=[6,2,-16,18,2,5,-11,2,13,11,-6,-41,6,5,-14,0,-3,-1,-2]
o5867=next(item for item in data['rootless_p0_frontier']['targets'] if item.get('fibre')==o5867_fibre)
o5867_fingerprint=hashlib.sha256(
 json.dumps(o5867_fibre,separators=(',',':')).encode()
).hexdigest()
assert o5867_fingerprint=='d676cab5918a08add2f743081cf76932d02adc9e18ab9c808dc05760fe0157dd'
assert o5867['fibre_fingerprint_sha256']==o5867_fingerprint
assert o5867['candidate_id']=={'q':12,'old_fibre_degree':2,'orbit_index':30357}
assert o5867['minimum_word_length_at_most_four']==4
assert o5867['best_words_by_exact_length']['4']['q4o164_parent_degree_sum']==8
assert o5867['best_words_by_exact_length']['4']['q4o164_parent_degree_max']==3
o5867_compiler=o5867['best_four_P0_word_by_parent_a_minus_b']
assert o5867_compiler['q4o164_parent_degree_sum']==8
assert o5867_compiler['q4o164_parent_degree_max']==3
assert o5867_compiler['q4o164_parent_a_minus_b_sum']==6
assert o5867_compiler['q4o164_parent_a_minus_b_max']==2
assert o5867['strict_parent_a_minus_b_improvements_with_more_P0_sections']['cutoff']==5
assert o5867['strict_parent_a_minus_b_improvements_with_more_P0_sections']['five_sections'] is None
assert o5867['strict_parent_a_minus_b_improvements_with_more_P0_sections']['six_sections'] is None
assert data['q12_o5867_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['q12_o5867_cert']['candidate_id']=={
 'q':12,'old_fibre_degree':2,'label':'q12o5867_after_q8o376'
}
assert data['q12_o5867_cert']['child']['root_data']==[0,0,1]
assert data['q12_o5867_route']['status']=='PASS_EXACT_Q323_FREE_Q4O1584_Q4O164_Q8O376_Q12O5867_PINNED_R17_ROUTE'
assert data['q12_o5867_route']['route'][-1]['candidate_id']==data['q12_o5867_cert']['candidate_id']
assert data['q12_expanded_pareto']['status']=='PASS_EXACT_Q12_EXPANDED_COMPILER_PARETO_AND_PINNED_R17_ALTERNATIVE'
assert data['q12_low_pole_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['q12_low_pole_route']['status']=='PASS_EXACT_Q323_FREE_Q4O1584_Q4O164_Q8O376_Q12_FP_M13_10_PINNED_R17_ROUTE'
assert data['q8_all_zero_known_horizontal']['status']=='PASS_BOUNDED_Q4O164_ALL_ZERO_Q8_NO_KNOWN_HORIZONTAL'
assert data['q8_all_zero_known_horizontal']['total_known_horizontal_candidates']==0
assert data['q8_odd_i4_rational_sum']['status']=='PASS_EXACT_QQ_Q4O164_ODD_I4_SECTION_LIFTS'
assert data['q8_odd_i4_rational_sum']['method']['rational_sum_scan_word_count_up_to_sign']==364
assert not data['q8_odd_i4_rational_sum']['exact_rational_signed_sums']
assert data['p1_exact_source']['status']=='PASS_EXACT_H92_P1'
assert data['p1_exact_source']['exact_weierstrass_square']
assert data['q8_horizontal_qq']['status']=='PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT'
assert data['q8_smooth_rr_qq']['status']=='PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO'
assert data['q8_smooth_rr_qq']['pointed_origin']['lattice_label']=='old_A11_component_6'
assert data['q8_smooth_rr_qq']['preferred_pointed_zero']['lattice_label']=='P1229'
assert data['q8_smooth_rr_qq']['preferred_pointed_zero']['source_exact_q8_degree']==1
assert data['q8_smooth_rr_qq']['preferred_pointed_zero']['quartic_sign']==-1
assert data['q12_o5867_degree1_mod131']['status']=='PASS_MOD131_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD'
assert data['q12_o5867_degree1_qq']['status']=='PASS_EXACT_QQ_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD'
assert data['q12_o5867_degree1_qq']['section']['degrees_x_y']==[4,6]
assert data['q12_o5867_degree1_qq']['section']['maximum_rational_bits']==800
assert data['q12_o5867_shell_classification']['status']=='PARTIAL_MODP_Q12O5867_AMBIGUOUS_SHELL_FOUR_WORD_FOUND'
replacement=data['q12_o5867_shell_classification']['q12o5867_realized_four_section_search']['Q2_intersection_refined_best_parent_degree_first_word']
assert replacement['lattice_class_indices']==[499,500,69,511]
assert replacement['q4o164_parent_degrees']==[4,2,1,5]
assert replacement['q4o164_parent_a_minus_b']==[4,2,1,4]
assert data['q12_o5867_replacement_seeds']['status']=='PASS_EXACT_QQ_Q12O5867_REPLACEMENT_CLASS499_CLASS511_SECTIONS'
assert data['q12_o5867_horizontal_qq']['status']=='PASS_EXACT_QQ_Q12O5867_TARGET_HORIZONTAL_SECTION'
assert data['q12_o5867_horizontal_qq']['target']['P_dot_O']==10
assert data['q12_o5867_smooth_rr_qq']['status']=='PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN'
assert data['q12_o5867_smooth_rr_qq']['smooth_RR']['ambient_dimension']==22
assert data['q12_o5867_smooth_rr_qq']['smooth_RR']['condition_rank']==20
assert data['q12_o5867_smooth_rr_qq']['smooth_RR']['h0']==2
assert data['q12_o5867_smooth_rr_qq']['child']['root_rank']==0
assert data['q12_o5867_endpoint_shell']['status']=='PASS_COMPLETE_MOD131_Q12O5867_ROOTLESS_P0_SHELL'
assert data['q12_o5867_endpoint_shell']['signed_section_count']==2622
assert data['q12_o5867_endpoint_selection']['status']=='PASS_MOD131_Q12O5867_ROOTLESS_REGULAR_SHORT_BASIS'
assert data['q12_o5867_endpoint_sections']['status']=='PASS_EXACT_QQ_Q12O5867_ROOTLESS_17_SELECTED_SECTIONS'
assert data['q12_o5867_endpoint_height']['status']=='PASS_EXACT_QQ_Q12O5867_ROOTLESS_RANK17_HEIGHT_BASIS_PINNED'
assert data['q12_o5867_endpoint_height']['height_gram_determinant']==948
assert sum(item['child_mw_rank']==17 for item in data['q12_rootless_neighbors']['neighbors'])==2
assert not [item for item in data['q8_direct_negative']['neighbors'] if item['child_mw_rank']==17]
assert not [item for item in data['q10_direct_negative']['neighbors'] if item['child_mw_rank']==17]
assert sum(item['child_mw_rank']==17 for item in data['q20_direct_neighbors']['neighbors'])==4
assert data['q20_direct_cost']['best_candidate']['candidate_id']=={'q':20,'old_fibre_degree':2,'orbit_index':2725}
assert data['q20_direct_cost']['best_candidate']['equation_cost_score']==65805
assert data['q6_mw15_cost']['best_candidate']['candidate_id']=={'q':6,'old_fibre_degree':2,'orbit_index':4013}
assert data['q6_mw15_zeros']['status']=='PASS_EXACT_PHYSICAL_AN_ALL_EFFECTIVE_ZERO_MARKINGS'
assert data['q6_mw15_p1229_terminal']['best_candidate']['equation_cost_score']==75592
assert data['q6_mw15_second_terminal']['best_candidate']['equation_cost_score']==75592
assert not [item for item in data['q6_mw15_horizontal_terminal']['neighbors'] if item['child_mw_rank']==17]
for key in ('q164_direct_high_norm_negative','q8o10712_q12_negative','degree1_negative','degree3_negative','degree4_negative'):
 assert not [item for item in data[key]['neighbors'] if item['child_mw_rank']==17]
payload={
 'schema':'elkies-k3.h3-q4o1584-route-optimization-handoff.v2',
 'status':'PASS_EXACT_Q4O1584_Q4O164_Q8O376_Q12_ROOTLESS_OPTIONS_HANDOFF_PROMOTED',
 'equation_frontier':'q12/o5867 and seventeen exact independent endpoint sections are complete over QQ; the active gates are Picard/discriminant, full-MW saturation promotion, and source identity',
 'best_equation_realizable_branch':[
  {'edge':'q4/orbit1584','source_zero':'old_A11_component_5','literal_divisor':'O(C5)+first_I6_affine+r3+2r4+r5','P_dot_O':0,'expected_RR_ambient':4,'child':'3A1+A3+D4/MW7','lattice_certificate':str(FILES['o1584_cert'].relative_to(ROOT)),'equation_status':'PASS exact QQ: literal 4-to-2 jet RR construction and global Jacobian','equation_artifact':str(FILES['o1584_qq'].relative_to(ROOT)),'equation_script':str(O1584_QQ_SCRIPT.relative_to(ROOT))},
  {'edge':'q4/orbit164','source_zero':'second_I6_affine_component','horizontal':'old_A11_component_0','literal_divisor':'O(second_I6_affine)+old_A11_component_0+2F+(F-r7)+(F-r9)','P_dot_O':0,'expected_RR_ambient':4,'child':'2A1+2A3/MW9','lattice_certificate':str(FILES['o164_cert'].relative_to(ROOT)),'equation_status':'PASS exact QQ: physical zero and C0 resolved, RR 4-to-2-to-2, global Jacobian','equation_marking_artifact':str(FILES['o164_equation_marking'].relative_to(ROOT)),'equation_artifact':str(FILES['o164_qq'].relative_to(ROOT)),'measured_runtime_seconds':185.65781034900283,'maximum_quartic_rational_bits':886429,'maximum_jacobian_AB_rational_bits':1329443},
 ],
 'physical_suffix_improvement':[
  {'edge':'corrected A3+2A2 to 5A1','source_zero':'old_A11_component_10','q':8,'old_fibre_degree':2,'P_dot_O':2,'child':'5A1/MW12','certificate':str(FILES['corrected_q8'].relative_to(ROOT))},
  {'edge':'5A1 to 4A1','source_zero':'old_A11_component_5','q':12,'old_fibre_degree':2,'P_dot_O':4,'child':'4A1/MW13','certificate':str(FILES['physical_q12'].relative_to(ROOT))},
  {'edge':'4A1 to A2','source_zero':'old_A11_component_0','q':4,'old_fibre_degree':2,'P_dot_O':1,'expected_RR_ambient':5,'child':'A2/MW15','certificate':str(FILES['a2_cert'].relative_to(ROOT)),'physical_child_marking':str(FILES['a2_marking'].relative_to(ROOT))},
 ],
 'promoted_physical_pinned_suffix':{
  'route':'current physical 3A3 --q4 corrected o323--> A3+2A2 --q8--> 5A1 --q12--> 4A1 --q4/o21633--> A2 --q8/o2102--> A1 --q8/o5165--> pinned rootless R17',
  'certificate':str(FILES['pinned_suffix'].relative_to(ROOT)),
  'endpoint_identification':'exact integral determinant-one isometry to pinned rank17_gram.txt',
  'late_edge_equation_costs':[
   {'edge':'A2 to A1 q8/o2102','P_dot_O':5,'expected_RR_ambient':13,'score':81176,'cost_artifact':str(FILES['a2_a1_cost'].relative_to(ROOT))},
   {'edge':'A1 to R17 q8/o5165','P_dot_O':6,'expected_RR_ambient':14,'score':90584,'cost_artifact':str(FILES['a1_r17_cost'].relative_to(ROOT))},
  ],
  'promotion_scope':'Promoted as the complete physical lattice lifting suffix. The last two edges require new horizontal sections and are not equation-lifted.',
 },
 'promoted_q323_free_pinned_route':{
  'route':'current physical 3A3 --q4/o1584--> D4+A3+3A1/MW7 --q4/o164--> 2A3+2A1/MW9 --q8/o376--> 4A1/MW13 --q12/o5867--> pinned rootless R17',
  'certificate':str(FILES['q12_o5867_route'].relative_to(ROOT)),
  'endpoint_identification':'exact determinant-minus-one integral isometry to pinned rank17_gram.txt, with full marked U and bidirectional unimodular NS transports',
  'compiler_frames':[
   {'edge':'q4/o1584','source_zero':'old_A11_component_5','horizontal':'first_I6_affine_component','literal_divisor':'O(C5)+first_I6_affine+r3+2r4+r5','P_dot_O':0,'expected_RR_ambient':4,'equation_status':'PASS exact QQ'},
   {'edge':'q4/o164','source_zero':'second_I6_affine_component','horizontal':'old_A11_component_0','literal_divisor':'O(second_I6_affine)+old_A11_component_0+2F+(F-r7)+(F-r9)','P_dot_O':0,'expected_RR_ambient':4,'equation_status':'PASS exact QQ','equation_marking_artifact':str(FILES['o164_equation_marking'].relative_to(ROOT)),'equation_artifact':str(FILES['o164_qq'].relative_to(ROOT))},
   {'edge':'q8/o376','source_zero':'old_A11_component_8','selected_child_zero':'P1229','horizontal_section':[5,1,-7,-2,-1,-8,-2,-9,-18,-12,-1,-10,1,6,-1,-1,3,1,1],'horizontal_in_already_explicit_subgroup':False,'unspanned_horizontal_rank_gap':5,'P_dot_O':4,'expected_RR_ambient':12,'explicit_degree_zero':['second_old_I6_I4_missing_component'],'explicit_degree_one':['P1229','first_old_I6_I4_missing_component','old_A11_component_0','old_A11_component_5','old_A11_component_6','old_A11_component_9','old_zero','second_I6_affine_component'],'score':28735,'cost_artifact':str(FILES['o164_o376_cost'].relative_to(ROOT)),'preferred_missing_direction_construction':{'strategy':'transport exact inherited P1 through the three certified q4 models, retain its unsplit degree-seven divisor, and Abel-reduce before QQ reconstruction','exact_source_artifact':str(FILES['p1_exact_source'].relative_to(ROOT)),'changed_zero_negative_audit':str(FILES['q8_all_zero_known_horizontal'].relative_to(ROOT)),'odd_I4_signed_sum_negative_audit':str(FILES['q8_odd_i4_rational_sum'].relative_to(ROOT)),'equation_artifact':str(FILES['q8_horizontal_qq'].relative_to(ROOT)),'resolved_rr_artifact':str(FILES['q8_smooth_rr_qq'].relative_to(ROOT)),'resolved_rr_dimensions':[12,4,2],'oriented_T0_origin':'old_A11_component_6','opposite_T0_sign':'old_A11_component_5','P1229_support':'25281/168246841','P1229_quartic_sign':-1,'pointing_identities':['81*A_P1229_pointed=A_child','729*B_P1229_pointed=B_child'],'status':'PASS exact QQ horizontal, resolved h0=2 pencil, quartic, 4A1 Jacobian, oriented T=0 marking, and direct P1229 pointing'}},
   {'edge':'q12/o5867','source_zero':'P1229','stable_candidate_identity':{'promoted_label':'q12o5867_after_q8o376','fibre':o5867_fibre,'fibre_fingerprint_sha256':o5867_fingerprint,'bounded_frontier_local_orbit_index':30357,'identity_note':'The promoted label and exact fibre fingerprint are stable; orbit 30357 is local to the bounded frontier enumeration.'},'horizontal_section':o5867['horizontal_section'],'P_dot_O':o5867['direct_P_dot_O'],'expected_RR_ambient':o5867['direct_expected_RR_ambient'],'score':56681,'lattice_certificate':str(FILES['q12_o5867_cert'].relative_to(ROOT)),'route_certificate':str(FILES['q12_o5867_route'].relative_to(ROOT)),'nominal_lattice_compiler':{'q4o164_parent_degrees':[3,2,1,2],'q4o164_parent_a_minus_b':[2,2,1,1],'status':'planning proxy; Q1 absent from complete shells at four good primes'},'equation_effective_compiler':{'equation_word':'499+500+69+511-489+933-913','new_lattice_classes':replacement['lattice_class_indices'],'q4o164_parent_degrees':replacement['q4o164_parent_degrees'],'q4o164_parent_a_minus_b':replacement['q4o164_parent_a_minus_b'],'classification_artifact':str(FILES['q12_o5867_shell_classification'].relative_to(ROOT)),'replacement_seed_artifact':str(FILES['q12_o5867_replacement_seeds'].relative_to(ROOT)),'horizontal_artifact':str(FILES['q12_o5867_horizontal_qq'].relative_to(ROOT)),'resolved_rr_artifact':str(FILES['q12_o5867_smooth_rr_qq'].relative_to(ROOT)),'resolved_rr_dimensions':[22,2],'child_fibres':'24I1','status':'PASS exact QQ horizontal, h0=2 smooth RR plane, quartic, and rootless Jacobian'},'endpoint_basis':{'mod131_shell_artifact':str(FILES['q12_o5867_endpoint_shell'].relative_to(ROOT)),'mod131_selected_basis_artifact':str(FILES['q12_o5867_endpoint_selection'].relative_to(ROOT)),'exact_QQ_sections_artifact':str(FILES['q12_o5867_endpoint_sections'].relative_to(ROOT)),'exact_height_artifact':str(FILES['q12_o5867_endpoint_height'].relative_to(ROOT)),'section_count':17,'height_determinant':948,'pinned_transport_determinant':-1,'status':'PASS exact QQ rank-at-least-17 height basis; full-MW saturation awaits Picard/discriminant gate'}},
  ],
  'promotion_scope':'q12/o5867 is the preferred optional final edge after q8/o376. Its edge and complete route are exactly lattice-certified to pinned R17; q12/o4484 is retained only as a certified fallback.',
  'equation_cost_comparison':{
   'new_edge_scores':[-1580,-1620,28735,56681],
   'new_total':82216,
   'previous_q323_free_total':149060,
   'saving_over_previous_q323_free':66844,
   'percent_saving_over_previous_q323_free':44.84368710586341,
   'superseded_physical_edge_scores':[29439,39489,56154,61493,81176,90584],
   'superseded_total':358335,
   'absolute_saving':276119,
   'percent_saving':77.05610671578272,
   'scope':'Same deterministic equation-cost scorer on each fully physical edge from current 3A3 to pinned R17; planning score, not measured QQ runtime.',
   'superseded_cost_artifacts':[str(FILES['old_q323_cost'].relative_to(ROOT)),str(FILES['old_q8_cost'].relative_to(ROOT)),str(FILES['old_q12_cost'].relative_to(ROOT))],
  },
  'superseded_five_edge_q323_free_fallback':{'route':'q4/o1584, q4/o164, q8/o376, q4/o24401, q8/o939','certificate':str(FILES['q323_free_route'].relative_to(ROOT)),'score':149060},
  'direct_rootless_norm_screen':{'source':'physical 4A1 with P1229 zero','sample_scope':'The original cross-norm screen used 10,000 MW quotient vectors. The q12 shell was subsequently expanded to 50,000 of 114,347,416 vectors; bounded search, not a completeness theorem.','q8_rootless_count':0,'q8_artifact':str(FILES['q8_direct_negative'].relative_to(ROOT)),'q10_rootless_count':0,'q10_artifact':str(FILES['q10_direct_negative'].relative_to(ROOT)),'q12_rootless_count_at_cap10000':2,'q12_rootless_count_at_cap50000':28,'q12_expanded_pareto_artifact':str(FILES['q12_expanded_pareto'].relative_to(ROOT)),'q14_q16_q18_rootless_count':0,'q20_rootless_count':4,'q20_best_score':65805,'q20_cost_disadvantage':9125,'q20_artifact':str(FILES['q20_direct_cost'].relative_to(ROOT))},
 },
 'q12o4484_fallback':{
  'route':'current physical 3A3 --q4/o1584--> D4+A3+3A1/MW7 --q4/o164--> 2A3+2A1/MW9 --q8/o376--> 4A1/MW13 --q12/o4484--> pinned rootless R17',
  'route_certificate':str(FILES['q12_direct_route'].relative_to(ROOT)),
  'edge_certificate':str(FILES['q12_rootless_cert'].relative_to(ROOT)),
  'cost_artifact':str(FILES['q12_rootless_cost'].relative_to(ROOT)),
  'section_word_artifact':str(FILES['q12_p0_word'].relative_to(ROOT)),
  'compiler_profile':{'individual_P_dot_O':[0,0,0,0],'q4o164_parent_degrees':[3,2,3,2],'parent_degree_sum':10,'parent_a_minus_b_sum':8},
  'retention_scope':'Fully lattice-certified fallback; superseded as the preferred optional compiler by q12/o5867.',
 },
 'q12_low_pole_alternative':{
  'route':'current physical 3A3 --q4/o1584--> D4+A3+3A1/MW7 --q4/o164--> 2A3+2A1/MW9 --q8/o376--> 4A1/MW13 --q12/fingerprint246597bd--> pinned rootless R17',
  'route_certificate':str(FILES['q12_low_pole_route'].relative_to(ROOT)),
  'edge_certificate':str(FILES['q12_low_pole_cert'].relative_to(ROOT)),
  'comparison_artifact':str(FILES['q12_expanded_pareto'].relative_to(ROOT)),
  'stable_fibre_identity':data['q12_expanded_pareto']['low_pole_alternative']['fibre'],
  'fibre_fingerprint_sha256':data['q12_expanded_pareto']['low_pole_alternative']['fibre_fingerprint_sha256'],
  'endpoint_identification':'exact determinant-one integral isometry to pinned rank17_gram.txt, with full marked U and bidirectional unimodular NS transports',
  'section_compiler':data['q12_expanded_pareto']['low_pole_alternative']['compiler'],
  'compiler_profile':data['q12_expanded_pareto']['low_pole_alternative']['compiler_word_profile'],
  'comparison_with_q12o5867':'Uses three new sections with current P.O=(0,0,1), parent a-b=(0,1,2), and parent degrees=(3,2,3). Relative to o5867 four-P0, the parent a-b sum drops 6 to 3 and branch count 4 to 3 at unchanged parent-degree sum 8; the named correction count rises 2 to 4.',
  'retention_scope':'Fully lattice-certified alternative. It uses one fewer branch and a smaller parent pole proxy, but introduces a P.O=1 branch and four corrections; it does not supersede the all-P.O=0 q12/o5867 compiler.',
 },
 'promoted_optional_q12o5867_endpoint':{
  'route':'current physical 3A3 --q4/o1584--> D4+A3+3A1/MW7 --q4/o164--> 2A3+2A1/MW9 --q8/o376--> 4A1/MW13 --q12/o5867--> pinned rootless R17',
  'route_certificate':str(FILES['q12_o5867_route'].relative_to(ROOT)),
  'edge_certificate':str(FILES['q12_o5867_cert'].relative_to(ROOT)),
  'endpoint_identification':'exact determinant-minus-one integral isometry to pinned rank17_gram.txt, with full marked U and bidirectional unimodular NS transports',
  'stable_candidate_identity':{
   'promoted_label':'q12o5867_after_q8o376',
   'fibre':o5867_fibre,
   'fibre_fingerprint_sha256':o5867_fingerprint,
   'bounded_frontier_local_orbit_index':30357,
   'identity_note':'The promoted label and exact fibre fingerprint are stable; orbit 30357 is local to the bounded frontier enumeration.',
  },
  'direct_horizontal':o5867['horizontal_section'],
  'direct_P_dot_O':o5867['direct_P_dot_O'],
  'direct_expected_RR_ambient':o5867['direct_expected_RR_ambient'],
  'preferred_section_compiler':{
   'strategy':'nominal four-P.O=0 lattice cost word; retained as a planning proxy',
   'artifact':str(FILES['rootless_p0_frontier'].relative_to(ROOT)),
   'new_section_count':4,
   'individual_P_dot_O':[0,0,0,0],
   'q4o164_parent_degrees':[item['q4o164_parent_degree'] for item in o5867_compiler['new_sections']],
   'q4o164_parent_a_minus_b':[item['q4o164_parent_a_minus_b'] for item in o5867_compiler['new_sections']],
   'parent_degree_sum':o5867_compiler['q4o164_parent_degree_sum'],
   'parent_degree_max':o5867_compiler['q4o164_parent_degree_max'],
   'parent_a_minus_b_sum':o5867_compiler['q4o164_parent_a_minus_b_sum'],
   'parent_a_minus_b_max':o5867_compiler['q4o164_parent_a_minus_b_max'],
   'known_section_correction':o5867_compiler['known_section_correction'],
   'lifted_branch_count':1,
   'degree1_branch':{
    'index_in_word':data['q12_o5867_degree1_qq']['compiler_branch']['index_in_promoted_four_section_word'],
    'exact_group_word':data['q12_o5867_degree1_qq']['exact_lattice_to_equation_word']['identity'],
    'mod131_artifact':str(FILES['q12_o5867_degree1_mod131'].relative_to(ROOT)),
    'QQ_artifact':str(FILES['q12_o5867_degree1_qq'].relative_to(ROOT)),
    'child_degrees':data['q12_o5867_degree1_qq']['section']['degrees_x_y'],
    'maximum_rational_bits':data['q12_o5867_degree1_qq']['section']['maximum_rational_bits'],
    'hensel_dimensions':[13,12,12],
   },
   'status':'nominal word not realized: its Q1 branch is absent from complete shells at four good primes',
  },
  'equation_effective_compiler':{
   'equation_word':'499+500+69+511-489+933-913',
   'new_lattice_classes':replacement['lattice_class_indices'],
   'q4o164_parent_degrees':replacement['q4o164_parent_degrees'],
   'q4o164_parent_a_minus_b':replacement['q4o164_parent_a_minus_b'],
   'classification_artifact':str(FILES['q12_o5867_shell_classification'].relative_to(ROOT)),
   'replacement_seed_artifact':str(FILES['q12_o5867_replacement_seeds'].relative_to(ROOT)),
   'horizontal_artifact':str(FILES['q12_o5867_horizontal_qq'].relative_to(ROOT)),
   'resolved_rr_artifact':str(FILES['q12_o5867_smooth_rr_qq'].relative_to(ROOT)),
   'resolved_rr_dimensions':[22,2],
   'jacobian_degrees':[8,12,24],
   'child_fibres':'24I1',
   'status':'PASS exact QQ rootless equation and seventeen-section determinant-948 height basis; full-MW saturation awaits Picard/discriminant gate',
  },
  'comparison_with_q12o4484':'Same q=12, P.O=10, RR22 and four P.O=0 branches. Under the parent-pole-first optimization, parent-degree sum drops 10 to 8 and parent a-b sum drops 8 to 6; both optimized words use two named corrections. Exhaustive length-five and length-six P.O=0 searches find no word with parent a-b total at most 5.',
  'mixed_three_section_fallback':o5867['mixed_low_pole_words']['three_P0_P0_P1']['best_word'],
  'promotion_scope':'PROMOTED as the preferred optional final edge and now equation-explicit over QQ. The nominal (3,2,1,2)/(2,2,1,1) word remains a planning proxy; the executed compiler uses (4,2,1,5)/(4,2,1,4). The endpoint marking to pinned R17 is exact.',
 },
 'negative_search_results':[
  {'hub':'physical A2/MW15','shell':'q4, old-fibre degree 2','result':'complete 97030-orbit shell contains no rootless/MW17 child','artifact':str(FILES['a2_q4_shell'].relative_to(ROOT))},
  {'hub':'physical q4/o164 2A3+2A1/MW9 with C8 zero','shell':'q10,q12,q14,q16, degree 2, cap 10000','result':'79793 screened orbits; no rootless child, maximum MW13','artifact':str(FILES['q164_direct_high_norm_negative'].relative_to(ROOT))},
  {'hub':'q8/o10712 2A2/MW13 provisional exact frame','shell':'q12, degree 2, cap 10000','result':'11874 screened orbits; no rootless child, maximum MW15','artifact':str(FILES['q8o10712_q12_negative'].relative_to(ROOT)),'boundary':'Negative already in the exact provisional frame, so no physical effective-zero promotion was needed.'},
  {'hub':'physical 4A1/MW13 with P1229 zero','shell':'q6 to physical 2A1/MW15, then q4,q6,q8 from all three effective zeros','result':'best entry score 48395; two terminal q8/rootless scores are 75592 and the new-horizontal zero has no rootless child; detour is dominated','entry_cost_artifact':str(FILES['q6_mw15_cost'].relative_to(ROOT)),'zero_artifact':str(FILES['q6_mw15_zeros'].relative_to(ROOT))},
  {'hub':'physical 4A1/MW13 with P1229 zero','shell':'degrees 1,3,4; q through 12,18,20 respectively; cap 10000','result':'145407 screened orbits total; no rootless child. Degree 1 preserves all four roots in every tested shell.','artifacts':[str(FILES[k].relative_to(ROOT)) for k in ('degree1_negative','degree3_negative','degree4_negative')]},
  {'hub':'fixed-corridor q4/o323 experiment toward q207','candidate':5887,'result':'Exact QQ P.O=10 section and RR 22-to-20-to-2 compile to the wrong child: 6I2+12I1 (6A1), while the good-fibre Abel trace is (18,17) instead of q207 target (50,50). Candidate5903 also fails its marked Abel trace; a widened 120-triple P.O=0-shell MITM produced 13 degree collisions and zero trace survivors.','measured_section_degrees_X_Y_Z':[24,34,10],'maximum_section_rational_bits':1604,'jacobian_degrees_A_B_Delta':[8,12,24],'maximum_jacobian_rational_bits':33977,'horizontal_artifact':str(FILES['q323_q207_candidate5887'].relative_to(ROOT)),'resolved_rr_artifact':str(FILES['q323_q207_candidate5887_rr'].relative_to(ROOT)),'status':'certified negative; no actual q207 compiler-selected zero'},
 ],
 'fixed_corridor_equation_correction':{
  'withdrawn_edge':'q4/o323 from current physical q4/o208 3A3 chamber',
  'withdrawal_artifact':str(FILES['q323_lifted_shell_correction'].relative_to(ROOT)),
  'reason':'The branch33 inverse is only a local target candidate and occurs in neither complete signed eleven-section graph solution; exact summand replay confirms exclusion.',
  'replacement_edge':'q4/o1599 to A3+2A2/MW10',
  'replacement_equation_artifact':str(FILES['q1599_qq'].relative_to(ROOT)),
  'replacement_lattice_certificate':str(FILES['q1599_cert'].relative_to(ROOT)),
  'horizontal_branch':16,
  'P_dot_O':1,
  'resolved_RR_dimensions':[5,3,2],
  'child_root_data':[7,24,36],
  'child_fibres':'I4+2I3+14I1; Euler 24; MW10 conditional on rho=19',
  'compiler_method':'exact lifted branch, resolved RR, quartic and Jacobian; no large Groebner calculation',
  'splice_gate':'Two inherited-anchor NS markings remain; one additional independent exact section anchor is required before composition with a pinned-R17 suffix.',
  'route_effect':'Does not alter the promoted q323-free q4/o1584--q4/o164--q8/o376--q12/o5867 route.',
 },
 'q323_compiler_obstruction':{
  'artifact':str(FILES['q323_mod131_obstruction'].relative_to(ROOT)),
  'prime':131,
  'polynomial_section_count':62,
  'polynomial_section_MW_rank':7,
  'ambient_MW_rank':8,
  'target_augments_rank':True,
  'conclusion':'The former q4/o323 horizontal is outside the complete tested polynomial-section subgroup by one MW direction; this independently reinforces the q323-free promotion.',
  'gamma_label_gate':'RESOLVED EXACTLY over QQ: C7 is the oriented component-one gamma branch, the opposite sign is second_I6_affine_component, and their sum is first_I6_affine_component. The q4/o164 equation lift uses the physical second-affine branch.',
 },
 'lifting_agent_notice':{
  'direct_current_3A3_track':'The older physical suffix begins with corrected q4/o323 and its hidden P.O=1 horizontal. It remains certified but is no longer preferred.',
  'q323_avoiding_track':'q4/o1584, q4/o164 and q8/o376 have the preferred optional q12/o5867 rootless continuation plus certified fingerprint246597bd and q12/o4484 alternatives; all endpoints are exactly pinned to R17.',
  'recommended_action':'Use the exact rank-17 section basis to prove a Picard/discriminant upper bound, promote full-MW saturation, and close the source-identity certificate.',
 },
 'promotion_gate':'The q12/o5867 four-edge route is fully lattice-certified to pinned R17 and promoted as the preferred optional final edge. Its corrected equation-effective word, horizontal, smooth h0=2 RR plane, quartic, rootless Jacobian, and seventeen exact independent sections with pinned determinant-948 height Gram are exact over QQ. Picard/discriminant, full-MW saturation promotion, and source identity remain. Fingerprint246597bd and q12/o4484 remain certified alternatives.',
 'inputs':{'paths':[str(p.relative_to(ROOT)) for p in HASH_FILES],'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in HASH_FILES}},
}
out=GEN/'elkies-k3-h3-q4o1584-route-optimization-handoff.json';out.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(out)
