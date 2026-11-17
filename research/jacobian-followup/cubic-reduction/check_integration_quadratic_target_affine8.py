#!/usr/bin/env python3
"""Audit complete affine-frame degree8 shear coverage and exact scope boundaries."""
from pathlib import Path
import json,hashlib,re,sys,subprocess,py_compile
B=Path(__file__).resolve().parent;O=B/'quadratic-target-affine8';Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in index['entries']};e=E['KCUB43'];vp=O/'validation.json';V=json.loads(vp.read_text())
 assert V['status']=='pass' and V['claim']=='KCUB43' and not V['goal_complete']
 assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS
 assert e['state']=='proved' and e['independent_replay'] and not e['formal_verification'] and not e['external_review']
 assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker']) and V['proof_sha256']==sha(Q/e['canonical_source'])
 assert set(e['dependencies'])==set(V['imports'])=={'KCUB41','KCUB42'} and e['narrows_problems']==['OP-GLDC'] and 'KCUB43' in E['OP-GLDC']['dependencies']
 for attack in E['OP-GLDC']['forbidden_attack_classes']:assert set(attack['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB41','KCUB42','KCUB43']:
  a=E[k];assert a['artifact_hash']=='sha256:'+sha(Q/a['checker'])
  for rel in [a['canonical_source'],a['checker']]+a['software_lock']:assert (Q/rel).is_file(),rel
 seen=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in seen:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  seen.add(k)
 visit('KCUB43',[])
 for k,pin in V['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 for key,rel in [('producer_sha256','quadratic_target_affine8.py'),('arithmetic_sha256','verify.py'),('input_sha256','quadratic-target-affine8/input.json'),('checkpoint_sha256','quadratic-target-affine8/checkpoint.json'),('rigidity_sha256','quadratic-target-affine8/rigidity.json')]:assert V[key]==sha(B/rel)
 C=json.loads((O/'checkpoint.json').read_text());assert C['status']=='complete' and C['source_sha256']==V['producer_sha256'] and C['input_sha256']==V['input_sha256'] and C['rigidity_sha256']==V['rigidity_sha256'] and C['matrix_entries']==483 and C['centre_parameters']==3 and C['shear_parameters']==2 and C['solver_queries']==C['candidate_grids']==0
 assert V['complete_shear_jet_matrix_entries']==483 and V['all_base_matrices_and_square_rows_reconstructed'] and V['exact_centre_dependent_output_equivalence_checked'] and V['five_formal_c2_separators_checked'] and V['all_alpha_zero_square_formulas_through6_checked'] and V['zero_branch_base_pivots_checked'] and V['full_generic_degree7_8_forms_checked']
 assert list(V['summed_rank_bounds'].values())==[13,14,14,15,15,16] and all(sum(V['exact_rank_profiles'][k])==v for k,v in V['summed_rank_bounds'].items())
 assert V['conditional_top_and_degree3_equality_coefficients_checked'] and V['universal_rank_and_UFD_equality_arguments_written_not_formalized'] and V['nonmonomial_perfect_square_has_at_least3_terms_written_char0_proof']
 assert [V[k] for k in ['minimum_nonlinear_shear_occurrences_all_affine_frames_lower','minimum_nonlinear_shear_occurrences_all_affine_frames_upper']]==[16,19] and V['sharp_nonlinear_family_minimum_not_established']
 assert V['all_affine_frames_of_complete_degree8_rational_target_family_minimum']==13 and V['equality13_requires_zero_shear_and_KCUB41_monomial_frames'] and V['no_raw13_improvement_in_entire_family'] and V['solver_queries']==V['candidate_grids']==0
 U=V['matching_existing_degree8_upper'];assert U['degree']==8 and U['counts']==[2,11,6] and U['determinant_one'] and not U['new_global_improvement']
 for k,rel in [('KCUB41','raw13-affine-jet/validation.json'),('KCUB42','raw13-target-field/validation.json')]:assert json.loads((B/rel).read_text())['proof_sha256']==sha(Q/E[k]['canonical_source'])
 parent=B/'raw13-target-field/integration-validation.json';PI=json.loads(parent.read_text());assert PI['status']=='pass' and PI['main_receipt_sha256']==sha(B/'raw13-target-field/validation.json') and PI['source_affine_changes_of_new_degree8_family_open_at_KCUB42_stage']
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB43' in l['claims']];research.validate_lessons(affected,index['entries'],check_reviews=True)
 a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC');assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims'])
 assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']) and all((Q/rel).is_file() for rel in a['sources'])
 assert 'frontier is empty' in a['next_step']  # positive source gate subsequently closed byKCUB44
 out=O/'integration-validation.json'
 for rel in ['QUADRATIC_TARGET_AFFINE_DEGREE8.md','quadratic-target-affine8/README.md','RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md','RAW13_AFFINE_ORBIT_MINIMUM.md','README.md']:
  f=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if '://' not in link:
    target=(f.parent/link.split('#')[0]).resolve();assert target==out or target.is_file(),link
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 for rel in ['quadratic_target_affine8.py','verify_quadratic_target_affine8.py','check_integration_quadratic_target_affine8.py']:
  f=B/rel;py_compile.compile(str(f),doraise=True);assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian-followup/README.md','research/jacobian-followup/cubic-reduction/README.md','research/jacobian-followup/cubic-reduction/EXPLORATION.md'];z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 result={'status':'pass','claim':'KCUB43','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'proof_sha256':V['proof_sha256'],'all_input_import_and_parent_receipt_pins_current':True,'complete_shear_entries':483,'all_six_parameter_centre_rank_branches':True,'both_generic_equality_obstructions_checked':True,'nonlinear_family_minimum_interval':[16,19],'nonlinear_family_sharp_minimum_established':False,'complete_degree8_target_affine_class_sharp_minimum':13,'equality_requires_zero_shear_and_old_monomial_frames':True,'no_solver_or_candidate_grid':True,'matching_control_is_not_global_improvement':True,'dependency_graph_acyclic':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'links_compilation_whitespace':True,'previous_goal_turn_classification':'progress:KCUB42','this_goal_turn_classification':'progress:KCUB43 complete joint-affine first-target-layer closure and nonlinear lower16','goal_complete':False,'global_positive_bounds_and_unrestricted_profiles_unchanged':True,'next_gate_at_KCUB43_stage':'Preflight a genuine positive-axis quadratic source shear and full target-field repair,with all coefficient-zero strata and universal certificate coverage.','navigation_attempts_current_parent_state':PI['navigation_attempts_current_parent_state'],'navigation_boundary':'Actual current-state attempts retained from the freshly rerun parent chain;two unrelated missing artifacts still block rendering. No generated status hand edit or artifact reconstruction.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS KCUB43 integration:483 full jets,all centre/shear branches and coupled equality premises;nonlinear[16,19],complete class sharp13;all pins,dependencies,reviews,work,links,compile,whitespace;global goal open')
if __name__=='__main__':main()
