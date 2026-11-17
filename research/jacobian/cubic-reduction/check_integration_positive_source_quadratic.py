#!/usr/bin/env python3
"""Audit the universal positive-source family proof and exact base-field9 certificate."""
from pathlib import Path
import json,hashlib,re,sys,subprocess,py_compile
B=Path(__file__).resolve().parent;O=B/'positive-source-quadratic';Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in index['entries']};e=E['KCUB44'];vp=O/'validation.json';V=json.loads(vp.read_text())
 assert V['status']=='pass' and V['claim']=='KCUB44' and not V['goal_complete']
 assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS
 assert e['state']=='proved' and e['independent_replay'] and not e['formal_verification'] and not e['external_review']
 assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker']) and V['proof_sha256']==sha(Q/e['canonical_source'])
 assert set(e['dependencies'])==set(V['imports'])=={'KCUB42','KCUB43'} and e['narrows_problems']==['OP-GLDC'] and 'KCUB44' in E['OP-GLDC']['dependencies']
 for attack in E['OP-GLDC']['forbidden_attack_classes']:assert set(attack['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB42','KCUB43','KCUB44']:
  a=E[k];assert a['artifact_hash']=='sha256:'+sha(Q/a['checker'])
  for rel in [a['canonical_source'],a['checker']]+a['software_lock']:assert (Q/rel).is_file(),rel
 seen=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in seen:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  seen.add(k)
 visit('KCUB44',[])
 for k,pin in V['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 for key,rel in [('producer_sha256','positive_source_quadratic.py'),('degree9_producer_sha256','raw13_target_field9.py'),('arithmetic_sha256','verify.py'),('companion_helper_sha256','verify_raw13_target_field.py'),('input_sha256','positive-source-quadratic/input.json'),('checkpoint_sha256','positive-source-quadratic/checkpoint.json'),('degree9_input_sha256','positive-source-quadratic/base-field9/input.json'),('degree9_checkpoint_sha256','positive-source-quadratic/base-field9/checkpoint.json'),('discovery_sha256','positive-source-quadratic/discovery.json')]:assert V[key]==sha(B/rel)
 C=json.loads((O/'checkpoint.json').read_text());D=json.loads((O/'input.json').read_text());assert C['status']=='complete' and C['source_sha256']==D['source_sha256']==V['producer_sha256'] and C['input_sha256']==V['input_sha256'] and C['parameters']==3 and C['obstruction_heads']==5 and C['solver_queries']==0 and C['formal_order_bound']==15
 assert D['source_weights']==[1,-1,-2] and D['parameter_weights']==[3,4,5] and D['obstruction_a_exponents']==[5,2,3,1,0] and D['obstruction_nonzero_constants']==[6,-108,144,24,12] and [s['order'] for s in D['leading_parameter_strata']]==[3,4,5]
 C9=json.loads((O/'base-field9/checkpoint.json').read_text());D9=json.loads((O/'base-field9/input.json').read_text());assert C9['status']=='complete' and C9['source_sha256']==D9['source_sha256']==V['degree9_producer_sha256'] and C9['input_sha256']==V['degree9_input_sha256'] and [C9[k] for k in ['columns','rank','points']]==[220,213,107]
 assert V['degree9_base_source_columns']==220 and V['degree9_base_field_kernel_dimension']==7 and V['degree9_independent_minor_size']==213 and V['exact_rational_specialization_points']==107
 assert V['prime']==1000000007 and V['prime_checked'] and V['minor_determinant_mod_prime']==116435059 and V['all_specialization_E_determinants_nonzero']
 assert V['all_source_quadratic_parameters_retained'] and V['exact_source_inverse_and_conjugacy_weights_checked'] and V['five_Taylor_obstruction_heads_checked'] and V['formal_leading_order_cases']==[{'leading_order':3,'last_order_used':9},{'leading_order':4,'last_order_used':12},{'leading_order':5,'last_order_used':15}]
 assert V['universal_formal_induction_and_weighted_initial_space_proof_written_not_formalized']
 assert V['source_degree7_intersection_also_dimension2'] is True
 assert V['every_nonzero_quadratic_low_source8_field_intersection_dimension']==2 and V['intersection_basis']==['1','H1 composed E_q'] and V['target_degree_cap'] is None and V['target_denominator_cap'] is None and V['no_Keller_rational_target_repair_through_source_degree8']
 assert V['extends_to_lower_affine_terms_and_affine_frames_after_the_source_shear'] and V['does_not_classify_nonlinear_source_shears_in_arbitrary_intermediate_frames'] and V['zero_quadratic_base_kernel_dimensions_degree7_8_9']==[4,7,7]
 assert V['transformed_source_degree9_not_classified'] and V['endpoint_degree_interval_for_nonzero_source_quadratic']==[9,10]
 U=V['matching_degree10_control'];assert U['row_degrees']==[7,9,10] and U['nonlinear_counts']==[8,14,15] and U['full_determinant_one'] and not U['global_improvement'] and V['all_parameter_collision_transport_checked']
 discovery=json.loads((O/'discovery.json').read_text());assert len(discovery['models'])==7 and V['pilot_parameter_triples_are_discovery_only'] and 'No universal field or parameter exclusion is inferred' in discovery['scope']
 for k,rel in [('KCUB42','raw13-target-field/validation.json'),('KCUB43','quadratic-target-affine8/validation.json')]:assert json.loads((B/rel).read_text())['proof_sha256']==sha(Q/E[k]['canonical_source'])
 parent=B/'quadratic-target-affine8/integration-validation.json';PI=json.loads(parent.read_text());assert PI['status']=='pass' and PI['main_receipt_sha256']==sha(B/'quadratic-target-affine8/validation.json')
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB44' in l['claims']];research.validate_lessons(affected,index['entries'],check_reviews=True)
 a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC');assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims'])
 assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']) and all((Q/rel).is_file() for rel in a['sources'])
 assert 'frontier is empty' in a['next_step'] and 'negative-axis source shear' in a['next_step'] and 'weights(-3,0,3)' in a['next_step'] and 'transformed degree9 theorem' in a['next_step']
 out=O/'integration-validation.json'
 for rel in ['POSITIVE_SOURCE_QUADRATIC_FIELD_OBSTRUCTION.md','positive-source-quadratic/README.md','RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md','QUADRATIC_TARGET_AFFINE_DEGREE8.md','README.md']:
  f=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if '://' not in link:
    target=(f.parent/link.split('#')[0]).resolve();assert target==out or target.is_file(),link
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 for rel in ['positive_source_quadratic.py','raw13_target_field9.py','verify_positive_source_quadratic.py','check_integration_positive_source_quadratic.py']:
  f=B/rel;py_compile.compile(str(f),doraise=True);assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian/README.md','research/jacobian/cubic-reduction/README.md','research/jacobian/cubic-reduction/EXPLORATION.md'];z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 result={'status':'pass','claim':'KCUB44','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'proof_sha256':V['proof_sha256'],'all_source_model_import_and_parent_pins_current':True,'base_field_degree9_columns_rank_kernel':[220,213,7],'exact_minor_mod_prime':[116435059,1000000007],'full_three_parameter_inverse_and_weight_conjugacy_checked':True,'all_five_Taylor_heads_and_three_leading_order_cases':True,'universal_weighted_degeneration_is_written_not_formalized':True,'nonzero_transformed_low8_field_dimension':2,'zero_quadratic_dimensions':[4,7,7],'no_rational_Keller_repair_through8':True,'transformed_degree9_remains_unclassified':True,'repair_endpoint_degree_interval':[9,10],'existing_degree10_control_is_not_global_improvement':True,'seven_fixed_parameter_pilots_retained_unused_as_universal_proof':True,'dependency_graph_acyclic':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'links_compilation_whitespace':True,'previous_goal_turn_classification':'progress:KCUB43','this_goal_turn_classification':'progress:KCUB44 all positive-axis quadratic source repairs through8 excluded by formal field obstruction','goal_complete':False,'global_positive_bounds_and_unrestricted_profiles_unchanged':True,'next_gate':'Negative-axis quadratic source shear with mixed weights(-3,0,3),including resonance1+3s=0;requires a new argument.','navigation_attempts_current_parent_state':PI['navigation_attempts_current_parent_state'],'navigation_boundary':'Actual current-scope renderer attempts from the freshly rerun parent chain fail on the same two unrelated missing artifacts;no generated status hand edit or reconstruction.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS KCUB44 integration:base9 rank213/220,all quadratic/zero branches and exact Taylor heads;universal low8 dimension2/no Keller repair;all pins,scope boundaries,parents,reviews,work,links,compile,whitespace;global goal open')
if __name__=='__main__':main()
