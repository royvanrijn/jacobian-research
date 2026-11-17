#!/usr/bin/env python3
"""Audit KCUB39 evidence, scopes and current work without research queries."""
from pathlib import Path
import hashlib,json,re,sys,subprocess,py_compile
B=Path(__file__).resolve().parent
Q=B.parents[1]
R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text())
 E={e['id']:e for e in index['entries']}
 e=E['KCUB39'];vp=B/'zero-weight/validation.json';V=json.loads(vp.read_text())
 assert V['status']=='pass' and V['claim']=='KCUB39' and not V['goal_complete']
 assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS
 assert e['state']=='proved' and e['independent_replay'] and not e['formal_verification'] and not e['external_review']
 assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker'])
 assert V['proof_sha256']==sha(Q/e['canonical_source'])
 assert set(e['dependencies'])==set(V['imports'])=={'KCUB38','KCUB23','KCUB36'}
 assert e['narrows_problems']==['OP-GLDC'] and 'KCUB39' in E['OP-GLDC']['dependencies']
 for attack in E['OP-GLDC']['forbidden_attack_classes']:
  assert set(attack['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB36','KCUB37','KCUB38','KCUB39']:
  a=E[k];assert a['artifact_hash']=='sha256:'+sha(Q/a['checker'])
  for rel in [a['canonical_source'],a['checker']]+a['software_lock']:
   assert (Q/rel).is_file(),rel
 seen=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in seen:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  seen.add(k)
 visit('KCUB39',[])
 for key,rel in [('producer_sha256','zero_weight.py'),('arithmetic_sha256','verify.py'),('checkpoint_sha256','zero-weight/checkpoint.json')]:
  assert V[key]==sha(B/rel)
 assert len(V['model_hashes'])==9
 for name,h in V['model_hashes'].items():assert h==sha(B/'zero-weight'/name)
 cp=json.loads((B/'zero-weight/checkpoint.json').read_text())
 assert cp['status']=='complete' and cp['source_sha256']==V['producer_sha256']
 for k,pin in V['imports'].items():
  assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 counts=['complete_degree7_models','native_atoms','native_coefficient_generators','boundary_removed_pure_x_atoms','boundary_reduced_atoms','principal_face_specializations','rank_one_leading_degree_cases','tail_leading_degree_cases']
 assert [V[k] for k in counts]==[9,240,438,162,78,12,27,16]
 assert len(V['generic_shear_inverse_controls'])==9 and all(c['both_generic_inverse_compositions'] for c in V['generic_shear_inverse_controls'])
 assert V['same_sign_matrix_unit_identities_checked'] and V['all_degree_one_zero_theorem_has_written_filtration_proof']
 assert V['minimum_collision_dimension_exactly_one_zero_form']==4 and V['dimension_lower_bound_is_written_all_degree_one_zero_and_dimension2_unit_proof']
 U=V['existing_four_dimensional_upper_control']
 assert U['dimension']==4 and U['weights']==[0,1,-1,-2] and U['degree']==7 and U['row_degrees']==[1,4,6,7]
 assert U['nonlinear_counts']==[0,2,5,6] and U['full_determinant_one'] and not U['new_witness']
 assert U['old_raw13_input_sha256']==sha(B/'raw13-specialization/input.json')
 assert len(V['remaining_diagonal_weight_raw12_type'])==1 and V['remaining_diagonal_weight_raw12_type'][0]['weight']==[1,-1,-1]
 assert V['strict_flag_below764_has_no_nontrivial_linear_Gm_symmetry'] and V['strict_corollary_uses_imported_KCUB23_and_known112_flag']
 assert V['old111_timeouts_untouched'] and V['no_new_solver_query'] and V['global_positive_bounds_and_unrestricted_profile_counts_unchanged']
 parent=B/'weighted-312-412/validation.json';P=json.loads(parent.read_text())
 assert V['parent_frontier_receipt_sha256']==sha(parent) and P['proof_sha256']==sha(Q/E['KCUB38']['canonical_source'])
 assert len(P['remaining_diagonal_weight_improvement_types'])==10  # frozen parent stage
 pi=json.loads((B/'weighted-312-412/integration-validation.json').read_text())
 assert pi['status']=='pass' and pi['main_receipt_sha256']==sha(parent)
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons']
 affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB39' in l['claims']]
 research.validate_lessons(affected,index['entries'],check_reviews=True)
 a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC')
 assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims'])
 assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims'])
 assert all((Q/rel).is_file() for rel in a['sources'])
 assert 'proposed next gate,not a completed zero-weight exclusion' not in a['next_step']  # one remaining type is the frozen KCUB39 stage
 out=B/'zero-weight/integration-validation.json'
 for rel in ['ZERO_WEIGHT_CLASSIFICATION.md','zero-weight/README.md','README.md','DISTINCT_NONZERO_WEIGHT_MINIMUM.md','WEIGHTED_124_125_EXCLUSION.md','WEIGHT_RESONANCE_CENSUS7.md']:
  f=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if '://' not in link:
    dest=(f.parent/link.split('#')[0]).resolve()
    assert dest==out or dest.is_file(),link
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 for rel in ['zero_weight.py','verify_zero_weight.py','check_integration_zero_weight.py']:
  f=B/rel;py_compile.compile(str(f),doraise=True)
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian/README.md','research/jacobian/cubic-reduction/README.md','research/jacobian/cubic-reduction/EXPLORATION.md']
 z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True)
 assert z.returncode==0,z.stdout+z.stderr
 old={rel:sha(Q/rel) for rel in E['KCUB23']['software_lock'] if rel.endswith(('checkpoint.json','run.log','degree7.log','input.json'))}
 result={'status':'pass','claim':'KCUB39','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'proof_sha256':V['proof_sha256'],'all_source_model_and_import_pins_current':True,'counts':{k:V[k] for k in counts},'minimum_dimension_exactly_one_zero_form':4,'existing_upper_control':U,'current_diagonal_weight_raw12_frontier_types':1,'only_remaining_weight':[1,-1,-1],'strict_linear_symmetry_corollary_imports_strict_KCUB23_only':True,'old111_artifact_hashes':old,'old_native_timeout_checkpoint_sha256':sha(B/'weight-resonance7/mixed-lift-checkpoint.json'),'old_reduced_timeout_checkpoint_sha256':sha(B/'weighted-312-412/strata-checkpoint.json'),'no_solver_query_or_old_timeout_mutation':True,'dependency_graph_acyclic':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'links_compilation_whitespace':True,'previous_goal_turn_classification':'progress:KCUB38','this_goal_turn_progress':'KCUB39 all-degree one-zero closure,sharp form dimension4,one-type frontier and strict linear-symmetry exclusion','goal_complete':False,'no_global_positive_bound_or_unrestricted_profile_change':True,'navigation_attempts_current_state':pi['navigation_attempts_from_current_parent_audit'],'navigation_boundary':'Actual same-state renderer attempts fail on pre-existing unrelated missing KDSQ6 document consumer and METHOD-EC-CENTRE-POLICY completion;no generated status hand edit or unrelated reconstruction.'}
 out.write_text(json.dumps(result,indent=2)+'\n')
 print('PASS KCUB39 integration:240 atoms/438 generators;12 faces;all pins,imports,dependency graph,method reviews,work,links,compilation and whitespace;sharp form dimension4;only full111 remains;global goal open')
if __name__=='__main__':main()
