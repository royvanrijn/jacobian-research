#!/usr/bin/env python3
"""Audit KCUB41 source pins,complete jet evidence and current repository scope."""
from pathlib import Path
import json,hashlib,re,sys,subprocess,py_compile
B=Path(__file__).resolve().parent;O=B/'raw13-affine-jet';Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in index['entries']};e=E['KCUB41'];vp=O/'validation.json';V=json.loads(vp.read_text())
 assert V['status']=='pass' and V['claim']=='KCUB41' and not V['goal_complete']
 assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS
 assert e['state']=='proved' and e['independent_replay'] and not e['formal_verification'] and not e['external_review']
 assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker']) and V['proof_sha256']==sha(Q/e['canonical_source'])
 assert set(e['dependencies'])==set(V['imports'])=={'KCUB40'} and e['narrows_problems']==['OP-GLDC'] and 'KCUB41' in E['OP-GLDC']['dependencies']
 for x in E['OP-GLDC']['forbidden_attack_classes']:assert set(x['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB40','KCUB41']:
  a=E[k];assert a['artifact_hash']=='sha256:'+sha(Q/a['checker'])
  for rel in [a['canonical_source'],a['checker']]+a['software_lock']:assert (Q/rel).is_file(),rel
 done=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in done:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  done.add(k)
 visit('KCUB41',[])
 for k,pin in V['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 for k,rel in [('producer_sha256','raw13_affine_jet.py'),('arithmetic_sha256','verify.py'),('input_sha256','raw13-affine-jet/input.json'),('checkpoint_sha256','raw13-affine-jet/checkpoint.json'),('discovery_sha256','raw13-affine-jet/discovery.json')]:assert V[k]==sha(B/rel)
 C=json.loads((O/'checkpoint.json').read_text());assert C['status']=='complete' and C['source_sha256']==V['producer_sha256'] and C['input_sha256']==V['input_sha256'] and C['discovery_sha256']==V['discovery_sha256']
 assert V['complete_homogeneous_degrees']==list(range(2,8)) and V['complete_jet_matrix_entries']==348 and V['independent_full_translation_not_binomial_entry_production'] and V['all_weight_covariance_monomials_checked']
 assert V['unit_multiplier_terms']==31 and V['max_multiplier_degree']==4 and V['selected_minor_units']=={'degree2':{'generators':4,'multiplier_terms':20,'identity_one':True},'degree3':{'generators':3,'multiplier_terms':11,'identity_one':True}}
 assert V['high_degree_and_zero_branch_minors_checked'] and V['exact_rank_profile_on_c_axis']==[2,3,3,2,2,1] and V['exact_rank_profile_off_c_axis']==[3,3,3,2,2,1]
 assert V['affine_nonlinear_occurrence_lower_bound']==13 and V['off_axis_translation_lower_bound']==14 and V['nonzero_c_axis_translation_source_columns']==18
 assert V['equality_source_centre_must_be_origin'] and V['equality_source_and_target_linear_parts_must_be_monomial_matrices'] and V['universal_rank_invariance_and_UFD_equality_proof_written_not_formalized']
 assert V['corollary_minimum13_in_every_affine_frame_of_every_degree7_diagonal_weight_collision']
 assert V['preliminary_primary_queries']==C['preliminary_primary_queries']==3 and V['all_preliminary_solver_answers_unused_as_proof'] and V['producer_and_checker_solver_queries']==0
 U=V['matching_old_upper'];assert U['source_input_sha256']==sha(B/'raw13-specialization/input.json') and U['counts']==[2,5,6] and U['degree']==7 and U['identity_linear_part'] and U['full_determinant_one'] and not U['new_witness']
 parent=B/'equal-negative-full7/validation.json';P=json.loads(parent.read_text());assert P['status']=='pass' and P['proof_sha256']==sha(Q/E['KCUB40']['canonical_source'])
 pi=json.loads((B/'equal-negative-full7/integration-validation.json').read_text());assert pi['status']=='pass' and pi['main_receipt_sha256']==sha(parent)
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB41' in l['claims']];research.validate_lessons(affected,index['entries'],check_reviews=True)
 a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC');assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims'])
 assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']) and all((Q/rel).is_file() for rel in a['sources'])
 assert 'frontier is empty' in a['next_step'] and 'Do not run an affine raw13 grid' in a['next_step']
 assert '12-variable cubic outside its completed blocks are still open' in a['next_step']
 out=O/'integration-validation.json'
 for rel in ['RAW13_AFFINE_ORBIT_MINIMUM.md','raw13-affine-jet/README.md','EQUAL_NEGATIVE_FULL_DEGREE7.md','README.md']:
  f=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if '://' not in link:
    target=(f.parent/link.split('#')[0]).resolve();assert target==out or target.is_file(),link
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines())
 for rel in ['raw13_affine_jet.py','verify_raw13_affine_jet.py','check_integration_raw13_affine_jet.py']:
  f=B/rel;py_compile.compile(str(f),doraise=True);assert not any(line.rstrip()!=line for line in f.read_text().splitlines())
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian/README.md','research/jacobian/cubic-reduction/README.md','research/jacobian/cubic-reduction/EXPLORATION.md'];z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 result={'status':'pass','claim':'KCUB41','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'proof_sha256':V['proof_sha256'],'all_models_imports_and_parent_receipt_pins_current':True,'complete_jet_entries':348,'unit_generators':7,'unit_multiplier_terms':31,'max_multiplier_degree':4,'exact_centre_rank_profiles_checked':True,'full_affine_orbit_minimum':13,'equality_original_centre_and_monomial_linear_frames':True,'existing_upper_control':U,'three_preliminary_small_solver_answers_preserved_unused':True,'no_solver_in_producer_or_checker':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'dependency_graph_acyclic':True,'links_compilation_whitespace':True,'previous_goal_turn_classification':'progress:KCUB40','this_goal_turn_classification':'progress:KCUB41 complete affine raw13 orbit and equality rigidity','goal_complete':False,'no_global_positive_bound_or_unrestricted_profile_change':True,'next_routes':'Genuinely nonlinear coordinated changes or a changed/unweighted construction. The12-variable cubic affine orbit remains outside this3D closure.','navigation_attempts_current_parent_state':pi['navigation_attempts'],'navigation_boundary':'Same final mathematical/metadata state as the freshly rerun KCUB40 audit;the two unrelated missing artifacts still block renderers. No generated-status hand edit or reconstruction.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS KCUB41 integration:all348 jets,two7-generator/31-term units,complete affine13 minimum and equality rigidity;metadata/imports/parent/work/reviews/links/compile/whitespace;global goal open')
if __name__=='__main__':main()
