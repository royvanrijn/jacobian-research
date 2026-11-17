#!/usr/bin/env python3
"""Current artifact,scope and metadata audit for the sharp restricted result."""
from pathlib import Path
import json,hashlib,sys,re,subprocess,py_compile
B=Path(__file__).resolve().parent;Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in index['entries']};e=E['KCUB38'];vp=B/'weighted-312-412/validation.json';V=json.loads(vp.read_text());assert V['status']=='pass' and V['claim']=='KCUB38' and not V['goal_complete'];assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS;assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker']) and V['proof_sha256']==sha(Q/e['canonical_source']);assert set(e['dependencies'])==set(V['imports']) and e['narrows_problems']==['OP-GLDC'] and 'KCUB38' in E['OP-GLDC']['dependencies']
 for attack in E['OP-GLDC']['forbidden_attack_classes']:assert set(attack['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB36','KCUB37','KCUB38']:
  entry=E[k];assert entry['artifact_hash']=='sha256:'+sha(Q/entry['checker'])
  for rel in [entry['canonical_source'],entry['checker']]+entry['software_lock']:assert (Q/rel).is_file(),rel
 done=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in done:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  done.add(k)
 visit('KCUB38',[])
 for k,rel in [('model_producer_sha256','weighted3412.py'),('strata_producer_sha256','weighted312_strata.py'),('boundary_producer_sha256','weighted312_boundary.py'),('arithmetic_sha256','verify.py'),('checkpoint_sha256','weighted-312-412/checkpoint.json'),('structural_identities_sha256','weighted-312-412/structural-identities.json')]:assert V[k]==sha(B/rel)
 for name,h in V['model_hashes'].items():assert h==sha(B/'weighted-312-412'/name)
 for name,h in V['stratum_model_hashes'].items():assert h==sha(B/'weighted-312-412'/(name+'.json'))
 for name,h in V['unit_hashes'].items():assert h==sha(B/'weighted-312-412'/(name+'-unit.json'))
 for k,pin in V['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 keys=['source_atoms','native_coefficient_generators','boundary_factorizations','p4_diagonal_identities','p4_active_support_cases','complete_normalization_strata','mu1_complete_boundary_branches','all_constant_pivots_verified','successful_unit_constant_pivots','exact_unit_identities','unit_generators','unit_multiplier_monomials','max_multiplier_degree'];assert [V[k] for k in keys]==[29,52,2,3,7,4,4,48,42,7,123,1096,7]
 assert len(V['changed_reduced_timeouts_preserved'])==1 and len(V['remaining_diagonal_weight_improvement_types'])==10 and V['p3_Lzero_product_derivative'];assert V['minimum_dimension_for_this_additional_form']==3 and V['dimension_lower_bound_is_written_rank_one_consequence'];U=V['matching_existing_upper_control'];assert U['degree']==7 and U['nonlinear_counts']==[2,5,6] and U['weight']==[1,-1,-2] and U['full_determinant_one'] and not U['new_witness'] and U['source_sha256']==sha(B/'raw13-specialization/input.json')
 p37=B/'weighted-124-125/validation.json';P=json.loads(p37.read_text());assert V['parent_frontier_receipt_sha256']==sha(p37) and P['proof_sha256']==sha(Q/E['KCUB37']['canonical_source']);p36=B/'weight-resonance7/validation.json';T=json.loads(p36.read_text());assert P['frontier_KCUB36_receipt_sha256']==sha(p36) and T['proof_sha256']==sha(Q/E['KCUB36']['canonical_source']);assert V['old_native_timeout_checkpoint_sha256']==P['old_four_timeout_checkpoint_sha256']==T['coefficient_checkpoint_sha256']==sha(B/'weight-resonance7/mixed-lift-checkpoint.json')
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB38' in l['claims']];research.validate_lessons(affected,index['entries'],check_reviews=True);a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC');assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims']);assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']);assert all((Q/p).is_file() for p in a['sources']);assert 'remaining15/14-slot' not in a['next_step']  # ten types is the frozen KCUB38 stage, not the current frontier
 for rel in ['DISTINCT_NONZERO_WEIGHT_MINIMUM.md','weighted-312-412/README.md','WEIGHTED_124_125_EXCLUSION.md','WEIGHT_RESONANCE_CENSUS7.md']:
  p=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
   if '://' not in link:assert (p.parent/link.split('#')[0]).is_file(),link
  assert not any(line.rstrip()!=line for line in p.read_text().splitlines())
 for rel in ['weighted3412.py','weighted312_strata.py','weighted312_boundary.py','verify_weighted3412.py','check_integration_weighted3412.py']:
  p=B/rel;py_compile.compile(str(p),doraise=True);assert not any(line.rstrip()!=line for line in p.read_text().splitlines()),rel
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian-followup/README.md','research/jacobian-followup/cubic-reduction/README.md','research/jacobian-followup/cubic-reduction/EXPLORATION.md'];z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 # Parent audits were rerun after their subsequent-result notices; preserve the old stage counts.
 parent_path=B/'weighted-124-125/integration-validation.json';parent=json.loads(parent_path.read_text());assert parent['status']=='pass' and parent['checker_receipt_sha256']==sha(p37) and parent['proof_sha256']==sha(Q/E['KCUB37']['canonical_source']);parent['subsequent_KCUB38_notice_refresh']={'main_checker_passed':True,'current_frontier_types':10,'sharp_additional_form_dimension':3,'affected_method_reviews':len(affected)};parent_path.write_text(json.dumps(parent,indent=2)+'\n')
 receipt={'status':'pass','claim':'KCUB38','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'checker_sha256':sha(Q/e['checker']),'proof_sha256':sha(Q/e['canonical_source']),'counts':{k:V[k] for k in keys},'matching_existing_upper_control':U,'minimum_dimension_for_this_additional_form':3,'written_dimension_lower_bound_via_KCUB36':True,'remaining_types':10,'all_main_model_unit_import_pins_current':True,'native_and_changed_reduced_timeouts_preserved':True,'changed_strata_checkpoint_sha256':sha(B/'weighted-312-412/strata-checkpoint.json'),'new_boundary_checkpoint_sha256':sha(B/'weighted-312-412/boundary-checkpoint.json'),'parent_stage_receipts_refreshed':True,'dependency_graph_acyclic':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'links_compilation_whitespace':True,'previous_goal_turn_classification':'progress:KCUB37 complete twelve-slot closures','this_goal_turn_progress':'KCUB38 final distinct-nonzero closures,sharp conditional raw13/degree7/dimension3,and ten-type broader frontier','goal_complete':False,'no_unrestricted_bound_or_profile_change':True,'navigation_attempts_from_current_parent_audit':parent['navigation_attempts'],'navigation_boundary':'Same current-state rendering attempts retained from the refreshed parent audit; no generated status hand edit or unrelated artifact reconstruction.'};(B/'weighted-312-412/integration-validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS KCUB38:current source/unit/import hashes;29 atoms,52 equations,48 pivots,7 units/1096 terms;matching13-term upper;written form-dimension3;10-type frontier;',len(affected),'method reviews;ledger/work/dependencies/links/compile/whitespace;goal open')
if __name__=='__main__':main()
