#!/usr/bin/env python3
"""Narrow metadata/provenance audit after the full mathematical replays."""
from pathlib import Path
import json,hashlib,sys,re,subprocess,py_compile
B=Path(__file__).resolve().parent;Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 idx=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in idx['entries']};n=E['KCUB37'];vp=B/'weighted-124-125/validation.json';V=json.loads(vp.read_text());pp=B/'weight-resonance7/validation.json';P=json.loads(pp.read_text());assert V['status']==P['status']=='pass' and V['claim']=='KCUB37' and not V['goal_complete'];assert render_status.REQUIRED_FIELDS<=set(n)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS;assert n['artifact_hash']=='sha256:'+sha(Q/n['checker']) and V['checker_sha256']==sha(Q/n['checker']) and V['proof_sha256']==sha(Q/n['canonical_source']);assert n['dependencies']==['KCUB36'] and n['narrows_problems']==['OP-GLDC'] and 'KCUB37' in E['OP-GLDC']['dependencies']
 for k in ['KCUB36','KCUB37']:
  e=E[k];assert e['artifact_hash']=='sha256:'+sha(Q/e['checker'])
  for rel in [e['canonical_source'],e['checker']]+e['software_lock']:assert (Q/rel).is_file(),rel
 for a in E['OP-GLDC']['forbidden_attack_classes']:assert set(a['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 visited=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:'),k;return
  assert k not in stack,(k,stack)
  if k in visited:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  visited.add(k)
 visit('KCUB37',[])
 for k,rel in [('producer_sha256','weighted124125.py'),('arithmetic_sha256','verify.py'),('checkpoint_sha256','weighted-124-125/checkpoint.json'),('structural_identities_sha256','weighted-124-125/structural-identities.json')]:assert V[k]==sha(B/rel)
 for name,h in V['model_hashes'].items():assert h==sha(B/'weighted-124-125'/name)
 keys=['source_atoms','native_coefficient_generators','complete_graph_transformations','univariate_coefficients','power_relation_numerators','integer_boundary_power_subtraction_leading_cases','graph_coefficients','affine_vectorfield_eigenvalues','constant_pivot_determinant','one_variable_unit_multiplier_terms'];assert [V[k] for k in keys]==[24,30,2,3,2,4,4,15,-2,4]
 assert V['weighted_square_shear_identity'] and V['quadratic_invariant_flow_checked'] and len(V['generic_inverse_charts'])==5 and V['old_timeouts_preserved'] and V['no_external_theorem_imported_for_two_weight_closures'];assert len(V['remaining_diagonal_weight_improvement_types'])==12 and {tuple(z['mixed_positive_negative_magnitudes']) for z in V['remaining_diagonal_weight_improvement_types'] if z['mixed_positive_negative_magnitudes']}=={(3,1,2),(4,1,2),(1,1,1)}
 assert V['frontier_KCUB36_receipt_sha256']==sha(pp) and P['proof_sha256']==sha(Q/E['KCUB36']['canonical_source']) and P['checker_sha256']==sha(Q/E['KCUB36']['checker']);assert V['old_four_timeout_checkpoint_sha256']==P['coefficient_checkpoint_sha256']==sha(B/'weight-resonance7/mixed-lift-checkpoint.json')
 for name,h in P['model_hashes'].items():assert h==sha(B/'weight-resonance7/mixed-cases'/(name+'.json'))
 for name,h in P['certificate_hashes'].items():assert h==sha(B/'weight-resonance7/mixed-cases'/(name+'-certificate.json'))
 for k,pin in P['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['canonical_sha256']==sha(Q/E[k]['canonical_source'])
 ls=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in ls if 'OP-GLDC' in l['claims'] or 'KCUB36' in l['claims'] or 'KCUB37' in l['claims']];research.validate_lessons(affected,idx['entries'],check_reviews=True)
 a=next(e for e in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if e['id']=='WORK-GLDC');assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims']);assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']);assert all((Q/p).is_file() for p in a['sources']);assert 'Prioritize(1,-2,-4)' not in a['next_step']
 for rel in ['WEIGHTED_124_125_EXCLUSION.md','WEIGHT_RESONANCE_CENSUS7.md','weighted-124-125/README.md','weight-resonance7/README.md']:
  p=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
   if '://' not in link:assert (p.parent/link.split('#')[0]).is_file(),link
  assert not any(line.rstrip()!=line for line in p.read_text().splitlines())
 for rel in ['weighted124125.py','verify_weighted124125.py','check_integration_weighted124125.py']:
  p=B/rel;py_compile.compile(str(p),doraise=True);assert not any(line.rstrip()!=line for line in p.read_text().splitlines()),rel
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian-followup/README.md','research/jacobian-followup/cubic-reduction/README.md','research/jacobian-followup/cubic-reduction/EXPLORATION.md'];z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 print('PASS ledger/hash/dependency acyclicity;24 atoms/30 equations;graph/univariate/unit/inverse receipts;parent447 certificate pins;four timeout checkpoint unchanged;12-type frontier;',len(affected),'method reviews;work fingerprints;links/compile/whitespace')
 pip=B/'weight-resonance7/integration-validation.json';parent=json.loads(pip.read_text());parent['checker_receipt_sha256']=sha(pp);parent['proof_sha256']=sha(Q/E['KCUB36']['canonical_source']);parent['subsequent_notice_refresh']={'claim':'KCUB37','parent_exact_checker_replayed':True,'all_parent_model_certificate_and_import_pins_checked':True,'old_four_native_timeouts_unchanged':True,'KCUB37_stage_frontier_types':12,'affected_method_reviews':len(affected)};pip.write_text(json.dumps(parent,indent=2)+'\n')
 nav=[]
 for cmd in [[str(R/'.venv/bin/python'),'research/scripts/render_status.py'],[str(R/'.venv/bin/python'),'research/scripts/research.py','render']]:
  z=subprocess.run(cmd,cwd=R,capture_output=True,text=True);out=(z.stdout+z.stderr).strip();nav.append({'command':cmd[1:],'exit_code':z.returncode,'output':out});print('NAV',cmd[1:],'exit',z.returncode,out.splitlines()[-1] if out else '(no output)')
 receipt={'status':'pass','claim':'KCUB37','integration_checker_sha256':sha(Path(__file__)),'checker_receipt_sha256':sha(vp),'proof_sha256':sha(Q/n['canonical_source']),'checker_sha256':sha(Q/n['checker']),'counts':{k:V[k] for k in keys},'generic_inverse_charts':5,'quadratic_invariant_flow_checked':True,'parent_receipt_and_all447_certificate_pins_current':True,'parent_integration_notice_hashes_refreshed':True,'old_native_timeout_checkpoint_unchanged_sha256':V['old_four_timeout_checkpoint_sha256'],'KCUB37_stage_frontier_types':12,'affected_method_reviews':len(affected),'ledger_fields_paths_hashes_dependencies_acyclic':True,'external_dependency_labels_treated_as_terminal':True,'work_fingerprints_sources_current':True,'links_compile_whitespace':True,'previous_goal_turn_classification':'progress:KCUB36 complete weight-ray census and fourteen-type stage','this_goal_turn_progress':'KCUB37 complete124/125 automorphism families and current twelve-type frontier','goal_complete':False,'no_positive_bound_or_unrestricted_profile_change':True,'navigation_attempts':nav,'navigation_boundary':'Actual command outcomes recorded; no generated status hand edit or unrelated artifact reconstruction.'};(B/'weighted-124-125/integration-validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
