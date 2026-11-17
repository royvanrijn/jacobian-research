#!/usr/bin/env python3
"""Current KCUB40 evidence/scope audit; no research calculations or solver calls."""
from pathlib import Path
import hashlib,json,re,sys,subprocess,py_compile
B=Path(__file__).resolve().parent;O=B/'equal-negative-full7';Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in index['entries']};e=E['KCUB40'];vp=O/'validation.json';V=json.loads(vp.read_text())
 assert V['status']=='pass' and V['claim']=='KCUB40' and not V['goal_complete']
 assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS
 assert e['state']=='proved' and e['independent_replay'] and not e['formal_verification'] and not e['external_review']
 assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker']) and V['proof_sha256']==sha(Q/e['canonical_source'])
 assert set(e['dependencies'])==set(V['imports'])=={'KCUB23','KCUB39','KCUB38','KCUB36'}
 assert e['narrows_problems']==['OP-GLDC'] and 'KCUB40' in E['OP-GLDC']['dependencies']
 for x in E['OP-GLDC']['forbidden_attack_classes']:assert set(x['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB23','KCUB36','KCUB37','KCUB38','KCUB39','KCUB40']:
  a=E[k];assert a['artifact_hash']=='sha256:'+sha(Q/a['checker'])
  for rel in [a['canonical_source'],a['checker']]+a['software_lock']:assert (Q/rel).is_file(),rel
 visited=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in visited:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  visited.add(k)
 visit('KCUB40',[])
 for k,pin in V['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 assert V['producer_sha256']==sha(B/'equal_negative_full7.py') and V['arithmetic_sha256']==sha(B/'verify.py') and V['checkpoint_sha256']==sha(O/'checkpoint.json')
 C=json.loads((O/'checkpoint.json').read_text());assert C['status']=='complete' and C['source_sha256']==V['producer_sha256'] and C['files']==V['model_hashes'] and C['new_solver_queries']==0
 for name,h in V['model_hashes'].items():assert sha(O/name)==h
 assert [V[k] for k in ['complete_native_atoms','post_leading_parameters','generic_opposite_power_controls','elliptic_normal_monomial_controls','elliptic_distinct_pole_controls']]==[33,21,8,15,66]
 assert V['prime_degree_pairs']==[[7,4],[7,5],[7,6]] and V['full_source_dictionary_checked']
 assert set(V['small_coefficient_obstructions'])=={'graph','hyperbola1','hyperbola2'} and all(v['all_coefficient_identities_exact'] and v['nonzero_constant_jacobian_impossible'] for v in V['small_coefficient_obstructions'].values())
 assert V['laurent_euler_jet_identity_checked'] and V['elliptic_differential_has_order_zero'] and V['vertical_affine_determinant_identity_checked'] and V['all_degree_curve_and_laurent_lemmas_are_written_not_formalized']
 baseline=O/'prior-evidence.json';assert V['old_evidence_sha256']==sha(baseline)
 for rel,h in json.loads(baseline.read_text())['immutable_old_KCUB23_files'].items():assert sha(Q/rel)==h
 old=Q/'jacobian/equal-negative-completion'
 for rel in ['checkpoint.json','cstar-general-II/checkpoint.json']:
  t=json.loads((old/rel).read_text());assert t['status']=='timeout-no-theorem' and t['time_budget_seconds']==40 and not t['unit_marker']
 assert V['old40second_timeouts_preserved'] and V['large_old_primary_certificate_unreplayed_and_unused'] and V['new_solver_queries']==0
 assert V['remaining_degree7_diagonal_weight_raw12_types']==0 and V['only_degree7_diagonal_weight_collision_type']==[1,-1,-2]
 assert [V[k] for k in ['restricted_minimum_degree','restricted_minimum_nonlinear_occurrences_through_degree7','restricted_minimum_collision_dimension_all_degrees']]==[7,13,3]
 U=V['matching_old_upper_control'];assert U['input_sha256']==sha(B/'raw13-specialization/input.json') and U['degree']==7 and U['counts']==[2,5,6] and U['full_determinant_one'] and not U['new_witness']
 # All parent stage proofs/receipts are current; their nonempty frontiers remain stage evidence.
 parents={'KCUB23':old/'validation.json','KCUB36':B/'weight-resonance7/validation.json','KCUB37':B/'weighted-124-125/validation.json','KCUB38':B/'weighted-312-412/validation.json','KCUB39':B/'zero-weight/validation.json'}
 for k,rel in parents.items():assert json.loads(rel.read_text())['proof_sha256']==sha(Q/E[k]['canonical_source'])
 pi=json.loads((B/'zero-weight/integration-validation.json').read_text());assert pi['status']=='pass' and pi['main_receipt_sha256']==sha(parents['KCUB39'])
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB40' in l['claims']];research.validate_lessons(affected,index['entries'],check_reviews=True)
 a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC')
 assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims'])
 assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']) and all((Q/rel).is_file() for rel in a['sources'])
 assert 'frontier is empty' in a['next_step'] and 'could still matter to raw sparsity' not in a['next_step']  # affine sparsification is subsequently closed byKCUB41
 out=O/'integration-validation.json'
 for rel in ['EQUAL_NEGATIVE_FULL_DEGREE7.md','equal-negative-full7/README.md','README.md','ZERO_WEIGHT_CLASSIFICATION.md','DISTINCT_NONZERO_WEIGHT_MINIMUM.md','WEIGHTED_124_125_EXCLUSION.md','WEIGHT_RESONANCE_CENSUS7.md','../equal-negative-completion/PROOF.md']:
  f=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if '://' not in link:
    target=(f.parent/link.split('#')[0]).resolve();assert target==out or target.is_file(),link
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 for rel in ['equal_negative_full7.py','verify_equal_negative_full7.py','check_integration_equal_negative_full7.py']:
  f=B/rel;py_compile.compile(str(f),doraise=True);assert not any(line.rstrip()!=line for line in f.read_text().splitlines())
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian/README.md','research/jacobian/cubic-reduction/README.md','research/jacobian/cubic-reduction/EXPLORATION.md']
 z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 nav=[]
 for args in [['research/scripts/render_status.py'],['research/scripts/research.py','render']]:
  z=subprocess.run([sys.executable]+args,cwd=R,capture_output=True,text=True);nav.append({'command':args,'exit_code':z.returncode,'output':(z.stdout+z.stderr)[-3000:]})
  if z.returncode:
   expected='KDSQ6: missing document consumer artifacts/generated-results/README.md' if len(args)==1 else 'METHOD-EC-CENTRE-POLICY: missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json'
   assert expected in z.stdout+z.stderr
 result={'status':'pass','claim':'KCUB40','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'proof_sha256':V['proof_sha256'],'model_and_import_pins_current':True,'native_atoms':33,'post_leading_parameters':21,'short_coefficient_obstructions':3,'all_old_inputs_scripts_logs_checkpoints_immutable':True,'old_timeout_statuses_unchanged':True,'large_old_primary_certificate_unreplayed_and_unused':True,'new_solver_queries':0,'parent_stage_receipts_current':True,'dependency_graph_acyclic':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'links_compilation_whitespace':True,'current_diagonal_weight_raw12_frontier_types':0,'only_degree7_diagonal_weight_collision_type':[1,-1,-2],'restricted_degree_support_dimension_minima':[7,13,3],'raw_count_scope':'Three-variable degree<=7 counted displayed diagonal weight coordinates;general source linear frames are outside this minimum.','matching_old_upper':U,'previous_goal_turn_classification':'progress:KCUB39','this_goal_turn_classification':'progress:KCUB40 full111 closure and completed diagonal-weight sharp frontier','global_goal_complete':False,'global_optima_still_open':['unrestricted three-variable degree within[5,7]','unrestricted raw nonlinear term minimum below existing13','normalized support within[10,48]','ordinary-cubic dimension within[4,12]','other special-form dimension and cubic/cube-linear support minima'],'no_global_positive_bound_or_unrestricted_profile_change':True,'navigation_attempts':nav,'navigation_boundary':'Actual final-state attempts fail on the two unrelated pre-existing missing artifacts. Generated status is not hand edited;no unrelated artifact reconstruction.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS KCUB40 integration:full33/outer21,three coefficient obstructions,all models/imports/immutable old evidence/parent receipts/metadata/work/links/compile/whitespace;weighted raw<=12 frontier empty;global goal open')
 for n in nav:print('NAV',n['command'],'exit',n['exit_code'],n['output'].splitlines()[-1] if n['output'] else '')
if __name__=='__main__':main()
