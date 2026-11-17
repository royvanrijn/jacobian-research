#!/usr/bin/env python3
"""Audit the full rational-target field and first nonlinear layer evidence/scopes."""
from pathlib import Path
import json,hashlib,re,sys,subprocess,py_compile
B=Path(__file__).resolve().parent;O=B/'raw13-target-field';Q=B.parents[1];R=Q.parent
sys.path.insert(0,str(Q/'scripts'))
import research,render_status
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 index=json.loads((Q/'MATH_STATUS.json').read_text());E={e['id']:e for e in index['entries']};e=E['KCUB42'];vp=O/'validation.json';V=json.loads(vp.read_text())
 assert V['status']=='pass' and V['claim']=='KCUB42' and not V['goal_complete']
 assert render_status.REQUIRED_FIELDS<=set(e)<=render_status.REQUIRED_FIELDS|render_status.OPTIONAL_FIELDS
 assert e['state']=='proved' and e['independent_replay'] and not e['formal_verification'] and not e['external_review']
 assert e['artifact_hash']=='sha256:'+sha(Q/e['checker']) and V['checker_sha256']==sha(Q/e['checker']) and V['proof_sha256']==sha(Q/e['canonical_source'])
 assert set(e['dependencies'])==set(V['imports'])=={'KCUB41'} and e['narrows_problems']==['OP-GLDC'] and 'KCUB42' in E['OP-GLDC']['dependencies']
 for attack in E['OP-GLDC']['forbidden_attack_classes']:assert set(attack['witnesses'])<=set(E['OP-GLDC']['dependencies'])
 for k in ['KCUB41','KCUB42']:
  a=E[k];assert a['artifact_hash']=='sha256:'+sha(Q/a['checker'])
  for rel in [a['canonical_source'],a['checker']]+a['software_lock']:assert (Q/rel).is_file(),rel
 seen=set()
 def visit(k,stack):
  if k not in E:assert k.startswith('external:');return
  assert k not in stack
  if k in seen:return
  for j in E[k]['dependencies']:visit(j,stack+[k])
  seen.add(k)
 visit('KCUB42',[])
 for k,pin in V['imports'].items():assert pin['scope_fingerprint']==research.claim_fingerprint(E[k]) and pin['proof_sha256']==sha(Q/E[k]['canonical_source'])
 for k,rel in [('producer_sha256','raw13_target_field.py'),('arithmetic_sha256','verify.py'),('input_sha256','raw13-target-field/input.json'),('checkpoint_sha256','raw13-target-field/checkpoint.json'),('discovery_sha256','raw13-target-field/discovery.json'),('degree8_input_sha256','raw13-target-field/degree8/input.json'),('degree8_checkpoint_sha256','raw13-target-field/degree8/checkpoint.json'),('known_degree8_sha256','raw13-target-field/known-degree8.json')]:assert V[k]==sha(B/rel)
 for prefix,producer,columns,rank,points in [('', 'raw13_target_field.py',120,116,58),('degree8/','raw13_target_field8.py',165,158,79)]:
  cp=O/(prefix+'checkpoint.json');C=json.loads(cp.read_text());D=json.loads((O/(prefix+'input.json')).read_text());assert C['status']=='complete' and C['source_sha256']==D['source_sha256']==sha(B/producer) and C['input_sha256']==sha(O/(prefix+'input.json'))
  assert C['columns']==columns and C['rank']==rank and C['points']==points and D['expected_rank']==rank
 assert V['full_generic_cubic_and_all_source_recovery_identities_checked'] and V['irreducibility_by_written_primitive_linear_target_argument'] and V['source_degree_caps']==[7,8]
 assert V['target_degree_cap'] is None and V['target_denominator_cap'] is None and V['source_monomial_columns']==120 and V['selected_minor_size']==116 and V['exact_rational_target_specialization_points']==58
 assert V['modulus']==1000000007 and V['prime_checked_by_trial_division'] and V['nonzero_minor_determinant_mod_prime']==866695135 and V['all_specialization_denominators_nonzero_mod_prime'] and V['kernel_dimension_exact']==4 and V['kernel_basis']==['1','H1','H2','H3']
 D8=V['degree8_field_rank'];assert [D8[k] for k in ['points','minor_size','columns','determinant_mod_prime']]==[79,158,165,224657276] and len(D8['E_determinants'])==79 and all(D8['E_determinants'])
 assert V['degree8_kernel_dimension']==7 and V['degree8_kernel_basis']==['1','H1','H2','H3','H1^2','S','T']
 assert V['all_rational_targets_with_source_degree_le7_are_affine'] and V['degree7_rigidity_persists_in_arbitrary_affine_source_target_frames'] and V['all_target_only_compositions_preserving_degree7_are_affine'] and V['minimum_source_degree_for_nonaffine_rational_target_with_polynomial_pullback']==8
 L=V['degree8_target_locus'];assert [L[k] for k in ['parameters','native_determinant_generators','linear_memberships','multiplier_terms','max_multiplier_degree']]==[9,27,7,17,1]
 assert L['normal_form']=='Y+(0,u1,u2)*Y1^2' and L['all_memberships_multiplied_exactly'] and L['generic_two_sided_inverse_checked'] and L['fixed_source_nonlinear_postcomposition_minimum']==19 and L['fixed_source_bound_uses13_old_and6_disjoint_new_source_columns']
 loc=O/'target-locus';assert L['input_sha256']==sha(loc/'input.json') and L['certificate_checkpoint_sha256']==sha(loc/'certificate-checkpoint.json') and L['certificate_log_sha256']==sha(loc/'certificate.log')
 for cpname,script,log,status in [('checkpoint.json','model.sing','run.log','primary-basis'),('certificate-checkpoint.json','certify.sing','certificate.log','complete-primary-certificates')]:
  C=json.loads((loc/cpname).read_text());assert C['status']==status and C['exit_code']==0 and C['script_sha256']==sha(loc/script) and C['log_sha256']==sha(loc/log) and C['input_sha256']==sha(loc/'input.json')
 U=V['matching_nonlinear_target_shear_control'];assert U['degree']==8 and U['nonlinear_counts']==[2,11,6] and U['full_determinant_one'] and not U['global_improvement']
 parent=B/'raw13-affine-jet/validation.json';P=json.loads(parent.read_text());assert P['status']=='pass' and P['proof_sha256']==sha(Q/E['KCUB41']['canonical_source'])
 pi=json.loads((B/'raw13-affine-jet/integration-validation.json').read_text());assert pi['status']=='pass' and pi['main_receipt_sha256']==sha(parent)
 lessons=json.loads((Q/'knowledge/lessons.json').read_text())['lessons'];affected=[l for l in lessons if 'OP-GLDC' in l['claims'] or 'KCUB42' in l['claims']];research.validate_lessons(affected,index['entries'],check_reviews=True)
 a=next(a for a in json.loads((Q/'knowledge/work_items.json').read_text())['actions'] if a['id']=='WORK-GLDC');assert len(a['claims'])==len(set(a['claims'])) and len(a['sources'])==len(set(a['sources'])) and set(a['reviewed_claims'])==set(a['claims'])
 assert all(a['reviewed_claims'][k]==research.claim_fingerprint(E[k]) for k in a['claims']) and all((Q/rel).is_file() for rel in a['sources'])
 assert 'frontier is empty' in a['next_step']  # degree8 affine-source gate is subsequently closed byKCUB43
 out=O/'integration-validation.json'
 for rel in ['RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md','raw13-target-field/README.md','RAW13_AFFINE_ORBIT_MINIMUM.md','README.md']:
  f=B/rel
  for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if '://' not in link:
    target=(f.parent/link.split('#')[0]).resolve();assert target==out or target.is_file(),link
  assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 for rel in ['raw13_target_field.py','raw13_target_field8.py','raw13_target_locus.py','verify_raw13_target_field.py','check_integration_raw13_target_field.py']:
  f=B/rel;py_compile.compile(str(f),doraise=True);assert not any(line.rstrip()!=line for line in f.read_text().splitlines()),rel
 paths=['research/MATH_STATUS.json','research/knowledge/lessons.json','research/knowledge/work_items.json','research/jacobian/README.md','research/jacobian/cubic-reduction/README.md','research/jacobian/cubic-reduction/EXPLORATION.md'];z=subprocess.run(['git','diff','--check','--']+paths,cwd=R,capture_output=True,text=True);assert z.returncode==0,z.stdout+z.stderr
 result={'status':'pass','claim':'KCUB42','integration_checker_sha256':sha(Path(__file__)),'main_receipt_sha256':sha(vp),'proof_sha256':V['proof_sha256'],'all_rank_locus_import_and_parent_pins_current':True,'source_degree7_field_kernel_dimension':4,'source_degree8_field_kernel_dimension':7,'independent_rank_minors':[[116,866695135],[158,224657276]],'prime':1000000007,'uncapped_rational_target_degree_and_denominators':True,'complete_target_locus':[9,27,7,17,1],'only_normalized_degree8_Keller_targets':'Y+(0,lambda,mu)*Y1^2','fixed_source_nonlinear_minimum':19,'source_affine_changes_of_new_degree8_family_open_at_KCUB42_stage':True,'matching_control_global_improvement':False,'primary_jobs_terminal_and_certificates_replayed':True,'dependency_graph_acyclic':True,'affected_method_reviews':len(affected),'work_fingerprints_current':True,'links_compilation_whitespace':True,'previous_goal_turn_classification':'progress:KCUB41','this_goal_turn_classification':'progress:KCUB42 uncapped full target-field rigidity and complete first nonlinear layer','goal_complete':False,'no_global_positive_bound_or_unrestricted_profile_change':True,'next_gate':'Joint affine source/output sparsification of H+(0,lambda,mu)*H1^2,with full parameter strata;or genuinely nonlinear source coordination/changed construction.','navigation_attempts_current_parent_state':pi['navigation_attempts_current_parent_state'],'navigation_boundary':'Same final mathematical/metadata scopes as the freshly rerun parent audits;two unrelated missing artifacts block rendering. No generated-status hand edit or reconstruction.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS KCUB42 integration:full rational target field116/158 minors,kernels4/7,complete degree8 locus27 generators/7 memberships17 terms;all pins,scopes,reviews,work,links,compile,whitespace;fixed-source19 only;global goal open')
if __name__=='__main__':main()
