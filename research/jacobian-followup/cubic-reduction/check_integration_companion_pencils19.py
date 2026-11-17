#!/usr/bin/env python3
"""Scoped conditional-tensor claims and unchanged-global-bound integration."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB21'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['proof_type']=='hybrid' and e['dependencies']==['KCUB3','KCUB5','KCUB20'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 vp=HERE/'companion-pencils19-validation.json';v=load(vp);assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_companion_pencils19.py') and v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py') and v['base_witness_sha256']==sha(HERE/'cubic12-integer-witness.json') and v['base_validation_sha256']==sha(HERE/'integer-validation12.json');assert (v['homogeneous_companion_dimension'],v['nonzero_cubic_rank'],v['shared_cube_lower_bound'],v['best_verified_cube_identity'],v['global_cubic_linear_upper_bound'])==(19,6,37,98,89);assert not v['global_bound_improved'] and not v['full_19_variable_determinant_expanded'];assert v['all_polynomial_and_rank_identities_replayed'] and v['polynomial_substitution_domain_guarded']
 for pencil,tag,producer,search,parameter,uniform in zip(v['pencils'],['companion-pencil19','companion-w19'],['prepare_companion_pencil19.py','prepare_companion_w19.py'],['search_companion_pencil19.py','search_companion_w19.py'],['-6','-3/2'],[100,101]):
  folder=HERE/tag;inp=load(folder/'input.json');assert inp['source_sha256']==sha(HERE/producer) and pencil['input_sha256']==sha(folder/'input.json');c=load(folder/'best-identity.json');assert c['source_sha256']==sha(HERE/search) and pencil['cube_identity_sha256']==sha(folder/'best-identity.json') and pencil['identity_parameter']==c['sigma']==parameter and c['size']==98;assert pencil['coefficient_matrix_rank']==19 and pencil['cube_identity_all_coefficients_checked'] and pencil['homogeneous_companion_collision_checked'];rank=pencil['derivative_rank'];assert rank['certificate_sha256']==sha(folder/'derivative-rank-certificate.json') and rank['rank_nonzero_parameter']==37 and rank['rank_at_zero']==35 and rank['unit_pivots']==37 and rank['full_residual_zero'];assert pencil['uniform_identity']['uniform_cube_count']==uniform
 for tag,script,forms,rank in [('companion-pencil19','search_companion_pencil19.py',277,165),('companion-pencil19-mixed','search_companion_pencil19_mixed.py',427,171),('companion-w19','search_companion_w19.py',271,163)]:
  folder=HERE/tag;d=load(folder/'input.json');q=load(folder/'search-checkpoint.json');assert len(d['forms'])==forms and d['rank']==rank and q['source_sha256']==sha(HERE/script) and q['input_sha256']==sha(folder/'input.json') and q['best_identity_sha256']==sha(folder/'best-identity.json') and q['status']=='bounded-no-improvement-below89' and q['best']==98 and q['time_budget_seconds']==40
 smt=load(HERE/'companion-pencil19/smt-checkpoint.json');assert smt['source_sha256']==sha(HERE/smt['executed_source_snapshot']) and smt['status']=='complete-with-unknown' and [(r['bound'],r['status']) for r in smt['tests']]==[(88,'unknown'),(97,'unknown')];assert smt['timeout_ms_per_query']==20000 and smt['input_sha256']==sha(HERE/'companion-pencil19/input.json')
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   if 'KCUB21' in item.get('reviewed_claims',{}):assert item['reviewed_claims']['KCUB21']==claim_fingerprint(e)
 for name in ['CHANGED_COMPANION_PENCILS.md','FULL_CUBIC_CARRIER_CLASSIFICATION.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB21'],validation_sha256=sha(vp),schema=True,source_hashes=True,links=True,reviewed_fingerprints=True,conditional_tensor_bounds_only=True,global_upper_bound_still89=True,dominated98_matrix_not_reported_as_improvement=True,three_bounded_pools_preserved=True,both_smt_unknowns_preserved=True,finite_pool_optimum_proved=False,global_minima_closed=False,umbrella_status='open');(HERE/'companion-pencils19-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
