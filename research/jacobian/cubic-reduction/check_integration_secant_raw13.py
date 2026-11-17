#!/usr/bin/env python3
"""Scoped authority and source checks for KCUB9 and KCUB10."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']}
 for identifier,dependency in [('KCUB9','KCUB8'),('KCUB10','BCRDB1')]:
  e=known[identifier];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
  assert e['state']=='proved' and e['dependencies']==[dependency] and e['programme_status']=='archived'
  assert e['independent_replay'] and not e['formal_verification'] and not e['external_review']
  assert not e['closes_problems'] and not e['narrows_problems'] and e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
  for p in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/p).is_file(),p
 out=HERE/'power-secant-two';v=load(out/'validation.json')
 assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_power_secant_two.py')
 assert v['primary_checkpoint_sha256']==sha(out/'checkpoint.json') and v['certification_sha256']==sha(out/'certification.json')
 assert v['fixed_span_validation_sha256']==sha(HERE/'power-span-split/validation.json')
 assert (v['source_cliques'],v['source_support_pairs'],v['coordinate_zero_pairs'],v['eligible_proportional_triples'],v['unit_identities'],v['max_multiplier_degree'])==(15,120,7,214,607,12)
 assert all(v[k] for k in ['all_nonedge_mixed_derivatives_zero','all_support_pairs_covered','all_membership_spaces_reconstructed','all_coordinate_maps_verified','all_unit_identities_multiplied'])
 for case in v['cases']:
  if 'certificate_sha256' in case:assert case['certificate_sha256']==sha(out/('pair-%02d-%02d-cert.json'%tuple(case['pair'])))
 raw=HERE/'raw13-rigidity';r=load(raw/'validation.json')
 assert r['status']=='pass' and r['source_sha256']==sha(HERE/'verify_raw13_rigidity.py')
 assert r['raw_input_sha256']==sha(HERE/'raw13-specialization/input.json') and r['normalized_input_sha256']==sha(HERE/'raw13-normalized/input.json')
 assert r['certification_sha256']==sha(raw/'certification-direct.json') and r['identities']==15 and r['baseline_total']==13 and r['baseline_nonlinear_counts']==[2,5,6]
 assert all(r[k] for k in ['universal_determinant_reconstructed','normalized_generators_reconstructed','baseline_collision_checked','all_polynomial_identities_multiplied','laurent_family_determinant_checked','laurent_family_collision_checked'])
 manifest=load(raw/'certification-direct.json')
 for tag,c in manifest['certificates'].items():assert c['sha256']==sha(raw/c['path'])
 old=load(HERE/'raw13-specialization/checkpoint.json');assert old['status']=='inconclusive' and sum(c['status']=='timeout' for c in old['systems'])==7
 for folder,script,status in [('power-three-rank2','power_three_rank2.py','finite-pool-no-improvement'),('power-three-rank1','power_three_rank1.py','finite-pool-no-improvement'),('power-tangent-search','power_tangent_search.py','completed-directions-no-improvement')]:
  data=load(HERE/folder/'checkpoint.json');assert data['source_sha256']==sha(HERE/script) and data['status']==status and data['best']==89
 assert load(HERE/'power-three-rank2/checkpoint.json')['triples_tested']==50138
 assert load(HERE/'power-three-rank1/checkpoint.json')['tasks']==123379
 assert load(HERE/'power-tangent-search/checkpoint.json')['directions']==930
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   for identifier in ['KCUB9','KCUB10']:
    if identifier in item.get('reviewed_claims',{}):assert item['reviewed_claims'][identifier]==claim_fingerprint(known[identifier])
 for name in ['TWO_CUBE_OBSTRUCTION.md','RAW13_PATTERN_RIGIDITY.md','README.md','SPAN_OPTIMALITY.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open'
 receipt=dict(status='pass',claims=['KCUB9','KCUB10'],secant_unit_identities=607,raw_pattern_identities=15,
  source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,global_minimum_claimed=False,umbrella_status='open',
  secant_validation_sha256=sha(out/'validation.json'),raw_validation_sha256=sha(raw/'validation.json'))
 (HERE/'secant-raw13-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
