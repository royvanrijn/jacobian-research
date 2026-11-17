#!/usr/bin/env python3
"""Scoped authority, proof-source and preserved-failure checks for KCUB11/12."""
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
 for identifier,dep in [('KCUB11','KCUB10'),('KCUB12','KCUB11')]:
  e=known[identifier];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
  assert e['state']=='proved' and e['dependencies']==[dep] and e['programme_status']=='archived'
  assert e['independent_replay'] and not e['formal_verification'] and not e['external_review']
  assert not e['closes_problems'] and not e['narrows_problems'] and e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
  for p in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/p).is_file(),p
 v=load(HERE/'weighted-tail7/validation.json');w=load(HERE/'weighted-newfirst6/validation.json')
 assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_tail7.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py')
 assert v['input_sha256']==sha(HERE/'weighted-tail7/input.json') and v['raw13_validation_sha256']==sha(HERE/'raw13-rigidity/validation.json')
 assert (v['nonlinear_atoms'],v['branch_unit_identities'],v['multiplier_monomials'],v['max_multiplier_degree'],v['unique_nonlinear_total'])==(19,3,15665,7,13)
 assert all(v[k] for k in ['full_three_variable_determinant_reconstructed','necessary_boundary_equations_reconstructed','linear_parameterizations_complete','all_branch_identities_multiplied','top_v_identities_checked','linear_g_coordinate_change_checked','old13_reduction_checked'])
 for case in v['cases']:assert case['certificate_sha256']==sha(HERE/'weighted-unit-certificates'/(case['tag']+'.json'))
 assert w['status']=='pass' and w['source_sha256']==sha(HERE/'verify_weighted_newfirst6.py') and w['utility_sha256']==sha(HERE/'weighted_verify_utils.py')
 assert w['input_sha256']==sha(HERE/'weighted-newfirst6/input.json') and w['tail_classification_sha256']==sha(HERE/'weighted-tail7/validation.json')
 assert w['parameter_certificate_sha256']==sha(HERE/'weighted-newfirst6/parameter-field/certificate.json')
 assert w['nonlinear_atoms']==17 and w['certificate_target']=='theta^7' and w['multiplier_monomials']==51
 assert all(w[k] for k in ['universal_three_variable_determinant_reconstructed','linear_parameterization_complete','cleared_identity_multiplied','exceptional_parameter_roots_exactly_zero','zero_slice_excluded_by_tail_classification'])
 for path,status in [('weighted-tail7/checkpoint.json','timeout'),('weighted-tail7/generic/checkpoint.json','timeout'),('weighted-tail7/nra/checkpoint.json','unknown'),('weighted-newfirst6/nra-checkpoint.json','unknown'),('weighted-newfirst6/groebner-checkpoint.json','timeout'),('weighted-newfirst7/nra-checkpoint.json','unknown')]:assert load(HERE/path)['status']==status
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   for identifier in ['KCUB11','KCUB12']:
    if identifier in item.get('reviewed_claims',{}):assert item['reviewed_claims'][identifier]==claim_fingerprint(known[identifier])
 for name in ['WEIGHTED_TAIL_CLASSIFICATION.md','WEIGHTED_MARKED_DEGREE.md','README.md','RAW13_PATTERN_RIGIDITY.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open'
 receipt=dict(status='pass',claims=['KCUB11','KCUB12'],tail_classification=True,marked_degree_lower_bound=7,marked_degree_attained=True,
  global_minima_closed=False,source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,unknown_timeouts_preserved=True,
  tail_validation_sha256=sha(HERE/'weighted-tail7/validation.json'),degree_validation_sha256=sha(HERE/'weighted-newfirst6/validation.json'),umbrella_status='open')
 (HERE/'weighted-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
