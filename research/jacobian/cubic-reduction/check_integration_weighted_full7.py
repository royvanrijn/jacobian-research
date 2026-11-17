#!/usr/bin/env python3
"""Narrow authority/evidence checks for KCUB13 and its retained failures."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={r['id']:r for r in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB13'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['dependencies']==['KCUB11'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for p in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/p).is_file(),p
 v=load(HERE/'weighted-newfirst7/validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_newfirst7.py') and v['input_sha256']==sha(HERE/'weighted-newfirst7/input.json');assert v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py') and v['tail_classification_sha256']==sha(HERE/'weighted-tail7/validation.json');assert (v['nonlinear_atoms'],v['coefficient_parameters'],v['free_parameters'],v['multiplier_monomials'],v['max_multiplier_degree'],v['unique_nonlinear_total'])==(23,21,19,1834,15,13);assert v['cleared_parameter_identities']==4 and v['exceptional_unit_identities']==1
 for key in ['universal_three_variable_determinant_reconstructed','quotient_jacobian_checked','complete_linear_parameterization_checked','completed_square_identity_checked','boundary_charts_and_wronskians_checked','boundary_constraints_reconstructed','all_used_certificates_multiplied','all_exceptional_branches_covered','new_first_coefficients_forced_zero','zero_slice_reduced_to_KCUB11']:assert v[key]
 for case in v['cases']:
  folder=HERE/'weighted-newfirst7/parameter-field'/case['mode'];assert case['certificate_sha256']==sha(folder/'certificate.json');ip=folder/('input-v0.json' if case['mode']=='parabolic-boundary' else 'input.json');assert case['input_sha256']==sha(ip)
 assert load(HERE/'weighted-newfirst7/nra-checkpoint.json')['status']=='unknown'
 for mode in ['cubicgraph','parabolic']:assert load(HERE/'weighted-newfirst7/parameter-field'/mode/'checkpoint.json')['status']=='timeout'
 assert load(HERE/'weighted-newfirst7/parameter-field/parabolic-boundary/error-v0.json')['status']=='renderer-error'
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   for identifier in ['KCUB12','KCUB13']:
    if identifier in item.get('reviewed_claims',{}):assert item['reviewed_claims'][identifier]==claim_fingerprint(known[identifier])
 for name in ['WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md','WEIGHTED_MARKED_DEGREE.md','WEIGHTED_TAIL_CLASSIFICATION.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open'
 receipt=dict(status='pass',claims=['KCUB13'],validation_sha256=sha(HERE/'weighted-newfirst7/validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,unknown_timeouts_preserved=True,geometric_proof_formalized=False,marked_class_minimum_terms=13,global_minima_closed=False,umbrella_status='open');(HERE/'weighted-full7-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
