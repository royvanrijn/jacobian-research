#!/usr/bin/env python3
"""Scoped ledger, evidence and retained-failure checks for KCUB14."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB14'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['dependencies']==['KCUB13','KCUB11'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 v=load(HERE/'weighted-unmarked7/validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_unmarked7.py') and v['input_sha256']==sha(HERE/'weighted-unmarked7/input.json');assert v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py') and v['tail_classification_sha256']==sha(HERE/'weighted-tail7/validation.json');assert (v['nonlinear_atoms'],v['coefficient_parameters'],v['assumed_collision_equations'],v['multiplier_monomials'],v['max_multiplier_degree'],v['unique_nonlinear_total'])==(23,21,0,30337,15,13);assert v['cleared_parameter_identities']==3 and v['full_ring_unit_identities']==5 and len(v['cases'])==8
 for key in ['universal_three_variable_determinant_reconstructed','quotient_jacobian_checked','all_boundary_constraints_reconstructed','linear_parameterizations_certified_with_allowed_monomial_pivots','all_used_certificates_multiplied','parallel_line_polynomial_identities_checked','all_exceptional_branches_covered','automatic_nonzero_collision_from_boundary_fixed_point','diagonal_normalization_nonzero_pair_checked']:assert v[key]
 for case in v['cases']:
  mode=case['mode'];folder=HERE/'weighted-unmarked7'/(mode if mode in ['graph','conic','cubicgraph'] else ('tails/'+mode if mode in ['Aone','Bzero'] else 'reduced/'+mode));assert case['certificate_sha256']==sha(folder/'certificate.json') and case['input_sha256']==sha(folder/'input.json')
 for path in ['linear/checkpoint.json','reduced/linear/checkpoint.json','reduced/parallel/lift-checkpoint.json','reduced/parallel-simple/lift-checkpoint.json','reduced/parallel-simple/liftstd-checkpoint.json']:assert load(HERE/'weighted-unmarked7'/path)['status']=='timeout'
 assert load(HERE/'weighted-unmarked7/parallel/lift-postprocessing-checkpoint.json')['status']=='postprocessing-stopped'
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   if 'KCUB14' in item.get('reviewed_claims',{}):assert item['reviewed_claims']['KCUB14']==claim_fingerprint(e)
 for name in ['WEIGHTED_UNMARKED_DEGREE7_CLASSIFICATION.md','WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md','WEIGHTED_TAIL_CLASSIFICATION.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB14'],validation_sha256=sha(HERE/'weighted-unmarked7/validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,unknown_timeouts_preserved=True,written_arguments_formalized=False,assumed_collisions=0,weighted_nonzero_coefficient_sector_minimum_terms=13,global_minima_closed=False,umbrella_status='open');(HERE/'weighted-unmarked7-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
