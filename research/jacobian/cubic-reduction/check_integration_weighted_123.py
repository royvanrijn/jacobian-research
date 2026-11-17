#!/usr/bin/env python3
"""Scoped ledger/evidence and retained-UNKNOWN checks for KCUB16."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB16'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['proof_type']=='hybrid' and e['dependencies']==['KCUB15'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 root=HERE/'weighted-sparse7';v=load(root/'1_-2_-3/validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_123.py') and v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py');assert (v['certificate_count'],v['multiplier_monomials'],v['max_multiplier_degree'])==(6,9017,13);assert v['complete_weights']==[1,-2,-3] and v['ordinary_degree_bound']==7 and v['all_maps_automorphisms'] and not v['other_two_weight_searches_closed']
 for key in ['all_three_universal_models_reconstructed','automatic_collision_gate_restrictions_checked','all_six_identities_multiplied','necessary_boundary_equations_reconstructed','balanced_degree23_control_checked','exceptional_parameter_values_exactly_resolved','univariate_integral_divisibility_algebra_checked']:assert v[key]
 for c in v['certificates']:
  folder=root/'1_-2_-3/slices'/c['case']/c['mode'];assert c['input_sha256']==sha(folder/'input.json') and c['certificate_sha256']==sha(folder/'certificate.json')
 for g in v['gates']:
  tag='_'.join(map(str,g['weights']));folder=root/tag;assert g['input_sha256']==sha(folder/'input.json') and g['negative_plane_triangular'];assert load(folder/'nra-checkpoint.json')['status']=='unknown' and load(folder/'groebner-checkpoint.json')['status']=='timeout'
 ip=HERE/'weighted-zero7/sources/plane-quartic-import.json';imp=load(ip);assert v['plane_theorem_import_sha256']==sha(ip) and imp['sha256']==sha(ip.with_name('nakai-baba-1977.pdf')) and not v['plane_theorem_reproved']
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   if 'KCUB16' in item.get('reviewed_claims',{}):assert item['reviewed_claims']['KCUB16']==claim_fingerprint(e)
 for name in ['WEIGHTED_123_DEGREE7_EXCLUSION.md','WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB16'],validation_sha256=sha(root/'1_-2_-3/validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,all_three_original_unknown_timeouts_preserved=True,complete_weight_exclusion=[1,-2,-3],original_two_other_queries_remain_unknown=True,mixed_class_subsequently_closed=("KCUB17" in known and known["KCUB17"]["state"]=="proved"),remaining_initial_weights_open=([] if "KCUB18" in known else ([[1,-1,-3]] if "KCUB17" in known else [[1,-1,-3],[2,-1,-3]])),third_class_subsequently_closed=("KCUB18" in known and known["KCUB18"]["state"]=="proved"),global_minima_closed=False,umbrella_status='open');(HERE/'weighted-123-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
