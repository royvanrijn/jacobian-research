#!/usr/bin/env python3
"""Scoped authority/evidence checks for KCUB18, with untouched solver failures."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB18'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['proof_type']=='hybrid' and e['dependencies']==['KCUB15'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 root=HERE/'weighted-113';v=load(root/'validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_113.py') and v['input_sha256']==sha(root/'input.json');assert v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py') and v['linear_audit_sha256']==sha(HERE/'verify_weighted_unmarked7.py');assert (v['nonlinear_slots'],v['certificate_count'],v['multiplier_monomials'],v['max_multiplier_degree'],v['plane_family_degree_bound'])==(20,7,5728,13,5);assert v['complete_weights']==[1,-1,-3] and v['ordinary_degree_bound']==7 and v['all_maps_automorphisms']
 for key in ['universal_three_variable_determinant_reconstructed','quotient_jacobian_checked','all_six_boundary_parameterizations_checked','all_seven_identities_multiplied','hard_quadratic_core_source_conditions_checked','core_laurent_coefficient_denominators_preserved','highest_transverse_coefficient_identities_checked','univariate_integral_identities_checked','plane_family_unimodular_linear_part_checked']:assert v[key]
 assert not v['written_case_split_formalized'] and not v['plane_theorem_reproved']
 for c in v['certificates']:
  folder=root/c['mode'];assert c['input_sha256']==sha(folder/'input.json') and c['certificate_sha256']==sha(folder/'certificate.json')
 assert len(v['branches'])==6
 for c in v['branches']:
  folder=root/'boundary'/(c['case']+'-'+c['degree']);assert c['input_sha256']==sha(folder/'input.json') and c['reduced_input_sha256']==sha(folder/'reduced/input.json')
 for tag in ['1_-1_-3','2_-1_-3','1_-2_-3']:assert load(HERE/'weighted-sparse7'/tag/'nra-checkpoint.json')['status']=='unknown' and load(HERE/'weighted-sparse7'/tag/'groebner-checkpoint.json')['status']=='timeout'
 for case in ['qzero','both','mu-zero','both-linear']:assert load(root/case/'checkpoint.json')['status']=='timeout'
 for branch in ['reduced','field','clean']:assert load(root/'boundary/both-quad'/branch/'checkpoint.json')['status']=='timeout'
 ip=root/'plane-degree5-import.json';imp=load(ip);assert v['plane_import_sha256']==sha(ip) and imp['primary_pdf_sha256']==sha(HERE/'weighted-zero7/sources/nakai-baba-1977.pdf') and imp['source_import_sha256']==sha(HERE/'weighted-zero7/sources/plane-quartic-import.json') and not imp['reproved_here']
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   for identifier in ['KCUB16','KCUB17','KCUB18']:
    if identifier in item.get('reviewed_claims',{}):assert item['reviewed_claims'][identifier]==claim_fingerprint(known[identifier])
 for name in ['WEIGHTED_113_DEGREE7_EXCLUSION.md','WEIGHTED_213_DEGREE7_EXCLUSION.md','WEIGHTED_123_DEGREE7_EXCLUSION.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB18'],validation_sha256=sha(root/'validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,all_original_unknown_timeouts_preserved=True,changed_hard_quadratic_timeouts_preserved=True,complete_weight_exclusion=[1,-1,-3],remaining_initial_weights_open=[],written_case_split_formalized=False,plane_theorem_reproved=False,global_minima_closed=False,umbrella_status='open');(HERE/'weighted-113-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
