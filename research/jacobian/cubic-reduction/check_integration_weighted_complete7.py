#!/usr/bin/env python3
"""Scoped status/evidence checks for KCUB15 and its attributed plane import."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB15'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['proof_type']=='hybrid' and e['dependencies']==['KCUB14'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 root=HERE/'weighted-zero7';v=load(root/'validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_zero7.py') and v['arithmetic_sha256']==sha(HERE/'verify.py') and v['linear_audit_sha256']==sha(HERE/'verify_weighted_unmarked7.py');assert v['prior_nonzero_sector_validation_sha256']==sha(HERE/'weighted-unmarked7/validation.json');assert (v['cleared_parameter_identities'],v['unit_identities'],v['multiplier_monomials'],v['max_multiplier_degree'],v['collision_equations_assumed'],v['weighted_counterexample_minimum_degree'],v['weighted_counterexample_minimum_terms_through_degree7'])==(9,4,12342,15,0,7,13);assert len(v['cases'])==13
 for key in ['universal_three_variable_determinants_reconstructed','general_mu_nu_completed_square_checked','all_boundary_constraints_reconstructed','linear_parameterizations_certified_with_allowed_monomial_pivots','all_certificates_multiplied','highest_w_identities_and_bounds_checked','written_univariate_lemma_source_bounds_and_algebra_checked','constant_first_quotient_reduces_to_plane_degree_at_most4']:assert v[key]
 for case,h in v['models'].items():assert h==sha(root/case/'input.json')
 for c in v['cases']:
  case,mode=c['case'],c['mode'];path='reduced/purecubic' if mode=='purecubic' else ('reduced/parabola' if case=='nu-zero' and mode=='parabola' and c['target']=='1' else ('high/generic' if mode=='high' else mode));folder=root/case/path;assert c['input_sha256']==sha(folder/'input.json') and c['certificate_sha256']==sha(folder/'certificate.json')
 ip=root/'sources/plane-quartic-import.json';imp=load(ip);assert v['plane_degree4_theorem_import_sha256']==sha(ip) and imp['sha256']==sha(root/'sources/nakai-baba-1977.pdf') and imp['doi']=='10.18910/8846' and imp['primary_statement_visually_checked'] and not imp['reproved_here'] and not v['plane_theorem_reproved'];assert load(root/'nu-zero/linear/checkpoint.json')['status']=='timeout'
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   if 'KCUB15' in item.get('reviewed_claims',{}):assert item['reviewed_claims']['KCUB15']==claim_fingerprint(e)
 for name in ['WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md','WEIGHTED_UNMARKED_DEGREE7_CLASSIFICATION.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB15'],validation_sha256=sha(root/'validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,timeouts_preserved=True,plane_degree4_import_attributed=True,plane_theorem_reproved=False,weighted_degree_minimum=7,weighted_term_minimum_through_degree7=13,global_minima_closed=False,umbrella_status='open');(HERE/'weighted-complete7-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
