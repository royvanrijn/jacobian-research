#!/usr/bin/env python3
"""Scoped authority, exact-evidence and exploratory-status seals for KCUB19."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB19'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['proof_type']=='hybrid' and e['dependencies']==['KCUB3'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 root=HERE/'carrier-reparam-wide';v=load(root/'validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_carrier_reparam12.py') and v['input_sha256']==sha(root/'input.json');assert v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py') and v['laurent_utility_sha256']==sha(HERE/'verify_weighted_unmarked7.py');assert (v['parameters'],v['high_equations'],v['exact_memberships'],v['multiplier_monomials'],v['max_multiplier_degree'])==(47,93,43,123,3);assert v['ordinary_degree']==3 and v['dimension']==12 and v['minimum_nonlinear_total']==36 and v['generic_nonlinear_total']==45;assert v['coefficient_factor_counts']=={'unit':33,'s':4,'t':5,'3p+2s':1,'pq+3t':2}
 for key in ['all_memberships_multiplied','all_source_substitutions_triangular','entire_carrier_composition_reconstructed','four_parameter_family_sufficient','minimum_family_diagonal_conjugacy_checked','all_degree_one_axis_obstruction']:assert v[key]
 assert len(v['one_axis_arbitrary_nonlinear_source_obstructions'])==12 and not v['written_elimination_formalized'];assert v['certificate_sha256']==sha(root/'classification-certificate.json') and v['family_sha256']==sha(root/'family.json');assert v['base_witness_sha256']==sha(HERE/'cubic12-integer-witness.json') and v['base_validation_sha256']==sha(HERE/'integer-validation12.json') and known['KCUB3']['state']=='proved'
 family=load(root/'family.json');assert family['source_sha256']==sha(HERE/'carrier_reparam_family.py') and family['input_sha256']==sha(root/'input.json');c=load(root/'classification-certificate.json');assert c['source_sha256']==sha(HERE/'carrier_reparam_lift.py') and c['input_sha256']==sha(root/'input.json') and c['singular_source_sha256']==sha(root/'classification-lift.sing')
 for tag,filename in [('quadratic-conjugacy12','quadratic_conjugacy12.py'),('carrier-family-graphs12','carrier_family_graph12.py'),('carrier-family-linear12','carrier_family_linear12.py')]:
  d=load(HERE/tag/'checkpoint.json');assert d['status']=='complete' and d['source_sha256']==sha(HERE/filename)
 lines=load(HERE/'quadratic-conjugacy12/checkpoint.json');assert lines['lines']==792 and len(lines['cases'])==792 and lines['best']==36
 graphs=load(HERE/'carrier-family-graphs12/checkpoint.json');assert len(graphs['cases'])==9 and all(c['status']=='unit-ideal-no-raw-cubic' for c in graphs['cases'])
 linear=load(HERE/'carrier-family-linear12/checkpoint.json');assert len(linear['cases'])==6 and [c['initial'] for c in linear['cases']]==[39,39,42,40,41,45] and [c['final'] for c in linear['cases']]==[39,36,39,40,36,40] and linear['global_best']==36
 audit=load(HERE/'substitution-domain-audit.json');assert audit['status']=='pass' and audit['source_sha256']==sha(HERE/'audit_substitution_domain.py') and len(audit['results'])==7 and not audit['domain_failures'] and not audit['checker_files_actually_written']
 for item in audit['results']:assert item['status']=='pass' and item['source_sha256']==sha(HERE/(item['module']+'.py'))
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   if 'KCUB19' in item.get('reviewed_claims',{}):assert item['reviewed_claims']['KCUB19']==claim_fingerprint(e)
 for name in ['CARRIER_REALIZATION_LIMIT.md','PROOF12.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB19'],validation_sha256=sha(root/'validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,scope='Sharp36 only in the displayed47-coefficient carrier realization, plus all-degree one-axis source shear obstruction in the fixed frame.',exploratory_runs_preserved=True,seven_module_substitution_domain_audit_pass=True,exploratory_runs_used_as_global_proof=False,positive_bounds_unchanged=True,global_minima_closed=False,umbrella_status='open');(HERE/'carrier-reparam-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
