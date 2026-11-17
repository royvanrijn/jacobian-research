#!/usr/bin/env python3
"""Scoped authority/evidence and finite-experiment status audit for KCUB20."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e['id']:e for e in load(ROOT/'MATH_STATUS.json')['entries']};e=known['KCUB20'];assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS;assert e['state']=='proved' and e['proof_type']=='hybrid' and e['dependencies']==['KCUB3','KCUB5'] and e['programme_status']=='archived';assert e['independent_replay'] and not e['formal_verification'] and not e['external_review'];assert not e['closes_problems'] and not e['narrows_problems'];assert e['artifact_hash']=='sha256:'+sha(ROOT/e['checker'])
 for path in [e['canonical_source'],e['checker'],*e['software_lock']]:assert (ROOT/path).is_file(),path
 root=HERE/'carrier-reparam-full';v=load(root/'validation.json');assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_carrier_full.py') and v['input_sha256']==sha(root/'input.json');assert v['arithmetic_sha256']==sha(HERE/'verify.py') and v['utility_sha256']==sha(HERE/'weighted_verify_utils.py') and v['laurent_utility_sha256']==sha(HERE/'verify_weighted_unmarked7.py');assert (v['parameters'],v['high_equations'],v['low_slots_before_reduction'],v['exact_memberships'],v['multiplier_monomials'],v['max_multiplier_degree'])==(135,540,184,127,410,5);assert (v['dimension'],v['ordinary_degree'],v['generic_nonlinear_total'],v['minimum_nonlinear_total'],v['next_possible_support'],v['nonlinear_output_rank'],v['minimum_cubic_output_rank'],v['minimum_standard_homogeneous_companion_dimension'])==(12,3,64,36,39,12,5,18)
 for key in ['all_memberships_multiplied','complete_weighted_degree23_triangular_atoms','all_source_substitutions_triangular','entire_composition_reconstructed','complete_eight_parameter_family_sufficiency_checked','diagonal_normalization_of_entire_family_checked','minimum_family_is_diagonal_KCUB3','remaining_rank_block_three_minors_checked','remaining_rank_block_coefficient_bezout_checked','no_affine_linear_output_component','rank5_full_linear_conjugacy_checked','polynomial_substitution_domain_guarded']:assert v[key]
 assert not v['written_elimination_and_counting_formalized'];assert v['certificate_sha256']==sha(root/'classification-certificate.json') and v['family_sha256']==sha(root/'family.json');assert v['base_witness_sha256']==sha(HERE/'cubic12-integer-witness.json') and v['base_validation_sha256']==sha(HERE/'integer-validation12.json')
 inp=load(root/'input.json');assert inp['source_sha256']==sha(HERE/'carrier_reparam_full.py') and inp['input_sha256']==sha(HERE/'cubic12-integer-witness.json');family=load(root/'family.json');assert family['source_sha256']==sha(HERE/'carrier_full_family.py') and family['input_sha256']==sha(root/'input.json');cert=load(root/'classification-certificate.json');assert cert['source_sha256']==sha(HERE/'carrier_full_lift.py') and cert['input_sha256']==sha(root/'input.json') and cert['singular_source_sha256']==sha(root/'classification-lift.sing')
 folders=sorted((HERE/'carrier-order-search').glob('*/checkpoint.json'));assert len(folders)==10;dims=[]
 for path in folders:
  c=load(path);assert c['status']=='primary-bases-complete' and c['source_sha256']==sha(HERE/'carrier_order_search.py') and c['input_sha256']==sha(path.with_name('input.json'));text=path.with_name('model.log').read_text();ids=re.findall(r'CARRIER:([0-9]+):DIM:(-?[0-9]+)',text);assert len(ids)==9 and all(dim=='-1' for row,dim in ids) and 'PROPER:' not in text;dims.append(int(re.search(r'MODEL_DIM:(-?[0-9]+)',text).group(1)));rc=load(path.with_name('rank-checkpoint.json'));assert rc['status']=='primary-bases-complete' and rc['source_sha256']==sha(HERE/'carrier_rank_search.py') and rc['input_sha256']==sha(path.with_name('input.json'));rt=path.with_name('rank-model.log').read_text();branches=re.findall(r'BRANCH:[0-9]+:PIVOT:[0-9]+:DIM:(-?[0-9]+)',rt);assert len(branches)==4 and all(dim=='-1' for dim in branches) and 'PROPER:' not in rt
 assert dims.count(-1)==2
 orth=load(HERE/'orthogonal-line12/checkpoint.json');assert orth['status']=='complete' and orth['source_sha256']==sha(HERE/'orthogonal_line12.py') and orth['input_sha256']==sha(HERE/'cubic12-integer-witness.json') and orth['best']==36 and orth['generic_support']==60 and [r['nonlinear_total'] for r in orth['scores']]==[36,54,59,59,59,59]
 sparse=load(root/'sparse-checkpoint.json');assert sparse['source_sha256']==sha(HERE/'carrier_full_sparse.py') and sparse['input_sha256']==sha(root/'family.json') and sparse['tests'][0]['bound']==35 and sparse['tests'][0]['status']=='unsat' and sparse['constant_nonlinear_occurrences']==31 and len(sparse['hyperplanes'])==16
 for filename,key in [('knowledge/lessons.json','lessons'),('knowledge/work_items.json','actions')]:
  for item in load(ROOT/filename)[key]:
   for identifier in ['KCUB19','KCUB20']:
    if identifier in item.get('reviewed_claims',{}):assert item['reviewed_claims'][identifier]==claim_fingerprint(known[identifier])
 for name in ['FULL_CUBIC_CARRIER_CLASSIFICATION.md','CARRIER_REALIZATION_LIMIT.md','PROOF12.md','README.md']:
  for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
   if not target.startswith(('http:','https:','#')):assert (HERE/target.split('#')[0]).is_file(),(name,target)
 assert known['OP-GLDC']['state']=='open';receipt=dict(status='pass',claims=['KCUB20'],validation_sha256=sha(root/'validation.json'),source_hashes=True,schema=True,links=True,reviewed_fingerprints=True,complete_original_order_atoms=135,exploratory_changed_orders=10,exploratory_runs_used_as_global_proof=False,primary_smt_unsat_not_proof_premise=True,affine_linear_component_deletion_excluded_in_class=True,minimum_standard_companion_in_class=18,rank5_locus_tensor_not_new=True,positive_bounds_unchanged=True,global_minima_closed=False,umbrella_status='open');(HERE/'carrier-full-integration.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
