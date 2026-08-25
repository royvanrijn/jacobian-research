#!/usr/bin/env python3
"""Build the durable q4/o1584 route-optimization handoff."""
import hashlib, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]; GEN=ROOT/'artifacts/generated-results'
FILES={
 'o1584_qq':ROOT/'artifacts/local/elkies-k3/q4o208-physical-q4o1584-rr-qq.json',
 'o1584_cert':GEN/'elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json',
 'o1584_zeros':GEN/'elkies-k3-h3-q4o208-physical-q4o1584-effective-zero-markings.json',
 'o164_cert':GEN/'elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json',
 'o164_zeros':GEN/'elkies-k3-h3-q4o208-q4o1584-q4o164-effective-zero-markings.json',
 'corrected_q8':GEN/'elkies-k3-h3-q4o208-corrected-a3-2a2-to-5a1-physical-q8-c10-certificate.json',
 'physical_q12':GEN/'elkies-k3-h3-physical-5a1-to-4a1-q12-certificate.json',
 'a2_cert':GEN/'elkies-k3-h3-physical-4a1-q4o21633-a2-certificate.json',
 'a2_marking':GEN/'elkies-k3-h3-physical-a2-effective-zero-marking.json',
 'a2_q4_shell':GEN/'elkies-k3-h3-physical-a2-q4d2-rootless-neighbors.json',
}
O1584_QQ_SCRIPT=ROOT/'elkies-k3/scripts/certify_h92_q4o208_physical_q4o1584_rr_qq.sage'
HASH_FILES=tuple(FILES.values())+(O1584_QQ_SCRIPT,)
data={k:json.loads(p.read_text()) for k,p in FILES.items()}
assert data['o1584_qq']['status']=='PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN'
assert data['o1584_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['o164_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['corrected_q8']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['physical_q12']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['a2_cert']['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
assert data['a2_marking']['status']=='PASS_EXACT_PHYSICAL_A2_EFFECTIVE_ZERO_MARKING'
assert not [n for n in data['a2_q4_shell']['neighbors'] if n['child_mw_rank']==17]
payload={
 'schema':'elkies-k3.h3-q4o1584-route-optimization-handoff.v1',
 'status':'PASS_EXACT_Q4O1584_ROUTE_OPTIMIZATION_HANDOFF_NOT_PROMOTED_TO_PINNED_R17',
 'equation_frontier':'q4o1584 physical D4+A3+3A1/MW7, exact QQ Jacobian',
 'best_equation_realizable_branch':[
  {'edge':'q4/orbit1584','source_zero':'old_A11_component_5','literal_divisor':'O(C5)+first_I6_affine+r3+2r4+r5','P_dot_O':0,'expected_RR_ambient':4,'child':'3A1+A3+D4/MW7','lattice_certificate':str(FILES['o1584_cert'].relative_to(ROOT)),'equation_status':'PASS exact QQ: literal 4-to-2 jet RR construction and global Jacobian','equation_artifact':str(FILES['o1584_qq'].relative_to(ROOT)),'equation_script':str(O1584_QQ_SCRIPT.relative_to(ROOT))},
  {'edge':'q4/orbit164','source_zero':'second_I6_affine_component','literal_divisor':'O(second_I6_affine)+old_A11_component_0+(F-r7)+(F-r9)','P_dot_O':0,'expected_RR_ambient':4,'child':'2A1+2A3/MW9','lattice_certificate':str(FILES['o164_cert'].relative_to(ROOT)),'equation_status':'recommended immediate next lifting target; continuation to pinned R17 is not yet certified'},
 ],
 'physical_suffix_improvement':[
  {'edge':'corrected A3+2A2 to 5A1','source_zero':'old_A11_component_10','q':8,'old_fibre_degree':2,'P_dot_O':2,'child':'5A1/MW12','certificate':str(FILES['corrected_q8'].relative_to(ROOT))},
  {'edge':'5A1 to 4A1','source_zero':'old_A11_component_5','q':12,'old_fibre_degree':2,'P_dot_O':4,'child':'4A1/MW13','certificate':str(FILES['physical_q12'].relative_to(ROOT))},
  {'edge':'4A1 to A2','source_zero':'old_A11_component_0','q':4,'old_fibre_degree':2,'P_dot_O':1,'expected_RR_ambient':5,'child':'A2/MW15','certificate':str(FILES['a2_cert'].relative_to(ROOT)),'physical_child_marking':str(FILES['a2_marking'].relative_to(ROOT))},
 ],
 'negative_search_results':[{'hub':'physical A2/MW15','shell':'q4, old-fibre degree 2','result':'complete 97030-orbit shell contains no rootless/MW17 child','artifact':str(FILES['a2_q4_shell'].relative_to(ROOT))}],
 'promotion_gate':'Do not switch the pinned-R17 lifting suffix yet: q4/o164 is fully lattice-certified and literal but not equation-lifted, and neither the o164 continuation nor the alternate physical A2/MW15 continuation has yet been fully certified to pinned R17.',
 'inputs':{'paths':[str(p.relative_to(ROOT)) for p in HASH_FILES],'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in HASH_FILES}},
}
out=GEN/'elkies-k3-h3-q4o1584-route-optimization-handoff.json';out.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(out)
