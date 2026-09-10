#!/usr/bin/env python3
"""Portable hash manifest for the class1 bridge and carrier proof bundle."""
import argparse,hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUTPUT=ART/'curve302_class1_bridge_manifest_v1.json'
SOURCES=[
 'curve302_class1_point_transport.sage','curve302_class1_fibre_carriers.sage',
 'curve302_class1_twist_sections.sage','curve302_class1_branch_veronese.sage',
 'curve302_class1_prescribed_core_glue.sage','verify_curve302_class1_bridge.sage',
 'verify_curve302_class1_bridge_bundle.py','package_curve302_class1_bridge.py',
 'icarm_curve302.py','half_lattice_pointed_sieve.py','alternate_quartic_covers.py',
 'research_runtime/finite_reduction.py','research_runtime/arithmetic.py',
 'research_runtime/binary.py','research_runtime/store.py','research_runtime/memory_store.py',
 'mod2_reduction_independence.py','certify_compact_r17_candidates.py','elliptic_candidate_record.py',
]
ARTIFACTS=[
 'curve302_class1_point_transport_protocol_v1.json','curve302_class1_point_transport_v1.json',
 'curve302_class1_fibre_carriers_v1.json','curve302_class1_bridge_verification_v1.json',
 'curve302_class1_twist_sections_v1.json','curve302_class1_branch_veronese_v1.json',
 'curve302_class1_prescribed_core_glue_obstruction_v1.json',
 'x1092_class1_realization_manifest_v1.json',
]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compute():
    paths=[ROOT/'elliptic-curves/cas'/n for n in SOURCES]+[ART/n for n in ARTIFACTS]
    paths += [ROOT/'elliptic-curves/notes/CURVE302_CLASS1_BRIDGE_AND_SHARED_CARRIERS_2026-09-10.md']
    paths += [ROOT/'elliptic-curves/ecsearch/q12o5867_specialization.py']
    bindings={}
    while paths:
        path=paths.pop();key=str(path.relative_to(ROOT))
        if key in bindings:continue
        bindings[key]=sha(path)
        if path.suffix!='.json':continue
        data=read(path)
        if not isinstance(data,dict):continue
        for field in ['bindings','inputs','input','proof_packets','software']:
            values=data.get(field,{})
            if not isinstance(values,dict):continue
            for name,digest in values.items():
                if not isinstance(digest,str)or not re.fullmatch('[0-9a-f]{64}',digest):continue
                p=ROOT/name;assert sha(p)==digest,(key,name);paths.append(p)
    protocol=read(ART/ARTIFACTS[0]);visibility=read(ART/'curve302_residual_visibility_geometry_v1.json')
    assert protocol['rows']==[{k:r[k]for k in ['id','public_word','strict_mod_2']}for r in visibility['directions']]
    assert protocol['integral_basis_determinant']==-1
    transport=read(ART/ARTIFACTS[1]);replay=read(ART/ARTIFACTS[3])
    assert transport['complete']and replay['status']=='PASS'
    assert replay['distinct_B_parameters']==replay['outside_generic_rational_span_count']==14
    assert len(replay['point_transport_and_rank_cases'])==14
    assert all(r['finite_group_proof']['rank_lower_bound']==18 for r in replay['point_transport_and_rank_cases'])
    core=read(ART/'curve302_class1_prescribed_core_glue_obstruction_v1.json')
    assert core['core_automorphisms']=='exactly +I and -I'
    assert core['refinement_class_counts']==[1,248,516]
    assert core['embeddings'][0]['unique_order_two_glue_bits']!=core['embeddings'][1]['unique_order_two_glue_bits']
    twist=read(ART/'curve302_class1_twist_sections_v1.json')
    assert len(twist['sections'])==14 and all(s['section_height']==8 and s['native_section_primitive']for s in twist['sections'])
    branch=read(ART/'curve302_class1_branch_veronese_v1.json')
    assert branch['matrix_rank']==5 and branch['projective_map_degree']==4
    return {'schema':'curve302.class1-bridge.manifest.v1','bindings':dict(sorted(bindings.items())),
            'software':{'Sage':'10.9','PARI':'2.17.3'},
            'results':{'prescribed_core_determinant':4100,'prescribed_core_integral_extension':'IMPOSSIBLE for the retained embeddings25/33, even after any core isometry',
                       'actual_realization_common_core_determinant':13104,'transported_directions':14,
                       'distinct_B_parameters':14,'B_fibres_with_certified_rank_at_least18':14,
                       'generic_B_rational_span_memberships':0,'native_carrier_genus':1,
                       'native_pair_count':91,'native_pair_joint_genus':5,
                       'smooth_B_fibre_branch_map':'degree4 Veronese embedding',
                       'primitive_native_twist_section_height':8,'geometric_twist_height_lower_bound':4,
                       'geometric_twist_rank_interval':[1,22],
                       'fibre_preserving_automorphism_action_mod_generic':'plus or minus the seed',
                       'second_independent_twist_section':'UNKNOWN','propagation_theorem':'UNKNOWN'},
            'replay_from_repository_root':'python3 research/elliptic-curves/cas/verify_curve302_class1_bridge_bundle.py',
            'foundation_replay':'timeout 180 sage -python research/elkies-k3/scripts/verify_x1092_class1_realization.sage',
            'limits':{'workers':1,'transport_case_seconds':60,'independent_case_seconds':60,
                      'complete_point_replay_seconds':180,'core_replay_seconds':120,'twist_replay_seconds':120,
                      'branch_replay_seconds':60,'point_searches':0,'new_fibrations':0,'translations_tried':0},
            'boundary':'Exact results for the specified overlap and existing class1 realization. Other core embeddings, other carriers, full twist ranks, class/Selmer upper bounds and seed propagation remain open.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();result=compute()
    if args.check:assert read(OUTPUT)==result
    else:
        with OUTPUT.open('x')as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print('PASS class1 bridge manifest;',len(result['bindings']),'file bindings')
