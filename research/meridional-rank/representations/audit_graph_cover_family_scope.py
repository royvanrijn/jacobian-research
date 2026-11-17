#!/usr/bin/env python3
"""Bind sixteen retained case certificates to a written classical family audit.

This checks finite provenance and uppers. Closed rigidity and BJW18 are
imported in the companion proof, not proved by this script.
"""
import hashlib
import json
from pathlib import Path
import time
from verify_double_transpositions import geometric_crossings,check_trace

HERE=Path(__file__).resolve().parent


def read(name):return json.loads((HERE/name).read_text())


def seals(record):
    for p,digest in record['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest,p


def main():
    started=time.monotonic();plan=read('graph_cover_family_scope_inputs.json')
    four=read('four_piece_cover_verification.json');three=read('three_piece_cover_verification.json');five=read('final_five_cover_verification.json')
    mutation=read('mutation_cover_verification.json');sources=read('graph_cover_family_source_audit.json')
    for record in (four,three,five,mutation,sources):seals(record)
    assert four['status']=='PASS independent four-piece peripheral and upper4 inputs'
    assert three['status']=='PASS independent three-root peripheral and upper4 inputs; five transfers unproved'
    assert five['status']=='PASS independent final-five source, peripheral, filling and upper4 inputs'
    # Historical source-catalogue counts are not mathematical replay inputs.
    panel={r['name']:r['pd'] for r in read('parabolic_panel.json')['cases']}
    joint={r['name']:r for r in read('joint_a5_results.json')['results']}
    four_roots={r['name']:r for r in four['roots']};cases={}
    for names in mutation['blocks'].values():
        root=names[0];data=four_roots[root];assert data['proposed_JSJ_pieces']==4
        assert data['fiber_intersections']==[1,1,1]
        for name in names:
            if name!=root:assert any(p['source']==root and p['target']==name for p in mutation['pairs'])
            cases[name]={'name':name,'model_source':root,'model_route':'four_piece','canonical_JSJ_pieces':4,
                'H1_order':data['source_H1_order'],'group_rank':3,'individual_upper_source':'four_piece_cover_verification.json'}
    assert set(cases)=={r['name'] for r in four['upper4']}
    for data in three['roots']:
        name=data['name'];assert name not in cases and data['JSJ_pieces']==3
        cases[name]={'name':name,'model_source':name,'model_route':'three_piece','canonical_JSJ_pieces':3,
            'H1_order':data['H1_order'],'group_rank':3,'individual_upper_source':'three_piece_cover_verification.json'}
    for data in five['targets']:
        name=data['name'];assert name not in cases and data['proposed_JSJ_pieces']==3
        assert data['nonabelian_filling_certificate']['image_order']==6
        cases[name]={'name':name,'model_source':name,'model_route':'final_five','canonical_JSJ_pieces':3,
            'H1_order':data['source_H1_order'],'group_rank':3,'individual_upper_source':'final_five_cover_verification.json'}
    assert set(cases)==set(panel)==set(plan['targets']) and len(cases)==16
    upper_by_receipt={key:{r['name']:r for r in record['upper4']} for key,record in
        [('four_piece_cover_verification.json',four),('three_piece_cover_verification.json',three),('final_five_cover_verification.json',five)]}
    for name in plan['targets']:
        case=cases[name];row=joint[name];upper=upper_by_receipt[case['individual_upper_source']][name]
        assert row['seeds']==upper['seeds'] and len(row['seeds'])==4
        check_trace(geometric_crossings(panel[name]),row['seeds'],row['seed_trace'])
        case.update({'upper4_seeds':row['seeds'],'upper4_steps':len(row['seed_trace']),
            'actual_cover_graph_manifold':'Established by companion written closed-rigidity argument',
            'existing_MRC_family':'BJW18 published Corollary1.3',
            'named_global_bridge_rank_priority':'UNESTABLISHED'})
    # Review and programme state belong to the compact topic metadata.
    # This replay binds the retained finite cases without asserting external review.
    files=['audit_graph_cover_family_scope.py', 'graph_cover_family_scope_inputs.json', 'graph_cover_family_source_audit.json', 'ORIGINAL_SIXTEEN_KNOWN_GRAPH_COVER_FAMILY.tex', 'four_piece_cover_verification.json', 'three_piece_cover_verification.json', 'final_five_cover_verification.json', 'mutation_cover_verification.json', 'parabolic_panel.json', 'joint_a5_results.json', 'verify_double_transpositions.py', 'parabolic_lifts.py']
    receipt={'status':'PASS sixteen exact case/source bindings for known graph-cover family',
        'cases':[cases[n] for n in plan['targets']],'existing_family_count':16,'new_MRC_family_cases':0,
        'retained_named_equalities':20,'other_sharp_seed_family_audit_names':['13a_3145','13a_3147','13a_3148','13a_3149'],
        'topology_proved_by_this_script':False,'full_MRC_proved':False,
        'assurance':'Finite source names, all16individual upper traces and retained exact lower/model receipt seals checked. Actual-cover homeomorphism uses AFW15 closed prime/lens classification in the companion written proof, and existing MRC-family scope uses BJW18 Cor1.3. Both are explicit classical imports, not formal or external verification. Exact named bridge/rank global priority and the other four sharp seeds family scopes remain open.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    assert time.monotonic()-started<plan['limits']['maximum_seconds_finite_binding']
    (HERE/'graph_cover_family_scope_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'existing_family_count':16,'new_MRC_family_cases':0,'individual_upper4':16}))


if __name__=='__main__':main()
