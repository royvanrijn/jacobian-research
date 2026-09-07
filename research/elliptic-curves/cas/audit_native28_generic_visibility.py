#!/usr/bin/env python3
"""Positive-control check of what direct generic-point visibility can diagnose."""
import argparse
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from search_observability import point_visibility
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000';OUT=ART/'native28_generic_visibility_v1.json'
def expected():
    v=cert.read(D/'verification.json');raw=cert.read(D/'result.json')
    if v['status']!='PASS' or v['completed_boxes']!=49 or v['rank_lower_bound']!=28 or len(raw['generic_points'])!=17:raise ArithmeticError('completed known28 recovery control required')
    observations=[];counts=Counter();seen=set();best=[None]*17;discrepancies=[]
    for i,c in enumerate(raw['charts']):
        r=c['search']
        if r['status']!='bounded_search_complete' or r['height_bound']!=125000:raise ArithmeticError('full control box required')
        for k,P in enumerate(raw['generic_points']):
            x,y=map(cert.F,P)
            for sign in (-1,1):
                a=point_visibility({**r,'completed_denominator':125000},(x,sign*y));row={'section':k,'sign':sign,'chart':i,'visibility':a};observations.append(digest(row));counts[a['status']]+=1
                if a['status']=='VISIBLE_NOT_RECORDED':discrepancies.append(row)
                if a.get('in_completed_box') or a['status']=='KNOWN_POINTED_ENDPOINT':seen.add(k)
                H=a.get('minimum_affine_height');key=(0 if a.get('at_parameter_infinity') or a['status']=='KNOWN_POINTED_ENDPOINT' else H if H is not None else 10**1000,i,sign)
                if best[k] is None or key<best[k][0]:best[k]=(key,row)
    if len(observations)!=1666:raise ArithmeticError('fixed49x17x2 control product differs')
    paths=[Path(__file__).resolve(),ROOT/'elliptic-curves/cas/search_observability.py',D/'verification.json',D/'result.json']
    return {'schema':'elliptic-curves.native28-generic-visibility.v1','status':'PASS_EXACT_OBSERVATIONS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'point_chart_observations':len(observations),'ordered_observation_digest':digest(observations),'status_counts':dict(counts),'generic_sections_with_a_visible_sign':len(seen),'section_minima':[b[1] for b in best],'discrepancies':discrepancies,'independently_certified_recovered_lower_bound':28,'claim_boundary':'Retrospective direct visibility of the original17 generic section representatives, both signs, on the same49 maps that independently recovered28 from only17 generic seeds. A generic point representative outside these boxes does not show that exceptional directions are invisible. This is an already known curve and a completed control, not a new rank result or a new point search.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('positive-control generic visibility differs')
    else:
        if OUT.exists():raise FileExistsError('preserve control visibility audit')
        checkpoint(OUT,d)
    print('KNOWN28 CONTROL GENERIC VISIBILITY',d['generic_sections_with_a_visible_sign'],'of17;',len(d['discrepancies']),'discrepancies',flush=True)
