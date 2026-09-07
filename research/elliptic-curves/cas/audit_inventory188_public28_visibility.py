#!/usr/bin/env python3
"""Exact retrospective visibility of a certified public direction missed by the old27 search."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from search_observability import point_visibility
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PROOF=ART/'inventory188_public28_reproduction_v1.json'
REPLAY=ART/'inventory188_public28_sage_replay_v1.json'
SOURCE=ART/'full11952_late64_r17_results_v1.json'
OUT=ART/'inventory188_public28_visibility_v1.json'


def expected():
    proof=cert.read(PROOF);replay=cert.read(REPLAY)
    if proof['status']!='PASS' or replay['status']!='PASS' or proof['rank_lower_bound']!=28:
        raise ArithmeticError('independent public28 proof required')
    if proof['independent_column_indices']!=list(range(27))+[53]:
        raise ArithmeticError('fixed certified exceptional representative required')
    for d in (proof,replay):
        if any(cert.hashed(ROOT/n)!=h for n,h in d['sources'].items()):raise ArithmeticError('public proof changed')
    source=cert.read(SOURCE);row=next(r for r in source['curves'] if r['parameter']==proof['parameter'])
    rawpath=ROOT/row['discovery_witness']['path']
    if cert.hashed(rawpath)!=row['discovery_witness']['sha256']:raise ArithmeticError('original transcript changed')
    raw=cert.read(rawpath)
    if raw['curve']!=proof['curve'] or row['rank_lower_bound']!=27 or row['completed_boxes']!=49 or len(raw['charts'])!=49:
        raise ArithmeticError('same completed49-chart old27 trial required')
    x,y=map(cert.F,proof['transported_public_points'][26]);observations=[]
    for i,chart in enumerate(raw['charts']):
        r=chart['search']
        if r['status']!='bounded_search_complete' or r['height_bound']!=125000:
            raise ArithmeticError('complete fixed125000 box required')
        for sign in (-1,1):
            observations.append({'chart':i,'sign':sign,'visibility':point_visibility({**r,'completed_denominator':125000},(x,sign*y))})
    finite=[r for r in observations if r['visibility'].get('minimum_affine_height') is not None]
    best=min(finite,key=lambda r:(r['visibility']['minimum_affine_height'],r['chart'],r['sign'])) if finite else None
    paths=[Path(__file__).resolve(),ROOT/'elliptic-curves/cas/search_observability.py',PROOF,REPLAY,SOURCE,rawpath]
    return {'schema':'elliptic-curves.inventory188-public28-visibility.v1','status':'PASS_EXACT_OBSERVATIONS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
            'local_id':proof['local_id'],'public_id':619,'public_point_index':26,
            'point':list(map(str,(x,y))),'observations':observations,'best_observation':best,
            'status_counts':dict(Counter(r['visibility']['status'] for r in observations)),
            'discrepancies':[r for r in observations if r['visibility']['status']=='VISIBLE_NOT_RECORDED'],
            'claim_boundary':'Retrospective exact coordinates for both signs of public point26, independently proved outside the old27 subgroup, in all49 original completed125000 charts. It is a public witness diagnostic, not prospective selection or a new point search. Being outside every box excludes this representative only, not all translates, all points of its quotient direction or the existence of a cheap recovery.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('public visibility replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve public visibility audit')
        checkpoint(OUT,d)
    print('PUBLIC28 REPRESENTATIVE VISIBILITY',d['status_counts'],'BEST',d['best_observation'])
