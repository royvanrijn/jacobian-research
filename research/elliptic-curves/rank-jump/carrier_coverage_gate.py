#!/usr/bin/env python3
"""Exact dimension gate: can the marked native dictionary cover retained gain?"""
import argparse
from collections import Counter
from pathlib import Path
import retrospective as r
from fibre_discrimination import hash_file

HERE=Path(__file__).resolve().parent
CENSUS=r.OUT/'rank_jump_fibre_discrimination_v1.json'
VERIFIED=r.OUT/'rank_jump_fibre_discrimination_verification_v1.json'
RELATIONS=r.OUT/'rank_jump_low_degree_triple_discrimination_v1.json'
OUTPUT=r.OUT/'rank_jump_marked_carrier_coverage_gate_v1.json'


def compute():
    census=r.read(CENSUS);verified=r.read(VERIFIED);relations=r.read(RELATIONS)
    for doc in (census,verified,relations):
        assert doc['status']=='PASS'
        for path,sha in doc['bindings'].items():assert hash_file(r.ROOT/path)==sha
    rr={x['block_key']:x for x in relations['rows']};best={};untested=[]
    for row in census['rows']:
        if row['block_key'] is None:untested.append(row['observation_id']);continue
        key=row['block_key']
        if key not in best or row['retained_quotient_rank']>best[key]['retained_quotient_rank']:best[key]=row
    rows=[]
    for key,old in sorted(best.items()):
        n=old['compatible_cover_count'];j=old['retained_quotient_rank'];g=old['generic_subgroup_rank']
        assert old['retained_rank_lower_bound']==g+j
        c=0
        if key in rr:
            x=rr[key];assert x['compatible_cover_count']==n and x['best_retained_quotient_rank']==j
            c=x['combined_relation_rank'];assert x['native_quotient_upper_bound']==n-c
        else:assert old['dictionary']=='11952' and n<=1
        assert 0<=c<=n
        cap=n-c;deficit=max(0,j-cap)
        rows.append({'block_key':key,'dictionary':old['dictionary'],'dictionary_size':old['dictionary_size'],
            'dictionary_coverage':old['dictionary_coverage'],'parameter':old['native_parameter'],
            'generic_subgroup_rank':g,'retained_quotient_rank':j,'compatible_cover_count':n,
            'certified_relation_rank':c,'marked_native_quotient_upper_bound':cap,
            'retained_directions_outside_marked_native_span_lower_bound':deficit,
            'can_cover_entire_retained_quotient':'NO' if deficit else 'NOT_EXCLUDED',
            'full_curve_rank':'UNKNOWN','retained_observation_id':old['observation_id']})
    groups=[]
    for name in ('published-R17','11952'):
        sub=[x for x in rows if x['dictionary']==name]
        groups.append({'dictionary':name,'addresses':len(sub),'positive_retained_gain_addresses':sum(x['retained_quotient_rank']>0 for x in sub),
            'provably_incomplete_marked_span_addresses':sum(x['can_cover_entire_retained_quotient']=='NO' for x in sub),
            'deficit_count_distribution':dict(sorted(Counter(str(x['retained_directions_outside_marked_native_span_lower_bound']) for x in sub).items()))})
    return {'schema':'rank-jump.marked-carrier-coverage-gate.v1','status':'PASS','rows':rows,'groups':groups,
        'untested_observation_ids':untested,
        'linear_algebra_lemma':'In V=(E_t(Q) tensor Q)/M_t let L be the retained j-dimensional span and N the span of n marked native points with c independent certified relation rows. Then dim(N)<=n-c and dim(L/(L intersection N))>=max(0,j-n+c). No containment of N in L is assumed.',
        'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (CENSUS,VERIFIED,RELATIONS,Path(__file__),HERE/'retrospective.py',HERE/'fibre_discrimination.py')},
        'boundary':'Retrospective coverage exclusion for the displayed native sections only. It does not bound all sections over the full native base change, rule out other dictionaries, or prove exact fibre ranks. Deficit zero does not prove coverage. Untested families remain untested; no point search or selector change.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    for g in d['groups']:print(g)
