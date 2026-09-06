#!/usr/bin/env python3
"""Retrospective exact visibility of older blind discoveries in completed metric boxes."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from search_observability import point_visibility
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
OLD=ART/'half_lattice_search_ablation_rank29_holdout_blind_v1.json'
PARENT=LOCAL/'native11952-metric49-control-v1/candidate-00/result.json'
TAIL=LOCAL/'native11952-metric49-tails-v1/candidate-00/result.json'
COVERAGE=ART/'native11952_tail_coverage_v1.json'
OUTPUT=ART/'native11952_metric_visibility_v1.json'


def expected():
    old=next(r for r in cert.read(OLD)['results'] if r['label']=='curve12-2024-rank29')
    parent,tail,coverage=cert.read(PARENT),cert.read(TAIL),cert.read(COVERAGE)
    if coverage['status']!='PASS' or coverage['complete_combined_boxes']!=49 or coverage['input_sha256']!=cert.hashed(TAIL) or coverage['parent_input_sha256']!=cert.hashed(PARENT):raise ArithmeticError('full metric coverage gate failed')
    arm=next(r for r in old['arms'] if r['id']=='generic-deepest43');roster=old['cover_records'];by_mask={r['mask']:r for r in roster}
    if set(arm['masks'])!={c['mask'] for c in parent['centres'][:43]}:raise ArithmeticError('original43 masks differ')
    comparisons=[]
    for c in parent['centres'][:43]:
        a=c['representative'];b=by_mask[c['mask']]['specialized_representative']
        if a!=b and a!=[-v for v in b]:raise ArithmeticError('representatives differ beyond sign')
        comparisons.append({'mask':c['mask'],'representative_sign':1 if a==b else -1})
    u=cert.F(6);original=tuple(map(cert.F,old['short_model']));model=tuple(map(cert.F,parent['curve']))
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6):raise ArithmeticError('old-control scale changed')
    if [[str(cert.F(x)/u**2),str(cert.F(y)/u**3)] for x,y in ((p['x'],p['y']) for p in old['generic_points'])]!=parent['generic_points']:raise ArithmeticError('old generic prefix changed')
    rows=[];counts=Counter()
    for index in arm['candidate_point_indices']:
        p=old['candidate_points'][index]['point'];x,y=cert.F(p['x'])/u**2,cert.F(p['y'])/u**3
        signs=[]
        for sign in (1,-1):
            observations=[]
            for j,c in enumerate(parent['charts']):
                v=point_visibility(c['search'],(x,sign*y))
                if v['status']=='UNSEARCHED_INTERVAL':v=point_visibility(tail['charts'][j]['search'],(x,sign*y))
                if v['status'] in ('UNSEARCHED_INTERVAL','VISIBLE_NOT_RECORDED'):raise ArithmeticError('metric completeness discrepancy')
                observations.append({'chart':j,'mask':c['centre']['mask'],**v});counts[v['status']]+=1
            signs.append({'sign':sign,'visible_charts':sum(v['status']=='VISIBLE_AND_RECORDED' for v in observations),'observations':observations})
        rows.append({'old_candidate_index':index,'point':[str(x),str(y)],'signs':signs,'visible_in_any_metric_box':any(r['visible_charts'] for r in signs)})
    return {'schema':'elliptic-curves.native11952-metric-visibility.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/search_observability.py',OLD,PARENT,TAIL,COVERAGE)},
        'same43_masks':True,'representatives_equal_up_to_sign':comparisons,'old_candidate_count':len(rows),'old_candidates_visible_in_metric_boxes':sum(r['visible_in_any_metric_box'] for r in rows),
        'observation_counts':dict(sorted(counts.items())),'rows':rows,
        'claim_boundary':'Retrospective comparison only, after the generic-point-only metric control and its fixed tails terminated. Same43 parity masks and specialized representatives up to sign, but the reduced finite boxes and implementations differ. Outside-box points do not exclude other representatives of their Mordell-Weil cosets or unseen points. No new curve or rank record.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args();data=expected()
    if a.check:
        if cert.read(OUTPUT)!=data:raise ArithmeticError('visibility proof differs')
    else:
        if OUTPUT.exists():raise FileExistsError('preserve visibility proof')
        checkpoint(OUTPUT,data)
    print('NATIVE11952 VISIBILITY',data['old_candidates_visible_in_metric_boxes'],'of',data['old_candidate_count'],data['observation_counts'],flush=True)
