#!/usr/bin/env sage-python
"""Retrospective omission states in the exact V2-completed rank31 basis."""
import argparse
import itertools
import json
import hashlib
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import numpy as np
from sage.all import matrix, vector, ZZ
from research_runtime.store import checkpoint

CAS=Path(__file__).resolve().parent
l=SourceFileLoader('terminal_landscapes',str(CAS/'v2_terminal_landscape.sage')).load_module()
b=l.b;D=l.D
a2=b.load('audit_curve302_exceptional_subgroup_landscape.sage')


def plan():
    if (D/'state-plan.json').exists():raise FileExistsError('preserve state plan')
    data=b.read(D/'coordinates.json');ambient=tuple(tuple(map(F,p)) for p in data['ambient_points']);model=tuple(map(F,data['model']))
    fixed=b.read(D/'fixed-M30.json');points=tuple(tuple(map(F,p)) for p in fixed['points'])+(tuple(map(F,data['missing_point'])),)
    rows=matrix(ZZ,data['M30_words']).stack(matrix.identity(ZZ,31)[data['missing_axis_index']:data['missing_axis_index']+1,:])
    assert abs(rows.det())==1
    h,asym=a2.rounded_metric(l.geo,model,ambient)
    table=b.read(a2.OUTPUT)
    digest=hashlib.sha256(json.dumps([list(map(int,row)) for row in h.rows()],separators=(',',':')).encode()).hexdigest()
    assert digest==table['fixed_basis']['rounded_metric_sha256']
    array=np.asarray(h.rows(),dtype=np.int64)
    states=[]
    for size in (1,2):
        for omitted in itertools.combinations(range(17,31),size):
            kept=[i for i in range(31) if i not in omitted]
            w=rows[kept,:];prepared=a2.prepare_subgroup(h,array,w.rows())
            g=prepared['gram'];scores=[]
            for i in omitted:
                target=rows.row(i);cross=w*h*target;hh=int(target*h*target)
                t=2*np.linalg.solve(prepared['gram_array'].astype(float),np.asarray(cross,dtype=float))
                seed=a2.nearest_plane(prepared['upper'],t)
                word,cost=a2.exact_coordinate_descent(g,cross,seed,hh)
                scores.append({'omitted':i,'fresh_common_metric_numerator':cost,'nearest_word':word})
            # A coordinate-subset state matches only if its rational span is
            # spanned by exactly 31-size ambient coordinate axes. Since rows
            # are a unimodular subset, equal rational spans then imply equality.
            span=w.row_space();axes=[i for i in range(17,31) if matrix.identity(ZZ,31).row(i) in span]
            match=sum(1<<(i-17) for i in axes) if len(axes)==14-size else None
            states.append({'id':'omit-'+'-'.join(str(i-17) for i in omitted),'omitted':list(omitted),'kept':kept,
                'rank':len(kept),'scores':scores,'score_max':max(r['fresh_common_metric_numerator'] for r in scores),
                'matching_agent2_state':match,'contained_agent2_axes':axes})
    rank29=sorted((s for s in states if s['rank']==29),key=lambda s:(s['score_max'],s['id']))
    selected={s['id'] for s in rank29[:12]}
    for state in (s for s in states if s['rank']==30):
        children=[s for s in rank29 if set(state['omitted'])<=set(s['omitted'])]
        selected.add(min(children,key=lambda s:(s['score_max'],s['id']))['id'])
    # Include the actually visited M29 predecessor even if the proxy ranks it poorly.
    selected.add('omit-12-13')
    payload={'B31_points':[list(map(str,p)) for p in points],'B31_ambient_words':[list(map(int,r)) for r in rows.rows()],
        'common_ambient_metric':[list(map(int,r)) for r in h.rows()], 'common_metric_sha256':digest,
        'states':states,'selected_rank29_ids':sorted(selected),
        'bounds':'Every one-point omission (14 including terminal M30) and every two-point omission (91) is scored in Agent2 common rounded metric. Vet all 13 alternative M30 states with exact k=32 unrestricted half-lattice enumeration per missing direction. For M29 vet top12 by worst missing-direction surrogate, best child of each M30, and actual V2 predecessor; k=16 per missing direction. No oracle-guided prospective execution.',
        'V2_M30_primitive_normal':list(map(int,matrix(ZZ,data['M30_words']).right_kernel().basis()[0]))}
    checkpoint(D/'state-plan.json',payload);print('PLANNED',len(states),'states; selected M29',len(selected),flush=True)


def run():
    plan=b.read(D/'state-plan.json');data=b.read(D/'coordinates.json');model=tuple(map(F,data['model']))
    selected=set(plan['selected_rank29_ids'])
    if (D/'exact-predecessor-plan.json').exists():selected.update(b.read(D/'exact-predecessor-plan.json')['selected_ids'])
    allpts=tuple(tuple(map(F,p)) for p in plan['B31_points'])
    for state in plan['states']:
        if state['id']=='omit-13':continue
        if state['rank']==29 and state['id'] not in selected:continue
        points=tuple(allpts[i] for i in state['kept'])
        for i in state['omitted']:
            result=l.vet(points,allpts[i],model,D/'states'/state['id']/('target-%d'%i),32 if state['rank']==30 else 16)
            print('STATE',state['id'],'target',i,'best',result['best']['height'],flush=True)


def exact_plan():
    path=D/'exact-predecessor-plan.json'
    if path.exists():raise FileExistsError('preserve exact predecessor plan')
    plan=b.read(D/'state-plan.json')
    rows=[b.read(D/'common-exact'/(s['id']+'.json')) for s in plan['states'] if s['rank']==29]
    rows.sort(key=lambda r:(max(x['minimum_numerator'] for x in r['targets'].values()),r['state']['id']))
    selected={r['state']['id'] for r in rows[:12]}
    for s in plan['states']:
        if s['rank']==30:selected.add(next(r['state']['id'] for r in rows if set(s['omitted'])<=set(r['state']['omitted'])))
    checkpoint(path,{'selected_ids':sorted(selected),'additional_ids':sorted(selected-set(plan['selected_rank29_ids'])),
        'rule':'After exact common-metric optimization of all91 predecessors, also vet top12 by maximum of the two exact missing-direction minima and best child of every M30. Keep all earlier surrogate-selected results. Still retrospective, not a V2 or V3 selection rule.'})


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['plan','run','exact-plan']);a=ap.parse_args()
    if a.action=='plan':plan()
    elif a.action=='exact-plan':exact_plan()
    else:run()
