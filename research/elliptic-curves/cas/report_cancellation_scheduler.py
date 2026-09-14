#!/usr/bin/env python3
"""Frozen paired comparison including independent rank certification costs."""
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from cancellation_scheduler_prepare import OUT, RAW, ARMS
from finite_cancellation_corpus import digest, pointkey, short, write


def summarize(rows):
    return {arm:{'cases':sum(r['arm']==arm for r in rows),
        'gains':sum(r['success'] for r in rows if r['arm']==arm),
        'cpu_seconds':sum(r['cpu_seconds'] for r in rows if r['arm']==arm),
        'calls':sum(r['calls'] for r in rows if r['arm']==arm),
        'literal_hit_cases':sum(r['literal_hit'] for r in rows if r['arm']==arm)} for arm in ARMS}


def compare(rows, candidate, reference):
    table={(r['case'],r['arm']):r for r in rows};groups=defaultdict(list)
    for r in rows:
        if r['arm']==reference:groups[r['family'],r['stratum']].append(r['case'])
    groups=[sorted(ids) for k,ids in sorted(groups.items())]
    totals=summarize(rows);a,b=totals[candidate],totals[reference]
    ratio=(a['gains']/a['cpu_seconds'])/(b['gains']/b['cpu_seconds']) if b['gains'] else None
    rng=np.random.default_rng(20260914);ratios=[];undefined=0
    for _ in range(10000):
        sample=[x for ids in groups for x in rng.choice(ids,len(ids),replace=True)]
        ca=sum(table[x,candidate]['cpu_seconds'] for x in sample);cb=sum(table[x,reference]['cpu_seconds'] for x in sample)
        ga=sum(table[x,candidate]['success'] for x in sample);gb=sum(table[x,reference]['success'] for x in sample)
        if gb:ratios.append((ga/ca)/(gb/cb))
        else:undefined+=1
    interval=list(map(float,np.quantile(ratios,[.0125,.9875]))) if ratios and not undefined else None
    gate=bool(a['gains']>=b['gains'] and ratio>=1.1 and interval[0]>1) if ratio is not None and interval else None
    return {'candidate':candidate,'reference':reference,'recoveries_per_cpu_ratio':ratio,
        'paired_stratified_ci975':interval,'undefined_bootstrap_samples':undefined,'gate_passed':gate,
        'candidate_only':sum(table[x,candidate]['success'] and not table[x,reference]['success'] for ids in groups for x in ids),
        'reference_only':sum(table[x,reference]['success'] and not table[x,candidate]['success'] for ids in groups for x in ids),
        'boundary':'Conditional paired bootstrap over the declared family/stratum roster; central 97.5 percent interval for each of two predeclared promotion comparisons. No universal policy or population theorem.'}


def main():
    plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
    supervision=json.loads((OUT/'supervision.json').read_text())
    if supervision['status']!='COMPLETE':raise ArithmeticError('incomplete benchmark')
    if len(supervision['records'])!=plan['cases']*len(ARMS):raise ArithmeticError('missing arms')
    inputs={c['id']:c for c in json.loads((OUT/'inputs.json').read_text())}
    oracle={o['id']:o for o in json.loads((OUT/'oracle.json').read_text())}
    rows=[];identities=set();unknowns=0
    for receipt in supervision['records']:
        ci,arm=receipt['case'],receipt['arm'];ident=(ci,arm)
        if ident in identities:raise ArithmeticError('duplicate arm')
        identities.add(ident)
        if receipt['status']!='COMPLETE' or receipt['protocol_sha256']!=ph:
            raise ArithmeticError('incomplete CPU or protocol binding')
        dest=RAW/'arms'/ci/arm;result=json.loads((dest/'result.json').read_text())
        verification=json.loads((dest/'independent-verification.json').read_text())
        if digest((dest/'result.json').read_bytes())!=receipt['result_sha256']:
            raise ArithmeticError('result changed')
        if digest((dest/'independent-verification.json').read_bytes())!=result['independent_verification_sha256']:
            raise ArithmeticError('independent verification changed')
        if verification['status']!='PASS' or result['success']!=result['candidate']:
            raise ArithmeticError('uncertified gain')
        unknowns+=verification['preparation_unknowns']
        case=inputs[ci];seed=case['seed'];_,base=short(seed['curve'],seed['points'])
        held={pointkey(p) for p in oracle[ci]['points']}-{pointkey(p) for p in base}
        hit=False
        for call in result['calls']:
            packet=json.loads((dest/call['file']).read_text())
            points=[[p['x'],p['y']] for p in packet['search']['finite_curve_points']]
            _,normalized=short(seed['curve'],points)
            hit=hit or bool({pointkey(p) for p in normalized}&held)
        rows.append({'case':ci,'arm':arm,'family':case['family'],'stratum':case['stratum'],
            'success':result['success'],'cpu_seconds':receipt['charged_cpu_seconds'],'calls':len(result['calls']),
            'initial_rank':result['initial_rank'],'rank_lower_bound':result['rank_lower_bound'],
            'literal_hit':hit,'independent_replay_cpu':result['proof_and_independent_replay_cpu_seconds'],
            'search_phase_cpu':result['search_phase_cpu_seconds'],
            'status':result['status']})
    scopes={}
    for name in ('all','shallow','deep'):
        selected=rows if name=='all' else [r for r in rows if r['stratum']==name]
        if not selected:continue
        scopes[name]={'totals':summarize(selected),'comparisons':[
            compare(selected,a,b) for a,b in [('adaptive_local','factor_free'),
                                            ('adaptive_local','adaptive_uniform'),
                                            ('adaptive_uniform','factor_free')]]}
    promoted=all(c['gate_passed'] is True for c in scopes['all']['comparisons'][:2]) and unknowns==0
    offline={n:json.loads((OUT/n).read_text())['cpu_seconds'] for n in ('fit.json','preflight.json','training-replay.json')}
    totals=scopes['all']['totals'];a=totals['adaptive_local']
    conservative={b:((a['gains']/(a['cpu_seconds']+sum(offline.values())))/
                    (totals[b]['gains']/totals[b]['cpu_seconds'])) if totals[b]['gains'] else None
                  for b in ('factor_free','adaptive_uniform')}
    with (OUT/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary={'status':'COMPLETE_FROZEN_ADAPTIVE_COMPARISON','scopes':scopes,'rows':rows,
        'promote_local_policy':promoted,'fresh_fibre_gate':'PASS' if promoted else 'NOT_PASSED',
        'preparation_unknowns':unknowns, 'offline_cpu_seconds':offline,
        'conservative_campaign_rate_ratios':conservative,
        'campaign_cost_boundary':'Sensitivity charges ALL recorded fitting and shared preflight CPU to the local candidate alone, in addition to its complete process CPU. Inherited historical corpus creation, development and retained landscape construction are outside this input interface and have no invented zero cost.',
        'protocol_sha256':ph,'files_sha256':{n:digest((OUT/n).read_bytes()) for n in ('protocol.json','fit.json','preflight.json','supervision.json','comparison.csv')},
        'rank32':'UNKNOWN','new_record':False,
        'boundary':'Known withheld controls only. Negative policy results do not imply that every cancellation policy fails. Rank upper bounds, saturation and fresh-fibre transfer remain unknown.'}
    write(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k in ('status','scopes','promote_local_policy','offline_cpu_seconds','conservative_campaign_rate_ratios')},indent=2),flush=True)


if __name__=='__main__':main()
