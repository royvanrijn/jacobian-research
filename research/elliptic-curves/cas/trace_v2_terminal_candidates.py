#!/usr/bin/env python3
"""Trace retrospectively vetted cosets through every immutable V2 landscape."""
import csv
import json
from pathlib import Path
import numpy as np
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'artifacts/local/elliptic-curves/v2-terminal-retrospective-v1'
V2=ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2'
def read(p):return json.loads(p.read_text())
def fp(word,columns):
    value=0
    for x,c in zip(word,columns):
        if int(x)%2:value^=c
    return value


def trace(witness):
    p=[int(x)%2 for x in witness['centre_word']]
    orbit=None
    with (ROOT/'artifacts/generated-results/elliptic-curves/curve302_parent_degree2_multisection_orbits_v1.tsv').open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            if [int(x)%2 for x in row['parent_MW17_w'].split()]==p[:17]:
                orbit={k:row[k] for k in ('orbit_mask','minimum_norm','category')};break
    stages=[]
    for wd in sorted((V2/'calibration302').glob('epoch-*')):
        s=read(wd/'selection.json');r=s['rank']
        record={'rank':r,'eligible_in_embedded_stage_parity_space':not any(p[r:]),'scored':False}
        matches=[(i,a) for i,a in enumerate(s['anchors']) if [x%2 for x in a['anchor']['representative'][:17]]==p[:17]]
        record['retained_base_anchor']=bool(matches)
        if matches and not any(p[r:]):
            ai,a=matches[0];ext=sum(p[17+i]<<i for i in range(r-17));key=fp(p[:r],s['fingerprint_columns'])
            data=np.load(wd/f'anchor-{ai:02d}-full.npz');norm=int(data['babai_norms'][ext])
            def better(data):
                norms=data['babai_norms'];ties=np.flatnonzero(norms==norm)
                return int(np.sum(norms>norm))+sum(fp(data['residues'][int(i)],s['fingerprint_columns'])<key for i in ties)
            local=better(data)+1;global_rank=sum(better(np.load(wd/f'anchor-{j:02d}-full.npz')) for j in range(32))+1
            refined=[row for row in a['refined'] if row['extension']==ext]
            finalists=[row for row in s['centres'] if row['fingerprint']==key]
            built=[str(path.relative_to(ROOT)) for path in sorted(wd.glob('chart-*.json')) if read(path)['centre']['fingerprint']==key]
            record.update(scored=True,extension_encoding=ext,babai_norm=norm,rank_within_anchor=local,
                rank_across_scored_landscape=global_rank,exact_CVP_shortlisted=bool(refined),
                exact_norm=refined[0]['metric_norm'] if refined else None,
                finalist=bool(finalists),built_same_coset_charts=built)
        stages.append(record)
    return {'witness':witness,'base_orbit':orbit,'parity':p,'stages':stages}


def main():
    trials=[read(p) for p in sorted((D/'terminal-landscape').glob('trial-*.json'))]
    finite=[w for r in trials for w in r.get('witnesses',[]) if int(w['coordinate'][1])!=0]
    winner=min(finite,key=lambda w:int(w['height']))
    bank=[read(p) for p in sorted((D/'retained-bank').glob('anchor-*.json'))]
    assert len(bank)==32,'wait for full retained-bank diagnostic'
    bank_winner=min((w for a in bank for c in a['candidates'] for w in c['witnesses'] if int(w['coordinate'][1])!=0),key=lambda w:int(w['height']))
    searched=read(D/'terminal-searched-charts/result.json')
    result={'globally_optimal_finite_height_witness':trace(winner),
        'best_vetted_retained_bank_witness':trace(bank_winner),
        'best_vetted_searched_chart':searched['best'],
        'global_height_lower_bound':1,'global_height_optimum_proved':int(winner['height'])==1,
        'terminal_chart_count':searched['count'],'terminal_affine_box_hits':searched['cheap']}
    checkpoint(D/'gate-trace.json',result)
    for name in ('globally_optimal_finite_height_witness','best_vetted_retained_bank_witness'):
        row=result[name];print(name,row['witness']['height'],row['base_orbit'],row['stages'][-1])


if __name__=='__main__':main()
