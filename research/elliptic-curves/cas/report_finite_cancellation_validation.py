#!/usr/bin/env python3
"""Frozen aggregate/per-stratum comparisons with paired case uncertainty."""
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from finite_cancellation_corpus import digest,write
from finite_cancellation_validation_audit import OUT

ARMS=('factor_free','real_q','real_q_g')


def totals(rows):
    return {arm:{'cases':sum(r['arm']==arm for r in rows),
        'recoveries':sum(bool(r['success']) for r in rows if r['arm']==arm),
        'literal_hit_cases':sum(bool(r['literal_withheld_points_recovered']) for r in rows if r['arm']==arm),
        'cpu_seconds':sum(r['charged_cpu_seconds'] for r in rows if r['arm']==arm),
        'point_calls':sum(r['point_calls'] or 0 for r in rows if r['arm']==arm),
        'unknown_worker_arms':sum(r['verification']!='PASS_EXACT_REPLAY' for r in rows if r['arm']==arm),
        'cpu_components':{key:sum(r.get('cpu_components',{}).get(key,0.) for r in rows if r['arm']==arm)
            for key in ['initial_verification_cpu_seconds','map_and_feature_cpu_seconds','point_backend_cpu_seconds','admission_cpu_seconds','rank_proof_cpu_seconds']}}
        for arm in ARMS}


def comparison(rows,candidate,reference):
    tab={(r['case'],r['arm']):r for r in rows};groups=defaultdict(list)
    for r in rows:
        if r['arm']==reference:groups[r['stratum'],r['family']].append(r['case'])
    groups=[sorted(v) for k,v in sorted(groups.items())];ids=sorted(c for g in groups for c in g)
    sums=totals(rows);a,b=sums[candidate],sums[reference]
    ratio=(a['recoveries']/a['cpu_seconds'])/(b['recoveries']/b['cpu_seconds']) if b['recoveries'] else None
    rng=np.random.default_rng(20260914);bs=[];cpu_bs=[];undefined=0
    for _ in range(10000):
        sample=[x for group in groups for x in rng.choice(group,len(group),replace=True)]
        ac=sum(tab[x,candidate]['charged_cpu_seconds'] for x in sample);bc=sum(tab[x,reference]['charged_cpu_seconds'] for x in sample)
        ag=sum(tab[x,candidate]['success'] for x in sample);bg=sum(tab[x,reference]['success'] for x in sample)
        cpu_bs.append(ac/bc)
        if bg:bs.append((ag/ac)/(bg/bc))
        else:undefined+=1
    ci=[float(x) for x in np.quantile(bs,[.025,.975])] if bs and not undefined else None
    finite_cost_ratio=a['cpu_seconds']/b['cpu_seconds']
    gate=None if ratio is None or ci is None else bool(a['recoveries']>=b['recoveries'] and ratio>=1.1 and ci[0]>1)
    return {'candidate':candidate,'reference':reference,'recovery_per_cpu_ratio':ratio,
        'paired_stratified_bootstrap_ci95':ci,'bootstrap_undefined_rate_samples':undefined,
        'cpu_ratio':finite_cost_ratio,'cpu_ratio_bootstrap_ci95':[float(x) for x in np.quantile(cpu_bs,[.025,.975])],
        'primary_gate_passed':gate,
        'original_style_time_gate_passed':bool(a['recoveries']>=b['recoveries'] and finite_cost_ratio<=.9),
        'paired_disagreements':{'candidate_only':sum(tab[x,candidate]['success'] and not tab[x,reference]['success'] for x in ids),
            'reference_only':sum(tab[x,reference]['success'] and not tab[x,candidate]['success'] for x in ids)},
        'cases_faster':sum(tab[x,candidate]['charged_cpu_seconds']<tab[x,reference]['charged_cpu_seconds'] for x in ids),
        'cases':len(ids),'scope':'Conditional on the retained family/stratum roster, 10000 paired case-bootstrap draws. No familywise correction or population-rank theorem.'}


def main():
    verified=json.loads((OUT/'verified.json').read_text());rows=verified['rows']
    scopes={}
    for name in ('all','shallow','deep'):
        rs=rows if name=='all' else [r for r in rows if r['stratum']==name]
        scopes[name]={'totals':totals(rs),'comparisons':{a+'_vs_'+b:comparison(rs,a,b)
            for a,b in [('real_q','factor_free'),('real_q_g','factor_free'),('real_q_g','real_q')]}}
    with (OUT/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['case','family','stratum','initial_rank','arm','success','literal_withheld_points_recovered','point_calls','charged_cpu_seconds','status'])
        writer.writeheader()
        for row in rows:writer.writerow({k:row[k] for k in writer.fieldnames})
    offline={}
    for filename in ('preflight.json','fit.json','training-verified.json','verified.json'):
        offline[filename]=json.loads((OUT/filename).read_text())['cpu_seconds']
    summary={'status':'COMPLETE_FROZEN_THREE_POLICY_VALIDATION','scopes':scopes,
        'cases':verified['cases'],'arms':len(rows),'point_calls':verified['point_calls'],
        'independent_gain_proofs':verified['independent_gain_proofs'],'residue_leaves':verified['residue_leaves'],
        'offline_cpu_seconds':offline,'source_sha256':digest(Path(__file__).read_bytes()),
        'files_sha256':{p.name:digest(p.read_bytes()) for p in sorted(OUT.iterdir()) if p.is_file() and p.name not in ('summary.json','comparison.png','comparison.pdf')},
        'boundary':'A separate clean ablation and retrospective validation. Historical V1 results and gate remain unchanged. Both fitted arms use identical pure-q and real features; only the full arm computes extra gcd features. Initial ranks M18-M19 and M24-M25, not M27+. Conditional successes and CPU costs are not a new rank record or population speedup theorem. Offline fitting/input audit/independent replay costs are shown separately from per-arm process-tree CPU.'}
    write(OUT/'summary.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12.5,4.1),layout='constrained')
    labels=['V3','Real + q','Real + q + gcd'];colors=['#64748b','#2563eb','#b45309']
    for ax,name in zip(axes,('all','shallow','deep')):
        ts=scopes[name]['totals'];values=[100*ts[a]['recoveries']/ts[a]['cpu_seconds'] for a in ARMS]
        bars=ax.bar(labels,values,color=colors,width=.68)
        for bar,a in zip(bars,ARMS):
            t=ts[a];ax.annotate(f"{t['recoveries']}/{t['cases']} gains\n{t['cpu_seconds']:.1f} CPU s",(bar.get_x()+bar.get_width()/2,bar.get_height()),xytext=(0,6),textcoords='offset points',ha='center',fontsize=9)
        ax.set_title({'all':'All 41 controls','shallow':'M18–M19 controls','deep':'M24–M25 controls'}[name])
        ax.set_ylim(0,max(values+[.1])*1.35);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Certified recoveries per 100 CPU seconds')
    fig.suptitle('Finite cancellation: separate q-only features and complete search costs',fontsize=14)
    fig.savefig(OUT/'comparison.png',dpi=180);fig.savefig(OUT/'comparison.pdf');plt.close(fig)
    summary['files_sha256'].update({n:digest((OUT/n).read_bytes()) for n in ('comparison.png','comparison.pdf')})
    write(OUT/'summary.json',summary)
    print(json.dumps({'status':summary['status'],'scopes':scopes,'offline_cpu_seconds':offline},indent=2))


if __name__=='__main__':main()
