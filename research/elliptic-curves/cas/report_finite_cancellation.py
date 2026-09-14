#!/usr/bin/env sage -python
"""Compact tables, exported figure and manifest from completed certificates."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from finite_cancellation_corpus import OUT, ROOT, digest, write

def main():
    read=lambda name:json.loads((OUT/name).read_text())
    a=read('analysis.json');c=read('census.json');cpu=read('cpu/verified.json');s=read('neighbour_states.json');q=read('q_ablation.json')
    bycase={}
    for r in cpu['rows']:bycase.setdefault(r['case'],{})[r['arm']]=r
    rows=[]
    for case,pair in bycase.items():
        b,n=pair['factor_free'],pair['finite_selector']
        rows.append({'case':case,'family':b['family'],'initial_rank':b['initial_rank'],
            'baseline_gain':int(b['success']),'selector_gain':int(n['success']),
            'baseline_literal_withheld_hits':b['literal_withheld_points_recovered'],'selector_literal_withheld_hits':n['literal_withheld_points_recovered'],
            'baseline_calls':b['point_calls'],'selector_calls':n['point_calls'],
            'baseline_cpu':b['charged_cpu_seconds'],'selector_cpu':n['charged_cpu_seconds']})
    with (OUT/'cpu/comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    rng=np.random.default_rng(20260914);times=np.array([[r['baseline_cpu'],r['selector_cpu']] for r in rows])
    boot=[]
    for _ in range(10000):
        sums=times[rng.integers(len(times),size=len(times))].sum(axis=0);boot.append(sums[1]/sums[0])
    interval=[float(x) for x in np.quantile(boot,[.025,.975])]
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.4),layout='constrained')
    plt.rcParams.update({'font.size':10})
    keys=['root_ball','outside_ball'];entries=[s['probability_neighbour_has_smaller_literal_H'][k] for k in keys]
    vals=[100*r['estimate'] for r in entries];err=np.array([[100*(r['estimate']-r['ci95'][0]) for r in entries],[100*(r['ci95'][1]-r['estimate']) for r in entries]])
    axes[0].bar(['Distinguished ball','Outside ball'],vals,color=['#247c78','#9da7b0'],yerr=err,capsize=3)
    axes[0].set(ylabel='Neighbour has smaller H (%)',ylim=(0,110),title='Exact finite-place distinction')
    for i,y in enumerate(vals):axes[0].text(i,y+3,f'{y:.1f}%',ha='center')
    entries=[a['outcomes']['factor_free']['H_squared_ratio'],a['outcomes']['ridge_real']['H_squared_ratio'],
             q['results']['five_j_folds']['real_plus_q_vs_factor_free'],a['outcomes']['ridge_real_and_finite']['H_squared_ratio']]
    axes[1].bar(['V3','Real','Real + q','Full finite'],[r['estimate'] for r in entries],color=['#9da7b0','#9da7b0','#7eacab','#247c78'])
    axes[1].set(ylabel='Geometric H² ratio to V3',title='Withheld model accessibility',ylim=(0,1.15))
    for i,r in enumerate(entries):axes[1].text(i,r['estimate']+.03,f'{r["estimate"]:.3f}',ha='center')
    b,n=cpu['totals']['factor_free'],cpu['totals']['finite_selector']
    axes[2].scatter([b['cpu_seconds'],n['cpu_seconds']],[b['recoveries'],n['recoveries']],s=90,c=['#7d8995','#247c78'])
    axes[2].annotate('Factor-free V3',(b['cpu_seconds'],b['recoveries']),xytext=(-8,-20),textcoords='offset points',ha='right')
    axes[2].annotate('Finite selector',(n['cpu_seconds'],n['recoveries']),xytext=(8,8),textcoords='offset points')
    axes[2].set(xlabel='Charged CPU seconds (lower is better)',ylabel='Next-direction recoveries / 12',xlim=(170,215),ylim=(8.5,10.7),yticks=[9,10],title='Completed V3 controls')
    for ax in axes:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    fig.suptitle('Finite cancellation predicts accessibility; the strict CPU gate remains unmet',fontsize=13)
    fig.savefig(OUT/'comparison.png',dpi=180);fig.savefig(OUT/'comparison.pdf');plt.close(fig)
    names=['protocol.json','census.json','corpus.json.gz','analysis.json','evaluation.json.gz','replay.json','prime_decomposition.json',
           'g_factorizations.json.gz','q_ablation.json','neighbour_states.json','cpu/protocol.json','cpu/inputs.json','cpu/supervision.json',
           'cpu/verified.json','cpu/comparison.csv','comparison.png','comparison.pdf']
    result={'status':'COMPLETE_BOUNDED_PREDICTOR_STUDY_STRICT_CPU_GATE_MISSED','curve_generic_banks':c['curve_generic_banks'],
        'j_groups':c['j_groups'],'exceptional_point_observations':c['point_observations'],'registry_endpoints_joined':c['registry_generic_joins'],
        'exact_signed_model_evaluations':a['exact_signed_model_evaluations'],'primary_accessibility_gate_passed':a['passed_predictors'],
        'cpu_totals':cpu['totals'],'cpu_ratio':cpu['cpu_ratio'],'strict_cpu_gate_passed':False,
        'cpu_ratio_case_bootstrap_ci95_exploratory':interval,'gain_per_cpu_ratio':n['recoveries']/n['cpu_seconds']/(b['recoveries']/b['cpu_seconds']),
        'literal_withheld_point_case_rate_per_cpu_ratio':n['cases_with_literal_withheld_recovery']/n['cpu_seconds']/(b['cases_with_literal_withheld_recovery']/b['cpu_seconds']),
        'files_sha256':{name:digest((OUT/name).read_bytes()) for name in names},
        'policy_sha256':digest((ROOT/'elliptic-curves/data/finite_cancellation_policy_v1.json').read_bytes()),
        'reporter_sha256':digest(Path(__file__).read_bytes()),
        'boundary':'Positive conditional predictive signal and a measured12-case recovery/CPU advantage;9.987859% CPU saving is below the frozen10% cutoff, and the exploratory case bootstrap includes1. Opt-in experimental policy, no default replacement or rank32 claim. No C/L optimization.'}
    write(OUT/'summary.json',result);print(json.dumps(result))

if __name__=='__main__':main()
