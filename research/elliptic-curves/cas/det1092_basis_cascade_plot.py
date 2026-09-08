#!/usr/bin/env python3
"""Static scientific plot from terminal, independently replayed cascade data."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];P=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
s=json.loads((P/'summary.json').read_text());assert all(r['replay_status']=='PASS_EXACT_REPLAY' for r in s['arms'])
fig,ax=plt.subplots(figsize=(9,4.8));styles={'original':('#0072B2','-','Original'),'recovered-P14':('#7B61A8','--','Recovered P14'),'recovered-P11':('#D55E00','-','Recovered P11'),'random-unimodular':('#009E73','-','Random unimodular')}
for r in s['arms']:
 x=[0]
 for c in r['charts_per_stage']:x.append(x[-1]+c)
 color,linestyle,label=styles[r['arm_id']]
 ax.step(x,r['rank_path'],where='post',color=color,linestyle=linestyle,linewidth=1.8,label=label)
 ax.plot(x[-1],r['rank_path'][-1],'o',color=color,markersize=4)
ax.set(xlabel='Charts executed',ylabel='Rank certified from retained points',title='Curve302: reconstruction cost from four MW17 bases')
ax.set_yticks(range(17,max(r['final_rank_lower_bound'] for r in s['arms'])+1));ax.set_ylim(16.7,max(r['final_rank_lower_bound'] for r in s['arms'])+.5);ax.set_xlim(left=0)
ax.grid(axis='y',color='#e4e4e4',linewidth=.6);ax.spines[['top','right']].set_visible(False);ax.legend(loc='lower right',frameon=False)
fig.text(.12,.045,'Same curve and initial integral subgroup; fresh searches; height bound125000.',fontsize=8)
fig.text(.12,.018,'Starting commit: '+s['starting_commit'],fontsize=6.5)
fig.tight_layout(rect=[0,.07,1,1])
for ext in ['png','svg']:
 target=P/('rank_path_comparison.'+ext)
 if target.exists():raise FileExistsError('preserve rendered figure')
 meta={'Software':'det1092_basis_cascade_plot.py','Starting commit':s['starting_commit']} if ext=='png' else {'Creator':'det1092_basis_cascade_plot.py','Description':'Starting commit: '+s['starting_commit'],'Date':None}
 fig.savefig(target,dpi=180,metadata=meta)
print('WROTE rank_path_comparison.png/svg')
