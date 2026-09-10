#!/usr/bin/env python3
"""Read-only score/cost/rank report, preserving unfinished and censored work."""
import argparse
import json
from pathlib import Path

def read(p):return json.loads(p.read_text())
def report(folder):
    plan=read(folder/'plan.json');root=Path(plan['root'])
    rows=[];windows=[]
    for w in sorted((root/'ordinary-search').glob('window-*')):
        if not (w/'scores.json').exists():continue
        scores=read(w/'scores.json');windows.append({'window':w.name,'scored':len(scores['rows']),
          'triage_cpu_seconds':scores['cold_cpu_seconds']+scores['overhead_cpu_seconds']+sum(r['score_direct_cpu_seconds'] for r in scores['rows'])})
        for case in sorted((w/'cases').glob('*')):
            selection=read(case/'selection.json');terminal=read(case/'terminal.json') if (case/'terminal.json').exists() else None
            costs=[read(p) for p in case.glob('batch-*/cost.json')]
            stages=[read(p) for p in case.glob('batch-*/result.json')]
            rank=max([r.get('rank_lower_bound') or 0 for r in stages],default=0) or None
            open_batches=sum(not (p/'cost.json').exists() for p in case.glob('batch-*'))
            rows.append({**selection,'rank_lower_bound':rank,'status':terminal['status'] if terminal else ('PAUSED' if (folder/'STOP').exists() and not open_batches else 'IN_PROGRESS'),
              'calls':sum(r.get('calls',0) for r in stages),
              'closed_batch_search_cpu_seconds':sum(r['total_cpu_seconds'] for r in costs),
              'open_batches':open_batches,
              'total_cpu_seconds':terminal['total_cpu_seconds'] if terminal else None})
    summary={}
    for name,control in [('ranked',False),('control',True)]:
        arm=[r for r in rows if r['control']==control]
        summary[name]={'dispatched':len(arm),'completed':sum(r['status']=='COMPLETE_BOUNDED' for r in arm),
          'censored_or_unresolved':sum(r['status'] not in ('COMPLETE_BOUNDED','IN_PROGRESS') for r in arm),
          'threshold_counts':{str(k):sum((r['rank_lower_bound'] or 0)>=k for r in arm) for k in (20,23,27,28,29,30,31,32)},
          'closed_batch_cpu_seconds':sum(r['closed_batch_search_cpu_seconds'] for r in arm)}
    return {'status':read(folder/'STATUS.json') if (folder/'STATUS.json').exists() else 'STARTING',
      'windows':windows,'arms':summary,'parameters':rows,'preparation_cost':read(folder/'preparation-cost.json'),
      'boundary':'Open-batch CPU is pending, not zero. Controller CPU includes its scoring work: do not add both twice. Scores are scheduling signals. Adaptive rank-based exposure and incomplete cohorts must be accounted for before claiming triage efficacy.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path)
    a=p.parse_args();result=report(a.folder);data=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if a.output:a.output.write_text(data)
    print(json.dumps({k:v for k,v in result.items() if k!='parameters'},indent=2))
