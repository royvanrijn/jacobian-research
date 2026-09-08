#!/usr/bin/env python3
"""Retrospective V3 universality panel over the 12 remaining curve-302 exceptional seeds.

The two strongest seeds, recovered-strict-02 and recovered-strict-03, already have
independently replayed 18->31 amplifier runs.  This panel asks the complementary
finite question: if each of the other known M31/M17 directions is supplied as the
single 18th generator, does the *unchanged* frozen V3 policy amplify it, and how far?

Oracle use is restricted to constructing each known rank-18 seed.  After the seed
input/proof/protocol are frozen, search uses the same numerical V3 policy as the
completed two-seed experiment.  This is retrospective mechanism analysis, not a
prospective selector or a new-rank claim.

Commands (standard Python): preflight | launch | resume | status | diagnose | stop | worker
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'
D = LOCAL/'curve302-seed-universality-panel-v1'
STATE = D/'state.json'
LOG = D/'worker.log'
SELF = Path(__file__).resolve()
CLOSURE = ART/'curve302_unlock_seed_closure_v2.json'
KNOWN_WINNERS = ('recovered-strict-02','recovered-strict-03')
EXPECTED_DIRECTIONS = (
    'recovered-local-01','recovered-local-02','recovered-local-03','recovered-local-04',
    'recovered-strict-01','recovered-strict-02','recovered-strict-03',
    'residual-strict-01','residual-strict-02','residual-strict-03','residual-strict-04',
    'residual-strict-05','residual-strict-06','residual-strict-07',
)


def read(path): return json.loads(Path(path).read_text())
def atomic(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name('.'+path.name+'.tmp-'+str(os.getpid()))
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n'); os.replace(tmp,path)

def process_alive(pid):
    if not isinstance(pid,int) or pid<=0:return False
    try: os.kill(pid,0)
    except (ProcessLookupError,PermissionError): return False
    try:
        stat=Path(f'/proc/{pid}/stat').read_text().split()
        return len(stat)>2 and stat[2]!='Z'
    except Exception:return True

def sage_launcher():
    for value in (os.environ.get('V3_SAGE'),shutil.which('sage'),str(Path.home()/'.local/bin/sage'),'/usr/bin/sage'):
        if value and Path(value).is_file() and os.access(value,os.X_OK):return str(Path(value).resolve())
    raise RuntimeError('Sage launcher not found; set V3_SAGE=/absolute/path/to/sage')

def tail(path,n=16000):
    path=Path(path)
    if not path.exists():return ''
    return path.read_bytes()[-n:].decode(errors='replace')

def update(status,**extra):
    old=read(STATE) if STATE.exists() else {}
    old.update(extra,status=status,updated_unix=time.time())
    if status in ('PREFLIGHTING','PREPARING','RUNNING_SEARCH','REPLAYING'):
        old['pid']=os.getpid()
    atomic(STATE,old)
    print('CURVE302_SEED_PANEL',status,json.dumps(extra,sort_keys=True),flush=True)


def panel_order():
    if not CLOSURE.exists():raise RuntimeError('unlock-seed closure v2 artifact is missing')
    data=read(CLOSURE)
    ids=tuple(data['direction_ids'])
    if set(ids)!=set(EXPECTED_DIRECTIONS) or len(ids)!=14:
        raise RuntimeError('frozen direction roster changed')
    ranking=data['single_seed_ranking']
    ranked=[row['seed_direction'] for row in ranking]
    if len(ranked)!=14 or set(ranked)!=set(EXPECTED_DIRECTIONS):
        raise RuntimeError('single-seed ranking is not a complete permutation')
    return tuple(seed for seed in ranked if seed not in KNOWN_WINNERS)


def load_base():
    path=CAS/'run_curve302_seeded_v3_amplifier.py'
    spec=importlib.util.spec_from_file_location('curve302_seed_panel_base',path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    # Rebind all evidence endpoints before calling any base operation.
    module.D=D; module.STATE=STATE; module.LOG=LOG; module.SELF=SELF
    module.SEEDS=panel_order()
    return module


def preflight_all(*,publish=True):
    base=load_base(); seeds=panel_order(); rows=[]
    if publish:update('PREFLIGHTING',case_count=len(seeds),charts=0)
    for i,seed in enumerate(seeds):
        base.construct_seed(seed)
        ctx=base.context(seed)
        if ctx.state.rank!=18:
            raise RuntimeError(f'{seed}: preflight rank {ctx.state.rank} != 18')
        rows.append({'seed_direction':seed,'rank':18,'index':i})
        print(f'CURVE302_SEED_PANEL_PREFLIGHT|seed={seed}|rank=18|index={i+1}/{len(seeds)}',flush=True)
    result={'status':'PASS_ZERO_CHART_PREFLIGHT','case_count':len(seeds),'charts':0,'order':list(seeds),'rows':rows}
    atomic(D/'preflight.json',result)
    if publish:update('PREFLIGHT_PASS',case_count=len(seeds),charts=0,order=list(seeds))
    print(f'CURVE302_SEED_PANEL_PREFLIGHT_COMPLETE|cases={len(seeds)}|charts=0|status=PASS',flush=True)
    return result


def worker():
    import det1092_v3_worker as searcher
    base=load_base(); seeds=panel_order()
    preflight_all(publish=False)
    results=[]
    for i,seed in enumerate(seeds):
        folder=D/seed
        verified=folder/'seeded-verified.json'
        if verified.exists():
            result=read(verified); results.append(result)
            print('CURVE302_SEED_PANEL_REUSE',seed,result['rank_lower_bound'],flush=True)
            continue
        ctx=base.context(seed)
        terminal=folder/'replay-M17/terminal.json'
        if not terminal.exists():
            update('RUNNING_SEARCH',case=seed,case_index=i+1,case_count=len(seeds),initial_rank=18)
            searcher.run_search(ctx)
        update('REPLAYING',case=seed,case_index=i+1,case_count=len(seeds))
        result=base.verify_case(ctx); results.append(result)
        update('CASE_VERIFIED',case=seed,case_index=i+1,case_count=len(seeds),
               rank_lower_bound=result['rank_lower_bound'],charts=result['charts'],stop_reason=result['stop_reason'])
    by_seed={row['seed_direction']:row for row in results}
    summary={
        'status':'COMPLETE_SEED_UNIVERSALITY_PANEL',
        'order':list(seeds),
        'known_winners_excluded':list(KNOWN_WINNERS),
        'results':results,
        'rank31_count':sum(row['rank_lower_bound']>=31 for row in results),
        'any_no_gain':any(row['rank_lower_bound']==18 for row in results),
        'minimum_terminal_rank':min(row['rank_lower_bound'] for row in results),
        'maximum_terminal_rank':max(row['rank_lower_bound'] for row in results),
        'total_charts':sum(row['charts'] for row in results),
        'all_cases_present':set(by_seed)==set(seeds),
        'claim_boundary':'Retrospective known-seed amplification panel; not a prospective selector or new rank claim.',
    }
    atomic(D/'summary.json',summary)
    update('COMPLETE_SEED_UNIVERSALITY_PANEL',rank31_count=summary['rank31_count'],
           total_charts=summary['total_charts'],minimum_terminal_rank=summary['minimum_terminal_rank'],
           maximum_terminal_rank=summary['maximum_terminal_rank'])


def launch(resume=False):
    D.mkdir(parents=True,exist_ok=True)
    old=read(STATE) if STATE.exists() else {}
    if process_alive(old.get('pid')):raise RuntimeError('seed universality worker is already alive')
    if old and not resume and old.get('status') not in (None,'NOT_LAUNCHED','PREFLIGHT_PASS'):
        raise RuntimeError('existing panel state; use resume')
    if not (D/'preflight.json').exists():
        raise RuntimeError('run foreground preflight before launching the long panel')
    sage=sage_launcher()
    with LOG.open('ab',buffering=0) as log:
        proc=subprocess.Popen([sage,'-python','-u',str(SELF),'worker'],cwd=ROOT,stdin=subprocess.DEVNULL,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
            env={**os.environ,'PYTHONUNBUFFERED':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    atomic(STATE,{'status':'LAUNCHED','pid':proc.pid,'sage':sage,'updated_unix':time.time(),'order':list(panel_order())})
    print('Launched curve302 seed universality panel pid='+str(proc.pid))

def status():
    row=read(STATE) if STATE.exists() else {'status':'NOT_LAUNCHED'}
    row['process_alive']=process_alive(row.get('pid'))
    cases={}
    for seed in panel_order():
        folder=D/seed
        if (folder/'seeded-verified.json').exists():
            v=read(folder/'seeded-verified.json'); cases[seed]={'status':'VERIFIED','rank':v['rank_lower_bound'],'charts':v['charts']}
        elif (folder/'replay-M17/terminal.json').exists():cases[seed]={'status':'SEARCH_SEALED_REPLAY_PENDING'}
        elif (folder/'protocol.json').exists():cases[seed]={'status':'SEED_PREPARED'}
        else:cases[seed]={'status':'PENDING'}
    row['cases']=cases; row['progress_tail']=tail(LOG).splitlines()[-12:]
    if row.get('status') in ('LAUNCHED','RUNNING_SEARCH','REPLAYING','PREFLIGHTING') and not row['process_alive']:
        row['effective_status']='STOPPED_REVIEW_REQUIRED'
    print(json.dumps(row,indent=2,sort_keys=True))

def diagnose():status();print(tail(LOG,40000))
def stop():
    row=read(STATE) if STATE.exists() else {}; pid=row.get('pid')
    if not process_alive(pid):print('No live seed universality worker');return
    os.killpg(pid,signal.SIGTERM);print('Stop requested; checkpoints retained')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('preflight','launch','resume','status','diagnose','stop','worker'))
    a=p.parse_args()
    if a.action=='preflight':preflight_all()
    elif a.action=='launch':launch(False)
    elif a.action=='resume':launch(True)
    elif a.action=='status':status()
    elif a.action=='diagnose':diagnose()
    elif a.action=='stop':stop()
    else:worker()

if __name__=='__main__':main()
