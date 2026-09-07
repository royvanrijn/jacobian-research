#!/usr/bin/env python3
"""Detached warm-start V3 transfer over the three frozen rank-27 11952 fibres.

This is a DIFFERENT hypothesis from generic-start transfer. The completed
native11952 control stayed 17->17 and remains immutable. This runner reuses
only the already-frozen warm seed/orbit/protocol inputs from v3-transfer-11952-v3
and asks whether V3 can exploit an already enlarged rank-27 subgroup.

No positive-control claim is made. Each warm case is run and independently
replayed under the existing supervisor. Operational failures are preserved,
surfaced by status/diagnose, and may be explicitly resumed after review.
"""
from __future__ import annotations
from importlib.machinery import SourceFileLoader
from pathlib import Path
import argparse, json, os, subprocess, sys, time, hashlib, shutil

SELF=Path(__file__).resolve(); CAS=SELF.parent
v3=SourceFileLoader('v3_warm_runtime',str(CAS/'v3_transfer_campaign_v3.sage')).load_module()
base=v3.base
ROOT=base.ROOT
D=ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v3'
AUTO=ROOT/'artifacts/local/elliptic-curves/v3-warm-start-overnight-v1'
STATE=AUTO/'state.json'; LOG=AUTO/'autorun.log'
WARM=('warm-11952-41','warm-11952-72','warm-11952-186')


def read(p): return json.loads(Path(p).read_text())
def atomic(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');tmp.replace(path)
def sha(p):
    h=hashlib.sha256();h.update(Path(p).read_bytes());return h.hexdigest()
def alive(pid):
    try:
        if not pid:return False
        os.kill(int(pid),0);return True
    except (ProcessLookupError,PermissionError,ValueError):return False

def tail(path,n=120):
    p=Path(path)
    if not p.exists():return None
    return '\n'.join(p.read_text(errors='replace').splitlines()[-n:])

def validate_authorization():
    control=D/'control-native11952/verified.json'
    if not control.exists(): raise RuntimeError('completed independently replayed native11952 control is required')
    c=read(control)
    if c.get('status')!='PASS_INDEPENDENT_TRANSFER_REPLAY' or c.get('initial_rank')!=17 or c.get('rank_lower_bound')!=17 or c.get('gain')!=0:
        raise RuntimeError('expected the preserved clean 17->17 generic control')
    for name,digest in c.get('bindings',{}).items():
        p=ROOT/name
        if not p.exists() or sha(p)!=digest: raise RuntimeError('control binding changed: '+name)
    roster=read(D/'roster.json'); ids={r['id']:r for r in roster['cases']}
    for case in WARM:
        if case not in ids: raise RuntimeError('missing frozen warm case '+case)
        folder=D/case
        if ids[case]['initial_rank']!=27: raise RuntimeError('warm seed is not frozen at rank 27: '+case)
        if sha(folder/'protocol.json')!=ids[case]['protocol_sha256']: raise RuntimeError('warm protocol changed: '+case)
        p=read(folder/'protocol.json')
        for name,digest in p['inputs'].items():
            q=ROOT/name
            if not q.exists() or sha(q)!=digest: raise RuntimeError('warm input changed: '+name)
    return {'schema':'v3-warm-start-transfer-authorization.v1','control':'17->17 finite no-gain',
            'hypothesis':'Can unchanged V3 exploit already enlarged independently certified rank-27 subgroups on parent 11952?',
            'cases':list(WARM),'runner_sha256':sha(SELF),'created_unix':time.time(),
            'claim_boundary':'Warm-start transfer only. The generic-start positive-control gate failed and is not overridden or relabeled as a pass.'}


def warm_bind(case):
    if case not in WARM: raise RuntimeError('warm runner refuses non-warm case')
    folder=D/case; roster=read(D/'roster.json'); entry=next(r for r in roster['cases'] if r['id']==case)
    if sha(D/'parent-bank.json')!=roster['parent_bank_sha256'] or sha(D/'gate.json')!=roster['gate_sha256'] or sha(D/'preparation.json')!=roster['preparation_sha256']:
        raise RuntimeError('frozen transfer preparation changed')
    if sha(folder/'protocol.json')!=entry['protocol_sha256']: raise RuntimeError('job protocol changed')
    engine=base.engine_for(folder); policy=read(folder/'protocol.json')
    base.check_bindings(ROOT,policy['inputs'])
    return engine,folder,policy


def bind_with_cert(case):
    engine,folder,policy=warm_bind(case)
    seed=read(folder/'seed-input.json'); proof=read(folder/'seed-proof.json')
    frozen_curve=v3._normal_curve(seed['curve']); frozen_points=v3._normal_points(seed['points'])
    ec,ep=engine.v1.seed(None); ec=v3._normal_curve(ec); ep=v3._normal_points(ep)
    if ec!=frozen_curve or v3._first_diff(ep,frozen_points) is not None: raise RuntimeError('engine/frozen warm seed mismatch')
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    from research_runtime.arithmetic import ArithmeticContext,CurveModel
    from research_runtime.mw_state import MWState
    primes=tuple(int(r['prime']) for r in proof['signatures']); torsion=int(proof['no_rational_2_torsion_prime'])
    cache=Cache(MemoryFactStore()); state=MWState.empty(ArithmeticContext.for_search(CurveModel(frozen_curve)),cache=cache,primes=primes,no_two_torsion_prime=torsion)
    for i,p in enumerate(frozen_points):
        before=state.rank; state=state.adjoin(p,cache=cache,extra_primes=())
        if state.rank!=before+1: raise RuntimeError(f'warm seed certificate prefix failed at {i}')
    if state.rank!=27 or tuple(state.basis)!=frozen_points: raise RuntimeError('warm marked rank-27 seed did not replay exactly')
    v3._active_seed=(ec,ep,state)
    return engine,folder,policy

base.bind_job=bind_with_cert
base.run_case.__globals__['bind_job']=bind_with_cert
base.run_case.__globals__['_initial_worker_state']=v3._initial_worker_state
base.replay_case.__globals__['bind_job']=bind_with_cert


def attempt_paths(case,action):
    stem='warm-'+action.replace('-worker','');folder=D/case
    return folder/(stem+'.supervisor.json'),folder/(stem+'.log')

def latest_failure():
    found=[]
    for c in WARM:
        for a in ('run-worker','replay-worker'):
            sup,log=attempt_paths(c,a)
            if sup.exists():
                r=read(sup)
                if r.get('outcome')!='completed':
                    found.append((sup.stat().st_mtime,c,a,r,log))
    if not found:return None
    _,c,a,r,log=max(found,key=lambda x:x[0])
    return {'case':c,'action':a,'supervisor':r,'log':str(log.relative_to(ROOT)),'log_tail':tail(log)}


def supervise(action,case):
    from research_runtime.supervisor import Limits,run
    sup,log=attempt_paths(case,action)
    if sup.exists():
        report=read(sup)
        if report.get('outcome')=='completed': return
        raise RuntimeError(f'preserved failed warm attempt exists for {case}: {sup}')
    seconds=base.LIMITS['case_wall_seconds'] if action=='run-worker' else base.LIMITS['replay_wall_seconds']
    command=[sys.executable,str(SELF),action,'--case',case]
    env={**os.environ,'V3_WARM_SUPERVISED':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    report=run(command,limits=Limits(seconds,base.LIMITS['rss_bytes']),cwd=ROOT,env=env,log_path=log,checkpoint_path=sup)
    if report['outcome']!='completed': raise RuntimeError(f'{case} {action} stopped: {report["outcome"]}')


def worker():
    try:
        auth=validate_authorization();AUTO.mkdir(parents=True,exist_ok=True)
        if not (AUTO/'authorization.json').exists(): atomic(AUTO/'authorization.json',auth)
        for case in WARM:
            folder=D/case
            if (folder/'warm-verified.json').exists(): continue
            atomic(STATE,{'status':'RUNNING_WARM_CASE','case':case,'pid':os.getpid(),'updated_unix':time.time()})
            if not (folder/'replay-M17/terminal.json').exists(): supervise('run-worker',case)
            supervise('replay-worker',case)
            verified=read(folder/'verified.json')
            out={'schema':'v3-warm-start-transfer-result.v1','case':case,'initial_rank':27,
                 'rank_lower_bound':verified['rank_lower_bound'],'gain':verified['rank_lower_bound']-27,
                 'stop_reason':verified['stop_reason'],'charts':verified['charts'],
                 'source_verified_sha256':sha(folder/'verified.json'),
                 'claim_boundary':'Warm-start transfer result; generic 17->17 control remained a clean no-gain.'}
            atomic(folder/'warm-verified.json',out)
        atomic(STATE,{'status':'COMPLETE_WARM_ROSTER','pid':os.getpid(),'updated_unix':time.time(),
                      'results':[read(D/c/'warm-verified.json') for c in WARM]})
    except Exception as exc:
        f=latest_failure()
        old=read(STATE) if STATE.exists() else {}
        atomic(STATE,{'status':'STOPPED_REVIEW_REQUIRED','case':old.get('case'),'pid':os.getpid(),
                      'error':repr(exc),'failure':f,'updated_unix':time.time()})
        raise


def launch():
    auth=validate_authorization();AUTO.mkdir(parents=True,exist_ok=True)
    if STATE.exists():
        s=read(STATE)
        if s.get('status')=='COMPLETE_WARM_ROSTER': print('Already complete');return
        if alive(s.get('pid')): raise RuntimeError(f'warm controller pid {s.get("pid")} is still alive')
        if s.get('status')=='STOPPED_REVIEW_REQUIRED':
            raise RuntimeError('preserved failed attempt exists; use diagnose, then resume after the underlying issue is understood/fixed')
    stream=LOG.open('ab',buffering=0)
    p=subprocess.Popen([sys.executable,str(SELF),'worker'],cwd=ROOT,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
    atomic(STATE,{'status':'LAUNCHED','pid':p.pid,'updated_unix':time.time(),'authorization':auth})
    print('Launched detached warm-start V3 roster pid='+str(p.pid))

def archive_failed(case,action):
    sup,log=attempt_paths(case,action)
    if not sup.exists():raise RuntimeError('no failed supervisor to archive')
    report=read(sup)
    if report.get('outcome')=='completed':raise RuntimeError('attempt completed; refusing archive-as-failure')
    stamp=time.strftime('%Y%m%dT%H%M%SZ',time.gmtime())+'-'+sha(sup)[:12]
    dst=D/case/'warm-failed-attempts'/stamp;dst.mkdir(parents=True,exist_ok=False)
    for p in (sup,log):
        if p.exists():shutil.move(str(p),str(dst/p.name))
    # Preserve an unsealed partial epoch, but never move sealed stages.
    replay=D/case/'replay-M17'
    if replay.exists():
        for ep in sorted(replay.glob('epoch-*')):
            if not (ep/'stage.json').exists() and any(ep.iterdir()):
                shutil.move(str(ep),str(dst/ep.name))
    return dst

def resume():
    if not STATE.exists():raise RuntimeError('no prior warm state')
    s=read(STATE)
    if alive(s.get('pid')):raise RuntimeError('controller still alive')
    f=latest_failure()
    if not f:raise RuntimeError('no preserved failed warm attempt found')
    dst=archive_failed(f['case'],f['action'])
    atomic(STATE,{'status':'RESUMING','archived':str(dst.relative_to(ROOT)),'updated_unix':time.time()})
    stream=LOG.open('ab',buffering=0)
    p=subprocess.Popen([sys.executable,str(SELF),'worker'],cwd=ROOT,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
    atomic(STATE,{'status':'RELAUNCHED','pid':p.pid,'updated_unix':time.time(),'archived':str(dst.relative_to(ROOT))})
    print('Relaunched detached warm-start roster pid='+str(p.pid))

def status():
    s=read(STATE) if STATE.exists() else None
    if s is None: print('NOT_LAUNCHED')
    else:
        x=dict(s);x['process_alive']=alive(x.get('pid'))
        if x.get('status') in ('RUNNING_WARM_CASE','LAUNCHED','RELAUNCHED') and not x['process_alive']:
            x['effective_status']='STOPPED_REVIEW_REQUIRED'
            x['failure']=latest_failure()
        print(json.dumps(x,indent=2,sort_keys=True))
    for c in WARM:
        p=D/c/'warm-verified.json'
        if p.exists(): print(c,json.dumps(read(p),sort_keys=True))
        else:
            stages=D/c/'replay-M17/stages.json'
            if stages.exists(): print(c,'stages',json.dumps(read(stages)[-1],sort_keys=True))
            else: print(c,'NOT_YET_COMPLETE')

def diagnose():
    f=latest_failure()
    print('NO_FAILED_WARM_ATTEMPT' if f is None else json.dumps(f,indent=2,sort_keys=True))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['launch','status','diagnose','resume','worker','run-worker','replay-worker']);ap.add_argument('--case',choices=WARM);a=ap.parse_args()
    if a.action=='launch': launch()
    elif a.action=='status': status()
    elif a.action=='diagnose': diagnose()
    elif a.action=='resume': resume()
    elif a.action=='worker': worker()
    else:
        if os.environ.get('V3_WARM_SUPERVISED')!='1': raise RuntimeError('worker action must be supervised')
        if a.action=='run-worker': base.run_case(a.case)
        else: base.replay_case(a.case)

if __name__=='__main__': main()
