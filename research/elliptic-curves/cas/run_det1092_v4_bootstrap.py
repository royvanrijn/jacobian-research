#!/usr/bin/env python3
"""Detached controller for the eight-case determinant-1092 V4 bootstrap run."""
from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid

import det1092_v4_bootstrap_contract as c
from v3_warm_support import atomic, bindings, process_info, read, require, same_process, sha, tail

SELF=Path(__file__).resolve(); ROOT=c.ROOT; AUTO=c.AUTO
STATE=AUTO/'state.json'; LOG=AUTO/'autorun.log'; LOCK=AUTO/'controller.lock'; WORKER=c.CAS/'det1092_v4_bootstrap_worker.py'


def sage_launcher():
    candidates=[os.environ.get('V4_SAGE'),shutil.which('sage'),str(Path.home()/'.local/bin/sage'),'/usr/bin/sage']
    for raw in candidates:
        if raw:
            p=Path(raw).expanduser()
            if p.is_file() and os.access(p,os.X_OK): return str(p.resolve())
    raise RuntimeError('Sage launcher not found; set V4_SAGE=/absolute/path/to/sage')


def controller_alive(state):
    return bool(state.get('start_token') and same_process(state.get('pid'),state.get('start_token')))


def take_lock():
    AUTO.mkdir(parents=True,exist_ok=True); fd=os.open(LOCK,os.O_CREAT|os.O_RDWR,0o600)
    try: fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB); return fd
    except BaseException:
        os.close(fd); raise RuntimeError('another V4 controller owns the lock')


def validate_result(case):
    path=c.D/case/'v4-verified.json'
    if not path.exists(): return None
    row=read(path)
    require(row.get('status')=='PASS_INDEPENDENT_DET1092_V4_REPLAY' and row.get('case')==case,
            'unrecognized V4 result for '+case)
    bindings(ROOT,row['bindings']); bindings(ROOT,row['sources']); return row


class Reporter:
    def __init__(self,session):
        info=process_info(os.getpid()); self.value={'pid':os.getpid(),'start_token':info['start_token'],
            'session':str(session.relative_to(ROOT))}; self.lock=threading.RLock(); self.done=threading.Event()
        self.thread=threading.Thread(target=self._beat,daemon=True)
    def update(self,status=None,**values):
        with self.lock:
            self.value.update(values)
            if status is not None:self.value['status']=status
            self.value['updated_unix']=time.time(); atomic(STATE,self.value)
            if status is not None: print('V4',status,json.dumps(values,sort_keys=True),flush=True)
    def _beat(self):
        while not self.done.wait(10):
            try:self.update()
            except Exception as exc:print('V4_HEARTBEAT_FAILED',repr(exc),flush=True)
    def __enter__(self):self.update('CONTROLLER_RUNNING');self.thread.start();return self
    def __exit__(self,*_):self.done.set();self.thread.join(timeout=2)


class ResourceStop(RuntimeError):pass


def run_step(action,case,session_path,session,deadline,report):
    from research_runtime.supervisor import Limits,run
    remaining=deadline-time.monotonic()
    if remaining<=0:raise ResourceStop('campaign wall budget exhausted')
    maximum=c.RESOURCE[action]
    attempt=session_path.parent/'attempts'/(case or '_prepare')/(action+'-'+uuid.uuid4().hex[:10]);attempt.mkdir(parents=True)
    command=[session['sage'],'-python','-u',str(WORKER),action]
    if case:command += ['--case',case]
    label={'prepare':'PREPARING','preflight':'PREFLIGHT','search':'SEARCHING','replay':'REPLAYING'}[action]
    report.update(label,case=case,worker_action=action,attempt=str(attempt.relative_to(ROOT)),
                  worker_log=str((attempt/'worker.log').relative_to(ROOT)))
    result=run(command,limits=Limits(min(maximum,remaining),c.RESOURCE['rss_bytes']),cwd=ROOT,
               env={**os.environ,'DET1092_V4_SUPERVISED':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},
               log_path=attempt/'worker.log',checkpoint_path=attempt/'supervisor.json')
    if result['outcome']!='completed':
        if result['outcome'] in ('strict_wall_timeout','strict_rss_limit'):
            raise ResourceStop((case or 'prepare')+' '+action+': '+result['outcome']+'; checkpoints retained')
        raise RuntimeError((case or 'prepare')+' '+action+': '+result['outcome']+'\n'+(tail(attempt/'worker.log') or ''))


def worker(session_path,lease_fd):
    require(os.fstat(lease_fd).st_ino==LOCK.stat().st_ino,'wrong V4 lock lease')
    signal.signal(signal.SIGHUP,signal.SIG_IGN)
    def interrupted(signum,frame):raise KeyboardInterrupt('V4 controller signal '+str(signum))
    signal.signal(signal.SIGTERM,interrupted)
    session=read(session_path); deadline=time.monotonic()+session['hours']*3600
    try:
        with Reporter(session_path) as report:
            try:
                require(session['sources']==c.own_sources(),'V4 sources changed after launch');bindings(ROOT,session['sources'])
                if not (c.D/'roster.json').exists():run_step('prepare',None,session_path,session,deadline,report)
                roster=c.validate_roster();results=[]
                for row in roster['cases']:
                    case=row['id'];old=validate_result(case)
                    if old is not None:results.append(old);continue
                    run_step('preflight',case,session_path,session,deadline,report)
                    if not (c.D/case/'bootstrap/terminal.json').exists():run_step('search',case,session_path,session,deadline,report)
                    run_step('replay',case,session_path,session,deadline,report)
                    result=validate_result(case);require(result is not None,'replay exited without V4 result')
                    results.append(result);report.update('CASE_VERIFIED',case=case,rank_lower_bound=result['rank_lower_bound'],
                                                         gain=result['gain'],charts=result['bootstrap_charts'])
                summary=c.result_summary(results);atomic(c.D/'summary.json',summary)
                report.update('COMPLETE_V4_ROSTER',results=results,summary=summary)
            except ResourceStop as exc:
                report.update('STOPPED_RESOURCE_BUDGET',error=str(exc))
            except BaseException as exc:
                report.update('STOPPED_REVIEW_REQUIRED',error=repr(exc));raise
    finally:os.close(lease_fd)


def launch(hours,resume=False):
    require(math.isfinite(hours) and 0<hours<=72,'--hours must be in (0,72]')
    fd=take_lock()
    try:
        old=read(STATE) if STATE.exists() else {}
        require(not controller_alive(old),'V4 controller is still alive')
        if old.get('status')=='COMPLETE_V4_ROSTER':print(json.dumps(old,indent=2));return
        require(resume or not old,'existing V4 state: use resume')
        AUTO.mkdir(parents=True,exist_ok=True); sage=sage_launcher();session_dir=AUTO/'sessions'/str(uuid.uuid4());session_dir.mkdir(parents=True)
        session_path=session_dir/'session.json';session={'schema':'det1092-v4-controller.v1','sources':c.own_sources(),
            'sage':sage,'hours':hours,'created_unix':time.time(),'action':'resume' if resume else 'launch','scope':c.CLAIM}
        atomic(session_path,session,immutable=True)
        if old:atomic(session_dir/'previous-state.json',old,immutable=True)
        with LOG.open('ab',buffering=0) as log:
            proc=subprocess.Popen([sys.executable,str(SELF),'worker','--session',str(session_path),'--lease-fd',str(fd)],
                cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,pass_fds=(fd,),
                env={**os.environ,'PYTHONUNBUFFERED':'1'})
        print('Launched V4 controller pid='+str(proc.pid));print('Status: python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py status')
    finally:os.close(fd)


def status():
    row=read(STATE) if STATE.exists() else {'status':'NOT_LAUNCHED'};row=dict(row);row['process_alive']=controller_alive(row)
    info=process_info(row.get('pid'));row['process_state']=info['state'] if info else None
    if row.get('worker_log'):
        recent=tail(ROOT/row['worker_log'],8192);row['progress_tail']=recent.splitlines()[-6:] if recent else []
    if row.get('results'):
        row['results']=[{k:r.get(k) for k in ('case','rank_lower_bound','gain','bootstrap_charts','stop_reason')} for r in row['results']]
    print(json.dumps(row,indent=2,sort_keys=True))
    if (c.D/'roster.json').exists():
        for item in read(c.D/'roster.json')['cases']:
            case=item['id'];p=c.D/case/'v4-verified.json';t=c.D/case/'bootstrap/terminal.json'
            if p.exists():
                r=read(p);print(case,'VERIFIED',json.dumps({k:r.get(k) for k in ('rank_lower_bound','gain','bootstrap_charts','stop_reason')},sort_keys=True))
            elif t.exists():print(case,'SEARCH_SEALED_REPLAY_PENDING',json.dumps(read(t),sort_keys=True))
            else:
                b=c.D/case/'bootstrap'; charts=len(list(b.glob('chart-*.json'))) if b.exists() else 0
                print(case,'PENDING','charts',charts)


def diagnose():
    row=read(STATE) if STATE.exists() else {};path=ROOT/row['worker_log'] if row.get('worker_log') else LOG
    print(tail(path) or 'NO_LOG')


def stop():
    row=read(STATE) if STATE.exists() else {};require(controller_alive(row),'no live owned V4 controller')
    pid=int(row['pid'])
    if hasattr(os,'pidfd_open'):
        fd=os.pidfd_open(pid)
        try:signal.pidfd_send_signal(fd,signal.SIGTERM)
        finally:os.close(fd)
    else:os.kill(pid,signal.SIGTERM)
    print('Sent SIGTERM to verified V4 controller',pid)


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['launch','resume','status','diagnose','stop','worker']);p.add_argument('--hours',type=float,default=24);p.add_argument('--session',type=Path);p.add_argument('--lease-fd',type=int);a=p.parse_args()
    if a.action=='launch':launch(a.hours,False)
    elif a.action=='resume':launch(a.hours,True)
    elif a.action=='status':status()
    elif a.action=='diagnose':diagnose()
    elif a.action=='stop':stop()
    else:worker(a.session,a.lease_fd)

if __name__=='__main__':main()
