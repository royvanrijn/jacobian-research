#!/usr/bin/env python3
"""Detached eight-fibre same-parent V3 pilot. Standard-library controller only.

launch/resume --hours 10 | status | diagnose | stop
The forty reserves are frozen but NEVER automatically searched by this script.
Each arithmetic service is a supervised Sage process with immutable inputs.
"""
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

import det1092_v3_contract as c
from v3_warm_support import (atomic, bindings, process_info, read, require,
                            same_process, sha, tail, terminal_structure)

SELF = Path(__file__).resolve()
STATE, LOG, LOCK = c.AUTO/'state.json', c.AUTO/'autorun.log', c.AUTO/'controller.lock'
ENGINE = c.CAS/'det1092_v3_worker.py'
ACTIVE = {'LAUNCHING','CONTROLLER_RUNNING','PREPARING','PREFLIGHT','SEARCHING','REPLAYING',
          'SEALED_SEARCH_REUSED','CASE_VERIFIED','CASE_RESOURCE_STOP'}


def sage_launcher():
    explicit = os.environ.get('V3_SAGE')
    candidates = [explicit] if explicit else [shutil.which('sage'), str(Path.home()/'.local/bin/sage'), '/usr/bin/sage']
    for value in candidates:
        if value:
            path = Path(value).expanduser()
            if path.is_file() and os.access(path,os.X_OK): return str(path.resolve())
    raise RuntimeError('Sage launcher not found; set V3_SAGE=/absolute/path/to/sage')


def worker_command(sage, action, case, session):
    require(action in ('prepare','preflight','search','replay'), 'unknown arithmetic action')
    require(action == 'prepare' or case is not None, 'case is required')
    argv = [str(sage),'-python','-u',str(ENGINE),action,'--session',str(session)]
    if case is not None: argv += ['--case',case]
    return argv


def alive(row):
    return same_process(row.get('pid'),row.get('start_token'))


def ensure_no_live_workers():
    for path in (c.AUTO/'sessions').glob('*/attempts/*/*/supervisor.json'):
        row = read(path)
        require(not same_process(row.get('pid'),row.get('start_token')), 'arithmetic worker still alive: '+str(path))


def take_lock():
    c.AUTO.mkdir(parents=True,exist_ok=True)
    fd = os.open(LOCK,os.O_CREAT|os.O_RDWR,0o600)
    try:
        fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB); return fd
    except BaseException:
        os.close(fd); raise RuntimeError('another pilot controller owns the lifetime lock')


def launch(hours, resume=False):
    require(math.isfinite(hours) and 0 < hours <= 24, 'hours must be positive and at most 24')
    fd = take_lock()
    try:
        old = read(STATE) if STATE.exists() else {}
        require(not alive(old), 'pilot controller is still alive')
        ensure_no_live_workers()
        require(resume or not old, 'existing pilot: use resume; no artifacts will be deleted')
        if old.get('status') == 'COMPLETE_EIGHT_PILOT':
            print(json.dumps(old,indent=2)); return
        # Early validation is outcome blind and fast. Rank checks stay in Sage.
        selected, path = c.selection_input(); c.choose(selected)
        for p, expected in ((c.REDUCED_PARENT,c.PARENT_BLOB),(c.ORBITS,c.ORBITS_BLOB),(c.LATTICE,c.LATTICE_BLOB)):
            require(c.git_blob(p) == expected, 'generic input differs: '+str(p))
        session_dir = c.AUTO/'sessions'/str(uuid.uuid4()); session_dir.mkdir(parents=True)
        session_path = session_dir/'session.json'
        inputs = {str(p.relative_to(c.ROOT)):sha(p) for p in
                  (path,c.REDUCED_PARENT,c.ORBITS,c.LATTICE,c.CALIBRATION/'protocol.json')}
        session = {'schema':'det1092-v3-session.v1','sources':c.own_sources(),'inputs':inputs,
                   'hours':hours,'sage':sage_launcher(),'created_unix':time.time(),
                   'action':'resume' if resume else 'launch','scope':c.CLAIM}
        atomic(session_path,session,immutable=True)
        if old: atomic(session_dir/'previous-controller-state.json',old,immutable=True)
        atomic(STATE,{'status':'LAUNCHING','session':str(session_path.relative_to(c.ROOT)),'updated_unix':time.time()})
        with LOG.open('ab',buffering=0) as log:
            proc = subprocess.Popen([sys.executable,str(SELF),'worker','--session',str(session_path),
                 '--lease-fd',str(fd)],cwd=c.ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                 start_new_session=True,pass_fds=(fd,),env={**os.environ,'PYTHONUNBUFFERED':'1'})
        print('Launched same-parent pilot pid='+str(proc.pid)+'; session='+str(session_path.relative_to(c.ROOT)))
        print('Status: python3 elliptic-curves/cas/run_det1092_v3_trial.py status')
    finally:
        os.close(fd)  # no LOCK_UN: child retains the inherited lease


class Reporter:
    def __init__(self,session):
        info = process_info(os.getpid())
        self.value = {'pid':os.getpid(),'start_token':info['start_token'],'session':str(session.relative_to(c.ROOT))}
        self.lock = threading.RLock(); self.done = threading.Event()
        self.thread = threading.Thread(target=self.heartbeat,daemon=True)

    def update(self,status=None,**values):
        with self.lock:
            self.value.update(values)
            if status is not None: self.value['status'] = status
            self.value['updated_unix'] = time.time(); atomic(STATE,self.value)
            if status is not None: print('DET1092',status,json.dumps(values,sort_keys=True),flush=True)

    def heartbeat(self):
        while not self.done.wait(10):
            try: self.update()
            except Exception as exc: print('HEARTBEAT_WRITE_FAILED',repr(exc),flush=True)

    def __enter__(self):
        self.update('CONTROLLER_RUNNING'); self.thread.start(); return self

    def __exit__(self,*_):
        self.done.set(); self.thread.join(timeout=2)


class ResourceStop(RuntimeError):
    pass


def verified(case):
    path = c.D/case/'trial-verified.json'
    if not path.exists(): return None
    r = read(path)
    require(r.get('status') == 'PASS_INDEPENDENT_DET1092_REPLAY' and r.get('case') == case and
            r.get('initial_rank') == 17 and r.get('gain') == r['rank_lower_bound']-17, 'invalid independent trial result')
    bindings(c.ROOT,r['sources']); bindings(c.ROOT,r['bindings'])
    return r


def run_step(action,case,session_path,session,deadline,report):
    from research_runtime.supervisor import Limits,run
    remaining = deadline-time.monotonic()
    if remaining <= 0: raise ResourceStop('session wall budget exhausted')
    attempt = session_path.parent/'attempts'/(case or 'preparation')/(action+'-'+uuid.uuid4().hex[:12])
    attempt.mkdir(parents=True)
    command = worker_command(session['sage'],action,case,session_path)
    report.update({'prepare':'PREPARING','preflight':'PREFLIGHT','search':'SEARCHING','replay':'REPLAYING'}[action],
        case=case,worker_action=action,attempt=str(attempt.relative_to(c.ROOT)),
        worker_log=str((attempt/'worker.log').relative_to(c.ROOT)))
    result = run(command,limits=Limits(min(c.RESOURCE[action],remaining),c.RESOURCE['rss_bytes']),
        cwd=c.ROOT,env={**os.environ,'DET1092_V3_SUPERVISED':'1','OPENBLAS_NUM_THREADS':'1',
                       'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},
        log_path=attempt/'worker.log',checkpoint_path=attempt/'supervisor.json')
    if result['outcome'] != 'completed':
        # A committed terminal can be INPUT to verification despite an abnormal
        # process exit. It is never automatically marked verified/completed.
        if action == 'search' and terminal_structure(c.D/case) is not None:
            report.update('SEALED_SEARCH_REUSED',case=case,original_outcome=result['outcome']); return
        if result['outcome'] in ('strict_wall_timeout','strict_rss_limit'):
            raise ResourceStop(f'{case or "preparation"} {action}: {result["outcome"]}')
        raise RuntimeError(f'{case} {action}: {result["outcome"]}\n'+(tail(attempt/'worker.log') or ''))
    if action == 'prepare': c.validate_roster()
    elif action == 'search': require(terminal_structure(c.D/case) is not None,'search exited without sealed terminal')
    elif action == 'replay': require(verified(case) is not None,'replay exited without independent certificate')


def worker(session_path,lease_fd):
    require(session_path is not None and lease_fd is not None,'launch through the controller')
    st = os.fstat(lease_fd); lock = LOCK.stat()
    require((st.st_dev,st.st_ino) == (lock.st_dev,lock.st_ino),'invalid inherited lifetime lease')
    signal.signal(signal.SIGHUP,signal.SIG_IGN)
    def interrupt(signum,frame): raise KeyboardInterrupt(f'controller signal {signum}')
    signal.signal(signal.SIGTERM,interrupt)
    session = read(session_path); deadline = time.monotonic()+session['hours']*3600
    try:
        with Reporter(session_path) as report:
            try:
                require(session['sources'] == c.own_sources(),'implementation changed after launch')
                bindings(c.ROOT,session['sources']); bindings(c.ROOT,session['inputs'])
                if not (c.D/'roster.json').exists():
                    run_step('prepare',None,session_path,session,deadline,report)
                roster = c.validate_roster()
                results, stopped = [], []
                for row in roster['pilot']:
                    case = row['id']; old = verified(case)
                    if old is not None: results.append(old); continue
                    try:
                        run_step('preflight',case,session_path,session,deadline,report)
                        terminal = terminal_structure(c.D/case)
                        if terminal is None: run_step('search',case,session_path,session,deadline,report)
                        else: report.update('SEALED_SEARCH_REUSED',case=case,charts=terminal['charts'])
                        run_step('replay',case,session_path,session,deadline,report)
                        result = verified(case); results.append(result)
                        report.update('CASE_VERIFIED',case=case,rank_lower_bound=result['rank_lower_bound'],gain=result['gain'])
                    except ResourceStop as exc:
                        stopped.append({'case':case,'reason':str(exc)})
                        report.update('CASE_RESOURCE_STOP',case=case,error=str(exc))
                        if time.monotonic() >= deadline: break
                summary = c.result_summary(results)
                report.update('COMPLETE_EIGHT_PILOT' if len(results) == 8 else 'STOPPED_RESOURCE_BUDGET',
                              results=results,resource_stops=stopped,summary=summary)
                atomic(c.D/'summary.json',{'results':results,'resource_stops':stopped,'summary':summary,
                       'roster_sha256':sha(c.D/'roster.json')})
            except ResourceStop as exc:
                report.update('STOPPED_RESOURCE_BUDGET',error=str(exc))
            except BaseException as exc:
                report.update('STOPPED_REVIEW_REQUIRED',error=repr(exc)); raise
    finally: os.close(lease_fd)


def status():
    r = dict(read(STATE)) if STATE.exists() else {'status':'NOT_LAUNCHED'}
    r['process_alive'] = alive(r)
    if r['status'] in ACTIVE and not r['process_alive']: r['effective_status'] = 'STOPPED_REVIEW_REQUIRED'
    if r.get('attempt'):
        path = c.ROOT/r['attempt']/'supervisor.json'
        if path.exists(): r['supervisor'] = read(path)
    if r.get('worker_log'):
        text = tail(c.ROOT/r['worker_log'],8192)
        r['progress_tail'] = text.splitlines()[-5:] if text else []
    if r.get('results'):
        r['results'] = [{k:x.get(k) for k in ('case','rank_lower_bound','gain','charts','stop_reason')} for x in r['results']]
    if r.get('error'): r['error'] = r['error'][-4000:]
    print(json.dumps(r,sort_keys=True,indent=2))
    if (c.D/'roster.json').exists():
        for row in read(c.D/'roster.json')['pilot']:
            folder = c.D/row['id']
            if (folder/'trial-verified.json').exists():
                out = read(folder/'trial-verified.json'); print(row['id'],'VERIFIED',out['rank_lower_bound'],out['stop_reason'])
            elif (folder/'replay-M17/terminal.json').exists(): print(row['id'],'SEARCH_SEALED_REPLAY_PENDING')
            else: print(row['id'],'NOT_YET_VERIFIED')


def diagnose():
    row = read(STATE) if STATE.exists() else {}
    status()
    if row.get('worker_log'): print(tail(c.ROOT/row['worker_log']) or 'No worker log yet')


def stop():
    row = read(STATE) if STATE.exists() else {}
    if not alive(row): print('No matching live controller'); return
    require(hasattr(os,'pidfd_open') and hasattr(signal,'pidfd_send_signal'),'safe pidfd signalling unavailable')
    fd = os.pidfd_open(row['pid'])
    try:
        require(alive(row),'controller identity changed before stop')
        signal.pidfd_send_signal(fd,signal.SIGTERM)
    finally: os.close(fd)
    print('Stop requested; checkpoints retained')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('launch','resume','status','diagnose','stop','worker'))
    p.add_argument('--hours',type=float,default=10); p.add_argument('--session',type=Path); p.add_argument('--lease-fd',type=int)
    a = p.parse_args()
    if a.action in ('launch','resume'): launch(a.hours,a.action == 'resume')
    elif a.action == 'status': status()
    elif a.action == 'diagnose': diagnose()
    elif a.action == 'stop': stop()
    else: worker(a.session,a.lease_fd)


if __name__ == '__main__': main()
