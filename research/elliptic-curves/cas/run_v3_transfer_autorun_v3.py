#!/usr/bin/env python3
"""Detached autorun for transfer wrapper v3, with diagnostics and safe repair/resume.

The controller is operational only. It surfaces nested worker failures, can
archive a failed supervised attempt for ordinary resume, and can archive the
entire prepared transfer when a source-level repair changes the frozen protocol.
No evidence is deleted and no failed budget is silently retried.
"""
from importlib.machinery import SourceFileLoader
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

SELF = Path(__file__).resolve()
CAS = SELF.parent
base = SourceFileLoader('v3_transfer_autorun_v1_preserved', str(CAS/'run_v3_transfer_autorun.py')).load_module()
base.TRANSFER = base.LOCAL/'v3-transfer-11952-v3'
base.AUTO = base.LOCAL/'v3-transfer-autorun-v3'
base.STATE = base.AUTO/'state.json'
base.LOG = base.AUTO/'autorun.log'
base.CAMPAIGN = CAS/'v3_transfer_campaign_v3.sage'
base.__file__ = str(SELF)


def tail(path, lines=100):
    path = Path(path)
    if not path.exists(): return None
    try: return '\n'.join(path.read_text(errors='replace').splitlines()[-lines:])
    except OSError as exc: return f'<cannot read {path}: {exc!r}>'


def supervisor(path):
    path=Path(path)
    if not path.exists(): return None
    try:return json.loads(path.read_text())
    except Exception as exc:return {'parse_error':repr(exc)}


def current_case():
    roster=base.TRANSFER/'roster.json'
    if not roster.exists(): return None
    for spec in base.CASE_SPECS:
        if not (base.TRANSFER/spec['id']/'verified.json').exists(): return spec
    return None


def failure_snapshot():
    candidates=[]
    prep=base.TRANSFER.parent/(base.TRANSFER.name+'-preparation')
    candidates.append(('prepare',prep/'prepare.supervisor.json',prep/'prepare.log'))
    for spec in base.CASE_SPECS:
        folder=base.TRANSFER/spec['id']
        candidates += [(spec['id']+':run',folder/'run.supervisor.json',folder/'run.log'),
                       (spec['id']+':replay',folder/'replay.supervisor.json',folder/'replay.log')]
    found=[]
    for label,sup,log in candidates:
        if not sup.exists() and not log.exists():continue
        mtime=max(p.stat().st_mtime for p in (sup,log) if p.exists())
        found.append((mtime,label,sup,log))
    if not found:return None
    _,label,sup,log=max(found)
    return {'label':label,'supervisor_path':str(sup.relative_to(base.ROOT)) if sup.exists() else None,
            'supervisor':supervisor(sup),'log_path':str(log.relative_to(base.ROOT)) if log.exists() else None,
            'log_tail':tail(log,120)}


def status():
    value=base.read_state()
    if value is None:print('NOT_LAUNCHED')
    else:
        pid=int(value.get('pid',0) or 0);value=dict(value)
        value['process_alive']=bool(pid and base.live(pid,value.get('proc_token')))
        print(json.dumps(value,indent=2,sort_keys=True))
    snap=failure_snapshot()
    if snap:
        print('\nNESTED_ATTEMPT');print(json.dumps(snap,indent=2,sort_keys=True))


def diagnose():
    snap=failure_snapshot()
    print('NO_SUPERVISED_ATTEMPT_FOUND' if snap is None else json.dumps(snap,indent=2,sort_keys=True))


def digest_paths(paths):
    h=hashlib.sha256()
    for p in paths:
        p=Path(p);h.update(str(p).encode())
        if p.is_file():h.update(p.read_bytes())
    return h.hexdigest()[:16]


def move_preserved(src,dst):
    src=Path(src);dst=Path(dst)
    if not src.exists():return
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():raise RuntimeError(f'archive destination already exists: {dst}')
    src.replace(dst)


def assert_controller_stopped():
    state=base.read_state() or {};pid=int(state.get('pid',0) or 0)
    if pid and base.live(pid,state.get('proc_token')):
        raise RuntimeError(f'controller pid {pid} is still alive; refusing concurrent recovery')


def recover_failed_attempt():
    """Preserve one failed worker/replay while keeping the same prepared protocol."""
    assert_controller_stopped();spec=current_case()
    if spec is None:raise RuntimeError('no unfinished prepared case to recover')
    folder=base.TRANSFER/spec['id'];terminal=folder/'replay-M17/terminal.json'
    action='replay' if terminal.exists() else 'run';sup=folder/(action+'.supervisor.json');log=folder/(action+'.log')
    if not sup.exists():raise RuntimeError(f'no preserved failed {action} supervisor for {spec["id"]}; run diagnose first')
    report=supervisor(sup) or {}
    if report.get('outcome')=='completed':raise RuntimeError(f'{action} supervisor says completed')
    stamp=time.strftime('%Y%m%dT%H%M%SZ',time.gmtime())+'-'+digest_paths([sup,log])
    archive=folder/'failed-attempts'/stamp;archive.mkdir(parents=True,exist_ok=False)
    for name in (action+'.supervisor.json',action+'.log',action+'.stderr.log',action+'.stdin'):
        move_preserved(folder/name,archive/name)
    if action=='run':
        replay=folder/'replay-M17'
        if replay.exists():
            for epoch in sorted(replay.glob('epoch-*')):
                if (epoch/'stage.json').exists():continue
                if any(epoch.iterdir()):move_preserved(epoch,archive/epoch.name)
    history=base.AUTO/'history';history.mkdir(parents=True,exist_ok=True)
    if base.STATE.exists():move_preserved(base.STATE,history/(stamp+'-state.json'))
    print(json.dumps({'recovered_case':spec['id'],'action':action,'archive':str(archive.relative_to(base.ROOT))},indent=2))


def archive_protocol_attempt():
    """Archive the whole V3 transfer after a source-level protocol repair."""
    assert_controller_stopped()
    prep=base.TRANSFER.parent/(base.TRANSFER.name+'-preparation')
    sources=[p for p in (base.TRANSFER,prep,base.STATE,base.LOG) if p.exists()]
    stamp=time.strftime('%Y%m%dT%H%M%SZ',time.gmtime())+'-'+digest_paths(sources)
    root=base.LOCAL/'v3-transfer-11952-v3-failed-attempts'/stamp
    root.mkdir(parents=True,exist_ok=False)
    move_preserved(base.TRANSFER,root/'campaign')
    move_preserved(prep,root/'preparation')
    if base.STATE.exists():move_preserved(base.STATE,root/'autorun-state.json')
    if base.LOG.exists():move_preserved(base.LOG,root/'autorun.log')
    print(json.dumps({'archived_protocol_attempt':str(root.relative_to(base.ROOT)),
                      'reason':'source-level seed-replay repair; new preparation required'},indent=2))


def spawn_worker():
    base.AUTO.mkdir(parents=True,exist_ok=True);stream=base.LOG.open('ab',buffering=0)
    proc=subprocess.Popen([sys.executable,str(SELF),'worker'],cwd=base.ROOT,stdin=subprocess.DEVNULL,
        stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True,
        env={**os.environ,'PYTHONUNBUFFERED':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    token=None
    for _ in range(20):
        token=base.proc_token(proc.pid)
        if token is not None:break
        time.sleep(.01)
    base.atomic_json(base.STATE,{'status':'RELAUNCHED','pid':proc.pid,'proc_token':token,
        'launched_at_unix':time.time(),'log':str(base.LOG.relative_to(base.ROOT))})
    print(f'Relaunched detached V3 transfer controller pid={proc.pid}')


def resume():recover_failed_attempt();spawn_worker()
def repair_resume():archive_protocol_attempt();spawn_worker()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['launch','status','diagnose','resume','repair-resume','worker'])
    a=ap.parse_args()
    if a.action=='launch':base.launch()
    elif a.action=='status':status()
    elif a.action=='diagnose':diagnose()
    elif a.action=='resume':resume()
    elif a.action=='repair-resume':repair_resume()
    else:base.worker()

if __name__=='__main__':main()
