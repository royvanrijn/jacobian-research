#!/usr/bin/env python3
"""Finite, checkpointed controller for the determinant1092 three-stage funnel."""
from __future__ import annotations

import argparse
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import det1092_funnel as f
from v3_warm_support import atomic, read, require, sha
from research_runtime.supervisor import run, Limits

SELF = Path(__file__).resolve()


def sage():
    path = shutil.which('sage')
    require(path is not None, 'Sage launcher required')
    return path


def work(folder):
    p = f.protocol(folder)
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        ledger_path = folder/'controller.json'
        ledger = read(ledger_path) if ledger_path.exists() else dict(stages=[], status='RUNNING')
        ledger.update(status='RUNNING', pid=os.getpid())
        atomic(ledger_path, ledger)
        resources = p['resources']

        def stage(name, command, seconds, terminal):
            old = [r for r in ledger['stages'] if r['name'] == name]
            if old and old[-1]['passed']:
                require(terminal.exists() and sha(terminal) == old[-1]['terminal_sha256'], 'completed stage output changed: '+name)
                return True
            reservations = ledger.setdefault('interrupted_reservations', [])
            spent = sum(r['supervision']['wall_seconds'] for r in ledger['stages'])+sum(r['seconds'] for r in reservations)
            remaining = resources['campaign_seconds']-spent
            # Repeated resumes share the declared total allocation for this job.
            job_remaining = seconds-sum(r['supervision']['wall_seconds'] for r in old)-sum(r['seconds'] for r in reservations if r['name'] == name)
            if min(remaining,job_remaining) < 1:
                ledger.update(status='BUDGET_EXHAUSTED', active_stage=name)
                deferred = ledger.setdefault('deferred_stages', [])
                if name not in deferred:
                    deferred.append(name)
                atomic(ledger_path,ledger); return False
            limits = Limits(min(remaining,job_remaining), resources['rss_bytes'])
            reservation = dict(name=name, seconds=limits.wall_seconds)
            reservations.append(reservation)
            ledger['active_stage'] = name; atomic(ledger_path,ledger)
            attempt = len(old)
            log = folder/'logs'/f'{name}-{attempt:02d}.log'
            record = run(command, limits=limits, log_path=log,
                         checkpoint_path=log.with_suffix('.supervisor.json'), cwd=f.ROOT,
                         env={**os.environ, 'PYTHONUNBUFFERED':'1', 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1'})
            passed = record['outcome'] == 'completed' and record['returncode'] == 0 and terminal.exists()
            reservations.remove(reservation)
            ledger['stages'].append(dict(name=name, passed=passed, supervision=record,
                terminal_sha256=sha(terminal) if passed else None))
            ledger['status'] = 'RUNNING' if passed else 'ERROR_OR_CENSORED'
            if passed and name in ledger.get('deferred_stages', []):
                ledger['deferred_stages'].remove(name)
            atomic(ledger_path,ledger)
            print('FUNNEL_STAGE',name,ledger['status'],flush=True)
            return passed

        scan_cmd = [sys.executable,str(f.CAS/'det1092_funnel.py')]
        sage_cmd = [sage(),'-python',str(f.CAS/'det1092_funnel_worker.py')]
        suffix = ['--directory',str(folder)]
        jobs = [('tables',sage_cmd+['tables']+suffix,120,folder/'tables-replay.json'),
                ('intake',scan_cmd+['scan']+suffix,resources['intake_seconds'],folder/'selection.json'),
                ('audit',scan_cmd+['audit']+suffix,resources['intake_seconds'],folder/'intake-replay.json')]
        for args in jobs:
            if not stage(*args):
                return
        if p['profile'] == 'smoke':
            if not stage('conic-control',sage_cmd+['conic-control']+suffix,120,
                         folder/'amplifiers/conic-control/preflight.json'):
                return
            if not stage('replay-conic-control',sage_cmd+['conic-control']+suffix,120,
                         folder/'amplifiers/conic-control/preflight.json'):
                return
        for row in read(folder/'selection.json')['seed_inputs']:
            case = row['id']; cf = folder/'seeds'/case
            if not stage('seed-'+case,sage_cmd+['seed']+suffix+['--case',case],resources['seed_seconds'],cf/'result.json'):
                continue
            if not stage('replay-seed-'+case,sage_cmd+['replay-seed']+suffix+['--case',case],resources['seed_seconds'],cf/'result.json'):
                continue
            if read(cf/'result.json')['status'] != 'CERTIFIED_M18':
                continue
            if p['profile'] == 'smoke':
                stage('preflight-'+case,sage_cmd+['preflight']+suffix+['--case',case],120,
                      folder/'amplifiers'/case/'preflight.json')
            else:
                stage('amplify-'+case,sage_cmd+['amplify']+suffix+['--case',case],
                      resources['amplifier_seconds']+resources['replay_seconds'],
                      folder/'amplifiers'/case/'queue.json')
        outcomes = {r['id']: read(folder/'seeds'/r['id']/'result.json')
                    for r in read(folder/'selection.json')['seed_inputs']
                    if (folder/'seeds'/r['id']/'result.json').exists()}
        queues = {path.parent.name:read(path) for path in (folder/'amplifiers').glob('*/queue.json')}
        unresolved = [r['name'] for r in ledger['stages'] if not r['passed'] and not any(
            q['name'] == r['name'] and q['passed'] for q in ledger['stages'])]
        unresolved = sorted(set(unresolved+ledger.get('deferred_stages', [])))
        require(len(outcomes) == len(read(folder/'selection.json')['seed_inputs']) or unresolved,
                'missing seed outcomes without a recorded unresolved stage')
        summary = dict(status='COMPLETE_FINITE_FUNNEL' if not unresolved else 'PARTIAL_FINITE_FUNNEL',
                       profile=p['profile'],parameters=read(folder/'selection.json')['parameters'],
                       seed_outcomes=outcomes,queues=queues,unresolved_stages=unresolved,
                       constructive_control_in_population_yield=False,
                       production_amplification_tested=p['profile'] == 'production' and bool(queues))
        atomic(folder/'summary.json',summary)
        ledger.update(status=summary['status'],active_stage=None); atomic(ledger_path,ledger)


def launch(folder):
    f.protocol(folder)
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        with (folder/'controller.log').open('ab',buffering=0) as log:
            proc = subprocess.Popen([sys.executable,str(SELF),'worker','--directory',str(folder)],
                                    stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                    start_new_session=True,cwd=f.ROOT,
                                    env={**os.environ,'PYTHONUNBUFFERED':'1'})
        atomic(folder/'launch.json',dict(pid=proc.pid,unix=time.time()))
    print('LAUNCHED',proc.pid,str(folder),flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['freeze','run','launch','resume','worker','status'])
    ap.add_argument('--directory',type=Path,required=True)
    ap.add_argument('--profile',choices=['smoke','production'],default='smoke')
    a = ap.parse_args(); folder = a.directory.resolve()
    require(folder.is_relative_to(f.ROOT/'artifacts/local/elliptic-curves'), 'use local artifact namespace')
    if a.action == 'freeze':
        f.freeze(folder,a.profile)
    elif a.action in ('launch','resume'):
        launch(folder)
    elif a.action in ('run','worker'):
        work(folder)
    else:
        for name in ('controller.json','progress.json','summary.json'):
            if (folder/name).exists():
                print(name, f.packed(read(folder/name)).decode())


if __name__ == '__main__':
    main()
