#!/usr/bin/env python3
"""Prepare/produce orbit8044 M18 seeds and drain a separate, bounded V3 queue.

The producer and V3 consumer have separate locks, ledgers and total budgets.
Resume reuses immutable stages; an interrupted unreceipted stage consumes its
reservation and stays censored. Queue publication is idempotent across crashes.
"""
import argparse
from collections import Counter
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys

import det1092_funnel as f
import orbit8044_seed_factory as factory
from v3_warm_support import atomic, bindings, read, require, sha
from research_runtime.supervisor import run, Limits

SELF = Path(__file__).resolve()
WORKER = f.CAS/'orbit8044_seed_factory_worker.py'
CHECKER = f.CAS/'verify_det1092_funnel_seed.sage'
ENV = {**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1', 'PYTHONUNBUFFERED':'1'}


class Stages:
    def __init__(self, folder, lane, total, rss):
        self.folder, self.lane, self.total, self.rss = folder, lane, total, rss
        self.path = folder/f'{lane}-controller.json'
        self.data = read(self.path) if self.path.exists() else dict(stages={}, cursor=0)
        self.data.update(status='RUNNING', pid=os.getpid(), active_stage=None)
        self.exhausted = False
        self.save()

    def save(self):
        atomic(self.path, self.data)

    def stage(self, name, command, seconds, terminal):
        old = self.data['stages'].get(name)
        receipt = self.folder/'logs'/f'{self.lane}-{name}.supervisor.json'
        if old and old['status'] == 'RUNNING':
            record = read(receipt) if receipt.exists() else None
            if record and record.get('outcome') == 'completed':
                require(record['command'] == command and record['limits']['wall_seconds'] == old['reserved_seconds'], 'recovered receipt changed')
                self.finish(name, record, terminal)
            else:
                old.update(status='CENSORED_INTERRUPTION', charged_seconds=old['reserved_seconds'])
                self.save()
            old = self.data['stages'][name]
        if old:
            if old['status'] == 'PASS':
                require(terminal.exists() and sha(terminal) == old['terminal_sha256'], 'completed stage output changed')
                return True
            return False
        spent = sum(r.get('charged_seconds', r['reserved_seconds']) for r in self.data['stages'].values())
        remaining = self.total-spent
        if remaining < 1:
            self.exhausted = True
            self.data.update(status='BUDGET_EXHAUSTED', active_stage=None)
            self.save()
            return False
        reserved = min(seconds, remaining)
        self.data['stages'][name] = dict(status='RUNNING', reserved_seconds=reserved, command=command)
        self.data['active_stage'] = name
        self.save()
        record = run(command, limits=Limits(reserved, self.rss),
            log_path=receipt.with_name(f'{self.lane}-{name}.log'),
            checkpoint_path=receipt, cwd=f.ROOT, env=ENV)
        return self.finish(name, record, terminal)

    def finish(self, name, record, terminal):
        passed = record['outcome'] == 'completed' and record['returncode'] == 0 and terminal.exists()
        old = self.data['stages'][name]
        old.update(status='PASS' if passed else 'ERROR_OR_CENSORED',
            charged_seconds=record['wall_seconds'], supervision=record,
            terminal_sha256=sha(terminal) if passed else None)
        self.save()
        return passed


def checked(folder):
    p = f.protocol(folder)
    require(p['profile'] == 'orbit8044_seed_factory', 'wrong factory protocol')
    require(sha(folder/'parameter-chart.json') == p['factory']['chart_sha256'], 'chart changed')
    require(sha(folder/'known-equations.json') == p['factory']['known_equations_sha256'], 'deduplication snapshot changed')
    bindings(f.ROOT, read(folder/'parameter-chart.json')['inputs'])
    return p


def freeze(parent, folder, max_height, maximum, exclude_runs):
    require(not (folder/'protocol.json').exists(), 'already frozen; use the same finite resume')
    list(factory.parameters(max_height, maximum))  # Validate finite bounds before invoking Sage.
    p = f.protocol(parent)
    sage = shutil.which('sage'); require(sage is not None, 'Sage unavailable')
    result = run([sage, '-python', str(WORKER), 'chart', '--output', str(folder/'parameter-chart.json')],
        limits=Limits(120, 1024**3), log_path=folder/'logs/chart.log',
        checkpoint_path=folder/'logs/chart.supervisor.json', cwd=f.ROOT, env=ENV)
    require(result['outcome'] == 'completed' and result['returncode'] == 0, 'chart preparation failed')
    known = []
    for excluded in exclude_runs:
        for path in sorted((excluded/'seeds').glob('*/m18.json')):
            packet = read(path)
            require(packet['status'] == 'CERTIFIED_M18', 'invalid existing seed packet')
            known.append(dict(case=excluded.name+'/'+path.parent.name, curve=packet['curve'],
                j=str(factory.rational_j(packet['curve'])), packet_path=str(path.relative_to(f.ROOT)),
                packet_sha256=sha(path)))
    atomic(folder/'known-equations.json', dict(status='FROZEN_EQUATION_ONLY_DEDUP_SNAPSHOT', equations=known), immutable=True)
    inputs = dict(p['inputs'])
    inputs.update({str(q.relative_to(f.ROOT)):sha(q) for q in [parent/'protocol.json', folder/'parameter-chart.json', folder/'known-equations.json']})
    p.update(schema='orbit8044-seed-factory.v1', profile='orbit8044_seed_factory',
        domain='orbit8044-parametrized-seed-factory-v1', population=0, maximum_draws=0,
        sources=f.source_bindings(), inputs=inputs,
        factory=dict(max_height=max_height, maximum_parameters=maximum,
            enumeration='Increasing(max(abs(a),b),b,a); b>0; gcd(a,b)=1; includes0/1; finite rational u only',
            chart_sha256=sha(folder/'parameter-chart.json'), known_equations_sha256=sha(folder/'known-equations.json'),
            seed_seconds=120, verification_seconds=120, seed_campaign_seconds=14400,
            seed_rss_bytes=1024**3, seed_workers=1, seed_quartic_search_charts=0,
            v3_seconds_per_fibre=14400, v3_campaign_seconds=43200, v3_rss_bytes=3*1024**3, v3_workers=1,
            maximum_attempts_per_stage=1, queue_rule='Every independently certified rational-isomorphism class; no score cutoff or seed quota',
            deduplication='Exact j buckets followed by exact rational isomorphism; twists are retained',
            missing_or_failed_independence='UNKNOWN; not a whole-curve exclusion'),
        scope='Height-ordered orbit8044 conic specializations, exact two-branch construction, exact generic-word filtering, independently certified M18 and durable deduplicated V3 queue. Producer and V3 consumer have separate explicit launches and finite budgets. Existing campaigns remain unchanged.')
    atomic(folder/'protocol.json', p, immutable=True)
    print('FROZEN_ORBIT8044_FACTORY', max_height, maximum, flush=True)


def produce(folder):
    p = checked(folder); plan = p['factory']
    with (folder/'factory.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        jobs = Stages(folder, 'factory', plan['seed_campaign_seconds'], plan['seed_rss_bytes'])
        registry = list(read(folder/'known-equations.json')['equations'])
        for address in factory.parameters(plan['max_height'], plan['maximum_parameters']):
            case = address['id']; seed = folder/'seeds'/case
            passed = jobs.stage('seed-'+case,
                [shutil.which('sage'), '-python', str(WORKER), 'confirm', '--directory', str(folder), '--case', case],
                plan['seed_seconds'], seed/'result.json')
            if jobs.exhausted:
                break
            if passed and read(seed/'result.json')['status'] == 'CERTIFIED_M18':
                verified = jobs.stage('verify-'+case,
                    [shutil.which('sage'), '-python', str(CHECKER), '--seed', str(seed), '--output', str(seed/'standalone-replay.json')],
                    plan['verification_seconds'], seed/'standalone-replay.json')
                if jobs.exhausted:
                    break
                if verified:
                    admission = factory.queue_admission(folder, case, registry)
                    if admission['status'] == 'QUEUED_M18':
                        registry.append(admission)
            jobs.data['cursor'] = address['index']+1
            jobs.save()
        jobs.data.update(active_stage=None, status='BUDGET_EXHAUSTED' if jobs.exhausted else
            ('COMPLETE_FINITE_FACTORY' if all(r['status'] == 'PASS' for r in jobs.data['stages'].values()) else 'PARTIAL_FINITE_FACTORY'))
        jobs.save()
        atomic(folder/'factory-summary.json', status(folder))


def drain(folder):
    p = checked(folder); plan = p['factory']
    with (folder/'v3.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        jobs = Stages(folder, 'v3', plan['v3_campaign_seconds'], plan['v3_rss_bytes'])
        for queue_path in sorted((folder/'v3-queue').glob('*.json')):
            queued = read(queue_path); case = queued['case']; seed = folder/'seeds'/case
            require(queued['protocol_sha256'] == sha(folder/'protocol.json') and
                    queued['seed_sha256'] == sha(seed/'m18.json') and
                    queued['proof_sha256'] == sha(seed/'standalone-replay.json'), 'V3 admission binding changed')
            jobs.stage(case, [shutil.which('sage'), '-python', str(f.CAS/'det1092_funnel_worker.py'),
                'amplify', '--directory', str(folder), '--case', case],
                plan['v3_seconds_per_fibre'], folder/'amplifiers'/case/'queue.json')
            if jobs.exhausted:
                break
        jobs.data.update(active_stage=None, status='BUDGET_EXHAUSTED' if jobs.exhausted else
            ('DRAINED_CURRENT_QUEUE_SNAPSHOT' if all(r['status'] == 'PASS' for r in jobs.data['stages'].values())
             else 'PARTIAL_V3_QUEUE_SNAPSHOT'))
        jobs.save()


def status(folder):
    p = read(folder/'protocol.json')
    counts = Counter(read(path)['status'] for path in (folder/'seeds').glob('*/result.json'))
    admissions = Counter(read(path)['status'] for path in (folder/'admissions').glob('*.json'))
    controllers = {}
    for lane in ['factory', 'v3']:
        path = folder/f'{lane}-controller.json'
        if path.exists():
            d = read(path)
            controllers[lane] = {k:d.get(k) for k in ['status', 'pid', 'cursor', 'active_stage']}
            controllers[lane]['charged_seconds'] = sum(r.get('charged_seconds', r['reserved_seconds']) for r in d['stages'].values())
    completed = len(list((folder/'amplifiers').glob('*/queue.json')))
    return dict(status='ORBIT8044_FACTORY_SNAPSHOT', controllers=controllers,
        maximum_parameters=p['factory']['maximum_parameters'], max_height=p['factory']['max_height'],
        seed_outcomes=dict(counts), admissions=dict(admissions),
        v3_replayed=completed, queued_without_terminal_replay=admissions.get('QUEUED_M18', 0)-completed,
        boundary='Finite certified subgroup factory; failed proofs and bounded misses are not rank upper bounds.')


def launch(folder, lane):
    checked(folder)
    with (folder/f'{lane}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with (folder/f'{lane}-controller.log').open('ab', buffering=0) as log:
            process = subprocess.Popen([sys.executable, str(SELF), 'produce' if lane == 'factory' else 'drain', '--directory', str(folder)],
                cwd=f.ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, env=ENV, start_new_session=True)
    print('LAUNCHED_ORBIT8044', lane, process.pid, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['freeze', 'produce', 'launch', 'drain', 'launch-v3', 'status'])
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--parent-run', type=Path)
    p.add_argument('--exclude-run', type=Path, action='append', default=[])
    p.add_argument('--max-height', type=int, default=64)
    p.add_argument('--maximum-parameters', type=int, default=4096)
    a = p.parse_args(); folder = a.directory.resolve()
    require(folder.is_relative_to(f.ROOT/'artifacts/local/elliptic-curves'), 'local evidence directory required')
    if a.action == 'freeze':
        require(a.parent_run is not None, 'parent protocol required')
        freeze(a.parent_run.resolve(), folder, a.max_height, a.maximum_parameters, [r.resolve() for r in a.exclude_run])
    elif a.action == 'status':
        import json
        print(json.dumps(status(folder), indent=2, sort_keys=True))
    elif a.action == 'produce':
        produce(folder)
    elif a.action == 'drain':
        drain(folder)
    else:
        launch(folder, 'factory' if a.action == 'launch' else 'v3')
