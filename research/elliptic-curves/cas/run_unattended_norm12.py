#!/usr/bin/env python3
"""Detached, token-free scheduling of sealed V3 cached continuations.

One worker at a time. Every batch searches, independently replays and exports.
Gains stop scheduling; complete-cloud reconciliation runs before reporting one.
Partial/failed batches are preserved and never automatically restarted.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
import zipfile

from export_productive_seed_result import export
from queue_verified_complement_pass import record as cursor_record
from research_runtime.supervisor import Limits, run as supervise

ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def process_token(pid):
    try:
        fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
        return None if fields[0] == 'Z' else fields[19]
    except (OSError, IndexError):
        return None


def next_action(initial_rank, terminal):
    if terminal['stop_reason'] == 'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION':
        return 'RECONCILE'
    if terminal['rank_lower_bound'] > initial_rank:
        return 'STOP_CERTIFIED_GAIN'
    if terminal['rank_lower_bound'] != initial_rank:
        raise ArithmeticError('rank decreased')
    if terminal['stop_reason'] == 'CHART_BUDGET_EXHAUSTED':
        return 'CONTINUE'
    if terminal['stop_reason'] == 'FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN':
        return 'EXHAUSTED'
    raise ArithmeticError('unrecognized terminal stop condition')


def status(folder):
    state = read(folder / 'state.json')
    token = process_token(state.get('pid', -1))
    state['process_alive'] = token is not None and token == state.get('process_token')
    if state['status'] == 'RUNNING' and not state['process_alive']:
        state['observed_status'] = 'INTERRUPTED_REQUIRES_REVIEW'
    return state


def worker(folder):
    lock = (folder / 'controller.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config_path = folder / 'config.json'
    config = read(config_path)
    config_hash = sha(config_path)
    state = read(folder / 'state.json')
    if state['status'] != 'CREATED':
        raise RuntimeError('preserve existing campaign; no automatic restart')
    state.update(status='RUNNING', pid=os.getpid(), process_token=process_token(os.getpid()),
                 started_at=time.time(), batches=[], jobs=config['jobs'])

    def checkpoint(**changes):
        state.update(changes, updated_at=time.time())
        save(folder / 'state.json', state)

    def phase(command, supervision, seconds):
        if supervision.exists():
            raise FileExistsError('preserve existing supervision')
        supervision.mkdir()
        checkpoint(phase=supervision.name)
        result = supervise(command, limits=Limits(wall_seconds=seconds, rss_bytes=3*1024**3),
                           log_path=supervision/'worker.log',
                           checkpoint_path=supervision/'supervisor.json')
        if result['outcome'] != 'completed' or result['returncode'] != 0:
            raise RuntimeError(f"bounded worker failed: {result['outcome']}")
        return result

    checkpoint()
    try:
        while len(state['batches']) < config['max_batches']:
            active = [j for j in state['jobs'] if j['status'] == 'READY']
            if not active:
                checkpoint(status='STOP_ALL_COVERAGE_EXHAUSTED')
                return
            for job in active:
                if (folder/'STOP').exists():
                    checkpoint(status='STOP_REQUESTED')
                    return
                if time.time()-state['started_at'] >= config['max_hours']*3600:
                    checkpoint(status='STOP_TIME_LIMIT')
                    return
                if len(state['batches']) >= config['max_batches']:
                    break
                if shutil.disk_usage(folder).free < config['min_free_gib']*1024**3:
                    checkpoint(status='STOP_DISK_RESERVE')
                    return
                if sha(config_path) != config_hash or any(sha(ROOT/n) != h for n,h in config['sources'].items()):
                    raise ArithmeticError('controller configuration or sources changed')
                parent = ROOT/job['parent']
                cursor_record(parent)  # sealed replay, budget, receipts, supervisors
                before = read(parent/'terminal.json')['rank_lower_bound']
                index = len(state['batches'])+1
                target = folder/f"batch-{index:04d}-{job['id']}"
                if target.exists():
                    raise FileExistsError('preserve partial batch')
                checkpoint(current_job=job['id'], current_batch=str(target.relative_to(ROOT)))
                runner = CAS/'run_complement_cached_seed_v3.py'
                for mode in ('search', 'replay'):
                    suffix = '-supervision' if mode == 'search' else '-replay-supervision'
                    phase(['sage','-python',str(runner),mode,'--folder',str(target),'--parent',str(parent)],
                          target.with_name(target.name+suffix), 1800)
                    if mode == 'search':
                        protocol = read(target/'protocol.json')
                        with zipfile.ZipFile(target/'frozen-sources.zip','w',zipfile.ZIP_DEFLATED) as archive:
                            for name, digest in protocol['sources'].items():
                                if sha(ROOT/name) != digest:
                                    raise ArithmeticError('search source changed')
                                archive.write(ROOT/name, name)
                result_path = target/'export.json'
                export(target, result_path)
                terminal = read(target/'terminal.json')
                action = next_action(before, terminal)
                if action == 'RECONCILE':
                    for mode in ('reconcile', 'reconcile-replay'):
                        phase([sys.executable,str(CAS/'reconcile_verified_v3_cloud.py'),
                               '--folder',str(target),'--output',str(target/'reconciled.json')],
                              target.with_name(target.name+'-'+mode+'-supervision'), 120)
                    action = 'STOP_CERTIFIED_GAIN'
                if action == 'CONTINUE':
                    cursor = cursor_record(target)
                    save(target/'suffix-queue.json', cursor)
                    if cursor['next_action'] == 'NO_UNVISITED_CENTRE_SUFFIX':
                        action = 'EXHAUSTED'
                entry = dict(job=job['id'], parent=str(parent.relative_to(ROOT)),
                             result=str(result_path.relative_to(ROOT)), result_sha256=sha(result_path),
                             rank_lower_bound=terminal['rank_lower_bound'],
                             cumulative_calls=terminal['charts'], action=action)
                if (target/'reconciled.json').exists():
                    entry['reconciled'] = str((target/'reconciled.json').relative_to(ROOT))
                    entry['rank_lower_bound'] = read(target/'reconciled.json')['rank_lower_bound']
                state['batches'].append(entry)
                job.update(parent=str(target.relative_to(ROOT)),
                           status='EXHAUSTED' if action == 'EXHAUSTED' else 'READY')
                checkpoint(phase='BATCH_REPLAYED')
                print(json.dumps(entry), flush=True)
                if action == 'STOP_CERTIFIED_GAIN':
                    checkpoint(status=action)
                    return
        checkpoint(status='STOP_BATCH_LIMIT')
    except BaseException:
        checkpoint(status='STOP_FAILURE_REQUIRES_REVIEW', error=traceback.format_exc())
        raise


def start(folder, parents, max_hours, max_batches, min_free_gib):
    if not 0 < max_hours <= 168 or not 0 < max_batches <= 1000 or min_free_gib < 1:
        raise ValueError('invalid declared campaign limits')
    jobs = []
    for parent in parents:
        parent = parent.resolve()
        cursor = cursor_record(parent)
        if cursor['next_action'] != 'CONSIDER_REMAINING_SUFFIX':
            raise ValueError('each initial job needs unvisited verified coverage')
        jobs.append(dict(id=f'job{len(jobs)+1}', parent=str(parent.relative_to(ROOT)), status='READY'))
    if not jobs or len({j['parent'] for j in jobs}) != len(jobs):
        raise ValueError('nonempty distinct parent list required')
    folder.mkdir(exist_ok=False)
    sources = [Path(__file__),CAS/'export_productive_seed_result.py',CAS/'queue_verified_complement_pass.py',
               CAS/'reconcile_verified_v3_cloud.py',CAS/'research_runtime/supervisor.py']
    save(folder/'config.json',dict(jobs=jobs,max_hours=max_hours,max_batches=max_batches,
                                  min_free_gib=min_free_gib,sources={str(p.relative_to(ROOT)):sha(p) for p in sources},
                                  policy='Single worker; round robin;100 calls per batch; replay before continuation; '
                                  'stop on gain/failure/exhaustion/limits; STOP file is honored between batches.'))
    save(folder/'state.json',dict(status='CREATED'))
    with (folder/'controller.log').open('ab',buffering=0) as log:
        proc = subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'worker','--folder',str(folder)],
                                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                start_new_session=True,close_fds=True)
    print(json.dumps(dict(pid=proc.pid,folder=str(folder))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['start','worker','status','stop'])
    parser.add_argument('--folder',required=True,type=Path)
    parser.add_argument('--parent',action='append',type=Path,default=[])
    parser.add_argument('--max-hours',type=float,default=24)
    parser.add_argument('--max-batches',type=int,default=100)
    parser.add_argument('--min-free-gib',type=int,default=10)
    a=parser.parse_args(); folder=a.folder.resolve()
    if a.mode=='start': start(folder,a.parent,a.max_hours,a.max_batches,a.min_free_gib)
    elif a.mode=='worker': worker(folder)
    elif a.mode=='stop': (folder/'STOP').touch(); print('Stop requested after the current batch and replay.')
    else: print(json.dumps(status(folder),indent=2))
