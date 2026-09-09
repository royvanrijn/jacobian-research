#!/usr/bin/env python3
"""Four-process scheduling for the frozen twelve-seed panel.

Only case scheduling changes. Each Sage subprocess calls the original panel's
context, search and independent replay, writing the original case checkpoints.
Dispatch follows the frozen seed order; completion order may differ. A separate
manifest binds this orchestration change without changing any V3 policy bytes.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import run_curve302_seed_universality_panel as panel
from v3_warm_support import atomic, read, require, sha

SELF = Path(__file__).resolve()
D = panel.D / 'parallel-v1'
MANIFEST = D / 'manifest.json'
STATE = D / 'state.json'
LOCK = D / 'controller.lock'


def acquire(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BaseException:
        os.close(fd)
        raise RuntimeError('already owned: ' + str(path))
    return fd


def identity(pid):
    try:
        raw = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
        if raw[0] == 'Z':
            return None
        return {'pid': pid, 'start_token': raw[19]}
    except (FileNotFoundError, ProcessLookupError):
        return None


def alive(record):
    return bool(record and identity(record['pid']) == record)


def check_bindings(manifest):
    for name, digest in manifest['bindings'].items():
        require(sha(panel.ROOT / name) == digest, 'changed frozen input/source: ' + name)
    require(list(panel.panel_order()) == manifest['order'], 'panel order changed')


def manifest(workers):
    preflight = read(panel.D / 'preflight.json')
    order = list(panel.panel_order())
    require(preflight['status'] == 'PASS_ZERO_CHART_PREFLIGHT' and
            preflight['charts'] == 0 and preflight['order'] == order,
            'complete original zero-chart preflight required')
    paths = {SELF, panel.SELF, panel.CLOSURE, panel.D / 'preflight.json'}
    for name in ('run_curve302_seeded_v3_amplifier.py', 'det1092_v3_worker.py',
                 'v3_warm_engine.py', 'v3_warm_replay.py', 'v3_warm_support.py'):
        paths.add(panel.CAS / name)
    for seed in order:
        folder = panel.D / seed
        paths.update(folder / name for name in ('protocol.json', 'seed-input.json', 'seed-proof.json'))
        policy = read(folder / 'protocol.json')
        for key in ('sources', 'seed_inputs'):
            for name, digest in policy.get(key, {}).items():
                path = panel.ROOT / name
                require(sha(path) == digest, 'original policy binding differs: ' + name)
                paths.add(path)
    return {'schema': 'curve302-seed-panel-parallel.v1', 'workers': workers,
            'order': order, 'bindings': {str(p.relative_to(panel.ROOT)): sha(p) for p in sorted(paths)},
            'change': 'Concurrent independent cases, dispatched in original order; unchanged per-case search and full replay.',
            'timing_boundary': 'Concurrent CPU scheduling can affect wall-time censoring; retain original per-chart limits and all censor records.'}


def verified(seed):
    folder = panel.D / seed
    path = folder / 'seeded-verified.json'
    if not path.exists():
        return None
    result = read(path)
    require(result['status'] == 'PASS_INDEPENDENT_SEEDED_V3_REPLAY' and
            result['seed_direction'] == seed and result['initial_rank'] == 18,
            'invalid replay receipt: ' + seed)
    require(result['terminal_sha256'] == sha(folder / 'replay-M17/terminal.json') and
            result['protocol_sha256'] == sha(folder / 'protocol.json'),
            'stale replay receipt: ' + seed)
    return result


def next_cases(order, completed, active, capacity):
    return [s for s in order if s not in completed and s not in active][:max(0, capacity)]


def run_case(seed):
    m = read(MANIFEST)
    check_bindings(m)
    require(seed in m['order'], 'seed outside frozen panel')
    fd = acquire(D / (seed + '.lock'))
    try:
        if verified(seed):
            return
        import det1092_v3_worker as searcher
        base = panel.load_base()
        ctx = base.context(seed)
        # panel.load_base retains the ORIGINAL panel.SELF for source binding.
        if not (ctx.folder / 'replay-M17/terminal.json').exists():
            print('PARALLEL_SEARCH', seed, flush=True)
            searcher.run_search(ctx)
        print('PARALLEL_REPLAY', seed, flush=True)
        base.verify_case(ctx)
        check_bindings(m)
        require(verified(seed) is not None, 'case did not publish a bound replay receipt')
    finally:
        os.close(fd)


def publish(state):
    atomic(STATE, state)
    # Original status/resume sees the live controller and cannot start a duplicate.
    atomic(panel.STATE, {'status': state['status'], 'pid': os.getpid(),
                        'updated_unix': time.time(), 'parallel_state': str(STATE),
                        'active_cases': list(state['active']), 'completed_cases': len(state['completed']),
                        'case_count': 12, 'order': state['order']})


def summarize(order):
    results = [verified(seed) for seed in order]
    require(all(results), 'incomplete panel cannot be summarized')
    summary = {'status': 'COMPLETE_SEED_UNIVERSALITY_PANEL', 'order': order,
               'known_winners_excluded': list(panel.KNOWN_WINNERS), 'results': results,
               'rank31_count': sum(r['rank_lower_bound'] >= 31 for r in results),
               'any_no_gain': any(r['rank_lower_bound'] == 18 for r in results),
               'minimum_terminal_rank': min(r['rank_lower_bound'] for r in results),
               'maximum_terminal_rank': max(r['rank_lower_bound'] for r in results),
               'total_charts': sum(r['charts'] for r in results), 'all_cases_present': True,
               'claim_boundary': 'Retrospective known-seed amplification panel; not a prospective selector or new rank claim.'}
    atomic(panel.D / 'summary.json', summary, immutable=True)


def controller(lock_fd):
    m = read(MANIFEST)
    check_bindings(m)
    state = {'status': 'PARALLEL_RUNNING', 'controller': identity(os.getpid()),
             'order': m['order'], 'workers': m['workers'], 'active': {}, 'completed': {},
             'started_unix': time.time()}
    procs = {}

    def interrupted(signum, frame):
        raise InterruptedError('parallel controller interrupted')

    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        while True:
            for seed in m['order']:
                if seed not in procs and seed not in state['completed']:
                    result = verified(seed)
                    if result:
                        state['completed'][seed] = {'rank': result['rank_lower_bound'], 'charts': result['charts']}
            for seed, proc in list(procs.items()):
                code = proc.poll()
                if code is None:
                    continue
                del procs[seed]
                del state['active'][seed]
                require(code == 0, f'{seed} worker failed with exit {code}; inspect its log')
                result = verified(seed)
                require(result is not None, seed + ' exited without replay receipt')
                state['completed'][seed] = {'rank': result['rank_lower_bound'], 'charts': result['charts']}
            if len(state['completed']) == len(m['order']):
                check_bindings(m)
                summarize(m['order'])
                state['status'] = 'COMPLETE_SEED_UNIVERSALITY_PANEL'
                publish(state)
                return
            for seed in next_cases(m['order'], state['completed'], procs, m['workers'] - len(procs)):
                check_bindings(m)
                with (D / (seed + '.log')).open('ab', buffering=0) as log:
                    proc = subprocess.Popen([panel.sage_launcher(), '-python', '-u', str(SELF),
                                             'case', '--seed', seed], cwd=panel.ROOT,
                                            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                            start_new_session=True,
                                            env={**os.environ, 'OPENBLAS_NUM_THREADS': '1',
                                                 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
                procs[seed] = proc
                state['active'][seed] = {'pid': proc.pid, 'started_unix': time.time(), 'log': str(D / (seed + '.log'))}
                print('PARALLEL_DISPATCH', seed, proc.pid, flush=True)
            state['updated_unix'] = time.time()
            publish(state)
            time.sleep(3)
    except BaseException as exc:
        # Every child owns a fresh session; stop only children we spawned.
        for proc in procs.values():
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
        for proc in procs.values():
            try:
                proc.wait(timeout=15)
            except subprocess.TimeoutExpired:
                # Keep controller ownership while a child remains alive.
                proc.wait()
        state.update(status='PARALLEL_STOPPED_REVIEW_REQUIRED', error=str(exc))
        publish(state)
        raise
    finally:
        os.close(lock_fd)


def retire_serial(old):
    pid = old.get('pid')
    if not pid or not identity(pid):
        return
    require(old.get('status') == 'REPLAYING', 'takeover requires a sealed search currently in replay')
    require((panel.D / old['case'] / 'replay-M17/terminal.json').exists(), 'current search is not sealed')
    token = identity(pid)
    cmd = Path(f'/proc/{pid}/cmdline').read_bytes().split(b'\0')
    require(os.fsencode(str(panel.SELF)) in cmd and b'worker' in cmd,
            'live PID is not the original panel worker')
    fd = os.pidfd_open(pid)
    try:
        require(identity(pid) == token, 'PID identity changed')
        signal.pidfd_send_signal(fd, signal.SIGTERM)
        deadline = time.monotonic() + 15
        while alive(token) and time.monotonic() < deadline:
            time.sleep(.1)
        require(not alive(token), 'original worker has not stopped; no parallel launch')
    finally:
        os.close(fd)


def launch(workers, takeover):
    require(1 <= workers <= 4, 'bounded panel concurrency is 1..4')
    fd = acquire(LOCK)
    try:
        if STATE.exists():
            saved = read(STATE)
            require(not alive(saved.get('controller')), 'parallel controller is alive')
            for child in saved.get('active', {}).values():
                require(not identity(child['pid']), 'a previous case worker may still be alive')
        m = manifest(workers)
        atomic(MANIFEST, m, immutable=True)
        old = read(panel.STATE)
        if old.get('pid') and identity(old['pid']):
            require(takeover, 'live serial worker requires --takeover')
            atomic(D / 'serial-state-before-takeover.json', old, immutable=True)
            retire_serial(old)
        with (D / 'controller.log').open('ab', buffering=0) as log:
            proc = subprocess.Popen([sys.executable, '-u', str(SELF), 'controller', '--lock-fd', str(fd)],
                                    cwd=panel.ROOT, stdin=subprocess.DEVNULL, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(fd,))
        atomic(panel.STATE, {'status': 'PARALLEL_LAUNCHED', 'pid': proc.pid,
                            'parallel_state': str(STATE), 'updated_unix': time.time()})
        print('PARALLEL_LAUNCHED', proc.pid, 'workers', workers, flush=True)
    finally:
        os.close(fd)


def status():
    s = read(STATE) if STATE.exists() else {'status': 'NOT_LAUNCHED', 'active': {}}
    s['controller_alive'] = alive(s.get('controller'))
    for seed, row in s['active'].items():
        row['process_alive'] = bool(identity(row['pid']))
        row['progress_tail'] = panel.tail(Path(row['log']), 5000).splitlines()[-3:]
        terminal = panel.D / seed / 'replay-M17/terminal.json'
        if terminal.exists():
            t = read(terminal)
            row['search_terminal'] = {k: t[k] for k in ('final_rank_lower_bound', 'charts', 'stop_reason')}
    print(json.dumps(s, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('launch', 'resume', 'status', 'controller', 'case'))
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--takeover', action='store_true')
    p.add_argument('--lock-fd', type=int)
    p.add_argument('--seed')
    a = p.parse_args()
    if a.action in ('launch', 'resume'):
        launch(a.workers, a.takeover)
    elif a.action == 'controller':
        require(a.lock_fd is not None, 'launch controller through the locked launcher')
        controller(a.lock_fd)
    elif a.action == 'case':
        run_case(a.seed)
    else:
        status()


if __name__ == '__main__':
    main()
