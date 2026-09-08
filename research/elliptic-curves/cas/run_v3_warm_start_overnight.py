#!/usr/bin/env python3
"""One maintained, Sage-free controller for the three authorized warm fibres.

  python3 elliptic-curves/cas/run_v3_warm_start_overnight.py resume --hours 10
  python3 elliptic-curves/cas/run_v3_warm_start_overnight.py status

The controller imports no Sage code. Exact search/replay services run ONLY as
`sage -python` subprocesses, under the existing wall/RSS supervisor. An advisory
lock prevents two controllers, PID start tokens detect stale/reused PIDs, and a
heartbeat reports SEARCHING versus REPLAYING. A sealed search is never rerun.
Old supervisor failures and all original checkpoints are left untouched.
Runbook: ../notes/V3_WARM_RUNBOOK.md.
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

from v3_warm_support import (WARM, atomic, bindings, process_info, read, require,
                             same_process, sha, tail, terminal_structure)

SELF = Path(__file__).resolve()
CAS = SELF.parent
ROOT = CAS.parents[1]
D = ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v3'
AUTO = ROOT/'artifacts/local/elliptic-curves/v3-warm-start-overnight-v1'
STATE, LOG, LOCK = AUTO/'state.json', AUTO/'autorun.log', AUTO/'controller.lock'
ENGINE = CAS/'v3_warm_engine.py'


def sources():
    return {str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in
            ('run_v3_warm_start_overnight.py', 'v3_warm_engine.py', 'v3_warm_replay.py', 'v3_warm_support.py')}


def sage_launcher():
    explicit = os.environ.get('V3_SAGE')
    candidates = [explicit] if explicit else [shutil.which('sage'), str(Path.home()/'.local/bin/sage'), '/usr/bin/sage']
    for candidate in candidates:
        if candidate:
            path = Path(candidate).expanduser()
            if path.is_file() and os.access(path, os.X_OK):
                return str(path.resolve())
    raise RuntimeError('Sage launcher not found. Set V3_SAGE=/absolute/path/to/sage')


def worker_command(sage, action, case, session):
    require(action in ('preflight', 'search', 'replay') and case in WARM, 'invalid worker action/case')
    return [str(sage), '-python', '-u', str(ENGINE), action, '--case', case, '--session', str(session)]


def validate_inputs():
    """Read-only, standard-library preflight; does not rebuild the parent bank."""
    control = read(D/'control-native11952/verified.json')
    require(control.get('status') == 'PASS_INDEPENDENT_TRANSFER_REPLAY' and
            control.get('initial_rank') == control.get('rank_lower_bound') == 17 and control.get('gain') == 0,
            'the completed generic-start no-gain control must remain intact')
    bindings(ROOT, control['bindings'])
    roster = read(D/'roster.json')
    input_hashes = {str((D/'roster.json').relative_to(ROOT)): sha(D/'roster.json')}
    for filename, key in (('parent-bank.json','parent_bank_sha256'), ('gate.json','gate_sha256'),
                           ('preparation.json','preparation_sha256')):
        require(sha(D/filename) == roster[key], f'changed {filename}')
        input_hashes[str((D/filename).relative_to(ROOT))] = roster[key]
    for case in WARM:
        entry = next(r for r in roster['cases'] if r['id'] == case)
        folder = D/case
        require(sha(folder/'protocol.json') == entry['protocol_sha256'], f'changed {case} protocol')
        policy = read(folder/'protocol.json')
        require(entry['initial_rank'] == policy['initial_rank'] == 27, 'warm initial rank changed')
        for key in ('inputs', 'sources', 'driver_sources'):
            bindings(ROOT, policy[key])
        input_hashes.update(policy['inputs'])
        input_hashes[str((folder/'protocol.json').relative_to(ROOT))] = entry['protocol_sha256']
    return input_hashes


def controller_alive(state):
    if state.get('start_token'):
        return same_process(state.get('pid'), state['start_token'])
    # Read legacy status safely, without treating a reused arbitrary PID as ours.
    row = process_info(state.get('pid'))
    if not row or row['state'] in ('Z', 'X'):
        return False
    try:
        argv = (Path('/proc')/str(row['pid'])/'cmdline').read_bytes().split(b'\0')
        return str(SELF).encode() in argv and b'worker' in argv
    except OSError:
        return False


def ensure_no_live_workers():
    paths = [p for case in WARM for p in (D/case).glob('warm-*.supervisor.json')]
    paths += list((AUTO/'sessions').glob('*/attempts/*/*/supervisor.json'))
    for path in paths:
        row = read(path)
        if same_process(row.get('pid'), row.get('start_token')):
            raise RuntimeError(f'supervised worker pid {row["pid"]} is still alive: {path}; do not launch concurrently')


def take_lock():
    AUTO.mkdir(parents=True, exist_ok=True)
    fd = os.open(LOCK, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return fd
    except BaseException:
        os.close(fd)
        raise RuntimeError('another warm controller owns the lock; inspect status, do not start a second one')


def launch(hours, *, resume=False):
    require(math.isfinite(hours) and 0 < hours <= 24, '--hours must be positive and at most 24')
    fd = take_lock()
    try:
        old = read(STATE) if STATE.exists() else {}
        require(not controller_alive(old), 'legacy controller is still alive; inspect it before stopping')
        ensure_no_live_workers()
        if old.get('status') == 'COMPLETE_WARM_ROSTER':
            print(json.dumps(old, indent=2))
            return
        require(resume or not old, 'existing attempt: use resume; nothing will be deleted')
        inputs = validate_inputs()
        sage = sage_launcher()
        session_dir = AUTO/'sessions'/str(uuid.uuid4())
        session_dir.mkdir(parents=True)
        session_path = session_dir/'session.json'
        session = {'schema': 'warm-controller-session.v2', 'sources': sources(), 'inputs': inputs,
                   'cases': list(WARM), 'sage': sage, 'hours': hours, 'created_unix': time.time(),
                   'action': 'resume' if resume else 'launch',
                   'scope': 'Original three warm fibres and frozen search budgets; no automatic retuning or new parameters.'}
        atomic(session_path, session, immutable=True)
        if old:
            atomic(session_dir/'previous-controller-state.json', old, immutable=True)
        atomic(STATE, {'status': 'LAUNCHING', 'session': str(session_path.relative_to(ROOT)),
                       'updated_unix': time.time()})
        with LOG.open('ab', buffering=0) as log:
            proc = subprocess.Popen([sys.executable, str(SELF), 'worker', '--session', str(session_path),
                                     '--lease-fd', str(fd)], cwd=ROOT, stdin=subprocess.DEVNULL,
                                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                    pass_fds=(fd,), env={**os.environ, 'PYTHONUNBUFFERED': '1'})
        # The inherited descriptor keeps the SAME flock alive. Do not LOCK_UN
        # in the parent, and do not overwrite the child's more advanced status.
        print(f'Launched controller pid={proc.pid}; session={session_path.relative_to(ROOT)}')
        print('Status: python3 elliptic-curves/cas/run_v3_warm_start_overnight.py status')
    finally:
        os.close(fd)


class Reporter:
    def __init__(self, session):
        self.lock = threading.RLock()
        info = process_info(os.getpid())
        self.value = {'pid': os.getpid(), 'start_token': info['start_token'],
                      'session': str(session.relative_to(ROOT))}
        self.done = threading.Event()
        self.thread = threading.Thread(target=self._heartbeat, daemon=True)

    def update(self, status=None, **values):
        with self.lock:
            self.value.update(values)
            if status is not None:
                self.value['status'] = status
            self.value['updated_unix'] = time.time()
            atomic(STATE, self.value)
            if status is not None:
                print('WARM', status, json.dumps(values, sort_keys=True), flush=True)

    def _heartbeat(self):
        while not self.done.wait(10):
            try:
                self.update()
            except Exception as exc:
                print('HEARTBEAT_WRITE_FAILED', repr(exc), flush=True)

    def __enter__(self):
        self.update('CONTROLLER_RUNNING')
        self.thread.start()
        return self

    def __exit__(self, *_):
        self.done.set()
        self.thread.join(timeout=2)


class ResourceStop(RuntimeError):
    pass


def verify_result(case):
    path = D/case/'warm-verified.json'
    if not path.exists():
        return None
    row = read(path)
    require(row.get('status') == 'PASS_INDEPENDENT_WARM_REPLAY' and row.get('case') == case,
            'unrecognized warm result; do not silently accept a legacy/partial certificate')
    bindings(ROOT, row['bindings'])
    bindings(ROOT, row['sources'])
    return row


def run_step(action, case, session_path, session, deadline, report):
    from research_runtime.supervisor import Limits, run
    from v3_transfer_contract import LIMITS
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise ResourceStop('campaign wall budget exhausted')
    maximum = 120 if action == 'preflight' else LIMITS['case_wall_seconds' if action == 'search' else 'replay_wall_seconds']
    attempt = session_path.parent/'attempts'/case/(action+'-'+uuid.uuid4().hex[:12])
    attempt.mkdir(parents=True)
    command = worker_command(session['sage'], action, case, session_path)
    report.update({'preflight': 'PREFLIGHT', 'search': 'SEARCHING', 'replay': 'REPLAYING'}[action],
                  case=case, worker_action=action, attempt=str(attempt.relative_to(ROOT)),
                  worker_log=str((attempt/'worker.log').relative_to(ROOT)))
    result = run(command, limits=Limits(min(maximum, remaining), LIMITS['rss_bytes']), cwd=ROOT,
                 env={**os.environ, 'V3_WARM_SUPERVISED': '1', 'OPENBLAS_NUM_THREADS': '1',
                      'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'},
                 log_path=attempt/'worker.log', checkpoint_path=attempt/'supervisor.json')
    if result['outcome'] != 'completed':
        # A committed terminal is usable INPUT to independent replay, not proof
        # that an anomalous process exit succeeded or that its rank is exact.
        if action == 'search' and terminal_structure(D/case) is not None:
            report.update('SEALED_SEARCH_PENDING_REPLAY', case=case, original_outcome=result['outcome'])
            return
        if result['outcome'] in ('strict_wall_timeout', 'strict_rss_limit'):
            raise ResourceStop(f'{case} {action}: {result["outcome"]}; checkpoints retained')
        raise RuntimeError(f'{case} {action}: {result["outcome"]}\n' + (tail(attempt/'worker.log') or ''))
    if action == 'search':
        require(terminal_structure(D/case) is not None, 'worker exited without a sealed terminal')
    if action == 'replay':
        require(verify_result(case) is not None, 'replay exited without independent certificate')


def worker(session_path, lease_fd):
    require(session_path is not None and lease_fd is not None, 'worker must be launched by controller')
    # An inherited flock is the lifetime lease, independent of mutable PID files.
    require(os.fstat(lease_fd).st_ino == LOCK.stat().st_ino, 'wrong inherited controller lease')
    signal.signal(signal.SIGHUP, signal.SIG_IGN)

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'controller received signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    session = read(session_path)
    deadline = time.monotonic() + session['hours']*3600
    try:
        with Reporter(session_path) as report:
            try:
                require(session['sources'] == sources(), 'controller sources changed after launch')
                bindings(ROOT, session['sources'])
                bindings(ROOT, session['inputs'])
                resource_stops = []
                results = []
                for case in WARM:
                    old = verify_result(case)
                    if old is not None:
                        results.append(old)
                        continue
                    try:
                        run_step('preflight', case, session_path, session, deadline, report)
                        if (D/case/'replay-M17/terminal.json').exists():
                            terminal = terminal_structure(D/case)
                            report.update('SEALED_SEARCH_REUSED', case=case, charts=terminal['charts'])
                        else:
                            run_step('search', case, session_path, session, deadline, report)
                        run_step('replay', case, session_path, session, deadline, report)
                        result = verify_result(case)
                        results.append(result)
                        report.update('CASE_VERIFIED', case=case, rank_lower_bound=result['rank_lower_bound'],
                                      gain=result['gain'])
                    except ResourceStop as exc:
                        resource_stops.append({'case': case, 'reason': str(exc)})
                        report.update('CASE_RESOURCE_STOP', case=case, error=str(exc))
                        if time.monotonic() >= deadline:
                            break
                        # Resource failures are case-local and not negative rank
                        # results. Proof/replay/integrity errors stop the roster.
                complete = len(results) == len(WARM)
                report.update('COMPLETE_WARM_ROSTER' if complete else 'STOPPED_RESOURCE_BUDGET',
                              results=results, resource_stops=resource_stops)
            except BaseException as exc:
                report.update('STOPPED_REVIEW_REQUIRED', error=repr(exc))
                raise
    finally:
        os.close(lease_fd)


def status():
    row = read(STATE) if STATE.exists() else {'status': 'NOT_LAUNCHED'}
    row = dict(row)
    row['process_alive'] = controller_alive(row)
    info = process_info(row.get('pid'))
    row['process_state'] = info['state'] if info else None
    if row['status'] in ('CONTROLLER_RUNNING','LAUNCHING','PREFLIGHT','SEARCHING','REPLAYING','RUNNING_WARM_CASE',
                          'LAUNCHED','RELAUNCHED','SEALED_SEARCH_PENDING_REPLAY','SEALED_SEARCH_REUSED') and not row['process_alive']:
        row['effective_status'] = 'STOPPED_REVIEW_REQUIRED'
    if row.get('attempt'):
        path = ROOT/row['attempt']/'supervisor.json'
        if path.exists():
            row['supervisor'] = read(path)
    if 'results' in row:
        row['results'] = [{k:r.get(k) for k in ('case','rank_lower_bound','gain','charts','stop_reason')} for r in row['results']]
    if row.get('error'):
        row['error'] = row['error'][-4000:]
    if row.get('worker_log'):
        recent = tail(ROOT/row['worker_log'], 8192)
        row['progress_tail'] = recent.splitlines()[-5:] if recent else []
    print(json.dumps(row, indent=2, sort_keys=True))
    for case in WARM:
        folder = D/case
        if (folder/'warm-verified.json').exists():
            result = read(folder/'warm-verified.json')
            print(case, 'VERIFIED', json.dumps({k:result.get(k) for k in ('rank_lower_bound','gain','charts','stop_reason')}, sort_keys=True))
        elif (folder/'replay-M17/terminal.json').exists():
            t = read(folder/'replay-M17/terminal.json')
            print(case, 'SEARCH_SEALED_REPLAY_PENDING', json.dumps({k:t[k] for k in ('charts','final_rank_lower_bound','stop_reason')}))
        else:
            stages = folder/'replay-M17/stages.json'
            print(case, 'SEARCH_INCOMPLETE', json.dumps(read(stages)[-1]) if stages.exists() and read(stages) else '')


def diagnose():
    status()
    row = read(STATE) if STATE.exists() else {}
    log = ROOT/row['worker_log'] if row.get('worker_log') else None
    if log is None:
        # Legacy failed replay is retained exactly where it was written.
        files = [p for case in WARM for p in (D/case).glob('warm-*.log')]
        log = max(files, key=lambda p:p.stat().st_mtime) if files else LOG
    print('\nWORKER_LOG', log)
    print(tail(log) or 'No worker log.')


def stop():
    row = read(STATE) if STATE.exists() else {}
    require(controller_alive(row), 'no live matching controller; nothing was killed')
    require(row.get('start_token') and hasattr(os, 'pidfd_open') and hasattr(signal, 'pidfd_send_signal'),
            'safe PID-handle stop unavailable for legacy controller; inspect its identity manually')
    fd = os.pidfd_open(int(row['pid']))
    try:
        require(same_process(row['pid'], row['start_token']), 'PID identity changed; nothing was killed')
        signal.pidfd_send_signal(fd, signal.SIGTERM)
    finally:
        os.close(fd)
    print('Sent SIGTERM to the verified controller; supervised child cleanup preserves checkpoints.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['launch', 'resume', 'status', 'diagnose', 'stop', 'worker'])
    parser.add_argument('--hours', type=float, default=10.0, help='total wall budget, not a completion estimate')
    parser.add_argument('--session', type=Path)
    parser.add_argument('--lease-fd', type=int)
    args = parser.parse_args()
    if args.action in ('launch','resume'):
        launch(args.hours, resume=args.action == 'resume')
    elif args.action == 'status':
        status()
    elif args.action == 'diagnose':
        diagnose()
    elif args.action == 'stop':
        stop()
    else:
        worker(args.session, args.lease_fd)


if __name__ == '__main__':
    main()
