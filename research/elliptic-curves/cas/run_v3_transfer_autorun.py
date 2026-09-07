#!/usr/bin/env python3
"""Detached fail-closed controller for the post-V3 transfer campaign.

One user-facing command:

    python3 elliptic-curves/cas/run_v3_transfer_autorun.py launch

Run it from research/.  ``launch`` starts one detached controller and returns
immediately.  The controller does not retune V3 or expand the frozen transfer
roster.  It waits for the already-running V3 detached supervisor to publish its
final replay/metric/package, prepares the fixed 11952 transfer once, and then
runs the gated four-case roster one case at a time through
``v3_transfer_campaign.sage next``.  Any preserved method/resource failure
stops the controller; there is no automatic retry or budget expansion.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]                 # research/
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
V3 = LOCAL / 'adaptive-visibility-cascade-v3'
TRANSFER = LOCAL / 'v3-transfer-11952-v1'
AUTO = LOCAL / 'v3-transfer-autorun-v1'
STATE = AUTO / 'state.json'
LOG = AUTO / 'autorun.log'
FINAL_REPLAY = V3 / 'replay-M17.json'
FINAL_METRIC = V3 / 'metric-replay-M17.json'
FINAL_PACKAGE = ART / 'adaptive_visibility_cascade_v3.json'
CAMPAIGN = CAS / 'v3_transfer_campaign.sage'
WAIT_SECONDS = 12 * 60 * 60
POLL_SECONDS = 15

sys.path.insert(0, str(CAS))
from v3_transfer_contract import CASE_SPECS, read_json  # noqa: E402


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.tmp-{os.getpid()}')
    temp.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    os.replace(temp, path)


def proc_token(pid: int):
    try:
        # Linux /proc field 22, robust to spaces/parentheses in comm.
        return Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return None


def live(pid: int, token: str | None) -> bool:
    return bool(token) and proc_token(pid) == token


def read_state():
    return json.loads(STATE.read_text()) if STATE.exists() else None


def update(status: str, **extra) -> None:
    old = read_state() or {}
    value = {**old, **extra, 'status': status, 'updated_at_unix': time.time()}
    atomic_json(STATE, value)
    print('AUTORUN', status, json.dumps(extra, sort_keys=True), flush=True)


def sage() -> str:
    candidates = [shutil.which('sage'), '/home/royvanrijn/.local/bin/sage',
                  '/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/sage']
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return str(Path(candidate).resolve())
    raise FileNotFoundError('sage executable not found; controller did not start transfer work')


def run_checked(argv) -> None:
    print('RUN', ' '.join(map(str, argv)), flush=True)
    subprocess.run(list(map(str, argv)), cwd=ROOT, check=True,
                   env={**os.environ, 'OPENBLAS_NUM_THREADS':'1',
                        'OMP_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1'})


def v3_ready() -> bool:
    if not (FINAL_REPLAY.exists() and FINAL_METRIC.exists() and FINAL_PACKAGE.exists()):
        return False
    replay = read_json(FINAL_REPLAY)
    if replay.get('rank_lower_bound') != 31:
        raise RuntimeError('final V3 replay exists but does not certify rank 31')
    package = read_json(FINAL_PACKAGE)
    # Do not depend on one historical packaging spelling; require an explicit
    # successful M17 target marker if present, otherwise the replay is the gate.
    for key in ('M17_target_achieved', 'm17_target_achieved'):
        if key in package and package[key] is not True:
            raise RuntimeError('V3 package explicitly says the M17 target was not achieved')
    return True


def wait_for_v3() -> None:
    deadline = time.monotonic() + WAIT_SECONDS
    update('WAITING_FOR_V3_FINAL_REPLAY', required=[str(FINAL_REPLAY.relative_to(ROOT)),
           str(FINAL_METRIC.relative_to(ROOT)), str(FINAL_PACKAGE.relative_to(ROOT))])
    while True:
        if v3_ready():
            update('V3_GATE_PASSED', v3_replay=str(FINAL_REPLAY.relative_to(ROOT)))
            return
        if time.monotonic() >= deadline:
            raise TimeoutError('V3 final replay/package did not appear within the fixed 12-hour controller wait')
        time.sleep(POLL_SECONDS)


def prepare_transfer(sage_bin: str) -> None:
    roster = TRANSFER / 'roster.json'
    if roster.exists():
        update('TRANSFER_ALREADY_PREPARED')
        return
    if TRANSFER.exists():
        raise RuntimeError('transfer directory exists without frozen roster; preserve and review partial preparation')
    update('PREPARING_TRANSFER')
    run_checked([sage_bin, '-python', CAMPAIGN, 'prepare', '--v3-replay', FINAL_REPLAY])
    if not roster.exists():
        raise RuntimeError('transfer preparation returned without frozen roster')
    update('TRANSFER_PREPARED')


def verified(case_id: str):
    path = TRANSFER / case_id / 'verified.json'
    return read_json(path) if path.exists() else None


def run_roster(sage_bin: str) -> None:
    for index, spec in enumerate(CASE_SPECS):
        report = verified(spec['id'])
        if report is None:
            update('RUNNING_CASE', case=spec['id'], case_index=index,
                   initial_rank=spec['initial_rank'])
            run_checked([sage_bin, '-python', CAMPAIGN, 'next'])
            report = verified(spec['id'])
            if report is None:
                raise RuntimeError(f"campaign next returned without verified result for {spec['id']}")
        update('CASE_VERIFIED', case=spec['id'], initial_rank=report['initial_rank'],
               rank_lower_bound=report['rank_lower_bound'], gain=report['gain'],
               stop_reason=report['stop_reason'])
        if spec['kind'] == 'positive-control' and report['gain'] <= 0:
            update('CONTROL_GATE_CLOSED', case=spec['id'],
                   reason='positive control produced no independently verified gain; warm jobs intentionally not released')
            return
    update('COMPLETE_FIXED_ROSTER', cases=len(CASE_SPECS),
           results={s['id']: verified(s['id'])['rank_lower_bound'] for s in CASE_SPECS})


def worker() -> None:
    signal.signal(signal.SIGHUP, signal.SIG_IGN)
    update('CONTROLLER_RUNNING', pid=os.getpid(), proc_token=proc_token(os.getpid()),
           controller_source=str(Path(__file__).relative_to(ROOT)))
    try:
        wait_for_v3()
        sage_bin = sage()
        update('SAGE_RESOLVED', sage=sage_bin)
        prepare_transfer(sage_bin)
        run_roster(sage_bin)
    except BaseException as exc:
        update('STOPPED_REVIEW_REQUIRED', error=repr(exc))
        raise


def launch() -> None:
    AUTO.mkdir(parents=True, exist_ok=True)
    old = read_state()
    if old:
        pid = int(old.get('pid', 0) or 0)
        if pid and live(pid, old.get('proc_token')):
            print(f"Already running pid={pid}; log={LOG.relative_to(ROOT)}")
            return
        if old.get('status') in ('COMPLETE_FIXED_ROSTER', 'CONTROL_GATE_CLOSED'):
            print(json.dumps(old, indent=2, sort_keys=True))
            return
        raise RuntimeError(f"previous autorun state is {old.get('status')!r}; preserve it and review before relaunch")
    stream = LOG.open('ab', buffering=0)
    proc = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), 'worker'],
        cwd=ROOT, stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
        start_new_session=True, close_fds=True,
        env={**os.environ, 'PYTHONUNBUFFERED':'1', 'OPENBLAS_NUM_THREADS':'1',
             'OMP_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1'})
    token = None
    for _ in range(20):
        token = proc_token(proc.pid)
        if token is not None:
            break
        time.sleep(0.01)
    atomic_json(STATE, {'status':'LAUNCHED', 'pid':proc.pid, 'proc_token':token,
        'launched_at_unix':time.time(), 'log':str(LOG.relative_to(ROOT))})
    print(f"Launched detached V3→transfer controller pid={proc.pid}")
    print(f"Status: python3 elliptic-curves/cas/run_v3_transfer_autorun.py status")
    print(f"Log: {LOG.relative_to(ROOT)}")


def status() -> None:
    value = read_state()
    if value is None:
        print('NOT_LAUNCHED')
        return
    pid = int(value.get('pid', 0) or 0)
    value['process_alive'] = bool(pid and live(pid, value.get('proc_token')))
    print(json.dumps(value, indent=2, sort_keys=True))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['launch', 'status', 'worker'])
    args = ap.parse_args()
    if args.action == 'launch':
        launch()
    elif args.action == 'status':
        status()
    else:
        worker()


if __name__ == '__main__':
    main()
