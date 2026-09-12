#!/usr/bin/env python3
"""Supervise the single authorized two-field skew-sieve extension."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import signal
import subprocess
import time
import run_two_class_relation_pilot as pilot
import two_class_skew_core as skew

OUT = pilot.OUT.with_name('two_class_skew_sieve_v1')
WORKER = Path(__file__).with_name('two_class_skew_sieve.sage')


def worker(row, plan):
    key = row['curve_key']
    dest = OUT/key
    if (dest/'supervisor.json').exists():
        return pilot.read(dest/'supervisor.json')
    dest.mkdir(parents=True, exist_ok=True)
    assert not list(dest.iterdir()), 'unresolved partial attempt; do not overwrite'
    command = [plan['runtime']['sage'], '-python', str(WORKER), '--input', str(pilot.OUT/'inputs'/f'{key}.json'),
               '--setup', str(pilot.OUT/key/'setup.json'), '--output', str(dest)]
    start = time.monotonic()
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        log, _ = proc.communicate(timeout=180)
        status = 'PASS_BOUNDED_WORKER' if proc.returncode == 0 and (dest/'result.json').exists() else 'UNKNOWN_WORKER_FAILURE'
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        log, _ = proc.communicate()
        status = 'UNKNOWN_TIMEOUT'
    (dest/'worker.log').write_text(log)
    record = {'status': status, 'command': command, 'wall_limit_seconds': 180,
              'wall_seconds': round(time.monotonic()-start, 6), 'returncode': proc.returncode,
              'curve_key': key, 't': row['t']}
    pilot.save(dest/'supervisor.json', record)
    return record


def main():
    old = pilot.check(pilot.OUT, complete=True)
    sources = {str(p): pilot.sha(p) for p in [Path(__file__), WORKER, Path(skew.__file__)]}
    frozen_old = pilot.historical.tree_hashes(pilot.OUT)
    plan = {'schema': 'elliptic-curves.two-class-skew-supervisor.v1', 'sources': sources,
            'old_pilot_tree': frozen_old, 'runtime': old['runtime'], 'rows': old['rows'],
            'limits': {'wall_seconds_per_field': 180, 'max_jobs': 2, 'A_times_B': 262144},
            'selection': 'same retrospective two-field engineering control, not a prospective rank prediction',
            'class_group_and_point_searches': False}
    pilot.save(OUT/'plan.json', plan)
    rows = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = [pool.submit(worker, r, plan) for r in old['rows']]
        for job in as_completed(jobs):
            record = job.result()
            rows.append(record)
            print('SKEW_SUPERVISOR', record['t'], record['status'], record['wall_seconds'], flush=True)
    pilot.check(pilot.OUT, complete=True)
    assert pilot.historical.tree_hashes(pilot.OUT) == frozen_old
    pilot.save(OUT/'verified.json', {'status': 'PASS_BOUNDED_SUPERVISION_AND_OLD_EVIDENCE_INTEGRITY',
               'rows': sorted(rows, key=lambda r: r['curve_key']), 'original_pilot_and_arithmetic_bytes_unchanged': True,
               'bindings': {str(p.relative_to(OUT)): pilot.sha(p) for p in sorted(OUT.rglob('*')) if p.is_file()}})


if __name__ == '__main__':
    main()
