#!/usr/bin/env python3
"""Supervise the frozen remaining-five R17 cache build and read-only replay once."""
import argparse
import sys
import time
from pathlib import Path
import encode_r17_remaining_joint_caches as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

D = batch.D / 'controller'


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve cache controller')
    batch.parent.protocol()
    paths = [Path(__file__).resolve(), Path(batch.__file__).resolve(), batch.parent.D / 'protocol.json',
             batch.CAS / 'research_runtime/store.py', batch.CAS / 'research_runtime/supervisor.py']
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.r17-remaining-cache-controller.v1',
        'sources': {str(p.relative_to(batch.ROOT)): cert.hashed(p) for p in paths},
        'wait_seconds': 15000,
        'jobs': [{'name': 'prepare', 'seconds': 120}, {'name': 'encode', 'seconds': 2400}, {'name': 'encoding-check', 'seconds': 2400}, {'name': 'score', 'seconds': 600}, {'name': 'score-check', 'seconds': 600}],
        'rss_bytes': 4294967296,
        'scope': 'Wait for complete full-cache construction and replay, then encode and replay all ten short/extended caches. Require exact short/extended score and good-count agreement on all4831 earlier outer scalar candidates,9662 component comparisons. No new point counting, parameter search, retries or point exposure.'})


def launch():
    p = cert.read(D / 'protocol.json')
    if any(cert.hashed(batch.ROOT / n) != h for n, h in p['sources'].items()):
        raise ArithmeticError('controller source binding changed')
    if (D / 'ledger.json').exists():
        raise FileExistsError('preserve cache controller attempt')
    ledger = {'status': 'WAITING_FOR_FULL_CACHES', 'rows': []}
    checkpoint(D / 'ledger.json', ledger)
    deadline = time.monotonic() + p['wait_seconds']
    while True:
        state = cert.read(batch.parent.D/'controller/ledger.json')['status']
        if state == 'PASS': break
        if state != 'RUNNING': raise ArithmeticError('upstream cache failed or censored')
        if time.monotonic() > deadline: raise TimeoutError('declared cache wait elapsed')
        time.sleep(5)
    ledger['status'] = 'RUNNING'; checkpoint(D/'ledger.json',ledger)
    for job in p['jobs']:
        name = job['name']
        s = run([sys.executable, str(Path(batch.__file__).resolve()), name],
                limits=Limits(job['seconds'], p['rss_bytes']), log_path=D / (name + '.log'),
                checkpoint_path=D / (name + '.supervisor.json'), cwd=batch.ROOT)
        ok = s['outcome'] == 'completed' and s['returncode'] == 0
        ledger['rows'].append({'name': name, 'status': 'PASS' if ok else 'FAILED_OR_CENSORED', 'supervision': s})
        checkpoint(D / 'ledger.json', ledger)
        print('REMAINING FIVE R17 BINARY CACHES', name, ledger['rows'][-1]['status'], flush=True)
        if not ok:
            ledger['status'] = 'FAILED_OR_CENSORED'
            checkpoint(D / 'ledger.json', ledger)
            raise ArithmeticError('cache stage failed or censored; no retry')
    ledger['status'] = 'PASS'
    checkpoint(D / 'ledger.json', ledger)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=['prepare', 'launch'])
    args = p.parse_args(); globals()[args.stage]()
