#!/usr/bin/env python3
"""Supervise the frozen remaining-five R17 cache build and read-only replay once."""
import argparse
import sys
from pathlib import Path
import build_r17_remaining_extended_projective_caches as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

D = batch.D / 'controller'


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve cache controller')
    batch.protocol()
    paths = [Path(__file__).resolve(), Path(batch.__file__).resolve(), batch.D / 'protocol.json',
             batch.CAS / 'research_runtime/store.py', batch.CAS / 'research_runtime/supervisor.py']
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.r17-remaining-cache-controller.v1',
        'sources': {str(p.relative_to(batch.ROOT)): cert.hashed(p) for p in paths},
        'jobs': [{'name': 'run', 'seconds': 7200}, {'name': 'check', 'seconds': 3600}],
        'rss_bytes': 4294967296,
        'scope': 'Supervise the fixed14740-table build and full replay for the five R17 families other than11952. Exactly15 benchmark tables reused once; at most14725 new GP calls with five workers,20seconds per call and80-case checkpoints. No retries, automatic parameter scan, score selection or point search.'})


def launch():
    p = cert.read(D / 'protocol.json')
    if any(cert.hashed(batch.ROOT / n) != h for n, h in p['sources'].items()):
        raise ArithmeticError('controller source binding changed')
    if (D / 'ledger.json').exists():
        raise FileExistsError('preserve cache controller attempt')
    ledger = {'status': 'RUNNING', 'rows': []}
    checkpoint(D / 'ledger.json', ledger)
    for job in p['jobs']:
        name = job['name']
        s = run([sys.executable, str(Path(batch.__file__).resolve()), name],
                limits=Limits(job['seconds'], p['rss_bytes']), log_path=D / (name + '.log'),
                checkpoint_path=D / (name + '.supervisor.json'), cwd=batch.ROOT)
        ok = s['outcome'] == 'completed' and s['returncode'] == 0
        ledger['rows'].append({'name': name, 'status': 'PASS' if ok else 'FAILED_OR_CENSORED', 'supervision': s})
        checkpoint(D / 'ledger.json', ledger)
        print('REMAINING FIVE R17 CACHES', name, ledger['rows'][-1]['status'], flush=True)
        if not ok:
            ledger['status'] = 'FAILED_OR_CENSORED'
            checkpoint(D / 'ledger.json', ledger)
            raise ArithmeticError('cache stage failed or censored; no retry')
    ledger['status'] = 'PASS'
    checkpoint(D / 'ledger.json', ledger)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=['prepare', 'launch'])
    args = p.parse_args(); globals()[args.stage]()
