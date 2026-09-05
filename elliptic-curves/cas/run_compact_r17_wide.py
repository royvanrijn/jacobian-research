#!/usr/bin/env python3
"""Run the fixed24-address compact R17 batch with four bounded workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import run, Limits
ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1'

def batch(directory):
    if (directory/'point-ledger.json').exists():
        raise FileExistsError('preserve fixed point batch')
    protocol = cert.read(directory/'point-protocol.json')
    for name, h in protocol['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('frozen point source differs')
    populations = {f: cert.read(directory/f/'population.json')['finalists'] for f in protocol['families']}
    rows = [{'family': f, 'index': i, 'parameter': populations[f][i]['parameter'], 'status': 'PENDING'}
            for i in range(4) for f in protocol['families']]
    if len(rows) != 24:
        raise ArithmeticError('fixed roster differs')
    checkpoint(directory/'point-launch.json', {'protocol_hash': digest(protocol),
        'runner_sha256': cert.hashed(Path(__file__).resolve()), 'rows': rows,
        'maximum_workers': 4, 'scope': 'Every address receives a new attempt. No catalogue input or replacement.'})
    ledger = {'status': 'RUNNING', 'rows': [dict(r) for r in rows]}
    checkpoint(directory/'point-ledger.json', ledger)
    def worker(row):
        folder = directory/row['family']/f"candidate-{row['index']:02}"
        path = folder/'result.json'
        s = run(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',
                 str(ROOT/'elliptic-curves/cas/search_compact_r17_wide.sage'), 'run',
                 '--directory', str(directory), '--family', row['family'], '--index', str(row['index'])],
                limits=Limits(300, 1610612736), log_path=folder/'worker.log',
                checkpoint_path=folder/'worker.supervisor.json', cwd=ROOT)
        result = {**row, 'status': 'FAILED_OR_CENSORED', 'supervision': s}
        if path.exists():
            data = cert.read(path)
            result.update(result_path=str(path.relative_to(ROOT)), result_sha256=cert.hashed(path),
                          rank_lower_bound=data['rank_lower_bound'], charts=len(data['charts']))
            if s['outcome'] == 'completed' and s['returncode'] == 0:
                result['status'] = data['status']
        return result
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(worker, row): i for i, row in enumerate(rows)}
        for future in as_completed(futures):
            row = future.result(); ledger['rows'][futures[future]] = row
            checkpoint(directory/'point-ledger.json', ledger)
            print('R17 H4096 POINT ATTEMPT', row['family'], row['parameter'], row['status'],
                  row.get('rank_lower_bound'), row.get('charts'), flush=True)
    ledger['status'] = 'COMPLETE_FIXED_BATCH_ATTEMPTS'
    checkpoint(directory/'point-ledger.json', ledger)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=DIRECTORY)
    batch(p.parse_args().directory.resolve())
