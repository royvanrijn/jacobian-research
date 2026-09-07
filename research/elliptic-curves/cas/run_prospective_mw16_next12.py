#!/usr/bin/env python3
"""Run all twelve frozen MW16 follow-on candidates under fixed resource limits."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits
ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/prospective-mw16-f5-next12-v1'

def batch(directory):
    if (directory/'point-ledger.json').exists():
        raise FileExistsError('preserve previous batch')
    protocol = cert.read(directory/'point-protocol.json')
    for name, h in protocol['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('frozen worker source changed')
    family = protocol['families'][0]
    population = cert.read(directory/family/'population.json')
    rows = [{'index':i, 'family':family, 'parameter':r['parameter'], 'status':'PENDING'}
            for i,r in enumerate(population['finalists'])]
    if len(rows) != 12:
        raise ArithmeticError('fixed roster differs')
    checkpoint(directory/'point-launch.json', {'protocol_hash':digest(protocol),
        'batch_source_sha256':cert.hashed(Path(__file__).resolve()), 'roster':rows,
        'maximum_workers':4, 'scope':'All twelve addresses; no catalogue input or replacement.'})
    ledger = {'status':'RUNNING', 'rows':[dict(r) for r in rows]}
    checkpoint(directory/'point-ledger.json', ledger)
    def run(row):
        i = row['index']; path = directory/family/f'candidate-{i:02}'/'result.json'
        try:
            r = capture(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',
                str(ROOT/'elliptic-curves/cas/search_prospective_mw16_next12.sage'), 'run',
                '--directory', str(directory), '--index', str(i)],
                limits=Limits(300, 1610612736), log_path=directory/family/f'worker-{i:02}.log')
            data = cert.read(path)
            result = {**row, 'status':data['status'], 'supervision':r.supervision}
        except Exception as e:
            result = {**row, 'status':'FAILED_OR_CENSORED', 'error':str(e)}
        if path.exists():
            data = cert.read(path)
            result.update(rank_lower_bound=data['rank_lower_bound'], charts=len(data['charts']),
                result_path=str(path.relative_to(ROOT)), result_sha256=cert.hashed(path))
        return result
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(run, row) for row in rows]):
            row = future.result(); ledger['rows'][row['index']] = row
            checkpoint(directory/'point-ledger.json', ledger)
            print('NEXT12 MW16', row['parameter'], row['status'], row.get('rank_lower_bound'), row.get('charts'), flush=True)
    ledger['status'] = 'COMPLETE_FIXED_BATCH_ATTEMPTS'
    checkpoint(directory/'point-ledger.json', ledger)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=DIRECTORY)
    batch(p.parse_args().directory.resolve())
