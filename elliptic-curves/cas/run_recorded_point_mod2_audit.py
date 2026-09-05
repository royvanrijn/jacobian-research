#!/usr/bin/env python3
"""Freeze and replay a bounded cohort of existing compact search point clouds."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import sys
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture, Limits

ROOT = Path(__file__).resolve().parents[2]
WORKER = ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py'

def run(directory, minimum, maximum):
    directory.mkdir()
    rows = []
    for folder in sorted((ROOT/'artifacts/local/elliptic-curves').iterdir()):
        if not folder.is_dir() or not any(s in folder.name for s in ('compact', 'prospective')):
            continue
        paths = set(folder.glob('*/result.json')) | set(folder.glob('*/*/result.json')) | set(folder.glob('*/*/*/result.json'))
        for path in sorted(paths):
            if 'previous' in path.parts:
                continue
            data = cert.read(path)
            if not minimum <= data.get('rank_lower_bound', -1) <= maximum or not data.get('charts') or not data.get('final_state'):
                continue
            basis = data['final_state']['state']['reductions']['points']
            seen = {(cert.F(p[0]), abs(cert.F(p[1]))) for p in basis}
            for chart in data['charts']:
                for p in chart['search']['finite_curve_points']:
                    seen.add((cert.F(p['x']), abs(cert.F(p['y']))))
            if len(seen) <= len(basis):
                continue
            rows.append({'index':len(rows), 'input':str(path.relative_to(ROOT)),
                'sha256':cert.hashed(path), 'original_rank_lower_bound':len(basis),
                'point_count':len(seen), 'input_status':data['status']})
    protocol = {'schema':'elliptic-curves.retained-mod2-cohort.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(), WORKER)},
        'maximum_workers':2, 'prime_bound':1000, 'build_seconds':300,
        'check_seconds':180, 'rss_bytes':1610612736, 'roster':rows,
        'rank_filter':[minimum, maximum],
        'scope':'All current retained compact/prospective transcripts in the rank filter with at least one additional point up to sign. Previous snapshots are omitted. Each input is hash-pinned. Lower bounds only; source searches can remain censored. No new search or public points.'}
    checkpoint(directory/'protocol.json', protocol)
    print('FROZEN POINT-CLOUD ROSTER', len(rows), flush=True)
    def audit(row):
        i = row['index']
        output = directory/f'{i:03}-certificate.json'
        stages = []
        try:
            if cert.hashed(WORKER) != protocol['sources'][str(WORKER.relative_to(ROOT))]:
                raise ArithmeticError('frozen worker source changed')
            commands = [('build', ['--input', str(ROOT/row['input']), '--input-sha256', row['sha256'], '--output', str(output)], 300),
                        ('check', ['--check', str(output)], 180)]
            for name, args, seconds in commands:
                r = capture([sys.executable, str(WORKER)]+args, limits=Limits(seconds, 1610612736), log_path=directory/f'{i:03}-{name}.log')
                stages.append({'stage':name, 'status':'PASS', 'supervision':r.supervision})
            proof = cert.read(output)
            result = {**row, 'status':'PASS', 'new_rank_lower_bound':proof['rank_lower_bound'],
                'output':str(output.relative_to(ROOT)), 'output_sha256':cert.hashed(output), 'stages':stages}
        except Exception as e:
            result = {**row, 'status':'FAILED', 'error':str(e), 'stages':stages}
        checkpoint(directory/f'{i:03}-summary.json', result)
        return result
    results = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(audit, row) for row in rows]):
            row = future.result()
            results.append(row)
            checkpoint(directory/'ledger.json', {'status':'RUNNING', 'rows':sorted(results, key=lambda r:r['index'])})
            print(row['index'], row['status'], row['original_rank_lower_bound'], '->', row.get('new_rank_lower_bound'), flush=True)
    checkpoint(directory/'ledger.json', {'status':'COMPLETE_DECLARED_AUDIT', 'rows':sorted(results, key=lambda r:r['index'])})

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--min-rank', type=int, required=True)
    p.add_argument('--max-rank', type=int, required=True)
    a = p.parse_args()
    run(a.directory.resolve(), a.min_rank, a.max_rank)
