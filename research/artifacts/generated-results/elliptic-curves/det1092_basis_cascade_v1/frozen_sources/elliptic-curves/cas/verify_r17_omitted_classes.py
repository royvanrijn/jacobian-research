#!/usr/bin/env python3
"""Replay all eight omitted-class attempts and all new recorded point clouds."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import sys
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/r17-omitted-generic-classes-v1'

def verify(directory):
    ledger = cert.read(directory/'ledger.json')
    if ledger['status'] != 'COMPLETE_DECLARED_ATTEMPTS' or len(ledger['rows']) != 8:
        raise ArithmeticError('fixed point batch unfinished')
    folder = directory/'verification'
    if (folder/'ledger.json').exists():
        raise FileExistsError('preserve verification attempts')
    checkpoint(folder/'protocol.json', {'point_ledger_sha256': cert.hashed(directory/'ledger.json'),
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (Path(__file__).resolve(),
            ROOT/'elliptic-curves/cas/replay_r17_omitted_classes.py', ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py')},
        'workers': 2, 'chart_replay_wall_seconds': 600, 'cloud_build_wall_seconds': 300,
        'cloud_check_wall_seconds': 180, 'rss_bytes': 1610612736,
        'scope': 'Every retained terminal prefix and complete recorded cloud, including censored searches. No catalogue access. Missing measurements and verification failures remain explicit.'})
    result = {'status': 'RUNNING', 'rows': []}
    checkpoint(folder/'ledger.json', result)
    def worker(row):
        if 'result_path' not in row:
            return {'family': row['family'], 'index': row['index'], 'status': 'NO_RETAINED_MEASUREMENT', 'input_row': row}
        path = ROOT/row['result_path']
        if cert.hashed(path) != row['result_sha256']:
            raise ArithmeticError('terminal worker checkpoint changed')
        label = row['family']+f"-{row['index']:02}"
        artifact = ART/f"r17_omitted_classes_recorded_mod2_{row['family']}_{row['index']:02}_v1.json"
        stages = [('charts', [sys.executable, str(ROOT/'elliptic-curves/cas/replay_r17_omitted_classes.py'), str(path)], 600),
            ('cloud-build', [sys.executable, str(ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py'),
                '--input', str(path), '--input-sha256', row['result_sha256'], '--output', str(artifact)], 300),
            ('cloud-check', [sys.executable, str(ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py'), '--check', str(artifact)], 180)]
        outcomes = []
        for stage, argv, seconds in stages:
            s = run(argv, limits=Limits(seconds, 1610612736), log_path=folder/(label+'-'+stage+'.log'),
                    checkpoint_path=folder/(label+'-'+stage+'.supervisor.json'), cwd=ROOT)
            outcomes.append({'stage': stage, 'status': 'PASS' if s['outcome'] == 'completed' and s['returncode'] == 0 else 'FAILED_OR_CENSORED', 'supervision': s})
            # A chart-replay failure does not prevent independent point proofs.
            if stage == 'cloud-build' and outcomes[-1]['status'] != 'PASS':
                break
        return {'family': row['family'], 'index': row['index'], 'parameter': row['parameter'],
                'rows': outcomes, 'status': 'PASS' if len(outcomes) == 3 and all(q['status'] == 'PASS' for q in outcomes) else 'FAILED_OR_CENSORED'}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(worker, row) for row in ledger['rows']]):
            row = future.result(); result['rows'].append(row)
            checkpoint(folder/'ledger.json', result)
            print('VERIFIED OMITTED R17 ATTEMPT', row['family'], row['index'], row['status'], flush=True)
    result['status'] = 'COMPLETE_DECLARED_VERIFICATION_ATTEMPTS'
    checkpoint(folder/'ledger.json', result)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=DIRECTORY)
    verify(p.parse_args().directory.resolve())
