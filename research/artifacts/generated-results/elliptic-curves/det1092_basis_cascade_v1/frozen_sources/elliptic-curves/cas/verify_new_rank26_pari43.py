#!/usr/bin/env python3
"""Verify the generic-only native11952 pari49 control."""
from pathlib import Path
import sys
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

ROOT = Path(__file__).resolve().parents[2]
D = ROOT/'artifacts/local/elliptic-curves/new-rank26-pari43-v1'
ART = ROOT/'artifacts/generated-results/elliptic-curves'


def verify():
    terminal = cert.read(D/'terminal.json')
    path = D/'candidate-00/result.json'
    if cert.hashed(path) != terminal['result_sha256']:
        raise ArithmeticError('terminal input changed')
    folder = D/'verification'
    if (folder/'protocol.json').exists():
        raise FileExistsError('preserve follow-up verification attempts')
    output = ART/'new_rank26_pari43_recorded_mod2_v1.json'
    checkpoint(folder/'protocol.json', {'input_sha256': cert.hashed(path),
        'terminal_sha256': cert.hashed(D/'terminal.json'),
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
            Path(__file__).resolve(), ROOT/'elliptic-curves/cas/replay_new_rank26_pari43.py',
            ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py')},
        'chart_replay_wall_seconds': 900, 'cloud_build_wall_seconds': 300,
        'cloud_check_wall_seconds': 180, 'rss_bytes': 1610612736,
        'scope': 'All retained charts and the complete raw point cloud, including any censored prefix. A replayed finite miss is not an upper rank bound.'})
    stages = [('charts', ['replay_new_rank26_pari43.py', str(path)], 900),
        ('cloud-build', ['audit_recorded_point_mod2_rank_v2.py', '--input', str(path),
                         '--input-sha256', cert.hashed(path), '--output', str(output)], 300),
        ('cloud-check', ['audit_recorded_point_mod2_rank_v2.py', '--check', str(output)], 180)]
    result = {'status': 'RUNNING', 'rows': []}
    checkpoint(folder/'result.json', result)
    for stage, args, seconds in stages:
        command = [sys.executable, str(ROOT/'elliptic-curves/cas'/args[0]), *args[1:]]
        r = run(command, limits=Limits(seconds, 1610612736), log_path=folder/(stage+'.log'),
                checkpoint_path=folder/(stage+'.supervisor.json'), cwd=ROOT)
        ok = r['outcome'] == 'completed' and r['returncode'] == 0
        result['rows'].append({'stage': stage, 'status': 'PASS' if ok else 'FAILED_OR_CENSORED', 'supervision': r})
        checkpoint(folder/'result.json', result)
        print('NEW RANK26 PARI43 VERIFICATION', stage, result['rows'][-1]['status'], flush=True)
        if stage == 'cloud-build' and not ok:
            break
    result['status'] = 'COMPLETE_DECLARED_REPLAY_ATTEMPTS'
    checkpoint(folder/'result.json', result)


if __name__ == '__main__':
    verify()
