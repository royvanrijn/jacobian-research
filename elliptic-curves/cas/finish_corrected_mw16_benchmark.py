#!/usr/bin/env python3
"""Pin proven cache identities and supervise one finite corrected-score benchmark."""
import argparse
import sys
from pathlib import Path
import benchmark_corrected_mw16_annulus as bench
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

D = bench.D/'controller'


def cache_gate():
    p = cert.read(bench.D/'protocol.json')
    expected = {}
    for r in cert.read(bench.base.OUT)['encoding']['rows']:
        expected[str(bench.base.cache(r['family'], r['label']).relative_to(bench.ROOT))] = r['sha256']
    for r in cert.read(bench.local.OUT)['rows']:
        expected[r['path']] = r['sha256']
    if len(expected) != 15 or p['cache_bindings'] != expected:
        raise ArithmeticError('benchmark caches differ from exact upstream encoding proofs')
    if any(cert.hashed(bench.ROOT/n) != h for n, h in expected.items()):
        raise ArithmeticError('proven cache bytes changed')
    return expected


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve corrected benchmark controller')
    paths = [Path(__file__).resolve(), Path(bench.__file__), bench.D/'protocol.json',
             bench.base.OUT, bench.local.OUT]
    checkpoint(D/'protocol.json', {'sources': {str(p.relative_to(bench.ROOT)): cert.hashed(p) for p in paths},
        'cache_bindings': cache_gate(), 'stages': [['run', 1800], ['check', 1200]],
        'scope': 'Exactly one frozen corrected-score benchmark, with all fifteen base/local cache '
                 'hashes required to match the upstream exact byte proofs. Run and read-only replay '
                 'are separately supervised. No retries, replacement slices, broader scan or points.'})


def launch():
    p = cert.read(D/'protocol.json')
    if any(cert.hashed(bench.ROOT/n) != h for n, h in p['sources'].items()) or p['cache_bindings'] != cache_gate():
        raise ArithmeticError('frozen corrected benchmark gate changed')
    if (D/'ledger.json').exists():
        raise FileExistsError('preserve corrected benchmark launch')
    d = {'status': 'RUNNING', 'rows': []}
    checkpoint(D/'ledger.json', d)
    try:
        for stage, seconds in p['stages']:
            s = run([sys.executable, str(Path(bench.__file__)), stage], limits=Limits(seconds, 2147483648),
                    log_path=D/(stage+'.log'), checkpoint_path=D/(stage+'.supervisor.json'), cwd=bench.ROOT)
            ok = s['outcome'] == 'completed' and s['returncode'] == 0
            d['rows'].append({'name': stage, 'status': 'PASS' if ok else 'FAILED_OR_CENSORED', 'supervision': s})
            checkpoint(D/'ledger.json', d)
            print('CORRECTED MW16 BENCHMARK', stage, d['rows'][-1]['status'], flush=True)
            if not ok:
                raise ArithmeticError('finite benchmark failed/censored; no retry')
        result = cert.read(bench.OUT)
        d['status'] = 'PASS' if result['cost_gate_passed'] else 'PASS_REPLAY_COST_GATE_FAILED'
        d['result_sha256'] = cert.hashed(bench.OUT)
        checkpoint(D/'ledger.json', d)
    except Exception as exc:
        d.update(status='FAILED_OR_CENSORED', reason=str(exc))
        checkpoint(D/'ledger.json', d)
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['prepare', 'launch'])
    a = p.parse_args()
    globals()[a.stage]()
