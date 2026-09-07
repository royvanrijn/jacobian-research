#!/usr/bin/env python3
"""Exact algebraic orbit audit on the already replayed initial43 new-rank26 charts."""
from pathlib import Path
import time
import certify_compact_r17_candidates as cert
from research_runtime.pointed_orbit_compression import compress
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1/07ca9/candidate-00/result.json'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/pointed_orbit_compression_audit_v1.json'


def audit():
    if OUT.exists():
        raise FileExistsError('preserve orbit audit')
    original = cert.read(INPUT)
    rows = []
    started = time.monotonic()
    for i, chart in enumerate(original['charts']):
        r = chart['search']; model = r['input']['curve']
        decode = lambda p: (p['x'], p['y'])
        subgroup = [decode(p) for p in r['input']['subgroup']]
        points = [decode(p) for p in r['finite_curve_points']]
        result = compress(model, subgroup, r['input']['centre']['coefficients'], points)
        if result['centre'] != list(decode(r['base_point'])):
            raise ArithmeticError('recomputed centre differs')
        rows.append({'chart': i, 'raw_points': len(points), **result})
    data = {'schema': 'elliptic-curves.pointed-orbit-compression-audit.v1',
        'status': 'PASS', 'input_path': str(INPUT.relative_to(ROOT)), 'input_sha256': cert.hashed(INPUT),
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
            Path(__file__).resolve(), ROOT/'elliptic-curves/cas/research_runtime/pointed_orbit_compression.py',
            ROOT/'elliptic-curves/cas/alternate_quartic_covers.py',
            ROOT/'elliptic-curves/cas/half_lattice_pointed_sieve.py',
            ROOT/'elliptic-curves/tests/test_pointed_orbit_compression.py')},
        'charts': rows, 'raw_admission_events': sum(r['raw_points'] for r in rows),
        'kept_admission_events': sum(len(r['kept_indices']) for r in rows),
        'verified_skip_relations': sum(len(r['skipped']) for r in rows),
        'wall_seconds': time.monotonic()-started,
        'claim_boundary': 'Every omitted point is exactly centre minus a kept point, and every centre is an integral word in its chart input subgroup. Thus compression preserves the generated integral subgroup. This is an algebraic audit of retained charts, not a new full search replay, timing speedup measurement, or proof that finite-admission records are identical. No running or frozen search worker uses this optional helper.'}
    checkpoint(OUT, data)
    print('EXACT ORBIT AUDIT', len(rows), 'charts;', data['raw_admission_events'], 'events ->', data['kept_admission_events'], flush=True)


if __name__ == '__main__':
    audit()
