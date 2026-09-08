#!/usr/bin/env python3
"""Independently certify every fresh initial measurement in the two full-score cohorts.

Scores and coverage are retained scheduling/worker metadata, not rank proofs.
The portable check needs neither Sage nor the discovery caches.
"""
import argparse
import json
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
COHORTS = ('compact-r17-top64-v1', 'compact-r17-h16384-v1')


def build(output):
    if output.exists():
        raise FileExistsError('use a new immutable output')
    rows, bindings = [], {}
    for cohort in COHORTS:
        directory = ROOT / 'artifacts/local/elliptic-curves' / cohort
        ledger, population = (cert.read(directory / name) for name in ('ledger.json', 'population.json'))
        if ledger['status'] != 'COMPLETE_DECLARED_BATCH':
            raise ArithmeticError('initial cohort unfinished')
        for name in ('ledger.json', 'population.json', 'protocol.json'):
            path = directory / name
            bindings[str(path.relative_to(ROOT))] = cert.hashed(path)
        for item in ledger['rows']:
            if item['status'] != 'COMPLETE':
                continue  # Reused measurements and known equations are reported separately.
            path = ROOT / item['path']; result = cert.read(path)
            model = tuple(map(cert.F, result['curve']))
            points = tuple(tuple(map(cert.F, p)) for p in result['final_state']['state']['reductions']['points'])
            proof = cert.checked_rank(model, points)
            scale = cert.family_check(result['parameter'], model, points)
            if len(points) != item['rank_lower_bound'] or result['status'] != 'COMPLETE':
                raise ArithmeticError('initial rank or completion assertion changed')
            candidate = population['finalists'][item['index']]
            if candidate['parameter'] != result['parameter']:
                raise ArithmeticError('population address changed')
            charts = [r['search'] for r in result['charts']]
            rows.append({'cohort': cohort, 'index': item['index'], 'parameter': result['parameter'],
                         'curve': result['curve'], 'points': [list(map(str, p)) for p in points],
                         'rank_certificate': proof, 'family_to_curve_scale_u': scale,
                         'score_units': candidate['score_units'],
                         'short_coefficient_bits': max(max(abs(q.numerator).bit_length(), q.denominator.bit_length()) for q in model),
                         'chart_count': len(charts),
                         'full_box_count': sum(r['completed_denominator'] >= r['denominator_end'] for r in charts),
                         'mean_denominator_coverage': sum(r['completed_denominator'] / r['denominator_end'] for r in charts) / len(charts),
                         'input': {'path': item['path'], 'sha256': cert.hashed(path)}})
            print('MEASUREMENT', cohort, result['parameter'], len(points), flush=True)
    sources = [Path(__file__).resolve(), Path(cert.__file__).resolve(), cert.MODEL, cert.SECTIONS,
               ROOT / 'elliptic-curves/cas/mod2_reduction_independence.py',
               ROOT / 'elliptic-curves/cas/elliptic_candidate_record.py']
    cert.write(output, {'schema': 'elliptic-curves.compact-r17-initial-measurements.v1',
                       'rows': rows, 'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in sources},
                       'discovery_bindings': bindings,
                       'claim_boundary': 'Exact lower bounds only. Scores and coverage are pinned discovery metadata, not independently rerun enumeration. Fresh means not reused within these cohorts; novelty is checked separately. No true-rank distribution or score enrichment claim.'})


def check(path):
    data = cert.read(path)
    for name, expected in data['sources'].items():
        if cert.hashed(ROOT / name) != expected:
            raise ArithmeticError('checker dependency changed: ' + name)
    counts = {name: Counter() for name in COHORTS}
    for row in data['rows']:
        model = tuple(map(cert.F, row['curve'])); points = tuple(tuple(map(cert.F, p)) for p in row['points'])
        proof = row['rank_certificate']
        actual = cert.checked_rank(model, points, [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
        if json.dumps(actual, sort_keys=True) != json.dumps(proof, sort_keys=True) or cert.family_check(row['parameter'], model, points) != row['family_to_curve_scale_u']:
            raise ArithmeticError('point proof or transport changed')
        counts[row['cohort']][len(points)] += 1
    print('REPLAYED INITIAL LOWER BOUNDS', {name: dict(sorted(count.items())) for name, count in counts.items()}, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else build(args.output)
