#!/usr/bin/env python3
"""Exact rank-increment targets against the retained conductor benchmark snapshot.

No factorization, point search, current-record assertion, or probability model.
"""
import argparse
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'artifacts/generated-results/elliptic-curves/new_curve_conductor_screen_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    if output.exists():
        raise FileExistsError('preserve payoff audit')
    summary, table = BASE/'summary.json', BASE/'exact_conductors.csv'
    benchmarks = json.loads(summary.read_text())['benchmarks']
    rows = []
    with table.open() as stream:
        for row in csv.DictReader(stream):
            rank, conductor = int(row['rank_lower_bound']), int(row['exact_conductor'])
            targets = [int(k) for k, v in benchmarks.items()
                       if int(k) > rank and conductor < int(v['recorded_conductor'])]
            if not targets:
                continue
            target = min(targets)
            benchmark = benchmarks[str(target)]
            rows.append({k: row[k] for k in ('inventory_id', 'family', 'parameter', 'exact_conductor')})
            rows[-1].update(current_lower_bound=rank, target_lower_bound=target,
                additional_directions=target-rank, benchmark_id=benchmark['id'],
                benchmark_conductor=benchmark['recorded_conductor'],
                conductor_ratio=str(Fraction(conductor, int(benchmark['recorded_conductor']))),
                missing_benchmark_conductors=benchmark['missing_conductor_ids'])
    rows.sort(key=lambda r: (r['additional_directions'], Fraction(r['conductor_ratio']), r['inventory_id']))
    result = {'status': 'PASS_EXACT_PINNED_BENCHMARK_COMPARISON',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in (summary, table)},
        'source_sha256': sha(Path(__file__)), 'point_searches': 0, 'rows': rows,
        'claim_boundary': 'Exact inequalities against a retained reported catalogue snapshot. '
            'Missing conductors remain missing. No current/universal record, rank attainability '
            'or discovery is asserted; refresh benchmarks after a certified gain.'}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(rows[0], sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output)
