#!/usr/bin/env python3
"""Exact content obstruction audit; written Hodge/gluing proof is an input."""
import argparse
import hashlib
import json
from collections import Counter
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / 'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
PLANNER = ROOT / 'artifacts/generated-results/elkies-k3-arithmetic-first-marked-t-foundry-v1.json'
OUTPUT = ROOT / 'artifacts/generated-results/elkies-k3-real-marking-content-v1.json'


def det3(g):
    a, b, c = g[0]
    d, e, f = g[1]
    h, i, j = g[2]
    return a*(e*j-f*i)-b*(d*j-f*h)+c*(d*i-e*h)


def build():
    catalogue = json.loads(INPUT.read_text())['surfaces']
    planner = json.loads(PLANNER.read_text())
    queue = {r['surface_id'] for r in planner['curve_identification_queue']}
    prior = {r['surface_id'] for r in planner['excluded_before_NS_or_equation_work']}
    positive = {r['surface_id'] for r in planner['positive_controls']}
    assert not (queue & prior or queue & positive or prior & positive)
    identifiers = {r['surface_id'] for r in catalogue}
    assert len(identifiers) == len(catalogue) == 827
    assert identifiers == queue | prior | positive
    rows = []
    for row in catalogue:
        g = row['literal_transcendental_gram']
        assert len(g) == 3 and all(len(r) == 3 for r in g)
        assert all(type(x) is int for r in g for x in r)
        assert all(g[i][j] == g[j][i] for i in range(3) for j in range(3))
        assert all(g[i][i] % 2 == 0 for i in range(3))
        assert det3(g) == -int(row['determinant'])
        m = gcd(*(x for r in g for x in r))
        assert m == row['similarity_normalization']['literal_content']
        rows.append({'surface_id': row['surface_id'], 'gram': g,
                     'determinant': -det3(g), 'content': m,
                     'full_real_marking_excluded': 2 % m != 0,
                     'previously_unresolved': row['surface_id'] in queue})
    excluded = {r['surface_id'] for r in rows if r['full_real_marking_excluded']}
    assert not positive & excluded
    survivors = [r for r in planner['curve_identification_queue']
                 if r['surface_id'] not in excluded]
    return {
        'schema': 'elkies-k3.real-marking-content.v1',
        'status': 'PASS_EXACT_INTEGER_AUDIT_WITH_WRITTEN_THEOREM_INPUT',
        'inputs': [{'path': str(p.relative_to(ROOT)),
                    'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                   for p in (INPUT, PLANNER)],
        'theorem_source': 'elkies-k3/REAL_MARKING_CONTENT_OBSTRUCTION_2026-09-14.md',
        'counts': {'catalogue': len(rows), 'excluded': len(excluded),
                   'already_excluded': len(excluded & prior),
                   'newly_excluded': len(excluded & queue),
                   'unresolved_remaining': len(survivors),
                   'positive_controls_excluded': len(excluded & positive)},
        'content_histogram': dict(sorted(Counter(str(r['content']) for r in rows).items())),
        'rows': sorted(rows, key=lambda r: r['surface_id']),
        'surviving_queue_ids': [r['surface_id'] for r in survivors],
        'surviving_low_genus_diagnostics': [r for r in planner['coarse_low_genus_diagnostic_shortlist']
                                         if r['surface_id'] not in excluded],
        'new_positive_handoff': [],
        'boundary': 'Content dividing two is necessary only. No full marked curve, rational point, rootless frame, or equation is certified by passing.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    payload = build()
    if args.check:
        assert json.loads(OUTPUT.read_text()) == payload, 'retained audit differs'
    else:
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    print(json.dumps(payload['counts'], sort_keys=True))
