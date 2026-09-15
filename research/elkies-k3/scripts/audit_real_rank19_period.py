#!/usr/bin/env sage-python
"""Replay the rank-19 real-place necessary-condition corollary."""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import QQ, QuadraticForm, hilbert_symbol, matrix

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / 'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
CONTENT = ROOT / 'artifacts/generated-results/elkies-k3-real-marking-content-v1.json'
PLANNER = ROOT / 'artifacts/generated-results/elkies-k3-arithmetic-first-marked-t-foundry-v1.json'
OUTPUT = ROOT / 'artifacts/generated-results/elkies-k3-real-rank19-period-v1.json'


def build():
    catalogue = json.loads(INPUT.read_text())['surfaces']
    content = json.loads(CONTENT.read_text())
    planner = json.loads(PLANNER.read_text())
    old_excluded = {r['surface_id'] for r in content['rows']
                    if r['full_real_marking_excluded']}
    survivors = set(content['surviving_queue_ids'])
    rows = []
    for row in catalogue:
        if row['surface_id'] in old_excluded:
            continue
        g = matrix(QQ, row['literal_transcendental_gram'])
        four_divisible = all(g[i, i] % 4 == 0 for i in range(3)) and all(
            g[i, j] % 2 == 0 for i in range(3) for j in range(3) if i != j)
        if not four_divisible or row['rational_isotropy']['isotropic']:
            continue
        assert g.det() == -row['determinant']
        q = QuadraticForm(QQ, g)
        assert q.signature_vector() == (2, 1, 0)
        diagonal = q.rational_diagonal_form()
        a, b, c = [diagonal[i, i] for i in range(3)]
        p = int(row['rational_isotropy']['pari_qfsolve_obstruction_prime'])
        assert p > 0
        symbol = int(hilbert_symbol(-a*b, -a*c, p))
        assert symbol == -1
        rows.append({'surface_id': row['surface_id'],
                     'gram': row['literal_transcendental_gram'],
                     'determinant': int(-g.det()),
                     'diagonal_quadratic_coefficients': [str(a), str(b), str(c)],
                     'clifford_parameters': [str(-a*b), str(-a*c)],
                     'obstruction_prime': p, 'hilbert_symbol': symbol,
                     'all_integral_squares_divisible_by_four': True,
                     'full_real_rank19_marking_excluded': True})
    new_excluded = {r['surface_id'] for r in rows}
    assert new_excluded <= survivors
    assert len(new_excluded) == len(rows) == 31
    assert not new_excluded & old_excluded
    positive = {r['surface_id'] for r in planner['positive_controls']}
    assert not positive & (old_excluded | new_excluded)
    return {'schema': 'elkies-k3.real-rank19-period.v1',
            'status': 'PASS_EXACT_LOCAL_REPLAY_WITH_WRITTEN_THEOREM',
            'inputs': [{'path': str(p.relative_to(ROOT)),
                        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                       for p in (INPUT, CONTENT, PLANNER)],
            'counts': {'catalogue': len(catalogue), 'additional_exclusions': len(rows),
                       'previous_content_exclusions': len(old_excluded),
                       'unresolved_remaining': len(survivors-new_excluded)},
            'rows': sorted(rows, key=lambda r: r['surface_id']),
            'surviving_queue_ids': [s for s in content['surviving_queue_ids']
                                    if s not in new_excluded],
            'surviving_low_genus_diagnostics': [r for r in content['surviving_low_genus_diagnostics']
                                               if r['surface_id'] not in new_excluded],
            'new_positive_handoff': [],
            'boundary': 'A necessary-condition pass does not certify a real or rational marking, a non-CM point, or an MW17 fibration.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.check:
        assert json.loads(OUTPUT.read_text()) == result, 'retained packet differs'
    else:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result['counts'], sort_keys=True))
