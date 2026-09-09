#!/usr/bin/env python3
"""Post-execution analysis only: no selection, point search, or extra budgets."""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from statistics import median

import audit_r17_60_panel as audit
from v3_warm_support import atomic, read, require, sha

panel = audit.panel
OUTPUT = panel.ART / 'r17_60_panel_analysis_v1.json'


def summarize(rows):
    seeded = [r for r in rows if r['first_seed']]
    return {
        'cases': len(rows), 'first_seeds': len(seeded),
        'first_call_seeds': sum(r['first_seed'] and r['seed_calls'] == 1 for r in rows),
        'rank_lower_bound_histogram': dict(sorted(Counter(str(r['rank_lower_bound']) for r in rows).items())),
        'mean_rank_lower_bound': str(Fraction(sum(r['rank_lower_bound'] for r in rows), len(rows))) if rows else None,
        'seed_calls': sum(r['seed_calls'] for r in rows),
        'median_calls_to_seed_conditional': str(median([Fraction(r['seed_calls']) for r in seeded])) if seeded else None,
        'complement_calls': sum(r['complement_calls'] for r in rows),
        'seed_cloud_added_directions': sum(r['seed_cloud_rank'] - 17 for r in rows),
        'complement_added_directions': sum(r['complement_added_directions'] for r in rows),
        'seeded_without_complement_gain': sum(r['complement_added_directions'] == 0 for r in seeded),
        'late_gain_cases': sum((r['last_complement_gain_call'] or 0) >= 75 for r in rows),
        'point_timeouts': sum(r['point_timeouts'] for r in rows),
        'map_timeouts': sum(r['map_timeouts'] for r in rows),
        'stop_reasons': dict(sorted(Counter(r['stop_reason'] for r in rows).items())),
    }


def build(complete=True, proofs=False):
    checked = audit.audit(complete=complete, proofs=proofs)
    rows, evidence = [], {}
    for item in checked['verified']:
        folder = panel.D / 'cases' / item['id']
        path = folder / 'result.json'
        r = read(path)
        cloud = read(folder / 'seed-reconciled.json')['rank_lower_bound']
        rows.append({
            'id': r['id'], 'family': r['family'], 'parameter': r['parameter'], 'stratum': r['stratum'],
            'first_seed': r['first_M18_found'], 'seed_calls': r['seed_calls'],
            'seed_cloud_rank': cloud, 'rank_lower_bound': r['rank_lower_bound'],
            'complement_calls': r['complement_calls'],
            'complement_added_directions': r['rank_lower_bound'] - cloud,
            'gaining_complement_calls': r['gaining_complement_calls'],
            'calls_per_added_direction': r['calls_per_added_direction'],
            'last_complement_gain_call': r['last_complement_gain_call'],
            'calls_since_last_complement_gain': r['complement_calls_since_last_gain'],
            'gain_timeline': r['gain_timeline'], 'stop_reason': r['stop_reason'],
            'point_timeouts': r['point_timeouts'], 'map_timeouts': r['map_timeouts'],
        })
        evidence[str(path.relative_to(panel.ROOT))] = sha(path)
    def key(row):
        return (-row['rank_lower_bound'], -row['gaining_complement_calls'],
                Fraction(row['calls_per_added_direction']) if row['calls_per_added_direction'] else Fraction(10**9),
                -(row['last_complement_gain_call'] or 0), row['id'])
    groups = {'all': summarize(rows)}
    for name in ('low', 'high', *panel.selection.FAMILIES):
        groups[name] = summarize([r for r in rows if r['stratum'] == name or r['family'] == name])
    return {
        'schema': 'r17-sixty-panel-analysis.v1',
        'status': 'COMPLETE_INDEPENDENTLY_VERIFIED_PANEL_ANALYSIS' if complete else 'PARTIAL_PANEL_ANALYSIS',
        'protocol_sha256': sha(panel.D / 'protocol.json'), 'point_searches': 0,
        'cases': rows, 'groups': groups, 'ranking': [r['id'] for r in sorted(rows, key=key)],
        'late_gain_candidates': [r['id'] for r in sorted(rows, key=key) if (r['last_complement_gain_call'] or 0) >= 75],
        'pending': checked['pending'], 'unresolved': checked['unresolved'],
        'bindings': evidence,
        'analysis_sources': {str(p.relative_to(panel.ROOT)): sha(p) for p in (Path(__file__), Path(audit.__file__))},
        'claim_boundary': 'Certified subgroup lower bounds and outcomes of one selected finite panel. Seed incidence is bounded detection, not arithmetic absence. Saved-cloud gains use original call positions. No exact ranks, public novelty, conductor records, general rank distribution, causal height effect, or performance beyond the frozen budget is established.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--partial', action='store_true')
    parser.add_argument('--proofs', action='store_true', help='run with Sage -python')
    args = parser.parse_args()
    result = build(complete=not args.partial, proofs=args.proofs)
    if args.write:
        atomic(OUTPUT, result, immutable=True)
    if args.check:
        require(read(OUTPUT) == result, 'analysis replay differs')
    print(json.dumps({'status': result['status'], 'groups': result['groups'],
                      'ranking': result['ranking'], 'late_gain_candidates': result['late_gain_candidates']}, indent=2))
