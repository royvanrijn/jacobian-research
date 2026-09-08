#!/usr/bin/env python3
"""Join the public28 control and retained norm preflight to parallel results.

This replays joins and finite rank evidence; it does not recompute the
earlier local-boundary or full class-group certificates.
"""
import argparse
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import json
from research_runtime.norm_ramification import cubic_data
from memory_rank_certificate import checked_rank

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'search_result_integration_v1.json'


def calculate():
    sources = {}
    def read(name):
        path = ART/(name+'.json')
        sources[str(path.relative_to(ROOT))] = sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    public = read('inventory188_public28_reproduction_v1')
    rows = read('rank_jump_fresh_strict_boundary_coordinate_comparison_v1')
    boundary = next(r for r in rows['rows'] if r['token'] == 'case-06')
    models = read('rank_jump_fresh_governing_panel_inputs_v1')
    model = next(r for r in models['cases'] if r['token'] == boundary['token'])
    assert boundary['status'] == public['status'] == 'PASS'
    assert boundary['parameter'] == public['parameter'] == '110314/102227'
    assert boundary['id'] == public['local_id']
    assert list(map(Q, model['model'])) == list(map(Q, public['curve']))
    proof = public['rank_certificate']
    replay = checked_rank(tuple(map(Q, public['curve'])),
                          [tuple(map(Q, p)) for p in public['independent_points']],
                          [r['prime'] for r in proof['signatures']],
                          proof['no_rational_2_torsion_prime'])
    assert json.dumps(replay, sort_keys=True) == json.dumps(proof, sort_keys=True)
    R = replay['rank_lower_bound']
    m, k, h, a = [boundary[key] for key in ['generic_rank', 'generic_strict_dimension',
                                          'boundary_upper_bound', 'additional_boundary_capacity_upper_bound']]
    assert a == h-(m-k) and (m,k,h,a,R) == (17,0,21,4,28)
    projection = read('retained_norm_preflight_inputs_v1')
    preflight = read('retained_norm_preflight_v1')
    independent = read('retained_norm_preflight_sage_v1')
    assert preflight['status'] == independent['status'] == 'PASS'
    panel = read('rank_jump_fresh_norm_projection_v1')
    batch = read('rank_jump_retained_norm_batch_capacity_inputs_v1')
    support = {r['token']: [c['place'] for c in r['local'] if c['place'] != 'infinity']
               for r in panel['rows']}
    support['retained-reference-296'] = batch['S_finite']
    counts = []
    for case in projection['cases']:
        _, forbidden = cubic_data(case['cubic_ascending'], case['elements'])
        # All recorded bad places are excluded from each isolation witness.
        assert all(forbidden % int(p) == 0 for p in support[case['id']])
        counts.append({'id': case['id'], 'bad_places_checked': len(support[case['id']])})
    sources[str(Path(__file__).resolve().relative_to(ROOT))] = sha256(Path(__file__).read_bytes()).hexdigest()
    return {'schema': 'elliptic-curves.search-result-integration.v1', 'status': 'PASS',
            'sources': sources,
            'public28_strict_consequence': {
                'local_id': public['local_id'], 'rank_lower_bound_replayed': R,
                'generic_dimension': m, 'generic_strict_dimension': k,
                'additional_boundary_capacity_upper_bound': a,
                'additional_strict_rational_dimension_lower_bound': max(0,R-m-a),
                'localized_class_dimension_lower_bound': max(k,R-h),
                'previous_rank27_strict_lower_bound': 6,
                'localized_class_dimension_upper_bound': 'UNKNOWN',
                'rank_upper_bound': 'UNKNOWN',
                'provenance': 'Public-point reproduction joined AFTER masked boundary arithmetic; never a prospective selection feature.'},
            'norm_preflight_bad_support_checks': counts,
            'norm_preflight_forced_zero_count': preflight['totals']['forced_zero'],
            'production_curve_exclusions': 0,
            'score_changes': 0,
            'following_parameter_scan': None,
            'boundary': 'Inherited certificates supply the boundary and bad sets; only the joins, support inclusion and public finite rank evidence are replayed here.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.check:
        assert result == json.loads(OUT.read_text())
    else:
        with OUT.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print('PASS: public28 forces at least seven additional strict rational directions; 428 preflight exclusions respect bad support.')
