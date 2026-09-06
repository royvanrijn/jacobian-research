#!/usr/bin/env python3
"""A frozen post-search comparison panel, independent of live search state."""
import argparse
from collections import Counter
import csv
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import subprocess
import retrospective as r

PROTOCOL = Path(__file__).with_name('COMPLETED_COHORT_PANEL_PROTOCOL.json')
INPUT = r.OUT/'rank_jump_completed_cohort_inputs_v1.json'
OUTPUT = r.OUT/'rank_jump_completed_cohort_panel_v1.json'
CSV = r.OUT/'rank_jump_completed_cohort_panel_v1.csv'
SOURCES = ['compact192_r17_results_v1.json', 'extended20_mw16_results_v1.json', 'new_mw16_followup_results_v1.json']


def height(parameter):
    q = F(parameter)
    return max(abs(q.numerator), q.denominator)


def export():
    rows, bindings = [], {}
    for name in SOURCES:
        path = r.OUT/name; data = r.read(path)
        bindings[str(path.relative_to(r.ROOT))] = r.digest(path.read_bytes())
        for old in data['curves']:
            m = len(old['generic_points'])
            assert m == (16 if old['family'].startswith('a1') else 17)
            assert old['points'][:m] == old['generic_points']
            adaptive = name == SOURCES[2]
            rows.append({'observation_id': name+':'+old['id'], 'cohort': name,
                         'family': old['family'], 'parameter': old['parameter'], 'source_id': old['id'],
                         'phase': 'adaptive' if adaptive else 'initial', 'model': old['curve'],
                         'generic_rank': m, 'points': old['points'],
                         'rank_lower_bound': old['rank_lower_bound'],
                         'rank_certificate': old['rank_certificate'],
                         'declared_initial_boxes': old['declared_charts'],
                         'completed_initial_boxes': old['completed_boxes'],
                         'completed_adaptive_boxes': old.get('completed_adaptive_boxes', 0),
                         'search_status': old['search_status'],
                         'adaptive_search_status': old.get('adaptive_search_status'),
                         'height_limit': 125000 if name == SOURCES[0] else 100000,
                         'prior_equation_matches': old['previous_matches'], 'catalogue_matches': old['icarm_matches']})
    assert len(rows) == 215
    r.write_new(INPUT, {'schema': 'rank-jump.completed-cohort-inputs.v1',
                       'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=r.ROOT, text=True).strip(),
                       'source_hashes': bindings, 'rows': rows,
                       'boundary': 'Copied proof inputs are immutable and replay without reading live search files. Rows with repeated equations are separate exposure observations.'})


def compute():
    inp = r.read(INPUT); rows = []
    for old in inp['rows']:
        model, points = r.short(old['model'], old['points'])
        assert model == old['model']
        m, n = old['generic_rank'], len(points)
        assert n == old['rank_lower_bound'] == old['rank_certificate']['rank_lower_bound']
        primes = [s['prime'] for s in old['rank_certificate']['signatures']]
        assert max(primes) <= 503
        blocks = [(p, r.roots_at(model[3], model[4], p)) for p in primes]
        assert all(roots for p, roots in blocks)
        sigs = [r.point_signature(model, P, blocks) for P in points]
        assert r.rank(sigs) == n and r.rank(sigs[:m]) == m
        # Equation-only cubic structure, with its own bounded prime dictionary.
        A, B = map(F, model[3:]); disc = -4*A**3-27*B*B
        assert disc
        square = disc > 0 and isqrt(disc.numerator)**2 == disc.numerator and isqrt(disc.denominator)**2 == disc.denominator
        irreducible = next((p for p in r.primes(503) if r.roots_at(model[3], model[4], p) == ()), None)
        assert irreducible is not None
        j = 1728*4*A**3/(4*A**3+27*B*B)
        assert old['completed_initial_boxes'] == old['declared_initial_boxes']
        if old['phase'] == 'adaptive': assert old['completed_adaptive_boxes'] == 301
        rows.append({k: old[k] for k in ['observation_id', 'cohort', 'family', 'parameter', 'source_id', 'phase',
                                        'completed_initial_boxes', 'completed_adaptive_boxes', 'height_limit']} | {
            'parameter_height': height(old['parameter']), 'j_invariant': str(j),
            'generic_subgroup_rank_exact': m, 'witness_subgroup_rank_exact': n,
            'witness_quotient_rank_exact': n-m, 'full_curve_rank': 'UNKNOWN',
            'full_quotient_rank': 'UNKNOWN', 'cubic_Galois_group': 'C3' if square else 'S3',
            'cubic_discriminant': str(disc), 'cubic_irreducibility_prime': irreducible,
            'recomputed_Kummer_fingerprints': sigs,
            'relative_joint_halving_degree': 4**(n-m),
            'halving_degree_base': 'Q(E[2], halves of the marked generic generators)',
            'common_rational_section_or_solubility_construction': 'UNKNOWN'})
    initial = [x for x in rows if x['phase'] == 'initial']
    assert len(initial) == 212
    # All initial j values are checked, rather than silently treating observations as curves.
    jgroups = {}
    for row in initial: jgroups.setdefault(row['j_invariant'], []).append(row['observation_id'])
    pairs = []
    compact = [x for x in initial if x['cohort'] == SOURCES[0]]
    for fam in sorted({x['family'] for x in compact if x['witness_subgroup_rank_exact'] == 26}):
        high = [x for x in compact if x['family'] == fam and x['witness_subgroup_rank_exact'] == 26]
        low = [x for x in compact if x['family'] == fam and x['witness_quotient_rank_exact'] == 0]
        def key(pair):
            a, b = pair; ah, bh = a['parameter_height'], b['parameter_height']
            return F(max(ah,bh),min(ah,bh)), abs(F(a['parameter'])-F(b['parameter'])), a['source_id'], b['source_id']
        a, b = min(((a,b) for a in high for b in low), key=key)
        pairs.append({'family': fam, 'high': a['observation_id'], 'low': b['observation_id'],
                      'parameter_height_ratio': str(key((a,b))[0]),
                      'absolute_parameter_distance': str(key((a,b))[1]),
                      'same_initial_box_count': a['completed_initial_boxes'] == b['completed_initial_boxes'],
                      'score_matched': False})
    counts = []
    for name in SOURCES:
        subset = [x for x in rows if x['cohort'] == name]
        counts.append({'cohort': name, 'observations': len(subset),
                       'quotient_lower_bound_counts': dict(sorted(Counter(str(x['witness_quotient_rank_exact']) for x in subset).items())),
                       'Galois_counts': dict(Counter(x['cubic_Galois_group'] for x in subset))})
    return {'schema': 'rank-jump.completed-cohort-panel.v1',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (Path(__file__), PROTOCOL, INPUT, Path(r.__file__))},
            'status': 'PASS', 'rows': rows, 'cohort_counts': counts, 'paired_cases': pairs,
            'initial_distinct_j_count': len(jgroups),
            'initial_equal_j_groups': [v for v in jgroups.values() if len(v)>1],
            'boundary': 'A measured lower-bound panel with complete declared initial boxes. Censored ranks remain UNKNOWN. Shared-field, CT and auxiliary-curve explanations beyond the proved halving-field statement remain uncomputed.'}


def csv_text(data):
    import io
    keys = ['observation_id', 'phase', 'family', 'parameter', 'parameter_height', 'generic_subgroup_rank_exact',
            'witness_subgroup_rank_exact', 'witness_quotient_rank_exact', 'full_curve_rank', 'cubic_Galois_group',
            'height_limit', 'completed_initial_boxes', 'completed_adaptive_boxes']
    s = io.StringIO(); w = csv.DictWriter(s, fieldnames=keys, lineterminator='\n'); w.writeheader()
    for row in data['rows']: w.writerow({k: row[k] for k in keys})
    return s.getvalue()


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['export', 'build', 'check']); a = p.parse_args()
    if a.mode == 'export': export()
    else:
        data = compute()
        if a.mode == 'check':
            assert r.read(OUTPUT) == data and CSV.read_text() == csv_text(data); print('PASS completed-cohort panel')
        else:
            r.write_new(OUTPUT, data)
            with CSV.open('x') as f: f.write(csv_text(data))
            print(data['status'], data['initial_distinct_j_count'], data['cohort_counts'], data['paired_cases'])
