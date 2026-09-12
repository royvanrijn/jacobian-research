#!/usr/bin/env python3
"""Read-only censoring/ramification analysis of the frozen prospective cohort.

This is retrospective analysis of a prospectively frozen cohort, not a new
out-of-sample prediction test. No unknown arithmetic is imputed.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse

import historical_external_arithmetic as historic
import arithmetic_frozen_snapshot as frozen
import wide_arithmetic_profile_core as core

ROOT = historic.ROOT
DEFAULT = ROOT / 'artifacts/generated-results/elliptic-curves/two_class_relation_pilot_v1'


def cohort(rows):
    if any((r['cohort'], r['selection_mode']) != historic.PROSPECTIVE for r in rows):
        raise ValueError('historical data must not enter prospective audit')
    return rows


def stats(rows):
    known = [r for r in rows if r.get('field_ramified_prime_count') is not None]
    local = [r for r in rows if r['local_status'] == 'PASS']
    high = sum(r['field_ramified_prime_count'] >= 9 for r in known)
    return {'count': len(rows), 'local_pass': len(local), 'local_unknown': len(rows)-len(local),
            'completion_rate': len(local)/len(rows) if rows else None,
            'median_base_log2_delta': core.median([r['log2_abs_minimal_discriminant'] for r in rows]),
            'median_base_log2_delta_pass': core.median([r['log2_abs_minimal_discriminant'] for r in local]),
            'median_base_log2_delta_unknown': core.median([r['log2_abs_minimal_discriminant'] for r in rows if r['local_status'] != 'PASS']),
            'ramification_known': len(known), 'ramified_ge9_known': high,
            'ramified_ge9_population_fraction_bounds': [high/len(rows), (high+len(rows)-len(known))/len(rows)] if rows else None,
            'median_ramified_primes': core.median([r['field_ramified_prime_count'] for r in known]),
            'tail_ge23_count': sum(r['final_rank_lower_bound'] >= 23 for r in rows),
            'median_rank_lower_bound': core.median([r['final_rank_lower_bound'] for r in rows])}


def grouped(rows, key):
    groups = defaultdict(list)
    for r in rows:
        groups[str(key(r))].append(r)
    return {k: stats(v) for k, v in sorted(groups.items())}


def centered_spearman(rows, metric, group_key):
    rows = [r for r in rows if r.get(metric) is not None]
    groups = defaultdict(list)
    for r in rows:
        groups[group_key(r)].append(r)
    xs, ys = [], []
    used = 0
    for group in groups.values():
        a = [r[metric] for r in group]
        b = [r['final_rank_lower_bound'] for r in group]
        if len(set(a)) < 2 or len(set(b)) < 2:
            continue
        # Within-stratum ranks, then demean; a descriptive stratified rank
        # association, not a partial-correlation causal identification claim.
        ar, br = core._average_ranks(a), core._average_ranks(b)
        ax, by = sum(ar)/len(ar), sum(br)/len(br)
        xs.extend(x-ax for x in ar)
        ys.extend(y-by for y in br)
        used += 1
    denominator = (sum(x*x for x in xs)*sum(y*y for y in ys))**0.5
    return {'centered_rank_association': sum(x*y for x, y in zip(xs, ys))/denominator if denominator else None,
            'pairs': len(xs), 'informative_strata': used}


def analyze(rows):
    cohort(rows)
    values = sorted(r['log2_abs_minimal_discriminant'] for r in rows)
    # Boundaries depend on BASE data available for every row, never LOCAL
    # completion or rank. Tied boundary values stay in the same band.
    cuts = [values[len(values)*i//5] for i in range(1, 5)]
    band = lambda r: sum(r['log2_abs_minimal_discriminant'] >= c for c in cuts) + 1
    rambin = lambda r: ('UNKNOWN' if r.get('field_ramified_prime_count') is None else
                       '0-4' if r['field_ramified_prime_count'] <= 4 else
                       '5-6' if r['field_ramified_prime_count'] <= 6 else
                       '7-8' if r['field_ramified_prime_count'] <= 8 else '9+')
    family11952 = [r for r in rows if r['family'] == '11952']
    complete = [r for r in rows if r.get('field_ramified_prime_count') is not None]
    correlation = lambda rs: {'n': len(rs), 'spearman_ramification_vs_rank_LB': core.spearman(
        [r['field_ramified_prime_count'] for r in rs], [r['final_rank_lower_bound'] for r in rs])}
    return {'schema': 'elliptic-curves.wide-censoring-ramification-audit.v1',
            'cohort': historic.PROSPECTIVE[0], 'selection_mode': historic.PROSPECTIVE[1],
            'analysis_mode': 'retrospective_analysis_of_frozen_prospective_population',
            'population': stats(rows), 'by_rank_stratum': grouped(rows, lambda r: r['rank_bucket']),
            'by_family': grouped(rows, lambda r: r['family']),
            'base_size_quintile_cuts': cuts, 'by_base_size_quintile': grouped(rows, band),
            'by_family_and_size_quintile': grouped(rows, lambda r: f"{r['family']}/Q{band(r)}"),
            'ramification_strata': grouped(rows, rambin),
            'same_family_11952': {'population': stats(family11952),
                'by_rank_stratum': grouped(family11952, lambda r: r['rank_bucket']),
                'ramification_strata': grouped(family11952, rambin)},
            'completed_case_ramification_association': correlation(complete),
            'completed_11952_ramification_association': correlation([r for r in complete if r['family'] == '11952']),
            'within_family_rank_association': centered_spearman(rows, 'field_ramified_prime_count', lambda r: r['family']),
            'within_family_and_size_rank_association': centered_spearman(rows, 'field_ramified_prime_count', lambda r: (r['family'], band(r))),
            'completion_vs_base_size_spearman': core.spearman(
                [r['log2_abs_minimal_discriminant'] for r in rows], [int(r['local_status']=='PASS') for r in rows]),
            'missing_data_boundary': 'Worst-case fraction intervals assign every unknown to low or high ramification. No missing-at-random assumption, inverse weighting or imputation.',
            'rank_boundary': 'All ranks are certified lower bounds; no exact ranks, p-values, predictive validation, Selmer dimensions or class ranks inferred.'}


def text_table(headers, rows):
    def cell(v):
        return ('UNKNOWN' if v is None else f'{v:.4f}' if isinstance(v, float) else str(v)).replace('|', r'\|')
    return '\n'.join(['| '+' | '.join(cell(h) for h in headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(cell(v) for v in row)+' |' for row in rows])


def markdown(data):
    lines = ['# Censoring and prospective-population ramification audit', '',
             'Retrospective analysis of the unchanged prospective 2,080-fibre population. Historical controls are excluded. No new arithmetic, searches or class groups.', '',
             '## LOCAL completion by BASE discriminant-size quintile', '',
             'Quintile boundaries use BASE information available on every fibre, without rank or LOCAL outcomes. Boundaries: '+str(data['base_size_quintile_cuts'])+'.', '',
             text_table(['Quintile', 'Fibres', 'LOCAL PASS', 'Completion', 'Median log2 |Delta|'],
                        [[k, v['count'], v['local_pass'], v['completion_rate'], v['median_base_log2_delta']] for k, v in data['by_base_size_quintile'].items()]), '',
             '## Completion by certified lower-bound stratum', '',
             text_table(['Rank LB', 'Fibres', 'LOCAL PASS', 'Completion', 'Known ramification >=9', 'Worst-case population fraction'],
                        [[k, v['count'], v['local_pass'], v['completion_rate'], v['ramified_ge9_known'], v['ramified_ge9_population_fraction_bounds']] for k, v in data['by_rank_stratum'].items()]), '',
             '## Parent/family completion', '',
             text_table(['Family', 'Fibres', 'LOCAL PASS', 'Completion', 'Median ramified primes (known)'],
                        [[k, v['count'], v['local_pass'], v['completion_rate'], v['median_ramified_primes']] for k, v in data['by_family'].items()]), '']
    for label, strata in [('All prospective', data['ramification_strata']),
                           ('Prospective 11952 only', data['same_family_11952']['ramification_strata'])]:
        lines += ['## '+label+': ramification strata', '',
                  text_table(['Ramified primes', 'Fibres', 'LB >=23', 'Tail fraction', 'Median rank LB', 'Median BASE log2 |Delta|'],
                    [[k, v['count'], v['tail_ge23_count'], v['tail_ge23_count']/v['count'], v['median_rank_lower_bound'], v['median_base_log2_delta']]
                     for k, v in strata.items()]), '']
    lines += ['## Descriptive diagnostics', '', core.stable_json({k: data[k] for k in (
        'completion_vs_base_size_spearman', 'completed_case_ramification_association',
        'completed_11952_ramification_association', 'within_family_rank_association',
        'within_family_and_size_rank_association')}), '', data['missing_data_boundary'], '', data['rank_boundary'], '',
        'A completion association is evidence of observed censoring structure, not proof of why a particular worker timed out. '
        'Conditioning on completed rows, family and coarse size bands does not remove unobserved selection bias. '
        'This audit cannot establish or refute a class-group mechanism.', '']
    return '\n'.join(lines)


def build(out):
    source = historic.DEFAULT_OUTPUT
    frozen.verify()
    rows = historic.read(source/'prospective_profiles.json')['rows']
    if len(rows) != 2080:
        raise ValueError('population fingerprint mismatch')
    data = analyze(rows)
    data['sources'] = {historic.worker.rel(p): historic.sha(p) for p in
                       (source/'prospective_profiles.json', source/'plan.json', Path(__file__), Path(core.__file__), Path(frozen.__file__))}
    return {'censoring_audit.json': core.stable_json(data), 'CENSORING.md': markdown(data)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['run', 'check'])
    ap.add_argument('--output', type=Path, default=DEFAULT)
    args = ap.parse_args()
    outputs = build(args.output)
    if args.command == 'run':
        args.output.mkdir(parents=True, exist_ok=True)
        for name, text in outputs.items():
            p = args.output/name
            if p.exists() and p.read_bytes() != text.encode():
                raise ValueError('refusing changed audit '+name)
            if not p.exists():
                p.write_bytes(text.encode())
    else:
        assert all((args.output/name).read_bytes() == text.encode() for name, text in outputs.items())
    print('CENSORING_AUDIT|PASS|prospective=2080|historical_excluded=PASS|no_new_arithmetic=PASS')


if __name__ == '__main__':
    main()
