#!/usr/bin/env python3
"""Deterministic read-only analysis of completed ancestry and skew-sieve rows."""
import argparse
from collections import Counter
import json
from pathlib import Path
import historical_external_arithmetic as storage

ROOT = storage.ROOT
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'rank_ancestry_principal28_v1'
SIEVE = ART/'two_class_skew_sieve_v1'
OLD = ART/'two_class_relation_pilot_v1'


def median_bits(histogram):
    target = sum(histogram.values())//2
    count = 0
    for k, v in sorted(histogram.items(), key=lambda kv: int(kv[0])):
        count += v
        if count > target:
            return int(k)


def report():
    geometry, verified, branches = [storage.read(OUT/name) for name in ['geometry.json', 'verified.json', 'branch-analysis.json']]
    assert verified['status'] == 'PASS_EXACT_14_VS_11_ANCESTRY_REPLAY'
    assert branches['status'] == 'PASS_EXACT_SELECTED_CARRIER_BRANCH_ANALYSIS'
    for payload in [verified, branches]:
        for rel, digest in payload['bindings'].items():
            assert storage.sha(ROOT/rel) == digest
    rows = []
    for cid in ['302', '11952']:
        chosen = [r for r in geometry['rows'] if r['fibre'] == cid]
        b = next(b for b in branches['rows'] if b['fibre'] == cid)
        rows.append({'fibre': cid, 'parameter': chosen[0]['parameter'], 'rank_lower_bound': 17+len(chosen),
                     'known_extra_directions': len(chosen), 'minimum_degree_each': 2,
                     'best_genus_each': sorted({r['best_genus_in_atlas'] for r in chosen}),
                     'global_minimum_genus': 'UNKNOWN', 'tested_covers': sum(len(r['covers']) for r in chosen),
                     'distinct_fixed_base_covers': sum(c['fibre'] == cid for c in geometry['cover_equivalence_classes']),
                     'shared_target_cover_count': sum(c['fibre'] == cid and len(c['targets']) > 1 for c in geometry['cover_equivalence_classes']),
                     'genus_histogram': dict(sorted(Counter(c['normalization_genus'] for r in chosen for c in r['covers']).items())),
                     'selected_carrier_pencils': 1, 'selected_distinct_carriers': len(chosen),
                     'pairwise_disjoint_selected_branch_divisors': b['pairwise_disjoint_geometric_branch_divisors'],
                     'selected_radical_squareclass_rank': b['quadratic_squareclass_rank_over_algebraic_closure_of_constants'],
                     'minimum_common_cover_degree_for_selected_extensions': b['compositum_degree'],
                     'same_generic_trace_in_pencil': b['same_exact_generic_trace']})
    sieve = []
    supervisor = storage.read(SIEVE/'verified.json')
    for rel, digest in supervisor['bindings'].items():
        assert storage.sha(SIEVE/rel) == digest
    for row in storage.read(SIEVE/'plan.json')['rows']:
        key = row['curve_key']
        s = storage.read(SIEVE/key/'supervisor.json')
        if s['status'] != 'PASS_BOUNDED_WORKER':
            sieve.append({'parameter': row['t'], 'status': s['status'], 'global_class_2rank_estimate': None})
            continue
        r = storage.read(SIEVE/key/'result.json')
        old = storage.read(OLD/key/'result.json')
        assert r['global_class_2rank_estimate'] is None and r['global_class_2rank_upper_bound'] is None
        sieve.append({'parameter': row['t'], 'status': r['status'], 'primitive_candidates': r['attempted'],
                      'noncanonical_relations': r['noncanonical_relations'], 'matrix_deficiency': r['final_matrix']['deficiency'],
                      'median_norm_numerator_bits_before': median_bits(old['norm_numerator_bit_histogram']),
                      'median_norm_numerator_bits_after': median_bits(r['norm_numerator_bit_histogram']),
                      'new_relation_pairs': [[v['a'], v['b']] for v in r['relations']],
                      'density_gate': r['density_gate'], 'global_class_2rank_estimate': None})
    result = {'status': 'PASS_PRINCIPAL28_BOUNDED_COMPARISON', 'principal_control': '11952 @ 110314/102227, rank >=28',
              'secondary_control': '11952 @ 921/653, rank >=25', 'comparison': rows, 'bounded_sieve': sieve,
              'interpretation': 'Similar existing genus-one pencils, but independent selected quadratic covers. No common-cover decomposition of the rank jumps is established.',
              'global_minimum_genus': 'UNKNOWN', 'no_new_point_searches': True,
              'sources': {str(p.relative_to(ROOT)): storage.sha(p) for p in [Path(__file__), OUT/'geometry.json', OUT/'verified.json', OUT/'branch-analysis.json', SIEVE/'verified.json']}}
    lines = ['# Principal-control ancestry: Curve302 +14 versus 11952 +11', '',
             'The principal control is **11952 at 110314/102227, rank at least28**. The 921/653 fibre remains a secondary rank-at-least25 control. All ranks here are certified lower bounds.', '',
             '| Frozen comparison | Curve302 at0 | 11952 at110314/102227 |',
             '|---|---:|---:|', '| Known extra directions over M17 | 14 | 11 |',
             '| Minimum covering degree, every target | 2 | 2 |', '| Best genus in this atlas, every target | 1 | 1 |',
             '| Global minimum genus | UNKNOWN | UNKNOWN |', '| Covers checked | 504 | 396 |',
             '| Distinct covers over the fixed base | 504 | 396 |', '| Covers shared by different chosen targets | 0 | 0 |',
             '| Existing genus-one pencils used | 1 | 1 |', '| Distinct genus-one carriers selected | 14 | 11 |',
             '| Independent selected quadratic radicals over Qbar(t) | 14 | 11 |',
             '| Degree needed to contain all selected quadratic extensions | 16,384 | 2,048 |', '',
             'These last degrees apply **only to the selected covers**, not to every possible ancestry of the points.', '',
             '## What is the same, and what is not established', '',
             'Both rosters lie on genus-one bisections from one existing alternate-fibration pencil: class1 on X1092 for302, and the old-R17 pencil on X948 for11952. Within each pencil the exact generic trace is constant, but the fourteen/eleven pencil parameters and branch quartics are distinct. All146 within-fibre pairwise branch gcds are1. Consequently these particular quadratic covers are independent even after extending the constants to Qbar; they do not collapse to a few common quadratic base changes.', '',
             'This supports similarity of the available low-degree carrier constructions, not a general mechanism for extreme rank jumps. There is no new qualitative difference in the extra three Curve302 directions in this atlas. Nor is this outcome C: simple ancestry exists. Rational bisections outside the atlas, full Neron-Severi orbit labels, other representative choices, and richer shared-cover families remain unresolved.', '',
             'The common atlas has36 covers per target: vertical x,34 signed generic-section chords, and one existing alternate-fibration carrier. The genus counts **per target** are identical within each fibre:', '',
             '| Genus | Curve302 | 11952 >=28 |', '|---|---:|---:|', '| 1 | 1 | 1 |', '| 3 | 30 | 10 |', '| 5 | 5 | 21 |', '| 7 | 0 | 4 |', '',
             'These higher-genus counts depend on the chosen parent presentation and atlas; they are not intrinsic rank invariants.', '',
             '## Every target', '',
             'Each linked record retains the exact original point, parent-model transport, all36 equations and maps, the monic branch divisor, constant twist and rational lift. Cover index35 (zero-based) is the unique genus-one winner in every row. Its branch divisor has four simple finite geometric points and no branch at infinity. Positive-genus carriers do not receive a smooth-rational-bisection orbit label.', '',
             '| Fibre | Target | Degree | Best genus | Carrier pencil | Exact record |', '|---|---|---:|---:|---|---|']
    for row in geometry['rows']:
        cid, j = row['fibre'], row['target']
        family = 'class1' if cid == '302' else 'old R17'
        lines.append(f'| {cid} | E{j} | 2 | 1 | {family} | [point, covers, branch maps](targets/{cid}-E{j}.json) |')
    lines += ['', 'The 11952 independent28 roster is the old27-point subgroup plus public union column53 (zero-based). Its first17 columns match the specialized generic basis in order, with recorded signs. E1–E10 are the old extra directions; E11 is the independent public addition. Curve302 retains exactly the previous fourteen representatives.', '',
              '## Certificates and minimum-degree argument', '',
              'Fresh finite-reduction certificates establish independence of31 and28 points. If a target lay on a degree-one multisection, that curve would be a generic section. Generic MW rank17 makes a nonzero multiple of it lie in the rational span of M17 (up to torsion); specialization at the smooth target fibre would contradict the certified independence. The constructed degree-two maps therefore attain the global minimum degree. They do not exclude genus0 within degree2.', '',
              '[Exact replay](verified.json):900 polynomial map identities,900 altered-ordinate negative controls, both generic prefixes and both rank certificates. All504 earlier Curve302 cover payloads are unchanged. [Branch analysis](branch-analysis.json) supplies every gcd and trace identity. Pairwise-disjoint branch support proves radical independence by valuations; Riemann–Hurwitz gives compositum genus212,993 or20,481. These are properties of these covers only.', '',
              '## One bounded mod-two sieve improvement', '',
              'The frozen small factor bases are unchanged. Equation coefficients select a dyadic skew rectangle minimizing an exact norm triangle bound, with A*B=2^18. Exact modular-root progression sieving removes all supported prime powers without a score cutoff. An independent primorial-GCD pass replays every primitive candidate, and Sage independently checks the mod-two matrix rank. Each field has a180-second cap; neither needed it.', '',
              '| Fibre | Primitive candidates | Noncanonical relations | Deficiency | Median norm bits, old -> skew |', '|---|---:|---:|---:|---:|']
    for r in sieve:
        if r['status'].startswith('PASS'):
            lines.append(f"| {r['parameter']} | {r['primitive_candidates']:,} | {r['noncanonical_relations']} | {r['matrix_deficiency']} | {r['median_norm_numerator_bits_before']} -> {r['median_norm_numerator_bits_after']} |")
        else:
            lines.append(f"| {r['parameter']} | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |")
    lines += ['', 'The sole historical relation is again(a,b)=(0,1): no new independent relation beyond the earlier pilot. Lower norm sizes did not produce serious relation density. Stop at this gate; no expanded class-group profile is justified by this experiment. This is a larger/skewed search region, not an equal-candidate causal estimate of the sieve change.', '',
              'Deficiency is **not an estimate of full class-group2-rank**: factor-base generation and omitted relations remain unproved. The distinctions between mod-two relation matrices, generation bounds and skew selection follow [Klagsbrun–Sherman–Weigandt, sections4–5](https://arxiv.org/abs/1606.07178). This is a small exact prototype, not their production NFS or certified Julia reduction.', '',
              'The existing consequences g(921/653)>=16 and g(110314/102227)>=18 still use known MW lower bounds and are not independent explanatory evidence. No BNF, new local-arithmetic correlations, or point searches were run. The original census, historical arithmetic and first relation-pilot directories remain byte-identical.', '']
    return result, '\n'.join(lines)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    result, text = report()
    if args.check:
        assert storage.read(OUT/'REPORT.json') == json.loads(json.dumps(result))
        assert (OUT/'SUMMARY.md').read_text() == text
        print('PASS byte-recomputed principal28 analysis')
    else:
        storage.save(OUT/'REPORT.json', result)
        if (OUT/'SUMMARY.md').exists():
            assert (OUT/'SUMMARY.md').read_text() == text
        else:
            (OUT/'SUMMARY.md').write_text(text)
        print('PASS principal28 report')
