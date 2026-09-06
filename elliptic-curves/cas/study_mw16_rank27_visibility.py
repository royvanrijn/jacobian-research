#!/usr/bin/env python3
"""Exact pointwise visibility certificate for MW16 at -1867/270.

The input preserves the chart geometry and raw square witnesses needed here.
It is retrospective: no claim about all representatives of a quotient class.
"""
import argparse
import json
from pathlib import Path

import certify_compact_r17_candidates as cert
import search_observability as obs
from alternate_quartic_covers import short_add
from half_lattice_pointed_sieve import linear_combination
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
INPUT = ART / 'mw16_rank27_visibility_input_v1.json'
OUTPUT = ART / 'mw16_rank27_visibility_v1.json'


def prepare():
    if INPUT.exists():
        raise FileExistsError('preserve visibility input')
    paths = [LOCAL / 'extended20-mw16-pari-v1/a1-fibration-01-052/result.json',
             LOCAL / 'mw16-new26-a1-fibration-01-052-adaptive-v1/result.json',
             ART / 'new_mw16_rank27_minimal_proof_v1.json']
    initial, adaptive, proof = [cert.read(p) for p in paths]
    point_proof = proof['curves'][0]
    arms = {}
    keys = ('base_point', 'short_model', 'short_model_x_shift', 'pointed_chart',
            'horizontal_matrix', 'ordinate_scale', 'coefficients', 'height_bound',
            'status', 'finite_curve_points', 'primitive_square_hits', 'infinity_checked')
    for name, data in [('initial', initial), ('adaptive', adaptive)]:
        arms[name] = []
        for row in data['charts']:
            r = row['search']
            chart = {k: r[k] for k in keys}
            chart['input'] = {'curve': r['input']['curve']}
            arms[name].append({'chart_number': row['index'] + 1,
                               'centre': row['centre'], 'search': chart})
    checkpoint(INPUT, {'schema': 'elliptic-curves.mw16-rank27-visibility-input.v1',
                       'original_sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
                       'point_proof': {k: point_proof[k] for k in
                                       ('discovery_curve', 'discovery_points', 'rank_certificate',
                                        'minimal_curve', 'points', 'parameter')},
                       'initial_points': initial['final_state']['state']['reductions']['points'],
                       'generic_points': initial['generic_points'], 'arms': arms})


def locate(record, point):
    # PARI completes a whole box or reports a timeout; it has no prefix field.
    # Use only exact geometry here, and check retained witnesses separately.
    geometry = {k: v for k, v in record.items() if k != 'height_bound'}
    result = obs.point_visibility(geometry, point)
    result.pop('claim_boundary')
    result['returned'] = any(obs.point(q) == point for q in record['finite_curve_points'])
    coordinate = result.get('coordinate')
    result['square_hit_recorded'] = coordinate is not None and any(
        list(map(str, hit[:2])) == coordinate for hit in record['primitive_square_hits'])
    h = result['minimum_affine_height']
    result['inside_declared_affine_box'] = h is not None and h <= record['height_bound']
    return result


def expected(data=None):
    data = data or cert.read(INPUT)
    proof = data['point_proof']
    model = tuple(map(cert.F, proof['discovery_curve']))
    points = [tuple(map(cert.F, p)) for p in proof['discovery_points']]
    initial = [tuple(map(cert.F, p)) for p in data['initial_points']]
    if proof['parameter'] != '-1867/270' or len(points) != 27 or initial != points[:26]:
        raise ArithmeticError('exact specialization and ordered 26-to-27 extension required')
    if data['generic_points'] != proof['discovery_points'][:16] or any(model[:3]):
        raise ArithmeticError('generic sixteen prefix and short model required')
    old = proof['rank_certificate']
    rank_proof = checked_rank(model, points, [s['prime'] for s in old['signatures']],
                              old['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(rank_proof)) != old:
        raise ArithmeticError('independent 27-point proof differs')
    inv = cert.weierstrass_invariants(tuple(map(cert.F, proof['minimal_curve'])))
    a1, _, a3, _, _ = map(cert.F, proof['minimal_curve'])
    for p, raw in zip(points, proof['points'], strict=True):
        x = p[0] - inv['b2']/12
        transported = (x, p[1] - (a1*x+a3)/2)
        if transported != tuple(map(cert.F, raw)) or not cert.is_on_weierstrass_curve(
                tuple(map(cert.F, proof['minimal_curve'])), transported):
            raise ArithmeticError('minimal model point transport differs')
    P = points[-1]
    audits, summaries = {}, {}
    for name, count in [('initial', 43), ('adaptive', 301)]:
        charts = data['arms'][name]
        if len(charts) != count or [r['chart_number'] for r in charts] != list(range(1, count+1)):
            raise ArithmeticError('complete declared roster required')
        rows = []
        for row in charts:
            record = row['search']
            if tuple(map(cert.F, record['input']['curve'])) != model or record['height_bound'] != 100000:
                raise ArithmeticError('curve or height differs')
            if record['status'] != 'bounded_search_complete':
                raise ArithmeticError('complete historical boxes required')
            c = row['centre']; word = c['representative']
            basis = initial if name == 'adaptive' else initial[:16]
            centre = linear_combination(model, basis, word)
            if centre != obs.point(record['base_point']):
                raise ArithmeticError('exact centre combination differs')
            parity = sum((int(v) % 2) << j for j, v in enumerate(word))
            declared = c['parity'] if name == 'adaptive' else c['mask']
            if parity != declared:
                raise ArithmeticError('centre parity differs')
            if name == 'adaptive' and parity != c['generic_mask'] | (c['quotient_word'] << 16):
                raise ArithmeticError('generic/quotient splitting differs')
            for sign in (1, -1):
                p = (P[0], sign*P[1])
                located = locate(record, p)
                partner = short_add(model, centre, (p[0], -p[1]))
                partner_location = locate(record, partner)
                if located['coordinate'] != partner_location['coordinate']:
                    raise ArithmeticError('pointed involution does not preserve coordinate')
                rows.append({'chart_number': row['chart_number'], 'sign': sign,
                             'parity': parity, **located,
                             'involution_partner': list(map(str, partner))})
        finite = [r for r in rows if r['minimum_affine_height'] is not None]
        minimum = min(finite, key=lambda r: r['minimum_affine_height'])
        visible = [r for r in rows if r['inside_declared_affine_box']]
        if any(not r['returned'] or not r['square_hit_recorded'] for r in visible):
            raise ArithmeticError('visible point lacks historical output witness')
        audits[name] = rows
        summaries[name] = {'chart_count': count, 'signed_point_evaluations': len(rows),
                           'minimum': {k: minimum[k] for k in
                                       ('chart_number', 'sign', 'minimum_affine_height', 'coordinate')},
                           'visible_signed_representatives': len(visible)}
    win = summaries['adaptive']['minimum']
    if win['chart_number'] != 86 or win['minimum_affine_height'] != 79466:
        raise ArithmeticError('historical discovery witness differs')
    return {'schema': 'elliptic-curves.mw16-rank27-visibility.v1', 'status': 'PASS',
            'input_sha256': cert.hashed(INPUT), 'rank_lower_bound': 27,
            'generic_rank': 16, 'displayed_quotient_rank': 11,
            'point': proof['discovery_points'][-1], 'minimal_model_point': proof['points'][-1],
            'summaries': summaries, 'signed_chart_audits': audits,
            'height_ratio': str(cert.F(summaries['initial']['minimum']['minimum_affine_height'],
                                      win['minimum_affine_height'])),
            'claim_boundary': 'Exact visibility of this point, its negative, and the chart-dependent '
            'partners C-P and C+P only. No exhaustive quotient-coset visibility, causal comparison '
            'of centre selection alone, new rank increase, upper bound, or global covering claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'build', 'check'])
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
        return
    result = expected()
    if args.stage == 'build':
        if OUTPUT.exists():
            raise FileExistsError('preserve visibility certificate')
        checkpoint(OUTPUT, result)
    elif result != cert.read(OUTPUT):
        raise ArithmeticError('visibility certificate differs')
    print(json.dumps(result['summaries'], sort_keys=True))


if __name__ == '__main__':
    main()
