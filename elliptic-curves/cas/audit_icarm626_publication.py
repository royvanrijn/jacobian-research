#!/usr/bin/env python3
"""Replay a frozen publication supplement; no network, search, or rank promotion."""
import argparse
from collections import Counter
from pathlib import Path

import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
INDEX = ART / 'new_high_rank_curve_index_v22.json'
PROOF = ART / 'small_conductor_rank22_proof_v1.json'
SNAPSHOT = ART / 'icarm626_publication_snapshot_v1.json'
OUT = ART / 'icarm626_publication_audit_v1.json'


def result():
    inventory, proof, snapshot = map(cert.read, (INDEX, PROOF, SNAPSHOT))
    curves = snapshot['curves']
    if len(curves) != snapshot['count'] or len({r['id'] for r in curves}) != 626:
        raise ArithmeticError('expected 626 distinct catalogue entries')
    if len(inventory['curves']) != 201:
        raise ArithmeticError('expected frozen 201-curve inventory')
    public = snapshot['submitted_curve']
    model = tuple(map(cert.F, proof['integral_model']))
    if (public['id'] != 626 or tuple(map(cert.F, public['ainvs'])) != model
            or public['rank_lower_bound'] != proof['rank_lower_bound']
            or int(public['conductor']) != int(proof['conductor'])):
        raise ArithmeticError('public entry differs from local certified result')
    points = {tuple(map(cert.F, p)) for p in public['points']}
    if points != {tuple(map(cert.F, p)) for p in proof['integral_points']}:
        raise ArithmeticError('submitted point set differs from certified point set')
    if not all(cert.is_on_weierstrass_curve(model, p) for p in points):
        raise ArithmeticError('public point off curve')
    for key in ('ainvs', 'rank_lower_bound', 'conductor'):
        if next(r for r in curves if r['id'] == 626)[key] != public[key]:
            raise ArithmeticError('single-curve and database downloads disagree')
    byj = {}
    for row in curves:
        inv = cert.weierstrass_invariants(tuple(map(cert.F, row['ainvs'])))
        byj.setdefault(inv['c4']**3 / inv['discriminant'], []).append(row)
    rows = []
    for row in inventory['curves']:
        inv = cert.weierstrass_invariants(tuple(map(cert.F, row['curve'])))
        j = inv['c4']**3 / inv['discriminant']
        if j != cert.F(row['j_invariant']):
            raise ArithmeticError('inventory invariant mismatch')
        matches = sorted(q['id'] for q in byj.get(j, [])
                         if cert.isomorphic(row['curve'], q['ainvs']))
        rows.append({k: row[k] for k in
                     ('id', 'family', 'parameter', 'rank_lower_bound', 'rank_provenance')}
                    | {'catalogue_matches': matches})
    unmatched = [r for r in rows if not r['catalogue_matches']]
    relevant = [r for r in curves if r['rank_lower_bound'] >= 22]
    paths = [Path(__file__).resolve(), Path(cert.__file__),
             ROOT / 'elliptic-curves/cas/elliptic_candidate_record.py', INDEX, PROOF, SNAPSHOT]
    return {
        'schema': 'elliptic-curves.icarm626-publication.v1', 'status': 'PASS',
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
        'catalogue_intake': snapshot['intake'], 'catalogue_count': len(curves),
        'inventory_count': len(rows), 'publication_status': rows,
        'unmatched_count': len(unmatched),
        'unmatched_by_rank_lower_bound': dict(sorted(Counter(
            str(r['rank_lower_bound']) for r in unmatched).items())),
        'optional_rank27_submission_shortlist': [r['id'] for r in unmatched
                                               if r['rank_lower_bound'] >= 27],
        'submitted': {'local_id': proof['curve_id'], 'icarm_id': 626,
                      'url': 'https://elliptic-rank.icarm.cloud/curve/626',
                      'equation_conductor_and_point_set_match': True},
        'smaller_recorded_rank_at_least22_conductor_ids': sorted(
            r['id'] for r in relevant if r['conductor'] is not None
            and int(r['conductor']) < int(proof['conductor'])),
        'missing_rank_at_least22_conductor_ids': sorted(
            r['id'] for r in relevant if r['conductor'] is None),
        'claim_boundary': 'Publication metadata and exact Q-isomorphism comparison only. '
            'V22 equations, points, ranks and frozen novelty claims are unchanged. '
            'The submitted point set is identical to the existing proof; no new rank '
            'certificate, independent implementation, universal novelty, record, or '
            'first-discovery priority is claimed. The shortlist is editorial advice.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    expected = result()
    if args.check:
        if cert.read(OUT) != expected:
            raise ArithmeticError('publication supplement differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve frozen publication supplement')
        cert.write(OUT, expected)
    print('ICARM626 PUBLICATION PASS:', expected['unmatched_count'], 'unmatched curves')
