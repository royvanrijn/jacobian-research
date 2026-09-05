#!/usr/bin/env python3
"""Replay the retained586-row comparison of every curve in the36-row inventory."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1'
INDEX = ART/'new_high_rank_curve_index_v2.json'
DATABASE = LOCAL/'database.json'
COMPARISON = LOCAL/'comparison.json'
CONDUCTOR = ART/'small_conductor_rank22_proof_v1.json'

def compute():
    database = cert.read(DATABASE)['curves']
    inventory = cert.read(INDEX)['curves']
    original = cert.read(COMPARISON)
    if len(database) != 586 or len(inventory) != 36 or original['catalogue_sha256'] != cert.hashed(DATABASE):
        raise ArithmeticError('fixed cohort binding differs')
    rows = []
    for row in inventory:
        matches = [r['id'] for r in database if cert.isomorphic(row['curve'], r['ainvs'])]
        rows.append({k: row[k] for k in ('family', 'id', 'parameter')})
        rows[-1]['matches'] = matches
    N = int(cert.read(CONDUCTOR)['conductor'])
    smaller = sorted([{'id': r['id'], 'conductor': r['conductor']} for r in database
                      if r['rank_lower_bound'] >= 22 and r['conductor'] and int(r['conductor']) < N],
                     key=lambda r: int(r['conductor']))
    missing = [r['id'] for r in database if r['rank_lower_bound'] >= 22 and not r['conductor']]
    expected = {'catalogue_sha256': cert.hashed(DATABASE), 'curves': rows,
                'smaller_rank22_conductors': smaller, 'rank22_missing_conductor': missing}
    if expected != original or any(r['matches'] for r in rows):
        raise ArithmeticError('fresh catalogue comparison differs')
    paths = (Path(__file__).resolve(), INDEX, DATABASE, COMPARISON, CONDUCTOR,
             Path(cert.__file__).resolve(), ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    return {'schema': 'elliptic-curves.refreshed-inventory-replay.v1', 'status': 'PASS_EXACT_Q_ISOMORPHISM_COMPARISON',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}, 'comparison': expected,
            'claim_boundary': 'No Q-isomorphism matches for these36 equations in the pinned586 snapshot. Exact point proofs remain in their independent certificates. This finite catalogue comparison does not establish universal novelty.'}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--output', type=Path)
    g.add_argument('--check', type=Path)
    a = p.parse_args()
    result = compute()
    if a.check:
        if cert.read(a.check) != result:
            raise ArithmeticError('fresh comparison replay artifact differs')
    else:
        if a.output.exists():
            raise FileExistsError('preserve comparison replay')
        cert.write(a.output, result)
    print('REPLAYED REFRESHED INVENTORY:36 curves;586 catalogue equations; no matches', flush=True)
