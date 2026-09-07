#!/usr/bin/env python3
"""Post-batch exact point proofs and pinned catalogue comparison for wide R17."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DATABASE = ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json'
PREVIOUS = ART/'prospective_mw16_next12_results_v1.json'

def sources():
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
        Path(__file__).resolve(), Path(cert.__file__).resolve(), Path(spec.__file__).resolve(), spec.ATLAS,
        ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',
        ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')}

def verify(row, families, projection, previous):
    model = tuple(map(cert.F, row['curve']))
    points = [tuple(map(cert.F, p)) for p in row['points']]
    proof = row['rank_certificate']
    actual = cert.checked_rank(model, points, [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
    if json.dumps(actual, sort_keys=True) != json.dumps(proof, sort_keys=True) or row['rank_lower_bound'] != len(points):
        raise ArithmeticError('independent point proof differs')
    original, generic = spec.specialize(families[row['family']], row['parameter'])
    u = cert.F(row['family_to_curve_scale_u'])
    if not u or model != (cert.F(0), cert.F(0), cert.F(0), original[3]/u**4, original[4]/u**6):
        raise ArithmeticError('family model transport differs')
    if tuple((x/u**2, y/u**3) for x, y in generic) != tuple(tuple(map(cert.F, p)) for p in row['generic_points']):
        raise ArithmeticError('generic point transport differs')
    matches = [r['id'] for r in projection if cert.isomorphic(model, r['ainvs'])]
    old = [r['address'] for r in previous if cert.isomorphic(model, r['curve'])]
    if matches != row['icarm_matches'] or old != row['previous_matches']:
        raise ArithmeticError('post-batch Q-isomorphism comparison differs')

def build(directory, output):
    if output.exists():
        raise FileExistsError('preserve wide R17 certificate')
    ledger = cert.read(directory/'point-ledger.json')
    if ledger['status'] != 'COMPLETE_FIXED_BATCH_ATTEMPTS' or len(ledger['rows']) != 24 or any(r['status'] == 'PENDING' for r in ledger['rows']):
        raise ArithmeticError('catalogue access requires terminal fixed24 batch')
    database = cert.read(DATABASE)
    projection = [{'id': r['id'], 'ainvs': r['ainvs']} for r in database['curves']]
    old = cert.read(PREVIOUS)
    previous = list(old['previous_equations']) + [
        {'address': PREVIOUS.name+':'+r['family']+':'+r['parameter'], 'curve': r['curve']} for r in old['curves']]
    families = {f['family']: f for f in cert.read(spec.ATLAS)['families']}
    rows, missing = [], []
    for entry in ledger['rows']:
        if 'result_path' not in entry:
            missing.append(entry); continue
        path = ROOT/entry['result_path']
        if cert.hashed(path) != entry['result_sha256']:
            raise ArithmeticError('terminal point input changed')
        data = cert.read(path)
        cloud_path = ART/f"compact_r17_wide_recorded_mod2_{entry['family']}_{entry['index']:02}_v1.json"
        cloud = cert.read(cloud_path)
        if cloud['status'] != 'COMPLETE_DECLARED_FINITE_AUDIT' or cloud['input_sha256'] != cert.hashed(path) or cloud['curve'] != data['curve']:
            raise ArithmeticError('complete-cloud certificate input differs')
        model = tuple(map(cert.F, data['curve']))
        row = {'family': data['family'], 'parameter': data['parameter'], 'curve': data['curve'],
            'points': cloud['independent_points'], 'generic_points': data['generic_points'],
            'family_to_curve_scale_u': data['family_to_curve_scale_u'], 'rank_certificate': cloud['rank_certificate'],
            'rank_lower_bound': cloud['rank_lower_bound'], 'icarm_matches': [r['id'] for r in projection if cert.isomorphic(model, r['ainvs'])],
            'previous_matches': [r['address'] for r in previous if cert.isomorphic(model, r['curve'])],
            'completed_charts': len(data['charts']), 'search_status': data['status'], 'supervision_status': entry['status'],
            'discovery_witness': {'path': str(path.relative_to(ROOT)), 'sha256': cert.hashed(path)},
            'complete_cloud_certificate': {'path': str(cloud_path.relative_to(ROOT)), 'sha256': cert.hashed(cloud_path)}}
        verify(row, families, projection, previous)
        rows.append(row)
        print('CERTIFIED R17 H4096', row['family'], row['parameter'], 'rank >=', row['rank_lower_bound'], 'catalogue', row['icarm_matches'], flush=True)
    pairs = [[j, i] for i, row in enumerate(rows) for j, other in enumerate(rows[:i]) if cert.isomorphic(row['curve'], other['curve'])]
    cert.write(output, {'schema': 'elliptic-curves.compact-r17-wide-results.v1', 'sources': sources(),
        'curves': rows, 'missing_measurements': missing, 'within_batch_isomorphic_pairs': pairs,
        'point_ledger_sha256': cert.hashed(directory/'point-ledger.json'),
        'selection_protocol_sha256': cert.hashed(directory/'protocol.json'), 'point_protocol_sha256': cert.hashed(directory/'point-protocol.json'),
        'previous_equations': previous, 'previous_source': {'path': str(PREVIOUS.relative_to(ROOT)), 'sha256': cert.hashed(PREVIOUS)},
        'catalogue': {'url': 'https://elliptic-rank.icarm.cloud/database.json', 'raw_sha256': cert.hashed(DATABASE),
                      'curve_count': len(projection), 'equations': projection, 'acknowledgement': 'ICARM, supported by NSF Grant DMS2425401'},
        'claim_boundary': 'Exact independent-point lower bounds on the fixed24-address R17 batch. Catalogue and prior-equation comparisons occur after terminal attempts. No exact rank, global minimality, conductor or universal novelty follows. All incomplete searches remain censored.'})

def check(path):
    data = cert.read(path)
    if data['sources'] != sources():
        raise ArithmeticError('wide R17 certificate source differs')
    families = {f['family']: f for f in cert.read(spec.ATLAS)['families']}
    for row in data['curves']:
        verify(row, families, data['catalogue']['equations'], data['previous_equations'])
        print('REPLAYED R17 H4096 CERTIFICATE', row['family'], row['parameter'], 'rank >=', row['rank_lower_bound'], flush=True)
    pairs = [[j, i] for i, row in enumerate(data['curves']) for j, other in enumerate(data['curves'][:i]) if cert.isomorphic(row['curve'], other['curve'])]
    if pairs != data['within_batch_isomorphic_pairs']:
        raise ArithmeticError('within-batch deduplication differs')

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--directory', type=Path)
    g.add_argument('--check', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    check(a.check) if a.check else build(a.directory.resolve(), a.output)
