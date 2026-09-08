#!/usr/bin/env python3
"""Versioned offline ICARM mirror and current publication view over frozen V22."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
CURRENT = ROOT / 'elliptic-curves/data/icarm_current.json'
INDEX = ART / 'new_high_rank_curve_index_v22.json'
BASELINE = ROOT / 'artifacts/local/elliptic-curves/inventory200-current-catalogue-v1/database.json'


def relative(path):
    return str(path.relative_to(ROOT))


def load_catalogue():
    """Load all public fields; their reported ranks are not local proof authority."""
    manifest = cert.read(CURRENT)
    raw = gzip.decompress((ROOT / manifest['snapshot']).read_bytes())
    if hashlib.sha256(raw).hexdigest() != manifest['raw_sha256']:
        raise ArithmeticError('current catalogue hash mismatch')
    return json.loads(raw)


def load_inventory():
    """Join current publication status without changing discovery points or ranks."""
    manifest = cert.read(CURRENT)
    for key, hashkey in [('inventory', 'inventory_sha256'), ('publication', 'publication_sha256')]:
        if cert.hashed(ROOT / manifest[key]) != manifest[hashkey]:
            raise ArithmeticError('current inventory/publication binding differs')
    old = cert.read(ROOT / manifest['inventory'])
    overlay = cert.read(ROOT / manifest['publication'])
    if overlay['catalogue_sha256'] != manifest['raw_sha256']:
        raise ArithmeticError('publication catalogue binding differs')
    updates = {r['id']: r for r in overlay['curves']}
    rows = []
    for r in old['curves']:
        updated = dict(r)
        updated['frozen_v22_publication'] = {
            k: r[k] for k in ('current_catalogue_matches', 'novelty_status')}
        updated.update(updates[r['id']])
        rows.append(updated)
    return rows


def build(intake_dir, check=False):
    raw_path, metadata_path = intake_dir / 'database.json', intake_dir / 'metadata.json'
    raw, metadata = raw_path.read_bytes(), cert.read(metadata_path)
    digest = hashlib.sha256(raw).hexdigest()
    if digest != metadata['sha256']:
        raise ArithmeticError('download metadata mismatch')
    db, old, inventory = json.loads(raw), cert.read(BASELINE), cert.read(INDEX)
    byid = {r['id']: r for r in db['curves']}
    if len(byid) != db['count'] or len(byid) != len(db['curves']):
        raise ArithmeticError('catalogue count or duplicate IDs')
    stem = f'icarm_catalogue_{db["count"]}_{digest[:12]}'
    snapshot, publication, delta = [ART / (stem + suffix) for suffix in
                                     ('.json.gz', '_publication.json', '_delta.json')]
    byj, point_count = {}, 0
    for r in db['curves']:
        model = tuple(map(cert.F, r['ainvs']))
        inv = cert.weierstrass_invariants(model)
        if not inv['discriminant'] or inv['discriminant'] != cert.F(r['discriminant']):
            raise ArithmeticError(f'bad discriminant: {r["id"]}')
        for p in r['points']:
            if not cert.is_on_weierstrass_curve(model, tuple(map(cert.F, p))):
                raise ArithmeticError(f'point off curve: {r["id"]}')
            point_count += 1
        byj.setdefault(inv['c4']**3/inv['discriminant'], []).append(r)
    rows = []
    for r in inventory['curves']:
        inv = cert.weierstrass_invariants(tuple(map(cert.F, r['curve'])))
        matches = sorted(q['id'] for q in byj.get(inv['c4']**3/inv['discriminant'], [])
                         if cert.isomorphic(r['curve'], q['ainvs']))
        rows.append({'id': r['id'], 'current_catalogue_matches': matches,
                     'novelty_status': 'CATALOGUE_MATCH' if matches else 'UNMATCHED_IN_CURRENT_SNAPSHOT',
                     'publication_catalogue_sha256': digest})
    unmatched = {r['id'] for r in rows if not r['current_catalogue_matches']}
    overlay = {'catalogue_sha256': digest, 'inventory_sha256': cert.hashed(INDEX),
               'curves': rows, 'unmatched_count': len(unmatched),
               'unmatched_by_local_rank': dict(sorted(Counter(
                   str(r['local_search_rank_lower_bound']) for r in inventory['curves']
                   if r['id'] in unmatched).items())),
               'claim_boundary': 'Publication metadata only. No point, rank or discovery provenance changes.'}
    before = {r['id']: r for r in old['curves']}
    changes = [{'id': i, 'changed_fields': sorted(k for k in set(before[i]) | set(byid[i])
                 if before[i].get(k) != byid[i].get(k))} for i in sorted(set(before) & set(byid))]
    changes = [r for r in changes if r['changed_fields']]
    benchmarks = []
    for rank in range(22, 28):
        def best(catalogue):
            candidates = [r for r in catalogue if r['rank_lower_bound'] >= rank and r['conductor']]
            q = min(candidates, key=lambda r: int(r['conductor']))
            return {'id': q['id'], 'recorded_conductor': q['conductor']}
        benchmarks.append({'rank_threshold': rank, 'previous': best(old['curves']),
            'current': best(db['curves']), 'missing_conductors': sorted(r['id'] for r in db['curves']
                if r['rank_lower_bound'] >= rank and r['conductor'] is None)})
    sources = [Path(__file__).resolve(), Path(cert.__file__),
               ROOT / 'elliptic-curves/cas/elliptic_candidate_record.py', INDEX, BASELINE, metadata_path]
    report = {'schema': 'elliptic-curves.catalogue-refresh.v1', 'status': 'PASS',
              'intake': metadata, 'baseline_count': old['count'], 'current_count': db['count'],
              'sources': {relative(p): cert.hashed(p) for p in sources},
              'added_ids': sorted(set(byid)-set(before)), 'removed_ids': sorted(set(before)-set(byid)),
              'existing_entry_changes': changes, 'point_memberships_checked': point_count,
              'equation_discriminants_checked': len(byid), 'conductor_benchmarks': benchmarks,
              'claim_boundary': 'Exact equation, point membership and Q-isomorphism audit. '
                  'Public ranks and conductors remain reported metadata unless separately certified. '
                  'No search, point-independence claim for the whole catalogue or universal record claim.'}
    manifest = {'schema': 'elliptic-curves.current-catalogue.v1', 'snapshot': relative(snapshot),
                'raw_sha256': digest, 'retrieved_at': metadata['retrieved_at'], 'count': db['count'],
                'inventory': relative(INDEX), 'inventory_sha256': cert.hashed(INDEX),
                'publication': relative(publication), 'delta': relative(delta),
                'intake_directory': relative(intake_dir),
                'usage': 'Post-discovery deduplication and explicitly retrospective controls. '
                    'Prospective selection continues to use its own frozen inputs.'}
    blob = (json.dumps(overlay, sort_keys=True, indent=2)+'\n').encode()
    manifest['publication_sha256'] = hashlib.sha256(blob).hexdigest()
    if check:
        if gzip.decompress(snapshot.read_bytes()) != raw:
            raise ArithmeticError('archived raw snapshot differs')
        if publication.read_bytes() != blob or cert.read(delta) != report or cert.read(CURRENT) != manifest:
            raise ArithmeticError('current refresh outputs differ')
    else:
        for p, content in [(snapshot, gzip.compress(raw, mtime=0)), (publication, blob),
                           (delta, (json.dumps(report, sort_keys=True, indent=2)+'\n').encode())]:
            if p.exists() and p.read_bytes() != content:
                raise FileExistsError('preserve frozen refresh: ' + str(p))
            if not p.exists():
                p.write_bytes(content)
        cert.write(CURRENT, manifest)
    joined = load_inventory()
    if len(joined) != len(inventory['curves']):
        raise ArithmeticError('inventory join lost a discovery')
    for original, current in zip(inventory['curves'], joined):
        if any(original[k] != current[k] for k in ('curve', 'points', 'rank_certificate',
                'rank_lower_bound', 'local_search_rank_lower_bound', 'rank_provenance')):
            raise ArithmeticError('refresh changed research proof data')
    print('CURRENT CATALOGUE PASS:', db['count'], 'curves;', point_count,
          'point memberships;', len(unmatched), 'unmatched inventory equations')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--intake', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    intake_dir = args.intake or ROOT / cert.read(CURRENT)['intake_directory']
    build(intake_dir.resolve(), args.check)
