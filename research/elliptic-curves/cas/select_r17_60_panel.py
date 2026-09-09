#!/usr/bin/env python3
"""Freeze a 5+5 height-stratified, six-fibration roster before point outcomes."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import certify_compact_r17_candidates as cert
from v3_warm_support import atomic, read, require, sha

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
FAMILIES = ('074d9', '07ca9', '08234', '08f72', '103b2', '11952')
POOLS = tuple(LOCAL / n / 'result.json' for n in (
    'r17-retained-extended-primes-v1', 'r17-discarded-shards-extended-v1',
    'r17-retention512-extended-v1'))
FIELDS = ('family', 'parameter', 'model', 'numerator', 'denominator',
          'combined_selection_units', 'combined_good')


def order_key(row):
    return (-row['combined_selection_units'], -row['combined_good'],
            row['denominator'], row['numerator'])


def j_key(model):
    inv = cert.weierstrass_invariants(tuple(map(F, model)))
    require(bool(inv['discriminant']), 'singular candidate')
    return inv['c4']**3 / inv['discriminant']


def select(snapshot, per_stratum=5):
    pool = snapshot['pool']
    require(len({(r['family'], str(F(r['parameter']))) for r in pool}) == len(pool), 'duplicate pool address')
    exclude = {}
    for model in snapshot['excluded_models']:
        exclude.setdefault(j_key(model), []).append(model)
    reserved = {tuple(a) for a in snapshot['reserved_addresses']}
    selected, skipped, counts = [], [], {}
    for family in FAMILIES:
        for stratum in ('low', 'high'):
            candidates = sorted((r for r in pool if r['family'] == family and
                                 ('low' if max(abs(r['numerator']), r['denominator']) <= 1024 else 'high') == stratum), key=order_key)
            accepted = []
            for pos, row in enumerate(candidates, 1):
                key = (family, str(F(row['parameter'])))
                model = row['model']; jk = j_key(model)
                if key in reserved or any(cert.isomorphic(tuple(map(F, model)), tuple(map(F, q))) for q in exclude.get(jk, [])):
                    skipped.append({'family': family, 'parameter': key[1], 'stratum': stratum})
                    continue
                if len(accepted) < per_stratum:
                    rec = {**row, 'id': f'{family}-{stratum}-{len(accepted)+1:02d}',
                           'stratum': stratum, 'stratum_score_position': pos}
                    accepted.append(rec)
                    exclude.setdefault(jk, []).append(model)
            require(len(accepted) == per_stratum, f'not enough fresh candidates: {family}/{stratum}')
            selected.extend(accepted)
            counts[f'{family}/{stratum}'] = len(candidates)
    # Interleave strata and fibrations; results cannot determine dispatch order.
    selected.sort(key=lambda r: (int(r['id'][-2:]), r['stratum'] != 'low', r['family']))
    return {'rows': selected, 'skipped': skipped, 'stratum_pool_counts': counts}


def freeze(folder):
    snapshot_path = folder / 'selection-input.json'
    if snapshot_path.exists():
        snapshot = read(snapshot_path)
    else:
        pool = []
        for path in POOLS:
            for raw in read(path)['rows']:
                row = {k: raw[k] for k in FIELDS}
                require(F(row['parameter']) == F(row['numerator'], row['denominator']) and row['denominator'] > 0,
                        'parameter coordinates differ')
                require(max(abs(row['numerator']), row['denominator']) <= 4096, 'wrong retained population')
                pool.append(row)
        require(len(pool) == 6144, 'fixed retained pool count differs')
        paths = sorted(set(list(ART.glob('*results*.json')) +
                           list(LOCAL.glob('*pari*/protocol.json')) + [ART / 'compact_atlas_endpoints_v2.json']))
        models, addresses = set(), set()
        def visit(value):
            if isinstance(value, dict):
                for key in ('curve', 'model'):
                    q = value.get(key)
                    if isinstance(q, list) and len(q) == 5:
                        try:
                            q = tuple(map(F, q))
                            if cert.weierstrass_invariants(q)['discriminant']:
                                models.add(tuple(map(str, q)))
                        except (ValueError, TypeError, ZeroDivisionError):
                            pass
                if value.get('family') in FAMILIES and 'parameter' in value:
                    try:
                        addresses.add((value['family'], str(F(value['parameter']))))
                    except (ValueError, TypeError, ZeroDivisionError):
                        pass
                for child in value.values():
                    visit(child)
            elif isinstance(value, list):
                for child in value:
                    visit(child)
        for path in paths:
            visit(read(path))
        snapshot = {'pool': pool, 'excluded_models': sorted(models),
                    'reserved_addresses': sorted(addresses),
                    'source_snapshots': {str(p.relative_to(ROOT)): sha(p) for p in [*POOLS, *paths]},
                    'boundary': 'Fresh relative to these equation/address snapshots; outcomes are projected out before selection. No literature-wide novelty claim.'}
        atomic(snapshot_path, snapshot, immutable=True)
    result = select(snapshot)
    result.update(status='PASS_FROZEN_R17_60_SELECTION', point_searches=0,
                  input_sha256=sha(snapshot_path), selector_sha256=sha(Path(__file__)))
    atomic(folder / 'roster.json', result, immutable=True)
    require(select(read(snapshot_path)) == {k: result[k] for k in ('rows', 'skipped', 'stratum_pool_counts')}, 'selection replay differs')
    reservation = {'status': 'FROZEN_R17_60_RESERVED', 'rows': [
        {k: r[k] for k in ('id', 'family', 'parameter', 'model')} for r in result['rows']],
        'roster_sha256': sha(folder / 'roster.json')}
    atomic(ART / 'r17_60_reserved_results_v1.json', reservation, immutable=True)
    print('R17_60_SELECTION_PASS|cases=60|low=30|high=30|point_searches=0', flush=True)
    return result
