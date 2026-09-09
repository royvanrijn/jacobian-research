#!/usr/bin/env python3
"""Fresh R17 selector for the autonomous rank hunter.

The selector never imposes a score cutoff. It freezes the full retained 6144-row
population plus an exact repository exclusion snapshot, then mixes score-frontier,
low-discrepancy and deterministic-hash exploration lanes.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import certify_compact_r17_candidates as cert
from autonomous_rank_hunter_policy import FAMILIES, height_stratum, lane_for_dispatch
from v3_warm_support import atomic, read, require, sha

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
POOLS = tuple(LOCAL / n / 'result.json' for n in (
    'r17-retained-extended-primes-v1',
    'r17-discarded-shards-extended-v1',
    'r17-retention512-extended-v1'))
FIELDS = ('family', 'parameter', 'model', 'numerator', 'denominator',
          'combined_selection_units', 'combined_good')


def order_key(row):
    return (-row['combined_selection_units'], -row['combined_good'],
            row['denominator'], row['numerator'])


def j_key(model):
    inv = cert.weierstrass_invariants(tuple(map(F, model)))
    require(bool(inv['discriminant']), 'singular candidate')
    return inv['c4'] ** 3 / inv['discriminant']


def _visit(value, models, addresses):
    if isinstance(value, dict):
        for key in ('curve', 'model', 'a_invariants'):
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
            _visit(child, models, addresses)
    elif isinstance(value, list):
        for child in value:
            _visit(child, models, addresses)


def freeze_snapshot(folder):
    """Freeze selection inputs once. Point outcomes never mutate this snapshot."""
    path = folder / 'selection-snapshot.json'
    if path.exists():
        return read(path)
    pool = []
    for source in POOLS:
        require(source.exists(), 'missing retained R17 pool: ' + str(source))
        for raw in read(source)['rows']:
            row = {k: raw[k] for k in FIELDS}
            require(row['family'] in FAMILIES, 'unexpected R17 family')
            require(F(row['parameter']) == F(row['numerator'], row['denominator']), 'parameter differs')
            require(row['denominator'] > 0 and max(abs(row['numerator']), row['denominator']) <= 4096,
                    'retained population boundary differs')
            pool.append(row)
    require(len(pool) == 6144, 'fixed retained pool count differs')
    require(len({(r['family'], str(F(r['parameter']))) for r in pool}) == len(pool), 'duplicate pool address')

    # Generated result packets + rendered inventory are the known-curve universe.
    # Do NOT scan retained candidate-pool result files here: that would classify the
    # unsearched 6144-row population itself as already known.
    paths = set(ART.glob('*.json'))
    db = ROOT / 'elliptic-curves/data/research_curves/database.json'
    if db.exists():
        paths.add(db)
    models, addresses = set(), set()
    for p in sorted(paths):
        try:
            _visit(read(p), models, addresses)
        except (OSError, json.JSONDecodeError):
            continue
    snapshot = {
        'schema': 'autonomous-r17-selection.v1',
        'pool': pool,
        'excluded_models': sorted(models),
        'reserved_addresses': sorted(addresses),
        'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in POOLS},
        'exclusion_source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)},
        'boundary': 'All 6144 retained H4096 rows; no score cutoff. Freshness is exact Q-isomorphism/address freshness relative to this frozen repository/local snapshot.'
    }
    atomic(path, snapshot, immutable=True)
    return snapshot


def _vdc(n):
    """Base-2 van der Corput point in [0,1), used only for broad deterministic coverage."""
    x, f = 0.0, 0.5
    n = int(n) + 1
    while n:
        x += f * (n & 1)
        n >>= 1
        f *= 0.5
    return x


def _is_excluded(row, model_index, reserved):
    key = (row['family'], str(F(row['parameter'])))
    if key in reserved:
        return True
    model = tuple(map(F, row['model']))
    jk = j_key(model)
    return any(cert.isomorphic(model, tuple(map(F, q))) for q in model_index.get(jk, ()))


def _family_order(attempted, dispatch_index):
    prior = {'074d9': 1.05, '07ca9': 1.00, '08234': 0.88,
             '08f72': 0.98, '103b2': 1.12, '11952': 1.08}
    rotation = int(dispatch_index) % len(FAMILIES)
    rotated = FAMILIES[rotation:] + FAMILIES[:rotation]
    return sorted(rotated, key=lambda f: (attempted.get(f, 0) / prior[f], rotated.index(f)))


def select_next(snapshot, campaign_rows, dispatch_index):
    lane = lane_for_dispatch(dispatch_index)
    requested_stratum, mode = lane.split('_', 1)
    attempted = {f: 0 for f in FAMILIES}
    used_addresses = set()
    used_models = []
    for row in campaign_rows:
        attempted[row['family']] = attempted.get(row['family'], 0) + 1
        used_addresses.add((row['family'], str(F(row['parameter']))))
        used_models.append(row['model'])

    model_index = {}
    for q in [*snapshot['excluded_models'], *used_models]:
        try:
            model_index.setdefault(j_key(q), []).append(q)
        except (ArithmeticError, ValueError, TypeError, ZeroDivisionError):
            pass
    reserved = set(map(tuple, snapshot['reserved_addresses'])) | used_addresses

    for family in _family_order(attempted, dispatch_index):
        candidates = [r for r in snapshot['pool'] if r['family'] == family and
                      height_stratum(r) == requested_stratum and not _is_excluded(r, model_index, reserved)]
        if not candidates:
            continue
        ordered = sorted(candidates, key=order_key)
        if mode == 'score':
            chosen = ordered[0]
        elif mode == 'spread':
            chosen = ordered[min(len(ordered)-1, int(_vdc(dispatch_index) * len(ordered)))]
        elif mode == 'hash':
            salt = f'autonomous-r17-v1/{dispatch_index}/{family}/{requested_stratum}/'
            chosen = min(ordered, key=lambda r: hashlib.sha256(
                (salt + str(F(r['parameter']))).encode()).digest())
        else:
            raise ValueError('unknown selection lane: ' + lane)
        result = dict(chosen)
        result.update(selection_lane=lane, selection_dispatch_index=int(dispatch_index),
                      stratum=requested_stratum,
                      score_position=ordered.index(chosen)+1,
                      eligible_in_family_stratum=len(ordered))
        return result

    alternate = 'high' if requested_stratum == 'low' else 'low'
    fake_index = next(i for i in range(dispatch_index + 1, dispatch_index + 25)
                      if lane_for_dispatch(i).startswith(alternate + '_'))
    return select_next(snapshot, campaign_rows, fake_index)
