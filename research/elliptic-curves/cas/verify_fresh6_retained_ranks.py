#!/usr/bin/env python3
"""Replay retained subgroup witnesses for the three six-fibre cohorts.

Select sufficient primes from saved binary rows, then recompute those finite
quotients and exact point identities. No prime search, CVP, point search,
factorization, persistent cache or original worker restart. This shares the
established finite-group implementation; it is not a new independent algorithm.
Full protocol/map/visibility claims require their separate historical replays.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from compact_atlas_specialization import ATLAS, specialize
from certify_compact_r17_candidates import isomorphic
from memory_rank_certificate import checked_rank

ROOT = Path(__file__).resolve().parents[2]
WITNESS = ROOT/'artifacts/generated-results/elliptic-curves/fresh6_retained_rank_witnesses_v1.json'
WITNESS_SHA256 = '75a3829bc7e6040cc88577446dfbfed792b7929a13c662dd3176177ccac1a5a4'
SCHEMA = 'elliptic-curves.fresh6-retained-rank-witnesses.v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.loads(json.dumps(value))


def select_signatures(signatures, width):
    """Greedily retain rank-increasing prime blocks; verify arithmetic later.

    At most width primes are needed. This does not claim minimum cardinality,
    and the saved rows serve only to propose a witness, never to prove it.
    """
    require(type(width) is int and 0 < width <= 32, 'invalid point count')
    pivots, selected, seen = {}, [], set()
    for signature in signatures:
        prime = signature['prime']
        require(type(prime) is int and 3 <= prime <= 1000 and prime % 2,
                'outside retained small-prime protocol')
        require(prime not in seen, 'duplicate prime block')
        seen.add(prime)
        grew = False
        for row in signature['rows']:
            require(len(row) == width and all(type(x) is int and x in (0, 1) for x in row),
                    'malformed binary row')
            bits = sum(x << i for i, x in enumerate(row))
            while bits:
                pivot = bits.bit_length()-1
                if pivot in pivots:
                    bits ^= pivots[pivot]
                else:
                    pivots[pivot] = bits
                    grew = True
                    break
        if grew:
            selected.append(signature)
    require(len(pivots) == width, 'saved rows do not span the point columns')
    return selected


def replay_row(row, family):
    packet = row['packet']
    model = tuple(map(Fraction, packet['curve']))
    points = tuple(tuple(map(Fraction, p)) for p in packet['points'])
    rank, proof = packet['rank_lower_bound'], packet['proof']
    require(rank == len(points) == proof['rank_lower_bound'], 'rank header differs')
    original, generic = specialize(family, row['parameter'])
    require(model == original and points[:17] == generic, 'native family or section prefix differs')
    selected = select_signatures(proof['signatures'], rank)
    require(selected == proof['signatures'], 'packet still contains redundant prime blocks')
    actual = checked_rank(model, points, [s['prime'] for s in selected],
                          proof['no_rational_2_torsion_prime'])
    require(canonical(actual) == proof, 'recomputed finite witness differs')
    return model, rank


def verify(data, atlas):
    require(data['schema'] == SCHEMA, 'unexpected witness schema')
    rows = data['rows']
    require(len(rows) == 37 and len({r['id'] for r in rows}) == 37,
            'incomplete or duplicate cohort endpoint roster')
    families = {f['family']: f for f in atlas['families']}
    cohorts = {c: [] for c in ('first', 'second', 'lowheight')}
    claims = set()
    results = []
    for row in rows:
        model, rank = replay_row(row, families[row['family']])
        claims.update(row['claims'])
        if row.get('seed_cohort'):
            cohort = cohorts[row['seed_cohort']]
            require(not any(isomorphic(model, old) for old in cohort),
                    'seed cohort contains rationally isomorphic curves')
            cohort.append(model)
        results.append({'id': row['id'], 'rank_lower_bound': rank,
                        'verified_primes': len(row['packet']['proof']['signatures'])})
    require(all(len(models) == 6 for models in cohorts.values()), 'incomplete seed cohort')
    require(claims == set(data['rank_claims']) and len(claims) == 17, 'rank-claim coverage differs')
    return results


def at_path(data, selector):
    for key in selector:
        data = data[key]
    return data


def check_sources(data):
    """Optional export/provenance check, without repeating arithmetic."""
    sources = {}
    for name, digest in data['source_inputs'].items():
        path = ROOT/name
        raw = path.read_bytes()  # Missing evidence never triggers reconstruction.
        require(hashlib.sha256(raw).hexdigest() == digest, 'original input changed: '+name)
        sources[name] = json.loads(raw)
    for row in data['rows']:
        raw = at_path(sources[row['source_file']], row['source_selector'])
        proof = {**raw['proof'], 'signatures': select_signatures(
            raw['proof']['signatures'], raw['rank_lower_bound'])}
        expected = {k: raw[k] for k in ('curve', 'points', 'rank_lower_bound')}
        expected['proof'] = proof
        require(row['packet'] == expected, 'retained rank projection differs: '+row['id'])
        identity = at_path(sources[row['identity_file']], row['identity_selector'])
        require(row['family'] == identity['family'] and row['parameter'] == identity['parameter'],
                'source identity differs: '+row['id'])
    return len(sources)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-sources', action='store_true',
                        help='Compare the projection with the pinned original exports; no arithmetic replay')
    args = parser.parse_args()
    raw = WITNESS.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == WITNESS_SHA256, 'rank witness bytes changed')
    data = json.loads(raw)
    if args.check_sources:
        print('PASS retained rank projection:', check_sources(data), 'original exports; no arithmetic')
        return
    atlas_raw = ATLAS.read_bytes()
    require(hashlib.sha256(atlas_raw).hexdigest() == data['atlas_sha256'], 'native atlas changed')
    results = verify(data, json.loads(atlas_raw))
    print('PASS', len(results), 'retained rank endpoints for 17 claims;',
          sum(r['verified_primes'] for r in results), 'small-prime blocks; no search')


if __name__ == '__main__':
    main()
