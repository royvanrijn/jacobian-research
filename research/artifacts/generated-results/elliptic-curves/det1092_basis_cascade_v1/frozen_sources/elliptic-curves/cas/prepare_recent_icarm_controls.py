#!/usr/bin/env python3
"""Bounded public-point proofs and separated seed/oracle recovery controls."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import certify_compact_r17_candidates as cert
import memory_rank_certificate as memory
from mod2_reduction_independence import (_primes_up_to,
    find_two_torsion_certificate_prime, combined_mod2_rank)
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from refresh_icarm_local_database import load_catalogue, CURRENT, ROOT, ART

OUT = ART / 'icarm_recent_controls_v1'
PROTOCOL = {
    'schema': 'elliptic-curves.public-control-protocol.v1',
    'public_ids': list(range(615, 626)), 'prime_bound': 997,
    'recovery_control_ids': list(range(615, 621)),
    'mask': 'Last public point withheld for 615,616,617,618,620. '
            '619 uses the historical local27 basis and public point26 (zero-based).',
    'purpose': 'Retrospective point-recovery calibration and smaller-conductor comparison.',
    'limits': 'One process; finite reduction primes at most997; checkpoint each public proof. '
              'No rational-point enumeration, factorization, descent or parameter scan.',
    'failure': 'Rank stays UNKNOWN unless full finite-column rank and torsion exclusion close.',
    'separation': 'seeds.json is the only worker input. oracle.json and public proofs are '
                  'validator-only data. No existing prospective protocol is changed.',
    'evaluation': 'Freeze a search algorithm and equal budgets before a trial; log every '
                  'completed/timeout box and runtime. Success requires an exactly independent '
                  '28th direction, which need not equal the withheld representative. '
                  'Report genuine-local619 separately from five synthetic deletions. '
                  'A miss does not prove rank27; a recovery is not a new curve or rank discovery.'}


def normalized(value):
    return json.loads(json.dumps(value))


def proof(model, points):
    cache = ReductionCache(MemoryFactStore())
    signatures, current_rank = [], 0
    for prime in _primes_up_to(997):
        if prime == 2:
            continue
        try:
            signature = cache.signature(model, points, prime)
        except ValueError:
            continue
        rank = combined_mod2_rank(signatures+[signature], len(points))
        if rank > current_rank:
            signatures.append(signature)
            current_rank = rank
        if rank == len(points):
            break
    rank = combined_mod2_rank(signatures, len(points))
    if rank != len(points):
        raise ArithmeticError(f'UNKNOWN: finite rank {rank} < {len(points)}')
    tp = find_two_torsion_certificate_prime(model, prime_bound=997)
    actual = memory.checked_rank(model, points, [s.prime for s in signatures], tp)
    if normalized(actual['signatures']) != normalized([asdict(s) for s in signatures]):
        raise ArithmeticError('finite signature replay disagrees')
    return actual


def put(name, value, check):
    path = OUT / name
    if path.exists():
        if cert.read(path) != normalized(value):
            raise ArithmeticError('frozen control output differs: ' + name)
    elif check:
        raise FileNotFoundError(path)
    else:
        cert.write(path, value)


def run(check=False):
    manifest = cert.read(CURRENT)
    if manifest['raw_sha256'] != '6a7ebb045f371832f54ba7509c7c802e484de14ef0fcbf0b2b10d0d019bc96d6':
        raise ArithmeticError('this control cohort is pinned to the626 snapshot')
    put('protocol.json', PROTOCOL, check)
    public = {r['id']: r for r in load_catalogue()['curves']}
    prepared = {}
    for identifier in PROTOCOL['public_ids']:
        r = public[identifier]
        model = tuple(map(cert.F, r['ainvs']))
        inv = cert.weierstrass_invariants(model)
        short = (cert.F(0), cert.F(0), cert.F(0), -inv['c4']/48, -inv['c6']/864)
        points = []
        for raw in r['points']:
            x, y = map(cert.F, raw)
            if not cert.is_on_weierstrass_curve(model, (x, y)):
                raise ArithmeticError('public point off curve')
            points.append((x+inv['b2']/12, y+(model[0]*x+model[2])/2))
        if len(points) != r['rank_lower_bound']:
            raise ArithmeticError('unexpected public point count')
        certificate = proof(short, points)
        data = {'public_id': identifier, 'public_submitter': r['submitter'],
                'public_commentary': r['commentary'], 'source_ainvs': r['ainvs'],
                'curve': list(map(str, short)), 'points': [list(map(str, p)) for p in points],
                'rank_certificate': certificate, 'proved_rank_lower_bound': len(points),
                'recorded_conductor': r['conductor'], 'conductor_proof_status': 'NOT_REPLAYED',
                'construction_status': 'PUBLIC_COMMENTARY_ONLY' if r['commentary'] else 'UNKNOWN',
                'catalogue_sha256': manifest['raw_sha256']}
        put(f'public_{identifier}.json', data, check)
        prepared[identifier] = (short, points)
        print(f'ICARM{identifier}: rank >= {len(points)} exact public-point proof PASS', flush=True)
    seeds, oracle = [], []
    old = next(r for r in cert.read(ART / 'new_high_rank_curve_index_v20.json')['curves']
               if r['id'] == 'new-20260906-188')
    for identifier in PROTOCOL['recovery_control_ids']:
        model, points = prepared[identifier]
        if identifier == 619:
            if tuple(map(cert.F, old['curve'])) != model:
                raise ArithmeticError('historical619 short model differs')
            seed = [tuple(map(cert.F, p)) for p in old['points']]
            withheld = points[26]
            kind = 'GENUINE_HISTORICAL_LOCAL27'
        else:
            seed, withheld = points[:27], points[27]
            kind = 'SYNTHETIC_PUBLIC_BASIS_DELETION'
        fullproof = proof(model, seed+[withheld])
        seedproof = memory.checked_rank(model, seed,
            [s['prime'] for s in fullproof['signatures']], fullproof['no_rational_2_torsion_prime'])
        case_id = f'icarm-{identifier}'
        seeds.append({'case_id': case_id, 'curve': list(map(str, model)),
                      'points': [list(map(str, p)) for p in seed], 'seed_rank_certificate': seedproof})
        oracle.append({'case_id': case_id, 'public_id': identifier, 'control_type': kind,
                       'withheld_point': list(map(str, withheld)), 'expected_minimum_rank': 28,
                       'full_28_point_certificate': fullproof,
                       'usage': 'Validator only; withheld coordinates must not enter worker selection.'})
    put('seeds.json', {'schema': 'elliptic-curves.point-recovery-seeds.v1', 'cases': seeds}, check)
    put('oracle.json', {'schema': 'elliptic-curves.point-recovery-oracle.v1', 'cases': oracle}, check)
    paths = [Path(__file__).resolve(), Path(cert.__file__), Path(memory.__file__),
             ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',
             ROOT/'elliptic-curves/cas/research_runtime/finite_reduction.py',
             ROOT/'elliptic-curves/cas/elliptic_candidate_record.py',
             ROOT/manifest['snapshot'], ART/'new_high_rank_curve_index_v20.json']
    put('replay.json', {'status': 'PASS', 'public_curves_certified': len(prepared),
        'public_rank_lower_bounds': {str(i): len(p[1]) for i, p in prepared.items()},
        'certified_recovery_controls': len(seeds), 'genuine_historical_controls': 1,
        'synthetic_deletion_controls': 5, 'searches_run': 0,
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
        'outputs': {p.name: cert.hashed(p) for p in sorted(OUT.glob('*.json'))
                    if p.name != 'replay.json'},
        'claim_boundary': 'Public-point reproduction, not new discoveries. Exact lower bounds '
            'and27-to28 seed/oracle separations are certified. No exact rank, conductor '
            'verification, family recognition, recovery success or search improvement yet.'}, check)
    print('SIX RECOVERY CONTROLS READY; NO POINT SEARCH RUN')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    run(args.check)
