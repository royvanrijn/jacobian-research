#!/usr/bin/env python3
"""One next-direction comparison on frozen pre-acquisition centre banks.

No lattice enumeration, parameter sweep, descent, or automatic bank expansion.
Use prepare_next_direction_benchmark.py for allowlisted input projection.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import sys
import time

from research_runtime.store import checkpoint, digest
from memory_rank_certificate import checked_rank
from v3_warm_engine import certified_state
from future_point_admission import FinitePointAdmission
from pointed_quartic_search import PointedQuarticSearch
from pointed_box_equivalence import box_key
from lean_preconditioned_map_receipts import obtain
from search_observability import multiply
import pari_pointed_backend as backend

IDENTITY = (1, 0, 0, 1)
SHEARS = (IDENTITY, (1, 1, 0, 1), (1, -1, 0, 1),
          (1, 0, 1, 1), (1, 0, -1, 1))
POLICIES = {
    'v3_dual_125k': {'maps': ['preconditioned_full', 'factor_free'], 'height': 125000, 'changes': [IDENTITY]},
    'factor_free_125k': {'maps': ['factor_free'], 'height': 125000, 'changes': [IDENTITY]},
    'v3_dual_500k': {'maps': ['preconditioned_full', 'factor_free'], 'height': 500000, 'changes': [IDENTITY]},
    'factor_free_shears_125k': {'maps': ['factor_free'], 'height': 125000, 'changes': SHEARS},
}


def read(p):
    return json.loads(Path(p).read_text())


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def homogeneous(coefficients, matrix):
    """Exact substitution for a binary form, retaining its declared degree."""
    n = len(coefficients)-1
    a, b, c, d = map(F, matrix)
    if a*d == b*c:
        raise ArithmeticError('singular coordinate change')
    out = [F(0)]*(n+1)
    for i, value in enumerate(map(F, coefficients)):
        for j in range(i+1):
            for k in range(n-i+1):
                out[j+k] += value*comb(i,j)*a**j*b**(i-j)*comb(n-i,k)*c**k*d**(n-i-k)
    return list(map(str, out))


def change_coordinates(mapping, matrix):
    result = deepcopy(mapping)
    result['second_matrix'] = list(map(str, multiply(tuple(map(F,mapping['second_matrix'])), matrix)))
    result['matrix'] = list(map(str, multiply(tuple(map(F,mapping['matrix'])), matrix)))
    for key in ('reduced_P', 'reduced_Q', 'discriminant_quartic'):
        result[key] = homogeneous(mapping[key], matrix)
    result['coordinate_policy'] = {'kind':'raw', 'matrix': result['matrix']}
    result['benchmark_change'] = list(matrix)
    return result


def guard_artifacts(folder):
    """Application read restriction, not an operating-system sandbox."""
    folder = folder.resolve()
    def hook(event, args):
        if event != 'open' or not isinstance(args[0], (str, bytes)):
            return
        p = Path(args[0]).resolve()
        code = Path(__file__).resolve().parent
        if 'artifacts' in p.parts and not (p.is_relative_to(folder) or p.is_relative_to(code)):
            raise PermissionError('benchmark worker cannot read sibling/historical artifacts: '+str(p))
    sys.addaudithook(hook)


def worker(folder, policy_name, seconds):
    started = time.monotonic()
    guard_artifacts(folder)
    seed = read(folder/'input.json')
    protocol = read(folder/'protocol.json')
    if seed['sha256'] != digest(seed['payload']):
        raise ArithmeticError('input projection changed')
    seed = seed['payload']
    if protocol['policies'] != json.loads(json.dumps(POLICIES)):
        raise ArithmeticError('frozen policies changed')
    model = tuple(map(F, seed['curve']))
    points = tuple(tuple(map(F,p)) for p in seed['points'])
    proof = checked_rank(model, points, seed['primes'], seed['torsion_prime'])
    checkpoint(folder/'initial-rank.json', proof)
    state = certified_state(model, points, proof)
    admission = FinitePointAdmission(model, points, prime_bound=1000)
    lane = POLICIES[policy_name]
    seen, rows = set(), []
    excluded = {tuple(row) for row in seed.get('completed_boxes', [])}
    status = 'NO_NEW_DIRECTION_IN_COMPLETED_PREFIX'
    for ci, centre in enumerate(seed['centres']):
        if time.monotonic()-started >= seconds-20:
            status = 'CENSORED_ARM_ALLOWANCE'
            break
        for key in lane['maps']:
            mapping, receipt, limited = obtain(folder, ci, key, model, points, centre, state, protocol)
            if limited:
                rows.append({'centre':ci, 'map':key, 'status':limited, 'receipt':receipt})
                continue
            for wi, W in enumerate(lane['changes']):
                candidate = change_coordinates(mapping, W)
                signature = (str(centre['representative']), *map(str, box_key(candidate['matrix'])), str(lane['height']))
                if signature in seen or signature in excluded:
                    continue
                seen.add(signature)
                if time.monotonic()-started >= seconds-20:
                    status = 'CENSORED_ARM_ALLOWANCE'
                    break
                search = PointedQuarticSearch(state=state, centre={'coefficients':centre['representative']},
                    coordinate_policy=candidate['coordinate_policy'])
                transcript, returned = backend.execute(search, candidate, lane['height'], 10, protocol['gp_sha256'])
                if backend.replay(search, candidate, transcript) != returned:
                    raise ArithmeticError('exact point transport replay changed')
                name = f'chart-{len(rows):04d}.json'
                checkpoint(folder/name, {'centre_index':ci, 'mapping':candidate, 'search':transcript})
                row = {'centre':ci, 'map':key, 'change':wi, 'chart':name, 'sha256':sha(folder/name), 'status':transcript['status']}
                rows.append(row)
                for p in returned:
                    if admission.consider(p)['status'] == 'INDEPENDENT_FINITE_COLUMN':
                        proof = checked_rank(model, admission.points, admission.primes, seed['torsion_prime'])
                        checkpoint(folder/'gain.json', {'curve':seed['curve'], 'points':[list(map(str,p)) for p in admission.points],
                            'proof':proof, 'chart':name, 'chart_sha256':row['sha256']})
                        status = 'NEW_DIRECTION_PENDING_INDEPENDENT_REPLAY'
                        break
                checkpoint(folder/'progress.json', {'status':status, 'centres_seen':ci+1,
                    'point_calls':sum('chart' in r for r in rows), 'wall_seconds':time.monotonic()-started})
                if (folder/'gain.json').exists():
                    break
            if (folder/'gain.json').exists() or status == 'CENSORED_ARM_ALLOWANCE':
                break
        if (folder/'gain.json').exists() or status == 'CENSORED_ARM_ALLOWANCE':
            break
    if status == 'NO_NEW_DIRECTION_IN_COMPLETED_PREFIX' and any(r['status'] != 'bounded_search_complete' for r in rows):
        status = 'NO_GAIN_WITH_INCOMPLETE_MAPS_OR_BOXES'
    checkpoint(folder/'result.json', {'status':status, 'input_sha256':sha(folder/'input.json'),
        'policy':policy_name, 'initial_rank':len(points), 'rows':rows,
        'worker_wall_seconds':time.monotonic()-started, 'independent_gain_verified':False,
        'boundary':'A fixed retained centre bank; historical landscape costs are not zero. No rank upper bound.'})


def replay(folder):
    guard_artifacts(folder)
    seed = read(folder/'input.json')['payload']
    model = tuple(map(F,seed['curve']))
    points = tuple(tuple(map(F,p)) for p in seed['points'])
    proof = checked_rank(model, points, seed['primes'], seed['torsion_prime'])
    # Digest canonicalization avoids treating tuple/list serialization as a failed certificate.
    if digest(proof) != digest(read(folder/'initial-rank.json')):
        raise ArithmeticError('initial finite certificate differs')
    state = certified_state(model, points, proof)
    result = read(folder/'result.json')
    returned_by_chart = {}
    for row in result['rows']:
        if 'chart' not in row:
            continue
        if sha(folder/row['chart']) != row['sha256']:
            raise ArithmeticError('chart hash changed')
        data = read(folder/row['chart'])
        centre = seed['centres'][row['centre']]
        search = PointedQuarticSearch(state=state, centre={'coefficients':centre['representative']},
            coordinate_policy=data['mapping']['coordinate_policy'])
        returned_by_chart[row['chart']] = backend.replay(search, data['mapping'], data['search'])
    gained = (folder/'gain.json').exists()
    if gained:
        gain = read(folder/'gain.json')
        new_points = tuple(tuple(map(F,p)) for p in gain['points'])
        if new_points[:-1] != points or new_points[-1] not in returned_by_chart[gain['chart']]:
            raise ArithmeticError('gain not produced by retained point map')
        if gain['chart_sha256'] != sha(folder/gain['chart']):
            raise ArithmeticError('gain chart changed')
        actual = checked_rank(model,new_points,[s['prime'] for s in gain['proof']['signatures']],seed['torsion_prime'])
        if digest(actual) != digest(gain['proof']):
            raise ArithmeticError('independent rank replay differs')
    checkpoint(folder/'verified.json', {'status':'PASS_POINT_AND_RANK_REPLAY',
        'independent_gain_verified':gained, 'rank_lower_bound':len(points)+int(gained),
        'result_sha256':sha(folder/'result.json'),
        'boundary':'Fresh-process exact map/square and finite-group replay; shared arithmetic implementation, not a new independent CAS implementation.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['worker','replay'])
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--policy', choices=POLICIES)
    parser.add_argument('--seconds', type=float, default=300)
    args = parser.parse_args()
    if args.mode == 'worker':
        if not args.policy:
            parser.error('--policy is required for a worker')
        worker(args.folder.resolve(),args.policy,args.seconds)
    else:
        replay(args.folder.resolve())
