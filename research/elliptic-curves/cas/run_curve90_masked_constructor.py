#!/usr/bin/env python3
"""Bounded M26 constructor control in sealed selector order, with cached admission.

No withheld coordinates or evaluation chart indices are read. Stop at the first
independently certified extra direction. This is recovery, not a new rank result.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

from research_runtime.store import checkpoint
from v3_warm_engine import certified_state, load
from future_point_admission import FinitePointAdmission
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend
from run_future_visibility_audit import guard, sources

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
PREP = LOCAL/'curve90-v3-preparation-v1'
MASKED = LOCAL/'future-v3-curve90-masked-v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def run(folder):
    if folder.exists():
        raise FileExistsError('preserve existing constructor attempt')
    seal, prior = read(MASKED/'sealed.json'), read(MASKED/'protocol.json')
    guard(prior)
    if (seal['status'] != 'SEALED_MASKED_M26_SELECTION'
            or seal['protocol_sha256'] != sha(MASKED/'protocol.json')
            or seal['selection_sha256'] != sha(MASKED/'landscape/selection.json')):
        raise ArithmeticError('unsealed masked selection')
    seed, selection = read(PREP/'seed-M26.json'), read(MASKED/'landscape/selection.json')
    if selection['basis'] != seed['points']:
        raise ArithmeticError('masked seed changed')
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    state = certified_state(model, basis, seed['proof'])
    old = load('constructor_v3_sources', CAS/'adaptive_visibility_cascade_v3.sage')
    inputs = [PREP/'seed-M26.json', MASKED/'landscape/selection.json',
              MASKED/'sealed.json', MASKED/'protocol.json']
    names = ('run_curve90_masked_constructor.py', 'future_point_admission.py',
             'memory_rank_certificate.py')
    protocol = {'schema': 'curve90-masked-constructor.v1', 'initial_rank': 26,
        'target_rank': 27, 'height': 125000, 'seconds_per_chart': 10,
        'maximum_charts': len(selection['centres']), 'prime_bound': 1000,
        'gp_sha256': sha(Path('/usr/bin/gp')),
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        'sources': {**sources(old), **{str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in names}},
        'scope': 'Retrospective masked recovery in sealed selector order; no oracle point '
                 'or visibility label read by search. Stop at first certified M27.'}
    folder.mkdir()
    checkpoint(folder/'protocol.json', protocol)
    admission = FinitePointAdmission(model, basis, prime_bound=1000)
    mapper = load('constructor_factor_free_mapper', CAS/'factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    previous = sha(folder/'protocol.json')
    start = time.monotonic()
    status, final_proof = 'BOUNDED_NO_CERTIFIED_EXTRA_DIRECTION', seed['proof']
    total = 0
    for index, centre in enumerate(selection['centres']):
        mapping = mapper.mapping(model, basis, centre)
        search = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                     coordinate_policy=mapping['coordinate_policy'])
        transcript, points = backend.execute(search, mapping, 125000, 10, protocol['gp_sha256'])
        if backend.replay(search, mapping, transcript) != points:
            raise ArithmeticError('backend witness replay differs')
        observations = []
        for point in points:
            observation = admission.consider(point)
            observations.append(observation)
            if len(admission.points) > len(basis):
                final_proof = checked_rank(model, admission.points, admission.primes,
                                           seed['proof']['no_rational_2_torsion_prime'])
                if final_proof['rank_lower_bound'] != 27:
                    raise ArithmeticError('standalone rank27 proof required')
                status = 'CERTIFIED_MASKED_M27_RECOVERY'
                break
        chart = {'index': index, 'previous_sha256': previous, 'centre': centre,
                 'mapping': mapping, 'search': transcript, 'observations': observations}
        path = folder/f'chart-{index:04d}.json'
        checkpoint(path, chart)
        previous, total = sha(path), index+1
        checkpoint(folder/'progress.json', {'charts': total, 'status': status,
            'rank_lower_bound': final_proof['rank_lower_bound'], 'last_chart_sha256': previous})
        print('MASKED_CONSTRUCTOR', total, status, flush=True)
        if status == 'CERTIFIED_MASKED_M27_RECOVERY':
            break
    guard(protocol)
    checkpoint(folder/'terminal.json', {'status': status, 'charts': total,
        'wall_seconds': time.monotonic()-start, 'last_chart_sha256': previous,
        'protocol_sha256': sha(folder/'protocol.json'),
        'curve': seed['curve'], 'points': [list(map(str, p)) for p in admission.points],
        'rank_certificate': final_proof, 'rank_lower_bound': final_proof['rank_lower_bound'],
        'claim_boundary': 'Retrospective recovery lower bound only; bounded misses remain UNKNOWN.'})


def replay(folder):
    protocol, terminal = read(folder/'protocol.json'), read(folder/'terminal.json')
    guard(protocol)
    seed = read(PREP/'seed-M26.json')
    selection = read(MASKED/'landscape/selection.json')
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    state = certified_state(model, basis, seed['proof'])
    previous = sha(folder/'protocol.json')
    returned = set(basis)
    for index in range(terminal['charts']):
        path = folder/f'chart-{index:04d}.json'
        chart = read(path)
        if (chart['index'] != index or chart['previous_sha256'] != previous
                or chart['centre'] != selection['centres'][index]):
            raise ArithmeticError('chart order/chain differs')
        search = PointedQuarticSearch(state=state,
            centre={'coefficients': chart['centre']['representative']},
            coordinate_policy=chart['mapping']['coordinate_policy'])
        if (chart['search']['height_bound'] != protocol['height']
                or chart['search']['timeout_seconds'] != protocol['seconds_per_chart']
                or chart['search']['gp_binary_sha256'] != protocol['gp_sha256']):
            raise ArithmeticError('constructor budget or backend differs')
        returned.update(backend.replay(search, chart['mapping'], chart['search']))
        previous = sha(path)
    points = tuple(tuple(map(F, p)) for p in terminal['points'])
    if (points[:26] != basis or any(p not in returned for p in points)
            or previous != terminal['last_chart_sha256']
            or terminal['protocol_sha256'] != sha(folder/'protocol.json')):
        raise ArithmeticError('terminal basis provenance differs')
    old = terminal['rank_certificate']
    proof = checked_rank(model, points, [s['prime'] for s in old['signatures']],
                         old['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(proof)) != old:
        raise ArithmeticError('independent terminal rank proof differs')
    guard(protocol)
    output = folder/'verified.json'
    if output.exists():
        raise FileExistsError('preserve replay')
    checkpoint(output, {'status': 'PASS_MASKED_CONSTRUCTOR_POINT_AND_RANK_REPLAY',
        'terminal_sha256': sha(folder/'terminal.json'), 'charts': terminal['charts'],
        'rank_lower_bound': proof['rank_lower_bound'],
        'claim_boundary': 'Exact returned-point provenance and lower bound; no new-rank discovery.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('run', 'replay'))
    p.add_argument('--folder', type=Path, required=True)
    args = p.parse_args()
    (run if args.mode == 'run' else replay)(args.folder)
