#!/usr/bin/env python3
"""Fixed-curve saturation audit and twelve-chart million-height experiment."""
import argparse
import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import audit_retained_cloud_modl as odd
import certify_compact_r17_candidates as cert
import certify_small_conductor_curve as original
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.memory_store import MemoryFactStore
from research_runtime.pointed_orbit_compression import compress
from research_runtime.preloaded_prime_state import preload
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.rotated_observation_state import rotate
from research_runtime.search_state import raw_state
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
D = LOCAL / 'small-conductor-targeted-study-v2'
PROOF = ART / 'small_conductor_rank22_proof_v1.json'
MAPS = LOCAL / 'small-conductor-adaptive-pari-maps-v1/maps.json'
ODD = ART / 'small_conductor_all_retained_modl_v1.json'
OUT = ART / 'small_conductor_targeted_study_v2.json'


def sources():
    paths = [Path(__file__).resolve(), Path(original.__file__),
             ROOT / 'elliptic-curves/cas/memory_rank_certificate.py']
    for module in ('search_state', 'cached_observation_state', 'memory_store',
                   'pointed_orbit_compression', 'preloaded_prime_state',
                   'quotient_only_reduction', 'rotated_observation_state', 'supervisor'):
        paths.append(ROOT / ('elliptic-curves/cas/research_runtime/' + module + '.py'))
    return {**backend.sources(), **odd.sources(),
            **{str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}}


def arithmetic():
    proof = cert.read(PROOF)
    original.verify(proof)
    model = tuple(map(cert.F, proof['short_model']))
    points = [tuple(map(cert.F, p)) for p in proof['short_points']]
    audits = []
    cloud = cert.read(ODD)
    if cloud['curve'] != proof['short_model'] or cloud['points'][:22] != proof['short_points']:
        raise ArithmeticError('same ordered subgroup required')
    for old in cloud['audits']:
        ell = old['modulus']
        signatures = [odd.finite.signature(model, points, r['prime'], ell)
                      for r in old['signatures']]
        rows = [row for sig in signatures for row in sig.rows]
        if odd.finite.pivots(rows, ell) != list(range(22)):
            raise ArithmeticError('saturation UNKNOWN')
        tp = old['no_rational_ell_torsion_prime']
        if not odd.ml.no_rational_l_torsion_reduction_certificate(model, tp, ell):
            raise ArithmeticError('torsion exclusion failed')
        audits.append({'modulus': ell, 'no_rational_ell_torsion_prime': tp,
                       'signatures': [asdict(sig) for sig in signatures]})
    if [a['modulus'] for a in audits] != [3, 5]:
        raise ArithmeticError('both odd moduli required')
    catalogue = cert.read(D / 'catalogue.json')
    comparison = []
    for rank in (22, 23, 24):
        rows = [r for r in catalogue['curves'] if r['rank_lower_bound'] >= rank]
        known = [r for r in rows if r.get('conductor')]
        best = min(known, key=lambda r: int(r['conductor']))
        comparison.append({'threshold': rank, 'minimum_id': best['id'],
                           'minimum_conductor': best['conductor'],
                           'minimum_over_target': str(cert.F(int(best['conductor']), int(proof['conductor']))),
                           'smaller_ids': sorted(r['id'] for r in known if int(r['conductor']) < int(proof['conductor'])),
                           'missing_ids': sorted(r['id'] for r in rows if not r.get('conductor'))})
    matches = [r['id'] for r in catalogue['curves']
               if cert.isomorphic(model, tuple(map(cert.F, r['ainvs'])))]
    primes = ['0'] + [p for p, e in proof['discriminant_factorization']]
    program = 'E=ellinit([' + ','.join(proof['integral_model']) + ']);print([' + ','.join('ellrootno(E,' + p + ')' for p in primes) + ']);quit\n'
    result = subprocess.run(['/usr/bin/gp', '-q'], input=program, text=True,
                            capture_output=True, timeout=15, check=True)
    if result.stderr:
        raise ArithmeticError(result.stderr)
    signs = json.loads(result.stdout)
    product = 1
    for sign in signs:
        if sign not in (-1, 1):
            raise ArithmeticError('invalid root number')
        product *= sign
    return {'saturated_at_primes': [2, 3, 5], 'odd_injectivity_certificates': audits,
            'mod2_certificate_path': str(PROOF.relative_to(ROOT)),
            'saturation_lemma': 'If ell Q=sum a_i P_i, independent finite quotient columns force ell to divide every a_i. Absence of rational ell-torsion gives Q=sum (a_i/ell)P_i. Thus this subgroup is ell-saturated, even without knowing the whole rank.',
            'catalogue_count': len(catalogue['curves']), 'catalogue_comparisons': comparison,
            'catalogue_q_isomorphism_matches': matches,
            'local_root_numbers': dict(zip(primes, signs)), 'global_root_number': product,
            'gp_program': program, 'gp_stdout': result.stdout,
            'boundary': 'Finite-prime saturation of the explicit22-point subgroup; no full saturation or upper rank bound. Root number controls analytic parity; algebraic parity is conjectural. Catalogue fields are external reported data.'}


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve frozen protocol')
    a = arithmetic()
    maps = cert.read(MAPS)
    if maps['status'] != 'COMPLETE_DECLARED_MAPS' or len(maps['rows']) != 301:
        raise ArithmeticError('complete original maps required')
    selected, words = [], set()
    for i, m in enumerate(maps['rows']):
        word = m['centre']['quotient_word']
        if word not in words:
            words.add(word)
            selected.append(i)
        if len(selected) == 12:
            break
    checkpoint(D / 'arithmetic.json', a)
    inputs = [PROOF, ODD, MAPS, D / 'catalogue.json', D / 'catalogue.metadata.json', D / 'arithmetic.json']
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.small-conductor-targeted-study.v1',
        'sources': sources(), 'inputs': {str(p.relative_to(ROOT)): cert.hashed(p) for p in inputs},
        'selected_indices': selected, 'height': 1000000, 'seconds_per_chart': 30,
        'worker_wall_seconds': 600, 'replay_wall_seconds': 300,
        'rss_bytes': 1610612736, 'maximum_workers': 1, 'target_rank': 23,
        'gp_sha256': cert.hashed(Path('/usr/bin/gp')),
        'selection': 'First twelve distinct nonzero six-bit quotient words in the existing ordered301-map roster, without using new point outcomes.',
        'gate': 'User-requested focused follow-up on3/17. Prior301 boxes at100000 are complete and still certify22. Reuse twelve exact maps at tenfold height to test a remaining coordinate-visibility limit; no new parameter or descent campaign.',
        'checkpoint_policy': 'Checkpoint each completed or timed-out chart and its exact admission history. Stop at provisional23 for independent replay. No automatic extra wave.',
        'boundary': 'At most twelve boxes on one curve. Timeout has no completed-box claim; a finite miss has no rank upper-bound implication.'})
    print('PREPARED', selected, 'SATURATED AT2,3,5; ROOT NUMBER', a['global_root_number'], flush=True)


def protocol():
    p = cert.read(D / 'protocol.json')
    if p['sources'] != sources() or any(cert.hashed(ROOT / name) != h for name, h in p['inputs'].items()):
        raise ArithmeticError('frozen input or source changed')
    return p


def initial():
    proof = cert.read(PROOF)
    cache = ReductionCache(MemoryFactStore())
    model = tuple(map(cert.F, proof['short_model']))
    points = [tuple(map(cert.F, q)) for q in proof['short_points']]
    state = raw_state(model, points, cache=cache, prime_bound=997)
    state = MWState.from_record(state.record(), cache=cache)
    state, bank = preload(state, cache, 997)
    if state.rank != 22 or tuple(tuple(map(cert.F, q)) for q in state.basis) != tuple(points):
        raise ArithmeticError('ordered certified22 seed required')
    return model, cache, state, bank


def experiment(replay=False):
    p = protocol()
    model, cache, state, bank = initial()
    maps = cert.read(MAPS)
    path = D / 'result.json'
    if not replay and path.exists():
        raise FileExistsError('preserve experiment')
    data = {'protocol_hash': digest(p), 'initial_state': state.record(), 'bank': bank,
            'status': 'RUNNING', 'charts': []}
    saved = cert.read(path) if replay else None
    if replay and any(saved[k] != data[k] for k in ('protocol_hash', 'initial_state', 'bank')):
        raise ArithmeticError('initial state differs')
    if not replay:
        checkpoint(path, data)
    indices = p['selected_indices'][:len(saved['charts'])] if replay else p['selected_indices']
    for offset, index in enumerate(indices):
        m = maps['rows'][index]
        state, archive = rotate(state)
        ap = D / 'states' / f'{index:03}.json'
        if replay:
            if cert.read(ap) != archive:
                raise ArithmeticError('archive differs')
        else:
            checkpoint(ap, archive)
        rep = m['centre']['representative'] + [0] * (state.rank - 22)
        search = PointedQuarticSearch(state=state, centre={'coefficients': rep}, coordinate_policy=m['coordinate_policy'])
        if replay:
            record = saved['charts'][offset]['search']
            if (record['height_bound'], record['timeout_seconds'], record['gp_binary_sha256']) != (p['height'], p['seconds_per_chart'], p['gp_sha256']):
                raise ArithmeticError('budget changed')
            points = backend.replay(search, m, record)
        else:
            record, points = backend.execute(search, m, p['height'], p['seconds_per_chart'], p['gp_sha256'])
        compression = compress(model, state.basis, rep, points)
        for j in compression['kept_indices']:
            state = state.adjoin(points[j], cache=cache)
        row = {'index': index, 'search': record, 'archive_sha256': cert.hashed(ap),
               'compression': compression, 'state_key': state.key,
               'observations': state.record()['state']['observations'], 'rank_lower_bound': state.rank}
        if replay and saved['charts'][offset] != row:
            raise ArithmeticError('exact admission replay differs')
        data['charts'].append(row)
        data.update(final_state=state.record(), rank_lower_bound=state.rank)
        if not replay:
            checkpoint(path, data)
        print('REPLAY' if replay else 'SEARCH', offset + 1, 'index', index, record['status'], 'rank >=', state.rank, flush=True)
        if state.rank >= p['target_rank']:
            break
    data['status'] = 'TARGET_REACHED' if state.rank >= p['target_rank'] else 'COMPLETE_DECLARED_ATTEMPTS'
    data['rank_certificate'] = checked_rank(model, state.basis, state.reductions.primes, state.no_two_torsion_prime)
    if replay:
        if json.loads(json.dumps(data)) != saved:
            raise ArithmeticError('terminal result differs')
    else:
        checkpoint(path, data)


def launch():
    p = protocol()
    if (D / 'ledger.json').exists():
        raise FileExistsError('preserve ledger')
    ledger = {'status': 'RUNNING', 'stages': []}
    checkpoint(D / 'ledger.json', ledger)
    for stage, seconds in [('worker', p['worker_wall_seconds']), ('replay', p['replay_wall_seconds'])]:
        result = run([sys.executable, str(Path(__file__).resolve()), stage],
                     limits=Limits(seconds, p['rss_bytes']), log_path=D / (stage + '.log'),
                     checkpoint_path=D / (stage + '.supervisor.json'), cwd=ROOT)
        ledger['stages'].append({'stage': stage, 'supervision': result})
        ledger['status'] = 'RUNNING' if result['outcome'] == 'completed' and result['returncode'] == 0 else 'FAILED_OR_CENSORED'
        checkpoint(D / 'ledger.json', ledger)
        if ledger['status'] != 'RUNNING':
            return
    ledger['status'] = 'PASS'
    checkpoint(D / 'ledger.json', ledger)


def summarize(check=False):
    p = protocol()
    a = arithmetic()
    if json.loads(json.dumps(a)) != cert.read(D / 'arithmetic.json'):
        raise ArithmeticError('arithmetic replay differs')
    ledger, result = cert.read(D / 'ledger.json'), cert.read(D / 'result.json')
    if ledger['status'] != 'PASS':
        raise ArithmeticError('completed exact replay required')
    data = {'schema': 'elliptic-curves.small-conductor-targeted-result.v1', 'status': 'PASS',
            'protocol': p, 'arithmetic': a, 'rank_lower_bound': result['rank_lower_bound'],
            'rank_certificate': result['rank_certificate'],
            'points': result['final_state']['state']['reductions']['points'],
            'chart_count': len(result['charts']),
            'complete_boxes': sum(r['search']['status'] == 'bounded_search_complete' for r in result['charts']),
            'chart_statuses': [r['search']['status'] for r in result['charts']],
            'point_occurrences': sum(len(r['search']['finite_curve_points']) for r in result['charts']),
            'bindings': {str(q.relative_to(ROOT)): cert.hashed(q) for q in [D / 'result.json', D / 'ledger.json']},
            'claim_boundary': 'A finite search outcome and certified subgroup saturation at2,3,5. No rank upper bound, whole-group basis, extra point unless independently certified, or absolute conductor record.'}
    if check:
        if json.loads(json.dumps(data)) != cert.read(OUT):
            raise ArithmeticError('summary differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve summary')
        checkpoint(OUT, data)
    print('STUDY PASS', data['complete_boxes'], '/', data['chart_count'], 'complete; rank >=', data['rank_lower_bound'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'launch', 'worker', 'replay', 'summarize', 'check'])
    stage = parser.parse_args().stage
    if stage in ('worker', 'replay'):
        experiment(stage == 'replay')
    elif stage in ('summarize', 'check'):
        summarize(stage == 'check')
    else:
        globals()[stage]()
