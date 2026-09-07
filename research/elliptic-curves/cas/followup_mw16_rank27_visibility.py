#!/usr/bin/env python3
"""Twelve frozen own27 charts at height one million; no adaptive refill."""
import argparse
import json
import sys
from pathlib import Path

import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
import audit_recorded_point_mod2_rank_v3 as mod2
import audit_retained_cloud_modl as odd
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import run, Limits
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState
from research_runtime.search_state import raw_state
from study_mw16_rank27_visibility import ROOT, ART, LOCAL, OUTPUT as VISIBILITY

D = LOCAL / 'mw16-rank27-visibility-followup-v1'
OLD = LOCAL / 'new27-specialized-parity-six-v1/new-20260906-90'
PROTOCOL = D / 'protocol.json'
RESULT = D / 'result.json'
MOD2 = ART / 'mw16_rank27_visibility_followup_mod2_v1.json'
ODD = ART / 'mw16_rank27_visibility_followup_modl_v1.json'
REPORT = ART / 'mw16_rank27_visibility_followup_v1.json'


def sources():
    return {**backend.sources(), **mod2.sources(), **odd.sources(),
            str(Path(__file__).resolve().relative_to(ROOT)): cert.hashed(Path(__file__).resolve())}


def prepare():
    if PROTOCOL.exists():
        raise FileExistsError('preserve frozen follow-up protocol')
    seed = cert.read(OLD / 'seed.json')
    maps = cert.read(OLD / 'maps.json')
    visibility = cert.read(VISIBILITY)
    if seed['parameter'] != '-1867/270' or len(seed['points']) != 27 or visibility['status'] != 'PASS':
        raise ArithmeticError('certified fixed curve and visibility gate required')
    if maps['status'] != 'COMPLETE_DECLARED_MAPS' or len(maps['rows']) != 49:
        raise ArithmeticError('original frozen specialized map roster required')
    candidates = sorted(maps['sample'], key=lambda c: (-c['metric_norm'], c['parity']))[:49]
    if maps['centres'] != candidates or [r['centre'] for r in maps['rows']] != candidates:
        raise ArithmeticError('original specialized ordering differs')
    old = seed['rank_certificate']
    checked_rank(tuple(map(cert.F, seed['curve'])), [tuple(map(cert.F, p)) for p in seed['points']],
                 [s['prime'] for s in old['signatures']], old['no_rational_2_torsion_prime'])
    checkpoint(D / 'seed.json', seed)
    checkpoint(D / 'maps.json', {'rows': maps['rows'][:12]})
    checkpoint(PROTOCOL, {
        'schema': 'elliptic-curves.mw16-rank27-visibility-followup-protocol.v1',
        'sources': sources(),
        'inputs': {str(p.relative_to(ROOT)): cert.hashed(p) for p in
                   (D / 'seed.json', D / 'maps.json', VISIBILITY)},
        'historical_sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in
                               (OLD / 'seed.json', OLD / 'maps.json')},
        'gp_sha256': cert.hashed(Path('/usr/bin/gp')),
        'charts': 12, 'height': 1000000, 'seconds_per_chart': 15,
        'worker_wall_seconds': 300, 'rss_bytes': 1610612736, 'maximum_workers': 1,
        'gate': 'The recovered 27th point has exact best signed height 79466 in its '
                'adaptive charts versus 1379753695791601 in the original generic charts. '
                'A separate 49-chart own27 policy already completed at 125000 without gain. '
                'Probe coordinate exposure above that height on twelve of those existing maps. '
                'This known-point comparison motivates a bounded test, not a sensitivity guarantee.',
        'selection': 'First twelve maps in the previously frozen 49-chart specialized '
                     'numerical-norm order on this curve own27 subgroup. No new geometry, '
                     'oracle visibility scoring, point-based reranking, extra parameters or refill.',
        'policy': 'Fixed initial 27-point subgroup for all twelve charts. Finish the declared '
                  'twelve boxes, checkpoint every completed/censored result; exact whole-cloud '
                  'finite-reduction audits at 2,3,5 determine any improved lower bound afterwards.',
        'claim_boundary': 'At most 180 seconds of point invocations and 300 seconds worker wall '
                          'time. Timeouts are censored with no prefix-completeness claim. '
                          'No automatic escalation, whole-curve upper bound or point absence.'})


def protocol():
    p = cert.read(PROTOCOL)
    if p['sources'] != sources() or any(cert.hashed(ROOT / n) != h for n, h in p['inputs'].items()):
        raise ArithmeticError('frozen source/input binding changed')
    return p


def initial():
    seed = cert.read(D / 'seed.json')
    cache = QuotientOnlyReductionCache(MemoryFactStore())
    model = tuple(map(cert.F, seed['curve']))
    points = [tuple(map(cert.F, q)) for q in seed['points']]
    raw = raw_state(model, points, cache=cache, prime_bound=997)
    state = CachedObservationMWState.from_record(raw.record(), cache=cache)
    if state.rank != 27 or list(map(list, state.basis)) != list(map(list, points)):
        raise ArithmeticError('exact own27 subgroup required')
    return seed, state


def worker():
    p = protocol()
    if RESULT.exists():
        raise FileExistsError('preserve point trial')
    seed, state = initial()
    data = {k: seed[k] for k in ('curve', 'family', 'parameter', 'generic_points')}
    data.update(protocol_hash=digest(p), status='RUNNING', charts=[],
                final_state=state.record(), rank_lower_bound=27)
    checkpoint(RESULT, data)
    for i, mapping in enumerate(cert.read(D / 'maps.json')['rows']):
        search = PointedQuarticSearch(state=state,
            centre={'coefficients': mapping['centre']['representative']},
            coordinate_policy=mapping['coordinate_policy'])
        record, _ = backend.execute(search, mapping, p['height'], p['seconds_per_chart'], p['gp_sha256'])
        data['charts'].append({'index': i, 'search': record})
        checkpoint(RESULT, data)
        print('MW16 VISIBILITY FOLLOWUP', i+1, record['status'],
              len(record['finite_curve_points']), 'points', flush=True)
    data['status'] = 'COMPLETE_DECLARED_ATTEMPT_PENDING_CLOUD_AUDIT'
    checkpoint(RESULT, data)


def replay():
    p = protocol()
    seed, state = initial()
    data = cert.read(RESULT)
    if data['protocol_hash'] != digest(p) or data['final_state'] != state.record():
        raise ArithmeticError('fixed subgroup/protocol changed')
    if any(data[k] != seed[k] for k in ('curve', 'family', 'parameter', 'generic_points')):
        raise ArithmeticError('trial curve changed')
    if len(data['charts']) != p['charts'] or data['status'] != 'COMPLETE_DECLARED_ATTEMPT_PENDING_CLOUD_AUDIT':
        raise ArithmeticError('declared attempt incomplete')
    for i, (row, mapping) in enumerate(zip(data['charts'], cert.read(D / 'maps.json')['rows'], strict=True)):
        r = row['search']
        if row['index'] != i or r['height_bound'] != p['height'] or r['timeout_seconds'] != p['seconds_per_chart'] or r['gp_binary_sha256'] != p['gp_sha256']:
            raise ArithmeticError('declared chart or budget changed')
        search = PointedQuarticSearch(state=state,
            centre={'coefficients': mapping['centre']['representative']},
            coordinate_policy=mapping['coordinate_policy'])
        backend.replay(search, mapping, r)
    print('REPLAYED TWELVE MW16 VISIBILITY FOLLOWUP CHARTS', flush=True)


def report():
    p = protocol()
    data, cloud, modl = [cert.read(path) for path in (RESULT, MOD2, ODD)]
    if cloud['input_sha256'] != cert.hashed(RESULT) or modl['input_sha256'] != cert.hashed(MOD2):
        raise ArithmeticError('cloud bindings differ')
    points = list(cert.read(D / 'seed.json')['points'])
    seen = {(cert.F(x), abs(cert.F(y))) for x, y in points}
    for row in data['charts']:
        for q in row['search']['finite_curve_points']:
            key = cert.F(q['x']), abs(cert.F(q['y']))
            if key not in seen:
                seen.add(key)
                points.append([q['x'], q['y']])
    if points != cloud['points'] or modl['points'] != points:
        raise ArithmeticError('complete retained cloud provenance differs')
    return {'schema': 'elliptic-curves.mw16-rank27-visibility-followup.v1', 'status': 'PASS',
            'bindings': {str(path.relative_to(ROOT)): cert.hashed(path) for path in
                         (PROTOCOL, RESULT, MOD2, ODD)},
            'attempted_charts': len(data['charts']), 'height': p['height'],
            'seconds_per_chart': p['seconds_per_chart'],
            'completed_boxes': sum(r['search']['status'] == 'bounded_search_complete' for r in data['charts']),
            'timeouts': sum(r['search']['status'] == 'bounded_search_timeout' for r in data['charts']),
            'backend_failures': sum(r['search']['status'] == 'backend_failure' for r in data['charts']),
            'retained_points_up_to_sign': len(points), 'mod2_lower_bound': cloud['rank_lower_bound'],
            'odd_prime_lower_bounds': {str(a['modulus']): a['finite_column_rank'] for a in modl['audits']},
            'claim_boundary': 'This fixed trial cloud includes the own27 seed and all retained '
            'output, not every historical point. No upper rank bound or unrecorded-point claim. '
            'Completed coverage trusts the pinned PARI execution; censored boxes remain unknown.'}


def audit():
    replay()
    mod2.build(RESULT, MOD2, 997, cert.hashed(RESULT))
    odd.build(MOD2, ODD)
    mod2.check(MOD2)
    odd.check(ODD)
    if REPORT.exists():
        raise FileExistsError('preserve report')
    checkpoint(REPORT, report())


def check():
    replay()
    mod2.check(MOD2)
    odd.check(ODD)
    if report() != cert.read(REPORT):
        raise ArithmeticError('follow-up report differs')
    print(json.dumps(cert.read(REPORT), sort_keys=True), flush=True)


def launch():
    p = protocol()
    for stage, seconds in [('worker', p['worker_wall_seconds']), ('audit', 180)]:
        supervisor = D / (stage + '.supervisor.json')
        if supervisor.exists():
            raise FileExistsError('preserve stage supervision')
        result = run([sys.executable, str(Path(__file__).resolve()), stage],
                     limits=Limits(seconds, p['rss_bytes']), log_path=D / (stage+'.log'),
                     checkpoint_path=supervisor, cwd=ROOT)
        print(stage, result['outcome'], result['returncode'], flush=True)
        if result['outcome'] != 'completed' or result['returncode'] != 0:
            raise ArithmeticError('stage failed or censored; retain checkpoint')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'worker', 'replay', 'launch', 'audit', 'check'])
    globals()[parser.parse_args().stage]()
