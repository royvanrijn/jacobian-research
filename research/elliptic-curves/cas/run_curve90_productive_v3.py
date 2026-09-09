#!/usr/bin/env python3
"""Bounded adaptive M27 discovery with productive anchors and two map policies.

New protocol: frozen historical anchor subset, V3 selection, integer CVP,
cached finite admission, immediate rebuild, and immutable per-chart receipts.
Existing runners and controls remain untouched. Run with Sage Python under
the shared supervisor. Replay uses original rational CVP and complete clouds.
"""
import argparse
import fcntl
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

from research_runtime.store import checkpoint
from v3_warm_engine import load, certified_state, _audit, _odd
from memory_rank_certificate import checked_rank
from future_point_admission import FinitePointAdmission
from visibility_productive import landscape
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
PREP = LOCAL/'curve90-v3-preparation-v1'
BANK = LOCAL/'curve90-productive-anchor-bank-v1'
MAPPERS = {'quartic_minimized': 'prepare_extended20_mw16_pari_batch.sage',
           'factor_free': 'factor_free_pari_mapping.sage'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def guard(protocol):
    for category in ('inputs', 'sources'):
        for p, h in protocol[category].items():
            require(sha(ROOT/p) == h, 'changed binding: '+p)
    require(sha(Path('/usr/bin/gp')) == protocol['gp_sha256'], 'GP executable changed')


def coordinate_key(mapping):
    m = tuple(map(F, mapping['matrix']))
    scale = next(x for x in m if x)
    return tuple(x/scale for x in m)


def freeze(folder):
    import sage.version
    from sage.all import pari
    import numpy
    old = load('productive_sources', CAS/'adaptive_visibility_cascade_v3.sage')
    bank = read(BANK/'anchor-bank.json')
    require(bank['status'] == 'COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'productive bank required')
    require(bank['protocol_sha256'] == sha(BANK/'protocol.json') and
            bank['gains_sha256'] == sha(BANK/'gains.json'), 'productive bank seal differs')
    seed = read(PREP/'seed-M27.json')
    prepared = read(PREP/'prepared.json')
    require(prepared['files']['seed-M27.json'] == sha(PREP/'seed-M27.json'), 'M27 seed changed')
    names = ('run_curve90_productive_v3.py', 'visibility_productive.py',
             'visibility_lattice_fast.py', 'future_point_admission.py',
             'memory_rank_certificate.py', 'v3_warm_engine.py', 'v3_warm_support.py',
             'prepare_extended20_mw16_pari_batch.sage')
    inputs = [PREP/'seed-M27.json', PREP/'prepared.json', BANK/'anchor-bank.json',
              BANK/'protocol.json', BANK/'gains.json']
    protocol = {'schema': 'curve90-productive-adaptive.v1', 'initial_rank': 27,
        'target_rank': 32, 'generic_rank': 16, 'scaled_shells': bank['shells'],
        'anchors_per_shell': 16, 'canonical_per_shell': 25,
        'exact_cvp_node_limit': 2000000, 'prime_bound': 1000,
        'height': 125000, 'seconds_per_chart': 10, 'max_charts': 4096,
        'max_epochs': 6, 'mapping_order': list(MAPPERS),
        'search_wall_limit_seconds': 7200, 'replay_wall_limit_seconds': 7200,
        'rss_limit_bytes': 3*1024**3, 'maximum_workers': 1,
        'software': {'sage': sage.version.version, 'pari': str(pari.version()), 'numpy': numpy.__version__},
        'gp_sha256': sha(Path('/usr/bin/gp')),
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        'sources': {**old.sources(), **{str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in names}},
        'selection': 'V3 diversified extension shortlist over exactly eight historical productive '
            'generic anchors. Full extension scoring, exact CVP. Factor-free chart profiles '
            'determine centre order. For each centre try quartic-minimized then factor-free '
            'coordinates, skipping identical projective coordinates. Rebuild immediately '
            'after the first standalone certified gain.',
        'scope': 'New directions beyond the certified M27 on curve90. Historical productive '
            'anchor selection is retrospective; no withheld/higher-rank point or visibility '
            'labels enter execution. The two masked pointwise nulls are not passed controls '
            'or rank upper bounds. This is a bounded discovery experiment, not a success prediction.'}
    folder.mkdir(exist_ok=False)
    checkpoint(folder/'seed.json', seed)
    checkpoint(folder/'bank.json', bank)
    protocol['inputs'].update({str((folder/n).relative_to(ROOT)): sha(folder/n)
                               for n in ('seed.json', 'bank.json')})
    checkpoint(folder/'protocol.json', protocol)
    return protocol


def normalized_selection(data):
    result = json.loads(json.dumps(data))
    for anchor in result['anchors']:
        anchor.pop('maps_sha256', None)
    return result


def search(folder):
    protocol = read(folder/'protocol.json') if folder.exists() else freeze(folder)
    guard(protocol)
    if (folder/'terminal.json').exists():
        print('SEALED_SEARCH_REUSED', flush=True)
        return
    seed, bank = read(folder/'seed.json'), read(folder/'bank.json')
    model = tuple(map(F, seed['curve']))
    points = tuple(tuple(map(F, p)) for p in seed['points'])
    proof = seed['proof']
    require(len(points) == 27, 'must start at M27')
    certified_state(model, points, proof)
    admission = FinitePointAdmission(model, points, prime_bound=1000)
    modules = {key: load('productive_'+key, CAS/name) for key, name in MAPPERS.items()}
    for module in modules.values():
        module.pari.allocatemem(256000000, silent=True)
    tested, total, stages = set(), 0, []
    stop = 'EPOCH_BUDGET_EXHAUSTED'
    for epoch in range(protocol['max_epochs']):
        wd = folder/f'epoch-{epoch:02d}'
        wd.mkdir(exist_ok=True)
        state = certified_state(model, points, proof)
        before = len(points)
        packet = {'curve': seed['curve'], 'points': [list(map(str, p)) for p in points], 'proof': proof}
        if (wd/'seed.json').exists():
            require(read(wd/'seed.json') == json.loads(json.dumps(packet)), 'resumed epoch seed differs')
        else:
            checkpoint(wd/'seed.json', packet)
        ld = wd/'landscape'
        if (ld/'selection.json').exists():
            selection = read(ld/'selection.json')
            require(selection['basis'] == packet['points'], 'resumed landscape basis differs')
        else:
            if ld.exists():
                ld.rename(wd/f'incomplete-landscape-{time.time_ns()}')
            selection = landscape(model, points, tested, ld, protocol, bank)
        previous = sha(ld/'selection.json')
        charts, skipped, grew = [], [], False
        for ci, centre in enumerate(selection['centres']):
            seen_coordinates = set()
            for key, module in modules.items():
                if total >= protocol['max_charts']:
                    break
                mapping = module.mapping(model, points, centre)
                coordinate = coordinate_key(mapping)
                if coordinate in seen_coordinates:
                    skipped.append({'centre_index': ci, 'policy': key, 'reason': 'SAME_PROJECTIVE_COORDINATE'})
                    continue
                seen_coordinates.add(coordinate)
                point_search = PointedQuarticSearch(state=state,
                    centre={'coefficients': centre['representative']},
                    coordinate_policy=mapping['coordinate_policy'])
                path = wd/f'chart-{len(charts):04d}.json'
                if path.exists():
                    chart = read(path)
                    require(chart['previous_sha256'] == previous and chart['centre_index'] == ci and
                            chart['policy'] == key and chart['mapping'] == json.loads(json.dumps(mapping)),
                            'resumed chart schedule differs')
                    returned = backend.replay(point_search, mapping, chart['search'])
                else:
                    transcript, returned = backend.execute(point_search, mapping, protocol['height'],
                        protocol['seconds_per_chart'], protocol['gp_sha256'])
                    require(backend.replay(point_search, mapping, transcript) == returned, 'point replay differs')
                    chart = {'index': len(charts), 'centre_index': ci, 'policy': key,
                        'previous_sha256': previous, 'mapping': mapping, 'search': transcript}
                    checkpoint(path, chart)
                charts.append(chart)
                previous = sha(path)
                total += 1
                for point in returned:
                    admission.consider(point)
                    if len(admission.points) > before:
                        points = tuple(admission.points)
                        proof = checked_rank(model, points, admission.primes,
                                             seed['proof']['no_rational_2_torsion_prime'])
                        require(proof['rank_lower_bound'] == before+1, 'standalone gain proof differs')
                        grew = True
                        checkpoint(wd/'gain.json', {'curve': seed['curve'],
                            'points': [list(map(str, p)) for p in points], 'proof': proof,
                            'chart_index': len(charts)-1, 'chart_sha256': previous})
                        break
                checkpoint(folder/'progress.json', {'phase': 'SEARCHING', 'epoch': epoch,
                    'charts': total, 'rank_lower_bound': len(points), 'last_chart_sha256': previous})
                print('PRODUCTIVE_SEARCH', epoch, total, 'rank', len(points), key, flush=True)
                if grew:
                    break
            tested.add(tuple(centre['point']))
            if grew or total >= protocol['max_charts']:
                break
        if not charts:
            stop = 'NO_FRESH_CENTRES'
            break
        # One complete cloud audit per epoch, never one cumulative copy per chart.
        snapshot = {'status': 'RETAINED_PRODUCTIVE_EPOCH_CLOUD', 'family': 'a1-fibration-01',
            'parameter': '-1867/270', 'curve': seed['curve'], 'charts': charts, 'final_state': state.record()}
        snap = wd/'cloud.json'
        if snap.exists():
            require(read(snap) == json.loads(json.dumps(snapshot)), 'resumed cloud differs')
        else:
            checkpoint(snap, snapshot)
        audit = _audit(snap, wd/'mod2.json')
        odd = _odd(wd/'mod2.json', wd/'modl.json')
        censored = sum(c['search']['status'] != 'bounded_search_complete' for c in charts)
        if audit['rank_lower_bound'] != len(points) or any(r != len(points) for r in odd.values()):
            stop = 'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION'
        elif len(points) >= protocol['target_rank']:
            stop = 'TARGET_LOWER_BOUND_REACHED'
        elif total >= protocol['max_charts']:
            stop = 'CHART_BUDGET_EXHAUSTED'
        elif grew:
            stop = 'REBUILD_AFTER_CERTIFIED_GAIN'
        else:
            stop = 'CENSORED_FINITE_EXPOSURE' if censored else 'FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN'
        stage = {'epoch': epoch, 'before': before, 'after': len(points), 'charts': len(charts),
            'skipped': skipped, 'last_chart_sha256': previous, 'stop_reason': stop,
            'selection_sha256': sha(ld/'selection.json'), 'mod2_sha256': sha(wd/'mod2.json'),
            'modl_sha256': sha(wd/'modl.json'), 'odd_ranks': odd, 'censored': censored}
        if (wd/'stage.json').exists():
            require(read(wd/'stage.json') == stage, 'resumed stage differs')
        else:
            checkpoint(wd/'stage.json', stage)
        stages.append(stage)
        guard(protocol)
        if stop != 'REBUILD_AFTER_CERTIFIED_GAIN':
            break
    checkpoint(folder/'terminal.json', {'status': 'SEARCH_FINISHED_PENDING_INDEPENDENT_REPLAY',
        'stop_reason': stop, 'charts': total, 'rank_lower_bound': len(points),
        'curve': seed['curve'], 'points': [list(map(str, p)) for p in points], 'proof': proof,
        'stages': stages, 'protocol_sha256': sha(folder/'protocol.json'),
        'claim_boundary': 'Subgroup lower bounds only. Unadmitted points and bounded misses remain UNKNOWN.'})


def replay(folder):
    protocol, terminal = read(folder/'protocol.json'), read(folder/'terminal.json')
    guard(protocol)
    require(terminal['protocol_sha256'] == sha(folder/'protocol.json'), 'terminal protocol differs')
    if (folder/'verified.json').exists():
        require(read(folder/'verified.json')['terminal_sha256'] == sha(folder/'terminal.json'), 'verified terminal changed')
        return
    seed, bank = read(folder/'seed.json'), read(folder/'bank.json')
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    expected_points = seed['points']
    all_returned = set(basis)
    tested, total = set(), 0
    modules = {k: load('replay_productive_'+k, CAS/n) for k, n in MAPPERS.items()}
    for module in modules.values():
        module.pari.allocatemem(256000000, silent=True)
    for stage in terminal['stages']:
        wd = folder/f"epoch-{stage['epoch']:02d}"
        packet = read(wd/'seed.json')
        require(packet['points'] == expected_points, 'adaptive prefix differs')
        basis = tuple(tuple(map(F, p)) for p in packet['points'])
        state = certified_state(model, basis, packet['proof'])
        rd = wd/'reference-landscape'
        if (rd/'selection.json').exists():
            regenerated = read(rd/'selection.json')
        else:
            if rd.exists():
                rd.rename(wd/f'incomplete-reference-{time.time_ns()}')
            regenerated = landscape(model, basis, tested, rd, protocol, bank, reference=True)
        saved = read(wd/'landscape/selection.json')
        require(normalized_selection(regenerated) == normalized_selection(saved), 'full reference selection differs')
        require(stage['selection_sha256'] == sha(wd/'landscape/selection.json'), 'selection hash differs')
        previous, index, skipped = sha(wd/'landscape/selection.json'), 0, []
        schedule = []
        for ci, centre in enumerate(saved['centres']):
            seen = set()
            for key, module in modules.items():
                mapping = module.mapping(model, basis, centre)
                coordinate = coordinate_key(mapping)
                if coordinate in seen:
                    skipped.append({'centre_index': ci, 'policy': key, 'reason': 'SAME_PROJECTIVE_COORDINATE'})
                    continue
                seen.add(coordinate)
                schedule.append((ci, centre, key, mapping))
        require(stage['charts'] <= len(schedule), 'more charts than declared schedule')
        if stage['stop_reason'] in ('FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN', 'CENSORED_FINITE_EXPOSURE'):
            require(stage['charts'] == len(schedule), 'incomplete exposure labelled exhausted')
        for ci, centre, key, mapping in schedule[:stage['charts']]:
                path = wd/f'chart-{index:04d}.json'
                chart = read(path)
                require(chart['index'] == index and chart['centre_index'] == ci and
                    chart['policy'] == key and chart['previous_sha256'] == previous and
                    chart['mapping'] == json.loads(json.dumps(mapping)), 'replayed chart order/map differs')
                search_object = PointedQuarticSearch(state=state,
                    centre={'coefficients': centre['representative']},
                    coordinate_policy=mapping['coordinate_policy'])
                require(chart['search']['height_bound'] == protocol['height'] and
                    chart['search']['timeout_seconds'] == protocol['seconds_per_chart'] and
                    chart['search']['gp_binary_sha256'] == protocol['gp_sha256'], 'chart budget differs')
                all_returned.update(backend.replay(search_object, mapping, chart['search']))
                previous, index, total = sha(path), index+1, total+1
                tested.add(tuple(centre['point']))
        require(previous == stage['last_chart_sha256'] and index == stage['charts'], 'chart chain differs')
        cloud = read(wd/'cloud.json')
        require(cloud['charts'] == [read(wd/f'chart-{j:04d}.json') for j in range(index)] and
                cloud['final_state'] == state.record() and cloud['curve'] == seed['curve'],
                'complete-cloud provenance differs')
        _audit(wd/'cloud.json', wd/'mod2.json')
        _odd(wd/'mod2.json', wd/'modl.json')
        require(sha(wd/'mod2.json') == stage['mod2_sha256'] and
                sha(wd/'modl.json') == stage['modl_sha256'], 'epoch proof differs')
        if stage['after'] > stage['before']:
            gain = read(wd/'gain.json')
            gp = tuple(tuple(map(F, p)) for p in gain['points'])
            require(gp[:len(basis)] == basis and all(p in all_returned for p in gp), 'gain point provenance differs')
            certified_state(model, gp, gain['proof'])
            expected_points = gain['points']
        checkpoint(wd/'reference-verified.json', {'status': 'PASS_FULL_PRODUCTIVE_EPOCH_REPLAY',
            'stage_sha256': sha(wd/'stage.json'), 'rank_lower_bound': stage['after']})
    require(total == terminal['charts'] and terminal['points'] == expected_points, 'terminal points/count differ')
    final_state = certified_state(model, tuple(tuple(map(F, p)) for p in terminal['points']), terminal['proof'])
    require(final_state.rank == terminal['rank_lower_bound'], 'terminal lower bound differs')
    guard(protocol)
    checkpoint(folder/'verified.json', {'status': 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY',
        'terminal_sha256': sha(folder/'terminal.json'), 'rank_lower_bound': terminal['rank_lower_bound'],
        'charts': total, 'stop_reason': terminal['stop_reason'],
        'claim_boundary': 'Exact subgroup lower bound and declared finite-policy replay, not an exact rank.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('search', 'replay'))
    parser.add_argument('--folder', required=True, type=Path)
    args = parser.parse_args()
    with (args.folder.parent/(args.folder.name+'.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        (search if args.mode == 'search' else replay)(args.folder)
