#!/usr/bin/env python3
"""Generic-only V3 comparison for the frozen two-cover recovery experiment.

Uses the maintained MW16 upper-shell selector, bounded map workers and PARI
point search unchanged. No previous point-search output is an input.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
sys.path.insert(0, str(CAS))
EXPERIMENT = ROOT/'artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1'
DEFAULT = ROOT/'artifacts/local/elliptic-curves/constructed-class-blind-v3-v2'
READS = set()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    from research_runtime.store import checkpoint
    checkpoint(path, data)


def data_guard(folder, *, atlas=False):
    allowed = {EXPERIMENT/'generic.json'}
    if atlas:
        allowed.add(ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json')

    def audit(event, args):
        if event != 'open' or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        path = Path(os.fsdecode(args[0])).resolve()
        if path.is_relative_to(ROOT/'artifacts'):
            if path not in allowed and not path.is_relative_to(folder):
                raise PermissionError('comparison data boundary: '+str(path))
            if not path.is_relative_to(folder):
                READS.add(str(path.relative_to(ROOT)))
    sys.addaudithook(audit)


def sources():
    # Freeze the shared implementation and this adapter before execution.
    names = ['visibility_future_upper.py', 'visibility_selection_v3.py',
        'visibility_generic_bank.py', 'visibility_lattice_fast.py', 'visibility_lattice_v2.py',
        'adaptive_visibility_cascade_v3.sage', 'adaptive_visibility_cascade.sage',
        'prospective_half_lattice_v3.sage', 'factor_free_pari_mapping.sage',
        'lean_preconditioned_map_receipts.py', 'lean_preconditioned_map_worker.py',
        'lean_preconditioned_full_pari_mapping.sage', 'lean_factor_free_pari_mapping.sage',
        'pari_pointed_backend.py', 'pointed_quartic_search.py', 'search_observability.py',
        'memory_rank_certificate.py', 'future_point_admission.py', 'v3_warm_engine.py',
        'v3_warm_support.py', 'compact_mw16_specialization.py',
        'research_runtime/supervisor.py']
    paths = [Path(__file__).resolve()] + [CAS/n for n in names]
    return {str(p.relative_to(ROOT)): sha(p) for p in paths}


def freeze(folder):
    folder.mkdir(parents=True, exist_ok=False)
    p = {'schema': 'constructed-class-generic-only-v3.v1',
        'family': 'a1-fibration-05', 'parameter': '3/17', 'generic_rank': 16,
        'scaled_shells': [16, 19, 20, 23], 'anchors_per_shell': 16,
        'canonical_per_shell': 25, 'exact_cvp_node_limit': 2000000,
        'generic_bank_node_limit': 20000000, 'generic_ellipsoid_bound': 23,
        'prime_bound': 1000, 'height': 125000, 'seconds_per_chart': 10,
        'max_charts': 100, 'seconds_per_map': 5, 'map_rss_bytes': 1024**3,
        'map_python': shutil.which('sage'), 'gp_sha256': sha('/usr/bin/gp'),
        'map_order': ['preconditioned_full', 'factor_free'],
        'max_epochs': 17, 'preparation_seconds': 60, 'landscape_seconds': 60,
        'total_computation_seconds': 1199.098360198026, 'rss_bytes': 2*1024**3,
        'prior_failed_attempt_seconds': 0.9016398019739427,
        'prior_failure': 'v1 rejected an attempted shared-cache open before point search; use explicit in-memory finite quotients. Original adapter and failure log retained in sibling constructed-class-blind-v3-v1.',
        'selection': 'Maintained visibility_future_upper V3: all MW16 upper shells, standard shell quotas, exact CVP, full extension shortlist, actual-centre deduplication. Immediate rebuild after a certified gain; stop on first no-gain epoch or explicit resource bound.',
        'inputs': {str((EXPERIMENT/'generic.json').relative_to(ROOT)): sha(EXPERIMENT/'generic.json'),
                   'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json':
                   sha(ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json')},
        'sources': sources(),
        'boundary': 'Generic-only execution on a user-selected historical control. Additional independent directions and fixed-cover recovery are separate endpoints. No rank upper bound or prospective fibre selection.'}
    write(folder/'protocol.json', p)
    return p


def guard(p, *, atlas=False):
    for name, h in p['sources'].items():
        assert sha(ROOT/name) == h, name
    for name, h in p['inputs'].items():
        if atlas or name.endswith('/generic.json'):
            assert sha(ROOT/name) == h, name
    assert sha('/usr/bin/gp') == p['gp_sha256']


def prepare(folder):
    from sage.all import ZZ, matrix, pari
    import sage.version
    from visibility_generic_bank import enumerate_bank
    from memory_rank_certificate import checked_rank
    from mod2_reduction_independence import find_two_torsion_certificate_prime
    from future_point_admission import FinitePointAdmission
    import compact_mw16_specialization as spec
    started = time.monotonic()
    p = read(folder/'protocol.json')
    data_guard(folder, atlas=True)
    guard(p, atlas=True)
    seed = read(EXPERIMENT/'generic.json')
    model = tuple(map(F, seed['curve']))
    points = tuple(tuple(map(F, pt)) for pt in seed['points'])
    family = next(f for f in read(spec.ATLAS)['families'] if f['fibration_id'] == p['family'])
    native_model, native_points = spec.specialize(family, p['parameter'])
    assert native_model[3] == 16*model[3] and native_model[4] == 64*model[4]
    signs = []
    for old, new in zip(native_points, points):
        assert old[0] == 4*new[0] and abs(old[1]) == 8*abs(new[1])
        signs.append(1 if old[1] == 8*new[1] else -1)
    G = matrix(ZZ, [[int(2*F(x)*signs[i]*signs[j]) for j, x in enumerate(row)]
                    for i, row in enumerate(family['generic_height_gram'])])
    assert G.is_positive_definite()
    U = matrix(ZZ, pari(G).qflllgram()).transpose()
    assert abs(U.det()) == 1
    write(folder/'generic-metric.json', {'gram': [list(map(int, r)) for r in G.rows()],
          'gram_scale': 2, 'signs_from_atlas': signs, 'native_scale': '2',
          'LLL': [list(map(int, r)) for r in U.rows()], 'determinant': str(G.det())})
    admission = FinitePointAdmission(model, points, prime_bound=p['prime_bound'])
    tp = find_two_torsion_certificate_prime(model)
    proof = checked_rank(model, points, admission.primes, tp)
    write(folder/'seed.json', {**seed, 'rank_lower_bound': 16, 'proof': proof})
    tick = time.monotonic()
    bank = enumerate_bank([list(map(int, r)) for r in (U*G*U.transpose()).rows()],
                          [list(map(int, r)) for r in U.rows()],
                          bound=p['generic_ellipsoid_bound'], shells=p['scaled_shells'],
                          node_limit=p['generic_bank_node_limit'])
    for row in bank['rows']:
        w = row['word']
        assert sum(w[i]*int(G[i,j])*w[j] for i in range(16) for j in range(16)) == row['norm']
    bank['generic_gram_scale'] = 2
    write(folder/'bank.json', bank)
    guard(p, atlas=True)
    write(folder/'prepared.json', {'status': 'PASS_GENERIC_ONLY_PREPARATION',
        'wall_seconds': time.monotonic()-started, 'bank_seconds': time.monotonic()-tick,
        'bank_rows': len(bank['rows']), 'bank_nodes': bank['nodes'],
        'sage': sage.version.version, 'pari': str(pari.version()),
        'input_reads': sorted(READS), 'protocol_sha256': sha(folder/'protocol.json'),
        'seed_sha256': sha(folder/'seed.json'), 'bank_sha256': sha(folder/'bank.json')})
    print('PREPARED', len(bank['rows']), time.monotonic()-started, flush=True)


def search(folder):
    from cysignals.alarm import alarm, cancel_alarm
    from cysignals.signals import AlarmInterrupt
    from v3_warm_engine import certified_state
    from visibility_future_upper import landscape
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from lean_preconditioned_map_receipts import obtain
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    started = time.monotonic()
    p = read(folder/'protocol.json')
    data_guard(folder)
    guard(p)
    prepared = read(folder/'prepared.json')
    assert prepared['seed_sha256'] == sha(folder/'seed.json')
    assert prepared['bank_sha256'] == sha(folder/'bank.json')
    seed, bank = read(folder/'seed.json'), read(folder/'bank.json')
    model = tuple(map(F, seed['curve']))
    points = tuple(tuple(map(F, pt)) for pt in seed['points'])
    proof = seed['proof']
    admission = FinitePointAdmission(model, points, prime_bound=p['prime_bound'])
    tested, stages, total, gain_times = set(), [], 0, []
    stop = 'EPOCH_BOUND'
    for epoch in range(p['max_epochs']):
        wd = folder/f'epoch-{epoch:02d}'
        wd.mkdir()
        state = certified_state(model, points, proof)
        write(wd/'seed.json', {'curve': seed['curve'], 'points': [list(map(str, pt)) for pt in points], 'proof': proof})
        tick = time.monotonic()
        alarm(p['landscape_seconds'])
        try:
            selected = landscape(model, points, tested, wd/'landscape', p, bank)
        except AlarmInterrupt:
            stop = 'CENSORED_LANDSCAPE_TIMEOUT'
            write(wd/'incomplete.json', {'status': stop, 'wall_seconds': time.monotonic()-tick})
            break
        finally:
            cancel_alarm()
        landscape_seconds = time.monotonic()-tick
        charts, attempts, grew = [], [], False
        for ci, centre in enumerate(selected['centres']):
            coordinates = set()
            for lane in p['map_order']:
                if total >= p['max_charts']:
                    break
                tick = time.monotonic()
                mapping, receipt, limit = obtain(wd, ci, lane, model, points, centre, state, p)
                attempt = {'centre_index': ci, 'policy': lane, 'receipt_sha256': receipt,
                           'wall_seconds': time.monotonic()-tick, 'limit': limit}
                attempts.append(attempt)
                if limit:
                    continue
                m = tuple(map(F, mapping['matrix']))
                scale = next(x for x in m if x)
                key = tuple(x/scale for x in m)
                if key in coordinates:
                    attempt['skip'] = 'SAME_PROJECTIVE_COORDINATE'
                    continue
                coordinates.add(key)
                ps = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                         coordinate_policy=mapping['coordinate_policy'])
                transcript, returned = backend.execute(ps, mapping, p['height'], p['seconds_per_chart'], p['gp_sha256'])
                assert backend.replay(ps, mapping, transcript) == returned
                chart = {'index': total, 'centre_index': ci, 'policy': lane,
                         'mapping': mapping, 'search': transcript}
                path = wd/f'chart-{len(charts):04d}.json'
                write(path, chart)
                charts.append(path.name)
                total += 1
                for point in returned:
                    result = admission.consider(point)
                    if result['status'] == 'INDEPENDENT_FINITE_COLUMN':
                        points = tuple(admission.points)
                        proof = checked_rank(model, points, admission.primes, seed['proof']['no_rational_2_torsion_prime'])
                        grew = True
                        gain_times.append({'rank': len(points), 'search_arm_seconds': time.monotonic()-started,
                                           'chart': total, 'point': list(map(str, point))})
                        write(wd/'gain.json', {'points': [list(map(str, pt)) for pt in points],
                            'proof': proof, 'chart': path.name, 'chart_sha256': sha(path)})
                        break
                write(folder/'progress.json', {'epoch': epoch, 'charts': total, 'rank_lower_bound': len(points),
                    'gain_times': gain_times, 'points': [list(map(str, pt)) for pt in points], 'proof': proof})
                print('V3_COMPARISON', epoch, total, 'rank', len(points), lane, transcript['status'], flush=True)
                if grew:
                    break
            tested.add(tuple(centre['point']))
            if grew or total >= p['max_charts']:
                break
        stages.append({'epoch': epoch, 'charts': charts, 'attempts': attempts,
                       'landscape_seconds': landscape_seconds, 'rank_after': len(points)})
        write(wd/'stage.json', stages[-1])
        if total >= p['max_charts']:
            stop = 'CHART_BOUND'
            break
        if not grew:
            stop = 'NO_GAIN_EPOCH'
            break
    guard(p)
    write(folder/'terminal.json', {'status': 'BOUNDED_V3_FINISHED_PENDING_REPLAY', 'stop_reason': stop,
        'initial_rank': 16, 'rank_lower_bound': len(points), 'charts': total,
        'curve': seed['curve'], 'points': [list(map(str, pt)) for pt in points], 'proof': proof,
        'stages': stages, 'gain_times': gain_times, 'wall_seconds': time.monotonic()-started,
        'input_reads': sorted(READS), 'protocol_sha256': sha(folder/'protocol.json')})


def verify(folder):
    from v3_warm_engine import certified_state
    from memory_rank_certificate import checked_rank
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    started = time.monotonic()
    p = read(folder/'protocol.json')
    data_guard(folder)
    guard(p)
    terminal = read(folder/'terminal.json')
    assert terminal['protocol_sha256'] == sha(folder/'protocol.json')
    initial = read(folder/'seed.json')
    model = tuple(map(F, initial['curve']))
    expected = initial['points']
    count = 0
    for stage in terminal['stages']:
        wd = folder/f"epoch-{stage['epoch']:02d}"
        seed = read(wd/'seed.json')
        assert seed['points'] == expected
        points = tuple(tuple(map(F, pt)) for pt in seed['points'])
        state = certified_state(model, points, seed['proof'])
        selection = read(wd/'landscape/selection.json')
        assert selection['basis'] == seed['points']
        returned_all = set()
        for name in stage['charts']:
            row = read(wd/name)
            assert row['index'] == count
            centre = selection['centres'][row['centre_index']]
            assert row['mapping']['centre'] == centre
            ps = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                     coordinate_policy=row['mapping']['coordinate_policy'])
            returned_all.update(backend.replay(ps, row['mapping'], row['search']))
            count += 1
        if (wd/'gain.json').exists():
            gain = read(wd/'gain.json')
            assert gain['points'][:-1] == expected and tuple(map(F, gain['points'][-1])) in returned_all
            assert gain['chart_sha256'] == sha(wd/gain['chart'])
            expected = gain['points']
        assert len(expected) == stage['rank_after']
    assert count == terminal['charts'] and expected == terminal['points']
    pts = tuple(tuple(map(F, pt)) for pt in expected)
    old = terminal['proof']
    fresh = checked_rank(model, pts, [s['prime'] for s in old['signatures']], old['no_rational_2_torsion_prime'])
    assert fresh == old
    write(folder/'verified.json', {'status': 'PASS_EXACT_POINT_MAP_AND_FINITE_RANK_REPLAY',
        'initial_rank': 16, 'rank_lower_bound': len(pts), 'charts': count,
        'terminal_sha256': sha(folder/'terminal.json'), 'wall_seconds': time.monotonic()-started,
        'scope': 'Every retained square, rational map, adaptive gain source and final finite independence proof. No re-enumeration of point boxes or independent full selector/CVP replay.'})
    print('VERIFIED', len(pts), count, flush=True)


def supervise(folder):
    from research_runtime.supervisor import run, Limits
    p = freeze(folder)
    elapsed = 0
    for mode, maximum in [('prepare', p['preparation_seconds']), ('search', p['total_computation_seconds']), ('verify', 60)]:
        budget = min(maximum, p['total_computation_seconds']-elapsed)
        if budget <= 0:
            break
        record = run([p['map_python'], '-python', str(Path(__file__).resolve()), mode, '--folder', str(folder)],
            limits=Limits(wall_seconds=budget, rss_bytes=p['rss_bytes']),
            log_path=folder/(mode+'.log'), checkpoint_path=folder/(mode+'-supervisor.json'))
        elapsed += record['wall_seconds']
        write(folder/'execution.json', {'last_stage': mode, 'outcome': record['outcome'], 'cumulative_wall_seconds': elapsed})
        print(mode, record['outcome'], record['wall_seconds'], flush=True)
        if record['outcome'] != 'completed':
            break


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['supervise', 'prepare', 'search', 'verify'])
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    args = parser.parse_args()
    globals()[args.mode](args.folder.resolve())
