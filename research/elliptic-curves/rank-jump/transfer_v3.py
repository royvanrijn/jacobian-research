"""Frozen-panel adapter of the maintained V3 upper-shell route.

Search and replay are retained from run_constructed_class_v3; only input
transport, cold generic admission and the outer metered deadline change.
"""
import argparse
from pathlib import Path
from fractions import Fraction as F
import time
from transfer_common import *


def prepare(folder):
    from sage.all import ZZ, matrix, pari
    from visibility_generic_bank import enumerate_bank
    import sage.version
    started=time.monotonic(); guard(folder)
    p=read(folder/'protocol.json'); packet=read(folder/'input.json')
    model,points=seed(packet)
    progress(folder,'generic_rank_certificate')
    a,proof=admission(model,points,p['prime_bound'])
    assert proof['rank_lower_bound']==16
    G=matrix(ZZ,[[int(2*F(x)) for x in row] for row in packet['generic_height_gram']])
    assert G.is_positive_definite()
    U=matrix(ZZ,pari(G).qflllgram()).transpose()
    assert abs(U.det())==1
    write(folder/'generic-metric.json',dict(gram=[list(map(int,row)) for row in G.rows()],gram_scale=2,LLL=[list(map(int,row)) for row in U.rows()]))
    write(folder/'seed.json',dict(curve=list(map(str,model)),points=points,rank_lower_bound=16,proof=proof))
    progress(folder,'generic_shell_bank')
    bank=enumerate_bank([list(map(int,row)) for row in (U*G*U.transpose()).rows()],
        [list(map(int,row)) for row in U.rows()],bound=p['generic_ellipsoid_bound'],
        shells=p['scaled_shells'],node_limit=p['generic_bank_node_limit'])
    for row in bank['rows']:
        w=row['word'];assert sum(w[i]*int(G[i,j])*w[j] for i in range(16) for j in range(16))==row['norm']
    bank['generic_gram_scale']=2;write(folder/'bank.json',bank)
    write(folder/'prepared.json',dict(status='PASS_GENERIC_ONLY_PREPARATION',wall_seconds=time.monotonic()-started,
        bank_rows=len(bank['rows']),bank_nodes=bank['nodes'],sage=sage.version.version,pari=str(pari.version()),
        seed_sha256=sha(folder/'seed.json'),bank_sha256=sha(folder/'bank.json')))


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
    guard(folder)
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
                if total >= p['max_charts'] or time.time() >= p['search_deadline_unix']:
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
            if grew or total >= p['max_charts'] or time.time() >= p['search_deadline_unix']:
                break
        stages.append({'epoch': epoch, 'charts': charts, 'attempts': attempts,
                       'landscape_seconds': landscape_seconds, 'rank_after': len(points)})
        write(wd/'stage.json', stages[-1])
        if total >= p['max_charts'] or time.time() >= p['search_deadline_unix']:
            stop = 'CHART_OR_TIME_BOUND'
            break
        if not grew:
            stop = 'NO_GAIN_EPOCH'
            break
    write(folder/'terminal.json', {'status': 'BOUNDED_V3_FINISHED_PENDING_REPLAY', 'stop_reason': stop,
        'initial_rank': 16, 'rank_lower_bound': len(points), 'charts': total,
        'curve': seed['curve'], 'points': [list(map(str, pt)) for pt in points], 'proof': proof,
        'stages': stages, 'gain_times': gain_times, 'wall_seconds': time.monotonic()-started,
        'protocol_sha256': sha(folder/'protocol.json')})


def verify(folder):
    from v3_warm_engine import certified_state
    from memory_rank_certificate import checked_rank
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    started = time.monotonic()
    p = read(folder/'protocol.json')
    guard(folder)
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



if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','search','verify']);ap.add_argument('--folder',type=Path,required=True)
    args=ap.parse_args();globals()[args.mode](args.folder.resolve())
