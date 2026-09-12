#!/usr/bin/env python3
"""Independent V3 witness replay; fixes only tuple/list certificate comparison.

The failed first checker and frozen search adapter are preserved unchanged.
"""
import argparse
import json
from pathlib import Path
import time
import run_constructed_class_v3 as r


def verify(folder):
    from v3_warm_engine import certified_state
    from memory_rank_certificate import checked_rank
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    started = time.monotonic()
    p = r.read(folder/'protocol.json')
    r.data_guard(folder)
    r.guard(p)
    terminal = r.read(folder/'terminal.json')
    assert terminal['protocol_sha256'] == r.sha(folder/'protocol.json')
    initial = r.read(folder/'seed.json')
    model = tuple(map(r.F, initial['curve']))
    expected = initial['points']
    count, complete, backend_seconds = 0, 0, 0.0
    for stage in terminal['stages']:
        wd = folder/f"epoch-{stage['epoch']:02d}"
        seed = r.read(wd/'seed.json')
        assert seed['points'] == expected
        points = tuple(tuple(map(r.F, pt)) for pt in seed['points'])
        state = certified_state(model, points, seed['proof'])
        selection = r.read(wd/'landscape/selection.json')
        assert selection['basis'] == seed['points']
        returned_all = set()
        for name in stage['charts']:
            row = r.read(wd/name)
            assert row['index'] == count
            centre = selection['centres'][row['centre_index']]
            assert row['mapping']['centre'] == centre
            ps = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                     coordinate_policy=row['mapping']['coordinate_policy'])
            returned_all.update(backend.replay(ps, row['mapping'], row['search']))
            count += 1
            complete += row['search']['status'] == 'bounded_search_complete'
            backend_seconds += row['search']['wall_seconds']
        if (wd/'gain.json').exists():
            gain = r.read(wd/'gain.json')
            assert gain['points'][:-1] == expected and tuple(map(r.F, gain['points'][-1])) in returned_all
            assert gain['chart_sha256'] == r.sha(wd/gain['chart'])
            expected = gain['points']
        assert len(expected) == stage['rank_after']
    assert count == terminal['charts'] and expected == terminal['points']
    pts = tuple(tuple(map(r.F, pt)) for pt in expected)
    old = terminal['proof']
    fresh = checked_rank(model, pts, [s['prime'] for s in old['signatures']], old['no_rational_2_torsion_prime'])
    assert json.loads(json.dumps(fresh)) == old
    r.write(folder/'verified-v2.json', {'status': 'PASS_EXACT_POINT_MAP_AND_FINITE_RANK_REPLAY',
        'initial_rank': 16, 'rank_lower_bound': len(pts), 'charts': count,
        'completed_charts': complete, 'backend_search_wall_seconds': backend_seconds,
        'terminal_sha256': r.sha(folder/'terminal.json'), 'wall_seconds': time.monotonic()-started,
        'replayer_sha256': r.sha(__file__),
        'retained_failed_checker': 'verify.log',
        'correction': 'Normalize dataclass tuple fields through JSON before comparing the independently recomputed certificate to its serialized representation. Search data and arithmetic checks unchanged.',
        'scope': 'Every retained square, rational map, adaptive gain source and final finite independence proof. No re-enumeration of point boxes or independent full selector/CVP replay.'})
    print('VERIFIED', len(pts), count, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--folder', type=Path, required=True)
    verify(p.parse_args().folder.resolve())
