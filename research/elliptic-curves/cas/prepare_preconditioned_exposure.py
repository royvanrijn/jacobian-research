#!/usr/bin/env python3
"""Freeze a short-pass landscape and compare its centres with named old runs.

No point search. The comparison is exact up to sign on the same equation;
it does not assert absence of overlaps with unlisted historical runs.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import run_preconditioned_seed_v3 as runner


def run(preparation, folder, histories):
    p = runner
    p.PREP = p.BANK = preparation
    protocol = runner.freeze(folder)
    protocol['sources'][str(Path(__file__).resolve().relative_to(p.ROOT))] = p.sha(Path(__file__))
    protocol['inputs'].update({str(path.relative_to(p.ROOT)):p.sha(path) for path in histories})
    p.checkpoint(folder/'protocol.json', protocol)
    seed = p.read(folder/'seed.json')
    epoch = folder/'epoch-00'
    epoch.mkdir()
    selection = p.landscape(tuple(map(F, seed['curve'])),
                            tuple(tuple(map(F, q)) for q in seed['points']), set(),
                            epoch/'landscape', protocol, p.read(folder/'bank.json'))
    key = lambda q: (F(q[0]), abs(F(q[1])))
    selected = {key(c['point']) for c in selection['centres']}
    rows = []
    old_all = set()
    for path in histories:
        data = p.read(path)
        p.require(data['curve'] == seed['curve'], 'comparison requires the identical equation')
        old = {key((c['search']['base_point']['x'], c['search']['base_point']['y']))
               for c in data['charts']}
        old_all.update(old)
        rows.append(dict(path=str(path.relative_to(p.ROOT)), sha256=p.sha(path),
                         distinct_centres=len(old), overlap_up_to_sign=len(old & selected)))
    p.guard(protocol)
    report = dict(status='PASS_EXACT_CENTRE_EXPOSURE_COMPARISON', point_searches=0,
                  selected_centres=len(selected), compared_old_centres=len(old_all),
                  overlap_up_to_sign=len(selected & old_all), histories=rows,
                  overlap_indices=[i for i, c in enumerate(selection['centres']) if key(c['point']) in old_all],
                  selection_sha256=p.sha(epoch/'landscape/selection.json'),
                  claim_boundary='Same-equation rational centre comparison, only against named histories. '
                  'No point absence, rank claim or independent CVP replay is asserted by this comparison.')
    p.checkpoint(folder/'exposure-preflight.json', report)
    print('EXPOSURE_PREFLIGHT', len(selected), len(old_all), report['overlap_up_to_sign'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preparation', required=True, type=Path)
    parser.add_argument('--folder', required=True, type=Path)
    parser.add_argument('--history', required=True, type=Path, action='append')
    args = parser.parse_args()
    run(args.preparation.resolve(), args.folder.resolve(), [p.resolve() for p in args.history])
