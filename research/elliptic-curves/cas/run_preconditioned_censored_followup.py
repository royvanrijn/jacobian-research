#!/usr/bin/env python3
"""Alternative maps at the censored centres of a verified, unchanged-basis pass.

Inherit the parent's verified CVPs; search only a frozen censored-centre roster.
Stop at the first certified new direction and hand that basis to later V3.
"""
import argparse
import fcntl
import json
from pathlib import Path
import time
from fractions import Fraction as F
import run_bounded_maps_seed_v3 as base
from preconditioned_map_receipts import obtain
from bounded_map_receipts import obtain as old_obtain

SOURCE = None
original_freeze = base.freeze


def context(source):
    read, sha, require = base.read, base.sha, base.require
    policy, terminal, verified = [read(source/n) for n in ('protocol.json','terminal.json','verified.json')]
    base.guard(policy)
    require(policy['schema'] == 'prepared-parent-bounded-maps.v1' and
            verified['status'] == 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' and
            verified['terminal_sha256'] == sha(source/'terminal.json') and
            verified['rank_lower_bound'] == terminal['rank_lower_bound'] and
            verified['charts'] == terminal['charts'] and
            terminal['protocol_sha256'] == sha(source/'protocol.json'), 'verified bounded-map parent required')
    require(len(terminal['stages']) == 1, 'this adapter requires a single unchanged-basis epoch')
    stage = terminal['stages'][0]
    require(stage['before'] == stage['after'] == terminal['rank_lower_bound'], 'parent basis changed')
    epoch = source/'epoch-00'
    seed, selection, reference = [read(epoch/n) for n in
                                 ('seed.json','landscape/selection.json','reference-landscape/selection.json')]
    require(seed['points'] == terminal['points'] == selection['basis'] and
            read(epoch/'stage.json') == stage and
            read(epoch/'reference-verified.json')['status'] == 'PASS_FULL_PRODUCTIVE_EPOCH_REPLAY' and
            read(epoch/'reference-verified.json')['stage_sha256'] == sha(epoch/'stage.json') and
            stage['selection_sha256'] == sha(epoch/'landscape/selection.json') and
            base.normalized_selection(selection) == base.normalized_selection(reference),
            'parent landscape proof differs')
    paths = [source/n for n in ('protocol.json','terminal.json','verified.json')]
    paths += [epoch/n for n in ('seed.json','stage.json','reference-verified.json',
                                'landscape/selection.json','reference-landscape/selection.json')]
    for suffix in ('-supervision','-replay-supervision'):
        path = source.with_name(source.name+suffix)/'supervisor.json'
        require(read(path)['outcome'] == 'completed', 'parent supervisor incomplete')
        paths.append(path)
    model = tuple(map(F,seed['curve']))
    points = tuple(tuple(map(F,p)) for p in seed['points'])
    state = base.certified_state(model,points,seed['proof'])
    indices = []
    for skipped in stage['skipped']:
        if skipped['policy'] != 'quartic_minimized' or skipped['reason'] not in ('strict_wall_timeout','strict_rss_limit'):
            continue
        ci = skipped['centre_index']
        mapping, receipt_sha, reason = old_obtain(epoch,ci,'quartic_minimized',model,points,
                                                  selection['centres'][ci],state,policy,replay=True)
        require(mapping is None and reason == skipped['reason'] and
                any(a['centre_index'] == ci and a['policy'] == 'quartic_minimized' and
                    a['receipt_sha256'] == receipt_sha for a in stage['map_attempts']), 'parent censor receipt differs')
        indices.append(ci)
        paths.append(epoch/f'map-{ci:04d}-quartic_minimized.json')
        paths += [p for p in (epoch/f'map-{ci:04d}-quartic_minimized').iterdir() if p.is_file()]
    require(0 < len(indices) <= 16 and len(set(indices)) == len(indices), 'bounded nonempty censor roster required')
    completed = {}
    previous = sha(epoch/'landscape/selection.json')
    for j in range(stage['charts']):
        path = epoch/f'chart-{j:04d}.json'
        chart = read(path)
        require(chart['index'] == j and chart['previous_sha256'] == previous, 'parent chart chain differs')
        previous = sha(path)
        paths.append(path)
        if chart['search']['status'] == 'bounded_search_complete':
            require(chart['search']['height_bound'] == 125000, 'parent box bound differs')
            key = tuple(chart['mapping']['centre']['point'])
            completed.setdefault(key,set()).add(base.coordinate_key(chart['mapping']))
    require(previous == stage['last_chart_sha256'], 'parent chart endpoint differs')
    return seed, selection, reference, indices, completed, paths


def freeze(folder):
    seed, selected, reference, indices, completed, paths = context(SOURCE)
    protocol = original_freeze(folder)
    prepared_seed = base.read(folder/'seed.json')
    base.require(prepared_seed['points'] == seed['points'] and prepared_seed['curve'] == seed['curve'],
                 'prepared seed differs from parent')
    protocol.update(schema='preconditioned-censored-centres.v1', parent_run=str(SOURCE.relative_to(base.ROOT)),
                    parent_centre_indices=indices, max_charts=len(indices), max_epochs=1,
                    target_rank=len(seed['points'])+1,
                    selection='Fixed resource-censored centres of the verified parent, in its order. '
                              'Use factor-free preconditioning then full PARI minimization. '
                              'Skip coordinate boxes already completed in the parent. Stop at first certified gain.',
                    scope='Alternative coordinates at previously censored centres, not recovery of the identical '
                          'unknown old map. Inherit parent CVP verification; verify all new maps and point witnesses. '
                          'A new direction is queued for later V3. No global minimality, exact rank or point absence.')
    protocol['inputs'].update({str(p.relative_to(base.ROOT)):base.sha(p) for p in paths})
    for name in ('run_preconditioned_censored_followup.py','preconditioned_map_receipts.py',
                 'preconditioned_map_worker.py','preconditioned_full_pari_mapping.sage'):
        protocol['sources'][str((base.CAS/name).relative_to(base.ROOT))] = base.sha(base.CAS/name)
    base.checkpoint(folder/'protocol.json',protocol)
    return protocol


def landscape(model, points, tested, wd, protocol, bank, *, reference=False):
    start = time.monotonic()
    seed, primary, ref, indices, _, _ = context(base.ROOT/protocol['parent_run'])
    base.require(not tested and list(map(list,points)) == [list(map(F,p)) for p in seed['points']] and
                 indices == protocol['parent_centre_indices'], 'fixed parent basis/roster differs')
    chosen = json.loads(json.dumps(ref if reference else primary))
    chosen['centres'] = [chosen['centres'][i] for i in indices]
    chosen['inherited_parent_centre_indices'] = indices
    wd.mkdir(parents=True,exist_ok=False)
    base.checkpoint(wd/'selection.json',chosen)
    base.checkpoint(wd/'timings.json',dict(landscape_wall_seconds=time.monotonic()-start,
                    reference_solver=False, inherited_independent_CVP_proof=True))
    return chosen


def mapping(epoch, ci, policy, model, points, centre, state, protocol, *, replay=False):
    result, receipt_sha, reason = obtain(epoch,ci,policy,model,points,centre,state,protocol,replay=replay)
    if result is not None:
        # Only six calls in the initial application; prefer rechecking to a mutable global cache.
        _, _, _, _, completed, _ = context(base.ROOT/protocol['parent_run'])
        if base.coordinate_key(result) in completed.get(tuple(centre['point']),set()):
            return None, receipt_sha, 'SAME_COORDINATE_AS_COMPLETED_PARENT'
    return result, receipt_sha, reason


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('search','replay'))
    parser.add_argument('--folder',type=Path,required=True)
    parser.add_argument('--parent',type=Path)
    parser.add_argument('--preparation',type=Path)
    args = parser.parse_args()
    SOURCE = args.parent.resolve() if args.parent else None
    base.PREP = base.BANK = args.preparation.resolve() if args.preparation else None
    base.MAPPERS = {'preconditioned_full':'preconditioned_full_pari_mapping.sage'}
    base.freeze, base.landscape, base.obtain_map = freeze, landscape, mapping
    folder = args.folder.resolve()
    base.require(folder.exists() or (SOURCE is not None and base.PREP is not None), 'new follow-up needs parent and preparation')
    with (folder.parent/(folder.name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        (base.search if args.mode == 'search' else base.replay)(folder)
