#!/usr/bin/env python3
"""Extend a verified V3 prefix by100 invocations, reusing sealed landscapes.

Creates a new output directory. Previous proofs and receipts remain unchanged.
Old point receipts are replayed; only the100 additional invocations are new
searches. Unchanged full reference landscapes inherit their prior independent
CVP verification through exact input and source bindings. A gain rebuilds as
usual, and new landscapes receive full independent rational-CVP replay.
"""
import argparse
import fcntl
from pathlib import Path
import shutil

import run_lean_preconditioned_seed_v3 as parent
from lean_preconditioned_map_receipts import obtain as fresh_map

SOURCE = None


def freeze(folder):
    source = SOURCE
    policy = parent.read(source/'protocol.json')
    parent.guard(policy)
    terminal = parent.read(source/'terminal.json')
    verified = parent.read(source/'verified.json')
    parent.require(policy['schema'] in ('prepared-parent-lean-preconditioned-maps.v1',
                                      'prepared-parent-lean-cached-continuation.v1'),
                   'verified short pass or cached continuation required')
    parent.require(terminal['protocol_sha256'] == parent.sha(source/'protocol.json') and
                   verified['status'] == 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' and
                   verified['terminal_sha256'] == parent.sha(source/'terminal.json') and
                   verified['rank_lower_bound'] == terminal['rank_lower_bound'] and
                   verified['charts'] == terminal['charts'], 'parent verification seal differs')
    parent.require(terminal['stop_reason'] == 'CHART_BUDGET_EXHAUSTED' and
                   terminal['charts'] == policy['max_charts'], 'complete budget prefix required')
    evidence = [source/n for n in ('protocol.json', 'terminal.json', 'verified.json')]
    for suffix in ('-supervision', '-replay-supervision'):
        path = source.with_name(source.name+suffix)/'supervisor.json'
        supervisor = parent.read(path)
        parent.require(supervisor['outcome'] == 'completed' and supervisor['returncode'] == 0,
                       'completed parent supervisors required')
        evidence.append(path)
    for stage in terminal['stages']:
        epoch = source/f"epoch-{stage['epoch']:02d}"
        parent.require(parent.read(epoch/'stage.json') == stage, 'parent stage differs')
        replay = parent.read(epoch/'reference-verified.json')
        parent.require(replay['status'] == 'PASS_FULL_PRODUCTIVE_EPOCH_REPLAY' and
                       replay['stage_sha256'] == parent.sha(epoch/'stage.json'),
                       'parent epoch replay seal differs')
        saved = epoch/'landscape/selection.json'
        reference = epoch/'reference-landscape/selection.json'
        parent.require(stage['selection_sha256'] == parent.sha(saved) and
                       parent.normalized_selection(parent.read(saved)) ==
                       parent.normalized_selection(parent.read(reference)),
                       'verified reference landscape differs')
        evidence.extend((epoch/'stage.json', epoch/'reference-verified.json', saved, reference))
    folder.mkdir(exist_ok=False)
    for name in ('seed.json', 'bank.json'):
        shutil.copyfile(source/name, folder/name)
    for stage in terminal['stages']:
        name = f"epoch-{stage['epoch']:02d}"
        shutil.copytree(source/name, folder/name)
    # Only this NEW directory's truncated epoch audits must be regenerated.
    last = folder/f"epoch-{terminal['stages'][-1]['epoch']:02d}"
    for name in ('stage.json', 'reference-verified.json', 'cloud.json', 'mod2.json', 'modl.json'):
        (last/name).unlink(missing_ok=True)
    map_locations = dict(policy.get('inherited_map_locations', {}))
    map_evidence = []
    for stage in terminal['stages']:
        epoch_name = f"epoch-{stage['epoch']:02d}"
        for attempt in stage['map_attempts']:
            map_name = f"map-{attempt['centre_index']:04d}-{attempt['policy']}"
            key = epoch_name + '/' + map_name
            origin = parent.ROOT/map_locations[key] if key in map_locations else source/epoch_name
            receipt = origin/(map_name+'.json')
            parent.require(parent.sha(receipt) == attempt['receipt_sha256'], 'inherited map receipt differs')
            map_locations[key] = str(origin.relative_to(parent.ROOT))
            map_evidence.append(receipt)
            map_evidence.extend(p for p in sorted((origin/map_name).rglob('*')) if p.is_file())
    evidence.extend(map_evidence)
    policy['inherited_map_locations'] = map_locations
    immutable = [folder/'seed.json', folder/'bank.json']
    for stage in terminal['stages']:
        epoch = folder/f"epoch-{stage['epoch']:02d}"
        immutable.append(epoch/'seed.json')
        immutable.extend(sorted(epoch.glob('chart-*.json')))
        immutable.extend(p for p in sorted(epoch.glob('map-*')) if p.is_file())
        immutable.extend(p for d in sorted(epoch.glob('map-*')) if d.is_dir() for p in sorted(d.rglob('*')) if p.is_file())
        for subdir in ('landscape', 'reference-landscape'):
            immutable.extend(p for p in sorted((epoch/subdir).rglob('*')) if p.is_file())
    # Bind both ends of every copied immutable file, not just a cached status.
    for path in immutable:
        original = source/path.relative_to(folder)
        parent.require(parent.sha(original) == parent.sha(path), 'copied evidence differs')
        evidence.append(original)
    policy.update(schema='prepared-parent-lean-cached-continuation.v1',
                  max_charts=terminal['charts']+100,
                  search_wall_limit_seconds=1800, replay_wall_limit_seconds=1800,
                  parent_run=str(source.relative_to(parent.ROOT)),
                  inherited_charts=terminal['charts'], additional_chart_budget=100)
    policy['inputs'].update({str(p.relative_to(parent.ROOT)): parent.sha(p)
                             for p in evidence+immutable})
    policy['sources'][str(Path(__file__).resolve().relative_to(parent.ROOT))] = parent.sha(Path(__file__))
    policy['scope'] += (' Cached continuation in a new directory:100 additional invocations. '
                        'Inherited reference landscapes were independently verified in the bound parent; '
                        'all retained point receipts and the extended clouds are replayed again. '
                        'New landscapes after a gain require full independent CVP replay. '
                        'Total chart counts include inherited exposure.')
    parent.checkpoint(folder/'protocol.json', policy)
    parent.guard(policy)
    return policy


def obtain_map(epoch, ci, policy, model, points, centre, state, protocol, *, replay=False):
    key = epoch.name + f'/map-{ci:04d}-{policy}'
    origin = protocol.get('inherited_map_locations', {}).get(key)
    if origin is not None:
        # Preserve the original command/path bindings and never restart an
        # inherited map worker, even if that map was resource-censored.
        return fresh_map(parent.ROOT/origin, ci, policy, model, points, centre,
                         state, protocol, replay=True)
    return fresh_map(epoch, ci, policy, model, points, centre, state, protocol, replay=replay)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'search', 'replay'))
    parser.add_argument('--folder', required=True, type=Path)
    parser.add_argument('--parent', type=Path)
    args = parser.parse_args()
    SOURCE = args.parent.resolve() if args.parent else None
    folder = args.folder.resolve()
    parent.require(folder.exists() or SOURCE is not None, 'new continuation needs --parent')
    parent.freeze = freeze
    parent.obtain_map = obtain_map
    if folder.exists():
        parent.require(parent.read(folder/'protocol.json')['schema'] ==
                       'prepared-parent-lean-cached-continuation.v1', 'not a cached continuation')
    with (folder.parent/(folder.name+'.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.mode == 'prepare':
            freeze(folder)
        else:
            (parent.search if args.mode == 'search' else parent.replay)(folder)
