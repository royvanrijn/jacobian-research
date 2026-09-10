#!/usr/bin/env python3
"""Create a smaller-worker continuation without mutating a sealed foundry.

The new campaign receives the complete scheduler state, generic inputs, every
constructed-parent equation and every live continuation packet.  Historical
job trees stay in the source campaign; their exported receipts are immutable.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from v3_warm_support import atomic, read, require, sha, same_process
import run_parent_foundry as foundry


def copy(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def handoff(source, destination, workers):
    source, destination = source.resolve(), destination.resolve()
    require(1 <= workers <= 4 and not destination.exists(), 'invalid continuation destination or worker count')
    config, root = foundry.guard(source)
    state = read(source/'state.json')
    require(state['status'] == 'STOPPED', 'source foundry must have drained before continuation')
    for role in ('guardian', 'controller'):
        lease = source/(role+'-lease.json')
        if lease.exists():
            row = read(lease)
            require(not same_process(row['pid'], row['token']), 'source '+role+' is still alive')
    destination.mkdir(parents=True)
    newroot = destination/'runtime/research'
    # Source code and generic-only inputs, omitting historical job receipts.
    shutil.copytree(root, newroot, ignore=shutil.ignore_patterns('parent-jobs', '__pycache__'))
    # Retain only packets required for future continuation and constructed
    # parents used as immutable scheduling inputs.
    needed = {p['path'] for p in state['parents'].values() if not p['baseline']}
    needed.update(f['packet'] for f in state['fibres'].values() if f.get('packet'))
    for relative in sorted(needed):
        old = root/relative
        require(old.is_file(), 'missing source continuation input: '+relative)
        copy(old, newroot/relative)
    for name in ('catalogue.json', 'prior-equations.json'):
        copy(source/name, destination/name)
    next_config = {**config, 'workers': workers, 'frozen_root': str(newroot),
                   'export': str(foundry.ART/destination.name),
                   'continued_from': str(source), 'continued_from_manifest_sha256': sha(source/'manifest.json'),
                   'continued_at': time.time()}
    atomic(destination/'config.json', next_config)
    next_state = json.loads(json.dumps(state))
    next_state.update(status='PREPARED', continuation_from=str(source), continued_at=time.time())
    atomic(destination/'state.json', next_state)
    manifest = {'files': {str(p.relative_to(newroot)): sha(p) for p in sorted(newroot.rglob('*')) if p.is_file()},
                'executables': config['sage'] and {str(p): sha(p) for p in (Path(config['sage']), Path('/usr/bin/gp'))},
                'config_sha256': sha(destination/'config.json')}
    atomic(destination/'manifest.json', manifest)
    export = Path(next_config['export']); export.mkdir(parents=True, exist_ok=False)
    for name in ('config.json', 'manifest.json'):
        copy(destination/name, export/name)
    atomic(destination/'handoff.json', {'status': 'PASS_DRAINED_CONTINUATION',
        'source': str(source), 'source_state_sha256': sha(source/'state.json'),
        'source_manifest_sha256': sha(source/'manifest.json'), 'workers_before': config['workers'],
        'workers_after': workers, 'retained_parent_inputs': sorted(needed),
        'boundary': 'Scheduler continuation only; no old job receipt is rewritten or reclassified.'})
    print(json.dumps({'status': 'PREPARED_CONTINUATION', 'folder': str(destination),
                      'workers': workers, 'inputs': len(needed)}))


def wait_for_drain(source, seconds):
    """Wait only for a previously requested graceful drain, never kill work."""
    while read(source/'state.json')['status'] != 'STOPPED':
        time.sleep(seconds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--workers', type=int, required=True)
    parser.add_argument('--wait-for-drain', action='store_true')
    parser.add_argument('--poll-seconds', type=float, default=15.0)
    parser.add_argument('--launch', action='store_true')
    args = parser.parse_args()
    require(args.poll_seconds > 0, 'poll interval must be positive')
    if args.wait_for_drain:
        wait_for_drain(args.source.resolve(), args.poll_seconds)
    handoff(args.source, args.destination, args.workers)
    if args.launch:
        foundry.launch(args.destination.resolve())
