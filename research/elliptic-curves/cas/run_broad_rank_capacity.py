#!/usr/bin/env python3
"""Resume a drained broad-rank campaign with a recorded scheduling capacity.

The original plan, manifest, arithmetic sources and per-fibre budgets stay frozen.
Only the controller's in-memory worker count changes. Each launch retains this
adapter and a hash-bound capacity receipt outside the original runtime snapshot.
"""
import argparse
import fcntl
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_controller(folder):
    """Check the snapshot before importing any of its code."""
    manifest = read(folder/'manifest.json')
    require(sha(folder/'plan.json') == manifest['plan_sha256'], 'plan changed')
    root = Path(read(folder/'plan.json')['root'])
    for name, digest in manifest['files'].items():
        require(sha(root/name) == digest, 'frozen input/source changed: '+name)
    for name, digest in manifest['executables'].items():
        require(sha(name) == digest, 'executable changed: '+name)
    cas = root/'elliptic-curves/cas'
    for name in ('broad_rank_policy', 'broad_rank_runtime'):
        loaded = sys.modules.get(name)
        require(loaded is None or Path(loaded.__file__).resolve() == (cas/(name+'.py')).resolve(),
                'capacity adapter needs a fresh process with frozen imports')
    sys.path.insert(0, str(cas))
    spec = importlib.util.spec_from_file_location('_broad_rank_capacity_controller', cas/'run_broad_rank_search.py')
    controller = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(controller)
    return controller


def run_at_capacity(controller, folder, workers):
    """Use the original controller/lock/STOP logic and original dispatch helpers."""
    require(type(workers) is int and 1 <= workers <= 8, 'workers must be in 1..8')
    original_guard = controller.guard

    def scheduling_guard(path):
        plan, root = original_guard(path)
        return {**plan, 'workers': workers}, root

    # run() reads its scheduling plan once through this module-level name.
    # Worker requests and dispatch hashes still read the unmodified plan file.
    controller.guard = scheduling_guard
    try:
        controller.run(folder)
    finally:
        controller.guard = original_guard


def validate_revision(folder, revision, expected_sha256):
    require(sha(revision/'capacity.json') == expected_sha256, 'capacity receipt changed')
    policy = read(revision/'capacity.json')
    require(policy['schema'] == 'broad-rank.capacity.v1', 'unknown capacity schema')
    require(Path(policy['folder']).resolve() == folder.resolve(), 'wrong campaign')
    require(type(policy['workers']) is int and 1 <= policy['workers'] <= 8, 'invalid capacity')
    for name in ('plan.json', 'manifest.json'):
        require(sha(folder/name) == policy[name+'_sha256'], name+' changed')
    require(sha(revision/'scheduler.py') == policy['adapter_sha256'], 'capacity adapter changed')
    return policy


def launch(folder, workers, reason):
    require(type(workers) is int and 1 <= workers <= 8, 'workers must be in 1..8')
    controller = load_controller(folder)
    plan, root = controller.guard(folder)
    with (folder/'launch.lock').open('a') as lock, (folder/'controller.lock').open('a') as controller_lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fcntl.flock(controller_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        require(not (folder/'STOP').exists(), 'STOP present; drain and inspect before removing it')
        require(not (folder/'COMPLETE.json').exists(), 'campaign already complete')
        for name in ('controller.json', 'launch.json'):
            require(not ((folder/name).exists() and controller.alive(read(folder/name))), name+' is still live')
        active = [p for p in root.glob('broad-cases/*/batch-*/driver-lease.json')
                  if controller.alive(read(p)) and not (p.parent/'seal.json').exists()]
        require(not active, 'bounded jobs still draining')
        revisions = folder/'capacity-revisions'
        revisions.mkdir(exist_ok=True)
        revision = Path(tempfile.mkdtemp(prefix='workers-'+str(workers)+'-', dir=revisions))
        shutil.copyfile(Path(__file__), revision/'scheduler.py')
        controller.write(revision/'capacity.json', {
            'schema': 'broad-rank.capacity.v1', 'folder': str(folder),
            'original_workers': plan['workers'], 'workers': workers, 'reason': reason,
            'plan.json_sha256': sha(folder/'plan.json'),
            'manifest.json_sha256': sha(folder/'manifest.json'),
            'adapter_sha256': sha(revision/'scheduler.py'),
            'scope': 'Scheduling capacity only; original plan and arithmetic snapshot remain unchanged.'})
        capacity_sha256 = sha(revision/'capacity.json')
        validate_revision(folder, revision, capacity_sha256)
        # The launch lock protects the gap before the child takes controller.lock.
        fcntl.flock(controller_lock, fcntl.LOCK_UN)
        with (folder/'controller.log').open('ab', buffering=0) as log:
            process = subprocess.Popen([sys.executable, str(revision/'scheduler.py'), 'run',
                '--folder', str(folder), '--revision', str(revision),
                '--capacity-sha256', capacity_sha256], cwd=root,
                env=controller.environment(plan['sage']), stdin=subprocess.DEVNULL,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        lease = {'pid': process.pid, 'token': controller.token(process.pid), 'workers': workers,
                 'revision': str(revision), 'capacity_sha256': sha(revision/'capacity.json')}
        controller.write(revision/'launch.json', lease)
        controller.write(folder/'launch.json', lease, False)
        controller.write(folder/'capacity.json', lease, False)
        print(json.dumps({'status': 'LAUNCHED', **lease}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('launch', 'run'))
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--workers', type=int)
    parser.add_argument('--revision', type=Path)
    parser.add_argument('--capacity-sha256')
    parser.add_argument('--reason', default='Explicit operator scheduling capacity override')
    args = parser.parse_args()
    folder = args.folder.resolve()
    if args.mode == 'launch':
        require(args.workers is not None, 'launch requires --workers')
        launch(folder, args.workers, args.reason)
    else:
        require(args.revision is not None and args.capacity_sha256,
                'run requires --revision and --capacity-sha256')
        revision = args.revision.resolve()
        require(Path(__file__).resolve() == revision/'scheduler.py', 'run the retained adapter copy')
        policy = validate_revision(folder, revision, args.capacity_sha256)
        run_at_capacity(load_controller(folder), folder, policy['workers'])


if __name__ == '__main__':
    main()
