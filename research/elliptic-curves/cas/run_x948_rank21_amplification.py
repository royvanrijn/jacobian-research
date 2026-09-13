"""One separately authorized, short amplification of the retained X948 rank21 seed.

Reuse the frozen foundry worker. This is a targeted follow-up, not Stage2 of the
unfinished matched experiment. No automatic restart, new parameter or retry.
"""
import argparse
import fcntl
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import time

from run_x948_seed_foundry import ROOT, Meter, guard, immutable, read, sha
from research_runtime.store import checkpoint
from research_runtime.supervisor import _start_token

DEFAULT = ROOT/'artifacts/local/elliptic-curves/x948-rank21-amplification-v1'
SOURCE = ROOT/'artifacts/local/elliptic-curves/x948-seed-foundry-v1'
FAMILY = 'x948-11952-a1-01510'
PACKET_HASH = '1149d043b38174742d15cf23397fa15a2e5663bf050d55c517419c41d6183f0d'
PARENT_HASH = '63b3b4eb23f2c9502d07a5c617be3ba0dc95d9be017782c2cf7b0b8ffd45cf86'


def prepare(folder):
    if folder.exists():
        raise FileExistsError('Preserve the existing run; no automatic restart.')
    started, own = time.monotonic(), time.process_time()
    old = read(SOURCE/'plan.json')
    oldroot = SOURCE/'runtime/research'
    oldjob = oldroot/f'jobs/stage1/035/{FAMILY}'
    assert sha(oldjob/'packet.json') == PACKET_HASH
    assert read(oldjob/'verified.json')['packet_sha256'] == PACKET_HASH
    assert read(oldjob/'result.json')['packet_sha256'] == PACKET_HASH
    assert sha(oldroot/f'parents/{FAMILY}.json') == PARENT_HASH
    root = folder/'runtime/research'
    root.mkdir(parents=True)
    sources = {}
    # Retained code, not a changing checkout or the old run's point-search tree.
    for relative, expected in old['sources'].items():
        source = SOURCE/relative
        assert sha(source) == expected, relative
        target = folder/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        sources[relative] = sha(target)
    script = root/'elliptic-curves/cas'/Path(__file__).name
    script.write_bytes(Path(__file__).read_bytes())
    sources[str(script.relative_to(folder))] = sha(script)
    inputs = {}
    copies = [(oldjob/'packet.json', root/'inputs/seed-M21.json'),
              (oldroot/f'parents/{FAMILY}.json', root/'parents/target.json')]
    copies += [(oldjob/name, root/'inputs/provenance'/name) for name in
               ('request.json', 'verified.json', 'result.json')]
    copies += [(oldjob/'bank-00/derivation.json', root/'inputs/provenance/old-bank.json')]
    for source, target in copies:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        inputs[str(target.relative_to(folder))] = sha(target)
    request = {'parent': 'parents/target.json', 'parent_sha256': PARENT_HASH,
        'packet': 'inputs/seed-M21.json', 'packet_sha256': PACKET_HASH,
        'parameter': '-16/15', 'phase': 'amplify', 'mapper': 'dual',
        'allowance': 128, 'bank_index': 1}
    job = root/'jobs/target'
    immutable(job/'request.json', request)
    inputs[str((job/'request.json').relative_to(folder))] = sha(job/'request.json')
    # Prove that the next generic parity proposals avoid every old bank proposal.
    seed = int(hashlib.sha256(('parent-banks-v1/'+FAMILY).encode()).hexdigest(), 16)
    modulus = 1 << 16
    multiplier, offset = (seed % modulus) | 1, (seed >> 16) % modulus
    masks = lambda a,b: {(multiplier*k+offset) % modulus for k in range(a,b)} - {0}
    previous, following = masks(0,64), masks(512,1024)
    assert not previous & following
    assert {r['mask'] for r in read(oldjob/'bank-00/derivation.json')['candidates']} == previous
    immutable(root/'inputs/parity-disjointness.json', {
        'multiplier': multiplier, 'offset': offset, 'modulus': modulus,
        'old_bank_indices': [0], 'new_bank_indices': list(range(8,16)),
        'old_masks': sorted(previous), 'new_masks': sorted(following),
        'intersection': [], 'boundary': 'Generic anchor classes are disjoint. The enlarged '
            'rank21 basis requires fresh landscapes; the existing runner also deduplicates '
            'identical height boxes under signed coordinate permutations.'})
    inputs[str((root/'inputs/parity-disjointness.json').relative_to(folder))] = sha(root/'inputs/parity-disjointness.json')
    plan = {k: old[k] for k in ('sage', 'sage_sha256', 'gp_sha256', 'rss_bytes')}
    plan.update(schema='x948-targeted-rank21-amplification.v1', created_unix=time.time(),
        family=FAMILY, parameter='-16/15', starting_certified_rank=21,
        maximum_elapsed_seconds=1200, maximum_workers=1, threads=1,
        input_replay_seconds=60, search_seconds=900, final_replay_seconds=60,
        maximum_point_calls=128, point_height=125000, point_seconds=10, map_seconds=5,
        sources=sources, frozen_inputs=inputs, source_plan_sha256=sha(SOURCE/'plan.json'),
        source_packet_sha256=PACKET_HASH, source_parent_sha256=PARENT_HASH,
        policy='Use the retained dual-map worker, fresh generic banks8..15 and the entire '
            'rank21 subgroup. Rebuild after independent gains. All costs and failed attempts count.',
        stop_rule='Stop at128 calls, rank32,900 search seconds, resource failure or STOP. '
            'Retain partial receipts. Final verification has a separate60-second reserve. '
            'No automatic restart, enlargement or additional fibre.',
        selection_boundary='Chosen retrospectively as the strongest fresh seed. This single '
            'targeted follow-up is separate from the incomplete matched Stage1/Stage2 experiment '
            'and cannot estimate amplification rates.',
        preparation_wall_seconds=time.monotonic()-started,
        preparation_cpu_seconds=time.process_time()-own)
    immutable(folder/'plan.json', plan)
    (folder/'status.sh').write_text('#!/bin/sh\nexec python3 "'+str(script)+'" status --folder "'+str(folder)+'"\n')
    (folder/'status.sh').chmod(0o755)
    print('PREPARED', folder, flush=True)


def verify_input(folder):
    from parent_foundry_worker import verify
    root = folder/'runtime/research'
    packet = read(root/'inputs/seed-M21.json')
    assert sha(root/'inputs/seed-M21.json') == PACKET_HASH
    assert sha(root/'parents/target.json') == PARENT_HASH
    assert packet['rank_lower_bound'] == len(packet['points']) == 21
    verify(packet, read(root/'parents/target.json'), '-16/15')
    immutable(folder/'input-verified.json', {'status': 'PASS_TWO_FINITE_IMPLEMENTATIONS',
        'rank_lower_bound': 21, 'packet_sha256': PACKET_HASH, 'parent_sha256': PARENT_HASH})
    print('INPUT_RANK21_VERIFIED', flush=True)


def controller(folder):
    plan = read(folder/'plan.json')
    if (folder/'state.json').exists():
        raise FileExistsError('A controller already started; no automatic restart.')
    state = {'status': 'VERIFYING_INPUT', 'stages': [], 'certified_rank_lower_bound': 21,
        'process_tree_cpu_seconds': plan['preparation_cpu_seconds'], 'started_unix': time.time()}
    meter = Meter(folder, plan, state)
    root = folder/'runtime/research'
    job = root/'jobs/target'
    worker = root/'elliptic-curves/cas/x948_seed_foundry_worker.py'
    try:
        with (folder/'controller.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
            row = meter.stage('input-replay', [plan['sage'], '-python', str(Path(__file__)),
                'verify-input', '--folder', str(folder)], plan['input_replay_seconds'], folder/'input-replay')
            if row['outcome'] != 'completed' or not (folder/'input-verified.json').exists():
                raise InterruptedError('INPUT_VERIFICATION_FAILED')
            state['status'] = 'AMPLIFYING'
            meter.save()
            row = meter.stage('amplification', [plan['sage'], '-python', str(worker),
                '--job', str(job)], plan['search_seconds'], folder/'search-supervision')
            if row['outcome'] != 'completed' or not (job/'result.json').exists():
                raise InterruptedError('SEARCH_STOPPED_PARTIAL_'+row['outcome'])
            state['status'] = 'VERIFYING_OUTPUT'
            row = meter.stage('output-replay', [plan['sage'], '-python', str(worker),
                '--job', str(job), '--verify-only'], plan['final_replay_seconds'], folder/'output-replay')
            if row['outcome'] != 'completed':
                raise InterruptedError('OUTPUT_VERIFICATION_FAILED')
            result = read(job/'result.json')
            assert result['packet_sha256'] == sha(job/'packet.json')
            assert result['initial_rank'] == 21 and result['calls'] <= 128
            state.update(status='FINISHED_NO_AUTOMATIC_CONTINUATION',
                certified_rank_lower_bound=result['rank_lower_bound'], result=result)
            immutable(folder/'seal.json', {'status': 'PASS_INDEPENDENT_OUTPUT_REPLAY',
                'result_sha256': sha(job/'result.json'), 'packet_sha256': sha(job/'packet.json'),
                'plan_sha256': sha(folder/'plan.json'),
                'rank_lower_bound': result['rank_lower_bound'], 'gain': result['rank_lower_bound']-21})
    except InterruptedError as error:
        state['status'] = str(error)
    except Exception as error:
        state.update(status='STOP_ERROR', error=repr(error))
        raise
    finally:
        state['active'] = None
        meter.save()


def launch(folder):
    plan = read(folder/'plan.json')
    guard(folder, plan)
    if (folder/'launch.json').exists():
        raise FileExistsError('No automatic restart.')
    script = folder/'runtime/research/elliptic-curves/cas'/Path(__file__).name
    command = [sys.executable, str(script), 'controller', '--folder', str(folder)]
    with (folder/'controller.log').open('w') as log:
        process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
    immutable(folder/'launch.json', {'pid': process.pid, 'token': _start_token(process.pid),
        'command': command, 'plan_sha256': sha(folder/'plan.json'), 'launched_unix': time.time()})
    print('DETACHED', process.pid, folder, flush=True)


def status(folder):
    state = read(folder/'state.json') if (folder/'state.json').exists() else {'status': 'PREPARED'}
    launch = read(folder/'launch.json') if (folder/'launch.json').exists() else None
    alive = bool(launch and launch['token'] and _start_token(launch['pid']) == launch['token'])
    print('Process:', 'running' if alive else 'not running', '|', state['status'])
    print('Target:', FAMILY, 't=-16/15 | independently sealed rank:', state.get('certified_rank_lower_bound',21))
    elapsed = time.time()-state['started_unix'] if alive and state.get('started_unix') else state.get('charged_wall_seconds',0)
    print('Elapsed seconds:', round(elapsed,1), '| completed process-tree CPU seconds:', round(state.get('process_tree_cpu_seconds',0),1))
    job = folder/'runtime/research/jobs/target'
    invocations = list((job/'point-invocations').glob('*/attempt-*.json'))
    print('Point calls started:', len(invocations), '/128')
    progress = sorted(job.glob('search-*/progress.json'))
    if progress:
        p = read(progress[-1])
        print('Latest online rank (pending final replay):', p['rank_lower_bound'], '|', progress[-1].parent.name)
    if state.get('result'):
        print('Certified gain:', state['result']['rank_lower_bound']-21,
              '| calls:', state['result']['calls'], '| censoring:', state['result']['point_or_map_censoring'])
    print('Details:', folder/'state.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['prepare','launch','controller','verify-input','status','stop'])
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    args = parser.parse_args()
    folder = args.folder.resolve()
    if args.mode == 'stop':
        immutable(folder/'STOP', {'requested_unix': time.time()})
        print('Stop requested; the current bounded stage drains. No automatic restart.')
    else:
        globals()[args.mode.replace('-','_')](folder)
