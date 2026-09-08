#!/usr/bin/env python3
"""Three fixed conic-split fibres: exact seed proofs, then unchanged V3.

One attempt per stage, no population refill and no quartic seed search.
Each proof has120s; each V3 including replay has14400s and3GiB RSS.
An interrupted attempt is preserved and never automatically relaunched.
"""
import argparse
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys

import det1092_funnel as f
from v3_warm_support import atomic, bindings, read, require, sha
from research_runtime.supervisor import run, Limits

SELF = Path(__file__).resolve()
CHECKER = f.CAS/'verify_det1092_funnel_seed.sage'


def freeze(parent, candidates, folder):
    require(not (folder/'protocol.json').exists(), 'already frozen')
    p = f.protocol(parent)
    extraction = read(candidates)
    rows = extraction['rows']
    require(extraction['status'] == 'PASS_EXACT_SPLIT_EXTRACTION', 'unverified extraction')
    require(len(rows) == 3 and len({r['id'] for r in rows}) == 3, 'expected three distinct split fibres')
    require(not extraction['already_in_original_seed_selection'], 'follow-up overlaps baseline')
    require(all(r['conic_splitting'] == 'SPLIT' for r in rows), 'nonsplit row')
    inputs = dict(p['inputs'])
    paths = [parent/'protocol.json', parent/'selection.json', parent/'intake-replay.json',
             candidates, candidates.parent/'protocol.json']
    inputs.update({str(path.relative_to(f.ROOT)): sha(path) for path in paths})
    p.update(schema='det1092-funnel-conic-split-seeds.v1', profile='exact_split_followup',
        domain='det1092-funnel-conic-split-seeds-v1', population=0, maximum_draws=0,
        inputs=inputs, sources=f.source_bindings(), seed_cap=3,
        followup=dict(cases=[r['id'] for r in rows], seed_seconds=120, verify_seconds=120,
                      amplifier_seconds=14400, rss_bytes=3*1024**3, workers=1,
                      maximum_attempts_per_stage=1, campaign_seconds=43920,
                      seed_quartic_charts=0),
        scope='All three exact conic-splitting fibres in the audited10million-address population, outside the original90. One conic-derived P18 per fibre; independently proved seeds enter unchanged V3. Separate follow-up yield; baseline selection unchanged.')
    atomic(folder/'protocol.json', p, immutable=True)
    atomic(folder/'selection.json', dict(status='SEALED_EXACT_SPLIT_FOLLOWUP',
        parameters=3, arithmetic_candidates=rows, seed_inputs=rows,
        protocol_digest=f.digest(f.packed(p))), immutable=True)


def confirm(folder, case):
    from det1092_funnel_worker import checked_protocol, prepare, conic_points, first_seed, seal_seed
    p = checked_protocol(folder)
    selection = read(folder/'selection.json')
    require(selection['protocol_digest'] == f.digest(f.packed(p)), 'selection binding differs')
    row = next(r for r in selection['seed_inputs'] if r['id'] == case)
    require(case in p['followup']['cases'], 'case outside fixed cohort')
    from fractions import Fraction
    s = Fraction(row['parameter'])
    arithmetic = f.Arithmetic(p)
    draw = int(case.removeprefix('funnel-'))
    expected = arithmetic.candidate(arithmetic.record(draw, s.numerator, s.denominator), row['role'])
    require(expected == row, 'equation-only candidate fails recomputation')
    seed = folder/'seeds'/case
    prepared = prepare(seed, row)
    if prepared is None:
        return
    model, base, proof = prepared
    points = conic_points(row)
    atomic(seed/'conic-points.json', [list(map(str, P)) for P in points], immutable=True)
    found = first_seed(model, base, proof, points)
    require(not list(seed.glob('chart-*.json')), 'quartic seed exposure forbidden')
    if found:
        seal_seed(seed, row, model, base, *found,
            evidence=dict(kind='exact_split_followup', sha256=sha(seed/'conic-points.json')))
    else:
        atomic(seed/'result.json', dict(status='CONIC_POINT_INDEPENDENCE_UNRESOLVED',
            parameter=row['parameter'], rank_lower_bound=17, quartic_seed_charts=0), immutable=True)
    print('CONIC_SEED', case, read(seed/'result.json')['status'], flush=True)


def worker(folder):
    p = f.protocol(folder)
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        require(not (folder/'controller.json').exists(), 'one attempt only; retain prior supervision')
        plan = p['followup']
        ledger = dict(status='RUNNING', pid=os.getpid(), stages=[], active_stage=None)
        atomic(folder/'controller.json', ledger)
        sage = shutil.which('sage')
        require(sage is not None, 'Sage unavailable')
        env = {**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1', 'PYTHONUNBUFFERED':'1'}

        def stage(name, command, seconds, terminal):
            ledger['active_stage'] = name
            atomic(folder/'controller.json', ledger)
            result = run(command, limits=Limits(seconds, plan['rss_bytes']),
                log_path=folder/'logs'/f'{name}.log',
                checkpoint_path=folder/'logs'/f'{name}.supervisor.json', cwd=f.ROOT, env=env)
            passed = result['outcome'] == 'completed' and result['returncode'] == 0 and terminal.exists()
            ledger['stages'].append(dict(name=name, passed=passed, supervision=result,
                terminal_sha256=sha(terminal) if passed else None))
            atomic(folder/'controller.json', ledger)
            print('SPLIT_STAGE', name, passed, flush=True)
            return passed

        certified = []
        for case in plan['cases']:
            seed = folder/'seeds'/case
            if not stage('confirm-'+case,
                [sage, '-python', str(SELF), 'confirm', '--directory', str(folder), '--case', case],
                plan['seed_seconds'], seed/'result.json'):
                continue
            if read(seed/'result.json')['status'] != 'CERTIFIED_M18':
                continue
            if stage('verify-'+case, [sage, '-python', str(CHECKER), '--seed', str(seed),
                '--output', str(seed/'standalone-replay.json')], plan['verify_seconds'], seed/'standalone-replay.json'):
                certified.append(case)
        for case in certified:
            stage('amplify-'+case,
                [sage, '-python', str(f.CAS/'det1092_funnel_worker.py'), 'amplify',
                 '--directory', str(folder), '--case', case],
                plan['amplifier_seconds'], folder/'amplifiers'/case/'queue.json')
        queues = {path.parent.name: read(path) for path in (folder/'amplifiers').glob('*/queue.json')}
        summary = dict(status='COMPLETE_FINITE_SPLIT_FOLLOWUP' if all(r['passed'] for r in ledger['stages'])
                       else 'PARTIAL_FINITE_SPLIT_FOLLOWUP',
            selected=plan['cases'], independently_certified_m18=certified, queues=queues,
            original_selection_changed=False, quartic_seed_charts=0,
            boundary='Rank lower bounds only; no finite miss is a rank upper bound.')
        atomic(folder/'summary.json', summary, immutable=True)
        ledger.update(status=summary['status'], active_stage=None)
        atomic(folder/'controller.json', ledger)


def launch(folder):
    p = f.protocol(folder)
    require(not (folder/'launch.json').exists(), 'already launched; no retry allocation')
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        atomic(folder/'launch.json', dict(protocol_sha256=sha(folder/'protocol.json'),
               bounds=p['followup'], sources=p['sources']), immutable=True)
        with (folder/'controller.log').open('ab', buffering=0) as log:
            proc = subprocess.Popen([sys.executable, str(SELF), 'worker', '--directory', str(folder)],
                cwd=f.ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        print('LAUNCHED_EXACT_SPLIT_FOLLOWUP', proc.pid, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze', 'launch', 'worker', 'confirm'])
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--parent-run', type=Path)
    parser.add_argument('--candidates', type=Path)
    parser.add_argument('--case')
    a = parser.parse_args()
    folder = a.directory.resolve()
    require(folder.is_relative_to(f.ROOT/'artifacts/local/elliptic-curves'), 'local evidence directory required')
    if a.action == 'freeze':
        require(a.parent_run is not None and a.candidates is not None, 'freeze requires parent and candidates')
        freeze(a.parent_run.resolve(), a.candidates.resolve(), folder)
    elif a.action == 'confirm':
        require(a.case is not None, 'case required')
        confirm(folder, a.case)
    elif a.action == 'worker':
        worker(folder)
    else:
        launch(folder)
