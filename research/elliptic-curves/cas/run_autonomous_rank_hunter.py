#!/usr/bin/env python3
"""Detached, non-AI explore/exploit record hunter for compact R17 elliptic curves.

The controller continuously alternates fresh-fibre seed acquisition with bounded
amplification of empirically hot cascades. Every arithmetic worker is supervised,
independently replayed and sealed before its result can affect scheduling.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from fractions import Fraction as F
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback

import autonomous_rank_hunter_policy as policy
import select_autonomous_r17 as selector
from research_runtime.supervisor import Limits, run as supervise
from v3_warm_support import atomic, read, require, sha

ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
D = LOCAL / 'autonomous-rank-hunter-v1'
SELF = Path(__file__).resolve()
ARITH = CAS / 'autonomous_rank_hunter_arithmetic.py'
SAGE = Path.home() / '.local/bin/sage'


def process(pid):
    try:
        fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
        if fields[0] == 'Z':
            return None
        return {'pid': int(pid), 'start_token': fields[19]}
    except (OSError, IndexError, ProcessLookupError):
        return None


def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BaseException:
        os.close(fd); raise RuntimeError('already running: ' + str(path))
    return fd


def save(path, value, immutable=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic(path, json.loads(json.dumps(value)), immutable=immutable)


def campaign_sources():
    names = (
        'run_autonomous_rank_hunter.py', 'autonomous_rank_hunter_policy.py',
        'select_autonomous_r17.py', 'autonomous_rank_hunter_arithmetic.py',
        'run_fresh6_seed_confirmation_v2.py', 'run_complement_seed_v3.py',
        'run_complement_cached_seed_v3.py', 'reconcile_verified_v3_cloud.py',
        'queue_verified_complement_pass.py', 'research_runtime/supervisor.py',
    )
    return [CAS / n for n in names]


def guard():
    p = read(D / 'protocol.json')
    require(sha(D / 'selection-snapshot.json') == p['selection_snapshot_sha256'], 'selection snapshot changed')
    for name, digest in p['sources'].items():
        require(sha(ROOT / name) == digest, 'campaign source changed: ' + name)
    require(sha(SAGE.resolve()) == p['sage_sha256'], 'Sage launcher changed')
    require(sha(Path('/usr/bin/gp')) == p['gp_sha256'], 'GP binary changed')
    return p


def initialize(workers, max_hours, max_cases, max_point_calls, min_free_gib):
    require(not D.exists(), 'campaign already exists; use status/resume')
    require(1 <= int(workers) <= 16 and 0 < float(max_hours) <= 720 and int(max_cases) > 0 and
            int(max_point_calls) > 0 and int(min_free_gib) >= 1, 'invalid campaign bounds')
    tmp = D.with_name(D.name + '.initializing')
    require(not tmp.exists(), 'stale initialization directory requires review: ' + str(tmp))
    tmp.mkdir(parents=True)
    snapshot = selector.freeze_snapshot(tmp)
    config = {
        'schema': 'autonomous-r17-rank-hunter.v1', 'workers': int(workers),
        'max_hours': float(max_hours), 'max_cases': int(max_cases),
        'max_point_calls': int(max_point_calls), 'min_free_gib': int(min_free_gib),
        'target_rank': 32,
        'policy': 'At least half of ordinary worker capacity explores fresh fibres. Exploitation is 100-call tranches selected by certified momentum/recency; stale fixed-bank work pivots representation or retires. Score/height/family are soft scheduling evidence, never rank exclusions.',
        'claim_boundary': 'Certified subgroup lower bounds only. Selection freshness is relative to the frozen repository/local snapshot. No exact-rank or live-world-record claim is made automatically.',
    }
    save(tmp / 'config.json', config, immutable=True)
    proto = {
        **config, 'selection_snapshot_sha256': sha(tmp / 'selection-snapshot.json'),
        'selection_pool_count': len(snapshot['pool']),
        'sources': {str(p.relative_to(ROOT)): sha(p) for p in campaign_sources()},
        'sage_sha256': sha(SAGE.resolve()), 'gp_sha256': sha(Path('/usr/bin/gp')),
    }
    save(tmp / 'protocol.json', proto, immutable=True)
    save(tmp / 'state.json', {'status': 'CREATED', 'dispatch_index': 0, 'cases': [],
                              'decisions': [], 'point_calls': 0})
    tmp.replace(D)


def phase(command, supervision, wall=1800, rss=3 * 1024**3):
    guard()
    require(not supervision.exists(), 'preserve existing supervision: ' + str(supervision))
    result = supervise(command, limits=Limits(wall_seconds=wall, rss_bytes=rss),
                       log_path=supervision / 'worker.log',
                       checkpoint_path=supervision / 'supervisor.json', cwd=ROOT,
                       env={**os.environ, 'PATH': str(SAGE.parent) + os.pathsep + os.environ.get('PATH', ''),
                            'PYTHONPATH': str(CAS), 'OPENBLAS_NUM_THREADS': '1',
                            'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    require(result['outcome'] == 'completed' and result['returncode'] == 0,
            f"bounded phase failed: {result['outcome']} ({supervision})")
    guard(); return result


def sage(*args):
    return [str(SAGE), '-python', '-u', *map(str, args)]


def _reconcile_v3(run, output):
    command = sage(CAS / 'reconcile_verified_v3_cloud.py', '--folder', run, '--output', output)
    phase(command, run.with_name(run.name + '-reconcile-supervision'), wall=600)
    phase(command, run.with_name(run.name + '-reconcile-replay-supervision'), wall=600)


def _reseed(preparation, packet, output):
    phase(sage(ARITH, 'reseed', '--case', output.parents[3], '--packet', packet,
               '--output', output, '--preparation', preparation),
          output.with_name(output.name + '-supervision'), wall=600)


def _fresh_segment(case, preparation, run, allowance):
    phase(sage(ARITH, 'v3-search', '--case', case, '--packet', preparation,
               '--output', run, '--max-calls', str(allowance)),
          run.with_name(run.name + '-supervision'))
    phase(sage(ARITH, 'v3-replay', '--case', case, '--packet', preparation,
               '--output', run, '--max-calls', str(allowance)),
          run.with_name(run.name + '-replay-supervision'))
    terminal = read(run / 'terminal.json')
    return terminal, terminal['charts']


def _cached_segment(case, parent_run, run):
    runner = CAS / 'run_complement_cached_seed_v3.py'
    phase(sage(runner, 'search', '--folder', run, '--parent', parent_run),
          run.with_name(run.name + '-supervision'))
    phase(sage(runner, 'replay', '--folder', run, '--parent', parent_run),
          run.with_name(run.name + '-replay-supervision'))
    terminal, protocol = read(run / 'terminal.json'), read(run / 'protocol.json')
    return terminal, terminal['charts'] - protocol['inherited_charts']


def _packet_rank(path):
    return int(read(path)['rank_lower_bound'])


def _gain_timeline(run, terminal, reconciliation, inherited, curve_offset):
    """Exact new-call positions from sealed V3 stages plus reconciled saved clouds."""
    offsets, subtotal, timeline = {}, 0, []
    for stage in terminal['stages']:
        offsets[int(stage['epoch'])] = subtotal
        subtotal += int(stage['charts'])
        if int(stage['after']) > int(stage['before']):
            call = subtotal
            if call > inherited:
                timeline.append({'phase': 'complement', 'call': curve_offset + call - inherited,
                                 'before': int(stage['before']), 'after': int(stage['after']),
                                 'source': str((run / f"epoch-{stage['epoch']:02d}/stage.json").relative_to(ROOT))})
    if reconciliation is not None:
        for gain in reconciliation['gains']:
            path = Path(gain['chart']); epoch = int(path.parts[0].split('-')[1])
            number = int(path.stem.split('-')[1]); call = offsets[epoch] + number + 1
            if call > inherited:
                timeline.append({'phase': 'complement-cloud', 'call': curve_offset + call - inherited,
                                 'before': int(gain['after']) - 1, 'after': int(gain['after']),
                                 'source': str((run / path).relative_to(ROOT))})
    return sorted(timeline, key=lambda x: (x['call'], x['after'], x['phase']))


def _run_tranche(case, generation, preparation, tranche, parent_run=None):
    """Spend at most100 new point calls, reconciling complete clouds before continuing."""
    bank = case / 'banks' / f'bank-{generation:02d}'
    runs = bank / 'runs'; runs.mkdir(parents=True, exist_ok=True)
    before_packet = next(preparation.glob('seed-M*.json'))
    before_rank = int(read(parent_run / 'terminal.json')['rank_lower_bound']) if parent_run else _packet_rank(before_packet)
    remaining, used, segment = 100, 0, 0
    current_prep = preparation
    latest_packet = parent_run / 'terminal.json' if parent_run else before_packet
    latest_run = parent_run
    cached_parent = parent_run
    curve_offset = sum(int(e['new_calls']) for e in _case_events(case))
    timeline = []
    while remaining > 0 and _packet_rank(latest_packet) < 32:
        run = runs / f'tranche-{tranche:04d}-segment-{segment:02d}'
        if cached_parent is not None and segment == 0:
            terminal, new_calls = _cached_segment(case, cached_parent, run)
            inherited = int(read(run / 'protocol.json')['inherited_charts'])
        else:
            terminal, new_calls = _fresh_segment(case, current_prep, run, remaining)
            inherited = 0
        require(0 <= new_calls <= remaining, 'tranche call accounting differs')
        used_before = used
        used += new_calls; remaining -= new_calls; latest_run = run
        latest_packet = run / 'terminal.json'; reconciliation = None
        if terminal['stop_reason'] == 'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION':
            reconciled = run / 'reconciled.json'; _reconcile_v3(run, reconciled)
            reconciliation = read(reconciled); latest_packet = reconciled
        timeline.extend(_gain_timeline(run, terminal, reconciliation, inherited, curve_offset + used_before))
        if reconciliation is not None and remaining and _packet_rank(latest_packet) < 32:
            rebuilt = bank / 'rebuilds' / f'tranche-{tranche:04d}-{segment:02d}'
            _reseed(current_prep, reconciled, rebuilt)
            current_prep = rebuilt; segment += 1; cached_parent = None; continue
        break
    after_rank = _packet_rank(latest_packet)
    require(len(timeline) == after_rank - before_rank, 'gain timeline/rank delta differs')
    event = {
        'kind': 'complement-tranche', 'bank_generation': generation, 'tranche': tranche,
        'before_rank': before_rank, 'after_rank': after_rank, 'new_calls': used,
        'directions': after_rank - before_rank, 'gain_timeline': timeline,
        'latest_run': str(latest_run.relative_to(ROOT)),
        'latest_packet': str(latest_packet.relative_to(ROOT)),
        'preparation': str(current_prep.relative_to(ROOT)),
        'finished_unix': time.time(),
    }
    ledger = case / 'ledger'; ledger.mkdir(exist_ok=True)
    save(ledger / f'bank-{generation:02d}-tranche-{tranche:04d}.json', event, immutable=True)
    return event


def _case_events(case):
    return [read(p) for p in sorted((case / 'ledger').glob('bank-*-tranche-*.json'))] if (case / 'ledger').exists() else []


def _latest_packet(case):
    events = _case_events(case)
    if events:
        return ROOT / events[-1]['latest_packet']
    return case / 'seed-reconciled.json'


def _current_generation(case):
    events = _case_events(case)
    return max((int(e['bank_generation']) for e in events), default=-1)


def _bank_preparation(case, generation):
    events = [e for e in _case_events(case) if int(e['bank_generation']) == generation]
    if events:
        return ROOT / events[-1]['preparation']
    return case / 'banks' / f'bank-{generation:02d}' / 'preparation'


def _next_tranche(case, generation):
    return sum(int(e['bank_generation']) == generation for e in _case_events(case))


def _queue_action(run):
    terminal = read(run / 'terminal.json')
    if terminal['stop_reason'] != 'CHART_BUDGET_EXHAUSTED':
        return 'NO_SUFFIX'
    try:
        from queue_verified_complement_pass import record
        return record(run)['next_action']
    except (ValueError, OSError, ArithmeticError):
        return 'NO_SUFFIX'


def build_summary(case):
    row = read(case / 'candidate.json')
    seed = read(case / 'seed-search/terminal.json') if (case / 'seed-search/terminal.json').exists() else None
    cloud = read(case / 'seed-reconciled.json') if (case / 'seed-reconciled.json').exists() else None
    events = _case_events(case)
    rank = int(cloud['rank_lower_bound']) if cloud else 17
    if events:
        rank = int(events[-1]['after_rank'])
    complement_calls = sum(int(e['new_calls']) for e in events)
    gain_calls = [int(g['call']) for e in events for g in e.get('gain_timeline', [])]
    generation = _current_generation(case)
    bank_events = [e for e in events if int(e['bank_generation']) == generation]
    bank_calls = sum(int(e['new_calls']) for e in bank_events)
    prior_calls = sum(int(e['new_calls']) for e in events if int(e['bank_generation']) < generation)
    bank_gains = [int(g['call']) - prior_calls for e in bank_events for g in e.get('gain_timeline', [])]
    bank_last = max(bank_gains, default=None)
    latest_run = ROOT / events[-1]['latest_run'] if events else None
    suffix = _queue_action(latest_run) if latest_run else 'NO_SUFFIX'
    seed_rank = int(cloud['rank_lower_bound']) if cloud else 17
    summary = {
        'id': row['id'], 'family': row['family'], 'parameter': row['parameter'],
        'selection_lane': row.get('selection_lane'), 'stratum': row.get('stratum'),
        'seed_calls': int(seed['charts']) if seed else 0,
        'seed_cloud_directions': max(0, seed_rank - 18) if seed and seed['rank_lower_bound'] >= 18 else 0,
        'rank_lower_bound': rank, 'complement_calls': complement_calls,
        'complement_directions': max(0, rank - seed_rank), 'gain_calls': gain_calls,
        'bank_generation': max(0, generation), 'bank_calls': bank_calls,
        'bank_last_gain_call': bank_last,
        'late_gain': bool(bank_last is not None and bank_last >= max(75, bank_calls - 25)),
        'has_suffix': suffix in ('CONSIDER_REMAINING_SUFFIX', 'REBUILD_FROM_CERTIFIED_GAIN'),
        'suffix_action': suffix,
        'latest_run': str(latest_run.relative_to(ROOT)) if latest_run else None,
        'latest_packet': str(_latest_packet(case).relative_to(ROOT)) if cloud else None,
        'events': events,
    }
    summary['decision'] = policy.continuation_decision(summary) if cloud and seed_rank >= 18 else {
        'action': 'retire', 'priority': -1, 'reason': 'no certified seed'}
    summary['conductor_screen'] = conductor_screen(row, rank)
    save(case / 'summary.json', summary)
    if rank >= 28:
        export_success(case, summary)
    return summary


def _integral_discriminant(model):
    from math import gcd
    q = list(map(F, model)); d = 1
    for x in q:
        d = d * x.denominator // gcd(d, x.denominator)
    weights = (1, 2, 3, 4, 6)
    a1, a2, a3, a4, a6 = [q[i] * d**weights[i] for i in range(5)]
    require(all(x.denominator == 1 for x in (a1, a2, a3, a4, a6)), 'integral scaling failed')
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    delta = -b2*b2*b8 - 8*b4*b4*b4 - 27*b6*b6 + 9*b2*b4*b6
    require(delta.denominator == 1 and delta, 'integral discriminant failed')
    return abs(int(delta))


def conductor_screen(row, rank):
    db = ROOT / 'elliptic-curves/data/research_curves/database.json'
    if rank < 18 or not db.exists():
        return {'status': 'NOT_SCREENED'}
    try:
        curves = read(db)['curves']
        known = [int(r['conductor']) for r in curves
                 if r.get('conductor') is not None and int(r.get('rank_lower_bound', 0)) >= rank]
        if not known:
            return {'status': 'NO_INTERNAL_BENCHMARK'}
        benchmark = min(known); delta = _integral_discriminant(row['model'])
        gap = math.log(delta) - math.log(benchmark)
        status = 'STRONG_CONDUCTOR_CANDIDATE' if gap <= 15 else 'NOT_COMPETITIVE_BY_DISCRIMINANT_SCREEN'
        return {'status': status, 'integral_discriminant_abs': str(delta),
                'benchmark_conductor': str(benchmark), 'log_discriminant_minus_log_benchmark': gap,
                'claim_boundary': 'Cheap integral-discriminant screen only; no conductor factorization or record claim.'}
    except Exception as exc:
        return {'status': 'SCREEN_UNRESOLVED', 'reason': str(exc)}


def export_success(case, summary):
    outdir = ART / 'autonomous_rank_hunter_v1'; outdir.mkdir(parents=True, exist_ok=True)
    packet = read(ROOT / summary['latest_packet'])
    value = {
        'status': 'CERTIFIED_AUTONOMOUS_RANK_HUNTER_RESULT',
        'id': summary['id'], 'family': summary['family'], 'parameter': summary['parameter'],
        'rank_lower_bound': summary['rank_lower_bound'], 'curve': packet['curve'],
        'points': packet['points'], 'proof': packet['proof'],
        'trajectory': summary['events'], 'conductor_screen': summary['conductor_screen'],
        'selection_snapshot_sha256': sha(D / 'selection-snapshot.json'),
        'claim_boundary': 'Exact subgroup certificate. Fresh relative to the frozen campaign exclusion snapshot; not an exact-rank or live external-catalogue novelty claim.'
    }
    path = outdir / f"{summary['id']}.json"
    if path.exists():
        old = read(path)
        if int(old['rank_lower_bound']) >= int(value['rank_lower_bound']):
            return
    save(path, value)


def explore(case):
    if not (case / 'generic/verified.json').exists():
        if not (case / 'generic').exists():
            phase(sage(ARITH, 'generic-prepare', '--case', case), case / 'generic-prepare-supervision', wall=600)
        phase(sage(ARITH, 'generic-replay', '--case', case), case / 'generic-replay-supervision', wall=600)
    seed = case / 'generic/candidate/seed-M17.json'; run = case / 'seed-search'
    runner = CAS / 'run_fresh6_seed_confirmation_v2.py'
    if not (run / 'protocol.json').exists():
        phase(sage(runner, 'prepare', '--folder', run, '--seed', seed), run.with_name('seed-prepare-supervision'), wall=600)
    if not (run / 'terminal.json').exists():
        phase(sage(runner, 'search', '--folder', run), run.with_name('seed-search-supervision'))
    if not (run / 'verified.json').exists():
        phase(sage(runner, 'replay', '--folder', run), run.with_name('seed-replay-supervision'))
    phase(sage(ARITH, 'seed-cloud', '--case', case), case / 'seed-cloud-supervision', wall=600)
    phase(sage(ARITH, 'seed-cloud', '--case', case), case / 'seed-cloud-replay-supervision', wall=600)
    cloud = read(case / 'seed-reconciled.json')
    if cloud['rank_lower_bound'] < 18:
        return build_summary(case)
    prep = case / 'banks/bank-00/preparation'
    if not prep.exists():
        phase(sage(ARITH, 'bank', '--case', case, '--packet', case / 'seed-reconciled.json',
                   '--output', prep, '--generation', '0'), prep.with_name('preparation-supervision'), wall=600)
    _run_tranche(case, 0, prep, 0)
    return build_summary(case)


def exploit(case, action):
    summary = build_summary(case); generation = int(summary['bank_generation'])
    if action == 'new_bank':
        generation += 1; packet = _latest_packet(case)
        prep = case / 'banks' / f'bank-{generation:02d}' / 'preparation'
        if not prep.exists():
            phase(sage(ARITH, 'bank', '--case', case, '--packet', packet, '--output', prep,
                       '--generation', str(generation)), prep.with_name('preparation-supervision'), wall=900)
        _run_tranche(case, generation, prep, _next_tranche(case, generation))
    elif action == 'continue_bank':
        events = [e for e in _case_events(case) if int(e['bank_generation']) == generation]
        require(events, 'continuation without bank history')
        latest = ROOT / events[-1]['latest_run']; q = _queue_action(latest)
        prep = _bank_preparation(case, generation)
        if q == 'REBUILD_FROM_CERTIFIED_GAIN':
            rebuilt = case / 'banks' / f'bank-{generation:02d}' / 'rebuilds' / f"continuation-{_next_tranche(case, generation):04d}"
            _reseed(prep, ROOT / events[-1]['latest_packet'], rebuilt); prep = rebuilt
            _run_tranche(case, generation, prep, _next_tranche(case, generation))
        elif q == 'CONSIDER_REMAINING_SUFFIX':
            _run_tranche(case, generation, prep, _next_tranche(case, generation), parent_run=latest)
        else:
            raise ArithmeticError('policy requested unavailable suffix')
    return build_summary(case)


def case_worker(case, action):
    fd = lock(case / 'case.lock')
    try:
        guard()
        result = explore(case) if action == 'explore' else exploit(case, action)
        print('AUTONOMOUS_CASE', case.name, result['rank_lower_bound'], result['decision'], flush=True)
    finally:
        os.close(fd)


def _dispatch(case, action):
    case_worker(case, action)
    return case.name, action


def _new_case(state):
    snapshot = read(D / 'selection-snapshot.json')
    rows = [read(D / 'cases' / cid / 'candidate.json') for cid in state['cases']]
    row = selector.select_next(snapshot, rows, state['dispatch_index'])
    cid = f"auto-{state['dispatch_index'] + 1:06d}"
    row = {**row, 'id': cid}
    case = D / 'cases' / cid; case.mkdir(parents=True, exist_ok=False)
    save(case / 'candidate.json', row, immutable=True)
    state['dispatch_index'] += 1; state['cases'].append(cid)
    return case


def _summaries(state):
    out = []
    for cid in state['cases']:
        p = D / 'cases' / cid / 'summary.json'
        if p.exists(): out.append(read(p))
    return out


def controller(fd):
    cfg = read(D / 'config.json'); state = read(D / 'state.json')
    resumable = ('CREATED', 'STOP_REQUESTED', 'STOP_TIME_LIMIT', 'STOP_DISK_RESERVE')
    require(state['status'] in resumable, 'state not resumable: ' + state['status'])
    started = state.get('started_unix', time.time()) if state['status'] == 'CREATED' else time.time()
    state.update(status='RUNNING', controller=process(os.getpid()), started_unix=started, active=[])
    save(D / 'state.json', state)
    try:
        with ThreadPoolExecutor(max_workers=cfg['workers']) as pool:
            active = {}
            while True:
                guard()
                if (D / 'STOP').exists(): stop_reason = 'STOP_REQUESTED'
                elif time.time() - started >= cfg['max_hours'] * 3600: stop_reason = 'STOP_TIME_LIMIT'
                elif state['point_calls'] >= cfg['max_point_calls']: stop_reason = 'STOP_POINT_LIMIT'
                elif len(state['cases']) >= cfg['max_cases']: stop_reason = 'STOP_CASE_LIMIT'
                elif shutil.disk_usage(D).free < cfg['min_free_gib'] * 1024**3: stop_reason = 'STOP_DISK_RESERVE'
                else: stop_reason = None
                summaries = _summaries(state)
                if any(int(s['rank_lower_bound']) >= 32 for s in summaries): stop_reason = 'STOP_TARGET32_FOUND'
                state['point_calls'] = sum(int(s.get('seed_calls', 0)) + int(s.get('complement_calls', 0)) for s in summaries)
                busy = {case.name for case, _ in active.values()}
                exploit = [(s['decision']['priority'], s) for s in summaries
                           if s['decision']['action'] in ('continue_bank', 'new_bank') and s['id'] not in busy
                           and not (D / 'cases' / s['id'] / 'QUARANTINED.json').exists()]
                exploit.sort(key=lambda x: (-x[0], -int(x[1]['rank_lower_bound']), x[1]['id']))
                explore_min = policy.exploration_slots(cfg['workers'], [s for _, s in exploit])

                while not stop_reason and len(active) < cfg['workers']:
                    running_explore = sum(task == 'explore' for _, task in active.values())
                    choose_explore = running_explore < explore_min
                    if choose_explore:
                        case = _new_case(state); task = 'explore'
                    elif exploit:
                        _, s = exploit.pop(0); case = D / 'cases' / s['id']; task = s['decision']['action']
                        state['decisions'].append({'unix': time.time(), 'id': s['id'], **s['decision']})
                    else:
                        case = _new_case(state); task = 'explore'
                    future = pool.submit(_dispatch, case, task); active[future] = (case, task)
                    save(D / 'state.json', {**state, 'active': [{'case': c.name, 'task': t} for c, t in active.values()]})

                if not active:
                    state.update(status=stop_reason or 'STOP_NO_WORK', active=[], finished_unix=time.time())
                    save(D / 'state.json', state); return
                done, _ = wait(active, timeout=5, return_when=FIRST_COMPLETED)
                for future in done:
                    case, task = active.pop(future)
                    try:
                        future.result(); build_summary(case)
                    except BaseException as exc:
                        failure = {'status': 'QUARANTINED_WORKER_FAILURE', 'case': case.name,
                                   'task': task, 'unix': time.time(), 'error': repr(exc),
                                   'claim_boundary': 'This case is removed from scheduling. No mathematical conclusion follows from the worker failure.'}
                        save(case / 'QUARANTINED.json', failure)
                        state['decisions'].append({'unix': time.time(), 'id': case.name,
                                                   'action': 'quarantine', 'reason': repr(exc)})
                save(D / 'state.json', {**state, 'active': [{'case': c.name, 'task': t} for c, t in active.values()]})
                if stop_reason and active:
                    continue
    except BaseException:
        state.update(status='STOP_FAILURE_REQUIRES_REVIEW', error=traceback.format_exc(), active=[])
        save(D / 'state.json', state); raise
    finally:
        os.close(fd)


def launch(resume=False):
    guard(); state = read(D / 'state.json')
    if resume:
        require(state['status'] in ('STOP_REQUESTED', 'STOP_TIME_LIMIT', 'STOP_DISK_RESERVE'),
                'state is not safely resumable: ' + state['status'])
        (D / 'STOP').unlink(missing_ok=True)
    fd = lock(D / 'controller.lock')
    try:
        with (D / 'controller.log').open('ab', buffering=0) as log:
            proc = subprocess.Popen([sys.executable, '-u', str(SELF), 'controller', '--lock-fd', str(fd)],
                cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                start_new_session=True, pass_fds=(fd,))
        print('AUTONOMOUS_RANK_HUNTER_LAUNCHED', proc.pid, str(D), flush=True)
    finally:
        os.close(fd)


def status():
    state = read(D / 'state.json')
    token = state.get('controller'); state['controller_alive'] = bool(token and process(token['pid']) == token)
    state['top'] = sorted(_summaries(state), key=lambda s: (-int(s['rank_lower_bound']), -float(s['decision'].get('priority', -1))))[:12]
    print(json.dumps(state, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('start', 'resume', 'controller', 'case-worker', 'status', 'stop'))
    p.add_argument('--workers', type=int, default=4); p.add_argument('--max-hours', type=float, default=168)
    p.add_argument('--max-cases', type=int, default=1000); p.add_argument('--max-point-calls', type=int, default=250000)
    p.add_argument('--min-free-gib', type=int, default=10); p.add_argument('--lock-fd', type=int)
    p.add_argument('--case', type=Path); p.add_argument('--task', choices=('explore', 'continue_bank', 'new_bank'))
    a = p.parse_args()
    if a.action == 'start': initialize(a.workers, a.max_hours, a.max_cases, a.max_point_calls, a.min_free_gib); launch()
    elif a.action == 'resume': launch(True)
    elif a.action == 'controller': controller(a.lock_fd)
    elif a.action == 'case-worker': case_worker(a.case.resolve(), a.task)
    elif a.action == 'status': status()
    else: (D / 'STOP').touch(); print('Stop requested; active work will finish replay before exit.')


if __name__ == '__main__':
    main()
