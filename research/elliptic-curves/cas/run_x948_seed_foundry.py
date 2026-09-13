#!/usr/bin/env python3
"""Frozen, detached X948 fibration diversification; one worker, ten hours.

Construction, control calibration, binary seed exposure and amplification are
separate gates. No replacement parameters, automatic expansion or restart.
"""
import argparse
import csv
import ctypes
import fcntl
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run, _start_token

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
DEFAULT = ROOT/'artifacts/local/elliptic-curves/x948-seed-foundry-v1'


def read(p):
    return json.loads(Path(p).read_text())


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def immutable(p, value):
    if p.exists():
        raise FileExistsError(p)
    checkpoint(p, value)


def prepare(folder):
    if folder.exists():
        raise FileExistsError('preserve the existing experiment')
    started, own = time.monotonic(), time.process_time()
    runtime = folder/'runtime/research'
    runtime.mkdir(parents=True)
    manifest = {}
    # Code only; no old clouds, exceptional points, scores or target catalogue.
    for directory in ('elliptic-curves/cas', 'elliptic-curves/ecsearch', 'elkies-k3/scripts'):
        for path in (ROOT/directory).rglob('*'):
            if path.is_file() and path.suffix in ('.py', '.sage', '.gp', '.c', '.cpp', '.h', '.sh') and '__pycache__' not in path.parts:
                dest = runtime/path.relative_to(ROOT)
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(path.read_bytes())
                manifest[str(dest.relative_to(folder))] = sha(dest)
    for path in (ROOT/'elliptic-curves').glob('*.py'):
        dest = runtime/path.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(path.read_bytes())
        manifest[str(dest.relative_to(folder))] = sha(dest)
    inputs = folder/'inputs'
    inputs.mkdir()
    art = ROOT/'artifacts/generated-results/elliptic-curves'
    native_path = art/'compact_six_r17_atlas_v1.json'
    native = next(p for p in read(native_path)['families'] if p['family'] == '11952')
    keys = ('family', 'A_coefficients_low_to_high', 'B_coefficients_low_to_high',
            'generic_height_gram', 'generic_rank_lower_bound', 'sections')
    source = {k: native[k] for k in keys}
    source.update(source_surface='X948', source_certificate_sha256=sha(native_path))
    immutable(inputs/'source.json', source)
    immutable(inputs/'native-atlas.json', {'families': [source]})
    table = ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv'
    (inputs/'norm8-table.tsv').write_bytes(table.read_bytes())
    proposals = []
    for row in csv.DictReader(table.open(), delimiter='\t'):
        if row['minimal_unoriented_count'] == '1':
            proposals.append({**{k: int(row[k]) for k in ('group_addition_upper_bound', 'support_count',
                'maximum_absolute_coefficient', 'coefficient_l1')}, 'priority': int(row['priority_rank']),
                'trace_word': list(map(int, row['section_basis_w'].split()))})
    immutable(inputs/'a1-proposals.json', {'rows': proposals, 'table_sha256': sha(table)})
    template_path = ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json'
    prior = [{'family': p['fibration_id'], 'A': p['pencil']['A_coefficients_low_to_high'],
              'B': p['pencil']['B_coefficients_low_to_high']} for p in read(template_path)['presentations']]
    history_path = ROOT/'artifacts/local/elliptic-curves/parent-foundry-v4/state.json'
    history = read(history_path)
    retained = {}
    bindings = {str(p.relative_to(ROOT)): sha(p) for p in (native_path, table, template_path, history_path)}
    for family in ('x948-11952-a1-00414', 'x948-11952-a1-00488', 'x948-11952-a1-01070'):
        record = history['parents'][family]
        paths = [ROOT/f'artifacts/local/elliptic-curves/parent-foundry-{v}/runtime/research'/record['path']
                 for v in ('v4', 'v3', 'v2')]
        path = next(p for p in paths if p.exists())
        if sha(path) != record['sha256']:
            raise ArithmeticError('historical parent binding changed')
        p = read(path)
        prior.append({'family': family, 'A': p['raw_A'], 'B': p['raw_B']})
        retained[family] = {k: p[k] for k in keys}
        bindings[str(path.relative_to(ROOT))] = sha(path)
    immutable(inputs/'prior-fibrations.json', {'rows': prior, 'bindings': bindings})
    compact_path = art/'compact_five_mw16_atlas_v1.json'
    base = next(p for p in read(compact_path)['families'] if p['fibration_id'] == 'a1-fibration-01')
    baseline = {k: base[k] for k in keys if k != 'family'}
    baseline['family'] = 'a1-fibration-01'
    for name, data in [('baseline', baseline), ('calibration-parent', retained['x948-11952-a1-00414'])]:
        immutable(runtime/f'parents/{name}.json', data)
    # The original constructor only needs these two generic files at its old paths.
    for src, dest in [(inputs/'native-atlas.json', runtime/native_path.relative_to(ROOT)),
                      (inputs/'norm8-table.tsv', runtime/table.relative_to(ROOT))]:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(src.read_bytes())
    parameters = sorted({Fraction(n,d) for d in range(1,33) for n in range(-32,33)
        if n and math.gcd(n,d) == 1 and 9 <= max(abs(n),d) <= 32},
        key=lambda q: hashlib.sha256(f'x948-fibration-seeds-v1/{q}'.encode()).digest())[:128]
    plan = {
        'schema': 'x948-fibration-seed-foundry.v1', 'created_unix': time.time(),
        'source_commit': subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        'maximum_workers': 1, 'maximum_elapsed_seconds': 36000,
        'rss_bytes': 3*1024**3, 'geometry_job_seconds': 180,
        'geometry_gate_seconds': 1800, 'seed_job_seconds': 300,
        'amplification_job_seconds': 900, 'replay_seconds': 60,
        'calibration_calls': 32, 'stage1_calls': 32, 'stage2_calls': 128,
        'point_height': 125000, 'point_seconds': 10, 'map_seconds': 5,
        'new_fibrations': 8, 'addresses_per_fibration': 128, 'seeds_to_amplify': 16,
        'baseline': 'a1-fibration-01', 'address_height_band': [9,32],
        'parameters': list(map(str, parameters)),
        'calibration': [
            {'id': 'old-mw16-01', 'parent': 'parents/baseline.json', 'parameter': '7/8'},
            {'id': 'old-new-a1-00414', 'parent': 'parents/calibration-parent.json', 'parameter': '1/8'}],
        'calibration_rule': 'Equation and generic basis only. Both old positive controls must independently seed. '
            'Compare factor-free and dual maps under 32 calls and equal whole-job caps. '
            'Among lanes passing both, freeze the lower aggregate complete CPU cost; ties use factor-free. '
            'No new fibration outcome may choose the detector. No automatic repair after a failed gate.',
        'geometry_rule': 'Complete retained 1266-class A1 stratum. Cheapest 64 remaining trace words; '
            '24-entry exact distance-diversified order; first eight inequivalent, rational, '
            'saturated MW16 parents with compact coefficient size at most 256 bits. '
            'Different separating invariants prove inequivalence. Equal invariants are excluded conservatively.',
        'stage1_rule': 'Same 128 addresses, interleaved by address then fibration. '
            'Certify the starting sixteen independently. Stop point searching on first gaining cloud; '
            'reconcile and verify its full return. No-seed, failed preparation and censored exposure stay distinct.',
        'stage2_rule': 'Only after all Stage1 jobs finish and the summary is sealed, freeze first sixteen '
            'certified seeds per fibration in address order (or all if fewer). '
            'Start with the sealed whole seed cloud; 128 further calls each. No selection by first jump size.',
        'stop_rule': 'Stop on fewer than eight admitted parents, failed positive calibration, '
            'source/verification failure, STOP request, disk reserve below 8 GiB or ten-hour ceiling. '
            'Never replace an address or retry a censored job automatically. '
            'Finish at most this first wave; winners require a separately frozen fresh-address launch.',
        'winner_rule': {'maximum_unresolved_fraction': 0.10, 'minimum_seed_count': 8,
            'seed_per_total_CPU_ratio': 1.5, 'minimum_matched_amp_seeds': 8,
            'amp_primary_tail': 23, 'amp_primary_tail_ratio': 2,
            'amp_primary_tail_extra_hits': 3,
            'boundary': 'Practical screening thresholds, not a population theorem. '
                'Amplification comparisons use the same address-ordered prefix size and complete matched exposure. '
                'Ranks 24..32 are descriptive secondary endpoints. No automatic scaling.'},
        'scope': 'Fibration diversification on X948 only. Rates describe fibration, frozen chart and detector jointly; '
            'neither causal intrinsic productivity nor a complete fibration classification is inferred.',
        'sources': manifest, 'input_bindings': bindings,
        'baseline_atlas_sha256': sha(compact_path),
        'preparation_wall_seconds': time.monotonic()-started,
        'preparation_cpu_seconds': time.process_time()-own,
        'sage': shutil.which('sage'), 'sage_sha256': sha(shutil.which('sage')),
        'gp_sha256': sha('/usr/bin/gp')}
    plan['frozen_inputs'] = {str(p.relative_to(folder)): sha(p) for p in inputs.iterdir() if p.is_file()}
    for p in (runtime/'parents').glob('*.json'):
        plan['frozen_inputs'][str(p.relative_to(folder))] = sha(p)
    immutable(folder/'plan.json', plan)
    (folder/'status.sh').write_text('#!/bin/sh\nexec python3 "'+str(folder/'runtime/research/elliptic-curves/cas'/Path(__file__).name)+'" status --folder "'+str(folder)+'"\n')
    (folder/'status.sh').chmod(0o755)
    print('PREPARED', folder, flush=True)


def guard(folder, plan):
    for key in ('sources', 'frozen_inputs'):
        for relative, value in plan[key].items():
            if sha(folder/relative) != value:
                raise ArithmeticError('frozen source/input changed: '+relative)
    if sha(plan['sage']) != plan['sage_sha256'] or sha('/usr/bin/gp') != plan['gp_sha256']:
        raise ArithmeticError('arithmetic executable changed')


class Meter:
    def __init__(self, folder, plan, state):
        self.folder, self.plan, self.state = folder, plan, state
        self.started = time.monotonic()
        self.previous_wall = state.get('charged_wall_seconds', plan['preparation_wall_seconds'])
        self.cpu = min(os.sched_getaffinity(0))
        os.sched_setaffinity(0, {self.cpu})
        os.environ.update({k: '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if ctypes.CDLL(None).prctl(36, 1, 0, 0, 0):
            raise RuntimeError('Linux child subreaper unavailable')

    def remaining(self):
        return self.plan['maximum_elapsed_seconds']-self.previous_wall-(time.monotonic()-self.started)

    def save(self):
        self.state.update(charged_wall_seconds=self.previous_wall+time.monotonic()-self.started, cpu_affinity=self.cpu)
        checkpoint(self.folder/'state.json', self.state)

    def stage(self, name, command, seconds, dest):
        if (self.folder/'STOP').exists():
            raise InterruptedError('USER_STOP')
        if shutil.disk_usage(self.folder).free < 8*1024**3:
            raise InterruptedError('DISK_RESERVE_STOP')
        if self.remaining() < seconds+2:
            raise InterruptedError('METERED_CEILING_STOP')
        guard(self.folder, self.plan)
        self.state['active'] = name
        self.save()
        dest.mkdir(parents=True, exist_ok=True)
        before, own, start = resource.getrusage(resource.RUSAGE_CHILDREN), time.process_time(), time.monotonic()
        record = run(command, limits=Limits(seconds, self.plan['rss_bytes']),
            log_path=dest/'worker.log', checkpoint_path=dest/'supervisor.json')
        while True:
            try:
                os.waitpid(-1, 0)
            except ChildProcessError:
                break
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime+time.process_time()-own
        row = {'name': name, 'wall_seconds': time.monotonic()-start, 'process_tree_cpu_seconds': cpu,
               'outcome': record['outcome'], 'returncode': record['returncode']}
        self.state['stages'].append(row)
        self.state['process_tree_cpu_seconds'] += cpu
        self.save()
        return row


def new_state(plan):
    return {'status': 'PREPARED', 'stages': [], 'admissions': [], 'exposures': [],
            'process_tree_cpu_seconds': plan['preparation_cpu_seconds']}


def admit(folder, plan, state, meter):
    if (folder/'roster.json').exists():
        return read(folder/'roster.json')
    script = folder/'runtime/research/elliptic-curves/cas/x948_seed_foundry_geometry.sage'
    base = [plan['sage'], '-python', str(script)]
    if not (folder/'proposal-order.json').exists():
        receipt = meter.stage('geometry-order', base+['order','--folder',str(folder)], 60, folder/'ordering')
        if receipt['outcome'] != 'completed' or receipt['returncode']:
            raise ArithmeticError('geometry proposal preflight failed')
    order = read(folder/'proposal-order.json')
    known = {p['invariant']['key'] for p in order['prior_invariants']}
    parents, begin = [], time.monotonic()
    for i, row in enumerate(order['rows']):
        if len(parents) == 8 or time.monotonic()-begin+plan['geometry_job_seconds'] > plan['geometry_gate_seconds']:
            break
        dest = folder/'constructions'/f"{i:02d}-{row['priority']:05d}"
        receipt = meter.stage('construct-'+str(row['priority']), base+['construct','--folder',str(folder),'--index',str(i)],
                              plan['geometry_job_seconds'], dest)
        item = {'priority': row['priority'], 'construction_index': i, 'outcome': 'UNRESOLVED_CONSTRUCTION',
                'cost': receipt}
        if receipt['outcome'] == 'completed' and receipt['returncode'] == 0:
            p, a = read(dest/'parent.json'), read(dest/'admission.json')
            if a['parent_sha256'] != sha(dest/'parent.json'):
                raise ArithmeticError('admission seal changed')
            key = p['novelty_invariant']['key']
            if key in known:
                item['outcome'] = 'DUPLICATE_OR_UNRESOLVED_INVARIANT_COLLISION'
            else:
                known.add(key)
                target = folder/'runtime/research/parents'/f"{p['family']}.json"
                target.write_bytes((dest/'parent.json').read_bytes())
                parents.append({'family': p['family'], 'path': str(target.relative_to(folder/'runtime/research')),
                    'sha256': sha(target), 'baseline': False, 'priority': row['priority'],
                    'coefficient_bits': p['compactification']['after_bits'],
                    'admission': str((dest/'admission.json').relative_to(folder)), 'admission_sha256': sha(dest/'admission.json'),
                    'invariant': key})
                item['outcome'] = 'ADMITTED_NEW_SATURATED_MW16'
        state['admissions'].append(item)
        meter.save()
    roster = {'status': 'PASS_EIGHT_NEW_FIBRATIONS' if len(parents) == 8 else 'STOP_INCOMPLETE_ADMISSION',
              'parents': parents, 'proposal_order_sha256': sha(folder/'proposal-order.json'),
              'boundary': 'Inequivalent to the frozen eight old X948 A1 fibrations and one another. '
                          'Marked integral embeddings retained; abstract U-orbit inequivalence not claimed.'}
    if len(parents) == 8:
        path = folder/'runtime/research/parents/baseline.json'
        roster['parents'].append({'family': plan['baseline'], 'path': 'parents/baseline.json',
                                  'sha256': sha(path), 'baseline': True})
    immutable(folder/'roster.json', roster)
    if roster['status'] != 'PASS_EIGHT_NEW_FIBRATIONS':
        raise InterruptedError('STOP_INCOMPLETE_ADMISSION')
    return roster


def exposure(folder, plan, state, meter, name, parent, parameter, phase, mapper, allowance, address=None, packet=None):
    root = folder/'runtime/research'
    job = root/'jobs'/name
    job.mkdir(parents=True, exist_ok=False)
    request = {'parent': parent['path'], 'parent_sha256': parent['sha256'], 'parameter': parameter,
               'phase': phase, 'mapper': mapper, 'allowance': allowance, 'bank_index': 0 if phase == 'seed' else 1}
    if packet:
        request.update(packet=packet['packet'], packet_sha256=packet['packet_sha256'])
    immutable(job/'request.json', request)
    cap = plan['seed_job_seconds'] if phase == 'seed' else plan['amplification_job_seconds']
    script = root/'elliptic-curves/cas/x948_seed_foundry_worker.py'
    receipt = meter.stage(name, [plan['sage'],'-python',str(script),'--job',str(job)], cap, job)
    result = read(job/'result.json') if (job/'result.json').exists() else {'status': 'CENSORED_JOB',
        'rank_lower_bound': None, 'binary_outcome_complete': False, 'calls': None}
    valid = receipt['outcome'] == 'completed' and receipt['returncode'] == 0
    if valid and result['status'] == 'PASS_CERTIFIED_EXPOSURE':
        replay = meter.stage(name+'/independent-replay', [plan['sage'],'-python',str(script),'--job',str(job),'--verify-only'],
                              plan['replay_seconds'], job/'independent-replay')
        valid = replay['outcome'] == 'completed' and replay['returncode'] == 0
        receipt = {**receipt, 'wall_seconds': receipt['wall_seconds']+replay['wall_seconds'],
                   'process_tree_cpu_seconds': receipt['process_tree_cpu_seconds']+replay['process_tree_cpu_seconds']}
        if not valid:
            raise ArithmeticError('independent point verification failed; stop scaling')
    elif not valid and receipt['outcome'] not in ('strict_wall_timeout','strict_rss_limit'):
        raise ArithmeticError('worker failed; preserve checkpoint and stop')
    row = {'name': name, 'family': parent.get('family'), 'parameter': parameter, 'phase': phase,
           'address_index': address, 'mapper': mapper, 'valid': valid, 'result': result, 'cost': receipt,
           'job': str(job.relative_to(folder))}
    state['exposures'].append(row)
    immutable(job/'seal.json', row)
    meter.save()
    return row


def calibration(folder, plan, state, meter):
    rows = []
    for control in plan['calibration']:
        p = read(folder/'runtime/research'/control['parent'])
        parent = {'family': p['family'], 'path': control['parent'], 'sha256': sha(folder/'runtime/research'/control['parent'])}
        for mapper in ('factor_free','dual'):
            rows.append(exposure(folder, plan, state, meter, f"calibration/{control['id']}/{mapper}",
                parent, control['parameter'], 'seed', mapper, plan['calibration_calls']))
    passed = []
    for mapper in ('factor_free','dual'):
        lane = [r for r in rows if r['mapper'] == mapper]
        if all(r['valid'] and (r['result']['rank_lower_bound'] or 0) > 16 for r in lane):
            passed.append((sum(r['cost']['process_tree_cpu_seconds'] for r in lane), mapper))
    result = {'status': 'PASS_POSITIVE_CALIBRATION' if passed else 'STOP_FAILED_POSITIVE_CALIBRATION',
              'mapper': min(passed)[1] if passed else None, 'rows': rows}
    immutable(folder/'calibration.json', result)
    if not passed:
        raise InterruptedError('STOP_FAILED_POSITIVE_CALIBRATION')
    return result['mapper']


def seed_summary(rows, parents):
    answer = []
    for parent in parents:
        data = [r for r in rows if r['family'] == parent['family'] and r['address_index'] is not None and r['phase'] == 'seed']
        seeds = [r for r in data if r['valid'] and (r['result'].get('rank_lower_bound') or 0) > 16]
        complete = sum(r['valid'] and r['result']['binary_outcome_complete'] for r in data)
        cpu = sum(r['cost']['process_tree_cpu_seconds'] for r in data)
        answer.append({'family': parent['family'], 'slots': len(data), 'seeded': len(seeds),
            'complete_binary_outcomes': complete, 'unresolved': len(data)-complete,
            'seed_fraction': len(seeds)/complete if complete else None,
            'CPU_seconds': cpu, 'CPU_per_seed': cpu/len(seeds) if seeds else None,
            'lower_all_slot_fraction': len(seeds)/len(data) if data else None,
            'upper_all_slot_fraction': (len(seeds)+len(data)-complete)/len(data) if data else None})
    return answer


def controller(folder, admission_only=False):
    plan = read(folder/'plan.json')
    state = read(folder/'state.json') if (folder/'state.json').exists() else new_state(plan)
    meter = Meter(folder, plan, state)
    try:
        with (folder/'controller.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
            state['status'] = 'ADMITTING'
            meter.save()
            roster = admit(folder, plan, state, meter)
            if admission_only:
                state['status'] = 'ADMITTED_AWAITING_LAUNCH'
                return
            state['status'] = 'CALIBRATING'
            meter.save()
            mapper = calibration(folder, plan, state, meter)
            state['status'], state['mapper'] = 'STAGE1', mapper
            meter.save()
            for i, t in enumerate(plan['parameters']):
                for parent in roster['parents']:
                    exposure(folder, plan, state, meter, f"stage1/{i:03d}/{parent['family']}",
                        parent, t, 'seed', mapper, plan['stage1_calls'], address=i)
            immutable(folder/'stage1-sealed.json', {'status': 'STAGE1_FINISHED',
                'summary': seed_summary(state['exposures'], roster['parents']),
                'job_seals': {r['job']: sha(folder/r['job']/'seal.json') for r in state['exposures'] if r['address_index'] is not None}})
            seeds = {}
            for p in roster['parents']:
                rows = [r for r in state['exposures'] if r['phase'] == 'seed' and r['address_index'] is not None and
                        r['family'] == p['family'] and r['valid'] and (r['result'].get('rank_lower_bound') or 0) > 16]
                seeds[p['family']] = sorted(rows, key=lambda r:r['address_index'])[:plan['seeds_to_amplify']]
            immutable(folder/'stage2-roster.json', {'stage1_sha256': sha(folder/'stage1-sealed.json'), 'seeds': seeds})
            state['status'] = 'STAGE2'
            meter.save()
            for k in range(plan['seeds_to_amplify']):
                for p in roster['parents']:
                    if k < len(seeds[p['family']]):
                        seed = seeds[p['family']][k]
                        exposure(folder, plan, state, meter, f"stage2/{k:02d}/{p['family']}", p, seed['parameter'],
                            'amplify', mapper, plan['stage2_calls'], address=seed['address_index'], packet=seed['result'])
            state['status'] = 'FIRST_WAVE_FINISHED_NO_AUTOMATIC_EXPANSION'
    except InterruptedError as error:
        state['status'] = str(error)
    except Exception as error:
        state.update(status='STOP_ERROR', error=repr(error))
        raise
    finally:
        state['active'] = None
        if (folder/'roster.json').exists():
            state['seed_summary'] = seed_summary(state['exposures'], read(folder/'roster.json')['parents'])
        meter.save()


def launch(folder):
    plan = read(folder/'plan.json')
    guard(folder, plan)
    if (folder/'launch.json').exists():
        raise FileExistsError('no automatic restart of this run')
    command = [sys.executable, str(folder/'runtime/research/elliptic-curves/cas'/Path(__file__).name),
               'controller', '--folder', str(folder)]
    with (folder/'controller.log').open('w') as log:
        process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
    immutable(folder/'launch.json', {'pid': process.pid, 'token': _start_token(process.pid),
              'command': command, 'plan_sha256': sha(folder/'plan.json'), 'launched_unix': time.time()})
    print('DETACHED', process.pid, folder, flush=True)


def status(folder):
    state = read(folder/'state.json') if (folder/'state.json').exists() else {'status':'PREPARED'}
    launch_record = read(folder/'launch.json') if (folder/'launch.json').exists() else None
    alive = bool(launch_record and _start_token(launch_record['pid']) == launch_record['token'])
    print('Process:', 'running' if alive else 'not running', '| stage:', state['status'])
    print('Active:', state.get('active'))
    print('Charged elapsed hours:', round(state.get('charged_wall_seconds',0)/3600,3),
          '| completed process-tree CPU hours:', round(state.get('process_tree_cpu_seconds',0)/3600,3))
    print('Admitted new fibrations:', sum(r['outcome']=='ADMITTED_NEW_SATURATED_MW16' for r in state.get('admissions',[])))
    if (folder/'roster.json').exists():
        for r in seed_summary(state.get('exposures',[]), read(folder/'roster.json')['parents']):
            print(r['family'], f"{r['seeded']}/{r['complete_binary_outcomes']} resolved seeds; {r['slots']}/128 slots; {r['unresolved']} unresolved")
    print('Details:', folder/'state.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['prepare','admit','launch','controller','status','stop'])
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    args = parser.parse_args()
    folder = args.folder.resolve()
    if args.mode == 'prepare': prepare(folder)
    elif args.mode == 'admit': controller(folder, admission_only=True)
    elif args.mode == 'launch': launch(folder)
    elif args.mode == 'controller': controller(folder)
    elif args.mode == 'status': status(folder)
    else:
        immutable(folder/'STOP', {'requested_unix': time.time()})
        print('Stop requested; the current bounded job drains. No new jobs will start.')
