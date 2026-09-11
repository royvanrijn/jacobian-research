#!/usr/bin/env python3
"""One offline runner: first-29 index, deficit-one basins, short-vector controls.

Consumes the completed short-vector-core experiment and its bound closure source.
Does not re-enumerate vectors, start V3, or change any existing experiment.
Commands: run | prepare | resume | status | check. All policies are frozen at
prepare, source/code snapshots are private, and failed stages are never reused.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time

from curve302_short_core_controls import (
    InvalidEvidence, QuotientMap, Vocabulary, band_for_rank, corrected_basins,
    first_direction_index, integer, panel_metrics, primitive, qnorm, require,
    simulate_panel, summarize_panels, track_observed,
)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
SELF = Path(__file__).resolve()
HELPER = CAS/'curve302_short_core_controls.py'
DEFAULT = ROOT/'artifacts/local/elliptic-curves/curve302-short-core-controls-v1'
NAMES = tuple([f'recovered-local-{i:02d}' for i in range(1,5)] +
              [f'recovered-strict-{i:02d}' for i in range(1,4)] +
              [f'residual-strict-{i:02d}' for i in range(1,8)])
AXES = ((NAMES[1], 1), (NAMES[0], 0), (NAMES[3], 3))
STAGES = ('index', 'basins', 'sampling')
SHORT_FILES = ('REPORT.json', 'enumeration.json', 'filtration.json', 'ranks-basins.json',
               'primitive-directions.tsv', 'plan.json')
STRUCTURE_FILES = ('REPORT.json', 'quotient-relations.json', 'trajectories.json')


def read(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'duplicate JSON key: '+key)
            result[key] = value
        return result
    def reject(value):
        raise InvalidEvidence('nonfinite JSON value: '+value)
    return json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=reject)


def sha(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='.'+path.name+'.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as out:
            json.dump(payload, out, sort_keys=True, indent=2, allow_nan=False)
            out.write('\n'); out.flush(); os.fsync(out.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def software():
    import sympy
    return {'python': platform.python_version(), 'sympy': sympy.__version__,
            'machine': platform.machine(), 'byteorder': sys.byteorder}


def discover_source(explicit):
    if explicit:
        return Path(explicit).resolve()
    parents = [ROOT/'artifacts/local/elliptic-curves', ROOT/'artifacts/generated-results/elliptic-curves']
    candidates = []
    for parent in parents:
        if not parent.exists():
            continue
        for file in parent.rglob('enumeration.json'):
            directory = file.parent
            if (directory/'REPORT.json').is_file():
                try:
                    if read(directory/'REPORT.json').get('status') == 'PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS':
                        candidates.append(directory.resolve())
                except (ValueError, OSError):
                    continue
    candidates = sorted(set(candidates))
    require(len(candidates) == 1, f'found {len(candidates)} completed short-vector bundles; pass --source explicitly')
    return candidates[0]


def validate_bundle(short, structure):
    short, structure = Path(short), Path(structure)
    report, source_plan = read(short/'REPORT.json'), read(short/'plan.json')
    require(report.get('status') == 'PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS', 'short-vector source is not complete')
    for stage in ('enumeration', 'filtration', 'ranks-basins'):
        require(report['outputs'][stage] == sha(short/f'{stage}.json'), 'short-vector output hash mismatch: '+stage)
    enumeration = read(short/'enumeration.json')
    require(enumeration.get('status') == 'PASS_COMPLETE_EXACT_ENUMERATION', 'enumeration is not complete')
    digest = sha(short/'primitive-directions.tsv')
    require(digest == enumeration['enumeration_sha256'] == report['enumeration_sha256'], 'vocabulary hash mismatch')
    structure_report = read(structure/'REPORT.json')
    require(structure_report.get('status') == 'PASS_THREE_CLOSURE_EXPERIMENTS', 'closure source is not complete')
    for filename in STRUCTURE_FILES:
        require(source_plan['source_hashes'][filename] == sha(structure/filename), 'closure source differs from enumerator binding: '+filename)
    for stage in ('quotient-relations', 'trajectories'):
        require(structure_report['outputs'][stage] == sha(structure/f'{stage}.json'), 'closure output hash mismatch: '+stage)
    rel, traj = read(structure/'quotient-relations.json'), read(structure/'trajectories.json')
    require(rel.get('status') == 'PASS_QUOTIENT_RELATION_ANALYSIS' and
            traj.get('status') == 'PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED', 'wrong source analysis status')
    require(tuple(rel['direction_ids']) == NAMES == tuple(traj['direction_ids']), 'direction roster/order changed')
    q = tuple(tuple(F(str(x)) for x in row) for row in rel['schur_quotient'])
    require(len(q) == 14 and all(len(r) == 14 for r in q), 'quotient form dimension changed')
    from sympy import Matrix
    qm = Matrix([[str(v) for v in row] for row in q])
    require(qm == qm.T and qm.is_positive_definite, 'quotient form is not symmetric positive definite')
    runs = []
    for raw in traj['runs']:
        name = raw['seed']
        require(name in NAMES, 'unknown seed')
        events = []
        for stage in raw['stages']:
            for item in stage['new']:
                require(item['integral'] is True and integer(item['denominator']) == 1, 'nonintegral trajectory event')
                word = tuple(integer(x) for x in item['quotient_word'])
                require(len(word) == 14 and any(word), 'bad quotient word')
                p = primitive(word)
                require(tuple(item['primitive_quotient_word']) == p, 'primitive word mismatch')
                events.append({'word': word, 'primitive': p, 'epoch': integer(stage['epoch'])})
        final = integer(raw['final_rank'])-17
        require(final == 1+len(events), 'trajectory length/final-rank mismatch')
        runs.append({'seed': name, 'seed_index': NAMES.index(name), 'events': events, 'final_dimension': final})
    require(len(runs) == 14 and len({r['seed'] for r in runs}) == 14 and sum(len(r['events']) for r in runs) == 180,
            'expected 14 distinct seeds and 180 acquisitions')
    require(sorted(r['final_dimension'] for r in runs) == [12]+[14]*13, 'endpoint census changed')
    require(next(r for r in runs if r['final_dimension'] == 12)['seed'] == NAMES[0], 'wrong rank-29 exception')
    filtration, basins = read(short/'filtration.json'), read(short/'ranks-basins.json')
    require(filtration.get('status') == 'PASS_INTRINSIC_AND_OBSERVED_FILTRATIONS' and
            basins.get('status') == 'PASS_COMPLETE_RANK_AND_BASIN_CENSUS', 'source stage status changed')
    cores = filtration['observed_common_integral_cores']
    require([integer(r['quotient_dimension']) for r in cores] == list(range(1,15)), 'incomplete/reordered common-core census')
    for row in cores:
        row['basis'] = [tuple(integer(x) for x in r) for r in row['basis']]
        require(all(len(r) == 14 for r in row['basis']) and len(row['basis']) == integer(row['rank']), 'malformed core basis')
    runs.sort(key=lambda r: r['seed_index'])
    return {'n': 14, 'form': q, 'runs': runs, 'cores': cores, 'basins': basins, 'enumeration': enumeration}


def validate_policy(policy):
    for name in ('panels', 'max_draws', 'rescue_after', 'rescue_scan_limit', 'stage_seconds', 'memory_bytes'):
        require(integer(policy[name]) > 0, 'nonpositive resource/policy value: '+name)
    require(policy['rescue_after'] <= policy['max_draws'], 'rescue begins after draw budget')
    require(policy['panels'] <= 10000, 'panel cap exceeds declared safe bound')
    require(policy['band_rule'] == 'decades_worst_rank_complete_ties_v1', 'unsupported band rule')


def prepare(folder, source=None, structure=None, *, panels=128, master='302-static-core-v1',
            max_draws=4096, rescue_after=64, rescue_scan_limit=4096, stage_seconds=3600, memory_bytes=4*1024**3):
    folder = Path(folder).resolve()
    require(not folder.exists(), 'output folder exists; use resume/check or a fresh --folder')
    short = discover_source(source)
    source_plan = read(short/'plan.json')
    structure = Path(structure).resolve() if structure else Path(source_plan['source']).expanduser().resolve()
    require(structure.is_dir(), 'original closure bundle moved; pass --structure with its byte-identical files')
    data = validate_bundle(short, structure)
    policy = {'panels': panels, 'master_seed': str(master), 'max_draws': max_draws, 'rescue_after': rescue_after,
              'rescue_scan_limit': rescue_scan_limit, 'stage_seconds': stage_seconds, 'memory_bytes': memory_bytes,
              'band_rule': 'decades_worst_rank_complete_ties_v1'}
    validate_policy(policy)
    folder.mkdir(parents=True, exist_ok=False)
    inputs = {}
    for label, origin, filenames in (('short', short, SHORT_FILES), ('structure', structure, STRUCTURE_FILES)):
        for filename in filenames:
            original = origin/filename
            before = sha(original)
            dest = folder/'inputs'/label/filename
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, dest)
            require(sha(dest) == before == sha(original), 'source changed while snapshotting: '+str(original))
            inputs[str(dest.relative_to(folder))] = before
    code = {}
    for src in (SELF, HELPER):
        dest = folder/'code'/src.name
        dest.parent.mkdir(exist_ok=True)
        shutil.copyfile(src, dest)
        code[str(dest.relative_to(folder))] = sha(dest)
        require(sha(src) == sha(dest), 'code changed while snapshotting')
    plan = {'schema': 'curve302-short-core-controls-plan.v1', 'policy': policy, 'inputs': inputs, 'code': code,
            'software': software(), 'source': str(short), 'structure': str(structure), 'stages': list(STAGES),
            'observed_direction_count': integer(data['enumeration']['direction_count']),
            'boundary': 'Retrospective controls on the displayed quotient. No V3, no new points, no fit to simulated outcomes, no propagation theorem.'}
    atomic(folder/'plan.json', plan)
    atomic(folder/'manifest.json', {'plan_sha256': sha(folder/'plan.json')})
    print('SHORT_CORE_CONTROL_PREPARED|panels={}|source={}'.format(panels, short), flush=True)


def guard(folder):
    folder = Path(folder).resolve()
    plan = read(folder/'plan.json')
    require(sha(folder/'plan.json') == read(folder/'manifest.json')['plan_sha256'], 'plan changed')
    require(plan.get('schema') == 'curve302-short-core-controls-plan.v1' and tuple(plan['stages']) == STAGES, 'wrong plan/stages')
    validate_policy(plan['policy'])
    for file in (SELF, HELPER):
        require(sha(file) == plan['code']['code/'+file.name],
                'controller version changed; use the frozen runner '+str(folder/'code'/SELF.name))
    require(software() == plan['software'], 'runtime version changed since prepare')
    for name, digest in {**plan['inputs'], **plan['code']}.items():
        file = (folder/name).resolve()
        require(file.is_relative_to(folder), 'snapshot path escapes output folder')
        require(file.is_file() and sha(file) == digest, 'snapshot changed: '+name)
    return plan


def load_data(folder, *, vocabulary=True):
    folder = Path(folder)
    data = validate_bundle(folder/'inputs/short', folder/'inputs/structure')
    maps, count = track_observed(data['runs'], 14, data['cores'])
    require(count == 194, 'prefix count changed')
    data['maps'] = maps
    if vocabulary:
        wanted = [e['primitive'] for r in data['runs'] for e in r['events']]
        enum = data['enumeration']
        data['vocab'] = Vocabulary.load(folder/'inputs/short/primitive-directions.tsv', 14,
                                        integer(enum['direction_count']), F(enum['bound']), wanted, data['form'])
        require(len(data['vocab']) >= 29, 'real experiment needs at least 29 directions')
        for run in data['runs']:
            for event in run['events']:
                index = data['vocab'].positions[event['primitive']]
                event['band'] = band_for_rank(data['vocab'], index)
        actual = max(qnorm(data['form'], v) for v in wanted)
        require(actual == F(enum['bound']), 'enumeration ceiling differs from actual acquisition maximum')
    return data


def verify_basin_compatibility(data):
    """Verify every source deficit against the actual integral prefixes."""
    by_dim = {r['quotient_dimension']: r for r in data['cores']}
    from curve302_short_core_controls import rational_rank
    for row in data['basins']['basins']:
        d = integer(row['landmark_dimension']); seed = row['seed']
        require(seed in data['maps'] and d-1 in data['maps'][seed], 'basin points to missing prefix')
        state, core = data['maps'][seed][d-1], by_dim[d]['basis']
        deficit = rational_rank([state.image(v) for v in core], 14-state.rank)
        require(deficit == integer(row['deficit_before']), 'source basin deficit changed')
        require(integer(row['core_rank']) == len(core), 'source basin core rank changed')


def stage_index(folder, destination, policy):
    data = load_data(folder)
    result = first_direction_index(data['vocab'])
    result.update({'schema': 'curve302-short-core-index.v1', 'observed_prefixes_rechecked': 194,
                   'source_direction_count': len(data['vocab'])})
    atomic(Path(destination)/'index.json', result)


def stage_basins(folder, destination, policy):
    data = load_data(folder, vocabulary=False)
    verify_basin_compatibility(data)
    result = corrected_basins(data['basins']['basins'])
    result['schema'] = 'curve302-deficit-one-basins.v1'
    atomic(Path(destination)/'basins.json', result)


def stage_sampling(folder, destination, policy):
    destination = Path(destination)
    data = load_data(folder)
    observed = panel_metrics(data['maps'], data['runs'], data['cores'], 14, AXES)
    assignments = [{'seed': r['seed'], 'final_dimension': r['final_dimension'],
                    'events': [{'step': i, 'primitive': e['primitive'], **e['band']} for i, e in enumerate(r['events'])]}
                   for r in data['runs']]
    atomic(destination/'bands.json', {'schema': 'curve302-frozen-height-bands.v1', 'rule': policy['band_rule'],
                                      'definition': 'Uniform primitive extension in observed static-rank decades 1..10, 11..100, 101..1000, etc; no norm ties split.',
                                      'runs': assignments})
    panels = []
    for index in range(policy['panels']):
        panel = simulate_panel(data['vocab'], data['runs'], data['cores'], 14, AXES,
                               policy['master_seed'], index, policy)
        atomic(destination/'panels'/f'{index:05d}.json', panel)
        panels.append(panel)
        print('SHORT_CORE_CONTROL_PANEL|index={}|complete={}'.format(index, panel['all_runs_complete']), flush=True)
    summary = summarize_panels(panels, observed)
    panel_hashes = {f'panels/{i:05d}.json': sha(destination/'panels'/f'{i:05d}.json') for i in range(policy['panels'])}
    result = {'schema': 'curve302-short-vector-sampling.v1',
              'status': 'PASS_REFERENCE_ENSEMBLE' if not summary['censored_runs'] else 'PASS_REFERENCE_ENSEMBLE_WITH_CENSORING',
              'policy': policy, 'observed': observed, 'summary': summary, 'bands_sha256': sha(destination/'bands.json'),
              'panel_hashes': panel_hashes, 'sampling_input_boundary':
              'The sampler reads seeds, frozen per-step height bands and the static vocabulary, not target-core information. Cores enter evaluation only. No outcome-based reselection or restart of censored panels.'}
    atomic(destination/'sampling.json', result)


FUNCTIONS = {'index': stage_index, 'basins': stage_basins, 'sampling': stage_sampling}


def stage_files(folder, stage):
    folder = Path(folder)
    files = [folder/f'{stage}.json']
    if stage == 'sampling':
        result = read(folder/'sampling.json')
        files += [folder/'bands.json']+[folder/name for name in result['panel_hashes']]
        require(sha(folder/'bands.json') == result['bands_sha256'], 'band output changed')
        for name, digest in result['panel_hashes'].items():
            require((folder/name).resolve().is_relative_to(folder.resolve()), 'panel path escape')
            require(sha(folder/name) == digest, 'panel output changed: '+name)
    return {str(f.relative_to(folder)): sha(f) for f in files}


@contextmanager
def locked(folder):
    import fcntl
    path = Path(folder)/'.controller.lock'
    with path.open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise InvalidEvidence('another controller owns this folder') from exc
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def supervise(command, logfile, seconds, memory_bytes):
    """Bound and terminate ONLY our child process group, never other campaigns."""
    import resource
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
    started = time.monotonic()
    with Path(logfile).open('w') as log:
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                start_new_session=True, preexec_fn=limits,
                                env={**os.environ, 'PYTHONUNBUFFERED':'1', 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1'})
        def terminate():
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try: proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL); proc.wait()
        old_handler = signal.getsignal(signal.SIGTERM)
        def interrupted(signum, frame):
            raise KeyboardInterrupt('controller terminated')
        signal.signal(signal.SIGTERM, interrupted)
        try:
            returncode = proc.wait(timeout=seconds)
            outcome = 'COMPLETED' if returncode == 0 else 'FAILED'
        except subprocess.TimeoutExpired:
            terminate(); returncode, outcome = proc.returncode, 'UNKNOWN_TIMEOUT'
        except BaseException:
            terminate(); raise
        finally:
            signal.signal(signal.SIGTERM, old_handler)
    return {'outcome': outcome, 'returncode': returncode, 'wall_seconds': time.monotonic()-started}


def run_stage(folder, stage, plan):
    folder = Path(folder)
    phase = folder/'phases'/stage
    if (phase/'seal.json').exists():
        seal = read(phase/'seal.json')
        require(seal['outputs'] == stage_files(folder, stage), 'sealed stage outputs changed')
        require(seal['receipt_sha256'] == sha(phase/'receipt.json'), 'supervisor receipt changed')
        require(read(phase/'receipt.json')['outcome'] == 'COMPLETED', 'sealed failed stage')
        return
    require(not (phase/'started.json').exists(), 'failed/interrupted unsealed stage; retain evidence and choose a fresh folder')
    atomic(phase/'started.json', {'stage': stage, 'pid': os.getpid(), 'time': time.time()})
    try:
        receipt = supervise([sys.executable, str(folder/'code'/SELF.name), '_stage', '--folder', str(folder), '--stage', stage],
                            phase/'worker.log', plan['policy']['stage_seconds'], plan['policy']['memory_bytes'])
    except BaseException:
        atomic(phase/'failure.json', {'status': 'UNKNOWN_INTERRUPTED'})
        raise
    atomic(phase/'receipt.json', receipt)
    require(receipt['outcome'] == 'COMPLETED', stage+' '+receipt['outcome']+'; inspect '+str(phase/'worker.log'))
    guard(folder)
    atomic(phase/'seal.json', {'receipt_sha256': sha(phase/'receipt.json'), 'outputs': stage_files(folder, stage)})


def summary_text(folder):
    folder = Path(folder)
    index, basins, sampling = (read(folder/f'{s}.json') for s in STAGES)
    literal = index['literal_first_29']
    lines = ['# Curve302 static short-core controls', '',
             f"First 29: rank {literal['rank']}; ambient index {literal['index_in_ambient']}; Smith invariants {literal['smith_invariants']}.",
             'First-29 integer reconstruction witnesses are in index.json. The full tied 29th shell index is reported separately.', '',
             f"Basins: {basins['overall']['nontrivial_cases']} deficit-one cases; {basins['overall']['already_contained_excluded']} automatic cases excluded.",
             'Pooled counts and equal-case averages at both bounds are in basins.json; they are not V3 probabilities.', '',
             f"Sampling: {sampling['summary']['panels']} frozen panels; {sampling['summary']['fully_completed_panels']} fully complete; {sampling['summary']['censored_runs']} censored runs.",
             'No unsuccessful panel was restarted or removed. Brackets below bound the finite empirical fraction under censoring, not sampling error.', '',
             '| Target | Observed all-run dimension | Simulated panels by that dimension |', '|---|---:|---|']
    for row in sampling['summary']['axes']:
        ref = row['all_runs_by_observed_dimension']
        lines.append(f"| {row['axis']} | {row['observed_all_runs_dimension']} | {ref['true']}/{ref['panels']}, unknown {ref['unknown']}; [{ref['empirical_fraction_lower']}, {ref['empirical_fraction_upper']}] |")
    lines += ['', '## Common integral intersections', '', '| Dimension | Cohort | Observed rank | Exact observed core reproduced |', '|---:|---:|---:|---|']
    for row in sampling['summary']['cores']:
        if not row['primary']:
            continue
        ref = row['equals_observed_core']
        lines.append(f"| {row['dimension']} | {row['participants']} | {row['observed_common_rank']} | {ref['true']}/{ref['panels']}; unknown {ref['unknown']} |")
    lines += ['', 'The simulator uses only the known static vocabulary, seeds, rank bands and primitive-extension rule. Cores are evaluated afterwards.',
              'This is retrospective model checking, not a formal significance test, a prospective rank selector, or a propagation theorem.',
              'The adaptive-predictor null and the failure of the intrinsic filtration to reproduce the entire common-core chain remain unchanged.',
              'This controller consumes the sealed exact enumeration; its check recomputes these three analyses, NOT the previous lattice enumeration.', '']
    return '\n'.join(lines)


def execute(folder):
    folder = Path(folder).resolve()
    with locked(folder):
        plan = guard(folder)
        try:
            for stage in STAGES:
                run_stage(folder, stage, plan)
            text = summary_text(folder)
            (folder/'SUMMARY.md').write_text(text)
            status = 'PASS_THREE_STATIC_CORE_CONTROLS'
            if read(folder/'sampling.json')['summary']['censored_runs']:
                status = 'PASS_THREE_STATIC_CORE_CONTROLS_WITH_CENSORED_SAMPLING'
            atomic(folder/'REPORT.json', {'status': status, 'plan_sha256': sha(folder/'plan.json'),
                                          'outputs': {s: sha(folder/f'{s}.json') for s in STAGES},
                                          'summary_sha256': sha(folder/'SUMMARY.md')})
            print(text, flush=True)
        except BaseException as exc:
            atomic(folder/'FAILURE.json', {'status':'UNKNOWN_INCOMPLETE', 'error_type': type(exc).__name__, 'message': str(exc)})
            raise


def check(folder):
    folder = Path(folder).resolve()
    with locked(folder):
        plan = guard(folder)
        report = read(folder/'REPORT.json')
        require(report['status'] in ('PASS_THREE_STATIC_CORE_CONTROLS', 'PASS_THREE_STATIC_CORE_CONTROLS_WITH_CENSORED_SAMPLING'), 'run incomplete')
        for stage in STAGES:
            phase = folder/'phases'/stage
            require((phase/'seal.json').is_file(), 'missing stage seal')
            seal = read(phase/'seal.json')
            require(seal['receipt_sha256'] == sha(phase/'receipt.json'), 'changed stage receipt')
            require(seal['outputs'] == stage_files(folder, stage), 'stage seal mismatch')
            require(report['outputs'][stage] == sha(folder/f'{stage}.json'), 'report hash mismatch')
        require(report['plan_sha256'] == sha(folder/'plan.json') and report['summary_sha256'] == sha(folder/'SUMMARY.md'), 'report binding mismatch')
        # Fresh outputs; a failed verification attempt is retained, not erased
        # by TemporaryDirectory cleanup. Successful scratch outputs are expendable.
        attempts = folder/'checks'
        attempts.mkdir(exist_ok=True)
        tmp = Path(tempfile.mkdtemp(prefix='attempt-', dir=attempts))
        try:
            for stage in STAGES:
                receipt = supervise([sys.executable, str(folder/'code'/SELF.name), '_stage', '--folder', str(folder),
                                     '--stage', stage, '--destination', str(tmp)], tmp/f'{stage}.log',
                                     plan['policy']['stage_seconds'], plan['policy']['memory_bytes'])
                atomic(tmp/f'{stage}-receipt.json', receipt)
                require(receipt['outcome'] == 'COMPLETED', 'deterministic '+stage+' failed; logs retained in '+str(tmp))
                require(stage_files(tmp, stage) == stage_files(folder, stage), 'deterministic output mismatch: '+stage)
            require(summary_text(tmp) == (folder/'SUMMARY.md').read_text(), 'summary recomputation mismatch')
            guard(folder)
        except BaseException as exc:
            atomic(tmp/'FAILURE.json', {'status':'UNKNOWN_VERIFICATION_FAILED', 'message':str(exc)})
            raise
        else:
            shutil.rmtree(tmp)
        atomic(folder/'CHECK.json', {'status':'PASS_DETERMINISTIC_THREE_STAGE_RECOMPUTATION',
                                    'plan_sha256':sha(folder/'plan.json'), 'outputs':report['outputs'],
                                    'boundary':'Recomputes these controls; does not repeat the source enumeration.'})
    print('SHORT_CORE_CONTROLS_CHECK|PASS|three_analyses_recomputed|no_reenumeration', flush=True)


def status(folder):
    folder = Path(folder)
    result = {'status': 'NOT_PREPARED', 'stages': []}
    if (folder/'plan.json').is_file():
        result['status'] = read(folder/'REPORT.json')['status'] if (folder/'REPORT.json').is_file() else 'INCOMPLETE'
        for stage in STAGES:
            phase = folder/'phases'/stage
            state = 'SEALED' if (phase/'seal.json').is_file() else 'STARTED_UNSEALED' if (phase/'started.json').is_file() else 'PENDING'
            result['stages'].append({'stage': stage, 'status': state})
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('run','prepare','resume','check','status','_stage'))
    parser.add_argument('--source', type=Path, help='completed short-vector-core result directory')
    parser.add_argument('--structure', type=Path, help='byte-identical closure-structure bundle if its original path moved')
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    parser.add_argument('--panels', type=int, default=128)
    parser.add_argument('--seed', default='302-static-core-v1')
    parser.add_argument('--max-draws', type=int, default=4096)
    parser.add_argument('--stage-seconds', type=int, default=3600)
    parser.add_argument('--memory-mib', type=int, default=4096)
    parser.add_argument('--stage', choices=STAGES, help=argparse.SUPPRESS)
    parser.add_argument('--destination', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.command == 'status': status(args.folder); return
    if args.command == 'check': check(args.folder); return
    if args.command == '_stage':
        plan = guard(args.folder)
        require(args.stage in STAGES, 'stage required')
        destination = args.destination or args.folder
        FUNCTIONS[args.stage](args.folder, destination, plan['policy'])
        return
    if args.command in ('run','prepare'):
        prepare(args.folder, args.source, args.structure, panels=args.panels, master=args.seed,
                max_draws=args.max_draws, rescue_after=min(64,args.max_draws),
                stage_seconds=args.stage_seconds, memory_bytes=args.memory_mib*1024**2)
        if args.command == 'prepare': return
    execute(args.folder)


if __name__ == '__main__':
    main()
