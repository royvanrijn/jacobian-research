#!/usr/bin/env python3
"""Detached four-worker, 60-fibre R17 seed/complement experiment.

freeze -> preflight -> launch; status and resume retain all prior evidence.
No deep suffix, additional candidates, or outcome-dependent budget escalation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from fractions import Fraction
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from v3_warm_support import atomic, read, require, sha
import select_r17_60_panel as selection
from research_runtime.supervisor import run as supervise, Limits

ROOT, CAS, ART = selection.ROOT, Path(__file__).resolve().parent, selection.ART
D = selection.LOCAL / 'r17-60-panel-v1'
SELF = Path(__file__).resolve()
SAGE = Path.home() / '.local/bin/sage'


def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BaseException:
        os.close(fd)
        raise RuntimeError('already running: ' + str(path))
    return fd


def process(pid):
    try:
        fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
        if fields[0] == 'Z':
            return None
        return {'pid': pid, 'start_token': fields[19]}
    except (FileNotFoundError, ProcessLookupError):
        return None


def check_protocol():
    p = read(D / 'protocol.json')
    for category in ('sources', 'inputs'):
        for name, digest in p[category].items():
            require(sha(ROOT/name) == digest, 'changed frozen ' + category + ': ' + name)
    require(sha(SAGE.resolve()) == p['sage_launcher_sha256'], 'Sage launcher changed')
    require(sha(Path('/usr/bin/gp')) == p['gp_sha256'], 'GP changed')
    return p


def freeze():
    require(not (D/'protocol.json').exists(), 'preserve frozen protocol; use preflight/status/resume')
    roster = selection.freeze(D)
    names = ('select_r17_60_panel.py', 'r17_60_arithmetic.py', 'run_r17_60_panel.py',
             'run_fresh6_seed_confirmation_v2.py', 'run_complement_seed_v3.py',
             'reconcile_verified_v3_cloud.py', 'prepare_reconciled_v3_continuation.py',
             'compact_atlas_specialization.py', 'certify_compact_r17_candidates.py',
             'visibility_complement_subset.py', 'visibility_lattice_fast.py',
             'visibility_lattice_v2.py', 'memory_rank_certificate.py', 'future_point_admission.py',
             'lean_preconditioned_map_receipts.py', 'lean_preconditioned_map_worker.py',
             'lean_preconditioned_full_pari_mapping.sage', 'lean_factor_free_pari_mapping.sage',
             'pointed_box_equivalence.py', 'prepare_extended20_mw16_pari_batch.sage',
             'research_runtime/supervisor.py')
    paths = {CAS/n for n in names}
    legacy = read(selection.LOCAL / 'adaptive-visibility-cascade-v3/protocol.json')
    paths.update(ROOT/n for n in legacy['sources'])
    import compact_atlas_specialization as atlas
    inputs = [D/'selection-input.json', D/'roster.json', atlas.ATLAS,
              ART/'r17_exact_maximum_parity_classes_v1.json']
    p = {'schema': 'r17-sixty-seed-complement.v1', 'cases': 60, 'workers': 4,
         'selection': '5 highest retained selection-score rows of height<=1024 and5 of1024<height<=4096 per fibration; exact equation/address exclusions and pairwise Q-isomorphism deduplication. No point outcomes or validation-band scores in selection.',
         'seed_budget': 'All43/49 generic maximum classes; at most86/98 actual point calls; stop on first independently certified M18.',
         'height': 125000, 'point_seconds': 10, 'map_seconds': 5,
         'complement_calls_total': 100, 'target_lower_bound': 32,
         'cloud_rule': 'Replay and reconcile the full winning seed cloud before complement preparation. After every full-cloud reconciliation stop, independently replay and reconcile before rebuilding; remaining complement allowance is100 minus all prior segment calls. No budget reset.',
         'parent_rule': '16 exact maximum generic parity classes outside the original single winning-mask span; ascending mask order extends the span first then fills. Both exact CVP solvers certify every candidate.',
         'stopping': 'First finite-policy exhaustion or target32 or cumulative100 complement calls. Reconciliation with a larger certified basis may rebuild within the remaining allowance. No suffix extension after a budget stop.',
         'reporting': 'Final certified rank, added directions, number of gaining calls, exact calls-per-direction, last gain call, calls since last gain, full gain timeline, and resource/censor outcomes. Report height strata and fibrations separately.',
         'scope': 'Selected finite cohort; no inference of a general rank distribution, exact rank, public novelty, conductor record or cross-family universality.',
         'phase_wall_seconds': 1800, 'phase_rss_bytes': 3*1024**3,
         'sage_launcher_sha256': sha(SAGE.resolve()), 'gp_sha256': sha(Path('/usr/bin/gp')),
         'sources': {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)},
         'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs}}
    atomic(D/'protocol.json', p, immutable=True)
    atomic(ART/'r17_60_panel_protocol_v1.json', p, immutable=True)
    atomic(ART/'r17_60_panel_roster_v1.json', roster, immutable=True)
    print('R17_60_FROZEN', sha(D/'protocol.json'), flush=True)


class ResourceStop(RuntimeError):
    pass


def step(case, phase, index=0, remaining=100):
    p = check_protocol()
    folder = D/'cases'/case if case else D
    name = f'{phase}-{index:02d}'
    atomic(folder/'progress.json', {'phase': phase, 'index': index, 'remaining_calls': remaining,
                                  'started_unix': time.time()})
    command = [str(SAGE), '-python', '-u', str(CAS/'r17_60_arithmetic.py'), phase]
    if case:
        command += ['--case', case, '--index', str(index), '--remaining', str(remaining)]
    record = supervise(command, limits=Limits(wall_seconds=1800 if phase in ('search','replay','seed-search','seed-replay') else 600,
                         rss_bytes=p['phase_rss_bytes']), log_path=folder/'logs'/(name+'.log'),
                       checkpoint_path=folder/'logs'/(name+'.supervisor.json'), cwd=ROOT,
                       env={**os.environ, 'PATH': str(SAGE.parent)+os.pathsep+os.environ['PATH'],
                            'PYTHONPATH': str(CAS), 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    if record['outcome'] in ('strict_wall_timeout', 'strict_rss_limit'):
        raise ResourceStop(phase + ': ' + record['outcome'])
    require(record['outcome'] == 'completed', f'{case}/{phase} failed; see {folder}/logs/{name}.log')
    check_protocol()


def preflight():
    p = check_protocol()
    require(selection.select(read(D/'selection-input.json')) ==
            {k: read(D/'roster.json')[k] for k in ('rows','skipped','stratum_pool_counts')}, 'roster replay differs')
    step(None, 'generic')
    require(len(read(D/'generic/prepared.json')['records']) == 60, 'generic roster incomplete')
    atomic(D/'preflight.json', {'status': 'PASS_SIXTY_GENERIC17_REPLAY', 'point_searches': 0,
                              'protocol_sha256': sha(D/'protocol.json'),
                              'generic_verified_sha256': sha(D/'generic/verified.json')}, immutable=True)
    print('R17_60_PREFLIGHT_COMPLETE|cases=60|rank=17|charts=0|status=PASS', flush=True)


def call_metrics(seed_calls, initial_rank, segment_data, seed_reconciliation):
    timeline = []
    if initial_rank >= 18:
        timeline.append({'phase': 'seed', 'call': seed_calls, 'before': 17, 'after': 18})
    for g in seed_reconciliation:
        timeline.append({'phase': 'seed-cloud', 'call': g['seed_call'], 'before': g['after']-1, 'after': g['after']})
    used = 0
    for terminal, rc in segment_data:
        offsets = {}; subtotal = 0
        for stage in terminal['stages']:
            offsets[stage['epoch']] = subtotal
            subtotal += stage['charts']
            if stage['after'] > stage['before']:
                timeline.append({'phase': 'complement', 'call': used+subtotal,
                                 'before': stage['before'], 'after': stage['after']})
        for g in rc['gains']:
            path = Path(g['chart']); epoch = int(path.parts[0].split('-')[1]); number = int(path.stem.split('-')[1])
            timeline.append({'phase': 'complement-cloud', 'call': used+offsets[epoch]+number+1,
                             'before': g['after']-1, 'after': g['after']})
        used += terminal['charts']
    return timeline, used


def case_run(case):
    fd = lock(D/'cases'/case/'case.lock')
    try:
        folder = D/'cases'/case
        if (folder/'result.json').exists():
            return read(folder/'result.json')
        row = next(r for r in read(D/'roster.json')['rows'] if r['id'] == case)
        if not (folder/'seed-search/protocol.json').exists():
            step(case, 'seed-prepare')
        if not (folder/'seed-search/terminal.json').exists():
            step(case, 'seed-search')
        if not (folder/'seed-search/verified.json').exists():
            step(case, 'seed-replay')
        seed = read(folder/'seed-search/terminal.json')
        # Always run twice: first constructs, second independently replays saved inputs.
        step(case, 'seed-cloud'); step(case, 'seed-cloud')
        packet = read(folder/'seed-reconciled.json')
        data = []; used = 0; reason = seed['status']
        if seed['rank_lower_bound'] == 18 and packet['rank_lower_bound'] < 32:
            if not (folder/'complement-preparation/prepared.json').exists():
                step(case, 'complement-prepare')
            for index in range(15):
                if used >= 100:
                    break
                runfolder = folder/f'complement-{index:02d}'
                if not (runfolder/'terminal.json').exists():
                    step(case, 'search', index, 100-used)
                if not (runfolder/'verified.json').exists():
                    step(case, 'replay', index, 100-used)
                step(case, 'reconcile', index); step(case, 'reconcile', index)
                terminal = read(runfolder/'terminal.json'); packet = read(folder/f'reconciled-{index:02d}.json')
                data.append((terminal, packet)); used += terminal['charts']
                require(used <= 100, 'complement allowance exceeded')
                reason = terminal['stop_reason']
                if used == 100 or packet['rank_lower_bound'] >= 32:
                    break
                if reason != 'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION':
                    break
                require(packet['rank_lower_bound'] > terminal['rank_lower_bound'], 'reconciliation did not resolve additional rank')
                if not (folder/f'preparation-{index+1:02d}/prepared.json').exists():
                    step(case, 'rebuild', index)
            else:
                raise ArithmeticError('unexpected reconciliation depth')
        timeline, check_used = call_metrics(seed['charts'], seed['rank_lower_bound'], data,
                                          read(folder/'seed-reconciled.json')['reconciliation_gains'])
        require(used == check_used, 'call accounting differs')
        complement_gains = [g for g in timeline if g['phase'].startswith('complement')]
        last = max((g['call'] for g in complement_gains), default=None)
        rank = packet['rank_lower_bound']; total = seed['charts']+used
        added = rank-17
        result = {'status': 'PASS_INDEPENDENT_R17_60_CASE', **{k: row[k] for k in ('id','family','parameter','stratum')},
                  'seed_calls': seed['charts'], 'first_M18_found': seed['rank_lower_bound'] == 18,
                  'complement_calls': used, 'total_calls': total, 'rank_lower_bound': rank,
                  'added_directions': added, 'gaining_complement_calls': len({g['call'] for g in complement_gains}),
                  'calls_per_added_direction': str(Fraction(total, added)) if added else None,
                  'last_complement_gain_call': last, 'complement_calls_since_last_gain': used-last if last else used,
                  'gain_timeline': timeline, 'stop_reason': reason, 'packet': packet,
                  'point_timeouts': seed['point_timeouts']+sum(s['censored'] for t,_ in data for s in t['stages']),
                  'map_timeouts': sum(a['skip_reason'] in ('strict_wall_timeout','strict_rss_limit') for a in seed.get('map_attempts',[]))
                                  +sum(s['censored_maps'] for t,_ in data for s in t['stages']),
                  'protocol_sha256': sha(D/'protocol.json'),
                  'bindings': {str(p.relative_to(ROOT)): sha(p) for p in folder.rglob('*.json')
                               if p.name in ('terminal.json','verified.json') or p.name.startswith('reconciled-') or p.name=='seed-reconciled.json'}}
        atomic(folder/'result.json', result, immutable=True)
        print('R17_60_CASE_COMPLETE', case, rank, total, flush=True)
        return result
    finally:
        os.close(fd)


def report():
    results = [read(D/'cases'/r['id']/'result.json') for r in read(D/'roster.json')['rows']
               if (D/'cases'/r['id']/'result.json').exists()]
    def key(r):
        if r['status'] != 'PASS_INDEPENDENT_R17_60_CASE':
            return (1, 0, 0, Fraction(10**9), 0, r['id'])
        return (0, -r['rank_lower_bound'], -r['gaining_complement_calls'],
                Fraction(r['calls_per_added_direction']) if r['calls_per_added_direction'] else Fraction(10**9),
                -(r['last_complement_gain_call'] or 0), r['id'])
    good = [r for r in results if r['status']=='PASS_INDEPENDENT_R17_60_CASE']
    grouped = {}
    for name in ('low','high',*selection.FAMILIES):
        group = [r for r in good if r['stratum']==name or r['family']==name]
        hist = {}
        for r in group:
            rank=str(r['rank_lower_bound']);hist[rank]=hist.get(rank,0)+1
        grouped[name]={'verified_cases':len(group),'first_seeds':sum(r['first_M18_found'] for r in group),
                       'final_rank_histogram':hist,'point_calls':sum(r['total_calls'] for r in group),
                       'point_timeouts':sum(r['point_timeouts'] for r in group),
                       'map_timeouts':sum(r['map_timeouts'] for r in group)}
    output = {'status': 'COMPLETE_BOUNDED_R17_60_PANEL' if len(results)==60 else 'PARTIAL_R17_60_PANEL',
              'completed': len(results), 'independently_verified': len(good), 'results': results,
              'groups':grouped,
              'ranking': [r['id'] for r in sorted(results,key=key)],
              'late_gain_candidates': [r['id'] for r in good if (r['last_complement_gain_call'] or 0)>=75],
              'protocol_sha256': sha(D/'protocol.json')}
    atomic(D/'summary.json', output)
    if len(results)==60:
        atomic(ART/'r17_60_panel_results_v1.json', output, immutable=True)
    return output


def worker(case):
    try:
        case_run(case)
    except ResourceStop as exc:
        atomic(D/'cases'/case/'result.json', {'status':'RESOURCE_STOP_UNRESOLVED','id':case,'reason':str(exc)}, immutable=True)


def dispatch(case):
    folder = D/'cases'/case
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/'worker.log').open('ab',buffering=0) as log:
        proc = subprocess.Popen([sys.executable,'-u',str(SELF),'worker','--case',case],
                                cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
        code=proc.wait()
    require(code==0,f'case worker failed: {case}; see {folder}/worker.log')
    return case


def controller(fd):
    state={'status':'RUNNING','controller':process(os.getpid()),'started_unix':time.time(),'active':[]}
    atomic(D/'state.json',state)
    try:
        order=[r['id'] for r in read(D/'roster.json')['rows']]
        pending=[c for c in order if not (D/'cases'/c/'result.json').exists()]
        with ThreadPoolExecutor(max_workers=4) as pool:
            active={}
            while pending or active:
                while pending and len(active)<4:
                    check_protocol();case=pending.pop(0);active[pool.submit(dispatch,case)]=case
                state.update(active=list(active.values()),updated_unix=time.time())
                atomic(D/'state.json',state)
                done,_=wait(active,timeout=5,return_when=FIRST_COMPLETED)
                for future in done:
                    future.result();del active[future];report()
        summary=report();state.update(status=summary['status'],active=[],finished_unix=time.time())
        atomic(D/'state.json',state)
    except BaseException as exc:
        state.update(status='STOPPED_REVIEW_REQUIRED',error=str(exc));atomic(D/'state.json',state)
        raise
    finally:
        os.close(fd)


def launch():
    check_protocol()
    require(read(D/'preflight.json')['protocol_sha256']==sha(D/'protocol.json'),'preflight missing/stale')
    fd=lock(D/'controller.lock')
    try:
        # A retained case lock prevents duplicate writes even following an interrupted controller.
        for case in [r['id'] for r in read(D/'roster.json')['rows']]:
            check=lock(D/'cases'/case/'case.lock');os.close(check)
        with (D/'controller.log').open('ab',buffering=0) as log:
            proc=subprocess.Popen([sys.executable,'-u',str(SELF),'controller','--lock-fd',str(fd)],
                cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                start_new_session=True,pass_fds=(fd,))
        print('R17_60_DETACHED_LAUNCHED',proc.pid,'workers=4',flush=True)
    finally:
        os.close(fd)


def status():
    state=read(D/'state.json') if (D/'state.json').exists() else {'status':'NOT_LAUNCHED','active':[]}
    token=state.get('controller');state['controller_alive']=bool(token and process(token['pid'])==token)
    state['cases']={}
    for row in read(D/'roster.json')['rows']:
        folder=D/'cases'/row['id'];r={}
        if (folder/'result.json').exists():
            result=read(folder/'result.json');r={k:result[k] for k in ('status','rank_lower_bound','total_calls','last_complement_gain_call') if k in result}
        elif row['id'] in state['active'] and (folder/'progress.json').exists():
            r=read(folder/'progress.json')
            stage = folder/('seed-search' if r['phase'].startswith('seed-') else f"complement-{r['index']:02d}")
            if (stage/'progress.json').exists():r['search_progress']=read(stage/'progress.json')
        if r:state['cases'][row['id']]=r
    print(json.dumps(state,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('freeze','preflight','launch','resume','status','report','worker','controller'))
    p.add_argument('--case');p.add_argument('--lock-fd',type=int)
    a=p.parse_args()
    if a.action=='freeze':freeze()
    elif a.action=='preflight':preflight()
    elif a.action in ('launch','resume'):launch()
    elif a.action=='status':status()
    elif a.action=='report':print(json.dumps(report(),indent=2))
    elif a.action=='worker':worker(a.case)
    else:controller(a.lock_fd)


if __name__=='__main__':
    main()
