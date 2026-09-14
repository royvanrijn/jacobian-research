"""Freeze, detach, meter and report the bounded three-arm height experiment."""
import argparse
import ctypes
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits, _start_token
from height_model_benchmark import ARMS, read, sha, require

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]


def freeze(folder):
    start, cpu = time.monotonic(), time.process_time()
    folder.mkdir(parents=True, exist_ok=False)
    old = ROOT/'artifacts/local/elliptic-curves/next-direction-benchmark-v1'
    plan = read(old/'plan.json')
    cases = []
    for rank in (27,29,30):
        case = next(c for c in plan['cases'] if c['id'] == f'curve302-M{rank}')
        require(sha(old/case['input']) == case['sha256'], 'changed retained input')
        target = folder/case['input']; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(old/case['input'], target)
        cases.append({k:case[k] for k in ('id','input','sha256','initial_rank','centres')})
    example = read(ROOT/'artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1/preconditioned_full-input.json')
    proof = example['prime_certificate']; require(sha(ROOT/proof['path']) == proof['sha256'], 'prime certificate changed')
    target = folder/'runtime'/proof['path']; target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ROOT/proof['path'], target)
    sources = {}
    for base in (CAS, CAS.parent/'ecsearch'):
        for path in sorted(base.rglob('*')):
            if not path.is_file() or path.suffix not in ('.py','.sage','.cpp','.h') or '__pycache__' in path.parts:
                continue
            target = folder/'runtime'/path.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, target)
            sources[str(target.relative_to(folder))] = sha(target)
    checkpoint(folder/'source-lock.json', sources)
    result = {'schema':'height-model-next-direction.v1', 'status':'FROZEN_NOT_STARTED',
        'cases':cases, 'arms':list(ARMS), 'arm_seconds':300, 'worker_seconds':240,
        'search_phase_seconds':210, 'support_seconds':30, 'diagnostic_seconds':300, 'aggregate_seconds':3060,
        'seconds_per_map':5, 'seconds_per_search':10, 'map_rss_bytes':1024**3,
        'rss_bytes':3*1024**3, 'maximum_workers':1, 'threads_per_worker':1,
        'map_python':shutil.which('sage'), 'sage_sha256':sha(shutil.which('sage')),
        'gp_sha256':sha('/usr/bin/gp'), 'source_lock_sha256':sha(folder/'source-lock.json'),
        'prime_certificate':proof, 'certified_primes':example['certified_primes'],
        'input_source_plan_sha256':sha(old/'plan.json'),
        'selection':'Original first256 centre order on Curve302 M27,M29,M30; no target inputs. '
            'Current cheaper factor-free V3 baseline at125000. Bounded arms construct full minimal '
            'model and first two original local neighbours, primes3..167 in fixed order, without enlargement. '
            'Only the previously verified distinct-prime portfolio formula is admitted; unsupported charts '
            'are explicitly incomplete. Uniform C/L picks the single. Joint max/min picks exactly two. '
            'Single H125000; pair least integer H matching that single guaranteed coverage.',
        'measurement':'One fresh worker per arm, no cross-arm map cache. Time starts before process import/input '
            'and includes starting rank, map construction/reduction, bounds and independent bound replay, '
            'search, exact map transport, duplicates, rank certification and fresh-process Sage finite-group '
            'replay. Linux RUSAGE_CHILDREN with subreaper meters complete descendant CPU plus controller CPU. '
            'Elapsed separately. Initial frozen Q-bank generation is outside supplied (E,M,Q) interface and '
            'historical cold landscape cost remains UNKNOWN. Bounded arms rebuild discriminant support '
            'from E with fresh factorization and unconditional prime proofs, no hints, in at most30 seconds. '
            'The retained prime certificate is available only to post-seal diagnostics; it is not used by search selection.',
        'stops':'First independently replayed quotient gain per arm; no retry or centre/box enlargement. '
            'Five minute total per arm, nine arms only. Stop on infrastructure/replay failure. '
            'No class construction, ancestry, production curves or rank32 attempt.',
        'diagnostics':'Only after all nine endpoints are sealed. Locate each literal winning point at its '
            'winning anchor in factor-free plus three minimal models; exact height, H^4/(D*Hx), '
            'guaranteed-radius membership and requested-box membership. Diagnostics are separately metered '
            'and cannot select models or feed point search. Different winners/anchors are not one common target.',
        'order':'Rotate arm order over the three cases; one repetition, no statistical speed claim.',
        'preparation_wall_seconds':time.monotonic()-start, 'preparation_cpu_seconds':time.process_time()-cpu}
    checkpoint(folder/'protocol.json', result)
    print('FROZEN', folder, 'nine five-minute arms; no point search yet', flush=True)


def validate(folder):
    protocol = read(folder/'protocol.json')
    require(sha(folder/'source-lock.json') == protocol['source_lock_sha256'], 'source lock changed')
    for relative, value in read(folder/'source-lock.json').items():
        require(sha(folder/relative) == value, 'frozen source changed: '+relative)
    require(sha(protocol['map_python']) == protocol['sage_sha256'], 'Sage launcher changed')
    require(sha('/usr/bin/gp') == protocol['gp_sha256'], 'GP binary changed')
    return protocol


def launch(folder):
    require(not (folder/'launch.json').exists(), 'preserve previous launch; no automatic retry')
    validate(folder)
    command = [sys.executable, str(folder/'runtime/elliptic-curves/cas'/Path(__file__).name),
               'controller','--folder',str(folder)]
    with (folder/'controller.log').open('w') as log:
        child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                                 stderr=subprocess.STDOUT, start_new_session=True)
    checkpoint(folder/'launch.json', {'pid':child.pid, 'token':_start_token(child.pid),
        'protocol_sha256':sha(folder/'protocol.json'), 'command':command, 'started_unix':time.time()})
    print('DETACHED', child.pid, folder, flush=True)


def controller(folder):
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        _controller(folder)


def _controller(folder):
    protocol=validate(folder)
    require(ctypes.CDLL(None).prctl(36,1,0,0,0)==0, 'subreaper unavailable')
    core=min(os.sched_getaffinity(0)); os.sched_setaffinity(0,{core})
    env={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')}
    os.environ.update(env)
    start=time.monotonic(); stages=[]; rows=[]
    state={'status':'RUNNING','cpu_affinity':core,'maximum_workers':1,
           'protocol_sha256':sha(folder/'protocol.json'),'rows':rows,'stages':stages}
    def save():
        state['elapsed_seconds']=time.monotonic()-start
        state['process_tree_cpu_seconds']=sum(s['process_tree_cpu_seconds'] for s in stages)
        checkpoint(folder/'state.json',state)
    def stage(dest, mode, seconds):
        state.update(active_arm=str(dest.relative_to(folder)),active_stage=mode);save()
        before=resource.getrusage(resource.RUSAGE_CHILDREN); own=time.process_time(); wall=time.monotonic()
        command=[protocol['map_python'],'-python',str(CAS/'height_model_benchmark.py'),mode,'--folder',str(dest)]
        outcome=run(command,limits=Limits(max(.1,seconds),protocol['rss_bytes']),
                    log_path=dest/(mode+'.log'),checkpoint_path=dest/(mode+'-supervisor.json'),env=env)
        while True:
            try: os.waitpid(-1,0)
            except ChildProcessError:break
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        cost=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime+time.process_time()-own
        row={'folder':str(dest.relative_to(folder)), 'mode':mode, 'wall_seconds':time.monotonic()-wall,
             'process_tree_cpu_seconds':cost,'outcome':outcome['outcome'],'returncode':outcome['returncode']}
        stages.append(row);save();return row
    try:
        for ci,case in enumerate(protocol['cases']):
            arms=list(ARMS[ci:]+ARMS[:ci])
            for arm in arms:
                if time.monotonic()-start+protocol['arm_seconds']>protocol['aggregate_seconds']:
                    state['status']='STOPPED_AGGREGATE_CEILING';save();return
                dest=folder/'runs'/case['id']/arm;dest.mkdir(parents=True,exist_ok=False)
                require(sha(folder/case['input'])==case['sha256'],'changed starting input')
                shutil.copyfile(folder/case['input'],dest/'input.json')
                checkpoint(dest/'protocol.json',{**protocol,'arm':arm})
                wall=time.monotonic(); runrow=stage(dest,'worker',protocol['worker_seconds'])
                verified=False; reprow=None
                normal=runrow['outcome']=='completed' and runrow['returncode']==0
                if normal and (remaining:=protocol['arm_seconds']-(time.monotonic()-wall))>1:
                    reprow=stage(dest,'replay',remaining)
                    verified=reprow['outcome']=='completed' and reprow['returncode']==0
                success=verified and read(dest/'verified.json')['independent_gain_verified']
                row={'case':case['id'],'arm':arm,'folder':str(dest.relative_to(folder)),
                    'success':success,'verified':verified,'wall_seconds':time.monotonic()-wall,
                    'cpu_seconds':runrow['process_tree_cpu_seconds']+(reprow['process_tree_cpu_seconds'] if reprow else 0),
                    'status':read(dest/'result.json')['status'] if (dest/'result.json').exists() else runrow['outcome'],
                    'allowance_seconds':protocol['arm_seconds']}
                checkpoint(dest/'seal.json',row);rows.append(row);save()
                print('SEALED',case['id'],arm,success,round(row['cpu_seconds'],3),flush=True)
                if (not normal and runrow['outcome'] not in ('strict_wall_timeout','strict_rss_limit')) or (normal and not verified):
                    state['status']='STOPPED_WORKER_OR_REPLAY_FAILURE';save();return
        checkpoint(folder/'all-arms-sealed.json', {'protocol_sha256':sha(folder/'protocol.json'),'rows':rows,
                   'seal_hashes':{r['folder']:sha(folder/r['folder']/'seal.json') for r in rows}})
        diagnostic=stage(folder,'diagnose',protocol['diagnostic_seconds'])
        state['status']='COMPLETE' if diagnostic['outcome']=='completed' and diagnostic['returncode']==0 else 'BENCHMARK_COMPLETE_DIAGNOSTICS_INCOMPLETE'
        checkpoint(folder/'summary.json', {'rows':rows,'diagnostic':diagnostic,
                   'boundary':'One rotated-order run per retained Curve302 transition. No cross-curve transfer or new rank record.'})
        save()
    except BaseException as error:
        state.update(status='STOPPED_EXCEPTION',error=f'{type(error).__name__}: {error}');save();raise


def status(folder):
    lease=read(folder/'launch.json') if (folder/'launch.json').exists() else None
    live=bool(lease and lease['token'] is not None and _start_token(lease['pid'])==lease['token'])
    state=read(folder/'state.json') if (folder/'state.json').exists() else {'status':'FROZEN_NOT_STARTED'}
    completion=read(folder/'completion.json') if (folder/'completion.json').exists() else None
    display={k:v for k,v in state.items() if k not in ('rows','stages')}
    if completion:
        display['controller_terminal_status']=display.get('status')
        display['status']=completion['status']
        display['diagnostic_cpu_seconds']=completion['diagnostic_cpu_seconds']
    print(json.dumps({'controller_live':live,**display,
                      'completed_arms':state.get('rows',[])},indent=2))
    if live and state.get('active_arm'):
        path=folder/state['active_arm']/'progress.json'
        if path.exists():
            data=read(path);print('Current progress:',json.dumps({k:v for k,v in data.items() if k!='rows'}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('freeze','launch','controller','status','stop'))
    p.add_argument('--folder',type=Path,required=True)
    a=p.parse_args();folder=a.folder.resolve()
    if a.mode=='stop':
        lease=read(folder/'launch.json')
        if lease['token'] is not None and _start_token(lease['pid'])==lease['token']:
            os.kill(lease['pid'],signal.SIGTERM);print('Stop sent to owned controller.')
        else:print('Controller not live; nothing signalled.')
    else:globals()[a.mode](folder)
