#!/usr/bin/env python3
"""Detached, metered comparison and gated follow-up; reuse the shared supervisor."""
import argparse
import ctypes
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

from next_direction_benchmark import read, sha
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits, _start_token


def launch(folder, hours):
    if not 0 < hours <= 8:
        raise ValueError('explicit ceiling must be in (0,8] worker-hours')
    if (folder/'launch.json').exists():
        raise FileExistsError('preserve prior run; no automatic restart')
    plan = read(folder/'plan.json')
    if sha(shutil.which('sage')) != plan['sage_launcher_sha256']:
        raise ArithmeticError('Sage launcher changed')
    for relative, value in plan['sources'].items():
        if sha(folder/relative) != value:
            raise ArithmeticError('source snapshot changed: '+relative)
    command = [sys.executable,str(folder/'runtime/elliptic-curves/cas'/Path(__file__).name),
        'controller','--folder',str(folder),'--hours',str(hours)]
    with (folder/'controller.log').open('w') as log:
        proc = subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    checkpoint(folder/'launch.json',{'pid':proc.pid,'token':_start_token(proc.pid),
        'hours':hours,'plan_sha256':sha(folder/'plan.json'),'command':command,'started_unix':time.time()})
    print('DETACHED',proc.pid,folder,flush=True)


def controller(folder, hours):
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        _controller(folder,hours)


def _controller(folder, hours):
    plan = read(folder/'plan.json')
    if ctypes.CDLL(None).prctl(36,1,0,0,0) != 0:
        raise RuntimeError('Linux child subreaper unavailable')
    cpu = min(os.sched_getaffinity(0)); os.sched_setaffinity(0,{cpu})
    environment = {k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']}
    os.environ.update(environment)
    started = time.monotonic()
    prep = plan['preparation_wall_seconds']
    state = {'status':'RUNNING','plan_sha256':sha(folder/'plan.json'),'cpu_affinity':cpu,
        'maximum_workers':1,'ceiling_seconds':hours*3600,'stages':[],
        'worker_tree_cpu_seconds':0.0,'preparation_wall_seconds':prep,
        'preparation_cpu_seconds':plan['preparation_cpu_seconds']}
    def remaining():
        return hours*3600-prep-(time.monotonic()-started)
    def save():
        state['charged_wall_seconds'] = prep+time.monotonic()-started
        checkpoint(folder/'state.json',state)
    def stage(arm, mode, seconds, policy):
        for relative, value in plan['sources'].items():
            if sha(folder/relative) != value:
                raise ArithmeticError('changed frozen implementation')
        state['active_arm'],state['active_stage'] = str(arm.relative_to(folder)),mode
        save()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        own, start = time.process_time(),time.monotonic()
        command = [shutil.which('sage'),'-python',str(folder/'runtime/elliptic-curves/cas/next_direction_benchmark.py'),
            mode,'--folder',str(arm),'--policy',policy,'--seconds',str(max(30,seconds-5))]
        result = run(command,limits=Limits(seconds,plan['rss_bytes']),log_path=arm/(mode+'.log'),
            checkpoint_path=arm/(mode+'-supervisor.json'),env=environment)
        while True:
            try: os.waitpid(-1,0)
            except ChildProcessError: break
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        used = after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime+time.process_time()-own
        row = {'arm':str(arm.relative_to(folder)),'mode':mode,'wall_seconds':time.monotonic()-start,
            'process_tree_cpu_seconds':used,'outcome':result['outcome'],'returncode':result['returncode']}
        state['stages'].append(row);state['worker_tree_cpu_seconds'] += used;save()
        return result['outcome']=='completed' and result['returncode']==0,row
    def arm(case, policy, seconds):
        dest = folder/'runs'/case['role']/case['id']/policy
        dest.mkdir(parents=True,exist_ok=False)
        if sha(folder/case['input']) != case['sha256']:
            raise ArithmeticError('input packet changed')
        shutil.copyfile(folder/case['input'],dest/'input.json')
        checkpoint(dest/'protocol.json',{k:plan[k] for k in ['policies','map_python','seconds_per_map','map_rss_bytes','gp_sha256']})
        good,first = stage(dest,'worker',seconds,policy)
        verified,second = False,{'process_tree_cpu_seconds':0.0,'wall_seconds':0.0}
        if good and remaining() > plan['replay_seconds']:
            verified,second = stage(dest,'replay',plan['replay_seconds'],policy)
        success = verified and read(dest/'verified.json')['independent_gain_verified']
        row = {'case':case['id'],'role':case['role'],'policy':policy,'success':success,
            'verified':verified,'cpu_seconds':first['process_tree_cpu_seconds']+second['process_tree_cpu_seconds'],
            'wall_seconds':first['wall_seconds']+second['wall_seconds'],'allowance_seconds':seconds,
            'penalized_cpu_seconds':first['process_tree_cpu_seconds']+second['process_tree_cpu_seconds'] if success else 2*(seconds+plan['replay_seconds'])}
        checkpoint(dest/'seal.json',row)
        # Infrastructure failures are a stop, not permission for a repair campaign.
        if not good and first['outcome'] not in ('strict_wall_timeout','strict_rss_limit'):
            raise RuntimeError('worker failed; preserve checkpoint and stop calibration')
        if good and not verified:
            raise RuntimeError('independent replay incomplete; stop calibration')
        return row

    rows = []
    policies = list(plan['policies'])
    development = [c for c in plan['cases'] if c['role']=='development']
    try:
        for i,case in enumerate(development):
            # Rotate lane order to distribute order effects without changing any lane's centres.
            for policy in policies[i:]+policies[:i]:
                if remaining() < plan['benchmark_arm_seconds']+plan['replay_seconds']+5:
                    state['status']='STOPPED_CEILING_DURING_CALIBRATION';save();return
                rows.append(arm(case,policy,plan['benchmark_arm_seconds']))
        scores = [{'policy':p,'successes':sum(r['success'] for r in rows if r['policy']==p),
            'penalized_cpu_seconds':sum(r['penalized_cpu_seconds'] for r in rows if r['policy']==p)} for p in policies]
        scores.sort(key=lambda s:(-s['successes'],s['penalized_cpu_seconds'],s['policy']))
        winner = scores[0]['policy']
        checkpoint(folder/'development-seal.json',{'rows':rows,'scores':scores,'winner':winner})
        if scores[0]['successes'] < 2:
            state['status']='STOPPED_NO_CALIBRATED_NEXT_DIRECTION_POLICY';save();return
        baseline = 'v3_dual_125k'
        validation_rows = []
        for case in [c for c in plan['cases'] if c['role']=='validation']:
            validated = []
            for policy in dict.fromkeys([baseline,winner]):
                if remaining() < plan['validation_arm_seconds']+plan['replay_seconds']+5:
                    state['status']='STOPPED_CEILING_BEFORE_VALIDATION';save();return
                validated.append(arm(case,policy,plan['validation_arm_seconds']))
            validation_rows.extend(validated)
            checkpoint(folder/'validation-seal.json',{'rows':validation_rows,'winner':winner})
            win = next(r for r in validated if r['policy']==winner)
            base = next(r for r in validated if r['policy']==baseline)
            if not win['success'] or (winner!=baseline and base['success'] and win['cpu_seconds'] > base['cpu_seconds']):
                state['status']='STOPPED_VALIDATION_GATE';save();return
        checkpoint(folder/'winner.json',{'policy':winner,'plan_sha256':sha(folder/'plan.json'),
            'development_seal_sha256':sha(folder/'development-seal.json'),
            'validation_seal_sha256':sha(folder/'validation-seal.json'),
            'scope':'Frozen fixed-bank comparison winner; conditional measured evidence, not a universal speed claim.'})
        production = [c for c in plan['cases'] if c['role']=='production']
        cap = min(7200,(remaining()-len(production)*(plan['replay_seconds']+5))/len(production))
        if cap < 300:
            state['status']='STOPPED_NO_PRODUCTION_ALLOWANCE';save();return
        state['production_arm_seconds']=cap;save()
        for case in production:
            rows.append(arm(case,winner,cap))
        state['status']='COMPLETED_FROZEN_NEXT_DIRECTION_PILOT'
        checkpoint(folder/'summary.json',{'rows':rows,'winner':winner,
            'fresh_next_direction_successes':sum(r['success'] for r in rows if r['role']=='production')})
        save()
    except BaseException as error:
        state['status']='STOPPED_FAILURE';state['error']=f'{type(error).__name__}: {error}';save()
        raise


def status(folder):
    launch = read(folder/'launch.json') if (folder/'launch.json').exists() else None
    alive = bool(launch and launch['token'] is not None and _start_token(launch['pid']) == launch['token'])
    state = read(folder/'state.json') if (folder/'state.json').exists() else {'status':'FROZEN_NOT_STARTED'}
    print(json.dumps({'status':state['status'],'controller_live':alive,
        **{k:state[k] for k in ['active_arm','active_stage','charged_wall_seconds','worker_tree_cpu_seconds','error'] if k in state}},indent=2))
    if 'active_arm' in state:
        p=folder/state['active_arm']/'progress.json'
        if p.exists(): print(json.dumps(read(p),indent=2))
        print('Log:',folder/state['active_arm']/(state['active_stage']+'.log'))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['launch','controller','status','stop'])
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--hours',type=float)
    a=p.parse_args();folder=a.folder.resolve()
    if a.mode in ('launch','controller'):
        if a.hours is None:p.error('an explicit --hours ceiling is required')
        if not 0<a.hours<=8:p.error('--hours must be in (0,8]')
        globals()[a.mode](folder,a.hours)
    elif a.mode=='status':status(folder)
    else:
        import signal
        lease=read(folder/'launch.json')
        if lease['token'] is not None and _start_token(lease['pid'])==lease['token']:
            os.kill(lease['pid'],signal.SIGTERM);print('Stop sent; owned workers are supervised; preserve all checkpoints.')
        else:print('Controller is not live; nothing signalled.')
