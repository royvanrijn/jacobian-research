"""Freeze, launch detached, follow and stop the user-authorized staged pilot.

One local arithmetic worker, one thread, sequential paired routes. There are
no retries, roster replacements, adaptive caps, BNF calls or agent launches.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import signal
import subprocess
import sys
import time

from transfer_common import *

PANEL=ROOT/'artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1'
DEFAULT=ROOT/'artifacts/local/elliptic-curves/fresh-constructor-transfer-v3'
PROGRAMMES=['transfer_common.py','transfer_constructor.py','transfer_lift.py','transfer_v3.py','transfer_compare.py','run_fresh_constructor_transfer.py']


def policy():
    return dict(schema='mw16-fresh-constructor-transfer-execution.v1',
        authorization='2026-09-12: metered commissioning, first two paired fresh fibres, conditional remaining six; maximum 36 worker-hours. Overnight unattended requested, not a minimum spend.',
        reference_cap_seconds=7200,aggregate_cap_seconds=129600,local_workers=1,local_threads=1,
        rss_bytes=2*1024**3,pari_stack_bytes=268435456,
        resource_measure='Charge elapsed route wall time including preparation, verification, failed attempts and remote requests; record local process-tree CPU and remote reported CPU separately. One fixed local worker. Remote hardware/threads unavailable; speed comparisons remain limited.',
        reference_order=['class','v3'],fresh_order=['class','v3','post_seal_union'],
        fresh_budget_formula='B = min(2*C, 7200), C includes full reference class cost and post-seal union verification',
        first_gate='No fresh fibre unless reference class-to-point and rank certificates independently pass within 7200 charged seconds.',
        second_gate='After first two fixed fresh fibres, continue remaining six only if at least one fresh class-to-point rank gain independently passes.',
        no_automatic_retry=True,lifting_class_limit=2,
        class_stage_schedule=dict(prepare='remaining construction allowance',construct='remaining allowance with min(1200, 25% of route cap) reserved for verification and lifting',verify_classes='remaining route allowance',compact='remaining route allowance',lifting='equal shares of the remaining lifting allowance across first two retained classes; no replacement',rank_verification_reserve='min(120 seconds, 5% of cap)',post_seal_union_reserve='min(60 seconds, 5% of cap)'),
        constructor=dict(smooth_bound=50000,anchor_prime_bound=1009,initial_box=512,
            strip_levels=[1024,2048,4096,8192],targets_per_level=512,monitor_interval=32,
            proof_prime_min=50001,proof_prime_max=55000,character_coordinates=96,
            method='Existing small representatives of unresolved ideal combinations, pivot columns eligible, generic-only protected unramified anchors, full prime blocks and support-disjoint characters. Cold empty atom pool; no formal quotient target.'),
        lifting=dict(quartic_bounds=[1000,100000],maximum_cover_coefficient_bits=20000,
            reduction='Degree 4 Minimise/Reduce, rational rank-three pencil/conic chart, degree 2 Minimise/Reduce; exact transport replay.',
            remote='Standard public Magma calculator: 60 CPU seconds and <=2 GiB per request, <=50000 input bytes, 65-second network timeout, no retries or service-limit workarounds.'),
        v3=dict(generic_rank=16,scaled_shells=[16,19,20,23],anchors_per_shell=16,canonical_per_shell=25,
            exact_cvp_node_limit=2000000,generic_bank_node_limit=20000000,generic_ellipsoid_bound=23,
            prime_bound=1000,height=125000,seconds_per_chart=10,max_charts=10000,seconds_per_map=5,
            map_rss_bytes=1024**3,map_order=['preconditioned_full','factor_free'],max_epochs=17,
            landscape_seconds=60,selection='Maintained upper-shell V3, immediate rebuild after certified gain, stop at first no-gain epoch, finite bank or metered endpoint.'),
        endpoints=['arithmetic preparation','independently certified strict class outside generic span','exact rational lift','independent rank certificate','post-seal union gain'],
        ordinary_ideal_class_endpoint='UNKNOWN unless separately certified; no Artin/half-ideal class search in this first pilot.',
        forbidden=['historical exceptional points as inputs','cached principal relations','old atom words or class columns','carrier work','fixed-word continuation','full class-group computation','roster replacements','automatic budget increases'])


def freeze(folder):
    if folder.exists():raise FileExistsError('preserve existing run; no automatic retry')
    selection=read(PANEL/'selection.json')
    # Verify the already frozen selection receipt and every sanitized input.
    expected={'family.json':'ae38533c564930af0c1cfdeeca6d72f5953a8bef759b8c00b1b8f38cb8fd4e22',
        'equation-aliases.json':'fcfe7df9146a8ad6b472fc3637497accb04612a1ed39fec1063ebee0748aa949',
        'selection.json':'aeb446e4478c2e75cdb53ca26894c9d51f08e79580bdb9ab7e36053bb4cee88d'}
    for name,h in expected.items():assert sha(PANEL/name)==h,name
    for name,h in selection['bindings'].items():assert sha(PANEL/name)==h,name
    sources=[PANEL/'commissioning-input.json']+[PANEL/'inputs'/f'fresh-{i:02d}.json' for i in range(1,9)]
    roster=[read(p) for p in sources]
    assert [r['parameter'] for r in roster]==['3/17','27/4','-29/4','18/17','-4/19','32/3','-17/14','-7/17','4/11']
    folder.mkdir(parents=True)
    implementation=folder/'implementation/research/elliptic-curves'
    source_manifest={}
    # Code only: snapshots survive unrelated concurrent repository cleanup.
    for area in ['cas','rank-jump','ecsearch']:
        base=ROOT/'elliptic-curves'/area
        for source in sorted(base.rglob('*')):
            if not source.is_file() or source.suffix not in {'.py','.sage','.cpp','.h'} or '__pycache__' in source.parts:continue
            relative=source.relative_to(ROOT/'elliptic-curves')
            destination=implementation/relative;destination.parent.mkdir(parents=True,exist_ok=True)
            data=source.read_bytes();destination.write_bytes(data)
            assert source.read_bytes()==data,'source changed while freezing'
            source_manifest[str(relative)]=sha(destination)
    for source,row in zip(sources,roster):
        write(folder/'inputs'/(row['id']+'.json'),row)
    p=policy();p.update(frozen_unix=time.time(),selection_sha256=expected['selection.json'],
        roster=[dict(id=r['id'],parameter=r['parameter'],input_sha256=sha(folder/'inputs'/(r['id']+'.json'))) for r in roster],
        source_manifest_sha256=None,
        hardware=dict(hostname=platform.node(),platform=platform.platform(),python=sys.version,available_cpus=sorted(os.sched_getaffinity(0))),
        cpu_affinity=min(os.sched_getaffinity(0)),sage=shutil.which('sage'),gp_sha256=sha('/usr/bin/gp'),
        data_boundary='Application-level read audit: per-arm equation/generic packet and own generated data only; post-seal evaluator gets sealed outputs. Code is snapshotted. Not OS isolation.')
    previous=ROOT/'artifacts/local/elliptic-curves/fresh-constructor-transfer-v1'
    p['prior_attempt_charges']={'class':0,'v3':0}
    if previous.exists() and previous!=folder:
        prior=read(previous/'state.json');assert prior['status']=='STOPPED_COMMISSIONING_GATE'
        p['prior_attempt_charges']={'class':prior['cases']['control']['class']['charged_seconds']+prior['cases']['control']['comparison']['charged_seconds'],
                                   'v3':prior['cases']['control']['v3']['charged_seconds']}
        p['prior_attempt']=dict(path=str(previous),state_sha256=sha(previous/'state.json'),reason='Snapshot omitted ecsearch import package; all routes stopped before field or point arithmetic. Packaging repaired, no change to mathematical policy; costs retained within original cap.')
    previous2=ROOT/'artifacts/local/elliptic-curves/fresh-constructor-transfer-v2'
    if previous2.exists() and previous2!=folder:
        prior=read(previous2/'state.json');assert prior['status']=='STOPPED_COMMISSIONING_GATE'
        assert read(previous2/'control/class/construction-terminal.json')['status']=='NO_NEW_CLASS_WITHIN_FROZEN_BANK'
        p['baseline_only_after_closed_class_gate']=True
        p['prior_attempt_charges']['class']+=prior['cases']['control']['class']['charged_seconds']+prior['cases']['control']['comparison']['charged_seconds']
        p['prior_attempt_charges']['v3']+=prior['cases']['control']['v3']['charged_seconds']
        p['native_packaging_repair']=dict(path=str(previous2),state_sha256=sha(previous2/'state.json'),reason='V3 native C++ helper omitted from code snapshot. Add native source and code-read permission; finish only the independent V3 reference within remaining original cap. The class bank is closed and is not rerun; no fresh fibres permitted.')
        shutil.copytree(previous2/'control/class',folder/'control/class')
    write(folder/'source-manifest.json',source_manifest);p['source_manifest_sha256']=sha(folder/'source-manifest.json')
    write(folder/'execution-protocol.json',p)
    # A compact public binding; preserve all historical design snapshots.
    public=PANEL/'execution-protocol-v3.json'
    if public.exists():raise FileExistsError(public)
    write(public,dict(**p,local_run=str(folder),protocol_sha256=sha(folder/'execution-protocol.json')))
    write(folder/'state.json',dict(status='FROZEN_NOT_LAUNCHED',charged_seconds=0,cases={},protocol_sha256=sha(folder/'execution-protocol.json')))
    return p


def token(pid):
    try:return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):return None


def launch(folder):
    freeze(folder)
    script=folder/'implementation/research/elliptic-curves/rank-jump/run_fresh_constructor_transfer.py'
    with (folder/'controller.log').open('x') as log:
        proc=subprocess.Popen([sys.executable,str(script),'run','--folder',str(folder)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
    write(folder/'launch.json',dict(pid=proc.pid,start_token=token(proc.pid),started_unix=time.time(),command=[sys.executable,str(script),'run','--folder',str(folder)],detached=True))
    print(json.dumps(dict(pid=proc.pid,folder=str(folder),status_command=f'python3 {Path(__file__).resolve()} status --folder {folder}'),indent=2))


class Controller:
    def __init__(self,folder):
        self.folder=folder;self.p=read(folder/'execution-protocol.json');self.state=read(folder/'state.json')
        self.code=folder/'implementation/research/elliptic-curves/rank-jump'
        self.manifest=read(folder/'source-manifest.json')
        self.case_cost=0;self.total=sum(self.p.get('prior_attempt_charges',{}).values())
        os.sched_setaffinity(0,{self.p['cpu_affinity']})
        self.env={name:'1' for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS']}
        self.env['PYTHONUNBUFFERED']='1';self.env['PYTHONDONTWRITEBYTECODE']='1'
        self.env['SAGE_NUM_THREADS']='1'
        self.started=time.monotonic()

    def update(self,**kw):
        self.state.update(kw,charged_seconds=self.total,updated_unix=time.time(),controller_elapsed_seconds=time.monotonic()-self.started)
        write(self.folder/'state.json',self.state)

    def stage(self,arm,script,mode,maximum,number=None):
        from research_runtime.supervisor import run,Limits
        if (self.folder/'STOP').exists():raise KeyboardInterrupt('user stop request')
        maximum=min(maximum,self.p['aggregate_cap_seconds']-self.total,
                    self.current_deadline-time.monotonic())
        if maximum<=1:return False
        for name,h in self.manifest.items():assert sha(self.code.parent/name)==h,name
        name=mode+(f'-{number}' if number is not None else '')
        command=[self.p['sage'],'-python',str(self.code/script)]
        if mode!='compare':command.append(mode)
        command+=['--folder',str(arm)]
        if number is not None:command+=['--number',str(number)]
        self.update(status='RUNNING',active_case=arm.parent.name if arm.name in ['class','v3','compare'] else arm.name,active_route=arm.name,active_stage=name,active_arm=str(arm),stage_started_unix=time.time(),stage_allowance_seconds=maximum)
        before=resource.getrusage(resource.RUSAGE_CHILDREN);own=time.process_time()
        record=run(command,limits=Limits(wall_seconds=max(0.1,maximum-1.5),rss_bytes=self.p['rss_bytes'],pari_stack_bytes=self.p['pari_stack_bytes'],poll_seconds=.5),log_path=arm/(name+'.log'),checkpoint_path=arm/(name+'-supervisor.json'),cwd=arm,env=self.env)
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        record.update(stage=name,route=arm.name,case=arm.parent.name,process_tree_cpu_seconds=(after.ru_utime+after.ru_stime)-(before.ru_utime+before.ru_stime),supervisor_cpu_seconds=time.process_time()-own,finished_unix=time.time())
        # The wait4 child totals include waited descendants of this sequential
        # worker. Detached escaped descendants are forbidden by the supervisor.
        self.case_cost+=record['wall_seconds'];self.total+=record['wall_seconds']
        with (self.folder/'meter.jsonl').open('a') as h:h.write(json.dumps(clean(record))+'\n')
        write(arm/(name+'-meter.json'),record)
        self.update(last_outcome=record['outcome'])
        return record['outcome']=='completed'

    def arm(self,case,route,cap):
        arm_started=time.monotonic()
        self.current_deadline=arm_started+cap-5
        arm=self.folder/case/route;arm.mkdir(parents=True)
        shutil.copyfile(self.folder/'inputs'/(case+'.json'),arm/'input.json')
        write(arm/'policy.json',self.p)
        self.case_cost=0;stages=[];reserve=min(120,.05*cap)
        success=False
        if route=='class':
            construction_cap=cap-min(1200,.25*cap)
            for mode in ['prepare','construct']:
                ok=self.stage(arm,'transfer_constructor.py',mode,construction_cap-self.case_cost)
                stages.append(dict(stage=mode,completed=ok))
                if not ok:break
            if all(s['completed'] for s in stages) and (arm/'classes.json').exists():
                for mode in ['verify_classes','compact']:
                    ok=self.stage(arm,'transfer_constructor.py',mode,cap-self.case_cost-reserve)
                    stages.append(dict(stage=mode,completed=ok))
                    if not ok:break
                if all(s['completed'] for s in stages):
                    count=min(2,len(read(arm/'classes.json')['classes']))
                    for i in range(count):
                        # Freeze the remaining equal share before this cover.
                        allowance=max(0,(cap-self.case_cost-reserve)/(count-i));start=self.case_cost
                        for mode in ['verify_compact','reduce_cover','search_cover','verify_lift']:
                            ok=self.stage(arm,'transfer_lift.py',mode,allowance-(self.case_cost-start),i)
                            stages.append(dict(stage=mode,number=i,completed=ok))
                            if not ok:break
                    self.stage(arm,'transfer_lift.py','verify_rank',cap-self.case_cost)
        else:
            p=dict(self.p['v3'],map_python=self.p['sage'],gp_sha256=self.p['gp_sha256'],search_deadline_unix=time.time()+cap-reserve-30)
            write(arm/'protocol.json',p)
            if self.stage(arm,'transfer_v3.py','prepare',cap-reserve):
                if self.stage(arm,'transfer_v3.py','search',cap-self.case_cost-reserve):
                    self.stage(arm,'transfer_v3.py','verify',cap-self.case_cost)
        verified=read(arm/'verified.json') if (arm/'verified.json').exists() else None
        if route=='class':success=bool(verified and verified.get('transfer_success'))
        hashes={p.name:sha(p) for p in arm.glob('*.json') if p.name!='seal.json'}
        complete_elapsed=time.monotonic()-arm_started
        overhead=max(0,complete_elapsed-self.case_cost)
        self.case_cost+=overhead;self.total+=overhead
        if self.case_cost>cap:success=False
        result=dict(case=case,route=route,status=verified['status'] if verified else 'INCOMPLETE_WITH_RETAINED_CHECKPOINTS',
            transfer_success=success,verified=bool(verified),charged_seconds=self.case_cost,cap_seconds=cap,files=hashes,
            orchestration_overhead_seconds=overhead,
            stage_evidence=dict(preparation='PASS' if (arm/('field.json' if route=='class' else 'prepared.json')).exists() else 'UNKNOWN',
                initial_rank='CERTIFIED_16' if (arm/'seed.json').exists() else 'UNKNOWN',
                class_construction='CERTIFIED' if (arm/'classes-verified.json').exists() else ('CANDIDATES_UNVERIFIED' if (arm/'classes.json').exists() else 'NO_CERTIFIED_NEW_CLASS'),
                rational_lifting='CERTIFIED' if success else 'UNKNOWN',rank_gain=(verified.get('rank_lower_bound',16)-16) if verified else 'UNKNOWN'),
            last_stage=stages[-1] if stages else None)
        write(arm/'seal.json',result)
        self.state['cases'].setdefault(case,{})[route]=dict(result,seal_sha256=sha(arm/'seal.json'))
        self.update()
        return result

    def compare(self,case,maximum):
        started=time.monotonic()
        self.current_deadline=started+maximum-3
        classarm=self.folder/case/'class';v3arm=self.folder/case/'v3'
        cs,vs=read(classarm/'seal.json'),read(v3arm/'seal.json')
        for arm,s in [(classarm,cs),(v3arm,vs)]:
            for name,h in s['files'].items():assert sha(arm/name)==h
        c=read(classarm/'verified.json') if cs['verified'] else {}
        v=read(v3arm/'terminal.json') if vs['verified'] else {}
        folder=self.folder/case/'compare';folder.mkdir()
        write(folder/'input.json',dict(packet=read(self.folder/'inputs'/(case+'.json')),class_verified=cs['verified'],v3_verified=vs['verified'],**{'class':c,'v3':v},seals=dict(class_sha256=sha(classarm/'seal.json'),v3_sha256=sha(v3arm/'seal.json'))))
        self.case_cost=0
        ok=self.stage(folder,'transfer_compare.py','compare',maximum)
        overhead=max(0,time.monotonic()-started-self.case_cost)
        self.case_cost+=overhead;self.total+=overhead
        result=read(folder/'comparison.json') if ok else dict(status='UNKNOWN_UNION_VERIFICATION_INCOMPLETE')
        result['charged_seconds']=self.case_cost
        self.state['cases'][case]['comparison']=result;self.update()
        return self.case_cost

    def paired(self,case,cap):
        union_reserve=min(60,.05*cap)
        prior=self.p.get('prior_attempt_charges',{}) if case=='control' else {}
        prior_class=prior.get('class',0);prior_v3=prior.get('v3',0)
        c=self.arm(case,'class',cap-union_reserve-prior_class)
        self.arm(case,'v3',cap-prior_v3)
        union=self.compare(case,min(union_reserve,cap-c['charged_seconds']-prior_class))
        return c,c['charged_seconds']+union+prior_class

    def execute(self):
        if self.p.get('baseline_only_after_closed_class_gate'):
            prior=self.p['prior_attempt_charges']
            self.state['cases']['control']={'class':read(self.folder/'control/class/seal.json')}
            self.arm('control','v3',self.p['reference_cap_seconds']-prior['v3'])
            self.compare('control',min(60,self.p['reference_cap_seconds']-prior['class']))
            self.update(status='STOPPED_COMMISSIONING_GATE',reason='Class bank closed without a new class. Independent V3 reference endpoint retained; no fresh fibres or larger bank launched.')
            return
        control,C=self.paired('control',self.p['reference_cap_seconds'])
        if not control['transfer_success'] or C>self.p['reference_cap_seconds']:
            self.update(status='STOPPED_COMMISSIONING_GATE',reason='Reference class route did not independently complete class construction, blind lifting and rank gain within the authorized allowance. No fresh route launched.')
            return
        B=min(2*C,self.p['reference_cap_seconds'])
        write(self.folder/'fresh-budget.json',dict(C_complete_class_route_seconds=C,B_per_route_seconds=B,control_seal_sha256=sha(self.folder/'control/class/seal.json')))
        self.update(C_seconds=C,B_seconds=B)
        successes=0
        for item in self.p['roster'][1:3]:
            c,_=self.paired(item['id'],B);successes+=int(c['transfer_success'])
        if not successes:
            self.update(status='STOPPED_FIRST_TWO_TRANSFER_GATE',reason='No certified fresh class-to-point success in the first two fixed fibres. Remaining six not started; no probability or Sha conclusion.')
            return
        self.update(remaining_six_gate='PASS_AT_LEAST_ONE_FRESH_END_TO_END_SUCCESS')
        for item in self.p['roster'][3:]:self.paired(item['id'],B)
        self.update(status='COMPLETED_FROZEN_PILOT',reason='All authorized stages complete. No campaign expansion or rank32 claim follows automatically.')


def run(folder):
    controller=Controller(folder)
    def stop_handler(signum,frame):raise KeyboardInterrupt('user stopped detached pilot')
    signal.signal(signal.SIGTERM,stop_handler);signal.signal(signal.SIGINT,stop_handler)
    try:controller.execute()
    except KeyboardInterrupt as exc:controller.update(status='STOPPED_BY_USER',reason=str(exc))
    except Exception as exc:
        import traceback
        traceback.print_exc()
        controller.update(status='STOPPED_CONTROLLER_ERROR',reason=type(exc).__name__+': '+str(exc))
        raise


def status(folder):
    state=read(folder/'state.json');launch=read(folder/'launch.json') if (folder/'launch.json').exists() else {}
    live=bool(launch.get('pid') and token(launch['pid'])==launch.get('start_token'))
    print(state['status'],'| controller',launch.get('pid'),'live' if live else 'not live')
    print('Charged route wall: %.2f h / 36 h maximum'%(state.get('charged_seconds',0)/3600))
    if state.get('reason'):print(state['reason'])
    if state.get('active_arm') and state['status']=='RUNNING':
        elapsed=time.time()-state['stage_started_unix']
        print('Active:',state['active_case'],state['active_route'],state['active_stage'],'| %.1f min / %.1f min stage allowance'%(elapsed/60,state['stage_allowance_seconds']/60))
        arm=Path(state['active_arm']);progress_file=arm/'progress.json'
        if progress_file.exists():
            progress=read(progress_file)
            summary={k:progress[k] for k in ['stage','epoch','charts','rank_lower_bound','atoms','dependencies','new_classes','matrix'] if k in progress}
            print('Latest progress:',json.dumps(summary))
        print('Log:',arm/(state['active_stage']+'.log'))
    for case,routes in state.get('cases',{}).items():
        for route,row in routes.items():
            if route=='comparison':print(case,'union:',row.get('union_rank_lower_bound','UNKNOWN'),'V3:',row.get('v3_rank_lower_bound','UNKNOWN'))
            else:print(case,route,row['status'],'%.1f min'%(row['charged_seconds']/60))


def stop(folder):
    (folder/'STOP').write_text('User stop requested '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'\n')
    launch=read(folder/'launch.json')
    if token(launch['pid'])==launch['start_token']:
        os.kill(launch['pid'],signal.SIGTERM);print('Stop sent to owned controller; active worker tree will be terminated and checkpoints retained.')
    else:print('Controller is already stopped; checkpoints retained.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['launch','run','status','stop']);ap.add_argument('--folder',type=Path,default=DEFAULT)
    args=ap.parse_args();globals()[args.mode](args.folder.resolve())
