#!/usr/bin/env python3
"""Detached fibration construction, fixed panels, and adaptive parent promotion.

No network/model calls or daily budgets. The guardian starts another bounded
controller run after both clean completion and recoverable failure. Explicit
stop drains existing workers. Frozen source bytes and exact evidence survive.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import fcntl
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
import zipfile

from v3_warm_support import read,sha,atomic,require,process_info,same_process
from research_runtime.supervisor import run as supervise,Limits
from parent_foundry_policy import (PANEL_SIZE,BASE_CALLS,parameter,tails,promote,
    utility,continuation,exploit_allowance,guardian_should_restart)

CAS=Path(__file__).resolve().parent;ROOT=CAS.parents[1]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DEFAULT=ROOT/'artifacts/local/elliptic-curves/parent-foundry-v2'


def env():
    return {**{k:v for k,v in os.environ.items() if not any(x in k.upper() for x in
        ('API_KEY','ACCESS_TOKEN','AUTH_TOKEN'))},'OPENBLAS_NUM_THREADS':'1',
        'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','PYTHONUNBUFFERED':'1'}


def lock(path):
    stream=path.open('a');fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB);return stream


def freeze(folder,workers):
    require(not folder.exists(),'preserve existing campaign; use launch or a new folder')
    frozen=folder/'runtime/research';frozen.mkdir(parents=True)
    paths=[]
    for directory in (ROOT/'elliptic-curves/cas',ROOT/'elliptic-curves/ecsearch',ROOT/'elkies-k3/scripts'):
        paths.extend(p for p in directory.rglob('*') if p.is_file() and
            p.suffix in ('.py','.sage','.gp','.c','.cpp','.h','.sh') and '__pycache__' not in p.parts)
    paths.extend((ROOT/'elliptic-curves').glob('*.py'))
    inputs=[ART/'compact_six_r17_atlas_v1.json',ART/'compact_five_mw16_atlas_v1.json',
        ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json',
        ART/'curve302_parent_degree2_multisection_orbits_v1.tsv',
        ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json',
        ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv']
    for p in paths+inputs:
        dest=frozen/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    # Catalogue only participates in post-result novelty and conductor screening.
    old=ROOT/'artifacts/local/elliptic-curves/high-rank-foundry-v3'
    previous=read(old/'config.json');oldroot=Path(previous['frozen_root'])
    catalogue=read(oldroot/previous['catalogue'])
    prior=[]
    for path in sorted(ART.glob('high-rank-foundry*/job-*.json')):
        r=read(path)
        if r.get('status')=='PASS_CERTIFIED_SEARCH':prior.append(r['packet']['curve'])
    prior.extend(r['ainvs'] for r in read(ROOT/'elliptic-curves/data/research_curves/database.json')['curves'])
    atomic(folder/'catalogue.json',catalogue,immutable=True)
    atomic(folder/'prior-equations.json',prior,immutable=True)
    config={'schema':'parent-foundry.configuration.v1','workers':workers,
        'sage':shutil.which('sage'),'frozen_root':str(frozen),'export':str(ART/folder.name),
        'panel_size':PANEL_SIZE,'panel_point_calls':BASE_CALLS,'height':125000,
        'construction_seconds':600,'evaluation_seconds':7200,'rss_bytes':3*1024**3,
        'min_free_gib':20,'run_jobs':64,'baseline_every':5,'exploit_share':0.30,
        'daily_budget':None,'created_at':time.time(),
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'scope':'Construct genuine A1 neighbours on X948 and X1092; exact MW16 gate; fixed panels; parent and fibre promotion. No artificial deletion cores or record-driven selection.'}
    atomic(folder/'config.json',config,immutable=True)
    data=frozen/'parent-inputs';data.mkdir()
    original_parent=ART/'curve302_recovered_mw17_parent_v1.json'
    atomic(data/'det1092-gram.json',{'generic_height_gram':read(original_parent)['generic_height_gram'],
        'source_sha256':sha(original_parent),'boundary':'Generic Gram only; no record equation, target parameter or public-group embedding.'},immutable=True)
    result=supervise([config['sage'],'-python',str(frozen/'elliptic-curves/cas/parent_foundry_inputs.sage'),
        '--output',str(data/'generic.json')],limits=Limits(300,3*1024**3),
        log_path=folder/'prepare.log',cwd=frozen,env=env())
    require(result['outcome']=='completed' and result['returncode']==0,'generic input freeze failed')
    generic=read(data/'generic.json');parents={}
    for row in generic['baselines']:
        row={**row,'baseline':True};path=data/(row['family']+'.json');atomic(path,row,immutable=True)
        parents[row['family']]={'path':str(path.relative_to(frozen)),'sha256':sha(path),
            'baseline':True,'source_surface':row['source_surface'],'generic_rank':row['generic_rank_lower_bound'],
            'target_slots':PANEL_SIZE,'dispatched':0,'panel':[],'promotions':0}
    state={'status':'PREPARED','parents':parents,'fibres':{},'jobs':{},'cursor':0,'next_job':0,
        'evaluation_dispatches':0,'runs':0,'epoch':0,'known_invariants':[r['invariant']['key'] for r in generic['known_a1']],
        'seconds':{'panel':0,'exploit':0,'construct':0},'created_at':time.time()}
    atomic(folder/'state.json',state)
    manifest={'files':{str(p.relative_to(frozen)):sha(p) for p in sorted(frozen.rglob('*')) if p.is_file()},
        'executables':{str(p):sha(p) for p in (Path(config['sage']),Path('/usr/bin/gp'))},
        'config_sha256':sha(folder/'config.json')}
    atomic(folder/'manifest.json',manifest,immutable=True)
    output=Path(config['export']);output.mkdir(parents=True,exist_ok=False)
    for name in ('config.json','manifest.json'):shutil.copyfile(folder/name,output/name)
    with zipfile.ZipFile(output/'frozen-runtime.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for name in manifest['files']:archive.write(frozen/name,'research/'+name)
    print(json.dumps({'status':'PREPARED','folder':str(folder),'proposals':len(generic['proposals'])}))


def guard(folder):
    config=read(folder/'config.json');manifest=read(folder/'manifest.json');root=Path(config['frozen_root'])
    require(sha(folder/'config.json')==manifest['config_sha256'],'configuration changed')
    for name,h in manifest['files'].items():require(sha(root/name)==h,'frozen input/source changed: '+name)
    for name,h in manifest['executables'].items():require(sha(Path(name))==h,'executable changed: '+name)
    return config,root


def schedule(folder,state,config,root,generic):
    active=[j for j in state['jobs'].values() if j['status']=='RUNNING']
    constructing=any(j['kind']=='construct' for j in active)
    backlog=sum(not p['baseline'] and p['dispatched']<p['target_slots'] for p in state['parents'].values())
    if not constructing and backlog<4 and state['cursor']<len(generic['proposals']):
        proposal=generic['proposals'][state['cursor']];state['cursor']+=1
        source=next(s for s in generic['sources'] if s['family']==proposal['source_family'])
        return {'kind':'construct',**proposal,'source':source}
    busy={j.get('fibre') for j in active}
    candidates=[(key,f) for key,f in state['fibres'].items() if key not in busy and continuation(f) and not f.get('quarantined')]
    seconds=state['seconds'];total=seconds['panel']+seconds['exploit']
    if candidates and (total==0 or seconds['exploit']<config['exploit_share']*total):
        key,f=max(candidates,key=lambda pair:(utility(pair[1]),pair[0]))
        p=state['parents'][f['family']]
        return {'kind':'exploit','fibre':key,'family':f['family'],'parent':p['path'],'parent_sha256':p['sha256'],
            'parameter':f['parameter'],'packet':f['packet'],'packet_sha256':f['packet_sha256'],
            'allowance':exploit_allowance(f),'height':config['height'],'bank_index':f['batches']}
    available=[(key,p) for key,p in state['parents'].items() if p['dispatched']<p['target_slots']]
    if not available:return None
    controls=[pair for pair in available if pair[1]['baseline']]
    new=[pair for pair in available if not pair[1]['baseline']]
    choices=controls if controls and (not new or state['evaluation_dispatches']%config['baseline_every']==0) else new or controls
    family,p=min(choices,key=lambda pair:(pair[1]['dispatched'],pair[0]))
    index=p['dispatched'];p['dispatched']+=1;state['evaluation_dispatches']+=1
    t=parameter(family,index)
    key=hashlib.sha256(f'{family}/{t}'.encode()).hexdigest()[:24]
    require(key not in state['fibres'],'duplicate panel parameter')
    state['fibres'][key]={'family':family,'parameter':t,'generic_rank':p['generic_rank'],'rank':None,
        'panel_index':index,'batches':0,'calls':0,'stale_calls':0}
    return {'kind':'panel','fibre':key,'family':family,'parent':p['path'],'parent_sha256':p['sha256'],
        'parameter':t,'panel_index':index,'allowance':BASE_CALLS,'height':config['height'],'bank_index':0}


def execute(job,request,config,root):
    if request['kind']=='construct':
        command=[config['sage'],'-python',str(root/'elliptic-curves/cas/parent_foundry_geometry.sage'),
                 '--request',str(job/'request.json'),'--output',str(job)]
        seconds=config['construction_seconds']
    else:
        command=[config['sage'],'-python',str(root/'elliptic-curves/cas/parent_foundry_worker.py'),'--job',str(job)]
        seconds=config['evaluation_seconds']
    return supervise(command,limits=Limits(seconds,config['rss_bytes']),log_path=job/'worker.log',
        checkpoint_path=job/'supervision.json',cwd=root,env=env())


def ingest(folder,state,jid,supervision,config,root):
    j=state['jobs'][jid];job=root/j['path'];request=read(job/'request.json')
    elapsed=supervision.get('wall_seconds',time.time()-j['started_at'])
    state['seconds'][request['kind']]+=elapsed
    j.update(status='DONE',finished_at=time.time(),supervision=supervision['outcome'])
    output=Path(config['export'])
    if request['kind']=='construct':
        if (job/'parent.json').exists() and (job/'verified.json').exists() and supervision['returncode']==0:
            p=read(job/'parent.json');v=read(job/'verified.json')
            require(v['parent_sha256']==sha(job/'parent.json'),'parent replay binding differs')
            key=p['novelty_invariant']['key'];j['family']=p['family']
            if key in state['known_invariants']:
                j['outcome']='DUPLICATE_OR_UNRESOLVED_INVARIANT_COLLISION'
            else:
                state['known_invariants'].append(key)
                state['parents'][p['family']]={'path':str((job/'parent.json').relative_to(root)),
                    'sha256':sha(job/'parent.json'),'baseline':False,'source_surface':p['source_surface'],
                    'generic_rank':16,'target_slots':PANEL_SIZE,'dispatched':0,'panel':[],'promotions':0,
                    'coefficient_bits':p['compactification']['after_bits'],'invariant':key}
                j['outcome']='ACCEPTED_NEW_FIBRATION'
            atomic(output/(jid+'-parent.json'),{'parent':p,'replay':v,'novelty_gate':j['outcome'],
                'relative_to':'Five pinned A1 families, seven rootless sources and previous accepted parents; not a universal novelty claim.'},immutable=True)
        elif (job/'not-a1.json').exists():j['outcome']='NOT_A1_STRATUM'
        else:j['outcome']='UNRESOLVED_CONSTRUCTION'
        return
    f=state['fibres'][request['fibre']];p=state['parents'][request['family']]
    if not (job/'result.json').exists() or supervision['returncode']!=0:
        result={'status':'UNRESOLVED_WORKER','rank_lower_bound':None,'exposure_complete':False,
            'calls':None,'job':j['path'],'outcome':supervision['outcome']}
        # One transient restart with the same immutable inputs and saved calls.
        j['status']='RETRY' if j.get('attempt',0)<1 else 'QUARANTINED'
        f['quarantined']=j['status']=='QUARANTINED'
        if j['status']=='RETRY':return
    else:
        result=read(job/'result.json');require(result['request_sha256']==sha(job/'request.json'),'request seal differs')
    if result.get('status')=='PASS_CERTIFIED_PARENT_EVALUATION':
        packet=read(root/result['packet']);v=read(job/'verified.json')
        require(sha(root/result['packet'])==result['packet_sha256']==v['packet_sha256'] and
            v['status']=='PASS_TWO_FINITE_IMPLEMENTATIONS','point certificate replay missing')
        require(len(packet['points'])==packet['rank_lower_bound']==result['rank_lower_bound'],'point rank header differs')
        old=f['rank'] if f['rank'] is not None else f['generic_rank']
        gain=result['rank_lower_bound']-old;require(gain>=0,'certified rank decreased')
        if gain:
            require(result['gain_timeline'],'gain has no chart provenance')
            last=max(e['call'] for e in result['gain_timeline'])
            stale=result['calls']-last
        else:stale=f['stale_calls']+result['calls']
        f.update(rank=result['rank_lower_bound'],calls=f['calls']+result['calls'],
            stale_calls=stale,last_gain=gain,batches=f['batches']+1,
            empty_batches=(f.get('empty_batches',0)+1 if result['calls']==0 else 0),
            packet=result['packet'],packet_sha256=result['packet_sha256'])
        from high_rank_foundry_intake import novelty,conductor_gate,matches,jkey
        catalogue=read(folder/'catalogue.json')
        result['catalogue_novelty']=novelty(packet['curve'],catalogue)
        prior=read(folder/'prior-equations.json')
        result['previous_repository_match']=bool(matches(packet['curve'],prior))
        result['conductor_gate']=conductor_gate(packet['curve'],result['rank_lower_bound'],catalogue)
        result['packet']=packet;result['replay']=v
        # Exact isomorphism check among earlier accepted fibres, after the fixed
        # panel is measured. Duplicates cannot become new curve claims.
        equations=state.setdefault('equations',{})
        jk=jkey(packet['curve']);result['other_foundry_matches']=[x['fibre'] for x in equations.get(jk,[])
            if x['fibre']!=request['fibre'] and matches(packet['curve'],[x['model']])]
        if not any(x['fibre']==request['fibre'] for x in equations.get(jk,[])):
            equations.setdefault(jk,[]).append({'fibre':request['fibre'],'model':packet['curve']})
    if request['kind']=='panel':
        summary={k:result.get(k) for k in ('status','rank_lower_bound','certified_jump_lower_bound','calls','exposure_complete')}
        summary.update(parameter=request['parameter'],panel_index=request['panel_index'],job=j['path'])
        # Replace a quarantined slot only through an explicit later retry receipt.
        p['panel'].append(summary)
        if len(p['panel'])==p['target_slots'] and not p['baseline']:
            if promote(p['panel'][-PANEL_SIZE:]):
                p['target_slots']+=PANEL_SIZE*(1+min(3,p['promotions']));p['promotions']+=1
    j['outcome']=result['status']
    atomic(output/(jid+'-result.json'),result,immutable=True)


def report(folder,state,config):
    rows=[]
    for family,p in state['parents'].items():
        rows.append({'family':family,'baseline':p['baseline'],'source_surface':p['source_surface'],
            'generic_rank':p['generic_rank'],'first_panel':tails(p['panel'][:PANEL_SIZE]),
            'all_panel_slots':len(p['panel']),'target_slots':p['target_slots'],'promotions':p['promotions']})
    best=sorted([{'id':key,**f} for key,f in state['fibres'].items() if f.get('rank') is not None],
        key=lambda f:(-f['rank'],f['calls']))[:12]
    data={'status':state['status'],'updated_at':time.time(),'runs':state['runs'],
        'construction_proposals_attempted':state['cursor'],
        'new_fibrations':sum(not p['baseline'] for p in state['parents'].values()),
        'running_jobs':[{k:j[k] for k in ('path','kind','started_at')} for j in state['jobs'].values() if j['status']=='RUNNING'],
        'job_outcomes':{x:sum(j.get('outcome')==x for j in state['jobs'].values()) for x in
            sorted({j.get('outcome','PENDING') for j in state['jobs'].values()})},
        'seconds':state['seconds'],'parents':rows,'best_fibres':best,
        'boundary':'Tail rates measure certified detection under fixed first-panel exposure. Unequal follow-up is separate. No exact ranks or world-record claims.'}
    output=Path(config['export']);atomic(output/'LIVE_STATUS.json',data)
    text=f"# Autonomous parent foundry\n\nStatus: {state['status']}. Controller runs: {state['runs']}.\n\n"
    text+=f"Constructed and accepted new fibrations: {data['new_fibrations']}; proposals attempted: {state['cursor']}.\n\n"
    text+='Every parent receives twelve fixed prospective slots, each with 64 point calls, bounded maps and full-cloud replay. No daily budget.\n\n'
    text+='| Parent | Generic rank | Slots | Δ≥3 | Δ≥5 | Δ≥8 | Incomplete |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for p in rows:
        t=p['first_panel']
        text+=f"| {p['family']} | {p['generic_rank']} | {t['slots']} | {t['3']['fraction'] or '—'} | {t['5']['fraction'] or '—'} | {t['8']['fraction'] or '—'} | {t['incomplete_exposures']} |\n"
    text+='\n'+data['boundary']+'\n'
    (output/'REPORT.md').write_text(text)


def controller(folder):
    lease=lock(folder/'controller.lock');config,root=guard(folder)
    state=read(folder/'state.json');generic=read(root/'parent-inputs/generic.json')
    state['runs']+=1;state['status']='RUNNING'
    for j in state['jobs'].values():
        if j['status']=='RUNNING':j['status']='RETRY'
    futures={};finished=0
    with ThreadPoolExecutor(max_workers=config['workers']) as pool:
        while True:
            stop=(folder/'STOP').exists();space=shutil.disk_usage(folder).free>=config['min_free_gib']*1024**3
            while not stop and space and len(futures)<config['workers'] and finished+len(futures)<config['run_jobs']:
                retry=next(((jid,j) for jid,j in state['jobs'].items() if j['status']=='RETRY'),None)
                if retry:
                    jid,j=retry;job=root/j['path'];request=read(job/'request.json');j['attempt']=j.get('attempt',0)+1
                else:
                    request=schedule(folder,state,config,root,generic)
                    if request is None:break
                    jid=f'job-{state["next_job"]:07d}';state['next_job']+=1
                    job=root/'parent-jobs'/jid;atomic(job/'request.json',request,immutable=True)
                    j={'path':str(job.relative_to(root)),'kind':request['kind'],'fibre':request.get('fibre'),'attempt':0}
                    state['jobs'][jid]=j
                j.update(status='RUNNING',started_at=time.time());atomic(folder/'state.json',state)
                futures[pool.submit(execute,job,request,config,root)]=jid
            state['status']='DRAINING' if stop else 'WAITING_FOR_DISK' if not space else 'RUNNING'
            report(folder,state,config);atomic(folder/'state.json',state)
            if not futures:
                if stop:state['status']='STOPPED';break
                if finished>=config['run_jobs']:state['status']='RUN_COMPLETE_RENEWING';break
                # Exhausted batches renew with larger panels; no stalled fibre
                # needs an AI decision to keep the outer programme moving.
                if state['cursor']>=len(generic['proposals']) and space:
                    state['epoch']+=1
                    for p in state['parents'].values():
                        if not p['baseline']:p['target_slots']+=PANEL_SIZE
                    state['status']='PARENT_POOL_EXHAUSTED_EXPANDING_PANELS'
                    atomic(folder/'state.json',state)
                time.sleep(5);continue
            done,_=wait(futures,timeout=2,return_when=FIRST_COMPLETED)
            for future in done:
                jid=futures.pop(future)
                try:supervision=future.result()
                except Exception:
                    supervision={'outcome':'SUPERVISOR_EXCEPTION','returncode':-1,'error':traceback.format_exc()}
                ingest(folder,state,jid,supervision,config,root);finished+=1
                atomic(folder/'state.json',state)
    report(folder,state,config);atomic(folder/'state.json',state)
    lease.close()


def guardian(folder):
    lease=lock(folder/'guardian.lock');config,root=guard(folder);crashes=0
    while not (folder/'STOP').exists():
        p=subprocess.Popen([sys.executable,str(root/'elliptic-curves/cas/run_parent_foundry.py'),
            'controller','--folder',str(folder)],cwd=root,env=env())
        atomic(folder/'controller-lease.json',{'pid':p.pid,'token':process_info(p.pid)['start_token']})
        rc=p.wait()
        if not guardian_should_restart((folder/'STOP').exists(),rc):break
        crashes=crashes+1 if rc else 0
        atomic(folder/'renewal.json',{'at':time.time(),'returncode':rc,'consecutive_crashes':crashes,
            'action':'START_NEXT_CONTROLLER_RUN'})
        time.sleep(min(60,2**min(crashes,6)))
    lease.close()


def launch(folder):
    config,root=guard(folder)
    if (folder/'guardian-lease.json').exists():
        p=read(folder/'guardian-lease.json')
        if same_process(p['pid'],p['token']):
            print(json.dumps({'status':'ALREADY_RUNNING',**p}));return
    (folder/'STOP').unlink(missing_ok=True)
    with (folder/'guardian.log').open('ab',buffering=0) as log:
        p=subprocess.Popen([sys.executable,str(root/'elliptic-curves/cas/run_parent_foundry.py'),
            'guardian','--folder',str(folder)],cwd=root,env=env(),stdin=subprocess.DEVNULL,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
    info=process_info(p.pid);require(info is not None,'guardian launch failed')
    atomic(folder/'guardian-lease.json',{'pid':p.pid,'token':info['start_token'],'launched_at':time.time()})
    print(json.dumps({'status':'LAUNCHED','pid':p.pid,'folder':str(folder)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('prepare','launch','controller','guardian','status','stop'))
    parser.add_argument('--folder',type=Path,default=DEFAULT);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();folder=args.folder.resolve()
    if args.mode=='prepare':
        require(1<=args.workers<=4,'worker count outside scope');freeze(folder,args.workers)
    elif args.mode=='launch':launch(folder)
    elif args.mode=='controller':controller(folder)
    elif args.mode=='guardian':guardian(folder)
    elif args.mode=='stop':atomic(folder/'STOP',{'requested_at':time.time()});print('Graceful stop requested.')
    else:
        config=read(folder/'config.json');p=read(folder/'guardian-lease.json') if (folder/'guardian-lease.json').exists() else None
        data=read(Path(config['export'])/'LIVE_STATUS.json');data['guardian_alive']=bool(p and same_process(p['pid'],p['token']))
        print(json.dumps(data,indent=2))
