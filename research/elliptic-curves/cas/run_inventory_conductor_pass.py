#!/usr/bin/env python3
"""Detached four-worker bounded pass over every missing ledger conductor."""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import os
from pathlib import Path
import subprocess
import sys
import time
from v3_warm_support import read, atomic, sha, require
from run_r17_60_panel import lock, process
from research_runtime.supervisor import run as supervise, Limits

CAS=Path(__file__).resolve().parent;ROOT=CAS.parents[1]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/inventory291-conductors-v1'
OUT=ART/'inventory291_conductors_v1'
SELF=Path(__file__).resolve();WORKER=CAS/'inventory_conductor_worker.py'
SAGE=Path.home()/'.local/bin/sage'


def freeze():
    require(not (D/'protocol.json').exists(),'preserve frozen pass')
    dbpath=ROOT/'elliptic-curves/data/research_curves/database.json';db=read(dbpath)
    missing=sorted([r for r in db['curves'] if r['conductor'] is None],key=lambda r:(-r['rank_lower_bound'],r['id']))
    paths=[]
    for row in missing:
        hints=[];known=list(row.get('known_bad_primes') or [])
        old=ROOT/'artifacts/local/elliptic-curves/conductor-inventory-continuation-v2/factors'/row['id']/'state.json'
        if old.exists(): hints=read(old)['factors']
        cert=ROOT/row['conductor_certificate'] if row.get('conductor_certificate') else None
        if cert and cert.exists():
            known += [r['prime'] for r in read(cert).get('local_data',[])]
        packet=dict(id=row['id'],curve=row['ainvs'],rank_lower_bound=row['rank_lower_bound'],
                    known_primes=sorted(set(known),key=int),factor_hints=hints,
                    previous_divisor=row.get('conductor_divisor'),previous_upper_bound=row.get('conductor_upper_bound'))
        path=D/'cases'/row['id']/'input.json';atomic(path,packet,immutable=True);paths.append(path)
    sources=[SELF,WORKER,CAS/'v3_warm_support.py',CAS/'research_runtime/supervisor.py',CAS/'run_r17_60_panel.py']
    protocol=dict(schema='inventory-missing-conductors.v1',inventory_count=db['count'],inventory_sha256=sha(dbpath),
                  ids=[r['id'] for r in missing],workers=4,build_seconds=90,replay_seconds=120,rss_bytes=2*1024**3,
                  sources={str(p.relative_to(ROOT)):sha(p) for p in sources},
                  inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},sage_launcher_sha256=sha(SAGE.resolve()),
                  scope='One bounded post-discovery pass on every missing ledger conductor. Reuse retained factors, checkpoint exact local data, independently replay without factoring. No point search, exact-rank claim, or unbounded factorization.')
    atomic(D/'protocol.json',protocol,immutable=True);atomic(OUT/'protocol.json',protocol,immutable=True)
    print('CONDUCTOR_PASS_FROZEN',len(missing),'targets',flush=True)


def guard():
    p=read(D/'protocol.json')
    for group in ('sources','inputs'):
        for name,digest in p[group].items():require(sha(ROOT/name)==digest,'frozen binding changed: '+name)
    require(sha(SAGE.resolve())==p['sage_launcher_sha256'],'Sage changed')
    return p


def case(identifier):
    f=D/'cases'/identifier;fd=lock(f/'case.lock')
    try:
        if (f/'result.json').exists():return
        p=guard();output=f/'conductor.json';build=f/'build.supervisor.json'
        argv=[str(SAGE),'-python','-u',str(WORKER),'--packet',str(f/'input.json'),'--output',str(output)]
        if not build.exists() or read(build)['outcome']=='running':
            require(not output.exists(),'interrupted build requires review; preserve partial evidence')
            rec=supervise(argv,limits=Limits(wall_seconds=p['build_seconds'],rss_bytes=p['rss_bytes']),
                          log_path=f/'build.log',checkpoint_path=build,cwd=ROOT)
        else:rec=read(build)
        require(rec['outcome'] in ('completed','strict_wall_timeout','strict_rss_limit'),'build engineering failure')
        if not output.exists():
            result={'id':identifier,'status':'RESOURCE_UNRESOLVED','build_outcome':rec['outcome']}
        else:
            replay=supervise(argv+['--check'],limits=Limits(wall_seconds=p['replay_seconds'],rss_bytes=p['rss_bytes']),
                             log_path=f/'replay.log',checkpoint_path=f/'replay.supervisor.json',cwd=ROOT)
            require(replay['outcome']=='completed','independent replay failed')
            result={'id':identifier,'status':'PASS_INDEPENDENT_CONDUCTOR_REPLAY','certificate':read(output),
                    'protocol_sha256':sha(D/'protocol.json'),'build_outcome':rec['outcome'],
                    'bindings':{str(x.relative_to(ROOT)):sha(x) for x in (f/'input.json',output,build,f/'replay.supervisor.json')}}
        guard();atomic(f/'result.json',result,immutable=True);atomic(OUT/(identifier+'.json'),result,immutable=True)
        print('CONDUCTOR_CASE_COMPLETE',identifier,result.get('certificate',{}).get('status',result['status']),flush=True)
    finally:os.close(fd)


def report():
    p=guard();rows=[read(D/'cases'/i/'result.json') for i in p['ids'] if (D/'cases'/i/'result.json').exists()]
    s=dict(status='COMPLETE' if len(rows)==len(p['ids']) else 'RUNNING',completed=len(rows),total=len(p['ids']),
           exact=sum(r.get('certificate',{}).get('status')=='EXACT' for r in rows),
           independently_replayed=sum(r['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY' for r in rows),
           certificates={r['id']:{'path':str((OUT/(r['id']+'.json')).relative_to(ROOT)),
                                 'sha256':sha(OUT/(r['id']+'.json'))} for r in rows},protocol_sha256=sha(D/'protocol.json'))
    atomic(D/'summary.json',s)
    if s['status']=='COMPLETE':atomic(OUT/'summary.json',s,immutable=True)
    return s


def dispatch(identifier):
    f=D/'cases'/identifier
    with (f/'worker.log').open('ab',buffering=0) as log:
        code=subprocess.call([sys.executable,'-u',str(SELF),'case','--id',identifier],cwd=ROOT,
                             stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
    require(code==0,'case failed: '+identifier)


def controller(fd):
    state=dict(status='RUNNING',controller=process(os.getpid()),started_unix=time.time(),active=[])
    try:
        pending=[i for i in guard()['ids'] if not (D/'cases'/i/'result.json').exists()]
        with ThreadPoolExecutor(max_workers=4) as pool:
            active={}
            while pending or active:
                while pending and len(active)<4:
                    i=pending.pop(0);active[pool.submit(dispatch,i)]=i
                state.update(active=list(active.values()));atomic(D/'state.json',state)
                done,_=wait(active,timeout=5,return_when=FIRST_COMPLETED)
                for future in done:future.result();del active[future];report()
        state.update(status='COMPLETE',active=[],finished_unix=time.time());atomic(D/'state.json',state)
    except BaseException as exc:
        state.update(status='STOPPED_REVIEW_REQUIRED',error=str(exc));atomic(D/'state.json',state);raise
    finally:os.close(fd)


def launch():
    p=guard();fd=lock(D/'controller.lock')
    try:
        for i in p['ids']:
            check=lock(D/'cases'/i/'case.lock');os.close(check)
        with (D/'controller.log').open('ab',buffering=0) as log:
            child=subprocess.Popen([sys.executable,'-u',str(SELF),'controller','--fd',str(fd)],cwd=ROOT,
                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,pass_fds=(fd,))
        print('CONDUCTOR_PASS_DETACHED',child.pid,flush=True)
    finally:os.close(fd)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('action',choices=['freeze','launch','status','controller','case'])
    a.add_argument('--id');a.add_argument('--fd',type=int);args=a.parse_args()
    if args.action=='freeze':freeze()
    elif args.action=='launch':launch()
    elif args.action=='case':case(args.id)
    elif args.action=='controller':controller(args.fd)
    else:
        s=read(D/'state.json') if (D/'state.json').exists() else {'status':'NOT_STARTED'}
        if 'controller' in s:s['controller_alive']=process(s['controller']['pid'])==s['controller']
        if (D/'summary.json').exists():s['summary']={k:v for k,v in read(D/'summary.json').items() if k!='certificates'}
        print(__import__('json').dumps(s,indent=2))
