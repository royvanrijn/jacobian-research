#!/usr/bin/env python3
"""Detached successor queue for the30 newly published foundry curves."""
import argparse
import json
import os
from pathlib import Path
import time
import run_inventory_conductor_long_v3 as parent
from v3_warm_support import atomic,read,require,sha

base=parent.base
ROOT,CAS,SAGE=parent.ROOT,parent.CAS,parent.SAGE
D=ROOT/'artifacts/local/elliptic-curves/foundry-conductors-queue-v1'
OUT=base.ART/'foundry_conductors_queue_v1'
SELF=Path(__file__).resolve()
WORKER=CAS/'inventory_conductor_worker_v2.py'
SELECTION=base.ART/'foundry_curve_ledger_snapshot_v1.json'


def freeze():
    require(not (D/'protocol.json').exists(),'preserve frozen queue')
    old=parent.guard();selection=read(SELECTION)
    database=read(ROOT/'elliptic-curves/data/research_curves/database.json')
    rows={r['id']:r for r in database['curves']}
    ids=sorted((i for i,r in selection['records'].items() if r['previous_rank_lower_bound'] is None),
               key=lambda i:(-rows[i]['rank_lower_bound'],i))
    require(len(ids)==30,'expected30 published additions')
    inputs={}
    for identifier in ids:
        row=rows[identifier];selected=selection['records'][identifier]['result']
        require(row['conductor'] is None and row['ainvs']==selected['packet']['curve']
                and row['rank_lower_bound']==selected['rank_lower_bound'], 'published input differs')
        source=dict(id=identifier,curve=row['ainvs'],rank_lower_bound=row['rank_lower_bound'],
                    known_primes=row.get('known_bad_primes') or [],factor_hints=[],
                    previous_divisor=row.get('conductor_divisor'),previous_upper_bound=row.get('conductor_upper_bound'))
        packet=parent.previous.integral_packet(source)
        path=D/'cases'/identifier/'input.json';atomic(path,packet,immutable=True)
        inputs[str(path.relative_to(ROOT))]=sha(path)
    sources=dict(old['sources']);sources[str(SELF.relative_to(ROOT))]=sha(SELF)
    value=dict(schema='foundry-conductor-successor-queue.v1',ids=ids,workers=4,
        build_seconds=1800,replay_seconds=120,rss_bytes=2*1024**3,
        sources=sources,inputs=inputs,sage_launcher_sha256=old['sage_launcher_sha256'],
        predecessor_protocol_sha256=sha(parent.D/'protocol.json'),selection_sha256=sha(SELECTION),
        scope='Only the30 newly published foundry curves. Start after the frozen99-case conductor pass completes and releases its controller lock. Four workers,1800-second build and120-second independent replay. No extra simultaneous conductor workers, no point searches, no auto-publication.')
    atomic(D/'protocol.json',value,immutable=True);atomic(OUT/'protocol.json',value,immutable=True)
    print('FOUNDRY_CONDUCTOR_QUEUE_FROZEN',len(ids),'cases',flush=True)


def guard():
    p=read(D/'protocol.json')
    require(sha(parent.D/'protocol.json')==p['predecessor_protocol_sha256'],'predecessor protocol changed')
    require(sha(SELECTION)==p['selection_sha256'],'selected discoveries changed')
    for group in ('sources','inputs'):
        for name,digest in p[group].items():require(sha(ROOT/name)==digest,'frozen binding changed: '+name)
    require(sha(SAGE.resolve())==p['sage_launcher_sha256'],'Sage changed')
    return p


def preflight():
    from inventory_conductor_worker_v2 import verify_transport
    p=guard()
    for identifier in p['ids']:verify_transport(read(D/'cases'/identifier/'input.json'))
    result=dict(status='PASS_QUEUED_CONDUCTOR_PREFLIGHT',cases=len(p['ids']),protocol_sha256=sha(D/'protocol.json'))
    atomic(D/'preflight.json',result,immutable=True);atomic(OUT/'preflight.json',result,immutable=True)
    print('FOUNDRY_CONDUCTOR_QUEUE_PREFLIGHT_PASS',len(p['ids']),'models; charts=0',flush=True)


def ready_to_start(state,completed,total):
    return state.get('status')=='COMPLETE' and not state.get('active') and completed==total


def predecessor_status():
    p=parent.guard();s=read(parent.D/'state.json')
    s['sealed_cases']=sum((parent.D/'cases'/i/'result.json').exists() for i in p['ids'])
    s['total']=len(p['ids'])
    if 'controller' in s:s['controller_alive']=base.process(s['controller']['pid'])==s['controller']
    return s


def controller(fd):
    waiting=dict(status='QUEUED_AFTER_CONDUCTOR_LONG_V3',controller=base.process(os.getpid()),
                 started_unix=time.time(),active=[])
    predecessor_lock=None
    try:
        while True:
            guard();s=predecessor_status()
            waiting.update(parent_status=s['status'],parent_completed=s['sealed_cases'],parent_total=s['total'])
            if ready_to_start(s,s['sealed_cases'],s['total']):
                try:predecessor_lock=base.lock(parent.D/'controller.lock')
                except RuntimeError:pass # Parent is finishing its last atomic write.
                else:break
            atomic(D/'state.json',waiting)
            time.sleep(15)
        # Keep the predecessor lock until this queue drains: an accidental
        # duplicate predecessor launch cannot add another four workers.
        base.controller(fd)
        fd=None
    except BaseException as error:
        waiting.update(status='STOPPED_REVIEW_REQUIRED',error=str(error))
        atomic(D/'state.json',waiting)
        raise
    finally:
        if fd is not None:
            try:os.close(fd)
            except OSError:pass # Frozen controller closes its own fd on failure.
        if predecessor_lock is not None:os.close(predecessor_lock)


def status():
    p=guard();s=read(D/'state.json') if (D/'state.json').exists() else dict(status='NOT_LAUNCHED')
    if 'controller' in s:s['controller_alive']=base.process(s['controller']['pid'])==s['controller']
    if s.get('controller_alive') is False and s['status'] in ('RUNNING','QUEUED_AFTER_CONDUCTOR_LONG_V3'):
        s['status']='INTERRUPTED_NEEDS_RELAUNCH'
    results=[read(D/'cases'/i/'result.json') for i in p['ids'] if (D/'cases'/i/'result.json').exists()]
    s['summary']=dict(completed=len(results),total=len(p['ids']),
        exact=sum(r.get('certificate',{}).get('status')=='EXACT' for r in results),
        independently_replayed=sum(r['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY' for r in results))
    parent_state=predecessor_status()
    s['predecessor']={k:parent_state.get(k) for k in ('status','controller_alive','sealed_cases','total')}
    print(json.dumps(s,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','preflight','launch','status','controller','case'])
    p.add_argument('--id');p.add_argument('--fd',type=int);a=p.parse_args()
    if a.action=='freeze':freeze()
    elif a.action=='preflight':preflight()
    elif a.action=='status':status()
    else:
        protocol=guard();v=read(D/'preflight.json')
        require(v['status']=='PASS_QUEUED_CONDUCTOR_PREFLIGHT' and v['protocol_sha256']==sha(D/'protocol.json'),'preflight required')
        base.D,base.OUT,base.SELF,base.WORKER=D,OUT,SELF,WORKER;base.guard=guard
        if a.action=='launch':base.launch()
        elif a.action=='controller':controller(a.fd)
        else:
            require(a.id in protocol['ids'],'case outside frozen queue');base.case(a.id)
