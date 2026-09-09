#!/usr/bin/env python3
"""Versioned continuation of the frozen conductor pass, without rerunning sealed cases."""
import argparse
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import run_inventory_conductor_pass as base
from v3_warm_support import atomic, read, require, sha

ROOT, CAS, SAGE = base.ROOT, base.CAS, base.SAGE
OLD, OLD_OUT = base.D, base.OUT
D = ROOT/'artifacts/local/elliptic-curves/inventory291-conductors-v2'
OUT = base.ART/'inventory291_conductors_v2'
SELF = Path(__file__).resolve()
WORKER = CAS/'inventory_conductor_worker_v2.py'


def integral_packet(source):
    curve = list(map(Fraction, source['curve']))
    scale = math.lcm(*(a.denominator for a in curve))
    integral = [a*scale**w for a,w in zip(curve, (1,2,3,4,6))]
    require(len(curve) == 5 and all(a.denominator == 1 for a in integral), 'invalid integral model')
    return dict(source, original_curve=source['curve'], integral_model_scale=str(scale),
                curve=[str(a) for a in integral],
                coordinate_map='x_new=scale^2*x_original; y_new=scale^3*y_original')


def freeze():
    require(not (D/'protocol.json').exists(), 'preserve frozen continuation')
    p=base.guard()
    fd=base.lock(OLD/'controller.lock')
    try:
        inputs={}; inherited={}; pending=[]; originals={}
        for identifier in p['ids']:
            f=OLD/'cases'/identifier
            lock=base.lock(f/'case.lock');os.close(lock)
            originals[str((f/'input.json').relative_to(ROOT))]=sha(f/'input.json')
            if (f/'result.json').exists():
                result=read(f/'result.json')
                require(result['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY', 'unverified retained result')
                require(read(OLD_OUT/(identifier+'.json'))==result, 'retained copies differ')
                for name,digest in result['bindings'].items():
                    require(sha(ROOT/name)==digest, 'retained evidence changed: '+name)
                    inherited[name]=digest
                for path in (f/'result.json', OLD_OUT/(identifier+'.json')):
                    inherited[str(path.relative_to(ROOT))]=sha(path)
            else:
                require(not (f/'conductor.json').exists(), 'unsealed old checkpoint needs review')
                path=D/'cases'/identifier/'input.json'
                atomic(path, integral_packet(read(f/'input.json')), immutable=True)
                inputs[str(path.relative_to(ROOT))]=sha(path);pending.append(identifier)
        sources=dict(p['sources'])
        for path in (SELF,WORKER):sources[str(path.relative_to(ROOT))]=sha(path)
        value=dict(p, schema='inventory-missing-conductors.continuation.v2', ids=pending,
                   sources=sources, inputs=inputs, original_inputs=originals, inherited=inherited,
                   parent_protocol_sha256=sha(OLD/'protocol.json'),
                   parent_total=len(p['ids']), retained_count=len(p['ids'])-len(pending),
                   scope='Resume only unsealed V1 cases. Exact rational denominator clearing with a checked Q-isomorphism; unchanged local conductor algorithms and90/120-second build/replay budgets. Preserve all sealed results and original models.')
        atomic(D/'protocol.json',value,immutable=True);atomic(OUT/'protocol.json',value,immutable=True)
        print('CONDUCTOR_V2_FROZEN',len(pending),'pending;',value['retained_count'],'retained')
    finally:os.close(fd)


def guard():
    p=read(D/'protocol.json')
    require(sha(OLD/'protocol.json')==p['parent_protocol_sha256'], 'parent protocol changed')
    for group in ('sources','inputs','original_inputs','inherited'):
        for name,digest in p[group].items():
            require(sha(ROOT/name)==digest, 'frozen binding changed: '+name)
    require(sha(SAGE.resolve())==p['sage_launcher_sha256'], 'Sage changed')
    return p


def preflight():
    from inventory_conductor_worker_v2 import verify_transport
    p=guard()
    for identifier in p['ids']:
        source=read(D/'cases'/identifier/'input.json')
        require(source==integral_packet(read(OLD/'cases'/identifier/'input.json')), 'input transport differs')
        verify_transport(source)
    value=dict(status='PASS_ZERO_SEARCH_MODEL_PREFLIGHT',cases=len(p['ids']),
               protocol_sha256=sha(D/'protocol.json'))
    atomic(D/'preflight.json',value,immutable=True);atomic(OUT/'preflight.json',value,immutable=True)
    print('CONDUCTOR_V2_PREFLIGHT_PASS',len(p['ids']),'models; charts=0')


def configure():
    # Reuse the byte-frozen supervision/locking code; child dispatch calls this
    # versioned entry point. All source hashes, including V1, remain frozen.
    base.D,base.OUT,base.SELF,base.WORKER=D,OUT,SELF,WORKER
    base.guard=guard


def status():
    p=guard()
    state=read(D/'state.json') if (D/'state.json').exists() else dict(status='NOT_STARTED')
    if 'controller' in state:
        state['controller_alive']=base.process(state['controller']['pid'])==state['controller']
        if state['status']=='RUNNING' and not state['controller_alive']:
            state['status']='INTERRUPTED_REVIEW_REQUIRED'
    rows=[read(OLD/'cases'/i/'result.json') for i in read(OLD/'protocol.json')['ids']
          if (OLD/'cases'/i/'result.json').exists()]
    new=[read(D/'cases'/i/'result.json') for i in p['ids'] if (D/'cases'/i/'result.json').exists()]
    rows+=new
    state['summary']=dict(completed=len(rows),total=p['parent_total'],
        retained_from_v1=p['retained_count'],completed_in_v2=len(new),
        exact=sum(r.get('certificate',{}).get('status')=='EXACT' for r in rows),
        independently_replayed=sum(r['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY' for r in rows))
    print(json.dumps(state,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('action',choices=['freeze','preflight','launch','status','controller','case'])
    p.add_argument('--id');p.add_argument('--fd',type=int);a=p.parse_args()
    if a.action=='freeze':freeze()
    elif a.action=='preflight':preflight()
    elif a.action=='status':status()
    else:
        protocol=guard();verified=read(D/'preflight.json')
        require(verified['status']=='PASS_ZERO_SEARCH_MODEL_PREFLIGHT'
                and verified['protocol_sha256']==sha(D/'protocol.json'), 'preflight required')
        configure()
        if a.action=='launch':base.launch()
        elif a.action=='controller':base.controller(a.fd)
        else:
            require(a.id in protocol['ids'], 'case outside frozen continuation')
            base.case(a.id)
