#!/usr/bin/env python3
"""Bounded longer factorization pass on the99 unresolved V1/V2 conductors."""
import argparse
import json
import math
import os
from pathlib import Path
import sys
import run_inventory_conductor_resume_v2 as previous
from v3_warm_support import atomic, read, require, sha

# Trusted local discriminant cofactors include more than20,000 decimal digits.
sys.set_int_max_str_digits(0)
base=previous.base
ROOT,CAS,SAGE=previous.ROOT,previous.CAS,previous.SAGE
D=ROOT/'artifacts/local/elliptic-curves/inventory291-conductors-long-v3'
OUT=base.ART/'inventory291_conductors_long_v3'
SELF=Path(__file__).resolve()
WORKER=CAS/'inventory_conductor_worker_v2.py'


def prior_records():
    rows={}
    for identifier in read(previous.OLD/'protocol.json')['ids']:
        parent=previous.D if (previous.D/'cases'/identifier/'result.json').exists() else previous.OLD
        folder=parent/'cases'/identifier
        result=read(folder/'result.json')
        require(result['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY','sealed prior replay required')
        rows[identifier]=(folder,result)
    return rows


def shared_splits(cofactors):
    hints={identifier:set() for identifier in cofactors}
    ids=sorted(cofactors)
    for j,left in enumerate(ids):
        for right in ids[j+1:]:
            g=math.gcd(cofactors[left],cofactors[right])
            for identifier in (left,right):
                if 1<g<cofactors[identifier]:hints[identifier].add(g)
    return hints


def freeze():
    require(not (D/'protocol.json').exists(),'preserve frozen long pass')
    p=previous.guard()
    require(read(previous.D/'state.json')['status']=='COMPLETE','finish previous continuation first')
    locks=[base.lock(parent/'controller.lock') for parent in (previous.OLD,previous.D)]
    try:
        rows=prior_records();inherited={};inputs={};origins={}
        cofactors={i:int(r['certificate']['remaining_cofactor']) for i,(f,r) in rows.items()
                   if r['certificate']['status']=='UNKNOWN'}
        hints=shared_splits(cofactors)
        ids=sorted(cofactors,key=lambda i:(len(str(cofactors[i])),-rows[i][1]['certificate']['rank_lower_bound'],i))
        for identifier,(folder,result) in rows.items():
            fd=base.lock(folder/'case.lock');os.close(fd)
            for name,digest in result['bindings'].items():
                require(sha(ROOT/name)==digest,'prior evidence changed: '+name);inherited[name]=digest
            origins[identifier]=str((folder/'result.json').relative_to(ROOT))
            inherited[origins[identifier]]=sha(folder/'result.json')
            if identifier not in cofactors:continue
            source=read(folder/'input.json');c=result['certificate']
            if 'original_curve' not in source:source=previous.integral_packet(source)
            source['known_primes']=sorted(set(source['known_primes'])|
                {r['prime'] for r in c['local_data']},key=int)
            source['factor_hints']=sorted(set(map(str,source['factor_hints']))|
                {str(g) for g in hints[identifier]}|{c['remaining_cofactor']},key=lambda n:(len(n),n))
            source['previous_divisor']=c['conductor_divisor']
            source['previous_upper_bound']=c['conductor_upper_bound']
            source['prior_result_sha256']=sha(folder/'result.json')
            path=D/'cases'/identifier/'input.json';atomic(path,source,immutable=True)
            inputs[str(path.relative_to(ROOT))]=sha(path)
        sources=dict(p['sources']);sources[str(SELF.relative_to(ROOT))]=sha(SELF)
        value=dict(schema='inventory-missing-conductors.long.v3',ids=ids,workers=4,
            build_seconds=1800,replay_seconds=120,rss_bytes=2*1024**3,
            sources=sources,inputs=inputs,inherited=inherited,origins=origins,
            sage_launcher_sha256=p['sage_launcher_sha256'],
            parents={str((q/'protocol.json').relative_to(ROOT)):sha(q/'protocol.json')
                     for q in (previous.OLD,previous.D)},
            previously_exact=62,original_pass_total=161,inventory_total=291,prepass_inventory_exact=192,
            shared_split_cases=sum(bool(v) for v in hints.values()),
            scope='One longer pass on all99 unresolved conductors;1800-second build and120-second replay per case, four workers. Schedule by residual cofactor digit length, then decreasing certified rank and ID. Reuse saved primes and pairwise-GCD splits, never alter old models or proofs. Timeouts remain UNKNOWN; no promise of complete factorization or new rank.')
        require(len(ids)==99 and len(rows)==161,'unexpected prior census')
        atomic(D/'protocol.json',value,immutable=True);atomic(OUT/'protocol.json',value,immutable=True)
        print('CONDUCTOR_LONG_V3_FROZEN',len(ids),'cases; shared-split cases',value['shared_split_cases'])
    finally:
        for fd in locks:os.close(fd)


def guard():
    p=read(D/'protocol.json')
    for group in ('sources','inputs','inherited','parents'):
        for name,digest in p[group].items():require(sha(ROOT/name)==digest,'frozen binding changed: '+name)
    require(sha(SAGE.resolve())==p['sage_launcher_sha256'],'Sage changed')
    return p


def preflight():
    from inventory_conductor_worker_v2 import verify_transport
    p=guard()
    for identifier in p['ids']:
        packet=read(D/'cases'/identifier/'input.json')
        old=read(ROOT/p['origins'][identifier])['certificate']
        require(packet['curve']==old['curve'],'conductor model changed')
        require({r['prime'] for r in old['local_data']}<=set(packet['known_primes']),'lost saved primes')
        require(old['remaining_cofactor'] in packet['factor_hints'],'lost residual cofactor')
        verify_transport(packet)
    value=dict(status='PASS_LONG_CONDUCTOR_PREFLIGHT',cases=len(p['ids']),protocol_sha256=sha(D/'protocol.json'))
    atomic(D/'preflight.json',value,immutable=True);atomic(OUT/'preflight.json',value,immutable=True)
    print('CONDUCTOR_LONG_V3_PREFLIGHT_PASS',len(p['ids']),'models; charts=0')


def status():
    p=guard();s=read(D/'state.json') if (D/'state.json').exists() else dict(status='NOT_STARTED')
    if 'controller' in s:
        s['controller_alive']=base.process(s['controller']['pid'])==s['controller']
        if s['status']=='RUNNING' and not s['controller_alive']:s['status']='INTERRUPTED_REVIEW_REQUIRED'
    rows=[read(D/'cases'/i/'result.json') for i in p['ids'] if (D/'cases'/i/'result.json').exists()]
    gain=sum(r.get('certificate',{}).get('status')=='EXACT' for r in rows)
    s['summary']=dict(completed=len(rows),total=len(p['ids']),new_exact=gain,
        inventory_exact_available=p['prepass_inventory_exact']+gain,inventory_unresolved=99-gain,
        independently_replayed=sum(r['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY' for r in rows),
        build_seconds=p['build_seconds'],workers=p['workers'])
    print(json.dumps(s,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','preflight','launch','status','case','controller'])
    p.add_argument('--id');p.add_argument('--fd',type=int);a=p.parse_args()
    if a.action=='freeze':freeze()
    elif a.action=='preflight':preflight()
    elif a.action=='status':status()
    else:
        protocol=guard();v=read(D/'preflight.json')
        require(v['status']=='PASS_LONG_CONDUCTOR_PREFLIGHT' and v['protocol_sha256']==sha(D/'protocol.json'),'preflight required')
        base.D,base.OUT,base.SELF,base.WORKER=D,OUT,SELF,WORKER;base.guard=guard
        if a.action=='launch':base.launch()
        elif a.action=='controller':base.controller(a.fd)
        else:
            require(a.id in protocol['ids'],'case outside frozen long pass');base.case(a.id)
