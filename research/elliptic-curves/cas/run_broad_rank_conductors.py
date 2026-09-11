#!/usr/bin/env python3
"""Bounded conductor intake for the thirty published broad-rank discoveries."""
import argparse
import json
import os
from pathlib import Path
import sys

import run_inventory_conductor_resume_v2 as transport
from publish_broad_rank_results import SNAPSHOT, REPLAY
from v3_warm_support import atomic, read, require, sha

base = transport.base
ROOT, CAS, SAGE = base.ROOT, base.CAS, base.SAGE
D = ROOT/'artifacts/local/elliptic-curves/broad-rank-conductors-v1'
OUT = base.ART/'broad_rank_conductors_v1'
SELF = Path(__file__).resolve()
WORKER = CAS/'inventory_conductor_worker_v2.py'


def freeze():
    require(not (D/'protocol.json').exists(), 'preserve frozen conductor queue')
    snapshot, replay = read(SNAPSHOT), read(REPLAY)
    require(replay['status'] == 'PASS_BROAD_LEDGER_FRESH_RANK_REPLAY'
            and replay['snapshot_sha256'] == sha(SNAPSHOT), 'replayed curve cohort required')
    records = sorted(snapshot['records'], key=lambda r:(-r['rank_lower_bound'],r['id']))
    require(len(records) == 30, 'expected thirty requested discoveries')
    inputs = {}
    for row in records:
        source = dict(id=row['id'],curve=row['packet']['curve'],rank_lower_bound=row['rank_lower_bound'],
                      known_primes=[],factor_hints=[],previous_divisor=None,previous_upper_bound=None)
        path = D/'cases'/row['id']/'input.json'
        atomic(path, transport.integral_packet(source), immutable=True)
        inputs[str(path.relative_to(ROOT))] = sha(path)
    sources = {SELF, WORKER, Path(base.__file__).resolve(), Path(transport.__file__).resolve(),
               CAS/'inventory_conductor_worker.py'}
    # Include all repository modules loaded by these existing controller helpers.
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.is_relative_to(ROOT) and path.is_file():
                sources.add(path)
    protocol = dict(schema='broad-rank.conductor-queue.v1',ids=[r['id'] for r in records],workers=4,
        build_seconds=1800,replay_seconds=120,rss_bytes=2*1024**3,
        sources={str(p.relative_to(ROOT)):sha(p) for p in sorted(sources)},inputs=inputs,
        curve_snapshot_sha256=sha(SNAPSHOT),curve_replay_sha256=sha(REPLAY),
        sage_launcher_sha256=sha(SAGE.resolve()),
        scope='One post-discovery pass over the thirty published broad-rank fibres above lower bound22. Four workers;1800-second factorization/build and120-second independent replay per curve;2GiB per process. Exact transport and existing Sage/PARI local arithmetic. Partial results remain UNKNOWN; no point searches or unbounded retries.')
    atomic(D/'protocol.json',protocol,immutable=True)
    atomic(OUT/'protocol.json',protocol,immutable=True)
    print('BROAD_CONDUCTORS_FROZEN',len(records),flush=True)


def guard():
    p = read(D/'protocol.json')
    require(p['curve_snapshot_sha256']==sha(SNAPSHOT) and p['curve_replay_sha256']==sha(REPLAY),
            'published curve binding changed')
    for group in ('sources','inputs'):
        for name,digest in p[group].items():
            require(sha(ROOT/name)==digest,'frozen conductor binding changed: '+name)
    require(sha(SAGE.resolve())==p['sage_launcher_sha256'],'Sage changed')
    return p


def preflight():
    from inventory_conductor_worker_v2 import verify_transport
    p = guard()
    for identifier in p['ids']:
        verify_transport(read(D/'cases'/identifier/'input.json'))
    receipt = dict(status='PASS_BROAD_CONDUCTOR_PREFLIGHT',cases=len(p['ids']),
                   protocol_sha256=sha(D/'protocol.json'))
    atomic(D/'preflight.json',receipt,immutable=True)
    atomic(OUT/'preflight.json',receipt,immutable=True)
    print('BROAD_CONDUCTOR_PREFLIGHT_PASS',len(p['ids']),flush=True)


def configure():
    base.D,base.OUT,base.SELF,base.WORKER = D,OUT,SELF,WORKER
    base.guard = guard


def status():
    p=guard();state=read(D/'state.json') if (D/'state.json').exists() else {'status':'NOT_STARTED'}
    if 'controller' in state:
        state['controller_alive']=base.process(state['controller']['pid'])==state['controller']
        if not state['controller_alive'] and state['status']=='RUNNING':
            state['status']='INTERRUPTED_REVIEW_REQUIRED'
    results=[read(D/'cases'/i/'result.json') for i in p['ids'] if (D/'cases'/i/'result.json').exists()]
    state['summary']={'total':len(p['ids']),'completed':len(results),
        'exact':sum(r.get('certificate',{}).get('status')=='EXACT' for r in results),
        'independently_replayed':sum(r['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY' for r in results)}
    print(json.dumps(state,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('freeze','preflight','launch','status','case','controller'))
    parser.add_argument('--id');parser.add_argument('--fd',type=int)
    args=parser.parse_args()
    if args.action=='freeze':freeze()
    elif args.action=='preflight':preflight()
    elif args.action=='status':status()
    else:
        protocol=guard();verified=read(D/'preflight.json')
        require(verified['status']=='PASS_BROAD_CONDUCTOR_PREFLIGHT'
                and verified['protocol_sha256']==sha(D/'protocol.json'),'preflight required')
        configure()
        if args.action=='launch':base.launch()
        elif args.action=='controller':base.controller(args.fd)
        else:
            require(args.id in protocol['ids'],'case outside requested cohort')
            base.case(args.id)
