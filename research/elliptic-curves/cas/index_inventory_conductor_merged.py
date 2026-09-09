#!/usr/bin/env python3
"""Publish portable, replayed snapshots across the preserved conductor passes."""
import argparse
from pathlib import Path
import tempfile
import run_inventory_conductor_long_v3 as long_pass
from v3_warm_support import atomic, read, require, sha
from index_inventory_conductor_results import index

ROOT=long_pass.ROOT
DEFAULT=ROOT/'artifacts/generated-results/elliptic-curves/inventory291_conductor_snapshot_v3.json'


def snapshot(path):
    long_pass.guard();require(not path.exists(),'preserve snapshot; use a new version')
    parents=[long_pass.previous.OLD,long_pass.previous.D,long_pass.D]
    # Freeze the list of sealed endpoints before reading payloads. Mutable live
    # conductor checkpoints never enter a publication snapshot.
    selected={}
    for parent in parents:
        for identifier in read(parent/'protocol.json')['ids']:
            folder=parent/'cases'/identifier
            if (folder/'result.json').exists():selected[identifier]=folder
    records={}
    for identifier,folder in sorted(selected.items()):
        result=read(folder/'result.json')
        require(result['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY','unverified sealed endpoint')
        for name,digest in result['bindings'].items():require(sha(ROOT/name)==digest,'case binding changed: '+name)
        records[identifier]=dict(input=read(folder/'input.json'),result=result,
                                origin=str((folder/'result.json').relative_to(ROOT)))
    require(len(records)==161,'complete original conductor roster required')
    value=dict(status='PASS_CHECKPOINTED_CONDUCTOR_SNAPSHOT',records=records,pass_total=161,
               exact_count=sum(r['result']['certificate']['status']=='EXACT' for r in records.values()),
               protocol_hashes={str((p/'protocol.json').relative_to(ROOT)):sha(p/'protocol.json') for p in parents},
               selection='Latest independently replayed sealed result per curve at snapshot cutoff; no live checkpoints.')
    atomic(path,value,immutable=True)
    print('MERGED_CONDUCTOR_SNAPSHOT',len(records),value['exact_count'],'exact',flush=True)


def check(path):
    from inventory_conductor_worker import run as integral_run
    from inventory_conductor_worker_v2 import run as transported_run
    value=read(path)
    require(value['status']=='PASS_CHECKPOINTED_CONDUCTOR_SNAPSHOT','invalid snapshot status')
    require(len(value['records'])==value['pass_total']==161,'snapshot roster count differs')
    with tempfile.TemporaryDirectory(prefix='merged-conductor-replay-') as directory:
        folder=Path(directory)
        for identifier,record in value['records'].items():
            result=record['result'];source=record['input'];certificate=result['certificate']
            require(result['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY','unverified record')
            require(source['id']==identifier==certificate['id'],'identity differs')
            atomic(folder/'input.json',source);atomic(folder/'certificate.json',certificate)
            worker=transported_run if 'original_curve' in source else integral_run
            worker(folder/'input.json',folder/'certificate.json',True)
    require(value['exact_count']==sum(r['result']['certificate']['status']=='EXACT'
                                     for r in value['records'].values()),'exact count differs')
    print('PASS_PORTABLE_MERGED_CONDUCTOR_SNAPSHOT',len(value['records']),value['exact_count'],'exact',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['write','check','index'])
    p.add_argument('--snapshot',type=Path,default=DEFAULT)
    p.add_argument('--claim',default='EC-INVENTORY291-CONDUCTOR-SNAPSHOT3-20260909')
    a=p.parse_args();path=a.snapshot.resolve()
    if a.action=='write':snapshot(path)
    elif a.action=='check':check(path)
    else:index(path,a.claim)
