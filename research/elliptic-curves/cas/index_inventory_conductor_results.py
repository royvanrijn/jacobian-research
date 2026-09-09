#!/usr/bin/env python3
"""Publish immutable, portable snapshots of independently replayed conductors."""
import argparse
from pathlib import Path
import tempfile
from v3_warm_support import read, atomic, sha, require
from run_inventory_conductor_pass import D, ART, ROOT, guard
from research_curve_refresh import MANIFEST


def snapshot(path):
    guard();require(not path.exists(),'preserve snapshot; choose a new version')
    records={}
    for identifier in read(D/'protocol.json')['ids']:
        folder=D/'cases'/identifier
        if not (folder/'result.json').exists():continue
        result=read(folder/'result.json')
        if result['status'] != 'PASS_INDEPENDENT_CONDUCTOR_REPLAY':continue
        for name,digest in result['bindings'].items():require(sha(ROOT/name)==digest,'case binding changed')
        records[identifier]={'input':read(folder/'input.json'),'result':result}
    require(records,'no verified conductors yet')
    value=dict(status='PASS_CHECKPOINTED_CONDUCTOR_SNAPSHOT',records=records,
               protocol_sha256=sha(D/'protocol.json'),
               exact_count=sum(r['result']['certificate']['status']=='EXACT' for r in records.values()),
               pass_total=len(read(D/'protocol.json')['ids']))
    atomic(path,value,immutable=True);print('CONDUCTOR_SNAPSHOT',len(records),value['exact_count'],'exact')


def check(path):
    from inventory_conductor_worker import run
    value=read(path)
    require(value['status']=='PASS_CHECKPOINTED_CONDUCTOR_SNAPSHOT','invalid snapshot status')
    with tempfile.TemporaryDirectory(prefix='conductor-snapshot-replay-') as directory:
        folder=Path(directory)
        for identifier,record in value['records'].items():
            require(record['result']['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY','unverified record')
            require(record['input']['id']==identifier==record['result']['certificate']['id'],'identity differs')
            atomic(folder/'input.json',record['input']);atomic(folder/'certificate.json',record['result']['certificate'])
            run(folder/'input.json',folder/'certificate.json',True)
    require(value['exact_count']==sum(r['result']['certificate']['status']=='EXACT' for r in value['records'].values()),'exact count differs')
    print('PASS_PORTABLE_CONDUCTOR_SNAPSHOT',len(value['records']),value['exact_count'],'exact')


def index(path,claim):
    value=read(path);manifest=read(MANIFEST);rel=str(path.relative_to(ROOT))
    claims={r['id']:r for r in read(ROOT/'MATH_STATUS.json')['entries']}
    require(claims[claim]['state']=='proved' and rel in claims[claim]['software_lock'],'proved bound ledger entry required')
    manifest['sources'][rel]=sha(path)
    for identifier in value['records']:
        entry=dict(id=identifier,source=rel,claim=claim,kind='checkpoint')
        if entry not in manifest['conductors']:manifest['conductors'].append(entry)
    atomic(MANIFEST,manifest)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['write','check','index'])
    p.add_argument('--snapshot',type=Path,default=ART/'inventory291_conductor_snapshot_v1.json')
    p.add_argument('--claim',default='EC-INVENTORY291-CONDUCTOR-SNAPSHOT-20260909');a=p.parse_args()
    if a.action=='write':snapshot(a.snapshot)
    elif a.action=='check':check(a.snapshot)
    else:index(a.snapshot,a.claim)
