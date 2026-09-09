#!/usr/bin/env python3
"""Merge the completed original and foundry-successor conductor queues."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile

from v3_warm_support import atomic, read, require, sha
import run_inventory_conductor_long_v3 as long_pass
from research_curve_refresh import ROOT, MANIFEST

ART = ROOT/'artifacts/generated-results/elliptic-curves'
V3 = ART/'inventory291_conductor_snapshot_v3.json'
LONG = ROOT/'artifacts/local/elliptic-curves/inventory291-conductors-long-v3'
QUEUE = ROOT/'artifacts/local/elliptic-curves/foundry-conductors-queue-v1'
DEFAULT = ART/'inventory321_conductor_snapshot_v4.json'


def sealed(parent, identifier):
    path = parent/'cases'/identifier/'result.json'
    require(path.exists(), 'missing sealed result: '+identifier)
    result = read(path)
    require(result['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY','unreplayed result: '+identifier)
    for name,digest in result['bindings'].items():
        require(sha(ROOT/name)==digest,'result binding changed: '+name)
    return dict(input=read(parent/'cases'/identifier/'input.json'), result=result,
                origin=str(path.relative_to(ROOT)))


def write(path):
    require(not path.exists(),'preserve immutable snapshot; choose a new version')
    previous = read(V3)
    records=dict(previous['records'])
    long_ids=read(LONG/'protocol.json')['ids']
    for identifier in long_ids: records[identifier]=sealed(LONG,identifier)
    queue_ids=read(QUEUE/'protocol.json')['ids']
    for identifier in queue_ids: records[identifier]=sealed(QUEUE,identifier)
    require(len(records)==191,'expected161 original plus30 foundry cases')
    value=dict(status='PASS_MERGED_CONDUCTOR_SNAPSHOT',records=records,
               record_count=len(records),exact_count=sum(
                   r['result']['certificate']['status']=='EXACT' for r in records.values()),
               source_snapshot_sha256=sha(V3),
               protocol_hashes={str((p/'protocol.json').relative_to(ROOT)):sha(p/'protocol.json')
                                for p in (LONG,QUEUE)},
               selection='Latest independently replayed sealed result per published curve at merge cutoff. Original V1/V2 certificates are retained through snapshot v3; completed long-pass and successor-queue results replace only their own records. No live checkpoints.')
    atomic(path,value,immutable=True)
    print('FINAL_CONDUCTOR_SNAPSHOT_WRITTEN',len(records),value['exact_count'],'exact',flush=True)


def check(path):
    from inventory_conductor_worker import run as integral_run
    from inventory_conductor_worker_v2 import run as transported_run
    value=read(path)
    require(value['status']=='PASS_MERGED_CONDUCTOR_SNAPSHOT','invalid final snapshot')
    require(len(value['records'])==value['record_count']==191,'record count differs')
    with tempfile.TemporaryDirectory(prefix='final-conductor-replay-') as directory:
        folder=Path(directory)
        for identifier,entry in sorted(value['records'].items()):
            result=entry['result'];source=entry['input'];certificate=result['certificate']
            require(source['id']==identifier==certificate['id'],'identity differs: '+identifier)
            require(result['status']=='PASS_INDEPENDENT_CONDUCTOR_REPLAY','unverified: '+identifier)
            atomic(folder/'input.json',source);atomic(folder/'certificate.json',certificate)
            (transported_run if 'original_curve' in source else integral_run)(
                folder/'input.json',folder/'certificate.json',True)
    require(value['exact_count']==sum(r['result']['certificate']['status']=='EXACT'
                                     for r in value['records'].values()),'exact count differs')
    print('PASS_FINAL_CONDUCTOR_SNAPSHOT',len(value['records']),value['exact_count'],'exact',flush=True)


def index(path, claim):
    value=read(path);manifest=read(MANIFEST);rel=str(path.relative_to(ROOT))
    claims={r['id']:r for r in read(ROOT/'MATH_STATUS.json')['entries']}
    require(claims[claim]['state']=='proved' and rel in claims[claim]['software_lock'],'proved claim binding required')
    entries={r['id']:r for r in manifest['conductors']}
    for identifier in value['records']:
        entries[identifier]=dict(id=identifier,source=rel,claim=claim,kind='checkpoint')
    manifest['conductors']=list(entries.values());manifest['sources'][rel]=sha(path)
    atomic(MANIFEST,manifest);print('FINAL_CONDUCTOR_INDEXED',len(value['records']))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['write','check','index'])
    parser.add_argument('--snapshot',type=Path,default=DEFAULT)
    parser.add_argument('--claim',default='EC-INVENTORY321-CONDUCTOR-SNAPSHOT4-20260909')
    args=parser.parse_args();path=args.snapshot.resolve()
    if args.action=='write':write(path)
    elif args.action=='check':check(path)
    else:index(path,args.claim)
