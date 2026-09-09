#!/usr/bin/env python3
"""Freeze and replay relevant foundry discoveries before adding them to the ledger."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import tempfile
import certify_compact_r17_candidates as cert
from research_curve_refresh import ROOT, MANIFEST
from v3_warm_support import atomic, read, require, sha

ART=ROOT/'artifacts/generated-results/elliptic-curves'
SNAPSHOT=ART/'foundry_curve_ledger_snapshot_v1.json'
CLAIM='EC-FOUNDRY-CURVE-LEDGER-SNAPSHOT-20260909'
DIRECTORIES=[f'high-rank-foundry-integration-v{i}' for i in range(1,5)]+['high-rank-foundry-v1','high-rank-foundry-v2']


def bucket(model):
    inv=cert.weierstrass_invariants(model)
    return inv['c4']**3/inv['discriminant']


def validate_record(record):
    require(record['status']=='PASS_CERTIFIED_SEARCH','certified search endpoint required')
    packet=record['packet'];receipt=record['certificate_replay']
    digest=hashlib.sha256((json.dumps(packet,indent=2,sort_keys=True)+'\n').encode()).hexdigest()
    require(digest==record['packet_sha256']==receipt['packet_sha256'],'portable packet hash differs')
    require(receipt['status']=='PASS_TWO_FINITE_IMPLEMENTATIONS','two independent replays required')
    require(record['rank_lower_bound']==packet['rank_lower_bound']==receipt['rank_lower_bound']==len(packet['points']),
            'certified rank headers differ')


def select(candidates,baseline):
    buckets={};selected={}
    for row in baseline:buckets.setdefault(bucket(row['ainvs']),[]).append(row)
    for origin,record in sorted(candidates,key=lambda x:(-x[1]['rank_lower_bound'],x[0])):
        validate_record(record)
        if record['rank_lower_bound']<22:continue
        model=record['packet']['curve'];peers=buckets.setdefault(bucket(model),[])
        old=next((r for r in peers if cert.isomorphic(model,r['ainvs'])),None)
        if old:
            if record['rank_lower_bound']<=old['rank_lower_bound']:continue
            identifier=old['id'];previous_rank=old['rank_lower_bound']
            # Candidates are descending rank; one strongest endpoint per curve.
            old['rank_lower_bound']=record['rank_lower_bound']
        else:
            address=record['family']+':'+str(F(record['parameter']))
            identifier='foundry-'+record['family']+'-'+hashlib.sha256(address.encode()).hexdigest()[:12]
            require(identifier not in selected,'curve ID collision')
            peers.append(dict(id=identifier,ainvs=model,rank_lower_bound=record['rank_lower_bound']))
            previous_rank=None
        selected[identifier]=dict(origin=origin,origin_sha256=sha(ROOT/origin),
            previous_rank_lower_bound=previous_rank,result=record)
    return selected


def freeze():
    require(not SNAPSHOT.exists(),'preserve snapshot; use a new version for later discoveries')
    path=ROOT/'elliptic-curves/data/research_curves/database.json'
    baseline=[{k:r[k] for k in ('id','ainvs','rank_lower_bound')} for r in read(path)['curves']]
    candidates=[]
    # Pin the list first: exports arriving during replay belong to a later update.
    paths=sorted(p for name in DIRECTORIES for p in (ART/name).glob('job-*.json'))
    for source in paths:
        r=read(source)
        if r.get('status')=='PASS_CERTIFIED_SEARCH' and r['rank_lower_bound']>=22:
            candidates.append((str(source.relative_to(ROOT)),r))
    records=select(candidates,json.loads(json.dumps(baseline)))
    value=dict(status='FROZEN_FOUNDRY_LEDGER_SELECTION',records=records,baseline=baseline,
               baseline_database_sha256=sha(path),threshold=22,
               new_count=sum(r['previous_rank_lower_bound'] is None for r in records.values()),
               strengthened_count=sum(r['previous_rank_lower_bound'] is not None for r in records.values()),
               rank_counts=dict(sorted(Counter(str(r['result']['rank_lower_bound']) for r in records.values()).items())),
               boundary='Distinct over Q relative to this baseline; not worldwide novelty. Strongest verified subgroup per curve, at a fixed publication cutoff. No new point search.')
    atomic(SNAPSHOT,value,immutable=True)
    print('FOUNDRY_LEDGER_FROZEN',value['new_count'],'new;',value['strengthened_count'],'strengthened;',value['rank_counts'],flush=True)


def check():
    from verify_high_rank_foundry_certificate import verify
    value=read(SNAPSHOT);baseline=value['baseline'];seen=[]
    with tempfile.TemporaryDirectory(prefix='foundry-ledger-replay-') as directory:
        for identifier,entry in value['records'].items():
            record=entry['result'];validate_record(record)
            require(record['rank_lower_bound']>=value['threshold'],'below editorial threshold')
            require(sha(ROOT/entry['origin'])==entry['origin_sha256'],'selected source changed')
            require(read(ROOT/entry['origin'])==record,'embedded source differs')
            path=Path(directory)/'certificate.json';atomic(path,record);verify(path)
            model=record['packet']['curve'];j=bucket(model)
            matches=[r for r in baseline if bucket(r['ainvs'])==j and cert.isomorphic(model,r['ainvs'])]
            if entry['previous_rank_lower_bound'] is None:require(not matches,'new curve duplicates baseline')
            else:
                require(len(matches)==1 and matches[0]['id']==identifier
                        and matches[0]['rank_lower_bound']==entry['previous_rank_lower_bound']<record['rank_lower_bound'],
                        'invalid subgroup improvement')
            require(not any(j==other_j and cert.isomorphic(model,other) for other_j,other in seen),'duplicate selected curve')
            seen.append((j,model))
    require(value['new_count']==sum(r['previous_rank_lower_bound'] is None for r in value['records'].values()),'new count differs')
    require(value['strengthened_count']==len(seen)-value['new_count'],'update count differs')
    require(value['rank_counts']==dict(Counter(str(r['result']['rank_lower_bound']) for r in value['records'].values())), 'rank counts differ')
    print('PASS_FOUNDRY_LEDGER_REPLAY',len(seen),'curves;',value['new_count'],'new;',value['strengthened_count'],'strengthened',flush=True)


def index():
    value=read(SNAPSHOT);manifest=read(MANIFEST);rel=str(SNAPSHOT.relative_to(ROOT))
    claim=next(r for r in read(ROOT/'MATH_STATUS.json')['entries'] if r['id']==CLAIM)
    require(claim['state']=='proved' and rel in claim['software_lock'],'proved snapshot claim required')
    entries={r['id']:r for r in manifest['curves']}
    for identifier,entry in value['records'].items():
        r=entry['result']
        entries[identifier]=dict(id=identifier,family=r['family'],parameter=r['parameter'],
            rank_lower_bound=r['rank_lower_bound'],kind='foundry',source=rel,record=identifier,claim=CLAIM)
    manifest['curves']=list(entries.values());manifest['sources'][rel]=sha(SNAPSHOT)
    atomic(MANIFEST,manifest);print('FOUNDRY_LEDGER_INDEXED',len(value['records']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','check','index']);a=p.parse_args()
    {'freeze':freeze,'check':check,'index':index}[a.action]()
