#!/usr/bin/env python3
"""Balanced H4096 extension, reusing only frozen generic censuses and trace tables."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import shutil
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h1024-v1'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1'
BINARY=ROOT/'artifacts/local/elliptic-curves/compact-r17-wide-v1/scanner'

def sources():
    paths=(Path(__file__).resolve(),spec.ATLAS,Path(spec.__file__).resolve(),Path(cert.__file__).resolve(),
           CAS/'newfamily/scan_rational_nagao_tables.cpp',PARENT/'protocol.json',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def bindings(directory):
    p=cert.read(directory/'protocol.json')
    if p['sources']!=sources() or p['scanner_binary_sha256']!=cert.hashed(BINARY):raise ArithmeticError('wide selection source changed')
    return p

def family_record(family):return next(r for r in cert.read(spec.ATLAS)['families'] if r['fibration_id']==family)

def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('wide protocol already frozen')
    parent=cert.read(PARENT/'protocol.json');ledger=cert.read(PARENT/'point-ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS':raise ArithmeticError('initial pilot unfinished')
    # A scheduling gate from prospective results, not a known-record label.
    if max((r.get('rank_lower_bound') or 0) for r in ledger['rows'])<25:raise ArithmeticError('declared positive-yield gate not met')
    censuses={f:cert.hashed(PARENT/f/'generic-census.json') for f in parent['families']}
    tables={f:{str(sign):{'path':str((PARENT/f/f'tables-{sign}.txt').relative_to(ROOT)),'sha256':cert.hashed(PARENT/f/f'tables-{sign}.txt')} for sign in (-1,1)} for f in parent['families']}
    p={'schema':'elliptic-curves.prospective-mw16-wide-selection.v1','sources':sources(),'scanner_binary_sha256':cert.hashed(BINARY),
       'families':parent['families'],'height':4096,'prime_bound':4093,'retained_per_family':128,'finalists_per_family':4,
       'generic_census_parent_protocol_hash':digest(parent),'generic_census_hashes':censuses,'trace_tables':tables,
       'score':parent['score'],'order':parent['order'],'population':'All20400078 signed primitive nonzero n/d per family with abs(n),d<=4096; five families,102000390 addresses. This contains the earlier H1024 box.',
       'mathematical_gate':'The completed compact-MW16 pilot certified a25-dimensional subgroup from16 generic sections. A balanced wider rational box tests candidate incidence on the same five compact families without choosing a family from its outcome. This does not predict ranks or remove the coordinate-visibility limitation.',
       'scanner_wall_seconds':180,'selection_worker_wall_seconds':400,'selection_worker_rss_bytes':1073741824,'maximum_selection_workers':2,
       'point_scope':'Separate fixed20-address protocol;43 charts, height100000,4seconds per chart,300seconds/1.5GiB per worker, at most4 concurrently. Exactly repeated initial addresses reuse immutable terminal measurements without refilling.',
       'target_free_boundary':parent['target_free_boundary'],'scope':'No new generic census or trace calculation: identical input bytes are reused and hashed. The new population and its point batch are separately frozen. No catalogue or public exceptional point is read during selection or execution.'}
    for f in parent['families']:
        destination=directory/f/'generic-census.json';destination.parent.mkdir(parents=True,exist_ok=True)
        if destination.exists():raise FileExistsError('census copy exists')
        shutil.copyfile(PARENT/f/'generic-census.json',destination)
        if cert.hashed(destination)!=censuses[f]:raise ArithmeticError('generic census copy changed')
    checkpoint(directory/'protocol.json',p);print('FROZEN BALANCED H4096 MW16 EXTENSION',flush=True)

def parse(text,sign):
    rows=[];summary=None
    for line in text.splitlines():
        v=line.split()
        if not v:continue
        if v[0]=='C':
            _,n,d,a,b,g,h=v;n=sign*int(n);d=int(d)
            rows.append({'numerator':n,'denominator':d,'parameter':str(F(n,d)),'score_units':int(F(a)*10**12),'good_primes':int(g)})
        elif v[0]=='S':summary=list(map(int,v[1:]))
    if summary is None:raise ArithmeticError('scanner did not finish')
    return rows,summary

def select(directory,family):
    p=bindings(directory);folder=directory/family
    if (folder/'population.json').exists():raise FileExistsError('wide population already retained')
    shards=[]
    for sign in (-1,1):
        source=p['trace_tables'][family][str(sign)];table=ROOT/source['path']
        if cert.hashed(table)!=source['sha256']:raise ArithmeticError('frozen trace table changed')
        result=capture([str(BINARY),str(table),'4096','4096','128','0','1'],limits=Limits(180,536870912),log_path=folder/f'scan-{sign}.log')
        rows,summary=parse(result.stdout,sign)
        shard={'rows':rows,'summary':summary,'supervision':result.supervision,'table_source':source,'protocol_hash':digest(p)}
        checkpoint(folder/f'scan-{sign}.json',shard);shards.append(shard)
    rows=[r for s in shards for r in s['rows']];rows.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
    if not all(rows[3]['score_units']>s['rows'][-1]['score_units'] for s in shards):raise ArithmeticError('finalist boundary requires tie-complete enumeration')
    count=sum(s['summary'][3] for s in shards)
    if count!=20400078:raise ArithmeticError('unexpected primitive population count')
    checkpoint(folder/'population.json',{'family':family,'protocol_hash':digest(p),'candidate_count':count,'retained_candidates':rows[:128],'finalists':rows[:4],
        'unused_H_band':'duplicate p5; unused and not validation','target_free':True})
    print('MW16 H4096 SELECTED',family,[(r['parameter'],r['score_units']/10**12) for r in rows[:4]],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','select']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family');a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else select(a.directory,a.family)
