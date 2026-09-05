#!/usr/bin/env sage-python
"""Finite height-16384 full-prime population, including the known rank-28 scale.

This is a separate disclosed expansion of the compact family experiment.
Retained smaller populations and their scores are never changed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
batch=SourceFileLoader('large_batch',str(CAS/'compact_r17_batch.sage')).load_module()
m=batch.m;full=batch.full
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
from certify_compact_r17_candidates import isomorphic,DATABASE
from elliptic_candidate_record import weierstrass_invariants
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-h16384-v1'


def prepare(directory):
    directory.mkdir(parents=True,exist_ok=True)
    parent=m.read(full.DIRECTORY/'protocol.json')
    binary=full.wide.DIRECTORY/'scanner'
    if m.hashed(binary)!=parent['scanner_binary_sha256']:raise ArithmeticError('scanner binary changed')
    tables={sign:full.DIRECTORY/f'tables-{sign}.txt' for sign in (-1,1)}
    for sign,path in tables.items():
        if m.hashed(path)!=parent['table_sha256'][str(sign)]:raise ArithmeticError('full-prime table changed')
    protocol={**m.read(batch.DIRECTORY/'protocol.json'),
        'schema':'elliptic-curves.compact-r17-h16384-protocol.v1','height':16384,
        'sources':m.provenance(),
        'selection_source':{str(p.relative_to(ROOT)):m.hashed(p) for p in
            (Path(__file__),CAS/'compact_r17_batch.sage',CAS/'compact_r17_fullscore.sage',CAS/'compact_r17_wide.sage',full.wide.SCANNER)},
        'parent_population_sha256':m.hashed(full.DIRECTORY/'population.json'),
        'table_sha256':parent['table_sha256'],'scanner_binary_sha256':m.hashed(binary),
        'seconds_per_shard':180,'max_parallel_scanners':2,'denominator_shards':4,
        'gate':'H4096 excludes the published rank-28 parameter of height 9529. Its completed fresh fibres reach rank 24; separately expand to H16384, scoring every prime through 4093 before selection. No statement that ranks 28 or 32 must occur.',
        'controls_after_population_freeze':['-2/377','-9529/5471'],
        'finalists':64,'first_survivors':1024}
    pp=directory/'protocol.json'
    if pp.exists() and m.read(pp)!=protocol:raise ArithmeticError('large protocol changed')
    checkpoint(pp,protocol)
    def shard(task):
        sign,k=task;path=directory/f'scan-{sign}-{k}.json'
        if path.exists():
            row=m.read(path)
            if row['protocol_hash']!=m.identity(protocol):raise ArithmeticError('retained shard changed')
            return row
        r=capture([str(binary),str(tables[sign]),'16384','16384','1024',str(k),'4'],
            limits=Limits(180,536870912),log_path=directory/f'scan-{sign}-{k}.log')
        rows,summary=full.wide.parse(r.stdout,sign)
        for row in rows:
            row['score_units']=row.pop('prefix_units');row.pop('extension_units');row.pop('extension_good')
        row={'protocol_hash':m.identity(protocol),'sign':sign,'shard':k,'rows':rows,'summary':summary,'supervision':r.supervision}
        checkpoint(path,row);print('LARGE SHARD',sign,k,summary[3],r.supervision['wall_seconds'],flush=True);return row
    with ThreadPoolExecutor(max_workers=2) as pool:
        shards=list(pool.map(shard,[(sign,k) for sign in (-1,1) for k in range(4)]))
    rows=[r for s in shards for r in s['rows']]
    if len(set(r['parameter'] for r in rows))!=len(rows):raise ArithmeticError('shard overlap')
    rows.sort(key=lambda r:(-r['score_units'],-r['prefix_good'],r['denominator'],r['numerator']))
    checkpoint(directory/'population.json',{'protocol_hash':m.identity(protocol),
        'candidate_count':sum(s['summary'][3] for s in shards),'retained_candidates':rows[:1024],
        'finalists':rows[:64],'public_points_or_record_equations_used_for_selection':False,
        'shard_hashes':{f'{s["sign"]}-{s["shard"]}':m.identity(s) for s in shards}})
    print('LARGE FROZEN',[(r['parameter'],r['score_units']/10**12) for r in rows[:16]],flush=True)


def roster(directory):
    protocol=m.read(directory/'protocol.json');pop=m.read(directory/'population.json')
    if pop['protocol_hash']!=m.identity(protocol):raise ArithmeticError('population changed')
    equations=[{'id':r['id'],'ainvs':r['ainvs']} for r in m.read(DATABASE)['curves']]
    byj={}
    for row in equations:
        inv=weierstrass_invariants(row['ainvs']);byj.setdefault(inv['c4']**3/inv['discriminant'],[]).append(row)
    prior={}
    for d in (m.DIRECTORY,full.wide.DIRECTORY,full.DIRECTORY,batch.DIRECTORY):
        for p in d.glob('*/result.json'):
            row=m.read(p)
            if row['status']=='COMPLETE':prior[row['parameter']]=p
    family=m.load_q12o5867_data(m.MODEL,m.SECTIONS);rows=[]
    for i,row in enumerate(pop['finalists']):
        t=m.F(row['parameter']);spec=m.evaluate_projective_specialization(family,t.numerator,t.denominator)
        inv=weierstrass_invariants(spec.model);j=inv['c4']**3/inv['discriminant']
        matches=[r['id'] for r in byj.get(j,[]) if isomorphic(spec.model,r['ainvs'])]
        old=prior.get(row['parameter'])
        rows.append({'index':i,'parameter':row['parameter'],'snapshot_matches':matches,
            'prior_complete_measurement':str(old.relative_to(ROOT)) if old else None,
            'prior_measurement_sha256':m.hashed(old) if old else None})
    checkpoint(directory/'roster.json',{'protocol_hash':m.identity(protocol),'database_sha256':m.hashed(DATABASE),'rows':rows})
    print('ROSTER',len(rows),'known',sum(bool(r['snapshot_matches']) for r in rows),'prior',sum(bool(r['prior_complete_measurement']) for r in rows),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','roster','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);a=p.parse_args()
    {'prepare':prepare,'roster':roster,'run':batch.run}[a.stage](a.directory)


if __name__=='__main__':main()
